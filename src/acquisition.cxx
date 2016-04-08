/*
* rtl_power_fftw, program for calculating power spectrum from rtl-sdr reciever.
* Copyright (C) 2015 Klemen Blokar <klemen.blokar@ad-vega.si>
*                    Andrej Lajovic <andrej.lajovic@ad-vega.si>
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>

#include "acquisition.h"
#include "dispatcher.h"
#include "device.h"
#include "interrupts.h"
#include "output.h"


template <typename T>
std::vector<T> read_inputfile(std::istream* stream) {
  // Parse baseline input line by line.
  std::vector<T> values;
  std::string line;
  while (std::getline(*stream, line)) {
    std::istringstream lineStream(line);

    if ((lineStream >> std::ws).peek() == '#') {
      // Commented lines don't count.
      continue;
    }
    // The strategy is: we read as much doubles from the line as we can,
    // and use the last one. Accomodates one column, two columns (the first
    // one being, for example, frequency), three columns, N columns;
    // anything goes.
    T value;
    unsigned int valuesRead = 0;
    while (lineStream >> value)
      valuesRead++;

    // We are being very relaxed here. No doubles in the line? Skip it.
    // As long as we end up with the right number of values, we're game.
    if (valuesRead > 0)
      values.push_back(value);
  }
  return values;
}

AuxData::AuxData(const Params& params) {
  std::istream* stream;
  std::ifstream fs;
  // Window function and baseline correction
  // If both are read from stdin, we read it all in one go,
  // and see if the total number of values adds up to what we need.
  // If so, the first half is window data, and the other half is
  // baseline data. Window data are floats, as our samples are only
  // 8-bit anyway, but baseline is double, as we can average a lot of
  // spectra and gain considerable precision in this manner.
  if (params.window && params.baseline && params.window_file == "-" && params.baseline_file == "-") {
    stream = &std::cin;
    Diagnostics() << "Reading baseline and window function from stdin.\n";
    std::vector<double> values = read_inputfile<double>(stream);
    if ((int)values.size() == 2*params.N) {
      std::size_t const half_size = window_values.size() / 2;
      for (std::size_t i=0; i < values.size(); i++) {
        if ( i < half_size )
          window_values.push_back((float)values[i]);
        else
          baseline_values.push_back(values[i]);
      }
      Diagnostics() << "Succesfully read " << window_values.size() << " window function points.\n"
                    << "Succesfully read " << baseline_values.size() << " baseline points.\n";
    }
    else {
      throw RPFexception(
        "Error reading window function and baseline from stdin. Expected "
        + std::to_string(2 * params.N) + " values, found "
        + std::to_string(values.size()) + ".",
      ReturnValue::InvalidInput);
    }
  }
  // In other scenarios we can safely read window function and
  // baseline data separately, as there is certainly no clash
  // on stdin anymore. The reason why we don't read all the data
  // separately in all cases is the ease with which we test for unsuitable input
  // data (too many points) if we read all of it and see if there is
  // just enough data points.
  else {
    if (params.window) {
      if (params.window_file == "-") {
        Diagnostics() << "Reading window function from stdin.\n";
        stream = &std::cin;
      }
      else {
        Diagnostics() << "Reading window function from file " << params.baseline_file << "\n";
        fs.open(params.window_file);
        if (!fs.good()) {
          throw RPFexception(
            "Could not open " + params.window_file + ". Quitting.",
            ReturnValue::InvalidInput);
        }
        stream = &fs;
      }
      window_values = read_inputfile<float>(stream);
      // Check for suitability.
      if ((int)window_values.size() == params.N) {
        Diagnostics() << "Succesfully read " << window_values.size() << " window function points.\n";
      }
      else {
        throw RPFexception(
          "Error reading window function. Expected " + std::to_string(params.N)
          + " values, found " + std::to_string(window_values.size()) + ".",
          ReturnValue::InvalidInput);
      }
    }
    if (params.baseline) {
      if (params.baseline_file == "-") {
        Diagnostics() << "Reading baseline from stdin.\n";
        stream = &std::cin;
      }
      else {
        Diagnostics() << "Reading baseline from file " << params.baseline_file << "\n";
        fs.open(params.baseline_file);
        if (!fs.good()) {
          throw RPFexception(
            "Could not open " + params.baseline_file + ". Quitting.",
            ReturnValue::InvalidInput);
        }
        stream = &fs;
      }
      baseline_values = read_inputfile<double>(stream);
      // Check for suitability.
      if ((int)baseline_values.size() == params.N) {
        Diagnostics() << "Succesfully read " << baseline_values.size() << " baseline points.\n";
      }
      else {
        throw RPFexception(
          "Error reading baseline. Expected " + std::to_string(params.N)
          + " values, found " + std::to_string(baseline_values.size()) + ".",
          ReturnValue::InvalidInput);
      }
    }
  }
}

Plan::Plan(Params& params_, int actual_samplerate_) :
  actual_samplerate(actual_samplerate_), params(params_)
{
  // Calculate the number of repeats according to the true sample rate.
  if (params.integration_time_isSet)
    params.repeats = ceil(actual_samplerate * params.integration_time / params.N);

  // Adjust buffer size
  if (!params.buf_length_isSet) {
    int64_t base_buf_multiplier = ceil((2.0 * params.N * params.repeats) / base_buf);
    // If less than approximately 1.6 MB of data is needed, make the buffer the
    // smallest possible while still keeping the size to a multiple of
    // base_buf. Otherwise, set it to 100 * base_buf.
    // If you know what should fit your purposes well, feel free to use the
    // command line options to override this.
    if (base_buf_multiplier <= default_buf_multiplier) {
      params.buf_length = base_buf * ((base_buf_multiplier == 0 ) ? 1 : base_buf_multiplier);
    }
  }

  // Make a plan of frequency hopping.
  // We're stuffing a vector full of frequencies that we wish to eventually tune to.
  if (params.freq_hopping_isSet) {
    double min_overhang = actual_samplerate*params.min_overlap/100;
    int hops = ceil((double(params.stopfreq - params.startfreq) - min_overhang) / (double(actual_samplerate) - min_overhang));
    if (hops > 1) {
      int overhang = (hops*actual_samplerate - (params.stopfreq - params.startfreq)) / (hops - 1);
      freqs_to_tune.push_back(params.startfreq + actual_samplerate/2.0);
      //Mmmm, thirsty? waah-waaah...
      for (int hop = 1; hop < hops; hop++) {
        freqs_to_tune.push_back(freqs_to_tune.back() + actual_samplerate - overhang);
      }
    }
    else
      freqs_to_tune.push_back((params.startfreq + params.stopfreq)/2);
  }
  // If there is only one hop, no problem.
  else {
    freqs_to_tune.push_back(params.cfreq);
  }
}

void Plan::print() const {
  Diagnostics diag;
  diag
    << "Number of bins: " << params.N << "\n"
    << "Total number of (complex) samples to collect: " << (int64_t)params.N*params.repeats << "\n"
    << "Buffer length: " << params.buf_length << "\n"
    << "Number of averaged spectra: " << params.repeats << "\n"
    << "Estimated time of measurements: " << (double)params.N * params.repeats / actual_samplerate << " seconds\n";

    if (params.strict_time)
      diag << "Acquisition will unconditionally terminate after " << params.integration_time << " seconds.\n";
}


Acquisition::Acquisition(const Params& params_,
                         AuxData& aux_,
                         Rtlsdr& rtldev_,
                         Dispatcher& dispatcher_,
                         int actual_samplerate_,
                         int freq_) :
  pwr(params_.N), actual_samplerate(actual_samplerate_), params(params_), aux(aux_), rtldev(rtldev_), dispatcher(dispatcher_),
  freq(freq_), queue_histogram(params.buffers+1, 0)
{ }


void Acquisition::run() {
  // Set center frequency.
  // There have been accounts of hardware being stubborn and refusing to
  // tune to the desired frequency on random occasions despite being able
  // to tune to that same frequency at other times. Such hiccups seem to
  // be rare. We handle them by a naive and stupid, but seemingly effective
  // method of persuasion.
  const int max_tune_tries = 3;
  bool success = false;
  for (int tune_try = 0; !success && tune_try < max_tune_tries; tune_try++)
  {
    Diagnostics() << "Tuning to " << freq << " Hz (try " << tune_try + 1 << ")\n";
    try {
      rtldev.set_frequency(freq);
      tuned_freq = rtldev.frequency();
      if (tuned_freq != 0)
        success = true;
    }
    catch (RPFexception) {}
  }

  // Check if the frequency was actually successfully set.
  if (!success) {
    //Warning: librtlsdr does not tell you of all cases when tuner cannot lock PLL, despite clearly writing so to the stderr!
    //TODO: Fix librtlsdr.
    throw TuneError(freq);
  }

  std::fill(pwr.begin(), pwr.end(), 0);

  // Record the start-of-acquisition timestamp.
  startAcqTimestamp = currentDateTime();

  Diagnostics() << "Device tuned to: " << tuned_freq << " Hz\n"
                << "Acquisition started at " << startAcqTimestamp << "\n";

  // Calculate the stop time. This will only be effective if --strict-time was given.
  using steady_clock = std::chrono::steady_clock;
  steady_clock::time_point stopTime = steady_clock::now() + std::chrono::milliseconds(int64_t(params.integration_time*1000));

  int64_t dataTotal = 2 * params.N * params.repeats;
  int64_t dataRead = 0;

  while (dataRead < dataTotal) {
    // Obtain an empty container from the dispatcher.
    size_t queueSize;
    DataContainer container = dispatcher.emptyContainers.get(queueSize);
    queue_histogram[queueSize]++;
    container.acquisition = this;
    RawBuffer& rawBuffer = *container.data;

    // Figure out how much data to read.
    int64_t dataNeeded = dataTotal - dataRead;
    if (dataNeeded >= params.buf_length)
      // More than one bufferful of data needed. Leave the rest for later.
      dataNeeded = params.buf_length;
    else {
      // Less than one whole buffer needed. Round the number of (real)
      // samples upwards to the next multiple of base_buf.
      dataNeeded = base_buf * ceil((double)dataNeeded / base_buf);
      if (dataNeeded > params.buf_length) {
        // Nope, too much. We'll still have to do this in two readouts.
        dataNeeded = params.buf_length;
      }
    }
    // Resize the buffer to match the needed amount of data.
    rawBuffer.resize(dataNeeded);

    bool read_success = rtldev.read(rawBuffer);
    deviceReadouts++;

    if (!read_success) {
      Diagnostics(LogLevel::Error) << "Error: dropped samples.\n";
      // There is effectively no data in this container - consider it empty.
      // Push the container to the front of the queue because it already has
      // a buffer of the correct size and we'll just pop it again on next
      // iteration.
      dispatcher.emptyContainers.push_front(container);
    }
    else {
      successfulReadouts++;
      dataRead += dataNeeded;
      dispatcher.occupiedContainers.push_back(container);
    }

    if (params.strict_time && (steady_clock::now() >= stopTime))
      break;

    // See if we have been instructed to conclude this measurement immediately.
    if (interrupts && checkInterrupt(InterruptState::FinishNow))
        break;
  }

  // Record the end-of-acquisition timestamp.
  endAcqTimestamp = currentDateTime();
  Diagnostics() << "Acquisition done at " << endAcqTimestamp << "\n";

  // Push a sentinel container into the queue to mark the end of acquisition.
  dispatcher.occupiedContainers.push_back({this, nullptr});
}

void Acquisition::printSummary() const {
  Diagnostics diag;
  diag << "Actual number of (complex) samples collected: "
       << (int64_t)params.N * repeatsProcessed << "\n"
       << "Actual number of device readouts: " << deviceReadouts << "\n"
       << "Number of successful readouts: " << successfulReadouts << "\n"
       << "Actual number of averaged spectra: " << repeatsProcessed << "\n"
       << "Effective integration time: "
       << (double)params.N * repeatsProcessed / actual_samplerate << " seconds\n";

  diag << "Buffer queue histogram: ";
  for (auto size : queue_histogram)
    diag << size << " ";
  diag << "\n";
}

void Acquisition::markResultsReady() {
  // Finalize the result: interpolate the central point to cancel DC bias.
  pwr[params.N/2] = (pwr[params.N/2 - 1] + pwr[params.N/2+1]) / 2;

  std::lock_guard<std::mutex> guard(mutex);
  resultsReady = true;
  event.notify_one();
}

void Acquisition::waitForResultsReady() const {
  std::unique_lock<std::mutex> lock(mutex);
  event.wait(lock, [this]{ return resultsReady; });
}

// Get current date/time, format is "YYYY-MM-DD HH:mm:ss UTC"
std::string Acquisition::currentDateTime() {
  time_t now = std::time(0);
  char buf[80];
  std::strftime(buf, sizeof(buf), "%Y-%m-%d %X UTC", std::gmtime(&now));
  return buf;
}
