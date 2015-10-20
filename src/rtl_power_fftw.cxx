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
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <mutex>
#include <string>
#include <thread>
#include <ctime>
#include <fstream>

#include <rtl-sdr.h>

#include "datastore.h"
#include "params.h"

static rtlsdr_dev_t *dev = nullptr;

// Get current date/time, format is "YYYY-MM-DD HH:mm:ss UTC"
std::string currentDateTime() {
  time_t now = std::time(0);
  char buf[80];
  std::strftime(buf, sizeof(buf), "%Y-%m-%d %X UTC", std::gmtime(&now));
  return buf;
}


int select_nearest_gain(int gain, const std::vector<int>& gain_table) {
  int dif = std::numeric_limits<int>::max();
  int selected = 0;
  for (const auto& trial_gain : gain_table) {
    int temp = abs(trial_gain - gain);
    if ( temp < dif ) {
      dif = temp;
      selected = trial_gain;
    }
  }
  return selected;
}


void print_gain_table(const std::vector<int>& gain_table) {
  std::cerr << "Available gains (in 1/10th of dB): ";
  for (unsigned int i = 0; i < gain_table.size(); i++) {
    if (i != 0)
      std::cerr << ", ";
    std::cerr << gain_table[i];
  }
  std::cerr << std::endl;
}

int read_rtlsdr(Buffer& buffer) {
  int n_read;
  rtlsdr_reset_buffer(dev);
  rtlsdr_read_sync(dev, buffer.data(), buffer.size(), &n_read);
  if (n_read != (signed)buffer.size()) {
    //fprintf(stderr, "Error: dropped samples.\n");
    return 1;
  }
  return 0;
}


int main(int argc, char **argv)
{
  int rtl_retval;
  std::vector<double> baseline_values;
  std::list<int> freqs_to_tune;

  Params params;
  ReturnValue retval = params.parse(argc, argv);
  if (retval != ReturnValue::Success)
    return (int)retval;

  // Convenient shortened names.
  int& N = params.N;
  int& buf_length = params.buf_length;
  bool& baseline = params.baseline;
  int64_t& repeats = params.repeats;

  //Baseline correction
  if (baseline) {
    std::istream* stream;
    std::ifstream fs;

    if (params.baseline_file == "-") {
      std::cerr << "Reading baseline from stdin." << std::endl;
      stream = &std::cin;
    }
    else {
      std::cerr << "Reading baseline from file " << params.baseline_file << std::endl;
      fs.open(params.baseline_file);
      stream = &fs;
    }

    // Parse baseline input line by line.
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
      double value;
      unsigned int valuesRead = 0;
      while (lineStream >> value)
        valuesRead++;

      // We are being very relaxed here. No doubles in the line? Skip it.
      // As long as we end up with the right number of values, we're game.
      if (valuesRead > 0)
        baseline_values.push_back(value);
    }

    // Check for suitability.
    if ((int)baseline_values.size() == N) {
      std::cerr << "Succesfully read " << baseline_values.size() << " baseline points." << std::endl;
    }
    else {
      std::cerr << "Error reading baseline. Expected " << N << " samples, found "
                << baseline_values.size() << "." << std::endl;
      std::cerr << "Ignoring baseline data." << std::endl;
      baseline = false;
    }
  }

  //Sanity checks
  //RTLSDR Device
  int num_of_rtls = rtlsdr_get_device_count();
  if (num_of_rtls == 0) {
    std::cerr << "Error: no RTL-SDR compatible devices found. Exiting." << std::endl;
    return (int)ReturnValue::NoDeviceFound;
  }
  if (params.dev_index >= num_of_rtls) {
    std::cerr << "Error: invalid device number. Only "<< num_of_rtls << " devices available. Exiting." << std::endl;
    return (int)ReturnValue::InvalidDeviceIndex;
  }
  rtl_retval = rtlsdr_open(&dev, (uint32_t)params.dev_index);
  if (rtl_retval < 0 ) {
    std::cerr << "Could not open rtl_sdr device " << params.dev_index << "." << std::endl;
    return rtl_retval;
  }

  //Available gains
  int number_of_gains = rtlsdr_get_tuner_gains(dev, nullptr);
  std::vector<int> gain_table(number_of_gains);
  rtlsdr_get_tuner_gains(dev, gain_table.data());
  print_gain_table(gain_table);
  int gain = select_nearest_gain(params.gain, gain_table);
  std::cerr << "Selected nearest available gain: " << gain
            << " (" << 0.1*gain << " dB)" << std::endl;
  rtlsdr_set_tuner_gain_mode(dev, 1);
  rtlsdr_set_tuner_gain(dev, gain);

  // Temporarily set the frequency to params.cfreq, just so that the device does not
  // complain upon setting the sample rate.
  rtl_retval = rtlsdr_set_center_freq(dev, (uint32_t)params.cfreq);
  std::this_thread::sleep_for(std::chrono::milliseconds(5));

  //Frequency correction
  if (params.ppm_error != 0) {
    rtl_retval = rtlsdr_set_freq_correction(dev, params.ppm_error);
    if (rtl_retval < 0)
      std::cerr << "Unable to set PPM error in rtl_sdr device." << std::endl;
    else
      std::cerr << "PPM error set to: "<< params.ppm_error << std::endl;
  }

  //Sample rate
  rtlsdr_set_sample_rate(dev, (uint32_t)params.sample_rate);
  int actual_samplerate = rtlsdr_get_sample_rate(dev);
  std::cerr << "Actual sample rate: " << actual_samplerate << " Hz" << std::endl;
  //It is only fair to calculate repeats with actual samplerate, not our wishes.
  if (params.integration_time_isSet)
    repeats = ceil((double)actual_samplerate * params.integration_time / N);

  //Frequency hopping
  //We're stuffing a vector full of frequencies that we wish to eventually tune to.
  if (params.freq_hopping_isSet) {
    int hops = ceil((params.stopfreq - params.startfreq) / actual_samplerate);
    int overhang = (hops*actual_samplerate - (params.stopfreq - params.startfreq)) / (hops + 1);
    freqs_to_tune.push_back(params.startfreq + actual_samplerate/2.0 - overhang);
    //Mmmm, thirsty? waah-waaah...
    for (int hop = 1; hop < hops; hop++) {
      freqs_to_tune.push_back( freqs_to_tune.back() + actual_samplerate - overhang);
    }
  }
  // If there is only one hop, no problem.
  else {
    freqs_to_tune.push_back(params.cfreq);
  }
  //Adjust buffer length in case of small sample batches.
  if (!params.buf_length_isSet) {
    int64_t base_buf_multiplier = ceil((2.0 * N * repeats) / base_buf);
    // Optimisation works like this: if we need to sample less than ~1.6MB, 
    // make the buffer the smallest possible optimized size that fits all
    // the data.
    // If it is longer, but not long enough to make it irrelevant how full
    // the last buffer is, try to optimize buffer size so that it should fit better,
    // but without overcomplicating the estimation.
    // If you know what should fit your purposes well, feel free to override this
    // simply by assigning buffer length yourself.
    if (base_buf_multiplier <= default_buf_multiplier) {
      buf_length = base_buf * ((base_buf_multiplier == 0 ) ? 1 : base_buf_multiplier);
    }
    else if (base_buf_multiplier <= default_buf_multiplier*default_buf_multiplier) {
      buf_length = base_buf * ceil(sqrt(base_buf_multiplier));
    }
    else {
      // Keep default length.
    }
  }

  int64_t readouts = ceil((2.0 * N * repeats) / buf_length);

  //Print info on capture time and associated specifics.
  std::cerr << "Number of bins: " << N << std::endl;
  std::cerr << "Total number of (complex) samples to collect: " << (int64_t)N*repeats << std::endl;
  std::cerr << "Buffer length: " << buf_length << std::endl;
  std::cerr << "Number of device readouts: " << readouts << std::endl;
  std::cerr << "Number of averaged spectra: " << repeats << std::endl;
  std::cerr << "Estimated time of measurements: " << readouts*0.5*(double)buf_length/actual_samplerate << " seconds" << std::endl;
  if (params.strict_time)
    std::cerr << "Acquisition will unconditionally terminate after " << params.integration_time << " seconds." << std::endl;

  //Begin the work: prepare data buffers
  Datastore data(N, buf_length, repeats, params.buffers);

  //Read from device and do FFT
  do {
    for (auto iter = freqs_to_tune.begin(); iter != freqs_to_tune.end();) {
      // Set center frequency.
      rtl_retval = rtlsdr_set_center_freq(dev, (uint32_t)*iter);
      int tuned_freq = rtlsdr_get_center_freq(dev);
      // This sleeping is inherited from other code. There have been hints of strange
      // behaviour if it was commented out, so we left it in. If you actually know
      // why this would be necessary (or, to the contrary, that it is complete
      // bullshit), you are most welcome to explain it here!
      std::this_thread::sleep_for(std::chrono::milliseconds(5));

      // Check if the frequency was actually successfully set.
      if ( rtl_retval < 0 || tuned_freq == 0 ) {
        //Warning: librtlsdr does not tell you of all cases when tuner cannot lock PLL, despite clearly writing so to the stderr!
        //TODO: Fix librtlsdr.
        std::cerr << "Unable to tune to " << *iter << ". Dropping from frequency list." << std::endl;
        iter = freqs_to_tune.erase(iter);
        continue;
      }
      else {
        ++iter;
      }

      std::cerr << "Device tuned to: " << tuned_freq << " Hz" << std::endl;
      std::fill(data.pwr.begin(), data.pwr.end(), 0);
      data.acquisition_finished = false;
      data.repeats_done = 0;

      std::thread t(&Datastore::fftThread, std::ref(data));

      // Record the start-of-acquisition timestamp.
      std::string startAcqTimestamp = currentDateTime();
      std::cerr << "Acquisition started at " << startAcqTimestamp << std::endl;

      // Calculate the stop time. This will only be effective if --strict-time was given.
      using steady_clock = std::chrono::steady_clock;
      steady_clock::time_point stopTime = steady_clock::now() + std::chrono::seconds(params.integration_time);

      std::unique_lock<std::mutex>
        status_lock(data.status_mutex, std::defer_lock);
      int64_t deviceReadouts = 0;
      int64_t successfulReadouts = 0;

      while (successfulReadouts < readouts) {
        // Wait until a buffer is empty
        status_lock.lock();
        data.queue_histogram[data.empty_buffers.size()]++;
        while (data.empty_buffers.empty())
          data.status_change.wait(status_lock);

        Buffer& buffer(*data.empty_buffers.front());
        data.empty_buffers.pop_front();
        status_lock.unlock();

        rtl_retval = read_rtlsdr(buffer);
        deviceReadouts++;

        if (rtl_retval) {
          fprintf(stderr, "Error: dropped samples.\n");
          // There is effectively no data in this buffer - consider it empty.
          status_lock.lock();
          data.empty_buffers.push_back(&buffer);
          status_lock.unlock();
          // No need to notify the worker thread in this case.
        }
        else {
          successfulReadouts++;
          status_lock.lock();
          data.occupied_buffers.push_back(&buffer);
          data.status_change.notify_all();
          status_lock.unlock();
        }

        if (params.strict_time && (steady_clock::now() >= stopTime))
          break;
      }

      // Record the end-of-acquisition timestamp.
      std::string endAcqTimestamp = currentDateTime();
      std::cerr << "Acquisition done at " << endAcqTimestamp << std::endl;

      status_lock.lock();
      data.acquisition_finished = true;
      data.status_change.notify_all();
      status_lock.unlock();
      t.join();

      // Print a summary.
      std::cerr << "Actual number of (complex) samples collected: "
        << (int64_t)N * data.repeats_done << std::endl;
      std::cerr << "Actual number of device readouts: " << deviceReadouts << std::endl;
      std::cerr << "Number of successful readouts: " << successfulReadouts << std::endl;
      std::cerr << "Actual number of averaged spectra: " << data.repeats_done << std::endl;
      std::cerr << "Effective integration time: " <<
        (double)N * data.repeats_done / actual_samplerate << " seconds" << std::endl;

      //Write out data.
      std::cout << "# rtl-power-fftw output" << std::endl;
      std::cout << "# Acquisition start: " << startAcqTimestamp << std::endl;
      std::cout << "# Acquisition end: " << endAcqTimestamp << std::endl;
      std::cout << "#" << std::endl;
      std::cout << "# frequency [Hz] power spectral density [dB/Hz]" << std::endl;

      //Interpolate the central point, to cancel DC bias.
      data.pwr[data.N/2] = (data.pwr[data.N/2 - 1] + data.pwr[data.N/2+1]) / 2;

      // Calculate the precision needed for displaying the frequency.
      const int extraDigitsFreq = 2;
      const int significantPlacesFreq =
        ceil(floor(log10(params.cfreq)) - log10(actual_samplerate/N) + 1 + extraDigitsFreq);
      const int significantPlacesPwr = 6;

      for (int i = 0; i < N; i++) {
        double freq = tuned_freq + (i - N/2.0) * actual_samplerate / N;
        double pwrdb = 10*log10(data.pwr[i] / data.repeats_done / N / actual_samplerate) - (baseline ? baseline_values[i] : 0);
        std::cout << std::setprecision(significantPlacesFreq)
                  << freq
                  << " "
                  << std::setprecision(significantPlacesPwr)
                  << pwrdb
                  << std::endl;
      }
      if (params.endless || params.freq_hopping_isSet) {
        // Separate different spectra with empty lines.
        std::cout << std::endl;
      }
      std::cout.flush();

      std::cerr << "Buffer queue histogram: ";
      for (auto size : data.queue_histogram)
        std::cerr << size << " ";
      std::cerr << std::endl;
    }
    if (params.endless) {
      // Separate measurement sets with two empty lines.
      std::cout << std::endl;
    }
  } while (params.endless);

  rtlsdr_close(dev);
  return 0;
}
