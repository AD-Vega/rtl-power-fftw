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
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <ctime>
#include <fstream>

#include <rtl-sdr.h>
#include <tclap/CmdLine.h>

#include "datastore.h"

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

int parse_frequency(std::string s) {
  std::istringstream ss(s);
  double f;
  std::string multiplier;
  ss >> f >> multiplier;
  if (multiplier == "k")
    f *= 1e3;
  else if (multiplier == "M")
    f *= 1e6;
  else if (multiplier == "G")
    f *= 1e9;
  else if (multiplier != "")
    return -1;
  return (int)f;
}

class NegativeArgException {
public:
  NegativeArgException(std::string msg_) : msg(msg_) {}
  std::string what() const { return msg; }
private:
  std::string msg;
};

template <typename T>
void ensure_positive_arg(std::list<TCLAP::ValueArg<T>*> list) {
  for (auto arg : list) {
    if (arg->isSet() && arg->getValue() < 0) {
      std::ostringstream message;
      message << "Argument to '" << arg->getName() << "' must be a positive number.";
      throw NegativeArgException(message.str());
    }
  }
}

int main(int argc, char **argv)
{
  int N = 512;
  int dev_index = 0;
  int gain = 372;
  int cfreq = 1420405752;
  int startfreq = 0; 
  int stopfreq = 0;
  int sample_rate = 2000000;
  int integration_time = 0;
  bool integration_time_isSet = false;
  int rtl_retval;
  int buffers = 5;
  const int base_buf = 16384;
  const int default_buf_multiplier = 100;
  int buf_length = base_buf * default_buf_multiplier;
  bool buf_length_isSet = false;
  int ppm_error = 0;
  bool endless = false;
  bool strict_time = false;
  bool baseline = false;
  bool freq_hopping_isSet = false;
  //bool clipped_output_isSet = false;
  std::vector<double> baseline_values;
  std::list<int> freqs_to_tune;
  //It is senseless to waste a full buffer of data unless instructed to do so.
  int64_t repeats = buf_length/(2*N);
  try {
    TCLAP::CmdLine cmd("Obtain power spectrum from RTL device using FFTW library.", ' ', "0.1");
    TCLAP::ValueArg<int> arg_buffers("","buffers","Number of read buffers (don't touch unless running out of memory).",false,buffers,"buffers");
    cmd.add( arg_buffers );
    TCLAP::ValueArg<int> arg_integration_time("t","time","Integration time in seconds (incompatible with -n).",false,integration_time,"seconds");
    cmd.add( arg_integration_time );
    TCLAP::SwitchArg arg_strict_time("T","strict-time","End measurement when the time set with --time option is up, regardless of gathered samples.",strict_time);
    cmd.add( arg_strict_time );
    TCLAP::ValueArg<int> arg_bufferlen("s","buffer-size","Size of read buffers (leave it unless you know what you are doing).", false, buf_length, "bytes");
    cmd.add( arg_bufferlen );
    TCLAP::ValueArg<int> arg_rate("r","rate","Sample rate of the receiver.",false,sample_rate,"samples/s");
    cmd.add( arg_rate );
    TCLAP::ValueArg<int> arg_ppm("p","ppm","Set custom ppm error in RTL-SDR device.", false, ppm_error, "ppm");
    cmd.add( arg_ppm );
    TCLAP::ValueArg<int64_t> arg_repeats("n","repeats","Number of scans for averaging (incompatible with -t).",false,repeats,"repeats");
    cmd.add( arg_repeats );
    TCLAP::ValueArg<int> arg_gain("g","gain","Receiver gain.",false, gain, "1/10th of dB");
    cmd.add( arg_gain );
    TCLAP::ValueArg<std::string> arg_freq("f","freq","Center frequency of the receiver or frequency range to scan.",false,"","Hz | Hz:Hz");
    cmd.add( arg_freq );
    TCLAP::ValueArg<int> arg_index("d","device","RTL-SDR device index.",false,dev_index,"device index");
    cmd.add( arg_index );
    TCLAP::SwitchArg arg_continue("c","continue","Repeat the same measurement endlessly.", endless);
    cmd.add( arg_continue );
    //Possible future option
    /*
    TCLAP::SwitchArg arg_clipped("C","clipped","Clip output of frequency range: include only points within bounds, in ascending order.", clipped_output_isSet);
    cmd.add( arg_clipped );
    */
    TCLAP::ValueArg<int> arg_bins("b","bins","Number of bins in FFT spectrum (must be even number)",false,N,"bins in FFT spectrum");
    cmd.add( arg_bins );
    TCLAP::ValueArg<std::string> arg_baseline("B","baseline","Subtract baseline, read baseline data from file or stdin.",false,"","file|-");
    cmd.add( arg_baseline );

    cmd.parse(argc, argv);

    try {
      // Ain't this C++11 f**** magic? Watch this:
      ensure_positive_arg<int>({&arg_bins, &arg_rate, &arg_gain, &arg_integration_time, &arg_index, &arg_buffers, &arg_bufferlen});
      ensure_positive_arg<int64_t>({&arg_repeats});
    }
    catch (NegativeArgException& e) {
      std::cerr << e.what() << std::endl;
      return 3;
    }

    dev_index = arg_index.getValue();
    N = arg_bins.getValue();
    //Number of bins should be even, to allow us a neat trick to get fftw output properly aligned.
    if (N % 2 != 0) {
      N++;
      std::cerr << "Number of bins should be even, changing to " << N << "." << std::endl;
    }
    gain = arg_gain.getValue();
    sample_rate = arg_rate.getValue();
    buffers = arg_buffers.getValue();
    buf_length = arg_bufferlen.getValue();
    endless = arg_continue.getValue();
    strict_time = arg_strict_time.getValue();
    //clipped_output_isSet = arg_clipped.getValue();

    // Due to USB specifics, buffer length for reading rtl_sdr device
    // must be a multiple of 16384. We have to keep it that way.
    // For performance reasons, the actual buffer length should be in the
    // MB range.
    if (buf_length % base_buf != 0) {
      buf_length = floor((double)buf_length/base_buf + 0.5) * base_buf;
      std::cerr << "Buffer length should be multiple of " << base_buf
                << ", changing to " << buf_length << "." << std::endl;
    }
    ppm_error = arg_ppm.getValue();
    if (arg_freq.isSet()) {
      std::string a_freq = arg_freq.getValue();
      std::size_t colon_position = a_freq.find(":");
      std::istringstream opt(a_freq);
      if (colon_position != std::string::npos) {
        std::string startFreqString, stopFreqString;
        if (getline(opt, startFreqString, ':') && getline(opt, stopFreqString)) {
          startfreq = parse_frequency(startFreqString);
          stopfreq = parse_frequency(stopFreqString);
          if (startfreq < 0 || stopfreq < 0 || stopfreq < startfreq) {
            std::cerr << "Invalid frequency range given to --freq: " 
                    << startfreq << ":" << stopfreq << ". "
                    << "Expecting positive numbers in ascending order, allowing the k,M,G multipliers. Exiting."
                    << std::endl;
            return 3;
          }
          else {
            freq_hopping_isSet = true;
            cfreq = (startfreq + stopfreq)/2;
          }
        }
        else {
          std::cerr << "Could not parse frequency range given to --freq: " 
                    << opt 
                    << ". Expecting form startfreq:stopfreq. Exiting." 
                    << std::endl;
          return 3;
        }
      }
      else {
        cfreq = parse_frequency(a_freq);
        if (cfreq < 0) {
          std::cerr << "Invalid frequency given to --freq: " 
                    << cfreq
                    << ". Expecting a positive number, allowing the k,M,G multipliers. Exiting."
                    << std::endl;
          return 3;
        }
      }
    }
    if (arg_repeats.isSet())
      repeats = arg_repeats.getValue();
    else
      repeats = buf_length/(2*N);
    if (arg_integration_time.isSet()) {
      integration_time = arg_integration_time.getValue();
      integration_time_isSet = true;
    }
    //Integration time
    if (arg_integration_time.isSet() + arg_repeats.isSet() > 1) {
      std::cerr << "Options -n and -t are mutually exclusive. Exiting." << std::endl;
      return 3;
    }
    if (arg_strict_time.isSet() && !arg_integration_time.isSet()) {
      std::cerr << "Warning: option --strict-time has no effect without --time." << std::endl;
      strict_time = false;
    }
    //Optimally adjust buffer length for small sample sizes only if buffer length is not user defined.
    if (arg_bufferlen.isSet()) {
      buf_length_isSet = true;
    }
    //Baseline correction
    if (arg_baseline.isSet()) {
      const std::string& fileName = arg_baseline.getValue();
      std::istream* stream;
      std::ifstream fs;

      if (fileName == "-") {
        std::cerr << "Reading baseline from stdin." << std::endl;
        stream = &std::cin;
      }
      else {
        std::cerr << "Reading baseline from file " << fileName << std::endl;
        fs.open(fileName);
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
        baseline = true;
      }
      else {
        std::cerr << "Error reading baseline. Expected " << N << " samples, found "
                  << baseline_values.size() << "." << std::endl;
        std::cerr << "Ignoring baseline data." << std::endl;
        baseline = false;
      }
    }
  }
  catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 4;
  }

  //Sanity checks
  //RTLSDR Device
  int num_of_rtls = rtlsdr_get_device_count();
  if (num_of_rtls == 0) {
    std::cerr << "Error: no RTL-SDR compatible devices found. Exiting." << std::endl;
    return 1;
  }
  if ( dev_index >= num_of_rtls) {
    std::cerr << "Error: invalid device number. Only "<< num_of_rtls << " devices available. Exiting." << std::endl;
    return 2;
  }
  rtl_retval = rtlsdr_open(&dev, (uint32_t)dev_index);
  if (rtl_retval < 0 ) {
    std::cerr << "Could not open rtl_sdr device " << dev_index << "." << std::endl;
    return rtl_retval;
  }

  //Available gains
  int number_of_gains = rtlsdr_get_tuner_gains(dev, nullptr);
  std::vector<int> gain_table(number_of_gains);
  rtlsdr_get_tuner_gains(dev, gain_table.data());
  print_gain_table(gain_table);
  gain = select_nearest_gain(gain, gain_table);
  std::cerr << "Selected nearest available gain: " << gain
            << " (" << 0.1*gain << " dB)" << std::endl;
  rtlsdr_set_tuner_gain_mode(dev, 1);
  rtlsdr_set_tuner_gain(dev, gain);

  // Temporarily set the frequency to cfreq, just so that the device does not
  // complain upon setting the sample rate.
  rtl_retval = rtlsdr_set_center_freq(dev, (uint32_t)cfreq);
  std::this_thread::sleep_for(std::chrono::milliseconds(5));

  //Frequency correction
  if (ppm_error != 0) {
    rtl_retval = rtlsdr_set_freq_correction(dev, ppm_error);
    if (rtl_retval < 0)
      std::cerr << "Unable to set PPM error in rtl_sdr device." << std::endl;
    else
      std::cerr << "PPM error set to: "<< ppm_error << std::endl;
  }

  //Sample rate
  rtlsdr_set_sample_rate(dev, (uint32_t)sample_rate);
  int actual_samplerate = rtlsdr_get_sample_rate(dev);
  std::cerr << "Actual sample rate: " << actual_samplerate << " Hz" << std::endl;
  //It is only fair to calculate repeats with actual samplerate, not our wishes.
  if (integration_time_isSet)
    repeats = ceil((double)actual_samplerate * integration_time / N);

  //Frequency hopping
  //We're stuffing a vector full of frequencies that we wish to eventually tune to.
  if (freq_hopping_isSet) {
    int hops = ceil((stopfreq - startfreq) / actual_samplerate);
    int overhang = (hops*actual_samplerate - (stopfreq - startfreq)) / (hops + 1);
    freqs_to_tune.push_back(startfreq + actual_samplerate/2.0 - overhang);
    //Mmmm, thirsty? waah-waaah...
    for (int hop = 1; hop < hops; hop++) {
      freqs_to_tune.push_back( freqs_to_tune.back() + actual_samplerate - overhang);
    }
  }
  // If there is only one hop, no problem.
  else {
    freqs_to_tune.push_back(cfreq);
  }
  //Adjust buffer length in case of small sample batches.
  if (!buf_length_isSet) {
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
  if (strict_time)
    std::cerr << "Acquisition will unconditionally terminate after " << integration_time << " seconds." << std::endl;

  //Begin the work: prepare data buffers
  Datastore data(N, buf_length, repeats, buffers);

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
      steady_clock::time_point stopTime = steady_clock::now() + std::chrono::seconds(integration_time);

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

        if (strict_time && (steady_clock::now() >= stopTime))
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
        ceil(floor(log10(cfreq)) - log10(actual_samplerate/N) + 1 + extraDigitsFreq);
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
      if (endless || freq_hopping_isSet) {
        // Separate different spectra with empty lines.
        std::cout << std::endl;
      }
      std::cout.flush();

      std::cerr << "Buffer queue histogram: ";
      for (auto size : data.queue_histogram)
        std::cerr << size << " ";
      std::cerr << std::endl;
    }
    if (endless) {
      // Separate measurement sets with two empty lines.
      std::cout << std::endl;
    }
  } while (endless);

  rtlsdr_close(dev);
  return 0;
}
