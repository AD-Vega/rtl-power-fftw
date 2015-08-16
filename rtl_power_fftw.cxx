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
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <complex>

#include <fftw3.h>
#include <rtl-sdr.h>
#include <tclap/CmdLine.h>

static rtlsdr_dev_t *dev = nullptr;

// Get current date/time, format is "YYYY-MM-DD HH:mm:ss UTC"
std::string currentDateTime() {
  time_t now = time(0);
  char buf[80];
  strftime(buf, sizeof(buf), "%Y-%m-%d %X UTC", gmtime(&now));
  return buf;
}

using Buffer = std::vector<uint8_t>;
using complex = std::complex<double>;

class Datastore {
  public:
    int N;
    int buffers;
    int64_t repeats;
    int64_t repeats_done = 0;

    std::mutex status_mutex;
    // Access to the following objects must be protected by locking
    // status_mutex.
    std::deque<Buffer*> empty_buffers;
    std::deque<Buffer*> occupied_buffers;
    bool acquisition_finished = false;
    std::condition_variable status_change;
    std::vector<int> queue_histogram;

    complex *inbuf, *outbuf;
    fftw_plan plan;
    std::vector<double> pwr;

    Datastore(int N, int buf_length, int64_t repeats, int buffers);
    ~Datastore();

    // Delete these so we don't accidentally mess anything up by copying
    // pointers to fftw_malloc'd buffers.
    Datastore(const Datastore&) = delete;
    Datastore(Datastore&&) = delete;
    Datastore& operator=(const Datastore&) = delete;
    Datastore& operator=(Datastore&&) = delete;
};

Datastore::Datastore(int N_, int buf_length, int64_t repeats_, int buffers_) :
  N(N_), buffers(buffers_), repeats(repeats_),
  queue_histogram(buffers_+1, 0), pwr(N)
{
  for (int i = 0; i < buffers; i++)
    empty_buffers.push_back(new Buffer(buf_length));

  inbuf = (complex*)fftw_alloc_complex(N);
  outbuf = (complex*)fftw_alloc_complex(N);
  plan = fftw_plan_dft_1d(N, (fftw_complex*)inbuf, (fftw_complex*)outbuf,
			  FFTW_FORWARD, FFTW_MEASURE);
}

Datastore::~Datastore() {
  for (auto& buffer : empty_buffers)
    delete buffer;

  for (auto& buffer : occupied_buffers)
    delete buffer;

  fftw_destroy_plan(plan);
  fftw_free(inbuf);
  fftw_free(outbuf);
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

void fft(Datastore& data) {
  std::unique_lock<std::mutex>
    status_lock(data.status_mutex, std::defer_lock);
  int fft_pointer = 0;
  while (true) {
    // Wait until we have a bufferful of data
    status_lock.lock();
    while (data.occupied_buffers.empty() && !data.acquisition_finished)
      data.status_change.wait(status_lock);
    if (data.occupied_buffers.empty()) {
      // acquisition finished
      break;
    }
    Buffer& buffer(*data.occupied_buffers.front());
    data.occupied_buffers.pop_front();
    status_lock.unlock();
    //A neat new loop to avoid having to have data buffer aligned with fft buffer.
    unsigned int buffer_pointer = 0;
    while (buffer_pointer < buffer.size() && data.repeats_done < data.repeats ) {
      while (fft_pointer < data.N && buffer_pointer < buffer.size()) {
        //The magic aligment happens here: we have to change the phase of each next complex sample
        //by pi - this means that even numbered samples stay the same while odd numbered samples
        //get multiplied by -1 (thus rotated by pi in complex plane).
        //This gets us output spectrum shifted by half its size - just what we need to get the output right.
        const double multiplier = (fft_pointer % 2 == 0 ? 1 : -1);
	complex bfr_val(buffer[buffer_pointer], buffer[buffer_pointer+1]);
	data.inbuf[fft_pointer] = (bfr_val - 127.0) * multiplier;
	buffer_pointer += 2;
        fft_pointer++;
      }
      if (fft_pointer == data.N) {
        fftw_execute(data.plan);
        for (int i = 0; i < data.N; i++) {
	  data.pwr[i] += std::norm(data.outbuf[i]);
        }
        data.repeats_done++;
        fft_pointer = 0;
      }
    }

    status_lock.lock();
    data.empty_buffers.push_back(&buffer);
    data.status_change.notify_all();
    status_lock.unlock();
  }
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
  int cfreq = 89300000;
  int sample_rate = 2000000;
  int integration_time = 0;
  int integration_time_isSet = 0;
  int rtl_retval;
  int buffers = 5;
  int buf_length = 16384*100;
  int ppm_error = 0;
  bool endless = false;
  //It is senseless to waste a full buffer of data unless instructed to do so.
  int64_t repeats = buf_length/(2*N);
  try {
    TCLAP::CmdLine cmd("Obtain power spectrum from RTL device using FFTW library.", ' ', "0.1");
    TCLAP::ValueArg<int> arg_integration_time("t","time","Integration time in seconds (incompatible with -n).",false,integration_time,"seconds");
    cmd.add( arg_integration_time );
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
    TCLAP::ValueArg<int> arg_freq("f","freq","Center frequency of the receiver.",false,cfreq,"Hz");
    cmd.add( arg_freq );
    TCLAP::ValueArg<int> arg_index("d","device","RTL-SDR device index.",false,dev_index,"device index");
    cmd.add( arg_index );
    TCLAP::SwitchArg arg_continue("c","continue","Repeat the same measurement endlessly.", endless);
    cmd.add( arg_continue );
    TCLAP::ValueArg<int> arg_bins("b","bins","Number of bins in FFT spectrum (must be even number)",false,N,"bins in FFT spectrum");
    cmd.add( arg_bins );
    TCLAP::ValueArg<int> arg_buffers("B","buffers","Number of read buffers (don't touch unless running out of memory).",false,buffers,"buffers");
    cmd.add( arg_buffers );

    cmd.parse(argc, argv);

    try {
      // Ain't this C++11 f**** magic? Watch this:
      ensure_positive_arg<int>({&arg_bins, &arg_freq, &arg_rate, &arg_gain, &arg_integration_time, &arg_index, &arg_buffers, &arg_bufferlen});
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
    cfreq = arg_freq.getValue();
    sample_rate = arg_rate.getValue();
    buffers = arg_buffers.getValue();
    buf_length = arg_bufferlen.getValue();
    endless = arg_continue.getValue();
    // Due to USB specifics, buffer length for reading rtl_sdr device
    // must be a multiple of 16384. We have to keep it that way.
    // For performance reasons, the actual buffer length should be in the
    // MB range.
    if (buf_length % 16384 != 0) {
      buf_length = floor((double)buf_length/16384.0 + 0.5)*16384;
      std::cerr << "Buffer length should be multiple of 16384, changing to " << buf_length << "." << std::endl;
    }
    ppm_error = arg_ppm.getValue();
    if (arg_repeats.isSet())
      repeats = arg_repeats.getValue();
    else
      repeats = buf_length/(2*N);
    if (arg_integration_time.isSet()) {
      integration_time = arg_integration_time.getValue();
      integration_time_isSet = 1;
    }
    //Integration time
    if (arg_integration_time.isSet() + arg_repeats.isSet() > 1) {
      std::cerr << "Options -n and -t are mutually exclusive. Exiting." << std::endl;
      return 3;
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

  //Center frequency
  rtlsdr_set_center_freq(dev, (uint32_t)cfreq);
  int tuned_freq = rtlsdr_get_center_freq(dev);
  std::cerr << "Device tuned to: " << tuned_freq << " Hz" << std::endl;
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
  if (integration_time_isSet == 1)
    repeats = ceil((double)actual_samplerate * integration_time / N);

  int64_t readouts = ceil((2.0 * N * repeats) / buf_length);

  //Print info on capture time and associated specifics.
  std::cerr << "Number of bins: " << N << std::endl;
  std::cerr << "Total number of (complex) samples to collect: " << (int64_t)N*repeats << std::endl;
  std::cerr << "Number of averaged spectra: " << repeats << std::endl;
  std::cerr << "Number of device readouts: " << readouts << std::endl;
  std::cerr << "Expected time of measurements: " << readouts*0.5*(double)buf_length/actual_samplerate << " seconds" << std::endl;

  //Begin the work: prepare data buffers
  Datastore data(N, buf_length, repeats, buffers);

  //Read from device and do FFT
  do {
    std::fill(data.pwr.begin(), data.pwr.end(), 0);
    data.acquisition_finished = false;
    data.repeats_done = 0;

    std::thread t(&fft, std::ref(data));

    // Record the start-of-acquisition timestamp.
    std::string startAcqTimestamp = currentDateTime();
    std::cerr << "Acquisition started at " << startAcqTimestamp << std::endl;

    std::unique_lock<std::mutex>
      status_lock(data.status_mutex, std::defer_lock);
    int64_t count = 0;
    while (count <= readouts) {
      // Wait until a buffer is empty
      status_lock.lock();
      data.queue_histogram[data.empty_buffers.size()]++;
      while (data.empty_buffers.empty())
        data.status_change.wait(status_lock);

      Buffer& buffer(*data.empty_buffers.front());
      data.empty_buffers.pop_front();
      status_lock.unlock();

      rtl_retval = read_rtlsdr(buffer);

      if (rtl_retval) {
        fprintf(stderr, "Error: dropped samples.\n");
        // There is effectively no data in this buffer - consider it empty.
        status_lock.lock();
        data.empty_buffers.push_back(&buffer);
        status_lock.unlock();
        // No need to notify the worker thread in this case.
      }
      else {
        count++;
        status_lock.lock();
        data.occupied_buffers.push_back(&buffer);
        data.status_change.notify_all();
        status_lock.unlock();
      }
    }
    // Record the end-of-acquisition timestamp.
    std::string endAcqTimestamp = currentDateTime();
    std::cerr << "Acquisition done at " << endAcqTimestamp << std::endl;

    status_lock.lock();
    data.acquisition_finished = true;
    data.status_change.notify_all();
    status_lock.unlock();
    t.join();

    //Write out.
    std::cout << "# rtl-power-fftw output" << std::endl;
    std::cout << "# Acquisition start: " << startAcqTimestamp << std::endl;
    std::cout << "# Acquisition end: " << endAcqTimestamp << std::endl;
    std::cout << "#" << std::endl;
    std::cout << "# frequency [Hz] power spectral density [dB/Hz]" << std::endl;

    //Interpolate the central point, to cancel DC bias.
    data.pwr[data.N/2] = (data.pwr[data.N/2 - 1] + data.pwr[data.N/2+1]) / 2;

    for (int i = 0; i < N; i++) {
      std::cout << tuned_freq + (i-N/2.0) * ( (double)actual_samplerate / ((double)N ) ) << " "
                << 10*log10(data.pwr[i]/ repeats) << std::endl;
    }
    if (endless) {
      // Separate measurement sets with empty lines.
      std::cout << std::endl;
    }
    std::cout.flush();

    std::cerr << "Buffer queue histogram: ";
    for (auto size : data.queue_histogram)
      std::cerr << size << " ";
    std::cerr << std::endl;
  } while (endless);

  rtlsdr_close(dev);
  return 0;
}
