/*
* rtl_power_fftw, program for calculating power spectrum from rtl-sdr reciever.
* Copyright (C) 2015 Klemen Blokar <klemen.blokar@ad-vega.si>
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
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <limits>
#include <thread>
#include <mutex>
#include <string>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <tclap/CmdLine.h>
//#include <libusb.h>
#include <rtl-sdr.h>
//#include <rtl-sdr/convenience/convenience.h>

// Indices of real and imaginary parts of complex numbers; for convenience.
#define RE 0
#define IM 1

#define BUFFERS 5

static rtlsdr_dev_t *dev = NULL;

// Greatest common denominator
int gcd(int a, int b) {
  if (b == 0)
      return a;
  return gcd(b, a % b);
}

// Lowest common multiplier
int lcm(int a, int b){
  return a*b / gcd(a, b);
}

class Buffer {
public:
  std::vector<uint8_t> array;
  bool ready = false;
  int usage = 0;
  std::mutex mutex;

  Buffer(size_t length) : array(length) {}
  uint8_t& operator[] (size_t n) { return array[n]; }
};

class Datastore {
  public:
    int N;
    int batches;
    int repeats;
    int repeats_done;
    std::atomic<bool> acquisition_done;
    std::vector<Buffer*> buffers;
    fftw_complex *inbuf, *outbuf;
    fftw_plan plan;
    std::vector<double> pwr;

    Datastore(int N, int buf_length, int batches, int repeats);
    ~Datastore();

    // Finds the first buffer with the specified status (readiness) and
    // returns a pointer to it. The selected buffer is returned with its
    // mutex locked and the caller is responsible for unlocking it.
    // Returns nullptr if no suitable buffer could be found.
    Buffer* find_buffer(bool ready);

    // Delete these so we don't accidentally mess anything up by copying
    // pointers to fftw_malloc'd buffers.
    Datastore(const Datastore&) = delete;
    Datastore(Datastore&&) = delete;
    Datastore& operator=(const Datastore&) = delete;
    Datastore& operator=(Datastore&&) = delete;
};

Datastore::Datastore(int N_, int buf_length, int batches_, int repeats_) :
  N(N_), batches(batches_), repeats(repeats_), repeats_done(0),
  acquisition_done(false), pwr(N)
{
  buffers.reserve(BUFFERS);
  for (int i = 0; i < BUFFERS; i++)
    buffers.push_back(new Buffer(buf_length));

  inbuf = fftw_alloc_complex(N);
  outbuf = fftw_alloc_complex(N);
  plan = fftw_plan_dft_1d(N, inbuf, outbuf, FFTW_FORWARD, FFTW_MEASURE);
}

Datastore::~Datastore() {
  for (auto& buffer : buffers)
    delete buffer;

  fftw_destroy_plan(plan);
  fftw_free(inbuf);
  fftw_free(outbuf);
}

Buffer* Datastore::find_buffer(bool ready) {
  for (auto& buffer : buffers) {
    if (buffer->mutex.try_lock()) {
      if (buffer->ready == ready)
        return buffer;
      else
        buffer->mutex.unlock();
    }
  }
  return nullptr;
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
  std::cerr << "Available gains: ";
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
  rtlsdr_read_sync(dev, buffer.array.data(), buffer.array.size(), &n_read);
  if (n_read != (signed)buffer.array.size()) {
    //fprintf(stderr, "Error: dropped samples.\n");
    return 1;
  }
  return 0;
}

void fft(Datastore& data) {
  bool do_exit = false;
  while (true) {
    // Try to find a ready buffer
    Buffer* bufptr = data.find_buffer(true);
    if (bufptr) {
      Buffer& buffer = *bufptr;
      std::lock_guard<std::mutex> lock(buffer.mutex, std::adopt_lock);

      for (int batch = 1;
           batch <= data.batches && data.repeats_done < data.repeats;
           batch++)
      {
        //std::cerr << "Processing repeat "<< data.repeats_done <<" of " << data.repeats << "in batch " << batch << "."<< std::endl;
        for (int i = 2*data.N*(batch-1), j = 0;
             i < 2*data.N*batch - 1;
             i += 4, j += 4) {
          //The magic aligment happens here: we have to change the phase of each next complex sample
          //by pi - this means that even numbered samples stay the same while odd numbered samples
          //get multiplied by -1 (thus rotated by pi in complex plane).
          //This gets us output spectrum shifted by half its size - just what we need to get the output right.
          data.inbuf[j/2][RE] = (double) buffer[i] - 127;
          data.inbuf[j/2][IM] = (double) buffer[i + 1] - 127;
          data.inbuf[j/2 + 1][RE] = ((double) buffer[i + 2] - 127) * -1;
          data.inbuf[j/2 + 1][IM] = ((double) buffer[i + 3] - 127) * -1;
        }
        fftw_execute(data.plan);
        for (int i = 0; i < data.N; i++) {
          data.pwr[i] += data.outbuf[i][RE] * data.outbuf[i][RE] + data.outbuf[i][IM] * data.outbuf[i][IM];
        }
        data.repeats_done++;
      }
      //Interpolate the central point, to cancel DC bias.
      data.pwr[data.N/2] = (data.pwr[data.N/2 - 1] + data.pwr[data.N/2+1]) / 2;
      buffer.ready = false;
    }
    else {
      if (do_exit) {
        break;
      }
      //usleep(500);
      if (data.acquisition_done) {
        do_exit = true;
      }
    }
  }
}

int main(int argc, char **argv)
{
  int N = 512;
  int repeats = 1;
  int dev_index = 0;
  int gain = 372;
  int cfreq = 89300000;
  int sample_rate = 2000000;
  int integration_time = 0;
  int rtl_retval;

  try {
    TCLAP::CmdLine cmd("Obtain power spectrum from RTL device using FFTW library.", ' ', "0.1");
    TCLAP::ValueArg<int> arg_bins("b","bins","Number of bins in FFT spectrum (must be multiple of 256)",false,512,"bins in FFT spectrum");
    cmd.add( arg_bins );
    TCLAP::ValueArg<int> arg_freq("f","freq","Center frequency of the receiver.",false,89300000,"Hz");
    cmd.add( arg_freq );
    TCLAP::ValueArg<int> arg_rate("r","rate","Sample rate of the receiver.",false,2000000,"samples/s");
    cmd.add( arg_rate );
    TCLAP::ValueArg<int> arg_gain("g","gain","Receiver gain.",false, 372, "1/10th of dB");
    cmd.add( arg_gain );
    TCLAP::ValueArg<int> arg_repeats("n","repeats","Number of scans for averaging (incompatible with -t).",false,1,"repeats");
    cmd.add( arg_repeats );
    TCLAP::ValueArg<int> arg_integration_time("t","time","Integration time in seconds (incompatible with -n).",false,0,"seconds");
    cmd.add( arg_integration_time );
    TCLAP::ValueArg<int> arg_index("d","device","RTL-SDR device index.",false,0,"device index");
    cmd.add( arg_index );

    cmd.parse(argc, argv);

    dev_index = arg_index.getValue();
    N = arg_bins.getValue();
    if (arg_repeats.isSet())
        repeats = arg_repeats.getValue();
    if (arg_integration_time.isSet())
        integration_time = arg_integration_time.getValue();
    //Integration time
    if (arg_integration_time.isSet() + arg_repeats.isSet() > 1) {
      std::cerr << "Options -n and -t are mutually exclusive. Exiting." << std::endl;
      return 3;
    }
    else if (arg_integration_time.isSet()) {
      repeats = (int)( (sample_rate * integration_time) / (double)N + 0.5 );
    }
    gain = arg_gain.getValue();
    cfreq = arg_freq.getValue();
    sample_rate = arg_rate.getValue();
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
  int number_of_gains = rtlsdr_get_tuner_gains(dev, NULL);
  std::vector<int> gain_table(number_of_gains);
  rtlsdr_get_tuner_gains(dev, gain_table.data());
  print_gain_table(gain_table);
  gain = select_nearest_gain(gain, gain_table);
  std::cerr << "Selected nearest available gain: " << gain << std::endl;
  rtlsdr_set_tuner_gain_mode(dev, 1);
  rtlsdr_set_tuner_gain(dev, gain);

  //Center frequency
  rtlsdr_set_center_freq(dev, (uint32_t)cfreq);
  int tuned_freq = rtlsdr_get_center_freq(dev);
  std::cerr << "Device tuned to: " << tuned_freq << " Hz." << std::endl;
  usleep(5000);

  //Sample rate
  rtlsdr_set_sample_rate(dev, (uint32_t)sample_rate);
  int actual_samplerate = rtlsdr_get_sample_rate(dev);

  //Print info on capture time
  std::cerr << "Number of averaged spectra: " << repeats << "." << std::endl;
  std::cerr << "Expected time of measurements: " << N*repeats/sample_rate << " seconds." << std::endl;

  //Number of bins should be even, to allow us a neat trick to get fftw output properly aligned.
  //rtl_sdr seems to be only able to read data from USB dongle in chunks of 256 (complex) data points.
  if (N % 256 != 0) {
    N = ((floor(N/256.0))+1)*256;
    std::cerr << "Number of bins should be multiple of 256, changing to " << N << "." << std::endl;
  }
  std::cerr << "Number of bins: " << N << std::endl;
  std::cerr << "Total number of (complex) samples to collect: " << N*repeats << std::endl;

  // Due to USB specifics, buffer length for reading rtl_sdr device
  // must be a multiple of 16384. We have to keep it that way.
  // For performance reasons, the actual buffer length should be in the
  // MB range.
  // TODO: allow user control over the buffer length.
  int buf_length = lcm(2*N, 16384*100);
  int readouts = ceil(2.0 * N * repeats / buf_length);
  int batches = buf_length / (2*N);
  std::cerr << "Data collection will proceed in " << readouts <<" readouts, each consisting of " << batches << " batches." << std::endl;

  //Begin the work: prepare data buffers
  Datastore data(N, buf_length, batches, repeats);
  for (int i = 0; i < N; i++) {
    data.pwr[i] = 0;
  }
  //Read from device and do FFT
  std::thread t(&fft, std::ref(data));
  int count = 0;
  while (count <= readouts) {
    // Try to find an empty buffer
    Buffer* bufptr = data.find_buffer(false);
    if (bufptr) {
      Buffer& buffer = *bufptr;
      std::lock_guard<std::mutex> lock(buffer.mutex, std::adopt_lock);
      rtl_retval = read_rtlsdr(buffer);
      if (rtl_retval) {
        fprintf(stderr, "Error: dropped samples.\n");
      }
      else {
        count++;
        buffer.ready = true;
      }
    }
  }
  std::cerr << "Acquisition_done." << std::endl;
  data.acquisition_done = true;
  t.join();
  //Write out.
  for (int i = 0; i < N; i++) {
    std::cout << i << " "
              << tuned_freq + (i-N/2.0) * ( (N-1) / (double)N  * (double)actual_samplerate / (double)N ) << " "
              << 10*log10(data.pwr[i]/ repeats) << std::endl;
  }
  std::cerr << "Buffer usage: ";
  for (const auto& buffer : data.buffers)
    std::cerr << buffer->usage << ", ";
  std::cerr << std::endl;
  rtlsdr_close(dev);
  return 0;
}
