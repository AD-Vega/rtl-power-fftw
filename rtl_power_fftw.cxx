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
#include <tclap/CmdLine.h>
//#include <libusb.h>
#include <rtl-sdr.h>
//#include <rtl-sdr/convenience/convenience.h>
#define RE 0
#define IM 1
#define BUFFERS 5

static rtlsdr_dev_t *dev = NULL;

int gcd(int a, int b) {
	if( b==0 ) return a;
	return gcd(b, a%b);
}

int lcm(int a, int b){
	return a*b / gcd(a,b);
}

class Datastore {
  public:
    int N;
    int buf_length;
    int batches;
    int repeats;
    int repeats_done;
    int acquisition_done;
    int buf_status[BUFFERS];
    uint8_t *buffer_arr[BUFFERS];
    std::mutex buf_mutex[BUFFERS];
    //uint8_t *buf8;
    fftw_complex *inbuf, *outbuf;
    fftw_plan plan;
    double *pwr;
    //std::mutex mu;
    Datastore(int a, int b, int c, int d);
};

Datastore::Datastore(int a, int b, int c, int d) {
  N = a;
  buf_length = b;
  batches = c;
  repeats = d;
  repeats_done = 0;
  acquisition_done = 0;
  for (int i = 0; i < BUFFERS; i++) {
    buffer_arr[i] = (uint8_t *) malloc (buf_length * sizeof(uint8_t));
    buf_status[i] = 0;
  }
  //buf8 = (uint8_t *) malloc (buf_length * sizeof(uint8_t));
  inbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  outbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  pwr = (double *) malloc (N*sizeof(double));
  plan = fftw_plan_dft_1d(N, inbuf, outbuf, FFTW_FORWARD, FFTW_MEASURE);
}

int select_nearest_gain(int gain, int size, int *gain_table) {
  int dif = std::numeric_limits<int>::max();
  int selected = 0;
  int temp = 0;
  for (int i=0; i < size; i++) {
    temp = abs(gain_table[i] - gain);
    if ( temp < dif ) {
      dif = temp;
      selected = gain_table[i];
    }
  }
  return selected;
}

void print_gain_table(int size, int *gain_table) {
  std::cerr << "Available gains: ";
  for (int i=0; i < size; i++) {
    if (i != 0) std::cerr << ", ";
    std::cerr << gain_table[i];
  }
  std::cerr << std::endl;
}

int read_rtlsdr(Datastore& data, int buf) {
  int n_read;
  rtlsdr_reset_buffer(dev);
  rtlsdr_read_sync(dev, data.buffer_arr[buf], data.buf_length, &n_read);
  if (n_read != data.buf_length) {
    //fprintf(stderr, "Error: dropped samples.\n");
    return 1;
  }
  return 0;
}

void fft(Datastore& data) {
  int canary = 0;
  while (true) {
    int k = 0;
    int buf_available = 0;
    while (k < BUFFERS) {
      if (data.buf_mutex[k].try_lock()) { 
	if (data.buf_status[k] == 1) {
	  buf_available = 1;
	  break;
	}
	else {
	  data.buf_mutex[k].unlock();
	}
      }
      k++;
    }
    if (buf_available == 1) {
      int i;
      int batch = 1;
      while (batch <= data.batches && data.repeats_done < data.repeats) {
	//std::cerr << "Processing repeat "<< data.repeats_done <<" of " << data.repeats << "in batch " << batch << "."<< std::endl;
	int j = 0;
	for (i = 2*data.N*(batch-1) ; i < 2*data.N*batch - 1 ; i += 4) {
	  //The magic aligment happens here: we have to change the phase of each next complex sample
	  //by pi - this means that even numbered samples stay the same while odd numbered samples
	  //get multiplied by -1 (thus rotated by pi in complex plane).
	  //This gets us output spectrum shifted by half its size - just what we need to get the output right.
	  data.inbuf[j/2][RE] = (double) data.buffer_arr[k][i] - 127;
	  data.inbuf[j/2][IM] = (double) data.buffer_arr[k][i + 1] - 127;
	  data.inbuf[j/2 + 1][RE] = ((double) data.buffer_arr[k][i+ 2] - 127) * -1;
	  data.inbuf[j/2 + 1][IM] = ((double) data.buffer_arr[k][i + 3] - 127) * -1;
	  j+=4;
	}
	fftw_execute(data.plan);
	for (i=0; i < data.N; i++) {
	  data.pwr[i] += data.outbuf[i][RE] * data.outbuf[i][RE] + data.outbuf[i][IM] * data.outbuf[i][IM];	
	}
	batch++;
	data.repeats_done++;
      }
      //Interpolate the central point, to cancel DC bias.
      data.pwr[data.N/2] = (data.pwr[data.N/2 - 1] + data.pwr[data.N/2+1]) / 2;
      data.buf_status[k] = 0;
      data.buf_mutex[k].unlock();
    }
    else {
      if (canary == 1) {
	break;
      }
      //usleep(500);
      if (data.acquisition_done == 1) {
	canary = 1;
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
  int r, num_of_gains, batches, overhang;
  int buf_length = 16384*100;
  double scans;
  int *gain_table;
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
    if (arg_repeats.isSet()) repeats = arg_repeats.getValue();
    if (arg_integration_time.isSet()) integration_time = arg_integration_time.getValue();
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
  r = rtlsdr_open(&dev, (uint32_t)dev_index);
  if (r < 0 ) {
    std::cerr << "Could not open rtl_sdr device " << dev_index << "." << std::endl;
    return r;
  }
  
  //Available gains
  num_of_gains = rtlsdr_get_tuner_gains(dev, NULL);
  gain_table = (int *) malloc(num_of_gains*sizeof(int));
  rtlsdr_get_tuner_gains(dev, gain_table);
  print_gain_table(num_of_gains, gain_table);
  gain = select_nearest_gain(gain, num_of_gains, gain_table);
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
  std::cerr << "Number of averaged samples: " << repeats << "." << std::endl;
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
  // should be a multiple of 16384. We have to keep it that way.
  scans = ceil((2*N*repeats)/buf_length);
  batches = buf_length/(2*N);
  overhang = buf_length % (2*N);
  if (overhang != 0) {
    buf_length = lcm(2*N, buf_length);
    scans = ceil((2*N*repeats)/buf_length);
    batches = buf_length/(2*N);
  }
  std::cerr << "Data collection will proceed in " << scans <<" scans, each consisting of " << batches << " batches." << std::endl;
  
  //Begin the work: prepare data buffers
  Datastore data(N, buf_length, batches, repeats);
  int i;
  for (i=0; i < N; i++) {
    data.pwr[i] = 0;
  }
  //Buffer usage stats
  int usage[BUFFERS] = {};
  //Read from device and do FFT
  std::thread t(&fft, std::ref(data));
  int count = 0;
  while (count <= scans) {
    int i = 0;
    int buf_available = 0;
    while (i < BUFFERS) {
      if (data.buf_mutex[i].try_lock()) { 
	if (data.buf_status[i] == 0) {
	  buf_available = 1;
          usage[i]++;
	  break;
	}
	else {
	  data.buf_mutex[i].unlock();
	}
      }
      i++;
    }
    if (buf_available == 1) {
      r = read_rtlsdr(data, i);
      if (r) {
	fprintf(stderr, "Error: dropped samples.\n");
      }
      else {
	count++;
	data.buf_status[i] = 1;
      }
      data.buf_mutex[i].unlock();
    }
  }
  std::cerr << "Acquisition_done." << std::endl;
  data.acquisition_done = 1;
  t.join();
  //Write out.
  for (i=0; i < N; i++) {
    std::cout << i << " " << tuned_freq + (i-N/2.0) * ( (N-1) / (double)N  * (double)actual_samplerate / (double)N ) << " " << 10*log10(data.pwr[i]/ repeats) << std::endl;
  }
  std::cerr << "Buffer usage: ";
  for (auto i : usage) std::cerr << i << ", ";
  std::cerr << std::endl;
  fftw_destroy_plan(data.plan);
  rtlsdr_close(dev);
  return 0;
}
