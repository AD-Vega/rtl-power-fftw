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
#include <pthread.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
//#include <libusb.h>
#include <rtl-sdr.h>
//#include <rtl-sdr/convenience/convenience.h>
#define RE 0
#define IM 1

static rtlsdr_dev_t *dev = NULL;

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


int main(int argc, char **argv)
{
  int N = 512;
  int repeats = 1;
  int dev_index = 0;
  int gain = 372;
  int cfreq = 89600000;
  int sample_rate = 2000000;
  int r, num_of_gains, buf_len, n_read;
  int *gain_table;
  try {
    TCLAP::CmdLine cmd("Obtain power spectrum from RTL device using FFTW library.", ' ', "0.1");
    TCLAP::ValueArg<int> arg_bins("b","bins","Number of bins in FFT spectrum",false,512,"even integer > 0");
    cmd.add( arg_bins );
    TCLAP::ValueArg<int> arg_freq("f","freq","Center frequency of the receiver.",false,89100000,"Hz");
    cmd.add( arg_freq );
    TCLAP::ValueArg<int> arg_rate("r","rate","Sample rate of the receiver.",false,2000000,"samples/s");
    cmd.add( arg_rate );
    TCLAP::ValueArg<int> arg_gain("g","gain","Receiver gain.",false, 372, "cB (1/10th of dB)");
    cmd.add( arg_gain );
    TCLAP::ValueArg<int> arg_repeats("n","repeats","Number of scans for averaging.",false,1,"repeats");
    cmd.add( arg_repeats );
    TCLAP::ValueArg<int> arg_index("d","device","RTL-SDR device index.",false,0,"device index");
    cmd.add( arg_index );
    
    cmd.parse(argc, argv);
    
    dev_index = arg_index.getValue();
    N = arg_bins.getValue();
    repeats = arg_repeats.getValue();
    gain = arg_gain.getValue();
    cfreq = arg_freq.getValue();
    sample_rate = arg_rate.getValue();
  }
  catch (TCLAP::ArgException &e) { 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
  }
  //Sanity checks
  int num_of_rtls = rtlsdr_get_device_count();
  if (num_of_rtls == 0) {
    std::cerr << "Error: no RTL-SDR compatible devices found. Exiting." << std::endl;
    return -1;
  }
  if ( dev_index >= num_of_rtls) {
    std::cerr << "Error: invalid device number. Only "<< num_of_rtls << " devices available. Exiting." << std::endl;
    return -2;
  } 
  r = rtlsdr_open(&dev, (uint32_t)dev_index);
  num_of_gains = rtlsdr_get_tuner_gains(dev, NULL);
  gain_table = (int *) malloc(num_of_gains*sizeof(int));
  rtlsdr_get_tuner_gains(dev, gain_table);
  print_gain_table(num_of_gains, gain_table);
  gain = select_nearest_gain(gain, num_of_gains, gain_table);
  std::cerr << "Selected nearest available gain: " << gain << std::endl;
  rtlsdr_set_tuner_gain_mode(dev, 1);
  rtlsdr_set_tuner_gain(dev, gain);
  rtlsdr_set_center_freq(dev, (uint32_t)cfreq);
  int tuned_freq = rtlsdr_get_center_freq(dev);
  std::cerr << "Device tuned to: " << tuned_freq << " Hz." << std::endl;
  usleep(5000);
  rtlsdr_set_sample_rate(dev, (uint32_t)sample_rate);
  int actual_samplerate = rtlsdr_get_sample_rate(dev);
  if (N % 2 != 0) {
    N++;
    std::cerr << "Number of bins should be an even number, changing to " << N << "." << std::endl;
  }
  std::cerr << "Number of bins: " << N << std::endl;
  
  buf_len = 2 * N;
  uint8_t *buf8;
  buf8 = (uint8_t *) malloc (buf_len * sizeof(uint8_t));
  fftw_complex *inbuf, *outbuf, *accumulator;
  double *pwr;
  inbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  outbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  pwr = (double *) malloc (N*sizeof(double));
  fftw_plan plan = fftw_plan_dft_1d(N, inbuf, outbuf, FFTW_FORWARD, FFTW_MEASURE);

  
  
  int count = 0;
  int i;
  for (i=0; i < N; i++) {
    pwr[i] = 0;
  }
  while (count <= repeats) {
    rtlsdr_reset_buffer(dev);
    rtlsdr_read_sync(dev, buf8, buf_len, &n_read);
    if (n_read != buf_len) {
      fprintf(stderr, "Error: dropped samples.\n");}
    else {
      for (i = 0 ; i < buf_len ; i += 4) {
	inbuf[i/2][RE] = (double) buf8[i] - 127;
	inbuf[i/2][IM] = (double) buf8[i + 1] - 127;
	inbuf[i/2 + 1][RE] = ((double) buf8[i+ 2] - 127) * -1;
	inbuf[i/2 + 1][IM] = ((double) buf8[i + 3] - 127) * -1;
	//printf("%i\t%g\t%g\t%g\n", i, inbuf[i][RE], inbuf[i][IM], w);
      }
      fftw_execute(plan);
      for (i=0; i < N; i++) {
	pwr[i] += sqrt(outbuf[i][RE] * outbuf[i][RE] + outbuf[i][IM] * outbuf[i][IM]);	
      }
      count++;
    }
  }
  pwr[N/2] = (pwr[N/2 - 1] + pwr[N/2+1]) / 2;
  for (i=0; i < N; i++) {
    //printf("%i\t%g\t%g\t%g\t%g\t%g\n", i, inbuf[i][RE], inbuf[i][IM], outbuf[i][RE], outbuf[i][IM], pwr[i]);
    std::cout << i << " " << tuned_freq + (i-N/2.0) * ( (N-1) / (double)N  * (double)actual_samplerate / (double)N ) << " " << pwr[i] << std::endl;
  }
  fftw_destroy_plan(plan);
  rtlsdr_close(dev);
  return 0;
}
