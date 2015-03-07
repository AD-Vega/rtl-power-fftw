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
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>

#include <pthread.h>
//#include <libusb.h>
#include <rtl-sdr.h>
//#include <rtl-sdr/convenience/convenience.h>
#define RE 0
#define IM 1

static rtlsdr_dev_t *dev = NULL;

int main(int argc, char **argv)
{
  int N = 1024;
  int dev_index = 0;
  int r;
  int gain = 290;
  int buf_len, n_read;
  buf_len = 2 * N;
  uint8_t *buf8;
  buf8 = (uint8_t *) malloc (buf_len * sizeof(uint8_t));
  //dev_index = verbose_device_search("0");
  r = rtlsdr_open(&dev, (uint32_t)dev_index);
  rtlsdr_set_tuner_gain(dev, gain);
  rtlsdr_set_center_freq(dev, (uint32_t)89300000);
  usleep(5000);
  rtlsdr_reset_buffer(dev);
  rtlsdr_set_sample_rate(dev, (uint32_t)1800000);
  rtlsdr_read_sync(dev, buf8, buf_len, &n_read);
  rtlsdr_close(dev);
  if (n_read != buf_len) {
    fprintf(stderr, "Error: dropped samples.\n");}
  fftw_complex *inbuf, *outbuf;
  inbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  outbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  fftw_plan plan = fftw_plan_dft_1d(N, inbuf, outbuf, FFTW_FORWARD, FFTW_MEASURE);
  int i;
  for (i = 0 ; i < buf_len ; i += 2) {
    inbuf[i/2][RE] = (double) buf8[i];
    inbuf[i/2][IM] = (double) buf8[i + 1];
    //printf("%i\t%g\t%g\t%g\n", i, inbuf[i][RE], inbuf[i][IM], w);
  }
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  for (i=0; i < N; i++) {
    printf("%i\t%g\t%g\t%g\t%g\t%g\n", i, inbuf[i][RE], inbuf[i][IM], outbuf[i][RE], outbuf[i][IM], sqrt(outbuf[i][RE] * outbuf[i][RE] + outbuf[i][IM] * outbuf[i][IM]));
  }
  
  return 0;
}
