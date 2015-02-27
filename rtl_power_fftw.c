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
#define RE 0
#define IM 1
int main(int argc, char **argv)
{
  int N = 32;
  fftw_complex *inbuf, *outbuf;
  inbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  outbuf = (fftw_complex *) malloc (N*sizeof(fftw_complex));
  fftw_plan plan = fftw_plan_dft_1d(N, inbuf, outbuf, FFTW_FORWARD, FFTW_MEASURE);
  int i; 
  double w;
  for (i=0; i < N; i++) {
    w = (double)i / (double)N * M_PI;
    inbuf[i][RE] = 2.0 * cos(6.0*w) + cos(15.0*w);
    inbuf[i][IM] = 2.0 * sin(6.0*w) + sin(15.0*w);
    printf("%i\t%g\t%g\t%g\n", i, inbuf[i][RE], inbuf[i][IM], w);
  }
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  for (i=0; i < N; i++) {
    printf("%i\t%g\t%g\t%g\t%g\t%g\n", i, inbuf[i][RE], inbuf[i][IM], outbuf[i][RE], outbuf[i][IM], sqrt(outbuf[i][RE] * outbuf[i][RE] + outbuf[i][IM] * outbuf[i][IM]));
  }  
  return 0;
}
