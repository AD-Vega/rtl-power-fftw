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

#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include "exceptions.h"

const int base_buf = 16384;
const int default_buf_multiplier = 100;

class Params {
public:
  Params(int argc, char **argv);

  int N = 512;
  int dev_index = 0;
  int gain = 372;
  int64_t cfreq = 1420405752;
  int64_t startfreq = 0; 
  int64_t stopfreq = 0;
  int sample_rate = 2000000;
  double integration_time = 0;
  bool integration_time_isSet = false;
  int buffers = 5;
  int buf_length = base_buf * default_buf_multiplier;
  bool buf_length_isSet = false;
  double min_overlap = 0;
  int ppm_error = 0;
  bool endless = false;
  bool strict_time = false;
  bool baseline = false;
  std::string baseline_file;
  bool window = false;
  std::string window_file;
  bool freq_hopping_isSet = false;
  //It is senseless to waste a full buffer of data unless instructed to do so.
  int64_t repeats = buf_length/(2*N);
  int threads = 1;
  std::string outputFile;
  int session_duration = -1;
  bool linear = false;
};

#endif // PARAMS_H
