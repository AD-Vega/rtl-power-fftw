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

#ifndef FFT_WORKER_H
#define FFT_WORKER_H

#include "dispatcher.h"

#include <complex>
#include <fftw3.h>
#include <mutex>
#include <thread>

using complex = std::complex<float>;

class FFTWorker {
public:
  FFTWorker(Dispatcher& dispatcher);
  ~FFTWorker();
  void startFFT();

  // Delete these so we don't accidentally mess anything up by copying
  // pointers to fftw_malloc'd buffers.
  FFTWorker(const FFTWorker&) = delete;
  FFTWorker(FFTWorker&&) = delete;
  FFTWorker& operator=(const FFTWorker&) = delete;
  FFTWorker& operator=(FFTWorker&&) = delete;

  Acquisition* acquisition = nullptr;
  complex *inbuf, *outbuf;

protected:
  Dispatcher& dispatcher;
  const int N;
  fftwf_plan plan;

  std::mutex controlMutex;
  std::condition_variable controlEvent;
  bool doFFT = false;
  bool doExit = false;

  static std::mutex fftwMutex;
  std::thread fftThread;
  // This function will be started in a separate thread.
  void fftOperation();
};

#endif // FFT_WORKER_H
