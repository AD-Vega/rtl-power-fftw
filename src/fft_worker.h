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

#include <complex>
#include <condition_variable>
#include <fftw3.h>
#include <mutex>
#include <thread>

class Acquisition;
class Dispatcher;

using complex = std::complex<float>;

// Each FFT worker spawns a separate thread for data processing. Data is passed
// in from the dispatching thread and is reaped by the accumulating thread of
// the Dispatcher.
class FFTWorker {
public:
  FFTWorker(Dispatcher& dispatcher);
  ~FFTWorker();
  // Trigger the processing of the data in the input buffer.
  void startFFT();

  // Delete these so we don't accidentally mess anything up by copying
  // pointers to fftw_malloc'd buffers.
  FFTWorker(const FFTWorker&) = delete;
  FFTWorker(FFTWorker&&) = delete;
  FFTWorker& operator=(const FFTWorker&) = delete;
  FFTWorker& operator=(FFTWorker&&) = delete;

  // A pointer to the acquisition to which the current data belongs.
  Acquisition* acquisition = nullptr;
  // Input buffer for FFTW.
  complex *inbuf;
  // Output buffer for FFTW.
  complex *outbuf;

protected:
  // This function will be started in a separate thread.
  void fftOperation();

  // A reference to the dispatcher (so we can notify it when we are done).
  Dispatcher& dispatcher;
  // Number of samples in the spectrum.
  const int N;
  // FFTW plan.
  fftwf_plan plan;

  // A mutex that guards the controlEvent variable and the doFFT and doExit flags.
  std::mutex controlMutex;
  // Used to signal an event (new data arrived, termination requested) to the worker.
  std::condition_variable controlEvent;
  // True if a FFT operation is requested; reset by the worker after the operation
  // is finished.
  bool doFFT = false;
  // True if the worker is requested to terminate.
  bool doExit = false;
  // The worker thread.
  std::thread fftThread;
  // Static mutex guarding the non-thread-safe functions in FFTW.
  static std::mutex fftwMutex;
};

#endif // FFT_WORKER_H
