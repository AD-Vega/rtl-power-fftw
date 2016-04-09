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

#include "fft_worker.h"

#include "acquisition.h"
#include "dispatcher.h"
#include "params.h"

std::mutex FFTWorker::fftwMutex;

FFTWorker::FFTWorker(Dispatcher& dispatcher_) :
  dispatcher(dispatcher_), N(dispatcher.params.N)
{
  {
    // From all the FFTW functions, fftw_execute() is the only one that is
    // guaranteed to be thread-safe. The allocations, plans etc. have to be
    // made only in one thread at once.
    std::lock_guard<std::mutex> guard(fftwMutex);

    inbuf = (complex*)fftwf_alloc_complex(N);
    outbuf = (complex*)fftwf_alloc_complex(N);
    plan = fftwf_plan_dft_1d(N, (fftwf_complex*)inbuf, (fftwf_complex*)outbuf,
                            FFTW_FORWARD, FFTW_MEASURE);
  }

  fftThread = std::thread(&FFTWorker::fftOperation, this);
}

FFTWorker::~FFTWorker() {
  {
    std::lock_guard<std::mutex> controlLock(controlMutex);
    doExit = true;
    controlEvent.notify_one();
  }
  fftThread.join();

  // Deallocate data.
  // Protect access to FFTW calls.
  std::lock_guard<std::mutex> guard(fftwMutex);
  fftwf_destroy_plan(plan);
  fftwf_free(inbuf);
  fftwf_free(outbuf);
}

void FFTWorker::startFFT() {
  std::lock_guard<std::mutex> controlLock(controlMutex);
  doFFT = true;
  controlEvent.notify_one();
}

void FFTWorker::fftOperation() {
  std::unique_lock<std::mutex> controlLock(controlMutex);
  while (!doExit) {
    while (!doFFT && !doExit)
      controlEvent.wait(controlLock);

    if (doExit)
      break;

    fftwf_execute(plan);
    doFFT = false;
    dispatcher.workerFinished(this);
  }
}
