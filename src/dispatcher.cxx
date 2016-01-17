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

#include "acquisition.h"
#include "dispatcher.h"
#include "fft_worker.h"

#include <iostream>
#include <bits/algorithmfwd.h>


Dispatcher::Dispatcher(const Params& params_, const AuxData& aux_, int numberOfFFTWorkers) :
  params(params_), aux(aux_)
{
  for (int i = 0; i < params.buffers; i++) {
    rawBuffers.emplace_back(params.buf_length);
    emptyContainers.push_back({nullptr, &rawBuffers.back()});
  }

  for (int i = 0; i < numberOfFFTWorkers; i++) {
    auto worker = new FFTWorker(*this);
    workers.push_back(worker);
    idleWorkers.push_back(worker);
  }

  dispatchingThread = std::thread(&Dispatcher::dispatchingOperation, this);
  accumulatingThread = std::thread(&Dispatcher::accumulatingOperation, this);
}

Dispatcher::~Dispatcher() {
  // Stop dispatching the data to workers. A nullptr in the buffer queue denotes
  // the end of the data stream and causes the dispatching thread to terminate.
  occupiedContainers.push_back({nullptr, nullptr});
  dispatchingThread.join();

  // Dismiss all the workers. If any of these are still busy, the call to delete
  // will block until they finish their respective operations.
  for (auto& worker : workers)
    delete worker;

  // Queue a nullptr to notify the accumulating thread that there are no workers
  // left.
  finishedWorkers.push_back(nullptr);
  accumulatingThread.join();
}

void Dispatcher::workerFinished(FFTWorker* worker) {
  finishedWorkers.push_back(worker);
}

void Dispatcher::dispatchingOperation() {
  DataContainer container;
  FFTWorker* worker = nullptr;

  unsigned int buffer_index = 0;
  int fft_index = 0;
  int64_t totalRepeats = 0;
  bool newAcquisition = true;

  while (true) {
    if (!container.acquisition) {
      container = occupiedContainers.get();

      if (container.acquisition == nullptr) {
        // This is a signal that data acquisition has stopped and that no
        // additional data will follow. We can terminate the dispatching
        // operation.
        return;
      }

      if (container.data == nullptr) {
        // We have received a sentinel denoting the end of data for this
        // particular acquisition. Check whether all the data has already
        // been processed. Use an atomic fetch_add operation to avoid racing
        // with the accumulating thread.
        auto repeatsToProcess = container.acquisition->repeatsToProcess.fetch_add(totalRepeats);
        repeatsToProcess += totalRepeats;
        if (repeatsToProcess == 0) {
          // All the data has been processed. Mark the results as ready.
          container.acquisition->markResultsReady();
        }
        else {
          // Data processing still running; the accumulating thread will mark
          // the result as ready when the time comes.
        }

        // Discard any partial data that could have been already copied to the
        // worker. We do this by simply returning the worker to the queue.
        if (worker) {
          idleWorkers.push_front(worker);
          worker = nullptr;
        }

        // Invalidate the current container and restart the loop in order to
        // fetch a new one.
        container.acquisition = nullptr;
        newAcquisition = true;
        break;
      }

      if (newAcquisition) {
        // This is a beginning of a new acquisition. Initialize counters.
        container.acquisition->repeatsProcessed = 0;
        container.acquisition->repeatsToProcess = 0;
        totalRepeats = 0;
        newAcquisition = false;
      }

      // We have a new container. Reset the buffer index.
      buffer_index = 0;
    }

    if (!worker) {
      worker = idleWorkers.get();
      worker->acquisition = container.acquisition;
      fft_index = 0;
    }

    auto& data = *container.data;
    while (fft_index < params.N && buffer_index < data.size()) {
      //The magic aligment happens here: we have to change the phase of each next complex sample
      //by pi - this means that even numbered samples stay the same while odd numbered samples
      //get multiplied by -1 (thus rotated by pi in complex plane).
      //This gets us output spectrum shifted by half its size - just what we need to get the output right.
      const float multiplier = (fft_index % 2 == 0 ? 1 : -1);
      complex bfr_val(data[buffer_index], data[buffer_index+1]);
      worker->inbuf[fft_index] = (bfr_val - complex(127.0, 127.0)) * multiplier;
      if (params.window)
        worker->inbuf[fft_index] *= aux.window_values[fft_index];
      buffer_index += 2;
      fft_index++;
    }

    if (fft_index == params.N) {
      totalRepeats++;
      worker->startFFT();
      worker = nullptr;
    }

    if (buffer_index == data.size()) {
      // We have used all the data from this buffer, so it can now be reused.
      container.acquisition = nullptr;
      emptyContainers.push_back(container);
    }
  }
}

void Dispatcher::accumulatingOperation()
{
  FFTWorker* worker;

  while ((worker = finishedWorkers.get()) != nullptr) {
    auto acquisition = worker->acquisition;
    auto outbuf = worker->outbuf;

    for (int i = 0; i < params.N; i++) {
      acquisition->pwr[i] += pow(outbuf[i].real(), 2) + pow(outbuf[i].imag(), 2);
    }

    idleWorkers.push_back(worker);
    acquisition->repeatsProcessed++;

    if (--(acquisition->repeatsToProcess) == 0) {
      // All the data has been processed. Mark the results as ready.
      acquisition->markResultsReady();
    }
    else {
      // The dispatching thread has not yet determined the total number of
      // repeats, so we don't know whether the processing for this acquisition
      // is finished. If it turns out that no more data follows, the dispatcher
      // will mark the results as ready.
    }
  }
}
