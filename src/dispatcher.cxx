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

#include "dispatcher.h"

#include "acquisition.h"
#include "fft_worker.h"
#include "params.h"

#include <iostream>


Dispatcher::Dispatcher(const Params& params_, const AuxData& aux_) :
  params(params_), aux(aux_)
{
  // Create a set of raw buffers.
  for (int i = 0; i < params.buffers; i++) {
    // Create a new buffer (as a member of a list).
    rawBuffers.emplace_back(params.buf_length);
    // Create an empty container with a pointer to the newly created buffer.
    emptyContainers.push_back({nullptr, &rawBuffers.back()});
  }

  // Create the FFT workers and add them to the queue of idle workers.
  for (int i = 0; i < params.threads; i++) {
    auto worker = new FFTWorker(*this);
    workers.push_back(worker);
    idleWorkers.push_back(worker);
  }

  // Run the dispatching and accumulating operations.
  dispatchingThread = std::thread(&Dispatcher::dispatchingOperation, this);
  accumulatingThread = std::thread(&Dispatcher::accumulatingOperation, this);
}

Dispatcher::~Dispatcher() {
  // Stop dispatching the data to workers. A container with a pointer to a null
  // acquisition denotes the end of the data stream and causes the dispatching
  // thread to terminate.
  occupiedContainers.push_back({nullptr, nullptr});
  dispatchingThread.join();

  // Dismiss all the workers. If any of these are still busy, the call to delete
  // will block until they finish their respective operations.
  for (auto& worker : workers)
    delete worker;

  // Queue a nullptr to notify the accumulating thread that there are no workers
  // left. The accumulating thread terminates upon receiving this message.
  finishedWorkers.push_back(nullptr);
  accumulatingThread.join();
}

void Dispatcher::workerFinished(FFTWorker* worker) {
  // A worker has finished its job. Queue it so that it will be dealt with by the
  // accumulating operation.
  finishedWorkers.push_back(worker);
}

void Dispatcher::dispatchingOperation() {
  DataContainer container;
  FFTWorker* worker = nullptr;

  // An index into the data buffer.
  unsigned int buffer_index = 0;
  // An index into the current worker's input buffer.
  int fft_index = 0;
  // The total number of repeats from the current acquisition that have been 
  // submitted for processing.
  int64_t totalRepeats = 0;
  bool newAcquisition = true;

  while (true) {
    // See if we have to fetch a new container. Otherwise we'll just continue
    // copying the data from the current one.
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
        // repeatsToProcess contains the value before addition.
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
        continue;
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
      // No worker assigned yet. Get an idle one.
      worker = idleWorkers.get();
      worker->acquisition = container.acquisition;
      fft_index = 0;
    }

    auto& data = *container.data;
    const auto dataSize = data.size();

    // Copy the data to the worker's input buffer until either the worker's input
    // buffer is full or we run out of data to copy.
    while (fft_index < params.N && buffer_index < dataSize) {
      // The "magic aligment" happens here: we change the phase of every second
      // complex sample by pi - this means that even-numbered samples stay the same
      // while odd-numbered samples get multiplied by -1 (thus rotated by pi in the
      // complex plane). This gets us an output spectrum shifted by half if its size
      // so that the DC component is in the middle.
      const float multiplier = (fft_index % 2 == 0 ? 1 : -1);
      complex bfr_val(data[buffer_index], data[buffer_index+1]);
      worker->inbuf[fft_index] = (bfr_val - complex(127.0, 127.0)) * multiplier;

      if (params.window)
        worker->inbuf[fft_index] *= aux.window_values[fft_index];

      buffer_index += 2; // Two bytes per complex sample.
      fft_index++;
    }

    if (fft_index == params.N) {
      // The worker's input buffer is full. Start the FFT.
      totalRepeats++;
      worker->startFFT();
      worker = nullptr;
    }

    if (buffer_index == dataSize) {
      // We have used all the data from this buffer, so it can now be reused.
      container.acquisition = nullptr;
      emptyContainers.push_back(container);
    }
  }
}

void Dispatcher::accumulatingOperation()
{
  FFTWorker* worker;

  // The loop runs until a nullptr is fetched instead of a pointer to a worker.
  while ((worker = finishedWorkers.get()) != nullptr) {
    auto acquisition = worker->acquisition;
    auto outbuf = worker->outbuf;

    // Add the results to the acquisition's accumulator. There is no need for
    // locking as this is the only thread that is allowed to access the
    // accumulator before the results are marked as ready.
    for (int i = 0; i < params.N; i++) {
      acquisition->pwr[i] += pow(outbuf[i].real(), 2) + pow(outbuf[i].imag(), 2);
    }

    idleWorkers.push_back(worker);
    acquisition->repeatsProcessed++;

    // When the dispatcher determines the total number of spectra passed off to
    // workers, it will atomically add that number to acquisition->repeatsToProcess
    // and check if the result is zero. If it is, it will mark the results as ready.
    // If it isn't, all the workers have not yet finished their job and this thread
    // is the one that will later mark the results as ready. It is important that
    // this marking only happens in one thread because once markResultsReady() is
    // called, the state of the Acquisition object immediately becomes indeterminate.
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
