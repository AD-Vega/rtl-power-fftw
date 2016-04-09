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

#ifndef DISPATCHER_H
#define DISPATCHER_H

#include "device.h"
#include "utilities.h"

#include <atomic>
#include <deque>
#include <list>
#include <thread>

class Params;
class AuxData;
class Acquisition;
class FFTWorker;


// DataContainer: a lightweight struct that relates a buffer of data to the
// Acquisition object to which the data belongs. Such data containers are
// used for passing the data between an Acquisition and the Dispatcher.
struct DataContainer {
  DataContainer() = default;
  DataContainer(Acquisition* acquisition_, RawBuffer* data_) :
    acquisition(acquisition_), data(data_)
  {}

  Acquisition* acquisition = nullptr;
  RawBuffer* data = nullptr;
};

// The dispatcher's role is to run two dedicated threads:
//
// * One thread receives the data that comes from the device (the Dispatcher itself
//   does not do any device-related work). The data is then dispatched to FFT
//   workers and the data buffers (once empty) queued to be reused by the currently
//   running acquisition.
//
// * The other thread is notified when a FFT worker has finished processing its data;
//   the thread takes the resulting spectrum and adds it to the accumulator of the
//   corresponding Acquisition.
//
class Dispatcher {
public:
  Dispatcher(const Params& params, const AuxData& aux);
  ~Dispatcher();

  // We have threads going on; don't allow copying the object.
  // If moving is ever needed, we'll implement it at that time.
  Dispatcher(const Dispatcher&) = delete;
  Dispatcher(Dispatcher&&) = delete;
  Dispatcher& operator=(const Dispatcher&) = delete;
  Dispatcher& operator=(const Dispatcher&&) = delete;

  // Any FFT worker that finishes its job calls this function and passes a pointer
  // to self as the argument.
  void workerFinished(FFTWorker* worker);

  // A reference to params. Declared public to allow access by FFT workers.
  const Params& params;
  // A reference to auxiliary data. Declared public to allow access by FFT workers.
  const AuxData& aux;
  // A queue holding empty containers. The Acquisition object will pop a container
  // from the queue, fill the associated buffer with data and push the container
  // to occupiedContainers.
  ConcurrentQueue<DataContainer> emptyContainers;
  // A queue holding occupied containers. A container will get popped by the
  // Dispatcher, the data distributed amongst the workers and the container
  // returned to emptyContainers to be reused.
  ConcurrentQueue<DataContainer> occupiedContainers;

protected:
  // The dispatching routine that runs in a separate thread.
  void dispatchingOperation();
  // The accumulating routine that runs in a separate thread.
  void accumulatingOperation();

  // A list of raw buffers. This list is used for easy allocation of the buffers
  // and raw buffer pointers are then used in real operation.
  std::list<RawBuffer> rawBuffers;
  // A list of pointers to FFT workers. This list is used to keep track of the
  // workers so that they can be terminated and destroyed when the Dispatcher's
  // destructor runs.
  std::list<FFTWorker*> workers;
  // A queue of idle workers.
  ConcurrentQueue<FFTWorker*> idleWorkers;
  // A queue of finished workers, i.e., those that have finished a FFT operation
  // and are waiting for the accumulating thread to reap the results.
  ConcurrentQueue<FFTWorker*> finishedWorkers;
  // The dispatching thread (see dispatchingOperation()).
  std::thread dispatchingThread;
  // The accumulating thread (see accumulatingOperation()).
  std::thread accumulatingThread;
};

#endif // DISPATCHER_H
