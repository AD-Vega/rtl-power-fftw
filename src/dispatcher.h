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

#include <atomic>
#include <complex>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <thread>
#include <vector>

#include "device.h"
#include "params.h"


template <typename T>
class ConcurrentQueue {
public:
  void push_back(T item) {
    std::lock_guard<std::mutex> queueLock(mutex);
    queue.push_back(item);
    event.notify_one();
  }

  void push_front(T item) {
    std::lock_guard<std::mutex> queueLock(mutex);
    queue.push_front(item);
    event.notify_one();
  }

  T get() {
    std::unique_lock<std::mutex> queueLock(mutex);
    return getItem(queueLock);
  }

  T get(size_t& queueSize) {
    std::unique_lock<std::mutex> queueLock(mutex);
    queueSize = queue.size();
    return getItem(queueLock);
  }

protected:
  T getItem(std::unique_lock<std::mutex>& queueLock) {
    while (queue.empty())
      event.wait(queueLock);
    T item = queue.front();
    queue.pop_front();
    return item;
  }

  std::mutex mutex;
  std::condition_variable event;
  std::deque<T> queue;
};


class Acquisition;

struct DataContainer {
  DataContainer() = default;
  DataContainer(Acquisition* acquisition_, RawBuffer* data_) :
    acquisition(acquisition_), data(data_)
  {}

  Acquisition* acquisition = nullptr;
  RawBuffer* data = nullptr;
};


class FFTWorker;

class Dispatcher {
public:
  Dispatcher(const Params& params, const AuxData& aux, int numberOfFFTWorkers);
  ~Dispatcher();

  // We have threads going on; don't allow copying the object.
  // If moving is ever needed, we'll implement it at that time.
  Dispatcher(const Dispatcher&) = delete;
  Dispatcher(Dispatcher&&) = delete;
  Dispatcher& operator=(const Dispatcher&) = delete;
  Dispatcher& operator=(const Dispatcher&&) = delete;

  void workerFinished(FFTWorker* worker);

  const Params& params;
  const AuxData& aux;
  ConcurrentQueue<DataContainer> emptyContainers;
  ConcurrentQueue<DataContainer> occupiedContainers;

protected:
  void dispatchingOperation();
  void accumulatingOperation();

  std::list<RawBuffer> rawBuffers;
  std::list<FFTWorker*> workers;
  ConcurrentQueue<FFTWorker*> idleWorkers;
  ConcurrentQueue<FFTWorker*> finishedWorkers;
  std::thread dispatchingThread;
  std::thread accumulatingThread;
};

#endif // DISPATCHER_H
