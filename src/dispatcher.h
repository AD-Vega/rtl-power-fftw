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
  Dispatcher(const Params& params, const AuxData& aux);
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
