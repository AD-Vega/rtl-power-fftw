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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "utilities.h"
#include <memory>
#include <thread>

class Acquisition;
class Params;
class AuxData;


class OutputStream {
public:
  virtual ~OutputStream() {};
  virtual void write(Acquisition& acquisition) = 0;
};


class TextStream : public OutputStream {
public:
  TextStream(const Params& params, AuxData& aux);
  void write(Acquisition& acquisition);

protected:
  const Params& params;
  const AuxData& aux;
};


class OutputWriter {
public:
  OutputWriter(OutputStream* stream);
  ~OutputWriter();

  // We have threads going on; don't allow copying the object.
  // If moving is ever needed, we'll implement it at that time.
  OutputWriter(const OutputWriter&) = delete;
  OutputWriter(OutputWriter&&) = delete;
  OutputWriter& operator=(const OutputWriter&) = delete;
  OutputWriter& operator=(const OutputWriter&&) = delete;

  ConcurrentQueue<std::shared_ptr<Acquisition>> queue;

protected:
  void run();

  OutputStream* stream;
  std::thread thread;
};

#endif // OUTPUT_H
