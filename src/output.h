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
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <thread>

class Acquisition;
class Params;
class AuxData;


class OutputStream {
public:
  virtual ~OutputStream() {};
  virtual void write(Acquisition& acquisition) = 0;
  virtual void writeDelimiter() = 0;
};


class TextStream : public OutputStream {
public:
  TextStream(const Params& params, AuxData& aux);
  void write(Acquisition& acquisition);
  void writeDelimiter();

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

  void queueData(std::shared_ptr<Acquisition> acquisition);
  void queueDelimiter();

protected:
  void run();

  // A message that can be passed to the writing thread.
  // It contains a command and, optionally, a shared ptr to Acquisition.
  struct Message {
    // Possible commands.
    enum class Command {
      WriteData,
      WriteDelimiter,
      Quit
    };

    Message(Command command_, std::shared_ptr<Acquisition> acquisition_ = nullptr) :
      command(command_), acquisition(acquisition_)
      {}

    // The message evaluates to true if the command is not Quit.
    operator bool() const {
      return command != Command::Quit;
    }

    Command command;
    std::shared_ptr<Acquisition> acquisition;
  };

  ConcurrentQueue<Message> queue;
  OutputStream* stream;
  std::thread thread;
};


enum class LogLevel {
  Debug,
  Operation,
  Info,
  Warning,
  Error
};


// This class is used to write diagnostic output to stderr. The writes are
// internally synchronized so that multiple threads do not interfere with
// each other. The usage pattern is as follows: create a Diagnostics() object
// (usually an anonymous temporary), then use the << operator to write to it
// as usual. The object will accumulate the messages into a buffer and when it
// is destroyed, it will send the whole content of the buffer at once to stderr.
// Only at this last point is the access to stderr synchronized so that the
// output stream remains locked for as little time as possible.
//
// In addition, a logging level can be supplied: the message will be printed
// only if the level exceeds a certain (global) threshold.

class Diagnostics {
public:
  Diagnostics(LogLevel level_ = LogLevel::Info) : level(level_) {}

  // Copying and moving is not allowed.
  Diagnostics(const Diagnostics&) = delete;
  Diagnostics(Diagnostics&&) = delete;
  Diagnostics& operator=(const Diagnostics&) = delete;
  Diagnostics& operator=(Diagnostics&&) = delete;

  // The destructor will actually trigger the writing to the output stream.
  ~Diagnostics();

  // Set the global logging threshold. A diagnostic message will only be
  // printed if its logging level is equal to or greater than the supplied
  // threshold.
  //
  // Warning: this function is not thread safe! It is meant to be called once
  // at the beginning of the program, before other threads are created.
  static void setThreshold(LogLevel threshold);

  // Operator << passes the content on to the buffer (or to a null stream).
  template <typename T>
  std::ostream& operator<<(const T& item) {
    if (level >= threshold)
      return buf << item;
    else
      return nullStream;
  }

protected:
  // The log level of the currently assembled message.
  LogLevel level;
  // The buffer for the currently assembled message.
  std::ostringstream buf;

  // For synchronizing access to the output stream.
  static std::mutex mutex;
  // Global log threshold.
  static LogLevel threshold;
  // A hack: we use an unopened std::ofstream as a null stream. Writing to this
  // object will simply ignore the input and will not produce any output at all.
  static std::ofstream nullStream;
};

#endif // OUTPUT_H
