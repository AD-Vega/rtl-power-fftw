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


// Class OutputStream is an interface through which new output formats can be
// added in a modular way. To use a particular format, implement the interface
// functions in a new class inheriting from OutputStream, create an instance of
// that class and pass a pointer to the instance to OutputWriter.
class OutputStream {
public:
  virtual ~OutputStream() {};

  // Called before a frequency scan begins.
  virtual void datasetBegin() = 0;
  // Called after a frequency scan ends.
  virtual void datasetEnd() = 0;
  // Called after every acquisition.
  virtual void write(Acquisition& acquisition) = 0;
};


// A text-based output stream that can write to either stdout or a file. The
// output contains a simple header preceding each spectrum. Consecutive spectra
// are separated by a single blank line and measurement sets are separated by
// two blank lines. This format is directly suitable as an input to gnuplot.
class TextStream : public OutputStream {
public:
  TextStream(const Params& params, AuxData& aux);
  void datasetBegin() { /* no-op */ }
  void datasetEnd();
  void write(Acquisition& acquisition);

protected:
  const Params& params;
  const AuxData& aux;
  // A pointer to the output stream (either stdout or a file).
  std::ostream* stream;
  // A unique_ptr to a file stream (if needed).
  std::unique_ptr<std::ofstream> fstream;
};


// An instance of this class runs a separate thread that receives finalized
// data from an acquisition and serializes it in a particular format.
class OutputWriter {
public:
  // The argument 'stream' can be a pointer to any object that implements the
  // OutputStream interface.
  OutputWriter(OutputStream* stream);
  ~OutputWriter();

  // We have threads going on; don't allow copying the object.
  // If moving is ever needed, we'll implement it at that time.
  OutputWriter(const OutputWriter&) = delete;
  OutputWriter(OutputWriter&&) = delete;
  OutputWriter& operator=(const OutputWriter&) = delete;
  OutputWriter& operator=(const OutputWriter&&) = delete;

  // Mark the beginning of a data set.
  void datasetBegin();
  // Mark the end of a data set.
  void datasetEnd();
  // Queue the data to be written.
  void queueData(std::shared_ptr<Acquisition> acquisition);

protected:
  // The function that runs in a separate thread.
  void run();

  // A message that can be passed to the writing thread.
  // It contains a command and, optionally, a shared ptr to Acquisition.
  struct Message {
    // Possible commands.
    enum class Command {
      DatasetBegin,
      DatasetEnd,
      WriteData,
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

  // The queue through which the writing thread receives notifications.
  ConcurrentQueue<Message> queue;
  // The output stream.
  OutputStream* stream;
  // The thread object in which run() gets invoked.
  std::thread thread;
};


// Loggin levels of diagnostic messages.
enum class LogLevel {
  Debug, // Debug messages.
  Operation, // Messages that get repeated with every acquisition.
  Info, // Messages that (more or less) appear once per program run.
  Warning, // For conditions that affect, but don't terminate, the operation.
  Error // For conditions that cause the program to terminate.
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
