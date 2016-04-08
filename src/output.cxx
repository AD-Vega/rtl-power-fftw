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

#include "output.h"
#include "acquisition.h"

#include <cmath>
#include <iomanip>
#include <iostream>


TextStream::TextStream(const Params& params_, AuxData& aux_) :
  params(params_), aux(aux_)
{
  if (params.outputFile.empty()) {
    stream = &std::cout;
  }
  else {
    fstream.reset(new std::ofstream(params.outputFile));
    if (*fstream) {
      stream = fstream.get();
    }
    else {
      throw RPFexception(
        "Could not open " + params.outputFile + " for writing.",
        ReturnValue::InvalidArgument);
    }
  }
}

void TextStream::write(Acquisition& acq) {
  std::string units = params.linear ? "power/Hz" : "dB/Hz";

  // Print the header
  *stream << "# rtl-power-fftw output\n"
          << "# Acquisition start: " << acq.startAcqTimestamp << "\n"
          << "# Acquisition end: " << acq.endAcqTimestamp << "\n"
          << "#\n"
          << "# frequency [Hz] power spectral density [" << units << "]\n";

  // Calculate the precision needed for displaying the frequency.
  const int extraDigitsFreq = 2;
  const int significantPlacesFreq = ceil(floor(log10(acq.tuned_freq))
    - log10(acq.actual_samplerate/params.N) + 1 + extraDigitsFreq);
  const int significantPlacesPwr = 6;

  for (int i = 0; i < params.N; i++) {
    double value = acq.pwr[i];
    if (!params.linear) {
      // Convert to dB.
      value = 10 * log10(value);
    }
    *stream << std::setprecision(significantPlacesFreq) << acq.frequency(i)
              << " "
              << std::setprecision(significantPlacesPwr) << value
              << "\n";
  }
  // Separate consecutive spectra with empty lines.
  *stream << std::endl;
}

void TextStream::writeDelimiter() {
  // Measurement sets are delimited with two empty lines. One of them was already
  // printed as a part of the last result.
  *stream << std::endl;
}


OutputWriter::OutputWriter(OutputStream* stream_) :
  stream(stream_)
{
  thread = std::thread(&OutputWriter::run, this);
}

OutputWriter::~OutputWriter() {
  queue.push_back(Message::Command::Quit);
  thread.join();
}

void OutputWriter::queueData(std::shared_ptr<Acquisition> acquisition) {
  queue.push_back({Message::Command::WriteData, acquisition});
}

void OutputWriter::queueDelimiter() {
  queue.push_back(Message::Command::WriteDelimiter);
}

void OutputWriter::run() {
  while (auto message = queue.get()) {
    if (message.command == Message::Command::WriteData) {
      message.acquisition->waitForResultsReady();
      message.acquisition->printSummary();
      stream->write(*message.acquisition);
    }
    if (message.command == Message::Command::WriteDelimiter) {
      stream->writeDelimiter();
    }
  }
}


std::mutex Diagnostics::mutex;
LogLevel Diagnostics::threshold = LogLevel::Operation;
std::ofstream Diagnostics::nullStream;

Diagnostics::~Diagnostics() {
  if (level < threshold)
    return;

  std::lock_guard<std::mutex> lock(mutex);
  std::cerr << buf.str() << std::flush;
}

void Diagnostics::setThreshold(LogLevel threshold_) {
  threshold = threshold_;
}
