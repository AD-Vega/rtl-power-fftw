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
{}

void TextStream::write(Acquisition& acq) {
  // A few shorthands to save us some typing.
  const auto& pwr = acq.pwr;
  const auto& tuned_freq = acq.tuned_freq;
  const auto& actual_samplerate = acq.actual_samplerate;

  // Print the header
  std::cout << "# rtl-power-fftw output" << std::endl;
  std::cout << "# Acquisition start: " << acq.startAcqTimestamp << std::endl;
  std::cout << "# Acquisition end: " << acq.endAcqTimestamp << std::endl;
  std::cout << "#" << std::endl;
  std::cout << "# frequency [Hz] power spectral density [dB/Hz]" << std::endl;

  // Calculate the precision needed for displaying the frequency.
  const int extraDigitsFreq = 2;
  const int significantPlacesFreq =
    ceil(floor(log10(tuned_freq)) - log10(actual_samplerate/params.N) + 1 + extraDigitsFreq);
  const int significantPlacesPwr = 6;

  for (int i = 0; i < params.N; i++) {
    double freq = tuned_freq + (i - params.N/2.0) * actual_samplerate / params.N;
    double pwrdb = 10*log10(pwr[i] / acq.repeatsProcessed / params.N / actual_samplerate)
                   - (params.baseline ? aux.baseline_values[i] : 0);
    std::cout << std::setprecision(significantPlacesFreq)
              << freq
              << " "
              << std::setprecision(significantPlacesPwr)
              << pwrdb
              << std::endl;
  }
  // Separate consecutive spectra with empty lines.
  std::cout << std::endl;
  std::cout.flush();
}


OutputWriter::OutputWriter(OutputStream* stream_) :
  stream(stream_)
{
  thread = std::thread(&OutputWriter::run, this);
}

OutputWriter::~OutputWriter() {
  queue.push_back(nullptr);
  thread.join();
}

void OutputWriter::run() {
  Acquisition* acquisition;
  while ((acquisition = queue.get()) != nullptr) {
    stream->write(*acquisition);
  }
}
