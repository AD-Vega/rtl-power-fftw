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

#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <exception>
#include <list>
#include <vector>
#include <rtl-sdr.h>

#include "params.h"
#include "datastore.h"
#include "device.h"

// This object contains data obtained from external sources.
class AuxData {
public:
  AuxData(const Params& params);

  // Values used for baseline correction.
  std::vector<double> baseline_values;
  // Values of the window function.
  std::vector<float> window_values;
};

// Plan of the measurement operation. Upon construction, this object adjusts
// several parameters according to the run time-data and creates a list of
// frequencies to tune to.
class Plan {
public:
  Plan(Params& params, int actual_samplerate);
  // Print a summary of the plan to stderr.
  void print() const;

  // List of frequencies to tune to.
  std::list<int> freqs_to_tune;
  // A cached value of the true sample rate.
  int actual_samplerate;

protected:
  Params& params;
};


// This exception gets thrown when the data acquisition procedure is unable to
// tune the hardware to the desired frequency despite several attempts. It is
// perfectly possible that this particular frequency is in a "dead" spot of the
// receiver, so TuneError must not be treated as a definite failure.
class TuneError : public std::exception {
public:
  TuneError(int freq_) : freq(freq_) {}
  const char* what() const noexcept {
    return "Could not tune to the given frequency.";
  }
  // The frequency that the receiver was unable to tune to.
  int frequency() { return freq; }

protected:
  int freq;
};


// This class performs data acquisition at a particular frequency.
class Acquisition {
public:
  Acquisition(const Params& params,
              AuxData& aux,
              Rtlsdr& rtldev,
              Datastore& data,
              int actual_samplerate,
              int freq);

  // Run the data acquisition.
  void run();
  // Print a summary of the acquisition (number of samples collected, number of
  // device readouts etc.) to stderr.
  void print_summary() const;
  // Write the gathered data to stdout.
  void write_data() const;

protected:
  // A helper function that returns the current date and time in the format
  // "YYYY-MM-DD HH:mm:ss UTC".
  static std::string currentDateTime();

  const Params& params;
  AuxData& aux;
  Rtlsdr& rtldev;
  Datastore& data;
  // A cached version of the actual sample rate.
  int actual_samplerate;
  // The frequency to tune to in this acquisition.
  int freq;
  // The actual frequency returned by the device.
  int tuned_freq;
  // Timestamp for the start of the acquisition.
  std::string startAcqTimestamp;
  // Timestamp for the end of the acquisition.
  std::string endAcqTimestamp;
  // Number of device readouts.
  int64_t deviceReadouts = 0;
  // Number of successful readouts (i.e., with no dropped samples).
  int64_t successfulReadouts = 0;
};

#endif // ACQUISITION_H
