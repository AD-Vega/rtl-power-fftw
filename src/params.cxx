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

#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <tclap/CmdLine.h>

#include "params.h"

int parse_frequency(std::string s) {
  std::istringstream ss(s);
  double f;
  std::string multiplier;
  ss >> f >> multiplier;
  if (multiplier == "k")
    f *= 1e3;
  else if (multiplier == "M")
    f *= 1e6;
  else if (multiplier == "G")
    f *= 1e9;
  else if (multiplier != "")
    return -1;
  return (int)f;
}

double parse_time(std::string s) {
  // The string should end with an unit. If no unit is present, we assume that
  // the last number is in seconds.
  std::string permitted_units = "dhms";
  if (permitted_units.find(s.back()) == std::string::npos)
    s.push_back('s');

  std::stringstream ss(s);
  double value;
  char unit;
  double t = 0;
  // We'll use these to prevent the same unit from being used twice (as it is
  // most likely a user mistake).
  bool days_consumed, hours_consumed, minutes_consumed, seconds_consumed;
  days_consumed = hours_consumed = minutes_consumed = seconds_consumed = false;

  while (ss >> value && ss.get(unit)) {
    if (unit == 'd' && !days_consumed) {
      t += value*86400;
      days_consumed = true;
    }
    else if (unit == 'h' && !hours_consumed) {
      t += value*3600;
      hours_consumed = true;
    }
    else if (unit == 'm' && !minutes_consumed) {
      t += value*60;
      minutes_consumed = true;
    }
    else if (unit == 's' && !seconds_consumed) {
      t += value;
      seconds_consumed = true;
    }
    else
      return -1;
  }

  if (ss.eof())
    return t;
  else {
    // Unconsumed input left - this indicates a parse error.
    return -1;
  }
}

class NegativeArgException {
public:
  NegativeArgException(std::string msg_) : msg(msg_) {}
  std::string what() const { return msg; }
private:
  std::string msg;
};

template <typename T>
void ensure_positive_arg(std::list<TCLAP::ValueArg<T>*> list) {
  for (auto arg : list) {
    if (arg->isSet() && arg->getValue() < 0) {
      std::ostringstream message;
      message << "Argument to '" << arg->getName() << "' must be a positive number.";
      throw NegativeArgException(message.str());
    }
  }
}

ReturnValue Params::parse(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Obtain power spectrum from RTL device using FFTW library.", ' ', "0.1");
    TCLAP::ValueArg<int> arg_buffers("","buffers","Number of read buffers (don't touch unless running out of memory).",false,buffers,"buffers");
    cmd.add( arg_buffers );
    TCLAP::ValueArg<std::string> arg_integration_time("t","time","Integration time (incompatible with -n).",false,"","seconds");
    cmd.add( arg_integration_time );
    TCLAP::SwitchArg arg_strict_time("T","strict-time","End measurement when the time set with --time option is up, regardless of gathered samples.",strict_time);
    cmd.add( arg_strict_time );
    TCLAP::ValueArg<int> arg_bufferlen("s","buffer-size","Size of read buffers (leave it unless you know what you are doing).", false, buf_length, "bytes");
    cmd.add( arg_bufferlen );
    TCLAP::ValueArg<int> arg_rate("r","rate","Sample rate of the receiver.",false,sample_rate,"samples/s");
    cmd.add( arg_rate );
    TCLAP::ValueArg<int> arg_ppm("p","ppm","Set custom ppm error in RTL-SDR device.", false, ppm_error, "ppm");
    cmd.add( arg_ppm );
    TCLAP::ValueArg<int64_t> arg_repeats("n","repeats","Number of scans for averaging (incompatible with -t).",false,repeats,"repeats");
    cmd.add( arg_repeats );
    TCLAP::ValueArg<int> arg_gain("g","gain","Receiver gain.",false, gain, "1/10th of dB");
    cmd.add( arg_gain );
    TCLAP::ValueArg<std::string> arg_freq("f","freq","Center frequency of the receiver or frequency range to scan.",false,"","Hz | Hz:Hz");
    cmd.add( arg_freq );
    TCLAP::ValueArg<int> arg_index("d","device","RTL-SDR device index.",false,dev_index,"device index");
    cmd.add( arg_index );
    TCLAP::SwitchArg arg_continue("c","continue","Repeat the same measurement endlessly.", endless);
    cmd.add( arg_continue );
    TCLAP::ValueArg<int> arg_bins("b","bins","Number of bins in FFT spectrum (must be even number)",false,N,"bins in FFT spectrum");
    cmd.add( arg_bins );
    TCLAP::ValueArg<std::string> arg_baseline("B","baseline","Subtract baseline, read baseline data from file or stdin.",false,"","file|-");
    cmd.add( arg_baseline );

    cmd.parse(argc, argv);

    try {
      // Ain't this C++11 f**** magic? Watch this:
      ensure_positive_arg<int>({&arg_bins, &arg_rate, &arg_gain, &arg_index, &arg_buffers, &arg_bufferlen});
      ensure_positive_arg<int64_t>({&arg_repeats});
    }
    catch (NegativeArgException& e) {
      std::cerr << e.what() << std::endl;
      return ReturnValue::InvalidArgument;
    }

    dev_index = arg_index.getValue();
    N = arg_bins.getValue();
    //Number of bins should be even, to allow us a neat trick to get fftw output properly aligned.
    if (N % 2 != 0) {
      N++;
      std::cerr << "Number of bins should be even, changing to " << N << "." << std::endl;
    }
    gain = arg_gain.getValue();
    sample_rate = arg_rate.getValue();
    buffers = arg_buffers.getValue();
    buf_length = arg_bufferlen.getValue();
    endless = arg_continue.getValue();
    strict_time = arg_strict_time.getValue();
    //clipped_output_isSet = arg_clipped.getValue();

    // Due to USB specifics, buffer length for reading rtl_sdr device
    // must be a multiple of 16384. We have to keep it that way.
    // For performance reasons, the actual buffer length should be in the
    // MB range.
    if (buf_length % base_buf != 0) {
      buf_length = floor((double)buf_length/base_buf + 0.5) * base_buf;
      std::cerr << "Buffer length should be multiple of " << base_buf
                << ", changing to " << buf_length << "." << std::endl;
    }
    ppm_error = arg_ppm.getValue();
    if (arg_freq.isSet()) {
      std::string a_freq = arg_freq.getValue();
      std::size_t colon_position = a_freq.find(":");
      std::istringstream opt(a_freq);
      if (colon_position != std::string::npos) {
        std::string startFreqString, stopFreqString;
        if (getline(opt, startFreqString, ':') && getline(opt, stopFreqString)) {
          startfreq = parse_frequency(startFreqString);
          stopfreq = parse_frequency(stopFreqString);
          if (startfreq < 0 || stopfreq < 0 || stopfreq < startfreq) {
            std::cerr << "Invalid frequency range given to --freq: " 
                    << startfreq << ":" << stopfreq << ". "
                    << "Expecting positive numbers in ascending order, allowing the k,M,G multipliers. Exiting."
                    << std::endl;
            return ReturnValue::InvalidArgument;
          }
          else {
            freq_hopping_isSet = true;
            cfreq = (startfreq + stopfreq)/2;
          }
        }
        else {
          std::cerr << "Could not parse frequency range given to --freq: " 
                    << opt 
                    << ". Expecting form startfreq:stopfreq. Exiting." 
                    << std::endl;
          return ReturnValue::InvalidArgument;
        }
      }
      else {
        cfreq = parse_frequency(a_freq);
        if (cfreq < 0) {
          std::cerr << "Invalid frequency given to --freq: " 
                    << cfreq
                    << ". Expecting a positive number, allowing the k,M,G multipliers. Exiting."
                    << std::endl;
          return ReturnValue::InvalidArgument;
        }
      }
    }
    if (arg_repeats.isSet())
      repeats = arg_repeats.getValue();
    else
      repeats = buf_length/(2*N);
    if (arg_integration_time.isSet()) {
      integration_time = parse_time(arg_integration_time.getValue());
      if (integration_time <= 0) {
        std::cerr << "Could not parse the value given to --time. "
                  << "Expecting format [WdXhYm]Z[s]. Exiting."
                  << std::endl;
        return ReturnValue::InvalidArgument;
      }
      integration_time_isSet = true;
    }
    //Integration time
    if (arg_integration_time.isSet() + arg_repeats.isSet() > 1) {
      std::cerr << "Options -n and -t are mutually exclusive. Exiting." << std::endl;
      return ReturnValue::InvalidArgument;
    }
    if (arg_strict_time.isSet() && !arg_integration_time.isSet()) {
      std::cerr << "Warning: option --strict-time has no effect without --time." << std::endl;
      strict_time = false;
    }
    //Optimally adjust buffer length for small sample sizes only if buffer length is not user defined.
    if (arg_bufferlen.isSet()) {
      buf_length_isSet = true;
    }
    baseline = arg_baseline.isSet();
    if (baseline)
      baseline_file = arg_baseline.getValue();
  }
  catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return ReturnValue::TCLAPerror;
  }

  return ReturnValue::Success;
}
