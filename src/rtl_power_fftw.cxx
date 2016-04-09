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

#include "acquisition.h"
#include "device.h"
#include "dispatcher.h"
#include "exceptions.h"
#include "interrupts.h"
#include "output.h"
#include "params.h"

#include <chrono>
#include <iostream>
#include <memory>
#include <rtl-sdr.h>
#include <thread>


int main(int argc, char **argv)
{
  ReturnValue final_retval = ReturnValue::Success;

  try {
    // Parse command line arguments.
    Params params(argc, argv);
    // Read auxiliary data for window function and baseline correction.
    AuxData auxData(params);

    // Set up RTL-SDR device.
    Rtlsdr rtldev(params.dev_index);

    // Print the available gains.
    auto gain_table = rtldev.gains();
    Diagnostics() << "Available gains (in 1/10th of dB): ";
    for (unsigned int i = 0; i < gain_table.size(); i++) {
      if (i != 0)
        Diagnostics() << ", ";
      Diagnostics() << gain_table[i];
    }
    Diagnostics() << "\n";

    // Select the gain nearest to the requested gain.
    int gain = rtldev.nearest_gain(params.gain);
    Diagnostics() << "Selected nearest available gain: " << gain
                  << " (" << 0.1*gain << " dB)\n";
    rtldev.set_gain(gain);

    // Temporarily set the frequency to params.cfreq, just so that the device does not
    // complain upon setting the sample rate. If this fails, it's not a big deal:
    // the sample rate will be read out just fine, but librtlsdr _might_ print an
    // error message to stderr.
    try {
      rtldev.set_frequency(params.cfreq);
    }
    catch (RPFexception&) {}

    // Set frequency correction
    if (params.ppm_error != 0) {
      rtldev.set_freq_correction(params.ppm_error);
      Diagnostics() << "PPM error set to: " << params.ppm_error << "\n";
    }

    // Set sample rate
    rtldev.set_sample_rate(params.sample_rate);
    int actual_samplerate = rtldev.sample_rate();
    Diagnostics() << "Actual sample rate: " << actual_samplerate << " Hz\n";

    // Create a plan of the operation. This will calculate the number of repeats,
    // adjust the buffer size for optimal performance and create a list of frequency
    // hops.
    Plan plan(params, actual_samplerate);
    // Print info on capture time and associated specifics.
    plan.print();

    Dispatcher dispatcher(params, auxData);
    TextStream stream(params, auxData);
    OutputWriter writer(&stream);

    // For limited session (scanning) duration.
    std::chrono::time_point<std::chrono::steady_clock> scanUntil =
      std::chrono::steady_clock::now() + std::chrono::seconds(params.session_duration);

    // Install a signal handler for detecting Ctrl+C.
    set_CtrlC_handler(true);

    //Read from device and do FFT
    do {
      writer.datasetBegin();
      for (auto iter = plan.freqs_to_tune.begin(); iter != plan.freqs_to_tune.end();) {
        // Begin a new data acquisition.
        auto acquisition = std::make_shared<Acquisition>(
          params, auxData, rtldev, dispatcher, actual_samplerate, *iter);
        try {
          // Read the required amount of data and process it.
          acquisition->run();
          iter++;
        }
        catch (TuneError &e) {
          // The receiver was unable to tune to this frequency. It might be just a "dead"
          // spot of the receiver. Remove this frequency from the list and continue.
          Diagnostics(LogLevel::Warning)
            << "Unable to tune to " << e.frequency() << ". Dropping from frequency list.\n";
          iter = plan.freqs_to_tune.erase(iter);
          continue;
        }

        // Pass the acquisition to the writer. The data will be output when it is ready.
        writer.queueData(acquisition);

        // Check for interrupts.
        if (checkInterrupt(InterruptState::FinishNow))
          break;
      }
      // Mark the end of a measurement set.
      writer.datasetEnd();

      // Loop if continuous mode or session duration is set, but quit if the frequency
      // list becomes empty, if the user interrupts the operation or if the session time
      // limit is exceeded.
    } while ((params.endless
              || (params.session_duration >= 0
                  && std::chrono::steady_clock::now() < scanUntil))
             && plan.freqs_to_tune.size()
             && !checkInterrupt(InterruptState::FinishPass)
            );

    if (plan.freqs_to_tune.size() == 0) {
      // No valid frequencies left in the list. This is certainly not OK.
      throw RPFexception("No valid frequencies left.", ReturnValue::AcquisitionError);
    }
  }
  catch (RPFexception &exception) {
    Diagnostics(LogLevel::Error) << exception.what() << "\n";
    final_retval = exception.returnValue();
  }

  return (int)final_retval;
}
