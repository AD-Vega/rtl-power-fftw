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
#include <rtl-sdr.h>

#include "acquisition.h"
#include "dispatcher.h"
#include "device.h"
#include "exceptions.h"
#include "interrupts.h"
#include "params.h"


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

    // Print the available gains and select the one nearest to the requested gain.
    rtldev.print_gains();
    int gain = rtldev.nearest_gain(params.gain);
    std::cerr << "Selected nearest available gain: " << gain
              << " (" << 0.1*gain << " dB)" << std::endl;
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
      std::cerr << "PPM error set to: " << params.ppm_error << std::endl;
    }

    // Set sample rate
    rtldev.set_sample_rate(params.sample_rate);
    int actual_samplerate = rtldev.sample_rate();
    std::cerr << "Actual sample rate: " << actual_samplerate << " Hz" << std::endl;

    // Create a plan of the operation. This will calculate the number of repeats,
    // adjust the buffer size for optimal performance and create a list of frequency
    // hops.
    Plan plan(params, actual_samplerate);
    // Print info on capture time and associated specifics.
    plan.print();

    Dispatcher dispatcher(params, auxData, 4);

    // Install a signal handler for detecting Ctrl+C.
    set_CtrlC_handler(true);

    //Read from device and do FFT
    do {
      for (auto iter = plan.freqs_to_tune.begin(); iter != plan.freqs_to_tune.end();) {
        // Begin a new data acquisition.
        Acquisition acquisition(params, auxData, rtldev, dispatcher, actual_samplerate, *iter);
        try {
          // Read the required amount of data and process it.
          acquisition.run();
          iter++;
        }
        catch (TuneError &e) {
          // The receiver was unable to tune to this frequency. It might be just a "dead"
          // spot of the receiver. Remove this frequency from the list and continue.
          std::cerr << "Unable to tune to " << e.frequency() << ". Dropping from frequency list." << std::endl;
          iter = plan.freqs_to_tune.erase(iter);
          continue;
        }

        acquisition.waitForResultsReady();
        // Print a summary (number of samples, readouts etc.) to stderr.
        acquisition.print_summary();
        // Write the gathered data to stdout.
        acquisition.write_data();
        // Print the histogram of the queue length to stderr.
        acquisition.printQueueHistogram();

        // Check for interrupts.
        if (checkInterrupt(InterruptState::FinishNow))
          break;
      }
      // Mark the end of a measurement set with an additional empty line
      // (one was already output as a terminator for the last data set).
      std::cout << std::endl;

      // Loop if continuous mode is set, but quit if the frequency list becomes
      // empty or if the user interrupts the operation.
    } while (params.endless && plan.freqs_to_tune.size()
             && !checkInterrupt(InterruptState::FinishPass));

    if (plan.freqs_to_tune.size() == 0) {
      // No valid frequencies left in the list. This is certainly not OK.
      throw RPFexception("No valid frequencies left.", ReturnValue::AcquisitionError);
    }
  }
  catch (RPFexception &exception) {
    std::cerr << exception.what() << std::endl;
    final_retval = exception.returnValue();
  }

  return (int)final_retval;
}
