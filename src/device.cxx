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

#include <chrono>
#include <iostream>
#include <limits>
#include <string>
#include <thread>

#include "device.h"
#include "exceptions.h"


Rtlsdr::Rtlsdr(int dev_index) {
  int num_of_rtls = rtlsdr_get_device_count();
  if (num_of_rtls == 0) {
    throw RPFexception(
      "No RTL-SDR compatible devices found.",
      ReturnValue::NoDeviceFound);
  }

  if (dev_index >= num_of_rtls) {
    throw RPFexception(
      "Invalid RTL device number. Only " + std::to_string(num_of_rtls) + " devices available.",
      ReturnValue::InvalidDeviceIndex);
  }

  int rtl_retval = rtlsdr_open(&dev, (uint32_t)dev_index);
  if (rtl_retval < 0 ) {
    throw RPFexception(
      "Could not open rtl_sdr device " + std::to_string(dev_index),
      ReturnValue::HardwareError);
  }
}

Rtlsdr::~Rtlsdr() {
  rtlsdr_close(dev);
}

std::vector<int> Rtlsdr::gains() const {
  int number_of_gains = rtlsdr_get_tuner_gains(dev, NULL);
  if (number_of_gains <= 0) {
    throw RPFexception(
      "RTL device: could not read the number of available gains.",
      ReturnValue::HardwareError);
  }
  std::vector<int> gains(number_of_gains);
  if (rtlsdr_get_tuner_gains(dev, gains.data()) <= 0) {
    throw RPFexception(
      "RTL device: could not retrieve gain values.",
      ReturnValue::HardwareError);
  }
  return gains;
}

uint32_t Rtlsdr::sample_rate() const {
  uint32_t sample_rate = rtlsdr_get_sample_rate(dev);
  if (sample_rate == 0) {
    throw RPFexception(
      "RTL device: could not read sample rate.",
      ReturnValue::HardwareError);
  }
  return sample_rate;
}

uint32_t Rtlsdr::frequency() const {
  uint32_t frequency = rtlsdr_get_center_freq(dev);
  if (frequency == 0) {
    throw RPFexception(
      "RTL device: could not read frequency.",
      ReturnValue::HardwareError);
  }
  return frequency;
}

bool Rtlsdr::read(Buffer& buffer) const {
  int n_read;
  rtlsdr_reset_buffer(dev);
  rtlsdr_read_sync(dev, buffer.data(), buffer.size(), &n_read);
  return n_read == (signed)buffer.size();
}

void Rtlsdr::set_gain(int gain) {
  int status = 0;
  status += rtlsdr_set_tuner_gain_mode(dev, 1);
  status += rtlsdr_set_tuner_gain(dev, gain);

  if (status != 0) {
    throw RPFexception(
      "RTL device: could not set gain.",
      ReturnValue::HardwareError);
  }
}

void Rtlsdr::set_frequency(uint32_t frequency) {
  if (rtlsdr_set_center_freq(dev, frequency) < 0) {
    throw RPFexception(
      "RTL device: could not set center frequency.",
      ReturnValue::HardwareError);
  }
  // This sleeping is inherited from other code. There have been hints of strange
  // behaviour if it was commented out, so we left it in. If you actually know
  // why this would be necessary (or, to the contrary, that it is complete
  // bullshit), you are most welcome to explain it here!
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
}

void Rtlsdr::set_freq_correction(int ppm_error) {
  if (rtlsdr_set_freq_correction(dev, ppm_error) < 0) {
    throw RPFexception(
      "RTL device: could not set frequency correction.",
      ReturnValue::HardwareError);
  }
}

void Rtlsdr::set_sample_rate(uint32_t sample_rate) {
  if (rtlsdr_set_sample_rate(dev, sample_rate)) {
    throw RPFexception(
      "RTL device: could not set sample rate.",
      ReturnValue::HardwareError);
  }
}

int Rtlsdr::nearest_gain(int gain) const {
  int dif = std::numeric_limits<int>::max();
  int selected = 0;
  for (const auto& trial_gain : gains()) {
    int temp = abs(trial_gain - gain);
    if ( temp < dif ) {
      dif = temp;
      selected = trial_gain;
    }
  }
  return selected;
}

void Rtlsdr::print_gains() const {
  auto gain_table = gains();

  std::cerr << "Available gains (in 1/10th of dB): ";
  for (unsigned int i = 0; i < gain_table.size(); i++) {
    if (i != 0)
      std::cerr << ", ";
    std::cerr << gain_table[i];
  }
  std::cerr << std::endl;
}
