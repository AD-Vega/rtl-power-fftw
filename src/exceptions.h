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

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>

enum class ReturnValue {
  Success = 0,
  NoDeviceFound = 1,
  InvalidDeviceIndex = 2,
  InvalidArgument = 3,
  TCLAPerror = 4,
  InvalidInput = 5,
  AcquisitionError = 6
};

class RPFexception : public std::runtime_error {
public:
    explicit RPFexception(const std::string& what, ReturnValue retval_) :
      runtime_error(what), retval(retval_) {}
    ReturnValue returnValue() const { return retval; }

private:
  ReturnValue retval;
};

#endif // EXCEPTIONS_H
