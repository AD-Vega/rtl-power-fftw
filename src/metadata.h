/*
* rtl_power_fftw, program for calculating power spectrum from rtl-sdr reciever.
* Copyright (C) 2015 Klemen Blokar <klemen.blokar@ad-vega.si>
*                    Andrej Lajovic <andrej.lajovic@ad-vega.si>
*
*                    Additions by: Mario Cannistra <mariocannistra@gmail.com>
*                                 (added this file to manage various metadata)
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

#ifndef METADATA_H
#define METADATA_H

#include <time.h>

extern int metaCols, metaRows;
extern int startFreq, endFreq, stepFreq;
extern time_t scanEnd, scanBeg;
extern float avgScanDur, sumScanDur;
extern std::string firstAcqTimestamp, lastAcqTimestamp;
extern int cntTimeStamps;

#endif // METADATA_H
