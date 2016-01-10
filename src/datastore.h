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

#ifndef DATASTORE_H
#define DATASTORE_H

#include <complex>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <vector>
#include <fftw3.h>

#include "params.h"

using Buffer = std::vector<uint8_t>;
using complex = std::complex<float>;

class Datastore {
  public:
    const Params& params;
    int64_t repeats_done = 0;

    std::mutex status_mutex;
    // Access to the following objects must be protected by locking
    // status_mutex.
    std::deque<Buffer*> empty_buffers;
    std::deque<Buffer*> occupied_buffers;
    bool acquisition_finished = false;
    std::condition_variable status_change;
    std::vector<int> queue_histogram;

    std::vector<float>& window_values;

    complex *inbuf, *outbuf;
    fftwf_plan plan;
    std::vector<double> pwr;

    Datastore(const Params& params, std::vector<float>& window_values);
    ~Datastore();

    // Delete these so we don't accidentally mess anything up by copying
    // pointers to fftw_malloc'd buffers.
    Datastore(const Datastore&) = delete;
    Datastore(Datastore&&) = delete;
    Datastore& operator=(const Datastore&) = delete;
    Datastore& operator=(Datastore&&) = delete;

    // This function will be started in a separate thread.
    void fftThread();
    void printQueueHistogram() const;
};

#endif // DATASTORE_H
