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

using Buffer = std::vector<uint8_t>;

#ifdef FLOAT_FFT
using fft_datatype = float;
using fftw_cdt = fftwf_complex;
#define FFTW(name) fftwf_##name
#else
using fft_datatype = double;
using fftw_cdt = fftw_complex;
#define FFTW(name) fftw_##name
#endif // FLOAT_FFT

using complex = std::complex<fft_datatype>;

class Datastore {
  public:
    int N;
    int buffers;
    int64_t repeats;
    int64_t repeats_done = 0;

    std::mutex status_mutex;
    // Access to the following objects must be protected by locking
    // status_mutex.
    std::deque<Buffer*> empty_buffers;
    std::deque<Buffer*> occupied_buffers;
    bool acquisition_finished = false;
    std::condition_variable status_change;
    std::vector<int> queue_histogram;

    complex *inbuf, *outbuf;
    FFTW(plan) plan;
    std::vector<double> pwr;

    Datastore(int N, int buf_length, int64_t repeats, int buffers);
    ~Datastore();

    // Delete these so we don't accidentally mess anything up by copying
    // pointers to fftw_malloc'd buffers.
    Datastore(const Datastore&) = delete;
    Datastore(Datastore&&) = delete;
    Datastore& operator=(const Datastore&) = delete;
    Datastore& operator=(Datastore&&) = delete;
};

void fft(Datastore& data);

#endif // DATASTORE_H
