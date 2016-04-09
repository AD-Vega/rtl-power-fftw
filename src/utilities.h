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

#ifndef UTILITIES_H
#define UTILITIES_H

#include <chrono>
#include <condition_variable>
#include <deque>

// A concurrent queue for an arbitrary data type. It is used as an interface for
// passing messages and data between the program threads.
template <typename T>
class ConcurrentQueue {
public:
  // Push an element to the back of the queue.
  void push_back(T item) {
    std::lock_guard<std::mutex> queueLock(mutex);
    queue.push_back(item);
    event.notify_one();
  }

  // Push an element to the front of the queue.
  void push_front(T item) {
    std::lock_guard<std::mutex> queueLock(mutex);
    queue.push_front(item);
    event.notify_one();
  }

  // Fetch an element from the front of the queue.
  T get() {
    std::unique_lock<std::mutex> queueLock(mutex);
    return getItem(queueLock);
  }

  // A variation on the get() method that also reports the initial queue length.
  T get(size_t& queueSize) {
    std::unique_lock<std::mutex> queueLock(mutex);
    queueSize = queue.size();
    return getItem(queueLock);
  }

protected:
  T getItem(std::unique_lock<std::mutex>& queueLock) {
    while (queue.empty())
      event.wait(queueLock);
    T item = queue.front();
    queue.pop_front();
    return item;
  }

  // A mutex that protects the access to the underlying deque.
  std::mutex mutex;
  // Used to notify the threads of a queue event (get, push).
  std::condition_variable event;
  // The actual queue.
  std::deque<T> queue;
};


// An abbreviation for the timestamp type.
using Timestamp = std::chrono::time_point<std::chrono::system_clock>;

// A helper function that returns the current date and time in the format
// "YYYY-MM-DD HH:mm:ss UTC".
std::string timeString(Timestamp& timestamp);

#endif // UTILITIES_H
