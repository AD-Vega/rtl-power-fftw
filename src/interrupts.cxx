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
#include <signal.h>

#include "interrupts.h"

std::atomic<int> interrupts(0);

// Signal handler that will get invoked on Ctrl+C
void CtrlC_handler(int signal) {
  if (++interrupts == static_cast<int>(InterruptState::FinishNow))
    set_CtrlC_handler(false);
}

void set_CtrlC_handler(bool install) {
#ifdef _WIN32
  signal(SIGINT, CtrlC_handler);
#else
  struct sigaction action;
  action.sa_handler = (install ? &CtrlC_handler : SIG_DFL);
  sigemptyset(&action.sa_mask);
  action.sa_flags = 0;
  sigaction(SIGINT, &action, nullptr);
#endif
}

bool checkInterrupt(InterruptState checkLevel) {
  static int reportedInterrupts = 0;
  const int currentInterrupts = interrupts;

  while (reportedInterrupts < currentInterrupts) {
    reportedInterrupts++;
    auto state = static_cast<InterruptState>(reportedInterrupts);

    if (state == InterruptState::FinishPass)
      std::cerr << "Interrupted, will try to finish this pass." << std::endl;
    else if (state == InterruptState::FinishNow)
      std::cerr << "Interrupted, finishing now." << std::endl;
  }

  return currentInterrupts >= static_cast<int>(checkLevel);
}
