/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2015 and later by Bastien Chevreux
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

// Adapted from
//   http://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c

#ifndef _util_timer_h
#define _util_timer_h

#include <chrono>

class HRTimer
{
private:
  std::chrono::time_point<std::chrono::high_resolution_clock> TI_start;
  typedef decltype(TI_start-TI_start) diff_t;

public:
  inline HRTimer() : TI_start(std::chrono::high_resolution_clock::now()) {}
  inline void reset() { TI_start = std::chrono::high_resolution_clock::now(); }
  inline diff_t diff() {
    return std::chrono::high_resolution_clock::now() - TI_start;
  }
  static inline uint64 toNano(diff_t & df) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(df).count();
  }
  static inline uint64 toMicro(diff_t & df) {
    return std::chrono::duration_cast<std::chrono::microseconds>(df).count();
  }
};


class LRTimer
{
private:
  std::clock_t TI_start;
  typedef decltype(TI_start-TI_start) diff_t;

public:
  inline LRTimer() : TI_start(std::clock()) {}
  inline void reset() { TI_start = std::clock(); }
  inline diff_t diff() {
    return std::clock() - TI_start;
  }
  static inline uint64 toNano(diff_t & df) {
    return df*1000;
  }
  static inline uint64 toMicro(diff_t & df) {
    return df;
  }
};

// conventional low resolution timer
// Linux 2.6.31:
// ~2.5 - 3x faster than HRTimer
// ~4x faster than LRTimer
class CLRTimer
{
private:
  timeval TI_start;
  timeval tmptv;
  typedef suseconds_t diff_t;

public:
  inline CLRTimer() {reset();}
  inline void reset() { gettimeofday(&TI_start,nullptr); }
  inline diff_t diff() {
    gettimeofday(&tmptv,nullptr);
    if(likely(tmptv.tv_sec==TI_start.tv_sec)) {
      return (tmptv.tv_usec-TI_start.tv_usec);
    }
    time_t sec=tmptv.tv_sec-TI_start.tv_sec-1;
    return sec*1000000+(1000000-TI_start.tv_usec)+tmptv.tv_usec;
  }
  inline auto diffAndReset() {
    auto ret=diff();
    reset();
    return ret;
  }
  static inline uint64 toNano(diff_t & df) {
    return df*1000;
  }
  static inline uint64 toMicro(diff_t & df) {
    return df;
  }
};

typedef CLRTimer MIRATimer;

#endif
