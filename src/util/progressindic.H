/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2000 and later by Bastien Chevreux
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

#ifndef _util_progressindic_h_
#define _util_progressindic_h_

#include "stdinc/defines.H"
#include <iostream>

template <class T>
class ProgressIndicator
{
private:
  T PI_lastval;
  T PI_offset;
  T PI_rangemax;

  uint32 PI_delaycounter;
  uint32 PI_delaytrigger;

  int8 PI_printeddot;
//  double PI_dlast;
//  double PI_doffset;
//  double PI_drangemax;

public:
  ProgressIndicator() {
    setDelayTrigger(1);
    reset(0,1);
  }

  ProgressIndicator(T min, T max) {
    setDelayTrigger(1);
    reset(min,max);
  }

  ProgressIndicator(T min, T max, uint32 dtrigger) {
    setDelayTrigger(dtrigger);
    reset(min,max);
  }

  void setDelayTrigger(uint32 dtrigger){
    PI_delaytrigger=dtrigger;
    if(PI_delaytrigger==0) PI_delaytrigger=1;
    PI_delaycounter=PI_delaytrigger-1;
  }

  bool delaytrigger() {
    if(++PI_delaycounter == PI_delaytrigger){
      PI_delaycounter=0;
      return true;
    }
    return false;
  }

  void reset(T min, T max) {
    if(max <= min) {
      max=min;
      ++max;
    }
    PI_lastval=min;
    --PI_lastval;
    PI_offset=min;
    PI_rangemax=max-min;

    PI_printeddot=-1;
  }
  inline void increaseprogress(T step=1, std::ostream & ostr = std::cout) {
    progress(PI_lastval+step,ostr);
  };
  inline void finishAtOnce(std::ostream & ostr = std::cout) {
    progress(PI_rangemax,ostr);
  };

  void progress(T act, std::ostream & ostr = std::cout) {
    //std::cout << "act1 " << act << endl;

    if(act <= PI_lastval) return;

    act-=PI_offset;

    //std::cout << "act2 " << act << endl;
    //std::cout << "PI_offset " << PI_offset << endl;
    //std::cout << "PI_lastval " << PI_lastval << endl;

    double ptodod=static_cast<double>(100.0*static_cast<double>(act)/static_cast<double>(PI_rangemax));
    if(ptodod<static_cast<double>(0)) return;
    int8 ptodo;
    if(ptodod>static_cast<double>(100)) {
      ptodo=100;
    }else{
      ptodo=static_cast<int8>(static_cast<double>(100.0*static_cast<double>(act)/static_cast<double>(PI_rangemax)));
    }
    int8 pdone=PI_printeddot;

    //std::cout << "todo " << ptodo << endl;
    //std::cout << "done " << pdone << endl;

    if(pdone<ptodo){
      while(pdone<ptodo){
	pdone++;
	if(pdone%10 == 0){
	  ostr << " [" << static_cast<uint16>(pdone) << "%] ";
	} else if(pdone%5 == 0){
	  ostr << "|";
	}else{
	  ostr << ".";
	}
      }
      PI_printeddot=pdone;
      ostr.flush();
    }

    PI_lastval=act;
  }
};

#endif
