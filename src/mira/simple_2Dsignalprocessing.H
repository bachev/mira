/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2013 and later by Bastien Chevreux
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


#ifndef _mira_simplesignalproc_h_
#define _mira_simplesignalproc_h_


template <class PI>
void erode(PI begin, PI end, uint32 thickness)
{
  auto counter=thickness;
  for(; begin!=end; ++begin){
    if(*begin){
      if(counter){
	*begin=0;
	--counter;
      }
    }else{
      counter=thickness;
    }
  }
}

template <class PI>
void dilate(PI begin, PI end, uint32 thickness)
{
  uint32 counter=0;
  for(; begin!=end; ++begin){
    if(*begin){
      counter=thickness;
    }else{
      if(counter){
	*begin=1;
	--counter;
      }
    }
  }
}


template <class PI>
void removeSmallPeaks(PI & peakindicator, uint32 minlen)
{
  auto start = peakindicator.begin();
  while(start != peakindicator.end()){
    if(*start){
      auto end=start;
      while(end != peakindicator.end() && *end) ++end;
      //cout << "peak is " << end-start << endl;
      if(end-start < minlen){
	//cout << "delete\n";
	for(; start != end; ++start){
	  *start=0;
	}
      }
      start=end;
    }else{
      ++start;
    }
  }
}


template <class PI>
int64 distToFirstThreshold(PI begin, PI end, uint64 threshold)
{
  auto counter=0;
  for(; begin!=end; ++begin, ++counter){
    if(*begin>=threshold) return counter;
  }
  return -1;
}

#endif
