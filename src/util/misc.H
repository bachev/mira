/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
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


#ifndef _util_misc_h
#define _util_misc_h


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <sys/time.h>

#include "stdinc/defines.H"

#include "stdinc/stlincludes.H"




//char * conditionalStrCpy(const char * src);

void dateStamp(std::ostream & ostr);

void ctout(char * ptr);
void ctinplace(char * ptr);

suseconds_t median_suseconds(std::vector<suseconds_t> & suv);
suseconds_t avg_suseconds(std::vector<suseconds_t> & suv);
suseconds_t diffsuseconds(timeval & oldtv);


void byteToHumanReadableSize(double b, std::ostream & o);

bool getSTDOUTFromCommand(const std::string & cmd, std::string & result);
bool checkRunabilityOfCommand(std::string cmd);

std::string toOwnBase32(uint64 val, uint32 bits);

void gff3Decode(const std::string & src, std::string & dst);
void gff3Code(const std::string & src, std::string & dst);
std::string htmlCode(const std::string & src);

inline const std::string & replaceEmptyString(const std::string & text){
  static const std::string string_minus("-");
  if(text.empty()) return string_minus;
  return text;
}


void dbgOpenWiggle(std::ofstream & ofs, const std::string & filename, const std::string & chrom, const std::string & descstr, uint32 maxview);
template <class CC> void dbgContainerToWiggle(CC & container,const std::string & chrom, const std::string & desc){
  if(container.empty()) return;
  auto maxv=container.front();
  for(auto & ce : container) maxv=std::max(maxv,ce);
  std::ofstream ofs;
  dbgOpenWiggle(ofs,desc+".wig",chrom,desc,maxv);
  for(auto & ce : container) ofs << static_cast<int64>(ce) << '\n';
};



#endif
