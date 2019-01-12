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

#include "io/phd.H"

using std::endl;

// Plain vanilla constructor
PHD::PHD()
{
  FUNCSTART("PHD::PHD()");

  FUNCEND();
}

PHD::~PHD()
{
  FUNCSTART("PHD::~PHD()");

  FUNCEND();
}

void PHD::discard()
{
  FUNCSTART("PHD::discard()");

  PHD_name="";
  PHD_sequence="";
  PHD_qualities.clear();
  PHD_peakindex.clear();

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//PHD::PHD(PHD const &other)
//{
//  FUNCSTART("PHD::PHD(PHD const &other)");
//
//  PHD_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, also needed by copy-constructor
//PHD const & PHD::operator=(PHD const & other)
//{
//  FUNCSTART("PHD const & PHD::operator=(PHD const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

std::ostream & operator<<(std::ostream &ostr, PHD const &thephd)
{
  //  FUNCSTART("friend std::ostream & PHD::operator<<(std::ostream &ostr, const  &???)");

  ostr << "BEGIN_SEQUENCE x\n\n";
  ostr << "BEGIN_DNA\n";
  for(uint32 i=0; i<thephd.PHD_sequence.size(); i++) {
    ostr << thephd.PHD_sequence[i] << " " << static_cast<uint16> (thephd.PHD_qualities[i]) << " " << thephd.PHD_peakindex[i] << endl;
  }

  ostr << "END_DNA\n\n";
  ostr << "END_SEQUENCE\n";

  //  FUNCEND();
  return ostr;
}


bool PHD::testIfEmpty()
{
  if(PHD_name.size()) return false;
  if(PHD_sequence.size()) return false;
  if(PHD_qualities.size()) return false;
  if(PHD_peakindex.size()) return false;

  return true;
}

void PHD::load(const char * filename)
{
  FUNCSTART("void PHD::load(const char * filename)");

  discard();

  std::ifstream fin(filename, std::ios::in|std::ios::ate);

  if(!fin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  if(!fin.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << filename);
  }
  fin.seekg(0, std::ios::beg);

  loadNextSeq(fin);

  FUNCEND();
}

void PHD::loadNextSeq(std::ifstream & phdfin)
{
  FUNCSTART("void PHD::loadNextSeq(std::ifstream & phdfin)");

  discard();

  std::string tmpstr;

  bool found=false;
  while(!phdfin.eof()){
    phdfin >> tmpstr;
    if(tmpstr.compare("BEGIN_SEQUENCE") == 0){
      found=true;
      phdfin >> PHD_name;
      break;
    }
  }
  if(!found) {
    MIRANOTIFY(Notify::FATAL, "File does not seem to be a PHD file?");
  }

  found=false;
  while(!phdfin.eof()){
    phdfin >> tmpstr;
    if(tmpstr.compare("BEGIN_DNA") == 0) {
	found=true;
	break;
    }
  }
  if(!found) {
    MIRANOTIFY(Notify::FATAL, "Missing BEGIN_DNA for sequence " << PHD_name);
  }

  //char base;
  uint16 quality;
  uint32 peakindex;
  found=false;
  while(!phdfin.eof()){
    phdfin >> tmpstr;
    if(tmpstr.compare("END_DNA") == 0) {
      found=true;
      break;
    }

    phdfin >> quality;
    phdfin >> peakindex;
    PHD_sequence+=static_cast<char>(toupper(tmpstr[0]));
    PHD_qualities.push_back(static_cast<uint8>(quality));
    PHD_peakindex.push_back(peakindex);
  }

  if(!found) {
    MIRANOTIFY(Notify::FATAL, "Missing END_DNA for sequence " << PHD_name);
  }

  FUNCEND();
}





#if 0

BEGIN_SEQUENCE pf84c05.s1
BEGIN_COMMENT
CHROMAT_FILE: pf84c05.s1
ABI_THUMBPRINT: 0
PHRED_VERSION: 0.980904.a
CALL_METHOD: phred
QUALITY_LEVELS: 99
TIME: Fri Oct 15 08:07:40 1999
TRACE_ARRAY_MIN_INDEX: 0
TRACE_ARRAY_MAX_INDEX: 9915
CHEM: prim
DYE: ET
END_COMMENT
BEGIN_DNA
c 6 5
g 6 11
c 8 26

#endif
