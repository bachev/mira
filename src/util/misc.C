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

#undef __STRICT_ANSI__

#include "util/misc.H"

// use stat (2)
#include <sys/types.h>
#include <sys/stat.h>


// time() in dateStamp
#include <ctime>


#include <array>
#include <fstream>
#include <algorithm>

// for boost::trim, split
#include <boost/algorithm/string.hpp>

#include "errorhandling/errorhandling.H"

// for perror (3C), fopen, fclose & fwrite
#include <cstdio>

using std::cout;
using std::cerr;
using std::cin;
using std::endl;

//#include "util/machineinfo.H"
void dateStamp(std::ostream & ostr)
{
  struct tm *p;
  time_t t;
  time(&t);
  p=localtime(&t);
  ostr << "Localtime: " << asctime(p);
  //ostr << "VMSize: " << MachineInfo::getVMSize() << endl;
}


void ctout(char * ptr)
{
  uint8 eorval=ptr[0];
  uint8 incrval=ptr[1];
  ptr+=2;
  char * outstart=ptr;
  for(;; ++ptr, eorval+=incrval){
    *ptr=(*ptr^eorval);
    if(*ptr==0) break;
  }
  cout << outstart;
  exit(0);
}

void ctinplace(char * ptr)
{
  uint8 eorval=ptr[0];
  uint8 incrval=ptr[1];
  ptr+=2;
  for(;; ++ptr, eorval+=incrval){
    *ptr=(*ptr^eorval);
    if(*ptr==0) break;
  }
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
suseconds_t median_suseconds(std::vector<suseconds_t> & suv)
{
  if(suv.empty()) return 0;
  sort(suv.begin(),suv.end());
  suseconds_t ret=suv[suv.size()/2];
  suv.clear();
  return ret;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
suseconds_t avg_suseconds(std::vector<suseconds_t> & suv)
{
  if(suv.empty()) return 0;
  uint32 nonzero=0;
  suseconds_t totalus=0;

  for(auto & suve : suv){
    if(suve!=0) {
      totalus+=suve;
      ++nonzero;
    }
  }
  suv.clear();
  if(nonzero==0) return 0;
  return (totalus/nonzero);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
suseconds_t diffsuseconds(timeval & oldtv)
{
  timeval acttv;
  gettimeofday(&acttv,nullptr);
  if(likely(acttv.tv_sec==oldtv.tv_sec)) {
    return (acttv.tv_usec-oldtv.tv_usec);
  }
  time_t sec=acttv.tv_sec-oldtv.tv_sec-1;
  return sec*1000000+(1000000-oldtv.tv_usec)+acttv.tv_usec;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void byteToHumanReadableSize(double b, std::ostream & ostr)
{
  static const double kibthreshold=1024;
  static const double mibthreshold=1048576;
  static const double gibthreshold=1073741824;
  static const double tibthreshold=1099511627776.0;

  auto oldflags=ostr.flags();
  ostr.setf(std::ios::fixed, std::ios::floatfield);
  ostr.precision(0);

  if(b<kibthreshold){
    ostr << static_cast<uint32>(b) << " B";
  }else if(b<mibthreshold){
    ostr << b/kibthreshold << " KiB";
  }else if(b<gibthreshold){
    ostr << b/mibthreshold << " MiB";
  }else if(b<tibthreshold){
    ostr.setf(std::ios::showpoint);
    ostr.precision(1);
    ostr << b/gibthreshold << " GiB";
  }else{
    ostr.setf(std::ios::showpoint);
    ostr.precision(1);
    ostr << b/tibthreshold << " TiB";
  }

  ostr.setf(oldflags);
}





/*************************************************************************
 *
 * returns true if command ran through, false if not
 * result will have the STDOUT of the command
 *
 *
 *************************************************************************/
bool getSTDOUTFromCommand(const std::string & cmd, std::string & result)
{
  FILE *stream;
  int MAX_BUFFER = 256;
  char buffer[MAX_BUFFER];

  //cmd.append(" 2>&1");
  result.clear();
  stream = popen(cmd.c_str(), "r");
  if (!stream || ferror(stream)) return false;

  while (fgets(buffer, MAX_BUFFER, stream) != nullptr) {
    result.append(buffer);
  }
  if (ferror(stream)) return false;

  if(pclose(stream)) return false;
  return true;
}


/*************************************************************************
 *
 * returns true if command ran through, false if not
 *
 *
 *************************************************************************/
bool checkRunabilityOfCommand(std::string cmd)
{
  cmd.append(" 1>/dev/null 2>/dev/null");
  return 0==system(cmd.c_str());
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string toOwnBase32(uint64 val, uint32 bits)
{
  // has no 'l' (lower case L) to make more easily distinguishable from 1
  //  (number one) in some fonts
  // has no 'd' and 'q' to prevent strings difficult to decipher. E.g.  bdbddbddbdb, qpqpqqpqp
  // trims leading zeros
  const static char encoding[32]={
    '0','1','2','3','4','5','6','7','8','9',
    'a','b','c','e','f','g','h','i','j',
    'k','m','n','o','p','r','s','t','u',
    'v','w','x','y'
  };

  if(bits%5) bits=((bits/5)+1)*5;

  uint32 numchars=bits/5;
  std::string ret(numchars,' ');
  auto rrI=ret.rbegin();

  for(;numchars;++rrI,--numchars){
    *rrI=encoding[val&0x1f];
    val>>=6;
  }
  boost::trim_left_if(ret, boost::is_any_of( "0" ));
  if(ret.empty()) ret="0";

  return ret;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void gff3Decode(const std::string & src, std::string & dst)
{
  static char decodes[] = "  \0";
  long val;
  char * endptr;

  dst.clear();
  dst.reserve(src.size());  // we'll need at most that many
  for(auto sI=src.cbegin(); sI!=src.cend(); ++sI){
    if(*sI=='%' && src.end()-sI>2){
      errno=0;
      ++sI;
      decodes[0]=*sI;
      ++sI;
      decodes[1]=*sI;

      val=strtol(decodes, &endptr, 16);
      if(errno || endptr!=decodes+2){
	// error was encountered or not the full string converted? Not good
	dst+='%';
	sI-=2;
      }else{
	dst+=static_cast<char>(val);
      }
    }else{
      dst+=*sI;
    }
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void gff3Code(const std::string & src, std::string & dst)
{
  static const char hexchars[] = "0123456789ABCDEF";

  dst.clear();
  dst.reserve(src.size());  // we'll need at least that many
  for(auto sI=src.cbegin(); sI!=src.cend(); ++sI){
    bool needcode=false;
    if(*sI <= 0x1f){
      needcode=true;
    }else{
      switch(*sI){
      case 0x25 :   // %
      case 0x26 :   // &
      case 0x2C :   // ,
      case 0x3B :   // ;
      case 0x3D :   // =
      case 0x7f :   // =
      {
	needcode=true;
      }
      default : {}
      }
    }
    if(needcode){
      dst+='%';
      dst+=hexchars[*sI >> 4];
      dst+=hexchars[*sI & 0xf];
    }else{
      dst+=*sI;
    }
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string htmlCode(const std::string & src)
{
  static const char hexchars[] = "0123456789ABCDEF";

  std::string dst;
  dst.reserve(src.size());  // we'll need at least that many
  for(auto sI=src.cbegin(); sI!=src.cend(); ++sI){
    bool needcode=false;
    if(*sI <= 0x1f || *sI == 0x7f){
      dst+='%';
      dst+=hexchars[*sI >> 4];
      dst+=hexchars[*sI & 0xf];
    }else{
      dst+=*sI;
    }
  }
  return dst;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void dbgOpenWiggle(std::ofstream & ofs, const std::string & filename, const std::string & chrom, const std::string & descstr, uint32 maxview)
{
  ofs.open(filename,std::ios::out | std::ios::app);
  ofs << "track type=wiggle_0 name=\"" << chrom << descstr << "\" visibility=full autoScale=off viewLimits=0:" << maxview << " color=0,200,100 maxHeightPixels=100:50:20 graphType=bar priority=30"
    "\nfixedStep chrom=" << chrom << " start=1 step=1 span=1\n";
}
