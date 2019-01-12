/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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


#include "io/fastq-mira.H"
#include "util/boostiostrutil.H"

//#include <boost/iostreams/filter/zlib.hpp>
//#include <boost/iostreams/filter/gzip.hpp>



const std::string FastQ::FQ_blanks=" \t";


// Plain vanilla constructor
FastQ::FastQ()
{
  FUNCSTART("FastQ::FastQ()");

  FQ_name.reserve(200);
  FQ_comment.reserve(200);
  FQ_seq.reserve(500);
  FQ_quals.reserve(500);
  FQ_dummy1.reserve(500);
  FQ_tracktellg=0;
  FQ_filesize=0;

  FUNCEND();
}

FastQ::~FastQ()
{
  FUNCSTART("FastQ::~FastQ()");

  discard();

  FUNCEND();
}


void FastQ::discard()
{
  FUNCSTART("FastQ::discard()");

  FUNCEND();
}




void FastQ::openFile(const std::string & filename, bool tracktellg)
{
  FUNCSTART("void FastQ::openFile(std::string & filename, bool tracktellg)");

  if(FQ_ifstream.is_open()){
    FQ_ifstream.close();
  }
  FQ_ifstream.clear();
  FQ_tracktellg=tracktellg;

  FQ_ifstream.open(filename,std::ios::in|std::ios::ate|std::ios::binary);
  if(!FQ_ifstream){
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }
  FQ_filesize=FQ_ifstream.tellg();
  if(FQ_filesize==0){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << filename);
  }
  FQ_ifstream.seekg(0, std::ios::beg);

  //FQ_fin.push(boost::iostreams::zlib_decompressor());
  //FQ_fin.push(boost::iostreams::gzip_decompressor());
  //FQ_fin.push(FQ_ifstream);
  setupStreamFilterByName(FQ_fin, FQ_ifstream, filename);

  FQ_linecount=0;
  FQ_tellgofread=0;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

int64_t FastQ::loadNext(std::streampos sp)
{
  FUNCSTART("int64_t FastQ::loadNext(std::streampos sp)");
  BUGIFTHROW(!FQ_ifstream.is_open(),"Tried to load a sequence while no file is opened.");
  FQ_fin.clear();
  FQ_fin.seekg(sp, std::ios::beg);
  FUNCEND();
  return loadNext();
}

/*************************************************************************
 *
 * returns size of read
 *  or -1 if end of file
 *
 *************************************************************************/

int64_t FastQ::loadNext()
{
  FUNCSTART("void FastQ::loadNext()");

  BUGIFTHROW(!FQ_ifstream.is_open(),"Tried to load a sequence while no file is opened.");

  FQ_name.clear();
  FQ_comment.clear();
  FQ_seq.clear();
  FQ_quals.clear();
  if(safeGetLineIntoName()){
    if(safeGetLineIntoSeq()){
      safeGetLineIntoQName();
      safeGetLineIntoQual();
    }
  }

  FUNCEND();
  if(FQ_name.empty() && FQ_fin.eof()) return -1;

  //cout << *this;
  return static_cast<int64>(FQ_seq.size());
}



bool FastQ::safeGetLineIntoName()
{
  FUNCSTART("void FastQ::safeGetLineIntName()");

  if(FQ_tracktellg) FQ_tellgofread=FQ_fin.tellg();

  if(FQ_fin.eof()) return false;
  if(static_cast<char>(FQ_fin.get())!='@'){
    if(!FQ_fin.eof()) MIRANOTIFY(Notify::FATAL,"First character of line " << FQ_linecount << " is not a '@'. That file is not in FASTQ format.");
  }
  if(FQ_fin.eof()) return false;
  while(getline(FQ_fin,FQ_name)){
    ++FQ_linecount;
    // get rid of '\r' from DOS
    while(!FQ_name.empty() && FQ_name[FQ_name.size()-1]=='\r') FQ_name.resize(FQ_name.size()-1);
    if(!FQ_name.empty()) break;
  }
  if(FQ_fin.eof()) return false;

  auto nameend=std::string::npos;
  nameend=FQ_name.find_first_of(FQ_blanks,0);
  if(nameend!=std::string::npos){
    FQ_comment=FQ_name.substr(nameend+1,FQ_name.size()-nameend);
    FQ_name.resize(nameend-1);
  }

  return true;
}

bool FastQ::safeGetLineIntoSeq()
{
  FUNCSTART("bool FastQ::safeGetLineIntoSeq()");

  while(getline(FQ_fin,FQ_seq)){
    ++FQ_linecount;
    // get rid of '\r' from DOS
    while(!FQ_seq.empty() && FQ_seq[FQ_seq.size()-1]=='\r') FQ_seq.resize(FQ_seq.size()-1);
    // allow for empty lines??
    if(!FQ_seq.empty()) break;
  }
  if(FQ_fin.eof()){
    MIRANOTIFY(Notify::FATAL,"At around line " << FQ_linecount << ": expected a sequence, none found.");
  }

  // continuation lines
  char peektmp=FQ_fin.peek();
  while(!FQ_fin.eof() && peektmp!='@' && peektmp!='+'){
    while(getline(FQ_fin,FQ_dummy1)){
      ++FQ_linecount;
      // get rid of '\r' from DOS
      while(!FQ_dummy1.empty() && FQ_dummy1[FQ_dummy1.size()-1]=='\r') FQ_dummy1.resize(FQ_dummy1.size()-1);
      // allow for empty lines??
      if(!FQ_dummy1.empty()) break;
    }
    FQ_seq.append(FQ_dummy1);
    peektmp=FQ_fin.peek();
  }
  if(FQ_fin.eof()){
    MIRANOTIFY(Notify::FATAL,"At around line " << FQ_linecount << ": expected a sequence, none found.");
  }

  return peektmp=='+';
}


void FastQ::safeGetLineIntoQName()
{
  FUNCSTART("void FastQ::safeGetLineIntQName()");

  if(static_cast<char>(FQ_fin.get())!='+'){
    MIRANOTIFY(Notify::FATAL,"First character of line " << FQ_linecount << " is not a '+'. We saw it moments ago, it should be ???");
  }
  if(static_cast<char>(FQ_fin.peek())=='\n'){
    ++FQ_linecount;
    FQ_fin.get();
  }else{
    while(getline(FQ_fin,FQ_dummy1)){
      ++FQ_linecount;
      // get rid of '\r' from DOS
      while(!FQ_dummy1.empty() && FQ_dummy1[FQ_dummy1.size()-1]=='\r') FQ_dummy1.resize(FQ_dummy1.size()-1);
      if(!FQ_dummy1.empty()) break;
    }

    if(FQ_dummy1.size()!=1){
      auto nameend=std::string::npos;
      nameend=FQ_dummy1.find_first_of(FQ_blanks,0);
      if(nameend!=std::string::npos){
	if(nameend){
	  FQ_dummy1.resize(nameend-1);
	}else{
	  FQ_dummy1.clear();
	}
      }
      if(FQ_dummy1!=FQ_name){
	MIRANOTIFY(Notify::FATAL,"Sequence name of read '" << FQ_name << "' is not equal to quality name '" << FQ_dummy1 << "'. Around line " << FQ_linecount);
      }
    }
  }

  return;
}


void FastQ::safeGetLineIntoQual()
{
  FUNCSTART("bool FastQ::safeGetLineIntoQual()");

  if(FQ_fin.eof()){
    MIRANOTIFY(Notify::FATAL,"At around line " << FQ_linecount << ": expected a sequence, none found.");
  }

  while(getline(FQ_fin,FQ_quals)){
    ++FQ_linecount;
    // get rid of '\r' from DOS
    while(!FQ_quals.empty() && FQ_quals[FQ_quals.size()-1]=='\r') FQ_quals.resize(FQ_quals.size()-1);
    // allow for empty lines??
    if(!FQ_quals.empty()) break;
  }
  if(FQ_fin.eof()){
    MIRANOTIFY(Notify::FATAL,"At around line " << FQ_linecount << ": expected a sequence, none found.");
  }

  // continuation lines
  while(!FQ_fin.eof() && FQ_quals.size() < FQ_seq.size()){
    while(getline(FQ_fin,FQ_dummy1)){
      ++FQ_linecount;
      // get rid of '\r' from DOS
      while(!FQ_dummy1.empty() && FQ_dummy1[FQ_dummy1.size()-1]=='\r') FQ_dummy1.resize(FQ_dummy1.size()-1);
      // allow for empty lines??
      if(!FQ_dummy1.empty()) break;
    }
    FQ_quals.append(FQ_dummy1);
  }

  if(FQ_seq.size()!=FQ_quals.size()){
    MIRANOTIFY(Notify::FATAL,"At around line " << FQ_linecount << ": sequence " << FQ_name << " has " << FQ_seq.size() << " bases, but it has " << FQ_quals.size() << " quality values ??? Here's what was read:\n" << *this);
  }

  return;
}


std::streampos FastQ::getStreamposOfRead() const {
  FUNCSTART("std::streampos FastQ::getStreamposOfRead()");
  BUGIFTHROW(!FQ_tracktellg,"Trying to get stream pos of a read, but file opened without tracking option");
  FUNCEND();
  return FQ_tellgofread;
}
