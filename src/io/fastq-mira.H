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

#ifndef _bas_fastq_mira_h_
#define _bas_fastq_mira_h_

#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"



class FastQ
{
  //Variables
private:
  static const std::string FQ_blanks;

  boost::iostreams::filtering_istream FQ_fin;
  std::ifstream    FQ_ifstream;
  std::string      FQ_linebuffer;
  size_t           FQ_linecount;
  std::streampos   FQ_filesize;
  std::streampos   FQ_tellgofread;
  bool             FQ_tracktellg;

  std::string FQ_name;
  std::string FQ_comment;
  std::string FQ_seq;
  std::string FQ_quals;

  // for splitting quality name line and comparing to FQ_name
  std::string FQ_dummy1;

public:


  //Functions
private:
  void foolCompiler();

  void splitNameLine(std::string & name,std::string & comment);
  size_t safeGetLineIntoLineBuffer();
  bool safeGetLineIntoName();
  void safeGetLineIntoQName();
  bool safeGetLineIntoSeq();
  void safeGetLineIntoQual();

public:
  FastQ();
  ~FastQ();

  FastQ(FastQ const &other) = delete;
  FastQ const & operator=(FastQ const & other) = delete;
  friend std::ostream & operator<<(std::ostream &ostr, FastQ const &fq){
    ostr << '@' << fq.FQ_name;
    if(!fq.FQ_comment.empty()) ostr << ' ' << fq.FQ_comment << '\n';
    if(!fq.FQ_quals.empty()){
      ostr << "+\n" << fq.FQ_quals << '\n';
    }
    return ostr;
  }

  void openFile(const std::string & filename, bool tracktellg=false);
  int64 loadNext();
  int64 loadNext(std::streampos sp);

  inline const std::string & getName() const {return FQ_name;}
  inline const std::string & getComment() const {return FQ_comment;}
  inline const std::string & getSequence() const {return FQ_seq;}
  inline const std::string & getQuals() const {return FQ_quals;}
  std::streampos getStreamposOfRead() const;
  void setStreampos(std::streampos sp) const;
  inline std::streampos getFileSize() const {return FQ_filesize;}
  inline size_t getLineCount() const {return FQ_linecount;}

  void discard();
};


#endif
