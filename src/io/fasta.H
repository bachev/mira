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


/*
 * FASTA load routines
 */

#ifndef _bas_fasta_h_
#define _bas_fasta_h_

#include <iostream>
#include <fstream>
#include <string>

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"
#include "errorhandling/errorhandling.H"



class FASTA
{
  //Variables
private:
  uint32 FA_valid;

  std::string         FA_fastaseqname;
  std::string         FA_qualseqname;
  std::string         FA_sequence;
  std::string         FA_comment;
  std::vector<int32>  FA_intvalues;


public:


  //Functions
private:
  void init();
  void zeroVars();
  void adjustNameNCBIHack(std::string & name, std::string * comment);

public:
  FASTA();
  FASTA(FASTA const &other);
  ~FASTA();

  FASTA const & operator=(FASTA const & other);
  friend std::ostream & operator<<(std::ostream &ostr, FASTA const & fas);

  const std::string & getSeqName() const { return FA_fastaseqname;};
  const std::string & getQualName() const { return FA_qualseqname;};
  const std::string & getSequence() const { return FA_sequence;};
  const std::string & getComment() const { return FA_comment;};
  const std::vector<base_quality_t> getQualities() const {
    std::vector<base_quality_t> bla;
    bla.reserve(FA_sequence.size()+1);
    for(auto & iv : FA_intvalues) bla.push_back(static_cast<base_quality_t>(iv));
    return bla;
  };

  const std::vector<int32> & getINTValues() const { return FA_intvalues;};

  void discard();
  void load(const char * fastain);
  void load(const char * fastain, const char * qualin);
  void loadQual(const char * qualin);
  void loadINT(const char * qualin, int32 maxvalue);

  void loadNext(std::ifstream & fastafin, std::ifstream & qualin);
  void loadNextSeq(std::ifstream & fastafin);
  void loadNextINTSeq(std::ifstream & qualfin, int32 maxvalue);

  void testIfSeqAndQualMatch();
  bool testIfEmpty();

  void dumpSequence(std::ostream & fout);
  void dumpQuality(std::ostream & fout);
};


#endif
