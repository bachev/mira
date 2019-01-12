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

// 	$Id$


#ifndef _exp_h_
#define _exp_h_

#include <cstdio>

#include <string>

#include <fstream>

#include "stdinc/stlincludes.H"

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "util/memusage.H"


#undef yyFlexLexer
#define yyFlexLexer EXPFlexLexer
#include <FlexLexer.h>


// class tag_t was a simple struct ... it has evolved :-/
// shouldn't be part of EXP file
// however, I don't know where to put elsewhere it right now

class tag_t {
public:
  uint32 from;                  // including
  uint32 to;                    // including
  char   strand;                // +, -, = means both
  std::string identifier;      /* Staden EXP doc says 4 byte id, but BNF
				   notation in TAG docs says string :( */
  std::string comment;

public:
  inline tag_t(){
    from=0;
    to=0;
    strand='?';
  };

  void serialiseComment(std::string & result) const;
  static void serialiseComment(const std::string & comment, std::string & result);

  bool extractGenBankKeyValueFromComment(const std::string & whatkex,
						std::string & result) const;
  static bool extractGenBankKeyValueFromComment(const std::string & comment,
						const std::string & whatkex,
						std::string & result);

  friend std::ostream & operator<<(std::ostream &ostr, const tag_t &t){
    ostr << "From: " << t.from <<  std::endl;
    ostr << "To: " << t.to << std::endl;
    ostr << "Strand: " << t.strand << std::endl;
    ostr << "Identifier: " << t.identifier << std::endl;
    ostr << "Comment: " << t.comment << std::endl;

    return ostr;
  }
};

// for Read
inline void emuoc_statsOnType(const tag_t & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{
  (void) c;
  (void) total_freecap;
  (void) total_lostba;

  total_size++;
  num_size=1;
  overhead=sizeof(c);
  bytes_size=overhead;
  iscontainer=true;

  lostbyalign=0;
  if(bytes_size%sizeof(void *)){
    lostbyalign=sizeof(void *)-bytes_size%sizeof(void *);
  }

  bytes_size+=estimateMemoryUsageOfContainer(c.identifier,false);
  bytes_size+=estimateMemoryUsageOfContainer(c.comment,false);


  return;
}



class EXP
{
private:
  // values taken from the exp-file, all strings 0 terminated.
  std::string EXP_bc;
  std::string EXP_cf;
  std::string EXP_cn;
  std::string EXP_cv;
  std::string EXP_dt;
  std::string EXP_en;
  std::string EXP_id;
  std::string EXP_le;
  std::string EXP_li;
  std::string EXP_ln;
  std::string EXP_lt;
  std::string EXP_ma;
  std::string EXP_mc;
  std::string EXP_mn;
  std::string EXP_mt;
  std::string EXP_ps;
  std::string EXP_sf;
  std::string EXP_sv;
  std::string EXP_ss;
  std::string EXP_tn;
  std::string EXP_op;
  std::string EXP_pn;

  std::string EXP_sq;

  int32 EXP_aq;
  int32 EXP_ch;
  int32 EXP_pr;
  int32 EXP_sc;
  int32 EXP_si_from;
  int32 EXP_si_to;
  int32 EXP_sp;
  /* left cutoffs: 10 meaning bases[0..9] cut , take [10..]
     right cutoffs: 100 meaning take bases [..99], cut [100..] */
  /* ergo: xL and xR make normal c delimiters like in
     for(int i=xL, i<xR, i++) ... */
  int32 EXP_cl;
  int32 EXP_cr;
  int32 EXP_ql;
  int32 EXP_qr;
  int32 EXP_sl;
  int32 EXP_sr;

  uint32 EXP_len_seq;

  std::vector<base_quality_t> EXP_qualities;
    //std::vector<int8> EXP_qualities;
  std::vector<tag_t> EXP_tags;

  std::vector<int32> EXP_on_ranges;


private:
  void   zeroVars();
  //char * conditionalStrCpy(const char * src);
  int32  gimmeAnInt(FlexLexer * lexer, const char * filename);
  void   getSequence(FlexLexer * lexer, const char * filename);

public:
  EXP();
  virtual ~EXP();

  void load(const char * filename, bool warning = false);
  void setSequenceFromFASTA(const char * seq);

  void discard();
  void dump();                           //dump an EXP to stdout
  //  int8 save(const char * filename);

  void swapTags(std::vector<tag_t> & tv) {
    EXP_tags.swap(tv);
  }

  inline int32 getSL() const { return EXP_sl; };
  inline int32 getSR() const { return EXP_sr; };
  inline int32 getQL() const { return EXP_ql; };
  inline int32 getQR() const { return EXP_qr; };
  inline int32 getCL() const { return EXP_cl; };
  inline int32 getCR() const { return EXP_cr; };
  inline int32 getAQ() const { return EXP_aq; };
  inline int32 getCH() const { return EXP_ch; };
  inline int32 getSC() const { return EXP_sc; };
  inline int32 getSIfrom() const { return EXP_si_from; };
  inline int32 getSIto() const { return EXP_si_to; };
  inline int32 getSP() const { return EXP_sp; };

  inline const std::string & getSequence() const { return EXP_sq; }
  inline const std::string & getSQ() const { return EXP_sq; }

  inline const std::vector<base_quality_t> & getQualities() const { return EXP_qualities;}
  inline const std::vector<base_quality_t> & getAV() const { return EXP_qualities;}
  inline const std::vector<int32> & getON() const { return EXP_on_ranges;}

  inline const std::string & getBC() const { return EXP_bc;}
  inline const std::string & getCF() const { return EXP_cf;}
  inline const std::string & getCN() const { return EXP_cn;}
  inline const std::string & getDT() const { return EXP_dt;}
  inline const std::string & getEN() const { return EXP_en;}
  inline const std::string & getID() const { return EXP_id;}
  inline const std::string & getLE() const { return EXP_le;}
  inline const std::string & getLI() const { return EXP_li;}
  inline const std::string & getLN() const { return EXP_ln;}
  inline const std::string & getLT() const { return EXP_lt;}
  inline const std::string & getMA() const { return EXP_ma;}
  inline const std::string & getMC() const { return EXP_mc;}
  inline const std::string & getMN() const { return EXP_mn;}
  inline const std::string & getMT() const { return EXP_mt;}
  inline const std::string & getPN() const { return EXP_pn;}
  inline const std::string & getPS() const { return EXP_ps;}
  inline const std::string & getSF() const { return EXP_sf;}
  inline const std::string & getSV() const { return EXP_sv;}
  inline const std::string & getTN() const { return EXP_tn;}
};

#endif
