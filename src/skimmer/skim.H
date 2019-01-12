/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 by Bastien Chevreux
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

#ifndef _bas_skim_h_
#define _bas_skim_h_

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include <iostream>
#include <math.h>

#include "errorhandling/errorhandling.H"
#include "mira/readpool.H"


typedef unsigned long vhash_t;

struct vhrap_t {
  vhash_t vhash;       // vhash,
  uint32 readid;       // readid
  uint16 hashpos;      // (hash)position

  //bool operator<(const vhrap_t & other) const {return vhash < other.vhash;};
};

struct matchwithsorter_t{
  uint32 otherid;
  int32  eoffset;
  int32  percent_in_overlap;
  uint32 numhashes;
  bool   taken;
};

struct matchwith_t{
  uint32 otherid;
  int32  eoffset;
};


struct readhashmatch_t{
  uint32 rid2;
  int32  eoffset;
  uint16 hashpos1;
  //uint16 hashpos2;
};

typedef multimap< int32, matchwith_t> possible_overlaps_t;
typedef possible_overlaps_t::value_type posoverlap_pair_t;


typedef vector<set<uint32> > bannedoverlappairs_t;

typedef map< uint32, uint32> simplematches_t;
typedef simplematches_t::value_type simplematch_pair_t;

class Skim
{
private:
  ReadPool * SKIM3_readpool;

  vector<vector<vhrap_t>::const_iterator > SKIM3_vashortcuts_begin;
  vector<vector<vhrap_t>::const_iterator > SKIM3_vashortcuts_end;

  vector<uint8> SKIM3_megahubs;

  vector<readhashmatch_t> SKIM3_readhashmatches;

  // counts how many times a read has excellent hits with other reads
  vector<uint32> * SKIM3_overlapcounter;
  vector<uint32> * SKIM3_rawhashitcounter;
  bannedoverlappairs_t * SKIM3_bannedoverlaps;
  
  ofstream SKIM3_posfmatchfout;
  ofstream SKIM3_poscmatchfout;


  uint32 SKIM3_possiblehits;
  uint32 SKIM3_acceptedhits;;


  uint8  SKIM3_basesperhash;
  uint8  SKIM3_hashsavestepping;
  int32  SKIM3_percentrequired;

  uint32 SKIM3_maxhitsperread;


  uint32 SKIM3_totalhitschosen;
  uint32 SKIM_totalpermbans;


  uint32 SKIM_partfirstreadid;
  uint32 SKIM_partlastreadid;

  ProgressIndicator * SKIM_progressindicator;
  uint32 SKIM_progressend;


  bool SKIM3_onlyagainstrails;


public:


  //Functions
private:
  void foolCompiler();

  void prepareSkim();

  uint32 transformSeqToVariableHash(const uint32 readid, 
				    const char * seq, 
				    uint32 slen,
				    const string & readname,
				    const uint8 basesperhash,
				    vector<vhrap_t>::iterator & vhraparrayI,
				    const bool countonly,
				    const uint8 hashsavestepping);
  
  void checkForHashes(const int8 direction);

  void checkForPotentialHits(const uint32 actreadid, 
			     vector<matchwithsorter_t> & tmpmatchwith);

  void makeVHRAPArrayShortcuts(vector<vhrap_t> & vhraparray,
			       const uint8 basesperhash);

//////////////////////////////////

  void init();
  
  uint32 computePartition(uint32 maxmemusage, 
			  bool computenumpartitions);


public:
  Skim();
  Skim(Skim const &other);
  ~Skim();

  Skim const & operator=(Skim const & other);
  friend ostream & operator<<(ostream &ostr, Skim const &theskim);

  void discard();


  void skimGo (ReadPool & rp, 
	       string               & posfmatchname,
	       string               & poscmatchname,
	       bannedoverlappairs_t & bannedoverlaps,
	       vector<uint32>       & overlapcounter,
	       vector<uint32>       & rawhashhitcounter,

	       uint32 maxmemusage = 15000000,
	       
	       bool onlyagainstrails=false,
	       
	       uint8  basesperhash    = 16,
	       uint8  hashsavestepping= 4,
	       int32  percentrequired = 50,
	       uint32 maxhitsperread  = 200);
  

};


#endif

