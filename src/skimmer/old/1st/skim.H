/*
 * Copyright (c) Deutsches Krebsforschungszentrum Heidelberg 1998.
 * All rights reserved.
 *
 *  routines
 *
 * Written by Bastien Chevreux (BaCh)
 *
 */

#ifndef _bas_skim_h_
#define _bas_skim_h_

#include <iostream.h>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"
#include "mira/readpool.H"


class Skim
{
  //Variables
private:
  uint32 SKIM_valid;

  uint16 * SKIM_forwardhashes;
  uint16 * SKIM_reversehashes;
  uint32 SKIM_forwardhashcount[65536];
  uint32 SKIM_reversehashcount[65536];
  uint32 SKIM_forwardreadptrfield[65536];
  uint32 SKIM_reversereadptrfield[65536];
  uint32 * SKIM_forwardrsph;             // read sort per hash
  uint32 * SKIM_reversersph;             // read sort per hash

  vector<uint32> SKIM_hashperseq;

public:


  //Functions
private:
  void foolCompiler();

  void prepareSkim(ReadPool & rp);
  uint16 * transformSeqToHash (const char * seq, 
			       uint32 slen,
			       const char * readname,
			       uint16 * hashp,
			       uint32 * hashcount);



public:
  Skim(ReadPool & rp);
  Skim(Skim const &other);
  ~Skim();

  Skim const & operator=(Skim const & other);
  friend ostream & operator<<(ostream &ostr, Skim const &theskim);

  void discard();


  void go ();

};


#endif

