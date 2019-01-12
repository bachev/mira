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

  struct rahp_t {        // read and hashposition
    uint32 read;
    uint16 hashpos;
  };

  //Variables
private:
  uint32 SKIM_valid;

  uint16 * SKIM_forwardhashes;
  uint16 * SKIM_reversehashes;
  uint8 * SKIM_forwardhashok;
  uint8 * SKIM_reversehashok;

  uint32   SKIM_forwardhashcount[65536];
  uint32   SKIM_reversehashcount[65536];
  uint32   SKIM_forwardreadptrfield[65536];
  uint32   SKIM_reversereadptrfield[65536];
  rahp_t * SKIM_forwardrahpsph;       // read and hashposition sort per hash
  rahp_t * SKIM_reverserahpsph;       // read and hashposition sort per hash

  vector<int32> SKIM_hashperseq;

public:


  //Functions
private:
  void foolCompiler();

  void prepareSkim(ReadPool & rp);
  void transformSeqToHash (const char * seq, 
			   uint32 slen,
			   const char * readname,
			   uint16 * hashp,
			   uint8 * hashokp,
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

