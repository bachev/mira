/*
 * Written by Bastien Chevreux (BaCh)
 *
 * See down below in this file for MRC copyright of SCF header
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ) and Bastien Chevreux
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

/*
 * SCF routines
 *
 * Written by Bastien Chevreux (BaCh)
 *
 * - Can load SCF V2 and V3 files.
 * - Can load (gzip|compress|pack) compressed files.
 * - Handles the deltaed SCF V3 sample traces.
 * - Can dump the data as readable information to stdout.
 *
 * - Can save SCF file as version 3 (unpacked)
 * - Can cutout portions of SCF files (and trash the rest)
 *
 * Data structures taken from the Staden Package and adapted to my
 *  conventions. See below for the Staden copyright.
 */

#ifndef _scf_h_
#define _scf_h_


#ifndef SOLARIS_CLUDGE
#include <netinet/in.h>
#else
#include <sys/byteorder.h>
#endif

#include <string>
#include <fstream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "stdinc/stlincludes.H"


// ---------------- begin included part

/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */


#define SCF_MAGIC (static_cast<uint32>('.')<<24)+(static_cast<uint32>('s')<<16)+(static_cast<uint32>('c')<<8)+(static_cast<uint32>('f'))


// define uncertainty code sets
#define CSET_DEFAULT 0  /* {A,C,G,T,-} */
#define CSET_STADEN  1
#define CSET_NC_IUB  2  /* Pharmacia A.L.F. */
#define CSET_ALF     3  /* extended NC_IUB */
#define CSET_ABI     4  /* {A,C,G,T,N} */
#define CSET_IBI     5  /* IBI/Pustell */
#define CSET_DNASTAR 6  /* DNA* */
#define CSET_DNASIS  7
#define CSET_PCGENE  8  /* IG/PC-Gene */
#define CSET_GENIE   9  /* MicroGenie */

// define samples to delta_delta values
#define DELTA_IT 1



/*
 * Type definition for the Header structure
 */
typedef struct {
    uint32 magic_number;       /* SCF_MAGIC */
    uint32 samples;            /* Number of elements in Samples matrix */
    uint32 samples_offset;     /* Byte offset from start of file */
    uint32 bases;              /* Number of bases in Bases matrix */
    uint32 bases_left_clip;    /* OBSOLETE: No. bases in left clip (vector) */
    uint32 bases_right_clip;   /* OBSOLETE: No. bases in right clip (qual) */
    uint32 bases_offset;       /* Byte offset from start of file */
    uint32 comments_size;      /* Number of bytes in Comment section */
    uint32 comments_offset;    /* Byte offset from start of file */
    char   version[4];         /* "version.revision" */
    uint32 sample_size;        /* precision of samples (in bytes) */
    uint32 code_set;           /* uncertainty codes used */
    uint32 private_size;       /* size of private data, 0 if none */
    uint32 private_offset;     /* Byte offset from start of file */
    uint32 spare[18];          /* Unused */
} SCF_Header;





/*
 * Header.sample_size == 1.
 */
typedef struct {
    uint8 sample_A;                    /* Sample for A trace */
    uint8 sample_C;                    /* Sample for C trace */
    uint8 sample_G;                    /* Sample for G trace */
    uint8 sample_T;                    /* Sample for T trace */
} SCF_Samples1;

/*
 * Header.sample_size == 2.
 */
typedef struct {
    uint16 sample_A;                    /* Sample for A trace */
    uint16 sample_C;                    /* Sample for C trace */
    uint16 sample_G;                    /* Sample for G trace */
    uint16 sample_T;                    /* Sample for T trace */
} SCF_Samples2;



/*
 * Type definition for the sequence data
 */
typedef struct {
  uint32 peak_index;        /* Index into Samples matrix for base position */
  uint8 prob_A;            /* Probability of it being an A */
  uint8 prob_C;            /* Probability of it being an C */
  uint8 prob_G;            /* Probability of it being an G */
  uint8 prob_T;            /* Probability of it being an T */
  char base;                /* Base called */
  //    uint8 spare[3];          /* Spare */
  // New for SCF V3.10
  uint8 prob_sub;
  uint8 prob_ins;
  uint8 prob_del;
} SCF_Bases;


/*
 * Type definition for the comments
 */
typedef char SCF_Comments;      /* Zero terminated list of \n separated entries */

/*
 * All of the above structs in a single scf format.
 */

// BaCh 1997
/* Nope, I'll take a class, thanks
typedef struct {
    SCF_Header header;
    union SCF_Samples {
        SCF_Samples1 *samples1;
        SCF_Samples2 *samples2;
    } samples;
    SCF_Bases *bases;
    SCF_Comments *comments;
    char *private_data;
} Scf;
*/

// -------------- End included part ---------------

typedef struct{
  union SCF_Samples {
    SCF_Samples1 *samples1;
    SCF_Samples2 *samples2;
  } samples;
  SCF_Bases *bases;
} scf_type2;


class SCF
{
  // Variables
private:
  static std::vector<std::string> SCF_suffixalternatives;

  char * SCF_tmpfnameforload;
  std::string SCF_lastloadedscfname;

protected:
  // Header, comments and private data are the same in all SCF versions
  SCF_Header     SCF_header;
  SCF_Comments * SCF_comments;
  char         * SCF_private_data;

  uint16 * SCF_samples_A;
  uint16 * SCF_samples_C;
  uint16 * SCF_samples_G;
  uint16 * SCF_samples_T;

  uint32 * SCF_peak_index;
  uint8  * SCF_prob_A;
  uint8  * SCF_prob_C;
  uint8  * SCF_prob_G;
  uint8  * SCF_prob_T;
  uint8  * SCF_prob_sub;
  uint8  * SCF_prob_ins;
  uint8  * SCF_prob_del;
  char   * SCF_bases;

  //if SCF version <3, then read or write the file with this temporary
  // structure. The struct itself is already here, the arrays must be
  // allocated and deallocated
  scf_type2 SCF_tmp_scf2;


  //class data
protected:
  uint32 SCF_valid_data;                // is the object initialised?


  //functions
private:
  // functions to convert byte order in architectures not
  //  working in the network byte order (e.g. Intel x86)
  //  on other machines, these functions should automatically
  //  be optimized away by the compiler
  inline void convert4ByteOrderToHost(uint32 *ptr, uint32 num);
  inline void convert2ByteOrderToHost(uint16 *ptr, uint32 num);
  inline void convertHeaderByteOrderToHost();

  inline void convert4ByteOrderToNetwork(uint32 *ptr, uint32 num);
  inline void convert2ByteOrderToNetwork(uint16 *ptr, uint32 num);
  inline void convertHeaderByteOrderToNetwork();

  void dumpSample(const uint16 * ptr);   //dump a Sample to stdout
  void undelta16(uint16 * ptr);            //undelta SCF-V3 Samples
  void delta16(uint16 * ptr);              //undelta SCF-V3 Samples
  void undelta8(uint8 * ptr);            //undelta SCF-V3 Samples
  void delta8(uint8 * ptr);              //undelta SCF-V3 Samples

  void recalcHeader();

  void zeroVars();

public:
  SCF();
  SCF(SCF const &other);
  virtual ~SCF();

  SCF const & operator=(SCF const & other);
//  friend ostream & operator<<(ostream &ostr, SCF const &scf);

  int8 load(const char * filename, std::string & emsg);
  const std::string & getNameOfLastLoadedFile() const { return SCF_lastloadedscfname;}
  int8 dump();                           //dump an SCF to stdout
  void save(const char * filename);
  void discard();                        // discard the SCF

  bool checkSCFDataOK();
  void correctPeakIndices();
  void transposeAmbiguityCodes();

  void cutBases(uint32 frombase, uint32 tobase);  // cut bases out of file

  // provide new bases to a file, replace old ones
  //  each base has only one probability, gotta write other functions for multiple
  //  probabilities.
//  void setBases(      uint32   numbases,
//		const char   * newbases,
//		const uint8  * newprobabilities,
//		const uint32 * newpeakpositions);
  void setBases(const std::string & newbases,
		const std::vector<uint8> & newprobabilities,
		const std::vector<uint32> & newpeakpositions);
  void addComment(const char * comment);

  uint32 getNumBases() const {return SCF_header.bases;};
  char getBase(uint32 index) const;
  uint8 getCalledBaseProb(uint32 index) const;
  uint8 getAProb(uint32 index) const;
  uint8 getCProb(uint32 index) const;
  uint8 getGProb(uint32 index) const;
  uint8 getTProb(uint32 index) const;

  uint32 getNumSamples() const { return SCF_header.samples; };
  uint16 getASample(uint32 samplepos) const;
  uint16 getCSample(uint32 samplepos) const;
  uint16 getGSample(uint32 samplepos) const;
  uint16 getTSample(uint32 samplepos) const;

  uint32 getPeakIndex(uint32 basepos) const;

  const SCF_Comments * getComments() const {return SCF_comments;};

};


#endif
