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


/*
 * Dynamic programming routines
 *
 * Given two sequences, this object computes the dynamic programming
 *  matrix for them. The default values make the score to an edit
 *  distance without terminal gap penalties.
 *
 * One will probably not want to use it as standalone! For the time
 *  being, it is inherited by the class Align.
 *
 * Can do:
 *   - compute dyn. matrix with and without terminal gap penalties
 *   - does not care about upper-lower case sequences; - aliases N
 *   - once acquired, the sequences are kept locally. So the source
 *      sequences do not need to be kept.
 *   - instance reusable, giving it new sequences will trash the
 *      old values.
 * TODO:
 *   - see .C
 */

#ifndef _dynamic_h_
#define _dynamic_h_

#include <stdinc/defines.H>

#include <errorhandling/errorhandling.H>
#include <util/dptools.H>

class MIRAParameters;

#define DYN_MATSIZE 128
#define DYN_BANDLIMIT 1<<30

class Dynamic
{
public:

  static uint64 DYN_alloccounts;
  static uint64 DYN_alloccountm;

  int32 DYN_maxscore;          // Max score within simmatrix
  int32 DYN_lastrc_maxscore;   /* Max score within last row and column
				 of simmatrix */


protected:
  MIRAParameters * DYN_miraparams;

  static int32 DYN_match_matrix[DYN_MATSIZE][DYN_MATSIZE];
  static int16 DYN_matvalid;

  char   * DYN_sequence1;         // sequence 1
  char   * DYN_sequence2;         // sequence 2

  uint32   DYN_s1size;            // size of seq1 array
  uint32   DYN_s2size;            // size of seq2 array

  uint32   DYN_len_seq1;          //  and its length are the ROWS
  uint32   DYN_len_seq2;          //  and its length are the COLUMNS

  int32    DYN_knowngaps;         /* number of gaps that need to be additionally
				     accounted for as these were seen in the
				     sequences */

  int32  * DYN_simmatrix;        /* Similarity matrix, always big enough to contain at least
				    (ROWS+1)*(COLUMNS+1) elements */
  uint32   DYN_smsize;           // size of the matrix;

  int32    DYN_leftbandx;        // x coordinate where left band diagonal cuts x axis (or DYN_BANDLIMIT if unused)
  int32    DYN_rightbandx;       // x coordinate where right band diagonal cuts x axis (or DYN_BANDLIMIT if unused)
  int32    DYN_bandsafety;       // how much safety buffer to have when looking for band hit
                                 //  atm: 1% of expected overlap, min 1, max 15

  int32    DYN_bandwidth;        // stored by computeBSimMatrix(): effective bandwidth used in last SW calc

  // deferred calculation of matrix (handled by Align class)
  bool   DYN_validseq;            // do we have valid sequences?
  bool   DYN_matrixcalculated;    // has the matrix for the sequences been calulated
  bool   DYN_calcwithoffset;      // calculate with offset?
  int32  DYN_eoffset;             // what offset

  suseconds_t DYN_timing_bswmatrix;  // time last bsw took
  suseconds_t DYN_timing_bswm_setup;
  suseconds_t DYN_timing_bswm_p1;
  suseconds_t DYN_timing_bswm_p2a;
  suseconds_t DYN_timing_bswm_p2b;
  suseconds_t DYN_timing_bswm_p3;
  suseconds_t DYN_timing_bswm_cleanband;
  suseconds_t DYN_timing_seqcopy;

private:
  void foolCompiler();

  int32 sequenceCopy(char * to, const char * from, uint32 len);
  void zeroVars();
  void matinit(char a, char b, int32 score);
  void init();

  void computeSimMatrix();
  void computeBSimMatrix();

public:
  Dynamic(MIRAParameters * params);
  Dynamic(MIRAParameters * params,
	  const char * seq1,
	  uint32 len1,
	  const char * seq2,
	  uint32 len2,
	  bool calcwithoffset=true,
	  int32 expectedoffset=0);
  virtual ~Dynamic();

  void discard();
  void setSequences(const char * seq1,
		    uint32 len1,
		    const char * seq2,
		    uint32 len2,
		    bool calcwithoffset=true,
		    int32 expectedoffset=0);
  void dump();
  void computeMatrix();

  void coutWhatWasGiven();
};


#endif
