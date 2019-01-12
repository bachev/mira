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

#ifndef _bas_align_h_
#define _bas_align_h_

#include <iostream>

#include <list>

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include "mira/parameters.H"
#include "mira/dynamic.H"

class AlignedDualSeq;


class Align : public Dynamic
{
public:
  static uint64 AL_alloccount;

private:
  MIRAParameters    * AL_miraparams;

  uint8  AL_valid;

  char * AL_alseq1;		/* ptr to the first aligned sequence
				   being build or having been build.
				   This is an array, mem allocated by
				   the instance  */
  char * AL_alseq2;		// dito for the second sequence

  uint32 AL_as12size;           // size of alseq1 / alseq2 array

  char * AL_seq1ptr;		// ptr for the actual position in seq 1
  char * AL_seq2ptr;		// dito for sequence 2

  int32  AL_id1;                // optional userdefined id for the seq.
  int32  AL_id2;                //  will be transfered to the ADS
  int8   AL_id1dir;             // optional direction for the seq (1/-1)
  int8   AL_id2dir;             //  will be transfered to the ADS

  uint32 AL_allen;		// length of the actual alignment consensus
  uint32 AL_align_maxlen;	/* maximum alowed length of it
				   =len(seq1)+len(seq2)  */

  uint8  AL_new_solution;

  bool   AL_error_hit_band;     /* true if traceback hit band in matrix,
				   bandwidth should be increased */

  uint32 AL_cutoff_counter;	// failsafe counter for stopping calc

  AlignedDualSeq       * AL_tmpads;
  std::list<AlignedDualSeq> * AL_adslist;

  // cached parameters from MIRAPARAMS so that rAlign doesn't need to
  //  get pointers to MIRAPARAMS every time
  // furthermore, we sometimes want to tweak these values just a bit
  //  from externally, but not directly in MIRAPARAMS
  // same type as originals in align_parameters structure of structs.H

  int32  AL_mpcache_dyn_score_multiplier;
  int32  AL_mpcache_dyn_score_gap;

  uint32 AL_mpcache_al_max_cutoff;
  uint32 AL_mpcache_al_min_score;
  uint32 AL_mpcache_al_min_overlap;
  uint32 AL_mpcache_al_min_relscore;

  // if these !=0, use these in above, else default from MIRAPARAMS
  uint32 AL_mpset_al_min_overlap;
  uint32 AL_mpset_al_min_relscore;

  // rAlign stores here the minimum distance it was away from left/right band
  // can be used to optimise alignment of very long sequences (PacBio)
  //  in later recalculations if this info is somehow stored
  // It's passed on to the AlignedDualSeq object for every alignment found
  int32 AL_minbanddistance;


  // RLE
  bool AL_userle;
  bool AL_rle_create_packed_align;
  bool AL_rle_create_stdleft_align;
  bool AL_rle_create_nonstdright_align;

  std::string         AL_rles1;    // space for packed seq1
  std::string         AL_rles2;    // space for packed seq2
  std::vector<uint32> AL_rlev1;    // rle values s1
  std::vector<uint32> AL_rlev2;    // rle values s2
  std::string         AL_unrles1;    // space for unpacked alignment of seq1
  std::string         AL_unrles2;    // space for unpacked alignment of seq2

  // Timing checks
  suseconds_t AL_timing_acquires;
  suseconds_t AL_timing_fullalign;
  suseconds_t AL_timing_prepalign;
  suseconds_t AL_timing_raligntot;
  suseconds_t AL_timing_ra_adsacquire;
  suseconds_t AL_timing_ra_adslist;

public:
  uint32 AL_no_solutions;       // number of alignment solutions found
  uint32 AL_no_diff_solutions;	// where no # these have different starting
  uint32 AL_max_relscore;

  bool AL_enforce_clean_ends;
  bool AL_affine_gap_scorees;

  //Functions
private:
  void foolCompiler();

  void init();
  void termAlign();
  void termAlignFancy();
  void rAlign(uint32 i, uint32 j, char lastdir, bool hadn);
  void prepareAlign(std::list<AlignedDualSeq> * adslist);
  void setRAlignParams();

  void pa_packSeqToRLE(const char * seq,
		       uint32 len,
		       std::string & rles,
		       std::vector<uint32> & rlev,
		       int32 & eoffset);
  void ra_expandRLEAlignments(const char * s1, const char * s2, bool shakedirstdleft);


public:
  Align(MIRAParameters * params);
  ~Align();

  void discard();
  inline void setUseRLE(bool b) {AL_userle=b;}
  inline void setAffineGapScore(bool b) {AL_affine_gap_scorees=b;}
  inline void setEnforceCleanEnds(bool b) {AL_enforce_clean_ends=b;}
  void acquireSequences(const char * seq1,
			uint32 len1,
			const char * seq2,
			uint32 len2,
			int32 id1,
			int32 id2,
			int8 id1dir,
			int8 id2dir,
			bool calcwithoffset=true,
			int32 expectedoffset=0);

  void useSpecialMinOverlap(uint32 a) {AL_mpset_al_min_overlap=a;}
  void useSpecialMinRelScore(uint32 a) {AL_mpset_al_min_relscore=a;}

  void simpleAlign(std::list<AlignedDualSeq> * adslist);
  void fullAlign(std::list<AlignedDualSeq> * adslist);

  bool wasBandHit() const {return AL_error_hit_band;}

  void MONITOR() {std::cout << "MONITOR: " << DYN_simmatrix << std::endl;}

  void coutWhatWasGiven();

  void dumpTimings();
  void resetTimings();
};


#endif
