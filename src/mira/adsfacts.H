/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2006 and later by Bastien Chevreux
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

#ifndef _adsfacts_h_
#define _adsfacts_h_

#include <fstream>
#include "stdinc/defines.H"



// Note that we need to save memory, therefore this class implicitly
//  limits the length of sequences that can be worked on to 2^15-1=32767
//  bases in length
//
// This is a crucial structure that will live a few million copies over
//  in memory for really large assemblies, so we need to keep it small
//
// This limitation is ok, an assembler should not try to align
//  sequences that come even close to this limit.

class AlignedDualSeqFacts
{
protected:
  readid_t ADSF_id1;
  readid_t ADSF_id2;

  uint16   ADSF_delta;            /* offset needed to get from *ADS_aligned_seq
				    to  *ADS_seq2;*/
  uint16   ADSF_id1_rightdelta;   // right offset for each id (total length
  uint16   ADSF_id2_rightdelta;   //  - endgaps)

  uint16   ADSF_total_len;            /* length of each aligned sequence and
					 consensus sequence (including
					 end gaps) */

  int8    ADSF_score_ratio;      // score to expected score in % from 0 to 100

  //uint16   ADSF_overlap_len;      // length of the overlap of both sequences
  //int8    ADSF_id1_direction;    // +1 forward || -1 complement
  //int8    ADSF_id2_direction;

  /* the direction of 1 & 2 compressed into one uint8 to save space
     0x0 = both reverse
     0x1 = id2 reverse, id1 forward
     0x2 = id2 forward, id1 reverse
     0x3 = both forward
   */
  uint8   ADSF_id1and2_directions:2;


   /*
      1 + acgtaaagggcccttt            1 + acgtaaagggcccttt         etc.pp
      2 +     aaaTggccctttacgt        2 -     aaaTggccctttacgt
      5pcm1 = 0;    3pcm1 = 8         5pcm1 = 0;    3pcm1 = 8
      5pcm2 = 3;    3pcm1 = 0         5pcm2 = 0;    3pcm1 = 3

      compressed version of ADS_... variables with 3 bits
      steps of 4: 0= 0-3, 1= 4-7, etc. ..., 7= >=28
      so, the above becomes:

      5pcm1 = 0;    3pcm1 = 2         5pcm1 = 0;    3pcm1 = 2
      5pcm2 = 0;    3pcm1 = 0         5pcm2 = 0;    3pcm1 = 0
  */
  uint16 ADSF_5pconmatch1:3;
  uint16 ADSF_3pconmatch1:3;
  uint16 ADSF_5pconmatch2:3;
  uint16 ADSF_3pconmatch2:3;

  uint16 ADSF_totalnonmatches:10;  // = ADS_nummismatches+ADS_numgaps; 1023 means >=1023

  //Functions
private:

public:
  friend std::ostream & operator<<(std::ostream &ostr, const AlignedDualSeqFacts &adsf);

  void serialiseOut(std::ostream & ostr);
  void serialiseIn(std::istream & ostr);

  inline readid_t getID1() const {return ADSF_id1;};
  inline readid_t getID2() const {return ADSF_id2;};
  readid_t getOtherID(readid_t id) const;
  int8 getSequenceDirection(readid_t id) const;

  uint32 getOffsetInAlignment(readid_t id) const;
  uint32 getRightOffsetInAlignment(readid_t id) const;

  uint32 get5pLenContiguousMatch(readid_t id) const;
  uint32 get3pLenContiguousMatch(readid_t id) const;
  uint32 getTotalNonMatches() const {return ADSF_totalnonmatches; };

  inline uint32 getTotalLen() const {return ADSF_total_len;};
  inline int8 getScoreRatio() const {return ADSF_score_ratio;};
  // added setScoreRatio() to enable fine tuning of
  //  troublemakers (get them out of the way for the pathfinder
  //  quick rules)
  inline void setScoreRatio(int8 sr) {ADSF_score_ratio=sr;};
  inline uint32 getOverlapLen() const {
    return ADSF_total_len-ADSF_delta-ADSF_id1_rightdelta-ADSF_id2_rightdelta;
  }

  void publicinit(readid_t id1,
		  readid_t id2,
		  uint16   delta,
		  uint16   id1_rightdelta,
		  uint16   id2_rightdelta,
		  uint16   totallen,
		  int8     id1_direction,
		  int8     id2_direction,
		  int8     score_ratio,
		  uint16 totalnonmatches, uint16 s5pcm1, uint16 s3pcm1, uint16 s5pcm2, uint16 s3pcm2);

  // these here are more for debugging: addRead() needs this for dumping help
  //  in the replay log
  inline uint16 getDelta() const { return ADSF_delta;}
  uint16 getRightDelta (readid_t) const;
};




#endif
