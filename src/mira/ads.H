/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1998-2000 by the German Cancer Research Center (Deutsches
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


#ifndef _bas_score_h_
#define _bas_score_h_

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include "mira/parameters.H"
#include "mira/adsfacts.H"

#define ADS_MATSIZE 128

// we need to save memory, therefore length variables only 16 bit
class AlignedDualSeq : public AlignedDualSeqFacts
{
private:
  const static bool ADS_initialisedstatics;

  static std::vector<double> ADS_powofstarcounter;

  static uint32 ADS_c_matrix_valid;

  static char   ADS_consensus_matrix[ADS_MATSIZE][ADS_MATSIZE];
  static char   ADS_consensus_gap_matrix[ADS_MATSIZE][ADS_MATSIZE];
  static int8   ADS_realscore_matrix[ADS_MATSIZE][ADS_MATSIZE];
  static int8   ADS_expectedscore_matrix[ADS_MATSIZE][ADS_MATSIZE];
  static uint32 ADS_s_matrix_valid;

private:
  static bool staticInitialiser();

  MIRAParameters    * ADS_miraparams;  // parameters

  char   * ADS_aligned_seq1;         // the aligned sequences as provided
  char   * ADS_aligned_seq2;         //  by alignement routines
  char   * ADS_consensus_seq;        // consensus sequence
  char   * ADS_consensus_gap_seq;    /* consensus sequence, but here
					a star in aligned sequences is
					shown as a star in the consensus */
  int32    ADS_cursize;              // size in bytes of each of the 4 arrays above


                                 // (pointing into ADS_aligned_seqX)!
  char   * ADS_seq1;             // begin of the 1st seq. without leading gaps
  char   * ADS_seq2;             // begin of 2nd seq. without leading gaps

  int32    ADS_score;            // score of this alignment
  int32    ADS_expected_score;   // maximum score one could have expected of it
  int32    ADS_weight;           // weight

  uint16   ADS_nummismatches;    /* num of base/base mismatches */
  uint16   ADS_numgaps;          /* num of gap/base mismatches,
				    gap/gap (or oldgaps) is a match */
  uint16   ADS_maxcontiguousgaps; // max num of contiguous gaps

  uint16   ADS_len1;             // length of sequence (without end-gaps)
  uint16   ADS_len2;             // length of sequence (without end-gaps)

  uint8   ADS_valid;                  // 0: invalid, !=0 valid

  uint8   ADS_contained;        /* 0 no containment
				   1 seq1 contains seq2
				   2 seq1 and seq2 contain each other
				*/
  bool    ADS_affinegapscore; /* if true, gaps malus are counted with affine
				 gap score */

  int32   ADS_minbanddistance; // how far Align::rAlign() was from a band in this solution
  int32   ADS_bandwidthused;   // effective bandwidth used by Dynamic::computeBSimMatrix()

  //Functions
private:
  void foolCompiler();
  int32 score(bool enforce_clean_ends,
	      bool affinegapscore
    );
  void consensus();
  void zeroVars();
  void init();
  void matinit(char * arr, char a, char b, char value);
  void calcEndsLenContiguousMatch();

public:
  AlignedDualSeq(MIRAParameters * params);
  AlignedDualSeq(MIRAParameters * params,
		 const char * seq1, const char * seq2,
		 readid_t id1,  readid_t id2,
		 int8 id1dir, int8 id2dir,
		 bool enforce_clean_ends,
		 bool affinegapscore
    );
  AlignedDualSeq(AlignedDualSeq const &other);
  ~AlignedDualSeq();

  AlignedDualSeq const & operator=(AlignedDualSeq const &other);
  friend std::ostream & operator<<(std::ostream &ostr, const AlignedDualSeq &ads);

  void  discard();
  void  saveMem(bool delete_seq=true, bool delete_consseq=true);
  void  acquireSequences(const char * seq1,
			 const char * seq2,
			 readid_t id1,
			 readid_t id2,
			 int8 id1dir,
			 int8 id2dir,
			 bool enforce_clean_ends,
			 bool affinegapscore,
			 int32 bandwidthused,
			 int32 minbanddistance);

  inline bool  isValid() const {return ADS_valid!=0;}
  inline bool  isContained() const {return ADS_contained>0;};
         bool  isContained(int32 id) const;
  inline int32 getScore() const {return ADS_score;};
  inline int32 getExpectedScore() const {return ADS_expected_score;};
  inline int32 getWeight() const {return ADS_weight;};

  inline uint16 getNumMismatches() const {return ADS_nummismatches;};
  inline uint16 getNumGaps() const {return ADS_numgaps;};

  const char * getAlignedSequence(readid_t id) const;
  uint32 getLenOfAlignedSequence(readid_t id) const;

  const char * getSequenceAtOverlapStart(readid_t id) const;

  inline const char * getConsensusSequence() const {return ADS_consensus_seq;};
  inline const char * getGapedConsensusSequence() const {return ADS_consensus_gap_seq;};

  // information about how the alignment was computed
  inline int32 getMinBandDistance() const { return ADS_minbanddistance;}
  inline int32 getBandwidthUsed() const { return ADS_bandwidthused;}


  bool clipper(uint32 winlen,
	       int32 numerr,
	       int32 & retgoodlen1,
	       int32 & retgoodlen2) const;

};




class ADSEstimator
{
private:
  readid_t ADSE_id1;
  readid_t ADSE_id2;

  uint32 ADSE_lexpandof1;
  uint32 ADSE_rexpandof1;
  uint32 ADSE_lexpandof2;
  uint32 ADSE_rexpandof2;

  uint32 ADSE_totallen;
  uint32 ADSE_overlaplen;

  uint32 ADSE_len1;
  uint32 ADSE_len2;
  int8   ADSE_dir1;
  int8   ADSE_dir2;

public:
  ADSEstimator();
  //ADSEstimator(int32 offsets1tos2,
  //	       uint32 len1, uint32 len2,
  //	       int32 id1, int32 id2,
  //	       int8 id1dir, int8 id2dir);

  void calcNewEstimate(int32 offsets1tos2,
		       uint32 len1, uint32 len2,
		       readid_t id1, readid_t id2,
		       int8 id1dir, int8 id2dir);
  void calcNewEstimateFromSkim(int32 offsets1tos2,
			       uint32 len1, uint32 len2,
			       readid_t id1, readid_t id2,
			       int8 id1dir, int8 id2dir);

  inline uint32 getEstimatedOverlap() const {return ADSE_overlaplen;}
  inline uint32 getEstimatedTotalLen() const {return ADSE_totallen;}

  uint32 getEstimatedLeftExpand(readid_t id) const;
  uint32 getEstimatedRightExpand(readid_t id) const;

  uint32 getLen(readid_t id) const;
  inline int32 getID1() const {return ADSE_id1;}
  inline int32 getID2() const {return ADSE_id2;}

  uint32 getContainmentLevel() const;  /* 0 == not contained
					  1 == one sequence is contained
					       in the other
					  2 == both sequences contain
					       each other
				       */
  int32 getIDOfContained() const;

  void getPositionsInForwardSequenceofAlignmentStart(
    readid_t id1, readid_t id2,
    int32 & offset1, int32 & offset2) const;


  friend std::ostream & operator<<(std::ostream &ostr, const ADSEstimator &adse){
    ostr << "adse: id1: " << adse.ADSE_id1
	 << " id2: " << adse.ADSE_id2
	 << " le1: " << adse.ADSE_lexpandof1
	 << " re1: " << adse.ADSE_rexpandof1
	 << " le2: " << adse.ADSE_lexpandof2
	 << " re2: " << adse.ADSE_rexpandof2
	 << " ln1: " << adse.ADSE_len1
	 << " ln2: " << adse.ADSE_len2
	 << " tl: " << adse.ADSE_totallen
	 << " ol: " << adse.ADSE_overlaplen
	 << " dr1: " << static_cast<int16>(adse.ADSE_dir1)
	 << " dr2: " << static_cast<int16>(adse.ADSE_dir2)
	 << std::endl;
    return ostr;
  }

};

#endif
