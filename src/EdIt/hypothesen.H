/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
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
 *
 * hypotheses generation - header file
 *
 * Written by Thomas Pfiserer
 * Version 0.2.0       03.09.1998
 */

#ifndef _edit_hypothesen_H
#define _edit_hypothesen_H


#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

//#include <vector>
//#include <map>
//#include <list>
//#include "algo.h"

#include "examine/scf_look.H"
#include "examine/buffer.H"
#include "examine/globals.H"

#include "io/scf.H"
#include "mira/read.H"
#include "mira/contig.H"
#include "util/progressindic.H"
#include "EdIt/parameters.H"


// 	$Id$

// Minimum length for a Contig to be edited
#define MIN_CONTIG_EDIT_LENGTH 200


struct editoperation {
  int16 sequenceId;
  int16 position;
  int16 operation;
  char  base_old;
  char  base_new;

  vector<afh_Info*> readEdits;      // FIXME V1
};



struct scfInfo {
  vector <int> scfInfo_position;
  const  Contig::contigread_t  *scfInfo_read;
  int32  scfInfo_readId;
  Contig *scfInfo_contig;
  char   *scfInfo_name;
  int32  scfInfo_lcutoff;
  int32  scfInfo_rcutoff;
  int32  scfInfo_length;
  int32  scfInfo_dbpos;
  int32  scfInfo_sequenceId;
  bool   scfInfo_confirmed;
  bool   scfInfo_hidden;
};


class sequenceInfo {
public:
  vector<char> si_sequence;
  int16        si_occurrences;    // how often is this sequence in the problem
  int16        si_lastEditPos;
  int16        si_insertCount;
};



class fault_region
{
private:
  static multitag_t fr_tagED_C;
  static multitag_t fr_tagED_D;
  static multitag_t fr_tagED_I;

  enum { fr_fillchar = '.', fr_gapchar = '*',
	 fr_eorchar  = '_', fr_maxsteps = 15 };
  enum { fr_insert = 1, fr_delete = 2, fr_alter = 3, fr_syntactic = 4  };

  cfh_Info scfEditHypotheses;
  int32  fr_longest;
  int32  fr_score;           // # pairwise mismatches
  int32  fr_steps;           // # edit operations
  int32  fr_solution_score;  // score of a solution
  bool   fr_has_eor;
  vector<sequenceInfo>   fr_sequenceList;
  vector<editoperation>  fr_edits;


  int32 resizeRegion(int32 size);
  int32 pairwiseMismatch();
  void  performEditOperationOnRegion(editoperation &edOp);
  int32 shufflePadInSequence(int32 read_id, int32 pos);
  int32 checkParameters(const int32 readId, const int32 readPos);
public:
  fault_region();
  ~fault_region();
  fault_region(fault_region const &other);
  const fault_region& operator= (const fault_region &other);
  void clear();
  vector<sequenceInfo>::iterator findSequence(const char *target);
  int32 addSequence(const char *s, scfInfo &aScfInfo);
  bool  possibleSolution() const;
  bool  validSolution() const;
  void  setColumnToBase(int32 col, char base);

  int32 pairwiseMismatchColumn(int32 col);

  void  performOperationsOnReads(vector<scfInfo> hp_scfInfo,
				const EDITParameters &p, int32 locomotion);
  int32 insertBase(const int32 readId, const int32 readPos, const char base);
  int32 deleteBase(const int32 readId, const int32 readPos);
  int32 alterBase(const int32 readId, const int32 readPos, const char base);
  char readBase(const int32 readId, const int32 readPos) const;
  int32 compare(fault_region &fr);
  int32 readSequenceCount() const;
  int32 readScore() const;
  int32 readSolutionScore() const;
  int32 calcSolutionScore();
  int32 editCount();
  bool  isBaseInColumn(int32 col, char aBase);
  const sequenceInfo& readSequence(int32 sequenceId);
  vector<editoperation>& readEditOperations();
  int32 size() const;
  int32 realSize() const;

  void  appendAfh(afh_Info *aAfh)
    { scfEditHypotheses.appendAfh(aAfh); }

  cfh_Info& getCfh()
    { return scfEditHypotheses; }

  bool  generateAllHypotheses(fault_region &target,
			      vector<scfInfo> &scf,
			      // cfh_Info &aHypotheses,
			      EDITParameters *p);

  afh_Info* createSingleHypotheses(scfInfo &theScfRead,
				   const sequenceInfo &s, int32 question,
				   char new_base, char old_base, int32 pos);

  void  shufflePads();

  friend ostream &operator<<(ostream &ostr, fault_region const &i);
};



class hypotheses_pool
{
private:
  enum { hp_usesize = 20, hp_poolsize = 30 };
  enum { fr_fillchar = '.', fr_gapchar = '*',
	 fr_eorchar  = '_' };

  vector<fault_region>  hp_ready;
  fault_region          hp_start;
  vector<fault_region>  hp_pool;
  vector<scfInfo>       hp_scfInfo;
  ofstream              logfile;
  EDITParameters        *hp_parameter;

public:
  hypotheses_pool(const fault_region startzustand, vector<scfInfo> &si);
  hypotheses_pool(const fault_region startzustand, vector<scfInfo> &si,
		  EDITParameters *p);
  ~hypotheses_pool();

  bool isDoubleStranded(bool& forward, bool &reverse);
  bool isDoublestrandedOrCoverage(bool& forward, bool &reverse, int32 c);
  bool usesDifferentChemistry();
  bool hasSufficientInformation(bool &forward, bool &reverse,
				bool useChemistry, int32 minCoverage);

  bool  makeDoubleStranded(Contig &aContig, int32 contigPosVon,
			  int32 contigPosBis,
			  bool forward, bool reverse,
			  const EDITParameters &p);
  void   pushStep(fault_region &fr);
  void   iterate(fault_region &fr);
  void   iteratePool();
  void   iteratePool(int32 anz);
  void   iterateFirstSolutionPlus(int32 anz);
  int32  minimalEdits();
  void   showResult();
  int32  locomotion(const fault_region &r) const;
  bool                  createCfhInfoList(EDITParameters &p);
  vector<fault_region>& getHypothesesList() { return hp_ready; }
  vector<scfInfo>&      getScfInfo() { return hp_scfInfo; }

  void   setFirstAndLastColumn(char first, char last);
  int32  checkPartialSolution(fault_region &fr);

  friend ostream &operator<<(ostream &ostr, hypotheses_pool const &i);

};


// find the consensus base
//char calledBase(const Contig::consensus_counts_t &cc);

// find next mismatch in a contig starting from the pos-th base
deque<Contig::consensus_counts_t>::iterator
findNextMismatch(const Contig &aContig, int32 &pos,
		  deque<Contig::consensus_counts_t> &theConsensus);
deque<Contig::consensus_counts_t>::iterator
findPrevMismatch(const Contig &aContig, int32 &pos,
		 deque<Contig::consensus_counts_t> &theConsensus);


// create a fault region for the given contig ranging from the base
// 'von' to the base 'bis'. returns an initialised hypotheses-pool object.
hypotheses_pool* createRegion(Contig &aContig, int32 von, int32 bis,
			      EDITParameters &p);

scfInfo*
faultRegionFromComplementRead(const Contig::contigread_t &R, int32 von, int32 bis);

scfInfo*
faultRegionFromRead(Contig &aContig, const Contig::contigread_t &R,
		    int32 readIndex, int32 von, int32 bis, char *sequence);


bool extendRead(const Contig &aContig, const Contig::contigread_t &aRead,
		int32 dbPos, int32 anzBases,
		int32 umgebung, int32 &offset_to_readpos);




int32 extendFaultRegionRight(const deque<Contig::consensus_counts_t>::iterator POS,
			   const deque<Contig::consensus_counts_t>::iterator END);

int32 extendFaultRegionLeft(const deque<Contig::consensus_counts_t>::iterator POS,
			  const deque<Contig::consensus_counts_t>::iterator START);

hypotheses_pool * createFaultRegion(Contig &aContig, int32 von, int32 bis,
				     EDITParameters &p);

bool examineFaultRegion(Contig &aContig, int32 von, int32 bis,
			hypotheses_pool *aPool, EDITParameters &p);

void findFaultRegionBack(Contig &aContig, EDITParameters &p);
void editContigBack(Contig &aContig, EDITParameters &p);
void editTowardsConsensusBack(Contig &aContig, EDITParameters &p);
void editSmallRegionsBack(Contig &aContig, EDITParameters &p,
			  int32 von, int32 bis);

int32 scfPosFromDBPos(const int32 pos, const Contig::contigread_t &aRead);
int32 scfInfo_contigPosFirstBase(scfInfo &si);


int32 dbPosFromReadPos(const int32 dbPos, const Contig::contigread_t &cRead);
int32 readPosFromContigPos(const int32 contigPos, const Contig::contigread_t &cRead);

#endif
