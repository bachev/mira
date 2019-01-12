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
 * faultRegion
 * There is no class for handling fault regions. This file is a collection
 * of the functions used to create fault regions from a Contig
 *
 * Written by Thomas Pfisterer
 * Version 0.80       16.03.1999
 *
 * May 2003 Small bugfixes by BaCh
 */


#include <stdlib.h>
#include <iostream>

#include "errorhandling/errorhandling.H"
#include "mira/readpool.H"
#include "mira/contig.H"
#include "mira/align.H"

#include "EdIt/hypothesen.H"
#include "caf/caf.H"




// *************************************************************************
// *
// *  Extend a fault region right
// *  returns the number of identical bases right of POS
// *
// *************************************************************************

int32 extendFaultRegionRight(const Contig::cccontainer_t::iterator POS,
	 		     const Contig::cccontainer_t::iterator END)
{
  Contig::cccontainer_t::iterator I = POS;

  if (I == END) return 0;

  try {
    char lastBase = calledBase(*(I+1));

    do {
      I++;
    } while (I != END && lastBase == calledBase(*I));
  }
  catch (Notify n) {
    n.handleError("EdIt Error: Can not extend read to the right (1)");
  }

  return I-POS-1;
}


// *************************************************************************
// *
// *  Extend a fault region left
// *  returns the number of identical bases left of POS
// *
// *************************************************************************

int32 extendFaultRegionLeft(const Contig::cccontainer_t::iterator POS,
			    const Contig::cccontainer_t::iterator START)
{
  Contig::cccontainer_t::iterator I = POS;

  if (I == START) return 0;

  try {
    char firstBase = calledBase(*(I-1));
    do {
      I--;
    } while (I != START && firstBase == calledBase(*I));
  }

  catch (Notify n) {
    n.handleError("EdIt Error: Can not extend read to the right (2)");
  }

  return I-POS+1;
}


// *******************************************************
//    checkColumnForMismatch:
//    compares the consensus bases with all reads
//    returns true if there is a mismatch, false otherwise
// *******************************************************


bool checkColumnForMismatch(const Contig &aContig, const int32 pos,
			    Contig::cccontainer_t &theConsensus)
{
  vector<char>::const_iterator seqIt;
  vector<int32> theActiveReads;
  const vector<Contig::contigread_t> &allReads = aContig.getContigReads();

  aContig.getCRIndicesAtContigPosition(theActiveReads, pos-1);
  vector<int32>::iterator I = theActiveReads.begin();

  Contig::cccontainer_t::iterator I_CONS;
  I_CONS = theConsensus.begin() + pos - 1;
  char consensusBase = calledBase(*I_CONS);

  // step through all active reads at this position
  while (I != theActiveReads.end()) {
    const Contig::contigread_t &r = allReads[*I];

    if (r.direction > 0) {
      seqIt = r.read.getClippedSeqIterator();
    } else {
      seqIt = r.read.getClippedComplementSeqIterator();
    }
    seqIt += I_CONS - theConsensus.begin() - r.offset;

    // search for a mismatch of the read with the consensus sequence

    if (isUndefinedBase(consensusBase)) {
      return true;
    }

    if (toupper(*seqIt) != consensusBase) {
      return true;
    };

    I++;
  }

  return false;
}


// *******************************************************
//    findPrevMismatch:
//    step backwards through a Contig base by base searching
//    for mismatches in the high quality parts starting with
//    position "pos". Returns an iterator to the fault
//     position and updates the "pos" variable
// *******************************************************


Contig::cccontainer_t::iterator
findPrevMismatch(const Contig &aContig, int32 &pos,
	         Contig::cccontainer_t &theConsensus)
{
  Contig::cccontainer_t::iterator I_CONS;

  try {
    I_CONS = theConsensus.begin() + pos - 1;
    // step through the contig base by base
    while (I_CONS != theConsensus.begin()) {
      if (checkColumnForMismatch(aContig, pos, theConsensus))
	 return I_CONS;

      I_CONS--;
      pos--;
    }
  }
  catch (Notify n) {
    throw Notify(Notify::WARNING, "EdIt Error: findPrevMismatch failed");
  }

  // no mismatch found: returns theConsensus.begin()
  return I_CONS;
}



// **********************************************************
//    createFaultRegion
//    creates and examines a region in the contig given
//    tries to double-strand single stranded fault regions
//    before hypotheses generation
// **********************************************************


hypotheses_pool* createFaultRegion(Contig &aContig, int32 von, int32 bis,
				   EDITParameters &p)
{
  hypotheses_pool *aPool;
  static int32 region_count = 0;

  if (p.evalAt(von) == false && p.evalAt(bis) == false) {
#ifndef RUN_ONLY
    if (p.isVerbose(7)) {
      cout << "FaultRegion out of edit range.... skipping." << endl;
    }
#endif
    return nullptr;
  }


#ifndef RUN_ONLY
  if (p.isVerbose()) {
    cout << "\nExamine Region " <<  von << " - " << bis
	 << "\t(" << ++region_count << ")" << endl;
  }
#endif

  aPool = createRegion(aContig, von, bis, p);

  if (aPool == nullptr) {
    // no region created - e.g. the region is to long

#ifndef RUN_ONLY
    if (p.isVerbose()) {
      cout << "Region not created" << endl;
    }
#endif
    return nullptr;
  }

#ifndef RUN_ONLY
  if (p.isVerbose(8)) {
    cout << "Region created." << endl;
    cout << *aPool << endl;
  }
#endif

  return aPool;
}



bool examineFaultRegion(Contig &aContig, int32 von, int32 bis,
			hypotheses_pool *aPool, EDITParameters &p)
{
  bool forward, reverse;

  if (aPool == nullptr) {
    return false;
  }

  // if not double-stranded and we try to make
  // it double stranded then try


  if (!aPool->isDoubleStranded(forward, reverse) &&
      p.makeDoubleStranded())
    {
      aPool->makeDoubleStranded(aContig, von, bis, forward, reverse, p);
    }

#ifndef RUN_ONLY
  if (p.isVerbose()) {
    cout << *aPool;
  }
#endif

  if (aPool->hasSufficientInformation(forward, reverse, true, p.getMinCoverage())) {

    //  if (aPool->isDoubleStrandedOrCoverage(forward, reverse,
    //					p.getMinCoverage())) {

    aPool->iterateFirstSolutionPlus(2);

    vector<fault_region>::iterator I;
    bool conf = false;
    int32 count = 0;

    if (aPool->createCfhInfoList(p) == false) {
      return false;
    }

    I = aPool->getHypothesesList().begin();

    if (p.isVerbose(4)) {
      aPool->showResult();
    }

    if (p.doEval() == true) {

      while (I != aPool->getHypothesesList().end()) {

	if (true == p.getEditLowQualityBases()) {
	  if (p.isVerbose()) {
	    cout << "Confirm Quality" << endl;
	  }
	  I->getCfh().confirmQuality(p.getLowQualityThreshold());
	}

#ifndef RUN_ONLY
	if (p.isVerbose(5)) {
	  cout << I->getCfh();
	  cout << *I;
	}
#endif

	conf = I->getCfh().eval(true, // confirm until rejected
				p.isVerbose(),
				p.getConfirmationThreshold(),
				true  // strict evalution
				);
#ifndef RUN_ONLY
	if (p.isVerbose(5)) {
          cout << "Evaluation result: " << conf << endl;
	}
#endif

	if (conf) break;

	I++;
      }
    }

    if (conf == true ) {
#ifndef RUN_ONLY
      if (p.isVerbose(3)) {
	cout << "\nConfirmed: " << endl;
	cout << I->getCfh() << endl;
      }
#endif

      if (true == p.doEdit()) {
	I->performOperationsOnReads(aPool->getScfInfo(), p,
				    aPool->locomotion(*I));
      }
    }
#ifndef RUN_ONLY
    else {
      if (p.isVerbose(3)) {
	cout << "Not strictly confirmed" << endl;
      }
    }
#endif



    // suche die beste unvollständige bestätigte hypothese....

    if (p.doEval() == true && conf == false ) {
      float    best_score = 0;
      fault_region *best_region = nullptr;

      I = aPool->getHypothesesList().begin();

#ifndef RUN_ONLY
      if (p.isVerbose(2)) {
	cout << "Searching best partial solution" << endl;
      }
#endif

      while (I != aPool->getHypothesesList().end()) {
	I->getCfh().eval(true,
			 p.isVerbose(8),
			 p.getConfirmationThreshold(),
			 false);

	if (p.isVerbose(3)) {
	  cout << "Score: " << I->getCfh().getRelScore() << endl;
	  cout << I->getCfh();
	}

	if (I->getCfh().getRelScore() > best_score &&
	    aPool->checkPartialSolution(*I) == 1) {

	  best_score = I->getCfh().getRelScore();
	  // BaCh: changed that for gcc3
	  //best_region = *(&I);
	  best_region = &(*I);
	}

	I++;
      }

      // den besten unvollständigen kandidaten editieren...
      if (conf == false && p.getStrictEvaluation() == false &&
	  best_region != nullptr && best_score > 0.8) {

#ifndef RUN_ONLY
	if (p.isVerbose(5)) {
	  cout << "Best partial solution:" << endl;
	  cout << best_region->getCfh();
	}
#endif

	(best_region->getCfh()).eval(false,
				     p.isVerbose(),
				     p.getConfirmationThreshold(),
				     false);

	(best_region->getCfh()).evalExtendedAlterOperations
	                             (p.getConfirmationThreshold());

	if (true == p.doEdit() && aPool->checkPartialSolution(*best_region)) {
	  best_region->performOperationsOnReads(aPool->getScfInfo(), p,
				      aPool->locomotion(*best_region));
	}
      }
    }

  } else {

#ifndef RUN_ONLY
  if (p.isVerbose()) {
    cout << "Hypothesis is single-stranded!" << endl << endl;
  }
#endif

  }

  delete aPool;
  return true;
}



// **********************************************************
//    findFaultRegionBack:
//    steps through a Contig via findPrevMismatch and creates
//    fault regions around them. Two mismatches are treated
//    as dependent (i.e. in different faultRegions) if they
//    are only seperated by a sequence of identical bases
// **********************************************************


void findFaultRegionBack(Contig &aContig, EDITParameters &p)
{
  findFaultRegionBack(aContig, p);
}



void editContigBack(Contig &aContig, EDITParameters &p)
{
  static int32 regionnr = 0;
  static int32 contignr = 0;
  Contig::cccontainer_t &theConsensus =
    const_cast<Contig::cccontainer_t&>(aContig.getConsensusCounts());
  Contig::cccontainer_t::iterator I_CONS = theConsensus.end();

  int32 pos = 0;
  int32 beginRegion = 0;
  int32 endRegion = 0;

  int32 oldContigSize;

  // set pos to the last base of the contig
  oldContigSize = aContig.getContigLength() - 1;
  pos = oldContigSize;

  ProgressIndicator<int32> pi(0, oldContigSize);

#ifndef RUN_ONLY
  if (p.isVerbose(0)) {
    cout << "\n========================= " << endl;
    cout << "New Contig (Nr. " << ++contignr << ")" << endl;
    cout << "Contig Length : " << pos << endl;
    cout << "========================= " << endl;

    if (p.isVerbose(9)) {
      cout << aContig;
    }
  }
#endif

  if (false == p.doContig(regionnr)) {
#ifndef RUN_ONLY
    if (p.isVerbose()) {
      cout << "Contig skipped... " << endl;
    }
  }
#endif


  if (pos < MIN_CONTIG_EDIT_LENGTH) {
#ifndef RUN_ONLY
    if (p.isVerbose()) {
      cout << "Contig skipped: not long enough" << endl;
    }
#endif

    return ;
  }


  while (I_CONS != theConsensus.begin() && p.doHypo()) {
    pos--;
    I_CONS = findPrevMismatch(aContig, pos, theConsensus);

    if (I_CONS == theConsensus.begin())
      break; // exit while-loop

    try {

#ifndef RUN_ONLY
      if (p.isVerbose()) {
	cout << "FR: Mismatch " << ++regionnr << " found at: " << pos << endl;
      }
#endif

      if (beginRegion == 0) {
	beginRegion =
	  pos + extendFaultRegionLeft(I_CONS, theConsensus.begin()) - 1;
	endRegion =
	  pos + extendFaultRegionRight(I_CONS, theConsensus.end()) - 1;
      } else {
	int32 rightBorder = pos +
	  extendFaultRegionRight(I_CONS, theConsensus.end()) - 1;

	if (rightBorder < beginRegion) {

	  /*
	    cout << "\nCALL SMALL REGIONS " << beginRegion << "\t"
	    << endRegion << endl;

	    editSmallRegionsBack(aContig, p, beginRegion, endRegion);
	    cout << "END SMALL REGIONS" << endl;
	  */
	  {
	    hypotheses_pool *aPool;
	    aPool = createFaultRegion(aContig, beginRegion, endRegion, p);

	    if (!examineFaultRegion(aContig, beginRegion, endRegion,
				    aPool, p)) {

	      editSmallRegionsBack(aContig, p, beginRegion, endRegion);
	    }
	  }

	  endRegion = pos +
	    extendFaultRegionRight(I_CONS, theConsensus.end()) - 1;

	  beginRegion= pos +
	    extendFaultRegionLeft(I_CONS, theConsensus.begin()) - 1;
      } else {

	  beginRegion=pos +
	    extendFaultRegionLeft(I_CONS,theConsensus.begin()) - 1;
	}
      }


    }
    catch (Notify n) {
      beginRegion = pos-1;
      endRegion = pos-1;
      n.handleError("EdIt Error: Unable to handle mismatch! (1)");
    }

    if (p.getShowProgress()) {
      pi.progress(oldContigSize - pos);
    }
  }

  // don't forget the last mismatch if there is one...
  if (beginRegion > 0) {
    try {
      hypotheses_pool *aPool;
      aPool = createFaultRegion(aContig, beginRegion, endRegion, p);

      examineFaultRegion(aContig, beginRegion, endRegion, aPool, p);
    }
    catch (Notify n) {
      n.handleError("EdIt Error: Unable to handle mismatch! (2)");
    }
  }

  if (p.getShowProgress()) {
    pi.progress(oldContigSize);
    cout << endl;
  }
}





/*********************************************
 *       editSmallRegionsBack                *
 *********************************************/


void editSmallRegionsBack(Contig &aContig, EDITParameters &p,
			  int32 von=-1, int32 bis=-1)
{
  static int32 regionnr = 0;
  int32 last_region_start = 0;
  int32 last_region_end   = 0;
  int32 ende = 0;

  Contig::cccontainer_t &theConsensus =
    const_cast<Contig::cccontainer_t&>(aContig.getConsensusCounts());
  Contig::cccontainer_t::iterator I_CONS = theConsensus.end();

  // set pos to the last base of the contig
  int32 pos = 0;

  if (bis < 0) {
    pos = aContig.getContigLength() - 1;
  } else {
    pos = bis;
  }

  if (von > 0) {
    ende = von;
  }


  while (I_CONS != theConsensus.begin() && pos > ende) {
    I_CONS = findPrevMismatch(aContig, pos, theConsensus);

    if (I_CONS == theConsensus.begin())
      break; // exit while-loop

    try {

#ifndef RUN_ONLY
      if (p.isVerbose()) {
	cout << "SR: Mismatch " << ++regionnr << " found at: " << pos << endl;
      }
#endif

      if (pos >= 2) {
	if (last_region_start == 0) {
	  last_region_start = pos;
	  last_region_end = pos;
	} else {
	  if (abs(pos - last_region_end) > 1) {
	    // Create region around lastpos
	    hypotheses_pool *aPool;

	    aPool = createFaultRegion(aContig, last_region_end-2,
				      last_region_start, p);

	    //	    examineFaultRegion(aContig, last_region_start-2,
	    //		       last_region_end, aPool, p);
	    examineFaultRegion(aContig, last_region_end-2,
			       last_region_start, aPool, p);

	    last_region_start = pos;
	    last_region_end = pos;
	  } else {
	    // extend the region
	    last_region_end = pos;
	  }
	}
      }

      if (p.isVerbose()) {
	cout << "   " << last_region_start << "\t" << last_region_end << endl;
      }
    }
    catch (Notify n) {
      pos--;
      n.handleError("EdIt Error: Unable to handle mismatch! (3)");
    }
    pos--;
  }

  if (last_region_end != 0 && pos >= ende) {
    hypotheses_pool *aPool;
    aPool = createFaultRegion(aContig, last_region_end-2,
			      last_region_start, p);

    examineFaultRegion(aContig, last_region_end-2,
		       last_region_start, aPool, p);
  }

}





/**********************************************
 *       editTowardsConsensusBack             *
 **********************************************/


void editTowardsConsensusBack(Contig &aContig, EDITParameters &p)
{
  static int32 regionnr = 0;

  Contig::cccontainer_t &theConsensus =
    const_cast<Contig::cccontainer_t&>(aContig.getConsensusCounts());
  Contig::cccontainer_t::iterator I_CONS = theConsensus.end();

  // set pos to the last base of the contig
  int32 pos = 0;
  pos = aContig.getContigLength() - 1;


  while (I_CONS != theConsensus.begin()) {
    I_CONS = findPrevMismatch(aContig, pos, theConsensus);

    if (I_CONS == theConsensus.begin())
      break; // exit while-loop

    try {

#ifndef RUN_ONLY
      if (p.isVerbose()) {
	cout << "C: Mismatch " << ++regionnr << " found at: " << pos << endl;
      }
#endif

      if (pos >= 2) {
	hypotheses_pool *aPool;
	aPool = createFaultRegion(aContig, pos-2, pos, p);

	// We only edit single columns, for several reasons (e.g. finding
	// the correct fault-class) we include adjacent columns. But we do
	// not edit the. A hypotheses with mismatches in these columns thus
	// can not be solved - we iron out these mismatches my brute force
        // in our hypotheses...

	aPool->setFirstAndLastColumn(calledBase(*(I_CONS-1)),
				     calledBase(*(I_CONS+1)));

	cout << "X5" << endl;
	examineFaultRegion(aContig, pos-2, pos, aPool, p);
      }

    }
    catch (Notify n) {
      n.handleError("EdIt Error: Unable to handle mismatch! (4)");
    }
    pos--;
  }
}



// *******************************************************
//    faultRegionFromRead:
//    create the scfInfo-structure for a fault region in
//    read 'R' (calling function must free the memory).
//    Copy the sequence of that fault region into
//    the paramter 'sequence'.
// *******************************************************


scfInfo* faultRegionFromRead(Contig &aContig, const Contig::contigread_t &R,
			     int32 readIndex, int32 von, int32 bis,
			     char *sequence)
{
  const vector<char>  &cseq  = R.read.getActualComplementSequence();
  const vector<char>  &seq   = R.read.getActualSequence();
  vector<char>::const_iterator I;
  //  vector<char>::const_iterator Ende;
  int32 readPos;
  bool  hiddenData = false;

  scfInfo *aScfInfo;

  readPos = readPosFromContigPos(von, R);

  if (R.direction > 0) {
    I = seq.begin() + readPos;
    //Ende = seq.end();
  } else {
    I = cseq.begin() + readPos;
    //Ende = cseq.end();
  }

  aScfInfo = new scfInfo;

  if (von > bis) {
    von ^= bis ^= von ^= bis;
  }

  aScfInfo->scfInfo_position.resize(bis - von + 3, -1);

  aScfInfo->scfInfo_dbpos = readPos;

  if ((von < R.offset || von >= R.offset + R.read.getLenClippedSeq()) &&
      (bis < R.offset || bis >= R.offset + R.read.getLenClippedSeq())) {
    hiddenData = true;
  }

  for (int32 loop = 0; loop <= bis - von; loop++) {
    if ((von + loop >= R.offset &&
	von + loop <  R.offset + R.read.getLenClippedSeq()) ||
	hiddenData) {
      sequence[loop] = toupper(*I);
      aScfInfo->scfInfo_position[loop] = scfPosFromDBPos(loop+readPos, R);
    } else {
      sequence[loop] = EOR_CHARACTER;
      aScfInfo->scfInfo_position[loop] = -1;
    }

    //if (von + loop <= R.offset + R.read.getLenClippedSeq()) {
    //if (I != Ende)
      I++;
    //}
  }

  sequence[bis - von + 1] = 0;

  aScfInfo->scfInfo_name = new char [R.read.getName().size() + 1];
  strcpy(aScfInfo->scfInfo_name, R.read.getName().c_str());
  aScfInfo->scfInfo_read = &R;
  aScfInfo->scfInfo_contig = &aContig;
  aScfInfo->scfInfo_readId = readIndex;
  aScfInfo->scfInfo_hidden = hiddenData;

  return aScfInfo;
}



// *******************************************************
//    createRegion:
//    create a faultRegion from 'von' to 'bis' in the
//    given Contig and initialise a hypotheses_pool with
//    this faultRegion. Include reads that have some valid
//    data between 'von' and 'bis'.
// *******************************************************


hypotheses_pool* createRegion(Contig &aContig, int32 von, int32 bis,
			      EDITParameters &p)
{
  scfInfo *aScfInfo;
  char sequence[100];  // FIXME 6
  const vector<Contig::contigread_t> &allReads = aContig.getContigReads();
  vector<int32> theActiveReads;
  vector<int32>::iterator I;
  vector<scfInfo> scfInfoList;
  fault_region    *aFaultRegion;
  hypotheses_pool *aPool;

  if (bis - von > 30) return nullptr;
  if (von < 0)  von = 0;
  if (bis >= aContig.getConsensusCounts().size()) { bis=aContig.getConsensusCounts().size()-1; }


  if (p.isVerbose(9)) {
    cout << "createRegion " << von << "-" << bis << endl;
  }

  {
    // create union of the reads at the begin and at the end of
    // the region in 'theActiveReads'
    vector<int32>  activeReads_von, activeReads_bis;
    less<int> vergleich;

    aContig.getCRIndicesAtContigPosition(activeReads_von, von);
    sort(activeReads_von.begin(), activeReads_von.end(), vergleich);
    aContig.getCRIndicesAtContigPosition(activeReads_bis, bis);
    sort(activeReads_bis.begin(), activeReads_bis.end(), vergleich);

    theActiveReads.reserve(activeReads_von.size()+activeReads_bis.size());
    set_union(activeReads_von.begin(), activeReads_von.end(),
	      activeReads_bis.begin(), activeReads_bis.end(),
	      inserter(theActiveReads, theActiveReads.begin()));
  }


  aFaultRegion = new fault_region;
  scfInfoList.clear();

  I = theActiveReads.begin();

  while (I != theActiveReads.end()) {
    const Contig::contigread_t &r = allReads[*I];
    if(r.read.usesAdjustments()){
      int32 readIndex = theActiveReads[I - theActiveReads.begin()];
      int32 readPos;

      aScfInfo = faultRegionFromRead(aContig, r, readIndex, von,
				     bis, sequence);

      aFaultRegion->addSequence(sequence, *aScfInfo);
      scfInfoList.push_back(*aScfInfo);

      delete aScfInfo;
    }
    I++;
  }

  aPool = new hypotheses_pool(*aFaultRegion, scfInfoList, &p);
  delete aFaultRegion;

  return aPool;
}




// *******************************************************
//    testADSforReadExtend:
//    examine a sequence alignment and look, if is can be
//    regarded as corresponing to the consensus with
//    sufficient quality
// *******************************************************


bool testADSforReadExtend(AlignedDualSeq &ads, int32 umgebung, int32 &offset)
{
  const char *s1 = ads.getAlignedSequence(1);
  const char *s2 = ads.getAlignedSequence(2);
  int32 verschiebung;
  int32 pos;

  // If we extend the read (we always extend from left to right!)
  // the begin should match quite good (without offset).

  if (ads.getOffsetInAlignment(1) != 0 ||
      ads.getOffsetInAlignment(2) != 0)
    {
      return false;
    }

  /*
  if (ads.getOffsetInAlignment(1) > umgebung ||
      ads.getOffsetInAlignment(2) > umgebung)
    {
      return false;
    }
  */

  if (ads.getRightOffsetInAlignment(1) == 0) {

    verschiebung = ads.getRightOffsetInAlignment(2);

#ifdef HYPOTHESES_VERBOSE
    cout << "Fall 1. Offset:" << offset << "  Verschiebung: "
    	 << verschiebung << endl;
#endif

    offset = - verschiebung; // - offset;
  } else {
    verschiebung = ads.getRightOffsetInAlignment(1);

#ifdef HYPOTHESES_VERBOSE
    cout << "Fall 2. Offset:" << offset << "  Verschiebung: "
	 << verschiebung << endl;
#endif

    offset = verschiebung - offset;
  }

  // the alignment should have a certain quality

  if (ads.getScoreRatio() < 85) {
    return false;
  }

  pos = ads.getTotalLen() - verschiebung - 1;

  //  the reads must not differ some bases before the region of interest

  for (int32 loop = 1; loop <= 2 * umgebung + 1; loop++) {
    if(s1[pos] != s2[pos] &&
       toupper(s1[pos]) != 'N' &&
       toupper(s2[pos] != 'N')) {
      return false;
    }
    pos--;
  }


  return true;
}




// *******************************************************
//    ExtendRead:
//    Examine a single read if it fits the consensus and
//    thus can be extended to the given number of bases
// *******************************************************


bool extendRead(const Contig &aContig, const Contig::contigread_t &aRead,
		int32 dbPos, int32 anzBases,
		int32 umgebung, int32 &offset_to_readpos)
{
  Contig::cccontainer_t& theConsensus =
    const_cast<Contig::cccontainer_t&>(aContig.getConsensusCounts());
  Contig::cccontainer_t::iterator I_CONS = theConsensus.begin();
  const Read &R = aRead.read;
  const vector<char> &csequence = R.getActualComplementSequence();
  const vector<char> &sequence =  R.getActualSequence();

  MIRAParameters mira_parameters;

  const_cast<align_parameters &>(mira_parameters.getAlignParams()).al_min_score=3;
  const_cast<align_parameters &>(mira_parameters.getAlignParams()).al_min_overlap=3;

  int32 readPos;
  int32 vonPos, bisPos;
  bool  result = false;


  // Find the positions in the read, where we have to look at and
  // initialize I_CONS to the corresponding base in the consensus

  readPos = readPosFromContigPos(dbPos, aRead);

#ifdef HYPOTHESES_VERBOSE
  cout << "Extend Read;  direction: " << aRead.direction << endl;
  cout << "dbpos:" << dbPos << "  readpos:" << readPos;
  cout << "  anzBases " << anzBases;
  cout << "  offset:" << aRead.offset << endl;
#endif

  if (dbPos > aRead.offset) {
    // We extend the read on the right
    I_CONS = I_CONS + dbPos - anzBases;
    vonPos = readPos - anzBases;
    bisPos = readPos - 1 + umgebung;
  } else {
    // We extend the read on the left
    I_CONS = I_CONS + dbPos +anzBases - 1;
    vonPos = readPos - umgebung;
    bisPos = readPos + anzBases - 1;
  }

  {
    // Copy the part of the read, we want to extend and the corresponding
    // part of the consensus (the base at dbPos is always on the rear)
    // and align them. We copy 'umgebung' bases more!

    Align myAlignment(&mira_parameters);
    char *read, *cons;
    list<AlignedDualSeq> adslist;
    list<AlignedDualSeq>::iterator I_ADS;

    read = new char [anzBases + umgebung + 1];
    cons = new char [anzBases + umgebung + 1];


    for (int32 loop = 0; loop < anzBases + umgebung; loop++) {
      if (dbPos > aRead.offset) {

        if (I_CONS != theConsensus.end()) {
          if (aRead.direction > 0) {
	    read[loop] = sequence[vonPos++];
	  } else {
	    read[loop] = csequence[vonPos++];
	  }
  	  cons[loop] = calledBase(*I_CONS);
	  I_CONS++;
        } else {
          read[loop] = 'N';
          cons[loop] = 'N';
        }

      } else {
	if (aRead.direction > 0) {
	  read[loop] = sequence[bisPos--];
	} else {
	  read[loop] = csequence[bisPos--];
	}
        cons[loop] = calledBase(*I_CONS);

        if (I_CONS != theConsensus.begin()) {
          I_CONS--;
        } else {
	  while (loop < anzBases + umgebung) {
            read[loop] = 'N';
            cons[loop] = 'N';
            loop++;
          }
        }
      }
    }


    read[anzBases + umgebung] = 0;
    cons[anzBases + umgebung] = 0;

    // Overlap is less than 5 bases long
    // we do not use the alignment routines for these cases

    if (anzBases + umgebung <= 5) {
      bool  ok = true;
      bool  undef = false;
      int32 l = 0;

      while (l < anzBases + umgebung && ok) {
	if (read[l] != cons[l]) {
	  if (undef) {
	    ok = false;
	  } else {
	    ok = (toupper(read[l]) == 'N' || toupper(cons[l]) == 'N');
	    undef = ok;
	  }
	}
	l++;
      }
      result = ok;

    } else {

      myAlignment.acquireSequences(read, anzBases + umgebung,
				   cons, anzBases + umgebung,
				   1, 2, // ReadIDs
				   1, 1, // Richtung
				   true, // calculate with offset
				   0     // expected offset
				   );

      myAlignment.fullAlign(&adslist,false,false);

      // we step through the list of possible alignments and search for
      // if one of them satisfies outr conditions

      I_ADS = adslist.begin();
      while (I_ADS != adslist.end() && !result) {
	result = testADSforReadExtend(*I_ADS, umgebung, offset_to_readpos);

	I_ADS++;
      }

      // turn the offset for reverse reads
      // if (aRead.direction < 0) offset_to_readpos = -offset_to_readpos;
    }

    delete [] read;
    delete [] cons;
  }
  return result;
}



// Bach: 30.07.2009
// check usesAdjustments()
int32 scfPosFromDBPos(const int32 pos, const Contig::contigread_t &aRead)
{
  const vector<int32> &adjust = aRead.read.getAdjustments();

  if (aRead.direction < 0) {
    if(aRead.read.usesAdjustments()) return adjust[aRead.read.getLenSeq() - 1 - pos];
    return aRead.read.getLenSeq() - 1 - pos;
  } else {
    if(aRead.read.usesAdjustments()) return adjust[pos];
    return pos;
  }
}



int32 scfInfo_contigPosFirstBase(scfInfo &si)
{
  return dbPosFromReadPos(si.scfInfo_dbpos, *(si.scfInfo_read));
}



// *******************************************************************
//  dbPosFromReadPos
//  We have a position in a read and we want to know the corresponding
//  position in the contig
//
//  readPosFromContigPos
//  We have a contig position and we want to know which basePos in the
//  read corresponds to this position
//
//
//  CONTIG                       Offset
//  ==============================+
//                 *--------------|------x--------+
//                 |<--leftClip-->|      readPos
//
//
//  CONTIG                          Offset
//  ================================+
//            *---------------------|-------x-------+
//            |<---len-rightClip--->|       readpos
//
//  ******************************************************************


int32 dbPosFromReadPos(const int32 readPos, const Contig::contigread_t &cRead)
{

  if (cRead.direction > 0) {
    return readPos + cRead.offset - cRead.read.getLeftClipoff();
  } else {
    return readPos + cRead.offset -
      (cRead.read.getLenSeq() -  cRead.read.getRightClipoff());
  }
}


int32 readPosFromContigPos(const int32 contigPos, const Contig::contigread_t &cRead)
{
  if (cRead.direction > 0) {
    return contigPos - cRead.offset + cRead.read.getLeftClipoff();
  } else {
    return contigPos - cRead.offset +
      (cRead.read.getLenSeq() - cRead.read.getRightClipoff());
  }
}
