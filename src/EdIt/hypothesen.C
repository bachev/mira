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
 * hypotheses generation
 *
 * Written by Thomas Pfisterer
 * Version 0.92       13.08.1999
 *         0.93       02.09.1999
 *         0.94       19.10.1999
 *
 */


#include "hypothesen.H"
#include "assert.h"


// TODO
// Was ist, wenn am Ende doch eine Read in einer Richtung fehlt (z.B. Quality
// war schlecht) und eigentlich sind welche zum Verlängern da.
// => Lösung(?): Nimm alle verlängerbaren Reads mit auf?
//
// Syntaktische Operationen evtl. mit eigenenm Operationstyp fr_insert_syntax
// und dann kein Aufruf von SCF_look und auch eigener penalty.
//
// Zwei Teilhypothesen werden nach dem einfachen Abstandskriterium getrennt.
// Hier sollten Repeats berücksichtigt werden
// ok. 07-09-97
//
// Mehrfaches Laden der gleichen SCF-Dateien bei mehrspaltigen Hypothesen
// vermeiden. Editieren auf den geladenen Hypothesen!
//
// 02.09 set tag in consensus if size has changed
// 19.10 friend funktionen zur ausgabe


//
//  Initialize new fault_region
//

fault_region::fault_region() {
  fr_sequenceList.clear();
  fr_sequenceList.reserve(5);

  fr_edits.clear();
  fr_edits.reserve(4);

  fr_steps = 0;
  fr_longest = 0;
  fr_score = 0;
  fr_solution_score = 0;
  fr_has_eor = false;
}



fault_region::fault_region(fault_region const &other)
{
  *this = other;
}



fault_region::~fault_region()
{
   clear();
}



void fault_region::clear()
{
  fr_sequenceList.clear();
  fr_edits.clear();

  fr_longest = 0;
  fr_score = 0;
  fr_steps = 0;
  fr_has_eor = false;
}



const fault_region& fault_region::operator= (const fault_region &other)
{
  if (this == &other) {
    return *this;
  }

  fr_sequenceList.clear();
  fr_edits.clear();

  // copy the object
  fr_longest = other.fr_longest;
  fr_score   = other.fr_score;
  fr_steps   = other.fr_steps;
  fr_edits   = other.fr_edits;
  fr_solution_score = other.fr_solution_score;
  fr_sequenceList = other.fr_sequenceList;
  fr_has_eor      = other.fr_has_eor;

  // FIXME
  //  scfEditHypotheses = ???

  return other;
}




char fault_region::readBase(const int32 readId, const int32 readPos) const
{
  return fr_sequenceList[readId].si_sequence[readPos];
}



int32 fault_region::readSequenceCount() const
{
  return fr_sequenceList.size();
}



int32 fault_region::readScore() const
{
  return fr_score;
}



vector<editoperation> &fault_region::readEditOperations()
{
  return fr_edits;
}



const sequenceInfo& fault_region::readSequence(int32 sequenceId)
{
  return fr_sequenceList[sequenceId];
}



inline int32 fault_region::size() const
{
  return fr_longest;
}


int32 fault_region::realSize() const
{
  int32 size = 0;

  for (int32 l=0; l < fr_longest; l++) {
    vector<sequenceInfo>::const_iterator I = fr_sequenceList.begin();
    while (I != fr_sequenceList.end()) {
      if (I->si_sequence[l] != fr_eorchar &&
	  I->si_sequence[l] != fr_fillchar) {
	size++;
	break;
      }
      I++;
    }
  }
  return size;
}


inline int32 fault_region::readSolutionScore() const
{
  return fr_solution_score;
}



void fault_region::setColumnToBase(int32 col, char base)
{
 vector<sequenceInfo>::iterator I = fr_sequenceList.begin();

 // cout << "SetColumnToBase " << col << " base: " << base << endl;
 while (I != fr_sequenceList.end()) {
   if (I->si_sequence[col] != fr_eorchar) {
     (I->si_sequence)[col] = base;
   }
   I++;
 }
}


// *************************************************************************
// *
// *  Search the fr_sequenceList if the sequence `target' is in the list
// *  The sequences may differ in fr_fillchar characters at the end.
// *
// *  return:  vector<sequenceInfo>::iterator to fr_sequenceList entry
// *
// *************************************************************************

vector<sequenceInfo>::iterator fault_region::findSequence(const char* target)
{
  vector<sequenceInfo>::iterator iter = fr_sequenceList.begin();
  int32  pos;

  if (target == nullptr) {
    throw Notify(Notify::FATAL, "findSequence: No valid Sequence (nullptr)");
  }

  while (iter != fr_sequenceList.end()) {
    vector<char>::iterator I = (*iter).si_sequence.begin();
    pos = 0;

    // skip over the identical part....

    while (I != (*iter).si_sequence.end() &&
	   (*I) == target[pos] &&
	   target[pos] != 0) {
      I++;
      pos++;
    }

    // they are really identical...
    // FIXME 1
    // if (*I == target[pos] && target[pos+1] == 0) return iter;
    //      pos+1 may access undefined memory

    if (I == (*iter).si_sequence.end() &&
	target[pos] == 0)
      {
	return iter;
      }

    // perhaps 'I' has only some additional fillchars...
    if (target[pos] == 0) {
      while (I != iter->si_sequence.end() && (*I) == fr_fillchar) I++;
      if (I == iter->si_sequence.end()) return iter;
    }

    iter++;
  }

  return iter;
}


// ***************************************************************************
// *
// *  Add sequence 's' to the fault region.
// *  Fill out scfInfo to link the fault region to the scf-file and the contig
// *
// *  return:     1  sequence found
// *              0  new sequence added to fault region
// *
// ***************************************************************************

int32 fault_region::addSequence(const char* s, scfInfo &aScfInfo)
{
  int32  newSize = strlen(s);

  // Search for sequence in the fault_region
  vector<sequenceInfo>::iterator I = findSequence(s);


  // the same sequence is already present in the fr_sequenceList
  // scfInfo must not be deleted in this case, it links the 'real' sequence
  // to the sequence in the fault_region

  if (I != fr_sequenceList.end()) {
    I->si_occurrences++;
    aScfInfo.scfInfo_sequenceId = I - fr_sequenceList.begin();
    return 1;
  }


  // create a new sequenceInfo and put it in the fr_sequenceList
  // link the scfInfo to the entry in fr_sequenceList
  {
    sequenceInfo seqInfo;
    seqInfo.si_occurrences = 1;
    seqInfo.si_insertCount = 0;
    seqInfo.si_lastEditPos = -1;
    seqInfo.si_sequence.reserve(newSize+2);   // to speed up inserts
    seqInfo.si_sequence.clear();

    for (int a = 0; a < newSize; a++) {
      fr_has_eor = fr_has_eor || (s[a] == fr_eorchar);
      seqInfo.si_sequence.push_back(toupper(s[a]));
    }

    fr_sequenceList.push_back(seqInfo);
    aScfInfo.scfInfo_sequenceId = fr_sequenceList.size() - 1;
    aScfInfo.scfInfo_confirmed = false;
  }


  // resize all entries in the list to the length of the longes sequence
  if (fr_sequenceList.size() == 0) {
    fr_longest = newSize;
  } else {
    if (newSize != fr_longest)
      resizeRegion(abs(newSize - fr_longest));
  }

  fr_score = pairwiseMismatch();

  return 0;
}



// ***************************************************************************
// *
// * Is this solution possible, given the scope of the current system?
// * Not used at the moment - all operations are possible! ;-)
// *
// ***************************************************************************


bool fault_region::possibleSolution() const
{
  return true;
}


// ***************************************************************************
// *
// * Is this solution valid. We do not allow solutions that insert a base
// * and delete the adjacent base afterwards.
// *
// ***************************************************************************


bool fault_region::validSolution() const
{

  vector<editoperation>::const_iterator I = fr_edits.begin();
  vector<editoperation>::const_iterator J;

  while (I != fr_edits.end()) {
    if (I->base_old == fr_fillchar && isRealBase(I->base_new)) {
      // insert a base...

      J = fr_edits.begin();
      while (J != fr_edits.end()) {
	if (J->base_new == fr_fillchar && isRealBase(J->base_old)) {
	  // and delete another base
	  if ( I->sequenceId == J->sequenceId &&
	      ( J->position   == I->position+1 ||
	        J->position+1 == I->position))
	    {
	      // in the same sequence and in adjacent positions
	      return false;
	    }
	}
	J++;
      }
    }
    I++;
  }

  if (fr_has_eor == false)
    return true;


  vector<sequenceInfo>::const_iterator S = fr_sequenceList.begin();
  vector<char>::const_iterator C;

  while (S != fr_sequenceList.end()) {
    C = S->si_sequence.begin();
    C++;  // next position...

    while (C+1 < S->si_sequence.end()) {
      if (*C == fr_eorchar &&
	  *(C-1) != fr_eorchar && *(C-1) != fr_fillchar &&
	  *(C+1) != fr_eorchar && *(C+1) != fr_fillchar)  {
	// cout << "Not valid! 2: " << *(C-1) << *C << *(C+1) << endl;
	return false;
      }

      // FIXME 4
      if (*C == fr_eorchar && *(C+1) == fr_gapchar) {
      	return false;
      }

      C++;
    }

    C = S->si_sequence.end() - 1;
    if (*C != fr_eorchar && *C != fr_fillchar && *(C-1) == fr_eorchar) {
      // cout << "Not valid! 3: " << *C << *(C-1) << endl;
       return false;
    }
    S++;
  }
  return true;
}





// ***************************************************************************
// *
// *  resize all sequences in the fault region
// *    increment 'fr_longest' by 'size'
// *    resize sequences to fr_longest
// *
// *   return:  -1 error / 0 ok
// ***************************************************************************


int32 fault_region::resizeRegion(int32 size)
{
  if (size < 1) {
    cerr << "resizeRegion: size<0 not implemented!\n";
    return -1;
  }

  fr_longest += size;
  for (uint32 l = 0; l < fr_sequenceList.size(); l++) {
    fr_sequenceList[l].si_sequence.resize(fr_longest, fr_fillchar);
  }
  return 0;
}


// ***************************************************************************
// *
// *  calculate number of pairwise mismatches in fault_region
// *  return:   number of pairwise mismatches
// *
// ***************************************************************************

int32 fault_region::pairwiseMismatch()
{
  int32 length_loop;
  int32 mismatches = 0;

  if (fr_sequenceList.size() < 2) return 0;

  for (length_loop = 1; length_loop < fr_longest; length_loop++) {
     mismatches += pairwiseMismatchColumn(length_loop);
  }

  return mismatches;
}


// ***************************************************************************
// *
// *  calculate number of pairwise mismatches in a single column
// *  return:   number of pairwise mismatches
// *
// ***************************************************************************


int32 fault_region::pairwiseMismatchColumn(int32 col) {
  int32 mismatches = 0;
  int32 s = fr_sequenceList.size();

  for (int32 l1=0; l1 < s; l1++) {
    vector<char> &seq1 = fr_sequenceList[l1].si_sequence;

    for (int32 l2=l1+1; l2 < s; l2++) {
      vector<char> &seq2 = fr_sequenceList[l2].si_sequence;
      if (seq1[col] != seq2[col] && seq1[col] != EOR_CHARACTER &&
	  seq2[col] != EOR_CHARACTER)
	mismatches++;
      //        mismatches += fr_sequenceList[l1].si_occurrences *
      //  fr_sequenceList[l2].si_occurrences;
    }
  }
  return mismatches;
}


// ***************************************************************************
// *
// *  to keep the search space compact we try to merge identical fault regions
// *  this is a criterion when regions are treated as identical
// *
// *  result  0  'identical'
// *          1  'different'
// ***************************************************************************

int32 fault_region::compare(fault_region &fr)
{
  // different score
  if (fr.readScore() != readScore())
    return 1;

  // different number of edits
  if (fr.fr_steps != fr_steps || fr.fr_steps >= fr_maxsteps)
    return 1;

  // different length ?????  -> fr_fillchar
  // if (fr.fr_longest != fr_longest)
  //   return 1;

  if (calcSolutionScore() != fr.calcSolutionScore())
    return 1;

  vector<sequenceInfo>::iterator iter1 = fr_sequenceList.begin();
  vector<sequenceInfo>::iterator iter2 = fr.fr_sequenceList.begin();

  while (iter1 != fr_sequenceList.end()) {
    vector<char> &seq1 = (*iter1).si_sequence;
    vector<char> &seq2 = (*iter2).si_sequence;
    vector<char>::iterator I = seq1.begin();
    vector<char>::iterator J = seq2.begin();

    // the first character is alway identical loop=0 => loop=1

    while(I != seq1.end() && J != seq2.end()) {
      if (*I != *J) return 1;
      I++;
      J++;
    }

    iter1++;
    iter2++;
  }
  return 0;
}



// ************************************************************
// *
// * check for a valid readId and readPos
// * return: 0 ok;
// *         1 nonsens operation
// *         throw an exception otherwise
// *
// ************************************************************


int32 fault_region::checkParameters(const int32 readId, const int32 readPos)
{
  if (readId > (int32)fr_sequenceList.size()) {
    throw Notify(Notify::WARNING,
		 "fault_region -- readId out of range");
  }

  if(readPos > fr_longest) {
    //    cout << "ReadPos   : " << readPos << endl;
    //    cout << "fr_longest: " << fr_longest << endl;
#ifdef LIGHTHEADED
    return 1;
#else
    throw Notify(Notify::WARNING,
    		 "fault_region -- readPos is out of range");
#endif
  }

  if (fr_sequenceList[readId].si_lastEditPos == readPos) {
    return 1;  // nonsens operation
  }

  return 0;
}



// ************************************************************
// *
// * insertBase in a given read
// *
// ************************************************************


int32 fault_region::insertBase(const int32 readId,
	  		       const int32 readPos, const char base)
{
  vector<char> &aSequence = fr_sequenceList[readId].si_sequence;

  if (checkParameters(readId, readPos) == 1) return 1;

  // last base is not a fillchar, we must extend the region
  if (aSequence[fr_longest-1] != fr_fillchar) {
    resizeRegion(1);
    aSequence = fr_sequenceList[readId].si_sequence;
  }


  // goto readpos, insert the base and delete resize the region
  // this is why we checked for the fillchar
  {
    vector<char>::iterator I =aSequence.begin();
    advance(I, readPos);
    aSequence.insert(I, base);
    aSequence.resize(fr_longest);
  }


  // note what we have edited in op and add it it to the list of
  // editoperations made
  {
    editoperation op;
    op.sequenceId = readId;
    op.position  = readPos;
    op.base_old  = fr_fillchar;
    op.base_new  = base;

    if (op.base_new == fr_gapchar) {
      op.operation = fr_syntactic;
    } else {
      op.operation = fr_insert;
    }

    fr_edits.push_back(op);
  }

  fr_score = pairwiseMismatch();
  fr_steps++;

  // mark where we have made the last edit
  fr_sequenceList[readId].si_lastEditPos = readPos;

  return 0;
}


// ************************************************************
// *
// * deleteBase in a given read
// *
// ************************************************************

int32 fault_region::deleteBase(const int readId, const int readPos)
{
  vector<char> &aSequence = fr_sequenceList[readId].si_sequence;

  if (checkParameters(readId, readPos) == 1) return 1;

  {
    // keep track of our edit-operations
    editoperation op;
    op.sequenceId = readId;
    op.position   = readPos;
    op.base_old   = aSequence[readPos];
    op.base_new   = fr_fillchar;

    if (op.base_old == fr_gapchar) {
      op.operation = fr_syntactic;
    } else {
      op.operation = fr_delete;
    }

    fr_edits.push_back(op);
  }

  {
    // goto readPos, delete the base, resize the region
    vector<char>::iterator I = aSequence.begin();
    advance(I, readPos);
    aSequence.erase(I);
    aSequence.resize(fr_longest, fr_fillchar);
  }

  fr_score = pairwiseMismatch();
  fr_steps++;

  fr_sequenceList[readId].si_lastEditPos = readPos;

  return 0;
}


// ************************************************************
// *
// * alterBase in a given read
// *
// ************************************************************


int32 fault_region::alterBase(const int32 readId,
			      const int32 readPos, const char base)
{
  vector<char> &aSequence = fr_sequenceList[readId].si_sequence;

  if (checkParameters(readId, readPos) == 1) return 1;

  // If we want to alter a gapchar into a base and there is the required
  // base on the left of some gaps it is possible to get the base by
  // shuffling the gaps to the right and the base to the left.
  // e.g.
  //      XXXX**AXXXX    ==>  XXXXA**XXXX
  //       A--^
  // if we need another A on the left we can alter a gap to an A later


  if (aSequence[readPos] == fr_gapchar) {
    int32 x = readPos;

    while (aSequence[x] == fr_gapchar && x+1 < fr_longest) { x++; }

    if (aSequence[x] == base && x < fr_longest) {

      //    cout << readId << " x " << x << " readPos" <<readPos<<" " << endl;
      while (x > readPos) {
	// We produce some new edit-operations which are only syntactic
	deleteBase(readId, x-1);
	insertBase(readId, x, fr_gapchar);
	x--;
      }

    return 0;
    }
  }


  // keep track of what we did...
  {
    editoperation op;
    op.sequenceId = readId;
    op.position  = readPos;
    op.base_old  = aSequence[readPos];
    op.base_new  = base;

    if (op.base_new == fr_gapchar) {
      if (op.base_old == fr_fillchar) {
	op.operation = fr_syntactic;
      } else {
	op.operation = fr_delete;
      }
    } else {
      if (op.base_old == fr_gapchar || op.base_old == fr_fillchar) {
	op.operation = fr_insert;
      } else {
	op.operation = fr_alter;
      }
    }
    fr_edits.push_back(op);
  }

  // do the edit
  aSequence[readPos] = base;

  // recalculate some parameters
  fr_score = pairwiseMismatch();
  fr_steps++;
  fr_sequenceList[readId].si_lastEditPos = readPos;

  return 0;
}


// ******************************************************************
// *
// * calculate number of edits to perform the edits on all scf-reads.
// *
// *******************************************************************

int32 fault_region::editCount()
{
  int32 editCount = 0;

  for (uint32 i = 0; i < fr_edits.size(); i++) {
    editoperation &op = fr_edits[i];

    if (op.operation == fr_delete || op.operation == fr_insert ||
	op.operation == fr_alter) {
      editCount += fr_sequenceList[op.sequenceId].si_occurrences;
    }
  }

  return editCount;
}


// ******************************************************************
// *
// * calculate the score of a solution.
// * all operations are rated by a certain penalty value. If the
// * solution changes the size of the fault region a small additional
// * 'unbalance_penalty' is added. This is to prevent unnecessary
// * changes in the fault_region size by syntactical operations.
// *
// *******************************************************************

int32 fault_region::calcSolutionScore()
{
  fr_solution_score = 0;

  const int32 delete_base_penalty = 12;
  const int32 alter_base_penalty  = 12;
  const int32 alter_n_penalty     = 7;
  const int32 insert_base_penalty = 12;
  const int32 insert_pad_penalty  = 2;
  const int32 unbalance_penalty   = 2;
  const int32 syntactic_penalty   = 1;
  int32 balance = 0;

  for (uint32 i = 0; i < fr_edits.size(); i++) {
    editoperation &op = fr_edits[i];
    int32 k = fr_sequenceList[op.sequenceId].si_occurrences;

    switch (op.operation) {
    case fr_delete:
      if (op.base_old != fr_gapchar) {
	fr_solution_score += (k * delete_base_penalty);
      } else {
	fr_solution_score -= insert_pad_penalty;
      }

      if (op.base_new == '.') {
	balance -= k;
      }

      break;

    case fr_insert:
      if (op.base_new != fr_gapchar) {
	fr_solution_score += (k * insert_base_penalty);
      } else {
        fr_solution_score += (k * insert_pad_penalty);
      }

      if (op.base_old == fr_fillchar) {
	balance += k;
      }
      break;

    case fr_alter:
      if (isUndefinedBase(op.base_old)) {
	fr_solution_score += (k * alter_base_penalty);
      } else {
	fr_solution_score += (k * alter_n_penalty);
      }
      break;

    case fr_syntactic:
      fr_solution_score += k * syntactic_penalty;
      if (op.base_old == fr_gapchar && op.base_new == fr_fillchar)
	balance -= k;
      if (op.base_old == fr_fillchar && op.base_new == fr_gapchar)
	balance += k;
      break;
    }
  }

  fr_solution_score += abs(balance * unbalance_penalty);
  return fr_solution_score;
}



// ******************************************************************
// *
// * Perform an editoperation on the fault region
// * This is necessary to do the edits step by step and look at the
// * intermediate states
// *
// *******************************************************************


void fault_region::performEditOperationOnRegion(editoperation &edOp)
{
  if (edOp.base_new == fr_fillchar) {
    deleteBase(edOp.sequenceId, edOp.position);
    fr_sequenceList[edOp.sequenceId].si_insertCount--;
  } else {
    if (edOp.base_old == fr_fillchar) {
      insertBase(edOp.sequenceId, edOp.position, edOp.base_new);
      fr_sequenceList[edOp.sequenceId].si_insertCount++;
    } else {
      alterBase(edOp.sequenceId, edOp.position, edOp.base_new);
    }
  }
}


// ******************************************************************
// *
// * PerformOperationsOnReads
// * The operations of the fault_region are performed on the reads.
// * The changes are tagged. Syntactic operations are als performed
// * but not tagged.
// *
// *******************************************************************


void fault_region::performOperationsOnReads(vector<scfInfo> aScfInfo,
					    const EDITParameters &p,
					    int32 locomotion)
{
  int32 seq_id;
  int32 ins_count = 0;   // count insert operations
  int32 del_count = 0;   // count delete operations
  int32 sum_loco = 0;

  vector<scfInfo>::iterator I;
  editoperation op;
  int32  contigpos = -1;
  char   s[200];
  bool   setTag;
  Contig *c = nullptr;

  FUNCSTART("void fault_region::perform_operations(vector<scfInfo> aScfInfo");

  // Step through all editoperations for the region

  if (fr_has_eor == true && locomotion != 0) {
    // FIXME
    // At the moment we have problems handling such cases
    return;
  }


  for (uint32 i = 0; i < fr_edits.size(); i++) {
    op = fr_edits[i];
    seq_id = op.sequenceId;

    vector<afh_Info*>::iterator XX = fr_edits[i].readEdits.begin();

    while (XX != fr_edits[i].readEdits.end()) {

      I = aScfInfo.begin();
      advance(I, (*XX)->scfInfoId);  //new

      setTag = true;

      if (I->scfInfo_sequenceId == seq_id &&
	  I->scfInfo_dbpos > 0) {

	contigpos = scfInfo_contigPosFirstBase(*I) + op.position;

#ifndef RUN_ONLY
	if (p.isVerbose()) {
	  cout << I->scfInfo_name << " : " << I->scfInfo_dbpos << "\t ";
	}
#endif
	// Is the edit-operation in the visible part of a read.
	// For insert operations the visible part also includes the position
	// before the first base, otherwise it starts with the first base

	bool visible;

	if (op.operation == fr_insert) {
	  visible=
	    (contigpos < (int32)((I->scfInfo_read)->offset +
				 I->scfInfo_read->read.getLenClippedSeq())
	     &&
	     contigpos >= (int32)(I->scfInfo_read)->offset - 1);

	} else {
	  visible =
	    (contigpos < (int32)((I->scfInfo_read)->offset +
				 I->scfInfo_read->read.getLenClippedSeq())
	     &&
	     contigpos >= (int32)(I->scfInfo_read)->offset);
	}


	//if (p.isVerbose()) {
	//  cout << " visible: " << visible;
	//  cout << " confirmed: " << (*XX)->isConfirmed();
	//}

	if (visible && ((*XX)->isConfirmed() ||
			(*XX)->getFaultClass() == fhc_SYNTACTIC)){

	  c = I->scfInfo_contig;
          op.base_new = tolower(op.base_new);

	  if (op.base_new == fr_fillchar && op.operation == fr_syntactic) {
	    // Syntactic insert operation - do not tag
	    op.operation = fr_delete;
	    setTag = false;

#ifndef RUN_ONLY
	    if (p.isVerbose()) {
	      cout << "Op changed to fr_delete" << endl;
	    }
#endif
	  }
	  if (op.base_old == fr_fillchar && op.operation == fr_syntactic) {
	    // Syntactic delete operation - do not tag
	    op.operation = fr_insert;
	    setTag = false;

#ifndef RUN_ONLY
	    if (p.isVerbose()) {
	      cout << "Op changed to fr_insert" << endl;
	    }
#endif
	  }

	  switch (op.operation) {

	  case fr_alter:
            if (op.base_old == fr_fillchar) {

	      sprintf(s, "Insert Base at %d %c=>%c %s %d", contigpos,
		      op.base_old, op.base_new, I->scfInfo_name,
		      I->scfInfo_position[op.position]);

#ifndef RUN_ONLY
	      if (p.isVerbose()) {
		cout << s << endl;
	      }
#endif

	      c->insertBaseInRead((char)op.base_new, contigpos,
				  I->scfInfo_readId);
	      ins_count++;
	    } else {

	      sprintf(s, "Alter Base at %d %c=>%c %s %d", contigpos,
		      op.base_old, op.base_new, I->scfInfo_name,
		      I->scfInfo_position[op.position]);

#ifndef RUN_ONLY
	      if (p.isVerbose()) {
		cout << s << endl;
	      }
#endif

	      c->changeBaseInRead((char)op.base_new,
				  contigpos, I->scfInfo_readId);
	    }

	    if (setTag && p.tagAlterOperations()) {
	      fr_tagED_C.comment=multitag_t::newComment(s);
	      try {
		c->addTagToRead(contigpos, contigpos, I->scfInfo_readId,
				fr_tagED_C);
	      }
	      catch (Notify n) {
	      }
	      p.logOperation(s);
	    }

	    break;

	  case fr_delete: {
	    int32 readpos = I->scfInfo_position[op.position];
	    if (op.base_new == fr_fillchar) {

	      sprintf(s, "Delete Base at %d %c=>%c %s %d", contigpos,
		      op.base_old, op.base_new,I->scfInfo_name,
		      I->scfInfo_position[op.position]);  // op.position-1

#ifndef RUN_ONLY
	      if (p.isVerbose()) {
		cout << s << endl;
	      }
#endif

	      //int32 readpos = I->scfInfo_position[op.position-1];

	      c->deleteBaseInRead(contigpos, I->scfInfo_readId);
	      del_count++;

	      fr_tagED_D.comment=multitag_t::newComment(s);
	      if (setTag && p.tagDeleteOperations() &&
		  contigpos-1 >= I->scfInfo_read->offset &&
		  contigpos-1 <= I->scfInfo_read->offset +
		                 I->scfInfo_read->read.getLenClippedSeq()) {

		c->addTagToRead(contigpos-1, contigpos-1, I->scfInfo_readId,
				fr_tagED_D);
	      }
	    } else {

	      sprintf(s, "Change Base at %d %c=>%c %s %d", contigpos,
		      op.base_old, op.base_new,I->scfInfo_name, readpos);

#ifndef RUN_ONLY
	      if (p.isVerbose()) {
		cout << s << endl;
	      }
#endif

	      c->changeBaseInRead('*', contigpos, I->scfInfo_readId);

	      if (setTag && p.tagDeleteOperations()) {
		//int32 readpos = I->scfInfo_position[op.position-1];
		//int32 readpos = I->scfInfo_position[op.position];
		int32 left;

		fr_tagED_D.comment=multitag_t::newComment(s);
                if (contigpos-1 >= I->scfInfo_read->offset &&
		    contigpos <= I->scfInfo_read->offset +
		                 I->scfInfo_read->read.getLenClippedSeq()) {

		  c->addTagToRead(contigpos-1, contigpos, I->scfInfo_readId,
				  fr_tagED_D);
		} else {
		  c->addTagToRead(contigpos, contigpos, I->scfInfo_readId,
				  fr_tagED_D);
		}
	      }
	    }
	    p.logOperation(s);
	    break;
	  }

	  case fr_insert:
	    if (op.base_old == fr_fillchar) {

	      sprintf(s, "Insert Base at %d %c=>%c %s %d", contigpos,
		      op.base_old, op.base_new,I->scfInfo_name,
		      I->scfInfo_position[op.position]);

#ifndef RUN_ONLY
	      if (p.isVerbose()) {
		cout << s << endl;
	      }
#endif

		c->insertBaseInRead((char)op.base_new, contigpos,
				    I->scfInfo_readId);

	      ins_count++;
	    } else {

	      sprintf(s, "Alter Base at %d %c=>%c %s %d", contigpos,
		      op.base_old, op.base_new,I->scfInfo_name,
		      I->scfInfo_position[op.position]);

#ifndef RUN_ONLY
	      if (p.isVerbose()) {
		cout << s << endl;
	      }
#endif

	      c->changeBaseInRead((char)op.base_new,
				  contigpos, I->scfInfo_readId);
	    }

	    if (setTag && p.tagInsertOperations()) {
	      fr_tagED_I.comment=multitag_t::newComment(s);
	      c->addTagToRead(contigpos, contigpos,  I->scfInfo_readId,
			      fr_tagED_I);
	    }
	    p.logOperation(s);
	    break;

          case fr_syntactic:
	    if (op.base_old == fr_eorchar && op.base_new == fr_gapchar)
	      del_count++;
	    if (op.base_old == fr_gapchar && op.base_new == fr_eorchar)
	      ins_count++;
	    break;

	  default:
            throw Notify(Notify::WARNING, "Illegal edit-operation");
	  }
	} else {
	  // Read is hidden or hypo is not confirmed
	}
      }
      // I++;
      XX++;
    }
  }


  // has the size of the fault region changed?
  // shift all reads that begin after the fault-region

  if (locomotion != 0 && c != nullptr &&
      contigpos < c->getConsensusCounts().size())
    {
      //      locomotion = (ins_count - del_count) / (int32)aScfInfo.size();

#ifndef RUN_ONLY
      if (p.isVerbose()) {
	cout << "AdjustReadOffsets (contigpos=" << contigpos
	     << ", loco=" << locomotion << ")" << endl;
      }
#endif

      // BaCh
      // cout << "Contig loco is:\n\n" << *c << endl;

#ifdef LIGHTHEADED

      // report the error and continue editing...
      try {
	//	c->adjustReadOffsets(contigpos, locomotion);
	c->adjustReadOffsets(scfInfo_contigPosFirstBase(*I)+1, locomotion);
      }
      catch(Notify n) {
	cout << "ERROR: Problems occurred while adjustingReadOffsets!" <<endl;
      }
#else
      //c->adjustReadOffsets(contigpos, locomotion);
      c->adjustReadOffsets(scfInfo_contigPosFirstBase(*I)+1, locomotion);
#endif

  }

  if (p.tagDeleteAsteriskColumns() && c != nullptr) {
    int32 oldsize = c->getContigLength();

    if (contigpos < c->getContigLength() - 1) {
#ifndef RUN_ONLY
      if (p.isVerbose()) {
	cout << "DeleteStarOnlyColumns between " << contigpos-10
	     << " - " << contigpos+1 << endl;
      }
#endif

      // BaCh
      // cout << "Contig is:\n\n" << *c << endl;

#ifdef LIGHTHEADED
      // report the error and continue editing...
      try {
	 c->deleteStarOnlyColumns(max(contigpos-10, 0), contigpos+1);
      }
      catch(Notify n) {
	n.handleError("Error while examining fault-region");
      }
#else
        c->deleteStarOnlyColumns(max(contigpos-10, 0), contigpos+1);
#endif
    }

    // BaCh
    //cout << "New contig is:\n\n" << *c << endl;
    //cout << "Done." << endl;

    sum_loco = locomotion - (oldsize - c->getContigLength());
  }

  if (p.tagConsensus() && sum_loco != 0) {
    if (sum_loco < 0) {
      c->addTagToConsensus(contigpos, contigpos, '=',
			   "ED_D", "EdIt: Consensus base deleted",
			   false);
    } else {
      c->addTagToConsensus(contigpos, contigpos, '=',
			   "ED_I", "EdIt: Consensus base inserted",
			   false);
    }
  }


  FUNCEND();
  return;
}




int32 fault_region::shufflePadInSequence(int read_id, int pos)
{
  int32 lauf = pos + 1;
  vector<char> &aSequence = fr_sequenceList[read_id].si_sequence;

  while (aSequence[lauf] == aSequence[pos]) lauf++;

  if (aSequence[lauf] == fr_gapchar) {
    deleteBase(read_id, lauf);
    insertBase(read_id, pos, fr_gapchar);
    return 1;
  }

  lauf = pos - 1;
  while (aSequence[lauf] == aSequence[pos]) lauf--;

  if (aSequence[lauf] == fr_gapchar) {
    insertBase(read_id, pos, fr_gapchar);
    deleteBase(read_id, lauf);
    return 1;
  }
  return 0;
}


void fault_region::shufflePads()
{
  int32 stern_count, other_count;

  for (int loop = 0; loop < fr_longest; loop++) {
    stern_count = 0;
    other_count = 0;

    // Count pads and bases in the column
    vector<sequenceInfo>::iterator I = fr_sequenceList.begin();
    while (I != fr_sequenceList.end()) {
      if (I->si_sequence[loop] == fr_gapchar) {
	stern_count += I->si_occurrences;
      } else {
	other_count += I->si_occurrences;
      }
      I++;
    }

    // if gaps are predominant try to shuffle the other reads
    if (stern_count > other_count) {
      vector<sequenceInfo>::iterator J = fr_sequenceList.begin();
      while (J != fr_sequenceList.end()) {
	if (J->si_sequence[loop] != fr_gapchar) {
	  shufflePadInSequence(J - fr_sequenceList.begin(), loop);
	}
	J++;
      }
    }
  }
}



bool fault_region::isBaseInColumn(int32 col, char aBase)
{
  vector<sequenceInfo>::iterator I = fr_sequenceList.begin();
  while (I != fr_sequenceList.end()) {
    if (I->si_sequence[col] == aBase) {
       return true;
    }
    I++;
  }
  return false;
}



ostream &operator<<(ostream &ostr, fault_region const &i)
{
  vector<editoperation>::const_iterator I = i.fr_edits.begin();

  ostr << endl;

  for (uint32 l1 = 0; l1 < i.fr_sequenceList.size(); l1++) {
    const vector<char> &seq = i.fr_sequenceList[l1].si_sequence;

    for (uint32 l2 = 0; l2 < seq.size(); l2++) {
      ostr << seq[l2];
    }

    ostr << "\t" << (int)i.fr_sequenceList[l1].si_occurrences;
    ostr << "  inserts: " << (int)i.fr_sequenceList[l1].si_insertCount;
    ostr << endl;

  }

  ostr << "Mismatches: " << i.fr_score   << "  Edits: " << i.fr_steps ;
  ostr << "  Score : "   << i.fr_solution_score
       << "  size: " << i.realSize() << endl;

  while (I != i.fr_edits.end()) {
    ostr << I->base_old   << " -> " << I->base_new << " Op.";
    ostr << (int)I->operation << "  ";
    ostr << (int)I->sequenceId << ":"    << (int)I->position
	 << endl;
    I++;
  }

  ostr << "------------------" << endl;

  return ostr;
}



/* --------------------------------------
        hypotheses_pool
   -------------------------------------- */

hypotheses_pool::hypotheses_pool(const fault_region startzustand,
				 vector<scfInfo> &si)
{
  hp_ready.clear();
  hp_ready.reserve(hp_poolsize);

  hp_pool.clear();
  hp_pool.reserve(hp_poolsize*2);

  hp_start = startzustand;
  hp_pool.push_back(hp_start);

  hp_scfInfo.clear();
  hp_scfInfo.swap(si);

  hp_parameter = nullptr; // new EDITParameters;
}


hypotheses_pool::hypotheses_pool(const fault_region startzustand,
				 vector<scfInfo> &si, EDITParameters *p)
{
  hp_ready.clear();
  hp_ready.reserve(hp_poolsize);

  hp_pool.clear();
  hp_pool.reserve(hp_poolsize*2);

  hp_start = startzustand;
  hp_pool.push_back(hp_start);

  hp_scfInfo.clear();
  hp_scfInfo.swap(si);

  hp_parameter = p;
}



hypotheses_pool::~hypotheses_pool()
{
  //  theBuffer.statistics();
  hp_ready.clear();
  hp_pool.clear();

  hp_start.clear();

  vector<scfInfo>::iterator I = hp_scfInfo.begin();
  while (I != hp_scfInfo.end()) {
    if (I->scfInfo_name != nullptr) {
      delete [] I->scfInfo_name;
    }
    I++;
  }

  hp_scfInfo.clear();
}



bool hypotheses_pool::isDoubleStranded(bool &forward, bool &reverse)
{
  forward = false;
  reverse = false;
  vector<scfInfo>::iterator I = hp_scfInfo.begin();

  while (I != hp_scfInfo.end()) {
    if (I->scfInfo_read->direction > 0) {
      forward = true;
    } else {
      reverse = true;
    }
    I++;
  }

  return (forward && reverse);
}


// Is doubleStranded _or_ has at a coverage of at least c

bool hypotheses_pool::isDoublestrandedOrCoverage(bool &forward,
						 bool &reverse, int32 c)
{
  int32   readCount = 0;
  forward = false;
  reverse = false;
  vector<scfInfo>::iterator I = hp_scfInfo.begin();

  while (I != hp_scfInfo.end()) {
    if (I->scfInfo_read->direction > 0) {
      forward = true;
    } else {
      reverse = true;
    }
    readCount++;
    I++;
  }

  return ((forward && reverse) || (c > 0 && readCount >= c));
}



// are different sequencing machines used in this fault region/ for this hypothesis?
// Evaluates the MACH flag in the comments of the SCF files. This is not always set
// correctly!
bool hypotheses_pool::usesDifferentChemistry()
{
  t_machinetype t = MT_undefined;
  vector<scfInfo>::iterator I = hp_scfInfo.begin();
  SCF_look *theTrace;

  if (hp_scfInfo.size() < 2) {
    return false;
  }


  while (I != hp_scfInfo.end()) {

    theTrace = ScfBuffer::bufferRead(I->scfInfo_read->read,
				     I->scfInfo_read->direction);

    if (t != MT_undefined) {
      if (t != theTrace->getMACHType()) {
	ScfBuffer::bufferDelete(theTrace);
	return true;
      }
    } else {
      t = theTrace->getMACHType();
      ScfBuffer::bufferDelete(theTrace);
    }
    I++;
  }

  return false;
}



// Test for sufficient information for the evaluation of a hypothesis
// - double stranded is sufficient
// - coverage of minCoverage' is sufficient (minCoverage > 0)
// - different machines used is sufficienty if useChemistry = true

bool hypotheses_pool::hasSufficientInformation(bool &forward, bool &reverse,
					       bool useChemistry, int32 minCoverage)
{
  if (isDoublestrandedOrCoverage(forward, reverse, minCoverage)) {
    return true;
  } else {
    if (useChemistry == true) {
      return usesDifferentChemistry();
    } else {
      return false;
    }
  }
}




void hypotheses_pool::pushStep(fault_region &fr)
{

  if (false == fr.validSolution()){
    // We don't want this solution...
    return;
  }

  if (fr.readScore() > 0) {
    vector<fault_region>::iterator iter = hp_pool.begin();
    while (iter != hp_pool.end() && iter->compare(fr) != 0) {
      ++iter;
    }

    if (iter == hp_pool.end()) {
      hp_pool.push_back(fr);
    }

  } else {
    vector<fault_region>::iterator iter = hp_ready.begin();

    while (iter != hp_ready.end() && iter->compare(fr) != 0) {
      ++iter;
    }

    if (iter == hp_ready.end()) {
      fr.calcSolutionScore();
      hp_ready.push_back(fr);
    }
  }
}


void hypotheses_pool::iterate(fault_region &fr)
{
  int32 firstFaultCol = 1;
  fault_region nextStep;


  while (firstFaultCol < fr.size() &&
	 0 == fr.pairwiseMismatchColumn(firstFaultCol)) {

    firstFaultCol++;
  }


  if (firstFaultCol == fr.size()) {
#ifdef HYPOTHESES_VERBOSE
    cout << "Nothing to correct!" << endl;
#endif
    return;
  }

  for (int32 read_loop=0;
       read_loop < fr.readSequenceCount();
       read_loop++)
    {

      if (fr.readBase(read_loop, firstFaultCol) != '_') {
	// e.g. _ -> . makes some problems...

	for (int8 loop = 0; loop < 5; loop++) {

	  if (fr.readBase(read_loop, firstFaultCol) != "*ACGT"[loop] &&
	      fr.isBaseInColumn(firstFaultCol, "*ACGT"[loop])) {

	    nextStep = fr;
	    if (0==nextStep.alterBase(read_loop,firstFaultCol,"*ACGT"[loop])) {
	      pushStep(nextStep);
	    }

	    nextStep = fr;
	    if(0==nextStep.insertBase(read_loop,firstFaultCol,"*ACGT"[loop])){
	      pushStep(nextStep);
	    }

	  }
	}
	nextStep = fr;

	if (0 == nextStep.deleteBase(read_loop, firstFaultCol))
	  pushStep(nextStep);
      }
    }

  return;
}



bool hp_compare(const fault_region &a, const fault_region &b)
{
  return (a.readScore() < b.readScore());
}



bool hp_compare_solution(const fault_region &a, const fault_region &b)
{
  return (a.readSolutionScore() < b.readSolutionScore());
}



void hypotheses_pool::iteratePool()
{
  vector<fault_region> hp_old;

  hp_old.clear();
  hp_old.reserve(hp_poolsize*2);
  hp_pool.swap(hp_old);

  vector<fault_region>::iterator iter = hp_old.begin();

  while (iter != hp_old.end()) {
    iterate(*iter);

    if (hp_pool.size() == 0) {
      hp_old.clear();
      return;
    }

    sort(hp_pool.begin(), hp_pool.end(), hp_compare);

    if (hp_pool.size() > hp_usesize) {
      hp_pool.erase(hp_pool.begin() + hp_usesize, hp_pool.end());
    }

    ++iter;
  }

  hp_old.clear();

}


void hypotheses_pool::iteratePool(int32 anz)
{
  if (anz < 1)  return;
  if (anz > 10) {
    cerr << "Endless iterations skipped...\n";
    cerr << "number of iterations set to 3\n";
    anz = 3;
  }

  for (int loop = 0; loop<anz; loop++) {
    iteratePool();
  }
}



void hypotheses_pool::iterateFirstSolutionPlus(int32 anz)
{
  int8 count = 0;

  if (anz > 5) {
    cerr << "Endless iterations skipped...\n";
    cerr << "Iterate until first solution.\n";
    anz = 0;
  }

  while (hp_ready.size() == 0 && count < 10) {
    iteratePool();
    count++;
  }
  iteratePool(anz);

}



ostream &operator<<(ostream &ostr, hypotheses_pool const &i)
{
  vector<fault_region>::const_iterator I = i.hp_pool.begin();

  ostr << "\nHypotheses-Pool" ;

  while (I != i.hp_pool.end()) {
    ostr << *I;
    ++I;
  }

  if (i.hp_parameter != nullptr && i.hp_parameter->isVerbose(4)) {
    vector<fault_region>::const_iterator iter = i.hp_ready.begin();

    for (uint32 l = 0; l < i.hp_scfInfo.size(); l++) {
      const scfInfo &aScf = i.hp_scfInfo[l];
      if (aScf.scfInfo_read == nullptr) {
	ostr << "  " << aScf.scfInfo_name << "\t";
      } else {
	ostr << "  " << aScf.scfInfo_read->read.getName() << "\t";
      }

      for (uint32 k = 0; k < aScf.scfInfo_position.size(); k++) {
	ostr << aScf.scfInfo_position[k] << "\t";
      }

      ostr << "\tc: " << aScf.scfInfo_confirmed
	   << "\th: " << aScf.scfInfo_hidden << endl;
    }
    ostr << "====================================" << endl;
  }
  return ostr;
}


void hypotheses_pool::showResult()
{
  vector<fault_region>::iterator iter = hp_ready.begin();

  cout << "Show result:\n";

  for (uint32 l = 0; l < hp_scfInfo.size(); l++) {
    const scfInfo &aScf = hp_scfInfo[l];
    if (aScf.scfInfo_read == nullptr) {
      cout << "  " << aScf.scfInfo_name << "\t";
    } else {
      cout << "  " << aScf.scfInfo_read->read.getName() << "\t";
    }

    for (uint32 k = 0; k < aScf.scfInfo_position.size(); k++) {
      cout << aScf.scfInfo_position[k] << "\t";
    }

    cout << "\tc: " << aScf.scfInfo_confirmed
	 << "\th: " << aScf.scfInfo_hidden << endl;
  }


  sort(hp_ready.begin(), hp_ready.end(), hp_compare_solution);

  while (iter != hp_ready.end()) {
    //if ((*iter).possibleSolution()) {
    cout << "Show Result" << endl;
    cout << *iter;
    cout << "Solution Score: " << iter->readSolutionScore() << endl;
    //   }
    ++iter;
  }

}



int32 hypotheses_pool::minimalEdits() {
  vector<fault_region>::iterator iter = hp_ready.begin();
  int32 minEds = 99999;

  while (iter != hp_ready.end()) {
     int32 i = iter->editCount();
     if (i < minEds) minEds = i;
     ++iter;
  }
  return minEds;
}




bool hypotheses_pool::createCfhInfoList(EDITParameters &p)
{
  vector<fault_region>::iterator iter = hp_ready.begin();
  int32 n = 1;             // count solutions checked
  int32 max_score = 0;     // max acceptable score for a solution

  hp_parameter = &p;

  if (iter == hp_ready.end()) {

    if (p.isVerbose()) {
      cout << "No solutions found!" << endl;
    }

    return false;
  }

  // Sort the solutions
  sort(hp_ready.begin(), hp_ready.end(), hp_compare_solution);

  // We generate the solutions - best first -
  max_score = ( (100 + p.getMaxScoreOverBest())
		 * iter->readSolutionScore())/100;

  while (iter != hp_ready.end() &&
        (iter->readSolutionScore() <= max_score ||
	 n <= p.getMinSolutionsChecked())) {

    fault_region editStart = hp_start;
    editStart.generateAllHypotheses(*iter, hp_scfInfo, hp_parameter);

    iter++;
    n++;
  }

  // Delete the remaining solutions with unsufficient score
  while(iter != hp_ready.end()) {
    hp_ready.erase(iter);
  }

  return true;
}



bool fault_region::generateAllHypotheses(fault_region &target,
					 vector<scfInfo> &scf,
					 EDITParameters *p)
{
  vector<editoperation>::iterator edOpLoop;
  vector<editoperation> &theEdits = target.readEditOperations();
  vector<sequenceInfo>::const_iterator I = fr_sequenceList.begin();

  int32 loop;
  int32 problem_offset = 0;   // necessary locomotions of the edit-posit
  int32 question = fhc_UNDEFINED;
  char  newBase  = GAP_BASE;
  char  oldBase  = GAP_BASE;
  bool  found;

  // check each sequence if it was editied at the given position


  edOpLoop = theEdits.begin();

  while (edOpLoop != theEdits.end()) {

    I = fr_sequenceList.begin();
    loop = 0;
    found = false;
    edOpLoop->readEdits.clear();
    problem_offset = 0;

    while (I!= fr_sequenceList.end()) {

      if (edOpLoop->sequenceId == loop) {

	oldBase = edOpLoop->base_old;
	newBase = edOpLoop->base_new;

	switch (edOpLoop->operation) {
	case fr_syntactic: {
	  question = fhc_SYNTACTIC;

	  found = true;
	  break;
	}
	case fr_delete: {
	  int32 x = edOpLoop->position-1;
	  int32 y = edOpLoop->position+1;

	  // TODO: (BaCh) fsck, this lets valgrind puke
	  //  first repair attempt led to errors in the edited contig
	  //
	  //  I'll have to live with it

	  while (I->si_sequence[x] == GAP_BASE) { x--; }
	  while (I->si_sequence[y] == GAP_BASE) { y++; }

	  if (I->si_sequence[edOpLoop->position] == I->si_sequence[x] ||
	      I->si_sequence[edOpLoop->position] == I->si_sequence[y]) {
	   question = fhc_OVERCALL;
	  } else {
	    question = fhc_ADDITIONAL;
	  }

	  found = true;
	  break;
	}
	case fr_alter: {
	  question = fhc_WRONG;
	  found = true;
	  break;
	}
	case fr_insert: {

	  if (newBase == I->si_sequence[edOpLoop->position-1] ||
	      newBase == I->si_sequence[edOpLoop->position]) {
	    question = fhc_UNDERCALL;
	    if (newBase == I->si_sequence[edOpLoop->position-1]) {
	      problem_offset = -1;
	    }
	  } else {
	    question = fhc_MISSING;
	  }
	  found = true;
	  break;
	}
	}

	vector<scfInfo>::iterator theScfRead = scf.begin();
	afh_Info *aProblem;

	while (theScfRead != scf.end() && found) {

	  if (theScfRead->scfInfo_sequenceId == loop) {
	    try {

	      //cout << "\n Edit position " << edOpLoop->position
	      //     << "  offset: " << problem_offset << endl;
	      //cout << "OldBase " << oldBase << " NewBase "
              //     << newBase << endl;
	      //cout << theScfRead->scfInfo_read->read.getName() << endl;

	      aProblem=createSingleHypotheses(*theScfRead, *I, question,
					      newBase, oldBase,
					      edOpLoop->position
					      +problem_offset);
	      if (aProblem != nullptr) {
		target.appendAfh(aProblem);
		aProblem->scfInfoId = theScfRead - scf.begin();
		edOpLoop->readEdits.push_back(aProblem);  // FIXME V1
	      } else {
		throw Notify(Notify::FATAL,
			     "Error creating hypotheses!");
	      }

	    }
	    catch(Notify n) {
	      cout << "Error evaluating hypotheses!" << endl;
	      throw Notify(Notify::WARNING, "Error evaluating hypotheses!");
	    }
	  }
	  theScfRead++;
	}

      }

      loop++;
      I++;
    }
    performEditOperationOnRegion(*edOpLoop);
    edOpLoop++;
  }
  return true;
}


afh_Info* fault_region::createSingleHypotheses(scfInfo &theScfRead,
					  const sequenceInfo &s,
					  int32 question, char new_base,
					  char old_base, int32 pos)
{
  const Contig::contigread_t *myRead = theScfRead.scfInfo_read;
  int32 ic = s.si_insertCount;
  char  oldDBBase;
  afh_Info *theProblem;

  // FIXME ???
  // Da Änderungen auf den Reads nachgezogen werden, sollten die
  // Positionen übereinstimmen ?!?
  // int32 dbPos  = theScfRead.scfInfo_dbpos + pos - ic;
  int32 dbPos  = theScfRead.scfInfo_dbpos + pos;
  int32 scfPos ;


  //cout << "start: " << theScfRead.scfInfo_dbpos  << "\tpos : " << pos
  //    << "\tic: " << ic
  //    << "\t" << theScfRead.scfInfo_read->read.getName() << endl;


  if (pos < ic) {
    scfPos = -1;
  } else {
    assert(pos-ic <= fr_longest);

    scfPos = theScfRead.scfInfo_position[pos- ic];
    //scfPos = theScfRead.scfInfo_position[pos];

    if (scfPos < 0 && question != fhc_MISSING &&
	question != fhc_UNDERCALL && question != fhc_SYNTACTIC){

      // sollte eigentlich eine Position dasein....
      // wir versuchen zu interpolieren

      if (abs(theScfRead.scfInfo_position[pos - ic + 1] -
	      theScfRead.scfInfo_position[pos - ic - 1]) == 2) {

	scfPos = (theScfRead.scfInfo_position[pos - ic + 1] +
		  theScfRead.scfInfo_position[pos - ic - 1]) / 2;
      } else {
	question = fhc_UNDEFINED;
      }
    }

  }


  oldDBBase = old_base;
  /*
  if (question == fhc_MISSING || question == fhc_UNDERCALL) {
    oldDBBase = GAP_BASE;
  } else {
    oldDBBase = old_base;
  }
  */


  const vector<base_quality_t> &q = myRead->read.getQualities();
  vector<base_quality_t>::const_iterator Q = q.begin();

  {
    int32 qual_pos;
    if (myRead->direction < 0) {
      qual_pos = myRead->read.getLenSeq() - 1 - dbPos;
      //  cout << "Adv -" << (int32)myRead->read.getLenSeq()-1-dbPos << endl;
    } else {
      qual_pos = dbPos;
      //  cout << "Adv + " << dbPos << endl;
    }
    if (qual_pos < 0) qual_pos = 0;
    advance(Q, qual_pos);
  }


  //cout << "Quality of " << myRead->read.getName() << " at " << dbPos
  //     <<  " ScfPos: " << scfPos << " is " << (int)(*Q) << endl;

  theProblem = new afh_Info(myRead, dbPos, scfPos, new_base,
			    oldDBBase, question);

  if (theProblem != nullptr && Q != q.end()) {
    theProblem->setOriginalBaseQuality(*Q);
  }

  return theProblem;
}



// *******************************************************
//    addExtendableReads:
//    Search the contig for reads in 'direction' where the
//    quality cutoff covers the position dbPos.
// *******************************************************


bool hypotheses_pool::makeDoubleStranded(Contig &aContig, int32 contigPosVon,
				 int32 contigPosBis, bool forward,
				 bool reverse, const EDITParameters &p)
{
  const vector<Contig::contigread_t> &allReads = aContig.getContigReads();
  vector<Contig::contigread_t>::const_iterator I= allReads.begin();
  scfInfo *aScfInfo;

  int32 umgebung = p.getDoublestrandUmgebung();
  int32 direction;
  bool  reads_added = false;


  if (!forward) {
    if (!reverse) direction = 0; else direction = 1;
  } else {
    direction = -1;
  }

  while (I != allReads.end()) {
    int32 relPosVon = contigPosVon - I->offset;
    int32 relPosBis = contigPosBis - I->offset;
    int32 offset = 0;
    bool  xtend  = false;
    const Read &R = I->read;


    if (I->direction == direction || direction == 0) {
      // the read is in the desired direction ?
      // >0 forward read    <0 reverse read  =0 very read

#ifdef DEBUG_MODE
                 cout << "->" <<R.getName() << " O: " << I->offset
      	   << " LC:" << R.getLeftClipoff() << " ("
      	   << R.getLSClipoff() << ") "
      	   << " RC:" << R.getRightClipoff() << " ("
      	   << R.getRSClipoff() << ") "
      	   << " L:" << R.getLenClippedSeq() << " ("
      	   << R.getLenSeq() << ") " <<  "  dbPos:" << contigPosVon
      	   << " Direction: " << I->direction <<endl;
#endif

      xtend = false;

      if (I->direction > 0) {
        // forward read in the leftClipOff part?

	if(relPosVon < 0 && relPosVon + R.getLeftClipoff() > 0) {
	  // can we extend the read?

	  if ((abs(relPosVon) + umgebung) < R.getLeftExtend() &&
	      (abs(relPosBis) + umgebung) < R.getLeftExtend()) {

	    xtend = extendRead(aContig, *I, contigPosVon, -relPosVon,
	    		       umgebung, offset);
	  }
	}

	// forward read in the rightClipoff part?


	if(relPosVon >= (int32)R.getLenClippedSeq() &&
	   relPosVon < (int32)(R.getLenSeq()-R.getLeftClipoff()) -  umgebung) {
	  // can we extend the read?

	  if (relPosBis < R.getLenClippedSeq() + R.getRightExtend()) {

	    xtend = extendRead(aContig, *I, contigPosVon,
			       relPosBis - R.getLenClippedSeq(),
			       umgebung, offset);
	  }
	}

      } else {
	if (relPosVon < 0 &&
	    relPosVon > (int32(-R.getLenSeq() + R.getRightClipoff()))) {
	  // reverse read in the left clipoff part?


	  if (-relPosVon + umgebung < R.getRightExtend()) {

	    // can we extend the read?  FIXME
	    xtend = extendRead(aContig, *I, contigPosVon,
			       -relPosVon, umgebung, offset);
	  }
	}

        if (relPosVon >= (int32)R.getLenClippedSeq() &&
	    relPosVon < (int32)(R.getLenClippedSeq()
				 + R.getLeftClipoff() + umgebung)) {
	  // reverse read in the right clipoff part?


	  if (relPosBis - R.getLenClippedSeq() + 1 < R.getLeftExtend()) {
	    // can we extend this read?
	    xtend = extendRead(aContig, *I, contigPosVon,
			       relPosVon - R.getLenClippedSeq(),
			       umgebung, offset);
	  }
	}
      }




      if (xtend) {
	// We found something to extend
	char *sequence;
 	// sequence will hold the null-terminated part of the sequence from
	// contigPosVon to ContigPosBis: [contigPonVon - ContigPosBis] + '\0'
	sequence = new char [contigPosBis - contigPosVon + 2];


#ifdef HYPOTHESES_VERBOSE
	cout << R.getName() << ": extended Fault Region: (Offset "
	     << offset << ")  Direction " << I->direction << endl;
#endif


	aScfInfo = faultRegionFromRead(aContig, *I, I - allReads.begin(),
				       contigPosVon - offset,
 				       contigPosBis - offset,
				       sequence);

	hp_start.addSequence(sequence, *aScfInfo);
	hp_scfInfo.push_back(*aScfInfo);
	hp_pool.clear();
	hp_pool.push_back(hp_start);

	delete aScfInfo;
	delete [] sequence;
	reads_added = true;

	// return reads_added;   // FIXME ... nur einen read aufdecken
      }
    }
    I++;
  }

  return reads_added;
}



int32 hypotheses_pool::locomotion(const fault_region &r) const
{
   return (r.realSize() - hp_start.realSize());
}



// checkPartialSolution
//    result  1:    all inserts/deletes confirmed; edit with all
//                  insert/delete operations - also syntactic!
//    result  0:    only non ins/del ops confirmed; edit without syntactic
//    result -1:    mixed confirmed...

int32 hypotheses_pool::checkPartialSolution(fault_region &r)
{
  // When is a partial solution valid? if it does not disturb the alignment!
  // Disturbing the alignment can be done by
  //    - insert base
  //    - delete base operations.
  // If we
  //  allow only solutions where all or none of the insert/delete
  //  operations are confirmed we should be safe.

  vector<scfInfo>::iterator I;
  editoperation op;
  vector<editoperation> edits = r.readEditOperations();

  bool all_confirmed = true;
  bool all_rejected  = true;


  for (uint32 i = 0; i < edits.size(); i++) {
    op = edits[i];

    if ((op.base_old == fr_fillchar) || (op.base_new == fr_fillchar)) {
      vector<afh_Info*>::iterator XX = op.readEdits.begin();

      while (XX != op.readEdits.end()) {
	all_confirmed = all_confirmed && (*XX)->isConfirmed();
	all_rejected  = all_rejected  && !(*XX)->isConfirmed();

	XX++;
      }
    }
  }

  if (hp_parameter->isVerbose(3)) {
    cout << "Checking partial solution: "
	 << "all confirmed " << all_confirmed << "  "
	 << "all rejected  " << all_rejected  << endl;
  }

  if (!(all_confirmed || all_rejected)) {
    return -1;
  } else {
    if (all_confirmed) return 1;
    return 0;
  }
}



void hypotheses_pool::setFirstAndLastColumn(char first, char last)
{
 hp_start.setColumnToBase(0, first);
 hp_start.setColumnToBase(hp_start.realSize()-1, last);

 if (hp_pool.size() > 0) {
   hp_pool.begin()->setColumnToBase(0, first);
   hp_pool.begin()->setColumnToBase(hp_start.realSize()-1, last);
 }
}
