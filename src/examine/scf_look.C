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
 * SCF examination methods
 *
 * 13.08.98  Now loading of scf-files using a contigread_t object to
 *           get the assignment to the scf-file possible
 * 04.11.98  SCF_buffer operating without instances
 * 17.04.2002 Changed load routine to honor full path names (BaCh)
 */


#include <math.h>
#include "examine/scf_look.H"
#include "io/scf.H"
#include "globals.H"
#include <stdlib.h>
#include <cstdio>
#include <assert.h>


int32 SCF_look::count = 0;


/* ------------------------ */
/*          scf_look        */
/* ------------------------ */


// Copy constructor
SCF_look::SCF_look(SCF_look const &other)
{

#ifdef VERBOSE
  count++;
  cout << "SCF_look count CopyOperator " << count;
#endif

  zeroVars();
  *this = other;
  return;
}


// Copy operator
const SCF_look& SCF_look::operator=(const SCF_look &other)
{
  if (this != &other) {
    discard();

    SCF::operator=(other);

    reversed   = other.reversed;
    db_samples = other.db_samples;
    theRead    = other.theRead;

    if (db_samples!=0) {
      db_peak_index = new uint32[db_samples];
      db_bases = new uint8[db_samples];

      memcpy(db_peak_index, other.db_peak_index, db_samples*sizeof(uint32));
      memcpy(db_bases, other.db_bases, db_samples*sizeof(uint8));
    }

    if (other.fName != nullptr) {
      fName = new char[strlen(other.fName)+1];
      strcpy(fName, other.fName);
    }

  }
  return *this;
}




SCF_look::SCF_look()
{
#ifdef VERBOSE
  count++;
  cout << "SCF_look count ^ " << count << endl;
#endif

  zeroVars();
}



SCF_look::~SCF_look()
{
#ifdef VERBOSE
  count--;
  cout << "SCF_look count v" << count << "  " << fName << endl;
#endif

  discard();
}



void SCF_look::zeroVars()
{
  // BaCh: Added these two inits, especially db_samples caused segfaults *grummel*
  reversed=0;
  db_samples=0;

  db_peak_index = nullptr;
  db_bases   = nullptr;
  theRead = nullptr;
  fName   = nullptr;
}


void SCF_look::discard()
{
  if (db_peak_index != nullptr) delete [] db_peak_index;
  if (db_bases != nullptr)      delete [] db_bases;
  //  if (theRead != nullptr)       delete [] theRead;
  if (fName != nullptr)         delete [] fName;

  zeroVars();
}



int8 SCF_look::setFilename(const char *filename)
{
  if (filename == nullptr) return 1;

  if (fName != nullptr) {
    delete [] fName;
    fName = nullptr;
  }

  fName = new char[strlen(filename) + 1];
  strcpy(fName, filename);
  return 0;
}



int8 SCF_look::load()
{
  if (fName == nullptr) return 1;
  return load(fName);
}



void SCF_look::interpolate() {
  int32 lastBase, step;
  int32 loop = 1;

  if (db_peak_index[0] == 0) {
    while (db_peak_index[loop] == 0) { loop++; }
    for (int32 l = 0; l < loop; l++) {
      db_peak_index[l] = (db_peak_index[loop] * l) / loop;
    }
  }

  for (; loop < db_samples; loop++) {
     if (db_peak_index[loop] == 0) {

      /* find the right step-width to interpolate between the bases */
      lastBase = loop;
      while (db_peak_index[lastBase] == 0) lastBase++;


      if (lastBase >= db_samples) {
	assert((db_samples - loop + 1) != 0);

        step = ((int32) SCF_header.samples - 1 -
		(int32) db_peak_index[loop-1]) / (db_samples - loop + 1);
      } else {
	assert((lastBase - loop + 1) != 0);

	step = ((int32) db_peak_index[lastBase] -
		(int32) db_peak_index[loop-1]) / (lastBase - loop + 1);
      }

      /* do the interpolation */
      while (db_peak_index[loop] == 0 && loop < db_samples) {
	db_peak_index[loop] = db_peak_index[loop-1] + step;
        // cout<<"Interpolate "<<loop << ": "<<db_peak_index[loop] << endl;
	loop++;
      }
    }
  }

}


/*
   Load an SCF-File and a *.ass file to connect the bases over the
   db-positions (instead of using the SCF-file positions.

*/
int8 SCF_look::load(const char *filename)
{
  FUNCSTART("int8 SCF_look::load(const char *filename)");

  ifstream assignment;
  char   *fNameAss;
  char   *fNameSCF;
  int32  index;

  if (setFilename(filename)) return 1;


  fNameSCF = new char[strlen(filename) + 5];
  strcpy(fNameSCF, filename);
  strcat(fNameSCF, "SCF");

  SCF::load(fNameSCF);
  transposeAmbiguityCodes();

  fNameAss = new char[strlen(filename) + 5];
  strcpy(fNameAss, filename);
  strcat(fNameAss, ".ass");
  assignment.open(fNameAss, ios::in);
  if(!assignment) {
    MIRANOTIFY(Notify::WARNING, "Assignmentfile not found for loading: " << fName);
  }

  assignment >> db_samples;

  db_peak_index = new uint32[db_samples];
  db_bases = new uint8[db_samples];

  {
     for (int32 loop=0; loop < db_samples; loop++) {
      assignment >> db_bases[loop] >> index;
      if (index != 0) {
	db_peak_index[loop] = SCF_peak_index[index-1];
      } else {
	db_peak_index[loop] = 0;
      }
    }
  }

  interpolate();

  if (db_peak_index[0] > db_peak_index[db_samples-1]) {
    reversed = 1;
  } else {
    reversed = 0;
  }

  if (fNameSCF != nullptr) {
    delete [] fNameSCF;
    fNameSCF = nullptr;
  }
  if (fNameAss != nullptr) {
    delete [] fNameAss;
    fNameAss = nullptr;
  }

  return 0;
}



int8 SCF_look::simpleLoad(const char *filename) {


  SCF::load(filename);
  transposeAmbiguityCodes();

  db_samples = getNumBases();

  db_peak_index = new uint32[db_samples];
  db_bases = new uint8[db_samples];


  for (int32 loop=0; loop < db_samples; loop++) {
    db_bases[loop] = SCF_bases[loop];
    db_peak_index[loop] = SCF_peak_index[loop];
    /*
    cout << loop << " :  " << db_bases[loop] << " "
	 << db_peak_index[loop] << endl;
    */
  }

  reversed = 0;

  return 0;
}


// TODO: BaCh 04.10.2008
//  cleanup this sickening mixture of old char * and string
int8 SCF_look::load(const Read &aRead, int32 richtung)
{
  FUNCSTART("int8 SCF_look::load(const Read &aRead, int32 richtung)");

  // BaCh 30.07.2009
  // handle reads that have no adjustments
  if(!aRead.usesAdjustments()){
    MIRANOTIFY(Notify::WARNING, "Read has no adjustments: " << aRead.getName());
  }

  char *fNameSCF = nullptr;
  const vector<char> &sequence = aRead.getActualSequence();
  const vector<char> &csequence = aRead.getActualComplementSequence();

  theRead = &aRead;

  if (setFilename(aRead.getName().c_str())) return 1;

  try {
    fNameSCF = new char[aRead.getSCFName().size() + 1];
    strcpy(fNameSCF, aRead.getSCFName().c_str());

    string dummy;
    aRead.getSCFFullPathName(dummy);
    SCF::load(dummy.c_str());
    transposeAmbiguityCodes();
  }
  catch (Notify n) {
    cerr << "Error loading SCF-File " << fNameSCF << endl;
    if (fNameSCF != nullptr) {
      delete [] fNameSCF;
      fNameSCF = nullptr;
    }
    discard();
    zeroVars();

    n.handleError("Error loading SCF-File");
    MIRANOTIFY(Notify::WARNING, "Error loading SCF-File: " << fNameSCF);
  }


  db_samples = sequence.size();

  db_peak_index = new uint32[db_samples];
  db_bases = new uint8[db_samples];


  const vector<int32> &alignToScf = aRead.getAdjustments();

  {
    vector<int32>::const_iterator I = alignToScf.begin();

    if (richtung < 0) {
      reversed = 1;
    } else {
      reversed = 0;
    }

    for(int32 loop=0; loop < db_samples; loop++) {
      int32 dummyPos;

      if (richtung > 0) {
        db_bases[loop] = sequence[loop];
	dummyPos = loop;
      } else {
	db_bases[loop] = csequence[loop];
	dummyPos = db_samples - loop - 1;
      }


      switch (toupper(db_bases[loop])) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case '*': break;
      default: db_bases[loop] = UNDEFINED_BASE;
      }

      if (*I < 0 || I == alignToScf.end()) {
	db_peak_index[dummyPos] = 0;
      } else {
	if (*I < (int32)SCF_header.bases) {
	  db_peak_index[dummyPos] = SCF_peak_index[*I];
	} else {
	  discard();
	  zeroVars();

	  MIRANOTIFY(Notify::WARNING, "Error: inconsistent trace data!");
	}
      }

      I++;
    }
  }

  interpolate();

  delete [] fNameSCF;

  return 0;
}



// =========================================================
// Some dirty functions
// hope you know, what you are using!!
//
// These funtions are made to simulate errors for generating negative
// examples in the lerning set. Bases are inserted, hypotheses are tested
// and the bases are deleted afterwards.
//
// the size of the array is not changed during insert/delete; thus
// the last base in the array is duplicated or corrupted!!
//
// We do *not* shift the base-quality values because we:
//     - expect they are not used
//     - regard delete/insert as __temporary__ operations that should be
//       corrected by the converse operation as soon as possible!
//

int32 SCF_look::deleteBaseDB(int32 pos)
{

  if (pos < 0 || pos >= db_samples) {
    return 1;
  }

  for (int32 loop = pos; loop < db_samples-1; loop++) {
    db_peak_index[loop] = db_peak_index[loop+1];
    db_bases[loop] = db_bases[loop+1];

    /*
    SCF_prob_A[loop] = SCF_prob_A[loop+1];
    SCF_prob_C[loop] = SCF_prob_C[loop+1];
    SCF_prob_G[loop] = SCF_prob_G[loop+1];
    SCF_prob_T[loop] = SCF_prob_T[loop+1];
    */
  }

  //db_samples--;  // do not change sample size
  return 0;
}


int32 SCF_look::insertBaseDB(int32 pos, int32 scfPos, char newBase)
{

  if (pos < 0 || pos > db_samples-1) {
    return -1;
  }

  for (int32 loop = db_samples-1; loop > pos; loop--) {
    db_peak_index[loop] = db_peak_index[loop-1];
    db_bases[loop]      = db_bases[loop-1];

    /*
    SCF_prob_A[loop] = SCF_prob_A[loop-1];
    SCF_prob_C[loop] = SCF_prob_C[loop-1];
    SCF_prob_G[loop] = SCF_prob_G[loop-1];
    SCF_prob_T[loop] = SCF_prob_T[loop-1];
    */
  }

  db_bases[pos] = newBase;

  {
    // find a Position for the new peak

    int32 newPeakPos = -1;
    int32 left;
    int32 right;

    if (scfPos >= 0) {
      newPeakPos = SCF_peak_index[scfPos];
    }

    if (pos > 0) {
      left = db_peak_index[pos-1];
    } else {
      left = db_peak_index[pos];
    }

    if (pos < db_samples - 1) {
      right = db_peak_index[pos+1];
    } else {
      right = db_peak_index[pos];
    }

    if (newPeakPos < left || newPeakPos > right) {
      if (pos > 0 && pos + 1 < db_samples) {
	db_peak_index[pos] = (left +  right) / 2;
      }
    } else {
      db_peak_index[pos] = newPeakPos;
    }
  }
  return 0;
}



int32 SCF_look::alterBaseDB(int32 pos, char newBase)
{
  if (pos < 1 || pos >= db_samples) {
    return -1;
  }

  db_bases[pos] = newBase;
  if (db_peak_index[pos] < 1) {
    db_peak_index[pos] = (db_peak_index[pos+1] + db_peak_index[pos-1]) / 2;
  }

  return 0;
}


#if 1

// BaCh: let's use the string class here
t_machinetype SCF_look::getMACHType() const
{
  string scfcomments=getComments();

  // use find ... MACH might be present more than once
  // lets hope the first one is the right one
  string::size_type startpos = scfcomments.find("MACH=");

  if (startpos != string::npos) {
    string::size_type endpos = scfcomments.find('\n',startpos);
    string::size_type foundpos;

    foundpos = scfcomments.find("ABI 3730",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_ABI_3700;;
    }
    foundpos = scfcomments.find("ABI 3700",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_ABI_3700;;
    }
    foundpos = scfcomments.find("ABI3730",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_ABI_3700;;
    }
    foundpos = scfcomments.find("ABI3700",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_ABI_3700;;
    }
    foundpos = scfcomments.find("ABI 373A",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_ABI_3700;;
    }
    foundpos = scfcomments.find("ABI373A",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_ABI_3700;;
    }
    foundpos = scfcomments.find("LI-COR",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_LICOR;
    }
    foundpos = scfcomments.find("MegaBACE",startpos);
    if (foundpos != string::npos && foundpos < endpos) {
      return MT_MegaBACE;
    }
  }

  return MT_undefined;
}

#else

// BaCh: valgrind pukes
t_machinetype SCF_look::getMACHType() const
{
  const SCF_Comments *comment;
  char token[] = "MACH=";
  char tokenvalue[30];

  int32   loop=0;
  int32   i=0;
  int32   j=0;
  bool    found = false;

  comment = getComments();

  while (comment[loop] != 0 && found == false) {
    if (comment[loop] == '\n' || loop == 0) {

      // Test for token
      i=0;
      if (loop > 0) loop++;  // skip the linebreak if not start of comment...

      while (comment[loop] != 0 && token[i] == comment[loop] && token[i] != 0) {
	i++; loop++;
      }

      // token found read tag....
      if (token[i] == 0) {
	j=0;
	while (j<20 && comment[loop] != 0 && comment[loop] != '\n') {
	  tokenvalue[j] = comment[loop];
	  j++; loop++;
	}
	tokenvalue[j] = 0;
	found =  j > 1;
      }
    }
    loop++;

  }

  if (strstr(tokenvalue, "ABI 373A") != nullptr) {
    //    cout << "ABI3700" << endl;
    // ABI 3700 !!!!
    return MT_ABI_3700;
  }
  if (strstr(tokenvalue, "LI-COR") != nullptr) {
    // LI-COR
    // cout << "LICOR" << endl;
    return MT_LICOR;
  }
  if (strstr(tokenvalue, "MegaBACE") != nullptr) {
    // MegaBACE
    // cout << "MagaBace" << endl;
    return MT_MegaBACE;
  }

  return MT_undefined;
}
#endif


// End of dirty functions
// ============================================================



char* SCF_look::getFileName() const
{
  return fName;
}


uint16* SCF_look::getTraceData(const char dbbase) const
{
 char aBase = toupper(dbbase);

 if (reversed != 0) aBase = invertBase(aBase);

 switch (aBase) {
    case 'A': return SCF_samples_A;
    case 'C': return SCF_samples_C;
    case 'G': return SCF_samples_G;
    case 'T': return SCF_samples_T;
    default: return nullptr;
  }
}



uint16* SCF_look::getTraceData(const char dbbase, int32 isReversed) const
{
 char aBase = dbbase;

 if (isReversed != 0) aBase = invertBase(aBase);
 switch (aBase) {
    case 'A': return SCF_samples_A;
    case 'C': return SCF_samples_C;
    case 'G': return SCF_samples_G;
    case 'T': return SCF_samples_T;
    default: return nullptr;
  }
}



void SCF_look::showDBArea(int32 pos, int32 circle, ostream &output ) {
  int32 start;
  int32 end;

  start = max(pos - circle, 0);
  end   = min(pos + circle, (int32)db_samples - 1);

  while (start <= end) {
    if (start == pos) output << "->";
    output << db_bases[start] << " \t";
    if (start == pos) output << "<- ";
    start++;
  }

  output << endl;
  start = max(pos - circle, 0);
  while (start <= end) {
    if (start == pos) output << "->";
    output << db_peak_index[start] << " \t";
    if (start == pos) output << "<- ";
    start++;
  }

  /*
  output << endl;
  start = max(pos - circle, 0);
  while (start <= end) {
    if (start == pos) output << "->";
    output << (int) getDBQuality(start) << "  \t";
    if (start == pos) output << "<- ";
    start++;
  }
  */

  output << endl;
  start = max(pos - circle, 0);
  while (start <= end) {
    if (start == pos) output << "->";
    output << start << "  \t";
    if (start == pos) output << "<- ";
    start++;
  }

}


// FIXME: scheint irgendwie nicht richtig zu funktionieren?!
/*
base_quality_t SCF_look::getDBQuality(int32 dbpos) const
{
    if (theRead != nullptr) {

    if (dbpos > 0 && dbpos < db_samples) {
      const vector<base_quality_t> &q = theRead->getQualities();
      vector<base_quality_t>::const_iterator I = q.begin();

      if (isReversed()) {
        advance(I, db_samples - dbpos - 1);
      } else {
	advance(I, dbpos);
      }

      return *I;
    }
  }

  return -1;

}
*/


void SCF_look::showArea(int32 pos, int32 circle, ostream &output) {
  int32 start;
  int32 end;

  start = max(pos - circle, 0);
  end   = min(pos + circle, (int)SCF_header.bases - 1);

  if (reversed == 0) {
    while (start <= end) {
      if (start == pos) output << "=>";
      output << SCF_bases[start] << " ";
      if (start == pos) output << "<= ";
      start++;
    }
  } else {
    while (start <= end) {
      if (end == pos) output << "=>";
      output << invertBase(SCF_bases[end]) << " ";
      if (end == pos) output << "<= ";
      end--;
    }
  }

}



SCF_distance* SCF_look::calcDBDistance(int32 dbPos)
{
  int32 posStart, posEnd;
  int32 count_steps, count_bases, radius = 5;
  int32 left, right, rep;
  int32 min = 999999;
  int32 max = 0;

  int32 repBases = 0;   // FIXME

  if (db_bases[dbPos] == '*') {
    posStart = dbPos;
    while (db_bases[posStart-1] == '*' && posStart > 0) posStart--;
    posEnd = dbPos;
    while (db_bases[posEnd+1] == '*' < db_samples-2) posEnd++;
    repBases = 0;
  } else {
    posStart = dbPos;
    repBases = findRepeatSize(posStart, posEnd);
  }

  // include distance to the neighbors

  if (posStart > 1) {
    do {
      posStart--;
      repBases++;
    } while (posStart > 1 && db_bases[posStart] == '*');
  }

  if (posEnd < db_samples-1) {
    do {
      posEnd++;
      repBases++;
    } while (posEnd < db_samples-1 && db_bases[posEnd] == '*');
  }

  if (posEnd > db_samples-1) {
    posEnd = db_samples-1;
    repBases = repBases - (posEnd - db_samples + 1);
  }

  count_steps = 0;
  count_bases = 0;
  while (count_bases < radius && count_steps < posStart) {
    if (GAP_BASE != db_bases[posStart - count_steps])
      { count_bases++; }
    count_steps++;
  }

  left = db_peak_index[posStart] - db_peak_index[posStart-count_steps];

  //cout << "left : " << left << " posStart " <<posStart << "  pi:"
  //   << db_peak_index[posStart] << "-" << db_peak_index[posStart-count_steps]
  //     << endl;

  count_steps = 0;
  while (count_bases < (2*radius) && (count_steps + posEnd) < db_samples-1) {
    if (GAP_BASE != db_bases[posEnd + count_steps])
      { count_bases++; }
    count_steps++;
  }
  right = db_peak_index[posEnd + count_steps] - db_peak_index[posEnd];

  //  cout << "right : " << right << " posEnd " <<posStart << "  pi:"
  //     << db_peak_index[posEnd] << "-" << db_peak_index[posEnd+count_steps]
  //     << endl;
  {
    // find minimum/maximum distance between two peaks
    // => skip gaps in the db_bases!!

    int32 loop = posStart;
    int32 step = 1;

    while (loop < posEnd) {
      while (db_bases[loop+step] == '*') step++;

      // BaCh:
      // gcc3 barfs (rightly) here with the abs()
      // dirty fix: dirty cast
      if (abs((int32) db_peak_index[loop+step] - (int32) db_peak_index[loop]) < min) {
	min = abs((int32) db_peak_index[loop+step] - (int32) db_peak_index[loop]);
      }
      if (abs((int32) db_peak_index[loop+step] - (int32) db_peak_index[loop]) > max) {
	max = abs((int32) db_peak_index[loop+step] - (int32) db_peak_index[loop]);
      }

      loop += step;
      step = 1;
    }
  }

  DEBUG_EDIT("Distances: min = " << min << "    max=" << max << endl);
  DEBUG_EDIT("peak  posStart = " << posStart << "  end=" << posEnd << endl);

  // BaCh:
  // gcc3 barfs (rightly) here with the abs()
  // dirty fix: dirty cast
  rep = abs((int32) db_peak_index[posStart] - (int32) db_peak_index[posEnd]);

  {
    SCF_distance *pd = new SCF_distance(dbPos);

    pd->setMinMaxProblemDistance(min, max);
    pd->setMeanVicDistance(left+right, count_bases);
    pd->setProblemDistance(rep, repBases-1);

#ifdef PARAMETER_VERBOSE
    cout << *pd;
#endif
    return pd;
  }
}


float SCF_look::singlePeakQuality(int32 dbpos, float sigma, char b)
{
  int32  minpos, maxpos;
  char   base;
  double x;

  if (b == '-') {
    base = db_bases[dbpos];
  } else {
    base = toupper(b);
  }


  if (isBase(base) == false) {
    // cout << "Not a base! " << endl;
    // not base: think about it!
    return -0.01;
  }

  if (dbpos < 1 || dbpos >= db_samples-1) {
    // at the fringe: think about it!
    // cout << "out of range" << endl;
    return -0.02;
  }

  wide_range(dbpos, dbpos, minpos, maxpos);

#ifndef RUNONLY
  assert(maxpos < (int32)SCF_header.samples);
#endif


  if (isRealBase(base)) {
    x = simpleNormedCorrelation(getTraceData(base), minpos, maxpos, sigma);
  } else {
    double dummy = 0;

    x = 0;
    for (int16 loop = 0; loop<4; loop++) {
      dummy += simpleNormedCorrelation(getTraceData("ACGT"[loop]),
				      minpos, maxpos, sigma);
      //if (dummy > x) {
      //	x = dummy;
      //}
    }
    x = dummy/4;
  }

  return x;
}



// einzelpeak:
// max        wird mit den maxima der base belegt
// lage       ein qualitätsmass für die lage des peaks
// rückgabe:  hoehenkriterium   2* peak/(min1 + min2)

float SCF_look::einzelpeak(int32 dbPos, char dbBase, int32 &max,
			   float &pos_rating)
{
  int32 min_left = 0, min_right = 0, max_pos;
  int32 von, bis;
  int32 pl, pr;        // Peak-Position links und rechts
  int32 fromBase, toBase;
  float hoehenkriterium = 0;
  uint16 *samples;

  max = 0;
  pos_rating = 1;  // we found no peak position of interest until now ;-)


  if (dbPos > 0) {
    // examine peak between two bases (dbPos and the real base on the
    // left). Use this if you want to find a missing peek in between.
    toBase = dbPos;
    while(0 == isBase(db_bases[toBase]) && toBase < db_samples-1) toBase++;
    fromBase = dbPos-1;
    while(0 == isBase(db_bases[fromBase]) && toBase > 1) fromBase--;

    if (fromBase < 1 || toBase >= db_samples-1) {
      return -1;
    }

    DEBUG_EDIT("Base Position: " << dbPos << " Base:" << dbBase << endl);
    DEBUG_EDIT("Between : " << db_bases[fromBase-1] << db_bases[fromBase] \
	       << "<->"<< db_bases[toBase] << db_bases[toBase+1] << endl);

    von = db_peak_index[fromBase];
    bis = db_peak_index[toBase];

    pl = von;
    pr = bis;


  } else {
    // examine peak around a base. Use this if you want to examine
    // an existing peak.
    dbPos = -dbPos;

    toBase = dbPos+1;
    while(false == isBase(db_bases[toBase]) &&
	  toBase < db_samples-1)
      toBase++;

    fromBase = dbPos-1;
    while(false == isBase(db_bases[fromBase]) &&
	  toBase > 1)
      fromBase--;

    if (fromBase < 1 || toBase >= db_samples-1) {
      return -1;
    }

    von = db_peak_index[fromBase];
    bis = db_peak_index[toBase];

    pl = bestPeakBetweenBases(getTraceData(db_bases[fromBase]), von-5, von+5);
    pr = bestPeakBetweenBases(getTraceData(db_bases[toBase]), bis-5, bis+5);

    if (pl < 1) { pl = von; }
    if (pr < 1) { pr = bis; }

    range(dbPos, von, bis);

    //cout << "pl = " << pl << "\t pr = " << pr << endl;
    //cout << "Range : " << von << "-" << bis << endl;
  }

  samples = getTraceData(dbBase);

  if (samples == nullptr) {
    pos_rating = 0.5;
    hoehenkriterium = 0.5;
    max = 0;
    return hoehenkriterium;
  }

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap
  if (pl  > pr ) { pl  ^= pr  ^=  pl ^= pr;  }  // swap

  min_left  = samples[von];
  min_right = samples[bis];

  max_pos = bestPeakBetweenBases(getTraceData(dbBase), von, bis) ;

  //cout << "Min_left " << min_left << "  min_right " << min_right << endl;
  //cout << "max_pos " << max_pos << endl;

  if (0 == max_pos) {
    DEBUG_EDIT("No Peak found between the Bases\n");
    pos_rating = 0.5;
    max = samples[max_pos];
    return -1;
  } else {
    max = samples[max_pos];
  }

  //cout << "Peak at: " << max_pos << "  max = " << max << endl;

  int32 min_abstand = min(max_pos - pl, pr - max_pos);

  if (min_abstand > 0) {
    pos_rating = 2.0 * ((float)min_abstand) / (float)(pr-pl);
  } else {
    pos_rating = 0.0;
  }

  // hoehenkriterium
  if (min_left == 0 && min_right == 0) {
    if (max > 0) {
      hoehenkriterium = 1.0;
    } else {
      hoehenkriterium = 0.0;
    }
  } else {
    assert ((min_left + min_right) != 0);
    hoehenkriterium = 0.5 * (float)(min_left + min_right) / max;

    if (hoehenkriterium > 1.0) hoehenkriterium = 1.0;
  }

  DEBUG_EDIT("Pos_rating : " << pos_rating <<
             " min_left: " << min_left << " min_right " << min_right <<
	     " Hoehenkrit: " << hoehenkriterium << endl);

  return hoehenkriterium;
}







int32 SCF_look::convergenceRegions(uint16 *samples, int32 von, int32 bis,
				   int32 &areaMin,  int32 &areaMax)
{
  int32 regStart = 0;
  int32 regEnd   = 0;
  int32 regArea  = 0;
  int32 regCount = 0;

  areaMax = 0;
  areaMin = 0;


  if (samples == nullptr) return -1;
  if (von < 1) von = 1;
  if (bis >= (int32)SCF_header.samples-1)
    bis = SCF_header.samples - 2;

  for (int32 loop = von+1; loop < bis; loop++) {
    if (samples[loop-1] + samples[loop+1] < 2*samples[loop]) {
      if (loop > regEnd+1) {  // gap or start
	if (regStart != 0) {
	  ++regCount;

	  //	  cout << "Region " << regCount << ". " << regStart <<" - "
	  //	<< regEnd << "  Area: " << regArea << endl;
	}

	if (regArea > areaMax) areaMax = regArea;
        if (areaMin == 0 || regArea < areaMin) areaMin = regArea;

	regStart = loop;
	regArea = 0;
      }
      regEnd = loop;
      regArea += samples[loop];
    }
  }

  if (regStart != 0 && regStart < bis-1) {
    ++regCount;
    //cout <<"Region " << regCount << ". " << regStart <<" - "
    //	 << regEnd << "  Area: " << regArea);
  }
  if (regArea > areaMax) areaMax = regArea;
  if (areaMin == 0 || regArea < areaMin) areaMin = regArea;


  DEBUG_EDIT("\nAreas : Anz: " << regCount << "  Max: " << areaMax
	     << "  Area Min: " << areaMin << endl);
  return regCount;
}



void SCF_look::peakRelations(int32 dbPos, char newDBBase, char oldDBBase,
			     float &newRatio, float &oldRatio,
			     float &oldNewRatio)
{
  int32 cRatio;
  int32 meanPeakSize[4];
  int32 von, bis;

  char newSCFBase;
  char oldSCFBase;

  float newMax, oldMax;

  if (isRealBase(newDBBase)) {
    newSCFBase = toSCFBase(newDBBase, isReversed());
  } else {
    int32 max[4], pos[4];
    findMaxDB(dbPos, max);
    indirectSort4(max, pos);

    newSCFBase = "ACGT"[pos[0]];
    newDBBase = toDBBase(newSCFBase, isReversed());
  }

  if (isRealBase(oldDBBase)) {
    oldSCFBase = toSCFBase(oldDBBase, isReversed());
  } else {
    int32 max[4], pos[4];
    findMaxDB(dbPos, max);
    indirectSort4(max, pos);

    if (newSCFBase != "ACGT"[pos[0]]) {
      oldSCFBase = "ACGT"[pos[0]];
    } else {
      oldSCFBase = "ACGT"[pos[1]];
    }
    oldDBBase = toDBBase(oldSCFBase, isReversed());
  }

  DEBUG_EDIT("Examine Peaks for " << oldDBBase << ":" << newDBBase << endl);

  if (dbPos >= db_samples || dbPos < 0) {
    cerr << "peakRelations: Examining base that is out of range!" << endl;
    newRatio    = 0.0;
    oldRatio    = 0.0;
    oldNewRatio = 0.0;
    return;
  }

  noiseRatioDB(dbPos, qualityRadius, cRatio, meanPeakSize);
  range(dbPos, von, bis);

  if ((von+bis) / 2 >= (int32) SCF_header.samples || (von+bis) < 0) {
    cerr << "Illegal position!!! " << endl;
  }

  uint16 *samples_old = getTraceData(oldDBBase);
  uint16 *samples_new = getTraceData(newDBBase);

  if (samples_old != nullptr) {
    int32 peakPosOld = bestPeakBetweenBases(samples_old, von, bis);

    if (peakPosOld > 0) {
      oldMax = samples_old[peakPosOld];
    } else {
      oldMax = samples_old[(von+bis)/2];
    }
  } else {
    oldMax = (SCF_samples_A[(von+bis)/2] + SCF_samples_C[(von+bis)/2] +
	      SCF_samples_G[(von+bis)/2] + SCF_samples_T[(von+bis)/2])/4;

	//findMaxDBHalf(dbPos, peakMax);
	//oldMax = valueForBase(peakMax, oldSCFBase);
  }


  if (samples_new != nullptr) {
    int32 peakPosNew = bestPeakBetweenBases(samples_new, von, bis);


    if (peakPosNew > 0) {
      newMax = samples_new[peakPosNew];
    } else {
      newMax = samples_new[(von+bis)/2];
    }
  } else {
    newMax = (SCF_samples_A[(von+bis)/2] + SCF_samples_C[(von+bis)/2] +
	      SCF_samples_G[(von+bis)/2] + SCF_samples_T[(von+bis)/2])/4;

    //findMaxDBHalf(dbPos, peakMax);
    //newMax = valueForBase(peakMax, newSCFBase);
  }

  //  int newMax  = valueForBase(peakMax, newSCFBase);
  float newMean = (float) valueForBase(meanPeakSize, newSCFBase);
  //int oldMax  = valueForBase(peakMax, oldSCFBase);
  float oldMean = (float) valueForBase(meanPeakSize, oldSCFBase);

  newRatio = oldRatio = 1.0;
  oldNewRatio = 0.0;

  // FIXME
  // Wenn es einen peak gibt, wäre es besser den peakwert zu nehmen!
  //  if (max_old == 0) { max_old = called_base; }
  //  if (max_new == 0) { max_new = hypo_base;}

  DEBUG_EDIT(newSCFBase <<" New Mean: "<<newMean<<" Max " << newMax << endl);
  DEBUG_EDIT(oldSCFBase << " Old Mean: "<<oldMean << " Max "<<oldMax << endl);

  if (newMean > 0.01) {
    newRatio = newMax / newMean;
  }
  if (oldMean > 0.01) {
    oldRatio = oldMax / oldMean;
  }
  if (oldMax > 0.01) {
    oldNewRatio = newMax/ oldMax;
  } else {
    oldNewRatio = 0.0;
  }

  DEBUG_EDIT("MaxOld : " << oldMax << " MaxNew: " << newMax << endl);

  return;
}



int32 SCF_look::findNextMax(const uint16 *samples, const int32 samplePos,
			    const bool forward)
{
  int32 i = samplePos;

  while (i > 0 && i < (int32)SCF_header.samples-1) {
    if (samples[i-1] < samples[i] && samples[i+1] <= samples[i]) {
      // Maximum
      return i;
    }
    if (forward) { i++; } else { i--; }
  }
  return i;
}



int32 SCF_look::findNextMin(const uint16 *samples, const int32 samplePos,
			    const bool forward)
{
  int32 i = samplePos;

  while (i > 0 && i < (int32)SCF_header.samples-1) {
    if (samples[i-1] >= samples[i] && samples[i+1] > samples[i]) {
      // Maximum
      return i;
    }

    if (forward) { i++; } else { i--; }
  }
  return i;
}


int32 SCF_look::findNextSlopeChange(const uint16 *samples,
				    const int32 samplePos,
				    const bool forward)
{
  int32 i = samplePos;

  while (i > 0 && i < (int32)SCF_header.samples-1) {
    if (slopeChangeTest(samples, i))
      return i;

    if (forward) { i++; } else { i--; }
  }
  return i;
}



float SCF_look::peakValleyValues(const int32 dbPos, const char base,
				 float result[6])
{
  bool  peak, conv, rest;  // what have we found
  int32 von, bis, mitte;   // area until halfway to the next base
  int32 EP, EV;            // Extrema-Peak and Extrema-Valley position
  int32 EC, ED;            // medium Peak - convergent  and divergent point
  int32 WC, WD;            // schwache Peaks - convergent/divergent
  int32 L, R;              // Positions to the left/right
  //int32 PeakPos, ValleyPos;// Positon of the actually used Peak/Valley
  float SP, MP, WP;        // Strong, medium and weak peak-values
  float SV, MV, WV;        // Strong, medium and weak valleys.
  uint16 *samples;         // samplesdate

  /*
  uint16 samples[15];


  samples[0] = 0;
  samples[1] = 39;
  samples[2] = 14;
  samples[3] = 47;
  samples[4] = 375;
  samples[5] = 404;
  samples[6] = 464;
  samples[7] = 704;
  samples[8] = 833;
  samples[9] = 702;
  samples[10] = 570;
  samples[11] = 1401;
  samples[12] = 1321;
  samples[13] = 0;
  samples[14] = 0;

  von = 1; bis = 12;
  */

  range(dbPos, von, bis);
  mitte = (von + bis + 1) / 2;

  samples = getTraceData(base);

  if (samples == nullptr) {
#ifdef VERBOSE_PARAMETER
    cout << "No Samples" << endl;
#endif
    return 0;
  }

  EP = bestPeakBetweenBases(samples, von, bis);
  EV = bestValleyBetweenBases(samples, von, bis);
  constSlopeBetweenBases(samples, von, bis, EC, ED);
  findBestWeakPeak(samples, von, bis, WC, WD);

#ifdef PARAMETER_VERBOSE
  cout << "===============================" << endl;
  cout << "Base : " << base << endl;
  cout << "Peaks    : EP " << EP << " EV " << EV << endl;
  cout << "Shoulders: EC " << EC << " ED " << ED << endl;
  cout << "Weak     : WC " << WC << " WD " << WD << endl;
#endif

  // Calculate Peak-Values
  peak = conv = rest = false;

  if (EP > 0) {
    int32 dp = abs(EP-mitte);
    if (EC > 0) {
      int32 dc = abs(EC-mitte);
      if (dp <= dc) {
	peak = true;
      } else {
	conv = true;
      }
    } else {
      peak = true;
    }
  } else {
    if (EC > 0) {
      conv = true;
    } else {
      rest = true;
    }
  }

  if (peak == true) {
    L = findNextMin(samples, EP, true);
    R = findNextMin(samples, EP, false);

    if (samples[EP] == 0) {
      SP = 0.0;
    } else {
      SP = 1.0 - (float)(samples[L] + samples[R]) / (2.0 * samples[EP]);
    }

    MP = 1.0 - SP;
    WP = 0.0;

    //    cout << "L: " << samples[L] << " R: " << samples[R] << " EP "
    // << samples[EP] << " ==> " << SP << endl;

    assert ((bis-von) != 0);

    SP = SP * (1.0 - (float)(abs(EP - mitte))/(float)(bis-von));
    MP = MP * (1.0 - (float)(abs(EP - mitte))/(bis-von));
  } else {
    SP = 0.0;
    MP = 0.0;
    WP = 1.0;
  }

  if (conv == true) {
    L = findNextSlopeChange(samples, EC-1, false);
    R = findNextSlopeChange(samples, EC+1, true);

    {
    int32 m = max(samples[L], samples[R]);

    SP = 0.0;
    if (m == 0) {
      WP = 1.0;
    } else {
      assert(m != 0);

      WP = ((float)m - samples[EC]) / m;
    }
    MP = 1.0 - WP;




    MP = MP * (1.0 - (float)abs(EC - mitte)/(float)(bis-von));
    WP = WP * (1.0 - (float)abs(EC - mitte)/(float)(bis-von));

    //    cout << "ConstSlope: "  << L << " " << EC << " " << R
    // << " ==> " << WP << endl;
    }
  }

  if (rest == true) {
     if (WC != 0) {
       assert((bis - von) != 0);

       WP = 1.0 - ((float)abs(WC - mitte)/(float)(bis-von));
    } else {
      WP = 0.0;
    }

    SP = 0.0;
    MP = 0.0;
  }


  // Now the same for valleys....
  peak = conv = rest = false;

  if (EV > 0) {
    int32 dv = abs(EV - mitte);
    if (ED > 0) {
      int32 dd = abs(ED - mitte);

      if (dv <= dd) {
	peak = true;
      } else {
	conv = true;
      }
    } else {
      peak = true;
    }
  } else {
    if (ED > 0) {
      conv = true;
    } else {
      rest = true;
    }
  }

  // Calculate Valley
  if (peak == true) {
    L = findNextMax(samples, EV, true);
    R = findNextMax(samples, EV, false);

    assert((samples[L] + samples[R]) != 0);
    assert((bis - von) != 0);

    SV = 1.0 - (float)samples[EV] / (float)(samples[L] + samples[R]);
    MV = 1.0 - SV;
    WV = 0.0;

    SV = SV * (1.0 - (float)(abs(EV - mitte))/(float)(bis-von));
    MV = MV * (1.0 - (float)(abs(EV - mitte))/(bis-von));

    //    cout << "L: " << L << " R: " << R << " EV " << EV
    //	 << " ==> " << SV << endl;
  } else {
    SV = 0.0;
    MV = 0.0;
    WV = 1.0;
  }


  if (conv == true) {
    L = findNextSlopeChange(samples, ED-1, false);
    R = findNextSlopeChange(samples, ED+1, true);

    if (samples[ED] == 0) {
      MV = 1.0;
    } else {
      MV = 1.0 - (float)(samples[ED] - min(samples[L], samples[R]))
	/ (float)samples[ED];
    }

    WV = 1.0 - MV;
    SV = 0;

    assert((bis - von) != 0);

    MV = MV * (1.0 - (float)abs(ED - mitte)/(float)(bis-von));
    WV = WV * (1.0 - (float)abs(ED - mitte)/(float)(bis-von));

    //    cout << "L: " << L << " R: " << R << " ED " << ED
    //	 << " ==> " << WV << endl;
  }

  if (rest == true) {
    if (WD != 0) {
      assert((bis - von) != 0);

      WV = 1.0 - ((float)abs(WD - mitte)/(float)(bis-von));

      //      cout << "WV " << WV << " mitte: " << mitte << " bis-von"
      //	   << bis-von <<endl;
    } else {
      WV = 0.0;
    }

    SV = 0.0;
    MV = 0.0;
  }


  result[0] = SP;
  result[1] = MP;
  result[2] = WP;
  result[3] = SV;
  result[4] = MV;
  result[5] = WV;

  return (SP + MP);
}



int32 SCF_look::countPeaksInTrace(uint16 *samples, int32 von, int32 bis,
				  int32 &peak, int16 &qual)
{
  int32 i = von;
  int32 anzPeak = 0;
  int32 peakSum = 0;
  int32 diffSum = 0;
  int32 lastPeak = 0;

  // Überlese absteigenden Teil
  while (samples[i] < samples[i-1] && i < bis)  i++;

  //1  lastPeak = samples[i];

  while (i < bis) {
    if (samples[i-1] <  samples[i] && samples[i+1] <= samples[i]) {
      // Maximum
      peakSum += samples[i];

      if (anzPeak != 0) {
	diffSum += abs(lastPeak - samples[i]);
      }

      //      cout << "Peak at : " << i << "\t" << samples[i-1]
      //	   << " " << samples[i] << " " << samples[i+1] << endl;
      anzPeak++;

      lastPeak = samples[i];
    }
    if (samples[i-1] > samples[i] && samples[i+1] >= samples[i]) {
      // Minimum
      peakSum += samples[i];

      if (anzPeak != 0) {
	diffSum += abs(lastPeak - samples[i]);
      }

      lastPeak = samples[i];
    }
    i++;
  }

  if (peakSum > 0 && diffSum > 0) {
    qual = (100 * diffSum) / peakSum;
  } else {
    if (peakSum > 0) {
      qual = 200 - (100 * (samples[von] + samples[bis])) / peakSum ;
    } else {
      qual = 0;
    }
  }

  peak = lastPeak;

  return anzPeak;
}


int32 SCF_look::calcConvergenceRegionQuality(int32 dbPos) {
  int32 fromBase, toBase;
  int32 von, bis, peak, q;

  // Skip gaps on both sides of the BOI
  toBase = dbPos;
  while(0 == isBase(db_bases[toBase]) && toBase < db_samples-1) toBase++;
  fromBase = dbPos;
  while(0 == isBase(db_bases[fromBase]) && toBase > 1) fromBase--;

  DEBUG_EDIT("  ConvergenceRegionQuality: " << endl);
  DEBUG_EDIT("  From Base " << fromBase << " toBase " << toBase << endl);

  if (fromBase < 1 || toBase >= db_samples-1) {
    return -1;
  }

  von = db_peak_index[fromBase];
  bis = db_peak_index[toBase];

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap

  if (von > 4) {
    von = von - 4;
  } else {
    von = 1;
  }

  if (bis + 5 < (int32)SCF_header.samples) {
    bis = bis + 4;
  } else {
    bis = SCF_header.samples - 2;
  }

  {
    uint16 *trace = getTraceData(db_bases[fromBase]);

    if (trace != nullptr) {
      peak = bestConvergenceRegion(trace, von, bis);
    } else {
      // e.g. if fromBase is an undefined base
      peak = von;  // not in the middle ->low quality
    }
  }

  q = 0;
  if (von != bis) {
    q = (200 * min(peak - von, bis - peak))/(bis-von);   // q in [0...100]
  }

  DEBUG_EDIT("  Quality [0..100] : " << q << endl);

  return q;
}


// Count the peaks in the trace for `base` between the base
// [basePos-1, basePos]

int32 SCF_look::countPeaksBetweenBases(char base, int32 basePos)
{
  int32  last_max_pos;
  int16  qual;         // peak quality measure
  uint16 *samples;

  samples = getTraceData(base);  // get the trace for this base
  if (samples == nullptr) {
    return -1;
  } else {
    return countPeaksInTrace(
	     samples,
	     db_peak_index[basePos-1] + 1,
	     db_peak_index[basePos] - 1,
	     last_max_pos, qual);
  }
}





int32 SCF_look::countConvergenceRegions(int32 dbBasePos, char base,
					int32 &minArea, int32 &maxArea)
{
  int32 dbBaseBis, dbBaseVon = dbBasePos;
  int32 minPos = 1;
  int32 maxPos = SCF_header.samples-1;
  int32 convRegs = 0;

  assert(isBase(base));

  if (isRealBase(base)) {
    findRepeatSize(dbBaseVon, dbBaseBis);
    wide_range(dbBaseVon, dbBaseBis, minPos, maxPos);
    convRegs = convergenceRegions(getTraceData(base), minPos, maxPos,
			minArea, maxArea);

    DEBUG_EDIT("Convergence Regions: " << convRegs << endl);
    return (abs(dbBaseBis - dbBaseVon) - convRegs);
  } else {
   int32 max[4], pos[4];
   int32 regCount1, regCount2;
   int32 areaMin1 = 0, areaMin2 = 0, areaMax1 = 0, areaMax2 = 0;
   char  base1, base2;

   findMaxDB(dbBasePos, max);
   indirectSort4(max, pos);

   base1 = "ACGT"[pos[0]];
   base2 = "ACGT"[pos[1]];

   wide_range(dbBasePos, dbBasePos, minPos, maxPos);

   regCount1 = convergenceRegions(getTraceData(base1),
 		            minPos, maxPos, areaMin1, areaMax1);

   regCount2 = convergenceRegions(getTraceData(base2),
			    minPos, maxPos, areaMin2, areaMax2);

   if (regCount1 == regCount2) {
     maxArea = (areaMax1 + areaMax2)/2;
     minArea = (areaMin1 + areaMin2)/2;
     return (regCount1 - 1);
   }
   if (regCount1 > regCount2) {
     minArea = areaMin1;
     maxArea = areaMax1;
     return (regCount1 - 1);
   } else {
     minArea = areaMin2;
     maxArea = areaMax2;
     return (regCount2 - 1);
   }
  }
}


int32 SCF_look::countRepeatPeaks(const int32 basePos, const char base,
				 int16 &anzRepeat, int16 &anzPeaks,
				 int16 &qual)
{
  FUNCSTART("int32 SCF_look::countRepeatPeaks(const int32 basePos, const char base, int16 &anzRepeat, int16 &anzPeaks, int16 &qual)");

  int32 minPos = 1;
  int32 maxPos = SCF_header.samples-1;
  int32 baseVon = basePos;
  int32 baseBis = basePos;
  int32 peakPos;
  uint16 *samples;

  // Ermittelt länge des repeats und setzt basePos auf die erste base
  // und baseBis auf die letzte Base des repeats. Ermittle dann die
  // zugehörigen sample positionen

  if (db_bases[basePos] != base) {
    anzRepeat = 0;
  } else {
    anzRepeat = findRepeatSize(baseVon, baseBis);
  }

  wide_range(baseVon, baseBis, minPos, maxPos);

  if (isBase(base)) {
    samples = getTraceData(base);
  } else {
    samples = getTraceData(db_bases[baseVon]);
  }

  if (samples == nullptr) {
    anzPeaks  = 0;
    qual      = 0;

    MIRANOTIFY(Notify::WARNING, "SCF_look::countRepeatPeaks: NoTrace");
  } else {
    anzPeaks = countPeaksInTrace(samples, minPos, maxPos, peakPos, qual);
    return anzRepeat - anzPeaks;
  }
}



int32 SCF_look::countDBPeaks(int32 basePos, int16 &anzRepeat,
			     int16 &anzPeaks, int16 &qual) {

  int32 minPos = 1;
  int32 maxPos = SCF_header.samples-1;
  uint16 * samples;

  anzRepeat = 1;
  anzPeaks  = 0;

  //  if (basePos < 0) basePos = problem->dbPos();

  {
    // Ermittelt länge des repeats und setzt basePos auf die erste base
    // und baseBis auf die letzte Base des repeats.
    int32 baseBis = 0;

    anzRepeat = findRepeatSize(basePos, baseBis);
    wide_range(basePos, baseBis, minPos, maxPos);
  }

  if (isRealBase(db_bases[basePos])) {
    int32 peakPos;

    samples = getTraceData(db_bases[basePos]);

    anzPeaks = countPeaksInTrace(samples, minPos, maxPos, peakPos, qual);

    DEBUG_EDIT("Extrema: " << anzPeaks << " Basen: " << anzRepeat<< endl);
  }

  return anzRepeat - anzPeaks;
}





// Returns the length of the repeat and sets the parameter pos_von
// to the first base and pos_bis to the last base of the repeat
int32 SCF_look::findRepeatSize(int32 &dbPosVon, int32 &dbPosBis)
{
  int32 base_count = 1;
  char  repeatBase = db_bases[dbPosVon];

  dbPosBis = dbPosVon;

  while ((repeatBase == db_bases[dbPosVon - 1] ||
	  db_bases[dbPosVon - 1] == GAP_BASE) && dbPosVon > 1) {

    dbPosVon--;
    if (db_bases[dbPosVon] != GAP_BASE) {
      base_count++;
    }
  }

  while ((repeatBase == db_bases[dbPosBis + 1] ||
	 db_bases[dbPosBis + 1] == GAP_BASE) && dbPosBis < db_samples-1) {
    dbPosBis ++;

    if (db_bases[dbPosBis] != GAP_BASE) {
      base_count++;
    }
  }

  return base_count;
}



//  Ermittelt für aufeinanderfolgende Peaks der gleichen Base das
//  prozentuale Verhältnisder mittleren Zwischentaltiefen
//  zur mittleren Peakhöhe in einer Umgebung um die angegebene Stelle

int16 SCF_look::peakQuality(int32 pos, int32 radius) {
  int32  minPos, maxPos;
  uint16 minimum, max1, max2;
  int32  von, bis;
  uint32 distance_sum = 0, distance_max = 0;
  uint16 *trace;
  char   base;

  minPos = max(2, pos - radius);
  maxPos = min(db_samples - 3, pos + radius);

  for (int32 loop = minPos + 1; loop <= maxPos; loop++){
    if (db_bases[loop] == db_bases[loop-1]) {
      base = db_bases[loop];

      if (isRealBase(base)) {
        trace = getTraceData(base);
	minimum=findMinBase(trace, db_peak_index[loop-1], db_peak_index[loop]);

	range(loop - 1, von, bis);
	max1 = findMaxBase(trace, von, bis) / 2;
	range(loop, von, bis);
	max2 = findMaxBase(trace, von, bis) / 2;

	distance_sum = distance_sum + max1 + max2 - minimum;
	distance_max = distance_max + max1 + max2;

	//cout << "Min " << min << "  Max " << max1 << ";" << max2;
	//cout << "  Sum " << distance_sum << "  Anz " << distance_max;
	//cout << " Base: " << db_bases[loop] << endl;
      }
    }
  }
  if (distance_max > 0) {
    return distance_sum * 100 / distance_max;
  } else {
    return 0;
  }
}



// Calculate the average relative-to-maximum noise of all three not-maximal traces in
// the given region (db positions)

float SCF_look::noiseRatioAll(uint32 min_pos, uint32 max_pos)
{
  int32 first_sample, last_sample;
  int32 max_array[4];
  float sum_noise = 0;

  if (min_pos > max_pos) {
    uint32 dummy = min_pos;
    min_pos = max_pos; max_pos = dummy;
  }

  for (uint32 loop = min_pos; loop <= max_pos; loop++) {
    int32 max_intensity = 0;
    int32 sum_intensity = 0;

    half_range(loop, first_sample, last_sample);
    max_intensity = findMaxDB(max_array, first_sample, last_sample);
    sum_intensity = -max_intensity; // subtract peak from the noise

    for (int16 j = 0; j < 4; j++) {
      sum_intensity += max_array[j];  // sum up the intensities
    }

    if (max_intensity > 0) {
      sum_noise += (float)sum_intensity / (float)max_intensity;
    } else {
      sum_noise += 1;
    }

  }

  return (sum_noise / float(max_pos - min_pos + 1));
}


// Calculate the mean of lowest intensity of all for traces for the given region
// (as db-positions)


float SCF_look::meanBaseline(uint32 min_pos, uint32 max_pos)
{
  int32 first_sample, last_sample;
  int32 min_array[4];
  int32 sum_intensity = 0;
  int32 min_intensity;


  if (min_pos > max_pos) {
    uint32 dummy = min_pos;
    min_pos = max_pos; max_pos = dummy;
  }

  // Step through all bases....
  for (uint32 loop = min_pos; loop <= max_pos; loop++) {
    range(loop, first_sample, last_sample);
    min_intensity = findMinDB(min_array, first_sample, last_sample);

    sum_intensity += min_intensity;
  }

  return (float)sum_intensity/(max_pos - min_pos + 1);
}



//  calculate for all bases the ratio called/best_non_called and
//  calculate the frequency of bases with a ratio less that 2
int32 SCF_look::calledNoncalledRatioDB(int32 pos, int32 radius)
{
  uint32 minPos, maxPos;
  int32 first_sample, last_sample;
  int32 max_array[4];
  int32 bad = 0;
  float value = 0.0;

  // find the region where we want to examine the trace...
  minPos = max(2, pos - radius);
  maxPos = min(db_samples - 2, pos + radius);


  for (uint32 loop = minPos; loop <= maxPos; loop++) {
    int32 m0 = 1, m1 = 1, m2 = 1;

    half_range(loop, first_sample, last_sample);
    m0 = findMaxDB(max_array, first_sample, last_sample);

    // max_array contains the highest intensity value for the four bases
    // in a small vicinity of the called base

    for (int16 j = 0; j < 4; j++) {
      if (max_array[j] > m1) {
	m2 = m1;
	m1 = max_array[j];
      } else {
	if (max_array[j] > m2) {
	  m2 = max_array[j];
	}
      }
    }

    /*
      if ((float(m1) / float(m2)) < min_value) {
        min_value = float(m1) /float(m2);
      }
    */

    assert (m2 != 0);

    if ((float(m1) / float(m2)) < 2.0) {
      bad++ ;
    }

    value += (float(m1) / float(m2));
  }

  assert((1 + maxPos - minPos) != 0);

  return ((100* bad) / (1 + maxPos - minPos));
}


/*
   Calculate the mean ratio between called bases and non-called bases
   in the area around pos. r_min is the worst ratio of all four bases
*/


int32 SCF_look::noiseRatioDB(int32 pos, int32 radius,
			     int32 &c_ratio, int32 meanPeakSize[4])
{
  int32 minPos = pos - radius;
  int32 maxPos = pos + radius;

  if (minPos < 1) { minPos = 1; }
  if (maxPos >= db_samples) { maxPos = db_samples-1; }

  if (db_bases[pos] == UNDEFINED_BASE) {
    return meanCalledLevel(minPos, maxPos, pos, c_ratio, meanPeakSize);
  } else {
    return noiseRatioBaseDB(minPos, maxPos, pos, c_ratio, meanPeakSize);
  }

}



/* Calculate the mean level of a called base. */
/* Returns 10 * signal / noise */

int32 SCF_look::meanCalledLevel(int32 minPos, int32 maxPos, int32 pos,
				int32 &relMax, int32 meanPeakSize[4])
{
  int32 summe = 0;
  int32 nc_summe = 0;
  int32 count = 0;
  int32 max, max_array[4], noise;
  int32 pos_max  = 0;  // base_max for the base at pos
  int32 left, right;

  {
    int32 dummy_pos = pos;

    while (dummy_pos < db_samples &&
	   db_bases[dummy_pos] == UNDEFINED_BASE)
      dummy_pos++;
    noiseRatioBaseDB(minPos, maxPos, dummy_pos, relMax, meanPeakSize);
  }

  half_range(pos, left, right);
  pos_max = findMaxDB(max_array, left+1, right-1);

  for (int32 loop = minPos; loop <= maxPos; loop++) {
    findMaxDB(loop, max_array);
    max = 0; noise = 0;

    if (isRealBase(db_bases[loop])) {
      for (int16 k = 0; k < 4; k++) {
	if (max < max_array[k]) {
          noise = max;
	  max = max_array[k];
	} else {
	  if (max_array[k] > noise) noise = max_array[k];
	}
      }
      summe += max;  // We assume the called base to be the highest
      nc_summe += noise;
      count++;
    }
  }


  if (count > 0) {
    int16 mean_call = summe / count;

    if (mean_call == 0) {
      relMax = 100;
    } else {
      relMax = (100 * pos_max) / mean_call;
    }

#ifdef PARAMETER_VERBOSE
    cout << "Mean of all called bases: "<< mean_call << endl;
    cout << "Max of base at position : "<< pos_max  << endl;
    cout << "Ratio average-call/call : "<< relMax   << endl;
    if (nc_summe != 0) {
      cout << "c_ratio                 : "<< (10*summe) / nc_summe << endl;
    } else {
      cout << "c_ratio --- no noise found" << endl;
    }
#endif

    if (nc_summe == 0) {
      // no noise at all - assume signal/noise = 5
      return 50;
    } else {
      return (10 * summe) / nc_summe;
    }

  } else {
    relMax = 100;
    return -1;
  }

  /*
  RelMax    ist der prozentuale Anteil, einer mittleren gecallten Base,
            den das Maximum an der untersuchten Stelle erreicht hat.
  Rückgabe  das mit Verhältnis 10 * signal/noise
  */
}



void SCF_look::meanTraceValues(int32 mean[4], int32 offset)
{
  mean[0] = mean[1] = mean[2] = mean[3] = 0;

  if ((int32)SCF_header.samples <= offset) {
    return;
  }

  for (uint32 loop = offset; loop < SCF_header.samples-1; loop++) {
    mean[0] += (int32)SCF_samples_A[loop];
    mean[1] += (int32)SCF_samples_C[loop];
    mean[2] += (int32)SCF_samples_G[loop];
    mean[3] += (int32)SCF_samples_T[loop];
  }

  for (int32 loop = 0; loop < 4; loop++) {
    assert(offset != (int32)SCF_header.samples);

    mean[loop] = mean[loop] / (SCF_header.samples - offset);
  }
}


bool SCF_look::localPeakOverMean(char base, int32 scfPos, int32 radius,
				 int32 &above_mean)
{
  uint16 *the_trace;
  int32  von, bis, sum, position;

  switch (base) {
  case 'A': the_trace = SCF_samples_A; break;
  case 'C': the_trace = SCF_samples_C; break;
  case 'G': the_trace = SCF_samples_G; break;
  case 'T': the_trace = SCF_samples_T; break;
  default: return false;
  }

  if (scfPos < 0 || scfPos >= (int32)SCF_header.bases) {
    return false;
  }

  position = SCF_peak_index[scfPos];
  von = 0;
  if (position > radius) {
    von = position - radius;
  } else {
    von = 0;
  }

  if (position + radius < (int32)SCF_header.samples) {
    bis = position + radius;
  } else {
    bis = SCF_header.samples-1;
  }

  sum = 0;

  sum = the_trace[(SCF_peak_index[scfPos-1] + position) / 2] +
        the_trace[(SCF_peak_index[scfPos+1] + position) / 2];

  /*
  for (int32 loop = von; loop < bis; loop++) {
    sum += (int32)the_trace[loop];
  }
  */

  if (bis > von) {
    above_mean = the_trace[position] - (sum / 2);
    return true;
  } else {
    return false;
  }
}




int32 SCF_look::noiseRatioBaseDB(int32 minPos, int32 maxPos, int32 pos,
				 int32 &c_ratio, int32 meanPeakSize[4]) {
  int32 max_array[4], num_called[4], num_ncalled[4], r_min = 0;
  int32 sum_called[4], sum_ncalled[4];

  c_ratio = -1;

  for (int16 i = 0; i < 4; i++) {
    num_called[i]  = 0;
    num_ncalled[i] = 0;
    sum_called[i]  = 0;
    sum_ncalled[i] = 0;
  }

  if (pos >= db_samples-1 || pos < 0) {
     return 0;
  }

  if (maxPos >= db_samples - 2) maxPos = db_samples - 2;
  if (minPos <=  1) minPos = 2;

  //  cout << "noiseRatioBaseDB " << endl
  //     << "minPos " << minPos << "  MaxPos " << maxPos << endl;
  /* wir betrachten die intervalle [first_sample; last_sample[ */

  for (int32 loop = minPos; loop <= maxPos; loop++) {

    findMaxDB(loop, max_array);

    for (int16 j = 0; j < 4; j++) {
      int32 jj = isReversed() ? 3-j : j;

      if (db_bases[loop] == "ACGT"[jj]) {
	//cout << "Peak :" << loop << " " << max_array[j] << endl;
	num_called[jj]++;
	sum_called[jj] += max_array[j];
      } else {
	if (db_bases[loop+1] != "ACGT"[jj] &&
	    db_bases[loop-1]!= "ACGT"[jj]){
	  num_ncalled[jj]++;
	  // cout << " NC:" << max_array[j];
	  sum_ncalled[jj] += max_array[j];
	}
      }
    }
  }

  {
    float ratio = 0;
    float ratio_called  = 0;  // mean height of the called peaks
    float ratio_ncalled = 0;  // mean max. height at non called positions

    for (int16 j=0; j<4; j++) {
      int32 jj = (isReversed() ? 3-j : j);

      if (num_called[jj] > 0) {
	ratio_called = (float)sum_called[jj] / (float)num_called[jj];
	// cout << "Ratio-Called: " << ratio_called ;

	if (db_bases[pos] == "ACGT"[j]) {
	  // calculate peak height
	  findMaxDB(pos, max_array);

	  if (ratio_called == 0) {
	    // no base in the vicinity found... assume ratio 1 (100%)
	    c_ratio = 100;
	  } else {
	    c_ratio = (int16)(100 * max_array[jj] / ratio_called);
	  }
	}
      }
      meanPeakSize[j] = (int32)ratio_called;


      if (num_ncalled[j] > 0) {
	ratio_ncalled = (float)sum_ncalled[j] / (float)num_ncalled[j];
	// cout<<" Ratio N-called : " << ratio_ncalled;
      }

      if (ratio_ncalled > 0) {
	ratio = min(1000, (int)(10 * ratio_called / ratio_ncalled));
	// cout << " Ratio:" <<ratio << endl;
	if ((ratio < r_min || r_min == 0) && ratio > 0){
	  r_min = (int)ratio;
	}
      }

    }
  }
  return r_min;
}



// scfBasePos ist aber die position der gesuchten scf-base in der DB!!
char SCF_look::getSCFBase(int32 scfBasePos) const
{
  FUNCSTART("char SCF_look::getSCFBase(int32 scfBasePos) const");
  if (scfBasePos < 0 || scfBasePos >= db_samples) {
    cout << "scf_look::getSCFBase SCF-Baseposition exceeds limit\n";
    cout << "File: " << getFileName() << "\tpos:" << scfBasePos << endl;

    MIRANOTIFY(Notify::WARNING, "SCF-Baseposition exceeds limit");
  }

  if (isReversed()) {
    return invertBase(db_bases[scfBasePos]);
  } else {
    return db_bases[scfBasePos];
  }
}



char SCF_look::getOriginalSCFBase(int32 scfBasePos) const
{
  FUNCSTART("char SCF_look::getOriginalSCFBase(int32 scfBasePos) const");
  if (!testSCFPos(scfBasePos)) {

    cout << "scf_look::getOriginalSCFBase scf-Baseposition exceeds limit\n";
    cout << "File: " << getFileName() << "\tpos:" << scfBasePos << endl;

    MIRANOTIFY(Notify::WARNING,"SCF-Baseposition exceeds limit");
  }

  if (isReversed()) {
    return invertBase(SCF_bases[scfBasePos]);
  } else {
    return SCF_bases[scfBasePos];
  }
}



char SCF_look::getDBBase(int32 dbbasePos) const
{
  FUNCSTART("char SCF_look::getDBBase(int32 dbbasePos) const");
  if (testDBPos(dbbasePos)) {
    return db_bases[dbbasePos];
  }

  cout << "scf_look::getDBBase DB-Baseposition exceeds limit\n";
  cout << "File: " << getFileName() << "\tpos:" << dbbasePos << endl;

  MIRANOTIFY(Notify::WARNING,"DB-Baseposition exceeds limit");
}




int32 SCF_look::getDBBasePeakVal(int32 dbBasePos) const
{
  FUNCSTART("int32 SCF_look::getDBBasePeakVal(int32 dbBasePos) const");
  uint16 *samples;

  if (testDBPos(dbBasePos)) {
    samples = getTraceData(db_bases[dbBasePos]);
    return samples[db_peak_index[dbBasePos]];
  }

  cout << "scf_look::getDBBasePeakVal DB-Baseposition exceeds limit\n";
  cout << "File: " << getFileName() << "\tpos:" << dbBasePos << endl;

  MIRANOTIFY(Notify::WARNING, "DB-Baseposition exceeds limit");
}




int32 SCF_look::calculateCompressionParameters(int32 dbPos, int32 radius,
					       float *inParameter)
{
  int32 used = 0;
  float b1, b2;

  inParameter[used++] = noiseRatioAll(dbPos - radius, dbPos - 1);
  inParameter[used++] = noiseRatioAll(dbPos + 1, dbPos + radius);

  b1 = meanBaseline(dbPos - radius, dbPos - 1);
  b2 = meanBaseline(dbPos + 1, dbPos + radius);

  if (b2 != 0) {
    inParameter[used++] = b1/b2;
  } else {
    inParameter[used++] = b1/1;
  }

  return used;
}






// ===================================================

peak_examination::peak_examination(SCF_look *aRead, int32 aDbPos) {
  theRead = aRead;
  dbPos   = aDbPos;
  theBase = toupper(theRead->getDBBase(dbPos));

  nonBorderRealBase =
    isBase(theBase) &&
    dbPos > 0 &&
    dbPos < theRead->getDBBases()-2;

  return;
}


bool peak_examination::peakHeightWithPrevMin(int32 &peak_height)
{
  int32  leftValleyPos  = dbPos-1;
  int32  rightValleyPos = dbPos+1;
  int32  von, bis, min;
  uint16 *samples;

  if (!nonBorderRealBase) return false;

  // find next and previous base which is different from the called base
  // We expect the trace to go down to the base level here
  while(leftValleyPos > 1 &&
	toupper(theRead->getDBBase(leftValleyPos)) == theBase) {
    leftValleyPos--;
  }

  while(rightValleyPos < theRead->getDBBases() - 2 &&
	toupper(theRead->getDBBase(rightValleyPos)) == theBase) {
    rightValleyPos++;
  }

  samples = theRead->getTraceData(theBase);
  if (samples == nullptr) return false;


  //  theRead->range(leftValleyPos, rightValleyPos, von, bis);
  //  min = theRead->findMinBase(samples, von, bis);

  {
    int32 min1, min2;
    theRead->range(leftValleyPos, von, bis);
    min1 = theRead->findMinBase(samples, von, bis);
    theRead->range(rightValleyPos, von, bis);
    min2 = theRead->findMinBase(samples, von, bis);
    min = (min1 + min2) / 2;
  }

  peak_height = theRead->getDBBasePeakVal(dbPos) - min;
  return true;
}

/*
float peak_examination::getPeakWeight()
{
  char  prev_base, this_base;
  int32 i = 0;
  float weight[17] = {
    0.69, 1.00, 0.70, 0.98,
    0.77, 1.00, 0.61, 0.99,
    1.00, 0.51, 0.56, 0.94,
    0.55, 1.00, 0.91, 0.60,
    1.00
  };

  if (false == nonBorderRealBase) return 1.0;

  if (theRead->isReversed()) {
    prev_base = toupper(theRead->getDBBase(dbPos-1));
  } else {
    prev_base = toupper(theRead->getDBBase(dbPos+1));
  }

  switch(prev_base) {
  case 'A': i = 0;  break;
  case 'C': i = 4;  break;
  case 'G': i = 8;  break;
  case 'T': i = 12; break;
  default:  i = 16;
  }

  switch(theBase) {
  case 'A': break;
  case 'C': i = i+1; break;
  case 'G': i = i+2; break;
  case 'T': i = i+3; break;
  default:  i = 16;
  }

  return weight[i];
}

*/






// ========================================
// ==========    SCF_distance   ===========
// ========================================

float SCF_distance::getMeanProblemDistance()
{
  if (problemBases > 1) {
    return (float)(problemDistanceSum / (problemBases - 1));
  } else {
    return (float)problemDistanceSum;
  }
}


float SCF_distance::getMeanVicDistance()
{
  if (vicDistanceSum != 0) {
    if (vicDistanceBases > 1) {
      return (float)(abs (vicDistanceSum) / vicDistanceBases);
    } else {
      return abs(vicDistanceSum);
    }
  } else {
    return 10.0;
    // Quit arbitrary value; should never happen...
  }
}


float SCF_distance::getRelativeProblemDistance()
{
  if (problemBases != 0 && vicDistanceSum != 0) {
    return (float)(problemDistanceSum * vicDistanceBases) /
      (float)(problemBases * vicDistanceSum);
  } else {
    return 1.0;
  }
}



float SCF_distance::getRelativeMinDistance()
{
  float m = getMeanVicDistance();

  if (m > 0.001) {
    return (float)minDistanceProblem / getMeanVicDistance();
  }
  return 1.0;
}



float SCF_distance::getRelativeMaxDistance()
{
  float m = getMeanVicDistance();

  if (m > 0.001) {
    return (float)maxDistanceProblem / getMeanVicDistance();
  }
  return 1.0;
}


float SCF_distance::getMinDistanceScore()
{
  float d1 = 50 * pow(getRelativeMinDistance() - 1, 4);

  return d1;
}

float SCF_distance::getMaxDistanceScore()
{
  float d2 = 50 * pow(getRelativeMaxDistance() - 1, 4);

  return d2;
}


float SCF_distance::getDistanceScore()
{
  float d3 = 50 * pow(getRelativeProblemDistance() - 1, 4);

  return d3;
}



ostream & operator<<(ostream &ostr, SCF_distance const &i) {
  ostr << "Peakpos:    " << i.peakPos << " Bases: " << i.problemBases << endl;
  ostr << "Problem Min " << i.minDistanceProblem
       << " Max: " << i.maxDistanceProblem << endl;
  ostr << "Vicinity Bases: " << i.vicDistanceBases << " Sum "
       << i.vicDistanceSum << endl;

  return ostr;
}
