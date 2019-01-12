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
 * SCF examination methods - misc.C
 *
 * Written by Thomas Pfisterer
 *
 * 06.07.99  Functions to handle/check bases etc. that do not belong
 *           to any class
 * 15.07.99  new functions isUndefinedBase/indirectSort4. isRealBase
 *           also with lowercase letters.
 */

#include "misc.H"
#include "assert.h"
#include "math.h"



// Sorts array "pos" such that references to "a" are sorted (descending)
// Is used to find the n-th-greatest value in array "a"

void indirectSort4(int32 a[4], int32 pos[4])
{
  for (int32 l=0; l<4; l++) pos[l] = l;

  for (int32 l1 = 0; l1 < 3; l1++) {
    for (int32 l2 = l1; l2 < 4; l2++) {
      if (a[pos[l1]] < a[pos[l2]]) {
	int32 d = pos[l1]; pos[l1] = pos[l2]; pos[l2] = d;
      }
    }
  }
}


// Given the base we want to alter into another base and the bases before
// and after the base this function calculates the corresponding fault class


int32 findFaultClass(const char prev, const char aBase,
		     const char next, const char newBase)
{

  if (newBase == GAP_BASE || newBase == '.') {
    if (aBase == GAP_BASE) return fhc_GAP;

    if (prev == aBase || next == aBase) {
      return fhc_OVERCALL;
    } else {
      return fhc_ADDITIONAL;
    }
  }

  if (aBase == GAP_BASE) {
    if (prev == newBase || next == newBase) {
      return fhc_UNDERCALL;
    } else {
      return fhc_MISSING;
    }
  }

  if (aBase == newBase) {
    return fhc_CORRECT;
  } else {
    return fhc_WRONG;
  }

}


// Find the fault class (and neural net) for the given fhc fault-class
// and the old base in the DB. Differences are e.g. altering an
// undefined base which has a seperate NNClass.

int32 faultClassToNNClass(const int32 faultClass, const char oldDBBase)
{
  int32 nnFaultClass;

  switch (faultClass) {
  case fhc_WRONG:
    if (isRealBase(oldDBBase)) {
      nnFaultClass = WRONGCALL_INDEX;
    } else {
      if (isUndefinedBase(oldDBBase)) {
	nnFaultClass = N_CALL_INDEX;
      } else {
	nnFaultClass = UNDEFINED_INDEX;
      }
    }
    break;
  case fhc_UNDERCALL:
    nnFaultClass = UNDERCALL_INDEX;
    break;
  case fhc_OVERCALL:
    if (isRealBase(oldDBBase)) {
      nnFaultClass = OVERCALL_INDEX;
    } else {
      if (isUndefinedBase(oldDBBase)) {
	nnFaultClass = N_PLUS_INDEX;
      } else {
	nnFaultClass = UNDEFINED_INDEX;
      }
    }
    break;
  case fhc_MISSING:
    nnFaultClass = MISSING_INDEX;
    break;

 case fhc_ADDITIONAL:
    if (isRealBase(oldDBBase)) {
      nnFaultClass = ADDITIONAL_INDEX;
    } else {
      if (isUndefinedBase(oldDBBase)) {
	nnFaultClass = N_PLUS_INDEX;
      } else {
	nnFaultClass = UNDEFINED_INDEX;
      }
    }

    break;
    // case fhc_GAP:
    // output << "nnFaultClass : Gap " << endl;
  default:
    nnFaultClass = UNDEFINED_INDEX;
  }

 return nnFaultClass;
}




int16 isRealBase(const char c)
{
  static char flags[128];
  static int8 is_init = 0;

  if (is_init == 0) {
    for (int32 loop = 0; loop < 128; loop++) { flags[loop] = 0; }
    flags['A'] =  flags['C'] = flags['G'] = flags['T'] = 1;
    flags['a'] =  flags['c'] = flags['g'] = flags['t'] = 1;
  }

  return flags[c];
}



int16 isBase(const char c)
{
  static char flags[128];
  static int8 is_init = 0;

  if (is_init == 0) {
    for (int32 loop = 0; loop < 128; loop++) { flags[loop] = 0; }
      flags['A'] = flags['C'] = flags['G'] = flags['T'] = 1;
      flags['a'] = flags['c'] = flags['g'] = flags['t'] = 1;
      flags['N'] = flags['n'] = flags['-'] = 1;
  }

  return flags[c];
}



int16 isUndefinedBase(const char c)
{
  return (c == '-' || c == 'N');
}



char invertBase(const char c)
{
  static char flags[128];
  static int8 is_init = 0;

  if (0 == is_init) {
    for (int32 loop = 0; loop < 128; loop++) { flags[loop] = 'N'; }
    flags['A'] = flags['a'] = 'T';
    flags['C'] = flags['c'] = 'G';
    flags['G'] = flags['g'] = 'C';
    flags['T'] = flags['t'] = 'A';
    flags['*'] = '*';
  }

 return flags[c];
}


int32 valueForBase(const int32 theArray[4], const char theSCFBase) {
  switch (theSCFBase) {
  case 'A':
  case 'a': return theArray[0];
  case 'C':
  case 'c': return theArray[1];
  case 'G':
  case 'g': return theArray[2];
  case 'T':
  case 't': return theArray[3];
  default:
    return ((theArray[0]+theArray[1]+theArray[2]+theArray[3])/4);
  }
}


char toDBBase(const char c, const int32 reversed)
{
  if (reversed == 0) return c;
  return invertBase(c);
}


char toSCFBase(const char c, const int32 reversed)
{
  if (reversed == 0) return c;
  return invertBase(c);
}


// *******************************************************
//    calledBase
//    simple method to calculate the consensus base.
//    actually the calculated base is not that important
//    ther is always a difference to one of the reads if
//    the reads don't agree unanimously
// *******************************************************


char calledBase(const Contig::consensus_counts_t &cc)
{
  if (cc.A > cc.C && cc.A > cc.G && cc.A > cc.T) { return 'A';};
  if (cc.C > cc.A && cc.C > cc.G && cc.C > cc.T) { return 'C';};
  if (cc.G > cc.A && cc.G > cc.C && cc.G > cc.T) { return 'G';};
  if (cc.T > cc.A && cc.T > cc.C && cc.T > cc.G) { return 'T';};
  return 'N';
}



void initializeGaussPattern(uint16* pattern, int32 size, double sigma)
{
  double offset = - ((double)size - 1) / 2;

  if (size < 1) {
    return;
  }

  for (int32 loop = 0; loop < size; loop++) {
    pattern[loop] = (uint16)(0.5 + 100 *
			     exp(- fabs((offset+loop)/sigma)));
  }

}



// Simple Korrelation of a small pattern in a larger pattern

double simpleNormedCorrelation(uint16 *sequence, int32 vonpos,
			      int32 bispos, double sigma)
{
  int32   patternsize = 9;
  //uint16  pattern[15];

  //  initializeGaussPattern(pattern, patternsize, sigma);

  // Dreieck 1
  //uint16  pattern[9] = { 0, 1, 2, 3, 4, 3, 2, 1, 0};

  // Dreieck 2
  // uint16  pattern[9] = { 1, 2, 3, 4, 5, 4, 3, 2, 1};

  // Gauss Sigma=1  (Gauss 1)
  // uint16  pattern[9] = { 2, 5, 14, 37, 100, 37, 14, 5, 2};

  // Gauss Sigma=4/3 (Gauss 11)
  //  uint16  pattern[9] = {11, 22, 47, 100, 47, 22, 11};

  // Gauss Sigma=2  (Gauss 2)
  // uint16  pattern[9] = { 14, 22, 37, 60, 100, 60, 37, 22, 14 };

  // Gauss Sigma=2.5  (Gauss 21)
  //  uint16  pattern[9] = { 20, 30, 45, 67, 100, 67, 45, 30, 20 };

  // Gaus Sigma=3   (Gauss 3)
  uint16 pattern[9] = {26, 37, 51, 72, 100, 72, 51, 37, 26};

  int32  maxpos = -1;
  int32  k, kmax = -1;
  int32  maximum = 1;

  if (sequence == nullptr) {
    cout << "No Trace!!" << endl;
    return 0;
  }


  for (int32 loop = vonpos; loop <= bispos; loop++) {
    if (sequence[loop] > maximum) { maximum = sequence[loop]; }
  }

  for (int32 loop = vonpos; loop+patternsize-1 <= bispos; loop++) {
    k = 0;
    for (int32 l = 0; l<patternsize; l++) {
      k = k + (pattern[l] * (int32)sequence[loop+l]);
    }
    if (k > kmax) {
      kmax = k;
      maxpos = loop;
    }
  }

  if (maxpos < 0) {
    return 0;
  }

  double ac1 = autoCorrelation(sequence, maxpos, maxpos+patternsize, 0);
  double ac2 = autoCorrelation(pattern, 0, patternsize-1, 0);

  // cout << "Autokorrelation Signal  : " << ac1 << endl;
  // cout << "Autokorrelation Pattern : " << ac2 << endl;
  double norm = sqrt(ac1 * ac2);

  if (norm <= 0.001) {
    return 0;
  }

  //#ifndef RUNONLY
  //if (norm < kmax) {
  //  cout << "Norm: " << norm << "\t kmax " << kmax << endl;
  //}

  if (norm <= 0.01 || kmax/norm > 1.01) {
    return 0;
  }
  // Norm should be in [0;1] (with a little epsilion)
  // assert(kmax/norm < 1.1);
  //#endif


  // Value betwen 0 and 100; the higher the value the better the korrelation

  return kmax/norm;
}




double autoCorrelation(uint16 *sequence, int32 vonpos, int32 bispos, int32 o)
{
  int32  offset = o;
  double s = 0;

  if (o < 0) offset = -o;

  if (bispos < vonpos) {
    int32 dummy = bispos; bispos = vonpos; vonpos = dummy;
  }

  //cout << "Autokorrelation " << vonpos << "-" << bispos << endl;

  for (int32 loop=vonpos; loop<=bispos; loop++) {
    s = s + (double)(sequence[loop]*sequence[loop+offset]);
  }

  return s;
}
