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
 */

#include "scf_look.H"



// area around a single base (halfway to adjacent peaks)
void SCF_look::range(const int32 basePos, int32 &von, int32 &bis)
{

  if (basePos <= 0) {
    von = 0;
  } else {
    if (basePos < db_samples-1) {
      // skip the gaps
      int32 i = 1;
      while (db_bases[basePos-i] == '*' && basePos > i) i++;

      von = (db_peak_index[basePos-i] + db_peak_index[basePos]) / 2;
    } else {
      von = db_peak_index[db_samples-1];
    }
  }

  if (basePos < db_samples - 1) {
    if (basePos > 0) {
      // skip the gaps
      int32 i = 1;
      while (db_bases[basePos+i] == '*' && basePos+i < db_samples-1) i++;

      bis = (db_peak_index[basePos+i] + db_peak_index[basePos]) / 2;
    } else {
      bis = db_peak_index[0];
    }
  } else {
    bis = db_peak_index[db_samples-1];
  }

#ifdef RANGECHECK_ON
  if (von < 0 ||  von > (int32)SCF_header.samples) {
    cerr << "VON ist out of range! 3 " << endl;
    //    cout << "VON ist out of range! 3 " << von << endl;
    throw Notify(Notify::WARNING, "range", "von out-of range");
  }
  if (bis < 0 ||  bis > (int32)SCF_header.samples) {
    cerr << "BIS ist out of range! 3 " << endl;
    //    cout << "BIS ist out of range! 3 " << bis << endl;
    //cerr << "Index: " << db_peak_index[basePos+1] << " ";
    //cerr <<  db_peak_index[basePos];
    //cerr << "BasePos: " << basePos << endl;
    //cout << "Bases in Read: " << db_samples << endl;
    throw Notify(Notify::WARNING, "range", "bis out of range");
  }
#endif

  if (von > bis) { von ^= bis ^= von ^= bis; }
}



// region around the bases (halfway to the adjacent peaks)
void SCF_look::range(const int32 dbPosVon, const int32 dbPosBis,
		     int32 &von, int32 &bis)
{

  if (dbPosVon < 1) {
    von = (db_peak_index[0]);
  } else {
    if (dbPosVon >= db_samples-1) {
      von = db_peak_index[db_samples-1];
    } else {
      // skip the gaps
      int32 i=1;
      while (db_bases[dbPosVon-i] == '*' && dbPosVon > i) i++;

      von = (db_peak_index[dbPosVon-i] + db_peak_index[dbPosVon]) / 2;
    }
  }

  if (dbPosBis >= (db_samples-1)) {
    bis = db_peak_index[db_samples-1];
  } else {
    if (dbPosBis < 0) {
      bis = db_peak_index[0];
    } else {
      // skip the gaps
      int32 i=1;
      while (db_bases[dbPosBis+i] == '*' && dbPosBis+i < db_samples-1) i++;

      bis = (db_peak_index[dbPosBis+i] + db_peak_index[dbPosBis]) / 2;
    }
  }

#ifdef RANGECHECK_ON
  if (von < 0 || bis < 0 ||
      von >= (int32)SCF_header.samples ||
      bis >= (int32)SCF_header.samples) {

    cerr << "basepos" << dbPosVon << "-" << dbPosBis << "  range "
	 << von << "-" << bis
	 << " #" << db_samples << endl;
    throw Notify(Notify::WARNING, "range.", "von/bis out of range");
  }
#endif

  if (von > bis)  { von ^= bis ^= von ^= bis; }
}


// narrow region around a base (one third to the next base)
void SCF_look::half_range(const int32 basePos, int32 &von, int32 &bis)
{
  int32 pos = basePos;

 if (basePos > (db_samples-1)) { pos = db_samples - 2; }

 if (pos < 1) {
    //    von = (0 + 2 * db_peak_index[1]) / 3;
    von = db_peak_index[0];
  } else {
    // skip the gap
    int32 i=1;
    while (db_bases[pos-i] == '*' && pos > i) i++;

    von = (db_peak_index[pos-i] + 2 * db_peak_index[pos]) / 3;
  }

  if (pos >= (db_samples-1)) {
    bis = db_peak_index[db_samples-1];
  } else {
    // skip the gap
    int32 i=1;
    while (db_bases[pos+i] == '*' && pos+i < db_samples-1) i++;

    bis = (db_peak_index[pos+i] + 2 * db_peak_index[pos]) / 3;
  }

#ifdef RANGECHECK_ON
  if (von < 0 || bis < 0 ||
      von >= (int32) SCF_header.samples ||
      bis >= (int32) SCF_header.samples) {

    cerr << "Out of Range: " << basePos << "half_range " << von << "-" << bis
	 << " #" << db_samples << endl;
    throw Notify(Notify::WARNING, "half_range", "von/bis out of range");
  }
#endif

  if (von > bis) { von ^= bis ^= von ^= bis; }
}


// wide region around the bases (two third of the distance to the
// adjacent peaks
void SCF_look::wide_range(const int32 dbPosVon, const int32 dbPosBis,
			  int32 &von, int32 &bis)
{
  if (dbPosVon < 1) {
    von = db_peak_index[0];
  } else {
    if (dbPosVon < (db_samples-1)) {
      // skip the gaps
      int32 i=1;
      while (db_bases[dbPosVon-i] == '*' && dbPosVon > i) i++;

      von = ((2* db_peak_index[dbPosVon-i]) + db_peak_index[dbPosVon]) / 3;
    } else {
      von = db_peak_index[db_samples-1];
    }
  }

  if (dbPosBis >= (db_samples-1)) {
    bis = db_peak_index[db_samples-1];
  } else {
    if (dbPosBis > 0) {
      // skip the gaps
      int32 i=1;
      while (db_bases[dbPosBis+i] == '*' && dbPosBis+i < db_samples-1) i++;

      bis = ((2* db_peak_index[dbPosBis+i]) + db_peak_index[dbPosBis]) / 3;
    } else {
      bis = db_peak_index[0];
    }
  }

#ifdef RANGECHECK_ON
  if (von < 0 || bis < 0 ||
      von >= (int32) SCF_header.samples ||
      bis >= (int32) SCF_header.samples) {

    cerr << "Out of range half_range " << von << "-" << bis
	 << " #" << db_samples << endl;
    throw Notify(Notify::WARNING, "half_range", "von/bis out of range");
  }
#endif

  if (von > bis) { von ^= bis ^= von ^= bis; }
}


// finds the maximum of each base at the base at "basePos"
// e.g. max[0] = highest value of base A around basePos
// returns: overall maximum
uint16 SCF_look::findMaxDB(int32 basePos, int32 max[4]) {
  int32 first_sample;
  int32 last_sample;

  range (basePos, first_sample, last_sample);
  return findMaxDB(max, first_sample, last_sample);

}



// same as findMaxDB with a smaller environment around the base
uint16 SCF_look::findMaxDBHalf(int32 basePos, int32 max[4]) {
  int32 first_sample;
  int32 last_sample;

  half_range (basePos, first_sample, last_sample);
  return findMaxDB(max, first_sample, last_sample);
}


// finds the maximum of each base in a range of sample values
// e.g. max[0] = highest value of "A"
// returns: overall maximum

uint16 SCF_look::findMaxDB(int32 max[4],
			   int32 first_sample, int32 last_sample) {
  uint16 max_all = 0;

  max[0] = max[1] = max[2] = max[3] = 0;

  max[0] = findMaxBase(SCF_samples_A, first_sample, last_sample);
  max[1] = findMaxBase(SCF_samples_C, first_sample, last_sample);
  max[2] = findMaxBase(SCF_samples_G, first_sample, last_sample);
  max[3] = findMaxBase(SCF_samples_T, first_sample, last_sample);


  for (int16 loop = 0; loop < 4; loop++) {
    if (max_all < max[loop]) max_all = max[loop];
  }

  return max_all;
}


// finds the minimum of each base in a range of sample values
// e.g. max[0] = smallest value of "A"
// returns: overall minimum

uint16 SCF_look::findMinDB(int32 min[4],
			   int32 first_sample, int32 last_sample) {
  uint16 min_all = 0;

  min[0] = min[1] = min[2] = min[3] = 0;

  min[0] = findMinBase(SCF_samples_A, first_sample, last_sample);
  min[1] = findMinBase(SCF_samples_C, first_sample, last_sample);
  min[2] = findMinBase(SCF_samples_G, first_sample, last_sample);
  min[3] = findMinBase(SCF_samples_T, first_sample, last_sample);


  for (int16 loop = 0; loop < 4; loop++) {
    if (min_all < min[loop]) min_all = min[loop];
  }

  return min_all;
}


/* Finds the maximum value of a base in a given range */
uint16 SCF_look::findMaxBase(uint16 *SCF_samples, int32 von, int32 bis) {
  uint16 max = 0;

  if (von > bis) { von ^= bis ^= von ^= bis; }

#ifdef RANGECHECK_ON
  if (von < 0 ||  von > (int32) SCF_header.samples) {
    cerr << "VON ist out of range! 2 " << endl;
    //cout << "VON ist out of range! 2 " << von << endl;
    //cout << "Von: " << von << " bis: " << bis << endl;
    throw Notify(Notify::WARNING, "range", "von out-of range");
  }
  if (bis < 0 ||  bis > (int32) SCF_header.samples) {
    cerr << "BIS ist out of range! 2 " << endl;
    //cout << "BIS ist out of range! 2 " << bis << endl;
    //cout << "Von: " << von << " bis: " << bis << endl;
    throw Notify(Notify::WARNING, "range", "bis out of range");
  }
#endif

  for (int32 loop = von; loop < bis; loop++) {
    if (SCF_samples[loop] > max) {
      max = SCF_samples[loop];
    }
  }
  return max;
}



/* Finds the minumum value of a base in a given range */
uint16 SCF_look::findMinBase(uint16 *SCF_samples, int32 von, int32 bis) {
  uint16 min = 10000;
  int16 step = 1;

  if (von > bis) { von ^= bis ^= von ^= bis; }

#ifdef RANGECHECK_ON
  if (von < 0 ||  von > (int32) SCF_header.samples) {
    cerr << "VON ist out of range! 1" << endl;
    //cout << "VON ist out of range! 1" << von << endl;

    throw Notify(Notify::WARNING, "range", "von out-of range");
  }
  if (bis < 0 ||  bis > (int32)SCF_header.samples) {
    cerr << "BIS ist out of range! 1" << endl;
    //cout << "BIS ist out of range! 1 " << bis << endl;
    throw Notify(Notify::WARNING, "range", "bis out of range");
  }
#endif

  for (int32 loop = von; loop < bis; loop += step) {
    if (SCF_samples[loop] < min) {
      min = SCF_samples[loop];
    }
  }
  return min;
}


// ==========================================================
// ==========================================================

// Searches the peak with the greatest distance to the borders of the
// searched intervall [von;bis] (sample positions)
int32 bestPeakBetweenBases(uint16 *samples, int32 von, int32 bis)
{
  int i;
  int bestPeak = 0;

  if (samples == nullptr) {
    DEBUG_EDIT("bestPeakBetweenBases: No trace!! " << endl);
    return 0;
  }

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap

  i = von+1;
  while (i < bis) {
    if (samples[i-1] < samples[i] && samples[i+1] <= samples[i]) {

      // Maximum
      if (bestPeak != 0) {
	if (min(bestPeak - von, bis - bestPeak) < min(i - von, bis - i)) {
	  bestPeak = i;
	}
      } else {
	bestPeak = i;
      }
    }
    i++;
  }

  return bestPeak;
}


// Searches the peak with the greatest distance to the borders of the
// searched intervall [von;bis] (sample positions)
int32 bestValleyBetweenBases(uint16 *samples, int32 von, int32 bis)
{
  int i;
  int bestValley = 0;

  if (samples == nullptr) {
    DEBUG_EDIT("bestValleyBetweenBase: No trace!! " << endl);
    return 0;
  }

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap

  i = von+1;
  while (i < bis) {
    if (samples[i-1] > samples[i] && samples[i+1] >= samples[i]) {
      // Minimum
      if (bestValley != 0) {
	if (min(bestValley - von, bis - bestValley) < min(i-von, bis-i)) {
	  bestValley = i;
	}
      } else {
	bestValley = i;
      }
    }
    i++;

  }
  return bestValley;
}



bool slopeChangeTest(const uint16 *samples, const int32 pos)
{
  float slopeDiff, meanSlope;

  slopeDiff = 2.0 * samples[pos] - samples[pos-1] - samples[pos+1];
  meanSlope = (float)(samples[pos+1] - samples[pos-1]) / 2.0;

  if (meanSlope == 0.0) {
    return !(slopeDiff == 0);
  }

  if ((slopeDiff / meanSlope) < 0) {
    return ((slopeDiff / meanSlope) < -0.75);
  } else {
    return ((slopeDiff / meanSlope) > 0.75);
  }
}



int32 findBestWeakPeak(const uint16 *samples, int32 von, int32 bis,
		       int32 &peak, int32 &valley)
{
  int32 i, mitte;
  int32 bestC, bestD;

  if (samples == nullptr) {
    DEBUG_EDIT("findBestWeakPea: No trace !!" << endl);
    return 0;
  }

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap
  i = von + 1;
  mitte = (von + bis + 1) / 2;
  bestC = bestD = 0;

  while (i < bis) {
    if(slopeChangeTest(samples, i)) {

	if ((samples[i-1] - 2*samples[i] + samples[i+1]) < 0) {
	  if (bestC == 0 || (abs(bestC - mitte) > abs(i - mitte))) {
	    bestC = i;
	  }
	} else {
	  if (bestD == 0 || (abs(bestD - mitte) > abs(i - mitte))) {
	    bestD = i;
	  }
	}
    }
    i++;
  }

  peak = bestC;
  valley = bestD;

  return 1;
}


int32 constSlopeBetweenBases(const uint16 *samples, int32 von,
			     int32 bis, int32 &bestC, int32 &bestD)
{
  int32 i;

  if (samples == nullptr) {
    DEBUG_EDIT("constSlopeBetweenBases: No trace!! " << endl);
    return 0;
  }

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap

  bestC = bestD = 0;
  i = von + 1;

  while (i < bis) {

    if (!slopeChangeTest(samples, i-1) &&
	slopeChangeTest(samples, i)  &&
	!slopeChangeTest(samples, i+1)) {


      if ((2 * samples[i] - samples[i-1] - samples[i+1]) > 0) {
	if (bestC != 0) {
	  if (min(bestC - von, bis - bestC) < min(i-von, bis-i)) {
	    bestC = i;
	  }
	} else {
	  bestC = i;
	}
      } else {
	if (bestD != 0) {
	  if (min(bestD - von, bis - bestD) < min(i-von, bis-i)) {
	    bestD = i;
	  }
	} else {
	  bestD = i;
	}
      }

    }
    i++;

  }
  return bestC;
}


int32 bestConvergenceRegion(uint16 *samples, int32 von, int32 bis)
{
  int32 regStart = -1;
  int32 regEnd   = -1;
  int32 bestCV = 0;
  int32 i;

  if (samples == nullptr) {
    DEBUG_EDIT("bestConvergenceRegion: No trace !!" << endl);
    return 0;
  }

  if (von > bis) { von ^= bis ^= von ^= bis; }  // swap

  i = von + 1;
  while (i < bis) {
    if (samples[i-1] + samples[i+1] < 2*samples[i]) {

      // Convergence
      if (i > regEnd+1) {  // gap or start
	if (regStart > 0) { // gap => old region finished

	  int32 mitte = (regStart + regEnd + 1) / 2;
	  if (bestCV == 0) bestCV = mitte;

	  if (min(bestCV - von, bis - bestCV) <
	      min(mitte - von, bis - mitte))
	    {
	      bestCV = mitte;
	    }
	}
	regStart = i;
      }
      regEnd = i;
    }
    i++;
  }

  if (regStart > 0 && regStart < bis-1) {
    if (bestCV == 0 ||
	min(bestCV - von, bis - bestCV) < min(i - von, bis - i)) {
      bestCV = (regStart + regEnd) / 2;
    }
  }

  return bestCV;
}
