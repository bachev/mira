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
 */


#ifndef _scflook_scflook_H_
#define _scflook_scflook_H_



// forward deklarationen
class SCF_buffer;
class afh_Info;
//struct contigread_t;

#include <string>

#include <math.h>
#include "ctype.h"

#include "io/scf.H"
#include "errorhandling/errorhandling.H"
#include "mira/read.H"
#include "mira/contig.H"

#include "buffer.H"
#include "prescan.H"
#include "globals.H"
#include "misc.H"


class SCF_distance {
private:
  int32 peakPos;

  int32 vicDistanceSum;
  int32 vicDistanceBases;

  int32 problemDistanceSum;
  int32 problemBases;
  int32 minDistanceProblem;
  int32 maxDistanceProblem;

public:
  SCF_distance(int32 pos) {peakPos = pos; }
  void setMinMaxProblemDistance(int32 min, int32 max)
    { minDistanceProblem = min; maxDistanceProblem = max; }
  void setMeanVicDistance(int32 distance, int32 bases)
    { vicDistanceSum = abs(distance); vicDistanceBases = bases; }
  void setProblemDistance(int32 distance, int32 bases)
    { problemDistanceSum = abs(distance); problemBases = bases; }

  int32 getProblemBases() { return problemBases; }
  float getMeanProblemDistance();
  float getMeanVicDistance();
  float getRelativeProblemDistance();
  float getRelativeMinDistance();
  float getRelativeMaxDistance();
  float getDistanceScore();
  float getMinDistanceScore();
  float getMaxDistanceScore();

  friend ostream & operator<<(ostream &ostr, SCF_distance const &i);
};



/* Create a class that can handle SCF-Files and more */
class SCF_look : public SCF {
  static int32 count;
private:
  enum { qualityRadius = 15 };
  int8   reversed;      /* read is the other way round in the db */
  int32  db_samples;
  const Read   * theRead;
  uint8  * db_bases;
  uint32 * db_peak_index;
  char   * fName;


  void  zeroVars();
  void  discard();

  void  wide_range(const int32 dbPosVon, const int32 dbPosBis,
		   int32 &von, int32 &bis);

  int32 meanCalledLevel(int32 minPos, int32 maxPos, int32 pos,
			int32 &relMax, int32 meanPeakSize[4]);
  int32 noiseRatioBaseDB(int32 minPos, int32 maxPos, int32 pos,
			 int32 &c_ratio, int32 meanPeakSize[4]);
  int32 countPeaksInTrace(uint16 *samples, int32 von, int32 bis,
			  int32 &peak, int16 &qual);
  void  interpolate();
public:
  SCF_look();
  ~SCF_look();
  SCF_look(SCF_look const &other); //: SCF(other);

  const SCF_look& operator=(const SCF_look &other);

  int8  load(const char *filename);
  int8  simpleLoad(const char *filename);
  int8  load(const Read &aRead, int32 richtung);
  int8  load();

  char   *getFileName() const;
  int8    setFilename(const char *filename);
  t_machinetype getMACHType() const;

  uint16 *getTraceData(const char dbbase) const;
  uint16 *getTraceData(const char dbbase, int32 isReversed) const;
  void    range(int32 basePos, int32 &von, int32 &bis);
  void    range(const int32 dbPosVon, const int32 dbPosBis,
	        int32 &von, int32 &bis);
  void   half_range(int32 basePos, int32 &von, int32 &bis);
  uint16 findMaxDB(int32 basePos, int32 max[4]);
  uint16 findMaxDB(int32 max[4], int32 first, int32 last);
  uint16 findMinDB(int32 min[4], int32 first, int32 last);
  uint16 findMaxDBHalf(int32 basePos, int32 max[4]);
  uint16 findMaxBase(uint16 *SCF_samples, int32 von, int32 bis);
  uint16 findMinBase(uint16 *SCF_samples, int32 von, int32 bis);


  void   meanTraceValues(int32 mean[4], int32 offset);
  bool   localPeakOverMean(char base, int32 scfPos, int32 radius,
			   int32 &above_mean);

  int32  findNextMax(const uint16 *samples, const int32 samplePos,
		     const bool forward);
  int32  findNextMin(const uint16 *samples, const int32 samplePos,
		     const bool forward);
  int32  findNextSlopeChange(const uint16 *samples, const int32 samplePos,
			     const bool forward);

  float  peakValleyValues(const int32 dbPos, const char newBase,
			  float result[6]);

  int32  findRepeatSize(int32 &dbPosVon, int32 &dbPosBis);

  SCF_distance* calcDBDistance(int32 dbPos);

  int32  calcConvergenceRegionQuality(int32 dbPos);
  int8   isReversed() const { return reversed; }



  int32  convergenceRegions(uint16 *samples, int32 von, int32 bis,
			    int32 &areaMin, int32 &areaMax);
  //  int32  countConvergenceRegions(int32 dbBasePos,
  //				 int32 &areaMin, int32 &areaMax);

  int32  countConvergenceRegions(int32 dbBasePos, char base,
				 int32 &areaMin, int32 &areaMax);

  float  einzelpeak(int32 dbPos, char dbBase, int32 &max,
		    float &pos_rating);

  float  singlePeakQuality(int32 dbpos, float sigma, char base = '-');

  void   peakRelations(int32 dbPos, char newDBBase, char oldDBBase,
		       float &newRatio, float &oldRatio, float &oldNewRatio);

  int32  countPeaksBetweenBases(char base, int32 basPos);


  int32 countRepeatPeaks(const int32 basePos, const char base,
			 int16 &anzRepeat, int16 &anzPeaks, int16 &qual);

  int32  countDBPeaks(int32 basePos, int16 &anzRepeat,
  		      int16 &anzPeak, int16 &qual);
  int16 peakQuality(int32 pos, int32 radius);

  float noiseRatioAll(uint32 von, uint32 bis);
  float meanBaseline(uint32 min_pos, uint32 max_pos);

  int32 calledNoncalledRatioDB(int32 pos, int32 radius);
  int32 noiseRatioDB(int32 pos, int32 radius, int32 &c_ratio,
		     int32 meanPeakSize[4]);
  void  showDBArea(int32 pos, int32 circle, ostream &output);
  void  showArea(int32 pos, int32 circle, ostream &output);

  char  getSCFBase(int32 scfBasePos) const;
  //base_quality_t getDBQuality(int32 dbpos) const;
  char  getOriginalSCFBase(int32 scfBasePos) const;
  char  getDBBase(int32 dbBasePos) const;
  int32 getDBBasePeakVal(int32 dbBasePos) const;
  int32 getDBBases() const { return db_samples; }

  inline bool testDBPos(int32 dbpos) const
    { return (dbpos >= 0 && dbpos < db_samples); }
  inline bool testSCFPos(int32 scfpos) const
    { return (scfpos >= 0 && scfpos < static_cast<int32>(SCF_header.bases)); }

  int32  deleteBaseDB(int32 pos);
  int32  insertBaseDB(int32 pos, int32 scfPos, char newBase);
  int32  alterBaseDB(int32 pos, char newBase);


  int32 calculateCompressionParameters(int32 dbPos, int32 radius, float *inParameter);

};



class peak_examination {
  SCF_look  *theRead;
  int32     dbPos;
  bool      nonBorderRealBase;
  char      theBase;

public:
  peak_examination(SCF_look *aRead, int32 aDbPos);

  bool  peakHeightWithPrevMin(int32 &peak_height);
  float getPeakWeight();

};





class afh_Info {
  enum { qualityRadius = 15 };
  afh_Info *next;

  SCF_look *theTrace;

  const Contig::contigread_t *contig_read;
  base_quality_t original_quality;

  int32 question;
  char  hypothesis_base;    // for inserts and wrong calls
  char  original_base;
  int32 status;
  float confirmed;
  int32 dbPosition;
  int32 scfPosition;



  void setDBPos(int32 newPos);
  void setSCFPos(int32 newPos);
  void shiftLeftAfterPos(int32 startPos);
  void shiftRightAfterPos(int32 startPos);

public:
  int32 scfInfoId;     // reference to a scfInfo in hypotheses generation

  afh_Info(afh_Info const &other);
  //afh_Info(SCF_look *r);
  afh_Info(SCF_look *r, const int32 dbPos, const int32 scfPos,
	   const char hypoBase, const char oldBase, const int32 fClass);
  afh_Info(const Contig::contigread_t *r, int32 dbPos, int32 scfPos,
	   char hypoBase, char oldBase, int32 fClass);
  afh_Info();
  ~afh_Info();
  afh_Info const & operator=(afh_Info const &other);
  afh_Info* nextAfh();
  int16 appendAfh(afh_Info *a);
  bool  loadRead();

  float eval(const bool verbose);

  SCF_look *performEditOnRead(afh_Info &a);
  SCF_look *getTraceData() { return theTrace; }

  bool  compareWith(afh_Info &a);
  int32 getFaultClass() const { return question; }
  bool  getIsExtendedAlterOperation();
  const char* getFaultClassName() const;
  const char* getFaultStatusName() const;
  const char* getReadName() const;
  const Contig::contigread_t* getContigRead() const { return contig_read; }
  int16 isSupportingHypotheses();
  bool  isNPlusHypotheses();
  bool  isNToBaseHypotheses();
  int16 isReversed();

  int32 scfPos() { return scfPosition; }
  int32 dbPos()  { return dbPosition;  }
  char  dbBase();
  char  getHypothesisBase()    { return hypothesis_base; }
  char  getOriginalBase()      { return original_base; }
  int16 getStatus() const   { return static_cast<int16>(status);   }
  int16 getQuestion() const { return static_cast<int16>(question); }
  void  confirm(const float confirmValue);
  void  confirm() { status = fhs_CONFIRMED; }
  void  reject(const float confirmValue);
  void  reject() { status = fhs_REJECTED; }
  bool  isConfirmed() const;

  base_quality_t getOriginalBaseQuality() { return original_quality;}
  void setOriginalBaseQuality(base_quality_t t) { original_quality = t; }

  int16 chemistryEvaluation(char before, char base, char after);

  friend ostream & operator<<(ostream &ostr, afh_Info const &i);
};



class cfh_Info {
private:
  afh_Info *first;
  afh_Info *last;

  int32 rejected_faults;
  int32 confirmed_positive;
  int32 confirmed_reverse;
  int32 rejected_n;
  int32 positive;
  int32 reversed;

  float rel_score;         // quality of confirmation

  void foolCompiler(void);
  void countInformation();
public:
  int32 cfh_id;

  cfh_Info(cfh_Info const &other);
  cfh_Info const & operator=(cfh_Info const &other);

  cfh_Info();
  cfh_Info(afh_Info *a);
  ~cfh_Info();

  afh_Info* firstAfh() const { return first; };
  bool  loadAfhTraces(afh_Info *a);
  float getRelScore() { return rel_score; };
  void  appendAfh(afh_Info *a);
  void  deleteAllBuffered();
  //int32 examineCfh();
  bool  examineCfhPartial(bool &isConfirmed);
  bool  needsInformationOfType(const int16 question, const int16 direction);
  bool  eval(const bool complete, const bool verbose,
	     const float threshold, const bool strict_N);
  bool  singleEval(const float threshold, int32 fhc_class = fhc_WRONG);
  bool  evalExtendedAlterOperations(const float threshold);
  bool  confirmQuality(const int32 threshold);
  afh_Info* readEditedBefore(afh_Info *afh_Info);

  friend ostream & operator<<(ostream &ostr, cfh_Info const &i);
};



int32 calculateParameters(int32 scfPos, int32 dbPos, int32 faultClass,
			char newBase, float *inParameter,
			SCF_look *aScfLook, ostream &output,
			bool &ko);

bool calculateThis(int32 calculation, int32 faultClass);


float evaluate(int32 scfPos, int32 dbPos, int32 faultClass,
	       char newBase, SCF_look *aScfLook, ostream &output);

float evaluate(int32 scfPos, int32 dbPos, int32 faultClass,
	       char newBase, SCF_look *aScfLook);


int32 bestPeakBetweenBases  (uint16 *samples, int32 von, int32 bis);
int32 bestValleyBetweenBases(uint16 *samples, int32 von, int32 bis);

bool  slopeChangeTest(const uint16 *samples, const int32 pos);

int32 constSlopeBetweenBases(const uint16 *samples, int32 von, int32 bis,
			     int32 &conv, int32 &div);

int32 findBestWeakPeak(const uint16 *samples, int32 von, int32 bis,
		       int32 &peak, int32 &valley);

int32 bestConvergenceRegion(uint16 *samples, int32 von, int32 bis);

#endif
