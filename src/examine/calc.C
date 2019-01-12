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


#include "examine/scf_look.H"
#include "examine/globals.H"

extern "C" {
  int add(float *in, float *out, int init);
  int missing(float *in, float *out, int init);
  int nplus(float *in, float *out, int init);
  int over(float *in, float *out, int init);
  int under(float *in, float *out, int init);
  int wrong(float *in, float *out, int init);
  int ncall(float *in, float *out, int init);
}

float evaluate(int32 scfPos, int32 dbPos, int32 faultClass,
	       char newBase, SCF_look *aScfLook)
{
  ofstream devNull;

  devNull.open("/dev/null", ios::out);

  return evaluate(scfPos, dbPos, faultClass, newBase, aScfLook, devNull);
}




float evaluate(int32 scfPos, int32 dbPos, int32 faultClass,
	       char newBase, SCF_look *aScfLook, ostream &output)
{
  float in_parameter[25];
  float out_parameter[1];
  int32 parameter_anz;
  int32 nnFaultClass;
  char  oldDBBase;
  bool  ko;


  if (!aScfLook->testDBPos(dbPos)) {
      // BACh: added the SILENT throw: if no filename, no point in warning
      //  (shouldn't this not be called anyway?)
      if(aScfLook->getFileName() != nullptr && strlen(aScfLook->getFileName())) {
	throw Notify(Notify::WARNING,
		     "evaluate: DB-Baseposition exceeds limit. Please check file: ",
		     aScfLook->getFileName());
      } else {
	throw Notify(Notify::SILENT, "");
      }
  }

  if (faultClass == fhc_SYNTACTIC) {
    if (newBase == '*') oldDBBase = '.';
    if (newBase == '.') oldDBBase = '*';

    output << "Evaluate " << aScfLook->getFileName() << " at: " << dbPos;
    aScfLook->showDBArea(dbPos, 6, output);
    return 1.0;

  } else {
    if (faultClass == fhc_UNDERCALL || faultClass == fhc_MISSING) {
      oldDBBase = '*';
    } else {
      oldDBBase = aScfLook->getDBBase(dbPos);
    }
  }


  output << "\nEvaluate " << aScfLook->getFileName() << " at: " << dbPos;
  output << " (" << scfPos << ")  " << oldDBBase
	 << " => " << newBase << "  faultClass " << faultClass << endl;

  if (dbPos < 0) {
    output << "DBPos < 0 found - abort examination" << endl;
    return 0.0;
  }

  if (dbPos >= aScfLook->getDBBases()) {
    output << "DBPos > readLen - abort examination" << endl;
    return 0.0;
  }


  /*
  if (scfPos < 0) {
    oldDBBase = '*';
  } else {
    //    oldDBBase  = aScfLook->getOriginalSCFBase(scfPos);
    oldDBBase = aScfLook->getDBBase(dbPos);
  }
  */



  if (scfPos > (int32)aScfLook->getNumBases()) {
    output << "Error in SCF-Position " << scfPos << " > "
	   << aScfLook->getNumBases() << endl;
    return 0.0;
  }

  nnFaultClass = faultClassToNNClass(faultClass, oldDBBase);

  if (nnFaultClass != UNDEFINED_INDEX) {
    try {
      parameter_anz = calculateParameters(scfPos, dbPos, nnFaultClass,
					  newBase, in_parameter,
					  aScfLook, output, ko);

      for (int32 loop = 0; loop < parameter_anz; loop++) {
	output << in_parameter[loop] << "  ";
      }
      output << endl;
    }

    catch(Notify n) {
      output << "Problems while calculating parameters.....\n" << endl;
      output << "Abort and reject evaluating. \n" << endl;
      return 0.0;
    }

    switch (nnFaultClass) {
    case WRONGCALL_INDEX:
      wrong(in_parameter, out_parameter, 0);
      break;
    case UNDERCALL_INDEX:
      under(in_parameter, out_parameter, 0);
      break;
    case MISSING_INDEX:
      missing(in_parameter, out_parameter, 0);
      break;
    case OVERCALL_INDEX:
      over(in_parameter, out_parameter, 0);
      break;
    case N_PLUS_INDEX:
      nplus(in_parameter, out_parameter, 0);
      break;
    case ADDITIONAL_INDEX:
      add(in_parameter, out_parameter, 0);
      break;
    case N_CALL_INDEX:
      ncall(in_parameter, out_parameter, 0);
      break;
    default: return 0.0;
    }

    output << "Result: " << out_parameter[0] << endl;

    if (ko) {
      output << "Apply k.o. rule" << endl;
      out_parameter[0] = 0;
    }

    return out_parameter[0];
  }
  return 1.0;
}



// Achtung!!!
// faultClass ist hier nicht die fhc klasse, sondern die hier
// verwendete Fehlerklasse!!
int32 calculateParameters(int32 scfPos, int32 dbPos, int32 faultClass,
			char newBase, float *inParameter,
			SCF_look *aScfLook, ostream &output, bool &ko)
{
 int32 areaMin, areaMax;
 int16 anzPeaks, anzRepeat, qual;
 int32 regCount;
 char  oldDBBase, oldSCFBase;
 int32 parameter_used = 0;
 int32 pos[4];

 int32 calcCode = 0;
 float meanVic_meanRep = 1.0;

 ko = false;

 // if (scfPos<0|| faultClass==UNDERCALL_INDEX || faultClass==MISSING_INDEX) {
 //  oldDBBase = GAP_BASE;
 //} else {
 oldDBBase  = toupper(aScfLook->getDBBase(dbPos));
 //}
 oldSCFBase = toupper(toSCFBase(oldDBBase, aScfLook->isReversed()));



 if (newBase == oldDBBase && faultClass == WRONGCALL_INDEX) {
   output << "Edited back to old Base!!!!" << endl;
   //return 0;
 }


 output << "\nOldBase was : " << oldDBBase << " (DB) " << oldSCFBase
	<< " (SCF)   newBase: " << newBase
	<< "  nnfaultClass " << faultClass
	<< "  scfPos: " << scfPos
        << "  dbPos: " << dbPos << endl;

 aScfLook->showDBArea(dbPos, 7, output);
 output << endl;


 switch (faultClass) {
 case WRONGCALL_INDEX:
   calcCode = WRONGCALL_OPT;
   break;
 case N_CALL_INDEX:
   calcCode = N_CALL_OPT;
   break;
 case UNDERCALL_INDEX:
   calcCode = UNDERCALL_OPT;
   break;
 case MISSING_INDEX:
   calcCode = MISSING_OPT;
   break;
 case OVERCALL_INDEX:
   calcCode = OVERCALL_OPT;
   break;
 case N_PLUS_INDEX:
   calcCode = N_PLUS_OPT;
   break;
 case ADDITIONAL_INDEX:
   calcCode = ADDITIONAL_OPT;
   break;
 default: return 0;
 }


 {
   int32 max[4];
   aScfLook->findMaxDB(dbPos, max);
   indirectSort4(max, pos);
 }


 if (DISTANCE_CODE & calcCode) {

   // --------------------------------------------------------------------
   // Calculates the following four parameters:
   // meanVic/meanRep:    distance in ROI relative to distances in vicinity
   // minRep/100:         smallest distance in ROI in relation to the mean
   //                     distance in the vicinity
   // maxRep/100:         same for maximum distance
   // repBases:           number of identical bases in the ROI
   //
   // Parameters:
   // the position of the ROI in the database
   // ---------------------------------------------------------------------

   output << "calculate distance parameter" << endl;


   SCF_distance *pd=nullptr;


   try {
     pd = aScfLook->calcDBDistance(dbPos);

     // BaCh: changed abs() to fabs() (seems what Thoams intended)
     if (fabs(pd->getRelativeProblemDistance()) < 0.01) {
       meanVic_meanRep = 1.0;
     } else {
       meanVic_meanRep = 1.0 / pd->getRelativeProblemDistance();
     }

     if (REL_DIST_PAR & calcCode) {
       inParameter[parameter_used++] = meanVic_meanRep;
     }


     if (REL_MIN_PAR & calcCode)
       inParameter[parameter_used++] = pd->getRelativeMinDistance();
     if (REL_MAX_PAR & calcCode)
       inParameter[parameter_used++] = pd->getRelativeMaxDistance();

     if (BASECOUNT_PAR & calcCode)
       inParameter[parameter_used++] = 1.0 / pd->getProblemBases();

     if (MIN_DIST_PAR & calcCode)
       inParameter[parameter_used++] = pd->getMinDistanceScore();
     if (MAX_DIST_PAR & calcCode)
       inParameter[parameter_used++] = pd->getMaxDistanceScore();
     if (MEAN_DIST_PAR & calcCode)
       inParameter[parameter_used++] = pd->getDistanceScore();

   }
   catch(Notify n) {
     output << "FEHLER: calcDBDistance " << dbPos << endl;

     if (REL_DIST_PAR & calcCode)
       inParameter[parameter_used++] = 0;
     if (REL_MIN_PAR & calcCode)
       inParameter[parameter_used++] = 0;
     if (REL_MAX_PAR & calcCode)
       inParameter[parameter_used++] = 0;

     if (BASECOUNT_PAR & calcCode)
       inParameter[parameter_used++] = 0;

     if (MIN_DIST_PAR & calcCode)
       inParameter[parameter_used++] = 0;
     if (MAX_DIST_PAR & calcCode)
       inParameter[parameter_used++] = 0;
     if (MEAN_DIST_PAR & calcCode)
       inParameter[parameter_used++] = 0;
   }

   if(pd!=nullptr) delete pd;

 }



 if (PEAKS_OLD_CODE & calcCode) {

   assert(isRealBase(oldDBBase));

   output << "calculate peaks_old " << oldDBBase << endl;
   aScfLook->countRepeatPeaks(dbPos, oldDBBase, anzRepeat, anzPeaks, qual);

   output << "Peaks : " << anzPeaks << "  Repeat: " << anzRepeat << endl;

   if (QUAL_OLD_PAR & calcCode)
     inParameter[parameter_used++] = float(qual) / 100;
   if (PEAKS_OLD_PAR & calcCode)
     inParameter[parameter_used++] = anzRepeat - anzPeaks;


   if (faultClass == OVERCALL_INDEX) {
     if (anzRepeat >= anzPeaks && meanVic_meanRep < 1.12) {
       output << "KO-Rule: anzRepeat == anzPeaks && meanVic_meanRep < 1.05\n";
       output << "meanVic_meanRep = " << meanVic_meanRep << endl;

       ko = true;
     }
   }
   if (faultClass == ADDITIONAL_INDEX) {
     if (anzPeaks == 0 && meanVic_meanRep < 1.2) {
       output << "KO-Rule: anzPeaks=0 && meanVic_meanRep < 1.10" << endl;
       ko = true;
     }
   }
 }



 if (PEAKS_NEW_CODE & calcCode) {
   assert(isRealBase(newBase));

   output << "calculate peaks for new " << newBase << endl;
   aScfLook->countRepeatPeaks(dbPos, newBase, anzRepeat, anzPeaks, qual);

   output << "Peaks : " << anzPeaks << "  Repeat: " << anzRepeat << endl;

   if (QUAL_NEW_PAR & calcCode)
     inParameter[parameter_used++] = float(qual)/100;
   if (PEAKS_NEW_PAR & calcCode)
     inParameter[parameter_used++] = anzRepeat - anzPeaks;


   if (faultClass == UNDERCALL_INDEX) {
     if (anzRepeat >= anzPeaks && meanVic_meanRep > .80) {
       output << "KO-Rule: anzRepeat == anzPeaks && meanVic_meanRep > 0.95\n";
       output << "meanVic_meanRep = " << meanVic_meanRep << endl;

       ko = true;
     }
   }
   if (faultClass == MISSING_INDEX) {
     if ((anzPeaks <= 0 && meanVic_meanRep > 0.75) || meanVic_meanRep > 0.90) {
       output << "KO-Rule: anzPeaks=0 && meanVic_meanRep> 0.67" << endl;
       ko = true;
     }
   }
 }


 if (CONV_OLD_CODE & calcCode) {
   // regCount        number of convergence regions (f"<0) - number of Bases
   // areaMin         Area under the smallest convergence region.
   // areaMax         Area under the largest convergence region.
   // These parameters make only sense in case of repeats (undercall or
   // overcall).

   output << "Calculate conv.regions for original base " << oldDBBase << endl;

   regCount = aScfLook->countConvergenceRegions(dbPos, oldDBBase,
						areaMin, areaMax);

   if (CONV_OLD_PAR & calcCode) {
     if (regCount != -1) {
       inParameter[parameter_used++] = 1/float(regCount+1);
     } else {
       inParameter[parameter_used++] = 1;
     }
   }

   if (AREA_OLD_PAR & calcCode) {
     if (areaMax < 0.01) {
       inParameter[parameter_used++] = 0.50;
     } else {
       inParameter[parameter_used++] = float(areaMin) / float(areaMax);
     }
   }
 }


 if (CONV_NEW_CODE & calcCode) {
   int32 regCount;

   output << "Calculate conv.regions for new base " << newBase << endl;

   regCount = aScfLook->countConvergenceRegions(dbPos, newBase,
						areaMin, areaMax);

   if (CONV_NEW_PAR & calcCode) {
     if (regCount != -1) {
       inParameter[parameter_used++] = 1/float(regCount+1);
     } else {
       inParameter[parameter_used++] = 1;
     }
   }

   if (AREA_NEW_PAR & calcCode) {
     if (areaMax < 0.01) {
       inParameter[parameter_used++] = 0.50;
     } else {
       inParameter[parameter_used++] = float(areaMin) / float(areaMax);
     }
   }
 }


 if (POSRAT_NEW_CODE & calcCode) {
   float new_pos_rating = 0;
   float hoehenkrit = 0;
   int32 max;


   if (oldDBBase == '*' || oldDBBase == '.') {
     output << "calculate posrating new (1) " << endl;

     hoehenkrit = aScfLook->einzelpeak(-dbPos, newBase, max, new_pos_rating);
   } else {
     output << "calculate posrating new (2) " << endl;
     hoehenkrit = aScfLook->einzelpeak(dbPos, newBase, max, new_pos_rating);
   }

   if (POS_NEW_PAR & calcCode)
     inParameter[parameter_used++] = new_pos_rating;
   if (SHAPE_NEW_PAR & calcCode)
     inParameter[parameter_used++] = hoehenkrit;

 }



 if (POSRAT_OLD_CODE & calcCode) {
   float old_pos_rating = 0;
   float hoehenkrit = 0;
   int32 max;


   if (newBase == '*' || newBase == '.') {
     output << "calculate posrating old (1)" << oldDBBase << endl;
     hoehenkrit = aScfLook->einzelpeak(-dbPos, oldDBBase, max,old_pos_rating);
   } else {
     output << "calculate posrating old (2)" << endl;
     hoehenkrit = aScfLook->einzelpeak(dbPos, oldDBBase, max,old_pos_rating);
   }

   if (POS_OLD_PAR & calcCode)
     inParameter[parameter_used++] = old_pos_rating;
   if (SHAPE_OLD_PAR & calcCode)
     inParameter[parameter_used++] = hoehenkrit;

 }


 if (PEAK1_REL_CODE & calcCode) {
   float newRatio, oldRatio, oldNewRatio;

   output << "calculate single peak relations" << endl;

   assert(isRealBase(newBase));

   aScfLook->peakRelations(dbPos, newBase, newBase,
			   newRatio, oldRatio, oldNewRatio);

   inParameter[parameter_used++] = newRatio;
 }



 if (PEAK_REL_CODE & calcCode) {
   float newRatio, oldRatio, oldNewRatio;

   output << "calculate peak relations" << endl;

   assert(isRealBase(newBase));

   aScfLook->peakRelations(dbPos, newBase, oldDBBase,
			   newRatio, oldRatio, oldNewRatio);

   if (INT_NEW_PAR & calcCode)
     inParameter[parameter_used++] = newRatio;
   if (INT_OLD_PAR & calcCode)
     inParameter[parameter_used++] = oldRatio;
   if (INT_OLDNEW_PAR & calcCode)
     inParameter[parameter_used++] = oldNewRatio;

   if (faultClass == WRONGCALL_INDEX || faultClass == N_CALL_INDEX){
     if (isRealBase(oldDBBase)) {
       if (newRatio < 0.12) {
	 output << "KO-Rule: newRatio < 0,12 " << newRatio << endl;
	 ko=true;
       }
       if (oldNewRatio < 0.12) {
	 output << "KO-rule: oldNewRatio < 0,12 " << oldNewRatio << endl;
	 ko =true;
       }
       if (newRatio + oldNewRatio < 0.35) {
	 output << "KO-rule: newRatio + oldNewRatio < 0,35" << endl;
	 output << "newRatio:"<<newRatio<<" oldNewRatio:"<<oldNewRatio<<endl;
	 ko = true;
       }
     } else {
       if ((newRatio < 0.10) && (oldNewRatio < 0.50)) {
	 output << "KO-rule: newRatio < 0.10 && oldNewRatio < 0.50 " << endl;
	 output << "newRatio: " << newRatio
		<< " oldNewRatio: " << oldNewRatio << endl;
	 ko = true;
       }
     }
   }
 }


 if (PEAK_VALLEY_CODE & calcCode) {
   float result[4][6];
   int32 oldPos = 0;
   int32 newPos = 1;

   output << "calculate peak-valley" << endl;

   if (isRealBase(oldDBBase)) {
     DEBUG_EDIT("Normal Case " << oldDBBase << " - " << newBase << endl);
     aScfLook->peakValleyValues(dbPos, newBase, result[oldPos]);
     aScfLook->peakValleyValues(dbPos, oldDBBase, result[newPos]);
   } else {
     if (newBase == GAP_BASE) {
       DEBUG_EDIT("N-Case (nplus) " << endl);

       aScfLook->peakValleyValues(dbPos, "ACGT"[pos[0]], result[oldPos]);
       aScfLook->peakValleyValues(dbPos, "ACGT"[pos[1]], result[newPos]);

     } else {
       DEBUG_EDIT("N-Case (alter) " << endl);
       float x;
       float max = 0;


       for (int32 i=0; i<4; i++) {
	 x = aScfLook->peakValleyValues(dbPos, "ACGT"[i], result[i]);

	 if ("ACGT"[i] == newBase) {
	   oldPos = i;
	 } else {
	   if (x > max) {
	     max = x;
	     newPos = i;
	   }
	 }
       }
     }
   }

   for (int32 i=0; i<6; i++) {
     inParameter[parameter_used++] = result[oldPos][i];
   }
 }


 if (COR_OLD_CODE & calcCode) {
   double sigma = 5.0;
   double old_value;

   output << "calculate correlation old base" << endl;

   if (!isBase(oldDBBase)) {
     output << "WARNING: expected a base !! " << endl;
     inParameter[parameter_used++] = 0.0;
     ko = true;
   } else {
     old_value = aScfLook->singlePeakQuality(dbPos, sigma, oldDBBase);
     inParameter[parameter_used++] = old_value;
   }
 }


 if (COR_NEW_CODE & calcCode) {
   double sigma = 5.0;
   double new_value;

   output << "calculate correlation new base" << endl;
   assert(isBase(newBase));
   new_value = aScfLook->singlePeakQuality(dbPos, sigma, newBase);
   inParameter[parameter_used++] = new_value;

   output << "New Base " << newBase << "\tvalue: " << new_value << endl;
 }

 /*
 if (calculateThis(CALC_TRIPLETT_RANKING, faultClass)) {
    prescan  aPrescan;
    int32    index;
    int32    ranking;

    output << "calculate triplett ranking" << endl;
    index = baseToIndex3(aScfLook->getOriginalSCFBase(scfPos-2),
			 aScfLook->getOriginalSCFBase(scfPos-1),
			 aScfLook->getOriginalSCFBase(scfPos));

    ranking = aPrescan.getRanking(index);

    output << aScfLook->getOriginalSCFBase(scfPos-2)
	   << aScfLook->getOriginalSCFBase(scfPos-1)
	   << aScfLook->getOriginalSCFBase(scfPos)
	   << "  index:  "  << index
	   << "  ranking: " << ranking << endl;

    inParameter[parameter_used++] = (float)ranking / 64.0;
 }
 */


 if (MACHINETYPE_CODE & calcCode) {
   int32 m = (int32)aScfLook->getMACHType();

   inParameter[parameter_used++] = (m == 1);
   inParameter[parameter_used++] = (m == 2);
   inParameter[parameter_used++] = (m == 3);

 }

 return parameter_used;

}




bool calculateThis(int32 calculation, int32 faultClass)
{
  int32 calculate[INDEX_COUNT];

  calculate[WRONGCALL_INDEX ] =
    CALC_POSRATING_OLD
    | CALC_POSRATING_NEW
    | CALC_PEAKRELATION
    //| CALC_PEAK_VALLEY
    //| MARK_UNDEFINED_BASE
    //| CALC_TRIPLETT_RANKING
    | CALC_CORRELATION_OLD
    | CALC_CORRELATION_NEW
    ;

  calculate[N_CALL_INDEX ] =
    CALC_POSRATING_OLD
    | CALC_POSRATING_NEW
    | CALC_PEAKRELATION
    //| CALC_PEAK_VALLEY
    //| MARK_UNDEFINED_BASE
    //| CALC_TRIPLETT_RANKING
    | CALC_CORRELATION_OLD
    | CALC_CORRELATION_NEW
    ;
  calculate[ADDITIONAL_INDEX] =
    CALC_DISTANCE
    | CALC_PEAKS_OLD
    | CALC_CONV_OLD
    | CALC_POSRATING_OLD
    // | CALC_PEAK_VALLEY
    // | CALC_POSRATING_NEW
    ;

  calculate[OVERCALL_INDEX  ] =
    CALC_DISTANCE
    | CALC_PEAKS_OLD
    | CALC_CONV_OLD
    | CALC_POSRATING_OLD
    ;

  calculate[N_PLUS_INDEX    ] =
    CALC_DISTANCE
    //| CALC_PEAKRELATION
    //| CALC_PEAK_VALLEY
    //| CALC_CONV_OLD
    ;

  calculate[UNDERCALL_INDEX ] =
    CALC_DISTANCE
    | CALC_PEAKS_NEW
    | CALC_CONV_NEW
    | CALC_POSRATING_NEW
    ;

  calculate[MISSING_INDEX   ] =
    CALC_DISTANCE
    | CALC_PEAKS_NEW
    | CALC_CONV_NEW
    | CALC_POSRATING_NEW
    ;


  return (calculate[faultClass] & calculation);
}
