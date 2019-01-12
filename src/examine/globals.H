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
 *
 * Global definitions for the status of a hypotheses and for the
 * possible fault-hypothesis classes
 *
 * Written by Thomas Pfisterer
 */

#ifndef __EDIT_GLOBALS__
#define __EDIT_GLOBALS__

#define TEST
//#define RUN_ONLY

//#define DEBUG_EDIT(blah) { cout << blah; }
#define DEBUG_EDIT(blah) {}
//#define VERBOSE

//#define RANGECHECK_ON
//#define DEBUG_MODE
//#define PARAMETER_VERBOSE
//#define HYPOTHESES_VERBOSE
//#define FAULT_HYPOTHESES_VERBOSE


#define LIGHTHEADED

#define UNDEFINED_BASE '-'
#define NO_BASE        '.'
#define GAP_BASE       '*'
// end-of-read-character: character to fill fault regions with if
// the read has ended or not startet yet-
#define EOR_CHARACTER '_'


#define UNDEFINED_INDEX -1
#define WRONGCALL_INDEX  0
#define ADDITIONAL_INDEX 1
#define OVERCALL_INDEX   2
#define N_PLUS_INDEX     3
#define UNDERCALL_INDEX  4
#define MISSING_INDEX    5
#define N_CALL_INDEX     6

#define INDEX_COUNT 7



#define CALC_DISTANCE           1
#define CALC_PEAKS_OLD          2
#define CALC_PEAKS_NEW          4
#define CALC_CONV_OLD           8
#define CALC_CONV_NEW           16
#define CALC_CORRELATION_OLD    32
#define CALC_CORRELATION_NEW    64
#define CALC_POSRATING_NEW     128
#define CALC_POSRATING_OLD     256
#define CALC_PEAKRELATION      512
#define CALC_PEAK_VALLEY      1024
#define CALC_TRIPLETT_RANKING 2048
#define MARK_UNDEFINED_BASE   4096



#define DISTANCE_CODE    0x0000007fl   // 00000000 00000000 00000000 01111111
#define REL_DIST_PAR     0x00000001l
#define REL_MIN_PAR      0x00000002l
#define REL_MAX_PAR      0x00000004l
#define BASECOUNT_PAR    0x00000008l
#define MIN_DIST_PAR     0x00000010l
#define MAX_DIST_PAR     0x00000020l
#define MEAN_DIST_PAR    0x00000040l

#define PEAKS_OLD_CODE   0x00000180l   // 00000000 00000000 00000001 10000000
#define QUAL_OLD_PAR     0x00000080l
#define PEAKS_OLD_PAR    0x00000100l

#define PEAKS_NEW_CODE   0x00000600l   // 00000000 00000000 00000110 00000000
#define QUAL_NEW_PAR     0x00000200l
#define PEAKS_NEW_PAR    0x00000400l

#define CONV_OLD_CODE    0x00001800l   // 00000000 00000000 00011000 00000000
#define CONV_OLD_PAR     0x00000800l
#define AREA_OLD_PAR     0x00001000l

#define CONV_NEW_CODE    0x00006000l   // 00000000 00000000 01100000 00000000
#define CONV_NEW_PAR     0x00002000l
#define AREA_NEW_PAR     0x00004000l

#define POSRAT_NEW_CODE  0x00018000l   // 00000000 00000001 10000000 00000000
#define POS_NEW_PAR      0x00008000l
#define SHAPE_NEW_PAR    0x00010000l

#define POSRAT_OLD_CODE  0x00060000l   // 00000000 00000110 00000000 00000000
#define POS_OLD_PAR      0x00020000l
#define SHAPE_OLD_PAR    0x00040000l

#define PEAK_REL_CODE    0x00380000l   // 00000000 00111000 00000000 00000000
#define INT_NEW_PAR      0x00080000l
#define INT_OLD_PAR      0x00100000l
#define INT_OLDNEW_PAR   0x00200000l

#define PEAK_VALLEY_CODE 0x0fc00000l   // 00001111 11000000 00000000 00000000

#define COR_OLD_CODE     0x10000000l   // 00010000 00000000 00000000 00000000
#define COR_NEW_CODE     0x20000000l   // 00100000 00010000 00000000 00000000

#define MACHINETYPE_CODE 0x40000000l   // 01000000 00000000 00000000 00000000
#define PEAK1_REL_CODE   0x80000000l


//#define WRONGCALL_OPT POSRAT_NEW_CODE | POS_OLD_PAR | PEAK_REL_CODE | CORREL_OLD_CODE | CORREL_NEW_CODE

#define WRONGCALL_OPT POSRAT_OLD_CODE|POSRAT_NEW_CODE|COR_OLD_CODE|COR_NEW_CODE|PEAK_REL_CODE|MACHINETYPE_CODE

#define N_CALL_OPT POSRAT_OLD_CODE|POSRAT_NEW_CODE|COR_OLD_CODE|COR_NEW_CODE|PEAK_REL_CODE

// Wrongcall (Optimal: 1,2,3,5,6,7,8,9)
// CALC_POSRATING_NEW  (1,2)
// CALC_POSRATING_OLD  (3,4)
// PEAK_RELATION (5,6,7)
// COR_OLD + NEW (8,9)
// MARK_UNDEFINED (10)

//#define ADDITIONAL_OPT REL_MAX_PAR|BASECOUNT_PAR|MEAN_DIST_PAR|PEAKS_OLD_CODE
#define ADDITIONAL_OPT DISTANCE_CODE|PEAKS_OLD_CODE|PEAKS_OLD_CODE|CONV_OLD_CODE|POSRAT_OLD_CODE|MACHINETYPE_CODE

// Additional (Optimal: 3,4,7,8,9)
// CALC_DISTANCE (1-7)
// CALC_PEAKS_OLD  (8,9)
// CALC_CONV_OLD (10,11)
// CALC_POSRATING_OLD (12,13)

//#define OVERCALL_OPT REL_DIST_PAR|REL_MIN_PAR|BASECOUNT_PAR|MAX_DIST_PAR|PEAKS_OLD_CODE|PEAKS_NEW_CODE|CONV_OLD_CODE|POSRAT_OLD_CODE

#define OVERCALL_OPT DISTANCE_CODE|PEAKS_OLD_CODE|PEAKS_OLD_CODE|CONV_OLD_CODE|POSRAT_OLD_CODE|MACHINETYPE_CODE

// Overcall (Optimal: ohne 3 evtl. 5,7)
// CALC_DISTANCE  (1-7)
// CALC_PEAKS_OLD (8,9)
// CALC_CONV_OLD  (10,11)
// CALC_POSRATING_OLD (12,13)

//#define UNDERCALL_OPT DISTANCE_CODE|PEAKS_NEW_CODE|POSRAT_NEW_CODE
#define UNDERCALL_OPT DISTANCE_CODE|PEAKS_NEW_CODE|CONV_NEW_CODE|POSRAT_NEW_CODE|MACHINETYPE_CODE

// Undercall (Optimal =
// CALC_DISTANCE  (1-7)
// CALC_PEAKS_NEW (8,9)
// CALC_CONV_NEW  (10,11)
// CALC_POSRATING_NEW (12,13)

//#define MISSING_OPT DISTANCE_CODE|QUAL_OLD_PAR|SHAPE_NEW_PAR
#define MISSING_OPT DISTANCE_CODE|PEAKS_NEW_CODE|CONV_NEW_CODE|POSRAT_NEW_CODE|PEAK1_REL_CODE|MACHINETYPE_CODE

// Missing (Optimal z.B. 1-7, 8,13)
// CALC_DISTANCE  (1-7)
// CALC_PEAKS_NEW (8,9)
// CALC_CONV_NEW  (10,11)
// CALC_POSRATING_NEW (12,13)

//#define N_PLUS_OPT REL_DIST_PAR|REL_MIN_PAR|MAX_DIST_PAR|MEAN_DIST_PAR
#define N_PLUS_OPT DISTANCE_CODE | MACHINETYPE_CODE

// NPlus (Optimal 1,2,6,7)
// CALC_DISTANCE (1-7)




/* List of possible states of a fault hypotheses  */
/* Achtung: Behandelter Status > fhs_UNHANDLED!!  */
enum fhStatus {fhs_UNDEFINED, fhs_UNHANDLED, fhs_CONFIRMED,
	       fhs_REJECTED, fhs_UNDECIDABLE};

/* List of possible atomic fault hypotheses classes (uncomplete) */
enum fhClass  {fhc_UNDEFINED, fhc_OVERCALL, fhc_ADDITIONAL, fhc_WRONG,
	       fhc_UNDERCALL, fhc_MISSING,  fhc_SYNTACTIC,
	       fhc_CORRECT,   fhc_GAP};

enum fhRating {rating_GOOD, rating_INDIFFERENT,
	       rating_BAD, rating_VERY_BAD};

#define IS_BETTER_THAN <


/*
   fhc_GAP and fhc_CORRECT are not fault-hypotheses; they state that the
   situation found is right but this is also a hypotheses than can be
   examined
*/

/*
  #define HYPOTHESES_ONLY
*/


int readParameter(char *name);

#define STAT_INFO_FILE "EdIt.stat"
#define PARAMETER_FILE "EdIt.par"



// Definition of machine types

enum t_machinetype {
  MT_undefined,
  MT_unknown,
  MT_ABI_3700,
  MT_LICOR,
  MT_MegaBACE
};


#endif
