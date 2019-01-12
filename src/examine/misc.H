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
 * misc functions - header file
 *
 * Written by Thomas Pfisterer
 * 09.10.1999
 */


#ifndef _edit_misc_H
#define _edit_misc_H


#include "scf_look.H"


//struct Contig::consensus_counts_t;


// get the fault class from what we want to change...
int32 findFaultClass(char prev, char aBase, char next, char newBase);


// Find the fault class (and neural net) for the given fhc fault-class
// and the old base in the DB. Differences are e.g. altering an
// undefined base which has a seperate NNClass.
int32 faultClassToNNClass(const int32 faultClass, const char oldDBBase);


// find the consensus base
char calledBase(const Contig::consensus_counts_t &cc);


// Sorts array "pos" such that references to "a" are sorted (descending)
// Is used to find the n-th-greatest value in array "a"
void  indirectSort4(int32 a[4], int32 pos[4]);



int16 isRealBase(const char c);
int16 isBase(const char c);
int16 isUndefinedBase(const char c);
char  invertBase(const char c);

int32 valueForBase(const int32 theArray[4], const char theBase);
char  toSCFBase(const char c, const int32 reversed);
char  toDBBase(const char c, const int32 reversed);

double simpleNormedCorrelation(uint16 *sequence, int32 vonpos,
			      int32 bispos, double sigma);
double autoCorrelation(uint16 *sequence, int32 vonpos, int32 bispos, int32 o);

#endif
