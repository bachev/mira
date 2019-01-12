/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2005 and later by Bastien Chevreux
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
 */


#ifndef _util_dptools_h
#define _util_dptools_h

#include <cstring>
#include <ctime>

#include <iostream>

#include <string>
#include <vector>

#include "stdinc/defines.H"

#include "stdinc/stlincludes.H"


class dptools {
private:
  static char nsTranslationTables[23][321];
  static char nsvalid_ACGTbases[256];
  static char nsvalid_bases[256];
  static char nscomplement_bases[256];
  static char nsvalidIUPAC_bases[256];
  static char nscomplementIUPAC_bases[256];
  static uint8 nsIUPAC_basebitmask[256];
  static uint8 nsIUPAC_basebitmaskN[256];
  static char nsIUPACbitmask[16];
  static uint8 nsbase2index_translationtable[256];

public:

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  inline static uint8 getIndexOfBase(const char c){
    return nsbase2index_translationtable[static_cast<uint8>(c)];
  }

  static uint8 getIndexOfCodon(const char b1, const char b2, const char b3);
  inline static uint8 getIndexOfCodon(const char * ptr) {return getIndexOfCodon(*ptr,*(ptr+1),*(ptr+2));}
  inline static uint8 getIndexOfCodon(const std::string & s) {return getIndexOfCodon(s.c_str());}

  inline static bool isValidACGTBase(const char c){
    return nsvalid_ACGTbases[static_cast<uint8>(c)];
  }

  inline static bool isValidACGTNBase(const char c){
    return nsvalid_ACGTbases[static_cast<uint8>(c)] || c=='n' || c=='N';
  }

  inline static bool isValidACGTStarBase(const char c){
    return nsvalid_ACGTbases[static_cast<uint8>(c)] || c=='*';
  }

  inline static bool isValidBase(const char c){
    return nsvalid_bases[static_cast<uint8>(c)];
  }

  inline static bool isValidIUPACBase(const char c){
    return nsvalidIUPAC_bases[static_cast<uint8>(c)];
  }

  inline static bool isValidStarBase(const char c){
    return nsvalid_bases[static_cast<uint8>(c)] || c=='*';
  }

  inline static bool isValidIUPACStarBase(const char c){
    return nsvalidIUPAC_bases[static_cast<uint8>(c)] || c=='*';
  }

  inline static char getComplementBase(const char c){
    return nscomplement_bases[static_cast<uint8>(c)];
  }

  inline static char getComplementIUPACBase(const char c){
    return nscomplementIUPAC_bases[static_cast<uint8>(c)];
  }

 /*************************************************************************
 *
 * both bases being a valid acgt or othe IUPAC code, check whether
 *  one base is contained in the other.
 * e.g.  a in a -> true
 *       a in c -> false
 *       a in m -> true  (a can be in base1 or base2, as can be m)
 *       m in r -> true  (they share the a)
 *
 * Exception: the 'n' base is NOT a catch all:
 * e.g.  a in n -> false
 *
 *************************************************************************/

//  inline static bool areBasesContained(const char base1, const char base2)
//  {
//    if(nsIUPAC_basebitmask[static_cast<uint8>(base1)]
//       && nsIUPAC_basebitmask[static_cast<uint8>(base2)]
//       && ((nsIUPAC_basebitmask[static_cast<uint8>(base1)] | nsIUPAC_basebitmask[static_cast<uint8>(base2)]) == nsIUPAC_basebitmask[static_cast<uint8>(base1)]
//	   || (nsIUPAC_basebitmask[static_cast<uint8>(base1)] | nsIUPAC_basebitmask[static_cast<uint8>(base2)]) == nsIUPAC_basebitmask[static_cast<uint8>(base2)]))
//      return true;
//
//    return false;
//  }

  inline static bool areBasesContained(const char base1, const char base2)
  {
    return (nsIUPAC_basebitmask[static_cast<uint8>(base1)]
	    & nsIUPAC_basebitmask[static_cast<uint8>(base2)]);
  }


 /*************************************************************************
 *
 * like the version without N, but there "is a in n" returns true
 *
 *************************************************************************/

  inline static bool areBasesContainedWithN(const char base1, const char base2)
  {
    return (nsIUPAC_basebitmaskN[static_cast<uint8>(base1)]
	    & nsIUPAC_basebitmaskN[static_cast<uint8>(base2)]);
  }



 /*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  inline static char calcIUPACConsensus(const char a, const char b)
  {
    return nsIUPACbitmask[nsIUPAC_basebitmask[static_cast<uint8>(a)]
  			  |nsIUPAC_basebitmask[static_cast<uint8>(b)]];
  }


 /*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
  inline static char calcIUPACConsensus(const std::vector<char> &v)
  {
    uint8 bm=0;

    for(auto & ve : v) {
      // not using bm|= ... as gcc4.3.2 somehow gives spurious warnings
      bm=static_cast<uint8>(bm|nsIUPAC_basebitmask[static_cast<uint8>(ve)]);
    }
    return nsIUPACbitmask[bm];
  }

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  inline static char calcIUPACConsensus(const std::string &v)
  {
    uint8 bm=0;

    for(auto & ve : v) {
      // not using bm|= ... as gcc4.3.2 somehow gives spurious warnings
      bm=static_cast<uint8>(bm|nsIUPAC_basebitmask[static_cast<uint8>(ve)]);
    }
    return nsIUPACbitmask[bm];
  }

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  inline static std::vector<char> getNucleicAcidFromIUPAC(const char base)
  {
    uint8 bm=nsIUPAC_basebitmaskN[static_cast<uint8>(base)];

    std::vector<char> v;
    if(bm & 8) v.push_back('A');
    if(bm & 4) v.push_back('C');
    if(bm & 2) v.push_back('G');
    if(bm & 1) v.push_back('T');
    return v;
  }

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  static bool hasNucleicAcidInIUPAC(char nucleicacid, const char base)
  {
    uint8 bm=nsIUPAC_basebitmaskN[static_cast<uint8>(base)];
    nucleicacid=static_cast<char>(toupper(nucleicacid));

    switch(nucleicacid){
    case 'A': {
      if(bm & 8) return true;
      break;
    }
    case 'C': {
      if(bm & 4) return true;
      break;
    }
    case 'G': {
      if(bm & 2) return true;
      break;
    }
    case 'T': {
      if(bm & 1) return true;
      break;
    }
    default: {}
    }

    return false;
  }


  static void codon2AminoAcids_wrapped(const uint8 transl_table,
				       char base1,
				       char base2,
				       char base3,
				       std::vector<std::string> & codonvariants,
				       std::vector<char> & aaresult,
				       std::vector<bool> & isstartresult);

  static void codon2AminoAcids(const uint8 transl_table,
			       const char base1,
			       const char base2,
			       const char base3,
			       std::vector<std::string> & codonvariants,
			       std::vector<char> & aaresult,
			       std::vector<bool> & isstartresult);

  static void dnaToProtein(const std::string & dna,
			   std::string & protresult,
			   std::string & dnaresult,
			   uint32 orfstart,
			   uint32 orfend,
			   const int8 orfdirection,
			   const uint8 transl_table,
			   const uint8 codonstart,
			   bool musttranslate);

  static void infoOnAAatDNAPos(const std::string & dna,
			       const uint32 whichdnapos,
			       const uint32 orfstart,
			       const int8 orfdirection,
			       const uint8 transl_table,
			       const uint8 codonstart,
			       std::string & codon,
			       std::vector<std::string> & codonvariants,
			       std::vector<char> & aminoacid,
			       std::vector<bool> & isstart,
			       int32 & aanumber,
			       int8 & posinaa);

  static bool isCodonStart(const uint8 transl_table,
			   const char base1,
			   const char base2,
			   const char base3);

  static double calcProteinIdentity(const std::string & p1,
				    const std::string & p2);
};

#endif
