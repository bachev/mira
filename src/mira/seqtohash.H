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

#ifndef _seqtohash_h
#define _seqtohash_h


#include "stdinc/defines.H"

class seqtohash
{
public:
  static const uint8 hashaddmatrix[256];
};


/*************************************************************************
 *
 * To use these two macros, one needs to define
 *
 * const uint8 * seq
 * uint64        slen
 * const uint32  basesperhash
 * const char *  namestr
 *
 * In the loop, two variables may be of interest:
 *   VHASHT acthash   contains the sequence hash
 *   uint64 seqi      position of this hash + basesperhash
 *
 *************************************************************************/


#define SEQTOHASH_LOOPSTART(VHASHT) \
  BUGIFTHROW(basesperhash>sizeof(VHASHT)*4,"SEQTOHASH_LOOPSTART basesperhash " << static_cast<uint16>(basesperhash) << " > allowed size for VHASHT ?"); \
{						\
  VHASHT acthash(0);							\
  VHASHT hashmask(1);							\
  /* *grml* undefined behaviour of left shift for 64 shifts in a 64 bit type makes this cludge necessary */ \
  /* the same for 32 shift in 32 bit types etc.pp */			\
  if(basesperhash==sizeof(VHASHT)*4){					\
    hashmask=0;								\
  }else{								\
    hashmask<<=(basesperhash*2);					\
  }									\
  --hashmask;								\
									\
  uint32  baseok=0;							\
									\
  for(uint64 seqi=0; likely(seqi<slen); ++seqi, ++seq){			\
    acthash<<=2;							\
    acthash&=hashmask;							\
    ++baseok;								\
									\
    if(likely(::seqtohash::hashaddmatrix[*seq])){			\
      acthash|=::seqtohash::hashaddmatrix[*seq]-1;			\
    }else{								\
      if(likely(dptools::isValidIUPACStarBase(*seq))) {			\
        /* the IUPAC bases are treated like N and X */			\
	/* break hash making (which is actually better than behaving */	\
	/*  like another character in case of multiple bases with */	\
	/*  IUPAC or '*') */						\
	acthash=0;							\
	baseok=0;							\
      } else {								\
	cout << "Illegal base '" << *seq << "' (ASCII " << static_cast<uint16>(*seq) << ") at position " << seqi << " in sequence "; \
	if(namestr!=nullptr) {						\
	  cout << namestr << endl;					\
	}else{								\
	  cout << "(no name given)" << endl;				\
	}								\
	exit(100);							\
      }									\
    }									\
									\
    if(likely(baseok >= basesperhash)) {


#define SEQTOHASH_LOOPEND }}}


#endif
