/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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


#ifndef _bas_vluint_h_
#define _bas_vluint_h_


/*
  This is a no-frills class for "large" integers as needed by MIRA, i.e., only
  a tiny subset of functions has been implemented.

  I had looked at
  1. https://github.com/calccrypto/uint256_t
  2. https://github.com/mjtryhorn/Boost.LargeInt
  3. http://src.chromium.org/svn/trunk/src/net/base/int128.h
  4. http://www.codef00.com/coding/uint128.h
  5. https://github.com/bitcoin/bitcoin/blob/master/src/uint256.cpp and *.h
  6. http://www.boost.org/doc/libs/1_55_0/libs/multiprecision/doc/html/index.html
  and the discussion at
    https://groups.google.com/forum/#!topic/boost-devel-archive/Maax-DSiE1M

  Mode of action:
  (1) uses a recursive mode (2*64 bit to define a 128, 2*128 to define a 256)
  (2) is a generalisation of (1), probably to "infinity"
  (3) & (4) use 2*64 bits
  (5) Uses an array of 8*32 bits
  (6) Arbitrary precision

  Drawbacks:
  (1) can do only up to 256 bits (and is "slow"); (2) is "slow"; (3) only 128 bits;
  (4) only 128 bits; (5) only 256 bits (and "slow" as based on 32 bits)
  (6) uses additional pointers (needs more memory)

  I decided for an own class which can do abitrary sizes without additional
  memory needed and relying on 64 bit operations. (5) was used as general
  storage setup, combined with 64 bit storage and use of standard C++ library
  calls for memory fills and comparisons where possible (they can be highly
  optimised and better than hand-rolled functions).

  TODO: I'm not sure the memory layout is the best possible. (5) uses a
  low-to-high scheme, i.e., the uint64 with the lowest bits is in payload[0]
  and the highest bits are in payload[n]. This is bad for comparison operators
  which need to go backwards through the array, leading to memory faults at
  every step! Doing it the other way round could be faster due to automatic
  cache prefetching.

  It would make sense for adding and subtracting (which are generally used
  more often than shifting), but as MIRA does not use these (and they're not
  implemented), turning the memory organisation around could make sense. TBD /
  tested when I have time.

  REMARKS:

  1) I am in awe. Somehow, a comparator involving operator&(uint64 &) is >50%
  slower than a comparator using operator&(Self &) when used in a std::sort()
  (e.g. mirabait -b bs168data ...) jumping from 15s to 22s runtime. I suspect
  the optimiser is able to change that into using the operator&= which does
  not involve creating a temporary copy. But still ...

  2) I am depressed. Using VLuint<64> instead of uint64_t leads to mirabait
  doubling in runtime. This does not feel right.

  3) I am scared because of the inexplicable segfault, see priv_compareLess()
  below.
 */


//#include <iostream>

#include "stdinc/defines.H"

// TODO: having extra namespace??? can we put that to VLuint without
//  using the same space every time???
namespace nsVLuint {
  // bit reverse table
  static const uint8 brtable[] =
  {
    0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0, 0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0, 0x70, 0xf0,
    0x08, 0x88, 0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8, 0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8,
    0x04, 0x84, 0x44, 0xc4, 0x24, 0xa4, 0x64, 0xe4, 0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4,
    0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac, 0x6c, 0xec, 0x1c, 0x9c, 0x5c, 0xdc, 0x3c, 0xbc, 0x7c, 0xfc,
    0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2, 0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2,
    0x0a, 0x8a, 0x4a, 0xca, 0x2a, 0xaa, 0x6a, 0xea, 0x1a, 0x9a, 0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa,
    0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6, 0x16, 0x96, 0x56, 0xd6, 0x36, 0xb6, 0x76, 0xf6,
    0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee, 0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe, 0x7e, 0xfe,
    0x01, 0x81, 0x41, 0xc1, 0x21, 0xa1, 0x61, 0xe1, 0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
    0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9, 0x19, 0x99, 0x59, 0xd9, 0x39, 0xb9, 0x79, 0xf9,
    0x05, 0x85, 0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5, 0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5,
    0x0d, 0x8d, 0x4d, 0xcd, 0x2d, 0xad, 0x6d, 0xed, 0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd,
    0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3, 0x63, 0xe3, 0x13, 0x93, 0x53, 0xd3, 0x33, 0xb3, 0x73, 0xf3,
    0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb, 0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb,
    0x07, 0x87, 0x47, 0xc7, 0x27, 0xa7, 0x67, 0xe7, 0x17, 0x97, 0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7,
    0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef, 0x1f, 0x9f, 0x5f, 0xdf, 0x3f, 0xbf, 0x7f, 0xff
  };

  // kmer reverse and complement table
  static const uint8 kractable[] =
  {
    0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f, 0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0xf,
    0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b, 0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0xb,
    0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27, 0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x7,
    0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23, 0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x3,
    0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e, 0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0xe,
    0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a, 0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0xa,
    0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26, 0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x6,
    0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22, 0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x2,
    0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d, 0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0xd,
    0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29, 0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x9,
    0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25, 0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x5,
    0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21, 0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x1,
    0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c, 0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0xc,
    0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28, 0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x8,
    0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24, 0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x4,
    0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20, 0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x0
  };
}

template<unsigned int VLBITS>
class VLuint
{
protected:
  typedef uint64 base_type;

  static_assert(std::is_pod<base_type>::value,"base_type must be a POD type in template class VLUint<>");
  static_assert(sizeof(base_type) >= 4, "base_type must be at least 32 bit wide in template class VLUint<>");
  static_assert(sizeof(base_type) == 8, "base_type must be 64 bit wide in template class VLUint<> for function reverseComplement(). If not, need to rewrite that function!");

  enum { WIDTH=VLBITS/(8*sizeof(base_type)) };
  base_type payload[WIDTH];

  typedef VLuint<VLBITS> Self;

  //Functions
private:
  // general comparison function
  // to use like in operator<()
  inline int priv_compareTo(const Self & other) const {
    for(int i = WIDTH-1; i >= 0; --i) {
      if(payload[i] < other.payload[i]) return -1;
      if(payload[i] > other.payload[i]) return 1;
    }
    return 0;
  }

  // so, here's one hell of a problem: the code below segfaults quickly in
  // mirabait (compiled with a VLuint<512>) on large k (e.g. -k 255)!
  //
  // This is *totally* the same code as the general priv_compareTo() function
  // above, except one if-line has been taken out. And it segfaults??? I have not
  // the slightest idea why ... I'm scared.
  //
  // seen 16.08.2014 on OSX, gcc 4.8.2, no optimisation or full -O3 + others
  //
  // inline int priv_compareLess(const Self & other) const {
  //   for(int i = WIDTH-1; i >= 0; --i) {
  //     if(payload[i] < other.payload[i]) return -1;
  //     // if(payload[i] > other.payload[i]) return 1;
  //   }
  //   return 0;
  // }

public:
  // STD. CONSTRUCTOR
  //
  // Note that no memory is cleared here, i.e., can contain random values upon
  // instantiation
  VLuint() {};

  // ASSIGNMENT CONSTRUCTOR
  //
  // I have them explicit to find/prevent unexpected conversions
  explicit VLuint(base_type x) {
    //cout << "ctor= int " << this << " ";
    payload[0]=x;
    std::fill(&payload[1], &payload[WIDTH], 0);
  };

  // Std. output operator (which is more for debugging than anything really useful)
  friend std::ostream & operator<<(std::ostream &ostr, Self const &vl) {
    //ostr << "W " << WIDTH << " ";
    for(uint32 i=WIDTH; i>0; ){
      ostr << std::hex << vl.payload[--i] << " ";
    }
    return ostr;
  }

  // ASSIGNMENT
  Self & operator=(const Self &other) {
    if(&other != this) {
      std::copy(&other.payload[0], &other.payload[WIDTH], &payload[0]);
    }
    return *this;
  }

  Self & operator=(const base_type x) {
    //cout << "opor= int " << this << " ";
    payload[0]=x;
    std::fill(&payload[1], &payload[WIDTH], 0);
    return *this;
  }

  // COMPARATORS
  inline bool operator==(const Self & other) const {
    return std::equal(&other.payload[0], &other.payload[WIDTH], &payload[0]);
  }
  inline bool operator!=(const Self & other) const {
    return !(*this==other);
  }

  inline bool operator<(const Self & other) const {
    return (priv_compareTo(other) < 0);
  }

  // Segfaults, see comments for priv_compareLess() above
  //inline bool operator<(const Self & other) const {
  //  return (priv_compareLess(other) < 0);
  //}


  // OPERATORS
  inline Self & operator--() {
    for(int i = 0; payload[i++]--==0 && i < WIDTH; ) {};
    return *this;
  }
  inline Self & operator++() {
    for(int i = 0; ++payload[i++]==0 && i < WIDTH; ) {};
    return *this;
  }

  inline Self operator&=(const Self & other) {
    for(int i = 0; i < WIDTH; ++i) {payload[i]&=other.payload[i];}
    return *this;
  }

  inline Self operator|=(const Self & other) {
    for(int i = 0; i < WIDTH; ++i) {payload[i]|=other.payload[i];}
    return *this;
  }

  inline Self operator~() const {
    Self ret;
    for(int i = 0; i < WIDTH; ++i) {ret.payload[i]=~payload[i];}
    return ret;
  }

  inline Self operator&(const Self & other) const {
    Self ret;
    for(int i = 0; i < WIDTH; ++i) {ret.payload[i]=payload[i] & other.payload[i];}
    return ret;
  }

// disabled, 50% slower in sort() than operator&(Self &) above
//  inline Self operator&(const base_type & other) const {
//    Self ret;
//    ret.payload[0]=payload[0] & other;
//    std::fill(&ret.payload[1], &ret.payload[WIDTH], 0);
//    return ret;
//  }

  inline Self operator|=(const base_type & other) {
    payload[0]|=other;
    return *this;
  }

  inline Self & operator>>=(base_type shift) {
    Self tmp(*this);
    std::fill(&payload[0], &payload[WIDTH], 0);
    int k = shift / (sizeof(base_type)*8);
    shift = shift % (sizeof(base_type)*8);
    for (int i = 0; i < WIDTH; ++i) {
      if (i-k-1 >= 0 && shift != 0) payload[i-k-1] |= (tmp.payload[i] << ((sizeof(base_type)*8)-shift));
      if (i-k >= 0) payload[i-k] |= (tmp.payload[i] >> shift);
    }
    return *this;
  }
  inline Self operator>>(base_type shift) const {
    return Self(*this) >>= shift;
  }

  inline Self & operator<<=(base_type shift) {
    Self tmp(*this);
    std::fill(&payload[0], &payload[WIDTH], 0);
    int k = shift / (sizeof(base_type)*8);
    shift = shift % (sizeof(base_type)*8);
    for (int i = 0; i < WIDTH; ++i) {
      if (i+k+1 < WIDTH && shift != 0) payload[i+k+1] |= (tmp.payload[i] >> ((sizeof(base_type)*8)-shift));
      if (i+k < WIDTH) payload[i+k] |= (tmp.payload[i] << shift);
    }
    return *this;
  }

  // CONVERTORS
  //
  // I have them explicit to find/prevent unexpected conversions
  // BaCh: 06.02.2015
  // Do not use size_t as cast target: on OSX, uint64 and size_t
  //  are NOT the same types! And I cannot make two conversion
  //  functions as under Linux (and others?), they are the same
  //  and the compiler complains about ambiguous conversions.
  // The base type of vluint has been defined as uint64 anyway.
  inline explicit operator uint64() {
    return payload[0];
  }
  inline explicit operator bool() {
    bool ret=false;
    for(int i = 0; i < WIDTH; ++i) {
      ret|=payload[i];
    };
    return ret;
  }

  // HASH FUNCTION
  //
  // to facilitate simple hashing where needed in C++ containers like
  // unordered_set/*_map. See vhash.H on how to do this.
  inline size_t hash() const {
    size_t ret=0;
    for (int hi = 0; hi < WIDTH; ++hi) {
      //boost::hash_combine(ret, payload[hi]);
      // To not to have to include boost here, do as BOOST does in hash_combine()
      ret ^= std::hash<base_type>()(payload[hi]) + 0x9e3779b9 + (ret << 6) + (ret >> 2);
    }
    // this was for testing, does not seem good
    //    ret ^= (ret>>32) | (ret<<32);
    return ret;
  }

  // Other
  // this function + the table it uses should somehow move to a vhash class!
  inline Self gnagna__this_is_so_badly_designed__kmerReverseComplement() const {
    Self ret;
    for(int pi = 0; pi < WIDTH; ++pi) {
      base_type oldpl=payload[pi];
      ret.payload[WIDTH-1-pi]=
	(static_cast<base_type>(nsVLuint::kractable[oldpl&0xFF])<<56)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >>  8) &0xFF])<<48)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >> 16) &0xFF])<<40)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >> 24) &0xFF])<<32)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >> 32) &0xFF])<<24)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >> 40) &0xFF])<<16)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >> 48) &0xFF])<<8)
	| (static_cast<base_type>(nsVLuint::kractable[(oldpl >> 56) &0xFF]))
	;
    }
    return ret;
  }

/*
  // Not sure if needed
  inline void reverseComplement() {
    Self tmp(*this);
    for(int pi = 0; pi < WIDTH; ++pi) {
      base_type oldpl=tmp.payload[pi];
      payload[WIDTH-1-pi]=
	(static_cast<base_type>(nsVLuint::bractable[oldpl&0xFF])<<56)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >>  8) &0xFF])<<48)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >> 16) &0xFF])<<40)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >> 24) &0xFF])<<32)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >> 32) &0xFF])<<24)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >> 40) &0xFF])<<16)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >> 48) &0xFF])<<8)
	| (static_cast<base_type>(nsVLuint::bractable[(oldpl >> 56) &0xFF]))
	;
    }
  }
*/

};


#endif
