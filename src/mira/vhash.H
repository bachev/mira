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

#ifndef _bas_vhash_h_
#define _bas_vhash_h_

#include "mira/vluint.H"


// TODO
// the typedef should somehow move to a vhash class which
// also implements reverseComplement
// the current code is really badly designed as VLuint should
// not have this kludge function gnagna__this_is_so_badly_designed__kmerReverseComplement()


#if 0
#define KMER_INTERNALTYPE
  typedef VLuint<64> vhash64_t;
  typedef VLuint<64> vhash128_t;
  typedef VLuint<64> vhash256_t;
  typedef VLuint<64> vhash512_t;

namespace std {
  template <> struct hash<VLuint<64> >
  {
    size_t operator()(const VLuint<64> & x) const
    {
      return x.hash();
    }
  };
};

#else

typedef uint64_t vhash64_t;        // faster
//  typedef VLuint<64> vhash64_t;
  typedef VLuint<128> vhash128_t;
  typedef VLuint<256> vhash256_t;
  typedef VLuint<512> vhash512_t;

// Specialisation of namespace std so that VLuint knows how to compute hashes of itself
// Needed for e.g., unordered_set or *_map
// See http://stackoverflow.com/questions/8157937/how-to-specialize-stdhashkeyoperator-for-user-defined-type-in-unordered?rq=1

// how to make this happen automatically for any instantiated VLuint size???

namespace std {
  template <> struct hash<VLuint<64> >
  {
    size_t operator()(const VLuint<64> & x) const
    {
      return x.hash();
    }
  };
};
namespace std {
  template <> struct hash<VLuint<128> >
  {
    size_t operator()(const VLuint<128> & x) const
    {
      return x.hash();
    }
  };
};
namespace std {
  template <> struct hash<VLuint<256> >
  {
    size_t operator()(const VLuint<256> & x) const
    {
      return x.hash();
    }
  };
};
namespace std {
  template <> struct hash<VLuint<512> >
  {
    size_t operator()(const VLuint<512> & x) const
    {
      return x.hash();
    }
  };
};
#endif

// just contains reverseComplement() atm ... should be somehow a derived class
//  of VLuint<> for every vhash<>_t
namespace nsvhash {
  // specialisation of reverseComplement() for uint64_t
  inline uint64 reverseComplement(uint64 vh, uint32 basesperhash) {
    auto oldpl=vh;
    uint64_t ret=
	(static_cast<uint64>(nsVLuint::kractable[oldpl&0xFF])<<56)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >>  8) &0xFF])<<48)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >> 16) &0xFF])<<40)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >> 24) &0xFF])<<32)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >> 32) &0xFF])<<24)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >> 40) &0xFF])<<16)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >> 48) &0xFF])<<8)
	| (static_cast<uint64>(nsVLuint::kractable[(oldpl >> 56) &0xFF]))
	;
    ret>>=(sizeof(uint64)*8-basesperhash*2);
    return ret;
  }
  // and template version for the rest
  template<class T>
  T reverseComplement(T & vh, uint32 basesperhash){
    return vh.gnagna__this_is_so_badly_designed__kmerReverseComplement()>>(sizeof(T)*8-basesperhash*2);
  }
}

#endif
