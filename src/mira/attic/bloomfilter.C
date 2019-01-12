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

#include "bloomfilter.H"
#include "util/dptools.H"

using namespace std;


// Plain vanilla constructor
template<class TVHASH_T>
BloomFilter<TVHASH_T>::BloomFilter(const uint8 bits, const uint32 numkeys)
{
  FUNCSTART("BloomFilter<TVHASH_T>::BloomFilter()");

  cout << "initialising Bloom filter with " << static_cast<uint16>(bits) << " bits (" << (0x1ull<<bits) << " elements), " << (0x1ull<<bits)/4 << " bytes\n";

  BUGIFTHROW(bits>64,"bits >64 ?");
  BUGIFTHROW(numkeys==0||numkeys>20,"numkeys == " << numkeys << " ???");
  BF_numkeys=numkeys;
  BF_bfaddressmask=(0x1uLL<<bits)-1;
  BF_bloomfield.resize((0x1ull<<bits)/4);
  reset();

  FUNCEND();
}


template<class TVHASH_T>
BloomFilter<TVHASH_T>::~BloomFilter()
{
  FUNCSTART("BloomFilter<TVHASH_T>::~BloomFilter()");

  discard();

  FUNCEND();
}


template<class TVHASH_T>
void BloomFilter<TVHASH_T>::discard()
{
  FUNCSTART("BloomFilter<TVHASH_T>::discard()");

  BF_bloomfield.clear();
  BF_level1count=0;
  BF_level2count=0;
  BF_numuniqkmers=0;
  BF_numkmerseenge2=0;
  BF_numkmerseenge3=0;

  FUNCEND();
}

template<class TVHASH_T>
void BloomFilter<TVHASH_T>::reset()
{
  discard();
  BF_bloomfield.resize(BF_bloomfield.capacity(),0);
}




// explicit template instantiations needed for the linker to create these for library files
template class BloomFilter<vhash64_t>;
#ifndef KMER_INTERNALTYPE
template class BloomFilter<vhash128_t>;
template class BloomFilter<vhash256_t>;
template class BloomFilter<vhash512_t>;
#endif
