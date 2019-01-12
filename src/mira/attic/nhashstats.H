/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 by Bastien Chevreux
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

#ifndef _bas_hashstats_h_
#define _bas_hashstats_h_

#include <unordered_map>

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/regex.hpp>

#include "stdinc/defines.H"

#include "mira/readpool.H"

template<class TVHASH_T> class BloomFilter;


#include <zlib.h>



template<typename TVHASH_T>
class NHashStatistics
{
#ifdef HSVHM_var
  static const TVHASH_T HS_MAXVHASHMASK;
#endif

public:
  // new hashstat

  typedef typename HashStatistics<TVHASH_T>::hscounts_t hscounts_t;
  typedef typename HashStatistics<TVHASH_T>::hashstat_t nhashstat_t;

private:
  // begin and end iterators
  struct hsvbendit_t {
    typename std::vector<nhashstat_t>::const_iterator b;
    typename std::vector<nhashstat_t>::const_iterator e;
  };

private:

  static uint32 HSN_hs_magic;

  BloomFilter<TVHASH_T> * HSN_bloomfilter;

  std::vector<nhashstat_t>                  HSN_hsv_hashstats;
  std::vector<hsvbendit_t>                  HSN_hsv_hsshortcuts;

  std::unordered_map<TVHASH_T,hscounts_t>    HSN_hsum_hashstats;

  std::vector<uint64>                                    HSN_hs_dist;

  uint8 HSN_hs_sortstatus;
  bool HSN_hs_needsconsolidation;

  uint32 HSN_basesperhash;

  std::string HSN_directory; // dir to write hashstat to

  uint16 HSN_needsteps;
  uint16 HSN_step;


private:
  void makeNHashStatArrayShortcuts(std::vector<nhashstat_t> & nhashstats,
				   const uint32 basesperhash,
				   std::vector<hsvbendit_t> & hsshortcuts);
  inline static bool sortHashStatComparatorByLow24bit(const nhashstat_t & a, const nhashstat_t & b){
    if((a.vhash & HS_MAXVHASHMASK) != (b.vhash & HS_MAXVHASHMASK)) {
      return (a.vhash & HS_MAXVHASHMASK) < (b.vhash & HS_MAXVHASHMASK);
    }
    return a.vhash < b.vhash;
  }
  inline static bool sortHashStatComparator(const nhashstat_t & a, const nhashstat_t & b){
    return a.vhash < b.vhash;
  }
  void sortLow24Bit(std::vector<nhashstat_t> & hashstats, uint8 & sortstatus) {
    std::sort(hashstats.begin(), hashstats.end(), sortHashStatComparatorByLow24bit);
    sortstatus=HashStatistics<TVHASH_T>::HSSS_LOW24BIT;
    HSN_hsv_hsshortcuts.clear();
  }
  void sortLexicographically(std::vector<nhashstat_t> & hashstats, uint8 & sortstatus) {
    std::sort(hashstats.begin(), hashstats.end(), sortHashStatComparator);
    sortstatus=HashStatistics<TVHASH_T>::HSSS_LEXIUP;
    HSN_hsv_hsshortcuts.clear();
  }
  uint64 calcHashDistrib(std::vector<uint64> & hsdist);

  void learnSequenceQuick1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse, bool lookuponly);
  void learnSequenceQuick2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void learnSequenceStep1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void learnSequenceStep2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void learnSequenceStep3(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);

  void saveHashVStatistics(std::ostream & ostr);
  void saveHashMStatistics(std::ostream & ostr);

  void trimHashVStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal);
  void trimHashMStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal);

public:
  NHashStatistics() : HSN_bloomfilter(nullptr), HSN_basesperhash(0), HSN_hs_sortstatus(HashStatistics<TVHASH_T>::HSSS_NOTSORTED) {};
  ~NHashStatistics();

  void setupNewAnalysis(const uint8  bfbits,
			const uint32 bfnumkeys,
			const uint32 basesperhash,
			uint16 numsteps);
  void analyseReadPool(ReadPool & rp);
  void deleteBloomFilter();

  //void prefetchSequence(const void * seqvoid, uint64 slen, const char * namestr);
  void learnSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void finaliseStep();

  const nhashstat_t * findVHash(const nhashstat_t & searchval);

  void trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal);

  std::vector<uint64> & getHashDistrib();
  void dumpHashDistrib(std::ostream & ostr);

  void dumpHealth(std::ostream & fout);
  void dumpHashCount(std::ostream & ostr);

  void saveHashStatistics(std::ostream & ostr);
  void loadHashStatistics(std::istream & istr);

  void saveHashStatistics(const std::string & filename, bool deleteoldfile);
  void loadHashStatistics(const std::string & filename);

  void sortLow24Bit(){
    sortLow24Bit(HSN_hsv_hashstats, HSN_hs_sortstatus);
  }
  void sortLexicographically(){
    sortLexicographically(HSN_hsv_hashstats, HSN_hs_sortstatus);
  }

  void moveStatCountMapToVector();


};


#endif
