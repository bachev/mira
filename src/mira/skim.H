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

#ifndef _bas_skim_h_
#define _bas_skim_h_

#include <boost/function.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

#include "stdinc/defines.H"

#include "mira/types_basic.H"
#include "mira/hashstats.H"

class ADSEstimator;
class ReadPool;
template <class T> class ProgressIndicator;


//#define SKIM3_MAXVHASHMASK 0xFFFFFFLL


#define SKIM3_MAXREADSIZEALLOWED MIRASRLEN

// 10G
#define SKIM3_SKIMMATCHFIRSTCHECK 10737418240LL
#define SKIM3_SKIMMATCHCHECKINCR 10737418240LL



typedef uint32 shash_t;


struct matchwithsorter_t{
  uint32 otherid;
  int32  eoffset;
  int32  percent_in_overlap;
  uint32 numhashes;                 // try to deprecate?

  uint32 estimscore;

  bool   taken;

  bool ol_stronggood:1;  // frequency: 2*bph-1 pos at 3, thereof bph-1 contiguous
  bool ol_weakgood:1;   // frequency: bph-1 positions contiguous at 3
  bool ol_belowavgfreq:1;   // frequency: bph-1 positions contiguous at <=3
  bool ol_norept:1;      // nothing >3 (but can contain 1 (single hashes == errors)
  bool ol_rept:1;      // bph-1 positions >=5 frequency

  // not saved to disk but used internally
  bool ol_fulllength:1;   // overlap is full length of one of the reads
  bool ol_fullencased:1;  /* got the status of fullencased:
			     100% match, no rept, len difference >=8
			  */

  friend std::ostream & operator<<(std::ostream &ostr, const matchwithsorter_t & m){
    ostr << "MWS:\t" << m.otherid
	 << '\t' << m.eoffset
	 << '\t' << m.percent_in_overlap
	 << '\t' << m.numhashes
	 << '\t' << m.estimscore
	 << "\ttak " << m.taken
	 << "\tsg  " << m.ol_stronggood
	 << "\twg  " << m.ol_weakgood
	 << "\tbaf " << m.ol_belowavgfreq
	 << "\tnrp " << m.ol_norept
	 << "\trep " << m.ol_rept
	 << "\tfll " << m.ol_fulllength
	 << "\tfle " << m.ol_fullencased
	 << '\n';
    return ostr;
  }
};


// only similar to matchwithsorter_t
struct skimhitforsave_t{
  uint32 rid1;
  uint32 rid2;
  int32  eoffset;
  int32  percent_in_overlap;
  uint32 numhashes;

  bool ol_stronggood:1;  // frequency: 2*bph-1 pos at 3, thereof bph-1 contiguous
  bool ol_weakgood:1;   // frequency: bph-1 positions contiguous at 3
  bool ol_belowavgfreq:1;   // frequency: bph-1 positions contiguous at <=3
  bool ol_norept:1;      // nothing >3 (but can contain 1 (single hashes == errors)
  bool ol_rept:1;      // bph-1 positions >=5 frequency

  // the following is used only by PacBio correction atm, should expand this
private:
  bool rid1dir:1;
  bool rid2dir:1;

public:
  inline int8 getRID1Dir() const { return (rid1dir ? 1 : -1);}
  inline int8 getRID2Dir() const { return (rid2dir ? 1 : -1);}
  inline void setRID1Dir(int8 d) { rid1dir=d>0;}
  inline void setRID2Dir(int8 d) { rid2dir=d>0;}

  friend std::ostream & operator<<(std::ostream &ostr, const skimhitforsave_t & e){
    ostr << "SHFS:\t" << e.rid1
	 << '\t' << e.rid2
	 << '\t' << e.rid1dir
	 << '\t' << e.rid2dir
	 << "\teo " << e.eoffset
	 << "\t% " << e.percent_in_overlap
	 << "\tnh " << e.numhashes
	 << "\tsg " << e.ol_stronggood
	 << "\twg " << e.ol_weakgood
	 << "\tbaf " << e.ol_belowavgfreq
	 << "\tnrp " << e.ol_norept
	 << "\trp " << e.ol_rept
	 << "\tr1d " << static_cast<int16>(e.getRID1Dir())
	 << "\tr2d " << static_cast<int16>(e.getRID2Dir())
	 << '\n';
    return ostr;
  }


  static bool stdSortCmp(const skimhitforsave_t & a, const skimhitforsave_t & b);
};



struct matchwith_t{
  uint32 otherid;
  int32  eoffset;
};


struct readhashmatch_t{
  uint32 rid2;
  int32  eoffset;
  uint16 hashpos1;
  uint16 hashpos2;
  Read::bhashstat_t bhashstats; // baseflags for this hash

  friend std::ostream & operator<<(std::ostream &ostr, const readhashmatch_t & rhm){
    ostr << "rid2: " << rhm.rid2
	 << "\teoffset: " << rhm.eoffset
	 << "\thp1: " << rhm.hashpos1
	 << "\thp2: " << rhm.hashpos2
	 << "\tbhs: " << rhm.bhashstats
	 << '\n';
    return ostr;
  }

};

typedef std::multimap< int32, matchwith_t> possible_overlaps_t;
typedef possible_overlaps_t::value_type posoverlap_pair_t;


// bop_t is organised as follows:
//  each ban of a pair of reads is stored in the vector of the read with the lower ID
// bop_t was initially std::vector<set<uint32> >, but these beasts are extremely expensive,
//  both in time and in memory (48 per set + 48 per element in set!)
// using a vector is way quicker (around 3 times), uses way less memory (24 + 4 per element)
//  and the growth strategy ensures that not too much memory is wasted (has a small hit on
//  run-time efficiency though as we might need to copy quite a couple of times)

struct bannedoverlappairs_t {
  typedef std::vector<std::vector<uint32> > bop_t;
  bop_t          bop;

  inline void nuke() {nukeSTLContainer(bop);};
  inline size_t size() {return bop.size();};
  inline void resize(size_t sz) {bop.resize(sz);};
  void insertBan(uint32 id1, uint32 id2) {
    FUNCSTART("void insertBan(uint32 id1, uint32 id2)");
    if(id1>id2) std::swap(id1,id2);
    BUGIFTHROW(id1>=bop.size(),"id1 " << id1 << " >= bop.size() " << bop.size())
    std::vector<uint32>::iterator iI=upper_bound(bop[id1].begin(),bop[id1].end(),id2);
    if(iI==bop[id1].end() || iI==bop[id1].begin() || (iI != bop[id1].begin() && (*(iI-1) != id2))){
      // must insert
      // moderate growth strategy if needed ... most reads will not have thousands of bans
      if(bop[id1].capacity()>=128
	 && bop[id1].capacity() == bop[id1].size()){
	// if we take over memory management ourselves, need to be careful:
	//  reserve() will invalidate our iI pointer (insert() handles that
	//  automatically)
	// So, just make sure we keep a valid pointer
	size_t distance=iI-bop[id1].begin();
	bop[id1].reserve(bop[id1].size()+128);
	iI=bop[id1].begin()+distance;
      }
      bop[id1].insert(iI,id2);
    }
  }
  inline bool checkIfBanned(uint32 id1, uint32 id2) {
    FUNCSTART("inline bool checkIfBanned(uint32 id1, uint32 id2)");
    if(id1>id2) std::swap(id1,id2);
    BUGIFTHROW(id1>=bop.size(),"id1 " << id1 << " >= bop.size() " << bop.size())
    return binary_search(bop[id1].begin(),bop[id1].end(),id2);
  }
  void getNumBans(size_t & banned, size_t & numsets) {
    banned=0;
    numsets=0;
    for(size_t i=0; i<bop.size(); i++){
      if(!bop[i].empty()){
	banned+=bop[i].size();
	++numsets;
      }
    }
  }
};

typedef std::map< uint32, uint32> simplematches_t;
typedef simplematches_t::value_type simplematch_pair_t;

template<typename TVHASH_T>
class Skim
{
  static const TVHASH_T SKIM3_MAXVHASHMASK;

  /****************************************************************************
   *
   * FARC
   *
   ****************************************************************************/
private:
  // each findAdaptorRightClip thread needs a couple of data structures ...
  struct farc_threaddata_t {
    std::vector<readhashmatch_t> readhashmatches;
    std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> singlereadvhraparray;
    std::vector<uint8> unused_tagmaskvector;
    std::vector<matchwithsorter_t> tmpmatchwith;
  };

  farc_threaddata_t  SKIM3_farcdata_fornonmultithread;
  ReadPool         * SKIM3_farc_searchpool;
  std::vector<int32>    * SKIM3_farc_results;
  uint32             SKIM3_farc_minhashes;
  uint8              SKIM3_farc_seqtype;
  std::vector<farc_threaddata_t> SKIM3_farcd_vector;

  // functions
  void farcThreadsDataInit(const uint32 threadnr);
  void farcThreadLoop(const uint32 threadnr);
  void checkForPotentialAdaptorHits(const int8 direction,
				    const uint32 actreadid,
				    Read & actread,
				    std::vector<matchwithsorter_t> & tmpmatchwith,
				    std::vector<readhashmatch_t> & readhashmatches);

  int32 findAdaptorRightClip_internal(Read & actread, uint32 minhashes, readid_t & ridofadapfound, farc_threaddata_t & farcd);

public:
  // These two calls need no preparation by the caller
  //int32 findAdaptorRightClip(Read & actread, uint32 minhashes);    // call for simple, non-multithreaded check of single read
  void findAdaptorRightClip(ReadPool & searchpool,
			    std::vector<int32> & results,
			    int8 seqtype,
			    uint32 minhashes,
			    uint32 numthreads); // call for multithreaded check of whole read pool

  // These calls for non-multithreaded and multi-threaded usage, caller decides everything
  void  prepareForMultithreadFarc(uint32 numthreads);
  int32 findAdaptorRightClip(Read & actread,
			     uint32 minhashes,
			     readid_t & ridofadapfound,
			     int32 threadid);  // call for multithreaded check of single read, caller decides in which thread structure ID it runs
                                               // threadid==-1 means: non-multithreaded version (same as above), needs no prepareForMultithreadFarc




  /****************************************************************************
   *
   * low BPH skim
   *
   ****************************************************************************/
private:

  struct srhashhits_t {
    uint32 numhashes;
    uint32 rid;

    // yeah, well, I know: this is ... errmmm ... well
    // but that function is really only needed for a high/low sort
    inline bool operator<(const srhashhits_t & other) const {
      return numhashes>other.numhashes;
    }
  };

  struct shashpos_t {
    shash_t shash;
    int32   hashpos;
  };

  struct lbphs_threaddata_t {
    std::vector<shash_t>  hashvec;
    std::vector<uint32>   hashhitperread;

    std::vector<srhashhits_t> hashhitsofoneread;

    std::vector<std::vector<int32> > r1posperhash;
    std::vector<shash_t>        r1pphlist;
    std::vector<shashpos_t> r2shashpos;

    std::vector<uint32>  r1r2poshistogram;
    std::vector<uint32>  r1r2phused;
    std::vector<uint32>  r1r2idstaken;

    std::vector<skimhitforsave_t> shfsvf;
    std::vector<skimhitforsave_t> shfsvr;
  };

  std::vector<lbphs_threaddata_t> SKIM3_lbphsd_vector;

  // 2ndtry
  std::vector<std::vector<uint32> >  SKIM3_lbphs_idsperhash;
  std::vector<uint32>           SKIM3_lbphs_numoverlapsperid;
  std::vector<uint32>           SKIM3_lbphs_maxoverlaphashesperid;

  uint64 totalphits;

  static shash_t SKIM3_lbphs_hashadd[128];


  // functions
  void lowBPHSkim();
  void lbphsThreadsDataInit(const uint32 numthreads);
  void lbphsThreadLoop(const uint32 threadnr);
  void lbphsLookAtRead(const uint32 actreadi, const uint32 threadnr, const int8 direction);

  void lbphsPrepareHashOverviewTable(uint32 & readi);
  bool lbphsPrepareOneSeqForHOT(const uint32 readi);
  void lbphsGetHashesOfOneSeq(const uint32 readi,
			      std::vector<shash_t> & hashvec,
			      const uint32 basesperhash,
			      const char * seq,
			      const uint32 slen
    );

  void lbphsSecondLineCheck(const uint32 actreadid, lbphs_threaddata_t & lbphsd, int8 direction);
  void lbphsSeqToPosOverviewTable(std::vector<std::vector<int32> > & posperhash,
				  std::vector<shash_t> & r1pphlist,
				  const uint32 basesperhash,
				  const char * seq,
				  const uint32 slen);
  void lbphsGetHashPosOfOneSeq(const uint32 readi,
			       std::vector<shashpos_t> & hashposvec,
			       const uint32 basesperhash,
			       const char * seq,
			       const uint32 slen);



  /****************************************************************************
   *
   * Core Skim
   *
   ****************************************************************************/
private:
  bool SKIM3_logflag_purgeunnecessaryhits;
  bool SKIM3_logflag_save2;

  ReadPool * SKIM3_readpool;

  typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> SKIM3_vhraparray;
  typename std::vector<typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::const_iterator > SKIM3_vashortcuts_begin;
  typename std::vector<typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::const_iterator > SKIM3_vashortcuts_end;

  // the "nullptr" replacement for the shortcut arrays above
  typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::const_iterator SKIM3_completevhraparray_end;

  // TODO: eventuall compress this into uint8 having 8 boolean
  std::vector<uint8> SKIM3_hasMNRr;
  std::vector<uint8> SKIM3_hasSRMr;
  std::vector<uint8> SKIM3_hasFpas;

  // counts how many times a read has excellent hits with other reads
  std::vector<uint32> * SKIM3_overlapcounter;

  // count how many times a read was saved as fully encased
  std::vector<uint8> SKIM3_fullencasedcounter;

  std::vector<uint32> * SKIM3_writtenhitsperid;
  std::vector<uint8> * SKIM3_megahubsptr;

  bannedoverlappairs_t * SKIM3_bannedoverlaps;

  std::ofstream SKIM3_posfmatchfout;
  std::ofstream SKIM3_poscmatchfout;
  std::string SKIM3_posfmatchfname;
  std::string SKIM3_poscmatchfname;
  uint64 SKIM3_posfmatchnextchecksize;
  uint64 SKIM3_poscmatchnextchecksize;

  uint32 SKIM3_megahubcap;

  uint32 SKIM3_numthreads;
  uint32 SKIM3_basesperhash;
  uint8  SKIM3_hashsavestepping;

  //int32  SKIM3_percentrequired;

  std::vector<int32>  SKIM3_percentrequired;
  std::vector<int32>  SKIM3_overlaplenrequired; // min len per sequencing type

  // if vector.size() != 0
  // skimGo() also hunts chimeras
  // outer vector size of readpool
  // inner vectors size of length of clipped sequences
  //  0 initially, set to 1 for consecutive hash matches
  //  chimeras show up as a 0 stretch within the 1s
  std::vector<std::vector<uint8> > SKIM3_chimerahunt;
  std::vector<int32> * SKIM3_chuntleftcut;
  std::vector<int32> * SKIM3_chuntrightcut;

  // overlap criterion levels for left and right extend
  // used to take overlaps only if they reach a given level
  // [0] is for norept overlaps, [1] for rept
  std::vector<std::vector<uint8>> * SKIM3_overlapcritlevelvl;
  std::vector<std::vector<uint8>> * SKIM3_overlapcritlevelvr;

  std::vector<uint32> SKIM3_largestencasementscoretodate;

  // for ultra-deeply connected reads: enough is enough
  // set to 1 for a read which whould not be looked at anymore
  //  by skim because it already has too many hits with
  //  overlap criterion 0 left/righ on disk
  // (currently only set for Solexa reads)
  std::vector<uint8>  SKIM3_nomorehitseie;

  uint32 SKIM3_maxhitsperread;
  bool SKIM3_forcetakestronggood;

  uint64 SKIM3_totalhitschosen;


  uint32 SKIM_partfirstreadid;
  uint32 SKIM_partlastreadid;

  ProgressIndicator<int64> * SKIM_progressindicator;
  int64 SKIM_progressend;

  bool SKIM3_onlyagainstrails;

  // each checkForHashes thread needs a couple of data structures ...
  struct cfh_threaddata_t {
    std::vector<readhashmatch_t> readhashmatches;
    std::vector<uint32> smallhist4repeats;
    std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> singlereadvhraparray;
    std::vector<matchwithsorter_t> tmpmatchwith;
    std::vector<uint8> tagmaskvector;
    std::vector<skimhitforsave_t> shfsv;
    std::ofstream * posmatchfout;
    // this vector is used to collect read ids which have a match
    //  and then quickly update SKIM3_writtenhitsperid inside a mutex
    std::vector<uint32> ridswithmatches;
    std::vector<std::vector<uint32>> uidswithnewcritlevelvl;  // 0 for norept, 1 for rept overlaps
    std::vector<std::vector<uint32>> uidswithnewcritlevelvr;
    std::vector<std::vector<uint8>> critlevellofnewuidsv;
    std::vector<std::vector<uint8>> critlevelrofnewuidsv;
  };

  std::vector<cfh_threaddata_t> SKIM3_cfhd_vector;


  struct threadworkercontrol_t {
    uint32 from;
    uint32 to;
    int8   direction;
    bool   flag_datavalid;
    bool   flag_endthread;
  };

  std::vector<threadworkercontrol_t> SKIM3_threadcontrol;

  boost::mutex SKIM3_coutmutex;
  boost::mutex SKIM3_resultfileoutmutex;
  boost::mutex SKIM3_globalclassdatamutex;
  //boost::mutex SKIM3_whpid_mutex;

  boost::mutex SKIM3_critlevelwrite_mutex;

  boost::mutex SKIM3_mutex;
  boost::condition SKIM3_master2slavesignal;
  boost::condition SKIM3_slave2mastersignal;

public:


  //Functions
private:
  void foolCompiler();

//  void prepareSkim(bool alsocheckreverse);
  void prepareSkim(uint32 fromid, uint32 toid, std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> & vhraparray, bool assemblychecks);
  void purgeMatchFileIfNeeded(int8 direction);
  void findPerfectRailMatchesInSkimFile(std::string & filename, const int8 rid2dir, std::vector<uint8> & prmatches);
  void purgeUnnecessaryHitsFromSkimFile(std::string & filename, const int8 rid2dir, std::vector<uint8> & prmatches);

  bool fillTagMaskVector(const uint32 readid, std::vector<uint8> & tagmaskvector);
  void fillTagStatusInfoOfReads();

  void startMultiThreading(const int8 direction,
			   const uint32 numthreads,
			   const uint32 readsperthread,
			   const uint32 firstid,
			   const uint32 lastid,
			   boost::function<void(uint32_t)> initfunc,
			   boost::function<void(uint32_t)> callfunc);
  void cfhThreadsDataInit(const uint32 numthreads);
  void cfhThreadLoop(const uint32 threadnr);
  void checkForHashes_fromto(const int8 direction,
			     const uint32 fromid,
			     const uint32 toid,
			     cfh_threaddata_t & cfhd);
  void checkForPotentialHits(const int8 direction,
			     const uint32 actreadid,
			     std::vector<matchwithsorter_t> & tmpmatchwith,
			     std::vector<readhashmatch_t> & readhashmatches,
			     std::vector<uint32> & smallhist4repeats);

  void selectPotentialHitsForSave2(const int8 direction,
				   const uint32 actreadid,
				   cfh_threaddata_t & cfhd);
  void updateCriterionLevels(const int8 direction,
			     const uint32 actreadid,
			     cfh_threaddata_t & cfhd);

  void chimeraHuntStoreOverlapCoverage(
    const int8 direction,
    const uint32 actreadid,
    const uint32 rid2,
    uint16 hp1min,
    uint16 hp1max,
    uint16 hp2min,
    uint16 hp2max);
  void chimeraHuntLocateChimeras();

  void makeVHRAPArrayShortcuts(std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> & vhraparray,
			       const uint32 basesperhash);


/*************************************************************************
 *
 * sorters VHRAPArray
 *
 *************************************************************************/

  // standard sort: from low to high, but lower 24bit grouped
  inline static bool sortVHRAPArrayElem_(const typename HashStatistics<TVHASH_T>::vhrap_t & a, const typename HashStatistics<TVHASH_T>::vhrap_t & b) {
    if((a.vhash & SKIM3_MAXVHASHMASK) != (b.vhash & SKIM3_MAXVHASHMASK)) {
      return (a.vhash & SKIM3_MAXVHASHMASK) < (b.vhash & SKIM3_MAXVHASHMASK);
    }
    return a.vhash < b.vhash;
  }

  // sorter to sort from low to high

  inline static bool compareVHRAPArrayElem_(const typename HashStatistics<TVHASH_T>::vhrap_t & one, const typename HashStatistics<TVHASH_T>::vhrap_t & other)
    {
      return one.vhash < other.vhash;
    };


/*************************************************************************
 *
 * sorters readhashmatches_t
 *
 *************************************************************************/

  // low to high
  inline static bool sortreadhashmatch_t_(const readhashmatch_t & a, const readhashmatch_t & b)
    {
      if(a.rid2 == b.rid2){
	if(a.eoffset == b.eoffset) return a.hashpos1 < b.hashpos1;
	return a.eoffset < b.eoffset;
      }
      return a.rid2 < b.rid2;
    }

/*************************************************************************
 *
 * sorters matchwithsorter_t
 *
 *************************************************************************/

  // high to low
  inline static bool sortMWByPercent_(const matchwithsorter_t & a, const matchwithsorter_t & b)
    {
      if(a.percent_in_overlap == b.percent_in_overlap) {
	return a.numhashes > b.numhashes;
      }
      return a.percent_in_overlap > b.percent_in_overlap;
    }

  // high to low
  inline static bool sortMWByEstimScore_(const matchwithsorter_t & a, const matchwithsorter_t & b)
    {
      if(a.estimscore == b.estimscore) {
	return a.percent_in_overlap > b.percent_in_overlap;
      }
      return a.estimscore > b.estimscore;
    }

  // high to low
  bool sortMWByNumHashes_(const matchwithsorter_t & a, const matchwithsorter_t & b)
    {
      if(a.numhashes == b.numhashes){
	return a.percent_in_overlap > b.percent_in_overlap;
      }
      return a.numhashes > b.numhashes;
    }




//////////////////////////////////

  void init();

  uint32 computePartition(uint32 maxmemusage,
			  bool computenumpartitions);

  void sFR_makeHashCounts(std::vector<uint32> & hashcounter, uint32 basesperhash);

////////////////////////////////////////

  void correctReadBaseStatisticsByRMB(ReadPool & rp,
				      const uint32 basesperhash);


public:
  Skim();
  ~Skim();

  Skim(Skim const &other) = delete;
  Skim const & operator=(Skim const & other) = delete;

  void discard();


  uint32 skimGo (ReadPool & rp,
		 std::string               & posfmatchname,
		 std::string               & poscmatchname,
		 bannedoverlappairs_t & bannedoverlaps,
		 std::vector<uint32>       & overlapcounter,
		 std::vector<uint32>       & writtenhitsperid,
		 std::vector<int32>        & chuntleftcut,
		 std::vector<int32>        & chuntrightcut,
		 std::vector<std::vector<uint8>>        & overlapcritlevell,
		 std::vector<std::vector<uint8>>        & overlapcritlevelr,
		 std::vector<uint8>        * meghubs,

		 uint32 numthreads,           //2
		 uint32 maxmemusage,          //  = 15000000

		 bool onlyagainstrails,       // =false
		 bool alsocheckreverse,       // =true

		 uint32 basesperhash,         // = 16
		 uint8  hashsavestepping,     // = 4
		 const std::vector<int32> & percentrequired,   // = 30 for each st
		 const std::vector<int32> & overlaplenrequired,   // = 30 for each st
		 uint32 maxhitsperread,
		 uint32 megahubcap,
		 bool forcetakestronggood);      // = 200

  void skimStreamPrepare(ReadPool & rp, uint32 basesperhash, uint8  hss, const char * additionalregexp=nullptr);


  static void getOverlapCriterionLevel(const uint32 actreadid,
				       const uint8 seqtype,
				       const ADSEstimator & adse,
				       const uint8 relscore,
				       uint8 & levell,
				       uint8 & levelr);

  static uint32 transformSeqToVariableHash(
    const uint32 readid,
    const Read & actread,
    const char * seq,
    uint32 slen,
    const uint32 basesperhash,
    typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::iterator & vhraparrayI,
    const bool countonly,
    const uint8 hashsavestepping,
    std::vector<uint8> & maskvector,
    const std::vector<Read::bposhashstat_t> & bposhashstats,
    int32 bfpos,
    const int32 bfposinc
    );

  void setExtendedLog(bool f) {
    SKIM3_logflag_purgeunnecessaryhits=f;
    SKIM3_logflag_save2=f;
  }

};


#endif
