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

#include <regex>

#include "stdinc/defines.H"

#include "mira/readpool.H"

#include <zlib.h>



template<typename TVHASH_T>
class HashStatistics
{

// GRIZZLYBEAR down below to perform timing experiments specific to a
// comparator used in std::sort()
// Result: I want to use HS_MAXVHASHMASK as a const static variable :-)
#if 1
#define HSVHM_var
  const static TVHASH_T HS_MAXVHASHMASK;
//#define GRIZZLYBEAR 0xffffffULL
#define GRIZZLYBEAR HS_MAXVHASHMASK
#else
#define HS_MAXVHASHMASK 0xffffffUL
#define GRIZZLYBEAR HS_MAXVHASHMASK
#endif


/*************************************************************************
 *
 *  hashstatistics public structs and vars
 *
 *************************************************************************/

public:
  struct vhrap_t {
    TVHASH_T vhash;       // vhash,
    uint32 readid;       // readid
    uint16 hashpos;      // (hash)position  (lowest pos: basesperhash-1)

    Read::bhashstat_t bhashstats; // baseflags for this hash

    //bool operator<(const vhrap_t & other) const {return vhash < other.vhash;};
    friend std::ostream & operator<<(std::ostream &ostr, const vhrap_t & v){
      ostr << "hash: " << std::hex << v.vhash
	   << "\trid: " <<std::dec << v.readid
	   << "\thpos: " << v.hashpos
	   << "\tbhs: " << v.bhashstats;
      return ostr;
    }
  };

  enum {HSSS_NOTSORTED=0, HSSS_LOW24BIT, HSSS_LEXIUP, HSSS_BYCOUNTUP, HSSS_BYCOUNTDOWN, HSSS_LEXBYCOUNT, HSSS_MASKUP};

  static void hash2string(TVHASH_T hash, uint32 basesperhash, std::string & str);
  static inline std::string hash2string(TVHASH_T hash, uint32 basesperhash) { std::string tmp; hash2string(hash,basesperhash,tmp); return tmp;}


#define MULTISEQTYPE 0xf

  // new hashstat
  class hscounts_t {
  private:
    uint8   lowposd4;      // lowest hash position/4 (lowest pos: 0, 255 means >=1020)
  public:
    uint32  fcount:24;
    bool    iskmerforkf:1;
    bool    iskmerforkr:1;
    bool    hasfwdrevthresholdok:1; // tempted, temporary solution, delete afterwards
    bool    b4:1;          // 1 bit, unused atm
  public:
    uint8   seqtype:4;     // sequencing type
    uint32  rcount:24;

  public:
    hscounts_t() {
      uint64 * ptr=reinterpret_cast<uint64 *>(&lowposd4);
      *ptr=0;
    }

    friend std::ostream & operator<<(std::ostream &ostr, const hscounts_t & v){
      ostr << "fcount: " << v.fcount
	   << "\trcount: " << v.rcount
	   << "\tlowposd4: " << static_cast<uint32>(v.lowposd4)
	   << "\tseqtype: " << static_cast<uint16>(v.seqtype)
	   << "\tkmf: " << v.iskmerforkf
	   << " kmr: " << v.iskmerforkr
	   << " fr: " << v.hasfwdrevthresholdok
	;
      return ostr;
    }

    inline void setLowPos(uint32 pos) {lowposd4=static_cast<uint8>(std::min(pos/4,static_cast<uint32>(255)));}
    inline uint32 getLowPos() const {return lowposd4*4;}
    inline uint32 getCount() const {return fcount+rcount;}
  };


  struct hashstat_t {
    hscounts_t hsc;
    TVHASH_T   vhash;       // vhash,

    hashstat_t() {
      vhash=0;
    }

    friend std::ostream & operator<<(std::ostream &ostr, const hashstat_t & v){
      ostr << "hash: " << std::hex << v.vhash << std::dec
	   << '\t' << v.hsc
	;
      return ostr;
    }
  };

  // Some avg frequency numbers for the hash statistics
  struct avg_freq_t {
    bool valid:1;
    uint16 padding:15;

    uint32 corrected;
    uint32 raw;
    uint32 taken;
    uint32 min;

    avg_freq_t() : valid(false), padding(0), corrected(0), raw(0), taken(0), min(0) {}
  };

  struct mhsheader_t {
    uint8  version=0;
    uint8  sortstatus=0;
    uint32 basesperhash=0;
    uint32 sizeofhash=0;      // in bytes, as given by sizeof(TVHASH_T)
    uint64 numelem=0;

    avg_freq_t freq;

    friend std::ostream & operator<<(std::ostream &ostr, const mhsheader_t & mhsh){
      ostr << "v: " << static_cast<uint16>(mhsh.version)
	   << "\tss: " << static_cast<uint16>(mhsh.sortstatus)
	   << "\tbph: " << mhsh.basesperhash
	   << "\tsoh: " << mhsh.sizeofhash
	   << "\tne: " << mhsh.numelem
	   << "\nh.v " << mhsh.freq.valid
	   << "\th.c " << mhsh.freq.corrected
	   << "\th.r " << mhsh.freq.raw
	   << "\th.t " << mhsh.freq.taken
	   << "\th.m " << mhsh.freq.min
	;
      return ostr;
    }
  };


/*************************************************************************
 *
 * hashstatistics_sdbg
 *
 * All of this is still a little bit brittle ...
 *
 *************************************************************************/

  // -------------------------------------------------------------------------------------
  // Structs and vars

public:
  struct dbgedits_t {
    uint32 pos=0;
    uint32 len=0;
    std::string replacement;

    dbgedits_t(uint32 p, uint32 l, std::string & r) : pos(p), len(l), replacement(r) {};
  };

  struct dbgseq_t {
    std::string * seqstrptr=nullptr;  // non-owning pointer into HS_hsv_seqcontainer elements. Allows fast sorting.
    uint32 medianhashcov=0;

    uint32 expectedcopynumber=0;

    int64 hsni_first=-1;  // hashstatnode_t index first (first kmer of sequence)
    int64 hsni_last=-1;  // hashstatnode_t index last (last kmer of sequence)
    uint32 graphid=-1;
    //dbgseq_t() : medianhashcov(0), expectedcopynumber(0), hsni_first(-1), hsni_last(-1), graphid(-1) {};
  };

private:
  struct hashstatnode_t {
    int64 next;
    int64 prev;

    uint32 graphid;

    hashstatnode_t() : next(-1), prev(-1), graphid(-1) {};

    friend std::ostream & operator<<(std::ostream &ostr, const hashstatnode_t & v){
      ostr << v.next
	   << '\t' << v.prev
	   << '\t' << v.graphid
	;
      return ostr;
    }
  };

  // for sdg
  std::vector<hashstatnode_t> HS_hsv_hashstatnodes;  // size of HS_hsv_hashstats
  std::vector<dbgseq_t> HS_hsv_dbgseqs;  //
  std::deque<std::string> HS_hsv_seqcontainer; // deque! we're storing pointers into elements during construction there! And vector might realocate.

  // -------------------------------------------------------------------------------------
  // Functions

private:
  void priv_sdbg_populateHashStatNodes();
  void priv_sdbg_assignGraphIDsMakeSeqsAndCollectStats();
  void priv_sdbg_calcExpectedCopyNumbers();

  void priv_sdbg_buildSequenceFromPath(int64 hsi1, int64 hsi2, std::string & retseq);
  bool priv_sdbg_checkPathDistance(int64 hsi1, int64 hsi2, int64 dist);

  inline static bool sortDBGComparatorLenDown(const dbgseq_t & a, const dbgseq_t & b){
    if(a.seqstrptr->size()==b.seqstrptr->size()){
      return a.medianhashcov>b.medianhashcov;
    }
    return (a.seqstrptr->size() > b.seqstrptr->size());
  }

  inline static bool sortDBGComparatorMHCDown(const dbgseq_t & a, const dbgseq_t & b){
    if(a.medianhashcov>b.medianhashcov){
       return a.seqstrptr->size()==b.seqstrptr->size();
    }
    return a.medianhashcov>b.medianhashcov;
  }


public:
  void buildSDBGraphs();
  void proposeSDBGEditsForSequence(const void * seqvoid, uint64 slen, const char * namestr, std::vector<dbgedits_t> & edits);
  bool checkSequenceForSDBGChimeras(const void * seqvoid, uint64 slen, const char * namestr);

  void sortDBGSeqsByLenDown();
  void sortDBGSeqsByMHCDown();
  void dumpDBGSeqs(std::ostream & ostr);


/*************************************************************************
 *
 * hashstatistics
 *
 *************************************************************************/


private:
  static uint32 HS_hsfilemagic;

  struct hsvbendit_t {
    typename std::vector<hashstat_t>::const_iterator b;
    typename std::vector<hashstat_t>::const_iterator e;
  };

  static size_t HS_numelementsperbuffer;

//  ReadPool * HS_readpoolptr;

  bool HS_abortall;

  std::string HS_tmpdirectorytodelete;

  double HS_freqest_minnormal=.4;
  double HS_freqest_maxnormal=1.6;
  double HS_freqest_repeat=2.0;
  double HS_freqest_heavyrepeat=8.0;
  double HS_freqest_crazyrepeat=20.0;
  uint32 HS_rarekmercount=2;
  uint32 HS_nastyrepeatratio=100;
  uint32 HS_nastyrepeatcoverage=0;

  uint32 HS_hs_basesperhash=0;

  std::vector<hashstat_t> HS_hsv_hashstats;
  std::vector<hsvbendit_t> HS_hsv_hsshortcuts;


  // for streaming old hash statistics
  std::string HS_hashstatfilename;
  std::vector<std::string> HS_hashfilenames;
  std::vector<size_t>      HS_elementsperfile;
  std::vector<gzFile>      HS_hashfiles;
  std::vector<std::vector<hashstat_t> > HS_hashfilebuffer;
  size_t HS_rightshift;

  //
  uint8 HS_hs_sortstatus; // not fully implemented yet!

  avg_freq_t HS_avg_freq;

  // size of ReadGroupLib::getNumSequencingTypes()
  // kmers having count < given value are marked as frequency 1 (or 0 ???) so that
  //  they can be masked later on in data processing
  //
  // assignReadBaseStatistics() will initialise it to {0, 0, ...} if not of correct size
  std::vector<uint32> HS_rarekmermasking;

  /*
    Baiting
   */

  std::vector<vhrap_t> HS_baiting_singlereadvhraparray;
  std::vector<uint8>   HS_baiting_tagmaskvector;

  bool HS_logflag_hashcount;

  /*
    Multithreading assign read base statistics
  */

  struct arbs_threadsharecontrol_t {
    boost::mutex accessmutex;
    int32 from;
    int32 to;
    int32 todo;
    int32 done;
    int32 stepping;

    // need to go via this as the boost:bind does not like a "ReadPool &" as parameter
    // and it also cannot have more than 9 parameters in total ... we'd have more with the below
    ReadPool * rpptr;
    size_t avghashcov;

    std::vector<hashstat_t> * hashstatsptr;
    std::vector<hsvbendit_t> * hsscptr;

    uint32 basesperhash;
    bool masknastyrepeats;
    bool truekmerforks;
  };

  //
  // Digital normalisation
  //

  std::vector<size_t> HS_diginorm_count;

  std::vector<uint8>  HS_diginorm_allow_s1;
  std::vector<size_t> HS_diginorm_vhashindexes_s1;

  //
  // A variable hash mask for the sortHashStatComparatorByMask
  // atm used by priv_calcKMerForks()
  //
  // Careful: due to static nature, cannot be used by multiple threads with different values!!!

  static TVHASH_T HS_vhashmask;

  /*
    Logging, debugging
  */


public:


  //Functions
private:
  // old sort & count kmer counting
  // -------------------------------------------------------------------------------------
  void priv_hashes2disk(ReadPool & rp,
			bool checkusedinassembly,
			bool alsorails,
			bool fwdandrev,
			uint32 basesperhash);

  size_t priv_writeCompressedHFB(size_t hfindex,
				 std::vector<hashstat_t> & hfb,
				 uint32 rarekmerearlykill,
				 gzFile & gzf,
				 bool force);
  void priv_compressHashStatBufferInPlace(std::vector<hashstat_t> & hsb,
					  uint32 rarekmerearlykill);
  void priv_calcFwdRevMinThresholdOKFlag(std::vector<hashstat_t> & hsb, uint32 fwdrevmin);
  size_t priv_createHashStatisticsFile(const std::string & hashstatfilename,
				       std::vector<std::string> & hashfilenames,
				       std::vector<size_t> & elementsperfile,
				       uint32 fwdrevmin,
				       uint32 rarekmerearlykill);
  // -------------------------------------------------------------------------------------
  // HashStatistics comparators
  inline static bool sortHashStatComparatorByLow24bit(const hashstat_t & a, const hashstat_t & b){
    if((a.vhash & GRIZZLYBEAR) != (b.vhash & GRIZZLYBEAR)) {
      return (a.vhash & GRIZZLYBEAR) < (b.vhash & GRIZZLYBEAR);
    }
    //if((GRIZZLYBEAR & a.vhash) != (GRIZZLYBEAR & b.vhash)) {
    //  return (GRIZZLYBEAR & a.vhash) < (GRIZZLYBEAR & b.vhash);
    //}
    return a.vhash < b.vhash;
  }

  inline static bool sortHashStatComparatorLexByCount(const hashstat_t & a, const hashstat_t & b){
    if(a.vhash==b.vhash){
      return a.hsc.fcount+a.hsc.rcount < b.hsc.fcount+b.hsc.rcount;
    }
    return a.vhash < b.vhash;
  }
  inline static bool sortHashStatComparatorLexicographicallyUp(const hashstat_t & a, const hashstat_t & b){
    return a.vhash < b.vhash;
  }
  inline static bool sortHashStatComparatorByCountUp(const hashstat_t & a, const hashstat_t & b){
    return a.hsc.fcount+a.hsc.rcount < b.hsc.fcount+b.hsc.rcount;
  }
  inline static bool sortHashStatComparatorByCountDown(const hashstat_t & a, const hashstat_t & b){
    return a.hsc.fcount+a.hsc.rcount > b.hsc.fcount+b.hsc.rcount;
  }
  inline static bool sortHashStatComparatorByMaskUp(const hashstat_t & a, const hashstat_t & b){
    if((a.vhash & HS_vhashmask) != (b.vhash & HS_vhashmask)) {
      return (a.vhash & HS_vhashmask) < (b.vhash & HS_vhashmask);
    }
    return a.vhash < b.vhash;
  }


  // -------------------------------------------------------------------------------------
  // HashStatistics sorter
  void priv_sorthelper(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr, uint8 finalstatus, const char * compname, bool (*ccallback)(const hashstat_t & a, const hashstat_t & b));

  void priv_sortLow24Bit(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr);
  void priv_sortLexicographicallyUp(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr);
  void priv_sortByCountUp(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr);
  void priv_sortByCountDown(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr);
  void priv_sortLexByCount(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr);
  void priv_sortMaskUp(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr);

  inline void priv_sortLow24Bit() {priv_sortLow24Bit(HS_hsv_hashstats,&HS_hs_sortstatus);}
  inline void priv_sortLexicographicallyUp() {priv_sortLexicographicallyUp(HS_hsv_hashstats,&HS_hs_sortstatus);}
  inline void priv_sortByCountUp() {priv_sortByCountUp(HS_hsv_hashstats,&HS_hs_sortstatus);}
  inline void priv_sortByCountDown() {priv_sortByCountDown(HS_hsv_hashstats,&HS_hs_sortstatus);}
  inline void priv_sortLexByCount() {priv_sortLexByCount(HS_hsv_hashstats,&HS_hs_sortstatus);}
  inline void priv_sortMaskUp() {priv_sortMaskUp(HS_hsv_hashstats,&HS_hs_sortstatus);}

  // -------------------------------------------------------------------------------------
  // ...
  void priv_makeHashStatArrayShortcuts();

  // -------------------------------------------------------------------------------------
  // ...
  void priv_calcAvgHashFreq(bool verbose=false);

  size_t priv_calcMidHashStatIndex(const std::vector<hashstat_t> & hashstats,
				   size_t dontcarepercent);

  void priv_ckmf_helper(TVHASH_T HashStatistics__vhashmask, uint32 mincount);
  void priv_ckmf_relaxed_helper(TVHASH_T HashStatistics__vhashmask, uint32 mincount, bool isfwd);

  // assignReadBaseStatistics()
  void priv_arb_thread(uint32 threadnum, arbs_threadsharecontrol_t * tscptr);
  void priv_arb_DoStuff(
    ReadPool & rp,
    size_t avgcov,
    std::vector<hashstat_t> & hashstats,
    const uint32 basesperhash,
    std::vector<hsvbendit_t> & hsshortcuts,
    bool masknastyrepeats,
    int32 fromid,
    int32 toid,
    bool truekmerforks
    );

  bool priv_dn_TestSingleSeq(Read & actread, std::vector<uint8> & dn_allow, std::vector<size_t> & dn_vhashindexes);

  void priv_saveHashVStatistics(gzFile & gzf);

  void priv_trimHashVStatsByFrequencyAND(uint32 minfwd, uint32 minrev, uint32 mintotal);
  void priv_trimHashVStatsByFrequencyANDOR(uint32 minfwd, uint32 minrev, uint32 mintotal);

  // prepareHashStatistics common part
  void priv_phsCommon(ReadPool * rpptr,
		      uint32 basesperhash,
		      int32 megabasesforbuffer,
		      const std::string & tmpdirectory);

  void priv_rde_helper1_set2zero(const std::vector<TVHASH_T> & kdv, TVHASH_T & andmask);

public:
  // TODO: back to non-static private once HSN has been merged
  static mhsheader_t priv_writeHashStatFileHeader(std::ostream & ostr,
						  uint32 basesperhash,
						  uint8 sortstatus,
						  uint64 numelem);

  mhsheader_t priv_writeHashStatFileHeader(gzFile & gzf,
					   uint32 basesperhash,
					   uint8 sortstatus,
					   uint64 numelem);

  HashStatistics();
  HashStatistics(HashStatistics const &other);
  ~HashStatistics() {discard();}

  HashStatistics const & operator=(HashStatistics const & other) = delete;

  static uint32 byteSizeOfHash(uint32 hashlen);

  size_t size() const { return HS_hsv_hashstats.size(); }

  void discard();
  bool hasStatistics() const { return !HS_hsv_hashstats.empty();}

  size_t getAvgHashFreqCorrected() const { return HS_avg_freq.corrected;};
  size_t getAvgHashFreqRaw() const { return HS_avg_freq.raw;};

  void setAvgHashFreqMinimum(size_t m) { HS_avg_freq.min=m;};

  void setHashFrequencyRatios(double freqest_minnormal,
			      double freqest_maxnormal,
			      double freqest_repeat,
			      double freqest_heavyrepeat,
			      double freqest_crazyrepeat,
			      uint32 rarekmercount,
			      uint32 nastyrepeatratio,
			      uint32 nastyrepeatcoverage);

  // ----------------- One step, complete readpool
  void computeHashStatistics(ReadPool & rp,
			     int32 memtouse,
			     bool checkusedinassembly,
			     bool alsorails,
			     bool fwdandrev,
			     uint32 fwdrevmin,
			     uint32 rarekmerearlykill,
			     uint32 basesperhash,
			     const std::string & hashstatfilename,
			     const std::string & tmpdirectory);
  // ----------------- One step, all files given
  uint64 computeHashStatistics(const std::list<std::string> & seqfiles,
			       int32 memtouse,
			       bool fwdandrev,
			       uint32 fwdrevmin,
			       uint32 rarekmerearlykill,
			       uint32 basesperhash,
			       const std::string & hashstatfilename,
			       const std::string & directory);
  // ----------------- Streaming
  void prepareStreamHashStatistics(uint32 basesperhash,
				   int32 memtouse,
				   const std::string & hashstatfilename,
				   const std::string & tmpdirectory);
  void prepareStreamAddNextSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void prepareStreamFinalise(uint32 fwdrevmin, uint32 rarekmerearlykill);
  // -----------------------------------------------------------------------

  void assignReadBaseStatistics_MultiThread(ReadPool & rp,
					    uint32 numthreads,
					    bool masknastyrepeats,
					    bool calckmerforks,
					    uint32 mincountkmerforks,
					    bool needfwdrev
    );

  //void showHashStatisticsInfo(std::string & filename);
  void showHashStatisticsInfo();

  void calcKMerForks(uint32 mincount, bool needfwdrev);

  void loadHashStatistics(const std::string & filename);
  void loadHashStatistics(gzFile & gzf);
  static const mhsheader_t loadHashStatisticsFileHeader(const std::string & filename);
  static const mhsheader_t loadHashStatisticsFileHeader(gzFile & gzf);
  void loadHashStatistics(const mhsheader_t & mhsh, gzFile & gzf);

  void saveHashStatistics(const std::string & filename);
  void saveHashStatistics(gzFile & gzf);


  uint32 getBasesPerHash() const {return HS_hs_basesperhash;}
  size_t getNumHashEntries() const {return HS_hsv_hashstats.size();}
  std::vector<hashstat_t> & getHashStats() {return HS_hsv_hashstats;}

  void dump(std::ostream & ostr);
  void dumpHSDebug(std::ostream & ostr);
  void dumpAsFASTA(std::ostream & ostr);

  void trimHashStatsByFrequencyAND(uint32 minfwd, uint32 minrev, uint32 mintotal);
  void trimHashStatsByFrequencyANDOR(uint32 minfwd, uint32 minrev, uint32 mintotal);
  void recalcFreqStats() {
    priv_calcAvgHashFreq();
  };

  void removeDustEndsRatio(uint8 replen=0, uint8 perc=67);
  void removeDustEndsFixed(uint8 replen, uint32 numbases);

  // use this one for multi-thread baiting on the same HashStatistics object
  //uint32 checkBaitHit(Read & actread, std::vector<vhrap_t> & baiting_singlereadvhraparray, std::vector<uint8> & baiting_tagmaskvector, bool changeseqcase);
  uint32 checkBaitHit(Read & actread, bool changeseqcase, char mask);

  // use this one for single thread baiting on the same HashStatistics object
  //inline uint32 checkBaitHit(Read & actread, bool changeseqcase) {return checkBaitHit(actread,HS_baiting_singlereadvhraparray,HS_baiting_tagmaskvector, changeseqcase);}

  const hashstat_t * findVHash(const hashstat_t & searchval);

  void digiNormReset() { HS_diginorm_count.clear();}
  bool digiNormTestRead(Read & actread, bool force);
  uint32 estimDigiNormCov(Read & actread);

  // for MiraDiff
  void inAAndB(HashStatistics & a, HashStatistics & b);
  void inANotB(HashStatistics & a, HashStatistics & b);

  inline void sortLow24Bit(){
    priv_sortLow24Bit(HS_hsv_hashstats, &HS_hs_sortstatus);
  }
  inline void sortLexicographicallyUp(){
    priv_sortLexicographicallyUp(HS_hsv_hashstats, &HS_hs_sortstatus);
  }
  inline void sortByCountUp(){
    priv_sortByCountUp(HS_hsv_hashstats, &HS_hs_sortstatus);
  }
  inline void sortByCountDown(){
    priv_sortByCountDown(HS_hsv_hashstats, &HS_hs_sortstatus);
  }

  void abortAll() {HS_abortall=true;}
};




#endif
