/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 and later by Bastien Chevreux
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
 *
 */


#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include <boost/filesystem.hpp>
#include "boost/format.hpp"

#include "errorhandling/errorhandling.H"

#include "util/machineinfo.H"
#include "util/dptools.H"
#include "util/fileanddisk.H"
#include "util/progressindic.H"

#include "mira/hashstats.H"

#include "mira/bloomfilter.H"
#include "mira/skim.H"
#include "mira/seqtohash.H"
#include "mira/readgrouplib.H"



using namespace std;
using boost::format;


//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif

#ifndef PUBLICQUIET
#define CLOCKSTEPS
#endif

#ifdef CLOCKSTEPS
#define TEBUG(bla)   {cout << bla; cout.flush();}
#else
#define TEBUG(bla)
#endif

//#define CEBUG(bla)   {cout << bla; cout.flush();}


// for timing a couple of things we need a "global" variable
// let's cheat and not put that into the class

#ifdef CLOCKSTEPS
timeval HS_CHEAT_tvfill;
#endif



#ifdef HSVHM_var
template<typename TVHASH_T>
const TVHASH_T NHashStatistics<TVHASH_T>::HS_MAXVHASHMASK(0xffffffUL);
#endif


template<typename TVHASH_T>
uint32 NHashStatistics<TVHASH_T>::HSN_hs_magic=0x4D4C6873;  // magic: "MLhs" MiraLibHashStat




/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
NHashStatistics<TVHASH_T>::~NHashStatistics()
{
  deleteBloomFilter();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::deleteBloomFilter()
{
  if(HSN_bloomfilter!=nullptr){
    delete HSN_bloomfilter;
    HSN_bloomfilter=nullptr;
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::analyseReadPool(ReadPool & rp)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::analsyeReadPool(ReadPool & rp)");

  ProgressIndicator<int64> pi(0,rp.size());

  for(uint32 ri=0; ri<rp.size(); ++ri){
    rp[ri].getClippedSeqAsChar();
    rp[ri].getClippedComplementSeqAsChar();
  }

  dateStamp(cout);

  for(uint32 step=1; step<=HSN_needsteps; ++step){
    pi.reset(0,rp.size());
    for(uint32 ri=0; ri<rp.size(); ++ri){
      pi.progress(ri);
      learnSequence(rp[ri].getClippedSeqAsChar(),
		    rp[ri].getLenClippedSeq(),
		    rp[ri].getName().c_str(),
		    0,
		    false);
      learnSequence(rp[ri].getClippedComplementSeqAsChar(),
		    rp[ri].getLenClippedSeq(),
		    rp[ri].getName().c_str(),
		    0,
		    true);
    }
    pi.finishAtOnce();
    cout << endl;
    finaliseStep();
    dateStamp(cout);
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::setupNewAnalysis(const uint8 bfbits, const uint32 bfnumkeys, const uint32 basesperhash, uint16 numsteps)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::setupNewAnalysis(const uint8 bfbits, const uint32 bfnumkeys,  const uint32 basesperhash, uint16 numsteps");
  BUGIFTHROW(HSN_bloomfilter!=nullptr,"HSN_bloomfilter!=nullptr ??");

  HSN_bloomfilter=new BloomFilter<TVHASH_T>(bfbits,bfnumkeys);
  HSN_basesperhash=basesperhash;

  BUGIFTHROW(numsteps==0,"numsteps==0 ???");
  BUGIFTHROW(numsteps>3,"numsteps " << numsteps << " not in 1,2,3.");

  HSN_needsteps=numsteps;

  if(numsteps==1){
    HSN_step=1001;
    cout << "Counting hashes (quick, slightly inaccurate, 1 pass): step 1" << endl;
  }else if(numsteps==2){
    HSN_step=2001;
    cout << "Counting hashes (quick, accurate, 2 pass): step 1" << endl;
  }else{
    HSN_step=3001;
    cout << "Counting hashes (accurate, savemem, 3 pass): step 1" << endl;
  }

  HSN_hs_sortstatus=HashStatistics<TVHASH_T>::HSSS_NOTSORTED;
  HSN_hs_needsconsolidation=false;

}

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::learnSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::learnSequenceStep(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  BUGIFTHROW(HSN_bloomfilter==nullptr,"HSN_bloomfilter==nullptr ???");

  switch(HSN_step){
  case 1001 : {
    learnSequenceQuick1(seqvoid,slen,namestr,seqtype,isreverse,false);
    break;
  }
  case 2001 : {
    learnSequenceStep1(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 2002 : {
    learnSequenceQuick1(seqvoid,slen,namestr,seqtype,isreverse,true);
    break;
  }
  case 3001 : {
    learnSequenceStep1(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 3002 : {
    learnSequenceStep2(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 3003 : {
    learnSequenceStep3(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 32678 : {
    BUGIFTHROW(true,"HSN_step 32678, nothing more to learn???");
    break;
  }
  default :{
    BUGIFTHROW(true,"unknown HSN_step "  << static_cast<int16>(HSN_step));
  }
  }
}


template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::finaliseStep()
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::finaliseStep1()");

  if(HSN_step==1){
    cout << "Counting hashes: finalised step 1, switching to step 2" << endl;
    HSN_hsv_hashstats.reserve(HSN_bloomfilter->getNumKMersSeenGE2());
    HSN_bloomfilter->reset();
  }else if(HSN_step==2){
    cout << "Counting hashes: finalising step 2 ..."; cout.flush();
    makeNHashStatArrayShortcuts(HSN_hsv_hashstats, HSN_basesperhash, HSN_hsv_hsshortcuts);
    cout << " done.\nCounting hashes: step 3" << endl;
  }else if(HSN_step==3){
    cout << "Trimming out erroneous single hashes ..."; cout.flush();
    auto oldsize=HSN_hsv_hashstats.size();
    trimHashStatsByFrequency(-1,-1,2);
    cout << " done. Trimmed " << oldsize-HSN_hsv_hashstats.size() << " hashes, " << HSN_hsv_hashstats.size() << " remaining" << endl;
    HSN_step=32767;
  }else if(HSN_step==1001){
    cout << "quick done\n";
    HSN_step=32767;
  }else if(HSN_step==2001){
    cout << "quick 2.1 done\n";
  }else if(HSN_step==2002){
    cout << "quick 2.2 done\n";
    HSN_step=32767;
  }else{
    BUGIFTHROW(true,"HSN_step is " << HSN_step << " ???");
  }

  ++HSN_step;
}


template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::learnSequenceQuick1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse, bool lookuponly)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::learnSequenceQuick(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  static hscounts_t tmphs;

  BUGIFTHROW(HSN_step!=1001 && HSN_step!=2002,"HSN_step!=1001 && HSN_step!=2002 ???");

  uint32 countinit=1;
  if(lookuponly) ++countinit;

  auto basesperhash=HSN_basesperhash;
  const uint8 * seq=static_cast<const uint8 *>(seqvoid);
  int bfres=0;
  SEQTOHASH_LOOPSTART(TVHASH_T){
    auto thispos=seqi-basesperhash;
    if(unlikely(isreverse)){
      thispos=slen-1-seqi;
    }
    auto umhsI=HSN_hsum_hashstats.find(acthash);
    if(umhsI!=HSN_hsum_hashstats.end()){
      if(thispos < umhsI->second.getLowPos()){
	umhsI->second.setLowPos(thispos);
      }
      if(umhsI->second.seqtype!=seqtype){
	seqtype=0xf;
      }
      if(unlikely(isreverse)){
	if(unlikely(++(umhsI->second.rcount)==0)) --(umhsI->second.rcount);
      }else{
	if(unlikely(++(umhsI->second.fcount)==0)) --(umhsI->second.fcount);
      }
    }else{
      if(lookuponly){
	bfres=HSN_bloomfilter->isNonUnique(acthash);
      }else{
	bfres=HSN_bloomfilter->addVHash(acthash);
      }
      if(bfres==1){
	// make new in unordered map!
	tmphs.setLowPos(thispos);
	tmphs.seqtype=seqtype;
	if(isreverse){
	  tmphs.fcount=0;
	  tmphs.rcount=countinit;
	}else{
	  tmphs.fcount=countinit;
	  tmphs.rcount=0;
	}
	HSN_hsum_hashstats[acthash]=tmphs;
      }
      // hmmm ... should not happen here
      BUGIFTHROW(bfres==2,"bfres==2 ???");
    }
  }SEQTOHASH_LOOPEND;
}


template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::learnSequenceQuick2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
}


template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::learnSequenceStep1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::learnSequenceStep1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  BUGIFTHROW(HSN_step!=1 && HSN_step!=2001,"HSN_step!=1 && HSN_step!=2001 ???");
//  HSN_bloomfilter->addSequenceToBloomfield(seqvoid, slen, HSN_basesperhash, namestr);

  auto basesperhash=HSN_basesperhash;
  for(uint32 xxi=0; xxi<2;++xxi){
    const uint8 * seq=static_cast<const uint8 *>(seqvoid);
    SEQTOHASH_LOOPSTART(TVHASH_T);
    if(xxi){
      (void) HSN_bloomfilter->addVHash(acthash);
    }else{
      HSN_bloomfilter->prefetchVHash(acthash);
    }
    SEQTOHASH_LOOPEND;
  }
}

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::learnSequenceStep2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::learnSequenceStep2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  static nhashstat_t tmphs;

  BUGIFTHROW(HSN_step!=2,"HSN_step!=2 ???");

  auto basesperhash=HSN_basesperhash;
  for(uint32 xxi=0; xxi<2;++xxi){
    const uint8 * seq=static_cast<const uint8 *>(seqvoid);
    SEQTOHASH_LOOPSTART(TVHASH_T);
    if(xxi){
      if(HSN_bloomfilter->addVHash(acthash)==1){
	BUGIFTHROW(HSN_hsv_hashstats.size()==HSN_hsv_hashstats.capacity(),"HSN_hsv_hashstats.size()==hstable.capacity() ???");
	tmphs.vhash=acthash;
	HSN_hsv_hashstats.push_back(tmphs);
      }
    }else{
      HSN_bloomfilter->prefetchVHash(acthash);
    }
    SEQTOHASH_LOOPEND;
  }
}


template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::learnSequenceStep3(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::learnSequenceStep3(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse))");

  BUGIFTHROW(HSN_step!=3,"HSN_step!=3 ???");

  nhashstat_t tmphs;

  // we expect to have much more kmers occurring more than once than single kmers
  // therefore, do not use a query to the bloomfilter before using find() on the
  //  vector<hashstat_t>, it's just a waste of time

  auto basesperhash=HSN_basesperhash;
  const uint8 * seq=static_cast<const uint8 *>(seqvoid);
  SEQTOHASH_LOOPSTART(TVHASH_T){
    tmphs.vhash=acthash;
    auto hptr=const_cast<nhashstat_t *>(findVHash(tmphs));
    if(likely(hptr!=nullptr)){
      auto thispos=seqi-basesperhash;
      if(unlikely(isreverse)){
	thispos=slen-1-seqi;
      }
      if(unlikely(hptr->hsc.fcount==0 && hptr->hsc.rcount==0)){
	hptr->hsc.seqtype=seqtype;
	hptr->hsc.setLowPos(thispos);
      }else if(thispos < hptr->hsc.getLowPos()){
	hptr->hsc.setLowPos(thispos);
      }

      if(unlikely(isreverse)){
	if(unlikely(++(hptr->hsc.rcount)==0)) --(hptr->hsc.rcount);
      }else{
	if(unlikely(++(hptr->hsc.fcount)==0)) --(hptr->hsc.fcount);
      }
    }
  }SEQTOHASH_LOOPEND;

}

#define prefetchrl(p)     __builtin_prefetch((p), 0, 3)

template<typename TVHASH_T>
const typename NHashStatistics<TVHASH_T>::nhashstat_t * NHashStatistics<TVHASH_T>::findVHash(const nhashstat_t & searchval)
{
  const nhashstat_t * ret=nullptr;

  // even if executed a couple of million times, this if takes virtually no time at all
  // so keep it
  if(unlikely(HSN_hsv_hsshortcuts.empty())) {
    makeNHashStatArrayShortcuts(HSN_hsv_hashstats, HSN_basesperhash, HSN_hsv_hsshortcuts);
  }

  auto hsindex=static_cast<size_t>(searchval.vhash & HS_MAXVHASHMASK);
  if(likely(!HSN_hsv_hashstats.empty())
     && HSN_hsv_hashstats.end() != HSN_hsv_hsshortcuts[hsindex].b){
    auto hsI=HSN_hsv_hsshortcuts[hsindex].b;
    // TODO: test with large & diverse data set effect of prefetch
    prefetchrl(&(*hsI));
    if(HSN_hsv_hsshortcuts[hsindex].e-hsI > 1){
      // with more than 12 bases in a hash, the array is subdivided
      // TODO: test with large & diverse data set whether this split in lower_bound
      //  vs. simple while loop is OK
      if(HSN_hsv_hsshortcuts[hsindex].e-hsI > 4){
	hsI=lower_bound(hsI,
			HSN_hsv_hsshortcuts[hsindex].e, // upperbound
			searchval,
			sortHashStatComparator);
      }else{
	while(hsI!=HSN_hsv_hsshortcuts[hsindex].e && hsI->vhash!=searchval.vhash){
	  ++hsI;
	}
      }
    }
    if(hsI != HSN_hsv_hashstats.end()
       && hsI->vhash == searchval.vhash) ret=&(*hsI);
  }

  return ret;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  if(HSN_hsum_hashstats.empty()){
    trimHashVStatsByFrequency(minfwd,minrev,mintotal);
  }else{
    trimHashMStatsByFrequency(minfwd,minrev,mintotal);
  }

}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::trimHashVStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::trimHashVStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  auto srcI=HSN_hsv_hashstats.begin();
  auto dstI=srcI;

  for(; srcI!=HSN_hsv_hashstats.end(); ++srcI){
    bool ok=true;
    if((minfwd>=0 && srcI->hsc.fcount<minfwd)
       || (minrev>=0 && srcI->hsc.rcount<minrev)
       || (mintotal>=0 && srcI->hsc.fcount+srcI->hsc.rcount < mintotal)){
      ok=false;
      CEBUG("rm\t");
    }else{
      CEBUG("keep\t");
    }
    CEBUG(srcI-HSN_hsv_hashstats.begin() << "\t" << *srcI << endl);
    *dstI=*srcI;
    if(ok)++dstI;
  }
  HSN_hsv_hashstats.resize(dstI-HSN_hsv_hashstats.begin());
  HSN_hsv_hsshortcuts.clear();
  HSN_hs_dist.clear();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::trimHashMStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::trimHashMStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  auto srcI=HSN_hsum_hashstats.begin();

  for(; srcI!=HSN_hsum_hashstats.end();){
    bool ok=true;
    if((minfwd>=0 && srcI->second.fcount<minfwd)
       || (minrev>=0 && srcI->second.rcount<minrev)
       || (mintotal>=0 && srcI->second.fcount+srcI->second.rcount < mintotal)){
      ok=false;
      CEBUG("rm\t");
    }else{
      CEBUG("keep\t");
    }
// no compile on gcc 4.6.1
//    CEBUG(srcI-HSN_hsum_hashstats.begin() << "\t" << *srcI << endl);
    if(ok) {
      ++srcI;
    }else{
      srcI=HSN_hsum_hashstats.erase(srcI);
    }
  }
}


/*************************************************************************
 *
 * needs:
 *  - hashstats filled with entries (can be unsorted, will be re-sorted
 *    anyway)
 *
 * returns:
 *  - hashstats array sorted by low 24 bit (low to high), then by vhash
 *  - elements .b and .e in hsshortcuts pointing to start and end of each
 *    low 24 bit group of same value
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::makeNHashStatArrayShortcuts(vector<nhashstat_t> & hashstats, const uint32 basesperhash, vector<hsvbendit_t> & hsshortcuts)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::makeNHashStatArrayShortcuts(vector<hashstat_t> & hashstats, const uint32 basesperhash, vector<hsvbendit_t> & hsshortcuts)");

  CEBUG("makeNHashStatArrayShortcuts: basesperhash: " << basesperhash << "\n");

  BUGIFTHROW(basesperhash==0, "basesperhash == 0 ???");

  for(size_t hsi=0; hsi<hashstats.size(); ++hsi){
    CEBUG(hashstats[hsi] << '\n');
  }

  sortLow24Bit(hashstats,HSN_hs_sortstatus);

  hsshortcuts.clear();
  auto hsI=hashstats.cbegin();

  {
    hsvbendit_t tmpb;
    tmpb.b=hashstats.end();
    tmpb.e=hashstats.end();
    hsshortcuts.resize(
      1<<(min(static_cast<uint32>(12),basesperhash)*2),
      tmpb
      );
  }

  if(hsI==hashstats.end()) return;

  CEBUG("hsshortcuts.size(): " << hsshortcuts.size() << endl);

  TVHASH_T acthash= (hsI->vhash & HS_MAXVHASHMASK);
  while(hsI != hashstats.end()){
    CEBUG("begin " << hex << acthash << dec << " is: " << *hsI << endl);
    hsshortcuts[static_cast<size_t>(acthash)].b=hsI;
    for(;(hsI != hashstats.end()) && ((hsI->vhash & HS_MAXVHASHMASK) == acthash); hsI++) {
      CEBUG("INC\n")
    }
    CEBUG("end " << hex << acthash << dec << " is: " << *hsI << endl);
    hsshortcuts[static_cast<size_t>(acthash)].e=hsI;
    //cout << "vhash: " << hex << acthash << "\t" << dec << hsshortcuts_end[acthash]-hsshortcuts_begin[acthash] << '\n';
    if(hsI != hashstats.end()) acthash= hsI->vhash & HS_MAXVHASHMASK;
  }

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::dumpHealth(ostream & fout)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::dumpHealth(ostream & fout)");
  BUGIFTHROW(HSN_bloomfilter==nullptr,"HSN_bloomfilter==nullptr ???");
  fout << *HSN_bloomfilter;
}



/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::moveStatCountMapToVector()
{
  if(HSN_hsum_hashstats.empty()) return;

  nhashstat_t tmphs;

  HSN_hsv_hashstats.clear();
  HSN_hsv_hashstats.reserve(HSN_hsum_hashstats.size());
  for(auto & hsme : HSN_hsum_hashstats){
    tmphs.vhash=hsme.first;
    tmphs.hsc=hsme.second;
    HSN_hsv_hashstats.push_back(tmphs);
  }
  HSN_hsum_hashstats.clear();
  makeNHashStatArrayShortcuts(HSN_hsv_hashstats, HSN_basesperhash, HSN_hsv_hsshortcuts);
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
uint64 NHashStatistics<TVHASH_T>::calcHashDistrib(vector<uint64> & hsdist)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::calcHashDistrib(vector<uint64> & hsdist)");

  hsdist.clear();

  if(HSN_hsv_hashstats.empty()) moveStatCountMapToVector();

  uint64 maxc=0;
  uint64 totalc=0;
  for(auto & hsve : HSN_hsv_hashstats){
    totalc+=static_cast<uint64>(hsve.hsc.fcount+hsve.hsc.rcount);
    maxc=max(maxc,static_cast<uint64>(hsve.hsc.fcount+hsve.hsc.rcount));
  }
  hsdist.resize(maxc+1,0);
  for(auto & hse : HSN_hsv_hashstats){
    ++hsdist[hse.hsc.fcount+hse.hsc.rcount];
  }

  return totalc;
}


/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
vector<uint64> & NHashStatistics<TVHASH_T>::getHashDistrib()
{
  if(HSN_hs_dist.empty()) calcHashDistrib(HSN_hs_dist);
  return HSN_hs_dist;
}


/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::dumpHashDistrib(ostream & ostr)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::dumpHashDistrib(vector<uint64> & hsdist, ostream & ostr) const");

  auto & hsdist=getHashDistrib();

  uint64 totalc=0;
  uint64 hsdi=0;
  for(auto hsde : hsdist){
    totalc+=hsde*hsdi;
    ++hsdi;
  }

  double dtotalc=static_cast<double>(totalc);
  uint64 cumc=0;
  hsdi=0;
  for(auto hsde : hsdist){
    cumc+=hsde*hsdi;
    double frac=static_cast<double>(cumc)/dtotalc;
    ostr << hsdi
	 << '\t' << hsde
	 << '\t' << hsde*hsdi
	 << '\t' << cumc
	 << '\t' << frac
	 << '\n';
    ++hsdi;
  }
}


/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::saveHashStatistics(const string & filename, bool deleteoldfile)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::saveHashStatistics(const string & filename, bool deleteoldfile)");

  ofstream fout;
  openFileForAppend(filename,fout,deleteoldfile);
  try{
    if(!fout){
      MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is the disk full? Are permissions set right?");
    }
    saveHashStatistics(fout);
  }
  catch(Notify n){
    cout << "Error for file " << filename << endl;
    n.handleError(THISFUNC);
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::loadHashStatistics(const string & filename)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::loadHashStatistics(const string & filename)");

  ifstream fin;
  try{
    fin.open(filename.c_str(),ios::in);
    if(!fin){
      MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is it present? Are permissions set right?");
    }
    loadHashStatistics(fin);
  }
  catch(Notify n){
    cout << "Error while loading file " << filename << endl;
    n.handleError(THISFUNC);
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::saveHashStatistics(ostream & ostr)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::saveHashStatistics(ostream & ostr)");

  if(HSN_hsum_hashstats.empty()){
    saveHashVStatistics(ostr);
  }else{
    saveHashMStatistics(ostr);
  }

}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::saveHashVStatistics(ostream & ostr)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::saveHashVStatistics(ostream & ostr)");

  HashStatistics<TVHASH_T>::priv_writeHashStatFileHeader(ostr,HSN_basesperhash,HSN_hs_sortstatus,HSN_hsv_hashstats.size());
  if(!HSN_hsv_hashstats.empty()){
    ostr.write(reinterpret_cast<const char *>(&HSN_hsv_hashstats[0]),
	       sizeof(nhashstat_t)*HSN_hsv_hashstats.size());
  }
  if(ostr.bad()){
    MIRANOTIFY(Notify::FATAL, "Could not save anymore the hash statistics (1). Disk full? Changed permissions?");
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::saveHashMStatistics(ostream & ostr)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::saveHashMStatistics(ostream & ostr)");

  HashStatistics<TVHASH_T>::priv_writeHashStatFileHeader(ostr,HSN_basesperhash,HashStatistics<TVHASH_T>::HSSS_NOTSORTED,HSN_hsum_hashstats.size());

  nhashstat_t tmphs;
  for(auto & hsume : HSN_hsum_hashstats){
    tmphs.vhash=hsume.first;
    tmphs.hsc=hsume.second;
    ostr.write(reinterpret_cast<const char *>(&tmphs), sizeof(nhashstat_t));
  }
  if(ostr.bad()){
    MIRANOTIFY(Notify::FATAL, "Could not save anymore the hash statistics. Disk full? Changed permissions?");
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::loadHashStatistics(istream & istr)
{
  FUNCSTART("bool NHashStatistics<TVHASH_T>::loadHashStatistics(istream & istr)");

  auto localmagic=HSN_hs_magic;
  istr.read(reinterpret_cast<char *>(&localmagic),4);
  if(istr.gcount()!=4
     || localmagic!=HSN_hs_magic) {
    MIRANOTIFY(Notify::FATAL,"No magic found?\n");
  }
  uint8 tmpbyte=0;
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);
  if(tmpbyte!=2) {
    MIRANOTIFY(Notify::FATAL,"Not version 2?\n");
  }
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);
  if(tmpbyte==0 || tmpbyte>32){
    MIRANOTIFY(Notify::FATAL,"Invalid kmer size " << static_cast<uint16>(tmpbyte) << " ???");
  }
  if(HSN_basesperhash!=0 && !HSN_hsv_hashstats.empty() && tmpbyte!= HSN_basesperhash){
    MIRANOTIFY(Notify::FATAL,"Current hashstat kmer size is " << HSN_basesperhash
	       << ", but kmer size in data to load is " << static_cast<uint16>(tmpbyte)
	       << " ???\n");
   }else{
    HSN_basesperhash=tmpbyte;
  }
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);
  if(!HSN_hsv_hashstats.empty()){
    HSN_hs_sortstatus=HashStatistics<TVHASH_T>::HSSS_NOTSORTED;
    HSN_hs_needsconsolidation=true;
    MIRANOTIFY(Notify::FATAL,"Appending to existing hashstat not implemented yet\n");
  }else{
    HSN_hs_sortstatus=tmpbyte;
  }
  // padd byte
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);

  uint64 numelem=0;
  istr.read(reinterpret_cast<char *>(&numelem),8);
  if(numelem){
    auto oldsize=HSN_hsv_hashstats.size();
    HSN_hsv_hashstats.resize(HSN_hsv_hashstats.size()+numelem);
    HSN_hsv_hsshortcuts.clear();
    istr.read(reinterpret_cast<char *>(&HSN_hsv_hashstats[oldsize]),numelem*sizeof(nhashstat_t));
    if(istr.gcount()!=numelem*sizeof(nhashstat_t)){
      MIRANOTIFY(Notify::FATAL,"Expected to read " << numelem*sizeof(nhashstat_t) << " bytes, but got " << istr.gcount() << endl);
    }
  }
}


/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void NHashStatistics<TVHASH_T>::dumpHashCount(ostream & ostr)
{
  FUNCSTART("void NHashStatistics<TVHASH_T>::dumpHashCount(ostream & ostr)");

  string tmpstr;
  for(auto & hse : HSN_hsv_hashstats){
    HashStatistics<TVHASH_T>::hash2string(hse.vhash,HSN_basesperhash,tmpstr);
    cout << tmpstr
	 << '\t' << hse.hsc.fcount
	 << '\t' << hse.hsc.rcount
	 << '\t' << hse.hsc.fcount + hse.hsc.rcount
	 << '\n';
  }
}


// explicit template instantiations needed for the linker to create these for library files
template class HashStatistics<vhash64_t>;
#ifndef KMER_INTERNALTYPE
template class HashStatistics<vhash128_t>;
template class HashStatistics<vhash256_t>;
template class HashStatistics<vhash512_t>;
#endif

template class NHashStatistics<vhash64_t>;
#ifndef KMER_INTERNALTYPE
template class NHashStatistics<vhash128_t>;
template class NHashStatistics<vhash256_t>;
template class NHashStatistics<vhash512_t>;
#endif
