/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2011 and later by Bastien Chevreux
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

#include "skim.H"

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include "errorhandling/errorhandling.H"

#include "util/dptools.H"

#include "util/stlimprove.H"


using std::cout;
using std::cerr;
using std::endl;


//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif


// TODO get 2nd line checks done and saved


template<class TVHASH_T>
shash_t Skim<TVHASH_T>::SKIM3_lbphs_hashadd[128];

/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<class TVHASH_T>
void Skim<TVHASH_T>::lowBPHSkim()
{
  FUNCSTART("void Skim<TVHASH_T>::lowBPHSkim()");

//  SKIM_progressindicator= new ProgressIndicator<int64>(0,SKIM3_readpool->size());
//
//  for(uint32 i=0; i<SKIM3_readpool->size();++i){
//    for(uint32 j=i+1; j<SKIM3_readpool->size();++j){
//      SKIM_progressindicator->progress(i);
//    }
//  }
//
//  delete SKIM_progressindicator;

  for(uint32 i=0;i < 128; ++i){
    SKIM3_lbphs_hashadd[i]=0xffff;
  }
  SKIM3_lbphs_hashadd['a']=0;
  SKIM3_lbphs_hashadd['A']=0;
  SKIM3_lbphs_hashadd['c']=1;
  SKIM3_lbphs_hashadd['C']=1;
  SKIM3_lbphs_hashadd['g']=2;
  SKIM3_lbphs_hashadd['G']=2;
  SKIM3_lbphs_hashadd['t']=3;
  SKIM3_lbphs_hashadd['T']=3;

  SKIM3_lbphs_idsperhash.clear();
  SKIM3_lbphs_idsperhash.resize(1<<(SKIM3_basesperhash*2));
  totalphits=0;

  SKIM3_lbphs_numoverlapsperid.clear();
  SKIM3_lbphs_numoverlapsperid.resize(SKIM3_readpool->size(),0);
  SKIM3_lbphs_maxoverlaphashesperid.clear();
  SKIM3_lbphs_maxoverlaphashesperid.resize(SKIM3_readpool->size(),0);

  cout << "Starting lbph skim" << endl;

  BUGIFTHROW(SKIM3_basesperhash>14,"bph is " << SKIM3_basesperhash << ". Using a bph >12 is a pretty bad idea memory wise.");

  for(uint32 readpart=0; readpart < SKIM3_readpool->size();){
    SKIM_partfirstreadid=readpart;
    lbphsPrepareHashOverviewTable(readpart);
    SKIM_partlastreadid=readpart;
    cout << "Prepared " << SKIM_partfirstreadid << " to " << SKIM_partlastreadid << endl;
    startMultiThreading(1,SKIM3_numthreads,1000,SKIM_partfirstreadid,SKIM3_readpool->size(),
    			boost::bind( &Skim<TVHASH_T>::lbphsThreadsDataInit, this, _1 ),
    			boost::bind( &Skim<TVHASH_T>::lbphsThreadLoop, this, _1 ));
  }

  cout << "Kill me now " << totalphits << endl;
  exit(0);

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsPrepareHashOverviewTable(uint32 & readi)
{
  for(auto & iphI=SKIM3_lbphs_idsperhash.begin(); iphI != SKIM3_lbphs_idsperhash.end(); ++iphI){
    iphI->clear();
  }

  bool loopok=true;
  for(;readi<SKIM3_readpool->size() && loopok; ++readi){
    loopok=lbphsPrepareOneSeqForHOT(readi);
  }

//  for(auto iphI=SKIM3_lbphs_idsperhash.begin();; iphI != SKIM3_lbphs_idsperhash.end(); ++iphI){
//    cout << "bhs" << readi << "\t" << iphI-SKIM3_lbphs_idsperhash.begin() << "\t" << iphI->size() << endl;
//  }

  // normalise overview table
  for(auto iphI=SKIM3_lbphs_idsperhash.begin(); iphI != SKIM3_lbphs_idsperhash.end(); ++iphI){
    if(!iphI->empty()){
      // TODO: parallel, non-parallel?
      // 31.10.2015 ... do I still use this code???
      mstd::ssort(*iphI);
      iphI->resize(mstd::unique(*iphI)-iphI->begin());
    }
  }

//  for(auto iphI=SKIM3_lbphs_idsperhash.begin();; iphI != SKIM3_lbphs_idsperhash.end(); ++iphI){
//    cout << "bhs" << readi << "\t" << iphI-SKIM3_lbphs_idsperhash.begin() << "\t" << iphI->size() << endl;
//  }

}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<class TVHASH_T>
bool Skim<TVHASH_T>::lbphsPrepareOneSeqForHOT(const uint32 readi)
{
  const char * seq=SKIM3_readpool->getRead(readi).getSeqAsChar();
  const uint32 slen=SKIM3_readpool->getRead(readi).getLenClippedSeq();
  const uint32 basesperhash=SKIM3_basesperhash;

  bool retvalue=true;

  shash_t acthash=0;
  shash_t hashmask=1;
  // note: bph kann never be >12 in these routines, hence no special handling of bph==16
  //  where the shift would be undefined
  hashmask<<=(basesperhash*2);
  --hashmask;

  uint32  baseok=0;

  for(int32 seqi=0; seqi<static_cast<int32>(slen); ++seqi, ++seq){
    ++baseok;

    acthash<<=2;
    acthash&=hashmask;
    acthash+=SKIM3_lbphs_hashadd[*seq];

    if(SKIM3_lbphs_hashadd[*seq]==0xffff){
      if(dptools::isValidIUPACStarBase(*seq)) {
	// the IUPAC bases are treated like N and X

	// break hash making (which is actually better than behaving
	//  like another character in case of multiple bases with
	//  IUPAC or '*')
	acthash=0;
	baseok=0;
      } else {
	cout << "Unknown base '" << *seq << "' (ASCII " << static_cast<uint16>(*seq) << ") at position " << seqi << " in _CLIPPED_ sequence " << readi << endl; // << actread.getName() << endl;
	exit(100);
      }
    }else if(baseok >= basesperhash){
      SKIM3_lbphs_idsperhash[acthash].push_back(readi);
      if(SKIM3_lbphs_idsperhash[acthash].size()>=512000) retvalue=false;
    }
  }

  if(readi>0 && readi%12000==0) retvalue=false;

  return retvalue;
}



template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsThreadsDataInit(const uint32 numthreads)
{
  FUNCSTART("void Skim<TVHASH_T>::lbphsThreadsDataInit(const uint32 numthreads)");

  SKIM3_lbphsd_vector.resize(numthreads);
  for(uint32 ti=0; ti<numthreads;++ti){
    SKIM3_lbphsd_vector[ti].hashhitperread.clear();
    SKIM3_lbphsd_vector[ti].hashhitperread.resize(SKIM3_readpool->size(),0);

    SKIM3_lbphsd_vector[ti].r1posperhash.clear();
    SKIM3_lbphsd_vector[ti].r1posperhash.resize(1<<(SKIM3_basesperhash*2));

    SKIM3_lbphsd_vector[ti].r1r2poshistogram.clear();
    SKIM3_lbphsd_vector[ti].r1r2poshistogram.resize(100000,0);    // FIXME

    SKIM3_lbphsd_vector[ti].r1r2phused.reserve(100000);
  }

  FUNCEND();
}

template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsThreadLoop(const uint32 threadnr)
{
  FUNCSTART("void Skim<TVHASH_T>::lbphsThreadloop(const uint32 threadnr)");

  // threads need their own try() catch() block

  try {
    CEBUG("Thread: " << threadnr << " starting.\n");

    BUGIFTHROW(threadnr>=SKIM3_lbphsd_vector.size(),"threadnr>=SKIM3_lbphsd_vector.size()???");

    // we'll jump out with a break;
    while(true){
      {
	boost::mutex::scoped_lock mylock(SKIM3_mutex);
	CEBUG("Thread " << threadnr << " waiting ...\n");
	while(!SKIM3_threadcontrol[threadnr].flag_datavalid
	      && ! SKIM3_threadcontrol[threadnr].flag_endthread){
	  SKIM3_master2slavesignal.wait(mylock);
	}
      }
      if(SKIM3_threadcontrol[threadnr].flag_datavalid){
	CEBUG("Thread " << threadnr << " working on " << SKIM3_threadcontrol[threadnr].from << " to " << SKIM3_threadcontrol[threadnr].to << "\n");

	for(uint32 readi=SKIM3_threadcontrol[threadnr].from; readi<SKIM3_threadcontrol[threadnr].to; ++readi){
	  lbphsLookAtRead(readi,threadnr,1);
	  lbphsLookAtRead(readi,threadnr,-1);
	  if(readi%1000==0) cout << "Doing " << readi << "\t" << totalphits << endl;
	  //if(actreadi==5000) exit(0);
	}

	boost::mutex::scoped_lock mylock(SKIM3_mutex);
	SKIM3_threadcontrol[threadnr].flag_datavalid=false;

	SKIM3_slave2mastersignal.notify_one();
      }else if(SKIM3_threadcontrol[threadnr].flag_endthread){
	CEBUG("Thread " << threadnr << "  exiting.\n");
	break;
      }
    }

  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  FUNCEND();
}


template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsLookAtRead(const uint32 actreadi, const uint32 threadnr, const int8 direction)
{
  FUNCSTART("void Skim<TVHASH_T>::lbphsLookAtRead(const uint32 readnr, const uint32 threadnr)");

  lbphs_threaddata_t & lbphsd=SKIM3_lbphsd_vector[threadnr];

  if(direction>0){
    lbphsGetHashesOfOneSeq(actreadi,
			   lbphsd.hashvec,
			   SKIM3_basesperhash,
			   SKIM3_readpool->getRead(actreadi).getClippedSeqAsChar(),
			   SKIM3_readpool->getRead(actreadi).getLenClippedSeq());
  }else{
    lbphsGetHashesOfOneSeq(actreadi,
			   lbphsd.hashvec,
			   SKIM3_basesperhash,
			   SKIM3_readpool->getRead(actreadi).getClippedComplementSeqAsChar(),
			   SKIM3_readpool->getRead(actreadi).getLenClippedSeq());
  }
  sort(lbphsd.hashvec.begin(),lbphsd.hashvec.end());
  lbphsd.hashvec.resize(unique(lbphsd.hashvec.begin(),lbphsd.hashvec.end())-lbphsd.hashvec.begin());

  for(auto & hve : lbphsd.hashvec){
    auto iphI=SKIM3_lbphs_idsperhash[hve].begin();
    uint32 loopcount=SKIM3_lbphs_idsperhash[hve].size();
    for(; loopcount!=0; --loopcount, ++iphI){
      lbphsd.hashhitperread[*iphI]+=1;
    }
  }

  // search for cutoff
  uint32 numhitschosen=0;
  uint32 threshold=0;
  if(1){
    uint32 maxhits=0;
    for(uint32 ri=SKIM_partfirstreadid; ri < SKIM_partlastreadid; ++ri){
      if(lbphsd.hashhitperread[ri]>maxhits) maxhits=lbphsd.hashhitperread[ri];
    }
    threshold=(maxhits*5)/10;
    uint32 diffstep=(maxhits-threshold)/2;
    uint32 numhits=0;
    for(; diffstep>0; diffstep/=2){
      numhits=0;
      for(uint32 ri=SKIM_partfirstreadid; ri < SKIM_partlastreadid; ++ri){
	if(lbphsd.hashhitperread[ri]>threshold) ++numhits;
      }
      if(numhits>25){
	threshold+=diffstep;
      }else if(numhits<10){
	threshold-=diffstep;
      }else{
	diffstep=0;
      }
    }
    numhitschosen=numhits;
  }


  lbphsd.hashhitsofoneread.clear();
  for(uint32 ri=SKIM_partfirstreadid; ri < SKIM_partlastreadid; ++ri){
    if(lbphsd.hashhitperread[ri]>threshold) {
      lbphsd.hashhitsofoneread.resize(lbphsd.hashhitsofoneread.size()+1);
      lbphsd.hashhitsofoneread.back().numhashes=lbphsd.hashhitperread[ri];
      lbphsd.hashhitsofoneread.back().rid=ri;
    }
  }

  sort(lbphsd.hashhitsofoneread.begin(), lbphsd.hashhitsofoneread.end());
  lbphsSecondLineCheck(actreadi, lbphsd, direction);


  numhitschosen=lbphsd.r1r2idstaken.size();
//    cout << "1: " << lbphsd.hashhitsofoneread[0].numhashes
//	 << "\n2: " << lbphsd.hashhitsofoneread[1].numhashes << endl;

  for(uint32 ri=SKIM_partfirstreadid; ri < SKIM_partlastreadid; ++ri){
    //cout << "makesure: " << actreadi << "\t" << ri << "\t# " << lbphsd.hashhitperread[ri] << endl;
    lbphsd.hashhitperread[ri]=0;
  }

  {
    boost::mutex::scoped_lock lock(SKIM3_coutmutex);
    //cout << "nhc\t" << lbphsd.r1r2idstaken.size() << "\n";
    totalphits+=numhitschosen;
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsGetHashesOfOneSeq(const uint32 readi, std::vector<shash_t> & hashvec, const uint32 basesperhash, const char * seq, const uint32 slen)
{
  //const char * seq=SKIM3_readpool->getRead(readi).getSeqAsChar();
  //const uint32 slen=SKIM3_readpool->getRead(readi).getLenClippedSeq();
  //const uint32 basesperhash=SKIM3_basesperhash;

  hashvec.clear();

  shash_t acthash=0;
  shash_t hashmask=1;
  // note: bph kann never be >12 in these routines, hence no special handling of bph==16
  //  where the shift would be undefined
  hashmask<<=(basesperhash*2);
  --hashmask;

  uint32  baseok=0;

  for(int32 seqi=0; seqi<static_cast<int32>(slen); ++seqi, ++seq){
    ++baseok;

    acthash<<=2;
    acthash&=hashmask;
    acthash+=SKIM3_lbphs_hashadd[*seq];

    if(SKIM3_lbphs_hashadd[*seq]==0xffff){
      acthash=0;
      baseok=0;
    }else if(baseok >= basesperhash){
      hashvec.push_back(acthash);
    }
  }

  return;
}




template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsSecondLineCheck(const uint32 actreadid, lbphs_threaddata_t & lbphsd, int8 direction)
{
//  //std::vector<std::vector<int32> > r1posperhash;

  lbphsd.r1r2poshistogram.clear();
  lbphsd.r1r2poshistogram.resize(100000,0);

  lbphsd.r1pphlist.clear();
  lbphsd.r2shashpos.clear();
  lbphsd.r1r2phused.clear();
  lbphsd.r1r2idstaken.clear();


  std::vector<shashpos_t> & r2shashpos=lbphsd.r2shashpos;

  if(direction>0){
    lbphsSeqToPosOverviewTable(lbphsd.r1posperhash,
			       lbphsd.r1pphlist,
			       SKIM3_basesperhash,
			       SKIM3_readpool->getRead(actreadid).getClippedSeqAsChar(),
			       SKIM3_readpool->getRead(actreadid).getLenClippedSeq());
  }else{
    lbphsSeqToPosOverviewTable(lbphsd.r1posperhash,
			       lbphsd.r1pphlist,
			       SKIM3_basesperhash,
			       SKIM3_readpool->getRead(actreadid).getClippedComplementSeqAsChar(),
			       SKIM3_readpool->getRead(actreadid).getLenClippedSeq());
  }

  uint32 dummy=0;

  auto srhhI=lbphsd.hashhitsofoneread.cbegin();
  for(; srhhI!= lbphsd.hashhitsofoneread.end(); ++srhhI){
    lbphsGetHashPosOfOneSeq(srhhI->rid, r2shashpos, SKIM3_basesperhash,
			    SKIM3_readpool->getRead(actreadid).getClippedSeqAsChar(),
			    SKIM3_readpool->getRead(actreadid).getLenClippedSeq());

    lbphsd.r1r2phused.clear();
    auto r2shpI=r2shashpos.cbegin();
    for(; r2shpI!=r2shashpos.end(); ++r2shpI){
      if(lbphsd.r1posperhash[r2shpI->shash].size()){
	auto r1pphI=lbphsd.r1posperhash[r2shpI->shash].cbegin();
	for(; r1pphI!=lbphsd.r1posperhash[r2shpI->shash].cend(); ++r1pphI){
	  int32 posdiff=*r1pphI-r2shpI->hashpos+50000;
	  if(lbphsd.r1r2poshistogram[posdiff]++==0) lbphsd.r1r2phused.push_back(posdiff);
	}
      }
    }

    uint32 maxseen=0;
    int32 indexmax=-1;
    for(auto & rpue : lbphsd.r1r2phused){
      if(lbphsd.r1r2poshistogram[rpue]>maxseen){
	maxseen=lbphsd.r1r2poshistogram[rpue];
	indexmax=rpue;
      }
    }

    if(indexmax>=0){
      if(indexmax>0) maxseen+=lbphsd.r1r2poshistogram[indexmax-1];
      if(indexmax<lbphsd.r1r2poshistogram.size()-2) maxseen+=lbphsd.r1r2poshistogram[indexmax+1];
    }

    //if(maxseen>0){
    //  boost::mutex::scoped_lock lock(SKIM3_coutmutex);
    //  cout << "ms\t" << maxseen << "\n";
    //}
    if(maxseen>=100) {
      ++dummy;
      //lbphsd.r1r2idstaken.push_back(srhhI->rid);
      std::vector<skimhitforsave_t> * actshfsv;
      if(direction>0){
	actshfsv=&lbphsd.shfsvf;
      }else{
	actshfsv=&lbphsd.shfsvr;
      }
      actshfsv->resize(actshfsv->size()+1);
      skimhitforsave_t & shfs=actshfsv->back();
      shfs.rid1=srhhI->rid;
      shfs.rid2=actreadid;
      shfs.eoffset=-(indexmax);
      shfs.percent_in_overlap=50;
      shfs.numhashes=maxseen;
      shfs.ol_stronggood  =false;
      shfs.ol_weakgood    =false;
      shfs.ol_belowavgfreq=false;
      shfs.ol_norept      =false;
      shfs.ol_rept        =false;
    }

    for(auto & rpue : lbphsd.r1r2phused){
      lbphsd.r1r2poshistogram[rpue]=0;
    }

    if(dummy==10) break;
  }

  for(auto & rple : lbphsd.r1pphlist) {
    lbphsd.r1posperhash[rple].clear();
  }
}


template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsGetHashPosOfOneSeq(const uint32 readi, std::vector<shashpos_t> & hashposvec, const uint32 basesperhash, const char * seq, const uint32 slen)
{
  hashposvec.clear();

  shash_t acthash=0;
  shash_t hashmask=1;
  // note: bph kann never be >12 in these routines, hence no special handling of bph==16
  //  where the shift would be undefined
  hashmask<<=(basesperhash*2);
  --hashmask;

  uint32  baseok=0;

  for(int32 seqi=0; seqi<static_cast<int32>(slen); ++seqi, ++seq){
    ++baseok;

    acthash<<=2;
    acthash&=hashmask;
    acthash+=SKIM3_lbphs_hashadd[*seq];

    if(SKIM3_lbphs_hashadd[*seq]==0xffff){
      acthash=0;
      baseok=0;
    }else if(baseok >= basesperhash){
      hashposvec.resize(hashposvec.size()+1);
      hashposvec.back().shash=acthash;
      hashposvec.back().hashpos=seqi;
    }
  }

  return;
}


template<class TVHASH_T>
void Skim<TVHASH_T>::lbphsSeqToPosOverviewTable(std::vector<std::vector<int32> > & posperhash, std::vector<shash_t> & r1pphlist, const uint32 basesperhash, const char * seq, const uint32 slen)
{
  r1pphlist.clear();

  shash_t acthash=0;
  shash_t hashmask=1;
  // note: bph kann never be >12 in these routines, hence no special handling of bph==16
  //  where the shift would be undefined
  hashmask<<=(basesperhash*2);
  --hashmask;

  uint32  baseok=0;

  for(int32 seqi=0; seqi<static_cast<int32>(slen); ++seqi, ++seq){
    ++baseok;

    acthash<<=2;
    acthash&=hashmask;
    acthash+=SKIM3_lbphs_hashadd[*seq];

    if(SKIM3_lbphs_hashadd[*seq]==0xffff){
      acthash=0;
      baseok=0;
    }else if(baseok >= basesperhash){
      posperhash[acthash].push_back(seqi);
      if(posperhash[acthash].size()==1) r1pphlist.push_back(acthash);
    }
  }

  return;
}
