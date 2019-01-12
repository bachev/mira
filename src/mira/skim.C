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

#include "mira/skim.H"

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include "errorhandling/errorhandling.H"

#include "util/fileanddisk.H"
#include "util/dptools.H"
#include "util/progressindic.H"

#include "util/stlimprove.H"

#include "mira/ads.H"
#include "mira/readpool.H"
#include "mira/vhash.H"

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


//#define CEBUG(bla)   {cout << bla; cout.flush();}


template<typename TVHASH_T>
const TVHASH_T Skim<TVHASH_T>::SKIM3_MAXVHASHMASK(0xffffffUL);


// sort id1 low to high,
//  on same id1 by skimweight high to low,
//  on same numhashes by id2 low to high
//  on same id2 by eoffset low to high
bool skimhitforsave_t::stdSortCmp(const skimhitforsave_t & a, const skimhitforsave_t & b)
{
  if(a.rid1 == b.rid1){
    if(a.numhashes == b.numhashes){
      if(a.rid2 == b.rid2){
	return a.eoffset < b.eoffset;
      }else{
	return a.rid2 < b.rid2;
      }
    }else{
      return a.numhashes > b.numhashes;
    }
  }
  return a.rid1 < b.rid1;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

// Plain vanilla constructor
template<typename TVHASH_T>
Skim<TVHASH_T>::Skim()
{
  FUNCSTART("Skim<TVHASH_T>::Skim(ReadPool & rp)");

  SKIM3_logflag_purgeunnecessaryhits=false;
  init();

  FUNCEND();
}


template<typename TVHASH_T>
void Skim<TVHASH_T>::init()
{
  FUNCSTART("Skim<TVHASH_T>::init()");


  SKIM3_numthreads=2;
  SKIM3_basesperhash=16;
  SKIM3_hashsavestepping=4;
  SKIM3_overlaplenrequired.clear();
  for(uint32 i=0;i<ReadGroupLib::getNumSequencingTypes(); i++){
    SKIM3_overlaplenrequired.push_back(20);
    SKIM3_percentrequired.push_back(50);
  }

  SKIM3_totalhitschosen=0;

  FUNCEND()
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
Skim<TVHASH_T>::~Skim()
{
  FUNCSTART("Skim<TVHASH_T>::~Skim()");
  //  ERROR("Not implemented yet.");

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void Skim<TVHASH_T>::discard()
{
  FUNCSTART("Skim<TVHASH_T>::discard()");

  FUNCEND();
}



/*************************************************************************
 *
 * returns number of megahubs
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
uint32 Skim<TVHASH_T>::skimGo(ReadPool & rp,
		    std::string               & posfmatchname,
		    std::string               & poscmatchname,
		    bannedoverlappairs_t & bannedoverlaps,
		    std::vector<uint32>       & overlapcounter,
		    std::vector<uint32>       & writtenhitsperid,
		    std::vector<int32>        & chuntleftcut,
		    std::vector<int32>        & chuntrightcut,
		    std::vector<std::vector<uint8>>        & overlapcritlevell,
		    std::vector<std::vector<uint8>>        & overlapcritlevelr,
		    std::vector<uint8>        * meghubsptr,
		    uint32 numthreads,
		    uint32 maxmemusage,
		    bool onlyagainstrails,
		    bool alsocheckreverse,
		    uint32 basesperhash,
		    uint8  hss,
		    const std::vector<int32> & percentrequired,
		    const std::vector<int32> & overlaplenrequired,
		    uint32 maxhitsperread,
		    uint32 megahubcap,
		    bool forcetakestronggood)
{
  FUNCSTART("uint32 Skim<TVHASH_T>::skimGo( ... )");

  BUGIFTHROW(posfmatchname.empty(),"posfmatchname.empty() ???");
  BUGIFTHROW(poscmatchname.empty(),"poscmatchname.empty() ???");

  dateStamp(cout);
  Read::setCoutType(Read::AS_CLIPPEDFASTA);

  init();

  SKIM3_readpool=&rp;

  SKIM3_numthreads=numthreads;
  if(SKIM3_numthreads<1) SKIM3_numthreads=1;
  if(SKIM3_numthreads>256) SKIM3_numthreads=256;

  if(megahubcap<150000) megahubcap=150000;
  SKIM3_megahubcap=megahubcap;

  // if extalso is true, the take not only the clipped sequence, but
  //  also the extends to it

  SKIM3_onlyagainstrails=onlyagainstrails;

  if(basesperhash > sizeof(TVHASH_T)*4){
    basesperhash= sizeof(TVHASH_T)*4;
  }
  SKIM3_basesperhash=basesperhash;
  SKIM3_hashsavestepping=hss;
  SKIM3_percentrequired=percentrequired;
  SKIM3_overlaplenrequired=overlaplenrequired;
  SKIM3_maxhitsperread=maxhitsperread;
  SKIM3_forcetakestronggood=forcetakestronggood;

  SKIM3_overlapcritlevelvl=&overlapcritlevell;
  SKIM3_overlapcritlevelvr=&overlapcritlevelr;

  SKIM3_megahubsptr=meghubsptr;

  // TODO: check whether these can be re-used between assembly passes
  //  would reduce skim files even a bit further in later passes
  //  Question: can an overlap criterion level become worse through editing???
  SKIM3_overlapcritlevelvl->clear();
  SKIM3_overlapcritlevelvl->resize(2);
  (*SKIM3_overlapcritlevelvl)[0].resize(SKIM3_readpool->size(),255);
  (*SKIM3_overlapcritlevelvl)[1].resize(SKIM3_readpool->size(),255);
  SKIM3_overlapcritlevelvr->clear();
  SKIM3_overlapcritlevelvr->resize(2);
  (*SKIM3_overlapcritlevelvr)[0].resize(SKIM3_readpool->size(),255);
  (*SKIM3_overlapcritlevelvr)[1].resize(SKIM3_readpool->size(),255);
  SKIM3_largestencasementscoretodate.clear();
  SKIM3_largestencasementscoretodate.resize(SKIM3_readpool->size(),0);
  SKIM3_nomorehitseie.clear();
  SKIM3_nomorehitseie.resize(SKIM3_readpool->size(),0);

  SKIM3_overlapcounter=&overlapcounter;
  SKIM3_overlapcounter->clear();
  SKIM3_overlapcounter->resize(SKIM3_readpool->size(),0);

  SKIM3_writtenhitsperid=&writtenhitsperid;
  SKIM3_writtenhitsperid->clear();

  SKIM3_posfmatchnextchecksize=SKIM3_SKIMMATCHFIRSTCHECK;
  SKIM3_poscmatchnextchecksize=SKIM3_posfmatchnextchecksize;

  SKIM3_bannedoverlaps=&bannedoverlaps;


  // for chimerahunt (if wished), but not if only against rails
  SKIM3_chimerahunt.clear();
  SKIM3_chuntleftcut=&chuntleftcut;
  SKIM3_chuntrightcut=&chuntrightcut;
  if(!onlyagainstrails && chuntleftcut.size()!=0){
    chuntleftcut.clear();
    chuntleftcut.resize(SKIM3_readpool->size(),0);
    chuntrightcut.clear();
    chuntrightcut.resize(SKIM3_readpool->size(),0);
    SKIM3_chimerahunt.resize(SKIM3_readpool->size());
    for(uint32 i=0; i<SKIM3_readpool->size(); i++){
      SKIM3_chimerahunt[i].resize(SKIM3_readpool->getRead(i).getLenClippedSeq(),0);
    }
  }


  SKIM_partfirstreadid=0;         // partition first read id
  SKIM_partlastreadid=0;         // partition last read id

  SKIM3_posfmatchfname=posfmatchname;
  SKIM3_posfmatchfout.open(posfmatchname, std::ios::out| std::ios::trunc | std::ios::binary);
  SKIM3_poscmatchfname=poscmatchname;
  SKIM3_poscmatchfout.open(poscmatchname, std::ios::out| std::ios::trunc | std::ios::binary);

  //std::ofstream mout;
  //mout.open(megahublogname, std::ios::out| std::ios::trunc);

  uint32 numpartitions=computePartition(maxmemusage*SKIM3_hashsavestepping,true);

  CEBUG("We will get " << numpartitions << " partitions.\n");

  //cout << "\nSKIMMER: Using maximum of " << maxmemusage << " hashes stored in memory, " << numpartitions << " partitions will be computed." << endl << endl;

  CEBUG("Progressend: " << SKIM_progressend << endl);

  if(SKIM3_megahubsptr!=nullptr){
    SKIM3_megahubsptr->clear();
    SKIM3_megahubsptr->resize(SKIM3_readpool->size(),0);
  }
  SKIM3_fullencasedcounter.resize(SKIM3_readpool->size(),0);

  fillTagStatusInfoOfReads();

  if(0){
  }else{
    cout << "Now running threaded and partitioned skimmer-" << SKIM3_basesperhash << " with " << numpartitions << " partitions in " << SKIM3_numthreads << " threads:" << endl;

    SKIM_progressindicator= new ProgressIndicator<int64>(0,SKIM_progressend);

    SKIM3_vhraparray.clear();
    for(uint32 actpartition=1; actpartition<=numpartitions; actpartition++){
      CEBUG("\nWorking on partition " << actpartition << "/" << numpartitions << endl);

      computePartition(maxmemusage*SKIM3_hashsavestepping,false);

      CEBUG("Will contain read IDs " << SKIM_partfirstreadid << " to " << SKIM_partlastreadid-1 << endl);

      prepareSkim(SKIM_partfirstreadid, SKIM_partlastreadid, SKIM3_vhraparray,true);
      if(!SKIM3_vhraparray.empty()){
	CEBUG("Checking forward hashes" << endl);
	startMultiThreading(1,
			    SKIM3_numthreads,
			    5000,
			    SKIM_partfirstreadid,
			    SKIM3_readpool->size(),
			    boost::bind( &Skim<TVHASH_T>::cfhThreadsDataInit, this, _1 ),
			    boost::bind( &Skim<TVHASH_T>::cfhThreadLoop, this, _1 ));
	purgeMatchFileIfNeeded(1);
	if(alsocheckreverse){
	  CEBUG("Checking reverse hashes" << endl);
	  startMultiThreading(-1,
			      SKIM3_numthreads,
			      5000,
			      SKIM_partfirstreadid,
			      SKIM3_readpool->size(),
			      boost::bind( &Skim<TVHASH_T>::cfhThreadsDataInit, this, _1 ),
			      boost::bind( &Skim<TVHASH_T>::cfhThreadLoop, this, _1 ));
	  purgeMatchFileIfNeeded(-1);
	}
	CEBUG("Done." << endl);
      }

      SKIM_partfirstreadid=SKIM_partlastreadid;
    }

    SKIM_progressindicator->finishAtOnce();
    delete SKIM_progressindicator;
    cout << " done.\n";
  }

  SKIM3_posfmatchfout.close();
  SKIM3_poscmatchfout.close();

  SKIM3_writtenhitsperid->resize(SKIM3_readpool->size(),0);

  {
    std::vector<uint8> perfectrailmatches;
    bool hasrails=false;
    for(uint32 ri=0; ri<SKIM3_readpool->size(); ++ri){
      if(SKIM3_readpool->getRead(ri).isRail()){
	hasrails=true;
	break;
      }
    }
    if(hasrails){
      perfectrailmatches.resize(SKIM3_readpool->size(),0);
      findPerfectRailMatchesInSkimFile(SKIM3_posfmatchfname,1,perfectrailmatches);
      findPerfectRailMatchesInSkimFile(SKIM3_poscmatchfname,-1,perfectrailmatches);
    }
    purgeUnnecessaryHitsFromSkimFile(SKIM3_posfmatchfname,1,perfectrailmatches);
    purgeUnnecessaryHitsFromSkimFile(SKIM3_poscmatchfname,-1,perfectrailmatches);
  }

  for(uint32 i=0; i<SKIM3_writtenhitsperid->size(); ++i){
    SKIM3_totalhitschosen+=(*SKIM3_writtenhitsperid)[i];
  }

  uint32 megahubs=0;

  if(SKIM3_megahubsptr!=nullptr){
    for(uint32 mhi=0; mhi<SKIM3_megahubsptr->size(); mhi++){
      if((*SKIM3_megahubsptr)[mhi]>0) {
	megahubs++;
      }
    }
  }
  //cout << "\nSkim summary:\n\taccepted: " << SKIM3_acceptedhits << "\n\tpossible: " << SKIM3_possiblehits  << "\n\tpermbans: " << SKIM3_totalpermbans;
  cout << "\n\nHits chosen: " << SKIM3_totalhitschosen << "\n\n";

//  mout.close();
  dateStamp(cout);

  cout << endl;

  if(SKIM3_chimerahunt.size()){
    chimeraHuntLocateChimeras();
  }



  FUNCEND();
  return megahubs;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void Skim<TVHASH_T>::fillTagStatusInfoOfReads()
{
  SKIM3_hasMNRr.clear();
  SKIM3_hasSRMr.clear();
  SKIM3_hasFpas.clear();
  SKIM3_hasMNRr.resize(SKIM3_readpool->size(),0);
  SKIM3_hasSRMr.resize(SKIM3_readpool->size(),0);
  SKIM3_hasFpas.resize(SKIM3_readpool->size(),0);

  for(uint32 actreadid=0; actreadid<SKIM3_readpool->size(); actreadid++){
    for(uint32 tn=0; tn<SKIM3_readpool->getRead(actreadid).getNumOfTags(); tn++){
      if(SKIM3_readpool->getRead(actreadid).getTag(tn).identifier==Read::REA_tagentry_idMNRr) SKIM3_hasMNRr[actreadid]=1;
      if(SKIM3_readpool->getRead(actreadid).getTag(tn).identifier==Read::REA_tagentry_idSRMr) SKIM3_hasSRMr[actreadid]=1;
      if(SKIM3_readpool->getRead(actreadid).getTag(tn).identifier==Read::REA_tagentry_idSOFApolyA_sequence) SKIM3_hasFpas[actreadid]=1;
    }
  }
}

/*************************************************************************
 *
 * computes either
 *  - starting from the first read id, the last read id to use for the
 *    next partition (SKIM_partfirstreadid and SKIM_partlastreadid)
 * or
 *  - starting from the first read id, the total number of partitions
 * to use in the next skim run.(returned)
 *
 * The maxmemusage is the dominating factor: each partition may not use
 *  (much) more than that.
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

template<typename TVHASH_T>
uint32 Skim<TVHASH_T>::computePartition(uint32 maxmemusage, bool computenumpartitions)
{
  FUNCSTART("uint32 Skim<TVHASH_T>::computePartition(uint32 maxmemusage, bool computenumpartitions)");

  uint32 numpartitions=0;
  uint32 totalseqlen=0;
  uint32 maxseqlen=0;

  SKIM_partlastreadid=SKIM_partfirstreadid;

  if(computenumpartitions) {
    SKIM_progressend=SKIM3_readpool->size()-SKIM_partlastreadid;
  }

  for(; SKIM_partlastreadid<SKIM3_readpool->size(); SKIM_partlastreadid++) {
    if(!SKIM3_readpool->getRead(SKIM_partlastreadid).hasValidData()
       || !SKIM3_readpool->getRead(SKIM_partlastreadid).isUsedInAssembly()) continue;

    if(SKIM3_readpool->getRead(SKIM_partlastreadid).getLenClippedSeq() > SKIM3_MAXREADSIZEALLOWED) {
      MIRANOTIFY(Notify::FATAL,"Read " << SKIM3_readpool->getRead(SKIM_partlastreadid).getName() << " is longer than SKIM3_MAXREADSIZEALLOWED (" << SKIM3_MAXREADSIZEALLOWED << ") bases. SKIM cannot handle this, aborting.\n");
    }

    maxseqlen=std::max(maxseqlen,SKIM3_readpool->getRead(SKIM_partlastreadid).getLenClippedSeq());

    totalseqlen+=SKIM3_readpool->getRead(SKIM_partlastreadid).getLenClippedSeq();
    //if(SKIM_takeextalso){
    //  totalseqlen+=SKIM3_readpool->getRead(SKIM_partlastreadid).getRightExtend();
    //}

    if(totalseqlen>maxmemusage) {
      if(computenumpartitions){
	totalseqlen=0;
	numpartitions++;

	SKIM_progressend+=SKIM3_readpool->size()-SKIM_partlastreadid;

      }else{
	SKIM_partlastreadid++;
	break;
      }
    }
  }

  if(computenumpartitions) {
    if(totalseqlen>0) {
      numpartitions++;
      SKIM_progressend+=SKIM3_readpool->size()-SKIM_partlastreadid;
    }
    SKIM_progressend*=2;
  }

  //// Compute SKIM_maxoffsets
  //if(maxseqlen>0){
  //  // this beauty is slow, but a cute hack to find out the number of
  //  //  the highest bit which power of two fits the given value
  //  // e.g. 1024 -> 10 -> 2^10 = 1024 fits 1024
  //  //      1025 -> 11 -> 2^11 = 2048 fits 1025
  //  SKIM_mo_shiftmultiplier=0;
  //  while((maxseqlen-1) >> ++SKIM_mo_shiftmultiplier);
  //}else{
  //  SKIM_mo_shiftmultiplier=4;
  //}
  //
  //// restrict shift multiplier to 11 (and therefore maxoffsets to 2048)
  //if(SKIM_mo_shiftmultiplier>11) SKIM_mo_shiftmultiplier=11;
  //SKIM_maxoffsets=(1<<SKIM_mo_shiftmultiplier);

  FUNCEND();

  return numpartitions;
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::prepareSkim(uint32 fromid, uint32 toid, std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> & vhraparray, bool assemblychecks)
{
  FUNCSTART("void Skim<TVHASH_T>::prepareSkim(bool alsocheckreverse)");

  vhraparray.clear();

  uint32 totalseqlen=0;
  uint32 totalseqs=0;

  for(uint32 seqnr=fromid; seqnr<toid; seqnr++) {
    if(!SKIM3_readpool->getRead(seqnr).hasValidData()) continue;
    if(assemblychecks
       && (!SKIM3_readpool->getRead(seqnr).isUsedInAssembly()
	   || (SKIM3_onlyagainstrails && !SKIM3_readpool->getRead(seqnr).isRail()))) continue;
    totalseqlen+=SKIM3_readpool->getRead(seqnr).getLenClippedSeq();
    totalseqs++;
    //if(SKIM_takeextalso) totalseqlen+=SKIM3_readpool->getRead(i).getRightExtend();
  }

  //dateStamp(cout);
  CEBUG("\nPreparing skim data: "  << fromid << " to " << toid << endl);
  CEBUG(totalseqs << " sequences to skim, totalling " << totalseqlen << " bases." << endl);


  // next loop:
  //  transform each read into a series of forward
  //   hashes, store them into the array along with info whether
  //   each hash position is valid or not,
  //  also store the info how many hashes each sequence produced

  uint32 totalhashes=0;

  if(totalseqlen>0){

    vhraparray.resize(totalseqlen/SKIM3_hashsavestepping);
    auto vhraparrayI=vhraparray.begin();
    std::vector<uint8> tagmaskvector;

    //ProgressIndicator P(partfirstreadid, partlastreadid);
    for(uint32 seqnr=fromid; seqnr < toid; seqnr++){
      //P.progress(seqnr);
      Read & actread= SKIM3_readpool->getRead(seqnr);
      if(!actread.hasValidData()) continue;
      if(assemblychecks
	 && (!SKIM3_readpool->getRead(seqnr).isUsedInAssembly()
	     || (SKIM3_onlyagainstrails && !SKIM3_readpool->getRead(seqnr).isRail()))) continue;

      uint32 slen=actread.getLenClippedSeq();
      //if(SKIM_takeextalso) slen+=actread.getRightExtend();

      const std::vector<Read::bposhashstat_t> & bposhashstats=actread.getBPosHashStats();
      int32 bfpos=actread.calcClippedPos2RawPos(0);
      int32 bfposinc=1;

      if(slen>=8) {
	fillTagMaskVector(seqnr, tagmaskvector);
	uint32 hashesmade= transformSeqToVariableHash(
	  seqnr,
	  actread,
	  actread.getClippedSeqAsChar(),
	  slen,
	  SKIM3_basesperhash,
	  vhraparrayI,
	  false,
	  SKIM3_hashsavestepping,
	  tagmaskvector,
	  bposhashstats,
	  bfpos,
	  bfposinc
	  );

	//CEBUG(seqnr << "\t" << totalhashes << "\t" << slen << endl);
	totalhashes+=hashesmade;
      }
    }

    //P.progress(partlastreadid);
    CEBUG("Totalseqlen " << totalseqlen << endl);
    CEBUG("Computed " << totalhashes << " linkpoints." << endl);

    if(totalhashes>0){
      CEBUG("Resizing array" << endl);
      vhraparray.resize(totalhashes);

      if(0){
	CEBUG("Partition unsorted:\n");
	auto vaI=vhraparray.cbegin();
	for(; vaI!=vhraparray.cend(); ++vaI){
	  cout << *vaI << '\n';
	}
	cout << "###########################" << endl;
      }

      CEBUG("Sorting array" << endl);
      mstd::psort(vhraparray, Skim<TVHASH_T>::sortVHRAPArrayElem_);

      if(0){
	CEBUG("Partition sorted:\n");
	auto vaI=vhraparray.cbegin();
	for(; vaI!=vhraparray.cend(); ++vaI){
	  cout << *vaI << '\n';
	}
	cout << "###########################" << endl;
      }

    }
  }

  CEBUG("Making shortcuts" << endl);
  makeVHRAPArrayShortcuts(vhraparray, SKIM3_basesperhash);

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void Skim<TVHASH_T>::purgeMatchFileIfNeeded(int8 direction)
{
  FUNCSTART("void Skim<TVHASH_T>::purgeMatchFileIfNeeded(int8 direction)");

  // never needed "during the run" when doing a mapping only, but towards
  //  end of skimGo() when purging using information on which reads have
  //  perfect hits.
  if(SKIM3_onlyagainstrails) return;

  std::ofstream * posmatchfout=nullptr;
  std::string * fname=nullptr;
  uint64 * nextchecksize=nullptr;
  if(direction>0){
    posmatchfout=&SKIM3_posfmatchfout;
    fname=&SKIM3_posfmatchfname;
    nextchecksize=&SKIM3_posfmatchnextchecksize;
  }else{
    posmatchfout=&SKIM3_poscmatchfout;
    fname=&SKIM3_poscmatchfname;
    nextchecksize=&SKIM3_poscmatchnextchecksize;
  }

  BUGIFTHROW(fname->empty(),"fname->empty() ???");

  if(posmatchfout->tellp() >= *nextchecksize){
    (*nextchecksize)+=SKIM3_SKIMMATCHCHECKINCR;

    posmatchfout->close();
    std::vector<uint8> dummy;
    purgeUnnecessaryHitsFromSkimFile(*fname,direction,dummy);

    posmatchfout->open(*fname, std::ios::out|std::ios::app);
    if(!posmatchfout){
      MIRANOTIFY(Notify::FATAL, "Could not reopen SKIM match file " << *fname);
    }
    if(posmatchfout->tellp() >= *nextchecksize){
      (*nextchecksize)=SKIM3_SKIMMATCHCHECKINCR+posmatchfout->tellp();
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 * Go through written skim hits on disk.
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void Skim<TVHASH_T>::findPerfectRailMatchesInSkimFile(std::string & filename, const int8 rid2dir, std::vector<uint8> & prmatches)
{
  FUNCSTART("void Skim<TVHASH_T>::findPerfectRailMatchesInSkimFile(std::string & filename, const int8 rid2dir, std::vector<uint8> & prmatches)");

  // temporary skim container
  std::vector<skimhitforsave_t> tsc;

  FILE * finfout;
  finfout = fopen(filename.c_str(),"r+");
  if(finfout == nullptr) {
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  //// which reads have a perfect rail match?
  //std::vector<uint8> prmatches;
  //prmatches.resize(SKIM3_readpool->size(),0);

  myFSeek(finfout, 0, SEEK_END);
  auto finsize=myFTell(finfout);
  rewind(finfout);

  uint64 lineno=0;

  ADSEstimator adse;

  long freadpos=0;
  long fwritepos=0;

  while(!feof(finfout)){
    tsc.resize(500000);
    myFSeek(finfout, freadpos, SEEK_SET);

    auto numread=myFRead(&tsc[0],sizeof(skimhitforsave_t),tsc.capacity(),finfout);

    if(numread==0) break;
    CEBUG("rsh4_pUHFNSF: read " << numread << endl;)
    lineno+=numread;

    freadpos=myFTell(finfout);
    CEBUG("new freadpos: " << freadpos << endl);

    if(numread<tsc.capacity()) tsc.resize(numread);

    for(const auto & tsce : tsc){
      // one of both must be rail
      if((*SKIM3_readpool)[tsce.rid1].isRail()
	 ||(*SKIM3_readpool)[tsce.rid2].isRail()){
	// but not both
	if(!((*SKIM3_readpool)[tsce.rid1].isRail() && (*SKIM3_readpool)[tsce.rid2].isRail())){
	    // we want a perfect match, i.e. 100%
	  if(tsce.percent_in_overlap == 100){
	    adse.calcNewEstimateFromSkim(
	      tsce.eoffset,
	      (*SKIM3_readpool)[tsce.rid1].getLenClippedSeq(),
	      (*SKIM3_readpool)[tsce.rid2].getLenClippedSeq(),
	      tsce.rid1,
	      tsce.rid2,
	      1,
	      rid2dir);
	    if(adse.getContainmentLevel()!=0){
	      uint32 railid=tsce.rid1;
	      uint32 nreadid=tsce.rid2;
	      if((*SKIM3_readpool)[tsce.rid2].isRail()){
		std::swap(railid,nreadid);
	      }
	      prmatches[nreadid]=1;
	    }
	  }
	}
      }
    }

  }
}

/*************************************************************************
 *
 * Go through written skim hits on disk. Compare saved hits to best level found
 * If both reads agree they have better partners at hand, throw out skim hit.
 * Do this only if both reads are not rails, else it might be that good, fully
 *  encased hits with 99% are thrown out because of seemingly better partially
 *  overlapping hits at 100% (and the part outside the overlap then wreaks
 *  havoc in reality with some non-identical repeat part)
 *
 * Careful: this rewrites and truncates the original file.
 * Length will awlays be <= initial length
 *
 * Careful: side effect SKIM3_writtenhitsperid gets filled with
 *  values ... this MUST happen!
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void Skim<TVHASH_T>::purgeUnnecessaryHitsFromSkimFile(std::string & filename, const int8 rid2dir, std::vector<uint8> & prmatches)
{
  FUNCSTART("void Skim<TVHASH_T>::purgeUnnecessaryHitsFromSkimFile(std::string & filename, const int8 rid2dir, std::vector<uint8> & prmatches)");

  BUGIFTHROW(filename.empty(),"filename.empty() ???");

  // temporary skim container
  std::vector<skimhitforsave_t> tsc;

  FILE * finfout;
  finfout = fopen(filename.c_str(),"r+");
  if(finfout == nullptr) {
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  std::ofstream logfout;
  if(SKIM3_logflag_purgeunnecessaryhits){
    std::string path,justfilename;
    splitFullPathAndFileName(filename,path,justfilename);
    std::string logfilename=path+"/elog.skim.puh."+justfilename;
    cout << "\nSkim: elog " << logfilename << endl;
    logfout.open(logfilename, std::ios::out|std::ios::trunc);
  }

  // some statistics: how many hits does each read have?
  std::vector<uint32> hitstats;
  hitstats.resize(SKIM3_readpool->size(),0);
  std::ofstream logcount;
  if(SKIM3_logflag_purgeunnecessaryhits){
    std::string path,justfilename;
    splitFullPathAndFileName(filename,path,justfilename);
    std::string logfilename=path+"/elog.skim.puh.count."+justfilename;
    cout << "\nSkim: elog count " << logfilename << endl;
    logcount.open(logfilename, std::ios::out|std::ios::trunc);
  }

  myFSeek(finfout, 0, SEEK_END);
  auto finsize=myFTell(finfout);
  rewind(finfout);

  cout << "\ntruncating " << filename << endl;

  uint64 lineno=0;
  uint32 bannedoverlapsfound=0;
  size_t totalhits=0;

  ADSEstimator adse;

  long freadpos=0;
  long fwritepos=0;

  while(!feof(finfout)){
    tsc.resize(500000);
    myFSeek(finfout, freadpos, SEEK_SET);

    auto numread=myFRead(&tsc[0],sizeof(skimhitforsave_t),tsc.capacity(),finfout);

    if(numread==0) break;
    CEBUG("rsh4_pUHFNSF: read " << numread << endl;)
    lineno+=numread;

    freadpos=myFTell(finfout);
    CEBUG("new freadpos: " << freadpos << endl);

    if(numread<tsc.capacity()) tsc.resize(numread);

    uint8 ocll=255;
    uint8 oclr=255;

    auto readI=tsc.cbegin();
    auto writeI=tsc.begin();
    for(; readI != tsc.cend(); ++readI){
      bool del1=false;
      bool del2=false;
      bool del3=false;
      bool del4=false;
      bool del5=false;
      if(!(*SKIM3_readpool)[readI->rid1].isRail()
	 && !(*SKIM3_readpool)[readI->rid2].isRail()){
	adse.calcNewEstimateFromSkim(
	  readI->eoffset,
	  (*SKIM3_readpool)[readI->rid1].getLenClippedSeq(),
	  (*SKIM3_readpool)[readI->rid2].getLenClippedSeq(),
	  readI->rid1,
	  readI->rid2,
	  1,
	  rid2dir);

	Skim<TVHASH_T>::getOverlapCriterionLevel(readI->rid1,
				       (*SKIM3_readpool)[readI->rid1].getSequencingType(),
				       adse,readI->percent_in_overlap,
				     ocll,oclr);
	uint32 ocvi=0;  // overlap criterion vector index 0 is for norept overlaps
	if(!readI->ol_norept) ocvi=1;  // index 1 is for all other rept overlaps (among them the rept)

	auto & oclvl = (*SKIM3_overlapcritlevelvl)[ocvi];
	auto & oclvr = (*SKIM3_overlapcritlevelvr)[ocvi];

	CEBUG("OCL: " << readI->rid1 << " " << static_cast<uint16>(ocll) << " " << static_cast<uint16>(oclr) << endl);

	// if it's a rail, then the read has no saying in the decision whether this skim should be deleted
	// if not, look at overlap criterion level left and right
	if((*SKIM3_readpool)[readI->rid1].isRail()
	   || (ocll>oclvl[readI->rid1]
	       && oclr>oclvr[readI->rid1])){
	  del1=true;
	}
	// if it's a Solexa and the best overlapcritlevel is not 0, we are probably in a under-coverage
	//  situation ... try to account for that by being less harsh
	// Values: in getOverlapCriterionLevel(), current overlapcritlevel for Solexa goes from
	//  0-29 for 100% matches
	//  30-59 for 99%
	//  ... etc up to 95% (including, == 149 max)
	// Therefore, if best overlapcritlevel != 0 && <= 29, then the best overlap is 100% albeit
	//  not as long as it could be when in high coverage situations
	// Therefore: low coverage
	// Therefore: we'll take all 100% matches for that read
	if(del1 && SKIM3_readpool->getRead(readI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	  if(ocll <= 29
	     && oclvl[readI->rid1]!= 0
	     && oclvl[readI->rid1] <= 29){
	    del1=false;
	  }else if(oclr <= 29
		   && oclvr[readI->rid1]!= 0
		   && oclvr[readI->rid1] <= 29){
	    del1=false;
	  }
	}

	Skim<TVHASH_T>::getOverlapCriterionLevel(readI->rid2,
				       (*SKIM3_readpool)[readI->rid2].getSequencingType(),
				       adse,readI->percent_in_overlap,
				       ocll,oclr);
	CEBUG("OCL: " << readI->rid2 << " " << static_cast<uint16>(ocll) << " " << static_cast<uint16>(oclr) << endl);
	// if it's a rail, then the read has no saying in the decision whether this skim should be deleted
	// if not, look at overlap criterion level left and right
	if((*SKIM3_readpool)[readI->rid2].isRail()
	   || (ocll>oclvl[readI->rid2]
	       && oclr>oclvr[readI->rid2])){
	  del2=true;
	}

	// if it's a Solexa and the best overlapcritlevel is not 0, we are probably in a under-coverage
	//  situation ... try to account for that by being less harsh
	if(del2 && SKIM3_readpool->getRead(readI->rid2).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	  if(ocll <= 29
	     && oclvl[readI->rid2]!= 0
	     && oclvl[readI->rid2] <= 29){
	    del2=false;
	  }else if(oclr <= 29
		   && oclvr[readI->rid2]!= 0
		   && oclvr[readI->rid2] <= 29){
	    del2=false;
	  }
	}

	// Special Solexa tests
	if(SKIM3_readpool->getRead(readI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	   && SKIM3_readpool->getRead(readI->rid2).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){

	  // Solexa elitist tests
	  // Solexa elitists (<=5 left/right) do not want to play with pariahs (both >+5 levels left/right)
	  if(oclvl[readI->rid1] <= 5
	     && oclvr[readI->rid1] <= 5){

	    if(oclvl[readI->rid2] > oclvl[readI->rid1]+5
	       && oclvr[readI->rid2] > oclvr[readI->rid1]+5){
	      del3=true;
	    }
	  }else if(oclvl[readI->rid2] <= 5
		   && oclvr[readI->rid2] <= 5){

	    if(oclvl[readI->rid1] > oclvl[readI->rid2]+5
	       && oclvr[readI->rid1] > oclvr[readI->rid2]+5){
	      del3=true;
	    }
	  }

	  // Solexa enough is enough tests
	  // but only if there is no strong good overlap

	  //bang here: take presence of HAF5-7 / MNRr into account. How?

	  if(readI->ol_stronggood && SKIM3_forcetakestronggood){
	    // forced take of strong good overlaps regardless of enough-is-enough
	    //  to ensure sg hits in highly repetitive reads get taken (should be on for de-novo genome only)
	    // TODO: do we need to do something or is this branch really empty? Think about it.
	  }else{
	    bool sdel1=false;
	    bool sdel2=false;
	    const uint8 ocritthresh=3;
	    if(hitstats[readI->rid1] >= 350){
	      SKIM3_nomorehitseie[readI->rid1]=1;
	      sdel1=true;
	    }else if(hitstats[readI->rid1] >= 50
		     && oclvl[readI->rid1]<=ocritthresh
		     && oclvr[readI->rid1]<=ocritthresh
	      ){
	      SKIM3_nomorehitseie[readI->rid1]=1;
	      sdel1=true;
	    }else if(hitstats[readI->rid1] >= 100
		     && (oclvl[readI->rid1]<=ocritthresh
			 || oclvr[readI->rid1]<=ocritthresh)){
	      SKIM3_nomorehitseie[readI->rid1]=1;
	      sdel1=true;
	    }

	    if(hitstats[readI->rid2] >= 350){
	      SKIM3_nomorehitseie[readI->rid2]=1;
	      sdel2=true;
	    }else if(hitstats[readI->rid2] >= 50
		     && oclvl[readI->rid2]<=ocritthresh
		     && oclvr[readI->rid2]<=ocritthresh
	      ){
	      SKIM3_nomorehitseie[readI->rid2]=1;
	      sdel2=true;
	    } else if(hitstats[readI->rid2] >= 100
		      && (oclvl[readI->rid2]<=ocritthresh
			  || oclvr[readI->rid2]<=ocritthresh)){
	      SKIM3_nomorehitseie[readI->rid2]=1;
	      sdel2=true;
	    }

	    if(sdel1&&sdel2) del4=true;
	  }
	}
      }

      if(readI->percent_in_overlap!=100
	 && !prmatches.empty()){
	if((*SKIM3_readpool)[readI->rid1].isRail() && prmatches[readI->rid2]){
	  del5=true;
	}else if((*SKIM3_readpool)[readI->rid2].isRail() && prmatches[readI->rid1]){
	  del5=true;
	}
      }

      // is there a template set and is it the same for both reads? keep that overlap
      if(SKIM3_readpool->getRead(readI->rid1).getTemplateID() >= 0
	 && SKIM3_readpool->getRead(readI->rid1).getTemplateID() == SKIM3_readpool->getRead(readI->rid2).getTemplateID()){
	del1=false;
	del2=false;
	del3=false;
	del4=false;
	del5=false;
      }

      // both reads PacBioLQ? currently don't purge
      if(SKIM3_readpool->getRead(readI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)
	 && SKIM3_readpool->getRead(readI->rid2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)){
	del1=false;
	del2=false;
	del3=false;
	del4=false;
	del5=false;
      }

      // strong good and rept? we're at a border, don't purge
      if(readI->ol_stronggood && readI->ol_rept){
	del1=false;
	del2=false;
	del3=false;
	del4=false;
	del5=false;
      }

//	del1=false;
//	del2=false;
//	del3=false;
//	del4=false;
//	del5=false;

      if(readI != writeI){
	*writeI=*readI;
      }
      ++writeI;
      CEBUG("DEL: " << del1 << " " << del2 << " " << del3 << " " << del4 << " " << del5 << endl);
      //CEBUG(readI->rid1 << ": " << static_cast<uint16>(oclvl[readI->rid1]) << " " << static_cast<uint16>(oclvr[readI->rid1])
      //	    << "\t\t" << readI->rid2 << ": " << static_cast<uint16>(oclvl[readI->rid2]) << " " << static_cast<uint16>(oclvr[readI->rid2]) << endl);

      if((del1 && del2) || del3 || del4 || del5){
	//if(0){
	--writeI;
	CEBUG("Purged: " << *readI);
	if(SKIM3_logflag_purgeunnecessaryhits){
	  logfout << "Purged:\t" << SKIM3_readpool->getRead(readI->rid1).getName()
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[0][readI->rid1])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[0][readI->rid1]) << ")"
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[1][readI->rid1])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[1][readI->rid1]) << ")"
		  << '\t' << SKIM3_readpool->getRead(readI->rid2).getName()
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[0][readI->rid2])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[0][readI->rid2]) << ")"
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[1][readI->rid2])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[1][readI->rid2]) << ")"
		  << "\t(" << static_cast<uint16>(ocll)
		  << "," << static_cast<uint16>(oclr) << ")"
		  << '\t' << *readI;
	}
      }else{
	if(std::min(readI->rid1,readI->rid2) == 0) {CEBUG("DINGO! ")};
	CEBUG("Kept: " << *readI);
	if(SKIM3_logflag_purgeunnecessaryhits){
	  logfout << "Kept:\t" << SKIM3_readpool->getRead(readI->rid1).getName()
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[0][readI->rid1])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[0][readI->rid1]) << ")"
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[1][readI->rid1])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[1][readI->rid1]) << ")"
		  << '\t' << SKIM3_readpool->getRead(readI->rid2).getName()
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[0][readI->rid2])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[0][readI->rid2]) << ")"
		  << " (" << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[1][readI->rid2])
		  << "," << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[1][readI->rid2]) << ")"
		  << "\t(" << static_cast<uint16>(ocll)
		  << "," << static_cast<uint16>(oclr) << ")"
		  << '\t' << *readI;
	}
	if(!SKIM3_writtenhitsperid->empty()){
	  (*SKIM3_writtenhitsperid)[readI->rid1]+=1;
	  (*SKIM3_writtenhitsperid)[readI->rid2]+=1;
	}

	++hitstats[readI->rid1];
	++hitstats[readI->rid2];
      }
    }

    // the resize thing is really not optimal ... one could write to file only a
    //  subset. However, at the moment just keep it for 100% safety
    CEBUG("Purge skim data. Old size: " << tsc.size() << endl);
    tsc.resize(tsc.size()-(readI-writeI));
    CEBUG("New size: " << tsc.size() << endl);

    if(!tsc.empty()){
      myFSeek(finfout, fwritepos, SEEK_SET);
      if(myFWrite(&tsc[0],
		  sizeof(skimhitforsave_t),
		  tsc.size(),
		  finfout) != tsc.size()){
	MIRANOTIFY(Notify::FATAL, "Could not write anymore to normalised skim file. Disk full? Changed permissions?");
      }
      fwritepos=myFTell(finfout);
      CEBUG("new fwritepos: " << fwritepos << endl);
    }
  }

  fclose(finfout);

  cout << "truncated " << filename << " from " << finsize << " to " << fwritepos << endl;
  if(truncate(filename.c_str(),fwritepos)){
    MIRANOTIFY(Notify::FATAL, "Could not truncate normalised skim file? Strange ...");
  }

  for(uint32 ri=0; ri<hitstats.size(); ++ri){
    logcount << SKIM3_readpool->getRead(ri).getName() << "\t" << hitstats[ri] << "\t"
	     << static_cast<uint32>((*SKIM3_overlapcritlevelvl)[0][ri]) << "\t"
	     << static_cast<uint32>((*SKIM3_overlapcritlevelvr)[0][ri]) << "\t"
	     << static_cast<uint32>((*SKIM3_overlapcritlevelvl)[1][ri]) << "\t"
	     << static_cast<uint32>((*SKIM3_overlapcritlevelvr)[1][ri]) << "\t";

    if(SKIM3_nomorehitseie[ri]) logcount << "enough";
    logcount << '\n';
  }
  logcount.close();

//  tsc.resize(500000);
//  finfout = fopen(filename.c_str(),"r+");
//  size_t numread=myFRead(&tsc[0],sizeof(skimhitforsave_t),500000,finfout);
//  CEBUG("Read anew: " << numread << endl);
//  tsc.resize(numread);
//  for(uint32 i=0;i < numread; i++){
//    CEBUG("dbgcheck: " << tsc[i]);
//  }

  FUNCEND();
  return;
}
//#define CEBUG(bla)






/*************************************************************************
 *
 * beware: SKIM3_vashortcuts_* arrays may be empty at the return of this
 *  function in case the vhraparray itself was empty! Account for that in
 *  the search functions!
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::makeVHRAPArrayShortcuts(std::vector<typename HashStatistics<TVHASH_T>::vhrap_t> & vhraparray, const uint32 basesperhash)
{
  //cout << "Making VHRAPArrayShortcuts" << endl;

  SKIM3_vashortcuts_begin.clear();
  SKIM3_vashortcuts_end.clear();
  SKIM3_completevhraparray_end=vhraparray.end();
  auto vaI=vhraparray.cbegin();
  if(vaI==vhraparray.end()) return;

  SKIM3_vashortcuts_begin.resize(
    1<<(std::min(static_cast<uint32>(12),basesperhash)*2),
    vhraparray.end()
    );

  SKIM3_vashortcuts_end.resize(
    1<<(std::min(static_cast<uint32>(12),basesperhash)*2),
    vhraparray.end()
    );

  TVHASH_T acthash= (vaI->vhash & SKIM3_MAXVHASHMASK);
  while(vaI != vhraparray.end()){
    auto vasi=static_cast<uint64>(acthash);
    SKIM3_vashortcuts_begin[vasi]=vaI;
    for(;vaI != vhraparray.end() && (vaI->vhash & SKIM3_MAXVHASHMASK) == acthash; vaI++) ;
    SKIM3_vashortcuts_end[vasi]=vaI;
    //cout << "vhash: " << hex << acthash << "\t" << dec << SKIM3_vashortcuts_end[acthash]-SKIM3_vashortcuts_begin[acthash] << '\n';
    if(vaI != vhraparray.end()) acthash= vaI->vhash & SKIM3_MAXVHASHMASK;
  }
}

//#define CEBUG(bla)
//#define CEBUGF(bla)




/*************************************************************************
 *
 * TODO: this is a mess, rewrite
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

template<typename TVHASH_T>
uint32 Skim<TVHASH_T>::transformSeqToVariableHash (const uint32 readid, const Read & actread, const char * seq, uint32 slen, const uint32 basesperhash, typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::iterator & vhraparrayI, const bool countonly, const uint8 hashsavestepping, std::vector<uint8> & tagmaskvector, const std::vector<Read::bposhashstat_t> & bposhashstats, int32 bfpos, const int32 bfposinc)
{
  FUNCSTART("void Skim<TVHASH_T>::transformSeqToVariableHash (...)");

  //BUGIFTHROW(basesperhash>32, "basesperhash > 32 ?");
  BUGIFTHROW(hashsavestepping<1, "hashsavestepping < 1 ?");
  if(slen>SKIM3_MAXREADSIZEALLOWED){
    MIRANOTIFY(Notify::FATAL,"Read " << actread.getName() << " is " << slen << " bp long and thus longer than SKIM3_MAXREADSIZEALLOWED (" << SKIM3_MAXREADSIZEALLOWED << ") bases. Skim cannot handle that, sorry.");
  }

  if(!tagmaskvector.empty()){
    BUGIFTHROW(tagmaskvector.size()!=slen,"tagmaskvector.size() " << tagmaskvector.size() << " != slen " << slen);
  }

//  Read::setCoutType(Read::AS_TEXT);
  CEBUG(actread);
  CEBUG("readid: " << readid << endl);
  CEBUG("seq: " << seq << endl);
  CEBUG("strlen(seq): " << strlen(seq) << endl);
  CEBUG("slen: " << slen << endl);
  CEBUG("bfpos: " << bfpos << endl);
  CEBUG("bfposinc: " << bfposinc << endl);

  auto bhsI=bposhashstats.cbegin();
  advance(bhsI,bfpos);

  TVHASH_T lasthash(0);
  TVHASH_T acthash(0);
  TVHASH_T hashmask(1);
  // *grml* undefined behaviour of left shift for 64 shifts in a 64 bit type makes this cludge necessary
  // VLuint would handle it by itself, but not a normal uint64
  if(basesperhash>=sizeof(TVHASH_T)*4){
    hashmask=0;
  }else{
    hashmask<<=(basesperhash*2);
  }
  --hashmask;

  // strictly speaking, the following two are no kmer hashes (TVHASH_T)
  // BUT: they need to be able to have one bit per baseperhash to follow the masking
  // making them of the same type as TVHASH_T helps to simplify coding
  TVHASH_T nonmaskedposbitvector(0);
  TVHASH_T nmpmask(1);
  // *grml* undefined behaviour of left shift for n shifts in a n-bit type makes this cludge necessary
  // VLuint would handle it by itself, but not a normal uint32/uint64
  if(basesperhash>=sizeof(TVHASH_T)*8){
    nmpmask=0;
  }else{
    nmpmask<<=basesperhash;
  }
  --nmpmask;

  CEBUG("sizeof TVHASH_T: " << sizeof(TVHASH_T) << '\n');
  CEBUG("bases per hash: " << static_cast<uint16>(basesperhash) << '\n');
  CEBUG("hashsavestepping: " << static_cast<uint16>(hashsavestepping) << '\n');
  CEBUG("hash mask: " << hex << hashmask << dec << '\n');
  CEBUG("nmpmask: " << hex << nmpmask << dec << '\n');

  // first hash made must also be saved
  uint8 hashsavecounter=1;

  uint32 goods=0;
  uint32 bads=0;
  uint32  baseok=0;
  auto initial_vaI=vhraparrayI;
  auto tmvI=tagmaskvector.cbegin();

  CEBUG("Hashing " << actread.getName() << '\t' << slen << '\t' << strlen(seq) << '\t' << actread.getLenClippedSeq() << "\n");


  char actbase;
  bool mustsavelasthash=false;
  uint16 lastposhashsaved=0;

  bool notright=false;
  for(uint32 seqi=0; seqi<slen; ++seqi, ++seq){
    // careful here: no continue statement until tmvI has been increased!!!
    // not doing this in loop increment as HashStatistics::checkBaitHit() calls us with empty vector
    lasthash=acthash;
    acthash<<=2;
    acthash&=hashmask;
    baseok++;

    actbase=static_cast<char>(toupper(*seq));

    CEBUG(seqi << '\t' << actbase << endl);

    switch (actbase) {
    case 'A' : break;
    case 'C' : {
      acthash|=1;
      break;
    }
    case 'G' : {
      acthash|=2;
      break;
    }
    case 'T' : {
      acthash|=3;
      break;
    }
    default : {
      if(dptools::isValidIUPACStarBase(actbase)) {
	// the IUPAC bases are treated like N and X

	mustsavelasthash=true;

	// break hash making (which is actually better than behaving
	//  like another character in case of multiple bases with
	//  IUPAC or '*')
	acthash=0;
	baseok=0;
      } else {
	cout << "Unknown base '" << *seq << "' (ASCII " << static_cast<uint16>(*seq) << ") at position " << seqi << " in _CLIPPED_ sequence " << actread.getName() << endl;
	exit(100);
      }
    }
    }

    // handling of masked positions
    bool lastposhadunmasked=false;
    if(static_cast<bool>(nonmaskedposbitvector)) lastposhadunmasked=true;
    nonmaskedposbitvector<<=1;
    nonmaskedposbitvector&=nmpmask;

    if(tagmaskvector.empty()){
      nonmaskedposbitvector|=1;
    }else{
      if(*tmvI==0){
	nonmaskedposbitvector|=1;
      }
      ++tmvI;
    }

    // using "continue" is now allowed (if needed)

    if(lastposhadunmasked && static_cast<bool>(nonmaskedposbitvector)==false
       && lastposhashsaved != seqi-1) {
      mustsavelasthash=true;
    }

    CEBUG(seqi << ' ' << *seq << ' ' << hex << nonmaskedposbitvector << dec << ' ' << mustsavelasthash << ' ');
    if(baseok >= basesperhash) {
      goods++;
      if(!countonly && static_cast<bool>(nonmaskedposbitvector)) {
	if(mustsavelasthash){
	  vhraparrayI->vhash=lasthash;
	  vhraparrayI->readid=readid;
	  vhraparrayI->hashpos=seqi-1;
	  vhraparrayI->bhashstats=(bhsI-bfposinc)->getBHashStat(-bfposinc);
	  // getBHashStat(-bfposinc) because while we're running "forward", we save
	  //  hashes only with a delay of 'basesperhash' and need to know the status
	  //  of the "past" bases ... and this info is readily available in
	  //  the BHashStat of the other strand

	  CEBUG("saved LG hash: " << *vhraparrayI << '\n');

	  lastposhashsaved=seqi-1;
	  vhraparrayI++;
	  // set hashsavecounter to 1 so that the next good hash
	  //  generated is saved!
	  hashsavecounter=1;
	  mustsavelasthash=false;
	} else if(--hashsavecounter == 0){
	  vhraparrayI->vhash=acthash;
	  vhraparrayI->readid=readid;
	  vhraparrayI->hashpos=seqi;
	  vhraparrayI->bhashstats=bhsI->getBHashStat(-bfposinc);
	  // getBHashStat(-bfposinc) because while we're running "forward", we save
	  //  hashes only with a delay of 'basesperhash' and need to know the status
	  //  of the "past" bases ... and this info is readily available in
	  //  the BHashStat of the other strand

	  CEBUG("saved hash: " << *vhraparrayI << '\n');

	  vhraparrayI++;
	  lastposhashsaved=seqi;
	  hashsavecounter=hashsavestepping;
	}
      }

    } else {
      //cout << "Missed hash" << endl;
      mustsavelasthash=false;
      if(seqi>=basesperhash) {
	bads++;
      }
    }
    CEBUG('\n');

    // this is a hack to make this routine work in -D_GLIBCXX_DEBUG mode
    // normally, this iterator should be handled by the for() statement,
    //  but in reverese cases where the read has no left clip, after the
    //  last loop the iterator would advance to the "-1" position where the
    //  STL debug containers are not happy with, even though the iterator
    //  wouldn't be used as the for() loops would stop right there
    // therefore, this cludge

    if(bfposinc<0){
      if(notright){
	MIRANOTIFY(Notify::FATAL, "Something's not right here.");
      }
      if(bhsI!=bposhashstats.begin()) {
	bhsI+=bfposinc;
      }else{
	notright=true;
      }
    }else{
      bhsI+=bfposinc;
    }
  }

  //for(uint32 i=0; i<basesperhash; i++, hashp++, hashokp++) {
  //  *hashp=0;
  //  *hashokp=0;
  //}

#ifdef CEBUGFLAG
  CEBUG("goods: " << goods << endl);
  CEBUG("bads: " << bads << endl);
#endif

  return (vhraparrayI-initial_vaI);
}

//#define CEBUG(bla)
//#define CEBUGF(bla)





//#define CEBUG(bla)   {cout << bla; cout.flush();}


//#define CEBUG(bla)   {boost::mutex::scoped_lock lock(SKIM3_coutmutex); cout << bla; cout.flush();}

// TODO: bad: direction should not be in this call, more of the called function
template<typename TVHASH_T>
void Skim<TVHASH_T>::startMultiThreading(const int8 direction, const uint32 numthreads, const uint32 readsperthread, const uint32 firstid, const uint32 lastid, boost::function<void(uint32_t)> initfunc, boost::function<void(uint32_t)> callfunc)
{
  // initialise task specific data by task specific init routine
  initfunc(numthreads);

  // initialise the data structure with which the master
  //  process (well, this process) will communicate with the
  //  worker threads
  // do this *before* creating the threads :-)
  {
    threadworkercontrol_t twc;
    twc.from=0;
    twc.to=0;
    twc.direction=direction;
    twc.flag_datavalid=false;
    twc.flag_endthread=false;

    SKIM3_threadcontrol.clear();
    SKIM3_threadcontrol.resize(numthreads,twc);
  }

  // create the number of worker threads we will use
  boost::thread_group workerthreads;
  for(uint32 i=0; i<numthreads;i++){
    workerthreads.create_thread(boost::bind(callfunc, i));
  }

  // main work distribution loop
  // gives each thread a part of the search space. If no
  //  thread is free, waits for a slave2master signal
  //  (which currently can only mean a thread has finished
  //  going through it's search space)

  uint32 startid=firstid;
  while(startid < lastid) {
    boost::mutex::scoped_lock mylock(SKIM3_mutex);

    // search thread that is idle
    uint32 tnr=0;
    for(; tnr<numthreads; tnr++){
      if(SKIM3_threadcontrol[tnr].flag_datavalid==false) break;
    }
    if(tnr==numthreads) {
      // no idle thread?
      //  well, wait for a slave2master signal
      SKIM3_slave2mastersignal.wait(mylock);
    }else{
      uint32 endid=startid+readsperthread;
      if(endid>lastid) endid=lastid;

      CEBUG("Giving " << startid << " to " << endid << " to thread " << tnr << "\n");
      SKIM3_threadcontrol[tnr].from=startid;
      SKIM3_threadcontrol[tnr].to=endid;
      SKIM3_threadcontrol[tnr].flag_datavalid=true;

      SKIM3_master2slavesignal.notify_all();

      startid=endid;
    }
  }

  // no more work to distribute
  // tell workerthreads to end as soon as they finished their
  //  current task
  CEBUG("Last packet given, flagging all threads the stop signal.\n");

  {
    boost::mutex::scoped_lock mylock(SKIM3_mutex);
    for(uint32 tnr=0; tnr<numthreads; tnr++){
      SKIM3_threadcontrol[tnr].flag_endthread=true;
    }
  }
  SKIM3_master2slavesignal.notify_all();

  // and wait for all threads of the threadgroup
  //  to return
  workerthreads.join_all();

}
//#define CEBUG(bla)

template<typename TVHASH_T>
void Skim<TVHASH_T>::cfhThreadsDataInit(const uint32 numthreads)
{
  FUNCSTART("void Skim<TVHASH_T>::cfhThreadsDataInit(const uint32 numthreads)");

  SKIM3_cfhd_vector.resize(numthreads);
  for(uint32 ti=0; ti<numthreads;++ti){
    SKIM3_cfhd_vector[ti].readhashmatches.clear();
    SKIM3_cfhd_vector[ti].readhashmatches.reserve(500000);
    SKIM3_cfhd_vector[ti].smallhist4repeats.clear();
    SKIM3_cfhd_vector[ti].smallhist4repeats.reserve(100);
    SKIM3_cfhd_vector[ti].singlereadvhraparray.clear();
    SKIM3_cfhd_vector[ti].singlereadvhraparray.reserve(5000);
    SKIM3_cfhd_vector[ti].tmpmatchwith.clear();
    SKIM3_cfhd_vector[ti].tmpmatchwith.reserve(2000);
    SKIM3_cfhd_vector[ti].tagmaskvector.clear();
    SKIM3_cfhd_vector[ti].tagmaskvector.reserve(2000);
    SKIM3_cfhd_vector[ti].shfsv.clear();
    SKIM3_cfhd_vector[ti].shfsv.reserve(100000);
    SKIM3_cfhd_vector[ti].ridswithmatches.clear();
    SKIM3_cfhd_vector[ti].ridswithmatches.reserve(10000);

    SKIM3_cfhd_vector[ti].uidswithnewcritlevelvl.clear();
    SKIM3_cfhd_vector[ti].uidswithnewcritlevelvl.resize(2);
    SKIM3_cfhd_vector[ti].uidswithnewcritlevelvr.clear();
    SKIM3_cfhd_vector[ti].uidswithnewcritlevelvr.resize(2);
    SKIM3_cfhd_vector[ti].critlevellofnewuidsv.clear();
    SKIM3_cfhd_vector[ti].critlevellofnewuidsv.resize(2);
    SKIM3_cfhd_vector[ti].critlevelrofnewuidsv.clear();
    SKIM3_cfhd_vector[ti].critlevelrofnewuidsv.resize(2);
    for(uint32 ocvi=0; ocvi<2; ++ocvi){
      SKIM3_cfhd_vector[ti].uidswithnewcritlevelvl[ocvi].reserve(10000);
      SKIM3_cfhd_vector[ti].uidswithnewcritlevelvr[ocvi].reserve(10000);
      SKIM3_cfhd_vector[ti].critlevellofnewuidsv[ocvi].reserve(10000);
      SKIM3_cfhd_vector[ti].critlevelrofnewuidsv[ocvi].reserve(10000);
    }
  }
  FUNCEND();
}

template<typename TVHASH_T>
void Skim<TVHASH_T>::cfhThreadLoop(const uint32 threadnr)
{
  FUNCSTART("void Skim<TVHASH_T>::threadloop(const uint32 threadnr)");

  // threads need their own try() catch() block

  try {
    CEBUG("Thread: " << threadnr << " starting.\n");

    BUGIFTHROW(threadnr>=SKIM3_cfhd_vector.size(),"threadnr>=SKIM3_cfhd_vector.size()???");
    cfh_threaddata_t & cfhd=SKIM3_cfhd_vector[threadnr];

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

	cfhd.posmatchfout=&SKIM3_posfmatchfout;
	if(SKIM3_threadcontrol[threadnr].direction<0) cfhd.posmatchfout=&SKIM3_poscmatchfout;
	checkForHashes_fromto(SKIM3_threadcontrol[threadnr].direction,
			      SKIM3_threadcontrol[threadnr].from,
			      SKIM3_threadcontrol[threadnr].to,
			      cfhd);

	boost::mutex::scoped_lock mylock(SKIM3_mutex);
	SKIM3_threadcontrol[threadnr].flag_datavalid=false;

	SKIM3_slave2mastersignal.notify_one();
      }else if(SKIM3_threadcontrol[threadnr].flag_endthread){
	CEBUG("Thread " << threadnr << "  exiting.\n");
	break;
      }
    }

    if(cfhd.shfsv.size()){
      boost::mutex::scoped_lock lock(SKIM3_resultfileoutmutex);
      cfhd.posmatchfout->write(reinterpret_cast<char*>(&cfhd.shfsv[0]),sizeof(skimhitforsave_t)*cfhd.shfsv.size());
      if(cfhd.posmatchfout->bad()){
	MIRANOTIFY(Notify::FATAL, "Could not write anymore to skimhit save6. Disk full? Changed permissions?");
      }
      cfhd.shfsv.clear();
    }

  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  FUNCEND();
}

//#define CEBUG(bla)




/*************************************************************************
 *
 *
 * TODO: counting fullencased is not thread safe atm! is it needed? should
 *        be approximative anyway
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::checkForHashes_fromto(const int8 direction, const uint32 fromid, const uint32 toid, cfh_threaddata_t & cfhd)
{
  FUNCSTART("void Skim<TVHASH_T>::checkForHashes_fromto(const int8 direction, const uint32 fromid, const uint32 toid, cfh_threaddata_t & cfhd)");

  // really?
  //BUGIFTHROW(Read::getNumSequencingTypes() >4, "Must be reworked for new sequencing types! (encasement shortcuts & others?");

  if(SKIM3_vashortcuts_begin.empty() || SKIM3_vashortcuts_end.empty()) return;

  cfhd.readhashmatches.clear();
  cfhd.singlereadvhraparray.clear();
  cfhd.tmpmatchwith.clear();
  cfhd.tagmaskvector.clear();;
  // do NOT clear std::vector<skimhitforsave_t> shfsv !!!
  cfhd.ridswithmatches.clear();
  for(uint32 ocvi=0; ocvi<2; ++ocvi){
    cfhd.uidswithnewcritlevelvl[ocvi].clear();
    cfhd.uidswithnewcritlevelvr[ocvi].clear();
    cfhd.critlevellofnewuidsv[ocvi].clear();
    cfhd.critlevelrofnewuidsv[ocvi].clear();
  }

  for(uint32 actreadid=fromid; actreadid<toid; actreadid++){
    //if(actreadid>100) return;

    // don't need to go through identified megahubs again
    if(SKIM3_megahubsptr!=nullptr && (*SKIM3_megahubsptr)[actreadid]>0) continue;

    // if this read has been fully encased by other reads, then also
    //  skip it
    if(SKIM3_fullencasedcounter[actreadid]) continue;

    Read & actread= SKIM3_readpool->getRead(actreadid);
    if(!actread.hasValidData()
      || !actread.isUsedInAssembly()) continue;

    uint32 slen=actread.getLenClippedSeq();

    if(slen<SKIM3_basesperhash) continue;

    cfhd.singlereadvhraparray.resize(slen);

    auto srvaI=cfhd.singlereadvhraparray.begin();
    const std::vector<Read::bposhashstat_t> & bposhashstats=actread.getBPosHashStats();

    uint32 hashesmade;
    if(direction>0) {
      fillTagMaskVector(actreadid, cfhd.tagmaskvector);
      int32 bfpos=actread.calcClippedPos2RawPos(0);
      int32 bfposinc=1;
      hashesmade=transformSeqToVariableHash(
	actreadid,
	actread,
	actread.getClippedSeqAsChar(),
	slen,
	SKIM3_basesperhash,
	srvaI,
	false,
	1,
	cfhd.tagmaskvector,
	bposhashstats,
	bfpos,
	bfposinc
	);
    }else{
      // TODO: first fill, then reverse is ... stupid
      // but do this only if some mask was really set (saves a bit of time)
      if(fillTagMaskVector(actreadid, cfhd.tagmaskvector)){
	mstd::reverse(cfhd.tagmaskvector);
      }
      int32 bfpos=actread.calcClippedComplPos2RawPos(0);
      int32 bfposinc=-1;
      hashesmade=transformSeqToVariableHash(
	actreadid,
	actread,
	actread.getClippedComplementSeqAsChar(),
	slen,
	SKIM3_basesperhash,
	srvaI,
	false,
	1,
	cfhd.tagmaskvector,
	bposhashstats,
	bfpos,
	bfposinc
	);
    }

    cfhd.singlereadvhraparray.resize(hashesmade);

    srvaI=cfhd.singlereadvhraparray.begin();
    uint32 truetestsm2hits=0;
    for(; srvaI != cfhd.singlereadvhraparray.end(); ++srvaI){
      auto lowerbound=SKIM3_vashortcuts_begin[static_cast<uint64>(srvaI->vhash & SKIM3_MAXVHASHMASK)];
      auto upperbound=SKIM3_vashortcuts_end[static_cast<uint64>(srvaI->vhash & SKIM3_MAXVHASHMASK)];

      // "SKIM3_empty_vector_vhrap_t.end()" is the "nullptr" replacement
      if(SKIM3_completevhraparray_end != lowerbound){
	if(SKIM3_basesperhash>12){
	  // with more than 12 bases in a hash, the vhrap array is
	  //  subdivided
	  auto p=equal_range(lowerbound,
			     upperbound,
			     *srvaI,
			     compareVHRAPArrayElem_);
	  lowerbound=p.first;
	  upperbound=p.second;
	}

	for(;lowerbound!=upperbound; lowerbound++){
	  truetestsm2hits++;

	  CEBUG("/// " << actreadid << '\t' << lowerbound->readid << '\n');

	  // hmmmm .....
	  // original: if(actreadid > lowerbound->readid){
	  // this fails spectacularly for mapping now that rails shifted to end of pool
	  // correct resolution would be adding
	  //
	  // but this might slow down the search quite a bit
	  if(actreadid > lowerbound->readid){
	    // NO! do not check this here ... terrible time penalty
	    // do that in checkForPotentialHits() !
	      //&& SKIM3_nomorehitseie[actreadid] == 0
	      //&& SKIM3_nomorehitseie[lowerbound->readid] == 0){
	    CEBUG("/// take!\n");
	    cfhd.readhashmatches.resize(cfhd.readhashmatches.size()+1);
	    cfhd.readhashmatches.back().rid2=lowerbound->readid;
	    cfhd.readhashmatches.back().hashpos1=srvaI->hashpos;
	    cfhd.readhashmatches.back().hashpos2=lowerbound->hashpos;
	    cfhd.readhashmatches.back().eoffset=srvaI->hashpos - lowerbound->hashpos;
	    cfhd.readhashmatches.back().bhashstats=srvaI->bhashstats;
	  }
	}
      }
    }

    if(actreadid % 1 == 0) {
      CEBUG("actreadid: " << actreadid << "\treadhashmatches.size(): " << cfhd.readhashmatches.size() << "\ttruetestsm2hits: " << truetestsm2hits << endl);
    }

    // Hmmm ... this does not represent the full truth, but without "if" there
    // is a partition effect in the data. Bad, but cannot be helped (except
    // going through all reads in all partitions, which is unnecessary for the
    // search itself and effectively.doubles the SKIM time)
    //if(SKIM_partfirstreadid==0)(*SKIM3_rawhashitcounter)[actreadid]+=truetestsm2hits;

    if(cfhd.readhashmatches.size()>0){
      bool ismegahub=false;
      if(!SKIM3_readpool->getRead(actreadid).isRail()
	 && !SKIM3_readpool->getRead(actreadid).isBackbone()
	 && SKIM3_megahubsptr!=nullptr && truetestsm2hits>SKIM3_megahubcap) {
	CEBUG("Potential megahub: " << actreadid << "\treadhashmatches.size(): " << cfhd.readhashmatches.size() << endl);

	// ok, potential megahub. To save the situation, iteratively throw out
	//  all hashes with a frequency > nlevel
	// Throws out first frequency 7, then 6, then 5 as last resort.
	// If the size of readhashmatches gets below the cap,
	//  then it's not treated as megahub

	for(uint32 nlevel=6; nlevel>=4 && cfhd.readhashmatches.size()>SKIM3_megahubcap; --nlevel){
	  auto dstI=cfhd.readhashmatches.begin();
	  for(auto srcI=dstI; srcI != cfhd.readhashmatches.end(); ++srcI){
	    *dstI=*srcI;
	    if(srcI->bhashstats.getFrequency()<=nlevel){
	      ++dstI;
	    }
	  }
	  cfhd.readhashmatches.resize(dstI-cfhd.readhashmatches.begin());
	}

	if(cfhd.readhashmatches.size() > SKIM3_megahubcap) {
	  CEBUG("Megahub confirmed: " << actreadid << "\treadhashmatches.size(): " << cfhd.readhashmatches.size() << endl);
	  (*SKIM3_megahubsptr)[actreadid]=1;
	  ismegahub=true;
	}
      }
      if(!ismegahub) {
	checkForPotentialHits(direction, actreadid, cfhd.tmpmatchwith, cfhd.readhashmatches, cfhd.smallhist4repeats);

	selectPotentialHitsForSave2(direction, actreadid,
				    cfhd);

      }
      cfhd.readhashmatches.clear();
    }
  }

  {
    boost::mutex::scoped_lock lock(SKIM3_coutmutex);
    SKIM_progressindicator->increaseprogress(toid-fromid);
  }

}

//#define CEBUG(bla)





/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::selectPotentialHitsForSave2(const int8 direction, const uint32 actreadid, cfh_threaddata_t & cfhd)
{
  FUNCSTART("void Skim<TVHASH_T>::selectPotentialHitsForSave2(const int8 direction, const uint32 actreadid, cfh_threaddata_t & cfhd)");

  cfhd.ridswithmatches.clear();
  if(cfhd.tmpmatchwith.empty()) return;
  CEBUG("start selectPotentialHitsForSave2()\n");

  updateCriterionLevels(direction,actreadid,
			cfhd);

  // use non-parallel sort as we are already in a multithreading environment here
  mstd::ssort(cfhd.tmpmatchwith, sortMWByEstimScore_);

  ADSEstimator adse;

  uint32 numleftext=SKIM3_maxhitsperread/2;
  uint32 numrightext=numleftext;
  uint8 ocll=255;
  uint8 oclr=255;

  bool takenrailfulllength=false;

  // take all which have the same or lower level (if it is not 255)
  // using iterator only because of CEBUG() below
  for(auto tmwI=cfhd.tmpmatchwith.begin(); tmwI!=cfhd.tmpmatchwith.end(); ++tmwI) {
    CEBUG("### " << tmwI-cfhd.tmpmatchwith.begin() << " " << (*SKIM3_readpool)[actreadid].getName() << " " << (*SKIM3_readpool)[tmwI->otherid].getName());
    CEBUG("\ntmwI: " << *tmwI << endl);

    adse.calcNewEstimateFromSkim(
      tmwI->eoffset,
      (*SKIM3_readpool)[actreadid].getLenClippedSeq(),
      (*SKIM3_readpool)[tmwI->otherid].getLenClippedSeq(),
      actreadid,
      tmwI->otherid,
      1,
      direction);
    CEBUG("ADSE: " << adse << endl);

    getOverlapCriterionLevel(actreadid,
			     SKIM3_readpool->getRead(actreadid).getSequencingType(),
			     adse,
			     static_cast<uint8>(tmwI->percent_in_overlap),
			     ocll,
			     oclr);

    uint32 ocvi=0;  // overlap criterion vector index 0 is for norept overlaps
    if(!tmwI->ol_norept) ocvi=1;  // index 1 is for all other rept overlaps (among them the rept)

    auto & oclvl = (*SKIM3_overlapcritlevelvl)[ocvi];
    auto & oclvr = (*SKIM3_overlapcritlevelvr)[ocvi];

    CEBUG("actreadid critlevel: " << static_cast<uint16>(oclvl[actreadid]) << " " << static_cast<uint16>(oclvr[actreadid]) << '\n');
    CEBUG("ocll: " << static_cast<uint16>(ocll) << " " << static_cast<uint16>(oclr) << '\n');

    bool regulartake=false;

    if((ocll!=255 && ocll<=oclvl[actreadid])
       || (oclr!=255 && oclr<=oclvr[actreadid])){
      regulartake=true;
    }else if(SKIM3_readpool->getRead(actreadid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
      // if it's a Solexa and the best overlapcritlevel is not 0, we are probably in a under-coverage
      //  situation ... try to account for that by being less harsh
      // see purgeUnnecessaryHitsFromSkimFile() for the value "29"
      if(ocll <= 29
	 && oclvl[actreadid]!= 0
	 && oclvl[actreadid] <= 29){
	regulartake=true;
      }else if(oclr <= 29
	       && oclvr[actreadid]!= 0
	       && oclvr[actreadid] <= 29){
	regulartake=true;
      }
    }

    if(regulartake){
      CEBUG("actreadid critlevel: " << static_cast<uint16>(oclvl[actreadid]) << " " << static_cast<uint16>(oclvr[actreadid]) << '\n');
      tmwI->taken=true;
      if(numleftext>0 && adse.getEstimatedLeftExpand(actreadid)>0) --numleftext;
      if(numrightext>0 && adse.getEstimatedRightExpand(actreadid)>0) --numrightext;
      CEBUG("+++++++++++++ take critlevel. nle " << numleftext << "\tnre: " << numrightext << '\n');
    } else if(tmwI->ol_fulllength){
      if((*SKIM3_readpool)[actreadid].isRail() || (*SKIM3_readpool)[tmwI->otherid].isRail()){
	// if one of the reads is a rail and there is a full-length overlap, we have to
	//  take that hit no matter what.
	// Reason: either it's a 100% hit, then it's obvious, or its a <100% hit and then
	//  the Smith-Waterman *needs* to have a look at that bugger to get the very best
	//  placement
	tmwI->taken=true;
	takenrailfulllength=true;
	CEBUG("++++++++++++ take rail\n");
      }else if((*SKIM3_readpool)[actreadid].getLenClippedSeq()>=(*SKIM3_readpool)[tmwI->otherid].getLenClippedSeq()){
	// if first read is larger and encase completely second read, take overlap
	//  if score is >= largest encasement score of other read seen to date
	if(tmwI->estimscore>=SKIM3_largestencasementscoretodate[tmwI->otherid]){
	  if(tmwI->estimscore>SKIM3_largestencasementscoretodate[tmwI->otherid]){
	    // is it worthwhile to make this thread safe???
	    SKIM3_largestencasementscoretodate[tmwI->otherid]=tmwI->estimscore;
	  }
	  tmwI->taken=true;
	  CEBUG("++++++++++++ take encasement\n");
	}
      }
    }
  }

  // if there are free capacities to any side, fill them with the best 3 extends
  //  but only if no full-length overlap was taken with a rail

  if(!takenrailfulllength){
    CEBUG("check freecap\n");
    if(numleftext>3) numleftext=3;
    if(numrightext>3) numrightext=3;
    for(auto & tmwe : cfhd.tmpmatchwith){
      if(numleftext==0 && numrightext==0) break;
      if(tmwe.taken==false){
	CEBUG("tmwI: " << tmwe << endl);
	adse.calcNewEstimateFromSkim(
	  tmwe.eoffset,
	  (*SKIM3_readpool)[actreadid].getLenClippedSeq(),
	  (*SKIM3_readpool)[tmwe.otherid].getLenClippedSeq(),
	  actreadid,
	  tmwe.otherid,
	  1,
	  direction);
	if(numleftext>0 && adse.getEstimatedLeftExpand(actreadid)>0){
	  tmwe.taken=true;
	  --numleftext;
	  CEBUG("************* take free capl " << numleftext << '\n');
	}
	if(!tmwe.taken && numrightext>0 && adse.getEstimatedRightExpand(actreadid)>0) {
	  tmwe.taken=true;
	  --numrightext;
	  CEBUG("************* take free capr " << numrightext << '\n');
	}
      }
    }
  }

  std::ofstream logfout;
  if(SKIM3_logflag_save2){
    std::string path,justfilename;
    if(direction>0){
      splitFullPathAndFileName(SKIM3_posfmatchfname,path,justfilename);
    }else{
      splitFullPathAndFileName(SKIM3_poscmatchfname,path,justfilename);
    }
    std::string logfilename=path+"/elog.skim.save2."+justfilename;
    logfout.open(logfilename, std::ios::out|std::ios::app|std::ios::ate);
  }

  // we made our choice, now save that
  for(auto & tmwe : cfhd.tmpmatchwith){
    if(tmwe.taken){
      if(SKIM3_logflag_save2){
	logfout << "taken:\t" << SKIM3_readpool->getRead(tmwe.otherid).getName()
		<< '\t' << SKIM3_readpool->getRead(actreadid).getName()
		<< '\t' << tmwe;
      }
      if(cfhd.shfsv.size()==cfhd.shfsv.capacity()){
	boost::mutex::scoped_lock lock(SKIM3_resultfileoutmutex);
	cfhd.posmatchfout->write(reinterpret_cast<char*>(&cfhd.shfsv[0]),sizeof(skimhitforsave_t)*cfhd.shfsv.size());
	if(cfhd.posmatchfout->bad()){
	  MIRANOTIFY(Notify::FATAL, "Could not write anymore to skimhit save5. Disk full? Changed permissions?");
	}
	cfhd.shfsv.clear();
      }
      cfhd.shfsv.resize(cfhd.shfsv.size()+1);
      skimhitforsave_t & shfs=cfhd.shfsv.back();
      shfs.rid1=tmwe.otherid;
      shfs.rid2=actreadid;
      shfs.eoffset=-(tmwe.eoffset);
      shfs.percent_in_overlap=tmwe.percent_in_overlap;
      shfs.numhashes=tmwe.numhashes;
      shfs.ol_stronggood  =tmwe.ol_stronggood;
      shfs.ol_weakgood    =tmwe.ol_weakgood;
      shfs.ol_belowavgfreq=tmwe.ol_belowavgfreq;
      shfs.ol_norept      =tmwe.ol_norept;
      shfs.ol_rept        =tmwe.ol_rept;

      CEBUG("save2: " << actreadid << "\n" << tmwe);

      cfhd.ridswithmatches.push_back(std::min(tmwe.otherid,actreadid));
    }else{
      if(SKIM3_logflag_save2){
	logfout << "dropped:\t" << SKIM3_readpool->getRead(tmwe.otherid).getName()
		<< '\t' << SKIM3_readpool->getRead(actreadid).getName()
		<< '\t' << tmwe;
      }
    }
  }

  cfhd.ridswithmatches.clear();

  CEBUG("end selectPotentialHitsForSave2()\n");

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Updates SKIM3_overlapcritlevell/r
 *
 * also finishes calculation of elements not initialised yet in the
 *  matchwithsorter_t vector (not done previously in checkForPotentialHits()
 *  to get down the number of calls to ADSEstimator calculations)
 *
 *  - stores the estimated score in estimscore
 *  - calculates ol_fulllength and ol_fullencased
 *
 * BaCh 13.06.2016: do NOT update critlevels for rail reads, that simply
 *   will lead to wrong mappings in some highly repetitive cases.
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void Skim<TVHASH_T>::updateCriterionLevels(const int8 direction, const uint32 actreadid, cfh_threaddata_t & cfhd)
{
  FUNCSTART("void Skim<TVHASH_T>::updateCriterionLevels(const int8 direction, const uint32 actreadid, std::vector<matchwithsorter_t> & tmpmatchwith)");

  CEBUG("start updateCriterionLevels() " << actreadid << '\n');

  ADSEstimator adse;

  for(uint32 ocvi=0; ocvi<2; ++ocvi){
    cfhd.uidswithnewcritlevelvl[ocvi].clear();
    cfhd.uidswithnewcritlevelvr[ocvi].clear();
    cfhd.critlevellofnewuidsv[ocvi].clear();
    cfhd.critlevelrofnewuidsv[ocvi].clear();
  }

  uint8 actreadcritlevell[2]={255,255};
  uint8 actreadcritlevelr[2]={255,255};
  uint8 tmpnewcll=255;
  uint8 tmpnewclr=255;

  for(auto & tmwe : cfhd.tmpmatchwith) {
    adse.calcNewEstimateFromSkim(
      tmwe.eoffset,
      (*SKIM3_readpool)[actreadid].getLenClippedSeq(),
      (*SKIM3_readpool)[tmwe.otherid].getLenClippedSeq(),
      actreadid,
      tmwe.otherid,
      1,
      direction);

    CEBUG("ADSE: " << adse << endl);

    CEBUG("tmwe before: " << tmwe << endl);

    //cout << "EstimO: " << adse.getEstimatedOverlap();
    tmwe.estimscore=adse.getEstimatedOverlap()*tmwe.percent_in_overlap*tmwe.percent_in_overlap;
    if(adse.getContainmentLevel()>0) {
      tmwe.ol_fulllength=true;
      if(tmwe.percent_in_overlap==100 && !tmwe.ol_rept
    	 && abs(static_cast<int32>(SKIM3_readpool->getRead(actreadid).getLenClippedSeq())
    		- static_cast<int32>(SKIM3_readpool->getRead(tmwe.otherid).getLenClippedSeq())) >= 8){
    	tmwe.ol_fullencased=true;
      }
    }

    CEBUG("tmwe after: " << tmwe << endl);

    uint32 ocvi=0;  // overlap criterion vector index 0 is for norept overlaps
    if(!tmwe.ol_norept) ocvi=1;  // index 1 is for all other rept overlaps (among them the rept)

    if(actreadcritlevell[ocvi]>0 || actreadcritlevelr[ocvi]>0){
      CEBUG("old: " << static_cast<uint16>(actreadcritlevell[ocvi]) << " " << static_cast<uint16>(actreadcritlevelr[ocvi]) << '\n');
      getOverlapCriterionLevel(actreadid,
			       SKIM3_readpool->getRead(actreadid).getSequencingType(),
			       adse,
			       static_cast<uint8>(tmwe.percent_in_overlap),
			       tmpnewcll,tmpnewclr);
      if(tmpnewcll<actreadcritlevell[ocvi]) actreadcritlevell[ocvi]=tmpnewcll;
      if(tmpnewclr<actreadcritlevelr[ocvi]) actreadcritlevelr[ocvi]=tmpnewclr;
      CEBUG("new: " << static_cast<uint16>(actreadcritlevell[ocvi]) << " " << static_cast<uint16>(actreadcritlevelr[ocvi]) << '\n');
    }

    // then for other id
    getOverlapCriterionLevel(tmwe.otherid,SKIM3_readpool->getRead(tmwe.otherid).getSequencingType(),
			     adse,
			     static_cast<uint8>(tmwe.percent_in_overlap),
			     tmpnewcll,tmpnewclr);
    if(tmpnewcll<(*SKIM3_overlapcritlevelvl)[ocvi][tmwe.otherid]){
      cfhd.uidswithnewcritlevelvl[ocvi].push_back(tmwe.otherid);
      cfhd.critlevellofnewuidsv[ocvi].push_back(tmpnewcll);
      CEBUG("tmpnewcll: pushback " << cfhd.uidswithnewcritlevelvl[ocvi].back() << "\t" << static_cast<uint16>(cfhd.critlevellofnewuidsv[ocvi].back()) << '\n');
    }
    if(tmpnewclr<(*SKIM3_overlapcritlevelvr)[ocvi][tmwe.otherid]){
      cfhd.uidswithnewcritlevelvr[ocvi].push_back(tmwe.otherid);
      cfhd.critlevelrofnewuidsv[ocvi].push_back(tmpnewclr);
      CEBUG("tmpnewclr: pushback " << cfhd.uidswithnewcritlevelvr[ocvi].back() << "\t" << static_cast<uint16>(cfhd.critlevelrofnewuidsv[ocvi].back()) << '\n');
    }
  }

  {
    boost::mutex::scoped_lock lock(SKIM3_critlevelwrite_mutex);
    if(!SKIM3_readpool->getRead(actreadid).isRail()){
      for(uint32 ocvi=0; ocvi<2; ++ocvi){
	if(actreadcritlevell[ocvi]<(*SKIM3_overlapcritlevelvl)[ocvi][actreadid]){
	  CEBUG("ari: update critlevell " << actreadid << " " << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[ocvi][actreadid]));
	  (*SKIM3_overlapcritlevelvl)[ocvi][actreadid]=actreadcritlevell[ocvi];
	  CEBUG(" to " << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[ocvi][actreadid]) << '\n');
	}
	if(actreadcritlevelr[ocvi]<(*SKIM3_overlapcritlevelvr)[ocvi][actreadid]){
	  CEBUG("ari: update critlevelr " << actreadid << " " << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[ocvi][actreadid]));
	  (*SKIM3_overlapcritlevelvr)[ocvi][actreadid]=actreadcritlevelr[ocvi];
	  CEBUG(" to " << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[ocvi][actreadid]) << '\n');
	}
      }
    }

    for(uint32 ocvi=0; ocvi<2; ++ocvi){
      auto nI=cfhd.uidswithnewcritlevelvl[ocvi].cbegin();
      auto cI=cfhd.critlevellofnewuidsv[ocvi].cbegin();
      for(; nI != cfhd.uidswithnewcritlevelvl[ocvi].cend(); ++nI, ++cI){
	// still check ... might have changed in the mean time by another thread
	if(!SKIM3_readpool->getRead(*nI).isRail()
	   && *cI < (*SKIM3_overlapcritlevelvl)[ocvi][*nI]){
	  CEBUG("ori: update critlevell " << *nI <<  " " << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[ocvi][*nI]));
	  (*SKIM3_overlapcritlevelvl)[ocvi][*nI]=*cI;
	  CEBUG(" to " << static_cast<uint16>((*SKIM3_overlapcritlevelvl)[ocvi][*nI]) << '\n');
	}
      }
      nI=cfhd.uidswithnewcritlevelvr[ocvi].cbegin();
      cI=cfhd.critlevelrofnewuidsv[ocvi].cbegin();
      for(; nI != cfhd.uidswithnewcritlevelvr[ocvi].cend(); ++nI, ++cI){
	// still check ... might have changed in the mean time by another thread
	if(!SKIM3_readpool->getRead(*nI).isRail()
	   && *cI < (*SKIM3_overlapcritlevelvr)[ocvi][*nI]){
	  CEBUG("ori: update critlevelr " << *nI <<  " " << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[ocvi][*nI]));
	  (*SKIM3_overlapcritlevelvr)[ocvi][*nI]=*cI;
	  CEBUG(" to " << static_cast<uint16>((*SKIM3_overlapcritlevelvr)[ocvi][*nI]) << '\n');
	}
      }
    }
  }

  CEBUG("end updateCriterionLevels()\n");

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * Returns criterium levels for overlaps extending left and right
 *
 *
 *   Sanger, 454, Ion, Text, PacBio
 *     crit0: 80% overlap
 *     crit1: 70% overlap
 *     crit2: 60% overlap
 *     crit3: 50% overlap
 *   Solexa:
 *     crit-level from 0 to 59, see code
 *
 * Special: if none of above, level = 240 if fully encased
 *          made for hits against rails, so that partial matches
 *          are not preferred to fully encased matches (which
 *          hits against a backbone should normally be)
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {if(adse.getID1()==2 && adse.getID2()==140820) {cout << bla; cout.flush();}}

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

template<typename TVHASH_T>
void Skim<TVHASH_T>::getOverlapCriterionLevel(const uint32 actreadid, const uint8 seqtype, const ADSEstimator & adse, const uint8 relscore, uint8 & levell, uint8 & levelr)
{
  FUNCSTART("void Skim<TVHASH_T>::getOverlapCriterionLevel(const uint32 actreadid, const uint8 seqtype, const ADSEstimator & adse, const uint8 relscore, uint8 & levell, uint8 & levelr)");

  uint32 overlapratiopc=100*adse.getEstimatedOverlap()/adse.getLen(actreadid);

  CEBUG(adse);
  CEBUG("gOCL: " << actreadid << " " << static_cast<uint16>(relscore) << " " << adse.getEstimatedOverlap() << " " << adse.getLen(actreadid) << " " << overlapratiopc);

  levell=255;
  levelr=255;
  switch(seqtype){
  case ReadGroupLib::SEQTYPE_SOLEXA : {
    if(relscore>=95){
      // level from 0 to 179
      uint8 startlevel=(100-relscore)*30;;
      startlevel+=29-(adse.getEstimatedOverlap()*29/adse.getLen(actreadid));
      BUGIFTHROW(startlevel>=200,"Startlevel>=200? " << static_cast<uint16>(startlevel) << "\tid1: " << adse.getID1() << " id2: " << adse.getID2() << "\n" << adse);

      // this is for reads longer than 100bp (MiSeq), because this is really "good enough"
      if(adse.getEstimatedOverlap()>=100 && relscore==100) startlevel=0;

      if(adse.getEstimatedLeftExpand(actreadid)>0){
	levell=startlevel;
      }
      if(adse.getEstimatedRightExpand(actreadid)>0){
	levelr=startlevel;
      }
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_SANGER :
  case ReadGroupLib::SEQTYPE_454GS20 :
  case ReadGroupLib::SEQTYPE_IONTORRENT :
  case ReadGroupLib::SEQTYPE_PACBIOHQ :
  case ReadGroupLib::SEQTYPE_PACBIOLQ :
  case ReadGroupLib::SEQTYPE_TEXT :{
    if(overlapratiopc>=80){
      if(adse.getEstimatedLeftExpand(actreadid)>0){
	levell=0;
      }
      if(adse.getEstimatedRightExpand(actreadid)>0){
	levelr=0;
      }
    }else if(overlapratiopc>=70){
      if(adse.getEstimatedLeftExpand(actreadid)>0){
	levell=1;
      }
      if(adse.getEstimatedRightExpand(actreadid)>0){
	levelr=1;
      }
    }else if(overlapratiopc>=60){
      if(adse.getEstimatedLeftExpand(actreadid)>0){
	levell=2;
      }
      if(adse.getEstimatedRightExpand(actreadid)>0){
	levelr=2;
      }
    }else if(overlapratiopc>=50){
      if(adse.getEstimatedLeftExpand(actreadid)>0){
	levell=3;
      }
      if(adse.getEstimatedRightExpand(actreadid)>0){
	levelr=3;
      }
    }
    break;
  }
  default : {
    BUGIFTHROW(true,"Unknown/unhandled seqtype " << static_cast<uint16>(seqtype));
  }
  }

  if(levell == 255 && levelr==255){
    uint32 clevel=adse.getContainmentLevel();
    if(clevel >1
       || (clevel==1 && adse.getIDOfContained()==actreadid)){
      levell=240;
      levelr=240;
    }
  }

  CEBUG("\tR: " << static_cast<uint16>(levell) << " " << static_cast<uint16>(levelr) << '\n');

  FUNCEND();
}
//#define CEBUG(bla)
//#define CEBUG2(bla)














/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG_extra_cFPH

template<typename TVHASH_T>
void Skim<TVHASH_T>::checkForPotentialHits(const int8 direction, const uint32 actreadid, std::vector<matchwithsorter_t> & tmpmatchwith, std::vector<readhashmatch_t> & readhashmatches, std::vector<uint32> & smallhist4repeats)
{
  //bool dodebug=false;

  tmpmatchwith.clear();

  // readhashmatches should not be empty ... normally.
  // but new method to deal with megahubs reduces this vector, keeping only
  //  'approximately normal' frequencies. Which in turn means: some vectors
  //  might be completely emptied
  // so, if it is empty, return immediately
  if(readhashmatches.empty()) return;

  // use non-parallel sort as we are already in a multithreading environment here
  mstd::ssort(readhashmatches, sortreadhashmatch_t_);

  bool actreadisrail=SKIM3_readpool->getRead(actreadid).isRail();
  bool actreadhasenough=SKIM3_nomorehitseie[actreadid]>0;

  uint32 possiblehits=0;
  uint32 acceptedhits=0;

  auto sI=readhashmatches.cbegin();
  uint32 countid=sI->rid2;
  while(sI != readhashmatches.cend()){
    uint32 rid2=sI->rid2;

    // disregard this potential match if
    //  1) both reads are rails
    if((actreadisrail && SKIM3_readpool->getRead(rid2).isRail())){
      for(;sI!=readhashmatches.cend() && sI->rid2==rid2; ++sI);
      if(sI!=readhashmatches.cend()) countid=sI->rid2;
      continue;
    }
    //  2) we scan only against rails and both reads are non-rails
    if((SKIM3_onlyagainstrails
	&& (!actreadisrail && !SKIM3_readpool->getRead(rid2).isRail()))){
      for(;sI!=readhashmatches.cend() && sI->rid2==rid2; ++sI);
      if(sI!=readhashmatches.cend()) countid=sI->rid2;
      continue;
    }
    // 3) both reads have the "enough-is-enough" status
    // nope, don't do that here! Do it when we know whether a strong match exists

    //  4) this read pair has been banned previously
    if((*SKIM3_bannedoverlaps).checkIfBanned(actreadid,rid2)){
      for(;sI!=readhashmatches.cend() && sI->rid2==rid2; ++sI);
      if(sI!=readhashmatches.cend()) countid=sI->rid2;
//      {
//	boost::mutex::scoped_lock lock(SKIM3_globalclassdatamutex);
//	SKIM3_totalpermbans++;
//      }
      continue;
    }

    ++possiblehits;
    if(possiblehits==1){
      CEBUG("Potential hits of " << actreadid << " (" << static_cast<int16>(direction) << '/' << SKIM3_readpool->getRead(actreadid).getLenClippedSeq() << ")\n----------------\n");
      CEBUG(SKIM3_readpool->getRead(actreadid) << endl);
      CEBUG("----------------\n");
    }

    uint16 oldhashpos=sI->hashpos1;
    uint16 hp1min=0xffff;
    uint16 hp1max=0;

    uint16 hp2min=0xffff;
    uint16 hp2max=0;

    int32  eoffsetmin=0x7fffffff;
    int32  eoffsetmax=0x80000000;
    int32  oldeoffset=sI->eoffset;

    int32  maxeoffsetjump=0;
    int32  weighteoffsetjumps=0;

    uint32 numhashes=0;

    bool flag_norept=true;
    size_t totalfreq3counter=0;
    size_t totalfreq5counter=0;
    size_t contiguousfreq3counter=0;
    size_t maxcontiguousfreq3counter=0;
    size_t contiguousfreq32counter=0;
    size_t maxcontiguousfreq32counter=0;


//    if((actreadid==52053 && rid2==208673)
//       || (rid2==52053 && actreadid==208673)) dodebug=true;

//    if(actreadid==0 || rid2==0) dodebug=true;
//#define CEBUG(bla)   {if(dodebug) cout << bla; cout.flush();}


    if(sI != readhashmatches.cend()){
      CEBUG("hhh " << SKIM3_readpool->getRead(actreadid).getName() << " " << SKIM3_readpool->getRead(sI->rid2).getName() << endl);
    }

    auto sIS=sI; // save current start, sIS used after after loop!
    for(;sI != readhashmatches.cend() && sI->rid2 == countid; ++sI){
      CEBUG(*sI);

      // this ensures that the eoffset between two following
      //  entries may not differ by too much (10 bases here)
      // IF they do, then this is treated like a different hit
      //  by breaking the loop

      // TODO: 100 only as tes tfor PacBio CLR
      if(abs(sI->eoffset - oldeoffset) > 100){
	CEBUG("BREAKER!\n");
	break;
      }
      ++numhashes;
      CEBUG("numhashes: " << numhashes << '\n');

      if(oldhashpos + SKIM3_hashsavestepping != sI->hashpos1){
	CEBUG("NOT CONTIGUOUS!\n");
	maxcontiguousfreq3counter=std::max(maxcontiguousfreq3counter,contiguousfreq3counter);
	maxcontiguousfreq32counter=std::max(maxcontiguousfreq32counter,contiguousfreq32counter);
	contiguousfreq3counter=0;
	contiguousfreq32counter=0;
	maxeoffsetjump=std::max(maxeoffsetjump,abs(abs(sI->eoffset)-abs(oldeoffset)));
	weighteoffsetjumps+=abs(abs(sI->eoffset)-abs(oldeoffset));
      }
      hp1min=std::min(hp1min,sI->hashpos1);
      hp1max=std::max(hp1max,sI->hashpos1);
      eoffsetmin=std::min(eoffsetmin,sI->eoffset);
      eoffsetmax=std::max(eoffsetmax,sI->eoffset);
      oldeoffset=sI->eoffset;

      hp2min=std::min(hp2min,sI->hashpos2);
      hp2max=std::max(hp2max,sI->hashpos2);

      if(sI->bhashstats.getFrequency() >= 5){
	totalfreq5counter++;
	flag_norept=false;
	contiguousfreq3counter=0;
	contiguousfreq32counter=0;
      }else if(sI->bhashstats.getFrequency() > 3){
	flag_norept=false;
	contiguousfreq3counter=0;
	contiguousfreq32counter=0;
      }else if(sI->bhashstats.getFrequency() == 3){
	contiguousfreq3counter++;
	maxcontiguousfreq3counter=std::max(maxcontiguousfreq3counter,contiguousfreq3counter);
	totalfreq3counter++;

	contiguousfreq32counter++;
	maxcontiguousfreq32counter=std::max(maxcontiguousfreq32counter,contiguousfreq32counter);
      }else if(sI->bhashstats.getFrequency() == 2){
	contiguousfreq32counter++;
	maxcontiguousfreq32counter=std::max(maxcontiguousfreq32counter,contiguousfreq32counter);
      }else{
	contiguousfreq3counter=0;
	contiguousfreq32counter=0;
      }


#ifdef CEBUG_extra_cFPH
      {
	boost::mutex::scoped_lock lock(SKIM3_coutmutex);
	CEBUG(sI->rid2
	      << "\t" << SKIM3_readpool->getRead(sI->rid2).getName()
	      << "\t" << SKIM3_readpool->getRead(sI->rid2).getLenClippedSeq()
	      << "\t" << sI->eoffset
	      << "\t" << sI->hashpos1
	      << "\t" << oldhashpos
	      << "\tfq: " << static_cast<uint16>(sI->bhashstats.getFrequency())
	      << "\t" << flag_norept
	      << ' ' << contiguousfreq3counter
	      << ' ' << maxcontiguousfreq3counter
	      << ' ' << contiguousfreq32counter
	      << ' ' << maxcontiguousfreq32counter
	      << ' ' << totalfreq3counter
	      << ' ' << totalfreq5counter
	      //<< "\t" << sI->hashpos2
	      << '\n');
      }
#endif

      oldhashpos=sI->hashpos1;
    }

    int32 maxoverlap;

    // adjust min positions for the hash length
    hp1min-=(SKIM3_basesperhash-1);
    hp2min-=(SKIM3_basesperhash-1);

    int32 eoffsetmean=eoffsetmin+(eoffsetmax-eoffsetmin)/2;

    // calc max overlap
    // currently only for one offset
    if(eoffsetmean<0){
      maxoverlap=std::min(SKIM3_readpool->getRead(rid2).getLenClippedSeq()+eoffsetmean,SKIM3_readpool->getRead(actreadid).getLenClippedSeq());
    }else{
      maxoverlap=std::min(SKIM3_readpool->getRead(actreadid).getLenClippedSeq()-eoffsetmean,SKIM3_readpool->getRead(rid2).getLenClippedSeq());
    }

    // correct the maxoverlap by the modulo of the hash steps as the
    //  border hashes will be found only in 1/(hash stepping) cases
    maxoverlap=maxoverlap-(maxoverlap%SKIM3_hashsavestepping);

    // hashe3soverlap is not the number of hashes in the overlap,
    // but the length of the overlap
    int32 hashesoverlap=hp1max-hp1min+1;

    int32 perc=100*hashesoverlap/maxoverlap;
    bool disregardperc=false;

    int32 minpercentrequired=std::min(
      SKIM3_percentrequired[SKIM3_readpool->getRead(actreadid).getSequencingType()],
      SKIM3_percentrequired[SKIM3_readpool->getRead(rid2).getSequencingType()]);

    uint32 maxnumhashes=((maxoverlap-SKIM3_basesperhash)/SKIM3_hashsavestepping)+1;

    CEBUG(static_cast<int16>(direction) << "\tmo: " << maxoverlap << "\tperc: " << perc << "\tari: " << actreadid << "\trid2: " << rid2 << "\tnumh: " << numhashes << "\tmnh: " << maxnumhashes << "\teom: " << eoffsetmean << "\teomin: " << eoffsetmin << "\teomax: " << eoffsetmax << "\tmej: " << maxeoffsetjump << endl);

    bool majorrecalc=false;

    if(maxnumhashes>0){

      if(SKIM3_readpool->getRead(actreadid).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)
	 || SKIM3_readpool->getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)){
	majorrecalc=true;
	// TODO: PacBio perhaps in recalc below, take also left & right of max
	//  (or monotonic down?) as numhashes
	// via flag?
      }else if(perc>=minpercentrequired
	 && numhashes>1){
	if((perc == 100 && maxeoffsetjump>=3)
	   || numhashes>maxnumhashes
	   // NO!!! this would be bad for microrepeats || weighteoffsetjumps>=3
	  ) {
	  majorrecalc=true;
	}
      }
      if(majorrecalc){
	CEBUG("part1\n");
	// find eoffset with most hashes (this is our new eoffset mean)
	//  and numhashes will now only concern that offset (disregarding frameshifts,
	//  but that cannot be helped)
	smallhist4repeats.clear();
	smallhist4repeats.resize(eoffsetmax-eoffsetmin+1,0);
	numhashes=0;
	for(auto rI=sIS; rI != sI; ++rI){
	  //cout << *rI;
	  //if(abs(rI->eoffset-eoffsetmin)>=smallhist4repeats.size()){
	  //  cout << "Eh, what?\n";
	  //  exit(0);
	  //}
	  ++smallhist4repeats[rI->eoffset-eoffsetmin];
	  if(smallhist4repeats[rI->eoffset-eoffsetmin] > numhashes){
	    numhashes=smallhist4repeats[rI->eoffset-eoffsetmin];
	  }
	}

//	cout << "hist\n";
//	for(uint32 ii=0; ii<smallhist4repeats.size(); ++ii){
//	  cout << ii << "\t" << smallhist4repeats[ii] << endl;
//	}

	int32 newmini=0;
	for(; newmini<smallhist4repeats.size(); ++newmini){
	  if(smallhist4repeats[newmini] == numhashes){
	    break;
	  }
	}
	int32 newmaxi=newmini;
	for(; newmaxi<smallhist4repeats.size(); ++newmaxi){
	  if(smallhist4repeats[newmaxi] != numhashes){
	    break;
	  }
	}
	--newmaxi;
	//cout << "mini: " << newmini << "\tmaxi: " << newmaxi << endl;
	newmini+=eoffsetmin;
	newmaxi+=eoffsetmin;
	//cout << "mini: " << newmini << "\tmaxi: " << newmaxi << endl;
	eoffsetmean=(newmaxi+newmini)/2;

	// recalc hp1min/max
	hp1min=0xffff;
	hp1max=0;

	for(auto rI=sIS;rI != sI; ++rI){
	  if(rI->eoffset == eoffsetmean){
	    hp1min=std::min(hp1min,rI->hashpos1);
	    hp1max=std::max(hp1max,rI->hashpos1);
	  }
	}
	// adjust min positions for the hash length
	hp1min-=(SKIM3_basesperhash-1);

	hashesoverlap=hp1max-hp1min+1;
	perc=100*hashesoverlap/maxoverlap;

	if(perc==100){
	  eoffsetmin=eoffsetmean;
	  eoffsetmax=eoffsetmean;
	}

	//cout << static_cast<int16>(direction) << "\tari: " << actreadid << "\trid2: " << rid2 << "\tnumh: " << numhashes << "\teom: " << eoffsetmean
	//     << "\tho: " << hashesoverlap << "\tper: " << perc
	//     << endl;

	// saver:
	// we might have been too harsh, very probably so in mapping asemblies for the following case:
	//
	//  C     ..............................
	//  R          ........*.....
	//
	// where only ~half of the hashes in R are counted. In that case, set percent to the
	//  minimum needed to not be thrown out.
	if(perc<minpercentrequired){
	  disregardperc=true;
	}

	CEBUG("recalc\n" << static_cast<int16>(direction) << "\tperc: " << perc << "\tari: " << actreadid << "\trid2: " << rid2 << "\tnumh: " << numhashes << "\tmnh: " << maxnumhashes << "\teom: " << eoffsetmean << "\teomin: " << eoffsetmin << "\teomax: " << eoffsetmax << endl);
      }

      if(SKIM3_readpool->getRead(actreadid).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)
	 && SKIM3_readpool->getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)){
	// TODO: PacBio make 20 configurable
	// PB / PB hits need at least 20 hashes
	if(numhashes<20){
	  numhashes=0;
	  maxnumhashes=0;
	  perc=0;
	  disregardperc=false;
	}
      }

      // look a bit closer at potential perfect matches
      if(perc == 100){
	if(eoffsetmin != eoffsetmax){
	  // this could not be: have at least two different expected offsets
	  //  and a 100% coverage. Side effects from intra-read repeats
	  //  or intre-read indel
	  // therefore, make sure this does not get through as a 100% match
	  //  by using the number of hashes as percentage
	  //  but do not fall below minimum required
	  //perc=99;
	  //cout << "dida\n";
	  perc=100*numhashes/maxnumhashes;
	  if(perc<minpercentrequired) disregardperc=true;
	}else if((numhashes-1)*SKIM3_hashsavestepping+SKIM3_basesperhash < maxoverlap){
	  // maxoverlap covers the whole potential overlap, but
	  //  there are not enough hashes supporting for 100% match
	  //  (base mismatch somewhere)
	  // reduce the percentage to show it's not a perfect match
	  //  but do not fall below minimum required
	  //perc=99;
	  //cout << "dide\n";
	  perc=100*numhashes/maxnumhashes;
	  if(perc<minpercentrequired) disregardperc=true;
	}
      }else if(eoffsetmin == eoffsetmax){
	if(perc>100) {
	  perc=100;
	}else{
	  if(perc>=minpercentrequired && numhashes==maxnumhashes){
	    CEBUG("maxnumhashes 100% saver: "  << perc << '\n');
	    perc=100;
	  }
	}
      }

      CEBUG("after closer look\n" << static_cast<int16>(direction) << "\tperc: " << perc << "\tari: " << actreadid << "\trid2: " << rid2 << "\tnumh: " << numhashes << "\tmnh: " << maxnumhashes << "\teom: " << eoffsetmean << "\teomin: " << eoffsetmin << "\teomax: " << eoffsetmax << endl);

      // if rail and only partial match -> reduce percentage
      if(!majorrecalc
	 && SKIM3_readpool->getRead(rid2).isRail()
	 && maxnumhashes+SKIM3_basesperhash*SKIM3_hashsavestepping<SKIM3_readpool->getRead(actreadid).getLenClippedSeq()){
	  if(perc>=minpercentrequired) disregardperc=true;
	  perc=100*numhashes/maxnumhashes;

	  CEBUG("partial match recalc\n" << static_cast<int16>(direction) << "\tperc: " << perc << "\tari: " << actreadid << "\trid2: " << rid2 << "\tnumh: " << numhashes << "\tmnh: " << maxnumhashes << "\teom: " << eoffsetmean << "\teomin: " << eoffsetmin << "\teomax: " << eoffsetmax << endl);

      }

      // reduce percentage by 3*largest gap we have encountered
      // and by number of gaps encountered
      if(perc>=minpercentrequired){
	perc-=3*maxeoffsetjump;
	if(!majorrecalc) perc-=weighteoffsetjumps;
	if(perc<=0) perc=1;
	if(perc<minpercentrequired) disregardperc=true;

	CEBUG("gap recalc\n" << static_cast<int16>(direction) << "\tperc: " << perc << "\tari: " << actreadid << "\trid2: " << rid2 << "\tnumh: " << numhashes << "\tmnh: " << maxnumhashes << "\teom: " << eoffsetmean << "\teomin: " << eoffsetmin << "\teomax: " << eoffsetmax << endl);
      }

      int32 minoverlaprequired=std::min(
	SKIM3_overlaplenrequired[SKIM3_readpool->getRead(actreadid).getSequencingType()],
	SKIM3_overlaplenrequired[SKIM3_readpool->getRead(rid2).getSequencingType()]);

#ifdef CEBUG_extra_cFPH
      {
	boost::mutex::scoped_lock lock(SKIM3_coutmutex);
	CEBUG("eomin: " << eoffsetmin << "\teomax: " << eoffsetmax
	      << "\tmor: " << minoverlaprequired
	      << "\tho: " << hashesoverlap
	      << "\t%: " << perc
	      << "\t%<: " << minpercentrequired << endl);
      }
#endif

      CEBUG("disr%: "  << disregardperc << '\n');
      CEBUG("ho: " << hashesoverlap << "\tmo required: " << minoverlaprequired << endl);

      if(perc>100) perc=100;

      // we take the hit if the overlap percentage is above threshold
      // NEW: or if both reads have MNRr tags
      // NEW: or we have a perfect overlap >= 17
      //   || (perc==100 && hashesoverlap>=17))
      // BaCh: 06.02.2015; nope, not "or >= 17" as this does not reflect the wish of the user!
      //       and it needlessly inflates the skim table for highly repetitive data sets
      if(hashesoverlap >= minoverlaprequired
	 && (perc>=minpercentrequired
	     || disregardperc
	     || (SKIM3_hasMNRr[actreadid] && SKIM3_hasMNRr[rid2]))) {

	CEBUG("accepted\n");
	acceptedhits++;

	// increase overlapcounter only for "real" reads,
	//  not for rails
	if(!SKIM3_readpool->getRead(actreadid).isRail()
	   && !SKIM3_readpool->getRead(rid2).isRail()){
	  boost::mutex::scoped_lock lock(SKIM3_globalclassdatamutex);
	  (*SKIM3_overlapcounter)[actreadid]+=1;
	  (*SKIM3_overlapcounter)[rid2]+=1;
	}

//#define CEBUG(bla)   {if(actreadid==273252 && rid2==273250) cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla; cout.flush();}
	matchwithsorter_t tmp;
	tmp.otherid=rid2;
	tmp.eoffset=eoffsetmean;

	tmp.percent_in_overlap=perc;
	tmp.numhashes=numhashes;
	tmp.estimscore=0;
	tmp.taken=false;

	// this was the standard, works good but sometimes too harsh
	// tmp.ol_belowavgfreq=maxcontiguousfreq32counter > (SKIM3_basesperhash/SKIM3_hashsavestepping-1);
	// tmp.ol_weakgood=maxcontiguousfreq3counter > (SKIM3_basesperhash/SKIM3_hashsavestepping-1);
	// this is worse in performance on lpla synthetic data
	//  tmp.ol_weakgood=maxcontiguousfreq3counter > 1;

	tmp.ol_belowavgfreq=false;
	if(maxcontiguousfreq32counter){
	  tmp.ol_belowavgfreq=(SKIM3_basesperhash+(maxcontiguousfreq32counter-1)*SKIM3_hashsavestepping) >= 26;
	}
	if(maxcontiguousfreq3counter){
	  tmp.ol_weakgood=(SKIM3_basesperhash+(maxcontiguousfreq3counter-1)*SKIM3_hashsavestepping) >= 20;
	  tmp.ol_stronggood=(SKIM3_basesperhash+(maxcontiguousfreq3counter-1)*SKIM3_hashsavestepping) >= SKIM3_basesperhash*2-1;
	}else{
	  tmp.ol_weakgood=false;
	  tmp.ol_stronggood=false;
	}

	// small hack for "very good reads" (atm Illumina / Moleculo or other "corrected" reads)
	if(perc==100 && SKIM3_basesperhash>60 && flag_norept){
	  if(SKIM3_basesperhash>90){
	    tmp.ol_stronggood=true;
	  }else{
	    tmp.ol_weakgood=true;
	  }
	}


	//tmp.ol_stronggood=tmp.ol_weakgood & (maxcontiguousfreq3counter-1 >= SKIM3_basesperhash*2/SKIM3_hashsavestepping);

	/* still too harsh sometimes

	// if no strong good and one of the reads is Solexa:
	//  extra rule for short (<64 bases) Solexa: strong good also for contiguous
	//  overlaps of >= 30
	if(!tmp.ol_stronggood){
	  bool acceptshort=false;
	  if(SKIM3_readpool->getRead(actreadid).getSequencingType() == ReadGroupLib::SEQTYPE_SOLEXA
	     && SKIM3_readpool->getRead(actreadid).getLenClippedSeq() < 64) {
	    acceptshort=true;
	  }else if(SKIM3_readpool->getRead(rid2).getSequencingType() == ReadGroupLib::SEQTYPE_SOLEXA
		   && SKIM3_readpool->getRead(rid2).getLenClippedSeq() < 64) {
	    acceptshort=true;
	  }
	  if(acceptshort){
	    tmp.ol_stronggood=(SKIM3_basesperhash+(maxcontiguousfreq3counter-1)*SKIM3_hashsavestepping) >= 30;
	    CEBUG("Should accept short: " << SKIM3_basesperhash+(maxcontiguousfreq3counter-1)*SKIM3_hashsavestepping << ' ' << tmp.ol_stronggood << '\n');
	  }
	}
	*/

	/* Nope, too lenient
	// all contiguous HAF33 overlaps >= 30 are stronggood
	if(!tmp.ol_stronggood){
	  tmp.ol_stronggood=(SKIM3_basesperhash+(maxcontiguousfreq3counter-1)*SKIM3_hashsavestepping) >= 30;
	}
	*/

	// 3) save only if there's a strong overlap or if both reads are not enough-is-enough together
	if(tmp.ol_stronggood || !(actreadhasenough && SKIM3_nomorehitseie[rid2]!=0)){

	  tmp.ol_norept=flag_norept;
	  tmp.ol_rept=(totalfreq5counter > 0);

	  // calc of these two would not be 100% accurate here, moved to updateCriterionLevels()
	  tmp.ol_fulllength=false;
	  tmp.ol_fullencased=false;

	  tmpmatchwith.push_back(tmp);

	  if(tmp.ol_rept) {CEBUG("\nREPT!!!\n")};
	  CEBUG("Pushing possible hit with offset: " << tmp.eoffset << endl
		<< rid2
		<< "\t" << actreadid
		<< "\t" << SKIM3_readpool->getRead(rid2).getLenClippedSeq()
		<< "\t" << hp1min
		<< "\t" << hp1max
		<< "\t" << eoffsetmin
		<< "\t" << eoffsetmax
		<< "\t" << maxoverlap
		<< "\t" << hashesoverlap
		<< "\t" << numhashes
		<< "\t" << minoverlaprequired
		<< "\t" << perc << '%'
		<< "\t" << maxcontiguousfreq3counter
		<< "\t" << maxcontiguousfreq32counter
		<< "\t" << totalfreq3counter
		<< "\t" << totalfreq5counter
		<< "\n" << tmp.ol_stronggood << ' ' << tmp.ol_weakgood << ' ' << tmp.ol_belowavgfreq << ' ' << tmp.ol_norept << ' ' << tmp.ol_rept
		<< '\n');
//#define CEBUG(bla)

	  if(SKIM3_chimerahunt.size()){
	    chimeraHuntStoreOverlapCoverage(direction, actreadid, rid2,
					    hp1min,hp1max,hp2min,hp2max);
	  }
	}
      }
    }

    if(sI!=readhashmatches.cend()) countid=sI->rid2;
  }

  if(possiblehits!=0){
    CEBUG("Numhits " << actreadid << "\t" << possiblehits << "\t" << acceptedhits << "\n\n");
  }

//  boost::mutex::scoped_lock lock(SKIM3_globalclassdatamutex);
//  SKIM3_possiblehits+=possiblehits;
//  SKIM3_acceptedhits+=acceptedhits;
}
//#define CEBUG(bla)
//#undef CEBUG_extra_cFPH


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::chimeraHuntStoreOverlapCoverage(const int8 direction, const uint32 actreadid, const uint32 rid2, uint16 hp1min, uint16 hp1max, uint16 hp2min, uint16 hp2max)
{
  bool cebug=false;
  //if(actreadid==0 || rid2==0) cebug=true;

//#define CEBUG(bla)   {if(cebug) cout << bla; cout.flush();}

  CEBUG("checkForChimeras: " << SKIM3_readpool->getRead(actreadid).getName()
	<< " (" << actreadid << ':' << static_cast<int16>(direction) << ") / "
	<< SKIM3_readpool->getRead(rid2).getName()
	<< " (" << rid2 << ":1)\n");
  CEBUG("hp1min: " << hp1min << '\n');
  CEBUG("hp1max: " << hp1max << '\n');
  CEBUG("hp2min: " << hp2min << '\n');
  CEBUG("hp2max: " << hp2max << '\n');

  // *sigh* Must also handle this:
  // instead of re-searching, get them passed by caller
  //
  // E0K6C4E01E3QB2 (-1) / E0K6C4E01C6VEK (1)
  // rid2: 8 eoffset: -49    hp1: 91 hp2: 140
  // rid2: 8 eoffset: -49    hp1: 93 hp2: 142
  // rid2: 8 eoffset: -49    hp1: 95 hp2: 144
  // rid2: 8 eoffset: -49    hp1: 97 hp2: 146
  // ...
  // rid2: 8 eoffset: -48    hp1: 206        hp2: 254
  // rid2: 8 eoffset: -48    hp1: 208        hp2: 256
  // rid2: 8 eoffset: -48    hp1: 210        hp2: 258
  // rid2: 8 eoffset: -46    hp1: 20 hp2: 66
  // rid2: 8 eoffset: -46    hp1: 22 hp2: 68
  // rid2: 8 eoffset: -46    hp1: 24 hp2: 70

  // ok, cut of 2 from each side

  std::vector<uint8> & id1hunt=SKIM3_chimerahunt[actreadid];
  std::vector<uint8> & id2hunt=SKIM3_chimerahunt[rid2];

  //id1hunt.clear();
  //id1hunt.resize(SKIM3_readpool->getRead(actreadid).getLenClippedSeq(),0);
  //id2hunt.clear();
  //id2hunt.resize(SKIM3_readpool->getRead(rid2).getLenClippedSeq(),0);

  if(hp2max-hp2min>4+2*SKIM3_hashsavestepping
     && hp1max-hp1min>4+2*SKIM3_hashsavestepping){

    // the +2 is "magic" ... at least it mixes MIRA not recognising
    //  some chimeras
    hp1min+=SKIM3_hashsavestepping+2;
    hp1max-=SKIM3_hashsavestepping+2;
    hp2min+=SKIM3_hashsavestepping+2;
    hp2max-=SKIM3_hashsavestepping+2;

    uint8 * ptr = &(id2hunt[hp2min]);
    for(uint16 i=0; i<hp2max-hp2min; i++, ptr++){
      *ptr=1;
    }

    if(direction>0){
      ptr = &(id1hunt[hp1min]);
      for(uint16 i=0; i<hp1max-hp1min; i++, ptr++){
	*ptr=1;
      }
    }else{
      ptr = &(id1hunt[SKIM3_readpool->getRead(actreadid).getLenClippedSeq()-1-hp1min]);
      for(uint32 i=0; i<hp1max-hp1min; i++, --ptr){
	*ptr=1;
      }
    }
  }

  if(cebug){
    cout << "id1hunt: " << SKIM3_readpool->getRead(actreadid).getName() << endl;
    uint32 conscounter=0;
    uint32 longest=0;
    for(uint32 i=0; i<id1hunt.size(); i++){
      cout << i << '\t' << static_cast<uint16>(id1hunt[i]) << '\n';
      if(id1hunt[i]){
	conscounter++;
	longest=std::max(longest,conscounter);
      }else{
	conscounter=0;
      }
    }
    cout << "id1longest: " << longest << endl;
    cout << "id2hunt: " << SKIM3_readpool->getRead(rid2).getName() << endl;
    conscounter=0;
    longest=0;
    for(uint32 i=0; i<id2hunt.size(); i++){
      cout << i << '\t' << static_cast<uint16>(id2hunt[i]) << '\n';
      if(id2hunt[i]){
	conscounter++;
	longest=std::max(longest,conscounter);
      }else{
	conscounter=0;
      }
    }
    cout << "id2longest: " << longest << endl;
  }


}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {if(actreadid==0) cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::chimeraHuntLocateChimeras()
{
  for(uint32 actreadid=0; actreadid<SKIM3_chimerahunt.size(); actreadid++){
    std::vector<uint8> & hunt = SKIM3_chimerahunt[actreadid];
    if(hunt.size()==0) continue;

    //if(actreadid==0){
    //  cout << "\nchimhunt\n";
    //  Read::setCoutType(Read::AS_TEXT);
    //  cout << SKIM3_readpool->getRead(actreadid);
    //  cout << "idhunt: " << SKIM3_readpool->getRead(actreadid).getName() << endl;
    //  for(uint32 i=0; i<hunt.size(); i++){
    //	cout << i << '\t' << static_cast<uint16>(hunt[i]) << '\n';
    //  }
    //}

    // we have this vector (with x=any of 0/1:
    // 000....xxx...000
    // now, to make thing easier in search, fill up leading
    // and traing 0s with 1
    // actually, if no fill up were used, this could also be used as
    //  automatic clipping ... but the analyseHashStats() in
    //  conjunction with the dataprocessing.C routine does it already
    //  quite well.
    //
    // BaCh: 09.02.2010
    // Actually, there are instances were the above fail, so let's try a fallback
    // If no chimera found, give back proposed right and left cuts as negative
    // values if these are longer than 4 bases.

    int32 proposedleft=0;
    for(int32 i=0; i<hunt.size() && hunt[i]==0; i++){
      hunt[i]=1;
      proposedleft--;
    };
    int32 proposedright=0;
    for(int32 i=hunt.size()-1; i>=0 && hunt[i]==0; --i){
      hunt[i]=1;
      proposedright--;
    };

    // not only look at holes in the alignments, but also make sure they correlate to
    //  rare or unique kmers by looking at the hash statistics
    auto & bposhashstats=SKIM3_readpool->getRead(actreadid).getBPosHashStats();
    uint32 bphsi=SKIM3_readpool->getRead(actreadid).getLeftClipoff();

    // search for holes with 0s >= 2*basesperhash-2
    int32 consecutivezeroes=0;
    int32 leftcut=0;
    int32 longestleftcut=0;
    int32 longestrightcut=0;

    int32 acti=0;
    bool foundcuts=false;
    for(; acti<hunt.size(); ++acti, ++bphsi){
      if(consecutivezeroes){
	if(hunt[acti]==0 && bposhashstats[bphsi].fwd.getFrequency()<2){
	  consecutivezeroes++;
	}else{
	  if(consecutivezeroes>=2*SKIM3_basesperhash-2){
	    foundcuts=true;
	    CEBUG("Chimera candidate "
		  << SKIM3_readpool->getRead(actreadid).getName()
		  << '\n');
	    CEBUG("consecutivezeroes: " << consecutivezeroes << '\n');
	    CEBUG("leftcut: " << leftcut << '\n');
	    CEBUG("acti: " << acti << '\n');
	    CEBUG("longestleftcut: " << longestleftcut << '\n');
	    CEBUG("longestrightcut: " << longestrightcut << '\n');
	    leftcut=acti;
	    consecutivezeroes=0;
	  }
	}
      }else{
	if(hunt[acti]==0 && bposhashstats[bphsi].fwd.getFrequency()<2){
	  consecutivezeroes++;
	}else{
	  if(acti-leftcut>longestrightcut-longestleftcut){
	    longestleftcut=leftcut;
	    longestrightcut=acti;
	  }
	}
      }
    }

    if(foundcuts){
      CEBUG("Chosen chimeric fragment " << SKIM3_readpool->getRead(actreadid).getName() << " : " << longestleftcut << '\t' << longestrightcut << '\n');
      (*SKIM3_chuntleftcut)[actreadid]=longestleftcut;
      (*SKIM3_chuntrightcut)[actreadid]=longestrightcut;
    }else{
      if(proposedright<-4) (*SKIM3_chuntrightcut)[actreadid]=proposedright;
      if(proposedleft<-4)(*SKIM3_chuntleftcut)[actreadid]=proposedleft;
    }
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * return false and empty tagmaskvector if nothing was set
 * return true and non-empty tagmaskvector if some mask was set
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
bool Skim<TVHASH_T>::fillTagMaskVector(const uint32 readid, std::vector<uint8> & tagmaskvector)
{
  tagmaskvector.clear();

  bool retvalue=false;

  if(SKIM3_hasMNRr[readid]
     || SKIM3_hasFpas[readid]){

    for(uint32 tn=0; tn<SKIM3_readpool->getRead(readid).getNumOfTags(); ++tn){
      if(SKIM3_readpool->getRead(readid).getTag(tn).identifier==Read::REA_tagentry_idMNRr
	 || SKIM3_readpool->getRead(readid).getTag(tn).identifier==Read::REA_tagentry_idSOFApolyA_sequence){
	CEBUG("MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAASK!\n");

	int32 from=SKIM3_readpool->getRead(readid).getTag(tn).from;
	from-=SKIM3_readpool->getRead(readid).getLeftClipoff();
	int32 to=SKIM3_readpool->getRead(readid).getTag(tn).to;
	to-=SKIM3_readpool->getRead(readid).getLeftClipoff();

	/* the masking routines will mask only if *ALL* positions
	   in a hash are masked. This allows for cases like this
	   (shown here: dot == masked)

	   ............A...............

	   to still make hashes that have the "A" in them

	   But for Fpas (poly A), we want stop dead at the beginning
	   of the poly A. Therefore, we need to expand the masked
	   area by (number of bases in a hash)-1
	*/
	if(SKIM3_readpool->getRead(readid).getTag(tn).identifier==Read::REA_tagentry_idSOFApolyA_sequence){
	  from-=(SKIM3_basesperhash-1);
	  to+=(SKIM3_basesperhash-1);
	}

	CEBUG("ftmv: " << from << " " << to << '\n');
	// cache lenclippedseq as read takes some time for this operation
	int32 lenclippedseq=static_cast<int32>(SKIM3_readpool->getRead(readid).getLenClippedSeq());
	for(int32 i=from; i<=to; ++i){
	  if(i>=0 && i<lenclippedseq) {
	    retvalue=true;
	    if(tagmaskvector.empty()){
	      tagmaskvector.resize(lenclippedseq,0);
	    }
	    tagmaskvector[i]=1;
	  }
	}
      }
    }
    if(!retvalue) tagmaskvector.clear();
  }

  //CEBUG("TMV " << SKIM3_readpool->getRead(readid).getName() << '\n');
  //for(uint32 i=0; i<tagmaskvector.size(); ++i){
  //  CEBUG(i << '\t' << static_cast<uint16>(tagmaskvector[i]) << '\n');
  //}

  return retvalue;
}


//#define CEBUG(bla)



/*************************************************************************
 *
 * Around each RMB base, set the base statitsics to "repeat"
 * More a test.
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::correctReadBaseStatisticsByRMB(ReadPool & rp, const uint32 basesperhash)
{
  ProgressIndicator<int32> P(0, rp.size());

  for(uint32 actreadid=0; actreadid<rp.size(); actreadid++){
    P.progress(actreadid);

    Read & actread= rp.getRead(actreadid);
    if(!actread.hasValidData()
       || !actread.isUsedInAssembly()
       || !actread.hasBaseHashStats()) continue;

    if(actread.hasTag(Read::REA_tagentry_idSRMr)
       || actread.hasTag(Read::REA_tagentry_idCRMr)){

      CEBUG("Read " << actread.getName() << endl);
      for(uint32 tagnr=0; tagnr<actread.getNumOfTags(); tagnr++){
	const multitag_t & acttag=actread.getTag(tagnr);
	if(acttag.identifier==Read::REA_tagentry_idSRMr
	   || acttag.identifier==Read::REA_tagentry_idCRMr){

	  CEBUG("Tag " << acttag << endl);
	  std::vector<Read::bposhashstat_t> & bposhashstats=const_cast<std::vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
	  if(!bposhashstats.empty()){
	    Read::setCoutType(Read::AS_TEXT);
	    CEBUG("Before\n" << actread);

	    uint32 bposfrom=acttag.from;
	    uint32 bposto=acttag.to;
	    if(bposto<bposfrom) std::swap(bposfrom,bposto);
	    if(bposfrom>=basesperhash-1){
	      bposfrom-=basesperhash-1;
	    }else{
	      bposfrom=0;
	    }

	    CEBUG("bposfrom: " << bposfrom << "\tbposto: " << bposto << endl);
	    for(; bposfrom<=bposto; bposfrom++){
	      if(bposhashstats[bposfrom].fwd.isValid()
		 && bposhashstats[bposfrom].fwd.getFrequency()<5){
		bposhashstats[bposfrom].fwd.setFrequency(5);
	      }
	      if(bposhashstats[bposfrom].rev.isValid()
		 && bposhashstats[bposfrom].rev.getFrequency()<5){
		bposhashstats[bposfrom].rev.setFrequency(5);
	      }
	    }

	    bposfrom=acttag.from;
	    bposto=acttag.to;
	    if(bposto<bposfrom) std::swap(bposfrom,bposto);
	    bposto+=basesperhash-1;
	    if(bposto>=bposhashstats.size()) bposto=bposhashstats.size()-1;

	    CEBUG("bposfrom: " << bposfrom << "\tbposto: " << bposto << endl);
	    for(; bposfrom<=bposto; bposfrom++){
	      if(bposhashstats[bposfrom].rev.isValid()
		 && bposhashstats[bposfrom].rev.getFrequency()<5){
		bposhashstats[bposfrom].rev.setFrequency(5);
	      }
	    }

	    CEBUG("After\n" << actread);

	  }
	}
      }
    }

  }
  P.finishAtOnce();
  cout << '\n';
}
//#define CEBUG(bla)






/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::skimStreamPrepare(ReadPool & rp, uint32 basesperhash, uint8 hss, const char * additionalregexp)
{
  FUNCSTART("void Skim<TVHASH_T>::skimStreamPrepare(ReadPool & rp, uint32 basesperhash, uint8 hss, const char * additionalregexp)");

  init();

  SKIM3_readpool=&rp;

  if(basesperhash>sizeof(TVHASH_T)*4){
    basesperhash=sizeof(TVHASH_T)*4;
  }
  SKIM3_basesperhash=basesperhash;
  SKIM3_hashsavestepping=hss;

  fillTagStatusInfoOfReads();

  prepareSkim(0, rp.size(), SKIM3_vhraparray,false);

  FUNCEND();
  return;
}
//#define CEBUG(bla)






/*************************************************************************
 *
 * seqtype -1 == all sequencing types, else only the given
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void Skim<TVHASH_T>::findAdaptorRightClip(ReadPool & searchpool, std::vector<int32> & results, int8 seqtype, uint32 minhashes, uint32 numthreads)
{
  FUNCSTART("void Skim<TVHASH_T>::checkForAdaptor(Read & actread, cfh_threaddata_t & cfhd)");

  results.clear();
  results.resize(searchpool.size(),-1);

  SKIM3_farc_searchpool=&searchpool;
  SKIM3_farc_results=&results;
  SKIM3_farc_minhashes=minhashes;
  SKIM3_farc_seqtype=seqtype;

  startMultiThreading(1,numthreads,10000,0,searchpool.size(),
		      boost::bind( &Skim<TVHASH_T>::farcThreadsDataInit, this, _1 ),
		      boost::bind( &Skim<TVHASH_T>::farcThreadLoop, this, _1 ));

  FUNCEND();
}

template<typename TVHASH_T>
int32 Skim<TVHASH_T>::findAdaptorRightClip(Read & actread, uint32 minhashes, readid_t & ridofadapfound, int32 threadid)
{
  FUNCSTART("int32 Skim<TVHASH_T>::findAdaptorRightClip(Read & actread, uint32 minhashes. int32 threadid)");

  int32 ret=0;
  if(threadid<0){
    ret=findAdaptorRightClip_internal(actread,minhashes,ridofadapfound,SKIM3_farcdata_fornonmultithread);
  }else{
    BUGIFTHROW(threadid>=SKIM3_farcd_vector.size(),"Oooooops, trying to use thread id " << threadid << ", but prepared only for " << SKIM3_farcd_vector.size() << " threads???");
    ret=findAdaptorRightClip_internal(actread,minhashes, ridofadapfound, SKIM3_farcd_vector[threadid]);
  }
  return ret;
}

#define CEBUG2(bla)
//#define CEBUG2(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
int32 Skim<TVHASH_T>::findAdaptorRightClip_internal(Read & actread, uint32 minhashes, readid_t & ridofadapfound, farc_threaddata_t & farcd)
{
  FUNCSTART("void Skim<TVHASH_T>::checkForAdaptor(Read & actread, cfh_threaddata_t & cfhd)");

  CEBUG("farc_i: " << actread.getName() << endl);

  if(SKIM3_vashortcuts_begin.empty() || SKIM3_vashortcuts_end.empty()) return -1;
  if(!actread.hasValidData()) return -1;
  uint32 slen=actread.getLenClippedSeq();
  if(slen<SKIM3_basesperhash) return -1;

  // don't really need to clear out this re-used vector
  //   - singlereadvhraparray needs to be big enough to be
  //     written into by transformSeqToVariableHash()
  if(farcd.singlereadvhraparray.size() < slen){
    farcd.singlereadvhraparray.resize(slen);
  }

  farcd.readhashmatches.clear();

  auto srvaI=farcd.singlereadvhraparray.begin();

  std::vector<Read::bposhashstat_t> & bposhashstats=const_cast<std::vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
  uint32 hashesmade;

  uint32 actreadid=0xffffffff;

  {
    int32 bfpos=0;
    int32 bfposinc=1;

    hashesmade=transformSeqToVariableHash(
      actreadid,
      actread,
      actread.getClippedSeqAsChar(),
      slen,
      SKIM3_basesperhash,
      srvaI,
      false,
      1,
      farcd.unused_tagmaskvector,
      bposhashstats,
      bfpos,
      bfposinc
      );
  }
  farcd.singlereadvhraparray.resize(hashesmade);

  CEBUG("hashesmade: " << hashesmade << endl);


  srvaI=farcd.singlereadvhraparray.begin();

  typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::const_iterator lowerbound;
  typename std::vector<typename HashStatistics<TVHASH_T>::vhrap_t>::const_iterator upperbound;
  for(; srvaI != farcd.singlereadvhraparray.end(); ++srvaI){
    lowerbound=SKIM3_vashortcuts_begin[static_cast<uint64>(srvaI->vhash & SKIM3_MAXVHASHMASK)];
    upperbound=SKIM3_vashortcuts_end[static_cast<uint64>(srvaI->vhash & SKIM3_MAXVHASHMASK)];

    // "SKIM3_empty_vector_vhrap_t.end()" is the "nullptr" replacement
    if(SKIM3_completevhraparray_end != lowerbound){
      if(SKIM3_basesperhash>12){
	// with more than 12 bases in a hash, the vhrap array is
	//  subdivided
	auto p=equal_range(lowerbound,
			   upperbound,
			   *srvaI,
			   compareVHRAPArrayElem_);
	lowerbound=p.first;
	upperbound=p.second;
      }

      for(;lowerbound!=upperbound; lowerbound++){
	//CEBUG("/// " << actreadid << '\t' << lowerbound->readid << '\n');
	//CEBUG("/// take!\n");
	farcd.readhashmatches.resize(farcd.readhashmatches.size()+1);
	farcd.readhashmatches.back().rid2=lowerbound->readid;
	farcd.readhashmatches.back().hashpos1=srvaI->hashpos;
	farcd.readhashmatches.back().hashpos2=lowerbound->hashpos;
	farcd.readhashmatches.back().eoffset=srvaI->hashpos - lowerbound->hashpos;
	farcd.readhashmatches.back().bhashstats=srvaI->bhashstats;

	CEBUG2("added: " << farcd.readhashmatches.back());
      }
    }
  }

  int32 retvalue=-1;

  if(!farcd.readhashmatches.empty()){
    checkForPotentialAdaptorHits(1, actreadid, actread, farcd.tmpmatchwith, farcd.readhashmatches);

    if(!farcd.tmpmatchwith.empty()){
      CEBUG2("Hits of: " << actread.getName() << endl);
      int32 leftmostpos=100000000;
      uint32 largesthash=0;
      auto leftmostokI=farcd.tmpmatchwith.cend();
      auto largestokI=leftmostokI;
      for(auto ssmwsI=farcd.tmpmatchwith.cbegin(); ssmwsI != farcd.tmpmatchwith.cend(); ++ssmwsI){
	CEBUG2(actread.getName() << " to " << SKIM3_readpool->getRead(ssmwsI->otherid).getName() << " : " << *ssmwsI);
	// have at least an estimated 50% identity, false positive adaptor identifications normally have well below 50%
	if(ssmwsI->numhashes>=minhashes && ssmwsI->percent_in_overlap>=50){
	  int32 numhashes=ssmwsI->numhashes;
	  if(ssmwsI->eoffset < 0) numhashes+=ssmwsI->eoffset;
	  if(numhashes>=minhashes){
	    if(ssmwsI->eoffset < leftmostpos){
	      CEBUG2("newleft!!!!!!\n");
	      leftmostpos=ssmwsI->eoffset;
	      leftmostokI=ssmwsI;
	    }
	    if(numhashes>largesthash){
	      CEBUG2("newlarge!!!!!!\n");
	      largesthash=numhashes;
	      largestokI=ssmwsI;
	    }
	  }
	}
      }
      if(largestokI!=farcd.tmpmatchwith.end()){
	//cout << "Chosen " << SKIM3_readpool->getRead(largestokI->otherid).getName() << ": " << *largestokI;
	ridofadapfound=largestokI->otherid;
	retvalue=largestokI->eoffset;
	if(retvalue<0) retvalue=0;
	if(leftmostokI!=largestokI){
	  //cout << "Leftmost != Biggest " << SKIM3_readpool->getRead(largestokI->otherid).getName() << ": " << *largestokI;
	}
      }
    }
  }

  return retvalue;
}

//#define CEBUG(bla)





/*************************************************************************
 *
 * a simplified version of checkForPotentialHits()
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG_extra_cFPH

template<typename TVHASH_T>
void Skim<TVHASH_T>::checkForPotentialAdaptorHits(const int8 direction, const uint32 actreadid, Read & actread, std::vector<matchwithsorter_t> & tmpmatchwith, std::vector<readhashmatch_t> & readhashmatches)
{
  //bool dodebug=false;

  CEBUG("Potential hits of " << actread.getName() << " (" << actreadid << ") (" << static_cast<int16>(direction) << '/' << actread.getLenClippedSeq() << ")\n----------------\n");
  CEBUG(actread << endl);
  CEBUG("----------------\n");

  tmpmatchwith.clear();

  // readhashmatches should not be empty ... normally.
  // but new method to deal with megahubs reduces this vector, keeping only
  //  'approximately normal' frequencies. Which in turn means: some vectors
  //  might be completely emptied
  // so, if it is empty, return immediately
  if(readhashmatches.empty()) return;

  // use non-parallel sort as we are already in a multithreading environment here
  mstd::ssort(readhashmatches, sortreadhashmatch_t_);

  auto sI=readhashmatches.cbegin();
  uint32 countid=sI->rid2;
  while(sI != readhashmatches.cend()){

    uint32 rid2=sI->rid2;
    uint16 oldhashpos=sI->hashpos1;
    uint16 hp1min=0xffff;
    uint16 hp1max=0;

    uint16 hp2min=0xffff;
    uint16 hp2max=0;

    int32  eoffsetmin=0x7fffffff;
    int32  eoffsetmax=0x80000000;
    int32 oldeoffset=sI->eoffset;

    uint32 numhashes=0;

    for(;sI != readhashmatches.cend() && sI->rid2 == countid; ++sI){
      CEBUG(*sI);

      // this ensures that the eoffset between two following
      //  entries may not differ by too much (2 bases here for adaptor search)
      // IF they do, then this is treated like a different hit
      //  by breaking the loop
      if(abs(sI->eoffset - oldeoffset) > 2){
	CEBUG("BREAKER!\n");
	break;
      }
      numhashes++;

      hp1min=std::min(hp1min,sI->hashpos1);
      hp1max=std::max(hp1max,sI->hashpos1);
      eoffsetmin=std::min(eoffsetmin,sI->eoffset);
      eoffsetmax=std::max(eoffsetmax,sI->eoffset);
      oldeoffset=sI->eoffset;

      hp2min=std::min(hp2min,sI->hashpos2);
      hp2max=std::max(hp2max,sI->hashpos2);

#ifdef CEBUG_extra_cFPH
      {
	boost::mutex::scoped_lock lock(SKIM3_coutmutex);
	CEBUG(sI->rid2
	      << "\t" << SKIM3_readpool->getRead(sI->rid2).getName()
	      << "\t" << SKIM3_readpool->getRead(sI->rid2).getLenClippedSeq()
	      << "\t" << sI->eoffset
	      << "\t" << sI->hashpos1
	      << "\t" << oldhashpos
	      << '\n');
      }
#endif

      oldhashpos=sI->hashpos1;
    }

    int32 maxoverlap;

    // adjust min positions for the hash length
    hp1min-=(SKIM3_basesperhash-1);
    hp2min-=(SKIM3_basesperhash-1);

    int32 eoffsetmean=eoffsetmin+(eoffsetmax-eoffsetmin)/2;

    // calc max overlap
    // currently only for one offset
    if(eoffsetmean<0){
      maxoverlap=std::min(SKIM3_readpool->getRead(rid2).getLenClippedSeq()+eoffsetmean,actread.getLenClippedSeq());
    }else{
      maxoverlap=std::min(actread.getLenClippedSeq()-eoffsetmean,SKIM3_readpool->getRead(rid2).getLenClippedSeq());
    }

    // correct the maxoverlap by the modulo of the hash steps as the
    //  border hashes will be found only in 1/(hash stepping) cases
    maxoverlap=maxoverlap-(maxoverlap%SKIM3_hashsavestepping);

    // hashe3soverlap is not the number of hashes in the overlap,
    // but the length of the overlap
    int32 hashesoverlap=hp1max-hp1min+1;

    int32 perc=100*hashesoverlap/maxoverlap;


    int32 minpercentrequired=0;

    // look a bit closer at potential perfect matches
    if(perc == 100){
      if(eoffsetmin != eoffsetmax){
	// this could not be: have at least two different expected offsets
	//  and a 100% coverage. Side effects from intra-read repeats
	// therefore, make sure this does not get through as a 100% match
	perc=99;
      }else if((numhashes-1)*SKIM3_hashsavestepping+SKIM3_basesperhash < maxoverlap){
	// maxoverlap covers the whole potential overlap, but
	//  there are not enough hashes supporting for 100% match
	//  (base mismatch somewhere)
	// reduce the percentage to show it's not a perfect match
	perc=99;
      }
    }else if(eoffsetmin == eoffsetmax){
      if(perc>100) {
	perc=100;
      }else{
	uint32 maxnumhashes=((maxoverlap-1-SKIM3_basesperhash)/SKIM3_hashsavestepping)+1;
	if(perc>=minpercentrequired && numhashes==maxnumhashes){
	  CEBUG("maxnumhashes 100% saver: "  << perc << '\n');
	  perc=100;
	}
      }
    }

    // 13 was a bit too small, still produces rare false positives even if percent_in_overlap >= 60
    // only workable countermeasure: 16 (and that is also more in-line with minimm overlaps one
    //  should use)
    int32 minoverlaprequired=16;

#ifdef CEBUG_extra_cFPH
    {
      boost::mutex::scoped_lock lock(SKIM3_coutmutex);
      CEBUG("eomin: " << eoffsetmin << "\teomax: " << eoffsetmax
	    << "\tmor: " << minoverlaprequired
	    << "\tho: " << hashesoverlap
	    << "\t%: " << perc
	    << "\t%<: " << minpercentrequired << endl);
    }
#endif

    // we take the hit if the overlap percentage is above threshold
    if(hashesoverlap >= minoverlaprequired
       && perc>=minpercentrequired){

//#define CEBUG(bla)   {if(actreadid==273252 && rid2==273250) cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla; cout.flush();}
      matchwithsorter_t tmp;
      tmp.otherid=rid2;
      tmp.eoffset=eoffsetmean;

      if(perc>100) perc=100;
      tmp.percent_in_overlap=perc;
      tmp.numhashes=numhashes;
      tmp.estimscore=0;
      tmp.taken=false;

      tmpmatchwith.push_back(tmp);

      CEBUG("Pushing possible hit with offset: " << tmp.eoffset << endl
	    << rid2
	    << "\t" << actreadid
	    << "\t" << SKIM3_readpool->getRead(rid2).getLenClippedSeq()
	    << "\t" << hp1min
	    << "\t" << hp1max
	    << "\t" << eoffsetmin
	    << "\t" << eoffsetmax
	    << "\t" << maxoverlap
	    << "\t" << hashesoverlap
	    << "\t" << numhashes
	    << "\t" << minoverlaprequired
	    << "\t" << perc << '%'
	    << '\n');
//#define CEBUG(bla)

    }
    if(sI!=readhashmatches.cend()) countid=sI->rid2;
  }

}
//#define CEBUG(bla)



template<typename TVHASH_T>
void Skim<TVHASH_T>::prepareForMultithreadFarc(uint32 numthreads)
{
  FUNCSTART("void Skim<TVHASH_T>::prepareForMultithreadFarc()");

  SKIM3_numthreads=numthreads;
  farcThreadsDataInit(SKIM3_numthreads);
}

template<typename TVHASH_T>
void Skim<TVHASH_T>::farcThreadsDataInit(const uint32 numthreads)
{
  FUNCSTART("void Skim<TVHASH_T>::cfhThreadsDataInit(const uint32 numthreads)");

  SKIM3_farcd_vector.resize(numthreads);
  for(uint32 ti=0; ti<numthreads;++ti){

    SKIM3_farcd_vector[ti].readhashmatches.clear();
    SKIM3_farcd_vector[ti].readhashmatches.reserve(2000);
    SKIM3_farcd_vector[ti].singlereadvhraparray.clear();
    SKIM3_farcd_vector[ti].singlereadvhraparray.reserve(2000);
    SKIM3_farcd_vector[ti].unused_tagmaskvector.clear();
    //SKIM3_farcd_vector[ti].unused_tagmaskvector.reserve(2000);
    SKIM3_farcd_vector[ti].tmpmatchwith.clear();
    SKIM3_farcd_vector[ti].tmpmatchwith.reserve(2000);

  }
  FUNCEND();
}

template<typename TVHASH_T>
void Skim<TVHASH_T>::farcThreadLoop(const uint32 threadnr)
{
  FUNCSTART("void Skim<TVHASH_T>::threadloop(const uint32 threadnr)");

  // threads need their own try() catch() block

  try {
    CEBUG("Thread: " << threadnr << " starting.\n");

    BUGIFTHROW(threadnr>=SKIM3_farcd_vector.size(),"threadnr>=SKIM3_farcd_vector.size()???");
    farc_threaddata_t & farcd=SKIM3_farcd_vector[threadnr];

    farcd.readhashmatches.clear();
    farcd.singlereadvhraparray.clear();
    farcd.unused_tagmaskvector.clear();
    farcd.tmpmatchwith.clear();

    readid_t dummy=0; // in this version, we do not give back the read id of the adaptor found, but need a variable to call the internal routine

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
	  if(SKIM3_farc_seqtype < 0
	     || SKIM3_farc_searchpool->getRead(readi).getSequencingType() == SKIM3_farc_seqtype){
	    int32 clip=findAdaptorRightClip_internal(SKIM3_farc_searchpool->getRead(readi),SKIM3_farc_minhashes,dummy, farcd);
	    if(clip>=0){
	      boost::mutex::scoped_lock lock(SKIM3_resultfileoutmutex);
	      (*SKIM3_farc_results)[readi]=clip;
	    }
	  }
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


// explicit template instantiations needed for the linker to create these for library files
template class Skim<vhash64_t>;
#ifndef KMER_INTERNALTYPE
template class Skim<vhash128_t>;
template class Skim<vhash256_t>;
template class Skim<vhash512_t>;
#endif
