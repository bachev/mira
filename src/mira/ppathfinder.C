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

#include "mira/ppathfinder.H"
#include "errorhandling/errorhandling.H"

#include "util/stlimprove.H"

// for time measurements
#include <sys/times.h>


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)

// cs1 for normal clocking ('user compatible' as is does not disturb user)
// cs2 is for extensive clocking output, more for analysis of MIRA behaviour
//  during development

#define CLOCK_STEPS1
#ifndef PUBLICQUIET
#define CLOCK_STEPS2
#endif

#define CLOCK_STEPS2


/*************************************************************************
 *
 * This should be a compile time check, but don't know how!
 *
 *************************************************************************/

bool PPathfinder::PPF_staticinit=PPathfinder::staticInit();
bool PPathfinder::staticInit()
{
  //cout << "QTG_END " << QTG_END << " QTE_END " << QTE_END << endl;
  if(QTG_END<QTE_END){
    cout << "PPathfinder problem: QTG_END " << QTG_END << " is < QTE_END " << QTE_END
	 << "\nThis is a problem which can only be fixed by changing the source code! " << endl;
    exit(101);
  }
  return true;
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

// Plain vanilla constructor
PPathfinder::PPathfinder(std::vector<MIRAParameters> * params, ReadPool * readpool, necontainer_t * overlap_edges, adsfcontainer_t * adsfacts, std::vector<Align> * aligncache, std::vector<int8> * used_ids, std::vector<uint8> * multicopies, std::vector<uint8> * hasmcoverlaps, std::vector<uint8> * hasreptoverlaps, std::vector<uint8> * hasnoreptoverlaps, std::vector<uint8> * istroublemaker, std::vector<uint8> * incorrectibleorchimera, std::vector<uint8> * maybespoilsport, std::vector<uint8> * wellconnected,std::vector<necontainer_t::iterator > * lowerbound_oedges_ptr, std::vector<Contig::templateguessinfo_t> * astemplateguesses)
{
  FUNCSTART("PPathfinder::PPathfinder()");

  PPF_actcontig_ptr=nullptr;

  PPF_miraparams_ptr=params;
  PPF_readpool_ptr=readpool;
  PPF_overlap_edges_ptr=overlap_edges;
  PPF_adsfacts_ptr=adsfacts;
  PPF_aligncache_ptr=aligncache;
  PPF_used_ids_ptr=used_ids;
  PPF_multicopies_ptr=multicopies;
  PPF_hasmcoverlap_ptr=hasmcoverlaps;
  PPF_hasreptoverlap_ptr=hasreptoverlaps;
  PPF_hasnoreptoverlap_ptr=hasnoreptoverlaps;
  PPF_istroublemaker_ptr=istroublemaker;
  PPF_incorrectibleorchimera_ptr=incorrectibleorchimera;
  PPF_maybespoilsport_ptr=maybespoilsport;
  PPF_wellconnected_ptr=wellconnected;
  PPF_lowerbound_oedges_ptr=lowerbound_oedges_ptr;
  PPF_astemplateguesses_ptr=astemplateguesses;

  //////
  if(PPF_hasnoreptoverlap_ptr->empty()){
    priv_ppFillNoRept();
  }
  if(PPF_hasreptoverlap_ptr->empty()){
    priv_ppFillRept();
  }
  if(PPF_maybespoilsport_ptr->empty()){
    priv_ppFillSpoilSport();
  }
  //////

  PPF_pafparams_ptr=&((*PPF_miraparams_ptr)[0].getPathfinderParams());

  PPF_ids_in_contig_list.reserve(readpool->size());
  PPF_ids_added_oltype.resize(readpool->size(),0);
  PPF_blacklisted_ids.resize(readpool->size(),0);
  PPF_tmparray.resize(readpool->size(),0);

  // give the small store iterators to banned overlaps a capacity of
  //  500k entries.
  // should be enough even for extremely deep RNASeq data
  PPF_overlapsbanned_smallstore.reserve(500000);

  PPF_wantscleanoverlapends=0;
  PPF_mintotalnonmatches=0;
  PPF_allowedseqtype=ReadGroupLib::SEQTYPE_END;

  FUNCEND();
}


PPathfinder::~PPathfinder()
{
  FUNCSTART("PPathfinder::~PPathfinder()");

  // do nothing???

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/
void PPathfinder::priv_ppFillNoRept()
{
  PPF_hasnoreptoverlap_ptr->clear();
  PPF_hasnoreptoverlap_ptr->resize(PPF_readpool_ptr->size(),false);
  for(auto & oee : *PPF_overlap_edges_ptr){
    if(oee.ol_norept){
      (*PPF_hasnoreptoverlap_ptr)[oee.rid1]=true;
      (*PPF_hasnoreptoverlap_ptr)[oee.linked_with]=true;
    }
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/
void PPathfinder::priv_ppFillRept()
{
  PPF_hasreptoverlap_ptr->clear();
  PPF_hasreptoverlap_ptr->resize(PPF_readpool_ptr->size(),false);
  for(auto & oee : *PPF_overlap_edges_ptr){
    if(oee.ol_rept){
      (*PPF_hasreptoverlap_ptr)[oee.rid1]=true;
      (*PPF_hasreptoverlap_ptr)[oee.linked_with]=true;
    }
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/
void PPathfinder::priv_ppFillSpoilSport()
{
  auto & rp=*PPF_readpool_ptr;
  auto & wc=*PPF_wellconnected_ptr;
  auto & mss=*PPF_maybespoilsport_ptr;

  mss.clear();
  mss.resize(rp.size(),false);

  for(size_t rpi=0; rpi<rp.size(); ++rpi){
    if(!wc[rpi]){
      auto & actread=rp[rpi];
      auto bhsI=actread.getBPosHashStats().cbegin();
      int32 xpos=actread.getLeftClipoff();
      if(xpos>=0){
	advance(bhsI,xpos);
	if(!bhsI->fwd.isValid()){
	  mss[rpi]=true;
	}
      }
      xpos=actread.getRightClipoff();
      if(xpos>0){
	bhsI=actread.getBPosHashStats().begin();
	advance(bhsI,xpos-1);
	if(!bhsI->rev.isValid()){
	  mss[rpi]=true;
	}
      }
    }
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::prepareForNewContig(Contig & con)
{
  FUNCSTART("void PPathfinder::prepareForNewContig(Contig & con)");

  PPF_actcontig_ptr=&con;

  for(auto & qu : PPF_queues){
    BUGIFTHROW(!qu.empty(),"Queue not empty?");
  }

  for(auto rid : PPF_ids_in_contig_list){
    PPF_ids_added_oltype[rid]=ADDED_NOTADDED;
  }
  PPF_ids_in_contig_list.clear();
  PPF_rails_in_contig_list.clear();

  std::vector<int8> & lr_used_ids = *PPF_used_ids_ptr;
  auto & cr=PPF_actcontig_ptr->getContigReads();
  for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
    if(pcrI.getORPID()>=0){
      PPF_ids_in_contig_list.push_back(pcrI.getORPID());
      lr_used_ids[pcrI.getORPID()]=1;
      if(pcrI->isRail() || pcrI->isBackbone()){
	PPF_ids_added_oltype[pcrI.getORPID()]=ADDED_BY_BACKBONE;
	if(pcrI->isRail()) PPF_rails_in_contig_list.push_back(pcrI.getORPID());
      }else{
	priv_insertRIDIntoDenovoQueues(pcrI.getORPID());
      }
    }
  }

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * almost like prepareForNewContig(), but here we trust the external
 *  used_ids and synch the contig according to it.
 * also, entries in PPF_ids_added_oltype do not get all cleared, just the ones
 *  of reads which "disappeared" between the time the pathfinder last saw
 *  the contig and now (due to cutdown by assembly etc.)
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::resyncContig()
{
  FUNCSTART("void PPathfinder::resyncContig()");

  BUGIFTHROW(PPF_actcontig_ptr==nullptr,"PPF_actcontig_ptr==nullptr ???");

  for(auto & qu : PPF_queues){
    BUGIFTHROW(!qu.empty(),"Queue not empty?");
  }

  std::vector<int8> & lr_used_ids = *PPF_used_ids_ptr;
  for(auto rid : PPF_ids_in_contig_list){
    if(!lr_used_ids[rid]){
      PPF_ids_added_oltype[rid]=0;
    }
  }
  PPF_ids_in_contig_list.clear();
  PPF_rails_in_contig_list.clear();

  auto & cr=PPF_actcontig_ptr->getContigReads();
  for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
    if(pcrI.getORPID()>=0){
      if(lr_used_ids[pcrI.getORPID()]){
	PPF_ids_in_contig_list.push_back(pcrI.getORPID());
	if(pcrI->isRail() || pcrI->isBackbone()){
	  if(pcrI->isRail()) PPF_rails_in_contig_list.push_back(pcrI.getORPID());
	}
	// nope ... both map() and denovo() should take care of that
	// priv_insertRIDIntoDenovoQueues(pcrI.getORPID());
      }else{
      }
    }
  }

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::denovo()
{
  FUNCSTART("void PPathfinder::denovo()");

  Contig::templateguessinfo_t tguess;

  if(PPF_ids_in_contig_list.size()
     &&(PPF_actcontig_ptr->getContigReads().size()==0
	|| PPF_actcontig_ptr->getContigLength()==0)){
    MIRANOTIFY(Notify::FATAL,"Ummm, forgot to call prepareForNewContig? " << PPF_ids_in_contig_list.size() << " " << PPF_actcontig_ptr->getContigReads().size() << " " << PPF_actcontig_ptr->getContigLength());
  }

  priv_basicSetup();

  readid_t startid=-1;
  bool alreadydone=false;
  size_t fdnmaxdist=0; // fillDenovoQueue maxdist
  if(PPF_actcontig_ptr->getContigLength()==0){
    startid=priv_getNextStartID();
    CEBUG("startid: " << startid << endl);
    if(startid<0) return;
    BUGIFTHROW(startid >= static_cast<readid_t>(PPF_used_ids_ptr->size()), "Starting with read id " << startid << " which is >= number of reads " << PPF_used_ids_ptr->size() << " ?");
    BUGIFTHROW((*PPF_used_ids_ptr)[startid], "Startid " << startid << " already used ???");
    CEBUG("Read: " << PPF_readpool_ptr->getRead(startid).getName() << endl);
    BUGIFTHROW(PPF_readpool_ptr->getRead(startid).getLenClippedSeq()==0,"PPF_readpool_ptr->getRead(" << startid << ").getLenClippedSeq()==0???\n" << PPF_readpool_ptr->getRead(startid).getName());

    // if we keep long repeats separated:
    //   if we start with a non-multicopy read, forbid
    //    multicopy/multicopy overlaps
    if((*PPF_multicopies_ptr)[startid]>0) {
#ifndef PUBLICQUIET
      cout << "\nStarted with multicopy.\n";
#endif
      if(PPF_pafparams_ptr->paf_use_genomic_algorithms) PPF_actcontig_ptr->setLongRepeatStatus(true);
    }

    ++PPF_readaddattempts;
    PPF_actcontig_ptr->addRead(*PPF_aligncache_ptr,
			       nullptr, startid, startid, 1,
			       (*PPF_multicopies_ptr)[startid],
			       0,
			       tguess,
			       PPF_contigerrstat);
    (*PPF_used_ids_ptr)[startid]=1;
    PPF_ids_in_contig_list.push_back(startid);
    PPF_ids_added_oltype[startid]=1;
    priv_showProgress();

    if(PPF_bsccontent==BSCC_SINGLETS){
      // singlets
      alreadydone=true;
      CEBUG("That's a singlet, PPF_bsccontent is " << static_cast<uint16>(PPF_bsccontent) << "\n");
    }
  }else{
    // we seem to continue construction of a contig ... prepare for that

    // first: contig may have lost reads, account for that and update some internal arrays
    // use PPF_tmparray as temporary for info "which reads are in the contig"
    for(auto index : PPF_tmparray_idxused) PPF_tmparray[index]=0;
    auto & cr=PPF_actcontig_ptr->getContigReads();
    for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
      if(pcrI.getORPID()>=0){
	PPF_tmparray[pcrI.getORPID()]=1;
	PPF_tmparray_idxused.push_back(pcrI.getORPID());
      }
    }
    // now adjust PPF_ids_added_oltype && PPF_ids_in_contig_list
    //  the blacklist was already taken care of by basicSetup
    // rewrite PPF_ids_in_contig_list in place, swapping unused ids to the end
    //  and then resizing the array afterwards
    BUGIFTHROW(PPF_ids_in_contig_list.empty(),"PPF_ids_in_contig_list.empty() ???");
    auto iicI=PPF_ids_in_contig_list.begin();
    auto lastusedI=PPF_ids_in_contig_list.end();
    while(iicI!=lastusedI){
      if(PPF_tmparray[*iicI]){
	++iicI;
      }else{
	BUGIFTHROW((*PPF_used_ids_ptr)[*iicI],"hmmmm, AS_used not reflecting reality? status for " << *iicI << " not 0!");
	PPF_ids_added_oltype[*iicI]=0;
	--lastusedI;
	if(iicI!=lastusedI) std::swap(*iicI,*lastusedI);
      }
    }
    PPF_ids_in_contig_list.resize(lastusedI-PPF_ids_in_contig_list.begin());
    // cleanup tmparray
    for(auto index : PPF_tmparray_idxused) PPF_tmparray[index]=0;
    PPF_tmparray_idxused.clear();

    // TODO:
    //  set this to llength of longest read or so?
    //  or else 1000
    fdnmaxdist=1000;

    BUGIFTHROW(PPF_ids_in_contig_list.size()!=PPF_actcontig_ptr->getContigReads().size(),"PPF_ids_in_contig_list.size() " << PPF_ids_in_contig_list.size() << " != PPF_actcontig_ptr->getContigReads().size()" << PPF_actcontig_ptr->getContigReads().size());
  }

  if(!alreadydone){
    priv_fillDenovoQueues(fdnmaxdist);
    priv_loopDenovo();
  }

  FUNCEND();
  return;
};
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_showProgress()
{
  const uint32 cpl=60;
  if(PPF_buildcontig_newlinecounter==0){
    cout << '[' << PPF_ids_in_contig_list.size() << "]\t";
    PPF_buildcontig_newlinecounter=cpl;
    PPF_timing_pathsearch=0;
    PPF_timing_connadd=0;
  }
#ifdef PUBLICQUIET
  PPF_contigerrstat.dumpStatus();
#else
  PPF_contigerrstat.dumpStatus(true);
#endif
  --PPF_buildcontig_newlinecounter;
  if(PPF_buildcontig_newlinecounter==0){
    cout << "   " << PPF_actcontig_ptr->getContigLength();
#ifdef CLOCK_STEPS1
    cout << "\tpft\t" << PPF_timing_pathsearch/cpl;
    cout << " / " << PPF_timing_connadd/cpl;
#endif
    cout << endl;
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_initialiseLowerBoundOEdges()
{
  // initialise the quick lowerbound_oedges lookup vector if not already done
  if(PPF_lowerbound_oedges_ptr->empty()){
#ifdef CLOCK_STEPS2
    timeval tv;
    gettimeofday(&tv,nullptr);
#endif
    auto & lowerbound_oedges=*PPF_lowerbound_oedges_ptr;
    lowerbound_oedges.resize(PPF_readpool_ptr->size(),PPF_overlap_edges_ptr->end());
    newedges_t tmp;
    for(uint32 li=0; li<lowerbound_oedges.size(); ++li) {
      tmp.rid1=li;
      lowerbound_oedges[li]=mstd::lower_bound(*PPF_overlap_edges_ptr,
					      tmp,
					      newedges_t::sortComparatorByRIDUp);
    }
#ifdef CLOCK_STEPS2
    cout << "Timing priv_initialiseLowerBoundOEdges: " << diffsuseconds(tv) << endl;
#endif
  }
  return;
}


/*************************************************************************
 *
 * called as basic init routine by denovo() or map()
 * take care not to overwrite/destroy structures/values which might be
 *  needed by iterative calls of assembly to map() or denovo()
 *  Those would rather be made by prepareForNewContig()
 *
 *************************************************************************/

void PPathfinder::priv_basicSetup()
{
  PPF_buildcontig_newlinecounter=0;
  PPF_timing_pathsearch=0;
  PPF_timing_connadd=0;
  PPF_readaddattempts=0;

  priv_initialiseLowerBoundOEdges();

  if(!PPF_overlapsbanned_smallstore.empty()){
    if(PPF_overlapsbanned_smallstore.size() < PPF_overlapsbanned_smallstore.capacity()){
#ifndef PUBLICQUIET
      cout << "Clear pf_banned quick: " << PPF_overlapsbanned_smallstore.size() << '\n';
#endif
      for(auto obI : PPF_overlapsbanned_smallstore){
	obI->pf_banned=false;
      }
    }else{
#ifndef PUBLICQUIET
      cout << "Clear pf_banned full\n";
#endif
      for(auto & ne : *PPF_overlap_edges_ptr){
	ne.pf_banned=false;
      }
    }
    PPF_overlapsbanned_smallstore.clear();
  }

  while(!PPF_blacklist_queues.empty()){
    for(auto rid : PPF_blacklist_queues.front()){
      PPF_blacklisted_ids[rid]=0;
    }
    PPF_blacklist_queues.pop();
  }

  return;
}


/*************************************************************************
 *
 * pilfered from Pathfinder::n4_searchBestStrongGoodStartEnvironment_sub
 *
 *************************************************************************/

void PPathfinder::priv_fdns_genome()
{
  FUNCSTART("void PPathfinder::priv_fdns_genome()");

#ifdef CLOCK_STEPS2
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  PPF_bsccontent=BSCC_GENOME_BESTQUAL;
  priv_fdns_g_subFillCache(true,true,true,true,true);

  if(PPF_beststartcache.empty()){
    PPF_bsccontent=BSCC_GENOME_MEDQUAL_TROUBLEMAKER_NOTSTRONG_NOTMULTICOPY_WELLCONNECTED_NOKMERFORK;
    priv_fdns_g_subFillCache(false,false,true,true,true);
    if(PPF_beststartcache.empty()){
      PPF_bsccontent=BSCC_GENOME_MEDQUAL_TROUBLEMAKER_NOTSTRONG_NOTMULTICOPY_WELLCONNECTED_NOKMERFORK;
      priv_fdns_g_subFillCache(false,false,true,true,true);
      if(PPF_beststartcache.empty()){
	PPF_bsccontent=BSCC_GENOME_MEDQUAL_TROUBLEMAKER_NOTSTRONG_MULTICOPY_WELLCONNECTED_NOKMERFORK;
	priv_fdns_g_subFillCache(false,false,false,true,true);
	if(PPF_beststartcache.empty()){
	  PPF_bsccontent=BSCC_GENOME_BADQUAL_NOTWELLCONNECTED;
	  priv_fdns_g_subFillCache(false,false,false,false,false);
	}
      }
    }
  }

  if(!PPF_beststartcache.empty()){
    if(PPF_beststartcache.size()<=65535){
      mstd::ssort(PPF_beststartcache, PPathfinder::beststartinfo_t::ltclustersize);
    }else{
      mstd::psort(PPF_beststartcache, PPathfinder::beststartinfo_t::ltclustersize);
    }
  }else{
    // found nothing ... now put singlets into list
    priv_fdns_fillSinglets();
  }

  if(!PPF_beststartcache.empty()){
    CEBUG("Startcache take " << PPF_beststartcache.back());
  }

#ifdef CLOCK_STEPS2
  cout << "Timing priv_fdns_genome total: " << diffsuseconds(tv)
       << "\n";
#endif


  FUNCEND();

  return;

  //  return 0;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_fdns_fillSinglets()
{
  FUNCSTART("void PPathfinder::priv_fdns_fillSinglets()");

#ifdef CLOCK_STEPS2
  timeval tvtmp;
  gettimeofday(&tvtmp,nullptr);
#endif
  PPF_bsccontent=BSCC_SINGLETS;
  {
    beststartinfo_t tmp;
    tmp.bsi_clustersize=0;
    tmp.bsi_numconnects=0;

    readid_t uid=0;
    for(auto tmpused : *PPF_used_ids_ptr){
      if(!tmpused){
	if(PPF_readpool_ptr->getRead(uid).isUsedInAssembly()
	   && PPF_readpool_ptr->getRead(uid).getLenClippedSeq() > 0) {
	  tmp.bsi_rid=uid;
	  PPF_beststartcache.push_back(tmp);
	}
      }
      ++uid;
    }
  }
  CEBUG("Startcache filled with singlets: " << PPF_beststartcache.size() << '\n');
#ifdef CLOCK_STEPS2
  cout << "Timing priv_fdns_fillSinglets: " << diffsuseconds(tvtmp) << endl;
#endif
}

/*************************************************************************
 *
 *  fills the cache for the best start sites
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::priv_fdns_g_subFillCache(bool wanttroublemakercheck, bool wantstronggoodcheck, bool wantmulticopycheck, bool wantwellconnectedcheck, bool wantnokmerfork)
{
  FUNCSTART("void PPathfinder::priv_fdns_g_subFillCache(bool wanttroublemakercheck, bool wantstronggoodcheck, bool wantmulticopycheck, bool wantwellconnectedcheck, bool wantnokmerfork)");

#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  const auto & PPF_used_ids = *PPF_used_ids_ptr;
  const auto & PPF_istroublemaker = *PPF_istroublemaker_ptr;
  const auto & PPF_multicopies = *PPF_multicopies_ptr;
  const auto & PPF_wellconnected = *PPF_wellconnected_ptr;

  std::vector<bool> uid_in_cluster(PPF_used_ids.size(),false);

  // in EST / RNASeq projects, things can get *very* ugly ... clusters of
  //  several million reads (rRNA etc). This takes ages, even on fast processors
  //  (E.g., 18 minutes for 10m 100bp where 40 to 50% may be rRNA).
  // Solution: restrict growth of the unlooked vector to (currently) 8000
  //  elements
  // Effect: in genome assemblies, may lead to several start sites per
  //  area. In EST / RNASeq assemblies, leads to several start sites per
  //  heavy multicopy area.
  //  Oh, and reduces time *a lot*, especially in conjunction with -PF:mscft
  //  (see end of function)

  std::vector<readid_t> unlooked;
  unlooked.reserve(8000);
  bool unlookedaccepts=true;
  bool unlookthreshreached=false;

  PPF_beststartcache.reserve(10000);

  bool checkcluster=true;
  // checking cluster ID does not make much sense when we're in the really
  //  bad part of an assembly
  // For assemblies with millions of reads, this leads to the start cache
  //  constantly running dry at the end of an assembly, needing constant
  //  refill and triggering small clsuter checks. That slows things down
  //  quite a bit.
  // switch it off in that case
  if(!wanttroublemakercheck
     && !wantstronggoodcheck
     && !wantmulticopycheck
     && !wantwellconnectedcheck) checkcluster=false;

  long int maxmusec=PPF_pafparams_ptr->paf_max_startcache_filltime*1000000;
  timeval starttv;
  timeval curtv;
  gettimeofday(&starttv,nullptr);

  for(readid_t actid=0; actid < static_cast<readid_t>(PPF_used_ids.size()); ++actid) {
    if(!PPF_readpool_ptr->getRead(actid).isUsedInAssembly()
       || PPF_readpool_ptr->getRead(actid).getLenClippedSeq() == 0) continue;

    CEBUG("aid: " << actid);
    CEBUG("\tused " << (int16) PPF_used_ids[actid]);
    CEBUG("\ttm " <<(uint16) PPF_istroublemaker[actid]);
    CEBUG("\n");
    if(PPF_used_ids[actid]!=0
       || (wanttroublemakercheck && PPF_istroublemaker[actid]!=0)
       || (wantmulticopycheck && PPF_multicopies[actid]!=0)
       || (wantwellconnectedcheck && !PPF_wellconnected[actid])
       || (wantnokmerfork && PPF_readpool_ptr->getRead(actid).hasKMerFork())
       || (checkcluster && uid_in_cluster[actid])) continue;

    uint32 goodclustersize=0;
    uint32 goodclusterbestid=0;
    uint32 maxconnects=0;
    unlooked.clear();
    unlooked.push_back(actid);

    unlookedaccepts=true;

    readid_t lookid;
    while(!unlooked.empty()){
      CEBUG("unlooked.size(): " << unlooked.size() << endl);
      lookid=unlooked.back();
      unlooked.pop_back();

      uid_in_cluster[lookid]=true;

      uint32 numconnects=0;
      auto oeI=(*PPF_lowerbound_oedges_ptr)[lookid];
      for(;oeI!=PPF_overlap_edges_ptr->end() && lookid==oeI->rid1;++oeI){
	CEBUG("  lid: " << oeI->linked_with);
	CEBUG("\tused " << (int16) PPF_used_ids[oeI->linked_with]);
	CEBUG("\ttm " <<(uint16) PPF_istroublemaker[oeI->linked_with]);
	CEBUG("\n");
	CEBUG("  " << *oeI);

	if((wantstronggoodcheck && !oeI->ol_stronggood)
	   || PPF_used_ids[oeI->linked_with]!=0
	   || (wanttroublemakercheck && PPF_istroublemaker[oeI->linked_with]!=0)
	   || (wantmulticopycheck && PPF_multicopies[oeI->linked_with]!=0)
	   || (wantwellconnectedcheck && !PPF_wellconnected[oeI->linked_with])
	   || uid_in_cluster[oeI->linked_with]) continue;
	++numconnects;

	if(unlookedaccepts) {
	  unlooked.push_back(oeI->linked_with);
	  // Check whether to stop looking for more
	  // Do not do that only when multicopies are also allowed as it may be that
	  //  people give projects *just* with heavy multicopy clusters and then
	  //  it would not work (as only a minority may then be multicopy with
	  //  regard to the already high copy number)
	  if(unlooked.size()==unlooked.capacity()){
	    unlookedaccepts=false;
	    unlookthreshreached=true;
	  }
	}
      }
      if(numconnects) ++goodclustersize;
      if(numconnects>maxconnects){
	goodclusterbestid=lookid;
	maxconnects=numconnects;
	CEBUG("goodclusterbestid: " << goodclusterbestid << endl);
	CEBUG("maxconnects: " << numconnects << endl);
      }
    }

    if(maxconnects){
      beststartinfo_t tmp;
      tmp.bsi_rid=goodclusterbestid;
      tmp.bsi_clustersize=goodclustersize;
      tmp.bsi_numconnects=maxconnects;
      PPF_beststartcache.push_back(tmp);
    }

    // see whether we need to care about time (-PF:mscft)
    if(maxmusec>=0 && actid%64==0 && !PPF_beststartcache.empty()) {
      if(diffsuseconds(starttv) >= maxmusec) {
	cout << "Non-deterministic behaviour of assembly likely: " << diffsuseconds(starttv) << " -PF:mscft threshold hit " << maxmusec << ".\n";
	break;
      }
    }
  }

#ifdef CLOCK_STEPS2
  cout << "Timing priv_fdns_subFillCache " << wanttroublemakercheck << " " << wantstronggoodcheck << " " << wantmulticopycheck << " " << wantwellconnectedcheck << " : " << diffsuseconds(tvtotal)
       << "\nStartcache size: " << PPF_beststartcache.size() << endl;
#endif

  if(unlookthreshreached){
    cout << "hit unlooked threshold\n";
  }

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_fdns_est()
{
  FUNCSTART("void PPathfinder::priv_fdns_est()");

  if(PPF_haflevel_max.empty()) priv_fillHAFLevelInfo();

#ifdef CLOCK_STEPS2
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  PPF_bsccontent=BSCC_EST_WELLCONNECTED_REPT6P;
  priv_fdns_e_subFillCache(true,6);
  if(PPF_beststartcache.empty()){
    PPF_bsccontent=BSCC_EST_WELLCONNECTED_REPT;
    priv_fdns_e_subFillCache(true,5);
    if(PPF_beststartcache.empty()){
      PPF_bsccontent=BSCC_EST_WELLCONNECTED_H4;
      priv_fdns_e_subFillCache(true,4);
      if(PPF_beststartcache.empty()){
	PPF_bsccontent=BSCC_EST_WELLCONNECTED_H3;
	priv_fdns_e_subFillCache(true,3);
	if(PPF_beststartcache.empty()){
	  PPF_bsccontent=BSCC_EST_WELLCONNECTED_H2;
	  priv_fdns_e_subFillCache(true,2);
	  if(PPF_beststartcache.empty()){
	    PPF_bsccontent=BSCC_EST_NOTWELLCONNECTED;
	    priv_fdns_e_subFillCache(false,0);
	  }
	}
      }
    }
  }

  if(!PPF_beststartcache.empty()){
    if(PPF_beststartcache.size()<=65535){
      mstd::ssort(PPF_beststartcache, PPathfinder::beststartinfo_t::ltclustersize);
    }else{
      mstd::psort(PPF_beststartcache, PPathfinder::beststartinfo_t::ltclustersize);
    }
  }else{
    // found nothing ... now put singlets into list
    priv_fdns_fillSinglets();
  }

  if(!PPF_beststartcache.empty()){
    CEBUG("Startcache take " << PPF_beststartcache.back());
  }

#ifdef CLOCK_STEPS2
  cout << "Timing priv_fdns_est total: " << diffsuseconds(tv)
       << "\n";
#endif
}

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::priv_fdns_e_subFillCache(bool wantwellconnectedcheck, uint8 minallowedfreq)
{
  FUNCSTART("void PPathfinder::priv_fdns_e_subFillCache(bool wantwellconnectedcheck, uint8 minallowedfreq)");

#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  const auto & PPF_used_ids = *PPF_used_ids_ptr;
  const auto & PPF_wellconnected = *PPF_wellconnected_ptr;

  std::vector<bool> uid_in_cluster(PPF_used_ids.size(),false);

  // in EST / RNASeq projects, things can get *very* ugly ... clusters of
  //  several million reads (rRNA etc). This takes ages, even on fast processors
  //  (E.g., 18 minutes for 10m 100bp where 40 to 50% may be rRNA).
  // Solution: restrict growth of the unlooked vector to (currently) 8000
  //  elements
  // Effect: in genome assemblies, may lead to several start sites per
  //  area. In EST / RNASeq assemblies, leads to several start sites per
  //  heavy multicopy area.
  //  Oh, and reduces time *a lot*, especially in conjunction with -PF:mscft
  //  (see end of function)

  std::vector<readid_t> unlooked;
  unlooked.reserve(8000);
  bool unlookedaccepts=true;
  bool unlookthreshreached=false;

  PPF_beststartcache.reserve(10000);

  long int maxmusec=PPF_pafparams_ptr->paf_max_startcache_filltime*1000000;
  timeval starttv;
  timeval curtv;
  gettimeofday(&starttv,nullptr);

  for(readid_t actid=0; actid < static_cast<readid_t>(PPF_used_ids.size()); ++actid) {
    CEBUG("aid: " << actid);
    CEBUG("\tused " << (int16) PPF_used_ids[actid]);
    CEBUG("\ttm " <<(uint16) (*PPF_istroublemaker_ptr)[actid]);
    CEBUG("\n");
    if(PPF_used_ids[actid]!=0
       || (PPF_haflevel_min[actid] < minallowedfreq)
       || (wantwellconnectedcheck && !PPF_wellconnected[actid])
       || uid_in_cluster[actid]) continue;

    uint32 goodclustersize=0;
    uint32 goodclusterbestid=0;
    uint32 maxconnects=0;
    unlooked.clear();
    unlooked.push_back(actid);

    unlookedaccepts=true;

    readid_t lookid;
    while(!unlooked.empty()){
      CEBUG("unlooked.size(): " << unlooked.size() << endl);
      lookid=unlooked.back();
      unlooked.pop_back();

      uid_in_cluster[lookid]=true;

      uint32 numconnects=0;
      auto oeI=(*PPF_lowerbound_oedges_ptr)[lookid];
      for(;oeI!=PPF_overlap_edges_ptr->end() && lookid==oeI->rid1;++oeI){
	CEBUG("  lid: " << oeI->linked_with);
	CEBUG("\tused " << (int16) PPF_used_ids[oeI->linked_with]);
	CEBUG("\ttm " <<(uint16) (*PPF_istroublemaker_ptr)[oeI->linked_with]);
	CEBUG("\n");
	CEBUG("  " << *oeI);

	if(PPF_used_ids[oeI->linked_with]!=0
	   || (PPF_haflevel_min[oeI->linked_with] < minallowedfreq)
	   || (wantwellconnectedcheck && !PPF_wellconnected[oeI->linked_with])
	   || uid_in_cluster[oeI->linked_with]) continue;
	++numconnects;

	if(unlookedaccepts) {
	  unlooked.push_back(oeI->linked_with);
	  // Check whether to stop looking for more
	  // Do not do that only when multicopies are also allowed as it may be that
	  //  people give projects *just* with heavy multicopy clusters and then
	  //  it would not work (as only a minority may then be multicopy with
	  //  regard to the already high copy number)
	  if(unlooked.size()==unlooked.capacity()){
	    unlookedaccepts=false;
	    unlookthreshreached=true;
	  }
	}
      }
      if(numconnects) ++goodclustersize;
      if(numconnects>maxconnects){
	goodclusterbestid=lookid;
	maxconnects=numconnects;
	CEBUG("goodclusterbestid: " << goodclusterbestid << endl);
	CEBUG("maxconnects: " << numconnects << endl);
      }
    }

    if(maxconnects){
      beststartinfo_t tmp;
      tmp.bsi_rid=goodclusterbestid;
      tmp.bsi_clustersize=goodclustersize;
      tmp.bsi_numconnects=maxconnects;
      PPF_beststartcache.push_back(tmp);
    }

    // see whether we need to care about time (-PF:mscft)
    if(maxmusec>=0 && actid%64==0 && !PPF_beststartcache.empty()) {
      if(diffsuseconds(starttv) >= maxmusec) {
	cout << "Non-deterministic behaviour of assembly likely: -PF:mscft threshold hit.\n";
	break;
      }
    }
  }

#ifdef CLOCK_STEPS2
  cout << "Timing priv_fdns_e_subFillCache " << wantwellconnectedcheck << " " << static_cast<uint16>(minallowedfreq) << ": " << diffsuseconds(tvtotal)
       << "\nStartcache size: " << PPF_beststartcache.size() << endl;
#endif

  if(unlookthreshreached){
    cout << "hit unlooked threshold\n";
  }

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_fillHAFLevelInfo()
{
  FUNCSTART("void PPathfinder::priv_fillHAFLevelInfo()");

  PPF_haflevel_max.clear();
  PPF_haflevel_max.resize(PPF_readpool_ptr->size(),0);
  PPF_haflevel_maxseen=0;
  PPF_haflevel_min.clear();
  PPF_haflevel_min.resize(PPF_readpool_ptr->size(),255);
  PPF_haflevel_minseen=255;
  for(readid_t ri=0; ri<PPF_readpool_ptr->size(); ++ri){
    for(auto & bphs : PPF_readpool_ptr->getRead(ri).getBPosHashStats()){
      if(bphs.fwd.isValid()){
	if(bphs.fwd.getFrequency() > PPF_haflevel_max[ri]) PPF_haflevel_max[ri]=bphs.fwd.getFrequency();
	if(bphs.fwd.getFrequency() > PPF_haflevel_maxseen) PPF_haflevel_maxseen=bphs.fwd.getFrequency();
	if(bphs.fwd.getFrequency() < PPF_haflevel_min[ri]) PPF_haflevel_min[ri]=bphs.fwd.getFrequency();
	if(bphs.fwd.getFrequency() < PPF_haflevel_minseen) PPF_haflevel_minseen=bphs.fwd.getFrequency();
      }
    }
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

readid_t PPathfinder::priv_getNextStartID()
{
#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  readid_t bestid=-1;
  PPF_bsrandry=false;
  if(!PPF_beststartcache.empty()){
    bestid=priv_gnsi_helper();
  }
  if(bestid<0 && PPF_beststartcache.empty()){
    priv_fillDenovoStartCache();
    bestid=priv_gnsi_helper();
  }
#ifdef CLOCK_STEPS2
  PPF_timing_pathsearch+=diffsuseconds(tv);
#endif
  return bestid;
}

readid_t PPathfinder::priv_gnsi_helper()
{
  const std::vector<int8> & lr_used_ids = *PPF_used_ids_ptr;

  // if there's something in the cache and
  //  the read it points to hasn't been included via other means
  //  (like jumping over a small repeat), return that cached entry
  while(!PPF_beststartcache.empty()){
    CEBUG("Startcache look " << PPF_beststartcache.back());
    auto bestid=PPF_beststartcache.back().bsi_rid;
    if(!lr_used_ids[bestid]){
      CEBUG("Startcache taken\n");
      return bestid;
    }
    PPF_beststartcache.pop_back();
  }

  return -1;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_fillDenovoQueues(size_t maxdist)
{
  auto & cr=PPF_actcontig_ptr->getContigReads();
  for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
    if(pcrI.getORPID()>=0){
      if(maxdist == 0
	 || PPF_actcontig_ptr->getContigLength() < maxdist
	 || pcrI.getReadStartOffset() <= maxdist
	 || pcrI.getReadStartOffset()+pcrI->getLenClippedSeq() >= PPF_actcontig_ptr->getContigLength()-maxdist){
	priv_insertRIDIntoDenovoQueues(pcrI.getORPID());
      }
    }
  }

  //cout << "DN Queues\n";
  //for(uint32 dqi=0; dqi<PPF_queues.size(); ++dqi){
  //  cout << dqi << "\t" << PPF_queues[dqi].size() << endl;
  //}
}


/*************************************************************************
 *
 * returns queue number where read was inserted (or QTG_END if not)
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
uint32 PPathfinder::priv_iridnq_genome(readid_t insertrid)
{
  FUNCSTART("void PPathfinder::priv_iridnq_genome(readid_t insertrid)");

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  BUGIFTHROW(insertrid<0, "insertrid " << insertrid << " < 0 ???");
  BUGIFTHROW(insertrid>=static_cast<readid_t>(PPF_readpool_ptr->size()), "insertrid " << insertrid << " >= rp size() " << PPF_readpool_ptr->size() << " ???");
  BUGIFTHROW(!PPF_ids_added_oltype[insertrid], "insertrid " << insertrid << " (" << PPF_readpool_ptr->getRead(insertrid).getName() << ") not in contig ???");

  CEBUG("Trying to insert " << PPF_readpool_ptr->getRead(insertrid).getName() << " (" << insertrid << ")...");

  // lr_ == local reference
  const std::vector<int8> &  lr_used_ids = *PPF_used_ids_ptr;
  const std::vector<uint8> &  lr_wellconnected = *PPF_wellconnected_ptr;
  //const std::vector<uint8> & lr_istroublemaker = *PPF_istroublemaker_ptr;

  auto oeI=(*PPF_lowerbound_oedges_ptr)[insertrid];
  auto bestoeI=oeI;
  uint32 bestoelevel=QTG_END;

  for(;oeI!=PPF_overlap_edges_ptr->end() && insertrid==oeI->rid1;++oeI){
    CEBUG("\nbest oe level " << bestoelevel);

    // don't bother looking at if overlap is banned
    CEBUG("\ncheck " << PPF_readpool_ptr->getRead(oeI->linked_with).getName() << " chkban...");
    if(oeI->pf_banned) continue;

    // don't bother looking at if read linked to is already used or temporarily blacklisted
    CEBUG(" chkuse of " << PPF_readpool_ptr->getRead(oeI->linked_with).getName() << " (" << oeI->linked_with << ")...");
    if(lr_used_ids[oeI->linked_with]) continue;

    // don't bother looking at if read linked to is temporarily blacklisted
    CEBUG(" chkblcklst ...");
    if(PPF_blacklisted_ids[oeI->linked_with]) continue;

    // of course, rails and backbones are not suited as new reads
    CEBUG(" chkrailbb ...");
    if(PPF_readpool_ptr->getRead(oeI->linked_with).isRail()
       || PPF_readpool_ptr->getRead(oeI->linked_with).isBackbone()) continue;

    CEBUG(" may take");

    readid_t linkedwithid=oeI->linked_with;
    readid_t linkedwith_partnerid=PPF_readpool_ptr->getRead(linkedwithid).getTemplatePartnerID();

    bool swb=oeI->ol_stronggood | oeI->ol_weakgood | oeI->ol_belowavgfreq;  // swb: overlap is Strong / Weak / Belowavg

    bool havedecision=false;

    if((*PPF_incorrectibleorchimera_ptr)[linkedwithid]){
      if(bestoelevel>QTG_INCORRECTIBLEORCHIMERA){
	bestoelevel=QTG_INCORRECTIBLEORCHIMERA;
	bestoeI=oeI;
      }
      havedecision=true;
    }else if((*PPF_maybespoilsport_ptr)[linkedwithid]){
      if(bestoelevel>QTG_MAYBESPOILSPORT){
	bestoelevel=QTG_MAYBESPOILSPORT;
	bestoeI=oeI;
      }
      havedecision=true;
    }

    if(!havedecision && swb && oeI->ol_norept){
      // 0-14
      if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOREPT && lr_wellconnected[linkedwithid]){
	// 0-2
	havedecision=true;
	if(oeI->ol_stronggood){
	  if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTSTRONG_WELLCONNECTED){
	    bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTSTRONG_WELLCONNECTED;
	    bestoeI=oeI;
	    break; // early termination of for loop, we won't find something better.
	  }
	}else if(oeI->ol_weakgood){
	  if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTWEAK_WELLCONNECTED){
	    bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTWEAK_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}else{
	  if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTBELOWAVG_WELLCONNECTED){
	    bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTBELOWAVG_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOREPT){
	// 3-5
	havedecision=true;
	if(oeI->ol_stronggood){
	  if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTSTRONG){
	    bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTSTRONG;
	    bestoeI=oeI;
	  }
	}else if(oeI->ol_weakgood){
	  if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTWEAK){
	    bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTWEAK;
	    bestoeI=oeI;
	  }
	}else{
	  if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTBELOWAVG){
	    bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTBELOWAVG;
	    bestoeI=oeI;
	  }
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOTREPT && lr_wellconnected[linkedwithid]){
	// 6-8
	havedecision=true;
	if(oeI->ol_stronggood){
	  if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTSTRONG_WELLCONNECTED){
	    bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTSTRONG_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}else if(oeI->ol_weakgood){
	  if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTWEAK_WELLCONNECTED){
	    bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTWEAK_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}else{
	  if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTBELOWAVG_WELLCONNECTED){
	    bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTBELOWAVG_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOTREPT){
	// 9-11
	havedecision=true;
	if(oeI->ol_stronggood){
	  if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTSTRONG){
	    bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTSTRONG;
	    bestoeI=oeI;
	  }
	}else if(oeI->ol_weakgood){
	  if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTWEAK){
	    bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTWEAK;
	    bestoeI=oeI;
	  }
	}else{
	  if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTBELOWAVG){
	    bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTBELOWAVG;
	    bestoeI=oeI;
	  }
	}
      }else if(swb && lr_wellconnected[linkedwithid]){
	// 12-14
	havedecision=true;
	if(oeI->ol_stronggood){
	  if(bestoelevel>QTG_OLNOREPTSTRONG_WELLCONNECTED){
	    bestoelevel=QTG_OLNOREPTSTRONG_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}else if(oeI->ol_weakgood){
	  if(bestoelevel>QTG_OLNOREPTWEAK_WELLCONNECTED){
	    bestoelevel=QTG_OLNOREPTWEAK_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}else{
	  if(bestoelevel>QTG_OLNOREPTBELOWAVG_WELLCONNECTED){
	    bestoelevel=QTG_OLNOREPTBELOWAVG_WELLCONNECTED;
	    bestoeI=oeI;
	  }
	}
      }
    }

    if(!havedecision
       && oeI->ol_stronggood
       && linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOREPT && lr_wellconnected[linkedwithid]){
      // 14a
      havedecision=true;
      if(bestoelevel>QTG_TPARTNERNOREPT_STRONG_WELLCONNECTED){
	bestoelevel=QTG_TPARTNERNOREPT_STRONG_WELLCONNECTED;
	bestoeI=oeI;
      }
    }

    if(!havedecision
       && swb && oeI->ol_norept){
      // 15-17
      havedecision=true;
      if(oeI->ol_stronggood){
	if(bestoelevel>QTG_OLNOREPTSTRONG){
	  bestoelevel=QTG_OLNOREPTSTRONG;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_weakgood){
	if(bestoelevel>QTG_OLNOREPTWEAK){
	  bestoelevel=QTG_OLNOREPTWEAK;
	  bestoeI=oeI;
	}
      }else{
	if(bestoelevel>QTG_OLNOREPTBELOWAVG){
	  bestoelevel=QTG_OLNOREPTBELOWAVG;
	  bestoeI=oeI;
	}
      }
    }

    if(!havedecision){
      // rest, 18-27

      if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOREPT && oeI->ol_norept){
	if(bestoelevel>QTG_TPARTNERNOREPT_OLNOREPTOTHER){
	  bestoelevel=QTG_TPARTNERNOREPT_OLNOREPTOTHER;
	  bestoeI=oeI;
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOREPT && !oeI->ol_norept && !oeI->ol_rept){
	if(bestoelevel>QTG_TPARTNERNOREPT_OLNOTREPTOTHER){
	  bestoelevel=QTG_TPARTNERNOREPT_OLNOTREPTOTHER;
	  bestoeI=oeI;
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOTREPT && oeI->ol_norept){
	if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOREPTOTHER){
	  bestoelevel=QTG_TPARTNERNOTREPT_OLNOREPTOTHER;
	  bestoeI=oeI;
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOTREPT && !oeI->ol_norept && !oeI->ol_norept && !oeI->ol_rept){
	if(bestoelevel>QTG_TPARTNERNOTREPT_OLNOTREPTOTHER){
	  bestoelevel=QTG_TPARTNERNOTREPT_OLNOTREPTOTHER;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_norept){
	// 22
	if(bestoelevel>QTG_OLNOREPTOTHER){
	  bestoelevel=QTG_OLNOREPTOTHER;
	  bestoeI=oeI;
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOREPT && oeI->ol_rept){
	if(bestoelevel>QTG_TPARTNERNOREPT_OLREPT){
	  bestoelevel=QTG_TPARTNERNOREPT_OLREPT;
	  bestoeI=oeI;
	}
      }else if(linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid]==ADDED_BY_NOTREPT && oeI->ol_rept){
	if(bestoelevel>QTG_TPARTNERNOTREPT_OLREPT){
	  bestoelevel=QTG_TPARTNERNOTREPT_OLREPT;
	  bestoeI=oeI;
	}
      }else if(!oeI->ol_norept & !oeI->ol_rept){
	if(bestoelevel>QTG_OLNOTREPT){
	  bestoelevel=QTG_OLNOTREPT;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept){
	if(bestoelevel>QTG_OLREPTSTRONG_WELLCONNECTED && lr_wellconnected[linkedwithid] && oeI->ol_stronggood){
	  bestoelevel=QTG_OLREPTSTRONG_WELLCONNECTED;
	  bestoeI=oeI;
	}else if(bestoelevel>QTG_OLREPTWEAK_WELLCONNECTED && lr_wellconnected[linkedwithid] && oeI->ol_weakgood){
	  bestoelevel=QTG_OLREPTWEAK_WELLCONNECTED;
	  bestoeI=oeI;
	}else if(bestoelevel>QTG_OLREPTSTRONG && oeI->ol_stronggood){
	  bestoelevel=QTG_OLREPTSTRONG;
	  bestoeI=oeI;
	}else if(bestoelevel>QTG_OLREPTWEAK && oeI->ol_weakgood){
	  bestoelevel=QTG_OLREPTWEAK;
	  bestoeI=oeI;
	}else if(bestoelevel>QTG_OLREPT){
	  bestoelevel=QTG_OLREPT;
	  bestoeI=oeI;
	}
      }else{
	if(bestoelevel>QTG_OLOTHER){
	  bestoelevel=QTG_OLOTHER;
	  bestoeI=oeI;
	}
      }
    }
  }


  if(bestoelevel!=QTG_END){
    PPF_queues[bestoelevel].push(ppfweightelem_t(bestoeI->best_weight,bestoeI));
    CEBUG("\nInserted " << PPF_readpool_ptr->getRead(bestoeI->rid1).getName() << " in queue " << bestoelevel << " with oe " << *bestoeI << endl);
  }else{
    CEBUG("\nNo insertion\n");
  }

#ifdef CLOCK_STEPS1
  PPF_timing_pathsearch+=diffsuseconds(tv);
#endif

  FUNCEND();
  return bestoelevel;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * returns queue number where read was inserted (or QTG_END if not)
 *  (really, qt_G_ and not qt_E_ !)
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
//#define CEBUG(bla)   {if(insertrid==10002) {cout << bla; cout.flush();} }
uint32 PPathfinder::priv_iridnq_est(readid_t insertrid)
{
  FUNCSTART("uint32 PPathfinder::priv_iridnq_est(readid_t insertrid)");

  CEBUG("priv_iridnq_est" << endl);

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  // TODO: for de-novo, this gets called anyway earlier
  //  but mapping not. See whether to put somewhere else
  //  for mapping
  if(PPF_haflevel_max.empty()) priv_fillHAFLevelInfo();

  BUGIFTHROW(insertrid<0, "insertrid " << insertrid << " < 0 ???");
  BUGIFTHROW(insertrid>=static_cast<readid_t>(PPF_readpool_ptr->size()), "insertrid " << insertrid << " >= rp size() " << PPF_readpool_ptr->size() << " ???");
  BUGIFTHROW(!PPF_ids_added_oltype[insertrid], "insertrid " << insertrid << " (" << PPF_readpool_ptr->getRead(insertrid).getName() << ") not in contig ???");

  CEBUG("Trying to insert " << PPF_readpool_ptr->getRead(insertrid).getName() << " (" << insertrid << ")...");

  // lr_ == local reference
  const std::vector<int8> &  lr_used_ids = *PPF_used_ids_ptr;
  const std::vector<uint8> &  lr_wellconnected = *PPF_wellconnected_ptr;
  //const std::vector<uint8> & lr_istroublemaker = *PPF_istroublemaker_ptr;

  bool has_tpartner;
  bool has_tpartnerwc;
  bool has_refwc;
  bool has_newwc;

  auto oeI=(*PPF_lowerbound_oedges_ptr)[insertrid];
  auto bestoeI=oeI;
  uint32 bestoelevel=QTG_END;
  for(;oeI!=PPF_overlap_edges_ptr->end() && insertrid==oeI->rid1;++oeI){
    // don't bother looking at if overlap is banned
    CEBUG("\ncheck " << PPF_readpool_ptr->getRead(oeI->linked_with).getName() << " chkban...");
    if(oeI->pf_banned) continue;

    // don't bother looking at if read linked to is already used or temporarily blacklisted
    CEBUG(" chkuse of " << PPF_readpool_ptr->getRead(oeI->linked_with).getName() << " (" << oeI->linked_with << ")...");
    if(lr_used_ids[oeI->linked_with]) continue;

    // don't bother looking at if read linked to is temporarily blacklisted
    CEBUG(" chkblcklst ...");
    if(PPF_blacklisted_ids[oeI->linked_with]) continue;

    // of course, rails and backbones are not suited as new reads
    CEBUG(" chkrailbb ...");
    if(PPF_readpool_ptr->getRead(oeI->linked_with).isRail()
       || PPF_readpool_ptr->getRead(oeI->linked_with).isBackbone()) continue;

    CEBUG(" may take");

    bool havedecision=false;
    readid_t linkedwithid=oeI->linked_with;

    if((*PPF_incorrectibleorchimera_ptr)[linkedwithid]){
      havedecision=true;
      if(bestoelevel>QTE_INCORRECTIBLEORCHIMERA){
	bestoelevel=QTE_INCORRECTIBLEORCHIMERA;
	bestoeI=oeI;
      }
    }

    if(!havedecision){
      readid_t linkedwith_partnerid=PPF_readpool_ptr->getRead(linkedwithid).getTemplatePartnerID();

      has_tpartner=linkedwith_partnerid>=0 && PPF_ids_added_oltype[linkedwith_partnerid];
      has_tpartnerwc=has_tpartner && lr_wellconnected[linkedwith_partnerid];
      has_refwc=lr_wellconnected[insertrid];
      has_newwc=lr_wellconnected[linkedwithid];

      if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	 && PPF_haflevel_min[linkedwith_partnerid]>=6
	 && PPF_haflevel_min[insertrid]>=6
	 && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_TPARTNERWCREPT6P_WCREPT6P_OLREPT_REPT6PWC){
	  bestoelevel=QTE_TPARTNERWCREPT6P_WCREPT6P_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	  break; // early termination of for loop, we won't find something better.
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[linkedwith_partnerid]>=5
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_TPARTNERWCREPT5_WCREPT5_OLREPT_REPT6PWC){
	  bestoelevel=QTE_TPARTNERWCREPT5_WCREPT5_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[linkedwith_partnerid]>=5
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=5
	){
	if(bestoelevel>QTE_TPARTNERWCREPT5_WCREPT5_OLREPT_REPT5WC){
	  bestoelevel=QTE_TPARTNERWCREPT5_WCREPT5_OLREPT_REPT5WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[insertrid]>=6
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_TPARTNERWC_WCREPT6P_OLREPT_REPT6PWC){
	  bestoelevel=QTE_TPARTNERWC_WCREPT6P_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_TPARTNERWC_WCREPT5_OLREPT_REPT6PWC){
	  bestoelevel=QTE_TPARTNERWC_WCREPT5_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=5
	){
	if(bestoelevel>QTE_TPARTNERWC_WCREPT5_OLREPT_REPT5WC){
	  bestoelevel=QTE_TPARTNERWC_WCREPT5_OLREPT_REPT5WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_TPARTNERWC_WC_OLREPT_REPT6PWC){
	  bestoelevel=QTE_TPARTNERWC_WC_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc
	       && PPF_haflevel_min[linkedwithid]>=5
	){
	if(bestoelevel>QTE_TPARTNERWC_WC_OLREPT_REPT5WC){
	  bestoelevel=QTE_TPARTNERWC_WC_OLREPT_REPT5WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartnerwc && has_refwc && has_newwc){
	if(bestoelevel>QTE_TPARTNERWC_WC_OLREPT_WC){
	  bestoelevel=QTE_TPARTNERWC_WC_OLREPT_WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_refwc && has_newwc
	       && PPF_haflevel_min[insertrid]>=6
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_WCREPT6P_OLREPT_REPT6PWC){
	  bestoelevel=QTE_WCREPT6P_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_refwc && has_newwc
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_WCREPT5_OLREPT_REPT6PWC){
	  bestoelevel=QTE_WCREPT5_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_refwc && has_newwc
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=5
	){
	if(bestoelevel>QTE_WCREPT5_OLREPT_REPT5WC){
	  bestoelevel=QTE_WCREPT5_OLREPT_REPT5WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_newwc
	       && PPF_haflevel_min[insertrid]>=6
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_REPT6P_OLREPT_REPT6PWC){
	  bestoelevel=QTE_REPT6P_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_newwc
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=6
	){
	if(bestoelevel>QTE_REPT5_OLREPT_REPT6PWC){
	  bestoelevel=QTE_REPT5_OLREPT_REPT6PWC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_newwc
	       && PPF_haflevel_min[insertrid]>=5
	       && PPF_haflevel_min[linkedwithid]>=5
	){
	if(bestoelevel>QTE_REPT5_OLREPT_REPT5WC){
	  bestoelevel=QTE_REPT5_OLREPT_REPT5WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartner && has_refwc && has_newwc){
	if(bestoelevel>QTE_TPARTNER_WC_OLREPT_WC){
	  bestoelevel=QTE_TPARTNER_WC_OLREPT_WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_tpartner && has_newwc){
	if(bestoelevel>QTE_TPARTNER_OLREPT_WC){
	  bestoelevel=QTE_TPARTNER_OLREPT_WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_refwc && has_newwc){
	if(bestoelevel>QTE_WC_OLREPT_WC){
	  bestoelevel=QTE_WC_OLREPT_WC;
	  bestoeI=oeI;
	}
      }else if(oeI->ol_rept && has_newwc){
	if(bestoelevel>QTE_OLREPT_WC){
	  bestoelevel=QTE_OLREPT_WC;
	  bestoeI=oeI;
	}
      }else if(has_tpartnerwc && has_refwc && has_newwc){
	if(bestoelevel>QTE_TPARTNERWC_WC_OL_WC){
	  bestoelevel=QTE_TPARTNERWC_WC_OL_WC;
	  bestoeI=oeI;
	}
      }else if(has_tpartner && has_refwc && has_newwc){
	if(bestoelevel>QTE_TPARTNER_WC_OL_WC){
	  bestoelevel=QTE_TPARTNER_WC_OL_WC;
	  bestoeI=oeI;
	}
      }else if(has_tpartner && has_newwc){
	if(bestoelevel>QTE_TPARTNER_OL_WC){
	  bestoelevel=QTE_TPARTNER_OL_WC;
	  bestoeI=oeI;
	}
//    }else if(has_tpartner){
//      if(bestoelevel>QTE_TPARTNER){
//	bestoelevel=QTE_TPARTNER;
//	bestoeI=oeI;
//      }
      }else if(has_refwc && has_newwc){
	if(bestoelevel>QTE_WC_OL_WC){
	  bestoelevel=QTE_WC_OL_WC;
	  bestoeI=oeI;
	}
      }else if(has_newwc){
	if(bestoelevel>QTE_OL_WC){
	  bestoelevel=QTE_OL_WC;
	  bestoeI=oeI;
	}
      }else if(PPF_haflevel_min[linkedwithid]>=6){
	if(bestoelevel>QTE_REPT6P){
	  bestoelevel=QTE_REPT6P;
	  bestoeI=oeI;
	}
      }else if(PPF_haflevel_min[linkedwithid]>=5){
	if(bestoelevel>QTE_REPT5){
	  bestoelevel=QTE_REPT5;
	  bestoeI=oeI;
	}
      }else if(PPF_haflevel_min[linkedwithid]>=4){
	if(bestoelevel>QTE_H4){
	  bestoelevel=QTE_H4;
	  bestoeI=oeI;
	}
      }else if(PPF_haflevel_min[linkedwithid]>=3){
	if(bestoelevel>QTE_H3){
	  bestoelevel=QTE_H3;
	  bestoeI=oeI;
	}
      }else{
	if(bestoelevel>QTE_OTHER){
	  bestoelevel=QTE_OTHER;
	  bestoeI=oeI;
	}
      }
    }
  }


  if(bestoelevel!=QTG_END){
    PPF_queues[bestoelevel].push(ppfweightelem_t(bestoeI->best_weight,bestoeI));
    CEBUG("\nInserted " << PPF_readpool_ptr->getRead(bestoeI->rid1).getName() << " in queue " << bestoelevel << endl);
  }else{
    CEBUG("\nNo insertion\n");
  }

#ifdef CLOCK_STEPS1
  PPF_timing_pathsearch+=diffsuseconds(tv);
#endif

  FUNCEND();
  return bestoelevel;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * return
 *   as value: in which queue it found something
 *   by caller: iterator, == PPF_overlap_edges_ptr->end() if not found
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
uint32 PPathfinder::priv_getNextOverlapFromDenovoQueue(necontainer_t::iterator & oeI)
{
  FUNCSTART("void PPathfinder::priv_getNextOverlapFromDenovoQueue()");

  CEBUG("priv_getNextOverlapFromDenovoQueue\n");

  oeI=PPF_overlap_edges_ptr->end();

  size_t qnum;

  // outer while to handle blacklisting
  while(true){
    for(qnum=0; qnum < PPF_queues.size() && oeI==PPF_overlap_edges_ptr->end(); ++qnum){
      CEBUG("Queue " << qnum << "\t" << PPF_queues[qnum].size() << endl);
      while(!PPF_queues[qnum].empty()){
	auto qe = PPF_queues[qnum].top();
	PPF_queues[qnum].pop();
	CEBUG("new qsize: " << PPF_queues[qnum].size() << endl);
	CEBUG("qe check " << PPF_readpool_ptr->getRead(qe.second->rid1).getName() << " "  << PPF_readpool_ptr->getRead(qe.second->linked_with).getName() << endl);

/*
  BaCh: 18.10.2014
  OK ... what on earth was I thinking???
  If a read was already added, we do not need to place the overlap edge into *any* queue
  BaCh: 26.11.2014
  Gaaaaahhhhhh ... so what was I thinking on 14.10. when I uncommented that block???
  Of course we need to reinsert rid1 into the denovo-queue if the current linked_with was
  already taken: there could still be *OTHER* reads. *bigsigh*
  If not, then low covered areas may have a premature stop lurking, even with 100% overlaps.
*/
	if(PPF_ids_added_oltype[qe.second->linked_with]){
	  CEBUG("going to insert " << qe.second->linked_with << "\t" << static_cast<uint16>(PPF_ids_added_oltype[qe.second->linked_with]) << " " << static_cast<uint16>((*PPF_used_ids_ptr)[qe.second->linked_with]) << endl);
	  BUGIFTHROW(((PPF_ids_added_oltype[qe.second->linked_with]>0)+(*PPF_used_ids_ptr)[qe.second->linked_with])==1,"Oooops, added by oltype and used ids do not agree? " << static_cast<uint16>(PPF_ids_added_oltype[qe.second->linked_with]) << " " << static_cast<uint16>((*PPF_used_ids_ptr)[qe.second->linked_with]) << endl);

	  // BaCh 04.03.2013
	  // what was I thinking when I had this?
	  //  || PPF_blacklisted_ids[qe.second->linked_with]){
	  // really a bad move as that may add blacklisted ids which are not in the contig!
	  size_t newqnum=priv_insertRIDIntoDenovoQueues(qe.second->rid1);
	  CEBUG("new qnum: " << qnum << " --> " << newqnum << endl);
	  if(newqnum<qnum) {
	    // if re-inserted in a higher-prio queue (i.e. due to a template partner having
	    //  been added to contig in the mean time), break handling of this queue
	    //  and have the outer for-loop restart at the higher-prio queue
	    qnum=newqnum-1; // -1 because of ++qnum in for-loop
	    break; // inner while
	  }
	}else if(!PPF_blacklisted_ids[qe.second->linked_with]){
	  oeI=qe.second; // this will stop the inner while
	  --qnum; // corrector: the for loop will increase qnum ("wrongly"), so correct for that
	  break;
	}

	if(!PPF_ids_added_oltype[qe.second->linked_with]
	   && !PPF_blacklisted_ids[qe.second->linked_with]){
	  oeI=qe.second; // this will stop the inner while
	  --qnum; // corrector: the for loop will increase qnum ("wrongly"), so correct for that
	  break;
	}

      }
    }
    if(oeI!=PPF_overlap_edges_ptr->end()) break;
    if(PPF_blacklist_queues.empty()) break;
    priv_munchBlacklist(true);
  }

  CEBUG("Returning from qnum " << qnum << ": " << *oeI << endl);
  FUNCEND();

  return static_cast<uint32>(qnum);
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * loop denovo
 *
 * Currently, this seems to handle genome and EST quite well
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::priv_ld_genome_and_est()
{
  FUNCSTART("void PPathfinder::priv_ld_genome()");

  Contig::templateguessinfo_t tguess;

#ifdef CLOCK_STEPS1
  timeval tv;

  timeval tvtotal;
  gettimeofday(&tvtotal,nullptr);
  suseconds_t ldtimear=0;
#endif

  struct tms mytms;
  times(&mytms);
  clock_t actclocks=mytms.tms_utime+mytms.tms_stime;
  clock_t maxallowedclocks=actclocks+PPF_pafparams_ptr->paf_maxcontigclockticks;


  std::vector<int8> &  lr_used_ids = *PPF_used_ids_ptr;

  nextreadtoadd_t nrta;
  decltype(PPF_overlap_edges_ptr->end()) oeI;

  // the forcegrow for addRead()
  //  0:  growth allowed
  //  >0: need to grow at least n bases (not used in newpathfinder)
  //  -1: no growth allowed (not used atm?)

  int32 forcegrow=0;

  bool buildprematurestop=false;

  // noaligncounter counts how many times the pathfinder itself skipped over calling Contig::addRead()
  //  when it saw it made no sense. Though every time a read is added to the contig, the counter is
  //  reset to 0
  // if noaligncounter reaches a threshold (2400 atm), we'll stop contig building
  //
  // This speeds up abortion a contig build process in hopeless situations,
  //  i.e., when both ends of a contig end in repeats which very probably
  //  cannot be crossed.

  uint32 noaligncounter=0;
  while(true){
#ifdef CLOCK_STEPS1
    gettimeofday(&tv,nullptr);
#endif
    nrta.foundqueuenum=priv_getNextOverlapFromDenovoQueue(oeI);
#ifdef CLOCK_STEPS1
    PPF_timing_pathsearch+=diffsuseconds(tv);
#endif
    CEBUG("nrta.foundqueuenum: " << static_cast<uint16>(nrta.foundqueuenum) << endl);
#ifndef PUBLICQUIET
    cout << "nrta.foundqueuenum: " << static_cast<uint16>(nrta.foundqueuenum) << endl;
#endif
    if(oeI==PPF_overlap_edges_ptr->end()) {
      CEBUG("breaking out\n");
      // queues are empty, this is a totally normal stopping of the build
      break;
    }
    nrta.refid=oeI->rid1;
    nrta.newid=oeI->linked_with;
    nrta.direction_newid=oeI->direction;
    nrta.ads_node=&(*PPF_adsfacts_ptr)[oeI->adsfindex];
    nrta.weight=oeI->best_weight;
    PPF_contigerrstat.reset();
    PPF_contigerrstat.code=Contig::ENOTCALLED;
    forcegrow=0;

#ifndef PUBLICQUIET
    //cout << "ADS Node: " << *(nrta.ads_node) << "\n";
#endif

    // for genomes, disallow growth of contig if ...
    bool doalign=true;
    if(PPF_pafparams_ptr->paf_use_genomic_algorithms && !PPF_actcontig_ptr->getLongRepeatStatus()){
      // this gets too dangerous, i.e., the overlap is repetitive and no
      //  paired-end in sight. BUT: only if this is a "regular" contig that
      //  was started from a non-multicopy read

      auto & newread = PPF_readpool_ptr->getRead(nrta.newid);
      auto tpid=newread.getTemplatePartnerID();
      if(newread.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){
	if(!(oeI->ol_stronggood || oeI->ol_weakgood)){
	  if(tpid>=0){
	    auto & tpread=PPF_readpool_ptr->getRead(tpid);
	    if(!PPF_ids_added_oltype[tpid]
	       && (tpread.hasFreqAvg() || !tpread.hasFreqRept())){
	      // so ... this Solexa overlap is not strong- or weakgood
	      //  and its template-partner has average frequency (or not
	      //  rept freq), but is not in the contig yet
	      // we could say "wait for a later timepoint to include this,
	      //  this may be safer" by simply doing
	      //    doalign=false;
	      //  but this can stop the building of a contig at places where,
	      //  for some reason or another, we have a "longer stretch with
	      //  less coverage than average. Longer means: > ~2x avg
	      //  read length.
	      //
	      // For this reason, we will take a detour:
	      // if overlap ratio=100 and overlap length>=35 and new read has no rept,
	      //  kmer we will still allow incorporation of read
	      // 35 maybe too small? On the other hand: we're already pretty low in the queues
	      //
	      // On the other hand, 35 may also be too large in projects with very
	      //  uneven coverage.
	      // Therefore, improved strategy:
	      //  - if tparther is already in contig: >=17bp, 100%
	      //  - if tparther not in contig
	      //     - overlap >= 35bp, accept
	      //     - if overlap <35bp and "normal" queues, relegate read (i.e.: pathfinder will try everything else before touching it again)
	      //     - if overlap <35bp and "relegated" queue, give it a try
	      if(!PPF_readpool_ptr->getRead(nrta.newid).hasFreqRept()       // maybe also check freq of ref?
		 && nrta.ads_node->getOverlapLen()>=17
		 && nrta.ads_node->getScoreRatio()==100){
		// this accepts all length>=17 where template partner is in contig
		// if not in contig, overlap must be >= 35
		if(PPF_ids_added_oltype[tpid] == 0
		   && nrta.ads_node->getOverlapLen()<35
		   && nrta.foundqueuenum<QTG_RELEGATEDBYPP){
		  doalign=false;
#ifndef PUBLICQUIET
		  cout << "FALSE1\n";
#endif
		}
	      }else{
		doalign=false;
#ifndef PUBLICQUIET
		cout << "FALSE2\n";
#endif
	      }
	    }else if(PPF_ids_added_oltype[tpid] != ADDED_BY_NOREPT){
	      forcegrow=-1;
	    }
	  }else{
	    if((!oeI->ol_norept || oeI->ol_rept)
	       && newread.hasKMerFork()){
	      forcegrow=-1;
	    }
	  }
	}
      }else if(nrta.foundqueuenum>=QTG_TPARTNERNOTREPT_OLREPT){
	if(!(oeI->ol_stronggood || oeI->ol_weakgood)){      // 05.11.2013: let's test that
	  forcegrow=-1;
	}else if(!oeI->ol_stronggood
		 && oeI->ol_weakgood
		 && newread.hasKMerFork()){
	  forcegrow=-1;
	}
      }

      // rule for rept(!) overlap with non-overlapping template partners
      // On 14.10.14, Added ol_stronggood to allow walking into rept area
      //  a bit more (for scaffolding). Hope this doesn't backfire.
      if(!oeI->ol_stronggood   // Added 14.10.14
	 && oeI->ol_rept
	 && tpid>=0
	 && tpid != nrta.refid){  // this is the non-pair overlap clause

	// if template partner not in contig (or added otherwise than by norept),
	// do not add
	if(PPF_ids_added_oltype[tpid]!=ADDED_BY_NOREPT){
	  doalign=false;
#ifndef PUBLICQUIET
	  cout << "FALSE3\n";
#endif
	}

/*
	// if template partner with a non-rept overlap exists
	//  but is not used yet, do not align!
	if((*PPF_hasnoreptoverlap_ptr)[tpid]
	   && !(*PPF_used_ids_ptr)[tpid]){
	  doalign=false;
	}else if((*PPF_hasnoreptoverlap_ptr)[tpid]
		 && (*PPF_used_ids_ptr)[tpid]
		 && PPF_ids_added_oltype[tpid]==ADDED_NOTADDED){
	  // if template partner with a non-rept overlap exists
	  //  but is not in this contig, do not align!
	  doalign=false;
	}else if(!(*PPF_hasnoreptoverlap_ptr)[tpid]
		 && !(*PPF_hasnoreptoverlap_ptr)[nrta.newid]){
	  // if both partner have no no-rept overlap (i.e., are completely
	  //  in a repeat), do not align
	  // TODO: 1. this of course leads to problems with PCR duplicates *sigh*
	  // TODO: 2. maybe better to blacklist that via ol_ flag?
	  doalign=false;
	}
//*/

      }
    }

#ifndef PUBLICQUIET
    cout << "doalign\t" << doalign << "\t" << PPF_readpool_ptr->getRead(nrta.newid).getName() << "\toechosen\tsg: " << oeI->ol_stronggood << " wg: " << oeI->ol_weakgood << " baf: " << oeI->ol_belowavgfreq << " nrp: " << oeI->ol_norept << " rep: " << oeI->ol_rept << endl;
#endif

    if(doalign){
#ifdef CLOCK_STEPS1
      gettimeofday(&tv,nullptr);
#endif
      ++PPF_readaddattempts;
      BUGIFTHROW(static_cast<uint16>(lr_used_ids[nrta.newid]),"PFcheck: newid already used??? " << nrta.newid << " " << static_cast<uint16>(lr_used_ids[nrta.newid]) << '\n');
      PPF_actcontig_ptr->addRead(*PPF_aligncache_ptr,
				 nrta.ads_node, nrta.refid, nrta.newid, nrta.direction_newid,
				 (*PPF_multicopies_ptr)[nrta.newid],
				 forcegrow,
				 tguess,
				 PPF_contigerrstat);
#ifdef CLOCK_STEPS1
      {
	auto tmpx=diffsuseconds(tv);
	PPF_timing_connadd+=tmpx;
	ldtimear+=tmpx;
      }
#endif
    }

    // in genome assemblies, some overlaps get one last chance
    if(!doalign
       && nrta.foundqueuenum<QTG_RELEGATEDBYPP
       && PPF_pafparams_ptr->paf_use_genomic_algorithms){
      PPF_queues[QTG_RELEGATEDBYPP].push(ppfweightelem_t(oeI->best_weight,oeI));
      doalign=true;
      PPF_contigerrstat.code=Contig::ERELEGATEDBYPP;
#ifndef PUBLICQUIET
      cout << "RELEGATED\n";
#endif
    }else{
      if(doalign && PPF_contigerrstat.code == Contig::ENOERROR){
	CEBUG("\nok, added" << endl);
	noaligncounter=0;
	PPF_ids_in_contig_list.push_back(nrta.newid);
	lr_used_ids[nrta.newid]=1;
	if(oeI->ol_norept){
	  PPF_ids_added_oltype[nrta.newid]=ADDED_BY_NOREPT;
	}else if(!oeI->ol_rept){
	  PPF_ids_added_oltype[nrta.newid]=ADDED_BY_NOTREPT;
	}else{
	  PPF_ids_added_oltype[nrta.newid]=ADDED_BY_OTHER;
	}
#ifdef CLOCK_STEPS1
	gettimeofday(&tv,nullptr);
#endif
	priv_insertRIDIntoDenovoQueues(oeI->linked_with);
#ifdef CLOCK_STEPS1
	PPF_timing_pathsearch+=diffsuseconds(tv);
#endif
	priv_munchBlacklist(false);
	//cout << "\nTGUESS " << tguess << endl;
	priv_storeTemplateGuess(nrta.newid,tguess);
      }else{
	CEBUG("\nnot added" << endl);
	priv_handleReadNotAligned(oeI,nrta);
	if(!doalign) ++noaligncounter;
      }
      CEBUG("\npost ar" << endl);
      // keep this here so that a possible ban above is taken into account
#ifdef CLOCK_STEPS1
      gettimeofday(&tv,nullptr);
#endif
      priv_insertRIDIntoDenovoQueues(oeI->rid1);
#ifdef CLOCK_STEPS1
      PPF_timing_pathsearch+=diffsuseconds(tv);
#endif
    }

    priv_showProgress();

    if(PPF_pafparams_ptr->paf_use_max_contig_buildtime && PPF_buildcontig_newlinecounter==0){
      times(&mytms);
      actclocks=mytms.tms_utime+mytms.tms_stime;
      if(actclocks>maxallowedclocks){
	cout << "\nMaximum build time for this contig reached, aborting build.\n";
	buildprematurestop=true;
	break;
      }
    }
    if(noaligncounter>=4800){
      cout << "\nProbable dead end, aborting build.\n";
      buildprematurestop=true;
      break;
    }
  }

  // in case we stopped the build prematurely
  // cleanup the priority queues ... we need to be tidy
  if(buildprematurestop){
    for(auto & pq : PPF_queues){
      while(!pq.empty()) pq.pop();
    }
  }


#ifdef CLOCK_STEPS1
  cout << "priv_ld addRead: " << ldtimear << endl;
  cout << "priv_ld total: " << diffsuseconds(tvtotal) << endl;
#endif

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::priv_handleReadNotAligned(necontainer_t::iterator oeI, nextreadtoadd_t const &nrta)
{
  FUNCSTART("void Pathfinder::priv_handleReadNotAligned(nextreadtoadd_t const &nrta)");

//  CEBUG("Banning " << oeI->rid1 << " (" << PPF_readpool_ptr->getRead(oeI->rid1).getName() << ")\t" << oeI->linked_with << " (" << PPF_readpool_ptr->getRead(oeI->linked_with).getName() << ")" << endl);

  // always, always ban the overlap which was not aligned
  // (or else endless loop possible in mapping)
  oeI->pf_banned=true;
  if(PPF_overlapsbanned_smallstore.size()<PPF_overlapsbanned_smallstore.capacity()){
    PPF_overlapsbanned_smallstore.push_back(oeI);
  }

  // Now look at the blocks of newedges_t whose reads are given
  //  by the vector and who overlap with nrta.newid
  auto raI=PPF_contigerrstat.reads_affected.begin();
  for(; raI != PPF_contigerrstat.reads_affected.end(); ++raI){
    //cout << "cesra: " << *raI << "\t" << PPF_readpool_ptr->getRead(*raI).getName() << endl;
    auto neI = (*PPF_lowerbound_oedges_ptr)[*raI];
    for(; neI != PPF_overlap_edges_ptr->end() && neI->rid1 == *raI; ++neI){
      if(!neI->pf_banned){
	if(neI->linked_with == nrta.newid){
	  //cout << "Banning2: " << neI->rid1 << '\t' << neI->linked_with << '\n';
	  neI->pf_banned=true;
	  if(PPF_overlapsbanned_smallstore.size()<PPF_overlapsbanned_smallstore.capacity()){
	    PPF_overlapsbanned_smallstore.push_back(neI);
	  }
	  //break;
	  // NOTE: should MIRA ever get multiple overlaps between
	  //  two reads, then the "break;" above must be removed
	  //  should be safe for now.
	  // BaCh 18.10.2014: this is now the case with overlap storage for small kmers
	}
      }
    }
  }

  if(PPF_blacklist_queues.empty()) PPF_blacklist_queues.push(std::vector<readid_t>());
  PPF_blacklist_queues.back().push_back(nrta.newid);
  PPF_blacklisted_ids[nrta.newid]=1;

  FUNCEND();
  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::priv_munchBlacklist(bool force)
{
  if((force && !PPF_blacklisted_ids.empty())
     || (PPF_ids_in_contig_list.size()%8 == 0
	 && PPF_blacklist_queues.size()>=10)){
    CEBUG("Munching start force("<<force<<") blacklist front with " << PPF_blacklist_queues.front().size() << " elements\n");
    for(auto rid : PPF_blacklist_queues.front()){
      PPF_blacklisted_ids[rid]=0;
    }
    uint64 dmok=0;
    uint64 dmnok=0;
    for(auto rid : PPF_blacklist_queues.front()){
      auto oeI=(*PPF_lowerbound_oedges_ptr)[rid];
      for(;oeI!=PPF_overlap_edges_ptr->end() && rid==oeI->rid1;++oeI){
	if(!oeI->pf_banned && PPF_ids_added_oltype[oeI->linked_with]){
	  auto rv=priv_insertRIDIntoDenovoQueues(oeI->linked_with);
	  if(rv!=QTG_END) {
	    CEBUG("demunch success " << *oeI << endl; ++dmok);
	  }else{
	    CEBUG("demunch fail " << *oeI << endl; ++dmnok);
	  }
	}
      }
    }
    PPF_blacklist_queues.pop();
    CEBUG("Munching end. ok: " << dmok << " nok: " << dmnok << endl);
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void PPathfinder::map()
{
  FUNCSTART("void PPathfinder::map()");

#ifdef CLOCK_STEPS1
  timeval tv;
#endif

  Contig::templateguessinfo_t tguess;

  cout << "Backbone assembly to " << PPF_actcontig_ptr->getContigName() << endl;
  for(auto qi=0; qi<PPF_queues.size(); ++qi){
    BUGIFTHROW(!PPF_queues[qi].empty(),"Queue chk 1 " << qi << " not empty?");
  }

  BUGIFTHROW(PPF_actcontig_ptr->getNumReadsInContig()==0,"No reads in contig " << PPF_actcontig_ptr->getContigName() << " ??? How should I map something to nothing?");

  priv_basicSetup();
  priv_prepareRailOverlapCache();

  nextreadtoadd_t nrta;
  nrta.newid=0;

  bool allowbbqmulticopies=false;
  bool allowbbqtroublemakers=false;
  bool allowbbqsmallhits=false;

  while(nrta.newid >= 0) {
    nrta.refid=-1;
    nrta.newid=-1;
    nrta.weight=0;
    nrta.direction_newid=0;
    nrta.ads_node=nullptr;

#ifdef CLOCK_STEPS1
    gettimeofday(&tv,nullptr);
#endif
    auto oeI=priv_findNextBackboneOverlapQuick(nrta,
					       allowbbqmulticopies,
					       allowbbqtroublemakers,
					       allowbbqsmallhits);

    if(nrta.newid < 0) {
      priv_prepareRailOverlapCache();
      oeI=priv_findNextBackboneOverlapQuick(nrta,
					    allowbbqmulticopies,
					    allowbbqtroublemakers,
					    allowbbqsmallhits);


      if(nrta.newid < 0) {
	CEBUG("allow everything\n");
	priv_prepareRailOverlapCache();
	allowbbqmulticopies=true;
	allowbbqtroublemakers=true;
	allowbbqsmallhits=true;
	oeI=priv_findNextBackboneOverlapQuick(nrta,
					      allowbbqmulticopies,
					      allowbbqtroublemakers,
					      allowbbqsmallhits);
      }
    }
#ifdef CLOCK_STEPS1
    PPF_timing_pathsearch+=diffsuseconds(tv);
    gettimeofday(&tv,nullptr);
#endif

    if(nrta.newid >= 0) {
      ++PPF_readaddattempts;
      PPF_actcontig_ptr->addRead(*PPF_aligncache_ptr,
				 nrta.ads_node, nrta.refid, nrta.newid, nrta.direction_newid,
				 (*PPF_multicopies_ptr)[nrta.newid],
				 0,
				 tguess,
				 PPF_contigerrstat);
      if(PPF_contigerrstat.code == Contig::ENOERROR) {
	CEBUG("\nok, added" << endl);
	//Contig::setCoutType(Contig::AS_TEXT);
	//cout << "nrnrnrnrnrnr\n" << *PPF_actcontig_ptr << endl;
	PPF_ids_in_contig_list.push_back(nrta.newid);
	(*PPF_used_ids_ptr)[nrta.newid]=1;
	PPF_ids_added_oltype[nrta.newid]=ADDED_BY_BACKBONE;
	priv_storeTemplateGuess(nrta.newid,tguess);
      }else{
	priv_handleReadNotAligned(oeI,nrta);
      }
#ifdef CLOCK_STEPS1
      PPF_timing_connadd+=diffsuseconds(tv);
#endif
      priv_showProgress();
    }
  }

  for(auto qi=0; qi<PPF_queues.size(); ++qi){
    BUGIFTHROW(!PPF_queues[qi].empty(),"Queue chk 2 " << qi << " not empty?");
  }

}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void PPathfinder::priv_prepareRailOverlapCache()
{
#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  // initialise PPF_railoverlapcache
  // each read that has an overlap with a rail is put into the
  //  vector once;
  PPF_railoverlapcache.clear();

  if(PPF_tmpproc_readalreadyrailed.empty()){
    PPF_tmpproc_readalreadyrailed.resize(PPF_used_ids_ptr->size(),0);
  }

  priv_prochelper1(ReadGroupLib::SEQTYPE_PACBIOHQ);
  priv_prochelper1(ReadGroupLib::SEQTYPE_SANGER);
  priv_prochelper1(ReadGroupLib::SEQTYPE_454GS20);
  priv_prochelper1(ReadGroupLib::SEQTYPE_SOLEXA);
  priv_prochelper1(ReadGroupLib::SEQTYPE_IONTORRENT);
  priv_prochelper1(ReadGroupLib::SEQTYPE_TEXT);
  // take the rest in case there should be some
  priv_prochelper1(ReadGroupLib::SEQTYPE_END);

  // done here, need to clean up our tmp data so that the next
  //  call finds it clean
  for(auto rid : PPF_railoverlapcache) PPF_tmpproc_readalreadyrailed[rid]=0;

#ifndef PUBLICQUIET
  cout << "Backbone overlap cache: " << PPF_railoverlapcache.size() << endl;
#endif
}

void PPathfinder::priv_prochelper1(uint8 seqtype)
{
  if(seqtype != ReadGroupLib::SEQTYPE_END
     && !ReadGroupLib::hasLibWithSeqType(seqtype)) return;

  for(auto rid : PPF_rails_in_contig_list){
    auto rcI=(*PPF_lowerbound_oedges_ptr)[rid];
    //if(!allowedrefids.empty() && !allowedrefids[rid]) continue;
    for(; rcI != PPF_overlap_edges_ptr->end() && rcI->rid1==rid; ++rcI){
      if(rcI->pf_banned) continue;
      if(seqtype != ReadGroupLib::SEQTYPE_END
	 && PPF_readpool_ptr->getRead(rcI->linked_with).getSequencingType()!=seqtype) continue;
      if((*PPF_used_ids_ptr)[rcI->linked_with]) continue;
      if(PPF_tmpproc_readalreadyrailed[rcI->linked_with]) continue;
      if(PPF_readpool_ptr->getRead(rcI->linked_with).isRail()) continue;

      // evil little rule for clean overlap ends ...
      if(PPF_wantscleanoverlapends > 0){
	CEBUG("PPFC check" << (*PPF_adsfacts_ptr)[rcI->adsfindex] << endl);
	if((*PPF_adsfacts_ptr)[rcI->adsfindex].get5pLenContiguousMatch(rcI->linked_with) <= PPF_wantscleanoverlapends
	   || (*PPF_adsfacts_ptr)[rcI->adsfindex].get3pLenContiguousMatch(rcI->linked_with) <= PPF_wantscleanoverlapends) {
	  CEBUG("PPFC denied\n");
	  continue;
	}
      }
      CEBUG("PPFC accepted\n");
      // ... and for minimum total non-matches
      if(PPF_mintotalnonmatches > 0
	 && (*PPF_adsfacts_ptr)[rcI->adsfindex].getTotalNonMatches() <= PPF_mintotalnonmatches){
	continue;
      }
      // ... and for allowed seqtype
      if(PPF_allowedseqtype!=ReadGroupLib::SEQTYPE_END
	 && PPF_readpool_ptr->getRead(rcI->linked_with).getSequencingType() != PPF_allowedseqtype) continue;

      PPF_tmpproc_readalreadyrailed[rcI->linked_with]=1;
      PPF_railoverlapcache.push_back(rcI->linked_with);
    }
  }
}

//#define CEBUG(bla)



/*************************************************************************
 *
 * Take one non-rail, find best place for a it in existing contig
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
necontainer_t::iterator PPathfinder::priv_findNextBackboneOverlapQuick(nextreadtoadd_t & resultread, bool allowmulticopies, bool allowtroublemakers, bool allowsmallhits)
{
  FUNCSTART("necontainer_t::iterator PPathfinder::priv_findNextBackboneOverlapQuick(nextreadtoadd_t & resultread, bool allowmulticopies, bool allowtroublemakers, bool allowsmallhits)");

  CEBUG("aaa: " << allowmulticopies << allowtroublemakers << allowsmallhits << endl);

  auto retoeI=PPF_overlap_edges_ptr->end();
  if(PPF_railoverlapcache.empty()){
    CEBUG("\troc empty\n");
    // found nothing in previous loop, exit
    return retoeI;
  }

  const std::vector<int8> & lr_used_ids = *PPF_used_ids_ptr;
  const std::vector<uint8> & lr_multicopies = *PPF_multicopies_ptr;
  const std::vector<uint8> & lr_istroublemaker = *PPF_istroublemaker_ptr;

  resultread.newid=-1;
  while(resultread.newid<0){
    bool continuesearch=true;
    int32 readid=-1;
    while(continuesearch && !PPF_railoverlapcache.empty()){
      readid=PPF_railoverlapcache.front();

      CEBUG("\nfnboq:\n");
      CEBUG("l: " << PPF_readpool_ptr->getRead(readid).getName());
      CEBUG("\tused: " << (int16) lr_used_ids[readid]);
      CEBUG("\tmc: " << (int16) lr_multicopies[readid]);
      CEBUG("\ttm: " << (int16) lr_istroublemaker[readid]);
      //if(!allowedrefids.empty()){
      //	CEBUG("\tar: " << (int16) allowedrefids[readid]);
      //}

      PPF_railoverlapcache.pop_front();

      if(lr_used_ids[readid]==0
	 && (allowmulticopies || lr_multicopies[readid]==0)
	 && (allowtroublemakers || lr_istroublemaker[readid]==0)){
	continuesearch=false;
      }
    }


    if(continuesearch){
      CEBUG("\tnot taken/found\n");
      // found nothing in previous loop, exit
      return retoeI;
    }
    BUGIFTHROW(readid<0,"readid " << readid << " < 0, shouldn't be here");

    CEBUG("\ttaken/found");

    // ok, found one. Now search the rail it fits best to *in this contig*
    auto oeI=(*PPF_lowerbound_oedges_ptr)[readid];
    for(; oeI!=PPF_overlap_edges_ptr->end() && oeI->rid1 == readid; ++oeI){
      // don't bother looking at if overlap is banned
      if(oeI->pf_banned) {CEBUG("\tbanned"); continue;}
      // must link to rail
      if(!PPF_readpool_ptr->getRead(oeI->linked_with).isRail()) {CEBUG("\tbanned"); continue;}
      // rail must be in this contig!
      if(!PPF_ids_added_oltype[oeI->linked_with]) {CEBUG("\tbanned"); continue;}
      //// rail must be allowed as refid
      //if((!allowedrefids.empty() && !allowedrefids[oeI->linked_with])) continue;

      // if necessary, take into account only rails that are
      //  - non-multicopies
      //  - non-troublemakers
      CEBUG("\tbasicok");
      if((allowmulticopies || lr_multicopies[oeI->linked_with]==0)
	 && (allowtroublemakers || lr_istroublemaker[oeI->linked_with] == 0)){
	CEBUG("\tmc&tm ok");
	//  - that have overlap length >= minim length (just to have
	//    good matches first)
	if(oeI->best_weight > resultread.weight){
	  CEBUG("\tbw ok");
	  if(allowsmallhits
	     || (*PPF_adsfacts_ptr)[oeI->adsfindex].getOverlapLen() >=
	     (*PPF_miraparams_ptr)[PPF_readpool_ptr->getRead(readid).getSequencingType()].getPathfinderParams().paf_bbquickoverlap_minlen) {
	    resultread.refid=oeI->linked_with;
	    resultread.newid=readid;
	    resultread.weight=oeI->best_weight;
	    resultread.direction_newid=oeI->direction;
	    resultread.ads_node=&(*PPF_adsfacts_ptr)[oeI->adsfindex];

	    // the edges are sorted by rid1, then by "bestweight: high to low"
	    // therefore, if we wound something, we do not need to check further
	    // edges of this read ... all remaining edges have lower or equal
	    // quality anyway. Therefore, break out of loop here
	    retoeI=oeI;
	    CEBUG("\treturned\n");
	    break;
	  }
	}
      }
    }
  }
  return retoeI;
}

//#define CEBUG(bla)

void PPathfinder::priv_storeTemplateGuess(readid_t newid, Contig::templateguessinfo_t & tguess)
{
  FUNCSTART("void PPathfinder::priv_storeTemplateGuess(readid_t refid, readid_t newid, Contig::templateguessinfo_t & tguess)");

  //cout << "STORE STORE STORE " << tguess << endl;

  if(PPF_astemplateguesses_ptr->empty()
     || tguess.rgid.isDefaultNonValidReadGroupID()) return;

  auto tpid=PPF_readpool_ptr->getRead(newid).getTemplatePartnerID();
  BUGIFTHROW(tpid == -1, "tpid == -1?");
  if( PPF_pafparams_ptr->paf_use_genomic_algorithms       // thorough check on repeats only for genome, not possible for EST
      && PPF_ids_added_oltype[newid]!=ADDED_BY_BACKBONE
      && (PPF_ids_added_oltype[tpid]!=ADDED_BY_NOREPT
	  || PPF_ids_added_oltype[newid]!=ADDED_BY_NOREPT)){
    return;
  }

  //cout << "REALLYSTORE\n";

  BUGIFTHROW(PPF_readpool_ptr->getRead(tpid).getTemplateID() != PPF_readpool_ptr->getRead(newid).getTemplateID(), "PPF_readpool_ptr->getRead(tpid).getTemplateID() " << PPF_readpool_ptr->getRead(tpid).getTemplateID() << " != " << PPF_readpool_ptr->getRead(newid).getTemplateID() << " PPF_readpool_ptr->getRead(newid).getTemplateID() ???");
  (*PPF_astemplateguesses_ptr)[PPF_readpool_ptr->getRead(newid).getTemplateID()]=tguess;
}
