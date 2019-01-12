/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2002 and later by Bastien Chevreux
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


#include "mira/newpathfinder.H"
#include "util/stlimprove.H"

#include <iostream>

// for time measurements
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include <unistd.h>

#include "errorhandling/errorhandling.H"



using namespace std;



//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}


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
 *
 *
 *
 *************************************************************************/

ostream & operator<<(ostream &ostr, nextreadtoadd_t const &nrta)
{
#ifndef PUBLICQUIET
  if(nrta.newid == -1) {
    ostr << "Invalid / empty nrta.\n";
    return ostr;
  }
#endif

  //ostr <<"RefID : " << nrta.refid << endl;
  //ostr <<"NewID : " << nrta.newid << endl;
  //ostr <<"Weight: " << nrta.weight << endl;
  //ostr <<"Dir: " << nrta.direction_newid << endl;
  //ostr <<"ADS: " << *(nrta.ads_node) << endl;

#ifndef PUBLICQUIET
  switch(nrta.foundmethod) {
  case Pathfinder::FOUND_NORMAL_AND_SRMB : {
    ostr << "FOUND_NORMAL_AND_SRMB ";
    break;
  }
  case Pathfinder::FOUND_NORMAL : {
    ostr << "FOUND_NORMAL ";
    break;
  }
  case Pathfinder::FOUND_SPEEDYMAP : {
    ostr << "FOUND_SPEEDYMAP ";
    break;
  }
  case Pathfinder::FOUND_NONMULTICOPY : {
    ostr << "FOUND_NONMULTICOPY ";
    break;
  }
  case Pathfinder::FOUND_TEMPLATE_AND_TPARTNER : {
    ostr << "FOUND_TEMPLATE_AND_TPARTNER ";
    break;
  }
  //case Pathfinder::FOUND_TEMPLATE_AND_OVERLAPPINGTPARTNER : {
  //  ostr << "FOUND_TEMPLATE_AND_OVERLAPPINGTPARTNER ";
  //  break;
  //}
  case Pathfinder::FOUND_TPARTNER_IN_CONTIG : {
    ostr << "FOUND_TPARTNER_IN_CONTIG ";
    break;
  }
  case Pathfinder::FOUND_OVERLAPPING_PARTNER : {
    ostr << "FOUND_OVERLAPPING_PARTNER ";
    break;
  }
  case Pathfinder::FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER : {
    ostr << "FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER ";
    break;
  }
  case Pathfinder::FOUND_RAIL_AND_NONMULTICOPY : {
    ostr << "FOUND_RAIL_AND_NONMULTICOPY ";
    break;
  }
  case Pathfinder::FOUND_RAIL_AND_MULTICOPY : {
    ostr << "FOUND_RAIL_AND_MULTICOPY ";
    break;
  }
  case Pathfinder::FOUND_STAGE1 : {
    ostr << "FOUND_STAGE1 ";
    break;
  }
  case Pathfinder::FOUND_STAGE2 : {
    ostr << "FOUND_STAGE2 ";
    break;
  }
  case Pathfinder::FOUND_STAGE100 : {
    ostr << "FOUND_STAGE100 ";
    break;
  }
  default : {
    ostr << "FOUND_??? unknown find method? ";
  }
  }

  ostr << nrta.weight << endl;

  ostr.flush();
#endif

  return ostr;
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Pathfinder::Pathfinder(vector<MIRAParameters> * params, ReadPool & readpool, vector<newedges_t> & forward_edges, vector<AlignedDualSeqFacts> & adsfacts)
{
  FUNCSTART("");

  PAF_miraparams=params;
  PAF_forward_edges=&forward_edges;
  PAF_readpool=&readpool;
  PAF_adsfacts=&adsfacts;

  // we don't know whether edges were banned ... play it safe
  //  and say there were so that the pf_banned flags get reset
  //  at first call of constructStepByStep()
  PAF_overlapsbanned=1;

  // give the small store iterators to banned overlaps a capacity of
  //  50k entries.
  PAF_overlapsbanned_smallstore.reserve(50000);

  PAF_uiitr_is_clean=false;
  PAF_usesbeststartcache=false;

  PAF_valid=1;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Pathfinder::~Pathfinder()
{
  FUNCSTART("Pathfinder::~Pathfinder()");

  discard();

  FUNCEND();
}

// Copy constructor
//  no discard needed as this object will be freshly created when
//  called through this constructor
//Pathfinder::Pathfinder(Pathfinder const &other)
//{
//  FUNCSTART("Pathfinder::Pathfinder(Pathfinder const &other)");
//
//  PAF_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}

// Copy operator, needed by copy-constructor
//Pathfinder const & Pathfinder::operator=(Pathfinder const & other)
//{
//  FUNCSTART("Pathfinder const & Pathfinder::operator=(Pathfinder const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, Pathfinder const &PAF)
//{
//  FUNCSTART("friend ostream & Pathfinder::operator<<(ostream &ostr, const  &PAF)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}

void Pathfinder::discard()
{
  FUNCSTART("Pathfinder::discard()");

  if(PAF_valid==0){
  }

  //PAF_lowerbound_oedges.resize(0);

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
bool Pathfinder__compareNewEdges_t_(const newedges_t & a,
				    const newedges_t & b);
bool Pathfinder__compareNewEdges_t_(const newedges_t & a, const newedges_t & b)
{
  return a.rid1 < b.rid1;
}

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Pathfinder::constructStepByStep(vector<Align> & aligncache, vector<int8> * used_ids, vector<int32> * ids_in_contig, vector<uint8> * multicopies, vector<uint8> * hasmcoverlaps, vector<uint8> * istroublemaker, vector<vector<newedges_t>::iterator> * lowerbound_oedges_ptr, Contig & con)
{
  FUNCSTART("void Pathfinder::constructFast()");

#ifndef PUBLICQUIET
  cout << "Old constructStepByStep().\n";
#endif

  PAF_pafparams=(*PAF_miraparams)[0].getPathfinderParams();

  PAF_used_ids_ptr=used_ids;
  PAF_ids_in_contig_ptr=ids_in_contig;
  PAF_multicopies_ptr=multicopies;
  PAF_hasmcoverlap_ptr=hasmcoverlaps,
  PAF_istroublemaker_ptr=istroublemaker;
  PAF_lowerbound_oedges_ptr=lowerbound_oedges_ptr;

  PAF_usesbeststartcache=false;

  PAF_railsincontig.clear();

  n3_basicCSBSSetup();


  if(con.getContigLength()==0){
#ifndef PUBLICQUIET
    cout << "\nStarted with empty contig.\n";
#endif

    CEBUG("csbs 3a\n");
    int32 startid=searchBestEnvironment();
    if(startid>=0) {
      CEBUG("csbs 3a1\n");

      // if we keep long repeats separated:
      //   if we start with a non-multicopy read, forbid
      //    multicopy/multicopy overlaps
      if((*PAF_multicopies_ptr)[startid]>0) {
#ifndef PUBLICQUIET
	cout << "\nStarted with multicopy.\n";
#endif
	con.setLongReapeatStatus(true);
      }else{
	cout << "\nStarted with non-multicopy.\n";
	if(PAF_musthonourmcmcflag){
#ifndef PUBLICQUIET
	  cout << "\nMay not add mc/mc overlaps.\n";
#endif
	  PAF_mayaddmcmcoverlaps=false;
	}
      }
#ifndef PUBLICQUIET
      cout << "mayaddmcmcoverlaps: " << PAF_mayaddmcmcoverlaps << endl;
#endif

      // BaCh. 05.12.2008
      // gcc version 4.1.2 20061115 (prerelease)
      // the above if/else is sometimes not executed (neither
      //  if nor else) when the buildContig() from below is above
      //  the if/else
      // Hunting that down s*cked big time.

      CEBUG("csbs 3a2\n");
      buildContig(aligncache,startid, con);
      CEBUG("csbs 3a3\n");

    } else {
      CEBUG("csbs 3a4\n");
      cout << "No more reads found to build a new contig!" << endl;
      FUNCEND();
      return;
    }
  } else {
#ifndef PUBLICQUIET
    cout << "\nStarted with existing contig.\n";
#endif
    CEBUG("csbs 3b\n");

    auto & crv=con.getContigReads();
    auto crI=crv.begin();
    for(; crI != crv.end(); ++crI){
      if(crI.getORPID() >= 0
	 && crI->isRail()){
	PAF_railsincontig.push_back(crI.getORPID());
      }
    }

    buildContig(aligncache, -1, con);
  }

  CEBUG("csbs 4\n");

  uint32 numreadspresent=PAF_used_ids_ptr->size();

  bool canloop=true;
  uint32 searchloop=1;
  for(;canloop && searchloop <30; searchloop++){
    cout << "\nRL" << searchloop << endl;
    zerovector(PAF_used_ids_in_this_run);
    PAF_blacklisted_ids.clear();
    PAF_blacklisted_ids.resize(numreadspresent,1);
    uint32 numreadsincontig=con.getNumReadsInContig();
    buildContig(aligncache, -1, con);
    if(numreadsincontig==con.getNumReadsInContig()) {
      cout << "\nThat's it for this contig." << endl;
      canloop=false;
    } else {
      cout << "\nGained " << con.getNumReadsInContig()-numreadsincontig << " reads.\n";
    }
  }

  if(canloop) {
    cout << "Stopped because reached maximum iterations allowed.\n";
  }

  FUNCEND();
}
//#define CEBUG(bla)







/*************************************************************************
 *
 * some basic setup for constructStepByStep()
 *
 * last remnants of n3_* routines ... see whether to rewrite
 *  or rename n4_*
 *
 *************************************************************************/


void Pathfinder::n3_basicCSBSSetup()
{
  FUNCSTART("void Pathfinder::constructFast()");

  PAF_pafparams=(*PAF_miraparams)[0].getPathfinderParams();

  uint32 numreadspresent=PAF_used_ids_ptr->size();
  PAF_used_ids_in_this_run.clear();
  PAF_used_ids_in_this_run.resize(numreadspresent,0);
  PAF_blacklisted_ids.clear();
  PAF_blacklisted_ids.resize(numreadspresent,1);
  PAF_readaddattempts=0;
  PAF_nonmulticopiesadded=0;

  PAF_musthonourmcmcflag=(*PAF_miraparams)[0].getAssemblyParams().as_keep_long_repeats_separated;
  PAF_mayaddmcmcoverlaps=true;

  CEBUG("csbs 1\n");

  // initialise all fields in the newedge_t structure of PAF_forward_edges
  //  used by the pathfinder
  // at the moment this is the pf_banned field (set when an overlap is
  //  rejected also set for all all reads that are already aligned
  //  in the contig at the potential insertion position

  if(PAF_overlapsbanned){
    if(!PAF_overlapsbanned_smallstore.empty()
       && PAF_overlapsbanned_smallstore.size() < PAF_overlapsbanned_smallstore.capacity()){
#ifndef PUBLICQUIET
      cout << "Clear pf_banned quick: " << PAF_overlapsbanned_smallstore.size() << '\n';
#endif
      vector<vector<newedges_t>::iterator>::iterator obssI=PAF_overlapsbanned_smallstore.begin();
      for(; obssI != PAF_overlapsbanned_smallstore.end(); obssI++){
	(*obssI)->pf_banned=false;
      }
    }else{
#ifndef PUBLICQUIET
      cout << "Clear pf_banned full: " << PAF_overlapsbanned << '\n';
#endif
      vector<newedges_t>::iterator neI=PAF_forward_edges->begin();
      for(; neI != PAF_forward_edges->end(); neI++){
	neI->pf_banned=false;
      }
    }
    PAF_overlapsbanned=0;
    PAF_overlapsbanned_smallstore.clear();
  }

  CEBUG("csbs 2\n");

  // initialise the quick lowerbound_oedges lookup vector if not already done
  if(PAF_lowerbound_oedges_ptr->empty()){
    vector<vector<newedges_t>::iterator> & lowerbound_oedges=*PAF_lowerbound_oedges_ptr;
    lowerbound_oedges.resize(PAF_readpool->size());
    newedges_t tmp;
    for(uint32 i=0; i<lowerbound_oedges.size(); i++) {
      //lowerbound_oedges[i]=PAF_forward_edges->lower_bound(i);
      tmp.rid1=i;
      lowerbound_oedges[i]=lower_bound(PAF_forward_edges->begin(),
				       PAF_forward_edges->end(),
				       tmp,
				       Pathfinder__compareNewEdges_t_);
    }
  }

  CEBUG("csbs 3\n");

  // initialise the has_multicopies_overlap if not already done
  if(PAF_hasmcoverlap_ptr->empty()){
    PAF_hasmcoverlap_ptr->resize(numreadspresent,0);
    vector<newedges_t>::const_iterator feI=PAF_forward_edges->begin();
    for(; feI != PAF_forward_edges->end(); feI++){
      if((*PAF_multicopies_ptr)[feI->rid1]
	 || (*PAF_multicopies_ptr)[feI->linked_with]){
	(*PAF_hasmcoverlap_ptr)[feI->rid1]=1;
	(*PAF_hasmcoverlap_ptr)[feI->linked_with]=1;
      }
    }
  }

  //PAF_idswithSRMBtags.clear();
  //PAF_idswithSRMBtags.resize(numreadspresent,0);
  //
  //for(uint32 i=0; i < numreadspresent; i++) {
  //  if(PAF_readpool->getRead(i).hasTag(Read::REA_tagSRMr)) PAF_idswithSRMBtags[i]=1;
  //}

}



/*************************************************************************
 *
 *
 *
 * last remnants of n3_* routines ... see whether to rewrite
 *  or rename n4_*
 *
 *************************************************************************/

void Pathfinder::n3_handleReadNotAligned(Contig::errorstatus_t & contigerrstat, nextreadtoadd_t const &nrta)
{
  FUNCSTART("void Pathfinder::n3_handleReadNotAligned(Contig::errorstatus_t & contigerrstat)");

#ifdef PUBLICQUIET
  contigerrstat.dumpStatus();
#else
  contigerrstat.dumpStatus(true,"n3_h ");
#endif

  // ban all alignments of that read with other reads that were
  //  at the same position

  //// debug: show all reads at that pos that the contig gave back
  //{
  //  cout << "\nReads affected:\n";
  //  vector<int32>::const_iterator raI;
  //  raI=contigerrstat.reads_affected.begin();
  //  for(; raI != contigerrstat.reads_affected.end(); raI++){
  //    cout << *raI << '\t' << (*PAF_readpool)[*raI].getName();
  //    if((*PAF_readpool)[*raI].isBackbone()){
  //      cout << "\tbackbone\n";
  //    } else if((*PAF_readpool)[*raI].isRail()){
  //      cout << "\trail\n";
  //    }else{
  //      cout << "\tseqtype: " << static_cast<uint16>((*PAF_readpool)[*raI].getSequencingType()) << '\n';
  //    }
  //  }
  //}


  {
    // First search look at the newedges_t block of nrta,newid overlaps
    //  and ban those overlaps whose reads are in the vector given back
    //  by the contig

    vector<newedges_t>::iterator neI = (*PAF_lowerbound_oedges_ptr)[nrta.newid];
    vector<int32>::const_iterator raI;
    for(; neI != PAF_forward_edges->end() && neI->rid1 == nrta.newid; neI++){
      if(!neI->pf_banned){
	raI=contigerrstat.reads_affected.begin();
	for(; raI != contigerrstat.reads_affected.end(); raI++){
	  if(*raI == neI->linked_with){
	    //cout << "Banning1: " << neI->rid1 << '\t' << neI->linked_with << '\n';
	    PAF_overlapsbanned++;
	    neI->pf_banned=true;
	    if(PAF_overlapsbanned_smallstore.size()<PAF_overlapsbanned_smallstore.capacity()){
	      PAF_overlapsbanned_smallstore.push_back(neI);
	    }
	    break;
	  }
	}
      }
    }

    // Now look at the blocks of newedges_t whose reads are given
    //  by by the vector and who overlap with nrta.newid
    raI=contigerrstat.reads_affected.begin();
    for(; raI != contigerrstat.reads_affected.end(); raI++){
      neI = (*PAF_lowerbound_oedges_ptr)[*raI];
      for(; neI != PAF_forward_edges->end() && neI->rid1 == *raI; neI++){
	if(!neI->pf_banned){
	  if(neI->linked_with == nrta.newid){
	    //cout << "Banning2: " << neI->rid1 << '\t' << neI->linked_with << '\n';
	    PAF_overlapsbanned++;
	    neI->pf_banned=true;
	    if(PAF_overlapsbanned_smallstore.size()<PAF_overlapsbanned_smallstore.capacity()){
	      PAF_overlapsbanned_smallstore.push_back(neI);
	    }
	    break;
	    // NOTE: should MIRA ever get multiple overlaps between
	    //  two reads, then the "break;" above must be removed
	    //  should be safe for now.
	  }
	}
      }
    }
  }

  FUNCEND();
  return;
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// returns ID of read with best matching env
// TODO: length >50?
uint32 Pathfinder::searchBestEnvironment()
{
  enum {SINGLE=0,
	NORMAL, NORMAL_NONMULTYCOPY,
	WITH_TEMPLATE, WITH_TEMPLATE_NONMULTYCOPY,
	WITH_TEMPLATE_PARTNER, WITH_TEMPLATE_PARTNER_NONMULTICOPIES,
	ENUMEND};

  FUNCSTART("void Pathfinder::searchBestEnviron()");

  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;
  const vector<uint8> & PAF_multicopies = *PAF_multicopies_ptr;
  const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;

  vector<int32> bestid;
  bestid.resize(ENUMEND,-1);
  bool foundid=false;


  CEBUG("sbe 1\n");

  {
    vector<uint32> bestweight;
    bestweight.resize(ENUMEND,0);

    uint32 ids=0;
    uint32 elements=0;

    for(int32 actid=0; actid < static_cast<int32>(PAF_used_ids.size()); actid++) {
      if(PAF_used_ids[actid]!=0) continue;
      if(PAF_istroublemaker[actid]!=0) continue;
      vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[actid];

      ids++;

      uint32 numids=0;
      vector<uint32> sumweight;
      sumweight.resize(ENUMEND,0);

      int32 actidtemplatepartnerid=PAF_readpool->getRead(actid).getTemplatePartnerID();

      for(;I!=PAF_forward_edges->end() && actid==I->rid1;I++){
	// don't count known troublemakers
	if(PAF_istroublemaker[I->linked_with]) continue;
	// if the linked read is already used, can't count that as
	//  plus for a new anchor point
	if(PAF_used_ids[I->linked_with]!=0) continue;

	numids++;
	// look if it's the single best overlap so far
	if(I->best_weight > sumweight[SINGLE]) {
	  sumweight[SINGLE]=I->best_weight;
	}

	// add weight to normal overlaps
	sumweight[NORMAL]+=I->best_weight;

	if(PAF_multicopies[actid]==0) {
	  sumweight[NORMAL_NONMULTYCOPY]+=I->best_weight;
	}

	// look if both reads have templatepartnerids
	if(actidtemplatepartnerid >= 0) {
	  int32 secondtemplatepartnerid=PAF_readpool->getRead(I->linked_with).getTemplatePartnerID();
	  if(secondtemplatepartnerid >= 0) {
	    sumweight[WITH_TEMPLATE]+=I->best_weight;

	    if(PAF_multicopies[actid]==0) {
	      sumweight[WITH_TEMPLATE_NONMULTYCOPY]+=I->best_weight;
	    }

	    // look if they're the same template
	    if(actid == secondtemplatepartnerid) {
	      sumweight[WITH_TEMPLATE_PARTNER]+=I->best_weight;

	      if(PAF_multicopies[actid]==0
		 && PAF_multicopies[actidtemplatepartnerid]==0){
		sumweight[WITH_TEMPLATE_PARTNER_NONMULTICOPIES]+=I->best_weight;
	      }
	    }
	  }
	}
      }
      for(uint32 i=0; i<ENUMEND; i++) {
	if(sumweight[i]>bestweight[i]){
	  foundid=true;
	  bestid[i]=actid;
	  bestweight[i]=sumweight[i];
	}
      }
      elements+=numids;
      CEBUG("ID: " << actid << " has " << numids << " elements. Summed weights: Single - " << sumweight[SINGLE] << "\tNormal - " << sumweight[NORMAL] << "\tTempl. - " << sumweight[WITH_TEMPLATE] << "\tTempl.+Partn. - " << sumweight[WITH_TEMPLATE_PARTNER] << "\n");
      //CEBUGF(" giving an averaged weight of " << (uint32) (sumweight/numids) << endl);
    }

    if(ids>0) {
      CEBUG("ID: Average elements per id: " << double(elements)/(double)ids << "\n");
    }
  }

  // Catch the reads which do not overlap with any other one to make single
  //  reads contigs out of them
  if(foundid==false){
    uint32 size=PAF_readpool->size();
    for(uint32 i=0; i<size; i++){
      if(PAF_used_ids[i]==0){
  	bestid[SINGLE]=i;
  	break;
      }
    }
  }

  int32 ret=-1;
  if(PAF_pafparams.paf_use_genomic_algorithms) {
    if(bestid[WITH_TEMPLATE_PARTNER_NONMULTICOPIES] >=0) {
      ret=bestid[WITH_TEMPLATE_PARTNER_NONMULTICOPIES];
      CEBUG("Using Best-ID (tempart nonMC/genomic): " << ret << "\n");
    } else if(bestid[WITH_TEMPLATE_NONMULTYCOPY] >=0) {
      ret=bestid[WITH_TEMPLATE_NONMULTYCOPY];
      CEBUG("Using Best-ID (temp nonMC/genomic): " << ret << "\n");
    } else if(bestid[NORMAL_NONMULTYCOPY] >=0) {
      ret=bestid[NORMAL_NONMULTYCOPY];
      CEBUG("Using Best-ID (normal nonMC/genomic): " << ret << "\n");
    } else if(bestid[NORMAL] >=0) {
      ret=bestid[NORMAL];
      CEBUG("Using Best-ID (normal/genomic): " << ret << "\n");
    } else if(bestid[SINGLE] >=0) {
      ret=bestid[SINGLE];
      CEBUG("Using Best-ID (single/genomic): " << ret << "\n");
    } else {
      CEBUG("No Best-ID genomic found.\n");
    }
  } else if(bestid[WITH_TEMPLATE_PARTNER] >=0) {
    ret=bestid[WITH_TEMPLATE_PARTNER];
    CEBUG("Using Best-ID with template partner: " << ret << "\n");
  } else if(bestid[WITH_TEMPLATE] >=0) {
    ret=bestid[WITH_TEMPLATE];
    CEBUG("Using Best-ID with template: " << ret << "\n");
  } else if(bestid[NORMAL] >=0) {
    ret=bestid[NORMAL];
    CEBUG("Using Best-ID (normal): " << ret << "\n");
  } else if(bestid[SINGLE] >=0) {
    ret=bestid[SINGLE];
    CEBUG("Using Best-ID (single): " << ret << "\n");
  } else {
    CEBUG("No Best-ID found.\n");
  }

  CEBUG("Returning: " << ret << '\n');

  FUNCEND();

  return ret;

  //  return 0;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Pathfinder::checkindex(int32 index)
{
  FUNCSTART("void checkindex(int32 index)");

  if(index <0 || index>=static_cast<int32>(PAF_used_ids_ptr->size())){
    cout << "Index given: " << index << "\tSize: " << PAF_used_ids_ptr->size() << "\n";
    MIRANOTIFY(Notify::INTERNAL, "Index overrun.");
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define FUNCTRACE(bla) { cout << THISFUNC << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Pathfinder::buildContig(vector<Align> & aligncache, int32 startid, Contig & con)
{
  FUNCSTART("void Pathfinder::buildContig(uint32 startid, Contig & con)");

  CEBUG("bc 1" << endl);

  uint32 lasttimingcout=0;

#ifndef PUBLICQUIET
  cout << "mayaddmcmcoverlaps: " << PAF_mayaddmcmcoverlaps << endl;
  cout << "musthonourmcmcflag: " << PAF_musthonourmcmcflag << endl;
#endif

  PAF_railoverlapcache.clear();

  // if we´re starting with an existing contig (e.g. when relooping
  //  or with backbones, the make sure to search whole contig first time)
  if(startid<0) {
#ifndef PUBLICQUIET
    cout << "startid<0.\n";
#endif
    PAF_skipwholecontigscan_counter=1;
  }else{
#ifndef PUBLICQUIET
    cout << "startid: " << startid << "\n";
#endif
    PAF_skipwholecontigscan_counter=PAF_pafparams.paf_skipwholecontigscan;
  }

  vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;

  Contig::errorstatus_t contigerrstat;
  // this should speed up things a bit, avoiding copying during extension
  //  although the 10000 might be not enough for really pathological cases
  //  with many SRMB regions and huge coverage
  contigerrstat.reads_affected.reserve(10000);

  // this is a kind of duplicate of PAF_used_ids_in_this_run, but I thought
  //  that it could give faster access as only the really used ids are in there

  //vector<int32> ids_in_contig;
  vector<int32> & ids_in_contig=*PAF_ids_in_contig_ptr;
  ids_in_contig.clear();
  ids_in_contig.reserve(PAF_used_ids.size()+10);

  PAF_buildcontig_newlinecounter=0;


  CEBUG("bc 2" << endl);

  bool isbackboneassembly=false;
  // do we start with a fresh contig?
  if(con.getContigLength()==0) {
#ifndef PUBLICQUIET
    cout << "empty contig\n";
#endif
    CEBUG("bc 3a 1" << endl);
    if(startid<0) {
      cerr << "Startid: " << startid << endl;
      MIRANOTIFY(Notify::INTERNAL, "Starting with non-existing read?");
    }

    if(startid >= static_cast<int32>(PAF_used_ids.size())) {
      cerr << "Startid: " << startid << endl;
      MIRANOTIFY(Notify::INTERNAL, "Starting with read id greater than number iof reads?\n");
    }

    if(PAF_used_ids[startid]) {
      cerr << "Startid: " << startid << endl;
      MIRANOTIFY(Notify::INTERNAL, "Starting with read that is already used?\n");
    }

    CEBUG("Startid: " << startid << endl);

#ifndef PUBLICQUIET
    cout << endl;
#endif

    ++PAF_readaddattempts;
    con.addRead(aligncache,
		nullptr, startid, startid, 1,
		(*PAF_multicopies_ptr)[startid],
		0,
		contigerrstat);

    CEBUG("bc 3a 2" << endl);

    PAF_used_ids[startid]=1;
    PAF_used_ids_in_this_run[startid]=1;
    ids_in_contig.push_back(startid);

#ifndef PUBLICQUIET
    cout << endl;
#else
    cout << "+";
#endif
  }else{
    // we start with an already made contig,
    // initialise needed values

#ifndef PUBLICQUIET
    cout << "existing contig\n";
#endif

    CEBUG("bc 3b" << endl);

    PAF_railsincontig.clear();
    auto & cr=con.getContigReads();

    for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
      if(pcrI.getORPID()>=0
	 && (pcrI->isRail()
	     || pcrI->isBackbone())){
	PAF_used_ids_in_this_run[pcrI.getORPID()]=1;
	ids_in_contig.push_back(pcrI.getORPID());
	//cout << "In contig:" << PAF_readpool->getRead(pcrI->id).getName() << endl;
	if(pcrI->isRail()) PAF_railsincontig.push_back(pcrI.getORPID());
      }
    }

    CEBUG("bc 3b 2" << endl);

//#ifndef PUBLICQUIET
//    Read::setCoutType(Read::AS_TEXTCLIPS);
//    cout << "Front read:\n" << cr.front().read << endl;
//#endif

    // check whether this this is a contig with a backbone (and therefore
    //  also rails). If yes, we'll use bb assembly later
    if(con.getNumBackbones()>0) {
      isbackboneassembly=true;
#ifndef PUBLICQUIET
      cout << "Backbone assembly to " << con.getContigName() << endl;
#endif

      prepareRailOverlapCache(con.getAllowedRefIDs());

    }
  }

  CEBUG("bc 4" << endl);

  struct tms mytms;
  times(&mytms);
  clock_t baseclocks=mytms.tms_utime+mytms.tms_stime;
  clock_t actclocks=baseclocks;
  clock_t maxallowedclocks=baseclocks+PAF_pafparams.paf_maxcontigclockticks;


#ifdef CLOCK_STEPS1
  vector<suseconds_t> us_pfsearch;
  vector<suseconds_t> us_conadd;
  vector<suseconds_t> us_overlapban;
#endif

  nextreadtoadd_t nrta;
  nrta.newid=0;

  uint32 statscounter=1;

  // findNextBackboneOverlapQuick() uses these
  bool allowbbqmulticopies=false;
  bool allowbbqtroublemakers=false;
  bool allowbbqsmallhits=false;

  while(nrta.newid >= 0) {
    cout.flush();

    PAF_skipwholecontigscan_counter--;

    if(PAF_skipwholecontigscan_counter<0) PAF_skipwholecontigscan_counter= PAF_pafparams.paf_skipwholecontigscan;
;

#ifndef PUBLICQUIET
    if(statscounter % 25 == 0) {
      cout << "rstatsm: len of c : " << con.getContigLength();
      cout << "\nrstatsm: ids in c : " << ids_in_contig.size();
      cout << "\nrstatsi: banned o : " << PAF_overlapsbanned;
#ifdef CLOCK_STEPS1
      cout << "\nrtimings: " << ids_in_contig.size();
      cout << " : " << avg_suseconds(us_pfsearch);
      cout << " / " << avg_suseconds(us_conadd);
      cout << " / " << avg_suseconds(us_overlapban);
#endif
      cout << '\n';
    }
    statscounter++;
#endif


    nrta.refid=-1;
    nrta.newid=-1;
    nrta.weight=0;
    nrta.direction_newid=0;
    nrta.ads_node=nullptr;

#ifdef CLOCK_STEPS1
    timeval us_start;
    gettimeofday(&us_start,nullptr);
#endif

#ifndef PUBLICQUIET
    cout << "Searching";
#endif
    if(nrta.newid < 0) {
      if(isbackboneassembly){
#ifndef PUBLICQUIET
	cout << " (backbone quick) ...";
	cout.flush();
#endif
	findNextBackboneOverlapQuick(nrta,
				     con.getAllowedRefIDs(),
				     allowbbqmulticopies,
				     allowbbqtroublemakers,
				     allowbbqsmallhits);

//	if(nrta.newid < 0) {
//#ifndef PUBLICQUIET
//	  cout << " (backbone normal) ...";
//	  cout.flush();
//#endif
//	  findNextBackboneOverlapNormal(ids_in_contig, nrta, con.getAllowedRefIDs());
//	}

	if(nrta.newid < 0) {
#ifndef PUBLICQUIET
	  cout << "\nRebuilding overlap cache:";
	  cout.flush();
#endif
	  prepareRailOverlapCache(con.getAllowedRefIDs());
#ifndef PUBLICQUIET
	  cout << " (backbone quick 2) ...";
	  cout.flush();
#endif
	  findNextBackboneOverlapQuick(nrta,
				       con.getAllowedRefIDs(),
				       allowbbqmulticopies,
				       allowbbqtroublemakers,
				       allowbbqsmallhits);


	  if(nrta.newid < 0) {
#ifndef PUBLICQUIET
	    cout << "\nNow allowing small hits:";
	    cout.flush();
#endif
	    prepareRailOverlapCache(con.getAllowedRefIDs());
#ifndef PUBLICQUIET
	    cout << " (backbone quick 3) ...";
	    cout.flush();
#endif
	    allowbbqmulticopies=true;
	    allowbbqtroublemakers=true;
	    allowbbqsmallhits=true;
	    findNextBackboneOverlapQuick(nrta,
					 con.getAllowedRefIDs(),
					 allowbbqmulticopies,
					 allowbbqtroublemakers,
					 allowbbqsmallhits);

	    if(nrta.newid < 0) {
#ifndef PUBLICQUIET
	      cout << "\nNow allowing multicopies and troublemakers:";
	      cout.flush();
#endif
	      prepareRailOverlapCache(con.getAllowedRefIDs());
#ifndef PUBLICQUIET
	      cout << " (backbone quick 4) ...";
	      cout.flush();
#endif
	    allowbbqmulticopies=true;
	    allowbbqtroublemakers=true;
	    findNextBackboneOverlapQuick(nrta,
					 con.getAllowedRefIDs(),
					 allowbbqmulticopies,
					 allowbbqtroublemakers,
					 allowbbqsmallhits);
	    }
	  }
	}



//	if(nrta.newid < 0) {
//#ifndef PUBLICQUIET
//	  cout << " (full next bb overlap) ...";
//	  cout.flush();
//#endif
//	  findNextOverlap(ids_in_contig, nrta);
//	}
      }else{
#ifndef PUBLICQUIET
	cout << " (full next overlap) ...";
	cout.flush();
#endif
	findNextOverlap(ids_in_contig, nrta);
      }
    }

#ifndef PUBLICQUIET
      cout << " done.\n";
      cout.flush();
#endif

#ifdef CLOCK_STEPS1
      us_pfsearch.push_back(diffsuseconds(us_start));
#endif


    if(nrta.newid >= 0) {
#ifndef PUBLICQUIET
      cout << nrta;
#endif

#ifdef CLOCK_STEPS1
      gettimeofday(&us_start,nullptr);
#endif
      ++PAF_readaddattempts;
      con.addRead(aligncache,
		  nrta.ads_node, nrta.refid, nrta.newid, nrta.direction_newid,
		  (*PAF_multicopies_ptr)[nrta.newid],
		  0,
		  contigerrstat);
#ifdef CLOCK_STEPS1
      us_conadd.push_back(diffsuseconds(us_start));
#endif

      cout.flush();

      //Contig::setCoutType(Contig::AS_TEXT);
      //cout << ",,,,,,\n" << con << "\n\n";


      if(PAF_buildcontig_newlinecounter==0){
	cout << "[" << ids_in_contig.size() << "] "; cout.flush();
      }
      if(contigerrstat.code == Contig::ENOERROR) {
	ids_in_contig.push_back(nrta.newid);
	PAF_used_ids[nrta.newid]=1;
	PAF_used_ids_in_this_run[nrta.newid]=1;

#ifndef PUBLICQUIET
	cout << "\t+" << endl;
#else
	cout << "+";
#endif

	FUNCTRACE("bC bla1\n");

	if((*PAF_multicopies_ptr)[nrta.newid] == 0) {
	  PAF_nonmulticopiesadded++;
	}

	FUNCTRACE("bC bla3\n");

	// if we added >=10 non-multicopies since the last reject
	//  de-blacklist the reads that the added read links to
	// TODO: check whether strategy can be improved with
	//  "minimum number of new bases in contig left or right"
	// TODO: make configurable?
	if(PAF_nonmulticopiesadded >= 10) {
	  vector<newedges_t>::iterator neI = (*PAF_lowerbound_oedges_ptr)[nrta.newid];
	  for(; neI != PAF_forward_edges->end() && neI->rid1 == nrta.newid; neI++){
	    if(PAF_blacklisted_ids[nrta.newid]==0) {
	      PAF_blacklisted_ids[nrta.newid]=1;
	    }
	  }
	}

	FUNCTRACE("bC bla4\n");

      } else {
	// countdown the blacklist counter (if blacklist is used)
	//  EXCEPT:
	//  - when the reference ID is a backbone rail
	//  - when the template partner is not a multicopy and in the contig
	if(PAF_pafparams.paf_use_emergency_blacklist){
	  bool doblacklist=true;
	  if(PAF_readpool->getRead(nrta.refid).isRail()){
	    doblacklist=false;
	  }else if(PAF_readpool->getRead(nrta.refid).getTemplatePartnerID()>=0){
	    if(! (*PAF_multicopies_ptr)[PAF_readpool->getRead(nrta.refid).getTemplatePartnerID()]
	       && PAF_used_ids_in_this_run[PAF_readpool->getRead(nrta.refid).getTemplatePartnerID()]){
	      doblacklist=false;
	    }
	  }

	  if(doblacklist && PAF_blacklisted_ids[nrta.newid]){
	    PAF_blacklisted_ids[nrta.newid]-=1;
	  }
	}

	// reset the non-multycopiy counter so that de-blacklisting
	//  can be performed once the counter will reach a threshold
	PAF_nonmulticopiesadded=0;

#ifdef PUBLICQUIET
	contigerrstat.dumpStatus();
#else
	contigerrstat.dumpStatus(true,"bCds ");
#endif

	{
	  // check if refid has been put into reads_affected by the contig
	  // if not, warn and put it in
	  bool isin=false;
	  for(auto raid : contigerrstat.reads_affected){
	    if(raid == nrta.refid){
	      isin=true;
	      break;
	    }
	  }
	  if(!isin){
	    cout << "\nWARNING: contig did not put refid into affected reads. Please file a bug report to\n\thttps://sourceforge.net/apps/trac/mira-assembler/\nand please send a mail to\n\tmira_talk@freelists.org\n";
	    cout << "current af size: " << contigerrstat.reads_affected.size() << endl;
	    contigerrstat.reads_affected.push_back(nrta.refid);
	    exit(0);
	  }
	}

	// ban all alignments of that read with other reads that were
	//  at the same position

	//// debug: show all reads at that pos that the contig gave back
	//{
	//  cout << "\nReads affected:\n";
	//  vector<int32>::const_iterator raI;
	//  raI=contigerrstat.reads_affected.begin();
	//  for(; raI != contigerrstat.reads_affected.end(); raI++){
	//    cout << *raI << '\t' << (*PAF_readpool)[*raI].getName();
	//    if((*PAF_readpool)[*raI].isBackbone()){
	//      cout << "\tbackbone\n";
	//    } else if((*PAF_readpool)[*raI].isRail()){
	//      cout << "\trail\n";
	//    }else{
	//      cout << "\tseqtype: " << static_cast<uint16>((*PAF_readpool)[*raI].getSequencingType()) << '\n';
	//    }
	//  }
	//}

#ifdef CLOCK_STEPS1
	gettimeofday(&us_start,nullptr);
#endif
	{
	  // First search look at the newedges_t block of nrta,newid overlaps
	  //  and ban those overlaps whose reads are in the vector given back
	  //  by the contig

	  vector<newedges_t>::iterator neI = (*PAF_lowerbound_oedges_ptr)[nrta.newid];
	  vector<int32>::const_iterator raI;
	  for(; neI != PAF_forward_edges->end() && neI->rid1 == nrta.newid; neI++){
	    if(!neI->pf_banned){
	      raI=contigerrstat.reads_affected.begin();
	      for(; raI != contigerrstat.reads_affected.end(); raI++){
		if(*raI == neI->linked_with){
		  //cout << "Banning1: " << neI->rid1 << '\t' << neI->linked_with << '\n';
		  PAF_overlapsbanned++;
		  neI->pf_banned=true;
		  if(PAF_overlapsbanned_smallstore.size()<PAF_overlapsbanned_smallstore.capacity()){
		    PAF_overlapsbanned_smallstore.push_back(neI);
		  }
		  break;
		}
	      }
	    }
	  }

	  // Now look at the blocks of newedges_t whose reads are given
	  //  by by the vector and who overlap with nrta.newid
	  raI=contigerrstat.reads_affected.begin();
	  for(; raI != contigerrstat.reads_affected.end(); raI++){
	    neI = (*PAF_lowerbound_oedges_ptr)[*raI];
	    for(; neI != PAF_forward_edges->end() && neI->rid1 == *raI; neI++){
	      if(!neI->pf_banned){
		if(neI->linked_with == nrta.newid){
		  //cout << "Banning2: " << neI->rid1 << '\t' << neI->linked_with << '\n';
		  PAF_overlapsbanned++;
		  neI->pf_banned=true;
		  if(PAF_overlapsbanned_smallstore.size()<PAF_overlapsbanned_smallstore.capacity()){
		    PAF_overlapsbanned_smallstore.push_back(neI);
		  }
		  break;
		  // NOTE: should MIRA ever get multiple overlaps between
		  //  two reads, then the "break;" above must be removed
		  //  should be safe for now.
		}
	      }
	    }
	  }
	}

#ifdef CLOCK_STEPS1
	us_overlapban.push_back(diffsuseconds(us_start));
#endif
      }

#ifndef PUBLICQUIET
#else
      PAF_buildcontig_newlinecounter++;
      if(PAF_buildcontig_newlinecounter==60){
	PAF_buildcontig_newlinecounter=0;
	cout << "   " << con.getContigLength();
#ifdef CLOCK_STEPS1
	cout << "\tpft\t" << avg_suseconds(us_pfsearch);
	cout << " / " << avg_suseconds(us_conadd);
	cout << " / " << avg_suseconds(us_overlapban);
#endif
	cout << endl;
	if(con.getNumReadsInContig() / 10000 > lasttimingcout){
	  con.coutAddReadTimings();
	  lasttimingcout=con.getNumReadsInContig()/10000;
	}
      }
#endif


    }

    //if(con.getNumReadsInContig()==1198){
    //  cout << "\npathfinder debugstop.\n";
    //  assout::saveAsCAF(con, "test.caf", true);
    //  abort();
    //}

    if(PAF_pafparams.paf_use_max_contig_buildtime){
      times(&mytms);
      actclocks=mytms.tms_utime+mytms.tms_stime;
      if(actclocks>maxallowedclocks){
	// the following should not be necessary, but just to be sure
	// (had a case where the disk got full and MIRA was stuck in
	//  endless loop ... calling times()
	// so I set the loop ending condition manually
	nrta.newid=-1;

	cout << "\nMaximum build time for this contig reached, aborting build.\n";
	PAF_railoverlapcache.clear();
	FUNCEND();
	return;
      }
    }
  }

  PAF_railoverlapcache.clear();

  FUNCEND();
}
//#define FUNCTRACE(bla)

//#define CEBUG(bla)






/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Pathfinder::prepareRailOverlapCache(const vector<bool> & allowedrefids)
{
  vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;

#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  // initialise PAF_railoverlapcache
  // each read that has an overlap with a rail is put into the
  //  vector once;
  PAF_railoverlapcache.clear();
  PAF_railoverlapcache.reserve(PAF_used_ids.size());

#ifdef CLOCK_STEPS2
  cout << "Timing pROC PAF_railoverlapcache: " << diffsuseconds(tv) << endl;
  gettimeofday(&tv,nullptr);
#endif

  vector<uint8> tmp_alreadyinvec(PAF_used_ids.size(),0);

#ifdef CLOCK_STEPS2
  cout << "Timing pROC tmp_alreadyinvec: " << diffsuseconds(tv) << endl;
  gettimeofday(&tv,nullptr);
#endif

  // put first non-sanger types in vector, then sanger types
  // as vector will be worked backwards, sanger will get first

  /* 139000 ticks (0.13s) for short contig (<10kb) & 8m reads
     242197 ticks for long contig (2.3m) & 8m reads
     333269 for 4.3m & 8m

  vector<newedges_t>::const_iterator rcI=PAF_forward_edges->begin();
  for(;rcI!=PAF_forward_edges->end(); rcI++){
    CEBUG('\n' << *rcI);
    CEBUG('1');

    if(PAF_used_ids[rcI->rid1]) continue;
    if(!PAF_used_ids_in_this_run[rcI->linked_with]) continue;

    if(!PAF_readpool->getRead(rcI->linked_with).isRail()) continue;
    if(PAF_readpool->getRead(rcI->rid1).isRail()) continue;
    if(tmp_alreadyinvec[rcI->rid1]) continue;
    if(rcI->pf_banned) continue;
    if(!allowedrefids.empty() && !allowedrefids[rcI->linked_with]) continue;
    // having this last means: if all SANGER mapped, hasSANGER will still be "false"
    if(PAF_readpool->getRead(rcI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_SANGER)) {
      hasSANGER=true;
      continue;
    }

    tmp_alreadyinvec[rcI->rid1]=1;
    PAF_railoverlapcache.push_back(rcI->rid1);
  }
  */

  // different strategy:
  // using the rails of this contig, go through their matches
  //
  //   ~7000 ticks for short contig (<10kb) & 8m reads (~15 to 20x faster than above)
  //   101570 ticks for long contig (2.3m) & 8m reads (~2.5x faster than above) in first map,
  //          then <20000 decreasing to 10000 ticks (~15x to 25x faster than above)
  //   333269 for 4.3m & 8m

  bool hasSANGER=false;
  bool has454=false;

  for(auto rid : PAF_railsincontig){
    auto rcI=(*PAF_lowerbound_oedges_ptr)[rid];
    if(!allowedrefids.empty() && !allowedrefids[rid]) continue;
    for(; rcI != PAF_forward_edges->end() && rcI->rid1==rid; ++rcI){
      if(PAF_used_ids[rcI->linked_with]) continue;
      if(tmp_alreadyinvec[rcI->linked_with]) continue;
      if(PAF_readpool->getRead(rcI->linked_with).isRail()) continue;
      if(rcI->pf_banned) continue;
      // having this last means: if all SANGER reads mapped, hasSANGER will still be "false"
      if(PAF_readpool->getRead(rcI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_SANGER)) {
	hasSANGER=true;
	continue;
      }
      if(PAF_readpool->getRead(rcI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_454GS20)) {
	has454=true;
	continue;
      }
      tmp_alreadyinvec[rcI->linked_with]=1;
      PAF_railoverlapcache.push_back(rcI->linked_with);
    }
  }

#ifdef CLOCK_STEPS2
  cout << "Timing pROC fill non-SANGER: " << diffsuseconds(tv) << endl;
  gettimeofday(&tv,nullptr);
#endif

#ifndef PUBLICQUIET
  cout << "Backbone overlap cache (non SANGER): " << PAF_railoverlapcache.size() << endl;
#endif

  if(has454){
    for(auto rid : PAF_railsincontig){
      auto rcI=(*PAF_lowerbound_oedges_ptr)[rid];
      if(!allowedrefids.empty() && !allowedrefids[rid]) continue;
      for(; rcI != PAF_forward_edges->end() && rcI->rid1==rid; ++rcI){
	// checking this first means that all other checks will not run in most cases in assemblies
	//  with majority of non-SANGER
	//
	// "!" is only diff to loop above (except placement at top of loop)
	if( ! PAF_readpool->getRead(rcI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_454GS20)) continue;
	//

	if(PAF_used_ids[rcI->linked_with]) continue;
	if(tmp_alreadyinvec[rcI->linked_with]) continue;
	if(PAF_readpool->getRead(rcI->linked_with).isRail()) continue;
	if(rcI->pf_banned) continue;
	tmp_alreadyinvec[rcI->linked_with]=1;
	PAF_railoverlapcache.push_back(rcI->linked_with);
      }
    }
  }

#ifdef CLOCK_STEPS2
  cout << "Timing pROC fill 454: " << diffsuseconds(tv) << endl;
  cout << "Timing pROC total: " << diffsuseconds(tvtotal) << endl;
#endif

  if(hasSANGER){
    for(auto rid : PAF_railsincontig){
      auto rcI=(*PAF_lowerbound_oedges_ptr)[rid];
      if(!allowedrefids.empty() && !allowedrefids[rid]) continue;
      for(; rcI != PAF_forward_edges->end() && rcI->rid1==rid; ++rcI){
	// checking this first means that all other checks will not run in most cases in assemblies
	//  with majority of non-SANGER
	//
	// "!" is only diff to loop above (except placement at top of loop)
	if( ! PAF_readpool->getRead(rcI->rid1).isSequencingType(ReadGroupLib::SEQTYPE_SANGER)) continue;
	//

	if(PAF_used_ids[rcI->linked_with]) continue;
	if(tmp_alreadyinvec[rcI->linked_with]) continue;
	if(PAF_readpool->getRead(rcI->linked_with).isRail()) continue;
	if(rcI->pf_banned) continue;
	tmp_alreadyinvec[rcI->linked_with]=1;
	PAF_railoverlapcache.push_back(rcI->linked_with);
      }
    }
  }

#ifdef CLOCK_STEPS2
  cout << "Timing pROC fill SANGER: " << diffsuseconds(tv) << endl;
  cout << "Timing pROC total: " << diffsuseconds(tvtotal) << endl;
#endif

#ifndef PUBLICQUIET
  cout << "Backbone overlap cache (non SANGER & SANGER): " << PAF_railoverlapcache.size() << endl;
#endif
}

//#define CEBUG(bla)


// Search for overlaps of reads that are in the contig
// expects resultread to contain illegal/empty reference
// version for backbone assembly

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Pathfinder::findNextBackboneOverlapNormal(const vector<int32> & readsincontig, nextreadtoadd_t & resultread, const vector<bool> & allowedrefids, bool allowmulticopies, bool allowtroublemakers)
{
  FUNCSTART("void Pathfinder::findNextBackboneOverlapNormal(vector<int32> & readsincontig, nextreadtoadd_t & resultread, const vector<bool> & allowedrefids, bool allowmulticopies, bool allowtroublemakers)");

  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;
  const vector<uint8> & PAF_multicopies = *PAF_multicopies_ptr;
  const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;

  BUGIFTHROW(readsincontig.size()==0, "No reads present in contig?");

  vector<int32>::const_iterator rcI=readsincontig.begin();
  //for(uint32 index=0; index <readsincontig.size() ; index++, rcI++) {
  //  int32 readid=readsincontig[index];
  for(; rcI != readsincontig.end() ; rcI++) {
    //// read in contig is blacklisted (==no more un-assembled partners)?
    //// -> next loop iteration
    // are backbones blacklisted???
    //if(PAF_blacklisted_ids[*rcI] == 0) continue;

    CEBUG("Check: " <<  PAF_readpool->getRead(*rcI).getName() << "\tRail? " << PAF_readpool->getRead(*rcI).isRail() << '\n');

    // read in contig is not rail? next loop iteration
    if(!PAF_readpool->getRead(*rcI).isRail()) continue;

    // contig has list of allowed refids and read is not allowed? next loop
    if(!allowedrefids.empty() && !allowedrefids[*rcI]) continue;

    vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[*rcI];

    for(;I!=PAF_forward_edges->end() && I->rid1 == *rcI;I++){
      // the edges are sorted by rid1, then by "bestweight: high to low"
      //  therefore, if the current bestweight is <= the one we have found already,
      //  we do not need to check further edges of this read ... all remaining
      //  edges have lower or equal quality anyway. Therefore, break out of loop
      //  here if possible
      if(I->best_weight <= resultread.weight) break;

      // don't bother looking at if overlap is banned
      if(I->pf_banned) continue;

      //int32 linkedwithid=I->linked_with;

      CEBUG("  with: " << PAF_readpool->getRead(I->linked_with).getName());

      // don't bother looking at if read linked to is already used
      // or read is blacklisted
      if(PAF_used_ids[I->linked_with] !=0) {
	CEBUG(" already used!\n");
	continue;
      }
      if(PAF_blacklisted_ids[I->linked_with]==0) {
	CEBUG(" blacklisted!\n");
	continue;
      }

      // if necessary, take into account only reads that are
      //  - non-multicopies
      //  - non-troublemakers
      if((allowmulticopies || PAF_multicopies[I->linked_with]==0)
	&& (allowtroublemakers || PAF_istroublemaker[I->linked_with] == 0)){

	CEBUG(" non-MC/TM. Ibw:" << I->best_weight << "\trrw: " << resultread.weight << "\toll: " <<(*PAF_adsfacts)[I->adsfindex].getOverlapLen());

	////  - that have overlap length >= min length for quick bb overlap (just to have a good match)
	//if(I->best_weight > resultread.weight
	//   && (*PAF_adsfacts)[I->adsfindex].getOverlapLen() >= (*PAF_miraparams)[PAF_readpool->getRead(I->linked_with).getSequencingType()].getPathfinderParams().paf_bbquickoverlap_minlen) {
	if(I->best_weight > resultread.weight){
	  CEBUG(" take it.");
	  resultread.refid=*rcI;
	  resultread.newid=I->linked_with;
	  resultread.weight=I->best_weight;
	  resultread.direction_newid=I->direction;
	  resultread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	  resultread.foundmethod=FOUND_RAIL_AND_NONMULTICOPY;
	}else{
	  CEBUG(" don't take.");
	}
	CEBUG('\n');
      }else{
	CEBUG(" MC or TM.\n");
      }
    }

  }
}

//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// reverse logic compared to findNextBackboneOverlapNormal()
// Take one non-rail, find best place for a it in existing contig
// Disregard multicopies and troublemakers

// expects resultread to contain illegal/empty reference
void Pathfinder::findNextBackboneOverlapQuick(nextreadtoadd_t & resultread, const vector<bool> & allowedrefids, bool allowmulticopies, bool allowtroublemakers, bool allowsmallhits)
{
  FUNCSTART("void Pathfinder::findNextBackboneOverlap(nextreadtoadd_t & resultread, const vector<bool> & allowedrefids, bool allowmulticopies, bool allowtroublemakers)");

#ifndef PUBLICQUIET
  cout << "roc: " << PAF_railoverlapcache.size() << ' ';
#endif

  if(PAF_railoverlapcache.empty()){
    CEBUG("\troc empty\n");
    // found nothing in previous loop, exit
    return;
  }

  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;
  const vector<uint8> & PAF_multicopies = *PAF_multicopies_ptr;
  const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;

  while(resultread.newid<0){
    bool continuesearch=true;
    int32 readid=-1;
    while(continuesearch && !PAF_railoverlapcache.empty()){
      readid=PAF_railoverlapcache.back();

      CEBUG("\nfnboq:\n");
      CEBUG("l: " << PAF_readpool->getRead(readid).getName());
      CEBUG("\tused: " << (int16) PAF_used_ids[readid]);
      CEBUG("\tmc: " << (int16) PAF_multicopies[readid]);
      CEBUG("\ttm: " << (int16) PAF_istroublemaker[readid]);
      if(!allowedrefids.empty()){
	CEBUG("\tar: " << (int16) allowedrefids[readid]);
      }

      PAF_railoverlapcache.pop_back();

      if(PAF_used_ids[readid]==0
	 && (allowmulticopies || PAF_multicopies[readid]==0)
	 && (allowtroublemakers || PAF_istroublemaker[readid]==0)){
	continuesearch=false;
//#ifndef PUBLICQUIET
//	cout << "stopsearch ";
//#endif
      }
//#ifndef PUBLICQUIET
//      else {
//	cout << "(* r *)";
//      }
//#endif
    }


    if(continuesearch){
      CEBUG("\tnot taken/found\n");
      // found nothing in previous loop, exit
      return;
    }

    CEBUG("\ttaken/found\n");

    // ok, found one. Now search the rail it fits best to *in this contig*
    vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[readid];

    for(;I!=PAF_forward_edges->end() && I->rid1 == readid;I++){
      // don't bother looking at if overlap is banned
      if(I->pf_banned) continue;
      // must link to rail
      if(!PAF_readpool->getRead(I->linked_with).isRail()) continue;
      // rail must be in this contig!
      if(!PAF_used_ids_in_this_run[I->linked_with]) continue;
      // rail must be allowed as refid
      if((!allowedrefids.empty() && !allowedrefids[I->linked_with])) continue;

      // if necessary, take into account only rails that are
      //  - non-multicopies
      //  - non-troublemakers
      if((allowmulticopies || PAF_multicopies[I->linked_with]==0)
	 && (allowtroublemakers || PAF_istroublemaker[I->linked_with] == 0)){
	//  - that have overlap length >= minim length (just to have
	//    good matches first)
	if(I->best_weight > resultread.weight){
	  if(allowsmallhits
	     || (*PAF_adsfacts)[I->adsfindex].getOverlapLen() >=
	     (*PAF_miraparams)[PAF_readpool->getRead(readid).getSequencingType()].getPathfinderParams().paf_bbquickoverlap_minlen) {
	    resultread.refid=I->linked_with;
	    resultread.newid=readid;
	    resultread.weight=I->best_weight;
	    resultread.direction_newid=I->direction;
	    resultread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    if(PAF_multicopies[readid]==0){
	      resultread.foundmethod=FOUND_RAIL_AND_NONMULTICOPY;
	    }else{
	      resultread.foundmethod=FOUND_RAIL_AND_MULTICOPY;
	    }
	    // the edges are sorted by rid1, then by "bestweight: high to low"
	    // therefore, if we wound something, we do not need to check further
	    // edges of this read ... all remaining edges have lower or equal
	    // quality anyway. Therefore, break out of loop here
	    break;
	  }
	}
      }
    }
  }
}

//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Search for overlaps of reads that are in the contig
// expects resultread to contain illegal/empty reference

//#define CEBUGF(bla)   {cout << bla; cout.flush();}

void Pathfinder::findNextOverlap(const vector<int32> & readsincontig, nextreadtoadd_t & resultread)
{
  FUNCSTART("void Pathfinder::findNextOverlap(vector<int32> & readsincontig, nextreadtoadd_t & resultread)");

  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;
  const vector<uint8> & PAF_multicopies = *PAF_multicopies_ptr;
  const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;

  BUGIFTHROW(readsincontig.size()==0, "No reads present in contig?");

  vector<nextreadtoadd_t> bestreads(FOUND_ENUMEND,resultread);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  //  partner not, then the partner must be in use (not necessarily
  //  the same contig)
  bool usestricttemplatepairrule=true;



  bool   usequickruleSAN=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrule454=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_use_quick_rule;
  bool   usequickruleION=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrulePBSHQ=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_PACBIOHQ].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrulePBSLQ=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_PACBIOLQ].getPathfinderParams().paf_use_quick_rule;
  bool   usequickruleTXT=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_TEXT].getPathfinderParams().paf_use_quick_rule;

  bool   usequickruleSHORTREADS=true;

  bool speedymap=false;

  // counter for how many "insanely high coverage" reads we found
  //uint32 indigestible=0;

  struct tms mytms;
  times(&mytms);
  clock_t baseclocks=mytms.tms_utime+mytms.tms_stime;
  clock_t actclocks=baseclocks;
  clock_t maxallowedclocks=baseclocks+PAF_pafparams.paf_nextread_maxcttoess;

  // going backwards through the contglist should be a tick faster than forward
  //  also helps to retain a bit of flexibility when emergency search stops occur
  //for(uint32 index=0; index < readsincontig.size() ; index++) {

  // like variable "notfoundyet" below, this speeds up search but
  //  makes order of checks important
  bool needsearchnormal=true;
  bool foundsolution=false;


  // make at most two loops:
  //  in first try find next read via shortcut (looking at last 50 reads
  //  added to contig)
  //  if not found, search all reads of contig for next read to add
  // start directly with second loop
  // TODO: this is ... less than beautiful.
  int8 searchloops=0;
  if(PAF_skipwholecontigscan_counter>1) searchloops=1;
  for(; searchloops>=0 && !foundsolution; searchloops--){

    // TODO: flexibilisieren? hits < länge oder relscore -> reloop?

    int32 dontlookatreads=0;

    int32 index=readsincontig.size()-1;

    // TODO: check whether 15 is reasonable. Make configurable.
    // (50 works well, but too slow for assemblies with 454 data)
    if(searchloops>0){
      if(PAF_skipwholecontigscan_counter>1) {
	dontlookatreads=readsincontig.size()-1-15;
	if(dontlookatreads<0) dontlookatreads=0;
      }
    }else{
      if(PAF_skipwholecontigscan_counter>1){
	index=readsincontig.size()-1-15-1;
	if(index<0) index=0;
      }
    }

    uint32 numreadschecked=0;
    for(; index >=dontlookatreads ; index--, numreadschecked++) {
      int32 readid=readsincontig[index];

      CEBUGF("fnO: new read " << readid << '\n');

      // read in contig is blacklisted (==no more un-assembled partners)?
      // -> next loop iteration
      if(PAF_blacklisted_ids[readid] == 0) continue;

      vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[readid];

      uint32 numpotentialpartners=0;

      for(;I!=PAF_forward_edges->end() && I->rid1 == readid;I++){

	CEBUGF("fNO: " <<  I->rid1 << '\t' << I->linked_with << '\n');

	// don't bother looking at if overlap is banned
	if(I->pf_banned) continue;

	CEBUGF("#1\n");
	int32 linkedwithid=I->linked_with;

	CEBUGF("#2\n");
	// don't bother looking at if read linked to is already used
	if(PAF_used_ids[linkedwithid] !=0) continue;

	CEBUGF("#3\n");

	// rails cannot be taken as linkedwithid
	if((*PAF_readpool)[linkedwithid].isRail()) continue;

	CEBUGF("#3a\n");

	// don't bother looking at if
	//  we may add mcmc, do not add non-mc/non-mc overlaps
	//  we may not add mcmc, do not add mc/mc overlaps
	if(PAF_musthonourmcmcflag){
	  if(PAF_mayaddmcmcoverlaps){
	    if(!(*PAF_multicopies_ptr)[readid]
	       && !(*PAF_multicopies_ptr)[linkedwithid]) continue;
	  }else{
	    if((*PAF_multicopies_ptr)[readid]
	       && (*PAF_multicopies_ptr)[linkedwithid]) continue;
	  }
	}

	CEBUGF("#4\n");

	// increase potential partners before looking at blacklisted ids
	//  as ids in the blacklist can be removed over time!
	numpotentialpartners++;

	if(PAF_blacklisted_ids[linkedwithid]==0) continue;

	CEBUGF("#5\n");

	// todo hier weiter erst reads deren template partner im contig ist!
	// has linked a SRMB tag?
	int32 otherread_partnerid=PAF_readpool->getRead(linkedwithid).getTemplatePartnerID();

	// do not build in multicopy reads with a paired end partner
	//  where the partner is not in use
	if(usestricttemplatepairrule
	   && otherread_partnerid>=0
	   && PAF_multicopies[linkedwithid]
	   && !PAF_multicopies[otherread_partnerid]
	   && !PAF_used_ids[otherread_partnerid]) continue;

	// IMPORTANT!
	// notfoundyet is used to save time
	// use of the "notfoundyet" variable implicits the order
	//  of evaluation at the end of the function to see which read
	//  to really take!
	bool notfoundyet=true;

	CEBUGF("#6\n");

	// if linkedwith is a troublemaker, don't care to look at it here!
	if(PAF_readpool->getRead(linkedwithid).hasTemplateInfo() && otherread_partnerid>=0){
	  CEBUGF("#7\n");
	  if(PAF_istroublemaker[linkedwithid] == 0) {
	    CEBUGF("#8\n");
	    // is linked with template partner already in the contig?
	    if(PAF_used_ids_in_this_run[otherread_partnerid] > 0) {
	      // yes, excellent
	      CEBUGF("#9\n");
	      if(I->best_weight > bestreads[FOUND_TPARTNER_IN_CONTIG].weight) {
		CEBUGF("#x1\n");
		bestreads[FOUND_TPARTNER_IN_CONTIG].refid=readid;
		bestreads[FOUND_TPARTNER_IN_CONTIG].newid=linkedwithid;
		bestreads[FOUND_TPARTNER_IN_CONTIG].weight=I->best_weight;
		bestreads[FOUND_TPARTNER_IN_CONTIG].direction_newid=I->direction;
		bestreads[FOUND_TPARTNER_IN_CONTIG].ads_node=&(*PAF_adsfacts)[I->adsfindex];
		notfoundyet=false;
		foundsolution=true;
		needsearchnormal=false;
	      }
	    } else {
	      CEBUGF("#a\n");
	      // nope, well, still good
	      if(I->best_weight > bestreads[FOUND_TEMPLATE_AND_TPARTNER].weight) {
		CEBUGF("#x2\n");
		bestreads[FOUND_TEMPLATE_AND_TPARTNER].refid=readid;
		bestreads[FOUND_TEMPLATE_AND_TPARTNER].newid=linkedwithid;
		bestreads[FOUND_TEMPLATE_AND_TPARTNER].weight=I->best_weight;
		bestreads[FOUND_TEMPLATE_AND_TPARTNER].direction_newid=I->direction;
		bestreads[FOUND_TEMPLATE_AND_TPARTNER].ads_node=&(*PAF_adsfacts)[I->adsfindex];

		notfoundyet=false;
		foundsolution=true;
		needsearchnormal=false;
		if(PAF_multicopies[linkedwithid]==0){
		  bestreads[FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER]=bestreads[FOUND_TEMPLATE_AND_TPARTNER];
		}
	      }
	    }
	  }
	}

	CEBUGF("#b\n");

	if(PAF_readpool->getRead(linkedwithid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)) {
	  CEBUGF("#c\n");
	  // temporarily backported and adapted from n4 functions
	  //
	  // special rule: if fully contained read with 100% match
	  //  take that one immediately and break the for loop
	  // immediately take that read
	  if((*PAF_adsfacts)[I->adsfindex].getScoreRatio()==100
	     && (*PAF_adsfacts)[I->adsfindex].getOverlapLen()== PAF_readpool->getRead(linkedwithid).getLenClippedSeq()){
	    CEBUGF("#x3\n");
	    bestreads[FOUND_NORMAL].refid=readid;
	    bestreads[FOUND_NORMAL].newid=linkedwithid;
	    bestreads[FOUND_NORMAL].weight=I->best_weight;
	    bestreads[FOUND_NORMAL].direction_newid=I->direction;
	    bestreads[FOUND_NORMAL].ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    notfoundyet=false;
	    foundsolution=true;
	    speedymap=true;
	    break;
	  }
	}

	CEBUGF("#d\n");
	if(notfoundyet && needsearchnormal) {
	  CEBUGF("#e\n");
	  if(I->best_weight > bestreads[FOUND_NORMAL].weight) {
	    CEBUGF("#x4\n");
	    bestreads[FOUND_NORMAL].refid=readid;
	    bestreads[FOUND_NORMAL].newid=linkedwithid;
	    bestreads[FOUND_NORMAL].weight=I->best_weight;
	    bestreads[FOUND_NORMAL].direction_newid=I->direction;
	    bestreads[FOUND_NORMAL].ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    notfoundyet=false;
	    foundsolution=true;

	    if(PAF_multicopies[linkedwithid]==0){
	      //cout << "^^^°°°";
	      bestreads[FOUND_NONMULTICOPY]=bestreads[FOUND_NORMAL];
	    }
	  }
	}
	if(notfoundyet && needsearchnormal) {
	  CEBUGF("#f\n");
	  bool linkedhasSRMB=PAF_readpool->getRead(linkedwithid).hasTag(Read::REA_tagentry_idSRMr);

	  if(linkedhasSRMB
	     && I->best_weight > bestreads[FOUND_NORMAL_AND_SRMB].weight) {
	    CEBUGF("#x5\n");
	    bestreads[FOUND_NORMAL_AND_SRMB].refid=readid;
	    bestreads[FOUND_NORMAL_AND_SRMB].newid=linkedwithid;
	    bestreads[FOUND_NORMAL_AND_SRMB].weight=I->best_weight;
	    bestreads[FOUND_NORMAL_AND_SRMB].direction_newid=I->direction;
	    bestreads[FOUND_NORMAL_AND_SRMB].ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    notfoundyet=false;
	    foundsolution=true;
	  }
	}
      }

      // If the read we looked at had absolutely no more (un-assembled)
      //  partner reads, blacklist it
      if(numpotentialpartners==0) {
	PAF_blacklisted_ids[readid]= 0;
      }


      if(speedymap){
#ifndef PUBLICQUIET
	cout << " speedymap (backported)";
#endif
	index=-1;
	break;
      }

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      // If wished, we'll use quick rules (premature stop of searching)
      //
      // TODO: hpow about usequickrule_other combining the others?
      if((usequickruleSHORTREADS || usequickrule454 || usequickruleION || usequickruleSAN || usequickrulePBSHQ || usequickrulePBSLQ || usequickruleTXT)
	 && (bestreads[FOUND_NONMULTICOPY].weight > 0
	     || bestreads[FOUND_TPARTNER_IN_CONTIG].weight >0 )){
	uint32 index2take=FOUND_NONMULTICOPY;
	if(bestreads[FOUND_TPARTNER_IN_CONTIG].weight >0){
	  index2take=FOUND_TPARTNER_IN_CONTIG;
	}

	bool mayusequick=false;
	uint8 seqtypenewread=PAF_readpool->getRead(bestreads[index2take].newid).getSequencingType();
	switch(seqtypenewread){
	case ReadGroupLib::SEQTYPE_SANGER : {
	  if(usequickruleSAN) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_454GS20  : {
	  if(usequickrule454) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_IONTORRENT : {
	  if(usequickruleION) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_SOLEXA :
	case ReadGroupLib::SEQTYPE_ABISOLID : {
	  if(usequickruleSHORTREADS) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_PACBIOLQ : {
	  if(usequickrulePBSLQ) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_PACBIOHQ : {
	  if(usequickrulePBSHQ) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_TEXT : {
	  if(usequickruleTXT) mayusequick=true;
	  break;
	}
	default : {
	  MIRANOTIFY(Notify::INTERNAL, "Oooops? seqtypenewread not handled by switch?" << static_cast<uint16>(seqtypenewread) << endl);
	}
	}

	if(mayusequick && n4_checkQuickRules(bestreads[index2take].ads_node,bestreads[index2take].newid)){
	  index=-1;
	  break;
	}


	// This thing is a BAD idea, test with GLL 14k were catastrophic
//	// every 5 reads checked:
//	// if the score ratio is 100% (regardless of the length) take it
//	if(numreadschecked>0 && numreadschecked%5 == 0
//	   && bestreads[FOUND_NONMULTICOPY].ads_node->getScoreRatio() == 100){
//#ifndef PUBLICQUIET
//	  cout << " 454 quick rule 2";
//#endif
//	  index=-1;
//	  break;
//	}


      }


      // the "used time" foolguard for extremely high coverage areas
      //  will not search for optimum then
      if(PAF_pafparams.paf_use_emergency_search_stop){
	times(&mytms);
	actclocks=mytms.tms_utime+mytms.tms_stime;
	//cout << "Actclocks: " << actclocks << "\tmaxallowedclocks: " << maxallowedclocks << endl;
	if(actclocks>maxallowedclocks){
	  // the following should not be necessary, but just to be sure
	  // (had a case where the disk got full and MIRA was stuck in
	  //  endless loop ... calling times()
	  // so I set the loop ending condition manually
	  index=-1;

#ifndef PUBLICQUIET
	  cout << " max time, ess activated ";
#else
	  cout << "*";
#endif
	  cout.flush();

	  break;
	}
      }

    }
  }

  // this order is implicitly defined by the order of checks above when
  //  using the "notfoundyet" bool variable
  if(PAF_pafparams.paf_use_genomic_algorithms) {

    CEBUGF("FOUND_RAIL_AND_NONMULTICOPY: ");
    CEBUGF(bestreads[FOUND_RAIL_AND_NONMULTICOPY].weight << endl);
    CEBUGF("FOUND_TPARTNER_IN_CONTIG: ");
    CEBUGF(bestreads[FOUND_TPARTNER_IN_CONTIG].weight << endl);
    CEBUGF("FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER: ");
    CEBUGF(bestreads[FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER].weight << endl);
    CEBUGF("FOUND_TEMPLATE_AND_TPARTNER: ");
    CEBUGF(bestreads[FOUND_TEMPLATE_AND_TPARTNER].weight << endl);
    CEBUGF("FOUND_NONMULTICOPY: ");
    CEBUGF(bestreads[FOUND_NONMULTICOPY].weight << endl);
    CEBUGF("FOUND_NORMAL: ");
    CEBUGF(bestreads[FOUND_NORMAL].weight << endl);
    CEBUGF("FOUND_NORMAL_AND_SRMB: ");
    CEBUGF(bestreads[FOUND_NORMAL_AND_SRMB].weight << endl);

    uint32 takewhich=FOUND_NORMAL_AND_SRMB;
    if(speedymap){
      resultread=bestreads[FOUND_NORMAL];
      resultread.foundmethod=FOUND_SPEEDYMAP;
    }else{
      if(bestreads[FOUND_NORMAL].refid >= 0) takewhich=FOUND_NORMAL;
      if(bestreads[FOUND_NONMULTICOPY].weight<<4 >= bestreads[takewhich].weight) takewhich=FOUND_NONMULTICOPY;
      if(bestreads[FOUND_TEMPLATE_AND_TPARTNER].weight<<4 >= bestreads[takewhich].weight) takewhich=FOUND_TEMPLATE_AND_TPARTNER;
      if(bestreads[FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER].weight<<4 >= bestreads[takewhich].weight) takewhich=FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER;
      if(bestreads[FOUND_TPARTNER_IN_CONTIG].weight<<4 >= bestreads[takewhich].weight) takewhich=FOUND_TPARTNER_IN_CONTIG;
      if(bestreads[FOUND_RAIL_AND_NONMULTICOPY].weight<<4 >= bestreads[takewhich].weight) takewhich=FOUND_RAIL_AND_NONMULTICOPY;
      resultread=bestreads[takewhich];
      resultread.foundmethod=takewhich;
    }

  } else {

    CEBUGF("FOUND_TPARTNER_IN_CONTIG: ");
    CEBUGF(bestreads[FOUND_TPARTNER_IN_CONTIG].weight << endl);
    CEBUGF("FOUND_TEMPLATE_AND_TPARTNER: ");
    CEBUGF(bestreads[FOUND_TEMPLATE_AND_TPARTNER].weight << endl);
    CEBUGF("FOUND_NORMAL: ");
    CEBUGF(bestreads[FOUND_NORMAL].weight << endl);
    CEBUGF("FOUND_NORMAL_AND_SRMB: ");
    CEBUGF(bestreads[FOUND_NORMAL_AND_SRMB].weight << endl);


    if(speedymap){
      resultread=bestreads[FOUND_NORMAL];
      resultread.foundmethod=FOUND_SPEEDYMAP;
    }else{
      if(bestreads[FOUND_TPARTNER_IN_CONTIG].refid >= 0){
	resultread=bestreads[FOUND_TPARTNER_IN_CONTIG];
	resultread.foundmethod=FOUND_TPARTNER_IN_CONTIG;
      }else if(bestreads[FOUND_TEMPLATE_AND_TPARTNER].refid >= 0){
	resultread=bestreads[FOUND_TEMPLATE_AND_TPARTNER];
	resultread.foundmethod=FOUND_TEMPLATE_AND_TPARTNER;
      }else if(bestreads[FOUND_NORMAL].refid >= 0){
	resultread=bestreads[FOUND_NORMAL];
	resultread.foundmethod=FOUND_NORMAL;
      }else if(bestreads[FOUND_NORMAL_AND_SRMB].refid >= 0){
	resultread=bestreads[FOUND_NORMAL_AND_SRMB];
	resultread.foundmethod=FOUND_NORMAL_AND_SRMB;
      }
    }
  }

  FUNCEND();

  return;
}

//#define CEBUGF(bla)




/*************************************************************************
 *
 * Returns true if next overlap of next read to add fullfils quickrule
 *  criteria
 *
 *************************************************************************/

bool Pathfinder::n4_checkQuickRules(const AlignedDualSeqFacts * ads_node, uint32 newreadid)
{
  // Change: minlen1 & 2 can now also be negative
  // if positive: absolute value of overlap
  // if negative: relative value of the overlap to the length of the newly added read
  // overlappercentage: percentage the overlap takes from the newly added read
  uint32 overlapperc=100*ads_node->getOverlapLen()/PAF_readpool->getRead(newreadid).getLenClippedSeq();
  uint8 seqtypenewread=PAF_readpool->getRead(newreadid).getSequencingType();

  if(ads_node->getScoreRatio() >= (*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minsim1){
    if((*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minlen1 >=0){
      if(ads_node->getOverlapLen() >= (*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minlen1){
#ifndef PUBLICQUIET
	cout << " quick rule 1 abs";
#endif
	return true;
      }
    }else{
      if(overlapperc >= -(*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minlen1){
#ifndef PUBLICQUIET
	cout << " quick rule 1 rel";
#endif
	return true;
      }
    }
  }else if(ads_node->getScoreRatio() >= (*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minsim2){
    if((*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minlen2 >=0){
      if(ads_node->getOverlapLen() >= (*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minlen2){
#ifndef PUBLICQUIET
	cout << " quick rule 2 abs";
#endif
	return true;
      }
    }else{
      if(overlapperc >= -(*PAF_miraparams)[seqtypenewread].getPathfinderParams().paf_quickrule_minlen2){
#ifndef PUBLICQUIET
	cout << " quick rule 2 rel";
#endif
	return true;
      }
    }
  }

  return false;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Pathfinder::n4_constructStepByStep(vector<Align> & aligncache, vector<int8> * used_ids, vector<int32> * ids_in_contig, vector<uint8> * multicopies, vector<uint8> * hasmcoverlaps, vector<uint8> * istroublemaker, vector<vector<newedges_t>::iterator> * lowerbound_oedges_ptr, vector<uint8> * wellconnected, Contig & con)
{
  FUNCSTART("void Pathfinder::n4_constructStepByStep(vector<Align> & aligncache, vector<int8> * used_ids, vector<uint8> * multicopies, vector<uint8> * hasmcoverlaps, vector<uint8> * istroublemaker, vector<vector<newedges_t>::iterator> * lowerbound_oedges_ptr, vector<int32> * proposedSRMrightclips, Contig & con)");

  PAF_pafparams=(*PAF_miraparams)[0].getPathfinderParams();

  PAF_used_ids_ptr=used_ids;
  PAF_ids_in_contig_ptr=ids_in_contig;
  PAF_multicopies_ptr=multicopies;
  PAF_wellconnected_ptr=wellconnected;
  PAF_hasmcoverlap_ptr=hasmcoverlaps,
  PAF_istroublemaker_ptr=istroublemaker;
  PAF_lowerbound_oedges_ptr=lowerbound_oedges_ptr;

  PAF_usesbeststartcache=false;

  n4_basicCSBSSetup();

  // we're adding now things to the contig. set
  PAF_uiitr_is_clean=false;
  // and use ids_in_contig before exiting to clean PAF_used_ids_in_this_run up
  // again

  if(con.getContigLength()==0){
    CEBUG("csbs 3a\n");
    int32 startid;
    if(PAF_pafparams.paf_use_genomic_algorithms) {
      startid=n4_searchBestStrongGoodStartEnvironmentGenome();
    }else{
      startid=n4_searchBestStrongGoodStartEnvironmentEST();
    }

    CEBUG("csbs 3a0\n");

    if(startid<0) return;

    if(startid >= static_cast<int32>(PAF_used_ids_ptr->size())) {
      cout << "Startid: " << startid << endl;
      MIRANOTIFY(Notify::INTERNAL, "n4 Starting with read id greater than number iof reads?\n");
    }

    if((*PAF_used_ids_ptr)[startid]) {
      cout << "Startid: " << startid << endl;
      MIRANOTIFY(Notify::INTERNAL, "n4 Starting with read that is already used?\n");
    }

    CEBUG("Startid: " << startid << endl);

    // if we keep long repeats separated:
    //   if we start with a non-multicopy read, forbid
    //    multicopy/multicopy overlaps
    if((*PAF_multicopies_ptr)[startid]>0) {
#ifndef PUBLICQUIET
      cout << "\nStarted with multicopy.\n";
#endif
      con.setLongReapeatStatus(true);
    }

    switch(PAF_bsccontent){
    case 0:
    case 1: {
      CEBUG("csbs 3a0\n");
      // when in genome mode, use best knowledge
      if(PAF_pafparams.paf_use_genomic_algorithms) {
	n4_buildContigStartWithStrongGood(aligncache,startid, con);
      }else{
	// in EST mode, use older routines
	// TODO: change ...
	n4_buildESTContigStartWithStrongGood(aligncache,startid, con);
	//buildContig(aligncache, startid, con);
      }
      break;
    }
//    case 1:
    case 2:
    case 3: {
      // TODO: make specific n4 routine for not strong good alignments
      CEBUG("csbs 3a1\n");
      buildContig(aligncache, startid, con);
      break;
    }
    case 4: {
      CEBUG("csbs 3a3\n");
      ++PAF_readaddattempts;
      Contig::errorstatus_t contigerrstat;
      con.addRead(aligncache,
		  nullptr, startid, startid, 1,
		  (*PAF_multicopies_ptr)[startid],
		  0,
		  contigerrstat);

      (*PAF_used_ids_ptr)[startid]=1;
      ids_in_contig->push_back(startid);
      break;
    }
    default: {
      MIRANOTIFY(Notify::INTERNAL, "Oooops? PAF_bsccontent in unknown state: " << PAF_bsccontent << endl);
    }
    }
  } else {
    // TODO: make specific n4 mapping routine
    //n4_mapToContig(aligncache, -1, con);
    CEBUG("csbs 3b\n");
    buildContig(aligncache, -1, con);
  }

  BUGIFTHROW(con.getContigLength() == 0, "Contig length == 0 ... hmmm.");

  // BaCh 28.01.2012: Urgh ... what is this zerovector doing here???
  //  this destroys the possibility for the assembly to know what has been
  //  going on in pathfinder!
  //zerovector(*PAF_ids_in_contig_ptr);

  PAF_uiitr_is_clean=true;

  FUNCEND();
}
//#define CEBUG(bla)





/*************************************************************************
 *
 * some basic setup for constructStepByStep()
 *
 *
 *************************************************************************/


void Pathfinder::n4_basicCSBSSetup()
{
  FUNCSTART("void Pathfinder::n4_basicCSBSSetup()");

  PAF_pafparams=(*PAF_miraparams)[0].getPathfinderParams();

  PAF_readaddattempts=0;

#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  //PAF_ids_in_contig_ptr->clear();

  uint32 numreadspresent=PAF_used_ids_ptr->size();

  // re-initialise PAF_used_ids_in_this_run if it is empty or
  //  probably not correctly initialised
  if(PAF_used_ids_in_this_run.empty() || ! PAF_uiitr_is_clean){
    PAF_used_ids_in_this_run.clear();
    PAF_used_ids_in_this_run.resize(numreadspresent,0);
    PAF_uiitr_is_clean=true;
  }

  // these are still needed at the moment if the traditional
  //  buildcontig routines are called (e.g. for mapping assemblies)

  // is now also needed by n4 for stage 100 mapping, but needs to get
  //  initialised a bit differently (larger number). Currently done when
  //  switching to stage 100, can be done here once old routines disappear
  //  completely
  PAF_blacklisted_ids.clear();
  PAF_blacklisted_ids.resize(numreadspresent,1);

#ifdef CLOCK_STEPS2
  cout << "Timing n4_basicCSBSSetup cleararrays: " << diffsuseconds(tv) << endl;
  gettimeofday(&tv,nullptr);
#endif

  PAF_nonmulticopiesadded=0;

  PAF_musthonourmcmcflag=(*PAF_miraparams)[0].getAssemblyParams().as_keep_long_repeats_separated;
  PAF_mayaddmcmcoverlaps=true;

  CEBUG("csbs 1\n");

  // initialise all fields in the newedge_t structure of PAF_forward_edges
  //  used by the pathfinder
  // at the moment this is the pf_banned field (set when an overlap is
  //  rejected also set for all all reads that are already aligned
  //  in the contig at the potential insertion position

  if(PAF_overlapsbanned){
    if(!PAF_overlapsbanned_smallstore.empty()
       && PAF_overlapsbanned_smallstore.size() < PAF_overlapsbanned_smallstore.capacity()){
#ifndef PUBLICQUIET
      cout << "Clear pf_banned quick: " << PAF_overlapsbanned_smallstore.size() << '\n';
#endif
      vector<vector<newedges_t>::iterator>::iterator obssI=PAF_overlapsbanned_smallstore.begin();
      for(; obssI != PAF_overlapsbanned_smallstore.end(); obssI++){
	(*obssI)->pf_banned=false;
      }
    }else{
#ifndef PUBLICQUIET
      cout << "Clear pf_banned full: " << PAF_overlapsbanned << '\n';
#endif
      vector<newedges_t>::iterator neI=PAF_forward_edges->begin();
      for(; neI != PAF_forward_edges->end(); neI++){
	neI->pf_banned=false;
      }
    }
    PAF_overlapsbanned=0;
    PAF_overlapsbanned_smallstore.clear();
  }

#ifdef CLOCK_STEPS2
  cout << "Timing n4_basicCSBSSetup init pf_banned: " << diffsuseconds(tv) << endl;
  gettimeofday(&tv,nullptr);
#endif

  CEBUG("csbs 2\n");

  // initialise the quick lowerbound_oedges lookup vector if not already done
  if(PAF_lowerbound_oedges_ptr->empty()){
    vector<vector<newedges_t>::iterator> & lowerbound_oedges=*PAF_lowerbound_oedges_ptr;
    lowerbound_oedges.resize(PAF_readpool->size(),PAF_forward_edges->end());
    newedges_t tmp;
    for(uint32 i=0; i<lowerbound_oedges.size(); i++) {
      //lowerbound_oedges[i]=PAF_forward_edges->lower_bound(i);
      tmp.rid1=i;
      lowerbound_oedges[i]=lower_bound(PAF_forward_edges->begin(),
				       PAF_forward_edges->end(),
				       tmp,
				       Pathfinder__compareNewEdges_t_);
    }
#ifdef CLOCK_STEPS2
    cout << "Timing n4_basicCSBSSetup lowerbound_oedges: " << diffsuseconds(tv) << endl;
#endif
  }

  CEBUG("csbs 3\n");

#ifdef CLOCK_STEPS2
  cout << "Timing n4_basicCSBSSetup total: " << diffsuseconds(tvtotal) << endl;
#endif

  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Pathfinder::n4_handleReadNotAligned(Contig::errorstatus_t & contigerrstat, nextreadtoadd_t const &nrta)
{
  FUNCSTART("void Pathfinder::n4_handleReadNotAligned(Contig::errorstatus_t & contigerrstat)");

#ifdef PUBLICQUIET
  contigerrstat.dumpStatus();
#else
  contigerrstat.dumpStatus(true,"n4_h ");
#endif

#ifndef PUBLICQUIET
  cout.flush();
#endif
  // ban all alignments of that read with other reads that were
  //  at the same position

  //// debug: show all reads at that pos that the contig gave back
  //{
  //  cout << "\nReads affected:\n";
  //  vector<int32>::const_iterator raI;
  //  raI=contigerrstat.reads_affected.begin();
  //  for(; raI != contigerrstat.reads_affected.end(); raI++){
  //    cout << *raI << '\t' << (*PAF_readpool)[*raI].getName();
  //    if((*PAF_readpool)[*raI].isBackbone()){
  //      cout << "\tbackbone\n";
  //    } else if((*PAF_readpool)[*raI].isRail()){
  //      cout << "\trail\n";
  //    }else{
  //      cout << "\tseqtype: " << static_cast<uint16>((*PAF_readpool)[*raI].getSequencingType()) << '\n';
  //    }
  //  }
  //}


  {
    // check if refid has been put into reads_affected by the contig
    // if not, warn and put it in
    bool isin=false;
    for(auto raid : contigerrstat.reads_affected){
      if(raid == nrta.refid){
	isin=true;
	break;
      }
    }
    if(!isin){
      cout << "\nWARNING: contig did not put refid into affected reads. Please file a bug report to\n\thttps://sourceforge.net/apps/trac/mira-assembler/\nand please send a mail to\n\tmira_talk@freelists.org\n";
      cout << "current af size: " << contigerrstat.reads_affected.size() << endl;
      contigerrstat.reads_affected.push_back(nrta.refid);
    }
  }


  {
    // First search look at the newedges_t block of nrta,newid overlaps
    //  and ban those overlaps whose reads are in the vector given back
    //  by the contig

    vector<newedges_t>::iterator neI = (*PAF_lowerbound_oedges_ptr)[nrta.newid];
    vector<int32>::const_iterator raI;
    for(; neI != PAF_forward_edges->end() && neI->rid1 == nrta.newid; neI++){
      if(!neI->pf_banned){
	raI=contigerrstat.reads_affected.begin();
	for(; raI != contigerrstat.reads_affected.end(); raI++){
	  if(*raI == neI->linked_with){
	    //cout << "Banning1: " << neI->rid1 << '\t' << neI->linked_with << '\n';
	    PAF_overlapsbanned++;
	    neI->pf_banned=true;
	    if(PAF_overlapsbanned_smallstore.size()<PAF_overlapsbanned_smallstore.capacity()){
	      PAF_overlapsbanned_smallstore.push_back(neI);
	    }
	    break;
	  }
	}
      }
    }

    // Now look at the blocks of newedges_t whose reads are given
    //  by by the vector and who overlap with nrta.newid
    raI=contigerrstat.reads_affected.begin();
    for(; raI != contigerrstat.reads_affected.end(); raI++){
      neI = (*PAF_lowerbound_oedges_ptr)[*raI];
      for(; neI != PAF_forward_edges->end() && neI->rid1 == *raI; neI++){
	if(!neI->pf_banned){
	  if(neI->linked_with == nrta.newid){
	    //cout << "Banning2: " << neI->rid1 << '\t' << neI->linked_with << '\n';
	    PAF_overlapsbanned++;
	    neI->pf_banned=true;
	    if(PAF_overlapsbanned_smallstore.size()<PAF_overlapsbanned_smallstore.capacity()){
	      PAF_overlapsbanned_smallstore.push_back(neI);
	    }
	    break;
	    // NOTE: should MIRA ever get multiple overlaps between
	    //  two reads, then the "break;" above must be removed
	    //  should be safe for now.
	  }
	}
      }
    }
  }

  FUNCEND();
  return;
}



///*************************************************************************
// *
// *  returns ID of read with best matching env
// *
// *
// *************************************************************************/
//
//
////#define CEBUG(bla)   {cout << bla; cout.flush(); }
//
//int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentEST()
//{
//  FUNCSTART("int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentEST()");
//
//  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;
//  const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;
//
//  int32 bestid=-1;
//  uint64 bestweight=0;
//  bool foundid=false;
//
//
//#ifdef CLOCK_STEPS2
//  timeval tv;
//  timeval tvtotal;
//  gettimeofday(&tv,nullptr);
//  tvtotal=tv;
//#endif
//
//  CEBUG("sbe 1\n");
//
//  for(int32 actid=0; actid < static_cast<int32>(PAF_used_ids.size()); actid++) {
//    CEBUG("aid: " << actid);
//    CEBUG("\tused " << (int16) PAF_used_ids[actid]);
//    CEBUG("\ttm " <<(uint16) PAF_istroublemaker[actid]);
//    CEBUG("\n");
//    if(PAF_used_ids[actid]!=0
//       || PAF_istroublemaker[actid]!=0) continue;
//
//    vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[actid];
//
//    uint32 sumweight=0;
//
//    for(;I!=PAF_forward_edges->end() && actid==I->rid1;I++){
//      CEBUG("  lid: " << I->linked_with);
//      CEBUG("\tused " << (int16) PAF_used_ids[I->linked_with]);
//      CEBUG("\ttm " <<(uint16) PAF_istroublemaker[I->linked_with]);
//      CEBUG("\n");
//      CEBUG("  " << *I);
//
//      if(!I->ol_stronggood
//	 || PAF_used_ids[I->linked_with]!=0
//	 || PAF_istroublemaker[I->linked_with]!=0) continue;
//
//      // add weight to normal overlaps
//      sumweight+=I->best_weight;
//    }
//
//    CEBUG("Sumweight: " << sumweight << '\n');
//
//    if(sumweight>bestweight){
//      foundid=true;
//      bestid=actid;
//      bestweight=sumweight;
//      CEBUG("New best weight: " << bestweight << '\n');
//    }
//  }
//
//  // we do not try to define singlets or troublemakers as Non-mc/non-mc
//  //  startpoint here. Let the caller decide what to do.
//
//#ifdef CLOCK_STEPS2
//  cout << "Timing n4_searchBestStrongGoodStartEnvironmentEST total: " << diffsuseconds(tvtotal) << endl;
//#endif
//
//  CEBUG("Returning: " << bestid << '\t' << bestweight << '\n');
//
//  FUNCEND();
//
//  return bestid;
//
//  //  return 0;
//}
////#define CEBUG(bla)   {cout << bla; cout.flush(); }



///*************************************************************************
// *
// *  returns ID of read with best matching env
// *
// *
// *************************************************************************/
//
//
////#define CEBUG(bla)   {cout << bla; cout.flush(); }
//
//int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentGenome()
//{
//  FUNCSTART("int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentGenome()");
//
//  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;
//  const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;
//
//  int32 bestid=-1;
//  uint64 bestweight=0;
//
//  uint32 largestclustersize=0;
//
//  bool foundid=false;
//
//
//#ifdef CLOCK_STEPS2
//  timeval tv;
//  timeval tvtotal;
//  gettimeofday(&tv,nullptr);
//  tvtotal=tv;
//#endif
//
//  vector<bool> uid_in_cluster(PAF_used_ids.size(),false);
//  vector<int32> unlooked;
//  unlooked.reserve(PAF_used_ids.size());
//
//  CEBUG("sbsgseg 1\n");
//
//  for(int32 actid=0; actid < static_cast<int32>(PAF_used_ids.size()); actid++) {
//    CEBUG("aid: " << actid);
//    CEBUG("\tused " << (int16) PAF_used_ids[actid]);
//    CEBUG("\ttm " <<(uint16) PAF_istroublemaker[actid]);
//    CEBUG("\n");
//    if(PAF_used_ids[actid]!=0
//       || PAF_istroublemaker[actid]!=0
//       || uid_in_cluster[actid]) continue;
//
//    uint32 goodclustersize=0;
//    uint64 goodclusterbestweight=0;
//    uint32 goodclusterbestid=0;
//    unlooked.clear();
//    unlooked.push_back(actid);
//
//    int32 lookid;
//    while(!unlooked.empty()){
//      CEBUG("unlooked.size(): " << unlooked.size() << endl);
//      lookid=unlooked.back();
//      unlooked.pop_back();
//
//      uid_in_cluster[lookid]=true;
//
//      vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[lookid];
//
//      uint64 sumweight=0;
//      for(;I!=PAF_forward_edges->end() && lookid==I->rid1;I++){
//	CEBUG("  lid: " << I->linked_with);
//	CEBUG("\tused " << (int16) PAF_used_ids[I->linked_with]);
//	CEBUG("\ttm " <<(uint16) PAF_istroublemaker[I->linked_with]);
//	CEBUG("\n");
//	CEBUG("  " << *I);
//
//	if(!I->ol_stronggood
//	   || PAF_used_ids[I->linked_with]!=0
//	   || PAF_istroublemaker[I->linked_with]!=0
//	   || uid_in_cluster[I->linked_with]) continue;
//	//if(goodclustersize==0) ++goodclustersize;
//	//++goodclustersize;
//	sumweight+=I->best_weight;
//	unlooked.push_back(I->linked_with);
//      }
//      if(sumweight>0) ++goodclustersize;
//
//      if(sumweight>goodclusterbestweight){
//	goodclusterbestid=lookid;
//	goodclusterbestweight=sumweight;
//	CEBUG("goodclusterbestid: " << goodclusterbestid << endl);
//	CEBUG("goodclusterbestweight: " << goodclusterbestweight << endl);
//      }
//    }
//    if(goodclusterbestweight>0) ++goodclustersize;
//
//    if(goodclustersize>largestclustersize){
//      foundid=true;
//      bestid=goodclusterbestid;
//      bestweight=goodclusterbestweight;
//      largestclustersize=goodclustersize;
//      CEBUG("New best id: " << bestid << endl);
//      CEBUG("New best weight: " << bestweight << endl);
//      CEBUG("New largest cluster: " << largestclustersize << '\n');
//    }
//  }
//
//  // we do not try to define singlets or troublemakers as Non-mc/non-mc
//  //  startpoint here. Let the caller decide what to do.
//
//#ifdef CLOCK_STEPS2
//  cout << "Timing n4_searchBestStrongGoodStartEnvironmentGenome total: " << diffsuseconds(tvtotal) << endl;
//#endif
//
//  CEBUG("Best id: " << bestid << endl);
//  CEBUG("Best weight: " << bestweight << endl);
//  CEBUG("Largest cluster: " << largestclustersize << '\n');
//  CEBUG("Returning: " << bestid << '\t' << bestweight << '\n');
//
//  FUNCEND();
//
//  return bestid;
//
//  //  return 0;
//}
////#define CEBUG(bla)




//#define CEBUG(bla)   {cout << bla; cout.flush(); }

bool Pathfinder::sortbeststartinfo_t_(const beststartinfo_t & a, const beststartinfo_t & b)
{
  return a.bsi_clustersize < b.bsi_clustersize;
}

/*************************************************************************
 *
 * returns ID of read with best matching env for genome/EST assemblies
 *
 * ! sets PAF_bsccontent to give assessement of what start sites
 *    are in the cache:
 *       0 == best quality, strong good overlaps found
 *       1 == medium quality, overlaps but not strong good, not multicopy
 *       2 == medium quality, overlaps but not strong good. multicopy
 *       3 == singlets
 *
 *
 * As the genome routines became more streamlined, same routines for genome
 *  and EST and use cache. Now also do not look at different status anymore,
 *  this has been merged for the time being
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }

int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentEST()
{
  FUNCSTART("int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentEST()");
  FUNCEND();

  return n4_searchBestStrongGoodStartEnvironment_sub();
}

int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentGenome()
{
  FUNCSTART("int32 Pathfinder::n4_searchBestStrongGoodStartEnvironmentGenome()");
  FUNCEND();
  return n4_searchBestStrongGoodStartEnvironment_sub();
}


/*************************************************************************
 *
 *  fills the cache for the best start sites
 *
 *************************************************************************/


void Pathfinder::n4_searchBestStrongGoodStartEnvironment_subFillCache(bool wanttroublemakercheck, bool wantstronggoodcheck, bool wantmulticopycheck, bool wantwellconnectedcheck)
{
  FUNCSTART("void Pathfinder::n4_searchBestStrongGoodStartEnvironment_subFillCache(bool trouble, bool stronggood)");

  const auto & PAF_used_ids = *PAF_used_ids_ptr;
  const auto & PAF_istroublemaker = *PAF_istroublemaker_ptr;
  const auto & PAF_multicopies = *PAF_multicopies_ptr;
  const auto & PAF_wellconnected = *PAF_wellconnected_ptr;

  vector<bool> uid_in_cluster(PAF_used_ids.size(),false);

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

  vector<int32> unlooked;
  unlooked.reserve(8000);
  bool unlookedaccepts=true;
  bool unlookthreshreached=false;

  PAF_beststartcache.reserve(10000);

  long int maxsec=PAF_pafparams.paf_max_startcache_filltime;
  timeval starttv;
  timeval curtv;
  gettimeofday(&starttv,nullptr);

  for(int32 actid=0; actid < static_cast<int32>(PAF_used_ids.size()); ++actid) {
    CEBUG("aid: " << actid);
    CEBUG("\tused " << (int16) PAF_used_ids[actid]);
    CEBUG("\ttm " <<(uint16) PAF_istroublemaker[actid]);
    CEBUG("\n");
    if(PAF_used_ids[actid]!=0
       || (wanttroublemakercheck && PAF_istroublemaker[actid]!=0)
       || (wantmulticopycheck && PAF_multicopies[actid]!=0)
       || (wantwellconnectedcheck && !PAF_wellconnected[actid])
       || uid_in_cluster[actid]) continue;

    uint32 goodclustersize=0;
    uint32 goodclusterbestid=0;
    uint32 maxconnects=0;
    unlooked.clear();
    unlooked.push_back(actid);

    unlookedaccepts=true;

    int32 lookid;
    while(!unlooked.empty()){
      CEBUG("unlooked.size(): " << unlooked.size() << endl);
      lookid=unlooked.back();
      unlooked.pop_back();

      uid_in_cluster[lookid]=true;

      vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[lookid];

      uint32 numconnects=0;
      for(;I!=PAF_forward_edges->end() && lookid==I->rid1;++I){
	CEBUG("  lid: " << I->linked_with);
	CEBUG("\tused " << (int16) PAF_used_ids[I->linked_with]);
	CEBUG("\ttm " <<(uint16) PAF_istroublemaker[I->linked_with]);
	CEBUG("\n");
	CEBUG("  " << *I);

	if((wantstronggoodcheck && !I->ol_stronggood)
	   || PAF_used_ids[I->linked_with]!=0
	   || (wanttroublemakercheck && PAF_istroublemaker[I->linked_with]!=0)
	   || (wantmulticopycheck && PAF_multicopies[I->linked_with]!=0)
	   || (wantwellconnectedcheck && !PAF_wellconnected[I->linked_with])
	   || uid_in_cluster[I->linked_with]) continue;
	++numconnects;

	if(unlookedaccepts) {
	  unlooked.push_back(I->linked_with);
	  // Check whether to stop looking for more
	  // Do not do that only when multicopies are also allowed as it may be that
	  //  people give projects *just* with heavy multicopy clusters and then
	  //  it would not work (as only a monirity may then be multicopy with
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
      PAF_beststartcache.push_back(tmp);
    }

    // see whether we need to care about time (-PF:mscft)
    if(maxsec>=0 && actid%50==0 && !PAF_beststartcache.empty()) {
      gettimeofday(&curtv,nullptr);
      if(curtv.tv_sec-starttv.tv_sec >= maxsec) {
	cout << "Non-deterministic behaviour of assembly likely: -PF:mscft threshold hit.\n";
	break;
      }
    }

  }

  if(unlookthreshreached){
    cout << "hit unlooked threshold\n";
  }

  FUNCEND();
}

/*************************************************************************
 *
 *  returns ID of read with best matching env
 *
 *************************************************************************/

int32 Pathfinder::n4_searchBestStrongGoodStartEnvironment_sub()
{
  FUNCSTART("int32 Pathfinder::n4_searchBestStrongGoodStartEnvironment_sub()");

  int32 bestid=-1;
  const vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;

  PAF_usesbeststartcache=true;

  // if there's something in the cache and
  //  the read it points to hasn't been included via other means
  //  (like jumping over a small repeat), return that cached entry
  while(!PAF_beststartcache.empty()){
    CEBUG("Startcache look " << PAF_beststartcache.back());
    bestid=PAF_beststartcache.back().bsi_rid;
    PAF_beststartcache.pop_back();

    if(!PAF_used_ids[bestid]){
      CEBUG("Startcache taken\n");
      return bestid;
    }
  }
  // else try to find / recompute something else
  bestid=-1;

  //const vector<uint8> & PAF_istroublemaker = *PAF_istroublemaker_ptr;

  uint32 largestclustersize=0;

#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif


  CEBUG("sbsgseg 1\n");

  PAF_bsccontent=0;
  n4_searchBestStrongGoodStartEnvironment_subFillCache(true,true,true,true);

#ifdef CLOCK_STEPS2
  cout << "Timing n4_searchBestStrongGoodStartEnvironment 1 1 1 1: " << diffsuseconds(tvtotal)
       << "\nStartcache size: " << PAF_beststartcache.size() << endl;
#endif

  if(PAF_beststartcache.empty()){
#ifdef CLOCK_STEPS2
    gettimeofday(&tv,nullptr);
#endif
    PAF_bsccontent=1;
    n4_searchBestStrongGoodStartEnvironment_subFillCache(false,false,true,true);
#ifdef CLOCK_STEPS2
    cout << "Timing n4_searchBestStrongGoodStartEnvironment 0 0 1 1: " << diffsuseconds(tv)
	 << "\nStartcache size: " << PAF_beststartcache.size() << endl;
#endif
    if(PAF_beststartcache.empty()){
#ifdef CLOCK_STEPS2
      gettimeofday(&tv,nullptr);
#endif
      PAF_bsccontent=2;
      n4_searchBestStrongGoodStartEnvironment_subFillCache(false,false,false,true);
#ifdef CLOCK_STEPS2
      cout << "Timing n4_searchBestStrongGoodStartEnvironment 0 0 0 1: " << diffsuseconds(tv)
	   << "\nStartcache size: " << PAF_beststartcache.size() << endl;
#endif
    }
    if(PAF_beststartcache.empty()){
#ifdef CLOCK_STEPS2
      gettimeofday(&tv,nullptr);
#endif
      PAF_bsccontent=3;
      n4_searchBestStrongGoodStartEnvironment_subFillCache(false,false,false,false);
#ifdef CLOCK_STEPS2
      cout << "Timing n4_searchBestStrongGoodStartEnvironment 0 0 0 0: " << diffsuseconds(tv)
	   << "\nStartcache size: " << PAF_beststartcache.size() << endl;
#endif
    }
  }

  if(!PAF_beststartcache.empty()){
#ifdef CLOCK_STEPS2
    gettimeofday(&tv,nullptr);
#endif
    sort(PAF_beststartcache.begin(),
	 PAF_beststartcache.end(),
	 Pathfinder::sortbeststartinfo_t_
      );
#ifdef CLOCK_STEPS2
    cout << "Timing n4_searchBestStrongGoodStartEnvironment sort: " << diffsuseconds(tv) << endl;
#endif
  }else{
    // found nothing ... now put singlets into list

#ifdef CLOCK_STEPS2
    gettimeofday(&tv,nullptr);
#endif

    PAF_bsccontent=4;

    beststartinfo_t tmp;
    tmp.bsi_clustersize=0;
    tmp.bsi_numconnects=0;

    uint32 size=PAF_readpool->size();
    for(uint32 i=0; i<size; i++){
      if(PAF_used_ids[i]==0){
	tmp.bsi_rid=i;
	PAF_beststartcache.push_back(tmp);
      }
    }
    CEBUG("Startcache filled with singlets: " << PAF_beststartcache.size() << '\n');
#ifdef CLOCK_STEPS2
    cout << "Timing n4_searchBestStrongGoodStartEnvironment singletfill: " << diffsuseconds(tv) << endl;
#endif
  }

  if(!PAF_beststartcache.empty()){
    CEBUG("Startcache take " << PAF_beststartcache.back());
    bestid=PAF_beststartcache.back().bsi_rid;
    largestclustersize=PAF_beststartcache.back().bsi_clustersize;
    PAF_beststartcache.pop_back();
  }

#ifdef CLOCK_STEPS2
  cout << "Timing n4_searchBestStrongGoodStartEnvironment total: " << diffsuseconds(tvtotal)
       << "\n";
#endif

  CEBUG("Best id: " << bestid << endl);
  CEBUG("Largest cluster: " << largestclustersize << '\n');

  FUNCEND();

  return bestid;

  //  return 0;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *  returns number of elements in start cache
 *  throw if startcache is not used!
 *
 *************************************************************************/

size_t Pathfinder::n4_getNumElementsInStartCache()
{
  FUNCSTART("size_t Pathfinder::n4_getNumElementsInStartCache()");
  //BUGIFTHROW(!PAF_usesbeststartcache,"Asked for elements in start cache while not using startcache?");

  // TODO: cludge, until all routines are migrated to n4* and use beststartcache
  if(!PAF_usesbeststartcache) return 1;

  FUNCEND();
  return PAF_beststartcache.size();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Pathfinder::n4_buildContigStartWithStrongGood(vector<Align> & aligncache, int32 startid, Contig & con)
{
  FUNCSTART("void Pathfinder::n4_buildContigStartWithStrongGood(vector<Align> & aligncache, int32 startid, Contig & con)");

#ifndef PUBLICQUIET
  cout << "n4_buildContigStartWithStrongGood\n";
#endif

  CEBUG("bc 1\n");

  vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;

  Contig::errorstatus_t contigerrstat;
  // this should speed up things a bit, avoiding copying during extension
  //  although the 10000 might be not enough for really pathological cases
  //  with many SRMB regions and huge coverage
  contigerrstat.reads_affected.reserve(10000);

  // this is a kind of duplicate of PAF_used_ids_in_this_run (which has
  //  the size of the readpool)
  // However, this array holds only those ids that have been built into
  //  the contig in this call, giving a faster access to the info
  vector<int32> & ids_in_contig = *PAF_ids_in_contig_ptr;
  ids_in_contig.clear();
  ids_in_contig.reserve(PAF_used_ids.size());

  // "duplicate" of ids_in_contig, but gets cleared after each
  //  stage1, stage2, ..., stageN round below.
  // Even faster access to the really newest additions to the contig
  //vector<int32> PAF_fresh_ids_in_contig;
  PAF_fresh_ids_in_contig.clear();
  PAF_fresh_ids_in_contig.reserve(PAF_used_ids.size());

  // "duplicate" of ids_in_contig, but contains ids that have overlaps
  //  with unassembled reads and therefore can contribute to growth
  // Once IDs have no open overlap anymore, they become dead ("cork")
  //  and are removed from the cambium vector (by the n4_find* functions)
  // On each change of stage, cambium_ids gets refilled with
  //  ids_in_contig to be sure that nothing gets lost
  //vector<int32> & PAF_cambium_ids;
  PAF_cambium_ids.clear();
  PAF_cambium_ids.reserve(PAF_used_ids.size());


  PAF_buildcontig_newlinecounter=0;


  CEBUG("bc 3a 1\n");
  if(startid<0) {
    cerr << "Startid: " << startid << endl;
    MIRANOTIFY(Notify::INTERNAL, "Starting with non-existing read?");
  }

  if(startid >= static_cast<int32>(PAF_used_ids.size())) {
    cerr << "Startid: " << startid << endl;
    MIRANOTIFY(Notify::INTERNAL, "Starting with read id greater than number iof reads?\n");
  }

  if(PAF_used_ids[startid]) {
    cerr << "Startid: " << startid << endl;
    MIRANOTIFY(Notify::INTERNAL, "Starting with read that is already used?\n");
  }

  CEBUG("Startid: " << startid << endl);

#ifndef PUBLICQUIET
  cout << endl;
#endif

  ++PAF_readaddattempts;
  con.addRead(aligncache,
	      nullptr, startid, startid, 1,
	      (*PAF_multicopies_ptr)[startid],
	      0,
	      contigerrstat);

  CEBUG("bc 3a 2\n");

  PAF_used_ids[startid]=1;
  PAF_used_ids_in_this_run[startid]=1;
  ids_in_contig.push_back(startid);
  PAF_fresh_ids_in_contig.push_back(startid);
  PAF_cambium_ids.push_back(startid);

#ifndef PUBLICQUIET
  cout << endl;
#else
  cout << "+";
#endif

  struct tms mytms;
  times(&mytms);
  clock_t baseclocks=mytms.tms_utime+mytms.tms_stime;
  clock_t actclocks=baseclocks;
  clock_t maxallowedclocks=baseclocks+PAF_pafparams.paf_maxcontigclockticks;

#ifdef CLOCK_STEPS1
  vector<suseconds_t> us_pfsearch;
  vector<suseconds_t> us_conadd;
  vector<suseconds_t> us_overlapban;
#endif

  nextreadtoadd_t nrta;
  nrta.newid=0;

  uint32 statscounter=1;

  // the forcegrow for addRead()
  //  will be set to 0 for adding NMC, 10-20? for MC
  int32 forcegrow=0;

  // stage 1 : add "strong good" overlaps
  // stage 2 : add paired-ends only
  // stage 3 : add overlaps that have frequencies below avg.

  uint32 stage=1;
  uint32 readsaddedinstages=0;

  /* forcenewstage is set by a stage if it explicitly wants to change
     to a new stage after the current found overlap has been added
     value of forcenewstage is the target stage */
  uint32 forcenewstage=0;

  while(stage < 1000) {
#ifndef PUBLICQUIET
    if(statscounter == 25) {
      cout << "rstatsm: len of c : " << con.getContigLength();
      cout << "\nrstatsm: ids in c : " << ids_in_contig.size();
      cout << "\nrstatsm: cids in c : " << PAF_cambium_ids.size();
      cout << "\nrstatsm: fids in c : " << PAF_fresh_ids_in_contig.size();
      cout << "\nrstatsi: banned o : " << PAF_overlapsbanned;
#ifdef CLOCK_STEPS1
      cout << "\nrstat timings: " << ids_in_contig.size();
      cout << " : " << avg_suseconds(us_pfsearch);
      cout << " / " << avg_suseconds(us_conadd);
      cout << " / " << avg_suseconds(us_overlapban);
#endif
      cout << '\n';
      statscounter=0;
    }
#endif


    nrta.refid=-1;
    nrta.newid=-1;
    nrta.weight=0;
    nrta.direction_newid=0;
    nrta.ads_node=nullptr;

#ifdef CLOCK_STEPS1
    timeval us_start;
    gettimeofday(&us_start,nullptr);
#endif

#ifndef PUBLICQUIET
    cout << "Searching";
#endif

    switch(stage){
    case 1 : {
#ifndef PUBLICQUIET
      cout << " (stage1 overlap) ...";
      cout.flush();
#endif
      forcegrow=0;
      n4_findNextGoodOverlap(true,false,true,PAF_cambium_ids, nrta, maxallowedclocks);

      break;
    }
    case 2 : {
#ifndef PUBLICQUIET
      cout << " (stage2 overlap) ...";
      cout.flush();
#endif
      forcegrow=0;
      n4_findNextPillarCounterpart(PAF_cambium_ids, nrta, maxallowedclocks);

      break;
    }
    case 3 : {
#ifndef PUBLICQUIET
      cout << " (stage3 overlap) ...";
      cout.flush();
#endif
      forcegrow=0;
      if(n4_findNextGoodOverlap(true,true,true,PAF_cambium_ids, nrta, maxallowedclocks)){
	// if function returns true, it has found that a read in the
	//  contig has a strong good overlap
	forcenewstage=1;
      }

      break;
    }
    case 100 : {
      // map to existing, allowing only overlaps with data
      //  that has a frequency >3
#ifndef PUBLICQUIET
      cout << " (stage100 overlap) ...";
      cout.flush();
#endif
      // disallow the contig to grow!
      forcegrow=-1;

      // At the moment we use the cambium ids anyway as allowed references
      //if(PAF_cambium_ids.empty()){
      //	n4_setupAllowedReferences(PAF_cambium_ids,ids_in_contig);
      //}

      n4_findNextAllowedRefOverlap(PAF_cambium_ids, nrta, maxallowedclocks);

      // nothing found? the clear allowed references
      if(nrta.newid < 0) {
	PAF_cambium_ids.clear();
      }
      break;
    }
    default : {
      cerr << "\nStage " << stage << '\n';
      MIRANOTIFY(Notify::INTERNAL, "stage not 1,2,3 ... this should not happen.\n");
    }
    }
#ifndef PUBLICQUIET
    cout << " done.\n";
    cout.flush();
#endif

#ifdef CLOCK_STEPS1
    us_pfsearch.push_back(diffsuseconds(us_start));
#endif


    if(nrta.newid < 0) {
      // nothing found
      // go to next stage and if we went through all stages,
      //  go back to stage 1 only if there had been a read added
      //  in the three stages. if not, this contigs is finished
#ifndef PUBLICQUIET
      cout << "\nWe were in stage " << stage << ", switching to stage " << stage+1 << endl;
      cout << "Reads added: " << readsaddedinstages << endl;
#endif
      ++stage;
      if(stage==4){
	if(readsaddedinstages){
	  // we added some reads in the stages. Perfect, restart
	  //  at stage 1 to see whether we could get further.
	  stage=1;
	  readsaddedinstages=0;
	}else{
	  // nope, so start the beef up loop now
	  stage=100;

	  uint32 numreadspresent=PAF_used_ids_ptr->size();
	  // is now also needed by n4 for stage 100 mapping, but needs to get
	  //  initialised a bit differently (larger number). Currently done here
	  //  switching to stage 100, can be done in basicSetup once old
	  //  routines disappear
	  PAF_blacklisted_ids.clear();
	  PAF_blacklisted_ids.resize(PAF_used_ids_ptr->size(),5);
	}
	PAF_fresh_ids_in_contig.clear();
      }else if(stage>=100){
	//stop the loop now
	stage=1000;
      }
      PAF_cambium_ids=ids_in_contig;
    }else{
#ifndef PUBLICQUIET
      cout << nrta;
#endif

#ifdef CLOCK_STEPS1
      gettimeofday(&us_start,nullptr);
#endif

      statscounter++;
      ++PAF_readaddattempts;
      con.addRead(aligncache,
		  nrta.ads_node, nrta.refid, nrta.newid, nrta.direction_newid,
		  (*PAF_multicopies_ptr)[nrta.newid],
		  forcegrow,
		  contigerrstat);
#ifdef CLOCK_STEPS1
      us_conadd.push_back(diffsuseconds(us_start));
#endif

      cout.flush();

      //Contig::setCoutType(Contig::AS_TEXT);
      //cout << ",,,,,,\n" << con << "\n\n";


      if(PAF_buildcontig_newlinecounter==0){
	cout << "[" << ids_in_contig.size() << "] "; cout.flush();
      }
      if(contigerrstat.code == Contig::ENOERROR) {
	ids_in_contig.push_back(nrta.newid);
	PAF_fresh_ids_in_contig.push_back(nrta.newid);
	PAF_cambium_ids.push_back(nrta.newid);
	++readsaddedinstages;
	PAF_used_ids[nrta.newid]=1;
	PAF_used_ids_in_this_run[nrta.newid]=1;

#ifndef PUBLICQUIET
	cout << "\t+\n";
#else
	cout << "+";
#endif
      }else{
#ifdef CLOCK_STEPS1
	gettimeofday(&us_start,nullptr);
#endif
	n4_handleReadNotAligned(contigerrstat,nrta);

	// handle blacklisting of failed alignment when refid was a short read
	//  and newid a long one
	if(stage==100
	   && PAF_readpool->getRead(nrta.refid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	   && !PAF_readpool->getRead(nrta.newid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	  if(PAF_blacklisted_ids[nrta.newid]) PAF_blacklisted_ids[nrta.newid]--;
	}
#ifdef CLOCK_STEPS1
	us_overlapban.push_back(diffsuseconds(us_start));
#endif
      }

      if(forcenewstage>0){
#ifndef PUBLICQUIET
	cout << "\nWe were in stage " << stage << ", being forced to stage " << stage+1 << endl;
	cout << "Reads added: " << readsaddedinstages << endl;
#endif
	readsaddedinstages=0;
	stage=forcenewstage;
	forcenewstage=0;
	PAF_fresh_ids_in_contig.clear();
	PAF_cambium_ids=ids_in_contig;
      }

#ifndef PUBLICQUIET
#else
      PAF_buildcontig_newlinecounter++;
      if(PAF_buildcontig_newlinecounter==60){
	PAF_buildcontig_newlinecounter=0;
	cout << "   " << con.getContigLength();
#ifdef CLOCK_STEPS1
	cout << "\t" << avg_suseconds(us_pfsearch);
	cout << " / " << avg_suseconds(us_conadd);
	cout << " / " << avg_suseconds(us_overlapban);
#endif
	cout << endl;
      }
#endif


    }

    //if(con.getNumReadsInContig()==1198){
    //  cout << "\npathfinder debugstop.\n";
    //  assout::saveAsCAF(con, "test.caf", true);
    //  abort();
    //}

    if(PAF_pafparams.paf_use_max_contig_buildtime){
      times(&mytms);
      actclocks=mytms.tms_utime+mytms.tms_stime;
      if(actclocks>maxallowedclocks){
	// the following should not be necessary, but just to be sure
	// (had a case where the disk got full and MIRA was stuck in
	//  endless loop ... calling times()
	// so I set the loop ending condition manually
	nrta.newid=-1;

	cout << "\nMaximum build time for this contig reached, aborting build.\n";
	FUNCEND();
	return;
      }
    }
  }

  FUNCEND();
}

//#define CEBUG(bla)





/*************************************************************************
 *
 * searches for two types of overlaps:
 *   allowbelowavgfreq == false -> only overlaps with a ol_strong (ol_weak)
 *                                 flag
 *                     == true -> overlaps with a ol_belowavgfreq
 *                                             or ol_norept flag
 *
 * returns whether the overlap chosen also has ol_strong (ol_weak) flag
 *    set
 *
 *************************************************************************/

//#define CEBUGF(bla)   {cout << bla; cout.flush(); }
bool Pathfinder::n4_findNextGoodOverlap(bool wantfreqcheck, bool allowbelowavgfreq, bool allowspeedypair, vector<int32> & cambiumidsincontig, nextreadtoadd_t & resultread, clock_t & maxallowedclocks)
{
  FUNCSTART("void Pathfinder::n4_findNextStrongGoodOverlap(bool allowbelowavgfreq, const vector<int32> & cambiumidsincontig, nextreadtoadd_t & resultread, clock_t & maxallowedclocks) ");

  if(cambiumidsincontig.size()==0) {
    return false;
  }

  bool returnvalue=false;

  // lr_ == local reference
  const vector<int8> &  lr_used_ids = *PAF_used_ids_ptr;
  const vector<uint8> & lr_istroublemaker = *PAF_istroublemaker_ptr;

  nextreadtoadd_t bestread=resultread;
  vector<newedges_t>::const_iterator bestreadI=PAF_forward_edges->end();
  bool bestreadispaired=false;
  uint32 peweighttobeat=0;


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif


  bool   usequickruleSAN=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrule454=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_use_quick_rule;
  bool   usequickruleION=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrulePBSHQ=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_PACBIOHQ].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrulePBSLQ=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_PACBIOLQ].getPathfinderParams().paf_use_quick_rule;
  bool   usequickruleTXT=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_TEXT].getPathfinderParams().paf_use_quick_rule;

  bool   usequickruleSHORTREADS=true;

  bool foundquick=false;

  // going backwards through the contglist should be a tick faster than forward
  //  also helps to retain a bit of flexibility when emergency search stops occur
  int32 rindex=cambiumidsincontig.size()-1;

  uint32 numreadschecked=0;


//#define CEBUGF(bla)   {if(numreadschecked>10000) cout << bla; cout.flush(); }


  for(; rindex >=0 ; --rindex, ++numreadschecked) {

    // TODO: not optimal: should perhaps set flag to sort cambium ids according to
    //  untaken overlaps? (every now and then, perhaps once cambium grows more
    //  than 1000 since last sort)
    // Test: stop runaway searches
    if(numreadschecked>3000
       && bestread.weight>0){
      rindex=-1;
      break;
    }

    int32 readid=cambiumidsincontig[rindex];

    // foundcandidate tells whether for this read id there are overlaps
    //  to unassembled reads
    // if not, the this read id is "dead" and is removed from the
    //  cambiumidsincontig further down (it becomes "cork" :-)
    bool foundcandidate=false;

    // speedymap: special quickrule for 100% contained reads that map
    //  100% but are shorter than the quickrule lengths (important for
    //  hybrid Solexa / other)
    bool speedymap=false;

    // speedypair: special quickrule for reads which have a partner in
    //  the contig
    bool speedypair=false;

    CEBUGF("fngO: new read " << readid << '\n');

    vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[readid];

    for(;I!=PAF_forward_edges->end() && I->rid1 == readid;I++){

      CEBUGF("edge: " << *I << endl);

      // don't bother looking at if overlap is banned
      if(I->pf_banned) continue;

      CEBUGF("#1\n");
      int32 linkedwithid=I->linked_with;

      CEBUGF("#2\n");
      // don't bother looking at if read linked to is already used
      if(lr_used_ids[linkedwithid] !=0) continue;

      CEBUGF("#3\n");

      bool isstrongmatch=false;

      int32 otherread_partnerid=PAF_readpool->getRead(linkedwithid).getTemplatePartnerID();

      // if the read has no template partner id or
      //  template partner is not in contig,, insist that
      //  overlap must have
      //   - either strong good flag (ol_weakgood sufficient at the moment)
      //   - or, if matches below avg freq are allowed and there is no
      //     strong match, no rept flag.
      //     (to prevent sequencing errors in rept regions to generate
      //     false matches)
      //
      // TODO: what about troublemakers? (PAF_istroublemaker)

      if(wantfreqcheck){
	if(otherread_partnerid <0
	   || !PAF_used_ids_in_this_run[otherread_partnerid]){

	  if(I->ol_weakgood) isstrongmatch=true;

	  if(!isstrongmatch){
	    if(allowbelowavgfreq){
	      // ol_rept is not good enough for checking
	      //   if(!I->ol_belowavgfreq || I->ol_rept) continue;
	      // this allows "assembly at sequencing error sites"
	      //  (r = repetitive, b = below avg because of error, * = error)
	      //
	      //  s0    rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	      //  s1    rrrrrrrrrrrrrrbbbbbb*bbbbbbrrrrrrrrrr
	      //  s2            rrrrrrbbbbbb*bbbbbbrrrrrrrrrrrrrrrrrrrrr
	      //
	      // therefore, if allowing joins in below avg regions, the reads
	      //  themselves may not have any rept

	      // BaCh 15.05.2011: some short overlaps have an avg freq (and not
	      //  below avg. freq). Therefore take decision on both ol_belowavgfreq
	      //  or ol_norept flags

	      if(!(I->ol_belowavgfreq || I->ol_norept)
		 || PAF_readpool->getRead(I->rid1).hasFreqRept()
		 || PAF_readpool->getRead(linkedwithid).hasFreqRept()) continue;

	    }else {
	      continue;
	    }
	  }
	}
      }

      CEBUGF("#4\n");

      foundcandidate=true;

      // has the potential new read a partner in the contig?
      // and is it allowed for quick overlap searches?
      if(otherread_partnerid>=0
	 && PAF_used_ids_in_this_run[otherread_partnerid]
	 && I->pf_allowquickoverlap) {
	CEBUGF("#4a\n");

	if(I->best_weight > bestread.weight) {
	  bestreadispaired=true;
	  bestreadI=I;
	  bestread.refid=readid;
	  bestread.newid=linkedwithid;
	  bestread.direction_newid=I->direction;
	  bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	  bestread.weight=I->best_weight;
	  peweighttobeat=I->best_weight + (I->best_weight>>3);
	  returnvalue=isstrongmatch;
	  if(allowspeedypair){
	    // yes, excellent, take it immediately by breaking from the for loop
	    CEBUGF("#4a1\n");
	    speedypair=true;
	    break;
	  }
	}
      } else {
	// nope, well, still good,
	CEBUGF("#4b\n");
	if(bestreadispaired){
	  CEBUGF("#4b1\n");
	  if(I->best_weight > peweighttobeat) {
	    CEBUGF("#4b2\n");
	    if(I->best_weight > bestread.weight) {
	      bestreadispaired=false;
	      bestreadI=I;
	      bestread.refid=readid;
	      bestread.newid=linkedwithid;
	      bestread.direction_newid=I->direction;
	      bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	      bestread.weight=I->best_weight;
	      returnvalue=isstrongmatch;
	      CEBUGF("#4b3\n");
	    }
	  }
	}else{
	  CEBUGF("#4c\n");
	  // special rule: if fully contained read with 100% match
	  //  take that one immediately and break the for loop
	  // set the speedymap flag so that quickrule down below also knows
	  //  about this (rationale: short solexas may be shorter than
	  //  the quickrule lengths)
	  if(PAF_readpool->getRead(linkedwithid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	     && !PAF_readpool->getRead(readid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	     && (*PAF_adsfacts)[I->adsfindex].getScoreRatio()==100
	     && (*PAF_adsfacts)[I->adsfindex].getOverlapLen()== PAF_readpool->getRead(linkedwithid).getLenClippedSeq()){
	    CEBUGF("#4d\n");
	    bestreadispaired=false;
	    bestreadI=I;
	    bestread.refid=readid;
	    bestread.newid=linkedwithid;
	    bestread.direction_newid=I->direction;
	    bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    bestread.weight=I->best_weight;
	    returnvalue=isstrongmatch;
	    speedymap=true;
	    break;
	  } else if(I->best_weight > bestread.weight) {
	    CEBUGF("#4e\n");
	    bestreadispaired=false;
	    bestreadI=I;
	    bestread.refid=readid;
	    bestread.newid=linkedwithid;
	    bestread.direction_newid=I->direction;
	    bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    bestread.weight=I->best_weight;
	    returnvalue=isstrongmatch;
	  }
	}
      }
    }

    CEBUGF("foundcandidate: " << foundcandidate << ' ' << bestread << endl);

    if((bestreadI!=PAF_forward_edges->end() && bestreadI->pf_allowquickoverlap)
       && (speedymap || speedypair)){
#ifndef PUBLICQUIET
      if(speedymap) cout << " speedymap";
      if(speedypair) cout << " speedypair";
#endif
      rindex=-1;
      foundquick=true;
      break;
    }

    if(foundcandidate){

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

      // If wished, we'll use quick rules (premature stop of searching)
      //
      if((bestreadI!=PAF_forward_edges->end() && bestreadI->pf_allowquickoverlap)
	 && (usequickruleSHORTREADS || usequickrule454 || usequickruleION || usequickruleSAN || usequickrulePBSHQ || usequickrulePBSLQ || usequickruleTXT)
	 && bestread.weight > 0){


	bool mayusequick=false;
	uint8 seqtypenewread=PAF_readpool->getRead(bestread.newid).getSequencingType();
	switch(seqtypenewread){
	case ReadGroupLib::SEQTYPE_SANGER : {
	  if(usequickruleSAN) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_454GS20  : {
	  if(usequickrule454) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_IONTORRENT : {
	  if(usequickruleION) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_SOLEXA :
	case ReadGroupLib::SEQTYPE_ABISOLID : {
	  if(usequickruleSHORTREADS) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_PACBIOHQ : {
	  if(usequickrulePBSHQ) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_PACBIOLQ : {
	  if(usequickrulePBSLQ) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_TEXT: {
	  if(usequickruleTXT) mayusequick=true;
	  break;
	}
	default : {
	  MIRANOTIFY(Notify::INTERNAL, "Oooops? seqtypenewread not handled by switch?" << static_cast<uint16>(seqtypenewread) << endl);
	}
	}

	if(mayusequick && n4_checkQuickRules(bestread.ads_node,bestread.newid)){
	  rindex=-1;
	  foundquick=true;
	  break;
	}

      }

    } else {
      // !foundcandidate
      // basically, this read id has become "cork"
      // remove it from the cambium ids by moving the last
      //  id in the vector to this place (can be the same!)
      //  and then reduce the vector size by 1

      cambiumidsincontig[rindex]=cambiumidsincontig[cambiumidsincontig.size()-1];
      cambiumidsincontig.resize(cambiumidsincontig.size()-1);
      CEBUGF("Cambiumcleanup " << cambiumidsincontig.size() << '\n');
    }


    // the "used time" foolguard for extremely high coverage areas
    //  will not search for optimum then
    if(PAF_pafparams.paf_use_emergency_search_stop){
      struct tms mytms;
      times(&mytms);
      clock_t actclocks=mytms.tms_utime+mytms.tms_stime;
      //cout << "Actclocks: " << actclocks << "\tmaxallowedclocks: " << maxallowedclocks << endl;
      if(actclocks>maxallowedclocks){
	// the following should not be necessary, but just to be sure
	// (had a case where the disk got full and MIRA was stuck in
	//  endless loop ... calling times()
	// so I set the loop ending condition manually
	rindex=-1;

#ifndef PUBLICQUIET
	cout << " max time, ess activated ";
#else
	cout << "*";
#endif
	cout.flush();

	break;
      }
    }
  }

  if(bestread.weight>0){
    resultread=bestread;

    resultread.foundmethod=FOUND_STAGE1;
  }else{
    // not really needed, the resultread itself is invalid anyway
    returnvalue=false;
  }

#ifndef PUBLICQUIET
  cout << "\nnrc_fngo: " << numreadschecked << "\tciic.size(): " << cambiumidsincontig.size();
#endif

  // cambiumshuffling
  // if too many reads were checked but the overlap was chosen via quickrule,
  //  we have too much cruft at the end of the cambium
  // move that cruft to the beginning of the cambium
  //  but not if the number of additional ids in cambium is less than 10
  if(foundquick
     && numreadschecked >500
     && numreadschecked < static_cast<int64>(cambiumidsincontig.size())-10) {
    numreadschecked--;

    PAF_tmpcambiumids.clear();
    PAF_tmpcambiumids.reserve(cambiumidsincontig.capacity());
    PAF_tmpcambiumids.insert(PAF_tmpcambiumids.end(),
		      cambiumidsincontig.end()-numreadschecked,
		      cambiumidsincontig.end());
    PAF_tmpcambiumids.insert(PAF_tmpcambiumids.end(),cambiumidsincontig.begin(),
		      cambiumidsincontig.end()-numreadschecked);
    PAF_tmpcambiumids.swap(cambiumidsincontig);
    PAF_tmpcambiumids.clear();
    numreadschecked++;

#ifndef PUBLICQUIET
    cout << "\tcambiumshuffling";
#endif
  }

#ifndef PUBLICQUIET
  cout << '\n';
#endif

  FUNCEND();

  return returnvalue;
}
//#define CEBUGF(bla)




/*************************************************************************
 *
 * allow only overlaps where the newid has a template partner in the
 *  contig
 *
 *************************************************************************/

//#define CEBUGF(bla)   {cout << bla; cout.flush(); }
void Pathfinder::n4_findNextPillarCounterpart(vector<int32> & cambiumidsincontig, nextreadtoadd_t & resultread, clock_t & maxallowedclocks)
{
  FUNCSTART("void Pathfinder::n4_findNextStrongGoodOverlap(const vector<int32> & cambiumidsincontig, nextreadtoadd_t & resultread, clock_t & maxallowedclocks) ");

  // lr_ == local reference
  const vector<int8> &  lr_used_ids = *PAF_used_ids_ptr;
  const vector<uint8> & lr_istroublemaker = *PAF_istroublemaker_ptr;

  if(cambiumidsincontig.size()==0) {
    return;
  }

  nextreadtoadd_t bestread=resultread;
  bool bestreadispaired=false;
  uint32 peweighttobeat=0;


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif



  bool   usequickruleSAN=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrule454=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_use_quick_rule;
  bool   usequickruleION=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrulePBSHQ=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_PACBIOHQ].getPathfinderParams().paf_use_quick_rule;
  bool   usequickrulePBSLQ=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_PACBIOLQ].getPathfinderParams().paf_use_quick_rule;
  bool   usequickruleTXT=(*PAF_miraparams)[ReadGroupLib::SEQTYPE_TEXT].getPathfinderParams().paf_use_quick_rule;

  bool   usequickruleSHORTREADS=true;

  bool foundquick=false;

  // going backwards through the contglist should be a tick faster than forward
  //  also helps to retain a bit of flexibility when emergency search stops occur
  int32 rindex=cambiumidsincontig.size()-1;

  uint32 numreadschecked=0;
  for(; rindex >=0 ; --rindex, ++numreadschecked) {
    int32 readid=cambiumidsincontig[rindex];

    // foundcandidate tells whether for this read id there are overlaps
    //  to unassembled reads
    // if not, the this read id is "dead" and is removed from the
    //  cambiumidsincontig further down (it becomes "cork" :-)
    bool foundcandidate=false;

    CEBUGF("fnpC: new read " << readid << '\n');

    vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[readid];

    for(;I!=PAF_forward_edges->end() && I->rid1 == readid;I++){

      // don't bother looking at if overlap is banned
      if(I->pf_banned) continue;

      CEBUGF("#1\n");
      int32 linkedwithid=I->linked_with;

      CEBUGF("#2\n");
      // don't bother looking at if read linked to is already used
      if(lr_used_ids[linkedwithid] !=0) continue;

      CEBUGF("#3\n");

      int32 otherread_partnerid=PAF_readpool->getRead(linkedwithid).getTemplatePartnerID();

      // don't bother if the read has no template partner id or
      //  template partner is not in contig,
      //
      // TODO: what about troublemakers? (PAF_istroublemaker)
      if(otherread_partnerid <0) continue;
      CEBUGF("#3a\n");
      if(!PAF_used_ids_in_this_run[otherread_partnerid]) continue;

      CEBUGF("#4\n");

      foundcandidate=true;

      // excellent, check whether it's the best overlap we have
      if(I->best_weight > bestread.weight) {
	CEBUGF("#6\n");

	bestreadispaired=true;
	bestread.refid=readid;
	bestread.newid=linkedwithid;
	bestread.direction_newid=I->direction;
	bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	bestread.weight=I->best_weight;
	peweighttobeat=I->best_weight + (I->best_weight>>3);
      }
    }

    // If wished, we'll use quick rules (premature stop of searching)
    //
    if(foundcandidate){

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

      if((usequickruleSHORTREADS || usequickrule454 || usequickruleION || usequickrulePBSHQ || usequickrulePBSLQ || usequickruleSAN || usequickruleTXT)
	 && bestread.weight > 0){


	bool mayusequick=false;
	uint8 seqtypenewread=PAF_readpool->getRead(bestread.newid).getSequencingType();
	switch(seqtypenewread){
	case ReadGroupLib::SEQTYPE_SANGER : {
	  if(usequickruleSAN) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_454GS20  : {
	  if(usequickrule454) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_IONTORRENT  : {
	  if(usequickruleION) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_SOLEXA :
	case ReadGroupLib::SEQTYPE_ABISOLID : {
	  if(usequickruleSHORTREADS) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_PACBIOHQ : {
	  if(usequickrulePBSHQ) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_PACBIOLQ : {
	  if(usequickrulePBSLQ) mayusequick=true;
	  break;
	}
	case ReadGroupLib::SEQTYPE_TEXT : {
	  if(usequickruleTXT) mayusequick=true;
	  break;
	}
	default : {
	  MIRANOTIFY(Notify::INTERNAL, "Oooops? seqtypenewread not handled by switch?" << static_cast<uint16>(seqtypenewread) << endl);
	}
	}

	if(mayusequick && n4_checkQuickRules(bestread.ads_node,bestread.newid)){
	  rindex=-1;
	  foundquick=true;
	  break;
	}


	// the "used time" foolguard for extremely high coverage areas
	//  will not search for optimum then
	if(PAF_pafparams.paf_use_emergency_search_stop){
	  struct tms mytms;
	  times(&mytms);
	  clock_t actclocks=mytms.tms_utime+mytms.tms_stime;
	  //cout << "Actclocks: " << actclocks << "\tmaxallowedclocks: " << maxallowedclocks << endl;
	  if(actclocks>maxallowedclocks){
	    // the following should not be necessary, but just to be sure
	    // (had a case where the disk got full and MIRA was stuck in
	    //  endless loop ... calling times()
	    // so I set the loop ending condition manually
	    rindex=-1;

#ifndef PUBLICQUIET
	    cout << " max time, ess activated ";
#else
	    cout << "*";
#endif
	    cout.flush();

	    break;
	  }
	}
      }
    }else{
      // !foundcandidate
      // basically, this read id has become "cork"
      // remove it from the cambium ids by moving the last
      //  id in the vector to this place (can be the same!)
      //  and then reduce the vector size by 1

      cambiumidsincontig[rindex]=cambiumidsincontig[cambiumidsincontig.size()-1];
      cambiumidsincontig.resize(cambiumidsincontig.size()-1);
    }
  }

  if(bestread.weight>0){
    resultread=bestread;

    resultread.foundmethod=FOUND_STAGE2;
  }

#ifndef PUBLICQUIET
  cout << "\nnrc_fnpc: " << numreadschecked << '\n';
#endif

  FUNCEND();

  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/



//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Pathfinder::n4_setupAllowedReferences(vector<int32> & allowed_references, vector<int32> & ids_in_contig)
{
  // temporary test, probably does not harm
  allowed_references=ids_in_contig;
}
//#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUGF(bla)   {cout << bla; cout.flush(); }

void Pathfinder::n4_findNextAllowedRefOverlap(vector<int32> & cambiumidsincontig, nextreadtoadd_t & resultread, clock_t & maxallowedclocks)
{
  FUNCSTART("");

  if(cambiumidsincontig.size()==0) {
    return;
  }

  bool usequickruleSHORTREADS=true;
  bool foundquick=false;


  // lr_ == local reference
  const vector<int8> &  lr_used_ids = *PAF_used_ids_ptr;

  nextreadtoadd_t bestread=resultread;
  bool bestreadispaired=false;
  uint32 peweighttobeat=0;


  int32 rindex=cambiumidsincontig.size()-1;

  uint32 numreadschecked=0;
  for(; rindex >=0 ; --rindex, ++numreadschecked) {
    int32 readid=cambiumidsincontig[rindex];

    // foundcandidate tells whether for this read id there are overlaps
    //  to unassembled reads
    // if not, the this read id is "dead" and is removed from the
    //  cambiumidsincontig further down (it becomes "cork" :-)
    bool foundcandidate=false;


    CEBUGF("fnarO: new read " << readid << '\n');
    CEBUGF("rindex " << rindex << '\n');

    vector<newedges_t>::const_iterator I=(*PAF_lowerbound_oedges_ptr)[readid];

    for(;I!=PAF_forward_edges->end() && I->rid1 == readid;I++){

      // don't bother looking at if overlap is banned
      //  or even whole linked read is banned
      if(I->pf_banned
	 || PAF_blacklisted_ids[I->linked_with]==0) continue;

      CEBUGF("#1\n");
      int32 linkedwithid=I->linked_with;

      CEBUGF("#2\n");
      // don't bother looking at if read linked to is already used
      if(lr_used_ids[linkedwithid] !=0) continue;

      CEBUGF("#3\n");

      int32 otherread_partnerid=PAF_readpool->getRead(linkedwithid).getTemplatePartnerID();

      // if the read has template partner id but
      //  template partner is not in contig, don't take
      //
      // TODO: what about troublemakers? (PAF_istroublemaker)
      if(otherread_partnerid >0
	 && !PAF_used_ids_in_this_run[otherread_partnerid]) continue;

      CEBUGF("#4\n");

      foundcandidate=true;

      // has the potential new read a partner in the contig?
      if(otherread_partnerid>=0
	 && PAF_used_ids_in_this_run[otherread_partnerid]) {
	// yes, excellent
	if(I->best_weight > bestread.weight) {
	  bestreadispaired=true;
	  bestread.refid=readid;
	  bestread.newid=linkedwithid;
	  bestread.direction_newid=I->direction;
	  bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	  bestread.weight=I->best_weight;
	  peweighttobeat=I->best_weight + (I->best_weight>>3);
	}
      } else {
	// nope, well, still good,
	if(bestreadispaired){
	  if(I->best_weight > peweighttobeat) {
	    if(I->best_weight > bestread.weight) {
	      bestreadispaired=false;
	      bestread.refid=readid;
	      bestread.newid=linkedwithid;
	      bestread.direction_newid=I->direction;
	      bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	      bestread.weight=I->best_weight;
	    }
	  }
	}else{
	  if(I->best_weight > bestread.weight) {
	    bestreadispaired=false;
	    bestread.refid=readid;
	    bestread.newid=linkedwithid;
	    bestread.direction_newid=I->direction;
	    bestread.ads_node=&(*PAF_adsfacts)[I->adsfindex];
	    bestread.weight=I->best_weight;
	  }
	}
      }
    }

    if(!foundcandidate){
      // basically, this read id has become "cork"
      // remove it from the cambium ids by moving the last
      //  id in the vector to this place (can be the same!)
      //  and then reduce the vector size by 1

      cambiumidsincontig[rindex]=cambiumidsincontig[cambiumidsincontig.size()-1];
      cambiumidsincontig.resize(cambiumidsincontig.size()-1);
    }else if(usequickruleSHORTREADS){
      // else, see if we can shorten things a bit for 454 & microreads

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      // TODO: what about PacBio ??? Not in yet
      if(PAF_readpool->getRead(bestread.newid).getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){

	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_SOLEXA].getPathfinderParams().paf_quickrule_minlen1
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_SOLEXA].getPathfinderParams().paf_quickrule_minsim1){
#ifndef PUBLICQUIET
	  cout << " Solexa quick rule 1";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}

	//  if they have an overlap >80 and expected score ratio >=90%,
	//  we just stop the search if they are NONMULTICOPY
	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_SOLEXA].getPathfinderParams().paf_quickrule_minlen2
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_SOLEXA].getPathfinderParams().paf_quickrule_minsim2){
#ifndef PUBLICQUIET
	  cout << " Solexa quick rule 2";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}
      }else if(PAF_readpool->getRead(bestread.newid).getSequencingType()==ReadGroupLib::SEQTYPE_454GS20){

	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_quickrule_minlen1
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_quickrule_minsim1){
#ifndef PUBLICQUIET
	  cout << " 454 quick rule 1";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}

	//  if they have an overlap >80 and expected score ratio >=90%,
	//  we just stop the search if they are NONMULTICOPY
	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_quickrule_minlen2
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getPathfinderParams().paf_quickrule_minsim2){
#ifndef PUBLICQUIET
	  cout << " 454 quick rule 2";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}
      }else if(PAF_readpool->getRead(bestread.newid).getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){

	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_quickrule_minlen1
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_quickrule_minsim1){
#ifndef PUBLICQUIET
	  cout << " ION quick rule 1";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}

	//  if they have an overlap >80 and expected score ratio >=90%,
	//  we just stop the search if they are NONMULTICOPY
	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_quickrule_minlen2
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getPathfinderParams().paf_quickrule_minsim2){
#ifndef PUBLICQUIET
	  cout << " ION quick rule 2";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}
      }else if(PAF_readpool->getRead(bestread.newid).getSequencingType()==ReadGroupLib::SEQTYPE_ABISOLID){

	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_ABISOLID].getPathfinderParams().paf_quickrule_minlen1
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_ABISOLID].getPathfinderParams().paf_quickrule_minsim1){
#ifndef PUBLICQUIET
	  cout << " ABISOLID quick rule 1";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}

	//  if they have an overlap >80 and expected score ratio >=90%,
	//  we just stop the search if they are NONMULTICOPY
	if(bestread.ads_node->getOverlapLen() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_ABISOLID].getPathfinderParams().paf_quickrule_minlen2
	   && bestread.ads_node->getScoreRatio() >= (*PAF_miraparams)[ReadGroupLib::SEQTYPE_ABISOLID].getPathfinderParams().paf_quickrule_minsim2){
#ifndef PUBLICQUIET
	  cout << " ABISOLID quick rule 2";
#endif
	  rindex=-1;
	  foundquick=true;
	  break;
	}
      }
    }

    // the "used time" foolguard for extremely high coverage areas
    //  will not search for optimum then
    if(PAF_pafparams.paf_use_emergency_search_stop){
      struct tms mytms;
      times(&mytms);
      clock_t actclocks=mytms.tms_utime+mytms.tms_stime;
      //cout << "Actclocks: " << actclocks << "\tmaxallowedclocks: " << maxallowedclocks << endl;
      if(actclocks>maxallowedclocks){
	// the following should not be necessary, but just to be sure
	// (had a case where the disk got full and MIRA was stuck in
	//  endless loop ... calling times()
	// so I set the loop ending condition manually
	rindex=-1;

#ifndef PUBLICQUIET
	cout << " max time, ess activated ";
#else
	cout << "*";
#endif
	cout.flush();

	break;
      }
    }


  }

  if(bestread.weight>0){
    resultread=bestread;

    resultread.foundmethod=FOUND_STAGE100;
  }

#ifndef PUBLICQUIET
  cout << "\nnrc_fnro: " << numreadschecked << '\n';
#endif

  FUNCEND();

  return;
}
//#define CEBUGF(bla)   {cout << bla; cout.flush(); }




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Pathfinder::n4_buildESTContigStartWithStrongGood(vector<Align> & aligncache, int32 startid, Contig & con)
{
  FUNCSTART("void Pathfinder::n4_buildESTContigStartWithStrongGood(vector<Align> & aligncache, int32 startid, Contig & con)");

#ifndef PUBLICQUIET
  cout << "n4_buildESTContigStartWithStrongGood\n";
#endif

  CEBUG("bc 1\n");

  vector<int8> & PAF_used_ids = *PAF_used_ids_ptr;

  Contig::errorstatus_t contigerrstat;
  // this should speed up things a bit, avoiding copying during extension
  //  although the 10000 might be not enough for really pathological cases
  //  with many SRMB regions and huge coverage
  contigerrstat.reads_affected.reserve(10000);

  // this is a kind of duplicate of PAF_used_ids_in_this_run (which has
  //  the size of the readpool)
  // However, this array holds only those ids that have been built into
  //  the contig in this call, giving a faster access to the info
  vector<int32> & ids_in_contig = *PAF_ids_in_contig_ptr;
  ids_in_contig.clear();
  ids_in_contig.reserve(PAF_used_ids.size());

  // "duplicate" of ids_in_contig, but gets cleared after each
  //  stage1, stage2, ..., stageN round below.
  // Even faster access to the really newest additions to the contig
  //vector<int32> PAF_fresh_ids_in_contig;
  PAF_fresh_ids_in_contig.clear();
  PAF_fresh_ids_in_contig.reserve(PAF_used_ids.size());

  // "duplicate" of ids_in_contig, but contains ids that have overlaps
  //  with unassembled reads and therefore can contribute to growth
  // Once IDs have no open overlap anymore, they become dead ("cork")
  //  and are removed from the cambium vector (by the n4_find* functions)
  // On each change of stage, cambium_ids gets refilled with
  //  ids_in_contig to be sure that nothing gets lost
  //vector<int32> & PAF_cambium_ids;
  PAF_cambium_ids.clear();
  PAF_cambium_ids.reserve(PAF_used_ids.size());


  PAF_buildcontig_newlinecounter=0;


  CEBUG("bc 3a 1\n");
  if(startid<0) {
    cerr << "Startid: " << startid << endl;
    MIRANOTIFY(Notify::INTERNAL, "Starting with non-existing read?");
  }

  if(startid >= static_cast<int32>(PAF_used_ids.size())) {
    cerr << "Startid: " << startid << endl;
    MIRANOTIFY(Notify::INTERNAL, "Starting with read id greater than number iof reads?\n");
  }

  if(PAF_used_ids[startid]) {
    cerr << "Startid: " << startid << endl;
    MIRANOTIFY(Notify::INTERNAL, "Starting with read that is already used?\n");
  }

  CEBUG("Startid: " << startid << endl);

#ifndef PUBLICQUIET
  cout << endl;
#endif

  ++PAF_readaddattempts;
  con.addRead(aligncache,
	      nullptr, startid, startid, 1,
	      (*PAF_multicopies_ptr)[startid],
	      0,
	      contigerrstat);

  CEBUG("bc 3a 2\n");

  PAF_used_ids[startid]=1;
  PAF_used_ids_in_this_run[startid]=1;
  ids_in_contig.push_back(startid);
  PAF_fresh_ids_in_contig.push_back(startid);
  PAF_cambium_ids.push_back(startid);

#ifndef PUBLICQUIET
  cout << endl;
#else
  cout << "+";
#endif

  struct tms mytms;
  times(&mytms);
  clock_t baseclocks=mytms.tms_utime+mytms.tms_stime;
  clock_t actclocks=baseclocks;
  clock_t maxallowedclocks=baseclocks+PAF_pafparams.paf_maxcontigclockticks;

#ifdef CLOCK_STEPS1
  vector<suseconds_t> us_pfsearch;
  vector<suseconds_t> us_conadd;
  vector<suseconds_t> us_overlapban;
#endif

  nextreadtoadd_t nrta;
  nrta.newid=0;

  uint32 statscounter=1;

  // the forcegrow for addRead()
  //  will be set to 0 for adding NMC, 10-20? for MC
  int32 forcegrow=0;

  // stage 1 : add best overlaps
  // stage 2 : add paired-ends only

  uint32 stage=1;
  uint32 readsaddedinstages=0;

  /* forcenewstage is set by a stage if it explicitly wants to change
     to a new stage after the current found overlap has been added
     value of forcenewstage is the target stage */
  uint32 forcenewstage=0;

  while(stage < 1000) {
#ifndef PUBLICQUIET
    if(statscounter == 25) {
      cout << "rstatsm: len of c : " << con.getContigLength();
      cout << "\nrstatsm: ids in c : " << ids_in_contig.size();
      cout << "\nrstatsm: cids in c : " << PAF_cambium_ids.size();
      cout << "\nrstatsm: fids in c : " << PAF_fresh_ids_in_contig.size();
      cout << "\nrstatsi: banned o : " << PAF_overlapsbanned;
#ifdef CLOCK_STEPS1
      cout << "\nrstat timings: " << ids_in_contig.size();
      cout << " : " << avg_suseconds(us_pfsearch);
      cout << " / " << avg_suseconds(us_conadd);
      cout << " / " << avg_suseconds(us_overlapban);
#endif
      cout << '\n';
      statscounter=0;
    }
#endif


    nrta.refid=-1;
    nrta.newid=-1;
    nrta.weight=0;
    nrta.direction_newid=0;
    nrta.ads_node=nullptr;

#ifdef CLOCK_STEPS1
    timeval us_start;
    gettimeofday(&us_start,nullptr);
#endif

#ifndef PUBLICQUIET
    cout << "Searching";
#endif

    switch(stage){
    case 1 : {
#ifndef PUBLICQUIET
      cout << " (stage1 overlap) ...";
      cout.flush();
#endif
      forcegrow=0;
      n4_findNextGoodOverlap(false,true,false,PAF_cambium_ids, nrta, maxallowedclocks);

      break;
    }
    case 2 : {
#ifndef PUBLICQUIET
      cout << " (stage2 overlap) ...";
      cout.flush();
#endif
      forcegrow=0;
      n4_findNextPillarCounterpart(PAF_cambium_ids, nrta, maxallowedclocks);

      break;
    }
    case 100 : {
      // map to existing, allowing only overlaps with data
      //  that has a frequency >3
#ifndef PUBLICQUIET
      cout << " (stage100 overlap) ...";
      cout.flush();
#endif
      // disallow the contig to grow!
      forcegrow=-1;

      // At the moment we use the cambium ids anyway as allowed references
      //if(PAF_cambium_ids.empty()){
      //	n4_setupAllowedReferences(PAF_cambium_ids,ids_in_contig);
      //}

      n4_findNextAllowedRefOverlap(PAF_cambium_ids, nrta, maxallowedclocks);

      // nothing found? the clear allowed references
      if(nrta.newid < 0) {
	PAF_cambium_ids.clear();
      }
      break;
    }
    default : {
      cerr << "\nStage " << stage << '\n';
      MIRANOTIFY(Notify::INTERNAL, "stage not 1,2,3 ... this should not happen.\n");
    }
    }
#ifndef PUBLICQUIET
    cout << " done.\n";
    cout.flush();
#endif

#ifdef CLOCK_STEPS1
    us_pfsearch.push_back(diffsuseconds(us_start));
#endif


    if(nrta.newid < 0) {
      // nothing found
      // go to next stage and if we went through all stages,
      //  go back to stage 1 only if there had been a read added
      //  in the three stages. if not, this contigs is finished
#ifndef PUBLICQUIET
      cout << "\nWe were in stage " << stage << ", switching to stage " << stage+1 << endl;
      cout << "Reads added: " << readsaddedinstages << endl;
#endif
      ++stage;
      if(stage==3){
	if(readsaddedinstages){
	  // we added some reads in the stages. Perfect, restart
	  //  at stage 1 to see whether we could get further.
	  stage=1;
	  readsaddedinstages=0;
	}else{
	  // nope, so start the beef up loop now
	  stage=100;

	  uint32 numreadspresent=PAF_used_ids_ptr->size();
	  // is now also needed by n4 for stage 100 mapping, but needs to get
	  //  initialised a bit differently (larger number). Currently done here
	  //  switching to stage 100, can be done in basicSetup once old
	  //  routines disappear
	  PAF_blacklisted_ids.clear();
	  PAF_blacklisted_ids.resize(PAF_used_ids_ptr->size(),5);
	}
	PAF_fresh_ids_in_contig.clear();
      }else if(stage>=100){
	//stop the loop now
	stage=1000;
      }
      PAF_cambium_ids=ids_in_contig;
    }else{
#ifndef PUBLICQUIET
      cout << nrta;
#endif

#ifdef CLOCK_STEPS1
      gettimeofday(&us_start,nullptr);
#endif

      statscounter++;
      ++PAF_readaddattempts;
      con.addRead(aligncache,
		  nrta.ads_node, nrta.refid, nrta.newid, nrta.direction_newid,
		  (*PAF_multicopies_ptr)[nrta.newid],
		  forcegrow,
		  contigerrstat);
#ifdef CLOCK_STEPS1
      us_conadd.push_back(diffsuseconds(us_start));
#endif

      cout.flush();

      //Contig::setCoutType(Contig::AS_TEXT);
      //cout << ",,,,,,\n" << con << "\n\n";


      if(PAF_buildcontig_newlinecounter==0){
	cout << "[" << ids_in_contig.size() << "] "; cout.flush();
      }
      if(contigerrstat.code == Contig::ENOERROR) {
	ids_in_contig.push_back(nrta.newid);
	PAF_fresh_ids_in_contig.push_back(nrta.newid);
	PAF_cambium_ids.push_back(nrta.newid);
	++readsaddedinstages;
	PAF_used_ids[nrta.newid]=1;
	PAF_used_ids_in_this_run[nrta.newid]=1;

#ifndef PUBLICQUIET
	cout << "\t+\n";
#else
	cout << "+";
#endif
      }else{
#ifdef CLOCK_STEPS1
	gettimeofday(&us_start,nullptr);
#endif
	n4_handleReadNotAligned(contigerrstat,nrta);

	// handle blacklisting of failed alignment when refid was a short read
	//  and newid a long one
	if(stage==100
	   && PAF_readpool->getRead(nrta.refid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	   && !PAF_readpool->getRead(nrta.newid).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	  if(PAF_blacklisted_ids[nrta.newid]) PAF_blacklisted_ids[nrta.newid]--;
	}
#ifdef CLOCK_STEPS1
	us_overlapban.push_back(diffsuseconds(us_start));
#endif
      }

      if(forcenewstage>0){
#ifndef PUBLICQUIET
	cout << "\nWe were in stage " << stage << ", being forced to stage " << stage+1 << endl;
	cout << "Reads added: " << readsaddedinstages << endl;
#endif
	readsaddedinstages=0;
	stage=forcenewstage;
	forcenewstage=0;
	PAF_fresh_ids_in_contig.clear();
	PAF_cambium_ids=ids_in_contig;
      }

#ifndef PUBLICQUIET
#else
      PAF_buildcontig_newlinecounter++;
      if(PAF_buildcontig_newlinecounter==60){
	PAF_buildcontig_newlinecounter=0;
	cout << "   " << con.getContigLength();
#ifdef CLOCK_STEPS1
	cout << "\t" << avg_suseconds(us_pfsearch);
	cout << " / " << avg_suseconds(us_conadd);
	cout << " / " << avg_suseconds(us_overlapban);
#endif
	cout << endl;
      }
#endif


    }

    //if(con.getNumReadsInContig()==1198){
    //  cout << "\npathfinder debugstop.\n";
    //  assout::saveAsCAF(con, "test.caf", true);
    //  abort();
    //}

    if(PAF_pafparams.paf_use_max_contig_buildtime){
      times(&mytms);
      actclocks=mytms.tms_utime+mytms.tms_stime;
      if(actclocks>maxallowedclocks){
	// the following should not be necessary, but just to be sure
	// (had a case where the disk got full and MIRA was stuck in
	//  endless loop ... calling times()
	// so I set the loop ending condition manually
	nrta.newid=-1;

	cout << "\nMaximum build time for this contig reached, aborting build.\n";
	FUNCEND();
	return;
      }
    }
  }

  FUNCEND();
}

//#define CEBUG(bla)
