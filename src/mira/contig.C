/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 * Copyright (C) 2000 and later by Bastien Chevreux
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


#include "contig.H"

#include "util/misc.H"
#include "assembly_output.H"

#include "mira/align.H"
#include "mira/ads.H"
#include "mira/readpool.H"

#include "util/stlimprove.H"

#include <unordered_set>


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)


//#define PARANOIABUGTRACKFLAG
#ifdef PARANOIABUGTRACKFLAG
#define paranoiaBUGSTAT(statement) { statement;}
#define paranoiaBUGIF(ifcond, statement) { if(ifcond) {statement;}}
#else
#define paranoiaBUGSTAT(statement)
#define paranoiaBUGIF(ifcond, statement)
#endif



enum { CON_CONSTARTSIZE=4000,
       CON_VECTORINCSIZE=4000,
       CON_READSSTARTSIZE=100,
       CON_READSINCSIZE=100};

bool Contig::CON_mastercebugflag=false;

bool Contig::CON_abortflag=false;

uint32 Contig::CON_id_counter=1;

uint32 Contig::CON_railcounter=0;

uint32 Contig::CON_cer_idcounter=0;





#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
const Contig::consensus_counts_t Contig::CON_concounts_zero=
{0,0,0,0,0,0,0, {0,0,0,0,0,0,0,0}, 0,0,0,
 0,{0},{0},{0},{0},{0},'@','@',0};
const Contig::consensus_counts_t Contig::CON_concounts_zero_nobb=
{0,0,0,0,0,0,0, {0,0,0,0,0,0,0,0}, 0,0,0};


uint8 Contig::CON_outtype=AS_TEXT;


//bool  Contig::CON_outputRails=true;
bool  Contig::CON_outputrails=false;


bool Contig::CON_static_ok=false;


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::setCEBUGFlag(int32 id1, int32 id2)
{
  CON_cebugflag=false;
  return;

  if(!CON_mastercebugflag) return;
  if(CON_reads.size()>=9736
     && CON_reads.size()<=18874) CON_cebugflag=true;
  return;

  if(CON_cebugnames.empty()){
//CON_cebugnames.insert("G8BUD:1800:2515");
//...
  }
  if(!CON_cebugnames.empty()){
    if(id1>=0){
      if(CON_cebugnames.find(CON_readpool->getRead(id1).getName())
	 != CON_cebugnames.end()) CON_cebugflag=true;
    }
    if(id2>=0){
      if(CON_cebugnames.find(CON_readpool->getRead(id2).getName())
	 != CON_cebugnames.end()) CON_cebugflag=true;
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


// Plain vanilla constructor
Contig::Contig(std::vector<MIRAParameters> * params, ReadPool & readpool) : CON_reads(readpool)
{
  FUNCSTART("Contig::Contig(MIRAParameters * params, ReadPool & readpool)");

  if(!CON_static_ok){
    // cannot use a real static initialiser because of the static
    //  initialisation fiasco which would happen otherwise with
    //  multitag_t string container
    CON_danger_zones_ids.push_back(multitag_t::newIdentifier("ALUS"));
    CON_danger_zones_ids.push_back(multitag_t::newIdentifier("REPT"));

    CON_baselock_ids.push_back(multitag_t::newIdentifier("SRMr"));
    CON_baselock_ids.push_back(multitag_t::newIdentifier("CRMr"));
    CON_baselock_ids.push_back(multitag_t::newIdentifier("WRMr"));

    CON_snplock_ids.push_back(multitag_t::newIdentifier("SIOr"));
    CON_static_ok=true;
  }

  CON_miraparams=params;

  //CON_aligncache = new Align(params, true);

  CON_readpool=&readpool;

  init();

  CON_id=CON_id_counter++;

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::init()
{
  zeroVars();

  CON_tmpcons_from_backbone=false;
  CON_mergenewsrreads=false;
  CON_hasforcemergeareas=false;

  CON_cebugflag=false;
  CON_verbose=true;     // contigs are a bit more verbose by default
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::discard()
{
  FUNCSTART("Contig::discard()");

  nukeSTLContainer(CON_2tmpcons);

  //CON_aligncache->discard();

  nukeSTLContainer(CON_counts);
  CON_reads.clear();
  nukeSTLContainer(CON_templates_present);
  nukeSTLContainer(CON_consensus_tags);
  nukeSTLContainer(CON_fixedconsseq);
  nukeSTLContainer(CON_fixedconsqual);
  nukeSTLContainer(CON_last_dangerous_overlaps);
  nukeSTLContainer(CON_targetcoverageperst);
//  nukeSTLContainer(CON_cheat_intelcons);
//  nukeSTLContainer(CON_cheat_intelconsqual);
  nukeSTLContainer(CON_allconsseq);
  nukeSTLContainer(CON_allconsqual);
  nukeSTLContainer(CON_alladjustments);
  nukeSTLContainer(CON_strainconsseq);
  nukeSTLContainer(CON_strainconsqual);
  nukeSTLContainer(CON_strainadjustments);
  nukeSTLContainer(CON_2tmpcons);
  nukeSTLContainer(CON_allowedrefids);

  zeroVars();

  FUNCEND();
}


void Contig::zeroVars()
{

  CON_readsperstrain.clear();
  CON_readsperreadgroup.clear();

  CON_markerpositions.clear();
  CON_mpindex_msrkceu_left=-1;
  CON_mpindex_msrkceu_right=-1;

  CON_longestreadseen=0;
  CON_longestrailseen=0;
  CON_longestnonbbreadseen=0;

  CON_conscalc_mincov=0;
  CON_conscalc_minqual=0;
  CON_conscalc_missingchar='X';
  CON_conscalc_assumediploid=false;
  CON_conscalc_allowiupac=true;
  CON_conscalc_addconstag=true;
  zeroStats();

  CON_specialsraddconditions=false;
  CON_ssrc_maxtotalerrors=-1;
  CON_ssrc_maxgaps=-1;
  CON_ssrc_maxmismatches=-1;

  for(uint32 i=0; i<NUMMERGESEQTYPES; i++){
    CON_nummergedreads_perseqtype[i]=0;
  }

  CON_nameprefix="stdname";
  CON_name.clear();
  CON_outtype=AS_TEXT;

  CON_contains_long_repeats_only=false;
  CON_contains_majority_digitallynormalised_reads=0;

  CON_us_steps.clear();
  CON_us_steps.resize(USCLO_END,0);
  CON_us_steps_iric.clear();
  CON_us_steps_iric.resize(USCLOIRIC_END,0);
  CON_us_steps_drfc.clear();
  CON_us_steps_drfc.resize(USCLODRFC_END,0);
  CON_us_steps_cons.clear();
  CON_us_steps_cons.resize(USCLOCONS_END,0);
  CON_track_numins=0;
  CON_track_numdels=0;

  CON_isbackbonecontig=false;

  definalise();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Contig::~Contig()
{
  FUNCSTART("Contig::~Contig()");

  discard();

  //delete CON_aligncache;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::definalise()
{
//  CON_cheat_intelcons.clear();
//  CON_cheat_intelconsqual.clear();

  CON_allconsseq.clear();
  CON_allconsqual.clear();
  CON_alladjustments.clear();
  CON_strainconsseq.clear();
  CON_strainconsqual.clear();
  CON_strainadjustments.clear();

  CON_stats.statsvalid=false;
  CON_finalised=false;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::saveMem()
{
  FUNCSTART("void Contig::saveMem()");

  //CON_aligncache->discard();

  nukeSTLContainer(CON_2tmpcons);

  definalise();
//  nukeSTLContainer(CON_cheat_intelcons);
//  nukeSTLContainer(CON_cheat_intelconsqual);

  nukeSTLContainer(CON_allconsseq);
  nukeSTLContainer(CON_allconsqual);
  nukeSTLContainer(CON_alladjustments);
  nukeSTLContainer(CON_strainconsseq);
  nukeSTLContainer(CON_strainconsqual);
  nukeSTLContainer(CON_strainadjustments);

  nukeSTLContainer(CON_allowedrefids);

  FUNCEND();
  return;
}


// Copy constructor
//  no discard needed as this object will be freshly created when
//  called through this constructor
Contig::Contig(Contig const &other) : CON_reads(*(other.CON_readpool))
{
  FUNCSTART("Contig::Contig(Contig const &other)");

  //CON_aligncache=new Align(CON_miraparams,true);

  init();

  *this=other;                               // call the copy operator

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Copy operator, needed by copy-constructor
Contig const & Contig::operator=(Contig const & other)
{
  FUNCSTART("Contig const & Contig::operator=(Contig const & other)");

  if(this != &other){
    discard();

    CON_cebugflag=other.CON_cebugflag;

    CON_miraparams=other.CON_miraparams;

    CON_id=other.CON_id;
    CON_nameprefix=other.CON_nameprefix;
    CON_name=other.CON_name;
    CON_finalised=other.CON_finalised;

    CON_readpool=other.CON_readpool;
    CON_reads=other.CON_reads;
    CON_counts=other.CON_counts;
    CON_templates_present=other.CON_templates_present;
    CON_consensus_tags=other.CON_consensus_tags;

    CON_targetcoverageperst=other.CON_targetcoverageperst;

    CON_2tmpcons=other.CON_2tmpcons;

    CON_tmpcons_from_backbone=other.CON_tmpcons_from_backbone;
    CON_mergenewsrreads=other.CON_mergenewsrreads;
    CON_hasforcemergeareas=other.CON_hasforcemergeareas;

    CON_specialsraddconditions=other.CON_specialsraddconditions;
    CON_ssrc_maxtotalerrors=other.CON_ssrc_maxtotalerrors;
    CON_ssrc_maxgaps=other.CON_ssrc_maxgaps;
    CON_ssrc_maxmismatches=other.CON_ssrc_maxmismatches;

    for(uint32 i=0; i<NUMMERGESEQTYPES; i++){
      CON_nummergedreads_perseqtype[i]=other.CON_nummergedreads_perseqtype[i];
    }

    CON_longestreadseen=other.CON_longestreadseen;
    CON_longestrailseen=other.CON_longestrailseen;
    CON_longestnonbbreadseen=other.CON_longestnonbbreadseen;

    CON_isbackbonecontig=other.CON_isbackbonecontig;

    CON_readsperstrain=other.CON_readsperstrain;
    CON_readsperreadgroup=other.CON_readsperreadgroup;

    CON_markerpositions=other.CON_markerpositions;
    CON_mpindex_msrkceu_left=other.CON_mpindex_msrkceu_left;
    CON_mpindex_msrkceu_right=other.CON_mpindex_msrkceu_right;

    CON_fixedconsseq=other.CON_fixedconsseq;
    CON_fixedconsqual=other.CON_fixedconsqual;

    //
    CON_conscalc_mincov=other.CON_conscalc_mincov;
    CON_conscalc_minqual=other.CON_conscalc_minqual;
    CON_conscalc_missingchar=other.CON_conscalc_missingchar;
    CON_conscalc_assumediploid=other.CON_conscalc_assumediploid;
    CON_conscalc_allowiupac=other.CON_conscalc_allowiupac;
    CON_conscalc_addconstag=other.CON_conscalc_addconstag;
    CON_allconsseq=other.CON_allconsseq;
    CON_allconsqual=other.CON_allconsqual;
    CON_alladjustments=other.CON_alladjustments;
    CON_strainconsseq=other.CON_strainconsseq;
    CON_strainconsqual=other.CON_strainconsqual;
    CON_strainadjustments=other.CON_strainadjustments;

    CON_allconsseq=other.CON_allconsseq;
    CON_allconsqual=other.CON_allconsqual;
    CON_strainconsseq=other.CON_strainconsseq;
    CON_strainconsqual=other.CON_strainconsqual;

    CON_last_dangerous_overlaps=other.CON_last_dangerous_overlaps;

    CON_contains_long_repeats_only=other.CON_contains_long_repeats_only;
    CON_contains_majority_digitallynormalised_reads=other.CON_contains_majority_digitallynormalised_reads;

    CON_stats=other.CON_stats;

  }

  FUNCEND();
  return *this;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Contig::setParams(std::vector<MIRAParameters> * params)
{
  FUNCSTART("Contig::setParams(MIRAParameters * params)");
  CON_miraparams=params;
  FUNCEND();
}



/*************************************************************************
 *
 *  TODO: rework for strain dependend data
 *
 *
 *************************************************************************/
int32 Contig::paddedPos2UnpaddedPos(uint32 padpos)
{
  FUNCSTART("int32 Contig::paddedPos2UnpaddedPos(uint32 padpos)");

  ensureConsensus(-1);

  BUGIFTHROW(CON_alladjustments.empty(),"CON_alladjustments.empty()?");

  if(padpos>=CON_alladjustments.size()) padpos=CON_alladjustments.size()-1;

  FUNCEND();
  return CON_alladjustments[padpos];
}





/*************************************************************************
 *
 * builds a vector which shows at which positions in the consensus there
 *  are reads that have a given tag
 *
 *   .........TTTTTTTT..................
 *      ...........TTTTTTT.................  <- reads
 *
 *   00000000011111111111100000000000000000  <- maskshadow
 *
 * returns true if anything set, else false.
 *
 * 24.06.2012 Had newMarkPossibleRepeats() routine once dump
 *  Marking possibly misassembled repeats: tcmalloc: large alloc 0 bytes == (nil) @
 *  Ouch, out of memory detected.
 * on me. Almost impossible, additional checks inserted.
 *
 *************************************************************************/


bool Contig::buildMaskShadow(std::vector<int8> & maskshadow, std::vector<multitag_t::mte_id_t> & masktagstrings, bool onlybackbone)
{
  FUNCSTART("void Contig::buildMaskShadow(std::vector<int8> & maskshadow, std::vector<string> & masktagstrings)");

  bool retvalue=false;

  //cout << "MASK #############" << endl;

  BUGIFTHROW(CON_counts.size()==0,"CON_counts.size()==0 ???");

  maskshadow.clear();
  maskshadow.resize(CON_counts.size(),0);
  auto pcrI=CON_reads.begin();
  for( ; pcrI!=CON_reads.end(); ++pcrI) {
    if(onlybackbone && !pcrI->isBackbone()) continue;
    const Read & actread = *pcrI;
    for(uint32 tagnum=0; tagnum<actread.getNumOfTags(); tagnum++) {
      const multitag_t & acttag=actread.getTag(tagnum);
      for(const auto & mtse : masktagstrings){
	if(acttag.identifier == mtse) {
	  // ok, this tag is one of those we search
	  // get the consensusposition from the read tag positions
	  //  (may be outside of bounds as tags can be in clipped parts!)
	  int32 cfrom=pcrI.unclippedReadPos2ContigPos(acttag.from);
	  int32 cto=pcrI.unclippedReadPos2ContigPos(acttag.to);

	  // positions may be reverse, swap then
	  if(cfrom > cto) std::swap(cfrom,cto);

	  //cout << actread.getName() << "\t" << acttag.from << "\t" << acttag.to << "\t" << cfrom << "\t" << cto << endl;

	  // if both positions are left, right out of bounds, next tag
	  if(cfrom <0 && cto<0) continue;
	  if(cfrom >= CON_counts.size() && cto>=CON_counts.size()) continue;
	  // those positions ot of bounds -> set to bounds
	  if(cfrom < 0) cfrom=0;
	  if(cto >= CON_counts.size()) cto=CON_counts.size()-1;
	  // Tag from .. to are INCLUDING!
	  for(int32 i=cfrom; i<=cto; ++i) {
	    maskshadow[i]=1;
	    retvalue=true;
	  }
	}
      }
    }
  }

  FUNCEND();
  return retvalue;
}






/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::zeroStats()
{
  CON_stats=constats_t();

  for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; ++i){
    CON_stats.readsperst[i]=0;
    CON_stats.totalbasesperst[i]=0;
    CON_stats.max_covperst[i]=0;
    CON_stats.avg_covperst[i]=0.0;
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::calcStats()
{
  FUNCSTART("void Contig::calcStats()");

  if(CON_stats.statsvalid) return;

  zeroStats();
  finalise();

  std::string consseq;
  std::vector<base_quality_t> consqual;

  newConsensusGet(consseq, consqual);

  CON_stats.conlength=CON_counts.size();

  // using a pointer to CON_counts is (a lot) faster than accessing it
  //  via [] each time
  auto ccI=CON_counts.begin();

  //cout << "tililili! " << CON_name << '\n';
  for(uint32 ci=0; ci<CON_counts.size(); ++ci, ++ccI){
    if(consseq[ci]!='*') ++CON_stats.conlength_nogap;
    //cout << ci << '\t' << consseq[ci] << '\n';
    switch(toupper(consseq[ci])){
    case 'A':{
      ++CON_stats.AinC;
      break;
    }
    case 'C':{
      ++CON_stats.CinC;
      break;
    }
    case 'G':{
      ++CON_stats.GinC;
      break;
    }
    case 'T':{
      ++CON_stats.TinC;
      break;
    }
    case '*':{
      ++CON_stats.starInC;
      break;
    }
    case '-':
    case 'N': {
      ++CON_stats.NinC;
      break;
    }
    case 'X':{
      ++CON_stats.XinC;
      break;
    }
    default: {
      if(dptools::isValidIUPACBase(toupper(consseq[ci]))){
	++CON_stats.IUPACinC;
      }else{
	++CON_stats.FunnyInC;
      }
    }
    }

    // in backbone assemblies, acount for positions not covered
    // this is only good for backbones with one read, but cannot
    //   be helped at the moment
    // alternative would be a pre-analysis: use array size of CON_counts,
    //   go thorugh all reads (non backbone, non rail) and set to 1
    //   for position covered
    if(ccI->getBBChar()!='@' && ccI->total_cov==1){
      ++CON_stats.numnocoverage;
    }

    CON_stats.max_coverage=std::max(CON_stats.max_coverage,
			       static_cast<uint32>(ccI->total_cov));
    for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; ++i){
      CON_stats.totalbasesperst[i]+=static_cast<uint64>(ccI->seqtype_cov[i]);
      CON_stats.max_covperst[i]=std::max(CON_stats.max_covperst[i],
				    static_cast<uint32>(ccI->seqtype_cov[i]));
    }

    CON_stats.starInR+=ccI->star;
    CON_stats.NinR+=ccI->N;
  }

  // ok, in the section above, rail reads would have been also counted in
  //  (in .coverage)
  // compute reads statistics and remove rail bases from base statistics

  {
    auto cI=CON_reads.begin();
    for(;cI!=CON_reads.end();++cI){
      if(cI->isRail()){
	CEBUG("Rail: " << cI->getName() << ": " << CON_stats.totalbasesperst[cI->getSequencingType()]);
	CON_stats.totalbasesperst[cI->getSequencingType()]-=cI->getLenClippedSeq();
	CEBUG(" -> " << CON_stats.totalbasesperst[cI->getSequencingType()] << endl);
      }else{

	// Nono, backbones must be counted in!
	// Scenario: assembly at 6x used as backbone, assembled further 1x
	//  If backbones were not counted, then it would look like all contigs
	//  had only 1x. Well, it'd be 1x new, but that's not what most people
	//  would want
	// && !cI->read.isBackbone()) {

	++CON_stats.total_reads;
	CON_stats.readsperst[cI->getSequencingType()]++;
	if(cI->hasQuality()){
	  CON_stats.numreads_withqual++;
	}else{
	  CON_stats.numreads_noqual++;
	}
      }
    }
  }

  CON_stats.avg_coverage=0.0;
  if(!CON_counts.empty()){
    uint64 totalbases=0;
    for(uint32 sti=0; sti<ReadGroupLib::SEQTYPE_END; ++sti){
      CON_stats.avg_covperst[sti]=static_cast<double>(CON_stats.totalbasesperst[sti])/CON_counts.size();
      totalbases+=CON_stats.totalbasesperst[sti];
    }
    CON_stats.avg_coverage=static_cast<double>(totalbases)/CON_counts.size();
  }

  {
    double qual=0.0;
    for(auto & cqe : consqual) qual+=static_cast<double>(cqe);
    CON_stats.avg_conqual=static_cast<uint16>((qual/consseq.size())+.5);
  }

  if(CON_stats.AinC+CON_stats.CinC+CON_stats.GinC+CON_stats.TinC>0.0){
    CON_stats.gccontent=100.0/(CON_stats.AinC+CON_stats.CinC+CON_stats.GinC+CON_stats.TinC)*(CON_stats.CinC+CON_stats.GinC);
  }

  // IUPACs are perhaps not set yet
  updateStatsFromConsensusTags(true, true, false, true, true);

  CON_stats.contains_long_repeats_only=CON_contains_long_repeats_only;
  CON_stats.islargecontig=0; // undefined, to be set by assembly_info

  CON_stats.statsvalid=true;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::updateStatsFromConsensusTags(bool countSRMcs, bool countWRMcs, bool countIUPACs, bool countSTMUs, bool countSTMSs)
{
  if(countSRMcs) CON_stats.numSRMc=0;
  if(countWRMcs) CON_stats.numWRMc=0;
  if(countIUPACs) CON_stats.IUPACinC=0;
  if(countSTMUs) CON_stats.numSTMU=0;
  if(countSTMSs) CON_stats.numSTMS=0;
  for(const auto & cte : CON_consensus_tags){
    if(countSRMcs && cte.identifier == CON_tagentry_idSRMc) ++CON_stats.numSRMc;
    if(countWRMcs && cte.identifier == CON_tagentry_idWRMc) ++CON_stats.numWRMc;
    if(countIUPACs && cte.identifier == CON_tagentry_idIUPc) ++CON_stats.IUPACinC;
    if(countSTMUs && cte.identifier == CON_tagentry_idSTMU) ++CON_stats.numSTMU;
    if(countSTMSs && cte.identifier == CON_tagentry_idSTMS) ++CON_stats.numSTMS;
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::setContigName(const std::string & name)
{
  FUNCSTART("void Contig::setContigName(const std::string & name)");

  CON_name=name;

  FUNCEND();
  return;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::string & Contig::getContigName() const
{
  FUNCSTART("std::string Contig::getContigName() const");

  //if(!CON_name.empty()) {
  //  FUNCEND();
  //  return CON_name;
  //}

  if(!CON_name.empty()) {
    return CON_name;
  }

  std::ostringstream ostr;

  if(!CON_nameprefix.empty()) {
    ostr << CON_nameprefix;
  }else{
    ostr << "contig";
  }

  if(CON_contains_majority_digitallynormalised_reads==0){
    const_cast<Contig *>(this)->priv_calcNumDigitallyNormalisedReads();
  }

  if(CON_contains_majority_digitallynormalised_reads>0){
    ostr << "_dn";
  }

  if(CON_reads.size()>1) {
    if(CON_contains_long_repeats_only){
      ostr << "_rep_c" << CON_id;
    }else{
      ostr << "_c" << CON_id;
    }
  }else{
    ostr << "_s" << CON_id;
  }

  //std::string cname=ostr.str();
  CON_name=ostr.str();

  FUNCEND();
  return CON_name;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_calcNumDigitallyNormalisedReads()
{
  CON_contains_majority_digitallynormalised_reads=-1;

  //if(CON_reads.size()>6){
  //  auto crI=CON_reads.begin();
  //  uint64 numdn=0;
  //  for(; crI!=CON_reads.end(); ++crI){
  //    if(crI->hasTag(Read::REA_tagentry_idDGNr)) ++numdn;
  //  }
  //  if(numdn*10000/CON_reads.size() > 3000){
  //    CON_contains_majority_digitallynormalised_reads=1;
  //  }
  //}

  auto crI=CON_reads.begin();
  uint64 numdn=0;
  for(; crI!=CON_reads.end(); ++crI){
    if(crI->getDigiNormMultiplier()>1){
      CON_contains_majority_digitallynormalised_reads=1;
      break;
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint32 Contig::getNumBackbones() const
{
  uint32 ret=0;
  for(uint32 rgi=0; rgi<CON_readsperreadgroup.size(); ++rgi){
    if(CON_readsperreadgroup[rgi]>0
       && ReadGroupLib::getReadGroupID(rgi).isBackbone()) ret+=CON_readsperreadgroup[rgi];
  }
  return ret;
}


/*************************************************************************
 *
 * Returns strain IDs of backbone reads
 * if no bb read in contig, an empty vector
 *
 *
 *************************************************************************/

void Contig::getStrainsOfBackbone(std::vector<int32> & sob) const
{
  sob.clear();
  for(uint32 rgi=0; rgi<CON_readsperreadgroup.size(); ++rgi){
    if(CON_readsperreadgroup[rgi]) {
      auto rgid=ReadGroupLib::getReadGroupID(rgi);
      if(rgid.isBackbone()){
	sob.push_back(rgid.getStrainID());
      }
    }
  }

  // multiple read groups may be present for the same strain ...
  if(sob.size()>1){
    // use non-parallel sort, this is a tiny thing
    mstd::ssort(sob);
    mstd::unique(sob);
  }

  return;
}

///*************************************************************************
// *
// *
// *
// *
// *************************************************************************/
//
//bool Contig::hasConsensusTag(const std::string & identifier) const
//{
//  FUNCSTART("bool Contig::hasTag(const std::string & identifier) const");
//
//  bool returnit=false;
//  for(uint32 i=0; i<CON_consensus_tags.size();i++){
//    if(CON_consensus_tags[i].identifier==identifier){
//      returnit=true;           // we found this id in the consensus tags
//      break;                  // we don't need to examine further
//    }
//  }
//
//  FUNCEND();
//
//  return returnit;
//}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Contig::hasConsensusTag(const multitag_t::mte_id_t identifier) const
{
  FUNCSTART("bool Contig::hasConsensusTag(const  multitag_t::mte_id_t identifier) const");

  bool returnit=false;
  for(uint32 i=0; i<CON_consensus_tags.size();i++){
    if(CON_consensus_tags[i].identifier==identifier){
      returnit=true;           // we found this id in the consensus tags
      break;                  // we don't need to examine further
    }
  }

  FUNCEND();

  return returnit;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::setSpecialSRAddConditions(const int32 maxtotalerrors, const int32 maxgaps, const int32 maxmismatches, const int32 cleanoverlapends)
{
  CON_ssrc_maxtotalerrors=maxtotalerrors;
  CON_ssrc_maxgaps=maxgaps;
  CON_ssrc_maxmismatches=maxmismatches;
  CON_ssrc_cleanoverlapends=cleanoverlapends;
  if(maxtotalerrors < 0 && maxgaps < 0 && maxmismatches < 0  && cleanoverlapends <= 0) {
    CON_specialsraddconditions=false;
  }else{
    CON_specialsraddconditions=true;
  }
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::addFirstRead(int32 id, int8 direction)
{
  FUNCSTART("void Contig::addFirstRead(uint32 id, int8 direction)");

  BUGIFTHROW(!CON_reads.empty(), "used addFirstRead while there are already reads in contig.");

  definalise();
  CON_last_dangerous_overlaps.clear();

  uint32 len=CON_readpool->getRead(id).getLenClippedSeq();
  CEBUG("Len of clipped read: " << len << endl);

  CON_counts.resize(len,CON_concounts_zero);

  //check in this first read
  CON_reads.placeRead(CON_readpool->getRead(id),id,0,direction);

  if(CON_readpool->getRead(id).getTemplateID()>=0){
    CON_templates_present.insert(CON_readpool->getRead(id).getTemplateID());
  }

  auto coveragemultiplier=CON_readpool->getRead(id).getDigiNormMultiplier();

  if(direction>0){
    updateCountVectors(0,
		       len,
		       CON_readpool->getRead(id).getClippedSeqIterator(),
		       CON_readpool->getRead(id).getSequencingType(),
		       true,
		       coveragemultiplier);
    // Put base locks in CON_counts that this read produces
    initialiseBaseLocks();
  }else{
    MIRANOTIFY(Notify::INTERNAL, "untested direction < 0");
  }

  //cout << "JJJJJJJ " << CON_readsperstrain.size() << " " << static_cast<uint32>(ReadGroupLib::getNumOfStrains()) << endl;
  if(CON_readsperstrain.size() < ReadGroupLib::getNumOfStrains()){
    CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);
  }
  CON_readsperstrain[CON_readpool->getRead(id).getStrainID()]=coveragemultiplier;
  //cout << "JJJJJJJ " << CON_readsperstrain.size() << " " << static_cast<uint32>(ReadGroupLib::getNumOfStrains()) << endl;

  if(CON_readsperreadgroup.size() < ReadGroupLib::getNumReadGroups()){
    CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);
  }
  CON_readsperreadgroup[CON_readpool->getRead(id).getReadGroupID().getLibId()]=coveragemultiplier;

  if(CON_readpool->getRead(id).getLenClippedSeq() > CON_longestreadseen){
    CON_longestreadseen=static_cast<int32>(CON_readpool->getRead(id).getLenClippedSeq());
  }
  if(!CON_readpool->getRead(id).isBackbone() && CON_readpool->getRead(id).getLenClippedSeq() > CON_longestnonbbreadseen){
    CON_longestnonbbreadseen=static_cast<int32>(CON_readpool->getRead(id).getLenClippedSeq());
  }
  // should never happen as rails should be added in a different manner ... but one never knows
  // maybe in future?
  if(CON_readpool->getRead(id).isRail() && CON_readpool->getRead(id).getLenClippedSeq() > CON_longestrailseen){
    CON_longestrailseen=static_cast<int32>(CON_readpool->getRead(id).getLenClippedSeq());
  }

  if(CON_fixedconsseq.size()){
    nukeSTLContainer(CON_fixedconsseq);
    nukeSTLContainer(CON_fixedconsqual);
  }

  FUNCEND();
}




/*************************************************************************
 *
 * addRead_wrapped might change some alignment parameters
 *
 * returns:
 *   errstat for communicating errors
 *   templateguess for communicating measure template info
 *
 * as it returns from a number of different points within the loop,
 *  setting back original this function takes care to set everything back
 *  to normal if needed.
 *
 *************************************************************************/

//#define BUGHUNT

void Contig::addRead(std::vector<Align> & aligncache, const AlignedDualSeqFacts * initialadsf, int32 refid, int32 newid, int32 direction_frnid, bool newid_ismulticopy, int32 forcegrow, templateguessinfo_t & templateguess, errorstatus_t & errstat)
{
  FUNCSTART("void Contig::addRead(std::vector<Align> & aligncache, const AlignedDualSeqFacts * initialadsf, int32 refid, int32 newid, int32 direction_frnid, bool newid_ismulticopy, int32 forcegrow, templateguess_t & templateguess, errorstatus_t & errstat)");

  //setCEBUGFlag(newid,refid);

  paranoiaBUGIF(aligncache.size()!=ReadGroupLib::SEQTYPE_END,
		MIRANOTIFY(Notify::INTERNAL, "aligncache is not size of available sequencing types???"));

  templateguess.rgid.resetLibId();

  // this is ugly: save the align params because
  //  addRead_wrapped() might change some *sigh*
  align_parameters oldalignparams=
    (*CON_miraparams)[CON_readpool->getRead(newid).getSequencingType()].getAlignParams();

  try {
#ifdef BUGHUNT
    {
      //rotate
      std::string system_rmdir = static_cast<std::string>("mv addRead_debug_pre.log addRead_debug_pre.log1");
      int tmp=system(system_rmdir.c_str());
      system_rmdir = static_cast<std::string>("mv addRead_debug_post.log addRead_debug_post.log1");
      tmp=system(system_rmdir.c_str());
      (void) tmp;

      std::ofstream eout("addRead_debug_pre.log", std::ios::out | std::ios::trunc);
      eout << "This file is for debugging and should never be seen by users!\n\n"
	"If you do see it, contact the author (bach@chevreux.org) immediately\n\n"
	"addRead_debug_pre.log!\n\n";
      priv_dumpReplay(eout,initialadsf,refid,newid,direction_frnid,newid_ismulticopy,forcegrow);
    }
#endif

    addRead_wrapped(aligncache,
		    initialadsf,
		    refid,
		    newid,
		    direction_frnid,
		    newid_ismulticopy,
		    forcegrow,
		    templateguess,
		    errstat);

    if(errstat.code!=ENOERROR){
      // remove an eventual guess for template placement
      templateguess.rgid.resetLibId();
      templateguess.splace_seen=ReadGroupLib::SPLACE_UNKNOWN;
      templateguess.tsize_seen=0;

      // check if refid has been put into reads_affected by the contig
      // if not, warn and dump
      bool isin=false;
      for(auto raid : errstat.reads_affected){
	if(raid == refid){
	  isin=true;
	  break;
	}
      }
      if(!isin){
	// TODO: revert back to throwing at some point
	cout << "\nWARNING/ERROR: ignored a missing refid" << endl;
	errstat.reads_affected.push_back(refid);
	//cout << "whooops, no refid? I see:\n";
	//for(auto raid : errstat.reads_affected){
	//  cout << CON_readpool->getRead(raid).getName() << endl;
	//}
	//BUGIFTHROW(true,"WARNING/ERROR: contig did not put refid into affected reads. Please file a bug report to\n\thttps://sourceforge.net/apps/trac/mira-assembler/\nand please send a mail to\n\tmira_talk@freelists.org\n");
      }
    }

#ifdef BUGHUNT
    {
      std::ofstream eout("addRead_debug_post.log", std::ios::out | std::ios::trunc);
      eout << "This file is for debugging and should never be seen by users!\n\n"
	"If you do see it, contact the author (bach@chevreux.org) immediately\n\n"
	"addRead_debug_post.log!\n\n";
      dumpReplay(eout,initialadsf,refid,newid,direction_frnid,newid_ismulticopy,forcegrow);
    }
#endif

  }
  catch(Notify n){
    std::ofstream eout("addReaderror_REPLAY_do_not_delete.log", std::ios::out | std::ios::trunc);
    priv_dumpReplay(eout,initialadsf,refid,newid,direction_frnid,newid_ismulticopy,forcegrow);
    cout << "\n\n\nA file named 'addReaderror_REPLAY_do_not_delete.log' has been written"
      "\nto the working directory. Please send it to the author to get that bug fixed!\n\n\n";
    n.handleError(THISFUNC);
  }

  //CON_miraparams->setAlignParams(oldalignparams);
  //(*CON_miraparams)[CON_readpool->getRead(newid).getSequencingType()].setAlignParams(oldalignparams);
  (*CON_miraparams)[CON_readpool->getRead(newid).getSequencingType()].getNonConstAlignParams()=oldalignparams;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUGFLAG
//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

//#define CEBUG(bla)   {if(docebug) {cout << bla; cout.flush();}}

#define CEBUG(bla)

void Contig::addRead_wrapped(std::vector<Align> & aligncache, const AlignedDualSeqFacts * initialadsf, int32 refid, int32 newid, int32 direction_frnid, bool newid_ismulticopy, int32 forcegrow, templateguessinfo_t & templateguess, errorstatus_t & errstat)
{
  FUNCSTART("void Contig::addRead_wrapped(std::vector<Align> & aligncache, const AlignedDualSeqFacts * initialadsf, int32 refid, int32 newid, int32 direction_frnid, bool newid_ismulticopy, int32 forcegrow, errorstatus_t & errstat)");

  //BUGIFTHROW(CON_cebugflag && CON_readpool->getRead(newid).getName()=="G8BUD:102:933","gna");

//#define ALIGNCHECK
//#define EXTRATRACEHELP

  //bool docebug=false;
  //if(CON_readpool->getRead(newid).getName()=="MISEQ:5:000000000-A1BMW:1:1101:25577:9445/1") {
  //  docebug=true;
  //}

  if(CON_readpool->getRead(newid).getLenClippedSeq()==0){
    // BaCh 29.05.2018
    // 'demoted' response from a fatal error to 'gracefully reject'
    // Reason: during the two-stage mapping, reads might be clipped after the first stage,
    //         e.g., in Contig::clipHomoPolyWithGapsAtEnds()
    //         This might bring reads to a length of zero and subsequent mapping
    //         triggers this response.
    //         Happened during a mapping of very old 36mers.
    //         Therefore, be graceful.
    errstat.code=EZEROLEN;
    return;
//    cout << "Gaaaaaah! Trying to add a read of length 0???\n";
//    Read::setCoutType(Read::AS_TEXT);
//    cout << CON_readpool->getRead(newid) << endl;
//    MIRANOTIFY(Notify::INTERNAL, "Gaaaaaah! Trying to add a read of length 0???");
  }

  timeval us_start;
  gettimeofday(&us_start,nullptr);

  if(getNumBackbones()>0) CON_tmpcons_from_backbone=true;

  errstat.code=ENOERROR;
  errstat.reads_affected.clear();

  uint8 newreadseqtype=CON_readpool->getRead(newid).getSequencingType();
  MIRAParameters & rt_params=(*CON_miraparams)[newreadseqtype];

  //CEBUG(*this);
  CEBUG("\nCON_counts.size(): " << CON_counts.size() << '\n');

#ifndef PUBLICQUIET
  cout << "###Adding: " <<CON_readpool->getRead(newid).getName() << " ("
       << ReadGroupLib::getShortNameOfSequencingType(newreadseqtype)
       << "," << newid << ")\t"
       << static_cast<uint16>(CON_readpool->getRead(newid).getStrainID()) << " (" << CON_readpool->getRead(newid).getStrainName() << ") ";
#endif

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  switch(newreadseqtype){
  case ReadGroupLib::SEQTYPE_SANGER :
  case ReadGroupLib::SEQTYPE_454GS20 :
  case ReadGroupLib::SEQTYPE_IONTORRENT :
  case ReadGroupLib::SEQTYPE_PACBIOHQ :
  case ReadGroupLib::SEQTYPE_PACBIOLQ :
  case ReadGroupLib::SEQTYPE_TEXT :
  case ReadGroupLib::SEQTYPE_SOLEXA : {
    break;
  }
  case ReadGroupLib::SEQTYPE_ABISOLID : {
    MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 4.");
    break;
  }
  default : {
    cout << "\nSeqtype: " << newreadseqtype << endl;
    MIRANOTIFY(Notify::INTERNAL, "Trying to add unknown seqtype? Needs more code 4.");
  }
  }


  CEBUG("\tNewid: " << newid << '\t' << CON_readpool->getRead(newid).getName() << endl);

#ifdef PARANOIABUGTRACKFLAG
  cout << "Warning: PARANOIABUGTRACKFLAG is compiled!";
  if(CON_readpool->getRead(newid).hasValidData()==false){
    cout << "Read has no valid data???";
    MIRANOTIFY(Notify::INTERNAL, "Read has no valid data???") ;
  }
  if(CON_readpool->getRead(newid).checkRead()){
    cout << "Read " << newid << " in readpool is trashed\n";
    cout << CON_readpool->getRead(newid);
    MIRANOTIFY(Notify::INTERNAL, "*sigh* read trashed") ;
  }
#endif

  definalise();
  CON_last_dangerous_overlaps.clear();

//  CEBUG("I've got this contig:\n"<<*this);
//  CEBUG("Must align this:\n");

  if(CON_reads.empty()){
	addFirstRead(refid, 1);
	CEBUG(*this);
	FUNCEND();
	return;
  }
  // should be set by addFirstRead
  BUGIFTHROW(CON_readsperstrain.size() < ReadGroupLib::getNumOfStrains(), "after addFirstRead: CON_readsperstrain.size() " << CON_readsperstrain.size() << " < ReadGroupLib::getNumOfStrains(() " << static_cast<uint16>(ReadGroupLib::getNumOfStrains()) << " ?");

  //// do not check readsperreadgroup
  //// Legal scenario: load two backbone, do a mapping, once first contig is done it will get CER reads
  ////  and contig 2 does not know anything of that
  //if(CON_readsperreadgroup.size() < ReadGroupLib::getNumReadGroups()){
  //  ReadGroupLib::debugDumpReadGroupInfo(cout);
  //  MIRANOTIFY(Notify::FATAL,"after addFirsRead: CON_readsperreadgroup.size() " << CON_readsperreadgroup.size() << " < CON_readpool->getNumReadGroups() " <<  ReadGroupLib::getNumReadGroups() << " ?");
  //}

#ifndef PUBLICQUIET
  cout << "From: " <<CON_readpool->getRead(refid).getName() << " ("
       << ReadGroupLib::getShortNameOfSequencingType(CON_readpool->getRead(refid).getSequencingType())
       << "," << refid << ")\t";

  cout << initialadsf->getOverlapLen() << '/';
  cout << static_cast<uint16>(initialadsf->getScoreRatio()) << "%\t";
#endif

  BUGIFTHROW(CON_readpool->getRead(newid).isRail(),"Tried to add a rail???");

  {
    bool strangeerror=false;
    if( !(initialadsf->getID1() == refid
	  || initialadsf->getID2() == refid)) {
      cout << "refid not in initialadsf???\n";
      strangeerror=true;
    }
    if( !(initialadsf->getID1() == newid
	  || initialadsf->getID2() == newid)) {
      cout << "newid not in initialadsf???\n";
      strangeerror=true;
    }
    if(strangeerror){
      MIRANOTIFY(Notify::INTERNAL,"Impossible error while adding read, investigate! This may be due faulty RAM.");
    }
  }


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  switch(CON_readpool->getRead(refid).getSequencingType()){
  case ReadGroupLib::SEQTYPE_SANGER :
  case ReadGroupLib::SEQTYPE_454GS20 :
  case ReadGroupLib::SEQTYPE_IONTORRENT :
  case ReadGroupLib::SEQTYPE_PACBIOHQ :
  case ReadGroupLib::SEQTYPE_PACBIOLQ :
  case ReadGroupLib::SEQTYPE_TEXT :
  case ReadGroupLib::SEQTYPE_SOLEXA : {
    break;
  }
  case ReadGroupLib::SEQTYPE_ABISOLID : {
    MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 4a.");
    break;
  }
  default : {
    cout << "\nSeqtype: " << CON_readpool->getRead(refid).getSequencingType() << endl;
    MIRANOTIFY(Notify::INTERNAL, "Having unknown seqtype as reference? Needs more code 4a.");
  }
  }

#ifdef PARANOIABUGTRACKFLAG
  // Test, if new read already in contig
  {
    auto crI= CON_reads.cbegin();
    for(uint32 i=0; i< CON_reads.size(); ++i, ++crI){
      //      cout << " Checking:" << crI->id;
      if(crI->orpid==newid) break;
    }
    if(i!=CON_reads.size()){
      dumpConReads();
      Read::setCoutType(Read::AS_TEXTSHORT);
      cout << CON_readpool->getRead(newid);
      MIRANOTIFY(Notify::INTERNAL, "new ID already in contig???");
    }
  }

  // Test, if ref read already in contig
  {
    auto crI= CON_reads.cbegin();
    for(uint32i=0; i< CON_reads.size(); ++i, ++crI){
      //      cout << " Checking:" << crI->id;
      if(crI->orpid==refid) break;
    }
    if(i==CON_reads.size()){
      dumpConReads();
      Read::setCoutType(Read::AS_TEXTSHORT);
      cout << CON_readpool->getRead(newid);
      MIRANOTIFY(Notify::INTERNAL, "ref ID not in contig???");
    }
  }
#endif

  CEBUG("Refid len:" << CON_readpool->getRead(refid).getLenClippedSeq() << '\n');
  CEBUG("Newid len:" << CON_readpool->getRead(newid).getLenClippedSeq() << '\n');

  CEBUG("Refid: " << refid << endl);
  CEBUG('\n'<< *initialadsf);


  CON_us_steps[USCLO_PRE]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  auto conrefreadI = CON_reads.getIteratorOfReadpoolID(refid);
  BUGIFTHROW(conrefreadI==CON_reads.end(),"RefID " << refid << " (" << CON_readpool->getRead(refid).getName() << ") not found on contig?");
  int32 offsetrefid=conrefreadI.getReadStartOffset();

  //gna if(errstat.reads_affected.empty())

  if(!CON_allowedrefids.empty()){
    if(!CON_allowedrefids[refid]){
      errstat.code=EREFIDNOTALLOWED;
      if(CON_readpool->getRead(refid).isRail()) {
	// we did not compute xcut and ycut yet.
	// to save time, we already stop here and take the positions
	//  of the refid rail as xcut and ycut
	priv_arw_getReadsAffected(refid, newid, errstat.reads_affected,
				  conrefreadI.getReadStartOffset(),
				  conrefreadI.getReadStartOffset()+conrefreadI->getLenClippedSeq());
      }else{
	errstat.reads_affected.push_back(refid);
      }
      FUNCEND();
      return;
    }
  }



  int32 direction_refid=conrefreadI.getReadDirection();


  CON_us_steps[USCLO_DIR]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);


  int32 direction_newid_incontig=direction_refid*direction_frnid;

  CEBUG("dir_frnid: " << direction_frnid);
  CEBUG("\tdir_refid: " << direction_refid<< endl;);
  CEBUG("Offset Refid:" << offsetrefid << " - " << offsetrefid+conrefreadI->getLenClippedSeq() <<endl;);

  //Align bla(CON_miraparams);
  std::list<AlignedDualSeq> madsl;
  decltype(madsl.cbegin()) adsI;


  // xcut and ycut are the posintions in the contig where a
  //  sequence must be made from for SW alignment
  // xcut and ycut are defined [...[  (xcut including, ycut excluding)

  // the ycut is less critical, but also influences slightly the
  //  speed of the banded SW ... there might be a little bit more to
  //  compute

  int32 xcut=offsetrefid; // will be adapted

  int32 xcutinc=1;
  int32 deltax;
  // same direction in contig and reference ads?
  if(direction_refid*initialadsf->getSequenceDirection(refid)>0){
    // yes
    if(initialadsf->getOffsetInAlignment(refid)==0){
      CEBUG("R1 offset.\n");
      deltax=initialadsf->getOffsetInAlignment(newid);

    }else{
      CEBUG("R2 offset.\n");
      deltax=initialadsf->getOffsetInAlignment(refid);
      xcutinc=-1;
    }
  }else{
    if(initialadsf->getRightOffsetInAlignment(refid)==0){
      CEBUG("R3 offset.\n");
      deltax=initialadsf->getRightOffsetInAlignment(newid);
    }else{
      CEBUG("R4 offset.\n");
      deltax=initialadsf->getRightOffsetInAlignment(refid);
      xcutinc=-1;
    }
  }

  CEBUG("deltax: " << deltax << '\n');
  CEBUG("xcutinc: " << xcutinc << '\n');


  int32 runlength=CON_readpool->getRead(newid).getLenClippedSeq();
  // BaCh 14.11.2012
  // to accomodate for read mappings which have a lot of inserts, runlength must be increased
  //  it may be that the length of the overlap is significantly larger than the initial read length
  //  e.g. Illumina 100bp with a clean 15bp insert makes for a 115bp overlap
  if(runlength < initialadsf->getOverlapLen()){
    runlength=initialadsf->getOverlapLen();
  }

  int32 ycut;

  // search for xcut: use refread sequence when possible,
  //  else the concount_t structure
  {
    int32 refpos=0;
    const char * refseq;
    {
      if(direction_refid>0){
	refseq=conrefreadI->getClippedSeqAsChar();
      }else{
	refseq=conrefreadI->getClippedComplementSeqAsChar();
      }
      auto ccI=CON_counts.begin();
      BOUNDCHECK(xcut, 0, CON_counts.size()+1);
      if(xcutinc<0 && xcut>0) xcut--;
      std::advance(ccI, xcut);
      int32 refgaps=0;
      bool refnotrail=!CON_readpool->getRead(refid).isRail();
      for(;deltax>=0; xcut+=xcutinc, refpos+=xcutinc, refseq+=xcutinc){
	CEBUG("deltax: " << deltax << "\txcut: " << xcut);
	if(refnotrail && refpos>=0 && refpos<conrefreadI->getLenClippedSeq()){
	  //if(0){
	  // BaCh 14.11.2012: hmmm ... I suppose I forgot to kill this branch (like in the similar loop below)
	  // BaCh 2013-03-08: NO! I did not forget to kill that branch, it should actually survive
	  //                  de-novo assemblies!!!
	  //                  There are therefore two conflicting requirements: de-novo which needs to be as
	  //                  exact as possible with regards to the reference sequence, and mapping which
	  //                  "in the second mapping round" should try to adhere to the intermediate
	  //                  consensus. Solution: take refread as long as possible if refread is not rail.
	  CEBUG("\tRead: " << *refseq << '\n');
	  if(*refseq!='*'){
	    deltax--;
	  }
	  ccI+=xcutinc;
	}else{
	  // BaCh 06.12.2012 ; error reported by "Kris" on 26.11.2012
	  //    "Trying to dereference an hditer pointing to end()???"
	  //    ->Thrown: inline TT & HDeque::hditer::dereference() const
	  //
	  // I am mystified: where did the "+1" below come from?
	  //     if(xcut>=0 && xcut<CON_counts.size()+1){
	  // even more mystifying: how did it survive so long after switching to hdeque???
	  //
	  // extremely vexing: see comment from 16.05.2010 "Why on earth did I have +1"
	  //  some 100 line below :-(((
	  if(xcut>=0 && xcut<CON_counts.size()){
	    // this part mimicks the makeTmpConsenus() behaviour of calling stars
	    if(CON_tmpcons_from_backbone
	       && ccI->getOriginalBBChar()!='N'
	       && ccI->getOriginalBBChar()!='X'
	       && ccI->getOriginalBBChar()!='@'){
	      if(ccI->getOriginalBBChar()!='*') {
		CEBUG("\toBB1: (base) ");
		deltax--;
	      }else{
		CEBUG("\tBB1: * ");
		++refgaps;
	      }
	      CEBUG(ccI->i_backbonecharorig << " " << ccI->i_backbonecharupdated << '\n');
	    }else{
	      ccctype_t maximum= std::max(ccI->A, std::max(ccI->C, std::max(ccI->G, ccI->T)));
	      if(unlikely(ccI->total_cov==0)){
		// BaCh 30.11.2012
		// should normally never happen, certainly not in de-novo
		// but the two-pass mapping may have this at the end of the contigs after first pass
		//  (should I decide not to go the chompFront() / chompBack() after 1st pass)
		//
		// treat it like a base (well, will be N)
		CEBUG("\tno read\n");
		deltax--;
	      }else if(maximum >0 && maximum > ccI->star) {
		//if(maximum/4 >= ccI->star) {
		if(maximum/4 >= (ccI->star)*2) {
		  deltax--;
		  CEBUG("\tCON: (base)\n");
		}else{
		  CEBUG("\tCON: *\n");
		  ++refgaps;
		}
	      }else{
		if(!((ccI->star >= ccI->X)
		     && (ccI->star >= ccI->X))){
		  deltax--;
		  CEBUG("\tCON: (base)\n");
		}else{
		  CEBUG("\tCON: *\n");
		  ++refgaps;
		}
	      }
	    }
	    ccI+=xcutinc;
	  }else{
	    // TODO: replace with one += and a break out of loop
	    CEBUG("\tccI out of contig bounds.\n");
	    deltax--;
	    runlength--;
	  }
	}
      }

      // Bach 14.11.2012
      // only important for mapping assemblies, and then only when ref sequence is from backbone, do not apply on de-novo!
      if(CON_tmpcons_from_backbone){
	// Step 1: have we landed in a gap? If yes, get the full length of the gap
	while(xcut>=0 && xcut<CON_counts.size()+1){
	  if(ccI->getOriginalBBChar()!='*') break;
	  ++refgaps;
	  xcut+=xcutinc;
	  ccI+=xcutinc;
	}
      }
    }

    CEBUG("xcut 1: " << xcut <<'\n');

    // xcut, and therefore the initial ycut, may be way negative, beware!
    ycut=xcut;
    if(ycut<0) ycut=0;

    {
      auto ccI=CON_counts.begin();
      if(ycut>0) std::advance(ccI, ycut);
      int32 refgaps=0;
      bool refnotrail=!CON_readpool->getRead(refid).isRail();
      for(;runlength>=0; ycut++, refpos++, refseq++){
	CEBUG("runlength: " << runlength << "\tycut: " << ycut);
	if(refnotrail && refpos>=0 && refpos<conrefreadI->getLenClippedSeq()){
	  // also see comment in loop above regarding refnotrail and the whole if statement
	  //if(0){
	  CEBUG("\tRead: " << *refseq << '\n');
	  if(*refseq!='*'){
	    runlength--;
	  }
	}else{
	  // Changed 16.05.2010:
	  // Found by using IndexedDeque segfaulting due to out of bounds
	  // Why on earth did I have +1
	  //  if(ycut>=0 && ycut<CON_counts.size()+1){
	  if(ycut>=0 && ycut<CON_counts.size()){
	    // this part mimicks the makeTmpConsenus() behaviour of calling stars
	    if(CON_tmpcons_from_backbone
	       && ccI->getBBChar()!='@'
	       && ccI->getBBChar()!='N'
	       && ccI->getBBChar()!='X'){
	      if(ccI->getBBChar()!='*') {
		CEBUG("\tBB2: (base) ");
		runlength--;
	      }else{
		CEBUG("\tBB2: * ");
		++refgaps;
	      }
	      CEBUG(ccI->i_backbonecharorig << " " << ccI->i_backbonecharupdated << '\n');
	    }else{
	      ccctype_t maximum= std::max(ccI->A, std::max(ccI->C, std::max(ccI->G, ccI->T)));
	      	      if(unlikely(ccI->total_cov==0)){
		// BaCh 30.11.2012
		// should normally never happen, certainly not in de-novo
		// but the two-pass mapping may have this at the end of the contigs after first pass
		//  (should I decide not to go the chompFront() / chompBack() after 1st pass)
		//
		// treat it like a base (well, will be N)
		CEBUG("\tno read\n");
		--runlength;
	      }else if(maximum >0 && maximum > ccI->star) {
		//if(maximum/4 >= ccI->star) {
		if(maximum/4 >= (ccI->star)*2) {
		  runlength--;
		  CEBUG("\tCON1: (base): " << *ccI << endl);
		}else{
		  CEBUG("\tCON1: *: " << *ccI << endl);
		++refgaps;
		}
	      }else{
		if(!((ccI->star >= ccI->X)
		     && (ccI->star >= ccI->X))){
		  runlength--;
		  CEBUG("\tCON2: (base): " << *ccI << endl);
		}else{
		  CEBUG("\tCON2: *: " << *ccI << endl);
		  ++refgaps;
		}
	      }
	    }
	    ccI++;
	  }else{
	    CEBUG("\tccI out of bounds.\n");
	    runlength--;
	  }
	}
      }

      // Revised: 16.01.2013: Only for mapping where sequence comes from backbone!
      // corrector for larger gaps where skim might not have told the whole truth by giving alignment
      //  offset for right part of the alignment ... xcut would be to far right, too.
      if(CON_tmpcons_from_backbone && refgaps){
	CEBUG("refgaps? " << direction_frnid << "\t" << initialadsf->get5pLenContiguousMatch(newid) << "\t" << initialadsf->get3pLenContiguousMatch(newid) << "\n");
	if(direction_frnid>0){
	  if(initialadsf->get5pLenContiguousMatch(newid)==0 && initialadsf->get3pLenContiguousMatch(newid)>0){
	    xcut-=refgaps;
	    CEBUG("refgaps xcut corrector f: " << refgaps << '\n');
	  }
	}else{
	  if(initialadsf->get5pLenContiguousMatch(newid)>0 && initialadsf->get3pLenContiguousMatch(newid)==0){
	    xcut-=refgaps;
	    CEBUG("refgaps xcut corrector r: " << refgaps << '\n');
	  }
	}
      }
    }
  }

  CEBUG("ycut 1: " << ycut <<'\n');

  ycut+=10; // add safety distance at the end

  // -2  as safety distance in front
  xcut-=2;

  CEBUG("xcut 2: " << xcut <<'\n');
  CEBUG("ycut 2: " << ycut <<'\n');


  // in some cases, xcut may be > size of contig (and ycut anyway)
  //  or ycut < 0 (and xcut anyway)
  // this can happen when a read is edited over and over again during
  //  contig assembly and then the expected offset in the adsfact
  //  is way off target
  // occurs with short matches at end of reads.
  //
  // only possibility to handle this at this stage: reject alignment

  if(xcut >= static_cast<int32>(CON_counts.size()) || ycut < 0){
#ifndef PUBLICQUIET
    cout << "rej: no align found (bounds 1)\t";
    cout.flush();
#endif
    errstat.code=ENOALIGN;
    errstat.reads_affected.push_back(refid);

    CON_us_steps[USCLO_XCUT]+=diffsuseconds(us_start);

    FUNCEND();
    return;
  }


  // but check that we're not hitting a gap base at the xcut position
  //  if yes, go back as far as needed to find a non-gap
  // Can happen because of the xcut-=2 above:
  if(xcut>0){
    auto ccI=CON_counts.begin();
    BOUNDCHECK(xcut, 0, CON_counts.size()+1);
    // we should not be out of bounds (see "bounds 1" check just above), but just in case
    if(xcut==CON_counts.size()) --xcut;
    std::advance(ccI, xcut);
    while(xcut>0){
      ccctype_t maximum= std::max(ccI->A, std::max(ccI->C, std::max(ccI->G, ccI->T)));
      if(maximum >0 && maximum > ccI->star) {
  	if(maximum/4 >= ccI->star) {
  	  // base
  	  break;
  	}
      }else{
  	if(!((ccI->star >= ccI->X)
  	     && (ccI->star >= ccI->X))){
  	  // base
  	  break;
  	}
      }
      --ccI;
      --xcut;
    }
  }

  CEBUG("xcut 3: " << xcut <<'\n');

  // especially in mapping alignments with many SNPs/indels and partial
  //  overlaps with the rail reads, we might have landed completely
  //  outside the reference rail. Darn.
  //
  // let's deal with that. It's a hack, and a bad one.

  // static_cast needed or gcc will convert RHS to unsigned, then LHS to unsigned ... and as xcut may be
  //  negative, hilarity ensues.
  if(xcut >= static_cast<int32>(offsetrefid + conrefreadI->getLenClippedSeq())){
    xcut=offsetrefid + conrefreadI->getLenClippedSeq()-10;
    if(xcut<0) xcut=0;
  }
  if(ycut <= offsetrefid){
    ycut=offsetrefid+10;
    if(ycut>CON_counts.size()) ycut=CON_counts.size();
  }


  CEBUG("xcut final: " << xcut <<'\n');
  CEBUG("ycut final: " << ycut <<'\n');

  BUGIFTHROW(xcut>ycut,"final: xcut " << xcut << " > ycut " << ycut);

  CON_us_steps[USCLO_XCUT]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);


  // Template handling1

  bool havematchingtemplatepartner=false;

  if(CON_readpool->getRead(newid).getTemplatePartnerID() != -1){
    // 1st: direction
#ifndef PUBLICQUIET
    cout << " tmplhand1 ";
    cout << CON_readpool->getRead(newid).getTemplatePartnerID() << " ";
    cout.flush();
#endif
    auto tppcrI=CON_reads.getIteratorOfReadpoolID(CON_readpool->getRead(newid).getTemplatePartnerID());
   // is partner already in contig?
    if(tppcrI==CON_reads.end()){
#ifndef PUBLICQUIET
      cout << "nic("
	   << CON_readpool->getRead(newid).getTemplatePartnerID();
      if(CON_readpool->getRead(newid).getTemplatePartnerID()>=0){
	cout << " - " << CON_readpool->getRead(CON_readpool->getRead(newid).getTemplatePartnerID()).getName();
      }
      cout << ") ";
      cout.flush();
#endif
    }else{
      // assume the template matches
      // if not, the checks below will correct for that
      havematchingtemplatepartner=true;

      // ** Check for direction **
      if(CON_readpool->getRead(newid).getReadGroupID().getSegmentPlacementCode() != ReadGroupLib::SPLACE_UNKNOWN){
#ifndef PUBLICQUIET
	cout << "dir";
	cout.flush();
#endif

	// ok, we found the read with the template which corresponds to the newly
	//  to insert read. Now check if they meet the constraints.
	//  If not, do not insert

#ifndef PUBLICQUIET
	cout << "\ndni: " << static_cast<int32>(direction_newid_incontig);
	cout << "\npdi: " << static_cast<int32>(tppcrI.getReadDirection());
	cout << "\ntbd: " << static_cast<int32>(CON_readpool->getRead(newid).getTemplateBuildDirection());
#endif
	if(direction_newid_incontig*tppcrI.getReadDirection() != CON_readpool->getRead(newid).getTemplateBuildDirection()){
	  // not direction wanted, not good
#ifndef PUBLICQUIET
	  cout << "templ in wrong dir";
#endif
	  havematchingtemplatepartner=false;
	  if(!CON_readpool->getRead(newid).getReadGroupID().getSPInfoOnly()){
	    errstat.code=ETEMPLATEDIRECTION;
	    if(xcut<0) xcut=0;
	    if(static_cast<uint32>(ycut)>CON_counts.size()){
	      ycut=CON_counts.size();
	    }
	    priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);
	    FUNCEND();
	    return;
	  }
	}
      }

      // for segment placement and template size check, we need the positions in the contig

      bool newreadisleft=true;

      // Note: play it safe while adding right extend ... take only
      //  part of it
      int32 leftrx=xcut;
      int32 leftry=ycut;
      if(direction_newid_incontig>0) {
	leftry+=CON_readpool->getRead(newid).getRightExtend()*2/3;
      }else{
	leftrx-=CON_readpool->getRead(newid).getRightExtend()*2/3;
      }

      int32 rightrx=tppcrI.getReadStartOffset();
      int32 rightry=tppcrI.getReadStartOffset()+tppcrI->getLenClippedSeq();
      if(direction_refid>0) {
	rightry+=tppcrI->getRightExtend()*2/3;
      }else{
	rightrx-=tppcrI->getRightExtend()*2/3;
      }

      if(rightrx<leftrx) {
	std::swap(rightrx,leftrx);
	std::swap(rightry,leftry);
	newreadisleft=false;  // well, yeah, it is on the right
      }

      int32 actinsertsize=rightry - leftrx;
      if(rightry<leftry) actinsertsize=leftry - leftrx;

      // we now have all info needed for storing measured template info
      bool guesstemplatesegplace=true;
      templateguess.tsize_seen=0; // just as default
      templateguess.splace_seen=ReadGroupLib::SPLACE_UNKNOWN; // just as default

      // Due to clipping at ends of reads,
      // E.g.:
      //         ----------->
      //             <----------
      // which can become
      //                ---->
      //             <----
      //
      // one should only guess on disjunct (non-overlapping) reads
      // or if overlapping, the right read must start >= 10 bp
      //  later than the left read
      if(leftry<rightrx){
	// if overlapping
	if(leftrx+10>=rightrx){
	  // and start of reads <= 10bp
	  guesstemplatesegplace=false;
	}
      }

#ifndef PUBLICQUIET
      cout << "GUESS TEMPLATE? " << guesstemplatesegplace << endl;
#endif

      if(guesstemplatesegplace){
	templateguess.rgid=CON_readpool->getRead(newid).getReadGroupID();  // this validates the entry (else its  rgid: 0)
	templateguess.tsize_seen=actinsertsize;
	// placement code is hardest
	if(direction_newid_incontig*tppcrI.getReadDirection() < 0){
#ifndef PUBLICQUIET
	  cout << "TEMPLATE FR/RF\n";
#endif
	  // FR or RF
	  if(newreadisleft){
	    if(direction_newid_incontig>0){
	      templateguess.splace_seen=ReadGroupLib::SPLACE_FR;
	    }else{
	      templateguess.splace_seen=ReadGroupLib::SPLACE_RF;
	    }
	  }else{
	    if(direction_newid_incontig>0){
	      templateguess.splace_seen=ReadGroupLib::SPLACE_RF;
	    }else{
	      templateguess.splace_seen=ReadGroupLib::SPLACE_FR;
	    }
	  }
	}else{
#ifndef PUBLICQUIET
	  cout << "TEMPLATE SAME DIR\n";
#endif
	  // same direction ...
	  //
	  // I'm sure that some boolean logic could get through this, but I'm not inclined to think about it now
	  templateguess.splace_seen=ReadGroupLib::SPLACE_SU; // just as default
	  if((direction_newid_incontig>0 && newreadisleft==true)
	     || (direction_newid_incontig<0 && newreadisleft==false)){
	    if(CON_readpool->getRead(newid).getTemplateSegment()==1){
	      templateguess.splace_seen=ReadGroupLib::SPLACE_SF;
	    }else{
	      templateguess.splace_seen=ReadGroupLib::SPLACE_SB;
	    }
	  }else{
	    if(CON_readpool->getRead(newid).getTemplateSegment()==1){
	      templateguess.splace_seen=ReadGroupLib::SPLACE_SB;
	    }else{
	      templateguess.splace_seen=ReadGroupLib::SPLACE_SF;
	    }
	  }
	}
      }
#ifndef PUBLICQUIET
      cout << " nril " << newreadisleft
	   << " " << leftrx << "-" << leftry
	   << " " << rightrx << "-" << rightry
	   << " " << templateguess;
      cout.flush();
#endif


      // ** Check for segment placement **
      // unknown does obviously need no check, "samedir unknown" already checked implicitly by direction check above
      {
	auto spc=CON_readpool->getRead(newid).getReadGroupID().getSegmentPlacementCode();
	if(spc != ReadGroupLib::SPLACE_UNKNOWN
	   && spc != ReadGroupLib::SPLACE_SU
	   && templateguess.splace_seen != ReadGroupLib::SPLACE_UNKNOWN
	   && templateguess.splace_seen != ReadGroupLib::SPLACE_SU){
#ifndef PUBLICQUIET
	  cout << " spl";
	  cout.flush();
#endif
	  // if reads overlap, stay cool and don't check as it could be that sequencing errors and clipping
	  //  muddy the picture
	  // That is: the +10 bp criterion of above is not taken into account here, only the non-overlapping part

	  bool hasplacementerror=false;
	  if(leftry<rightrx){
#ifndef PUBLICQUIET
	    cout << " chk";
	    cout.flush();
#endif
	    if(spc!=templateguess.splace_seen) hasplacementerror=true;
	  }

	  if(hasplacementerror){
	    errstat.code=ESEGMENTPLACEMENT;

	    if(xcut<0) xcut=0;
	    if(static_cast<uint32>(ycut)>CON_counts.size()){
	      ycut=CON_counts.size();
	    }
	    priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);
	    FUNCEND();
	    return;
	  }
	}
      }

      // ** Check for template size **
      if(CON_readpool->getRead(newid).getInsizeFrom() >= 0
	 || CON_readpool->getRead(newid).getInsizeTo() >= 0){
#ifndef PUBLICQUIET
	cout << " dist";
	cout.flush();
#endif

	// allow 15% error in calculated insert size
	int32 aisp10=actinsertsize+actinsertsize*15/100;
	int32 aism10=actinsertsize-actinsertsize*15/100;

	int32 tif=tppcrI->getInsizeFrom();
	int32 tit=tppcrI->getInsizeTo();

	if(tppcrI->getInsizeFrom() >=0 && aisp10 < tif){
	// distance too small
#ifndef PUBLICQUIET
	  cout << "templ too small: " << tif << " min allowed, got " << aism10 << "-" << aisp10;
#endif
	  havematchingtemplatepartner=false;
	  if(!CON_readpool->getRead(newid).getReadGroupID().getTSInfoOnly()){
	    errstat.code=ETEMPLATESIZELT;
	  }
	}
	if(errstat.code == ENOERROR && tppcrI->getInsizeTo() >=0 && aism10 > tit){
	  // distance too big
#ifndef PUBLICQUIET
	  cout << "templ too big: " << tit << " max allowed, got " << aism10 << "-" << aisp10;
#endif
	  havematchingtemplatepartner=false;
	  if(!CON_readpool->getRead(newid).getReadGroupID().getTSInfoOnly()){
	    errstat.code=ETEMPLATESIZEGT;
	  }
	}
	if(errstat.code != ENOERROR){
	  if(xcut<0) xcut=0;
	  if(static_cast<uint32>(ycut)>CON_counts.size()){
	    ycut=CON_counts.size();
	  }
	  priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);
	  FUNCEND();
	  return;
	}
      }
    }
#ifndef PUBLICQUIET
    cout << " done\n";
    cout.flush();
#endif
  }

  CON_us_steps[USCLO_TEMPL1]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);


  bool doneit;
  bool maymapthisread;
  uint32 foolsafe=0;


#ifdef ALIGNCHECK
  Align checkbla(CON_miraparams);
#endif

  do{
    CEBUG("xcut: " << xcut << "\t");
    CEBUG("ycut: " << ycut << endl);
    CEBUG("CON_counts size: " << CON_counts.size() << endl);

    int32 eoffset=0;
    if(xcut<0){
      eoffset=xcut;
      xcut=0;
    }
    if(static_cast<uint32>(ycut)>CON_counts.size()){
      ycut=CON_counts.size()+1;
    }

    CEBUG("xcut: " << xcut << "\t");
    CEBUG("ycut: " << ycut << endl);
    CEBUG("eoffset: " << eoffset << "\t");

    // The following should happen only very, very rarely
    // Normally almost impossible, but maybe triggered in projects where
    //  reads get edited over and over again
    if(xcut> static_cast<int32>(CON_counts.size())
       || ycut < 0
       || ycut <= xcut){
#ifndef PUBLICQUIET
      cout << "rej: no align found (bounds 2)\t";
      cout.flush();
#endif
      errstat.code=ENOALIGN;
      errstat.reads_affected.push_back(refid);

      CON_us_steps[USCLO_XCUT]+=diffsuseconds(us_start);

      FUNCEND();
      return;
    }

    // if makeTmpConsenus() returns true, a N or X was encountered in the consensus
    //  and we may not map this read
    maymapthisread=!(makeTmpConsensus(xcut, ycut, CON_tmpcons_from_backbone, true));

    if(CON_2tmpcons.size()==0){
      // This should never happen here (the above checks should have made sure of that)
      // If it does ... we'll simply make it easy:
      //  dump out error to log and continue, ignoring it the best we can
      cout << "Sheeesh, length of temporary consensus is 0? Error, but continuing as this probably will not affect contig-building anyway."
	   << "\nxcut: " << xcut
	   << "\nycut: " << ycut
	   << "\neoffset: " << eoffset
	   << endl;
      madsl.clear();

      errstat.code=EUNSPECIFIED;
      errstat.reads_affected.push_back(refid);

      // rather not do this here ... we really do not know what xcut and ycut are
      //
      //// this can throw ... if it does, things are really, really, really botched
      //try {
      //	if(CON_readpool->getRead(refid).isRail())  {
      //	  getRailsAsReadsAffected(refid, errstat.reads_affected, xcut, ycut);
      //	}else{
      //	  getReadORPIDsAtContigPosition(errstat.reads_affected, xcut, ycut);
      //	}
      //}
      //catch (Notify n) {
      //	n.gravity=Notify::WARNING;
      //	n.handleError("internal to addRead_wrapped()");
      //	errstat.reads_affected.clear();
      //	errstat.reads_affected.push_back(refid);
      //}

      CON_us_steps[USCLO_XCUT]+=diffsuseconds(us_start);

      FUNCEND();
      return;
    }

    try{
      if(direction_frnid>=0){
	if(direction_refid>=0){
	  CEBUG("C1 align.\n");
#ifdef EXTRATRACEHELP
	  cout << "Read sequence:\n";
	  for(uint32 qwer=0; qwer<CON_readpool->getRead(newid).getLenClippedSeq();qwer++){
	    cout << CON_readpool->getRead(newid).getClippedSeqAsChar()[qwer];
	  }
#endif
	  aligncache[newreadseqtype].acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, 1, true, eoffset);
#ifdef ALIGNCHECK
	  checkbla.acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, 1);
#endif
	}else{
	  CEBUG("C2 align.\n");
#ifdef EXTRATRACEHELP
	  cout << "Read sequence:\n";
	  for(uint32 qwer=0; qwer<CON_readpool->getRead(newid).getLenClippedSeq();qwer++){
	    cout << CON_readpool->getRead(newid).getClippedComplementSeqAsChar()[qwer];
	  }
#endif
	  aligncache[newreadseqtype].acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedComplementSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, -1, true, eoffset);

	  //if(xcut == 468 && ycut == 1251) {
	  //  cout << "Read sequence:\n";
	  //  for(uint32 qwer=0; qwer<CON_readpool->getRead(newid).getLenClippedSeq();qwer++){
	  //    cout << CON_readpool->getRead(newid).getClippedComplementSeqAsChar()[qwer];
	  //  }
	  //  cout << endl;
	  //}


#ifdef ALIGNCHECK
	  checkbla.acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedComplementSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, -1);
#endif
	}
      }else{
	if(direction_refid>=0){
	  CEBUG("C3 align.\n");
#ifdef EXTRATRACEHELP
	  cout << "Read sequence:\n";
	  for(uint32 qwer=0; qwer<CON_readpool->getRead(newid).getLenClippedSeq();qwer++){
	    cout << CON_readpool->getRead(newid).getClippedSeqAsChar()[qwer];
	  }
#endif
	  aligncache[newreadseqtype].acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedComplementSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, -1, true, eoffset);
#ifdef ALIGNCHECK
	  checkbla.acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedComplementSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, -1);
#endif
	}else{
	  CEBUG("C4 align.\n");
#ifdef EXTRATRACEHELP
	  cout << "Read sequence:\n";
	  for(uint32 qwer=0; qwer<CON_readpool->getRead(newid).getLenClippedSeq();qwer++){
	    cout << CON_readpool->getRead(newid).getClippedSeqAsChar()[qwer];
	  }
#endif
	  aligncache[newreadseqtype].acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, 1, true, eoffset);
#ifdef ALIGNCHECK
	  checkbla.acquireSequences(
	    CON_2tmpcons.c_str(),
	    CON_2tmpcons.size(),
	    CON_readpool->getRead(newid).getClippedSeqAsChar(),
	    CON_readpool->getRead(newid).getLenClippedSeq(),
	    -1, newid, 1, 1);
#endif
	}
      }

#ifdef EXTRATRACEHELP
      cout << "\ntmpcons sequence:\n";
      for(uint32 qwer=0, cpl=0; qwer<CON_2tmpcons.size(); qwer++, cpl++){
	if(cpl==60) {
	  cout << endl;
	  cpl=0;
	}
	cout << CON_2tmpcons[qwer];
      }
      cout << endl;
#endif

      CEBUG("Done acquiring.\n");
      madsl.clear();
      CEBUG("madsl cleared.\n");

      bool enforce_clean_ends=(*CON_miraparams)[CON_readpool->getRead(newid).getSequencingType()].getAlignParams().ads_enforce_clean_ends;
      // if the reference read is a rail or backbone, do not
      //  enforce clean ends! (to find SNPs!)
      if(CON_readpool->getRead(refid).isBackbone()
	 || CON_readpool->getRead(refid).isRail()) enforce_clean_ends=false;

      // and override with specialsradd condition
      if (CON_ssrc_cleanoverlapends > 0 ) enforce_clean_ends = true;

      /* BaCh 29.03.2009
	  enforcing clean ends here is ... dangerous.
          reads is a contig that contain a sequencing error
          otherwise prevent the correct extension of the contig

                  |    .    |    .    |    .    |
          ID1:-1  GGGCATTGTCTGCCACCTCTAACTCTACCTAA
          ID2:500 GGGCATTGTCTGCCAGCTCTAAC
          480-                   X

	 where the -1 contig is a single read with a sequencing error
	 near the end (C instead of G)

	 In general, enforce_clean_ends should be off in de-novo (genome, EST),
	 but on if one wants to do a clustering of transcript contigs.

	 USE AT OWN RISK!
      */

      bool affinegapscore=false;
      if(CON_readpool->getRead(refid).isBackbone()
	 || CON_readpool->getRead(refid).isRail()) affinegapscore = true;

      aligncache[newreadseqtype].setEnforceCleanEnds(enforce_clean_ends);
      aligncache[newreadseqtype].setAffineGapScore(affinegapscore);
      aligncache[newreadseqtype].fullAlign(&madsl);
    }
    catch(Notify n){
      cout << std::dec << "Ouch ... error in alignment detected.\n";
      cout << "xcut: " << xcut << "\t";
      cout << "ycut: " << ycut << endl;
      cout << "eoffset: " << eoffset << "\n";

      cout << "dir_frnid: " << direction_frnid;
      cout << "\tdir_refid: " << direction_refid<< endl;
      cout << "Offset Refid:" << offsetrefid<<endl;

      Read::setCoutType(Read::AS_TEXTSHORT);
      cout << CON_readpool->getRead(newid);

      aligncache[newreadseqtype].coutWhatWasGiven();

      // throw again so that addRead() catches and we get info to replay
      throw Notify(n);
    }

    CEBUG("Done full align.\n");


//#ifdef EXTRATRACEHELP
//    if(CLASS_debug_counter==112){
//      setCoutType(AS_TEXT);
//      cout << *this;
//      setCoutType(AS_DEBUG);
//      cout << *this;
//    }
//#endif


#ifdef ALIGNCHECK
    CEBUG("ALC");
    if(madsl.empty()){

      // scrap this, AlignedDualSeqFacts does not drag this along
      //  (saving memory)
      //if(initialadsf->getOverlap()>50) {
      //	cout << "Missed badly?\n";
      //	cout << *initialadsf;
      //	cout << "Tmp cons.: \n" << CON_tmpcons << endl;
      //	cout << *this;
      //}

      std::list<AlignedDualSeq> tadsl;
      checkbla.fullAlign(&tadsl);
      if(!tadsl.empty()){
	CEBUG("Waaaaah! BSW failed!" << endl);

	//CEBUG("I had this contig:\n"<<*this);
	CEBUG("Had to align this:\n");
	  CEBUG("Refid: " << refid);
	  CEBUG("\tNewid: " << newid);
	  CEBUG("dir_frnid: " << direction_frnid);
	  CEBUG("\tdir_refid: " << direction_refid<< endl;);
	  CEBUG("Offset Refid:" << offsetrefid<<endl;);
	  CEBUG("xcut: " << xcut << "\t");
	  CEBUG("ycut: " << ycut << endl);
	  CEBUG("eoffset: " << eoffset << "\t");

#ifdef CEBUGFLAG
	  CEBUG(" ----------------------------------------------------- \n");
	  CEBUG("# solutions found: "<< tadsl.size() << endl);

	  for(auto & adse : tadsl){
	    CEBUG(adse<<endl<<endl)
	  }

	  CEBUG(" ----------------------------------------------------- \n");
#endif
      }
    }
#endif

#ifdef CEBUGFLAG
    CEBUG(" ----------------------------------------------------- \n");
    CEBUG("# solutions found: "<< madsl.size() << endl);


    for(auto & adse : madsl){
      CEBUG(adse);
    }

    CEBUG(" ----------------------------------------------------- \n");
#endif

    if(madsl.empty()){
      //    throw Notify(Notify::INTERNAL, THISFUNC, "No solution?!? At least one expected.\n");

      CON_us_steps[USCLO_SWALIGN]+=diffsuseconds(us_start);

#ifndef PUBLICQUIET
      cout << "rej: no align found";
      if(aligncache[newreadseqtype].wasBandHit()) cout << " BH";
      cout << '\t';
      //cout << "Ov: " << initialadsf->getOverlap() << "\t";
      //cout << "ESR: " << static_cast<uint16>(initialadsf->getScoreRatio()) << "\t";
      //cout << "S: " << initialadsf->getScore() << "\t";
      //cout << "ES: " << initialadsf->getExpectedScore() << "\t";
      cout.flush();
#endif
      errstat.code=ENOALIGN;
      priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);

      FUNCEND();
      return;
    }

    // TODO: was tun?
    // z.Zt. erst einmal den besten nehmen. wird wahrsch. der richtige sein.
    // darauf achten, dass, wenn möglich, offset==0 ist
    {
      int32 bestratio=0;
      int32 bestweight=0;
      int32 offset=10000;

      auto Ibest=madsl.cbegin();
      for(adsI=Ibest; adsI!=madsl.cend(); ++adsI){
	if(adsI->getScoreRatio()>=bestratio){
	  CEBUG(offset << "\t" << adsI->getOffsetInAlignment(newid));
	  if(adsI->getScoreRatio()>=bestratio){
	    //|| (offset>0 && I->getOffsetInAlignment(newid)<offset)){
	    if(adsI->getWeight() > bestweight) {
	      offset=adsI->getOffsetInAlignment(newid);
	      bestratio=adsI->getScoreRatio();
	      bestweight=adsI->getWeight();
	      Ibest=adsI;
	    }
	  }
	}
      }
      adsI=Ibest;
    }

    CEBUG(" ----------------------------------------------------- \n");
    CEBUG("# solution chosen: " << endl);
    CEBUG(*adsI);
    CEBUG(" ----------------------------------------------------- \n");

    doneit=true;

    ++foolsafe;

    // should the alignment not be perfect, calculate new xcut (ycut)
    //  coordinates
    // IF we're in the last foolsaferedone loop, make sure we still go another
    //  loop
    if(xcut >0){
      if(adsI->getOffsetInAlignment(-1)!=0){
	CEBUG("Redo because of wrong left cut. Too far right.\n");
	doneit=false;
	if(foolsafe==5) foolsafe--;
	xcut-=adsI->getOffsetInAlignment(-1)+7;
      }
      //else if(I->getOffsetInAlignment(newid) > 5){
      //	// TODO: check whether we could accept this regardless of
      //	//  offset. check: right offsets.
      //	CEBUG("Redo because of wrong left cut. Too far left.\n");
      //	doneit=false;
      //	if(foolsafe<5) xcut+=I->getOffsetInAlignment(newid);
      //}
    }

    if(adsI->getRightOffsetInAlignment(-1)!=0
       && ycut < CON_counts.size()+1){
      CEBUG("Redo because of wrong right cut. Not taken enough cons.\n");
      doneit=false;
      if(foolsafe==5) foolsafe--;
      ycut+=adsI->getRightOffsetInAlignment(-1)+10;
    }
    if(aligncache[newreadseqtype].wasBandHit()){
      CEBUG("Redo because of band hit in banded SW.\n");
      const_cast<align_parameters &>(rt_params.getAlignParams()).al_kmin=200;
      const_cast<align_parameters &>(rt_params.getAlignParams()).al_kmin=400;
      const_cast<align_parameters &>(rt_params.getAlignParams()).al_kpercent=80;
      doneit=false;
    }

    if(foolsafe>1){
#ifndef PUBLICQUIET
      cout << "fsc (" << foolsafe << ")\t";
      cout.flush();
#endif
    }
    if(foolsafe>=5){
#ifndef PUBLICQUIET
      cout << "already 5 iterations to find optimum ads and still not found, using suboptimum.\n";
      cout.flush();

      CEBUG(*adsI);
#endif
      doneit=true;
    }

  }while(doneit==false);

  CON_us_steps[USCLO_SWALIGN]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  // need this, because if foolsafe kicks in, xcut may be <0 and ycut > size
  if(xcut<0){
    xcut=0;
  }
  if(static_cast<uint32>(ycut)>CON_counts.size()){
    ycut=CON_counts.size();
  }


  CEBUG("CON_specialsraddconditions: " << CON_specialsraddconditions << endl);
  CEBUG("newreadseqtype: " << ReadGroupLib::getNameOfSequencingType(newreadseqtype) << endl);

#ifndef PUBLICQUIET
  cout << "RRRR 5p: " << adsI->get5pLenContiguousMatch(newid) << "\t3p: " << adsI->get3pLenContiguousMatch(newid) << endl;
#endif

  // check for short reads whether we are violating against special rules for mismatch numbers
  if(CON_specialsraddconditions
    && (newreadseqtype==ReadGroupLib::SEQTYPE_SOLEXA
	|| newreadseqtype==ReadGroupLib::SEQTYPE_ABISOLID)){
    CEBUG("CON_ssrc_maxtotalerrors: " << CON_ssrc_maxtotalerrors << endl);
    if(CON_ssrc_maxtotalerrors >=0){
      if(static_cast<int32>(adsI->getNumMismatches())
	 + static_cast<int32>(adsI->getNumGaps()) > CON_ssrc_maxtotalerrors){
#ifndef PUBLICQUIET
	cout << "saf:mte " << (static_cast<int32>(adsI->getNumMismatches()) + static_cast<int32>(adsI->getNumGaps())) << ">" << CON_ssrc_maxtotalerrors;
#endif
	errstat.code=ESPECIALSRADDFAIL;
      }
    }
    CEBUG("CON_ssrc_maxmismatches: " << CON_ssrc_maxmismatches << endl);
    if(CON_ssrc_maxmismatches >=0){
      if(static_cast<int32>(adsI->getNumMismatches()) > CON_ssrc_maxmismatches){
#ifndef PUBLICQUIET
	cout << " saf:mnm " << adsI->getNumMismatches() << ">" << CON_ssrc_maxmismatches;
#endif
	errstat.code=ESPECIALSRADDFAIL;
      }
    }
    CEBUG("CON_ssrc_maxgaps: " << CON_ssrc_maxgaps << endl);
    if(CON_ssrc_maxgaps >=0){
      if(static_cast<int32>(adsI->getNumGaps()) > CON_ssrc_maxgaps){
#ifndef PUBLICQUIET
	cout << " saf:mng " << static_cast<int32>(adsI->getNumGaps()) << ">" << CON_ssrc_maxgaps;
#endif
	errstat.code=ESPECIALSRADDFAIL;
      }
    }
    CEBUG("CON_ssrc_cleanoverlapends: " << CON_ssrc_cleanoverlapends << endl);
    if(CON_ssrc_cleanoverlapends >=0){
      if(static_cast<int32>(adsI->get5pLenContiguousMatch(newid)) < CON_ssrc_cleanoverlapends
	 || static_cast<int32>(adsI->get3pLenContiguousMatch(newid)) < CON_ssrc_cleanoverlapends
	){
#ifndef PUBLICQUIET
	cout << " saf:coe " << static_cast<int32>(adsI->get5pLenContiguousMatch(newid))
	     << "," << static_cast<int32>(adsI->get3pLenContiguousMatch(newid)) << "<" << CON_ssrc_cleanoverlapends;
#endif
	errstat.code=ESPECIALSRADDFAIL;
      }
    }
    if(errstat.code==ESPECIALSRADDFAIL){
#ifndef PUBLICQUIET
      cout << " specialsradd failed\t";
      //cout << *adsI;
#endif

      priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);

      FUNCEND();
      return;
    }
//    else{
//      cout << "accepting this: " << endl;
//      cout << "maxtotal: " << CON_ssrc_maxtotalerrors;
//      cout << "maxmis: " << CON_ssrc_maxmismatches;
//      cout << "maxgap: " << CON_ssrc_maxgaps;
//      cout << *adsI;
//    }
  }

#ifndef PUBLICQUIET
    cout << "\tASR: " << static_cast<uint16>(adsI->getScoreRatio()) << '\t';
    //cout << *I;
#endif

  // First check on rodirs: relaxed parameters (rodirs*2)
  // Reason: we haven't checked template restriction yet
  //  if we lateron see that we do not have a matching template,
  //  then we'll use the 'stricter' rodirs and cmrs value
  if(initialadsf->getScoreRatio() >
     adsI->getScoreRatio()+(2*rt_params.getContigParams().con_reject_on_drop_in_relscore)){
    // REMOVEME
#ifndef PUBLICQUIET
    cout << "Dead end (even lax)\t";
    cout << "ESR: " << static_cast<uint16>(initialadsf->getScoreRatio());
    cout << "\tASR: " << static_cast<uint16>(adsI->getScoreRatio()) << '\t';
    //cout << *adsI;
#endif
    //    cout << *adsI;

    errstat.code=EDROPINRELSCORE;
    //errstat.reads_affected.push_back(refid);

    // Diese entscheidung hier ist ... gefährlich, wenn nicht gar falsch,
    //  denn es kann sein, dass ein potentiell guter match auch dieselben
    //  Fehler hatte wie dieses Align.
    // falsch!? ist legitim? diese bans werden ja nicht auf permbans uebertragen,
    //  weshalb es fuer _diesen_ Contig ja stimmen wuerde.

    priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);

    FUNCEND();
    return;
  }


  // compute whether the contig will grow (left or right)

  // check left side
  int32 expect_growleft=0;
  if(adsI->getOffsetInAlignment(newid)==0 && adsI->getOffsetInAlignment(-1)!=0){
    expect_growleft=-(xcut-adsI->getOffsetInAlignment(-1));
    if(expect_growleft<0) expect_growleft=0;
  }
  int32 expect_growright=0;
  if(ycut==CON_counts.size()
     && adsI->getRightOffsetInAlignment(newid)==0
     && adsI->getRightOffsetInAlignment(-1)!=0){
    expect_growright=adsI->getRightOffsetInAlignment(-1);
  }

  if(forcegrow>0
     && std::max(expect_growleft, expect_growright) < forcegrow){
    errstat.code=EFORCEDGROWTHNOTREACHED;

    priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);

    FUNCEND();
    return;
  }

  // look whether the contig may grow
  // no matter what the pathfinder says: if the new read has a mate in this contig and passed
  //  all template tests, then the contig will allow growth
  if(!havematchingtemplatepartner
     && forcegrow<0
     && (expect_growleft > 0
	 || expect_growright > 0)){
    errstat.code=EGROWTHNOTALLOWED;
    priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);
    FUNCEND();
    return;
  }

  // now check whether we're aligning in an area where we already
  //  reached max coverage
  //
  // New 14.12.2008: if contig grows, skip that as we will very
  //  probably be in a multicopy environment and if we add
  //  so much coverage that maxcoverage allowed is reached,
  //  the building stops (which is bad).
  // TODO: do better, perhaps before the alignment?
  //  furthermore, change reduceSKim() and Pathfinder to
  //  select links where multicopy is involved to have at least
  //  20 to 30 bases extension!
  // New 03.01.2009: scrap the above, leads to heavy overcompression
  // New 29.04.2011:
  // taken out unconditional check:
  //  old routine always gave back "true" anyway except if coverage
  //  reached 16384
  //  with new 32 bit counters, will never be really reached anyway

  if(newid_ismulticopy
     && !havematchingtemplatepartner
//     && expect_growleft==0
//     && expect_growright==0
     && (*CON_miraparams)[0].getAssemblyParams().as_uniform_read_distribution
     && !checkFreeCoverageForAddingRead(newreadseqtype, xcut, ycut)){
    errstat.code=EMAXCOVERAGEREACHED;
    priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);
    FUNCEND();
    return;
  }

  CON_us_steps[USCLO_PREINSCHK]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  auto coveragemultiplier=CON_readpool->getRead(newid).getDigiNormMultiplier();

  if(CON_readpool->getRead(refid).isRail()
     && CON_mergenewsrreads
     && newreadseqtype==ReadGroupLib::SEQTYPE_SOLEXA
     && rt_params.getContigParams().con_mergeshortreads){
    // forcemerge: for reads not mapping 100%
    bool forcemerge=false;
    // first force merge: reads < 100% but <= -CO:msrme errors
    if(CON_hasforcemergeareas){
      //// for testing
      //forcemerge=true;

      // TODO: is using xcut & ycut OK? Probably yes.
      //  (especially ycut is not 100% on the end of the alignment!)
      auto ccI=CON_counts.begin();
      std::advance(ccI,xcut);
      for(uint32 ii=xcut; ii<ycut; ++ii, ++ccI){
	if(ccI->forcemergearea) {
	  forcemerge=true;
	  break;
	}
      }
    }

    if(adsI->getNumMismatches()+adsI->getNumGaps() > rt_params.getContigParams().con_msr_maxerrors){
      maymapthisread=false;
    }

    // now check whether we can merge this read
    // if we want to map, we need 100% score ratio
    CEBUG("CON_mergenewsrreads: " << CON_mergenewsrreads << endl);
    CEBUG("maymapthisread: " << maymapthisread << endl);
    CEBUG("adsI->getScoreRatio() == 100: " << (adsI->getScoreRatio() == 100) << endl);

    if(forcemerge
       || (maymapthisread && adsI->getScoreRatio() == 100)){

      CEBUG("Try map.\n");

      bool canmap=true;

      // sometimes the 100% matches have some weirdnesses against N's
      // ADS should nowadays take care of this, but just to be sure:
      if(!forcemerge){
	if(canmap &&
	   (adsI->getNumMismatches()>0
	    || adsI->getNumGaps()>0)){
#ifndef PUBLICQUIET
	  cout << " has mismatches";
#endif
	  canmap=false;
	}
      }

      // after this point, "canmap" is the only deciding variable left, "forcemerge" not looked
      //  at anymore!

      // furthermore, we may not be extending the contig left or right

      if(expect_growleft>0) {
	canmap=false;
#ifndef PUBLICQUIET
	cout << " growl";
#endif
      }
      if(expect_growright>0) {
	canmap=false;
#ifndef PUBLICQUIET
	cout << " growr";
#endif
      }

// TODO: check whether can be taken out as should be taken care of by makeTmpConsensus()
//
//      // ok, consensus does not grow ... but it might be in a part
//      //  that initially was not in the contig.
//      if(canmap &&
//	 (CON_counts[xcut].getBBChar()=='@'
//	  || CON_counts[ycut-1].getBBChar()=='@')) {
//#ifndef PUBLICQUIET
////      cout << "xcut: " << xcut << "\tycut: " << ycut << endl;
////      dumpAsDebug(cout);
//	cout << " ingrown";
//#endif
//	canmap=false;
//      }

      //cout << "\ncanmap2: " << canmap << endl;

      // maybe we'd like to have reads at contig ends not mapped (e.g.
      //  for scaffolding)
      if(canmap
	 && rt_params.getContigParams().con_msr_keependsunmapped!=0
	 && CON_mpindex_msrkceu_left>=0  // are we using markerpositions anyway?
	 && CON_readpool->getRead(newid).getTemplatePartnerID() >= 0){

	// this method is susceptible to be less than optimal for data with lots of gap columns
	//  (high coverage 454 & Ion)
	// TODO: see whether I need to take countermeasures like if in 2*dist, then count
	//  gap-columns to get a much better true distance (left and right)

	int32 dist=rt_params.getContigParams().con_msr_keependsunmapped;
	if(dist<0){
	  dist=CON_readpool->getRead(newid).getInsizeTo();
	}

	//cout << "mpl: " << CON_markerpositions[CON_mpindex_msrkceu_left]
	//     << "\tmpr: " << CON_markerpositions[CON_mpindex_msrkceu_right]
	//     << "\nxcut: " << xcut
	//     << "\tycut: " << ycut
	//     << "\tdist: " << dist
	//     << "\nl: " << dist+CON_markerpositions[CON_mpindex_msrkceu_left]
	//     << "\tl: " << CON_markerpositions[CON_mpindex_msrkceu_right]-dist
	//  ;

	// if no distance or when read has no paired-end partner, we could map
	//  anyway and would not need to check further
	if(dist>0){
	  // check whether we are in boundaries, if yes, do not map
	  if(xcut <= dist+CON_markerpositions[CON_mpindex_msrkceu_left]
	     || ycut >= CON_markerpositions[CON_mpindex_msrkceu_right]-dist){
	    canmap=false;
	  }
	}

	//cout << "\ncanmap2: " << canmap << endl;
      }

      // map if possible, else the read will be normally added
      if(canmap){
#ifndef PUBLICQUIET
	cout << " mapping ...";
	cout.flush();
#endif

	bool wasmapped=true;
	try{
	  wasmapped=insertMappedReadInContig(*adsI, newreadseqtype, xcut,
					     direction_frnid, direction_newid_incontig,
					     coveragemultiplier, forcemerge);
	}
	catch (Notify n){
	  cout << "Uh oh ... not good\n";
	  cout << "forcemerge " << forcemerge << endl;
	  cout << "CON_mergenewsrreads " << CON_mergenewsrreads << endl;
	  cout << "maymapthisread " << maymapthisread << endl;
	  cout << "newreadseqtype " << static_cast<uint16>(newreadseqtype);
	  cout << "rt_params.getContigParams().con_mergeshortreads " << rt_params.getContigParams().con_mergeshortreads << endl;
	  cout << "adsI->getScoreRatio() " << static_cast<uint16>(adsI->getScoreRatio()) << endl;
	  cout << "canmap " << canmap << endl;
	  cout << *adsI;
	  // throw again so that addRead() catches and we get info to replay
	  throw Notify(n);
	}

	CON_us_steps[USCLO_INSCONM]+=diffsuseconds(us_start);

	if(wasmapped){
	  // count this in the readsperstrain statistics
	  CON_readsperstrain[CON_readpool->getRead(newid).getStrainID()]+=coveragemultiplier;

	  if(CON_readsperreadgroup.size() < ReadGroupLib::getNumReadGroups()){
	    CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);
	  }
	  CON_readsperreadgroup[CON_readpool->getRead(newid).getReadGroupID().getLibId()]+=coveragemultiplier;
	  return;
	}else{
#ifndef PUBLICQUIET
	  cout << " failedmap";
	  cout.flush();
#endif
	}
      }else{
#ifndef PUBLICQUIET
	cout << " cantmap";
	cout.flush();
#endif
      }
    }
  }

  //cout << *I;

  // Template handling, part 2

  if(!havematchingtemplatepartner){
#ifndef PUBLICQUIET
    cout << " tmplhand2...";
    cout.flush();
#endif

    // Second check on rodirs and cmrs: strict parameters (rodirs)
    // Reason: if there was no matching tpartner, apply strict parameter

    auto checkmrs=rt_params.getContigParams().con_min_relscore;
    if(checkmrs<0) checkmrs=rt_params.getAlignParams().al_min_relscore;
    if(initialadsf->getScoreRatio() > adsI->getScoreRatio()+rt_params.getContigParams().con_reject_on_drop_in_relscore
       || initialadsf->getScoreRatio() < checkmrs){
      // REMOVEME

#ifndef PUBLICQUIET
      cout << "Dead end (strict)\t";
      cout << "ESR: " << static_cast<uint16>(initialadsf->getScoreRatio());
      cout << "\tASR: " << static_cast<uint16>(adsI->getScoreRatio()) << '\t';
      // cout << *I;
#endif

      errstat.code=EDROPINRELSCORE;
      if(initialadsf->getScoreRatio() < checkmrs){
	errstat.code=EALIGNREJECTMRS;
      }

      //errstat.reads_affected.push_back(refid);

      // Diese entscheidung hier ist ... gefährlich, wenn nicht gar falsch,
      //  denn es kann sein, dass ein potentiell guter match auch dieselben
      //  Fehler hatte wie dieses Align.
      // falsch!? ist legitim? diese bans werden ja nicht auf permbans uebertragen,
      //  weshalb es fuer _diesen_ Contig ja stimmen wuerde.

      priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);

      FUNCEND();
      return;
    }
  }




#ifndef PUBLICQUIET
  cout << " inscon...";
  cout.flush();
#endif

  gettimeofday(&us_start,nullptr);
  auto nprI=insertReadInContig(*adsI, xcut,direction_frnid, direction_refid, coveragemultiplier);
  CON_us_steps[USCLO_INSCON]+=diffsuseconds(us_start);

  // count this in the readsperstrain statistics
  CON_readsperstrain[CON_readpool->getRead(newid).getStrainID()]+=coveragemultiplier;

  if(CON_readsperreadgroup.size() < ReadGroupLib::getNumReadGroups()){
    CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);
  }
  CON_readsperreadgroup[CON_readpool->getRead(newid).getReadGroupID().getLibId()]+=coveragemultiplier;


#ifndef PUBLICQUIET
  cout << "done";
  cout.flush();
#endif


#ifndef PUBLICQUIET
  cout << " updbasloc...";
  cout.flush();
#endif

  // Put base locks in CON_counts that this read produces
  gettimeofday(&us_start,nullptr);
  updateBaseLocks(nprI, true);
  CON_us_steps[USCLO_UPDBLOCKS]+=diffsuseconds(us_start);

#ifndef PUBLICQUIET
  cout << "done";
  cout.flush();
#endif

#ifndef PUBLICQUIET
  cout << " chkcon...";
  cout.flush();
#endif
  paranoiaBUGSTAT(checkContig());
#ifndef PUBLICQUIET
  cout << "done";
  cout.flush();
#endif


#ifndef PUBLICQUIET
  cout << " ansrmbzone...";
  cout.flush();
#endif

  {
    gettimeofday(&us_start,nullptr);
    auto arz=analyseRMBZones(nprI);
    CON_us_steps[USCLO_ANRMBZ]+=diffsuseconds(us_start);

    if(arz == true){
#ifndef PUBLICQUIET
      cout << "rej: srmb zone\t";
#endif
      //if(CON_reads.back().read.getName()=="GBIBI26TF") {
      //  cout << "\nMUSTSAVE!\n";
      //  saveAsGAP4DA("error_out.gap4da", cout);
      //  exit(0);
      //}
      // remove the read from the contig as we don't want it

      gettimeofday(&us_start,nullptr);
      deleteRead(nprI);
      CON_us_steps[USCLO_DELREAD]+=diffsuseconds(us_start);

      errstat.code=ESRMBMISMATCH;
      priv_arw_getReadsAffected(refid, newid, errstat.reads_affected, xcut, ycut);

      FUNCEND();
      return;
    }
  }

#ifndef PUBLICQUIET
  cout << "done";
  cout.flush();
#endif

#ifndef PUBLICQUIET
  cout << " andngrzone...";
  cout.flush();
#endif

  {
////++++//// Skip this as currently not really used by anyone
////++++////    BUGIFTHROW(true,"need redo 4 for PlacedContigReads");
////++++////    if(analyseDangerZones(CON_reads.back()) == true){
////++++////      // TODO: in ALUS & co strengere Kriterien.
////++++////      // remove the read from the contig as we don't want it
////++++////
////++++////#ifndef PUBLICQUIET
////++++////      cout << "rej: danger zone\t";
////++++////#endif
////++++////
////++++////#ifdef CLOCK_STEPS
////++++////      gettimeofday(&us_start,nullptr);
////++++////#endif
////++++////      deleteRead(newid);
////++++////#ifdef CLOCK_STEPS
////++++////      CON_us_steps[USCLO_DELREAD]=diffsuseconds(us_start);
////++++////#endif
////++++////
////++++////      errstat.code=EDANGERZONE;
////++++////
////++++////#ifdef CLOCK_STEPS
////++++////      gettimeofday(&us_start,nullptr);
////++++////#endif
////++++////      if(CON_readpool->getRead(refid).isRail())  {
////++++////	getRailsAsReadsAffected(refid, errstat.reads_affected, xcut, ycut);
////++++////      }else{
////++++////	getReadORPIDsAtContigPosition(errstat.reads_affected, xcut, ycut);
////++++////      }
////++++////#ifdef CLOCK_STEPS
////++++////      CON_us_steps[USCLO_GRACP]=diffsuseconds(us_start);
////++++////#endif
////++++////
////++++////      FUNCEND();
////++++////      return;
////++++////    }
  }

#ifndef PUBLICQUIET
  cout << "done";
  cout.flush();
#endif

  if(CON_fixedconsseq.size()){
    nukeSTLContainer(CON_fixedconsseq);
    nukeSTLContainer(CON_fixedconsqual);
  }

  CON_contains_majority_digitallynormalised_reads=0;

  FUNCEND();
  return;
}
#define CEBUG(bla)
#define CEBUGF(bla)
#undef CEBUGFLAG

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_arw_getReadsAffected(const readid_t refid, const readid_t newid, std::vector<readid_t> & reads_affected, const int32 xcut, const int32 ycut)
{
  timeval us_start;
  gettimeofday(&us_start,nullptr);
  if(CON_readpool->getRead(refid).isRail())  {
    // in mapping this can get quite slow for large contigs with hundreds of thousands of reads
    // workaround: for "small" reads just push the refid.
    // This may lead to a second / third alignment attempt if the pathfinder has also hits to neighbouring
    //  rails, but that is still a lot faster that the call to getRailsAsReadsAffected()
    if(CON_readpool->getRead(newid).getLenClippedSeq()<600)  {
      reads_affected.push_back(refid);
    }else{
      getRailsAsReadsAffected(refid, reads_affected, xcut, ycut);
    }
  }else{
    getReadORPIDsAtContigPosition(reads_affected, xcut, ycut);
  }
  if(CON_us_steps.size()>USCLO_GRACP) CON_us_steps[USCLO_GRACP]+=diffsuseconds(us_start);
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::coutAddReadTimings()
{
  if(CON_us_steps.size()){
    cout << "\nccon timings: "
	 << "\ncct pre\t" << CON_us_steps[USCLO_PRE]
	 << "\ncct dir\t" << CON_us_steps[USCLO_DIR]
	 << "\ncct xcu\t" << CON_us_steps[USCLO_XCUT]
	 << "\ncct tmp\t" << CON_us_steps[USCLO_TEMPL1]
	 << "\ncct sw \t" << CON_us_steps[USCLO_SWALIGN]
	 << "\ncct pic\t" << CON_us_steps[USCLO_PREINSCHK]
	 << "\ncct ico\t" << CON_us_steps[USCLO_INSCON]
	 << "\ncct icm\t" << CON_us_steps[USCLO_INSCONM]
	 << "\ncct upd\t" << CON_us_steps[USCLO_UPDBLOCKS]
	 << "\ncct del\t" << CON_us_steps[USCLO_DELREAD]
	 << "\ncct rmz\t" << CON_us_steps[USCLO_ANRMBZ]
	 << "\ncct gcp\t" << CON_us_steps[USCLO_GRACP]
	 << '\n';
  }
  if(CON_us_steps_iric.size()){
    cout << "\nccon i timings (" << CON_track_numins << "): "
	 << "\nccit insglcc\t" << CON_us_steps_iric[USCLOIRIC_INSGLCC]
	 << "\nccit insglaro\t" << CON_us_steps_iric[USCLOIRIC_INSGLARO]
	 << "\nccit insglact\t" << CON_us_steps_iric[USCLOIRIC_INSGLACT]
	 << "\nccit insgltot\t" << CON_us_steps_iric[USCLOIRIC_INSGLTOT]
	 << "\nccit pr\t\t" << CON_us_steps_iric[USCLOIRIC_PR]
	 << "\nccit templ\t" << CON_us_steps_iric[USCLOIRIC_TEMPL]
	 << "\nccit index\t" << CON_us_steps_iric[USCLOIRIC_INDEX]
	 << "\nccit biglccins\t" << CON_us_steps_iric[USCLOIRIC_BIGLCCINS]
	 << "\nccit biglinterpol\t" << CON_us_steps_iric[USCLOIRIC_BIGLINTERPOL]
	 << "\nccit biglupdtags\t" << CON_us_steps_iric[USCLOIRIC_BIGLUPDTAGS]
	 << "\nccit biglfpcri\t" << CON_us_steps_iric[USCLOIRIC_BIGLFPCRI]
	 << "\nccit bigllgap\t" << CON_us_steps_iric[USCLOIRIC_BIGLLIGAPINREADS]
	 << "\nccit biglshiftread\t" << CON_us_steps_iric[USCLOIRIC_BIGLSHIFTREADS]
	 << "\nccit bigltot\t" << CON_us_steps_iric[USCLOIRIC_BIGLTOT]
	 << "\nccit insgr\t" << CON_us_steps_iric[USCLOIRIC_INSGR]
	 << "\nccit ucv\t" << CON_us_steps_iric[USCLOIRIC_UCV]
	 << "\n";
  }
  if(CON_us_steps_drfc.size()){
    cout << "\nccon d timings (" << CON_track_numdels << "): "
	 << "\nccdt ubl\t"   << std::setw(14) << CON_us_steps_drfc[USCLODRFC_UBL]
	 << "\nccdt ucv\t"   << std::setw(14) << CON_us_steps_drfc[USCLODRFC_UCV]
	 << "\nccdt itf\t"   << std::setw(14) << CON_us_steps_drfc[USCLODRFC_ITF]
	 << "\nccdt itb\t"   << std::setw(14) << CON_us_steps_drfc[USCLODRFC_ITB]
	 << "\nccdt ccef\t"  << std::setw(14) << CON_us_steps_drfc[USCLODRFC_CCEF]
	 << "\nccdt cceb\t"  << std::setw(14) << CON_us_steps_drfc[USCLODRFC_CCEB]
	 << "\nccdt sdt\t"   << std::setw(14) << CON_us_steps_drfc[USCLODRFC_SDT]
	 << "\nccdt sr\t\t"  << std::setw(14) << CON_us_steps_drfc[USCLODRFC_SR]
	 << "\nccdt smp\t"   << std::setw(14) << CON_us_steps_drfc[USCLODRFC_SMP]
	 << "\nccdt dt\t\t"  << std::setw(14) << CON_us_steps_drfc[USCLODRFC_DT]
	 << "\nccdt rr\t\t"  << std::setw(14) << CON_us_steps_drfc[USCLODRFC_RR]
	 << "\nccdt dsoc\t"  << std::setw(14) << CON_us_steps_drfc[USCLODRFC_DSOC]
	 << "\nccdt total\t" << std::setw(14) << CON_us_steps_drfc[USCLODRFC_TOTAL]
	 << "\n";
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::setContigCoverageTarget(std::vector<uint32> covtarget)
{
  FUNCSTART("void Contig::setContigCoverageTarget(std::vector<uint32)");

  CON_targetcoverageperst.clear();
  if(covtarget.empty()) return;

  BUGIFTHROW(covtarget.size()!=ReadGroupLib::SEQTYPE_END, "covtarget.size()!=ReadGroupLib::SEQTYPE_END ??");
  for(uint32 i=0; i<covtarget.size(); i++){
    CON_targetcoverageperst.push_back(static_cast<uint16>(covtarget[i]));
  }

  FUNCEND();
}


/*************************************************************************
 *
 * looks whether adding a read between the given positions would not
 *  exceed the maximum allowed coverage
 *
 * return true for "it's ok to add"
 * false for "maximum already reached"
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool Contig::checkFreeCoverageForAddingRead(const uint8 newreadseqtype, const int32 xcut, const int32 ycut)
{
  FUNCSTART("void Contig::bool Contig::checkFreeCoverageForAddingRead(const uint8 newreadseqtype, const int32 xcut, const int32 ycut)");

  CEBUG("\ncovcheck: st:" << static_cast<uint16>(newreadseqtype) << "\txcut: " << xcut << "\tycut: " << ycut << '\n');

  ccctype_t maxcovallowed=0;
  if(!CON_targetcoverageperst.empty()){
    ccctype_t tmp=static_cast<ccctype_t>(CON_targetcoverageperst[newreadseqtype]) * ((*CON_miraparams)[newreadseqtype]).getAssemblyParams().as_urd_cutoffmultiplier;
    maxcovallowed=std::max(tmp,maxcovallowed);
  }

  auto ccI=CON_counts.begin();
  std::advance(ccI,xcut);

  if(maxcovallowed<10) maxcovallowed=10;

  // The 16383 number stems from the 16 bit base counters in CON_concounts[]
  //  which use increments by 4
  // 1073741823 is for 32 bit counters (== 2^32 / 4 -1)
  for(int32 ci=xcut; ci<ycut && ccI!=CON_counts.end(); ++ccI, ++ci){
    CEBUG("ci: " << ci << '\t' << ccI->seqtype_cov[0] << '\t' << ccI->seqtype_cov[1] << '\t' << ccI->seqtype_cov[2] << '\t' << ccI->seqtype_cov[3] << "\tc: " << ccI->total_cov);
    if(ccI->seqtype_cov[newreadseqtype]>maxcovallowed
      || ccI->total_cov == 1073741823) {
      return false;
    }
    CEBUG('\n');
  }

  FUNCEND();
  return true;
}
//#define CEBUG(bla)




/*************************************************************************
 *
 * Meant for Solexa data
 *
 * looks at ADS for short reads, looks how far it can cut back the read
 * cutback rule: base quality <6, mismatches not farther than 3 apart
 *
 *************************************************************************/

//#define CEBUG(bla)   {if(CON_reads.size()>=9 && CON_reads.size()<=11 ){cout << bla; cout.flush();}}
//#define CEBUG(bla)   {if(CON_counts.size()>=5900 && CON_counts.size()<=6100 ){cout << bla; cout.flush();}}

//#define CEBUG(bla)   {cout << bla; cout.flush();}

int32 Contig::analyseADSForCuttingBackCERMap(const AlignedDualSeq & ads, int32 direction_frnid)
{
  FUNCSTART("void Contig::analyseADSForCuttingBackCERMap(const AlignedDualSeq & ads, int32 direction_frnid)");

  //BUGIFTHROW(ads.getOffsetInAlignment(-1)>0, "ads.getOffsetInAlignment(-1) > 0 ?");
  //BUGIFTHROW(ads.getRightOffsetInAlignment(-1)>0, "ads.getRightOffsetInAlignment(-1) > 0 ?");

  CEBUG(ads);

  // This is a read that claims it wants to be mapped but has too many errors
  int32 id=ads.getOtherID(-1);

  if(CON_readpool->getRead(id).getSequencingType() != ReadGroupLib::SEQTYPE_SOLEXA
     && CON_readpool->getRead(id).getSequencingType() != ReadGroupLib::SEQTYPE_ABISOLID){
    MIRANOTIFY(Notify::INTERNAL, "Trying to cut back other than Solexa / SOLiD ?");
  }

  // get pointer to the aligned contig
  const char * contigptr;
  const char * readptr;
  //int32 indexincontig;
  //int32 indexinread;

  // ... and advance it to the beginning of the read
  // we made sure above that (ads.getOffsetInAlignment(-1) == 0) is true

  /*
    in this routine, we will run always in read forward direction, i.e., the two
    forward sequences

    ID1:-1  GTTCCGGATTATGGCTTCACGCCCGCTACGCCCGATAATGTCGCCAAGA
    ID2:83   TTCCGGATTATGGCTTCACGCCCGCTCCGCACGTTA
    0-                                 X   X  X

    will run read 83 (and the contig) from left to right and if the read is
    reversed (as is 84)

    ID1:-1  GATCATATCCCGCCGGTCCCAGGCATGGTTGCCCAGCGTAATCACATCG
    ID2:84   AACAAACCCCGCCGGTCCCAGGCATGGTTGCCCAGC
    0-        X  X X

    will run from right to left.
   */

  int32 ptrincr=1;
  if(direction_frnid>0){
    contigptr= ads.getAlignedSequence(-1)+ads.getOffsetInAlignment(id);
    readptr= ads.getAlignedSequence(id);
    //indexinread=0;
    //indexincontig= offsetnewread;
  }else{
    contigptr= ads.getAlignedSequence(-1)+ads.getOffsetInAlignment(id)+ads.getOverlapLen()-1;
    readptr= ads.getAlignedSequence(id)+ads.getOverlapLen()-1;
    ptrincr=-1;
  }

  auto qv=CON_readpool->getRead(id).getQualities();
  auto qvI=qv.cbegin();
  std::advance(qvI,CON_readpool->getRead(id).getLeftClipoff());

  int32 canstartcutpos=-1;
  int32 lastcutpos=-1;
  int32 currentreadpos=CON_readpool->getRead(id).getLeftClipoff();

  for(uint32 ioverlap=0; ioverlap<ads.getOverlapLen(); ioverlap++, contigptr+=ptrincr, readptr+=ptrincr, qvI++, currentreadpos++){
    //CEBUG("io: " << ioverlap << "\ti_c: " << indexincontig);
    //CEBUG("\ti_r: " << indexinread);
    CEBUG("\nio: " << ioverlap);
    CEBUG("\t*rptr: " << *readptr);
    CEBUG("\t*qvI: " << static_cast<uint16>(*qvI));
    CEBUG("\t*cptr: " << *contigptr);
    CEBUG("\tcrp: " << currentreadpos);

    BUGIFTHROW(qvI==qv.end(),"quality out of bounds? 2");

    if(*readptr=='*') {
      currentreadpos--;
      // the following works only because we do not have (I hope)
      //  gaps at the very end of sequences in an ADS
      qvI--;
    }

    if(*contigptr=='#' && *readptr=='*') continue;
    if(*contigptr!=*readptr
       && *contigptr != 'N'
       && *readptr != 'N'){
      if(*readptr=='*'){
	// the following works only because we do not have (I hope)
	//  gaps at the very end of sequences in an ADS
	if((*qvI+*(qvI+1))/2<6) {
	  if(canstartcutpos<0) canstartcutpos=currentreadpos;
	  lastcutpos=currentreadpos;
	}
      }else{
	if(*qvI<6) {
	  if(canstartcutpos<0) canstartcutpos=currentreadpos;
	  lastcutpos=currentreadpos;
	}
      }
    }else{
      if((currentreadpos-lastcutpos)>3) {
	canstartcutpos=-1;
	lastcutpos=-1;
      }
    }
    CEBUG("\tcscp: " << canstartcutpos << "\tlcp: " << lastcutpos);
  }
  CEBUG("\ncutback result:\ncscp: " << canstartcutpos << "\tlcp: " << lastcutpos);

  canstartcutpos-=2;
  if(canstartcutpos<20) canstartcutpos=-1;

  return canstartcutpos;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * "maps" a read, return true if successful
 *
 * if read or contig contain IUPACs that differ, mapping cannot be done
 *  (return false)
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {if(CON_reads.size()>=9 && CON_reads.size()<=11 ){cout << bla; cout.flush();}}
//#define CEBUG(bla)   {if(CON_counts.size()>=5900 && CON_counts.size()<=6100 ){cout << bla; cout.flush();}}

//#define CEBUG(bla)   {cout << bla; cout.flush();}

bool Contig::insertMappedReadInContig(const AlignedDualSeq & ads, const uint8 newreadseqtype, const uint32 coffset, const int32 direction_frnid, const int32 direction_newid_incontig, int32 coveragemultiplier, bool forcemerge)
{
  //hier weiter: gaps anders behandeln (mitzaehlen) desgleichen auch in insertRead()

  FUNCSTART("bool Contig::insertMappedReadInContig(const AlignedDualSeq & ads, const uint8 newreadseqtype, const uint32 coffset, const int32 direction_frnid, int32 coveragemultiplier, bool forcemerge)");

  CEBUG("CON_reads:" << CON_reads.size() << endl);
  CEBUG("CON_counts:" << CON_counts.size() << endl);

  // This is a read that claims it can be mapped, so better be sure
  //  the consensus is completely covering it!
  BUGIFTHROW(ads.getOffsetInAlignment(-1)>0, "ads.getOffsetInAlignment(-1) > 0 ?");
  BUGIFTHROW(ads.getRightOffsetInAlignment(-1)>0, "ads.getRightOffsetInAlignment(-1) > 0 ?");

  // TODO: do we need that?
  definalise();

  // precalculate the offset into CON_counts.bbcounts[] according to
  //  sequencing type
  unsigned int sr_seqtypeoffset=0;
  if(newreadseqtype == ReadGroupLib::SEQTYPE_SOLEXA){
    sr_seqtypeoffset=0;
  }else{
    MIRANOTIFY(Notify::INTERNAL, "Trying to map other than Solexa?");
  }

  int32 id=ads.getOtherID(-1);

  // precalculate the bbstrains bitmask
  bbstrainmask_t strainmask=0;
  if(newreadseqtype==ReadGroupLib::SEQTYPE_SOLEXA){
    strainmask=getBBStrainMask(CON_readpool->getRead(id).getStrainID());
  }

  //CEBUG(*this);

  CEBUG("coffset: " << coffset << endl);

  // compute the offset of this read in the contig.
  uint32 offsetnewread=coffset+ads.getOffsetInAlignment(id);


  // mapping is done in two steps
  // first loop is simulation to determine if read can really be mapped,
  // second loop the mapping is really done
  // if we force the merge ... well, no first loop

  for(uint32 loopi=0; loopi<2; loopi++){
    if(forcemerge) loopi=1;

    CEBUG("loopi: " << loopi << endl);

    // get pointer to the aligned contig
    const char * contigptr;
    const char * readptr;
    int32 indexincontig;
    int32 indexinread;

    // ... and advance it to the beginning of the overlap
    // we made sure above that (ads.getOffsetInAlignment(-1) == 0) is true

    contigptr= ads.getAlignedSequence(-1)+ads.getOffsetInAlignment(id);
    readptr= ads.getAlignedSequence(id);
    indexincontig= offsetnewread;
    indexinread=0;

    auto ccI=CON_counts.begin();
    std::advance(ccI,offsetnewread);

    // get the qualities for this read and position a pointer on the first
    auto qv=CON_readpool->getRead(id).getQualities();
    auto qvI=qv.cbegin();
    if(direction_frnid>0){
      CEBUG("advancef " << CON_readpool->getRead(id).getLeftClipoff() << endl);
      std::advance(qvI,CON_readpool->getRead(id).getLeftClipoff());
    }else{
      CEBUG("advancer " << CON_readpool->getRead(id).getRightClipoff()-1 << endl);
      std::advance(qvI,CON_readpool->getRead(id).getRightClipoff()-1);
    }

    for(uint32 ioverlap=0; ioverlap<ads.getOverlapLen(); ++ioverlap, ++contigptr, ++readptr, ++indexincontig, ++indexinread, ++ccI, qvI+=direction_frnid){
      CEBUG("io: " << ioverlap << "\ti_c: " << indexincontig);
      CEBUG("\ti_r: " << indexinread);
      CEBUG("\t*rptr: " << *readptr);
      CEBUG("\t*qvI: " << (uint16) *qvI);
      CEBUG("\t*cptr: " << *contigptr);

//      /* WARNING!!!!!!!!!!!!!!!!!!!!!!!!!
//
//	 These two checks below which look at the boundaries of the quality
//	 iterator only work because at the moment MIRA enforces a left clip of 1 on
//	 all short reads to hide the tag in GAP4
//
//	 rethink that strategy should this change one day!
//      */
//      if(forcemerge){
//	// on forced merge, do everything so that we don't stop
//	if(qvI==qv.begin() || qvI==qv.end()) qvI-=direction_frnid;
//      }else{
//	BUGIFTHROW(qvI==qv.begin(),"quality out of bounds? 1");
//	BUGIFTHROW(qvI==qv.end(),"quality out of bounds? 2");
//      }

      if(indexincontig<0 || indexincontig >= CON_counts.size()+1) {
	// force it ... somehow ... to continue the for loop;
	if(forcemerge) continue;

	cout << "indexincontig: " << indexincontig << endl;
	cout << "\nCON_counts.size()+1: " << CON_counts.size()+1;

	cout << "\n\n" << ads << endl;

	setCoutType(AS_CAF);
	cout << *this;
	setCoutType(AS_TEXT);
	cout << *this;
	MIRANOTIFY(Notify::INTERNAL, "BOUNDCHECK error");
      }

      if(loopi){
	ccI->seqtype_cov[newreadseqtype]+=coveragemultiplier;
	ccI->total_cov+=coveragemultiplier;
      }

      BUGIFTHROW(!forcemerge && *contigptr=='*',"!forcemerge && *contigptr=='*' ?");

      // we may have a gap in the read, this is ok as long as the
      //   consensus has the code for "old gap" ('#')
      if(*readptr=='*'){
	if(!forcemerge && *contigptr!='#'){
	  // return only if we do not try to force the merge
	  // if we force ... do nothing

	  //cout << "\n\n" << ads;
	  //cout << "Alert: read has " << static_cast<char>(*readptr);
	  //cout << "\tcontig has " << static_cast<char>(*contigptr);
	  //cout << "\nio: " << ioverlap << "\ti_c: " << indexincontig;
	  //cout << "\ti_r: " << indexinread;
	  //cout << "\n";
	  //throw Notify(Notify::INTERNAL, THISFUNC, "logical error");

	  return false;
	}else{
	  if(loopi){
	    // set the strain bitmask
	    ccI->bbstrains[sr_seqtypeoffset]|=strainmask;

	    BUGIFTHROW(qvI-qv.begin()<0 || qvI-qv.begin() >= qv.size(),"chk 1 qvI " << qvI-qv.begin() << " out of bounds wrt " << qv.size());
	    // see whether we need to adapt the gap quality
	    uint16 gapqual=*qvI;
	    // move the quality pointer back by one to counterbalance
	    //  the for loop (illegally) munching up one quality value
	    //  as this was a *new* gap in the read
	    qvI-=direction_frnid;
	    BUGIFTHROW(qvI-qv.begin()<0 || qvI-qv.begin() >= qv.size(),"chk 2 qvI " << qvI-qv.begin() << " out of bounds wrt " << qv.size());
	    if(qvI!=qv.end()) gapqual=(gapqual+(*qvI))/2;

	    // increase the gap count
	    if (direction_newid_incontig > 0) {
	      ccI->bbcountsf[sr_seqtypeoffset]+=coveragemultiplier;
	      if(gapqual>ccI->bbbestqualsf[sr_seqtypeoffset]){
		ccI->bbbestqualsf[sr_seqtypeoffset]=static_cast<base_quality_t>(gapqual);
		BUGIFTHROW(ccI->bbbestqualsf[sr_seqtypeoffset]>100,"qualchk 1f >100");
	      }
	    } else {
	      ccI->bbcountsr[sr_seqtypeoffset]+=coveragemultiplier;
	      if(gapqual>ccI->bbbestqualsr[sr_seqtypeoffset]){
		ccI->bbbestqualsr[sr_seqtypeoffset]=static_cast<base_quality_t>(gapqual);
		BUGIFTHROW(ccI->bbbestqualsr[sr_seqtypeoffset]>100,"qualchk 1f >100");
	      }
	    }

	  }else{
	    // in the mean time, the qvI iterator must be decreased/increased as
	    //  no quality value exists at the moment in the read and we must counter
	    //  the increase/decrease in the for loop
	    qvI-=direction_frnid;

	  }
	}
      }else{
	if(forcemerge
	   || toupper(*readptr) == ccI->getOriginalBBChar()
	   || toupper(*readptr) == 'N'){            // in dubio pro reo
	  if(loopi){
	    ccI->bbstrains[sr_seqtypeoffset]|=strainmask;
	    BUGIFTHROW(qvI-qv.begin()<0 || qvI-qv.begin() >= qv.size(),"chk 3 qvI " << qvI-qv.begin() << " out of bounds wrt " << qv.size());
	    if (direction_newid_incontig > 0) {
	      ccI->bbcountsf[sr_seqtypeoffset]+=coveragemultiplier;
	      if(*qvI > ccI->bbbestqualsf[sr_seqtypeoffset]){
		ccI->bbbestqualsf[sr_seqtypeoffset]=*qvI;
		BUGIFTHROW(ccI->bbbestqualsf[sr_seqtypeoffset]>100,"qualchk 2 >100");
	      }
	    } else {
	      ccI->bbcountsr[sr_seqtypeoffset]+=coveragemultiplier;
	      if(*qvI > ccI->bbbestqualsr[sr_seqtypeoffset]){
		ccI->bbbestqualsr[sr_seqtypeoffset]=*qvI;
		BUGIFTHROW(ccI->bbbestqualsr[sr_seqtypeoffset]>100,"qualchk 2 >100");
	      }
	    }
	  }
	}else if(!forcemerge && *contigptr!='#'){
	  //cout << "\n\n" << ads;
	  //
	  //cout << "Alert2: read has " << static_cast<char>(*readptr);
	  //cout << "\tcontig has " << static_cast<char>(*contigptr);
	  //cout << "\nio: " << ioverlap << "\ti_c: " << indexincontig;
	  //cout << "\ti_r: " << indexinread;
	  //cout << "\nbackbone has: " << static_cast<char>(ccI->getBBChar()) << endl;
	  //throw Notify(Notify::INTERNAL, THISFUNC, "logical error 2");

	  // return only if we do not try to force the merge
	  // if we force ... do nothing
	  return false;
	}
      }
    }
  }

  CON_nummergedreads_perseqtype[sr_seqtypeoffset]+=coveragemultiplier;

  return true;
}

//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

//#define IRICTRACE

// coffset:  the offset of the contig part in the contig used in the ads
PlacedContigReads::const_iterator Contig::insertReadInContig(const AlignedDualSeq & ads, uint32 coffset, int32 direction_frnid, int32 direction_refid, int32 coveragemultiplier)
{
  FUNCSTART("void Contig::insertReadInContig(AlignedDualSeq & ads, uint32 coffset, int32 direction_frnid, int32 direction_refid, int32 coveragemultiplier)");

#ifdef CLOCK_STEPS_IRIC
  timeval us_start;
  timeval us_tmpstart;
  gettimeofday(&us_start,nullptr);
  uint32 numgapsincon=0;
  uint32 numreadinserts=0;
#endif

  ++CON_track_numins;

  CEBUG("CON_reads:" << CON_reads.size() << endl);
  CEBUG("CON_counts:" << CON_counts.size() << endl);

#ifdef IRICTRACE
  cout << "predef "; cout.flush();
#endif

  definalise();

#ifdef IRICTRACE
  cout << "postdef "; cout.flush();
#endif


  int32 newrid=ads.getOtherID(-1);

  //if(CON_readpool->getRead(id).getName()=="125691_2930_3355"){
  //  cout << "Dingodong\n";
  //  //cout << *this;
  //}

  //CEBUG(*this);

  CEBUG("coffset: " << coffset << endl);

  // compute the offset of this read in the contig.
  //  it might change later if the contig has grown to the left, but
  //  basically it's like this:
  int32 offsetnewread;
  offsetnewread=coffset+ads.getOffsetInAlignment(newrid);

  BUGIFTHROW(coffset > CON_counts.size(),"Something is wrong with coffset: " << coffset);
  BUGIFTHROW(ads.getOffsetInAlignment(newrid) > CON_counts.size(),"Something is wrong with that ADS!\n" << ads);
  BUGIFTHROW(offsetnewread>CON_counts.size(),"offsetnewread>CON_counts.size() ???");

  // Has the consensus grown to the left?
  int32 grownleft=0;
  if(ads.getOffsetInAlignment(newrid)==0 && ads.getOffsetInAlignment(-1)!=0){
    grownleft=-(coffset-ads.getOffsetInAlignment(-1));
    CEBUG("Grown left: " << grownleft << endl);

    // WRONG: BUGIFTHROW(grownleft<0, "grownleft <0 ???");
    // if grownleft <0 -> routine was not called with 'optimal' leftaligned ads
    // This can happen when the iteration in addRead() did not find the optimum
    //  and stoppped.
    // No big deal, this ain't an error.

    if(grownleft>0){
      // Yes, it has. So insert the needed positions into the
      //  contig and ...

      gettimeofday(&us_tmpstart,nullptr);
      CON_counts.insert(CON_counts.begin(), grownleft, CON_concounts_zero);
      CON_us_steps_iric[USCLOIRIC_INSGLCC]+=diffsuseconds(us_tmpstart);
      gettimeofday(&us_tmpstart,nullptr);

      // ... this new read is the (an) anchor in the contig, so offset=0 ...
      offsetnewread=0;

      //// ... and adjust the beginning of the reads which are behind.
      //for(uint32 k=0; k<CON_reads.size(); k++){
      //	CON_reads[k].offset+=grownleft;
      //	CEBUG("========" << CON_reads[k].offset);
      //}
      CON_reads.shiftReads(0,grownleft);
      shiftMarkerPositions(0,grownleft);

      CON_us_steps_iric[USCLOIRIC_INSGLARO]+=diffsuseconds(us_tmpstart);
      gettimeofday(&us_tmpstart,nullptr);

      // ... and adjust the consensus tags which are behind.
      for(uint32 k=0; k<CON_consensus_tags.size(); k++){
	CON_consensus_tags[k].from+=grownleft;
	CON_consensus_tags[k].to+=grownleft;
      }

      CON_us_steps_iric[USCLOIRIC_INSGLACT]+=diffsuseconds(us_tmpstart);
      gettimeofday(&us_tmpstart,nullptr);

    }else{
      grownleft=0;
    }
  }

#ifdef IRICTRACE
  cout << "postl1 "; cout.flush();
#endif

  CON_us_steps_iric[USCLOIRIC_INSGLTOT]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  // now place the read into the PlacedContigReads of the contig
  auto nprI=CON_reads.placeRead(CON_readpool->getRead(newrid),
			       newrid,
			       offsetnewread,
			       direction_frnid*direction_refid);

  CON_us_steps_iric[USCLOIRIC_PR]+=diffsuseconds(us_start);

#ifdef IRICTRACE
  cout << "postreadplace "; cout.flush();
#endif


  // add template if any
  if(CON_readpool->getRead(newrid).getTemplateID()>=0){
    gettimeofday(&us_start,nullptr);
    CON_templates_present.insert(CON_readpool->getRead(newrid).getTemplateID());
    CON_us_steps_iric[USCLOIRIC_TEMPL]+=diffsuseconds(us_start);
  }

#ifdef IRICTRACE
  cout << "postti "; cout.flush();
#endif

  gettimeofday(&us_start,nullptr);

  CEBUG("Offset of new read: " << offsetnewread << endl);
  CEBUG("Direction: " << direction_frnid*direction_refid << endl);

  // get pointer to the aligned contig
  const char * contigptr;
  const char * readptr;
  int32 indexincontig;
  int32 indexinread;

  // ... and advance it to the beginning of the overlap
  if(ads.getOffsetInAlignment(-1) == 0){
    contigptr= ads.getAlignedSequence(-1)+ads.getOffsetInAlignment(newrid);
    readptr= ads.getAlignedSequence(newrid);
    indexincontig= offsetnewread;
    indexinread=0;
  }else{
    contigptr= ads.getAlignedSequence(-1);
    readptr= ads.getAlignedSequence(newrid)+ads.getOffsetInAlignment(-1);
    //WRONG! causes Heisenbug
    // indexincontig=newread.offset+ads.getOffsetInAlignment(-1);

    //also wrong
    //indexincontig=newread.offset;

    if(grownleft>0) {
      indexincontig=offsetnewread+ads.getOffsetInAlignment(-1);
    }else{
      indexincontig=offsetnewread;
    }

    indexinread=ads.getOffsetInAlignment(-1);
  }

  CEBUG("Index in read: " << indexinread << endl);
  CEBUG("Index in contig: " << indexincontig << endl);
  BUGIFTHROW(indexinread<0, "indexinread < 0 ?");
  BUGIFTHROW(indexincontig<0, "indexincontig < 0 ?");

  CON_us_steps_iric[USCLOIRIC_INDEX]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  // Now go through the aligned contig and read base per base
  //  and see whether there must be inserted a gap in either one

  // TODO: in case this loop is still a major time eater, think of possible
  //  ways around
  // 1) either collect from the ADS all edits and then make the inserts in CON_reads
  //    looping through only once
  // 2) probably better: track "longest read" in contig, then start CON_reads loop
  //    not from begin(), but from read at nrpI.getReadStartOffset()-"len(longestread)"
  //    (eventually treat backbones differently???)
  // 3) combine 1+2 ???

  {
    // will be initialised the first time a gap must be inserted in
    //  reads. Once initialised, needs not to be searched again for
    auto quickfirstpcrI(CON_reads.end());

    base_quality_t gapqual;

    auto ccI=CON_counts.begin();
    for(int32 olli=0; olli<ads.getOverlapLen(); ++contigptr, ++readptr, ++indexincontig, ++indexinread, ++olli){
      CEBUG("olli: " << olli << "\ti_c: " << indexincontig);
      CEBUG("\ti_r: " << indexinread);
      CEBUG("\t*rptr: " << *readptr);
      CEBUG("\t*cptr: " << *contigptr);


      //BOUNDCHECK(indexincontig, 0, CON_counts.size()+1);

      if(indexincontig<0 || indexincontig >= CON_counts.size()+1) {
	cout << "indexincontig: " << indexincontig << endl;
	cout << "\nCON_counts.size()+1: " << CON_counts.size()+1;

	cout << "\n\n" << ads << endl;

	setCoutType(AS_CAF);
	cout << *this;
	setCoutType(AS_TEXT);
	cout << *this;
	MIRANOTIFY(Notify::INTERNAL, "BOUNDCHECK error");
      }
//      CEBUG(" CC: " << CON_counts[indexincontig].A << "\t"
//	    <<CON_counts[indexincontig].C << "\t"
//	    <<CON_counts[indexincontig].G << "\t"
//	    <<CON_counts[indexincontig].T << "\t"
//	    <<CON_counts[indexincontig].star << "\t"
//	    <<CON_counts[indexincontig].coverage << endl);


      // is there a gap in the contig?
      // TODO: wrong?
      if(*contigptr=='*'){
	// yes, is there also a gap in the read?
	BUGIFTHROW(*readptr=='*',"Both strands in ads have a '*' ?");
	CEBUG("\t* in con\n");

#ifdef CLOCK_STEPS_IRIC
	++numgapsincon;
	gettimeofday(&us_tmpstart,nullptr);
#endif

	// luckily (well, this should really be the normal case), only
	//  the contig has a *
	// so we must insert this into the consensus counts
	ccI=CON_counts.begin();

	CEBUG("Increase CON_counts: " << CON_counts.size());
	BUGIFTHROW(indexincontig>CON_counts.size(),"indexincontig (" << indexincontig << ") > CON_counts.size() (" << CON_counts.size() << ") ???");
	std::advance(ccI, indexincontig);
	ccI=CON_counts.insert(ccI, CON_concounts_zero);

	CON_us_steps_iric[USCLOIRIC_BIGLCCINS]+=diffsuseconds(us_tmpstart);
	gettimeofday(&us_tmpstart,nullptr);

	// CON_concounts_zero has '@' as standard backbone char, change to *
	ccI->i_backbonecharorig='*';
	ccI->i_backbonecharupdated='*';
	interpolateSRMValuesInCONcounts(ccI);
	CON_us_steps_iric[USCLOIRIC_BIGLINTERPOL]+=diffsuseconds(us_tmpstart);
	gettimeofday(&us_tmpstart,nullptr);

	CEBUG("\t" << CON_counts.size() << endl);

#ifdef IRICTRACE
	cout << "preuTBI "; cout.flush();
#endif
	// push up the consensus tags
	updateTagBaseInserted(indexincontig);
	CON_us_steps_iric[USCLOIRIC_BIGLUPDTAGS]+=diffsuseconds(us_tmpstart);
	gettimeofday(&us_tmpstart,nullptr);

#ifdef IRICTRACE
	cout << "postuTBI "; cout.flush();
#endif

	// now, we must insert a * in all the reads at that are at this
	//  position, and update 'star' and 'coverage' counters
	// for that, go through all the reads (except this new inserted one,
	//  that's what the -1 is for as the new read has been inserted at
	//  the back) ...

	auto lastpostocheck=nprI.getReadStartOffset()+nprI->getLenClippedSeq();
	if(quickfirstpcrI==CON_reads.end()){
	  quickfirstpcrI=getFirstPCRIForReadsCoveringPosition(indexincontig);
	  CON_us_steps_iric[USCLOIRIC_BIGLFPCRI]+=diffsuseconds(us_tmpstart);
	  gettimeofday(&us_tmpstart,nullptr);
	}
	// as PlacedContigReads::iterator is sorted by offset, we can stop the loop
	//  once we left the area covered by the newly inserted read
	for(auto pcrI=quickfirstpcrI; pcrI!=CON_reads.end() && pcrI.getReadStartOffset() <= lastpostocheck; ++pcrI) {
	  // we don't want to work on the newly inserted read;
	  if(pcrI.getORPID() == newrid) continue;
	  // ... and insert a * in those reads covering this position ...
	  if(indexincontig > pcrI.getReadStartOffset() &&
	     indexincontig < pcrI.getReadStartOffset()+pcrI->getLenClippedSeq()){ // checked: <!, not <=
#ifdef CLOCK_STEPS_IRIC
	    ++numreadinserts;
#endif
	    CEBUG("Read id: " << pcrI.getORPID());
	    CEBUG("\tfrom: " << pcrI.getReadStartOffset());
	    CEBUG("\tto: " << pcrI.getReadStartOffset()+pcrI->getLenClippedSeq());
	    CEBUG("\tin!");
	    CEBUG(ccI->A << "\t" <<ccI->C << "\t" <<ccI->G << "\t" <<ccI->T << "\t" <<ccI->star << "\t" <<ccI->total_cov << endl);
	    ccI->star+=coveragemultiplier;
	    ccI->total_cov+=coveragemultiplier;
	    ccI->seqtype_cov[pcrI->getSequencingType()]+=coveragemultiplier;
	    CEBUG(ccI->A << "\t" <<ccI->C << "\t" <<ccI->G << "\t" <<ccI->T << "\t" <<ccI->star << "\t" <<ccI->total_cov << endl);

	    int32 posinread;
	    posinread= indexincontig - pcrI.getReadStartOffset();
	    CEBUG("posinread: " << posinread);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
	    // TODO: what are gapquals for PacBio??? different depending HQ /LQ ...
	    if(pcrI.getReadDirection()>0){
	      if(pcrI->isSequencingType(ReadGroupLib::SEQTYPE_454GS20)
		 || pcrI->isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)) {
		gapqual=1;
	      }else{
		CEBUG("preqAQICS1 ");
		gapqual=pcrI->queryAverageQualInClippedSequence(posinread-1, posinread, false, true);
		CEBUG("postqAQICS2 ");
	      }
	      CEBUG("preqiBICS ");
	      const_cast<Read &>(*pcrI).insertBaseInClippedSequence('*', gapqual, posinread, true);
	      CEBUG("postiBICS ");
	    }else{
	      if(pcrI->isSequencingType(ReadGroupLib::SEQTYPE_454GS20)
		 || pcrI->isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)) {
		gapqual=1;
	      }else{
		CEBUG("preqAQICS2 ");
		gapqual=pcrI->queryAverageQualInClippedComplementSequence(posinread-1, posinread, false, true);
		CEBUG("postqAQICS2 ");
	      }
	      CEBUG("preqiBICCS ");
	      const_cast<Read &>(*pcrI).insertBaseInClippedComplementSequence('*', gapqual, posinread, true);
	      CEBUG("postqiBICCS ");
	    }
	    if(pcrI->getLenClippedSeq() > CON_longestreadseen) CON_longestreadseen=pcrI->getLenClippedSeq();
	    if(!pcrI->isBackbone() && pcrI->getLenClippedSeq() > CON_longestnonbbreadseen) CON_longestnonbbreadseen=pcrI->getLenClippedSeq();
	    // should never happen as rails should be added in a different manner ... but one never knows
	    // maybe in future?
	    if(pcrI->isRail() && pcrI->getLenClippedSeq() > CON_longestrailseen){
	      CON_longestrailseen=static_cast<int32>(pcrI->getLenClippedSeq());
	    }
	  }
	  CEBUG(endl);
	}
	CON_us_steps_iric[USCLOIRIC_BIGLLIGAPINREADS]+=diffsuseconds(us_tmpstart);
	gettimeofday(&us_tmpstart,nullptr);

	// ... and push the reads being at or behind this position one up
	CON_reads.shiftReads(indexincontig,1);
	// same thing for contig marker positions
	shiftMarkerPositions(indexincontig,1);
	CON_us_steps_iric[USCLOIRIC_BIGLSHIFTREADS]+=diffsuseconds(us_tmpstart);
	gettimeofday(&us_tmpstart,nullptr);

#ifdef IRICTRACE
	cout << "postfl1 "; cout.flush();
#endif

	CEBUG(ccI->A << "\t" <<ccI->C << "\t" <<ccI->G << "\t" <<ccI->T << "\t" <<ccI->star << "\t" <<ccI->total_cov << endl);
      }else{
	// No gap in the contig. Could there be one in the read?
	// If not, do nothing ...

	if(*readptr=='*'){
	  CEBUG("\t* in rea");
	  // if yes, we must insert a gap exactly in this read
	  int32 posinread;
	  posinread= indexincontig - nprI.getReadStartOffset();

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

	  // TODO: what are gapquals for PacBio??? different depending HQ / LQ
	  if(nprI.getReadDirection()>0){
	    if(nprI->isSequencingType(ReadGroupLib::SEQTYPE_454GS20)
	      || nprI->isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)) {
	      gapqual=1;
	    }else{
	      gapqual=nprI->queryAverageQualInClippedSequence(posinread-1, posinread, false, true);
	    }
	    const_cast<Read &>(*nprI).insertBaseInClippedSequence('*', gapqual, posinread, true);
	  }else{
	    if(nprI->isSequencingType(ReadGroupLib::SEQTYPE_454GS20)
	       || nprI->isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)) {
	      gapqual=1;
	    }else{
	      gapqual=nprI->queryAverageQualInClippedComplementSequence(posinread-1, posinread, false, true);
	    }
	    const_cast<Read &>(*nprI).insertBaseInClippedComplementSequence('*', gapqual, posinread, true);
	  }
	}
      }
      CEBUG("\n");
    }
  }

  CON_us_steps_iric[USCLOIRIC_BIGLTOT]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

#ifdef IRICTRACE
  cout << "prex1 "; cout.flush();
#endif
  // good, now see if we have to add something to the right of the contig
  //  (happens when the read grows over the right of the contig)
  CEBUG("nprI.getReadStartOffset(): " << nprI.getReadStartOffset());
  CEBUG("\tnprI->getLenClippedSeq(): " << nprI->getLenClippedSeq());
  CEBUG("\tCON_counts.size(): " << CON_counts.size() << endl);

  int32 offset=nprI.getReadStartOffset();
  int32 nrlen=static_cast<int32>(nprI->getLenClippedSeq());
  int32 grownright=offset+nrlen-CON_counts.size();
  if(grownright>0){
    CEBUG("Grown right: " << grownright << endl);
    CON_counts.resize(CON_counts.size()+grownright,CON_concounts_zero);
  }
  CON_us_steps_iric[USCLOIRIC_INSGR]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  if(nrlen > CON_longestreadseen){
    CON_longestreadseen=nrlen;
  }
  if(!nprI->isBackbone() && (nrlen > CON_longestnonbbreadseen)){
    CON_longestnonbbreadseen=nrlen;
  }
  // should never happen as rails should be added in a different manner ... but one never knows
  // maybe in future?
  if(nprI->isRail() && (nrlen > CON_longestrailseen)){
    CON_longestrailseen=nrlen;
  }

#ifdef IRICTRACE
  cout << "prex2 "; cout.flush();
#endif

  decltype(nprI->getClippedSeqIterator()) sI;
  if(nprI.getReadDirection() > 0){
    sI=nprI->getClippedSeqIterator();
  }else{
    sI=nprI->getClippedComplementSeqIterator();
  }

#ifdef IRICTRACE
  cout << "preuCV "; cout.flush();
#endif

  updateCountVectors(offset,
		     nrlen,
		     sI,
		     nprI->getSequencingType(),
		     true,
		     coveragemultiplier);

  CON_us_steps_iric[USCLOIRIC_UCV]+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);


#ifdef IRICTRACE
  cout << "postuCV "; cout.flush();
#endif

  return nprI;
}
//#define CEBUG(bla)
//#define CEBUG(bla)

//#undef CLOCK_STEPS



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Contig::interpolateSRMValuesInCONcounts(cccontainer_t::iterator ccI)
{
  FUNCSTART("void Contig::interpolateSRMValuesInCONcounts(cccontainer_t::iterator ccI)");
  // We need to adjust the counts of short read mappings:
  //  take the mean of the positions left and right (if possible)
  // Also, the bbstrains must be guessed

  for(unsigned int actsrtype=0; actsrtype<NUMMERGESEQTYPES; actsrtype++){
    uint16 poslooked=0;
    ccctype_t bbcountnewposf=0;
    base_quality_t qualf=0;
    ccctype_t bbcountnewposr=0;
    base_quality_t qualr=0;
    // theoretically, the "if" should not be needed
    if(ccI != CON_counts.begin()){
      bbcountnewposf=(ccI-1)->bbcountsf[actsrtype];
      qualf=(ccI-1)->bbbestqualsf[actsrtype];
      bbcountnewposr=(ccI-1)->bbcountsr[actsrtype];
      qualr=(ccI-1)->bbbestqualsr[actsrtype];
      poslooked++;
    }
    // theoretically, the "if" should not be needed
    if((ccI+1) != CON_counts.end()){
      bbcountnewposf+=(ccI+1)->bbcountsf[actsrtype];
      qualf+=(ccI+1)->bbbestqualsf[actsrtype];
      bbcountnewposr+=(ccI+1)->bbcountsr[actsrtype];
      qualr+=(ccI+1)->bbbestqualsr[actsrtype];
      poslooked++;
    }
    if(poslooked){
      ccI->bbcountsf[actsrtype]=bbcountnewposf/poslooked;
      ccI->bbbestqualsf[actsrtype]=qualf/poslooked;
      ccI->bbcountsr[actsrtype]=bbcountnewposr/poslooked;
      ccI->bbbestqualsr[actsrtype]=qualr/poslooked;
    }else{
      // should never happen (reads of size 1???), but anyway
      ccI->bbcountsf[actsrtype]=0;
      ccI->bbcountsr[actsrtype]=0;
      MIRANOTIFY(Notify::INTERNAL, "I'm on a branch I shouldn't be. Really!");
    }


    // theoretically, the "if" should not be needed
    if(ccI != CON_counts.begin()
       && (ccI+1) != CON_counts.end()){
      ccI->bbstrains[actsrtype]=((ccI-1)->bbstrains[actsrtype])
	& ((ccI+1)->bbstrains[actsrtype]);
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

void Contig::initialiseBaseLocks()
{
  FUNCSTART("void Contig::initialiseBaseLocks()");

  for(auto & cce : CON_counts){
    cce.baselock=0;
    cce.snplock=0;
  }

  for(auto pcrI=CON_reads.begin();pcrI != CON_reads.end(); ++pcrI){
    updateBaseLocks(pcrI,true);
  }

  FUNCEND();
}


/*************************************************************************
 *
 * When adding: transfers base locks from read tag into the CON_counts[].baselock
 * When removing: decreases CON_counts[].baselock where baselock tags are present in read
 *
 *************************************************************************/

//#define CEBUG(bla)   {if(dodebug) {cout << bla; cout.flush();}}
void Contig::updateBaseLocks(PlacedContigReads::const_iterator pcrI, bool addiftrue)
{

  FUNCSTART("void Contig::updateBaseLocks(PlacedContigReads::const_iterator pcrI, bool addiftrue)");

  //bool dodebug=pcrI->getName()=="...";

  // one = 1 seems ... oh well
  // but might want to change the icrement of baselocks in the future
  uint16 one=1;
  if(addiftrue){
    CEBUG("Adding base locks in " << pcrI->getName() << endl);
  }else{
    CEBUG("Removing base locks in " << pcrI->getName() << endl);
    one=-one;
  }

  for(int32 tagi=0; tagi<pcrI->getNumOfTags(); ++tagi){
    const multitag_t & acttag=pcrI->getTag(tagi);
    bool baselock=false;
    bool snplock=false;
    for(int32 j=0; j<CON_baselock_ids.size(); ++j){
      if(acttag.identifier==CON_baselock_ids[j]){
	baselock=true;
	break;
      }
    }
    for(int32 j=0; j<CON_snplock_ids.size(); ++j){
      if(acttag.identifier==CON_snplock_ids[j]){
	snplock=true;
	break;
      }
    }

    if(baselock || snplock){
      CEBUG("Tag #" << tagi << "\tfrom: " << acttag.from << "\tto: " << acttag.to << "\t" << acttag.getIdentifierStr());
      int32 contigpos=pcrI.unclippedReadPos2ContigPos(acttag.from);
      CEBUG("\tconpos: " << contigpos << endl);
      // no, bad, it _might_ really be outside the contig
      //  e.g. if tags are in clipped part
      //BOUNDCHECK(contigpos, 0, CON_counts.size());
      if(contigpos>=0 && contigpos<CON_counts.size()){
	auto ccI=CON_counts.begin();
	std::advance(ccI, contigpos);
	auto readstart=pcrI.getReadStartOffset();
	auto readend=readstart+pcrI->getLenClippedSeq();
	CEBUG("readstart: " << readstart << "\treadend: " << readend << endl);
	for(int32 j=0; j<=acttag.to-acttag.from; ++contigpos, ++ccI, ++j){
	  CEBUG("cpos: " << contigpos << "\t" << *ccI << endl);
	  // We mind only for the parts of the lock in the used read
	  if(contigpos>=readstart && contigpos < readend){
	    CEBUG("Incr!\n");
	    if(baselock) ccI->baselock+=one;
	    if(snplock) ccI->snplock+=one;
	  }else{
	    CEBUG("nop\n");
	  }
	}
      }
    }
  }

  CEBUG("done updbl\n");

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * update the count vectors
 * the CON_counts vector mus have enough space to contain the updated seq
 * from    is a position in the contig
 * len     is the length of the updated portion
 * updateI is the sequence which updates into the contig
 * seqtype is the seqtype of the read that is added
 *
 * coveragemultiplier>0 : adding to contig
 * coveragemultiplier<0 : removing to contig
 *
 * Does NOT update the baselock! This is done in updateBaseLocks
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Contig::updateCountVectors(const int32 from, const int32 len, std::vector<char>::const_iterator updateI, const uint32 seqtype, const bool addiftrue, int32 coveragemultiplier)
{

  FUNCSTART("void Contig::updateCountVectors(const int32 from, const int32 len, std::vector<char>::const_iterator updateI, const uint32 seqtype, const bool addiftrue, int32 coveragemultiplier)");

  CEBUG("From: "<<from<<"\tLen: "<<len<<"\tCON_counts.size():" << CON_counts.size() << endl);

  BUGIFTHROW(from<0, "from " << from << " < 0 ?");

  if(from+len>CON_counts.size()){
    cout << "Error:\n";
    cout << "from: " << from << endl;
    cout << "len: " << len << endl;
    cout << "size of contig: " << CON_counts.size() << endl;
    MIRANOTIFY(Notify::INTERNAL, "from + len > size of contig?");
  }

  auto ccI=CON_counts.begin();
  BOUNDCHECK(from, 0, CON_counts.size());
  std::advance(ccI, from);
  CEBUG("ccI " << ccI << endl);

  int32 one=1*coveragemultiplier;
  int32 two=2*coveragemultiplier;
  int32 four=4*coveragemultiplier;

  for(int32 i=0; i<len; ++ccI, ++updateI, ++i){
    CEBUG("before i " << i << "\tcp " << from+i << "\tccI: " << *ccI << endl);
    char thechar=*updateI;
    CEBUG("char: " << thechar << endl);
    switch(toupper(thechar)){
    case '-':
    case 'N': {
      // BaCh 13.08.2018
      // I think I now do not need to account for Ns like this anymore
      //ccI->A+=one;
      //ccI->C+=one;
      //ccI->G+=one;
      //ccI->T+=one;
      ccI->N+=one;
      break;
    }
    case 'X': {ccI->X+=one; break;}
    case 'A': {ccI->A+=four; break;}
    case 'C': {ccI->C+=four; break;}
    case 'G': {ccI->G+=four; break;}
    case 'T': {ccI->T+=four; break;}

    case 'M': {ccI->A+=two; ccI->C+=two; break;}
    case 'R': {ccI->A+=two; ccI->G+=two; break;}
    case 'W': {ccI->A+=two; ccI->T+=two; break;}
    case 'S': {ccI->C+=two; ccI->G+=two; break;}
    case 'Y': {ccI->C+=two; ccI->T+=two; break;}
    case 'K': {ccI->G+=two; ccI->T+=two; break;}

    case 'V': {ccI->A+=one; ccI->C+=one; ccI->G+=one; break;}
    case 'H': {ccI->A+=one; ccI->C+=one; ccI->T+=one; break;}
    case 'D': {ccI->A+=one; ccI->G+=one; ccI->T+=one; break;}
    case 'B': {ccI->C+=one; ccI->G+=one; ccI->T+=one; break;}

    case '*': {ccI->star+=one; break;}
    default: {
      cout << "WHY? Illegal char: " << (uint16) thechar << " >>" << thechar << "<<\n";
      MIRANOTIFY(Notify::FATAL, "Unexpected base.");
    }
    }
    ccI->total_cov+=one;
    ccI->seqtype_cov[seqtype]+=one;
    CEBUG("after i " << i << "\tcp " << from+i << "\tccI: " << *ccI << endl);
  }

  FUNCEND();
}
//#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::finalise()
{
  FUNCSTART("void Contig::finalise()");

  //checkContig();

  if(CON_finalised==false){
    makeTmpConsensus(0, CON_counts.size(), CON_tmpcons_from_backbone, true);

    CON_finalised=true;
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

PlacedContigReads::const_iterator Contig::deleteRead(PlacedContigReads::const_iterator pcrI)
{
  //TODO: What happens, when a contig breaks? (no more coverage in the middle)

  FUNCSTART("void Contig::deleteRead(uint32 id)");

  CEBUG("before deleteRead():\n"; setCoutType(AS_TEXT); cout << *this << endl);

#ifdef CLOCK_STEPS_DRFC
  timeval us_tv;
  gettimeofday(&us_tv,nullptr);
  timeval us_tvtotal=us_tv;
#endif

  ++CON_track_numdels;

  if(CON_fixedconsseq.size()){
    nukeSTLContainer(CON_fixedconsseq);
    nukeSTLContainer(CON_fixedconsqual);
  }

  CEBUG("Deleting: " << id << "\nContig before:\n");
  CEBUG(*this);

  definalise();

  paranoiaBUGSTAT(checkContig());

  auto coveragemultiplier=pcrI->getDigiNormMultiplier();

  // update of internal contig statistics
  CON_readsperstrain[pcrI->getStrainID()]-=coveragemultiplier;
  CON_readsperreadgroup[pcrI->getReadGroupID().getLibId()]-=coveragemultiplier;

  // Remove the base locks in CON_counts that this read produced
  gettimeofday(&us_tv,nullptr);
  updateBaseLocks(pcrI, false);
  CON_us_steps_drfc[USCLODRFC_UBL]+=diffsuseconds(us_tv);

  gettimeofday(&us_tv,nullptr);

  // Remove the read from the CON_counts
  auto offset=pcrI.getReadStartOffset();
  int32 rlen=pcrI->getLenClippedSeq();
  decltype(pcrI->getClippedSeqIterator()) sI;
  if(pcrI.getReadDirection() > 0){
    sI=pcrI->getClippedSeqIterator();
  }else{
    sI=pcrI->getClippedComplementSeqIterator();
  }

  CEBUG("Offset: " << offset);
  CEBUG("\nrlen: " << rlen);

  updateCountVectors(offset,
		     rlen,
		     sI,
		     pcrI->getSequencingType(),
		     false,
		     -coveragemultiplier);
  CON_us_steps_drfc[USCLODRFC_UCV]+=diffsuseconds(us_tv);

  CEBUG(*this);


  int32 frontdeletions=0;
  int32 enddeletions=0;
  {
    // Readjust the ends of the contig if the removed read was there

    // Readjust front
    //   readjust the offsets if needed
    int32 maxchecklen=rlen;
    BUGIFTHROW(maxchecklen<0, "front: maxchecklen < 0 ?");
    if(maxchecklen>0 && offset==0){
      frontdeletions=chompFront(maxchecklen);
    }

    BUGIFTHROW(frontdeletions<0, "frontdeletions < 0 ?");

    if(offset>0){
      offset=offset-frontdeletions;
      BUGIFTHROW(offset<0, "offset >0 and offset-frontdeletions < 0 ?");
    }

    BUGIFTHROW(offset<0, "offset < 0 ?");

    // Readjust end
    {
      maxchecklen=rlen-frontdeletions;
      BUGIFTHROW(maxchecklen<0, "end: maxchecklen< 0 ?");
      if(maxchecklen>0){
	enddeletions=chompBack(maxchecklen);
      }
    }
  }

  BUGIFTHROW(enddeletions<0, "enddeletions < 0 ?");

  // delete template if any
  {
    int32 tid=pcrI->getTemplateID();
    if(tid>=0){
      gettimeofday(&us_tv,nullptr);
      auto msI=CON_templates_present.find(tid);
      if(msI==CON_templates_present.end()){
	MIRANOTIFY(Notify::INTERNAL, "Template not present in list though read has one?");
      }
      CON_templates_present.erase(msI);
      CON_us_steps_drfc[USCLODRFC_DT]+=diffsuseconds(us_tv);
    }
  }

  // Remove read
  gettimeofday(&us_tv,nullptr);
  auto retI=CON_reads.removeRead(pcrI);
  CON_us_steps_drfc[USCLODRFC_RR]+=diffsuseconds(us_tv);

  CEBUG("after removeRead:\n" << *this << endl);

#if 1
  CEBUG("deleteStarOnlyColumns(" << offset << ", " << offset+rlen-frontdeletions-enddeletions << ")\n");
  gettimeofday(&us_tv,nullptr);
  deleteStarOnlyColumns(offset, offset+rlen-frontdeletions-enddeletions);
  CON_us_steps_drfc[USCLODRFC_DSOC]+=diffsuseconds(us_tv);
#endif

  CON_us_steps_drfc[USCLODRFC_TOTAL]+=diffsuseconds(us_tvtotal);

  paranoiaBUGSTAT(checkContig());

  CEBUG("Contig after:\n" << *this);

  FUNCEND();

  return retI;
}






/*************************************************************************
 *
 * "doshiftreads==false" is a hack for trimMapOverhang() which needs
 *  chompFront() functionality WITHOUT shifting reads
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
int32 Contig::chompFront(int32 maxchecklen, bool doshiftreads)
{
  FUNCSTART("int32 Contig::chompFront(int32 maxchecklen, bool doshiftreads)");

#ifdef CLOCK_STEPS_DRFC
  timeval us_tv;
#endif

  gettimeofday(&us_tv,nullptr);

  if(maxchecklen<0) maxchecklen=static_cast<int32>(CON_counts.size());
  int32 frontdeletions=0;
  auto ccI= CON_counts.begin();
  while(ccI!=CON_counts.end()
	&& ccI->total_cov==0
	&& frontdeletions<maxchecklen){
    ++frontdeletions;
    ++ccI;
  }
  if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_ITF]+=diffsuseconds(us_tv);

  if(frontdeletions>0){
    CEBUG("CONTIG-- delread at beginning!\n");
    // ok, we've got to delete bases at the front of the consensus
    gettimeofday(&us_tv,nullptr);
    CON_counts.erase(CON_counts.begin(), ccI);
    if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_CCEF]+=diffsuseconds(us_tv);
    CEBUG("Deleted << " << frontdeletions << " at front.\n");

    // push consensus tags down, deleting all those which are completely
    //  in the deleted part, adjusting others.
    gettimeofday(&us_tv,nullptr);
    for(auto ctI=CON_consensus_tags.begin(); ctI!=CON_consensus_tags.end();){
      ctI->from-=frontdeletions;
      ctI->to-=frontdeletions;
      // static casts are mean trick as the multitag variables are uint32!
      if(static_cast<int32>(ctI->to) < 0){
	ctI=CON_consensus_tags.erase(ctI);
      }else{
	if(static_cast<int32>(ctI->from) < 0) ctI->from=0;
	++ctI;
      }
    }
    if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_SDT]+=diffsuseconds(us_tv);

    // push down the offsets of the other reads
    if(doshiftreads){
      gettimeofday(&us_tv,nullptr);
      CON_reads.shiftReads(1,-frontdeletions);
      if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_SR]+=diffsuseconds(us_tv);
    }
    // same thing for contig marker positions
    gettimeofday(&us_tv,nullptr);
    shiftMarkerPositions(1,-frontdeletions);
    if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_SMP]+=diffsuseconds(us_tv);
  }

  if(frontdeletions>0) definalise();

  return frontdeletions;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
int32 Contig::chompBack(int32 maxchecklen)
{
  FUNCSTART("int32 Contig::chompBack(int32 maxchecklen)");

  if(CON_counts.begin()==CON_counts.end()) return 0;

#ifdef CLOCK_STEPS_DRFC
  timeval us_tv;
#endif

  gettimeofday(&us_tv,nullptr);

  if(maxchecklen<0) maxchecklen=static_cast<int32>(CON_counts.size());
  int32 enddeletions=0;
  auto ccI= CON_counts.end();
  --ccI;
  if(ccI->total_cov==0){
    CEBUG("CONTIG-- delread at end!\n");
    gettimeofday(&us_tv,nullptr);
    while(ccI->total_cov==0 &&
	  enddeletions<maxchecklen){
      CEBUG(".");
      ++enddeletions;
      if(ccI==CON_counts.begin()) break;
      --ccI;
    }
    if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_ITB]+=diffsuseconds(us_tv);
    if(enddeletions>0){
      ++ccI;
      CEBUG("Counts.size() before: " << CON_counts.size() << endl);
      gettimeofday(&us_tv,nullptr);
      CON_counts.erase(ccI,CON_counts.end());
      if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_CCEB]+=diffsuseconds(us_tv);
      CEBUG("Counts.size() after: " << CON_counts.size() << endl);
      CEBUG("Deleted << " << enddeletions << " at end.\n");

      // and delete consensus tags there
      gettimeofday(&us_tv,nullptr);
      int32 newsize=CON_counts.size()-1;
      auto ctI=CON_consensus_tags.begin();
      while(ctI!=CON_consensus_tags.end()){
	if((ctI->from > newsize)){
	  ctI=CON_consensus_tags.erase(ctI);
	}else{
	  if(ctI->to > newsize) ctI->to=newsize;
	  ++ctI;
	}
      }
      if(CON_us_steps_drfc.size()) CON_us_steps_drfc[USCLODRFC_SDT]+=diffsuseconds(us_tv);
    }
  }

  if(enddeletions>0) definalise();

  return enddeletions;
}
#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Contig::trimMapOverhang()
{
  FUNCSTART("void Contig::trimMapOverhang()");

  //cout << "MMMMMMMMMM\n";
  //CON_reads.debugDump(false);

  if(getNumBackbones()==0) return;

  int32 firstknownbbpos=-1;
  int32 lastknownbbpos=-1;
  for(auto crI= CON_reads.begin(); crI!=CON_reads.end(); ++crI){
    if(crI->isBackbone()){
      if(static_cast<int32>(crI.getReadStartOffset()+crI->getLenClippedSeq()) > lastknownbbpos){
	lastknownbbpos=crI.getReadStartOffset()+crI->getLenClippedSeq();
      }
      if(firstknownbbpos<0){
	firstknownbbpos=crI.getReadStartOffset();
      }
    }
  }

  CEBUG("fkbb: " << firstknownbbpos << endl);
  CEBUG("lkbb: " << lastknownbbpos << endl);
  CEBUG("gcl: " << getContigLength() << endl);

  BUGIFTHROW(firstknownbbpos<0,"firstknownbbpos " << firstknownbbpos << " <0 ???");
  BUGIFTHROW(lastknownbbpos<0,"lastknownbbpos " << lastknownbbpos << " <0 ???");

  // check: are there reads completely outside the backbone? if yes: no-can-do
  for(auto crI= CON_reads.begin(); crI!=CON_reads.end(); ++crI){
    if(!crI->isBackbone()){
      if(crI.getReadStartOffset() < firstknownbbpos
	 && crI.getReadStartOffset()+crI->getLenClippedSeq() <= firstknownbbpos){
	// oooops, impossible at front
	return;
      }
      if(crI.getReadStartOffset() >= lastknownbbpos){
	// oooops, impossible at end
	return;
      }
    }
  }

  bool docorrect=false;
  // we're fine, now get things rolling
  //
  // first, the back of the contig if needed
  if(lastknownbbpos < getContigLength()){
    docorrect=true;
    // zero out CON_counts
    auto ccI=CON_counts.end(); // cannot use rbegin, not implemented in HDeque
    CEBUG("numdelsteps: " << getContigLength()-lastknownbbpos << endl);
    for(auto numdelsteps=getContigLength()-lastknownbbpos; numdelsteps!=0; --numdelsteps){
      --ccI;
      *ccI=CON_concounts_zero;
    }

    // shorten overhanging reads
    for(auto crI= CON_reads.begin(); crI!=CON_reads.end(); ++crI){
      int32 shorten=static_cast<int32>(crI.getReadStartOffset()+crI->getLenClippedSeq())-lastknownbbpos;
      CEBUG(crI->getName() << " shortenb " << shorten << endl);
      if(shorten>0){
	// very hacky, but we're forcing our hand on the reads of the PCR container here
	Read & ncr=const_cast<Read &>(*crI);
	if(crI.getReadDirection()>0){
	  ncr.setRQClipoff(crI->getRightClipoff()-shorten);
	}else{
	  ncr.setLQClipoff(crI->getLeftClipoff()+shorten);
	}
      }
    }
  }

  // now the front of the contig if needed
  if(firstknownbbpos>0){
    docorrect=true;
    // zero out CON_counts
    auto ccI=CON_counts.begin();
    for(auto numdelsteps=firstknownbbpos; numdelsteps!=0; --numdelsteps, ++ccI){
      *ccI=CON_concounts_zero;
    }

    // shorten overhanging reads
    for(auto crI= CON_reads.begin(); crI!=CON_reads.end(); ++crI){
      int32 shorten=firstknownbbpos-static_cast<int32>(crI.getReadStartOffset());
      CEBUG(crI->getName() << "\tfkbbp: " << firstknownbbpos << "\trso: " << crI.getReadStartOffset() <<"\tshortenf " << shorten << endl);
      if(shorten>0){
	// very hacky, but we're forcing our hand on the reads of the PCR container here
	Read & ncr=const_cast<Read &>(*crI);
	if(crI.getReadDirection()>0){
	  ncr.setLQClipoff(crI->getLeftClipoff()+shorten);
	}else{
	  ncr.setRQClipoff(crI->getRightClipoff()-shorten);
	}
      }
    }
  }

  if(docorrect){
    chompBack(-1);

    if(firstknownbbpos>0){
      // hacky thing
      chompFront(-1,false);
      // hacky thing
      CON_reads.shiftReadsBounceZero(0,-firstknownbbpos);
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

void Contig::checkContig()
{
#ifdef PARANOIABUGTRACK
  FUNCSTART("void Contig::checkContig()");

  // TODO: ich muss, da ansonsten nicht als CAF abspeicherbar
  //         (erlaube ich 0 längenreads oder nicht?)

  auto crI= CON_reads.begin();
  for(uint32 i=0; i<CON_reads.size(); i++, crI++){
    if(crI->offset<0){
      MIRANOTIFY(Notify::INTERNAL, "offset < 0?");
    }
    if(!CON_counts.empty()){
      if(crI->offset>=static_cast<int32>(CON_counts.size())){
	cout << "Error. Contig_id: " << CON_id << endl;
	cout << "First read: " << CON_reads[0].read;
	MIRANOTIFY(Notify::INTERNAL, "offset >= contigsize?");
      }
      if(crI->offset+crI->read.getLenClippedSeq()>CON_counts.size()){
	MIRANOTIFY(Notify::INTERNAL, "offset+readlen > contigsize?");
      }
    }
  }

  if(CON_miraparams->size() != Contig::SEQTYPE_END){
    MIRANOTIFY(Notify::INTERNAL, "CON_miraparams->size() != Contig::SEQTYPE_END");
  }

#if 1
  if(!CON_counts.empty() &&
     CON_counts.front().coverage == 0
     ){
    MIRANOTIFY(Notify::INTERNAL, "0-0-0-0-0 N (zero coverage) at front?");
  }
  if(!CON_counts.empty() &&
     CON_counts.back().coverage == 0
     ){
    MIRANOTIFY(Notify::INTERNAL, "0-0-0-0-0 N (zero coverage) at end?");
  }
#endif

#ifdef CEBUGFLAG
  for(auto ccI=CON_counts.begin(); ccI != CON_counts.end(); ++ccI){
    if(ccI->A + ccI->C + ccI->G + ccI->T + (ccI->star)*4 != (ccI->coverage)*4){
      MIRANOTIFY(Notify::INTERNAL, "Gna, error in Con_counts :/.");
    }
  }
#endif

  FUNCEND();
#endif
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO Fehler, wenn clips sich gegenüber dem original verschieben

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Contig::initialiseContig(const std::list<contig_init_read_t> & rlist,
			      std::vector<multitag_t> & tags,
			      const std::string & contigname,
			      std::string & fixedseq,
			      std::vector<base_quality_t> & fixedqual)
{
  FUNCSTART("void Contig::initialiseContig(const std::list<contig_init_read_t> & rlist, std::vector<multitag_t> & tags, const std::string & contigname, std::string & fixedseq, std::vector<base_quality_t> & fixedqual)");

  discard();

  CON_name=contigname;
  CON_fixedconsseq=fixedseq;
  CON_fixedconsqual=fixedqual;

  CEBUG("::initC1 fixed " << CON_fixedconsseq.size() << endl);

  CON_counts.clear();

  //cout << "KKKKKKKKK " << CON_readpool->getNumOfStrainInReadpool() << endl;
  // ... and strain tracking
  BUGIFTHROW(ReadGroupLib::getNumOfStrains()==0,"CON_readpool->getNumOfStrainInReadpool()==0 ???");
  CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);

  BUGIFTHROW(ReadGroupLib::getNumReadGroups()==0,"CON_readpool->getNumReadGroups()==0 ???");
  CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);

  // Correct for weird things:
  // 1) joining contigs in gap5 inserts gap bases having a quality of 0, leading to the possibility
  //    that a single wrong bases overrules hundred of valid gaps.
  //    Counterstrategy: in reads, re-assign quality to gap bases with quality 0
  // 2) no "two" atm

  CEBUG("Fixing zero-gap qualities\n");
  for(auto & rle : rlist){
    CON_readpool->getRead(rle.id).fixZeroGapQuals();
  }


  //
  CEBUG("Going through list (size " << rlist.size() << "):\n");
  uint32 computedclen=0;
  for(auto & rle : rlist){
    CEBUG("New element:\t");
    CEBUG("Offset: " << rle.offset_in_contig << "\tID: " << rle.id);
    CEBUG("\tDir: " << rle.direction << endl);
    CEBUG("Read: " << CON_readpool->getRead(rle.id).getName() << endl);

    auto pcrI=CON_reads.placeRead(CON_readpool->getRead(rle.id),
				  rle.id,
				  rle.offset_in_contig,
				  rle.direction);
    computedclen=std::max(computedclen,static_cast<uint32>(pcrI.getReadStartOffset())+static_cast<uint32>(pcrI->getLenClippedSeq()));

    auto coveragemultiplier=CON_readpool->getRead(rle.id).getDigiNormMultiplier();
    CON_readsperstrain[pcrI->getStrainID()]+=coveragemultiplier;
    CON_readsperreadgroup[pcrI->getReadGroupID().getLibId()]+=coveragemultiplier;
    if(pcrI->getLenClippedSeq() > CON_longestreadseen) CON_longestreadseen=pcrI->getLenClippedSeq();
    if(pcrI->isBackbone()) CON_isbackbonecontig=true;
    if(!pcrI->isBackbone() && pcrI->getLenClippedSeq() > CON_longestnonbbreadseen) CON_longestnonbbreadseen=pcrI->getLenClippedSeq();
    if(pcrI->isRail() && pcrI->getLenClippedSeq() > CON_longestrailseen) CON_longestrailseen=pcrI->getLenClippedSeq();

    // add template if any
    if(pcrI->getTemplateID()>=0){
      CON_templates_present.insert(pcrI->getTemplateID());
    }
  }

  CEBUG("end of list" << endl);

  // CAF as written by gap5 may start at another offset than 0
  // MIRA expects the first read to have an offset of 0, therefore correct for that
  int32 shiftoffset=-CON_reads.begin().getReadStartOffset();
  if(shiftoffset != 0){
    computedclen+=shiftoffset;
    CON_reads.shiftReads(0,shiftoffset);
  }

  CEBUG("Computed contig length: " << computedclen << endl);

  BUGIFTHROW(computedclen<=0,"Contig has size " << computedclen << " which is <= 0 ???");
  CON_counts.resize(computedclen, CON_concounts_zero);

  for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
    auto coveragemultiplier=pcrI->getDigiNormMultiplier();
    try{
      if(pcrI.getReadDirection() > 0){
	updateCountVectors(pcrI.getReadStartOffset(),
			   pcrI->getLenClippedSeq(),
			   pcrI->getClippedSeqIterator(),
			   pcrI->getSequencingType(),
			   true,
			   coveragemultiplier);
      }else{
	updateCountVectors(pcrI.getReadStartOffset(),
			   pcrI->getLenClippedSeq(),
			   pcrI->getClippedComplementSeqIterator(),
			   pcrI->getSequencingType(),
			   true,
			   coveragemultiplier);
      }
    }
    catch(Notify n){
      cout << "Readpool id: " << pcrI.getORPID() << endl;
      cout << "Error while adding this read (1):\n" << *pcrI;
      cout << "In this contig (output might crash):" << endl << *this << endl;
      n.handleError(THISFUNC);
    }


    // initialise backbone and short read mapping structures if this read was used as a
    //  backbone
    if(pcrI->isBackbone()){
      auto ccI=CON_counts.begin();
      auto rso=pcrI.getReadStartOffset();
      std::advance(ccI,rso);
      const char * consb=pcrI->getSeqAsChar();
      const base_quality_t * qptr=&(pcrI->getQualities()[rso]);
	for(uint32 ci=0; ci < pcrI->getLenClippedSeq(); ++ccI, ++consb, ++ci, ++qptr){
	BUGIFTHROW(ccI == CON_counts.end(), "ccI == CON_counts.end() ???");
	BUGIFTHROW(*consb == 0, "*consb == 0 ?");
	// Don't init bbcounts[], bbbestquals[] and bbstrains[] as they
	//  might already have been initialised by CER reads ... well,
	//  will be when this functionality is implemented.
	// So, only initialise backbone chars
	ccI->i_backbonecharorig=toupper(*consb);
	ccI->i_backbonecharupdated='@';
	ccI->i_backbonequalorig=*qptr;
      }
    }
  }

  CEBUG("CON_counts before trimming: " << CON_counts.size() << endl);

  // trim end of contig if needed
  uint64 numpops=0;
  while(!CON_counts.empty() &&
	CON_counts.back().total_cov == 0
	){
    CON_counts.pop_back();
    ++numpops;
  }
  if(numpops){
    CON_fixedconsseq.resize(CON_fixedconsseq.size()-numpops);
    CON_fixedconsqual.resize(CON_fixedconsseq.size());
  }
  CEBUG("Numpops front: " << numpops << endl);

  // trim start of contig if needed
  numpops=0;
  while(!CON_counts.empty() &&
	CON_counts.front().total_cov == 0
	){
    CON_counts.pop_front();
    ++numpops;
  }
  if(numpops){
    CON_fixedconsseq.erase(CON_fixedconsseq.begin(),CON_fixedconsseq.begin()+numpops);
    CON_fixedconsqual.erase(CON_fixedconsqual.begin(),CON_fixedconsqual.begin()+numpops);
  }
  CEBUG("Numpops back: " << numpops << endl);

  BUGIFTHROW(CON_counts.empty(),"After trimming, contig has zero length?");

  // must copy each tag as we have to "upgrade" a tag_t to a consensustag_t
  {
    CON_consensus_tags.clear();
    CON_consensus_tags.resize(tags.size());
    for(uint32 i=0;i!=tags.size();++i){
      CON_consensus_tags[i]=tags[i];
      CON_consensus_tags[i].from+=shiftoffset;
      CON_consensus_tags[i].to+=shiftoffset;
    }
  }

  //CEBUG(*this);

  try{
    checkContig();
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  definalise();

  CEBUG("CON_fixedconsseq.size(): " << CON_fixedconsseq.size() << " CON_counts.size(): " << CON_counts.size() << endl);
  if(!CON_fixedconsseq.empty()
     && CON_fixedconsseq.size() != CON_counts.size()) {
    // broken initialisation??? (seen in MAF of contig with error in align)
    // force recalc
    CON_fixedconsseq.clear();
    CON_fixedconsqual.clear();
    cout << "\nMinor problem with contig " << getContigName() << ", fixing.\n";
  }

  if(CON_fixedconsseq.size() != CON_fixedconsqual.size()
     || CON_fixedconsseq.empty()
     || CON_fixedconsqual.empty()){
    std::string dummy1;
    std::vector<base_quality_t> dummy2;
    newConsensusGet(dummy1, dummy2);
  }

  CEBUG("::initC2 fixed " << CON_fixedconsqual.size() << endl);
  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::vector<int32> &
Contig::getCRIndicesAtContigPosition(std::vector<int32> & vec, int32 pos1, int32 pos2) const
{
  FUNCSTART("std::vector<int32> & Contig::getCRIndicesAtContigPosition(std::vector<int32> & vec, int32 pos1, int32 pos2) const");

  BUGIFTHROW(true,"Redo for PlacedContigReads");

////++++////  BUGIFTHROW(pos1 < 0, "pos1 < 0?");
////++++////  BUGIFTHROW(pos1 >= static_cast<int32>(CON_counts.size()), "pos1 > size of contig?");
////++++////
////++++////  if(pos2==-1) pos2=pos1;
////++++////
////++++////  BUGIFTHROW(pos2 < pos1, "pos2 < pos1");
////++++////
////++++////  if(pos2 > static_cast<int32>(CON_counts.size())) pos2=static_cast<int32>(CON_counts.size()-1);
////++++////
////++++////  vec.clear();
////++++////  vec.reserve(100);
////++++////
////++++////  int32 i=0;
////++++////  std::vector<contigread_t>::const_iterator I= CON_reads.begin();
////++++////
////++++////  for(I= CON_reads.begin(); I != CON_reads.end(); I++, i++){
////++++////    if(I->offset > pos2
////++++////       || static_cast<int32>(I->offset+I->read.getLenClippedSeq()) <= pos1) continue;
////++++////    vec.push_back(i);
////++++////  }

  FUNCEND();

  return vec;
}


/*************************************************************************
 *
 * TODO: make intelligent version which does not start from begin()
 *   (needs contig to know length of longest read)
 *
 *************************************************************************/

void Contig::getPCRIteratorsAtContigPosition(std::vector<PlacedContigReads::const_iterator> & vec, int32 pos1, int32 pos2) const
{
  FUNCSTART("void Contig::getPCRIteratorsAtContigPosition(std::vector<PlacedContigReads::const_iterator> & vec, int32 pos1, int32 pos2) const");

  BUGIFTHROW(pos1 < 0, "pos1 < 0?");
  BUGIFTHROW(pos1 >= static_cast<int32>(CON_counts.size()), "pos1 > size of contig?");

  if(pos2==-1) pos2=pos1;

  BUGIFTHROW(pos2 < pos1, "pos2 < pos1");

  if(pos2 > static_cast<int32>(CON_counts.size())) pos2=static_cast<int32>(CON_counts.size()-1);

  vec.clear();

  auto pcrI(getFirstPCRIForReadsCoveringPosition(pos1));
  for(; pcrI != CON_reads.end() && pcrI.getReadStartOffset() <= pos2; ++pcrI){
    int32 endpoint=static_cast<int32>(pcrI.getReadStartOffset()+pcrI->getLenClippedSeq())-1;
    if(endpoint >= pos1) vec.push_back(pcrI);
  }

  FUNCEND();

  return;
}




/*************************************************************************
 *
 * Returns the original readpool-IDs of the reads between contig positions
 *  pos1 and pos2
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Contig::getReadORPIDsAtContigPosition(std::vector<int32> & vec, int32 pos1, int32 pos2) const
{
  FUNCSTART("void Contig::getReadORPIDsAtContigPosition(std::vector<int32> & vec, int32 pos1, int32 pos2) const");

  BUGIFTHROW(pos1 < 0, "pos1 (" << pos1 << ") < 0?");
  BUGIFTHROW(pos1 >= static_cast<int32>(CON_counts.size()), "pos1 (" << pos1 << ") > size of contig (" << CON_counts.size() << ") ?");

  if(pos2==-1) pos2=pos1;

  BUGIFTHROW(pos2 < pos1, "pos2 (" << pos2 << ") < pos1 (" << pos1 << ")");

  if(pos2 > static_cast<int32>(CON_counts.size())) pos2=static_cast<int32>(CON_counts.size()-1);

  vec.clear();

  CEBUG("gROACP: from " << pos1 << '\t' << pos2 <<endl);

  auto pcrI(getFirstPCRIForReadsCoveringPosition(pos1));
  for(; pcrI != CON_reads.end() && pcrI.getReadStartOffset() <= pos2; ++pcrI){
    int32 endpoint=static_cast<int32>(pcrI.getReadStartOffset()+pcrI->getLenClippedSeq())-1;
    CEBUG(pcrI->getName() << '\t' << pcrI.getReadStartOffset() << '\t' << pcrI->getLenClippedSeq() << '\t' << endpoint << '\t');
    if(endpoint >= pos1){
      CEBUG("taken\n");
      vec.push_back(pcrI.getORPID());
    }else{
      CEBUG("dropped\n");
    }
  }

  return;
}
//#define CEBUG(bla)




/*************************************************************************
 *
 * Starting from a reference id that is a rail, takes reads in readpool
 *  before and after that are rails and that fall into the same xcut,
 *  ycut region as given and pushed them into the vector
 * xcut, ycut =   [...[
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Contig::getRailsAsReadsAffected(const int32 refid, std::vector<int32> & reads_affected, const int32 xcut, const int32 ycut)
{
  FUNCSTART("void Contig::getRailsAsReadsAffected(const int32 refid, std::vector<int32> & reads_affected)");

  BUGIFTHROW(!CON_readpool->getRead(refid).isRail(),"!CON_readpool->getRead(refid).isRail() ???");
  BUGIFTHROW(CON_longestrailseen==0,"CON_longestrailseen==0???");

  CEBUG("grara xcut " << xcut << "\tycut " << ycut << endl);

  bool refisin=false;
  reads_affected.clear();
  auto pcrI=getFirstPCRIForRailsCoveringPosition(xcut);
  CEBUG("grara 1st: " << pcrI << endl);
  for(; pcrI != CON_reads.end() && pcrI.getReadStartOffset() < ycut; ++pcrI){
    CEBUG("grara check " << pcrI << endl);
    if(pcrI.getReadStartOffset()+pcrI->getLenClippedSeq() < xcut) continue;
    if(pcrI->isRail()){
      CEBUG("pushed back" << endl);
      reads_affected.push_back(pcrI.getORPID());
      if(pcrI.getORPID()==refid) refisin=true;
    }
  }
  CEBUG("grara rasize: " << reads_affected.size() << endl);

  if(!refisin){
    cout << "BUG BUG BUG BUG BUG\n";
    cout << "Searched for range " << xcut << " " << ycut << endl;
    cout << "CON_longestrailseen: " << CON_longestrailseen << endl;
    cout << "First PCRI: " << getFirstPCRIForRailsCoveringPosition(xcut) << endl;
    for(auto tpcrI= CON_reads.begin(); tpcrI != CON_reads.end(); ++tpcrI){
      cout << tpcrI << endl;
    }
    MIRANOTIFY(Notify::FATAL, "*sigh*");
  }

  FUNCEND();
  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::deleteBaseInRead(int32 contigposition, int32 readindex)
{
  FUNCSTART("void Contig::deleteBaseInRead(int32 contigposition, int32 readindex)");

  BUGIFTHROW(true,"Redo for PlacedContigReads");

////++++////  paranoiaBUGSTAT(checkContig());
////++++////
////++++////  BUGIFTHROW((readindex<0 || readindex>=CON_reads.size()),
////++++////	     "Readindex <0 || >= CON_reads.size() -> Read not in this contig!");
////++++////
////++++////  contigread_t & theread=CON_reads[readindex];
////++++////
////++++////  BUGIFTHROW((contigposition < theread.offset ||
////++++////	      contigposition >= theread.offset+theread.read.getLenClippedSeq()),
////++++////	     "Read does not cover given contigposition");
////++++////
////++++////  definalise();
////++++////
////++++////  int32 posinread=contigposition-theread.offset;
////++++////  int32 updatelen=theread.read.getLenClippedSeq()-posinread;
////++++////
////++++////  CEBUG("contigpos: " << contigposition);
////++++////  CEBUG("\nPPPosinread: " << posinread);
////++++////  CEBUG("\nupdatelen: " << updatelen<<endl);
////++++////
////++++////  bool deletein=false;
////++++////  bool deleteend=false;
////++++////
////++++////  // now only clip contig when size is affected, never ever at all move
////++++////  // other reads around ---> this may lead to a contig break :(
////++++////  if(theread.offset+theread.read.getLenClippedSeq()==CON_counts.size()
////++++////     && CON_counts.back().total_cov==1) deleteend=true;
////++++////
////++++////  CEBUG("deletein: " << deletein << endl);
////++++////  CEBUG("deleteend: " << deleteend << endl);
////++++////
////++++////  std::vector<char>::const_iterator I;
////++++////
////++++////  if(theread.direction>0){
////++++////    I=theread.read.getClippedSeqIterator();
////++++////    BOUNDCHECK(posinread, 0, theread.read.getLenClippedSeq());
////++++////    std::advance(I, posinread);
////++++////    updateCountVectors(contigposition,
////++++////		       updatelen,
////++++////		       I,
////++++////		       theread.read.getSequencingType(),
////++++////		       false);
////++++////    theread.read.deleteBaseFromClippedSequence(posinread);
////++++////    I=theread.read.getClippedSeqIterator();
////++++////  }else{
////++++////    I=theread.read.getClippedComplementSeqIterator();
////++++////    BOUNDCHECK(posinread, 0, theread.read.getLenClippedSeq());
////++++////    std::advance(I, posinread);
////++++////    updateCountVectors(contigposition,
////++++////		       updatelen,
////++++////		       I,
////++++////		       theread.read.getSequencingType(),
////++++////		       false);
////++++////    theread.read.deleteBaseFromClippedComplementSequence(posinread);
////++++////    I=theread.read.getClippedComplementSeqIterator();
////++++////  }
////++++////
////++++////  if(deletein==true){
////++++////    auto cI=CON_counts.begin();
////++++////    BOUNDCHECK(contigposition, 0, CON_counts.size());
////++++////    std::advance(cI, contigposition);
////++++////    CON_counts.erase(cI);
////++++////    for(uint32 i=0; i<CON_reads.size(); i++){
////++++////      if(CON_reads[i].offset>contigposition) CON_reads[i].offset--;
////++++////    }
////++++////    // adjust tags in consensus
////++++////    // TODO: untested (ok?)
////++++////    updateTagBaseDeleted(contigposition);
////++++////  }
////++++////  if(deleteend==true){
////++++////    CON_counts.pop_back();
////++++////    // adjust tags in consensus
////++++////    // TODO: untested (ok?)
////++++////    updateTagBaseDeleted(contigposition);
////++++////  }
////++++////
////++++////  // BOUNDS checked above
////++++////  std::advance(I, posinread);
////++++////  updateCountVectors(contigposition,
////++++////		     updatelen-1,
////++++////		     I,
////++++////		     theread.read.getSequencingType(),
////++++////		     true);
////++++////
////++++////  paranoiaBUGSTAT(checkContig());
////++++////
////++++//////  if(contigposition==4711){
////++++//////    cout << *this;
////++++//////  }

  FUNCEND();

}



/*************************************************************************
 *
 * Adds/Subtracts the offset to the offset of the reads in the contig
 *  which occur _after_ contigposition
 *
 *************************************************************************/

void Contig::adjustReadOffsets(int32 contigposition, int32 offset)
{
  FUNCSTART("void Contig::adjustReadOffsets(int32 contigposition, int32 offset)");

  paranoiaBUGSTAT(checkContig());

  BUGIFTHROW(contigposition<0,"contigposition<0 ?");

  if(offset!=0){
    CON_reads.shiftReads(contigposition,offset);
    shiftMarkerPositions(contigposition,offset);
    priv_rebuildConCounts();
  }

  paranoiaBUGSTAT(checkContig());

  FUNCEND();

}


void Contig::shiftMarkerPositions(int32 contigposition, int32 offset)
{
  for(auto & mp : CON_markerpositions){
    if(mp>=contigposition) mp+=offset;
  }
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_rebuildConCounts()
{
  FUNCSTART("void Contig::rebuildConsensus()");

  // do NOT clear the entire element! ... just what's needed: the counters and locks
  for(auto ccI=CON_counts.begin(); ccI != CON_counts.end(); ++ccI){
    ccI->A=0;             // ACGT are extended counters
    ccI->C=0;             // ccctype_t is enough for a coverage of 16k reads.
    ccI->G=0;
    ccI->T=0;
    ccI->N=0;             // N, X, star and coverage are normal counters.
    ccI->X=0;
    ccI->star=0;
    ccI->total_cov=0;
    ccI->baselock=0;     // if > 0 then only one of A, C, T, G may be set, otherwise it's a misassembly
    ccI->snplock=0;
    for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; ++i){
      ccI->seqtype_cov[i]=0;
    };
  }

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); pcrI++){
    auto coveragemultiplier=pcrI->getDigiNormMultiplier();
    if(pcrI.getReadDirection()>0){
      updateCountVectors(pcrI.getReadStartOffset(),
			 pcrI->getLenClippedSeq(),
			 pcrI->getClippedSeqIterator(),
			 pcrI->getSequencingType(),
			 true,
			 coveragemultiplier);
    }else{
      updateCountVectors(pcrI.getReadStartOffset(),
			 pcrI->getLenClippedSeq(),
			 pcrI->getClippedComplementSeqIterator(),
			 pcrI->getSequencingType(),
			 true,
			 coveragemultiplier);
    }
  }

  for(auto ccI=CON_counts.begin(); ccI != CON_counts.end(); ++ccI){
    auto addcount = ccI->bbcountsf[0]+ccI->bbcountsr[0];
    ccI->seqtype_cov[ReadGroupLib::SEQTYPE_SOLEXA] += addcount;
    ccI->total_cov += addcount;
  }

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO: base and quality check? Offsets after?
void Contig::insertBaseInRead(char base,
			      int32 contigposition,
			      int32 readindex,
			      base_quality_t quality)
{
//  if(contigposition==4709){
//    cout << *this;
//  }


  FUNCSTART("void Contig::insertBaseInRead(char base, base_quality_t quality, int32 contigposition, int32 readindex)");

  BUGIFTHROW(true,"Redo for PlacedContigReads");

////++++////  paranoiaBUGSTAT(checkContig());
////++++////
////++++////  BUGIFTHROW((readindex<0 || readindex>=CON_reads.size()),
////++++////	     "Readindex <0 || >= CON_reads.size() -> Read not in this contig!");
////++++////
////++++////  contigread_t & theread=CON_reads[readindex];
////++++////
////++++////  BUGIFTHROW((contigposition < theread.offset ||
////++++////	      contigposition >= theread.offset+theread.read.getLenClippedSeq()),
////++++////	     "Read does not cover given contigposition");
////++++////
////++++////  definalise();
////++++////
////++++////  int32 posinread=contigposition-theread.offset;
////++++////  int32 updatelen=theread.read.getLenClippedSeq()-posinread;
////++++////
////++++////
////++++////  bool insertin=false;
////++++////  bool insertend=false;
////++++////
////++++////  // now only extend contig when size is affected, never ever at all move
////++++////  // other reads around
////++++////  if(theread.offset+theread.read.getLenClippedSeq()==CON_counts.size()
////++++////     //    FALSCH!!!  && CON_counts.back().coverage==1)
////++++////     ) insertend=true;
////++++////
////++++////
////++++////  CEBUG("insertin: " << insertin << endl);
////++++////  CEBUG("insertend: " << insertend << endl);
////++++////
////++++////  base_quality_t rqual=0;
////++++////
////++++////  std::vector<char>::const_iterator I;
////++++////
////++++////  if(theread.direction>0){
////++++////    if(quality==BQ_UNDEFINED || quality>100){
////++++////      quality=theread.read.queryAverageQualInClippedSequence(posinread-1, posinread, false, true);
////++++////    }
////++++////    rqual=quality;
////++++////
////++++////    I=theread.read.getClippedSeqIterator();
////++++////    BOUNDCHECK(posinread, 0, theread.read.getLenClippedSeq()+1);
////++++////    std::advance(I, posinread);
////++++////    updateCountVectors(contigposition,
////++++////		       updatelen,
////++++////		       I,
////++++////		       theread.read.getSequencingType(),
////++++////		       false);
////++++////    theread.read.insertBaseInClippedSequence(base, rqual, posinread,true);
////++++////    I=theread.read.getClippedSeqIterator();
////++++////  }else{
////++++////    if(quality==BQ_UNDEFINED || quality>100){
////++++////      quality=theread.read.queryAverageQualInClippedComplementSequence(posinread-1, posinread, false, true);
////++++////    }
////++++////    rqual=quality;
////++++////
////++++////    I=theread.read.getClippedComplementSeqIterator();
////++++////    BOUNDCHECK(posinread, 0, theread.read.getLenClippedSeq()+1);
////++++////    std::advance(I, posinread);
////++++////    updateCountVectors(contigposition,
////++++////		       updatelen,
////++++////		       I,
////++++////		       theread.read.getSequencingType(),
////++++////		       false);
////++++////    theread.read.insertBaseInClippedComplementSequence(base, rqual, posinread,true);
////++++////    I=theread.read.getClippedComplementSeqIterator();
////++++////  }
////++++////
////++++////  if(insertin==true){
////++++////    auto ccI=CON_counts.begin();
////++++////    BOUNDCHECK(contigposition, 0, CON_counts.size()+1);
////++++////    std::advance(ccI, contigposition);
////++++////    ccI=CON_counts.insert(ccI,CON_concounts_zero);
////++++////    ccI->getBBChar()='*';
////++++////    interpolateSRMValuesInCONcounts(ccI);
////++++////
////++++////    for(uint32 i=0; i<CON_reads.size(); i++){
////++++////      if(CON_reads[i].offset>contigposition) CON_reads[i].offset++;
////++++////    }
////++++////    // adjust tags in consensus
////++++////    // TODO: untested (ok?)
////++++////    updateTagBaseInserted(contigposition);
////++++////  }
////++++////  if(insertend==true){
////++++////    CON_counts.push_back(CON_concounts_zero);
////++++////  }
////++++////
////++++////  std::advance(I, posinread);
////++++////  updateCountVectors(contigposition,
////++++////		     updatelen+1,
////++++////		     I,
////++++////		     theread.read.getSequencingType(),
////++++////		     true);
////++++////
////++++////  paranoiaBUGSTAT(checkContig());
////++++////
////++++//////  if(contigposition==4709){
////++++//////    cout << *this;
////++++//////  }
////++++////
////++++////

  FUNCEND();

  return;
}








/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// untested, should work

char Contig::getBaseInRead(const int32 contigposition, const PlacedContigReads::const_iterator & pcrI) const
{
  FUNCSTART("char Contig::getBaseInRead(const int32 contigposition, const PlacedContigReads::const_iterator & pcrI) const");

  BUGIFTHROW((contigposition < pcrI.getReadStartOffset() ||
	      contigposition >= pcrI.getReadStartOffset() + pcrI->getLenClippedSeq()),
	     "Read does not cover given contigposition");

  int32 readpos=pcrI.contigPos2UnclippedReadPos(contigposition);
  FUNCEND();
  if(pcrI.getReadDirection() > 0){
    return toupper(pcrI->getBaseInSequence(readpos));
  }else{
    return toupper(pcrI->getBaseInComplementSequence(readpos));
  }
  return '$';
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// untested, should work

base_quality_t Contig::getQualityInRead(int32 contigposition, const PlacedContigReads::const_iterator & pcrI) const
{
  FUNCSTART("char Contig::getQualityInRead(int32 contigposition, int32 readindex)");

  BUGIFTHROW((contigposition < pcrI.getReadStartOffset() ||
	      contigposition >= pcrI.getReadStartOffset() + pcrI->getLenClippedSeq()),
	     "Read does not cover given contigposition");

  int32 readpos=pcrI.contigPos2UnclippedReadPos(contigposition);
  FUNCEND();
  if(pcrI.getReadDirection() > 0){
    return pcrI->getQualityInSequence(readpos);
  }

  return pcrI->getQualityInComplementSequence(readpos);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO: check base and qualities?
void Contig::changeBaseInRead(char base,
			      int32 contigposition,
			      int32 readindex,
			      base_quality_t quality)
{
  FUNCSTART("void Contig::changeBaseInRead(char base, base_quality_t quality, int32 contigposition, int32 readindex)");

  BUGIFTHROW(true,"Redo for PlacedContigReads");

////++++////  BUGIFTHROW((readindex<0 || readindex>=CON_reads.size()),
////++++////	     "Readindex <0 || >= CON_reads.size() -> Read not in this contig!");
////++++////
////++++////  //  Read & theread=CON_reads[readindex].read;
////++++////  contigread_t & theread=CON_reads[readindex];
////++++////
////++++////  BUGIFTHROW((contigposition < theread.offset ||
////++++////	      contigposition >= theread.offset+theread.read.getLenClippedSeq()),
////++++////	     "Read does not cover given contigposition");
////++++////
////++++////
////++++////  definalise();
////++++////
////++++////  base_quality_t rqual=0;
////++++////
////++++////  int32 posinread=contigposition-CON_reads[readindex].offset;
////++++////  std::vector<char>::const_iterator I;
////++++////
////++++////  if(theread.direction>0){
////++++////    if(quality==BQ_UNDEFINED || quality>100){
////++++////      quality=theread.read.queryAverageQualInClippedSequence(posinread-1, posinread+1, true, true);
////++++////    }
////++++////    rqual=quality;
////++++////
////++++////    I=theread.read.getClippedSeqIterator();
////++++////    BOUNDCHECK(posinread, 0, theread.read.getLenClippedSeq());
////++++////    std::advance(I, posinread);
////++++////    updateCountVectors(contigposition,
////++++////		       1 ,
////++++////		       I,
////++++////		       theread.read.getSequencingType(),
////++++////		       false);
////++++////    theread.read.changeBaseInClippedSequence(base, rqual, posinread);
////++++////  }else{
////++++////    if(quality==BQ_UNDEFINED || quality>100){
////++++////      quality=theread.read.queryAverageQualInClippedComplementSequence(posinread-1, posinread+1, true, true);
////++++////    }
////++++////    rqual=quality;
////++++////
////++++////    I=theread.read.getClippedComplementSeqIterator();
////++++////    BOUNDCHECK(posinread, 0, theread.read.getLenClippedSeq());
////++++////    std::advance(I, posinread);
////++++////    updateCountVectors(contigposition,
////++++////		       1,
////++++////		       I,
////++++////		       theread.read.getSequencingType(),
////++++////		       false);
////++++////    theread.read.changeBaseInClippedComplementSequence(base, rqual, posinread);
////++++////  }
////++++////
////++++////  updateCountVectors(contigposition,
////++++////		     1,
////++++////		     I,
////++++////		     theread.read.getSequencingType(),
////++++////		     true);
////++++////
  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::addTagToRead(uint32 contigpositionfrom,
			  uint32 contigpositionto,
			  int32 readindex,
			  multitag_t & mt)
{
  FUNCSTART("void Contig::addTagToRead(int32 contigpositionfrom, int32 contigpositionto, int32 readindex, const multitag & mt)");

  BUGIFTHROW(true,"Redo for PlacedContigReads");

////++++////  BUGIFTHROW((readindex<0 || readindex>=CON_reads.size()),
////++++////	     "Readindex <0 || >= CON_reads.size() -> Read not in this contig!");
////++++////
////++++////  //  Read & theread=CON_reads[readindex].read;
////++++////  contigread_t & theread=CON_reads[readindex];
////++++////
////++++////  BUGIFTHROW((contigpositionfrom < theread.offset ||
////++++////	      contigpositionfrom > theread.offset+theread.read.getLenClippedSeq()),
////++++////	     "Read does not cover given contigposition (from)");
////++++////  BUGIFTHROW((contigpositionto < static_cast<uint32>(theread.offset) ||
////++++////	      contigpositionto > theread.offset+theread.read.getLenClippedSeq()),
////++++////	     "Read does not cover given contigposition (to)");
////++++////
////++++////
////++++////  definalise();
////++++////
////++++////  if(contigpositionto < contigpositionfrom) std::swap(contigpositionfrom,
////++++////						 contigpositionto);
////++++////
////++++////  //TODO: use contigPos2UnclippedReadPos or similar???
////++++////
////++++////  uint32 posinreadfrom;
////++++////  uint32 posinreadto;
////++++////  uint32 len=(contigpositionto-contigpositionfrom);
////++++////  if(CON_reads[readindex].direction>0){
////++++////    posinreadfrom=contigpositionfrom-theread.offset+theread.read.getLeftClipoff();
////++++////    posinreadto=posinreadfrom+len;
////++++////  }else{
////++++////    posinreadto=theread.read.getLenClippedSeq()-(contigpositionfrom-theread.offset)+theread.read.getLeftClipoff()-1;
////++++////    posinreadfrom=posinreadto-len;
////++++////  }
////++++////
////++++////  mt.from=posinreadfrom;
////++++////  mt.to=posinreadto;
////++++////
////++++////  theread.read.addTagO(mt);

  FUNCEND();
}




// TODO: write version with multitag_t::mte_id_t??

/*************************************************************************
 *
 * Adds a tag to the consensus if this tag is not already present
 *
 *
 *************************************************************************/

Contig::consensustag_t & Contig::addTagToConsensus(const uint32 contigposfrom,
						   const uint32 contigposto,
						   const char strand,
						   const char * identifier,
						   const char * comment,
						   const bool doublecheck
  )
{
  FUNCSTART("void Contig::addTagToConsensus(int32 contigposfrom, int32 contigposto, const char * identifier, const char * comment)");

  BUGIFTHROW((contigposto < contigposfrom),
	     "contigposto < contigposfrom not allowed");
  //BUGIFTHROW((contigposfrom < 0),
  //	     "contigposfrom < 0 not allowed");
  BUGIFTHROW((contigposfrom >= CON_counts.size()),
	     "contigposfrom > size of contig");
  //BUGIFTHROW((contigposto < 0),
  //	     "contigposto < 0 not allowed");
  BUGIFTHROW((contigposto >= CON_counts.size()),
	     "contigposto > size of contig");

  BUGIFTHROW(!(strand=='+' || strand=='-' || strand=='='),
	     "strand is not +, - or =");

  // TODO: commented out 22.04.2007, check whether ok
  // not needed? and disturbs setting tags in makeIntelligentConsensus()
  // definalise();

  consensustag_t tmptag;
  tmptag.identifier=multitag_t::newIdentifier(identifier);
  tmptag.comment=multitag_t::newComment(comment);
  tmptag.from=contigposfrom;
  tmptag.to=contigposto;
  tmptag.setStrand(strand);
  tmptag.additionalinfo_initialised=false;

  auto ctI = CON_consensus_tags.begin();

  if(doublecheck){
    for(; ctI!=CON_consensus_tags.end(); ++ctI){
      if(ctI->from==tmptag.from
	 && ctI->to==tmptag.to
	 && ctI->getStrand()==tmptag.getStrand()
	 && ctI->identifier==tmptag.identifier
	) {
	break;
      }
    }
  }else{
    ctI = CON_consensus_tags.end();
  }

  if(ctI==CON_consensus_tags.end()) {
    CON_consensus_tags.push_back(tmptag);
    return CON_consensus_tags.back();
  }

  *ctI=tmptag;

  FUNCEND();
  return *ctI;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::updateTagBaseDeleted(uint32 contigpos)
{
  FUNCSTART("void Contig::updateTagBaseDeleted(uint32 contigpos)");
  //BUGIFTHROW((contigpos < 0),
  //	     "contigpos < 0 not allowed");
  BUGIFTHROW((contigpos >= CON_counts.size()),
	     "contigpos > size of contig");

  for(auto ctI=CON_consensus_tags.begin(); ctI!=CON_consensus_tags.end();){
    if(ctI->to == ctI->from && ctI->to == contigpos){
      ctI=CON_consensus_tags.erase(ctI);
    }else{
      if(contigpos <= ctI->to) ctI->to-=1;
      if(contigpos < ctI->from) ctI->from-=1;
      ++ctI;
    }
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::updateTagBaseInserted(uint32 contigpos)
{
  FUNCSTART("void Contig::updateTagInsert(uint32 contigpos)");

  //BUGIFTHROW((contigpos < 0),
  //	     "contigpos < 0 not allowed");
  BUGIFTHROW((contigpos >= CON_counts.size()),
	     "contigpos > size of contig");

  for(auto & te : CON_consensus_tags){
    if(te.from >= contigpos) te.from+=1;
    if(te.to >= contigpos) te.to+=1;
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


int32 Contig::getRealReadPos(const uint32 contigpos, const PlacedContigReads::const_iterator & pcrI) const
{
  if(pcrI.getReadDirection() > 0){
    return pcrI.contigPos2UnclippedReadPos(contigpos);
  }
  return pcrI->calcComplPos(pcrI.contigPos2UnclippedReadPos(contigpos));
}





/*************************************************************************
 *
 * Transfers SRMr tags in reads to a ORMB tag in the contig
 *
 *
 *************************************************************************/
void Contig::transposeReadSRMTagsToContig()
{
  FUNCSTART("void Contig::transposeReadSRMTagsToContig()");

  BUGIFTHROW(true,"Redo for PlacedContigReads");

////++++////  //cout << "We are in contig: " << getContigName() << endl;
////++++////
////++++////  auto rI=CON_reads.cbegin();
////++++////  for(;rI!=CON_reads.cend();rI++){
////++++////
////++++////    const contigread_t & ric =*rI;
////++++////
////++++////    //cout << "We are in read: " << ric.read.getName() << endl;
////++++////
////++++////    for(uint32 i=0; i<ric.read.getNumOfTags(); i++){
////++++////      multitag_t acttag=ric.read.getTag(i);
////++++////      //cout << "Tag #" << i << "\tfrom: " << acttag.from << "\tto: " << acttag.to << "\t" << acttag.identifier[0] << acttag.identifier[1] << acttag.identifier[2] << acttag.identifier[3] << endl;
////++++////      if(acttag.identifier==Read::REA_tagentry_idSRMr
////++++////	 || acttag.identifier==Read::REA_tagentry_idWRMr ) {
////++++////
////++++////	int32 tleftbound;
////++++////	int32 trightbound;
////++++////
////++++////	if(ric.direction>0){
////++++////	  tleftbound=acttag.from-ric.read.getLeftClipoff();
////++++////	  trightbound=acttag.to-ric.read.getLeftClipoff();
////++++////	}else{
////++++////	  tleftbound=ric.read.getRightClipoff()-acttag.to-1;
////++++////	  trightbound=ric.read.getRightClipoff()-acttag.from-1;
////++++////	}
////++++////	if(tleftbound<0) tleftbound=0;
////++++////	if(trightbound<0) trightbound=0;
////++++////
////++++////	int32 lpos=ric.offset+tleftbound;
////++++////	int32 rpos=ric.offset+trightbound;
////++++////	//cout << "Tagging bounds " << ric.offset+tleftbound << "\t" << ric.offset+trightbound << endl;
////++++////	// this is needed as some RMB tags might be in presently clipped
////++++////	//  parts of a read
////++++////	if(lpos < static_cast<int32>(CON_counts.size())
////++++////	   && rpos>=0 ) {
////++++////	  // adjust for partly present tags
////++++////	  if(lpos<0) lpos=0;
////++++////	  if(rpos >= static_cast<int32>(CON_counts.size())) rpos=CON_counts.size()-1;
////++++////	  //cout << "Tagging bounds cleaned" << lpos << "\t" << rpos << endl;
////++++////	  addTagToConsensus(lpos,
////++++////			    rpos,
////++++////			    '=',
////++++////			    "ORMB",
////++++////			    "Old repeat marker base",
////++++////			    true);
////++++////	}
////++++////      }
////++++////    }
////++++////  }

  FUNCEND();
}



/*************************************************************************
 *
 * actualises template ids in the contig reads from the readpool
 * also refills the template ids present vector of the contig
 * also refills the CON_readsperstrain vector
 *
 * (routine needed for backbone assembly. data loaded for backbone
 *  would conflict with data for assembly)
 *
 *************************************************************************/

void Contig::recalcTemplateIDsAndStrainPresent()
{
  FUNCSTART("void Contig::recalcTemplateIDsAndStrainPresent()");

  CON_templates_present.clear();
  CON_readsperstrain.clear();
  CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);

  try{
    auto pcrI=CON_reads.begin();
    for(; pcrI != CON_reads.end(); ++pcrI){
      if(pcrI.getOriginalReadPoolID() >= 0){
	Read & actread = const_cast<Read &>(*pcrI);
	actread.setTemplateID(CON_readpool->getRead(pcrI.getOriginalReadPoolID()).getTemplateID());
	actread.setTemplatePartnerID(CON_readpool->getRead(pcrI.getOriginalReadPoolID()).getTemplatePartnerID());
	actread.setTemplate(CON_readpool->getRead(pcrI.getOriginalReadPoolID()).getTemplate());
	if(actread.getTemplateID()>=0){
	  CON_templates_present.insert(actread.getTemplateID());
	}

	++CON_readsperstrain[CON_readpool->getRead(pcrI.getOriginalReadPoolID()).getStrainID()];
      }
    }
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  FUNCEND();
  return;
}


/*************************************************************************
 *
 * build the "rails" for a backbone (except for singlets)
 *   the rails are pseudo reads of a given length that are
 *    deduced from the consensus of the backbone contigs
 *   they are inserted into the contig to cover it without
 *    gaps. Eventually with overlaps.
 * Schematic drawing (a-f reads, 1-4 rails), overlap 2.
 *
 *   aaaaaaaaaaaaaaaaaaaaa
 *            bbbbbbbbbbbbbbbbbbbbb
 *	          cccccccccccccccccccccc
 *		   ddddddddddddddddddddddddddddddddd
 *		      eeeeeeeeeeeeeeeeeeeeee
 *		              ffffffffffffffffffffffffffffff
 *   1111111111111111
 *                 22222222222222222
 *                                333333333333333
 *                                             4444444444444
 *
 *   Aims:
 *    - allow backbones at all :-)
 *    - reduce the work for the assembler a lot
 *    - hold back automatic editor to insert or delete whole
 *      columns
 *
 * Also inserts rails as reads into the readpool (not when simulateonly
 *  is true)
 *
 * Returns number of rails inserted (or simulated to insert)
 *
 * Also initialise values in CON_counts structure for backbone info
 *
 * Additional functionality later: add the reads as "backbonelets"
 *  to circumvent gap4 maximum read length for directed assembly
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

size_t Contig::addRails(const uint32 raillength, const uint32 railoverlap, const std::string & straintxt, const bool forcestrainset, const std::string & railfromstrain, const bool simulateonly)
{
  std::vector<base_quality_t> tmp;

  if(CON_verbose && !simulateonly){
    cout << "Adding rails: length " << raillength << " and overlap " << railoverlap << endl;
  }

  return addSubsequences(raillength,
			 railoverlap,
			 straintxt,
			 forcestrainset,
			 railfromstrain,
			 false, "", tmp, "", true,
			 simulateonly);
}


size_t Contig::addSubsequences(const uint32 raillength, const uint32 railoverlap, const std::string & straintxt, const bool forcestrainset, const std::string & railfromstrain, const bool asbackbone, const std::string & bbseq, const std::vector<base_quality_t> & bbqualvec, const std::string & backbonename, const bool initccbbvalues, const bool simulateonly)
{
  FUNCSTART("void Contig::addRails()");

  size_t numrailsinserted=0;

  // find a readgroup which is (backbone|rail) and has the chosen strain
  // if there is none, create one

  ReadGroupLib::ReadGroupID rgid;
  {
    bool foundrgid=false;
    for(size_t ri=0;ri<ReadGroupLib::getNumReadGroups();++ri){
      rgid=ReadGroupLib::getReadGroupID(ri);
      if((asbackbone && rgid.isBackbone())
	 || (!asbackbone && rgid.isRail())){
	if(rgid.getStrainName()==straintxt){
	  foundrgid=true;
	  break;
	}
      }
    }
    if(!foundrgid) {
      CEBUG("rgid not found, creating new one\n");
      rgid=ReadGroupLib::newReadGroup();
      rgid.setStrainName(straintxt);
      if(asbackbone){
	rgid.setBackbone(true);
      }else{
	rgid.setRail(true);
      }
      rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);
      CEBUG("RGID is now " << rgid << endl);
    }
  }

  // This can happen with multiple backbone sequences! First backbone might create a new readgroup
  //  which the following ones do not know (yet).
  // In that case, expand!
  if(ReadGroupLib::getNumReadGroups() >= CON_readsperreadgroup.size()){
    CEBUG("expanding CON_readsperreadgroup to " << ReadGroupLib::getNumReadGroups() << endl);
    CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);
  }

  // make sure we get the latest consensus
  definalise();

  CEBUG("railfromstrain: " << railfromstrain << endl);
  CEBUG("len: " << raillength << endl);
  CEBUG("overlap: " << railoverlap << endl);


  std::string consseq;
  std::vector<base_quality_t> consqual;

  if(asbackbone) {
    consseq=bbseq;
    consqual=bbqualvec;
  }else{
    //getConsensus(consseq, consqual, true,0,0);
    if(!railfromstrain.empty()){
      int32 strainidtotake=0;
      if(!ReadGroupLib::getStrainIDOfStrain(railfromstrain,strainidtotake)){
	cout.flush();
	cout << "Error: did not find strain \"" << railfromstrain << "\" in readpool." << endl;
	MIRANOTIFY(Notify::FATAL, "Have been asked to use non-existing strain for building rails.");
      }
      //std::string strainseq;
      //std::vector<base_quality_t> strainqual;

      CEBUG("Strainidtotake: " << strainidtotake << endl);

      CEBUG(">>>> strainidtotake: " << strainidtotake << endl);

      newConsensusGet(consseq, consqual, strainidtotake);

      // TODO: eventually "merge" strainseq with consseq
      //hier weiter: strainseq und consseq mergen nach consseq

    }else{
      newConsensusGet(consseq, consqual);
    }

    //// fill quality vector with qual (if given)
    //{
    //  std::vector<base_quality_t>::iterator bI=consqual.begin();
    //  for(; bI!=consqual.end(); bI++) *bI=bbquals;
    //}
  }

  CEBUG(">>>> " << consseq.substr(0,200) << endl);

  // init the CON_counts backbone values if needed
  if(initccbbvalues && !simulateonly){
    const char * consb=&consseq[0];
    for(auto ccI=CON_counts.begin(); ccI != CON_counts.end(); ++consb, ++ccI){
      ccI->bbcountsf[0]=0;
      ccI->bbbestqualsf[0]=0;
      ccI->bbcountsr[0]=0;
      ccI->bbbestqualsr[0]=0;
      ccI->bbstrains[0]=0;
      ccI->i_backbonecharorig=toupper(*consb);
      ccI->i_backbonecharupdated='@';
    }
  }

  // keep readname here so that at the end we can easily print out the name of the last rail
  std::string readname;

  // make the railreads
  {
    std::string cutseq;
    cutseq.reserve(raillength+10);
    std::vector<base_quality_t> cutqual;
    cutqual.reserve(raillength+10);

    uint32 transi=0;
    uint32 actoffset=0;
    char tmpcons;
    for(;transi<consseq.size(); ++transi){
      tmpcons=consseq[transi];
      if(tmpcons=='@') tmpcons='N';
      cutseq+=tmpcons;
      cutqual.push_back(consqual[transi]);

      // check whether one rail is ready
      // be sure that last rail is not < ... bases
      if((cutseq.size()==raillength && consseq.size()-transi>raillength)
	 || cutseq.size()>=29900
	 || transi==consseq.size()-1){

	numrailsinserted++;

	if(!simulateonly){
	  CON_railcounter++;
	  std::ostringstream ostr;

	  if(asbackbone) {
	    ostr << "bb_" << CON_railcounter << "_" << backbonename;
	  } else {
	    ostr << "rr_####" << CON_readpool->size() << "####";
	  }

	  readname=ostr.str();

	  if(CON_verbose && numrailsinserted==1 && !simulateonly){
	    cout << getContigName() << " first rail: " << readname << endl;
	  }

	  Read & newread=CON_readpool->getRead(CON_readpool->provideEmptyRead());
	  CEBUG("C::aSs 1"<< endl);
	  newread.setReadGroupID(rgid);
	  CEBUG("C::aSs 2"<< endl);

	  try{
	    //newread.setFileNamesFromFASTAName(readname);
	    CEBUG("C::aSs 3"<< endl);
	    newread.setName(readname);
	    CEBUG("C::aSs 4"<< endl);
	    newread.setSequenceFromString(cutseq);
	    CEBUG("C::aSs 5"<< endl);
	    newread.setQualities(cutqual);
	  }
	  catch(Notify n){
	    cout << "\nOuch ... error while adding subsequence.\n";
	    cout << "cs.size: " << consseq.size() << endl;
	    cout << "transi: " << transi << endl;
	    cout << "actoffset: " << actoffset << endl;
	    cout << "readname: " << readname << endl;
	    cout << "strain: " << straintxt << endl;
	    cout << "seq: " << cutseq << endl;
	    n.handleError(THISFUNC);
	  }
	  // set CON_longestrailseen as last rail may be longer!
	  if(!asbackbone
	     && newread.getLenClippedSeq() > CON_longestrailseen) CON_longestrailseen=newread.getLenClippedSeq();

	  //CON_reads.debugDump(false);
	  CEBUG("C::aSs 6 "<< CON_readpool->size()-1 << " ao: " << actoffset << endl);
	  CON_reads.placeRead(newread,CON_readpool->size()-1,actoffset,1);
	  CEBUG("C::aSs 7"<< endl);
	  ++CON_readsperreadgroup[rgid.getLibId()];
	  CEBUG("C::aSs 8"<< endl);
	}
	cutseq.clear();
	cutqual.clear();
	if(transi!=consseq.size()-1) transi-=railoverlap;
	actoffset=transi+1;
      }
    }
  }

  if(CON_verbose && !simulateonly){
    cout << getContigName() << " last rail: " << readname << endl;
  }

  definalise();

  FUNCEND();
  return numrailsinserted;
}


//#define CEBUG(bla)



/*************************************************************************
 *
 * Remove all rails from a contig
 *
 *************************************************************************/

void Contig::removeRails()
{
  FUNCSTART("void Contig::removeRails()");

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end();){
    if(pcrI->isRail()){
      pcrI=CON_reads.removeRead(pcrI);
    }else{
      ++pcrI;
    }
  }

  priv_rebuildConCounts();
  definalise();

  FUNCEND();
  return;
}


/*************************************************************************
 *
 * Remove all non-backbone, non-rails from contig
 *
 *************************************************************************/

void Contig::stripToBackbone()
{
  FUNCSTART("void Contig::stripToBackbone()");

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end();){
    if(!pcrI->isRail() && !pcrI->isBackbone()){
      pcrI=CON_reads.removeRead(pcrI);
    }else{
      ++pcrI;
    }
  }

  // BaCh 20.8.2015: do NOT use this!
  //   priv_rebuildConCounts();
  // it's not clearing out enough the CON_counts structure
  //  and leads to double CER read coverage afterwards

  // OK, reads were deleted, now re-init the CON_counts structure
  // easiest way (and probably quickest, too):
  // memset every element (saving i_backbonecharupdated and i_backbonecharorig)
  for(auto & cce : CON_counts){
    auto tmpibbco=cce.i_backbonecharorig;
    auto tmpibbcu=cce.i_backbonecharupdated;
    memset(&cce,0,sizeof(consensus_counts_t));
    cce.i_backbonecharorig=tmpibbco;
    cce.i_backbonecharupdated=tmpibbcu;
  }

  // now update CON_count with the backbones/rails
  for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
    auto coveragemultiplier=pcrI->getDigiNormMultiplier();
    try{
      if(pcrI.getReadDirection() > 0){
	updateCountVectors(pcrI.getReadStartOffset(),
			   pcrI->getLenClippedSeq(),
			   pcrI->getClippedSeqIterator(),
			   pcrI->getSequencingType(),
			   true,
			   coveragemultiplier);
      }else{
	updateCountVectors(pcrI.getReadStartOffset(),
			   pcrI->getLenClippedSeq(),
			   pcrI->getClippedComplementSeqIterator(),
			   pcrI->getSequencingType(),
			   true,
			   coveragemultiplier);
      }
    }
    catch(Notify n){
      cout << "Readpool id: " << pcrI.getORPID() << endl;
      cout << "Error while updating this read (3):\n" << *pcrI;
      cout << "In this contig (output might crash):" << endl << *this << endl;
      n.handleError(THISFUNC);
    }
  }

  definalise();

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Contig::transformCERMappingsToCoverageReads()
{
  FUNCSTART("void Contig::transformCERMappingsToCoverageReads()");

  CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);

  // go through all short read types
  // srtype==0 Solexa

  std::string srmseq;
  srmseq.reserve(CON_counts.size());

  std::vector<base_quality_t> srmqual;
  srmqual.reserve(CON_counts.size());

  std::vector<int32> strainsinsrm;

  Read dummyread;

  for(unsigned int actsrtype=0; actsrtype<NUMMERGESEQTYPES; actsrtype++){
    // go through complete CON_counts, looking for areas which have a
    //  a short read coverage of the actsrtype

    bool needanotherpass=true;
    uint32 passnum=0;

    int32 firstrangestart = 0;
    int32 lastrangeend = CON_counts.size();

    while(needanotherpass){
      /*
      {
	cout << "DBG CON_counts:\n";
	uint32 pos=0;
	for(auto & cce : CON_counts){
	  cout << pos << ":\t" << cce << '\n';
	  ++pos;
	  if(pos==150) break;
	}
      }
      */
      passnum++;
      if ( passnum % 50 == 0 ){
	cout << "CERpass: " << passnum << " / " << CON_stats.max_coverage << " (max)" << endl;
      }
      needanotherpass=false;
      srmseq.clear();
      srmqual.clear();
      int32 rangestart=-1;       // -1, not inrange
      int32 rangeendseen=-1;         // end of the last range seen in this pass
      int32 ccpos=firstrangestart;
      bbstrainmask_t thisstrainmask=0;
      bool hasnongap=false;
      bool needsave=true;
      uint32 rangesthispass=0;
      auto ccI=CON_counts.begin();
      ccI.advance(firstrangestart);
      CEBUG("frs " << firstrangestart << "\nlre " << lastrangeend << "\n");
      for(; ccpos < lastrangeend; ++ccpos, ++ccI){
	CEBUG("ccpos: " << ccpos << '\t' << *ccI << '\n');
	if(rangestart>=0){
	  // do we need to finish?
	  if(ccI->bbstrains[actsrtype]!= thisstrainmask
	     || ccI->bbcountsf[actsrtype] + ccI->bbcountsr[actsrtype] == 0
	     || ccI->getOriginalBBChar()=='@'
	     || ccpos+1==CON_counts.size()){
	    // yes, look whether stretch is long enough or we need to
	    //  save this part
	    // the first line (hasnongap) is a fix to prevent almost empty reads
	    //  (only gaps) being added (with no strain): e.g. "X*******"
	    if(hasnongap
	       && ((needsave && ccpos-rangestart>0) || ccpos-rangestart>30)){
	      CEBUG("Would add: " << rangestart << "\t" <<ccpos);
	      CEBUG("\n" << srmseq << "\n");

	      std::ostringstream ostr;
	      if(actsrtype==0){
		ostr << "_cer_sxa_" << CON_cer_idcounter << '_';
	      }else{
		MIRANOTIFY(Notify::FATAL,"SOLiD not supported atm");
		ostr << "_cer_sid_" << CON_cer_idcounter << '_';
	      }

	      definalise();

//	      CON_reads.resize(CON_reads.size()+1);
//	      CON_reads.back().offset=-1;
//	      CON_reads.back().orpid=-1;
//	      CON_reads.back().direction=1;
//	      CON_reads.back().offset=rangestart;
//  	      CEBUG("Adding:\n" << ostr.str() << '\n');
//	      CEBUG("Orpid: " << CON_reads.back().orpid);
//	      CEBUG("\nOffset: " << CON_reads.back().offset);
//	      CEBUG("\nDirection: " << CON_reads.back().direction);
//	      Read & newread=CON_reads.back().read;

	      Read & newread=const_cast<Read &>(*(CON_reads.placeRead(dummyread,-1,rangestart,1)));

	      // we won't need adjustments, let's save some memory
	      newread.disallowAdjustments();

	      newread.setName(ostr.str());
	      newread.setSequenceFromString(srmseq);
	      newread.setLQClipoff(1);

	      if(srmqual.empty()){
		newread.setQualities(1);
		newread.setQualityFlag(false);
	      }else{
		newread.setQualities(srmqual);
	      }

	      CEBUG("\nLen: " << newread.getLenClippedSeq());
	      CEBUG(endl);

	      // set a strain, if possible
	      getMappedBBStrainIDsFromMask(strainsinsrm,thisstrainmask);

	      // find a readgroup which is CoverageEquivalent and has the chosen strain
	      // if there is none, create one

	      ReadGroupLib::ReadGroupID rgid;
	      {
		std::string nameofsis0;
		if(!strainsinsrm.empty()){
		  nameofsis0=ReadGroupLib::getStrainOfStrainID(strainsinsrm[0]);
		}
		bool foundrgid=false;
		for(size_t ri=0;ri<ReadGroupLib::getNumReadGroups();++ri){
		  rgid=ReadGroupLib::getReadGroupID(ri);
		  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA
		     && actsrtype!=0) continue;
		  if(rgid.isCoverageEquivalentRead()){
		    if(nameofsis0.empty() && rgid.getStrainName().empty()){
		      foundrgid=true;
		      break;
		    }else if(!nameofsis0.empty() && rgid.getStrainName()==nameofsis0){
		      foundrgid=true;
		      break;
		    }
		  }
		}
		if(!foundrgid) {
		  rgid=ReadGroupLib::newReadGroup();
		  rgid.setCoverageEquivalentRead(true);
		  if(!strainsinsrm.empty()) rgid.setStrainName(nameofsis0);

		  if(actsrtype==0){
		    rgid.setSequencingType(ReadGroupLib::SEQTYPE_SOLEXA);
		  }else{
		    MIRANOTIFY(Notify::INTERNAL, "Unknown actsrtype? Do we have more than Solexa and SOLiD?.");
		  }
		  CON_readsperreadgroup.resize(ReadGroupLib::getNumReadGroups(),0);
		}
	      }

	      newread.setReadGroupID(rgid);
	      ++CON_readsperreadgroup[rgid.getLibId()];

	      //assout::saveAsCAF(*this, "transfsrmdebug.caf", true);

	      if ( rangesthispass == 0 ) {
		firstrangestart = rangestart;
	      }
	      ++rangesthispass;

	      rangeendseen = ccpos + 1;
	      ++CON_cer_idcounter;
	      needanotherpass=true;
	    }
	    rangestart=-1;
	    --ccI;
	    --ccpos;
	    needsave=true;
	    srmseq.clear();
	    srmqual.clear();
	    hasnongap=false;
	  }else{
	    //// no, just append current bb char to sequence
	    if(ccI->getOriginalBBChar() != '*') hasnongap=true;
	    srmseq.push_back(tolower(ccI->getOriginalBBChar()));
	    BUGIFTHROW(ccI->bbbestqualsf[actsrtype]>=100,"Whooops 1f, qual " << static_cast<uint16>(ccI->bbbestqualsf[actsrtype]) << " is >= 100. Position " << ccI-CON_counts.begin());
	    BUGIFTHROW(ccI->bbbestqualsr[actsrtype]>=100,"Whooops 1r, qual " << static_cast<uint16>(ccI->bbbestqualsr[actsrtype]) << " is >= 100. Position " << ccI-CON_counts.begin());
	    if(passnum==1) {
	      auto newqual = static_cast<uint32>(ccI->bbbestqualsf[actsrtype])
		+ static_cast<uint32>(ccI->bbbestqualsr[actsrtype]);
	      if (newqual > 90) newqual = 90;
	      srmqual.push_back(static_cast<base_quality_t>(newqual));
	    }

	    // while thisstrainmask is set by the first position of the CER, it might be that we
	    //  started in a gap region (very unlikely, but not impossible)
	    //  Therefore, OR the current mask to thisstrainmask
	    thisstrainmask|=ccI->bbstrains[actsrtype];

	    if(ccI->bbcountsf[actsrtype]>0) {
	      ccI->bbcountsf[actsrtype]--;
	    } else if(ccI->bbcountsr[actsrtype]>0) {
	      ccI->bbcountsr[actsrtype]--;
	    } else {
	      MIRANOTIFY(Notify::INTERNAL, "I'm on a branch I shouldn't be 1. Really!");
	    }
	  }
	}else{
	  if(ccI->bbcountsf[actsrtype] + ccI->bbcountsr[actsrtype] > 0){
	    thisstrainmask=ccI->bbstrains[actsrtype];
	    rangestart=ccpos;
	    srmseq.clear();
	    srmqual.clear();
	    srmseq.push_back('X');
	    srmseq.push_back(tolower(ccI->getOriginalBBChar()));
	    //srmseq.push_back('n');
	    if(passnum==1) {
	      srmqual.push_back(0);
	      BUGIFTHROW(ccI->bbbestqualsf[actsrtype]>=100,"Whooops 2f, qual " << static_cast<uint16>(ccI->bbbestqualsf[actsrtype]) << " is >= 100. Position " << ccI-CON_counts.begin());
	      BUGIFTHROW(ccI->bbbestqualsr[actsrtype]>=100,"Whooops 2r, qual " << static_cast<uint16>(ccI->bbbestqualsr[actsrtype]) << " is >= 100. Position " << ccI-CON_counts.begin());
	      auto newqual = static_cast<uint32>(ccI->bbbestqualsf[actsrtype])
		+ static_cast<uint32>(ccI->bbbestqualsr[actsrtype]);
	      if (newqual > 90) newqual = 90;
	      srmqual.push_back(static_cast<base_quality_t>(newqual));
	    }
	  }else{
	    //if(rangesthispass>0)
	    needsave=false;
	  }
	  if(ccI->bbcountsf[actsrtype]>0) {
	    ccI->bbcountsf[actsrtype]--;
	  } else if(ccI->bbcountsr[actsrtype]>0) {
	    ccI->bbcountsr[actsrtype]--;
//	  } else {
//	    MIRANOTIFY(Notify::INTERNAL, "I'm on a branch I shouldn't be 2. Really!");
	  }
	}
      }
      if (rangeendseen >= 0 && rangeendseen < lastrangeend) {
	lastrangeend = rangeendseen;
      }
    }
  }

  /*
  {
    cout << "DBG CON_counts end:\n";
    uint32 pos=0;
    for(auto & cce : CON_counts){
      cout << pos << ":\t" << cce << '\n';
      ++pos;
      if(pos==150) break;
    }
  }
  */

  FUNCEND();
  return;
}

//#define CEBUG(bla)



/*************************************************************************
 *
 * should the readpool outside have been reordered, the contig must
 *  know about this
 * If reversemap is empty, it is recalculated ... basically the caller
 *  does not need to care about it
 *
 *************************************************************************/

void Contig::exchangeReadIDs(std::vector<uint32> & newids, std::vector<uint32> & reversemap)
{
  FUNCSTART("void Contig::exchangeReadIDs(std::vector<uint32> & newids)");

  definalise();
  CON_templates_present.clear();
  CON_allowedrefids.clear();

  if(reversemap.empty()){
    reversemap.resize(newids.size(),-1);
    for(uint32 ni=0; ni<newids.size(); ++ni){
      //cout << "newids[" <<i << "]:\t" << newids[i] << endl;
      BUGIFTHROW(newids[ni]>= reversemap.size(),"newids[ni]>= reversemap.size() ???");
      reversemap[newids[ni]]=ni;
    }
  }

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
    if(pcrI.getORPID() >= 0) {
      BUGIFTHROW(pcrI.getORPID() >= reversemap.size(),"crI->orpid >= newids.size() ???");
      pcrI.setORPID(reversemap[pcrI.getORPID()]);
    }
    if(pcrI->getTemplatePartnerID()>=0){
      const_cast<Read &>(*pcrI).setTemplatePartnerID(reversemap[pcrI->getTemplatePartnerID()]);
    }
    if(pcrI->getTemplateID()>=0){
      CON_templates_present.insert(pcrI->getTemplateID());
    }
    //cout << "exchrid: a\t" << *crI;
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::getMappedBBStrainIDsFromMask(std::vector<int32> & strains, bbstrainmask_t mask)
{
  strains.clear();
  for(uint8 i=0; i<sizeof(mask)*8; ++i){
    if(mask & 1) strains.push_back(i);
    mask=mask>>1;
  }
  return;
}



/*************************************************************************
 *
 * Helper functions to iterate through the contig in forward direction
 *  and always know which reads cover which part (fast)
 *
 * Warning: after initialising the rcci, never ever perform contig
 *  changing operations, this will invalidate the rcci!
 *
 *************************************************************************/

void Contig::rcci_t::init(const std::vector<int32> * allowedstrainidsptr, const std::vector<uint8> * allowedreadtypesptr, bool takerails, bool takebackbones, bool takereadswithoutreadpool)
{
  FUNCSTART("void Contig::rcci_t::init(const std::vector<int32> * allowedstrainidsptr, const std::vector<uint8> * allowedreadtypesptr, bool takerails, bool takebackbones, bool takereadswithoutreadpool)");

  RCCI_actcontigpos=0;
  RCCI_allowedstrainids.clear();
  RCCI_allowedreadtypes.clear();
  if(allowedstrainidsptr != nullptr) RCCI_allowedstrainids=*allowedstrainidsptr;
  if(allowedreadtypesptr != nullptr) RCCI_allowedreadtypes=*allowedreadtypesptr;
  RCCI_takerails=takerails;
  RCCI_takebackbones=takebackbones;
  RCCI_takereadswithoutreadpool=takereadswithoutreadpool;

  // a vector (read_ids_in_col) keeps the ids of the reads that are
  //  covering a specific position of the contig
  // reserving 1000 position should be enough for 99.9% of all cases,
  //  but is automatically extended by STL if needed.
  RCCI_pcrIs_in_col.clear();
  RCCI_pcrIs_in_col.reserve(1000);
  RCCI_newPCRIsonlastupdate=RCCI_pcrIs_in_col.end();

  RCCI_mpcrI = RCCI_contig->CON_reads.begin();

  update();

  FUNCEND();
}

void Contig::rcci_t::update()
{
  FUNCSTART("void Contig::readColContigIteratorUpdate(rcci_t & rcci)");

  CEBUG("rccit pos " << RCCI_actcontigpos << endl);

  // update the vector that
  //  keeps track of the reads that are
  //  covering a specific position of the contig
  // works quite simple: reads in the vector that are ending get
  //  thrown out, new reads starting are entered.

  // first delete those who fall out at this new position
  auto Ifrom=RCCI_pcrIs_in_col.begin();
  auto Ito=Ifrom;
  for(;Ifrom != RCCI_pcrIs_in_col.end(); ++Ifrom){
    *Ito=*Ifrom;
    if((*Ifrom).getReadStartOffset() + (*Ifrom)->getLenClippedSeq() > RCCI_actcontigpos) {
      ++Ito;
    }else{
      CEBUG("rccit upd thrown out " << *Ito << endl);
    }
  }
  if(Ito != Ifrom) {
    RCCI_pcrIs_in_col.resize(Ito-RCCI_pcrIs_in_col.begin(),
			     RCCI_contig->CON_reads.end());      /* the resize() will always reduce the vector,
								    but the template needs a default value anyway */

  }

  auto oldnumel=RCCI_pcrIs_in_col.size();
  // now insert ids of reads that have newly started at this position
  // Don't take railreads or backbones on demand
  for(;RCCI_mpcrI != RCCI_contig->CON_reads.end() && RCCI_mpcrI.getReadStartOffset() == RCCI_actcontigpos; ++RCCI_mpcrI){
    bool takeit=true;
    if(RCCI_mpcrI->isRail()){
      if(!RCCI_takerails) takeit=false;
    }else if(RCCI_mpcrI->isBackbone()){
      if(!RCCI_takebackbones) takeit=false;
    }else if(!RCCI_allowedstrainids.empty()){
      takeit=false;
      for(auto asie : RCCI_allowedstrainids){
	if(RCCI_mpcrI->getStrainID() == asie){
	  takeit=true;
	  break;
	}
      }
    }else if(!RCCI_allowedreadtypes.empty()){
      takeit=false;
      for(auto arte : RCCI_allowedreadtypes){
	if(RCCI_mpcrI->getSequencingType() == arte){
	  takeit=true;
	  break;
	}
      }
    }
    if(takeit){
      RCCI_pcrIs_in_col.push_back(RCCI_mpcrI);
      CEBUG("rccit upd taken " << *Ito << endl);
    }
  }

  // compute the iterator to the newly inserted elements
  RCCI_newPCRIsonlastupdate=RCCI_pcrIs_in_col.begin();
  std::advance(RCCI_newPCRIsonlastupdate,oldnumel);

  FUNCEND();
  return;
}




/*************************************************************************
 *
 * ercci is like a rcci, but instead of holding only the read ids of
 *  allowed strains & sequencing types, it holds all read ids
 *  differentiated per sequencing type per strain
 *
 *************************************************************************/

void Contig::ercci_t::init(bool takerails, bool takebackbones, uint32 numstrains)
{
  FUNCSTART("void Contig::ercci_t::init(bool takerails, bool takebackbones, uint32 numstrains)");

  BUGIFTHROW(numstrains==0,"numstrains==0???");

  ERCCI_actcontigpos=0;
  ERCCI_takebackbones=takebackbones;
  ERCCI_takerails=takerails;
  ERCCI_mpcrI = ERCCI_contig->CON_reads.begin();

  // per sequencing type and per strain ids,
  //  a vector keeps the ids of the reads that are
  //  covering a specific position of the contig
  // reserving 1000 position should be enough 99.9% of all cases,
  //  is automatically extended by STL if needed.

  ERCCI_pcrI_st_st.clear();
  ERCCI_pcrI_st_st.resize(ReadGroupLib::getNumSequencingTypes());
  for(uint32 seqtype=0; seqtype<ReadGroupLib::getNumSequencingTypes(); seqtype++){
    ERCCI_pcrI_st_st[seqtype].resize(numstrains);
    for(uint32 strainid=0; strainid<numstrains; ++strainid){
      ERCCI_pcrI_st_st[seqtype][strainid].clear();
      ERCCI_pcrI_st_st[seqtype][strainid].reserve(1000);
    }
  }

  update();

  FUNCEND();
  return;
}

void Contig::ercci_t::update()
{
  FUNCSTART("void Contig::ercci_t::update()");

  // update the vector that
  //  keeps track of the reads that are
  //  covering a specific position of the contig
  // works quite simple: reads in the vector that are ending get
  //  thrown out, new reads starting are entered.

  // first delete those who fall out at this new position

  CEBUG("rccite pos " << ERCCI_actcontigpos << endl);

  for(uint32 seqtype=0; seqtype<ReadGroupLib::getNumSequencingTypes(); ++seqtype){
    for(uint32 strainid=0; strainid<ERCCI_pcrI_st_st[seqtype].size(); ++strainid){
      auto Ifrom=ERCCI_pcrI_st_st[seqtype][strainid].begin();
      auto Ito=Ifrom;
      for(;Ifrom != ERCCI_pcrI_st_st[seqtype][strainid].end(); ++Ifrom){
	*Ito=*Ifrom;
	if((*Ifrom).getReadStartOffset()+(*Ifrom)->getLenClippedSeq() > ERCCI_actcontigpos) {
	  ++Ito;
	}else{
	  CEBUG("rccite upd thrown out " << *Ito << endl);
	}
      }
      if(Ito != Ifrom) {
	ERCCI_pcrI_st_st[seqtype][strainid].resize(Ito-ERCCI_pcrI_st_st[seqtype][strainid].begin(),
						   ERCCI_contig->CON_reads.end());    /* the resize() will always reduce the vector,
											 but the template needs a default value anyway */
      }
    }
  }

  // add in the ones which are arriving anew
  for(;ERCCI_mpcrI != ERCCI_contig->CON_reads.end() && ERCCI_mpcrI.getReadStartOffset() == ERCCI_actcontigpos; ++ERCCI_mpcrI){
    if(ERCCI_mpcrI->isRail() && !ERCCI_takerails) continue;
    if(ERCCI_mpcrI->isBackbone() && !ERCCI_takebackbones) continue;
    uint32 seqtype=ERCCI_mpcrI->getSequencingType();
    uint32 strainid=static_cast<uint32>(ERCCI_mpcrI->getStrainID());
    // better safe than sorry: the following two BUGIFTHROW should not trigger, never ever at all
    //  but due to an upstream error with strain ids ... this seg fault took two freaking weeks to pin down :-(((((((((
    BUGIFTHROW(seqtype>=ERCCI_pcrI_st_st.size(),"seqtype " << seqtype << " >=ERCCI_pcrI_st_st.size() " << ERCCI_pcrI_st_st.size() << " ???");
    BUGIFTHROW(strainid>=ERCCI_pcrI_st_st[seqtype].size(),"strainid " << strainid << " >=ERCCI_pcrI_st_st[seqtype].size() " << ERCCI_pcrI_st_st[seqtype].size() << " ???");
    ERCCI_pcrI_st_st[seqtype][strainid].push_back(ERCCI_mpcrI);
    CEBUG("rccite upd taken " << ERCCI_mpcrI << endl);
  }

  FUNCEND();
  return;
}


void Contig::ercci_t::advance()
{
  ++ERCCI_actcontigpos;
  update();
}





/*************************************************************************
 *
 * Simply go through all reads and delet tags with specified tagid
 *
 *
 *************************************************************************/

void Contig::deleteTagsInReads(const multitag_t::mte_id_t identifier)
{
  for(auto & pcre : CON_reads){
    const_cast<Read &>(pcre).deleteTag(identifier);
  }

  return;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::setupAsBackBoneContig()
{
  FUNCSTART("void Contig::setupAsBackBoneContig()");

  // first, trash all consensus caches so that MIRA uses newly computed consensus
  //  sequences later on. Especially important when working with several strains
  trashConsensusCache(false);

  std::vector<bool> readgroupspresent(ReadGroupLib::getNumReadGroups(),false);
  std::vector<bool> bbreadgroupspresent(readgroupspresent);
  std::unordered_set<std::string> strainnamesinbb;

  CON_isbackbonecontig=true;

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
    if(pcrI.getORPID()==-1
       || CON_readpool->getRead(pcrI.getORPID()).isRail()
       || CON_readpool->getRead(pcrI.getORPID()).isCoverageEquivalentRead()) continue;
    readgroupspresent[static_cast<uint32>(CON_readpool->getRead(pcrI.getORPID()).getReadGroupID().getLibId())]=true;
    readgroupspresent[static_cast<uint32>(pcrI->getReadGroupID().getLibId())]=true;

    if(!pcrI->getStrainName().empty()){
      strainnamesinbb.insert(pcrI->getStrainName());
    }
  }

  if(strainnamesinbb.empty()){
    MIRANOTIFY(Notify::FATAL,"Contig " << getContigName() << " has no strainnames given it its data. Not good, make sure one is present (try putting 'strainname=' in manifest file).");
  }

  for(auto & sne : strainnamesinbb){
    cout << "Contig " << getContigName() << " has strain " << sne << endl;
  }

//  if(strainnamesinbb.size()>1){
//    MIRANOTIFY(Notify::FATAL,"Contig " << getContigName() << " has " << strainnamesinbb.size() << " strainnames. Not good, make sure only one is present (put 'strainname=' in manifest file).");
//  }

  const std::string & strainname=*strainnamesinbb.begin();

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
    if(pcrI.getORPID()==-1
       || CON_readpool->getRead(pcrI.getORPID()).isRail()
       || CON_readpool->getRead(pcrI.getORPID()).isCoverageEquivalentRead()) continue;
    CON_readpool->getRead(pcrI.getORPID()).getReadGroupID().setBackbone(true);
    pcrI->getReadGroupID().setBackbone(true);

    if(pcrI->getStrainName().empty()){
      CON_readpool->getRead(pcrI.getORPID()).getReadGroupID().setStrainName(strainname);
      pcrI->getReadGroupID().setStrainName(strainname);
    }
  }


  // setup quick lookup for -CO:msrkceu (keep contig ends unmerged)
  CON_mpindex_msrkceu_left=CON_markerpositions.size();
  CON_markerpositions.push_back(0);
  CON_mpindex_msrkceu_right=CON_markerpositions.size();
  CON_markerpositions.push_back(getContigLength());


  // setup data for areas in genome where we want to force merge
  // how memory-inefficient having it done by buildMaskShadow()
  //  instead of writing directly to CON_counts
  // but I don't care atm
  CON_hasforcemergeareas=false;
  {
    for(auto & cce : CON_counts){
      cce.forcemergearea=0;
    }

    std::vector<int8> maskshadow;
    std::vector<multitag_t::mte_id_t> masktagstrings;
    masktagstrings.push_back(Read::REA_tagentry_idMFSM);

    buildMaskShadow(maskshadow,masktagstrings,true);
    auto mI=maskshadow.cbegin();
    for(auto ccI=CON_counts.begin(); ccI != CON_counts.end(); ++mI, ++ccI){
      ccI->forcemergearea=*mI;
      CON_hasforcemergeareas|=*mI;
    }
  }

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Contig::errorstatus_t::dumpStatus(bool longmsg, const char * additionalmsg)
{
  FUNCSTART("void Contig::errorstatus_t::dumpStatus(const char * additionalmsg)");

  if(additionalmsg != nullptr) cout << additionalmsg;

  switch(code) {
  case Contig::ENOERROR : {
    if(longmsg){
      cout << "\t+\n";
    }else{
      cout << '+';
    }
    break;
  }
  case Contig::EZEROLEN : {
    if(longmsg){
      cout << "\t-\tzero len\n";
    }else{
      cout << '0';
    }
    break;
  }
  case Contig::ENOTCALLED : {
    if(longmsg){
      cout << "\t-\tnot called\n";
    }else{
      cout << '*';
    }
    break;
  }
  case Contig::ENOALIGN : {
    if(longmsg){
      cout << "\t-\tno align found\n";
    }else{
      cout << 'a';
    }
    break;
  }
  case Contig::EDROPINRELSCORE : {
    if(longmsg){
      cout << "\t-\tdrop in relscore too high\n";
    }else{
      cout << 'd';
    }
    break;
  }
  case Contig::EALIGNREJECTMRS : {
    if(longmsg){
      cout << "\t-\talign rejected by mrs\n";
    }else{
      cout << 'A';
    }
    break;
  }
  case Contig::ETEMPLATEDIRECTION : {
    if(longmsg){
      cout << "\t-\ttemplate in wrong direction\n";
    }else{
      cout << 'T';
    }
    break;
  }
  case Contig::ESEGMENTPLACEMENT : {
    if(longmsg){
      cout << "\t-\tplacement of segment wrong\n";
    }else{
      cout << 'P';
    }
    break;
  }
  case Contig::ETEMPLATESIZELT : {
    if(longmsg){
      cout << "\t-\tmismatch in template size (<)\n";
    }else{
      cout << '<';
    }
    break;
  }
  case Contig::ETEMPLATESIZEGT : {
    if(longmsg){
      cout << "\t-\tmismatch in template size (>)\n";
    }else{
      cout << '>';
    }
    break;
  }
  case Contig::ESRMBMISMATCH : {
    if(longmsg){
      cout << "\t-\tmismatch in SRMB zone\n";
    }else{
      cout << 'R';
    }
    break;
  }
  case Contig::ESPECIALSRADDFAIL : {
    if(longmsg){
      cout << "\t-\tfailed special SR add rules\n";
    }else{
      cout << 'c';
    }
    break;
  }
  case Contig::EREFIDNOTALLOWED : {
    if(longmsg){
      cout << "\t-\trefid not allowed\n";
    }else{
      cout << 'r';
    }
    break;
  }
  case Contig::EMAXCOVERAGEREACHED : {
    if(longmsg){
      cout << "\t-\tmaxcoverage reached\n";
    }else{
      cout << 'x';
    }
    break;
  }
  case Contig::EDANGERZONE : {
    if(longmsg){
      cout << "\t-\ttoo many mismatches in danger zone(s)\n";
    }else{
      cout << 'z';
    }
    break;
  }
  case Contig::EFORCEDGROWTHNOTREACHED : {
    if(longmsg){
      cout << "\t-\tforced growth not reached\n";
    }else{
      cout << 'g';
    }
    break;
  }
  case Contig::EGROWTHNOTALLOWED : {
    if(longmsg){
      cout << "\t-\tgrowth not allowed\n";
    }else{
      cout << 'G';
    }
    break;
  }
  case Contig::ERELEGATEDBYPP : { // actually set by pathfinder
    if(longmsg){
      cout << "\t-\trelegated\n";
    }else{
      cout << '\\';
    }
    break;
  }
  case Contig::EUNSPECIFIED : {
    if(longmsg){
      cout << "\t-\tunspecified reject\n";
    }else{
      cout << '?';
    }
    break;
  }
  default : {
    MIRANOTIFY(Notify::INTERNAL, "Unknown errorcode from the contig object.");
  }
  }
}
