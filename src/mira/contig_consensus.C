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

#include "util/dptools.H"
#include "util/timer.H"

#include "util/stlimprove.H"

#include <boost/lexical_cast.hpp>


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


// atm this slows down execution a bit
//#define CLOCK_STEPS_CONS
// atm this slows down execution ... a lot
//#define CLOCK_STEPS_CONSSUB


#define VCOUT(bla)   {if(CON_verbose) {cout << bla;} }

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {if(CON_cebugflag) cout << bla;}
//#define CEBUG(bla) {cout << bla;}


// new ... calculate main consensus and cache it
// consensus for strains: moved as "on-demand" calulation to getConsensus()

//#define CEBUG(bla) {cout << bla;}
void Contig::calcConsensi(int32 mincoverage, base_quality_t minqual, char missingcoveragechar)
{
  FUNCSTART("void Contig::calcConsensi(int32 mincoverage, base_quality_t minqual, char missingcoveragechar)");

  //bool assumediploid=true;
  bool assumediploid=false;
  CON_conscalc_assumediploid=assumediploid;
  bool allowiupac=!(*CON_miraparams)[0].getContigParams().con_force_nonIUPACconsensus;
  CON_conscalc_allowiupac=allowiupac;

  //CON_cebugflag=true;

  CEBUG("calcConsensi(). mincov: " << mincoverage << "\tminqual: " << static_cast<uint16>(minqual) << "\tmissingcchar: " << missingcoveragechar << endl);

  // loading from different files (backbones reads etc) make the number of strain change over time
  // i.e., contigs loaded earlier may have a smaller CON_readsperstrain.size() than is good
  if(CON_readsperstrain.size() < ReadGroupLib::getNumOfStrains()){
    CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);
  }

  CEBUG("Rebuild for " << CON_readsperstrain.size() << " strains in readpool\n");

  if(CON_conscalc_mincov!=mincoverage
     || CON_conscalc_minqual!=minqual
     || CON_conscalc_missingchar!=missingcoveragechar
     || CON_strainconsseq.size()==0
     || CON_readsperstrain.size() > CON_strainconsseq.size()){    // BaCh 17.09.2013: bugfix ">=" into ">", else always needless recalcs

    CEBUG("need recalc.:\n");
    CEBUG("old values: " << CON_conscalc_mincov << " " << static_cast<uint16>(CON_conscalc_minqual) << " " << CON_conscalc_missingchar << " " << CON_strainconsseq.size() << " " << CON_readsperstrain.size() << endl);
    CON_conscalc_mincov=mincoverage;
    CON_conscalc_minqual=minqual;
    CON_conscalc_missingchar=missingcoveragechar;
    CEBUG("new values: " << CON_conscalc_mincov << " " << static_cast<uint16>(CON_conscalc_minqual) << " " << CON_conscalc_missingchar << " " << CON_strainconsseq.size() << " " << CON_readsperstrain.size() << endl);

    // make a consensus for every strain
    // cache them in CON_allconsseq and *qual
    makeIntelligentConsensus(CON_allconsseq,
			     CON_allconsqual,
			     &CON_alladjustments,
			     0,
			     CON_counts.size(),
			     -1,                   // strainid
			     CON_conscalc_mincov,
			     CON_conscalc_minqual,
			     missingcoveragechar,
			     assumediploid,
			     allowiupac,
			     false             // no addconstag, as strainid == -1 anyway
      );

    CEBUG("CON_readsperstrain.size(): " << CON_readsperstrain.size() << endl);

    CON_strainconsseq.clear();
    CON_strainconsqual.clear();
    CON_strainadjustments.clear();
    CON_strainconsseq.resize(CON_readsperstrain.size());
    CON_strainconsqual.resize(CON_readsperstrain.size());
    CON_strainadjustments.resize(CON_readsperstrain.size());

    uint32 numstrains=0;
    for(uint32 si=0; si<CON_readsperstrain.size(); ++si){
      if(CON_readsperstrain[si]>0) ++numstrains;
    }

    // now for all strains.
    // strains not present in contig will have sequence, quality and adjustments
    //  pre-filled with default values (@, 0 and -1)
    // strains present will have empty seq+qual+adj ... to be calculated on demand

    // Change: no precalculated strains, only pre-filled for clear cases!
    // Should getConsensus() ask for them later, they
    //  will be calculated on demand (using CON_conscalc_* values)
    for(uint32 si=0; si<CON_readsperstrain.size(); ++si){
      if(CON_readsperstrain[si]==0) {
	CON_strainconsseq[si].resize(CON_allconsseq.size(),missingcoveragechar);
	CON_strainconsqual[si].resize(CON_allconsqual.size(),0);
	CON_strainadjustments[si].resize(CON_allconsqual.size(),-1);
      }else{
	if(numstrains==1 && mincoverage==0 && minqual==0){
	  // take over from allstrains in this very special case
	  CON_strainconsseq[si]=CON_allconsseq;
	  CON_strainconsqual[si]=CON_allconsqual;
	  CON_strainadjustments[si]=CON_alladjustments;
	}
      }
    }
  }else{
  }

  FUNCEND();

  return;
}
//#define CEBUG(bla)


// makes sure consensus and all adjoining structures are calculated and valid
// works like newConsensusGet() below, but does not return the consensus nor qualities
//
// trick: calling newConsensusGet() with the src==target will make sure things are not
//  unnecessarily copied
void Contig::ensureConsensus(int32 strainidtotake)
{
  FUNCSTART("void Contig::ensureConsensus(int32 strainidtotake)");

  CEBUG("ensureConsensus(): " << strainidtotake << endl);
  if(CON_allconsseq.empty() || CON_strainconsseq.empty()){
    CEBUG("ensureConsensus(): cons empty, need recalc" << endl);
    calcConsensi(0,0,'X');
  }
  if(strainidtotake<0){
    newConsensusGet(CON_allconsseq, CON_allconsqual, strainidtotake);
  }else{
    newConsensusGet(CON_strainconsseq[strainidtotake], CON_strainconsqual[strainidtotake], strainidtotake);
  }
  FUNCEND();
}

void Contig::newConsensusGet(std::string & target, std::vector<base_quality_t> & qual, int32 strainidtotake)
{
  FUNCSTART("void Contig::newConsensusGet(std::string & target, std::vector<base_quality_t> & qual, int32 strainidtotake)");

//  if(CON_abortflag) {
//    uint16 * bombme=nullptr;
//    *bombme=0xdead;
//    cout << bombme;
//  }

  CEBUG("newConsensusGet(): gimme strain " << strainidtotake << endl);

  // loading from different files (backbones reads etc) make the number of strain change over time
  // i.e., contigs loaded earlier may have a smaller CON_readsperstrain.size() than is good
  if(CON_readsperstrain.size() < ReadGroupLib::getNumOfStrains()){
    CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);
  }

  BUGIFTHROW(strainidtotake>=static_cast<int32>(CON_readsperstrain.size()),
	     "strainidtotake " << strainidtotake << " >= CON_readsperstrain.size() "
	     << static_cast<int32>(CON_readsperstrain.size()) << " ?");

  if(CON_fixedconsseq.size() && CON_fixedconsqual.size() && strainidtotake<0){
    CEBUG("newConsensusGet(): get fixed" << endl);
    target=CON_fixedconsseq;
    qual=CON_fixedconsqual;
  }else{
    CEBUG("CON_allconsseq.size(): " << CON_allconsseq.size() << endl);
    CEBUG("CON_strainconsseq.size(): " << CON_strainconsseq.size() << endl);
    if(CON_allconsseq.empty() || CON_strainconsseq.empty()
       || (strainidtotake>=0 && strainidtotake>=CON_strainconsseq.size())){
      CEBUG("something's empty, need recalc" << endl);
      calcConsensi(0,0,'X');
    }
    if(strainidtotake<0){
      CEBUG("newConsensusGet(): get allcons" << endl);
      CEBUG("CON_allconsseq.size(): " << CON_allconsseq.size() << endl);
      CEBUG("CON_allconsqual.size(): " << CON_allconsqual.size() << endl);
      target=CON_allconsseq;
      qual=CON_allconsqual;
    }else{
      CEBUG("newConsensusGet(): get strain " << strainidtotake << endl);
      BUGIFTHROW(strainidtotake>=CON_strainconsseq.size(),"something's utterly wrong: strainidtotake>=CON_strainconsseq.size() ???");
      // on demand calculation
      if(CON_strainconsseq[strainidtotake].empty()){
	CEBUG("check on demand calculation\n");
	uint32 numstrains=0;
	for(uint32 si=0; si<CON_readsperstrain.size(); si++){
	  if(CON_readsperstrain[si]>0) numstrains++;
	}
	if(numstrains==1){
	  CEBUG("only 1 strain, can take main consensus\n");
	  CON_strainconsseq[strainidtotake]=CON_allconsseq;
	  CON_strainconsqual[strainidtotake]=CON_allconsqual;
	  CON_strainadjustments[strainidtotake]=CON_alladjustments;
	}else{
	  CEBUG("must do calculation\n");
	  makeIntelligentConsensus(CON_strainconsseq[strainidtotake],
				   CON_strainconsqual[strainidtotake],
				   &CON_strainadjustments[strainidtotake],
				   0,
				   CON_counts.size(),
				   strainidtotake,
				   CON_conscalc_mincov,
				   CON_conscalc_minqual,
				   CON_conscalc_missingchar,
				   CON_conscalc_assumediploid,
				   CON_conscalc_allowiupac,
				   CON_conscalc_addconstag);
	}
      }else{
	CEBUG("take cached\n");
      }
      target=CON_strainconsseq[strainidtotake];
      qual=CON_strainconsqual[strainidtotake];
    }
  }

  BUGIFTHROW(target.size() != CON_counts.size(),"strainidtotake " << strainidtotake << ": target.size() " << target.size() << " != CON_counts.size() " << CON_counts.size());

  FUNCEND();
}
//#define CEBUG(bla)


//void Contig::OLDgetConsensus1(std::string & target, std::vector<base_quality_t> & qual, bool markspecials, int32 mincoverage, base_quality_t minqual, int32 strainidtotake, char missingcoveragechar, ostream * ostr, bool contagsintcs)
//{
//  FUNCSTART("void Contig::getConsensus(std::string & target, std::vector<base_quality_t> & qual, bool markspecials, int32 mincoverage, int32 strainidtotake, ostream * ostr, bool contagsintcs)");
//
//  CEBUG("getCons()\n");
//
//  CON_cebugflag=true;
//
//  if(CON_cheat_intelcons.empty()
//     || CON_cheat_intelcons_markspecials!=markspecials
//     || CON_cheat_intelcons_mincov!=mincoverage
//     || CON_cheat_intelcons_minqual!=minqual
//     || CON_cheat_intelcons_strainidtotake!=strainidtotake
//     || ostr != nullptr) {
//    bool mustcompute=true;
//    if(!CON_cheat_intelcons.empty()
//       && CON_cheat_intelcons_mincov==mincoverage
//       && CON_cheat_intelcons_minqual==minqual
//       && CON_cheat_intelcons_strainidtotake == strainidtotake) {
//      if(markspecials) {
//	if(CON_cheat_intelcons_markspecials) {
//	  mustcompute=false;
//	}
//      } else {
//	mustcompute=false;
//	if(CON_cheat_intelcons_markspecials) {
//	  // just make the sequence make uppercase
//	  for(uint32 i=0; i<CON_cheat_intelcons.size(); i++) {
//	    CON_cheat_intelcons[i]=toupper(CON_cheat_intelcons[i]);
//	  }
//	}
//      }
//    }
//    if(mustcompute || ostr != nullptr) {
//      makeIntelligentConsensus(CON_cheat_intelcons,
//			       CON_cheat_intelconsqual,
//			       0,
//			       CON_counts.size(),
//			       markspecials,
//			       mincoverage,
//			       minqual,
//			       strainidtotake,
//			       missingcoveragechar,
//			       ostr,
//			       contagsintcs);
//    }
//    CON_cheat_intelcons_markspecials=markspecials;
//    CON_cheat_intelcons_mincov=mincoverage;
//    CON_cheat_intelcons_minqual=minqual;
//    CON_cheat_intelcons_strainidtotake=strainidtotake;
//  }
//
//  target=CON_cheat_intelcons;
//  qual=CON_cheat_intelconsqual;
//
//  FUNCEND();
//
//  return;
//}






/*************************************************************************
 *
 * if the routines to decide for a base (helper2 routines) could not get
 *  clear base but the user wants one, this is doing a shootout based on
 *  majority vote.
 *
 * groups with forward and reverse count double in the read count
 *
 * if majority vote still does not work (all the same), then the last
 *  one wins. I.e., * takes precedence over T, this over G, over C, over A
 *
 *************************************************************************/

// MIC all CEBUG start
//#define CEBUG(bla) {cout << bla;}

void Contig::makeIntelligentConsensus_helper3(char & thisbase, base_quality_t & thisqual, const std::vector<nngroups_t> & groups, const std::vector<char> & IUPACbasegroups)
{
  uint32 maxcount=0;

  for(uint32 actgroup=0; actgroup<groups.size(); actgroup++){
    for(uint32 actbase=0; actbase<IUPACbasegroups.size(); actbase++){
      if(IUPACbasegroups[actbase]==groups[actgroup].base){
	uint32 groupcount=groups[actgroup].urdids.size();
	if(groups[actgroup].forwarddircounter>0
	   && groups[actgroup].complementdircounter>0){
	  groupcount*=2;
	}
	if(groupcount>=maxcount){
	  maxcount=groupcount;
	  thisbase=groups[actgroup].base;
	  thisqual=groups[actgroup].groupquality;
	}
      }
    }
  }

}


/*************************************************************************
 *
 * Return by value:
 *  num valid groups
 * Return by variable:
 *  thisbase/qual: chosen base/gap (one of ACGT*) and its qual if single choice were to be made
 *  iupacbasegroups: contains all bases (any combination of ACGT*) which could be valid
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}
void Contig::makeIntelligentConsensus_helper2_calcSOLEXA(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups, const bbstrainmask_t strainmask)
{
  FUNCSTART("uint8 Contig::makeIntelligentConsensus_helper2_calcSOLEXA(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups, const bbstrainmask_t strainmask)");

  thisbase=' ';
  thisqual=0;
  IUPACbasegroups.clear();

  bool hasmergedbases=false;
  if(CON_isbackbonecontig){
    // Ok, check whether the merged bases belong to this strain
    //  if not, well then no merged bases exist
    hasmergedbases=(ccI->getOriginalBBChar()!='@') & (ccI->bbcountsf[0] + ccI->bbcountsr[0] > 0);
    CEBUG("obc: " << ccI->getOriginalBBChar()
	  << "\tbbcf[0]: " << ccI->bbcountsf[0]
	  << "\tbbcr[0]: " << ccI->bbcountsr[0]
	  << '\n');
    CEBUG("hmb1: " << hasmergedbases << '\n');
    if(!(ccI->bbstrains[0] & strainmask)) hasmergedbases=false;
    CEBUG("hmb2: " << hasmergedbases << '\n');
  }

  CEBUG("bbchar: " << ccI->getOriginalBBChar() << "\tbbcounts: " << ccI->bbcountsf[0] + ccI->bbcountsr[0] << "\tbbbestquals: " << static_cast<uint16>(ccI->bbbestquals[0]) << "\tbbstrains: " << std::hex << static_cast<uint16>(ccI->bbstrains[0]) << std::dec << "\tHasmergedb: " << hasmergedbases <<'\n');
  CEBUG("Strainmask: " << static_cast<uint64>(strainmask) << '\n');

  bool groupschosen[groups.size()];  // init in for loop below
  uint32 groupcounts[groups.size()]; // init in for loop below
  uint32 maxcount=0;

  for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
    uint32 groupcount=groups[actgroup].urdids.size();
    if(hasmergedbases && ccI->getOriginalBBChar() == groups[actgroup].base) groupcount+=ccI->bbcountsf[0] + ccI->bbcountsr[0];
    if(groups[actgroup].forwarddircounter>=1
       && groups[actgroup].complementdircounter>=1){
      if(groups[actgroup].forwarddircounter>=2
	 && groups[actgroup].complementdircounter>=2){
	groupcount*=2;
      }else{
	groupcount*=3;
	groupcount/=2;
      }
    }
    if(groupcount>=maxcount){
      maxcount=groupcount;
    }
    groupcounts[actgroup]=groupcount;
    groupschosen[actgroup]=false;
  }

  CEBUG("GROUPCOUNTS: ");
  for(uint32 i=0;i<5;++i){
    CEBUG(groupcounts[i] << " ");
  }
  CEBUG(endl);

  for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
    if(groupcounts[actgroup]>=maxcount/2
       || (groups[actgroup].forwarddircounter>=2
	   && groups[actgroup].complementdircounter>=2)){
      IUPACbasegroups.push_back(groups[actgroup].base);
      groupschosen[actgroup]=true;
    }
  }

  CEBUG("GROUPSCHOSEN: ");
  for(uint32 i=0;i<5;++i){
    CEBUG(groupschosen[i] << " ");
  }
  CEBUG(endl);

  BUGIFTHROW(IUPACbasegroups.size()==0,"acp: " << actcontigpos << ": num groups ==0 ???");

  for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
    if(groupcounts[actgroup]==maxcount){
      thisbase=groups[actgroup].base;
      thisqual=groups[actgroup].groupquality;
      break;
    }
  }

  BUGIFTHROW(IUPACbasegroups.size()==0,"IUPACbasegroups.size()==0");

  if(IUPACbasegroups.size()==1){
    // need to do something?
  }else{
    uint32 tqual=0;
    CEBUG("singleding\n");
    for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
      CEBUG("Newgroup\n" << groups[actgroup]);
      if(groupschosen[actgroup]==true){
	CEBUG("CHOSEN!\n\n");
	tqual+=groups[actgroup].groupquality;
      }
    }
    tqual/=IUPACbasegroups.size();
    if(tqual>90) tqual=90;
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

void Contig::makeIntelligentConsensus_helper2_calcPACBIOHQ(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups)
{

  // no info. atm, same as 454 without 40:60 rule

  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?
  // TODO: write PacBioLQ/HQ routine

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=std::max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_PACBIOHQ].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
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

void Contig::makeIntelligentConsensus_helper2_calc454GS20(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups)
{
  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

  // if max count is gap, but there are other bases
  if(maxsize_i==4 && runnerup>0){
    // apply 40:60 rule
    if(100*runnerup/(maxsize+runnerup) >= 40){
      std::swap(maxsize,runnerup);
      std::swap(maxsize_i,runnerup_i);
      std::swap(maxqual,runnerupqual);
    }
  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=std::max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }
  }

  //// TODO: testing
  //if(maxsize_i>=0 && runnerup_i>=0){
  //  if(maxsize+runnerup == 0){
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454","DOH?");
  //  }else if(100*runnerup/(maxsize+runnerup) >= 30){
  //    //ostringstream ostr;
  //    //ostr << static_cast<char>(groups[maxsize_i].base) << ": " << maxsize;
  //    //ostr << " " << static_cast<char>(groups[runnerup_i].base) << ": " << runnerup;
  //    //ostr << "  -  " << 100*runnerup/(maxsize+runnerup) << "%";
  //    //
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454",ostr.str().c_str());
  //  }
  //}

  return;
}



/*************************************************************************
 *
 * atm a pure copy of 454
 *
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calcIonTorrent(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups)
{
  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

  // if max count is gap, but there are other bases
  // Test, maybe IonTorrent is a bit different.
  // Actually, it's almost the best rule there is ...
  //  35 perhaps better, but needs to be verified.
  // Astonishing
  if(maxsize_i==4 && runnerup>0){
    // apply 40:60 rule
    if(100*runnerup/(maxsize+runnerup) >= 40){
      std::swap(maxsize,runnerup);
      std::swap(maxsize_i,runnerup_i);
      std::swap(maxqual,runnerupqual);
    }
  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=std::max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }
  }

  //// TODO: testing
  //if(maxsize_i>=0 && runnerup_i>=0){
  //  if(maxsize+runnerup == 0){
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454","DOH?");
  //  }else if(100*runnerup/(maxsize+runnerup) >= 30){
  //    //ostringstream ostr;
  //    //ostr << static_cast<char>(groups[maxsize_i].base) << ": " << maxsize;
  //    //ostr << " " << static_cast<char>(groups[runnerup_i].base) << ": " << runnerup;
  //    //ostr << "  -  " << 100*runnerup/(maxsize+runnerup) << "%";
  //    //
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454",ostr.str().c_str());
  //  }
  //}

  return;
}




/*************************************************************************
 *
 * adapted copy of new solexa routine
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calcText(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups)
{
  FUNCSTART("void Contig::makeIntelligentConsensus_helper2_calcText(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups)");

  // Idea: pure coverage. Every group >= maximum_cov/2  is considered.

  thisbase=' ';
  thisqual=0;
  IUPACbasegroups.clear();

  bool groupschosen[groups.size()];  // init in for loop below
  uint32 groupcounts[groups.size()]; // init in for loop below
  uint32 maxcount=0;

  for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
    uint32 groupcount=groups[actgroup].urdids.size();
    if(groupcount>=maxcount){
      maxcount=groupcount;
    }
    groupcounts[actgroup]=groupcount;
    groupschosen[actgroup]=false;
  }

  CEBUG("GROUPCOUNTS: ");
  for(uint32 i=0;i<5;++i){
    CEBUG(groupcounts[i] << " ");
  }
  CEBUG(endl);

  for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
    if(groupcounts[actgroup]>=maxcount/2){
      IUPACbasegroups.push_back(groups[actgroup].base);
      groupschosen[actgroup]=true;
    }
  }

  CEBUG("GROUPSCHOSEN: ");
  for(uint32 i=0;i<5;++i){
    CEBUG(groupschosen[i] << " ");
  }
  CEBUG(endl);

  BUGIFTHROW(IUPACbasegroups.size()==0,"acp: " << actcontigpos << ": num groups ==0 ???");

  for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
    if(groupcounts[actgroup]==maxcount){
      thisbase=groups[actgroup].base;
      thisqual=groups[actgroup].groupquality;
      break;
    }
  }

  BUGIFTHROW(IUPACbasegroups.size()==0,"IUPACbasegroups.size()==0");

  if(IUPACbasegroups.size()==1){
    // need to do something?
  }else{
    uint32 tqual=0;
    CEBUG("singleding\n");
    for(uint32 actgroup=0; actgroup<groups.size(); ++actgroup){
      CEBUG("Newgroup\n" << groups[actgroup]);
      if(groupschosen[actgroup]==true){
	CEBUG("CHOSEN!\n\n");
	tqual+=groups[actgroup].groupquality;
      }
    }
    tqual/=IUPACbasegroups.size();
    if(tqual>90) tqual=90;
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

void Contig::makeIntelligentConsensus_helper2_calcSangerQual(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const std::vector<nngroups_t> & groups, std::vector<char> & IUPACbasegroups, const base_quality_t maxqual, const uint32 maxcount)
{
  // Idea: all groups with qual >= 30
  //    or groups with qual within (X) (with X==5 now) of maxqual
  //    or groups with count >= 3 and quality >=20 and forward/reverse

  IUPACbasegroups.clear();
  int32 avgqual=0;
  int32 groupstaken=0;
  int32 verygoodgroups=0;

  // New: if force non-IUPAC consensus is wished
  // decrease the level with which we look at the groups
  // gradually until we have to find something (qual == 0)
  base_quality_t goodgroupqual=35;
  base_quality_t lessergroupqual=goodgroupqual;
  if(lessergroupqual>=10) lessergroupqual-=10;

  do{
    if(goodgroupqual>=5) goodgroupqual-=5;
    if(lessergroupqual>=5) lessergroupqual-=5;
    for(uint32 i=0; i<groups.size(); i++){
      if(groups[i].urdids.size()
	 &&(groups[i].groupquality >= goodgroupqual
	    || groups[i].groupquality+5 >= maxqual
	    || (groups[i].urdids.size() >= 3          // TODO: check 2 or 3 (was 3, testing with 2)
		&& groups[i].groupquality >= lessergroupqual
		&& groups[i].forwarddircounter>0
		&& groups[i].complementdircounter>0))) {
	avgqual+=groups[i].groupquality;
	IUPACbasegroups.push_back(groups[i].base);
	groupstaken++;

	if(groups[i].groupquality >= goodgroupqual
	   && groups[i].forwarddircounter>0
	   && groups[i].complementdircounter>0) {
	  verygoodgroups++;
	}
      }
    }
  } while((*CON_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getContigParams().con_force_nonIUPACconsensus_perseqtype
	  && groupstaken==0 && goodgroupqual>0);

  CEBUG("groups taken: " << groupstaken << "\t");
#ifdef CEBUGFLAG
  {
    for(uint32 i=0; i< IUPACbasegroups.size(); i++){
      CEBUG(IUPACbasegroups[i]);
    }
  }
#endif
  CEBUG(endl);
  if(groupstaken == 0) {
    // this can still happen if we do not force a non-IUPAC
    // or if it's a column entirely made of N (or X) (or both)

    // well, calculate the base from all columnbases
    thisbase=' ';
    for(auto & ge : groups){
      if(ge.urdids.size() && ge.base != '*') thisbase=dptools::calcIUPACConsensus(thisbase,ge.base);
    }

    avgqual=0;
    int32 n=0;
    uint32 totalids=0;
    for(uint32 i=0; i<groups.size()-1; i++){
      if(groups[i].urdids.size()) {
	totalids+=groups[i].urdids.size();
	avgqual+=groups[i].groupquality;
	n++;
      }
    }
    // test if there are more stars than other bases
    if(groups[groups.size()-1].urdids.size() > totalids) {
      // oh well, more stars than bases, make it a star
      thisbase='*';
      thisqual=groups[groups.size()-1].groupquality;
    } else if(n) {
      thisqual=avgqual/n;
    } else {
      thisqual=0;
    }
  } else if(groupstaken == 1) {
    // only one group taken (perfect)
    thisbase=IUPACbasegroups.front();
    int32 bahqual=maxqual;
    for(uint32 i=0; i<groups.size(); i++){
      if(groups[i].base!=thisbase) {
	bahqual-=groups[i].groupquality;
      }
    }
    if(bahqual<0) bahqual=0;
    thisqual=bahqual;
  } else {
    // ouch, more than one group was taken as good

    bool needIUPACconsensus=true;

    CEBUG("cpos: " << actcontigpos << "\tpossible IUPAC\tvgg: " << verygoodgroups << endl);
    //for(uint32 i=0; i< groups.size(); i++){
    //  cout << groups[i] << endl;
    //}

    // First, check whether we have a very good group
    //  as searched for above
    if(verygoodgroups==1){
      // if yes, rebuild IUPACbasegroups with very good groups
      // then, if only one group remains, perfect. Else it'll be
      //  a IUPAC

      CEBUG("recheck group: ");

      std::vector<char> newIUPACbasegroups;
      int32 newavgqual=0;
      int32 newgroupstaken=0;
      for(uint32 i=0; i<groups.size(); i++){
	if(groups[i].urdids.size()
	   && (groups[i].groupquality >= goodgroupqual
	       && groups[i].forwarddircounter>0
	       && groups[i].complementdircounter>0)){
	  newavgqual+=groups[i].groupquality;
	  newIUPACbasegroups.push_back(groups[i].base);
	  newgroupstaken++;
	}
      }

      if(newgroupstaken == 1) {
	CEBUG("success!");
	needIUPACconsensus=false;

	IUPACbasegroups=newIUPACbasegroups;
	avgqual=newavgqual;
	groupstaken=newgroupstaken;

	// only one group taken (perfect)
	thisbase=IUPACbasegroups.front();
	int32 bahqual=maxqual;
	for(uint32 i=0; i<groups.size(); i++){
	  if(groups[i].base!=thisbase) {
	    bahqual-=groups[i].groupquality;
	  }
	}
	if(bahqual<0) bahqual=0;
	thisqual=bahqual;
      }
    }
    CEBUG(endl);

    if(needIUPACconsensus){
      if((*CON_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getContigParams().con_force_nonIUPACconsensus_perseqtype){
	// well, can't get a good decision, but user wants a non-IUPAC
	// make a majority-vote shootout among the chosen
	//   groups
	makeIntelligentConsensus_helper3(thisbase,
					 thisqual,
					 groups,
					 IUPACbasegroups);
      }else{
	// we'll have to make
	//  a IUPAC consensus here. Alas, if a star (*) is part of this group,
	//  it will 'go under' (the good bases always win). So check if the star
	//  has max quality.
	bool isstar=false;
	if(groups[groups.size()-1].groupquality == maxqual) {
	  // yep, now check if one of the bases has the same qual
	  isstar=true;
	  for(uint32 i=0; i<groups.size()-1; i++){
	    if(groups[i].groupquality == maxqual) {
	      // in dubio pro reo: it could be a base
	      isstar=false;
	    }
	  }
	}
	CEBUG("isstar: " << isstar << endl);
	if(isstar) {
	  thisbase='*';
	  thisqual=groups[groups.size()-1].groupquality;
	} else {
	  // really several different good bases.
	  thisqual=avgqual/groupstaken;
	  thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
	}
      }
    }
  }
}



/*************************************************************************
 *
 * For a given readtype, return thisbase and thisqual as result for this
 *  position in the contig
 * Also return the potential number of solutions (bases), the one which
 *  was chosen plus the ones which were not
 *
 * all other arguments are passovers from the main function (saving STL
 *  setup and memory allocation time)
 *
 * Side effects: manifold. Most notably: all the passovers from the
 *  main function which are not const will be overwritten with values
 *  reflecting the calculation of base and quality of this read type
 *
 * E.g.: the bases that were considered to make the consensus are in
 *  IUPACbasegroups, the base-groups in groups etc.
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}

void Contig::makeIntelligentConsensus_helper1(char & thisbase, base_quality_t & thisqual, std::vector<char> & IUPACbasegroups, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const int32 mincoverage, std::vector<nngroups_t> & groups, std::vector<nngroups_t> & maskedshadowgroups, const std::vector<rpicsocache_t> & read_pcrIs_in_col, std::vector<int8> & maskshadow, uint8 actreadtype, int32 strainidtotake, const bbstrainmask_t strainmask, char missingcoveragechar)
{
  FUNCSTART("void Contig::makeIntelligentConsensus_helper1(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const int32 mincoverage, std::vector<nngroups_t> & groups, std::vector<nngroups_t> & maskedshadowgroups, std::vector<char> & IUPACbasegroups, const std::vector<int32> & read_pcrIs_in_col, std::vector<int8> & maskshadow, uint8 actreadtype, char missingcoveragechar)");

#if 0
  if(read_pcrIs_in_col.empty()){
    thisbase='!';
    FUNCEND();
    return;
  }
  if(CON_counts[actcontigpos].total_cov < mincoverage){
    thisbase='N';
    FUNCEND();
    return;
  }

#else

  size_t nummapped=0;
  if(actreadtype == ReadGroupLib::SEQTYPE_SOLEXA && CON_isbackbonecontig){
    // Ok, check whether the merged bases belong to this strain
    //  if not, well then no merged bases exist
    if(ccI->bbcountsf[0] + ccI->bbcountsr[0]){
      if(ccI->bbstrains[0] & strainmask){
	nummapped = ccI->bbcountsf[0] + ccI->bbcountsr[0];
      }
    }
  }
  if(read_pcrIs_in_col.size()+nummapped == 0
     || read_pcrIs_in_col.size()+nummapped < mincoverage){

    // maybe we can get out quick if there's no backbone read present
    // so check for backbone ...
    bool hasbb=false;
    for(auto & rpice : read_pcrIs_in_col){
      if (rpice.pcrI->isBackbone()) {
	hasbb=true;
	break;
      }
    }
    // ... and get out if there is none
    if (!hasbb){
      thisbase='N';
      thisqual=0;
      FUNCEND();
      return;
    }
  }
#endif

  const contig_parameters & con_rt_params = (*CON_miraparams)[actreadtype].getContigParams();

  CEBUG("conpos: " << actcontigpos << "\tnumreads: " << read_pcrIs_in_col.size() << endl);

  /*
     slowwwwwwwwwwwwwwwwwwwwwww

     on a 170MB CAF file, caf2fasta goes from 1:46 to 1:40 when
     replacing these two lines with the loop below (from 49s for
     output down to 43s)

     groups=emptygroups;
     maskedshadowgroups=emptygroups;

  */

  bool xonly=true;
  bool hasmaskedset=false;

  if(nummapped>0 && ccI->getOriginalBBChar() != 'X') xonly=false;

  uint32 numcerreads=0;
  uint32 cergroupindex=0;

  // taking these out of the for loop apparently makes it slightly faster
  char           base;
  base_quality_t qual;
  int32 realreadpos;
  base_quality_t maxqual;
  uint32 maxcount;

  for(auto & rpice : read_pcrIs_in_col){

    //if(ric.read.isShortReadMapping()) continue;

    // yup, looks like this is better than using "pcrI->" all the time
    //  it is at least faster
    // maybe too hard for GCC to optimise across whole function

#ifdef CLOCK_STEPS_CONSSUB
    timersub.reset();
#endif

    auto & actread=*(rpice.pcrI);

    // this is not too slow when called hundreds of times per read
    //     int32 readpos=rpice.pcrI.contigPos2UnclippedReadPos(actcontigpos);
    // as trial: replicate contigPos2UnclippedReadPos() here, but using
    //  the cached value of the read start position
    if (rpice.dir > 0) {
      realreadpos = actcontigpos - rpice.readstartpos + rpice.leftclip;
      base=toupper(actread.nocheckGetBaseInSequence(realreadpos));
      qual=actread.nocheckGetQualityInSequence(realreadpos);
    } else {
      int32 readpos = rpice.lenseq - rpice.rightclip + actcontigpos - rpice.readstartpos;
      base=toupper(actread.nocheckGetBaseInComplementSequence(readpos));

      //realreadpos=actread.calcComplPos(readpos);
      realreadpos=rpice.lenseq-1-readpos;
      qual=actread.nocheckGetQualityInSequence(realreadpos);

    }

#ifdef CLOCK_STEPS_CONSSUB
    CON_us_steps_cons[USCLOCONS_H1_PL_2]+=timersub.diffAndReset();
#endif

    // this one gives a problem
    //   if(base != 'X' && base != '*') xonly=false;
    // when one retrieves the consensus of a single strain in multiple
    //  strain assemblies: it may well be only a star!
    // preliminary fix: back to checking of X only

    if(base != 'X') xonly=false;

    CEBUG("\t" << base << "\t" << (uint16) qual << endl);

    if(rpice.iscer){
      ++numcerreads;
    }

    // bases near start/end of reads might be dangerous because
    //  of possible vector leftovers or still bad quality
    // so, if there is more than one read in this column,
    // ...
    // - if within 10 bases of vector clip, max quality = 2*distance to svec
    // - if not, lower the quality values of the bases that are near ends
    //   (con_endreadmarkexclusionarea bases) of (eventually clipped) reads
    // must be regular read (i.e. not backbone, rail or shortreadmapping

    if(likely(!rpice.israil
	      && !rpice.isbackbone
	      && !rpice.iscer)){
      if(likely(read_pcrIs_in_col.size() >1)) {
	if((actread.getLSClipoff() >0
	    && realreadpos < actread.getLSClipoff()+10)
	   || (actread.getRSClipoff() < rpice.lenseq
	       && realreadpos > actread.getRSClipoff()-10)){
	  CEBUG(actread.getName()<< ": near seq vec in read, lowering the quality.\n");
	  int32 distance;
	  if(realreadpos < actread.getLSClipoff()+10){
	    distance=realreadpos-actread.getLSClipoff();
	  }else{
	    distance=actread.getRSClipoff()-realreadpos;
	  }
	  if(qual > distance*2){
	    qual=distance*2;
	  }
	} else if(realreadpos <= rpice.leftclip+con_rt_params.con_endreadmarkexclusionarea
		  || realreadpos >= rpice.rightclip-con_rt_params.con_endreadmarkexclusionarea) {
	  CEBUG(actread.getName()<< ": near end of read, lowering the quality.\n");

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

	  if(rpice.seqtype==ReadGroupLib::SEQTYPE_SOLEXA){
	    // decrease by one
	    if(likely(qual)) --qual;
	  }else{
	    // all other atm
	    if(qual>=10) {
	      qual-=10;
	    } else {
	      qual=5;
	    }
	    if(unlikely(qual<5)) qual=5;
	  }
	}
      }
    }

#ifdef CLOCK_STEPS_CONSSUB
    CON_us_steps_cons[USCLOCONS_H1_PL_3]+=timersub.diffAndReset();
#endif


    // No railreads at all, have been taken out earlier!
    //// Quality of bases from railreads are set to 0 so as not
    ////  to be counted twice (as bases are the same as in backbone)
    //if(rpice.israil) qual=0;



    bool maskedset=false;


    // TODO: rework the if(...hasTag(Read::REA_tagFpas,realreadpos))
    //  to check for all mask taks in a masktagstring vector

//    if(maskshadow[actcontigpos]) {
//      // remember that the readpos computing routine doesn't take care
//      //  of direction, so we have to complement that position in reverse cases
//      int32 realreadpos=readpos;
//      if(rpice.dir<0){
//	realreadpos=actread.calcComplPos(readpos);
//      }
//
//      CEBUG("MASKED: " << actcontigpos << "\t");
//      CEBUG(actread.getName() << "\t" << realreadpos << "\t");
//      if(actread.hasTag(Read::REA_tagFpas,realreadpos)) {
//	CEBUG("in" << endl);
//	for(uint32 bindex=0; bindex<maskedshadowgroups.size(); bindex++) {
//	  if(maskedshadowgroups[bindex].base==base) {
//	    maskedshadowgroups[bindex].urdids.push_back(actreadid);
//	    maskedshadowgroups[bindex].quals.push_back(qual);
//	    maskedshadowgroups[bindex].directions.push_back(rpice.dir);
//	    if(rpice.dir>0){
//	      maskedshadowgroups[bindex].hasforwarddir=true;
//	    }else{
//	      maskedshadowgroups[bindex].hascomplementdir=true;
//	    }
//	    // special case: treat short read mapping as both forward and reverse
//	    if(rpice.iscer){
//	      maskedshadowgroups[bindex].hasforwarddir=true;
//	      maskedshadowgroups[bindex].hascomplementdir=true;
//	    }
//	    maskedset=true;
//	    hasmaskedset=true;
//	    break;
//	  }
//	}
//      }
//    }


    if(likely(maskedset==false)) {
      CEBUG("out" << endl);
      //for(uint32 bindex=0; bindex<groups.size(); bindex++) {
      uint32 tmpgi=0;
      for(auto & groupe : groups){
	if(groupe.base==base) {
	  groupe.urdids.push_back(rpice.urdid);
	  groupe.quals.push_back(qual);
	  groupe.directions.push_back(rpice.dir);
	  if(rpice.dir>0){
	    groupe.forwarddircounter++;
	  }else{
	    groupe.complementdircounter++;
	  }

	  // special case: treat short read mapping as both forward and reverse
	  if(rpice.iscer){
	    // we already counted forward
	    //groupe.forwarddircounter++;
	    groupe.complementdircounter++;
	    cergroupindex=tmpgi;
	  }
	  break;
	}
	++tmpgi;
      }

#ifdef CLOCK_STEPS_CONSSUB
      CON_us_steps_cons[USCLOCONS_H1_PL_4]+=timersub.diff();
#endif
    }
  }

  // basically, the loop above counted bases for CER reads twice: fwd and reverse
  // now correct the counters to reflect true coverage
  // if number of CER reads is uneven, forward direction gets one more than complement direction
  if(numcerreads>0){
    groups[cergroupindex].complementdircounter-=numcerreads/2;
    groups[cergroupindex].forwarddircounter-=numcerreads-numcerreads/2;
  }

#ifdef CLOCK_STEPS_CONS
  CON_us_steps_cons[USCLOCONS_H1_countPCRI]+=read_pcrIs_in_col.size();
  CON_us_steps_cons[USCLOCONS_H1_PCRI]+=timeract.diffAndReset();
#endif

  if(hasmaskedset){
    // check whether there is anything in groups
    //  and maskedshadowgroups selected
    // if groups empty and masked not, copy masked to groups
    bool groupsempty=true;
    bool maskedempty=true;
    for(uint32 bindex=0; bindex<groups.size(); ++bindex) {
      if( groups[bindex].urdids.size()) groupsempty=false;
      if( maskedshadowgroups[bindex].urdids.size()) maskedempty=false;
    }
    if(groupsempty && maskedempty==false) {
      CEBUG("EMPTY: " << actcontigpos << endl);
      groups=maskedshadowgroups;
    }
#ifdef CLOCK_STEPS_CONS
    CON_us_steps_cons[USCLOCONS_H1_EGROUP]+=timeract.diffAndReset();
#endif
  }

  if(xonly) {
    thisbase='X';
  } else {
    maxqual=0;
    maxcount=0;
    for (auto & groupe : groups){
      calcGroupQual(groupe);
      if(groupe.groupquality>maxqual){
	maxqual=groupe.groupquality;
      }
      if(groupe.urdids.size()>maxcount){
	maxcount=groupe.urdids.size();
      }
      //if(actcontigpos>830 && actcontigpos <920) {
      CEBUG(actcontigpos << "\tb: " << groupe.base << "\tgq: " << (uint16) groupe.groupquality << "\ts: " << groupe.urdids.size() << endl);
      //}
    }

#ifdef CLOCK_STEPS_CONS
    CON_us_steps_cons[USCLOCONS_H1_GQUAL]+=timeract.diffAndReset();
#endif

    //if(actcontigpos>830 && actcontigpos <920) {
    CEBUG(actcontigpos << "maxqual: " << (uint16) maxqual << "\t" << "maxcount: " << maxcount << endl);
    //}

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

    // replacing switch() with if-else cascade showed no improvement
    switch(actreadtype) {
    case ReadGroupLib::SEQTYPE_SANGER : {
      makeIntelligentConsensus_helper2_calcSangerQual(thisbase,
						      thisqual,
						      actcontigpos,
						      groups,
						      IUPACbasegroups,
						      maxqual,
						      maxcount);
      break;
    }
    case ReadGroupLib::SEQTYPE_454GS20 : {
      makeIntelligentConsensus_helper2_calc454GS20(thisbase,
						   thisqual,
						   actcontigpos,
						   groups,
						   IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_IONTORRENT : {
      makeIntelligentConsensus_helper2_calcIonTorrent(thisbase,
						      thisqual,
						      actcontigpos,
						      groups,
						      IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_PACBIOHQ :
    case ReadGroupLib::SEQTYPE_PACBIOLQ : {
      // TODO: eventually two different for hq and lq ?
      makeIntelligentConsensus_helper2_calcPACBIOHQ(thisbase,
						    thisqual,
						    actcontigpos,
						    groups,
						    IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_TEXT : {
      makeIntelligentConsensus_helper2_calcText(thisbase,
						thisqual,
						actcontigpos,
						groups,
						IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_SOLEXA : {
      makeIntelligentConsensus_helper2_calcSOLEXA(
	thisbase,
	thisqual,
	actcontigpos,
	ccI,
	groups,
	IUPACbasegroups,
	strainmask
	);
      break;
    }
    case ReadGroupLib::SEQTYPE_ABISOLID : {
      cout << "Actreadtype: " << static_cast<uint16>(actreadtype) << endl;
      MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support.");
      break;
    }
    default: {
      cout << "Actreadtype: " << static_cast<uint16>(actreadtype) << endl;
      MIRANOTIFY(Notify::INTERNAL, "Unknown read type.");
    }
    }

#ifdef CLOCK_STEPS_CONS
    CON_us_steps_cons[USCLOCONS_H1_CALLH2]+=timeract.diffAndReset();
#endif

    if(!dptools::isValidIUPACStarBase(thisbase)){
      CEBUG("ALERT! This is not a valid base: '" << thisbase << "'\t: " << static_cast<uint16>(thisbase) << '\n');
    }

  }

#ifdef CLOCK_STEPS_CONS
  CON_us_steps_cons[USCLOCONS_H1_TOTAL]+=timerh1total.diff();
#endif


  FUNCEND();
}





///*************************************************************************
// *
// *
// *
// *************************************************************************/
//
//void Contig::priv_mic_helper_chooseDiploidBase(char & thisbase, base_quality_t & thisqual,const std::vector<uint8> & allstpossiblebases, const std::vector<std::vector<nngroups_t>> & groups, uint8 numpossiblebases)
//{
//  FUNCSTART("void Contig::priv_mic_helper_chooseDiploidBase(char & thisbase, base_quality_t & thisqual,const std::vector<uint8> & allstpossiblebases, const std::vector<std::vector<nngroups_t>> & groups, uint8 numpossiblebases)");
//
//  if(numpossiblebases==2){
//    if(allstpossiblebases['A']) thisbase='A';
//    thisbase=dptools::calcIUPACConsensus(thisbase,allstpossiblebases['C']);
//    thisbase=dptools::calcIUPACConsensus(thisbase,allstpossiblebases['G']);
//    thisbase=dptools::calcIUPACConsensus(thisbase,allstpossiblebases['T']);
//    thisqual=22; // TODO: continue here
//  }else{
//    // now what???
//    BUGIFTHROW(true,"need to cont here");
//  }
//  if(allstpossiblebases['*']) {
//    thisbase=tolower(thisbase);
//  }
//}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::priv_mic_helper_chooseOneBase(char & thisbase, base_quality_t & thisqual,const std::vector<uint32> & allstpossiblebases, const std::vector<std::vector<nngroups_t>> & groups)
{
  if (!dptools::isValidACGTStarBase(thisbase)) {
    if(allstpossiblebases['A']){
      thisbase='A';
    }else if(allstpossiblebases['C']){
      thisbase='C';
    }else if(allstpossiblebases['G']){
      thisbase='G';
    }else if(allstpossiblebases['T']){
      thisbase='T';
    }else{
      thisbase='*';
    }
  }
  uint32 totalqual=0;
  for(auto & stgroup : groups){
    for(auto & groupe : stgroup){
      if(groupe.base==thisbase && groupe.urdids.size()) {
	totalqual+=groupe.groupquality;
	thisqual=std::max(thisqual,groupe.groupquality);
	CEBUG("coB: " << groupe.base << " " << thisbase << " " << static_cast<uint16>(thisqual) << " " << totalqual << "\n");
      }
    }
  }
  // quality = 75% of sum of group quals ... or maxqual if higher
  totalqual=totalqual*75/100;
  thisqual=static_cast<uint8>(std::max(static_cast<uint32>(thisqual),totalqual));
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}
void Contig::priv_mic_helper_chooseMultiBase(char & thisbase, base_quality_t & thisqual,const std::vector<uint32> & allstpossiblebases, const std::vector<std::vector<nngroups_t>> & groups)
{
  FUNCSTART("void Contig::priv_mic_helper_chooseMultiBase(char & thisbase, base_quality_t & thisqual,const std::vector<uint32> & allstpossiblebases, const std::vector<std::vector<nngroups_t>> & groups)");

  thisbase='!';
  uint8 seenbases=0;
  if(allstpossiblebases['A']) {
    thisbase='A';
    seenbases=1;
  }
  if(allstpossiblebases['C']) {
    thisbase=dptools::calcIUPACConsensus(thisbase,'C');
    ++seenbases;
  }
  if(allstpossiblebases['G']) {
    thisbase=dptools::calcIUPACConsensus(thisbase,'G');
    ++seenbases;
  }
  if(allstpossiblebases['T']) {
    thisbase=dptools::calcIUPACConsensus(thisbase,'T');
    ++seenbases;
  }
  if(allstpossiblebases['*']) {
    ++seenbases;
  }

  if(seenbases==0){
    thisbase='N';
    thisqual=0;
  }else{
    uint32 totalqual=0;
    seenbases=0;
    for(auto & stgroup : groups){
      for(auto & groupe : stgroup){
	if(groupe.urdids.size()) {
	  CEBUG("abc: " << groupe.base << " " << thisbase << " -> " << dptools::areBasesContained(groupe.base,thisbase) << endl);
	  if(dptools::areBasesContainedWithN(groupe.base,thisbase) || groupe.base=='*'){
	    totalqual+=groupe.groupquality;
	    ++seenbases;
	  }
	}
      }
    }
    BUGIFTHROW(seenbases==0,"seenbases==0??? No group set? Strange.");
    thisqual=totalqual/seenbases;
  }
}
//#define CEBUG(bla)

/*************************************************************************
 *
 * Calculate the 'true' consensus and gives back a string with consensus
 *  and a vector with the base quality of each base
 *
 * strainidtotake: <0 means "all", >=0 means "exactly those reads with that id"
 *
 * If the ostream parameter is != nullptr, also writes a .tcs live file
 *  to it
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}

void Contig::makeIntelligentConsensus(std::string & target, std::vector<base_quality_t> & qual, std::vector<int32> * targetadjustments, int32 from, int32 to, int32 strainidtotake, int32 mincoverage, base_quality_t minqual, char missingcoveragechar, bool assumediploid, bool allowiupac, bool addconstag)//, ostream * ostr, bool contagsintcs)
{
  FUNCSTART("void Contig::makeIntelligentConsensus(std::string & target, std::vector<base_quality_t> & qual, int32 from, int32 to, int32 mincoverage, base_quality_t minqual, int32 strainidtotake)");//, ostream * ostr, bool contagsintcs)");

  static const char acgtgapbases[]="ACGT*";

  VCOUT("makeIntelligentConsensus() complete calc .. "; cout.flush());

  //CON_cebugflag=true;

  CEBUG("MIC\n");
  CEBUG("from " << from << endl);
  CEBUG("to " << to << endl);
  CEBUG("mincoverage " << mincoverage << endl);
  CEBUG("minqual " << static_cast<uint16>(minqual) << endl);
  CEBUG("missingcovchar " << missingcoveragechar << endl);
  CEBUG("strainidtotake " << strainidtotake << endl);
  CEBUG("addconstag " << addconstag << endl);


  BUGIFTHROW(from>to,"from>to?");
  BUGIFTHROW(from<0, "from < 0 ?");

  suseconds_t mict_fin=0;
  suseconds_t mict_pre=0;
  suseconds_t mict_shadow=0;
  suseconds_t mict_fallout=0;
  suseconds_t mict_newin=0;
  suseconds_t mict_helper1=0;
  suseconds_t mict_restofloop=0;
  suseconds_t mict_totalloop=0;

  timeval us_start;
  gettimeofday(&us_start,nullptr);

  bool strainisreference=false;
  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(rgid.getStrainID()==strainidtotake
       && rgid.isBackbone()) strainisreference=true;
  }
  CEBUG("strainisreference " << strainisreference << endl);


  // Finalising the contig initialises the output order structure vector
  finalise();
  mict_fin+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  //const contig_parameters & con_params = CON_miraparams->getContigParams();

  if( to > static_cast<int32>(CON_counts.size())) to=CON_counts.size();
  int32 len_target=to-from;

  //target.resize(len_target);
  target.clear();
  target.reserve(len_target+10);
  qual.clear();
  qual.reserve(len_target+10);

  // for calculating the adjustments, only do this when whole
  //  consensus is calculated
  if(targetadjustments != nullptr && from==0 && to==CON_counts.size()) {
    targetadjustments->clear();
    targetadjustments->reserve(len_target+10);
  }
  int32 unpaddedposcounter=0;


  nngroups_t emptygroup;
  emptygroup.base='!';
  emptygroup.valid=false;
  // TODO: make also use of bothdirectionspresent in this routine?
  emptygroup.forwarddircounter=0;
  emptygroup.complementdircounter=0;
  emptygroup.groupquality=0;
  std::vector<nngroups_t> emptygroups;
  for(uint32 i=0; i<5; i++) {
    emptygroups.push_back(emptygroup);
    emptygroups[i].base= acgtgapbases[i];
  }

/*
  // if there's a stream, we're dumping TCS
  // initialise a quick lookup vector to point at the positions in the
  //  consensus that have a tag
  // a for any tag
  // d for dangerous tag
  std::vector<bool> tcs_aconstagpositions;
  std::vector<bool> tcs_dconstagpositions;
  if(ostr != nullptr && contagsintcs){
    tcs_aconstagpositions.resize(CON_counts.size(),false);
    tcs_dconstagpositions.resize(CON_counts.size(),false);
    auto I=CON_consensus_tags.cbegin();
    for(; I!=CON_consensus_tags.end(); I++) {
      for(uint32 i=I->from; i<=I->to; i++) tcs_aconstagpositions[i]=true;
      if(I->identifier == CON_tagentry_idSRMc
	 || I->identifier == CON_tagentry_idWRMc
	 || I->identifier == CON_tagentry_idDGPc
	 || I->identifier == CON_tagentry_idUNSc
	 || I->identifier == CON_tagentry_idIUPc){
	for(uint32 i=I->from; i<=I->to; i++) tcs_dconstagpositions[i]=true;
      }
    }
  }
  // temporary vectors for TCS
  std::vector<int32> tcs_totalgroupcount;
  tcs_totalgroupcount.reserve(emptygroups.size());
  std::vector<int32> tcs_totalgroupqual;
  tcs_totalgroupqual.reserve(emptygroups.size());
*/


  std::vector<char> seqtypepicks;

  std::vector<nngroups_t> maskedshadowgroups=emptygroups;

  // for each read type, we will predict a base and quality and such
  std::vector<char> pred_baseperst(ReadGroupLib::getNumSequencingTypes());
  std::vector<base_quality_t> pred_qualperst(ReadGroupLib::getNumSequencingTypes());

  // as well as the possible alternative bases
  std::vector<std::vector<char>> possiblebases(ReadGroupLib::getNumSequencingTypes());
  std::vector<uint32> allstpossiblebases(256,0); // used flag wheter ACGT* was seen and count occurrences
  std::vector<uint32> allstpossiblebasesf(256,0); // same, forward only
  std::vector<uint32> allstpossiblebasesr(256,0); // same, reverse only

  std::vector<uint32> countsofpossiblebases; // just counts, for sorting later


  // as well as a goodness level estimate of the prediction
  std::vector<std::vector<int8>> predictlevels(ReadGroupLib::getNumSequencingTypes());

  // for each read type a std::vector<nngroups_t> to hold all estimates
  //  regarding bases and qualities
  std::vector<std::vector<nngroups_t>> groupsvec;
  groupsvec.reserve(ReadGroupLib::getNumSequencingTypes());

  // a vector (read_ids_in_col) keeps the ids of the reads that are
  //  covering a specific position of the contig
  // reserving 1000 position should be enough 99.9% of all cases,
  //  is automatically extended by STL if needed.

  // now with different read types (sanger, 454 etc), we need a vector
  //  of iterators to the PCR. Was formerly a vector of read-ids, but
  //  not possible anymore

  // and now we cache the read start offsets, thus using a struct

  std::vector<std::vector<rpicsocache_t> > read_pcrIs_in_col;
  read_pcrIs_in_col.resize(ReadGroupLib::getNumSequencingTypes());

  // ok, fill in some
  for(auto & rpice : read_pcrIs_in_col){
    rpice.reserve(1000);
    groupsvec.push_back(emptygroups);
  }

  // fill the vector in case we are starting within the contig
  if(from>0) {
    // TODO: this is wrong when going only after specific strainids!
    // correct that ASAN (as soon as needed)

    MIRANOTIFY(Notify::INTERNAL, "starting with from > 0 not available yet.");

    //getReadIDsAtContigPosition(read_pcrIs_in_col,from);
  }


  mict_pre=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  std::vector<int8> maskshadow;
  std::vector<multitag_t::mte_id_t> masktagstrings;

  // Bach: 17.08.2008
  // the new strategy of tagging poly-AT sites and keeping them in
  //  the read (no clipping) makes it necessary to keep sequence
  //  under Fpas tags as full valid member of consensus somputation
  // Therefore, Fpas may NOT be put into the masktagstrings anymore!
  // masktagstrings.push_back(Read::REA_tagFpas);

  if(!buildMaskShadow(maskshadow,masktagstrings,false)) maskshadow.clear();

  mict_shadow=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  std::vector<uint8> cached_contigseqtypes(ReadGroupLib::getNumSequencingTypes(),0);
  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); ++actseqtype){
    cached_contigseqtypes[actseqtype]=hasSeqTypeData(actseqtype);
  }

  auto pcrI=CON_reads.begin();
  auto ccI=CON_counts.cbegin();
  std::advance(ccI,from);

  // this is the loop that updates the vector that
  //  keeps track of the reads that are
  //  covering a specific position of the contig
  // works quite simple: reads in the vector that are ending get
  //  thrown out, new reads starting are entered.

  CEBUG("CON_counts.size(): " << CON_counts.size() << endl);

  bbstrainmask_t strainmask=-1;
  if(CON_isbackbonecontig && strainidtotake>=0) strainmask=getBBStrainMask(strainidtotake);

  timeval us_loop;

  for(uint32 actcontigpos=from; actcontigpos<to; ++ccI, ++actcontigpos){
    gettimeofday(&us_loop,nullptr);
    // updating the pcrIs of the reads at that position
    CEBUG("cc acp: " << actcontigpos << endl);
    // first delete those who fall out at this new position
    for(uint32 actseqtype=0; actseqtype<cached_contigseqtypes.size(); ++actseqtype){
      if(cached_contigseqtypes[actseqtype]){
	CEBUG("cc art " << actseqtype << " rpic " << read_pcrIs_in_col[actseqtype].size() << endl);
	auto Ifrom=read_pcrIs_in_col[actseqtype].begin();
	auto Ito=Ifrom;
	// pcrIs are expensive on copy: by not blindly copying every iterator but only when needed,
	//  we save between 30% and 40% time (maybe 5% total in a total cons. calc, but still)
	bool needcopy=false;
	for(;Ifrom != read_pcrIs_in_col[actseqtype].end(); ++Ifrom){
	  if(needcopy) *Ito=*Ifrom;
	  if(Ifrom->readstartpos+Ifrom->lenclippedseq > actcontigpos) {
	    ++Ito;
	  }else{
	    needcopy=true;
	    CEBUG("cc thrown out " << *Ito << endl);
	  }
	}
	//if(Ito != Ifrom) {
	//  read_pcrIs_in_col[actseqtype].resize(Ito-read_pcrIs_in_col[actseqtype].begin());
	//}
	// PlacedContigReads::const_iterator has no default constructor,
	//  i.e., resize cannot be used because it wants a default constructor
	//  though we always reduce the vector here
	// but we can use erase()
	CEBUG("loopend" << endl);
	if(needcopy){
	  CEBUG("cc popping " << Ito - Ifrom << " elements\n");
	  read_pcrIs_in_col[actseqtype].erase(Ito,Ifrom);
	}
      }
    }
    mict_fallout+=diffsuseconds(us_loop);
    gettimeofday(&us_loop,nullptr);

    // now insert ids of reads that have newly started at this position
    // Don't take railreads, backbone is there.
    for(;pcrI != CON_reads.end() && pcrI.getReadStartOffset() == actcontigpos; ++pcrI){
      if( ! pcrI->isRail()) {

	// new! one can also make a choice which strains to take into account
	// either all (strainidtotake<0 or exactly the one chosen)
	if(strainidtotake < 0
	   || pcrI->getStrainID() == strainidtotake) {
	  read_pcrIs_in_col[pcrI->getSequencingType()].emplace_back(rpicsocache_t(pcrI,
										  pcrI.getReadStartOffset(),
										  pcrI.getURDID(),
										  pcrI.getReadDirection(),
										  pcrI->getLeftClipoff(),
										  pcrI->getRightClipoff(),
										  pcrI->getLenSeq(),
										  pcrI->getLenClippedSeq(),
										  pcrI->getSequencingType(),
										  pcrI->isRail(),
										  pcrI->isBackbone(),
										  pcrI->isCoverageEquivalentRead()));
	  CEBUG("cc taken " << pcrI << endl);
	}
      }
    }
    mict_newin+=diffsuseconds(us_loop);
    gettimeofday(&us_loop,nullptr);

    // for each read type, predict a base and a quality,
    //  but also get back the potential solutions and number thereof

    //std::vector<uint32> numsolutionsvec(ReadGroupLib::getNumSequencingTypes(),0);

    allstpossiblebases['A']=0;
    allstpossiblebases['C']=0;
    allstpossiblebases['G']=0;
    allstpossiblebases['T']=0;
    allstpossiblebases['*']=0;
    allstpossiblebasesf['A']=0;
    allstpossiblebasesf['C']=0;
    allstpossiblebasesf['G']=0;
    allstpossiblebasesf['T']=0;
    allstpossiblebasesf['*']=0;
    allstpossiblebasesr['A']=0;
    allstpossiblebasesr['C']=0;
    allstpossiblebasesr['G']=0;
    allstpossiblebasesr['T']=0;
    allstpossiblebasesr['*']=0;
    uint32 numreadtypeswithsolution=0;
    for(uint32 actseqtype=0; actseqtype<cached_contigseqtypes.size(); ++actseqtype){
      if(pred_baseperst[actseqtype]!='!') {
#ifdef CLOCK_STEPS_CONS
	MIRATimer timerh1total;
	MIRATimer timeract;
#endif
#ifdef CLOCK_STEPS_CONSSUB
	MIRATimer timersub;
#endif
	for(auto & ge : groupsvec[actseqtype]){
	  ge.reset();
	}
	if(!maskshadow.empty()){
	  for(auto & ge : maskedshadowgroups){
	    ge.reset();
	  }
	}
#ifdef CLOCK_STEPS_CONS
	CON_us_steps_cons[USCLOCONS_H1_MGROUPS]+=timeract.diffAndReset();
#endif
      }
      pred_baseperst[actseqtype]='!';
      pred_qualperst[actseqtype]=0;
      possiblebases[actseqtype].clear();
      if(!cached_contigseqtypes[actseqtype]) continue;
      if(unlikely(read_pcrIs_in_col[actseqtype].empty())) continue;

      makeIntelligentConsensus_helper1(pred_baseperst[actseqtype],
				       pred_qualperst[actseqtype],
				       possiblebases[actseqtype],
				       actcontigpos,
				       ccI,
				       mincoverage,
				       groupsvec[actseqtype],
				       maskedshadowgroups,
				       read_pcrIs_in_col[actseqtype],
				       maskshadow,
				       actseqtype,
				       strainidtotake,
				       strainmask,
				       missingcoveragechar
	);

      // Count num occurrences for ACGT*
      if(!possiblebases[actseqtype].empty()){
	++numreadtypeswithsolution;
	for(auto & pbe : possiblebases[actseqtype]){
	  for(auto & groupe : groupsvec[actseqtype]){
	    if(groupe.base==pbe){
	      allstpossiblebasesf[pbe]+=groupe.forwarddircounter;
	      allstpossiblebasesr[pbe]+=groupe.complementdircounter;
	      allstpossiblebases[pbe]+=groupe.forwarddircounter+groupe.complementdircounter;
	      if(actseqtype==ReadGroupLib::SEQTYPE_SOLEXA
		 && CON_isbackbonecontig
		 && toupper(ccI->getOriginalBBChar())==pbe){
		allstpossiblebases[pbe] += ccI->bbcountsf[0] + ccI->bbcountsr[0];
	      }
	      break; // found, we can break innner loop
	    }
	  }
	}
      }
    }

    CEBUG("PB-A: " << allstpossiblebases['A'] << " PB-C: " << allstpossiblebases['C'] << " PB-G: " << allstpossiblebases['G'] << " PB-T: " << allstpossiblebases['T'] << " PB-*: " << allstpossiblebases['*'] << endl);

#ifndef CLOCK_STEPS_CONS
    mict_helper1+=diffsuseconds(us_loop);
#endif
    gettimeofday(&us_loop,nullptr);


    /* after the previous loop, we have the following important things
       per read type (other things too, less important at the moment):
       - the picked base and quality in pred_baseperst[] and
         pred_qualperst[]
       - the same for iupac (base, qual, and whether gap is involved)
       - possible alternatives (including picked one) in possiblebases[]
       - all group calculations for all read types in groupsvec[]
       - count of possible bases in allstpossiblebases['ACGT*']
     */

    // see how many seq types made a prediction
    // what a lazy, lazy way with a vector, push_back, sort and unique ...
    seqtypepicks.clear();
    CEBUG("pred_b:        ");
    for(auto & pb : pred_baseperst) {
      CEBUG(pb);
      if(dptools::isValidACGTStarBase(pb)) {
	seqtypepicks.push_back(pb);
      }
    }
    if(seqtypepicks.size()>1){
      // use non-parallel sort, this is a tiny thing
      mstd::ssort(seqtypepicks);
      mstd::unique(seqtypepicks);
    }
    CEBUG('\n');

    // Party time: let's find out which base we take
    char thisbase='@';
    base_quality_t thisqual=0;

    uint8 numpossiblebases=(allstpossiblebases['A']>0)+(allstpossiblebases['C']>0)+(allstpossiblebases['G']>0)+(allstpossiblebases['T']>0)+(allstpossiblebases['*']>0);

    CEBUG("numpossiblebases 1: " << static_cast<uint16>(numpossiblebases) << endl);

    int8 bestpredictlevel=127;
    if(numpossiblebases>1){
      bestpredictlevel=priv_mic_helper_rateGoodnessLevelOfAllGroups(ccI,groupsvec,predictlevels,possiblebases);

      // careful: recalcs allstpossiblebases and numpossiblebases!
      // also sets predictlevels of untake groups to 127
      priv_mic_helper_filterGoodnessLevels(ccI,
					   bestpredictlevel+2,
					   predictlevels,
					   groupsvec,
					   allstpossiblebases,
					   numpossiblebases,
					   false);

      CEBUG("numpossiblebases 2: " << static_cast<uint16>(numpossiblebases) << endl);

      // uh oh .. maybe bad solexas?
      // Redo, this time with filtering of solexa by solexa group quality
      if(numpossiblebases>1){
	// careful: recalcs allstpossiblebases and numpossiblebases!
	// also sets predictlevels of untake groups to 127
	priv_mic_helper_filterGoodnessLevels(ccI,
					     bestpredictlevel+2,
					     predictlevels,
					     groupsvec,
					     allstpossiblebases,
					     numpossiblebases,
					     true);

	CEBUG("numpossiblebases 3: " << static_cast<uint16>(numpossiblebases) << endl);
      }
    }

    // if still more than one picked base and when IUPACs are allowed (but not diploid),
    //  be a bit stricter and filter at +1 level
    if(numpossiblebases>1
       && allowiupac
       && !assumediploid){
      // careful: recalcs allstpossiblebases and numpossiblebases!
      priv_mic_helper_filterGoodnessLevels(ccI,
					   bestpredictlevel+1,
					   predictlevels,
					   groupsvec,
					   allstpossiblebases,
					   numpossiblebases,
					   true);
      CEBUG("numpossiblebases 4: " << static_cast<uint16>(numpossiblebases) << endl);
    }

    CEBUG("final numpossiblebases: " << static_cast<uint16>(numpossiblebases) << endl);

    if(numpossiblebases==1){
      priv_mic_helper_chooseOneBase(thisbase, thisqual, allstpossiblebases, groupsvec);
    }else if(unlikely(numpossiblebases==0)){
      if(strainisreference){
	// if we are in the reference strain, and could not get a good base
	//   back, it's probably a IUPAC reference base we could not resolve
	// in that case, simply take the reference base we have
	thisbase=ccI->getOriginalBBChar();
	thisqual=ccI->getOriginalBBQual();
      }else{
	if(ccI->N>0){
	  // did we count any N? If yes, probably a column completely made out of N
	  thisbase='N';
	}else{
	  // if not, tough luck. probably an uncovered base.
	  thisbase='X';
	}
      }
    }else{
      if(assumediploid || allowiupac){
	priv_mic_helper_chooseMultiBase(thisbase, thisqual, allstpossiblebases, groupsvec);
      }else{
	// careful: recalcs allstpossiblebases and numpossiblebases!
	priv_mic_helper_filterGoodnessLevels(ccI,
					     bestpredictlevel,
					     predictlevels,
					     groupsvec,
					     allstpossiblebases,
					     numpossiblebases,
					     true);
	// For just two possible bases
	// if we take all strains and there is a valid acgt in the backbone, and that
	//  was also found in the valid groups, then take the backbone reference
	if( numpossiblebases == 2
	    && strainidtotake == -1
	    && dptools::isValidACGTBase(ccI->getOriginalBBChar())
	    && allstpossiblebases[toupper(ccI->getOriginalBBChar())] ) {
	  CEBUG("possible base is also in backbone, taking that\n");
	  thisbase=ccI->getOriginalBBChar();
	  numpossiblebases=1;
	}
	if(numpossiblebases>1){
	  CEBUG("seqtypepicks.size() " << seqtypepicks.size() << endl);
	  if(seqtypepicks.size()==1) {
	    // pick the one from the seqtype prediction
	    for(uint8 xi=0; xi<5; ++xi){
	      if(acgtgapbases[xi]!=seqtypepicks[0]) allstpossiblebases[acgtgapbases[xi]]=0;
	    }
	  }else{
	    // what's too much is too much: simply take the one with the most calls
	    // reverse loop to prefer * over T over G, C, A
	    uint32 tmpmax=0;
	    for(int8 xi=4; xi>=0; --xi){
	      if(allstpossiblebases[acgtgapbases[xi]]>tmpmax){
		tmpmax=allstpossiblebases[acgtgapbases[xi]];
	      }else{
		allstpossiblebases[acgtgapbases[xi]]=0;
	      }
	    }
	  }
	  numpossiblebases=1; // we WILL take only one at any rate
	}
	priv_mic_helper_chooseOneBase(thisbase, thisqual, allstpossiblebases, groupsvec);
      }
    }


/*
    TODO what about STMU / STMS tags?
    rewrite calc for that
    addTagToConsensus(actcontigpos,
		      actcontigpos,
		      '=',
		      multitag_t::getIdentifierStr(CON_tagentry_idSTMS).c_str(),
		      tagstr.str().c_str(),
		      true);
*/


//    if(hassomething){
//      BUGIFTHROW(thisbase=='@',"Coverage present but base is '@' ... this is not healthy.\n");
//    }else{
//      BUGIFTHROW(thisbase!='@',"No coverage present but base is not '@' ... this is not healthy.\n");
//      thisbase=missingcoveragechar;
//      thisqual=0;
//    }

    if(numpossiblebases>1 && allstpossiblebases['*']){
      thisbase=tolower(thisbase);
      CEBUG("gapasgroup, tolower: " << actcontigpos << endl);
    }else if(thisqual<=30){
      thisbase=tolower(thisbase);
      CEBUG("qual <= 30, tolower: " << actcontigpos << endl);
    }

    CEBUG("thisbase " << thisbase << "\tthisqual " << (uint16) thisqual << endl);

    if(thisqual<minqual && thisbase!=missingcoveragechar){
      thisbase='N';
      thisqual=0;
      CEBUG("Minqual not reached, changed to: thisbase " << thisbase << "\tthisqual " << (uint16) thisqual << endl);
    }

    if(thisqual>90) thisqual=90;

    if(addconstag && strainidtotake >= 0
       && numpossiblebases>1){

      countsofpossiblebases.clear();
      for(auto & labase : acgtgapbases){
	if(allstpossiblebasesf[labase] >= 7
	   && allstpossiblebasesr[labase] >= 7) {
	  countsofpossiblebases.push_back(allstpossiblebases[labase]);
	}
      }
      // use non-parallel sort, this is a tiny/smallish thing
      mstd::ssort(countsofpossiblebases,std::greater<uint32>());
      auto sumcounts=mstd::accumulate(countsofpossiblebases);

      //cout << "ACP: " << actcontigpos << "\t" << sumcounts << "\t" << countsofpossiblebases[0] << "\t" << countsofpossiblebases[1] << "\t" << 100.0/sumcounts*countsofpossiblebases[1] << endl;

      // TODO: change that 0.0 to a configurable number
      if(countsofpossiblebases.size() > 1
	 && 100.0/sumcounts*countsofpossiblebases[1] >= 0.0){
	//std::string consastr("Variants within "+ReadGroupLib::getStrainOfStrainID(strainidtotake)+": ");
	std::string consastr("Hetlevel "
			     + boost::lexical_cast<std::string>(countsofpossiblebases.size())
			     + " within "+ReadGroupLib::getStrainOfStrainID(strainidtotake)+": ");
	bool firstentry=true;
	for(auto & labase : acgtgapbases){
	  if(allstpossiblebases[labase]){
	    if(!firstentry) consastr+=" / ";
	    firstentry=false;
	    consastr+=labase;
	    consastr+=" (";
	    consastr+=boost::lexical_cast<std::string>(allstpossiblebases[labase]);
	    consastr+=')';
	  }
	}
	auto & newtag = addTagToConsensus(actcontigpos,
					  actcontigpos,
					  '=',
					  multitag_t::getIdentifierStr(CON_tagentry_idSAOc).c_str(),
					  consastr.c_str(),
					  true);
	newtag.addAdditionalInfo(allstpossiblebasesf['A'],
				 allstpossiblebasesf['C'],
				 allstpossiblebasesf['G'],
				 allstpossiblebasesf['T'],
				 allstpossiblebasesf['*'],
				 allstpossiblebasesr['A'],
				 allstpossiblebasesr['C'],
				 allstpossiblebasesr['G'],
				 allstpossiblebasesr['T'],
				 allstpossiblebasesr['*'],
				 0,0,0,0,0
	  );
      }
    }

    target+=thisbase;
    qual.push_back(thisqual);

    // calc the adjustments
    if(targetadjustments!= nullptr && from==0 && to==CON_counts.size()) {
      if(thisbase=='*') {
	targetadjustments->push_back(-1);
      }else{
	targetadjustments->push_back(unpaddedposcounter);
	unpaddedposcounter++;
      }

      BUGIFTHROW(targetadjustments->size()!=target.size(), "gna1");
      BUGIFTHROW(targetadjustments->size()!=qual.size(), "gna2");

    }

/*
    // dump out .tcs file to ostr if given
    if(ostr!=nullptr) {
      //throw Notify(Notify::INTERNAL, THISFUNC, "must be adapted to multiple read types");

      // name and position
      // don't call paddedPos2UnpaddedPos as this would lead to recursion!
      *ostr << setw(20) << left << getContigName()
	    << setw(9) << right << actcontigpos
	    << setw(9) << CON_adjustments.back()
	    << " | " << thisbase
	    << setw(3) << static_cast<uint16>(thisqual)
	    << " |" << setw(5) << read_pcrIs_in_col.size();

      tcs_totalgroupcount.clear();
      tcs_totalgroupcount.resize(emptygroups.size(),0);
      tcs_totalgroupqual.clear();
      tcs_totalgroupqual.resize(emptygroups.size(),0);
      for(uint32 actseqtype=0; actseqtype<groupsvec.size(); actseqtype++){
	for(uint32 grpi=0; grpi<groupsvec[actseqtype].size(); grpi++){
	  tcs_totalgroupcount[grpi]+=groupsvec[actseqtype][grpi].urdids.size();
	  tcs_totalgroupqual[grpi]+=groupsvec[actseqtype][grpi].groupquality;
	  if(tcs_totalgroupqual[grpi]>90) tcs_totalgroupqual[grpi]=90;
	}
      }
      for(uint32 grpi=0; grpi<tcs_totalgroupcount.size(); grpi++){
	*ostr << setw(5) << tcs_totalgroupcount[grpi];
      }

      *ostr << " |";
      uint32 numgroupswithqual=0;
      for(uint32 grpi=0; grpi<tcs_totalgroupcount.size(); grpi++){
	if(tcs_totalgroupcount[grpi]) {
	  *ostr << setw(3) << tcs_totalgroupqual[grpi];
	  numgroupswithqual++;
	}else{
	  *ostr << " --";
	}
      }

      // TODO: different characters for different cases?
      char bstatus='?';
      if(!dptools::isValidBase(thisbase)
	 && dptools::isValidIUPACBase(thisbase)){
	bstatus='M';
      }else if(thisbase=='*'){
	if(thisqual>=40){
	  bstatus=':';
	}else if(thisqual<30){
	  bstatus='m';
	}else if(numgroupswithqual>2){
	  bstatus='m';
	}else{
	  bstatus=':';
	}
      }else if(thisqual<30 && numgroupswithqual>1){
	bstatus='m';
      }else{
	bstatus=':';
      }
      if(!tcs_aconstagpositions.empty() && tcs_dconstagpositions[actcontigpos]){
	bstatus='$';
      }
      if(bstatus!=':'){
	*ostr << " | !" << bstatus << " |";
      }else{
	*ostr << " |  : |";
      }

      if(!tcs_aconstagpositions.empty() && tcs_aconstagpositions[actcontigpos]){
        auto I=CON_consensus_tags.cbegin();
	bool doneoutput=false;
	for(; I!=CON_consensus_tags.cend(); I++) {
	  if((I->from <= actcontigpos) && ((I->to) >= actcontigpos)) {
	    if(doneoutput){
	      *ostr << ' ' << I->getIdentifierStr();
	    }else{
	      *ostr << " \"" << I->getIdentifierStr();
	      doneoutput=true;
	    }
	  }
	}
	if(doneoutput) *ostr << "\"";

      }
      *ostr << "\n";
    }
*/

    mict_restofloop+=diffsuseconds(us_loop);
  }
  mict_totalloop=diffsuseconds(us_start);

  if(CON_verbose){
    cout << "mict_fin        " << mict_fin << endl;
    cout << "mict_pre        " << mict_pre << endl;
    cout << "mict_shadow     " << mict_shadow << endl;
    cout << "mict_fallout    " << mict_fallout << endl;
    cout << "mict_newin      " << mict_newin << endl;
#ifndef CLOCK_STEPS_CONS
    cout << "mict_helper1    " << mict_helper1 << endl;
#else
    cout << "mict_helper1    " << CON_us_steps_cons[USCLOCONS_H1_TOTAL] << endl;
    cout << "  H1_MGROUPS    " << CON_us_steps_cons[USCLOCONS_H1_MGROUPS] << endl;
    cout << "  H1_PCRI       " << CON_us_steps_cons[USCLOCONS_H1_PCRI] << endl;
    cout << "    H1_PL_1     " << CON_us_steps_cons[USCLOCONS_H1_PL_1] << endl;
    cout << "    H1_PL_2     " << CON_us_steps_cons[USCLOCONS_H1_PL_2] << endl;
    cout << "    H1_PL_3     " << CON_us_steps_cons[USCLOCONS_H1_PL_3] << endl;
    cout << "    H1_PL_4     " << CON_us_steps_cons[USCLOCONS_H1_PL_4] << endl;
    cout << "    count       " << CON_us_steps_cons[USCLOCONS_H1_countPCRI] << endl;
    cout << "  H1_EGROUP     " << CON_us_steps_cons[USCLOCONS_H1_EGROUP] << endl;
    cout << "  H1_GQUAL      " << CON_us_steps_cons[USCLOCONS_H1_GQUAL] << endl;
    cout << "  H1_CALLH2     " << CON_us_steps_cons[USCLOCONS_H1_CALLH2] << endl;
#endif
    cout << "mict_newin      " <<   mict_newin << endl;
    cout << "mict_restofloop " <<   mict_restofloop << endl;
    cout << "mict_totalloop  " <<   mict_totalloop << endl;
  }

  // seconds on contig with ~4m reads
  //mict_pre                  55
  //mict_shadow           536072
  //mict_fallout        5 855515
  //mict_newin            330100
  //mict_helper1       28 934855
  //mict_restofloop       288773
  //mict_totalloop     35 875607


  CEBUG("target (first 100): " << target.substr(0,100) << endl);

  VCOUT("done." << endl);

  FUNCEND();

  return;
}

//#define CEBUG(bla)


///*************************************************************************
// *
// *
// *
// *************************************************************************/
/*
void Contig::priv_mic_helper_filterGoodnessLevels(cccontainer_t::const_iterator ccI, int8 maxacceptlevel, std::vector<std::vector<int8>> & predictlevels, std::vector<std::vector<nngroups_t>> & allgroups, std::vector<uint32> & allstpossiblebases, uint8 & numpossiblebases)
{
  allstpossiblebases['A']=0;
  allstpossiblebases['C']=0;
  allstpossiblebases['G']=0;
  allstpossiblebases['T']=0;
  allstpossiblebases['*']=0;
  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); ++actseqtype){
    for(uint32 gi=0; gi<predictlevels[actseqtype].size(); ++gi){
      if(predictlevels[actseqtype][gi]>maxacceptlevel) {
	predictlevels[actseqtype][gi]=127;
      }else{
	allstpossiblebases["ACGT*"[gi]]+=allgroups[actseqtype][gi].urdids.size();
	if(actseqtype==ReadGroupLib::SEQTYPE_SOLEXA
	   && CON_isbackbonecontig
	   && toupper(ccI->getOriginalBBChar())==allgroups[actseqtype][gi].base){
	  allstpossiblebases["ACGT*"[gi]]+=ccI->bbcounts[0];
	}
      }
    }
  }
  numpossiblebases=(allstpossiblebases['A']>0)+(allstpossiblebases['C']>0)+(allstpossiblebases['G']>0)+(allstpossiblebases['T']>0)+(allstpossiblebases['*']>0);

  CEBUG("PB-A: " << static_cast<uint16>(allstpossiblebases['A']) << " PB-C: " << static_cast<uint16>(allstpossiblebases['C']) << " PB-G: " << static_cast<uint16>(allstpossiblebases['G']) << " PB-T: " << static_cast<uint16>(allstpossiblebases['T']) << " PB-*: " << static_cast<uint16>(allstpossiblebases['*']) << endl);
  CEBUG("numpossiblebases after filter " << static_cast<int16>(maxacceptlevel) << ": " << static_cast<uint16>(numpossiblebases) << endl);

  return;
}
//#define CEBUG(bla)
*/


/*************************************************************************
 *
 * filtersxa: necessary for bad MiSeq data I've seen
 *   If set, take highest group sxa qual and filter out all that have
 *   >30 difference (1000x less likely)
 *
 *************************************************************************/

void Contig::priv_mic_helper_filterGoodnessLevels(cccontainer_t::const_iterator ccI, int8 maxacceptlevel, std::vector<std::vector<int8>> & predictlevels, std::vector<std::vector<nngroups_t>> & allgroups, std::vector<uint32> & allstpossiblebases, uint8 & numpossiblebases, bool filtersxa)
{
  allstpossiblebases['A']=0;
  allstpossiblebases['C']=0;
  allstpossiblebases['G']=0;
  allstpossiblebases['T']=0;
  allstpossiblebases['*']=0;

  base_quality_t maxsxaqual = 0;
  base_quality_t takesxaqual = 0;
  if ( filtersxa ) {
    const uint32 sxast=ReadGroupLib::SEQTYPE_SOLEXA;

    for(uint32 gi=0; gi<predictlevels[sxast].size(); ++gi){
      maxsxaqual = std::max(maxsxaqual, allgroups[sxast][gi].groupquality);
      CEBUG("mq " << static_cast<int16>(maxsxaqual) << "\t" << static_cast<int16>(allgroups[sxast][gi].groupquality)
	    << "\t" << "ACGT*"[gi] << '\n');
    }
    if ( maxsxaqual > 30 ) takesxaqual = maxsxaqual - 30;
  }

  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); ++actseqtype){
    for(uint32 gi=0; gi<predictlevels[actseqtype].size(); ++gi){
      if(predictlevels[actseqtype][gi]>maxacceptlevel) {
	predictlevels[actseqtype][gi]=127;
      }else{
	if ( actseqtype==ReadGroupLib::SEQTYPE_SOLEXA ){
	  if ( allgroups[actseqtype][gi].groupquality < takesxaqual ){
	    predictlevels[actseqtype][gi]=127;
	  } else {
	    allstpossiblebases["ACGT*"[gi]]+=allgroups[actseqtype][gi].urdids.size();
	    if(CON_isbackbonecontig
	       && toupper(ccI->getOriginalBBChar())==allgroups[actseqtype][gi].base){
	      allstpossiblebases["ACGT*"[gi]] += ccI->bbcountsf[0] + ccI->bbcountsr[0];
	    }
	  }
	} else {
	  allstpossiblebases["ACGT*"[gi]]+=allgroups[actseqtype][gi].urdids.size();
	}
      }
    }
  }
  numpossiblebases=(allstpossiblebases['A']>0)+(allstpossiblebases['C']>0)+(allstpossiblebases['G']>0)+(allstpossiblebases['T']>0)+(allstpossiblebases['*']>0);

  CEBUG("PB-A: " << static_cast<uint16>(allstpossiblebases['A']) << " PB-C: " << static_cast<uint16>(allstpossiblebases['C']) << " PB-G: " << static_cast<uint16>(allstpossiblebases['G']) << " PB-T: " << static_cast<uint16>(allstpossiblebases['T']) << " PB-*: " << static_cast<uint16>(allstpossiblebases['*']) << endl);
  CEBUG("numpossiblebases after filter lvl " << static_cast<int16>(maxacceptlevel)
	<< "tsxa " << filtersxa << ' ' << static_cast<int16>(takesxaqual)
	<< ": " << static_cast<uint16>(numpossiblebases) << endl);

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * returns by value: best level encountered
 * implicitly: predictlevels
 *
 *************************************************************************/

int8 Contig::priv_mic_helper_rateGoodnessLevelOfAllGroups(cccontainer_t::const_iterator ccI, const std::vector<std::vector<nngroups_t>> & allgroups, std::vector<std::vector<int8>> & predictlevels, const std::vector<std::vector<char>> & possiblebases)
{
  FUNCSTART("int8 Contig::priv_mic_helper_rateGoodnessLevelOfAllGroups(cccontainer_t::const_iterator ccI, const std::vector<std::vector<nngroups_t>> & allgroups, std::vector<std::vector<int8>> & predictlevels, const std::vector<std::vector<char>> & possiblebases)");

  int8 minlevel=127;
  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); ++actseqtype){
    predictlevels[actseqtype].clear();
    predictlevels[actseqtype].resize(5,127);
    for(uint32 gi=0; gi<allgroups[actseqtype].size(); ++gi){
      predictlevels[actseqtype][gi]=priv_mic_helper_rateGoodnessLevelOfGroup(
	ccI,
	allgroups[actseqtype][gi],
	possiblebases[actseqtype].size(),
	actseqtype);
      BUGIFTHROW(predictlevels[actseqtype][gi]<0, "problem in hybrid cons calc: level < 0");
      minlevel=std::min(minlevel,predictlevels[actseqtype][gi]);
    }
  }
  return minlevel;
}

/*************************************************************************
 *
 *
 *
 *
 *
 *************************************************************************/

int8 Contig::priv_mic_helper_rateGoodnessLevelOfGroup(cccontainer_t::const_iterator ccI, const nngroups_t & group, uint32 numpossiblebases, uint8 seqtype)
{
  FUNCSTART("int8 Contig::rateGoodnessLevelOfGroup(cccontainer_t::const_iterator ccI, nngroups_t & group, uint32 numpossiblebases, uint8 seqtype)");

  // empty group? quick dispatch
  if(group.urdids.empty()) return 127;

  // basically, we'll prefer read types that have forward and
  //  complement direction present.

  CEBUG("predicting level for " << group.base << " seqtype " << static_cast<uint16>(seqtype) << '\n' << group);
  int8 level=127;

  if(seqtype==ReadGroupLib::SEQTYPE_ABISOLID){
    MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 3.");
  }


  bool hasfr=(group.forwarddircounter>0) && (group.complementdircounter>0);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  switch(seqtype) {
  case ReadGroupLib::SEQTYPE_SANGER :{
    if(hasfr &&
       group.groupquality >= 35){
      level=0;
    }else if(hasfr
	     && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=3
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_454GS20 :{
    if(hasfr
       && group.urdids.size() >= 12
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 8
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=6
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_IONTORRENT :{
    if(hasfr
       && group.urdids.size() >= 8
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 6
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=4
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_PACBIOHQ :
  case ReadGroupLib::SEQTYPE_PACBIOLQ :{
    // No info. atm use same as 454
    // TODO: change that for pacbio hq and lq
    if(hasfr
       && group.urdids.size() >= 12
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 8
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=6
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_TEXT :{
    // TODO: that's just a first shot at atm, would need to refine
    if(hasfr
       && group.urdids.size() >= 8
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 6
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=4
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_SOLEXA :{
    // Solexa special case: we need to take the merged reads
    //  into account
    uint32 realsize=group.urdids.size();
    if(ccI->getOriginalBBChar()==group.base){
      realsize += ccI->bbcountsf[0] + ccI->bbcountsr[0];
      if(ccI->bbcountsf[0] && ccI->bbcountsr[0]){
	hasfr=true;
      }
    }
    if(hasfr
       && realsize >= 12
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && realsize >= 8
	      && group.groupquality >= 30){
      level=1;
    } else if(hasfr
	      && realsize >= 6
	      && group.groupquality >= 25){
      level=2;
    } else if(group.groupquality >= 20
	      && realsize >= 4
	      && numpossiblebases == 1){
      level=3;
    } else if(group.groupquality >= 10){
      level=4;
    } else {
      level=5;
    }
    break;
  }
  default : {
    cout << "seqtype: " << seqtype << endl;
    MIRANOTIFY(Notify::INTERNAL, "Unknown seqtype to rate?");
  }
  }

  FUNCEND();

  CEBUG("Level predicted: " << static_cast<uint16>(level) << '\n');

  return level;
}

// MIC all CEBUG end
//#define CEBUG(bla)




/*************************************************************************
 *
 * makes a quick consensus calculation of a portion of the contig
 *  stored in CON_2tmpcons (std::string)
 *  String can have ACGTN. But also 1234. '1'=a, 2='c', 3='g', 4='t'
 *    Used to denote "that base or a gap"
 *
 * from and to positions in the consensus
 * from inclusive
 * to exclusive
 * usenumgaps: use 1/2/3/4 to denote a/c/g/t-or-gap
 *
 * returns bool: whether a base not defined by *original* backbone or, if backbone,
 *  non-ACGT* was put into consensus
 *  Used to define whether read can be mapped or not.
 *
 *************************************************************************/

bool Contig::makeTmpConsensus(int32 from, int32 to, bool tmpconsfrombackbone, bool usenumgaps)
{
//#define CEBUG(bla)   {cout << bla; cout.flush();}

  FUNCSTART("void Contig::makeTmpConsensus(int32 from, int32 to)");

  definalise();
  CEBUG("\nFrom: " << from << "\tTo: " << to<<"\tCON_counts.size(): "<<CON_counts.size());

  BUGIFTHROW(from>to,"from>to?");
  BUGIFTHROW(from<0, "from < 0 ?");

  if(to>CON_counts.size()) to=CON_counts.size();
  uint32 len_tmpcons=to-from;

  if(CON_2tmpcons.capacity()<2000 || CON_2tmpcons.capacity()<len_tmpcons){
    CON_2tmpcons.reserve(std::max(len_tmpcons,static_cast<uint32>(2000)));
  }
  CON_2tmpcons.resize(len_tmpcons);

  bool hasNonBBMappable=!tmpconsfrombackbone;
  {
    auto toptr=CON_2tmpcons.begin();
    auto ccI=CON_counts.cbegin();
    BOUNDCHECK(from, 0, CON_counts.size()+1);
    std::advance(ccI, from);

    for(uint32 i=from; i<to ;++i, ++toptr, ++ccI){
      *toptr=0;
      if(tmpconsfrombackbone && ccI->getOriginalBBChar()!='@'){
	*toptr=ccI->getBBChar();
	hasNonBBMappable |= !dptools::isValidACGTStarBase(ccI->getOriginalBBChar());
	if(ccI->i_backbonecharupdated!='@'
	   && ccI->i_backbonecharorig != ccI->i_backbonecharupdated) hasNonBBMappable=true;
      }else{
	hasNonBBMappable=true;

	ccctype_t maximum= std::max(ccI->A, std::max(ccI->C, std::max(ccI->G, ccI->T)));
	uint8 counts=0;
	//CEBUGF(ccI->A << "\t" << ccI->C << "\t" << ccI->G << "\t" << ccI->T << "\t" << ccI->N << "\t" << ccI->star << "\n");

	// is any ACGT set?
	if(maximum >0 && maximum > ccI->star) {
	  if(ccI->A==maximum){
	    counts++;
	    *toptr='A';
	  }
	  if(ccI->C==maximum){
	    counts++;
	    *toptr='C';
	  }
	  if(ccI->G==maximum){
	    counts++;
	    *toptr='G';
	  }
	  if(ccI->T==maximum){
	    counts++;
	    *toptr='T';
	  }
	  if(counts>1) {
	    *toptr='N';
	  }

	  //// can be somewhat problematic with 454 data
	  //// calls the base until 50/50, then the gap
	  //if(maximum/4 < ccI->star) *toptr='*';

	  // this prefers to call gaps
	  // calls the base until 1/3 base, 2/3 gap, then the gap
	  // this should help, together with the "expected gap" # in
	  //  alignments, to further reduce to a maximum this kind of
	  //  base jiggling in homopolymers
	  //
	  //     ...*AAAAAAAAA...
	  //     ...*AAAAAAAAA...
	  //     ...AAAAAAAAA*...
	  //     ...AAAAAAAAA*...
	  //     ...AAAAAAAAA*...
	  //     ...*AAAAAAAAA...
	  //     ...*AAAAAAAAA...
	  //     ...AAAAAAAAA*...
	  //     ...*AAAAAAAAA...

	  if(usenumgaps && maximum/4 < (ccI->star)*2) {
	    switch(*toptr){
	    case 'A': {
	      *toptr='1';
	      break;
	    }
	    case 'C': {
	      *toptr='2';
	      break;
	    }
	    case 'G': {
	      *toptr='3';
	      break;
	    }
	    case 'T': {
	      *toptr='4';
	      break;
	    }
	    default: {
	      *toptr='*';
	    }
	    }
	  }
	} else {
	  if(unlikely(ccI->total_cov==0)){
	    // BaCh 30.11.2012
	    // should normally never happen, certainly not in de-novo
	    // but the two-pass mapping may have this at the end of the contigs after first pass
	    //  (should I decide not to go the chompFront() / chompBack() after 1st pass)
	    //
	    // treat it like a base (well, will be X)
	    *toptr='N';
	  }else if((ccI->star >= ccI->X)
	     && (ccI->star >= ccI->N)){
	    *toptr='*';
	  } else if(ccI->N){
	    *toptr='N';
	  }else{
	    *toptr='X';
	  }
	}
      }
      if(*toptr==0){
	MIRANOTIFY(Notify::INTERNAL,"Ooooops? makeTmpConsensus encountered the unexpected situation of an uncalled base? Please contact the author immediately.");
      }
    }
  }
  CEBUG("Tmp_cons: >>>" <<CON_2tmpcons << "<<<" << endl);

  FUNCEND();

  return hasNonBBMappable;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::updateBackboneConsensus()
{
  makeTmpConsensus(0,CON_counts.size(),false,true);
  auto ccI=CON_counts.begin();
  auto c2tc=CON_2tmpcons.c_str();

  // TODO: check for coverage, check for fwd,rev
  rcci_t rcci(this,
	      nullptr,  // all strainids
	      nullptr,  // all seqtypes
	      false,            // don't take rails
	      false,           // nor backbones
	      false);   // nor reads without readpool-reads

  for(; ccI != CON_counts.end(); ++ccI, ++c2tc, rcci.advance()){
    if(rcci.getPCRIsInCol().size()>=6){
      uint32 fwd=0;
      uint32 rev=0;
      for(auto & pcrI : rcci.getPCRIsInCol()){
	if(pcrI.getReadDirection()>0){
	  ++fwd;
	}else{
	  ++rev;
	}
      }
      if(fwd>=2 && rev>=2){
	if(ccI->i_backbonecharorig!='@'){
	  ccI->i_backbonecharupdated=*c2tc;
	}
      }
    }
  }
}
