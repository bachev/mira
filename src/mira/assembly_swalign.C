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


#include "mira/assembly.H"

// BOOST
#include <boost/algorithm/string.hpp>

#include "errorhandling/errorhandling.H"
#include "util/progressindic.H"
#include "util/dptools.H"
#include "util/fileanddisk.H"
#include "caf/caf.H"

#include "mira/align.H"
#include "mira/ads.H"

#include "util/stlimprove.H"


#if 0
#include <valgrind/memcheck.h>
#define VALGRIND_LEAKCHECK
#endif


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)

// cs1 for normal clocking ('user compatible' as is does not disturb)
//  cs2 for extensive clocking output, more for analysis of MIRA behaviour
//  during development

#ifndef PUBLICQUIET
#define CLOCK_STEPS2
#endif
#define CLOCK_STEPS2


//#define TRACKMEMUSAGE 1
#define TRACKMEMUSAGE 0





/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::setupAlignCache(std::vector<Align> & aligncache)
{
  for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++) {
    Align a(&AS_miraparams[i]);
    aligncache.push_back(a);
  }
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}
//#define ALIGNCHECK

// testing
#define DEBUGEND_L '\n'

void Assembly::makeAlignmentsFromPosMatchFile(const std::string & filename, const int32 version, const int8 direction, const bool trans100percent, bool (* checkfunction)(Assembly & as,int32,int32) , std::ofstream & matchfout, std::ofstream & rejectfout)
{
  FUNCSTART("void Assembly::makeAlignmentsFromPosMatchFile(const std::string & filename)");

  BUGIFTHROW(filename.empty(),"filename.empty() ???");
  CEBUG("makeAlignmentsFromPosMatchFile: " << filename << "\ttrans100percent=" << trans100percent << endl);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt! Also adapt if other aligns are used."
#endif

  std::vector<Align> chkalign;
  setupAlignCache(chkalign);

  std::list<AlignedDualSeq> madsl;

  uint32 potentialalignments=0;
  uint32 totalseqsaligned=0;
  uint32 permbansevaded=0;
  uint32 checkfunrejected=0;
  uint32 trans100saved=0;

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  std::ifstream posffin(filename, std::ios::in|std::ios::ate|std::ios::binary);

  if(!posffin){
    MIRANOTIFY(Notify::FATAL, "File not found. This should have been written earlier by MIRA: " << filename);
  }


  ProgressIndicator<std::streamsize> P (0, posffin.tellg(),2000);

  posffin.seekg(0, std::ios::beg);

#ifdef ALIGNCHECK
  Align checkbla(&AS_miraparams[ReadGroupLib::SEQTYPE_SANGER]);
  // true for using memcache
  //Align checkbla(AS_miraparams, true);
#endif

//  struct matchwithsorter_t{
//    uint32 rid1;
//    uint32 rid2;
//    int32  eoffset;
//    int32  percent_in_overlap;
//    uint32 numhashes;
//  } posmatch;

  skimhitforsave_t posmatch;

  while(!posffin.eof()){
    //CEBUG("pindic: " << pindic<<endl);

    posffin.read(reinterpret_cast<char *>(&posmatch),sizeof(posmatch));

    if(posffin.eof()) break;

    if(P.delaytrigger()) P.progress(posffin.tellg());

    potentialalignments++;

    CEBUG("Looking: " << posmatch.rid1 << " " << posmatch.rid2 << "\t" <<  AS_readpool.getRead(posmatch.rid1).getName() << "\t" <<  AS_readpool.getRead(posmatch.rid2).getName() << '\n');

    if(AS_permanent_overlap_bans.checkIfBanned(posmatch.rid1,posmatch.rid2) > 0) {
      CEBUG("PermBan for: " << posmatch.rid1 << " " << posmatch.rid2<<"\tskipping\n");
      permbansevaded++;
      continue;
    }

    if(AS_readpool.getRead(posmatch.rid1).isRail()
       && AS_readpool.getRead(posmatch.rid2).isRail()) {
      CEBUG("Both are rails: " << posmatch.rid1 << " " << posmatch.rid2<<"\tskipping\n");
      continue;
    }

    // version 0 == pre-assembly pass for vector clipping and/or
    //  read extension
    if((version >0 && version <as_fixparams.as_startbackboneusage_inpass)
       && (AS_readpool.getRead(posmatch.rid2).isRail()
	   || AS_readpool.getRead(posmatch.rid2).isRail())){
      CEBUG("One is rail and pass < startbackboneusage: " << posmatch.rid1 << " " << posmatch.rid2<<"\tskipping\n");
      continue;
    }

    // normally the sequences should have a length >0
    // but due to some clipping being done after SKIM (chimera etc.), it
    //  may happen they are 0 now. If that's the case, discard this possible match
    if(AS_readpool[posmatch.rid1].getLenClippedSeq() == 0
       || AS_readpool[posmatch.rid2].getLenClippedSeq() == 0) continue;

    if(!checkfunction(*this,posmatch.rid1,posmatch.rid2)){
      CEBUG("Read combination rejected by check function.\n");
      checkfunrejected++;
      continue;
    }

    // don't use the 100% transfer rule if
    //  - not 100% (d'oh)
    //  - a SRMr tag present or in each read a CRMr
    bool canuse100perctrans=trans100percent;
    if(posmatch.percent_in_overlap != 100
       || AS_readpool.getRead(posmatch.rid1).hasTag(Read::REA_tagentry_idSRMr)
       || AS_readpool.getRead(posmatch.rid2).hasTag(Read::REA_tagentry_idSRMr)){
      canuse100perctrans=false;
    }else if(AS_readpool.getRead(posmatch.rid1).hasTag(Read::REA_tagentry_idCRMr)
	     && AS_readpool.getRead(posmatch.rid2).hasTag(Read::REA_tagentry_idCRMr)){
      canuse100perctrans=false;
    }

    if(canuse100perctrans){
      CEBUG("100% trans rule.\n");
      if(matchfout.is_open()){
	matchfout << AS_readpool.getRead(posmatch.rid1).getName() << "\t" <<AS_readpool.getRead(posmatch.rid2).getName() << DEBUGEND_L;
      }

      //AS_CUMADSLofstream << madsl.begin()->getWeight() << '\t'
      //			 << static_cast<int16>(direction) << '\t';
      //madsl.begin()->serialiseOut(AS_CUMADSLofstream);
      //AS_CUMADSLofstream << '\n';

      bool swapped=false;
      if(posmatch.eoffset<0){
	// swap if this
	//   1             --------
	//   2      -----------
	std::swap(posmatch.rid1, posmatch.rid2);
	posmatch.eoffset=-posmatch.eoffset;
	swapped=true;
      }

      int32 overlaplen;
      int32 totallen;

      int32 rdls=AS_readpool.getRead(posmatch.rid1).getLenClippedSeq()-posmatch.eoffset-AS_readpool.getRead(posmatch.rid2).getLenClippedSeq();
      // only two cases left due to swapping of ids above
      if(rdls>=0) {
	//   1      -----------------
	//   2           ----------
	overlaplen=AS_readpool.getRead(posmatch.rid2).getLenClippedSeq();
	totallen=AS_readpool.getRead(posmatch.rid1).getLenClippedSeq();
      }else{
	//   1      -----------------
	//   2           --------------
	overlaplen=AS_readpool.getRead(posmatch.rid2).getLenClippedSeq()+rdls;
	totallen=posmatch.eoffset+AS_readpool.getRead(posmatch.rid2).getLenClippedSeq();
      }

      // overlap must be >= smallest allowed minimal overlap
      // ... or simply >= 17
      // BaCh: 06.02.2015; nope, not "or >= 17" as this does not reflect the wish of the user!
      if(overlaplen >= AS_miraparams[AS_readpool[posmatch.rid1].getSequencingType()].getAlignParams().al_min_overlap
	 || overlaplen >= AS_miraparams[AS_readpool[posmatch.rid2].getSequencingType()].getAlignParams().al_min_overlap){

	AS_CUMADSLofstream << overlaplen*10000 << '\t'
			   << static_cast<int16>(direction) << '\t'
			   << posmatch.ol_stronggood << '\t'
			   << posmatch.ol_weakgood << '\t'
			   << posmatch.ol_belowavgfreq << '\t'
			   << posmatch.ol_norept << '\t'
			   << posmatch.ol_rept << '\t';

	AS_CUMADSLofstream << posmatch.rid1
			   << '\t' << posmatch.rid2;

	if(swapped){
	  AS_CUMADSLofstream << '\t' << static_cast<int16>(direction);
	  AS_CUMADSLofstream << "\t1";
	}else{
	  AS_CUMADSLofstream << "\t1";
	  AS_CUMADSLofstream << '\t' << static_cast<int16>(direction);
	}
	AS_CUMADSLofstream << '\t' << posmatch.eoffset;

	if(rdls>=0) {
	  //   1      -----------------
	  //   2           ----------
	  AS_CUMADSLofstream << '\t' << 0;
	  AS_CUMADSLofstream << '\t' << rdls;
	}else{
	  //   1      -----------------
	  //   2           --------------
	  AS_CUMADSLofstream << '\t' << -rdls;
	  AS_CUMADSLofstream << '\t' << 0;
	}
	AS_CUMADSLofstream << '\t' << overlaplen
			   << '\t' << totallen
			   << '\t' << static_cast<uint16>(posmatch.percent_in_overlap);

	if(rdls>=0) {
	  //   1      -----------------
	  //   2           ----------
	  AS_CUMADSLofstream << "\t0\t0\t0\t7\t7";
	}else{
	  //   1      -----------------
	  //   2           --------------
	  if(direction>0){
	    AS_CUMADSLofstream << "\t0\t0\t7\t7\t0";
	  }else if(swapped){
	    AS_CUMADSLofstream << "\t0\t7\t0\t7\t0";
	  }else{
	    AS_CUMADSLofstream << "\t0\t0\t7\t0\t7";
	  }
	}

	AS_CUMADSLofstream << '\n';

	if(AS_CUMADSLofstream.bad()){
	  MIRANOTIFY(Notify::FATAL, "Could not write anymore to disk (at 100%trans). Disk full? Changed permissions?");
	}

	AS_numADSFacts_fromalignments++;
	trans100saved++;
      }else{
	// smaller, reject
	// well, do nothing for now, perhaps increase a counter later
      }
    }else{

      computeSWAlign(madsl, posmatch.rid1, posmatch.rid2, posmatch. eoffset, direction, chkalign, -1);

      totalseqsaligned++;

      CEBUG("Solutions found: " << madsl.size() << '\n');

      //if(madsl.size()) cout << madsl.front();

#ifdef ALIGNCHECK
      CEBUG("Alignment: " << AS_readpool.getRead(posmatch.rid1).getName() << " and " << AS_readpool.getRead(posmatch.rid2).getName());
      if(madsl.size()>0){
	CEBUG(" found\n");
	cout <<" ----------------------------------------------------- \n";
	cout <<"# solutions found: "<< madsl.size() << endl;
	for(const auto & adse : madsl){
	  cout <<*adse;
	}
	cout << " ----------------------------------------------------- \n";
      }else{
	CEBUG(" missed\n");

	std::list<AlignedDualSeq> tadsl;
	if(direction>0){
	  checkbla.acquireSequences(
	    static_cast<const char *>(AS_readpool.getRead(posmatch.rid1).getClippedSeqAsChar()),
	    AS_readpool.getRead(posmatch.rid1).getLenClippedSeq(),
	    static_cast<const char *>(AS_readpool.getRead(posmatch.rid2).getClippedSeqAsChar()),
	    AS_readpool.getRead(posmatch.rid2).getLenClippedSeq(),
	    posmatch.rid1,
	    posmatch.rid2,
	    1,
	    1);
	}else{
	  checkbla.acquireSequences(
	    static_cast<const char *> (AS_readpool.getRead(posmatch.rid1).getClippedSeqAsChar()),
	    AS_readpool.getRead(posmatch.rid1).getLenClippedSeq(),
	    static_cast<const char *> (AS_readpool.getRead(posmatch.rid2).getClippedComplementSeqAsChar()),
	    AS_readpool.getRead(posmatch.rid2).getLenClippedSeq(),
	    posmatch.rid1,
	    posmatch.rid2,
	    1,
	    -1);
	}
	checkbla.fullAlign(&tadsl,false,true);

	if(tadsl.size()!=0){
	  cout << "Dammit, Offset-BSW lost a solution!\n";
	  cout << "predicted offset: " << posmatch.eoffset << endl;
	  cout <<" ----------------------------------------------------- \n";
	  cout <<"# solutions found: "<< tadsl.size() << endl;
	  for(const auto & adse : tadsl){
	    cout <<*adse;
	  }
	  cout << " ----------------------------------------------------- \n";
	}
      }
#endif
      if(as_fixparams.as_tmpf_ads.size()!=0){
	if(madsl.size()!=0){
	  //matchfout << posmatch.rid1 << " " << posmatch.rid2 << "\t" << I->second.eoffset << endl;
	  if(matchfout.is_open()){
	    matchfout << AS_readpool.getRead(posmatch.rid1).getName() << "\t" <<AS_readpool.getRead(posmatch.rid2).getName() << DEBUGEND_L;
	  }
	}else{
	  if(rejectfout.is_open()){
	    rejectfout << AS_readpool.getRead(posmatch.rid1).getName() << "\t" << static_cast<int16>(direction) << "\t" <<AS_readpool.getRead(posmatch.rid2).getName() << DEBUGEND_L;
	  }
	}
      }

      cleanupMADSL(madsl, posmatch.rid1, posmatch.rid2, direction,
		   posmatch.ol_stronggood, posmatch.ol_weakgood, posmatch.ol_belowavgfreq,
		   posmatch.ol_norept, posmatch.ol_rept);

      //if(madsl.size()>0 && posmatch.percent_in_overlap==100
      //	&& madsl.front().getScoreRatio() == 99){
      //	cout <<" ----------------------------------------------------- \n";
      //	cout <<"# dingdong found: "<< madsl.size() << endl;
      //	{
      //	  std::list<AlignedDualSeq>::const_iterator Itmp=madsl.begin();
      //	  while(Itmp!=madsl.end()){
      //	    cout <<*Itmp; Itmp++;
      //	  }
      //	}
      //
      //	cout << " ----------------------------------------------------- \n";
      //}


    }
  }
  P.finishAtOnce();
  cout << "\nAlignment stats:";
  cout << "\nPotential:   " << potentialalignments;
  cout << "\nCalculated:  " << totalseqsaligned;
  cout << "\nEvaded (PB): " << permbansevaded;
  cout << "\nRejected (checkfun): " << checkfunrejected;
  cout << "\nTrans 100 saved: " << trans100saved;


  // count banned pairs
  {
    size_t banned=0;
    size_t numsets=0;

    AS_permanent_overlap_bans.getNumBans(banned,numsets);

    cout << "\n\nBanned overlap pairs: " << banned
	 << "\tin " << numsets << " sets.";
  }

  cout << "\n\n";
}
//#define CEBUG(bla)
//#define CEBUGF(bla)


/*************************************************************************
 *
 * TODO: handling of hintbandwidth is a cludge atm and will not work for
 *  parallel Smith-Watermans (for that Align, Dynamic and AlignedDualSeq
 *  will need to use not pointer to MIRAParameters, but pointers to their
 *  own params which can be changed on the fly)
 *
 * But we're not there yet.
 *
 *************************************************************************/

void Assembly::computeSWAlign(std::list<AlignedDualSeq> & madsl, uint32 rid1, uint32 rid2, int32 eoffset, int8 direction, std::vector<Align> & chkalign, int32 hintbandwidth)
{
  FUNCSTART("void Assembly::computeSWAlign(std::list<AlignedDualSeq> & madsl, uint32 rid1, uint32 rid2, int32 eoffset, int8 direction, std::vector<Align> & chkalign)");

  BUGIFTHROW(AS_needalloverlaps.size()!=AS_readpool.size(),"AS_needalloverlaps.size()!=AS_readpool.size() ???");

  CEBUG("Acquiring: " << rid1 << " "
	<< rid2<< "\teofset: " << eoffset
	<< "\tdirection: " << static_cast<int16>(direction) << '\n');


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt! Also adapt if other aligns are used."
#endif
  // use the sanger align as default
  uint8 usealign=ReadGroupLib::SEQTYPE_SANGER;

  // if any of the reads is PacBioLQ, use the PacBioLQ align
  // else if any of the reads is SOLEXA, use the SOLEXA align
  // else if any of the reads is PacBioHQ, use the PacBioHQ align
  // else if any of the reads is 454, use the 454 align
  // else if any of the reads is Text, use the Text align
  if(AS_readpool.getRead(rid1).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)
     || AS_readpool.getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)){
    usealign=ReadGroupLib::SEQTYPE_PACBIOLQ;
  } else if(AS_readpool.getRead(rid1).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	    || AS_readpool.getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
    usealign=ReadGroupLib::SEQTYPE_SOLEXA;
  }else if(AS_readpool.getRead(rid1).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOHQ)
	   || AS_readpool.getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOHQ)){
    usealign=ReadGroupLib::SEQTYPE_PACBIOHQ;
  }else if(AS_readpool.getRead(rid1).isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)
	   || AS_readpool.getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)){
    usealign=ReadGroupLib::SEQTYPE_IONTORRENT;
  }else if(AS_readpool.getRead(rid1).isSequencingType(ReadGroupLib::SEQTYPE_454GS20)
	   || AS_readpool.getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_454GS20)){
    usealign=ReadGroupLib::SEQTYPE_454GS20;
  }else if(AS_readpool.getRead(rid1).isSequencingType(ReadGroupLib::SEQTYPE_TEXT)
	   || AS_readpool.getRead(rid2).isSequencingType(ReadGroupLib::SEQTYPE_TEXT)){
    usealign=ReadGroupLib::SEQTYPE_TEXT;
  }

  bool enforce_clean_ends=AS_miraparams[usealign].getAlignParams().ads_enforce_clean_ends;
  // if any read is a rail or backbone, do not use the clean ends
  //  requirement. This is to align reads that contain true SNP in
  //  the end positions
  if(AS_readpool.getRead(rid1).isBackbone()
     || AS_readpool.getRead(rid1).isRail()
     || AS_readpool.getRead(rid2).isBackbone()
     || AS_readpool.getRead(rid2).isRail()){
    enforce_clean_ends=false;
  }

  auto & aparams=AS_miraparams[usealign].getNonConstAlignParams();
  int32 savealkmin=aparams.al_kmin;
  int32 savealkmax=aparams.al_kmax;
  if(hintbandwidth>0){
    aparams.al_kmin=hintbandwidth;
    aparams.al_kmax=hintbandwidth;
  }

  try{
    if(direction>0){
      chkalign[usealign].acquireSequences(
	static_cast<const char *> (AS_readpool.getRead(rid1).getClippedSeqAsChar()),
	AS_readpool.getRead(rid1).getLenClippedSeq(),
	static_cast<const char *> (AS_readpool.getRead(rid2).getClippedSeqAsChar()),
	AS_readpool.getRead(rid2).getLenClippedSeq(),
	rid1,
	rid2,
	1,
	1,
	true,
	eoffset);
    }else{
      chkalign[usealign].acquireSequences(
	static_cast<const char *> (AS_readpool.getRead(rid1).getClippedSeqAsChar()),
	AS_readpool.getRead(rid1).getLenClippedSeq(),
	static_cast<const char *> (AS_readpool.getRead(rid2).getClippedComplementSeqAsChar()),
	AS_readpool.getRead(rid2).getLenClippedSeq(),
	rid1,
	rid2,
	1,
	-1,
	true,
	eoffset);
    }
  }
  catch(Notify n) {
    Read::setCoutType(Read::AS_TEXT);
    cout << "Ouch, having a problem here. Tried to acquire the following reads:\n"
	 << AS_readpool.getRead(rid1)
	 << endl
	 << AS_readpool.getRead(rid2)
	 << endl << "with posmatch:\n"
	 << rid1
	 << "\t" << rid2
	 << "\t" << eoffset
	 << endl;
    n.handleError(THISFUNC);
  }
  madsl.clear();

  CEBUG("usealign: " << static_cast<uint16>(usealign) << endl);
  //chkalign[usealign].coutWhatWasGiven();
  if(AS_needalloverlaps[rid1] || AS_needalloverlaps[rid2]){
    CEBUG("go down with requirements\n");
    chkalign[usealign].useSpecialMinRelScore(50);
    chkalign[usealign].setEnforceCleanEnds(false);
    chkalign[usealign].setAffineGapScore(false);
    chkalign[usealign].fullAlign(&madsl);
    chkalign[usealign].useSpecialMinRelScore(0);
  }else{
    CEBUG("normal req align");
    chkalign[usealign].setEnforceCleanEnds(enforce_clean_ends);
    chkalign[usealign].setAffineGapScore(false);
    chkalign[usealign].fullAlign(&madsl);
  }

  if(hintbandwidth>0){
    aparams.al_kmin=savealkmin;
    aparams.al_kmax=savealkmax;
  }

  FUNCEND();
}



/*************************************************************************
 *
 * Check the generated match-table. Pick the possible matching
 *  candidates and check them with the Smith-Waterman.
 * Save the results into these variables:
 *   AS_global_adslist     a list of all dual aligned sequences which
 *                         have been retained by the Align
 *   AS_edges_forward      for forward reads.
 *                         a multimap. key is the number of the read
 *                         in the sequencepool, value is a struct
 *                         showing the destination read and 1 or more
 *                         ADS objects. Best score ratio and pointer to
 *                         ADS element in global_adslist is also stored.
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Assembly::ma_takeall(Assembly & as, int32 rid1, int32 rid2)
{
  // function is called indirectly, tell compiler not to worry
  //  about unused variables in this case
  (void) as;
  (void) rid1;
  (void) rid2;

  return true;
}

bool Assembly::ma_needRRFlag(Assembly & as, int32 rid1, int32 rid2)
{
  if(as.AS_readsforrepeatresolve[rid1]
     || as.AS_readsforrepeatresolve[rid2]){
    return true;
  }
  return false;
}

bool Assembly::ma_needRRFlagAndBothCRMr(Assembly & as, int32 rid1, int32 rid2)
{
  if(as.AS_readsforrepeatresolve[rid1]
     || as.AS_readsforrepeatresolve[rid2]){
    if(as.AS_readpool.getRead(rid1).hasTag(Read::REA_tagentry_idCRMr)
       && as.AS_readpool.getRead(rid2).hasTag(Read::REA_tagentry_idCRMr)) return true;
  }
  return false;
}

// is this one used anymore since the introduction of
//  AS_readsforrepeatresolve???
bool Assembly::ma_needSRMrOrTwoCRMr(Assembly & as, int32 rid1, int32 rid2)
{
  if(as.AS_readpool.getRead(rid1).hasTag(Read::REA_tagentry_idSRMr)) return true;
  if(as.AS_readpool.getRead(rid2).hasTag(Read::REA_tagentry_idSRMr)) return true;
  if(as.AS_readpool.getRead(rid1).hasTag(Read::REA_tagentry_idCRMr)
     && as.AS_readpool.getRead(rid2).hasTag(Read::REA_tagentry_idCRMr)) return true;
  return false;
}

void Assembly::makeAlignments(bool (* checkfunction)(Assembly & as,int32,int32), bool takefullskimfilenames, const bool trans100percent, int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)
{
  FUNCSTART("void Assembly::makeAlignments()");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  cout << "Making alignments.\n";

  std::string adsfacts_fn;
  if(tmpfname.size()){
    adsfacts_fn=buildFileName(version, prefix, postfix, tmpfname, ".adsfacts");
  }else{
    adsfacts_fn=buildFileName(version, prefix, postfix,
			      as_fixparams.as_tmpf_ads,
			      ".adsfacts");
  }
  std::string adsf_fn;
  if(tmpfname.size()){
    adsf_fn=buildFileName(version, prefix, postfix, tmpfname, ".forward");
  }else{
    adsf_fn=buildFileName(version, prefix, postfix,
			  "elog.ads",
			  ".forward");
  }
  std::string adsc_fn;
  if(tmpfname.size()){
    adsc_fn=buildFileName(version, prefix, postfix, tmpfname, ".complement");
  }else{
    adsc_fn=buildFileName(version, prefix, postfix,
			  "elog.ads",
			  ".complement");
  }
  std::string adsr_fn;
  if(tmpfname.size()){
    adsr_fn=buildFileName(version, prefix, postfix, tmpfname, ".reject");
  }else{
    adsr_fn=buildFileName(version, prefix, postfix,
			  "elog.ads",
			  ".reject");
  }

  AS_CUMADSLofstream.open(adsfacts_fn, std::ios::out|std::ios::trunc);
  AS_CUMADSLofstream.close();

  //nukeSTLContainer(AS_adsfacts);
  //nukeSTLContainer(AS_confirmed_edges);
  AS_adsfacts.clear();
  AS_confirmed_edges.clear();


  //directory_parameters const & dir_params= AS_miraparams->getDirectoryParams();

  nukeSTLContainer(AS_readhitmiss);
  nukeSTLContainer(AS_readhmcovered);
  nukeSTLContainer(AS_count_rhm);

  bool cpvwanted=false;
  for(uint32 i=0; i < ReadGroupLib::SEQTYPE_END; i++){
    if(AS_miraparams[i].getAssemblyParams().as_clip_possible_vectors
       && AS_seqtypespresent[i]) cpvwanted=true;
  }

  if(cpvwanted){
    CEBUG("Want cpv\n");
    if(!AS_steps[ASVECTORSCLIPPED]){
      CEBUG("No vectors clipped so far.\n");

#if TRACKMEMUSAGE
      cout << "\ndmi ma 00\n";
      dumpMemInfo();
#endif

      AS_readhitmiss.resize(AS_readpool.size());
      AS_readhmcovered.resize(AS_readpool.size());
      AS_count_rhm.resize(AS_readpool.size());

#if TRACKMEMUSAGE
      cout << "\ndmi ma 10\n";
      dumpMemInfo();
#endif
      for(uint32 i=0; i<AS_readpool.size(); i++) {
	AS_readhitmiss[i].resize(AS_readpool[i].getLenClippedSeq(),0);
	AS_readhmcovered[i].resize(AS_readpool[i].getLenClippedSeq(),0);
      }

#if TRACKMEMUSAGE
  cout << "\ndmi ma 20\n";
  dumpMemInfo();
#endif

    }
  }else{
    CEBUG("Don't want cpv\n");
  }

  // to find troublemakers
  // do not clear AS_istroublemaker! (carry over inbetween each call,
  //  once a troublemaker, always a troublemaker)
  // TODO: check whether this is ok
  AS_istroublemaker.resize(AS_readpool.size(),0);

  AS_allrmbsok.resize(AS_readpool.size(),0);
  AS_weakrmbsnotok.resize(AS_readpool.size(),0);
  AS_probablermbsnotok.resize(AS_readpool.size(),0);

  std::ofstream ffout;
  std::ofstream cfout;
  std::ofstream rfout;
  if(AS_logflag_adsdump){
    // forward, complement, reject
    ffout.open(adsf_fn, std::ios::out);
    cfout.open(adsc_fn, std::ios::out);
    rfout.open(adsr_fn, std::ios::out);
  }

  AS_CUMADSLofstream.open(adsfacts_fn, std::ios::out|std::ios::ate);
  AS_numADSFacts_fromalignments=AS_numADSFacts_fromshreds;

#if TRACKMEMUSAGE
  cout << "\ndmi ma 30\n";
  dumpMemInfo();
#endif

  try{
    {
      if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      cout << "\nAligning possible forward matches:\n";

      std::string ourskimfilename=AS_posfmatch_filename;
      if(takefullskimfilenames) ourskimfilename=AS_posfmatch_full_filename;

      makeAlignmentsFromPosMatchFile(
	ourskimfilename.c_str(),
	version,
	1,
	trans100percent,
	checkfunction,
	ffout,
	rfout);
    }

#if TRACKMEMUSAGE
  cout << "\ndmi ma 40\n";
  dumpMemInfo();
#endif

    //dump454OvercallsArrays();
    //exit(100);

    {
      if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      cout << "\nAligning possible complement matches:\n";

      std::string ourskimfilename=AS_poscmatch_filename;
      if(takefullskimfilenames) ourskimfilename=AS_poscmatch_full_filename;

      makeAlignmentsFromPosMatchFile(
	ourskimfilename.c_str(),
	version,
	-1,
	trans100percent,
	checkfunction,
	cfout,
	rfout);
    }

#if TRACKMEMUSAGE
    cout << "\ndmi ma 50\n";
    dumpMemInfo();
#endif

    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

    // calculate possible vectorclips and store them,
    //  release big chunk of memory that is otherwise wasted
    if(AS_readhitmiss.size()) {
      calcPossibleSeqVectorClipoffs(version, "", "_pass");
    }
    nukeSTLContainer(AS_readhitmiss);
    nukeSTLContainer(AS_readhmcovered);
    nukeSTLContainer(AS_count_rhm);

#if TRACKMEMUSAGE
    cout << "\ndmi ma 60\n";
    dumpMemInfo();
#endif



    // Find troublemakers
    // Later on, when hits are loaded again, reduce the weight of overlaps
    //  of troublemakers that are in in the overlap graph
    // Troublemakers are reads that have more overlaps with mismatches
    //  in WRMr tagged positions with other reads than overlaps with
    //  no problems.
    // huntSpoilSports also tags reads at contig ends as troublemakers

    //cout << "WRMB statistics\n";

    for(uint32 i=0; i<AS_allrmbsok.size(); i++){
      //cout << AS_allrmbsok[i] << "\t" << AS_weakrmbsnotok[i] << "\t" << AS_probablermbsnotok[i];
      if(AS_weakrmbsnotok[i] > AS_allrmbsok[i]
    	 || AS_probablermbsnotok[i] > AS_allrmbsok[i]) {
    	// ok, it's a troublemaker
    	AS_istroublemaker[i]=1;
      }
    }

    nukeSTLContainer(AS_allrmbsok);
    nukeSTLContainer(AS_weakrmbsnotok);
    nukeSTLContainer(AS_probablermbsnotok);
    // do not clear AS_istroublemaker!

#if TRACKMEMUSAGE
  cout << "\ndmi ma 70\n";
  dumpMemInfo();
#endif

    cout << '\n';
#if 0
    {
      for(uint32 i=0; i<AS_readpool.size(); i++){
	cout << i<< ":" << AS_readpool[i].getName();
	cout << '\n';
	uint32 counts=AS_edges_forward.count(i);
	cout << i << " is " << counts << " times as key in the forward_edge_mmap.\n";
	//    cout << " \t" << "forward: " << map.count(i);
	cout << '\n';
	cout.flush();
      }
    }
#endif

#if 0
    // output of RHM vector
    cout << "RHM vectors:\n";
    for(uint32 i=0; i<AS_readhitmiss.size(); i++) {
      cout << "Id: " << i << "\t" << AS_readpool[i].getName();
      cout << "\t# ADS: " << AS_count_rhm[i] << endl;
      for(uint32 j=0; j<AS_readhitmiss[i].size(); j++) {
	cout << AS_readhitmiss[i][j] << " ";
      }
      cout << endl;
      for(uint32 j=0; j<AS_readhmcovered[i].size(); j++) {
	cout << AS_readhmcovered[i][j] << " ";
      }
      cout << endl;
    }

#endif

    AS_CUMADSLofstream.close();
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  AS_steps[ASADSLISTOK]=1;

  FUNCEND();

  //exit(0);
}

//#define CEBUG(bla)
//#define CEBUGF(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::priv_loadAlignmentsFromFile(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)
{
  FUNCSTART("void Assembly::loadAlignmentsFromFile(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  std::string adsfacts_fn;
  if(tmpfname.size()){
    adsfacts_fn=buildFileName(version, prefix, postfix, tmpfname, ".adsfacts");
  }else{
    adsfacts_fn=buildFileName(version, prefix, postfix,
			      as_fixparams.as_tmpf_ads,
			      ".adsfacts");
  }

  std::string fnpovl(buildDefaultCheckpointFileName(AS_miraparams[0].getFileParams().chkpt_persistentoverlaps));
  if(!fileExists(fnpovl)){
    fnpovl.clear();
  }

  try{
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    cout << "Counting number of alignments in files ...";
    cout .flush();
//TODO add ProgressIndic
    uint64 totaladsfacts=countLinesInFile(adsfacts_fn);
    if(!fnpovl.empty()){
      totaladsfacts+=countLinesInFile(fnpovl);
    }
    cout << " done.\nExpecting " << totaladsfacts << " alignments.\n";


    if(as_fixparams.as_dateoutput) dateStamp(cout);
    cout << "Loading confirmed " << totaladsfacts << " overlaps from disk (will need approximately ";
    byteToHumanReadableSize(
      static_cast<double>(sizeof(AlignedDualSeqFacts))*totaladsfacts
      +static_cast<double>(sizeof(newedges_t))*totaladsfacts*2,
      cout );
    cout << " RAM):" << endl;

    // Reduce memory fragmentation, keep these two buggers all time
    //  reserved in memory (do NOT nuke under normal circumstances!)
    AS_adsfacts.clear();
    AS_confirmed_edges.clear();

#ifdef MIRA_ADSSTORE_CAPACITY
    if(AS_adsfacts.capacity() < totaladsfacts){
      // if current capacity is not enough, nuke and reserve with
      //  15% additional capacity
      // Idea is that in the first pass we will get a number of ADS that
      //  should be pretty close to the number we get in subsequent passes.
      // Therefore if we stay within the capacity of the vectors, no
      //  new block need to be allocated

      //nukeSTLContainer(AS_adsfacts);
      //nukeSTLContainer(AS_confirmed_edges);

#if TRACKMEMUSAGE
      cout << "\n\n\nOMG OMG OMG ... we must reserve anew!\n\n\n";
      cout << "\ndmi laff  omg 00\n";
      dumpMemInfo();
#endif

      AS_adsfacts.reserve(totaladsfacts + (totaladsfacts/100*15));

#if TRACKMEMUSAGE
      cout << "\ndmi laff  omg 10\n";
      dumpMemInfo();
#endif

      AS_confirmed_edges.reserve(totaladsfacts*2 + (totaladsfacts*2/100*15));

#if TRACKMEMUSAGE
      cout << "\ndmi laff  omg 20\n";
      dumpMemInfo();
#endif
    }
#endif //MIRA_ADSSTORE_CAPACITY

#if TRACKMEMUSAGE
    cout << "\ndmi laff 00\n";
    dumpMemInfo();
#endif

    AS_adsfacts.resize(totaladsfacts);
    AS_confirmed_edges.resize(totaladsfacts*2);

#if TRACKMEMUSAGE
  cout << "\ndmi laff 10\n";
  dumpMemInfo();
#endif

    // TODO: *sigh* try to save even more space .... somehow?

    uint64 runningADSFactnumber=0;
    uint64 numbannedADSFacts=0;

    ProgressIndicator<int64> P(0, AS_adsfacts.size());
    priv_laffhelper(adsfacts_fn,
		    P,runningADSFactnumber,numbannedADSFacts);
    if(!fnpovl.empty()){
      priv_laffhelper(fnpovl,
		      P,runningADSFactnumber,numbannedADSFacts);
    }
    P.finishAtOnce();

    if(runningADSFactnumber+numbannedADSFacts != AS_adsfacts.size()) {
      MIRANOTIFY(Notify::INTERNAL, "Error while loading adsfacts, less facts in file than calculated earlier:" << AS_adsfacts.size() << " vs. " << runningADSFactnumber << " + " << numbannedADSFacts);
    }
    if(numbannedADSFacts){
      cout << "Found " << numbannedADSFacts << " banned facts.\n";
    }
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    if(runningADSFactnumber < AS_adsfacts.size() ){
      cout << "Resizing pool to " << runningADSFactnumber << " overlaps.\n";
      if(runningADSFactnumber){
	AS_adsfacts.resize(runningADSFactnumber);
	AS_confirmed_edges.resize((runningADSFactnumber-1)*2);
      }else{
	AS_adsfacts.clear();
	AS_confirmed_edges.clear();
      }
    }
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  // Apply malus to overlaps we do not want to be taken early
  for(auto & cee : AS_confirmed_edges){
    uint32 malus=getOverlapMalusDivider(cee.rid1, cee.linked_with);
//      cout << "Malus\t" << AS_readpool[cee.rid1].getName()
//	   << "\t" << AS_readpool[cee.linked_with].getName()
//	   << "\t" << malus << endl;
    if(malus>1){
      cee.best_weight/=malus;
    }
  }

  // Make sure there's no overlap with a best_weight of "0"
  // (pathfinder is not prepared for this)
  for(auto & cee : AS_confirmed_edges){
    if(cee.best_weight==0) cee.best_weight=1;
  }

  // Sort overlaps

  cout << "\n\nSorting confirmed overlaps (this may take a while) ... ";
  cout.flush();
  if(AS_confirmed_edges.size()<=65535){
    mstd::ssort(AS_confirmed_edges, newedges_t::sortComparatorByRIDUpByWeightDown);
  }else{
    mstd::psort(AS_confirmed_edges, newedges_t::sortComparatorByRIDUpByWeightDown);
  }
  cout << "done.\n" << endl;

  if(as_fixparams.as_dateoutput) dateStamp(cout);

#if 0
  for(const auto & cee : AS_confirmed_edges){
    cout << AS_readpool[cee.rid1].getName() << '\t' << AS_readpool[cee.linked_with].getName() << '\t' << cee;
  }
#endif


  // output of possible clusters

  if(as_fixparams.as_tmpf_ads.size()!=0){
    cout << "Generating clusters:\n";

    std::string pclusters_fn;
    if(tmpfname.size()){
      pclusters_fn=buildFileName(version, prefix, postfix, tmpfname, ".adsfacts.pclusters");
    }else{
      pclusters_fn=buildFileName(version, prefix, postfix,
				as_fixparams.as_tmpf_ads,
				".adsfacts.pclusters");
    }

    std::ofstream pfout;
    pfout.open(pclusters_fn, std::ios::out);

#if 0
#else
    // fast, but memory intensive
    std::vector<int32> clusteridperread;
    std::vector<std::list<int32> > readinclusterlist;
    std::vector<int8> dummy; // empty, so that all reads get clustered

    clusterUnassembledReads(clusteridperread,readinclusterlist, dummy);

    {
      // sort is not needed, but might be easier to understand what is
      //  happening when looking at the file. And it's quick enough.

      cout << "\nSorting clusters:\n";

      size_t riclsize=readinclusterlist.size();
      ProgressIndicator<int32> P(0,riclsize);
      for(size_t i=0; i<riclsize; i++){
	P.increaseprogress(1);
	if(!readinclusterlist[i].empty()){
	  // these are lists, use list sorting
	  readinclusterlist[i].sort();
	}
      }
      P.finishAtOnce();
    }

    if(as_fixparams.as_dateoutput) dateStamp(cout);
    {
      cout << "\nWriting clusters:\n";

      uint32 clustercount=0;
      size_t riclsize=readinclusterlist.size();
      ProgressIndicator<int32> P(0,riclsize);
      for(size_t i=0; i<riclsize; ++i){
	//cout << "cc "<<i << '\t' << clustercount;
	P.increaseprogress(1);
	if(!readinclusterlist[i].empty()){
	  for(const auto & ricle : readinclusterlist[i]){
	    pfout << clustercount << '\t' << AS_readpool[ricle].getName() << '\n';
	  }
	  ++clustercount;
	}
      }
      P.finishAtOnce();
    }

#endif

    // write reads that did not cluster
    {
      cout << "\nWriting unclustered reads ... ";
      cout.flush();
      for(size_t j=0; j<AS_readpool.size(); j++) {
	if(clusteridperread[j]==-1) {
	  pfout << "-1\t" << AS_readpool[j].getName() << '\n';
	}
      }
      cout << "done.\n";
    }
    pfout.close();
  }

  if(AS_logflag_loadedoverlaps){
    std::string elfn(buildFileName(version, prefix, postfix,
			      "elog.overlapsinmem",
			      ".txt"));
    cout << "Logging overlaps in mem to " << elfn << endl;
    std::ofstream eout(elfn, std::ios::out);

    for(auto & cee : AS_confirmed_edges){
      eout << AS_readpool[cee.rid1].getName()
	   << "\t" << AS_readpool[cee.linked_with].getName()
	   << "\t" << cee
	   << "\n";
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::priv_laffhelper(const std::string & fn, ProgressIndicator<int64> & P, uint64 & runningADSFactnumber, uint64 & numbannedADSFacts)
{
  FUNCSTART("void Assembly::priv_laffhelper(const std::string & fn, ProgressIndicator<int64> & P, uint64 & runningADSFactnumber, uint64 & numbannedADSFacts)");

  std::ifstream finfin(fn, std::ios::in);
  if(!finfin){
    MIRANOTIFY(Notify::FATAL, "File not found? MIRA read it a few moments ago, it MUST exist: " << fn);
  }

  while(!finfin.eof()){
    uint32 bestweight;
    int16 direction;
    finfin >> bestweight;
    if(finfin.eof()) break;
    finfin >> direction;
    if(runningADSFactnumber >= AS_adsfacts.size()) {
      MIRANOTIFY(Notify::INTERNAL, "Error while loading adsfacts, more facts in file than calculated earlier. Calc: " << AS_adsfacts.size() << "\tNow at " << runningADSFactnumber);
    }
    bool flag_stronggood;
    bool flag_weakgood;
    bool flag_belowavgfreq;
    bool flag_norept;
    bool flag_rept;

    finfin >> flag_stronggood;
    finfin >> flag_weakgood;
    finfin >> flag_belowavgfreq;
    finfin >> flag_norept;
    finfin >> flag_rept;

    AS_adsfacts[runningADSFactnumber].serialiseIn(finfin);

    auto rid1=AS_adsfacts[runningADSFactnumber].getID1();
    auto rid2=AS_adsfacts[runningADSFactnumber].getID2();
    // reduce the hit weights of troublemakers
    // AND
    // set the score ratio to below the one for pathfinder
    //  quickrules
    if(AS_istroublemaker[rid1]
       || AS_istroublemaker[rid2]){
      if(bestweight>=100){
	bestweight/=100;
      }else if(bestweight>=10){
	bestweight/=10;
      }else{
	bestweight/=2;
      }
      uint8 st=AS_readpool[rid1].getSequencingType();
      int8 minsr=AS_miraparams[st].getPathfinderParams().paf_quickrule_minsim1;
      minsr=std::min(minsr,AS_miraparams[st].getPathfinderParams().paf_quickrule_minsim2);
      st=AS_readpool[rid2].getSequencingType();
      minsr=std::min(minsr,AS_miraparams[st].getPathfinderParams().paf_quickrule_minsim1);
      minsr=std::min(minsr,AS_miraparams[st].getPathfinderParams().paf_quickrule_minsim2);
      if(minsr>0) --minsr;
      AS_adsfacts[runningADSFactnumber].setScoreRatio(minsr);
    }

    // insert only ADSFacts where the reads are not permanently banned
    //  from overlapping
    // BaCh 14.11.2014: OR, this is new, check here whether that read might have fallen into
    //  disgrace. E.g., a chimera search in a later pass kicking away a read which has
    //  saved short overlaps
    // Let's hope the isUsedInAssembly() flag of a read tells the truth here and we did not
    //  forget to set it right.
    bool banned=!(AS_readpool[rid1].isUsedInAssembly() & AS_readpool[rid2].isUsedInAssembly());
    // well, and as I do not trust myself: check if one of the reads has length 0 ...
    //   ... a clear indicator that it has been sorted out earlier but not flagged as such
    // at the latest the contig would moan if it had to add a read of length 0
    if(likely(!banned)){
      if(AS_readpool[rid1].getLenClippedSeq()==0){
	banned=true;
	AS_readpool[rid1].setUsedInAssembly(false);
      }
      if(AS_readpool[rid2].getLenClippedSeq()==0){
	banned=true;
	AS_readpool[rid2].setUsedInAssembly(false);
      }
    }
    if(likely(!banned)){
      banned|=AS_permanent_overlap_bans.checkIfBanned(rid1,rid2);
    }
    if(likely(!banned)){
      AS_confirmed_edges[runningADSFactnumber*2].rid1=rid1;
      AS_confirmed_edges[runningADSFactnumber*2].linked_with=rid2;
      AS_confirmed_edges[runningADSFactnumber*2].best_weight=bestweight;
      AS_confirmed_edges[runningADSFactnumber*2].adsfindex=runningADSFactnumber;
      AS_confirmed_edges[runningADSFactnumber*2].direction=direction;
      AS_confirmed_edges[runningADSFactnumber*2].ol_stronggood=flag_stronggood;
      AS_confirmed_edges[runningADSFactnumber*2].ol_weakgood=flag_weakgood;
      AS_confirmed_edges[runningADSFactnumber*2].ol_belowavgfreq=flag_belowavgfreq;
      AS_confirmed_edges[runningADSFactnumber*2].ol_norept=flag_norept;
      AS_confirmed_edges[runningADSFactnumber*2].ol_rept=flag_rept;

      AS_confirmed_edges[runningADSFactnumber*2+1]=AS_confirmed_edges[runningADSFactnumber*2];
      AS_confirmed_edges[runningADSFactnumber*2+1].rid1=rid2;
      AS_confirmed_edges[runningADSFactnumber*2+1].linked_with=rid1;

      ++runningADSFactnumber;
    }else{
      ++numbannedADSFacts;
    }
    P.increaseprogress();
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

uint32 Assembly::getOverlapMalusDivider(int32 id1, int32 id2)
{
  FUNCSTART("uint32 Assembly::getOverlapMalusDivider(int32 id1, int32 id2)");

  if(AS_overlapcritlevelvl.empty()) return 1;

  BUGIFTHROW(id1>=AS_overlapcritlevelvl[1].size(),"id1>=AS_overlapcritlevell.size() ?");
  BUGIFTHROW(id2>=AS_overlapcritlevelvl[1].size(),"id2>=AS_overlapcritlevell.size() ?");

  auto ocll1=std::min(AS_overlapcritlevelvl[0][id1],AS_overlapcritlevelvl[1][id1]);
  auto oclr1=std::min(AS_overlapcritlevelvr[0][id1],AS_overlapcritlevelvr[1][id1]);
  auto ocll2=std::min(AS_overlapcritlevelvl[0][id2],AS_overlapcritlevelvl[1][id2]);
  auto oclr2=std::min(AS_overlapcritlevelvr[0][id2],AS_overlapcritlevelvr[1][id2]);

  if(ocll1==0
     && oclr1==0
     && ocll2==0
     && oclr2==0) return 1;

  if(ocll1==255
     || oclr1==255
     || ocll2==255
     || oclr2==255) return 1000;

  if(ocll1==240
     || oclr1==240
     || ocll2==240
     || oclr2==240) return 100;

  FUNCEND();
  return 10;
}

/*************************************************************************
 *
 * if usedids is empty: cluster all reads
 * if not empty: cluster all reads which are marked unused (0) in usedids
 *
 * returns:
 *  - in clusteridperread the clister-id for each read (-1 if
 *    unclustered)
 *  - readinclusterlist vector of lists: for each cluster, the read-ids
 *     contained within. Also contains real singlets, but not orphans (see
 *     below).
 *    Warning: if vector usedids was used, then orphans (singlets with
 *      with overlaps only to reads already used) will not appear as
 *      own cluster
 *    Warning: will contain empty clusters with no associated reads.
 *      Getting them out would mean recalc, and is not needed if callers
 *      know that this may happen.
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::clusterUnassembledReads(std::vector<int32> & clusteridperread, std::vector<std::list<int32> > & readinclusterlist, const std::vector<int8> & usedids)
{
  clusteridperread.clear();
  clusteridperread.resize(AS_readpool.size(),-1);
  readinclusterlist.clear();
  readinclusterlist.reserve(AS_readpool.size()/16);

  uint32 clustercount=0;

  //
  //cout << "CUR: " << usedids.size() << endl;
  //for(uint32 uid=0; uid<usedids.size(); ++uid){
  //  cout << "uid: " << uid << "\t" << static_cast<uint16>(usedids[uid])<<endl;
  //}
  //cout.flush();
  //cout << "AS_confirmed_edges.size(): " << AS_confirmed_edges.size() << endl;

  {
    ProgressIndicator<int32> P(0,
			       static_cast<int32>(AS_confirmed_edges.end()-AS_confirmed_edges.begin()));
    auto I=AS_confirmed_edges.cbegin();
    for(;I != AS_confirmed_edges.cend(); I++) {
      if(!usedids.empty() && (usedids[I->rid1] || usedids[I->linked_with])) continue;
      int32 cnum1=clusteridperread[I->rid1];
      int32 cnum2=clusteridperread[I->linked_with];
      //cout << "link: " << I->rid1 << " (" << cnum1 << ")\t" << I->linked_with << " (" << cnum2 << ")\t";
      if(cnum1==-1 && cnum2==-1) {
	//cout << "new cluster: " << clustercount;
	clusteridperread[I->rid1]=clustercount;
	clusteridperread[I->linked_with]=clustercount;
	readinclusterlist.resize(clustercount+1);
	readinclusterlist[clustercount].push_back(I->rid1);
	readinclusterlist[clustercount].push_back(I->linked_with);
	clustercount++;
      } else if(cnum1==-1) {
	//cout << "link " << I->rid1 << "\tinto cluster " << cluster[I->linked_with];
	clusteridperread[I->rid1]=clusteridperread[I->linked_with];
	readinclusterlist[clusteridperread[I->linked_with]].push_back(I->rid1);
      } else if(cnum2==-1) {
	//cout << "link " << I->linked_with << "\tinto cluster " << clusteridperread[I->rid1];
	clusteridperread[I->linked_with]=clusteridperread[I->rid1];
	readinclusterlist[clusteridperread[I->rid1]].push_back(I->linked_with);
      } else {
	if (cnum1 != cnum2) {
	  // uh oh ... we have to merge both these clusters

	  int32 killed=std::max(cnum1,cnum2);
	  int32 survive=std::min(cnum1,cnum2);

	  //cout << "merge. Kill cluster " << killed << "\tsurvive " << survive;
	  //
	  //cout << "\nSize(killed): " << readinclusterlist[killed].size();
	  //cout << "\nSize(survive): " << readinclusterlist[survive].size();
	  //cout << "\nMoving to cluster " << survive << " the reads:";
	  for(const auto & ricle : readinclusterlist[killed]){
	    clusteridperread[ricle]=survive;
	  }
	  readinclusterlist[survive].splice(readinclusterlist[survive].end(),
					    readinclusterlist[killed]);
	  //cout << "\nSize(killed): " << readinclusterlist[killed].size();
	  //cout << "\nSize(survive): " << readinclusterlist[survive].size();
	}
      }
      //cout << '\n';
      P.increaseprogress(1);
    }
    P.finishAtOnce();
  }

//  {
//    uint32 numclu=0;
//    for(size_t ricli=0; ricli<readinclusterlist.size(); ricli++){
//      if(!readinclusterlist[ricli].empty()){
//	++numclu;
//      }
//    }
//    cout << "Seeing " << numclu << " clusters.\n";
//  }

  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::minimiseMADSL(std::list<AlignedDualSeq> & madsl)
{
  CEBUG("minimiseMADSL() 1 madsl.size(): " << madsl.size() << endl);

  int32 bestweight=0;
  for(auto adsI= madsl.cbegin(); adsI!=madsl.cend(); ){
    //CEBUG("MADSL entry:\n" << *adsI << endl);
    //cout << "MADSL entry:\n" << *adsI << endl;
    if(adsI->isValid()==false){
      adsI=madsl.erase(adsI);
    }else{
      if(adsI->getWeight()>bestweight) bestweight=adsI->getWeight();
      ++adsI;
    }
  }
  CEBUG("minimseMADSL() 2 madsl.size(): " << madsl.size() << endl);

  // Save Memory, take only the best
  for(auto adsI= madsl.cbegin(); adsI!=madsl.cend(); ){
    if(adsI->getWeight() != bestweight){
      adsI=madsl.erase(adsI);
    } else {
      ++adsI;
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

//lpwcfs1_9079
//#define CEBUG(bla)   {if(rid1==9078 || rid2==9078) cout << bla; cout.flush(); }

void Assembly::cleanupMADSL(std::list<AlignedDualSeq> & madsl, uint32 rid1, uint32 rid2, int8 direction, bool flag_stronggood, bool flag_weakgood, bool flag_belowavgfreq, bool flag_norept, bool flag_rept)
{
  FUNCSTART("void Assembly::cleanupMADSL(std::list<AlignedDualSeq> & madsl, uint32 rid1, uint32 rid2, int8 direction, bool flag_stronggood, bool flag_weakgood, bool flag_belowavgfreq, bool flag_norept, bool flag_rept)");

  CEBUG("calling minimise\n");
  minimiseMADSL(madsl);

  CEBUG("cleanupMADSL() 3 madsl.size(): " << madsl.size() << endl);

  if(madsl.size()==0){
    // put both ids in permanent overlap banlist so that they
    //  won't make it through skim the next pass
    CEBUG("Banning " << AS_readpool[rid1].getName() << " & " << AS_readpool[rid2].getName() << endl);
    AS_permanent_overlap_bans.insertBan(rid1,rid2);
  }else{
    //cout << madsl.front();
    //CEBUG(madsl.front());

    // From here on, we'll take only the first solution anyway!

    // ok, check for mismatches in RMB regions
    if(AS_readpool.getRead(rid1).getNumOfTags() >0
       || AS_readpool.getRead(rid2).getNumOfTags() >0) {

      int32 adscheck=checkADSForRepeatMismatches(*madsl.begin());
      if (adscheck<0) {

	// ban all except template partners
	if(static_cast<uint32>(AS_readpool.getRead(rid1).getTemplatePartnerID())!=rid2) {
	  CEBUG("ADS RMB banning: " << rid1 << ": " << AS_readpool.getRead(rid1).getName() << "\t" << rid2 << ": " << AS_readpool.getRead(rid2).getName() << endl);
	  //cout << "ADS RMB banning: " << rid1 << ": " << AS_readpool.getRead(rid1).getName() << "\t" << rid2 << ": " << AS_readpool.getRead(rid2).getName() << endl;

	  // put both ids in permanent overlap banlist
	  AS_permanent_overlap_bans.insertBan(rid1,rid2);

	  AS_probablermbsnotok[rid1]++;
	  AS_probablermbsnotok[rid2]++;
	  // and get out of here
	  // TODO: modify caller to log these into ads_reject_prmb file
	  return;
	}
      }else if(adscheck>0) {
	AS_weakrmbsnotok[rid1]++;
	AS_weakrmbsnotok[rid2]++;
      }else{
	AS_allrmbsok[rid1]++;
	AS_allrmbsok[rid2]++;
      }
    }


    if(AS_miraparams[0].getAssemblyParams().as_clip_possible_vectors
       && ! AS_steps[ASVECTORSCLIPPED]) {
      CEBUG("transcribing hits\n");
      transcribeHits(*madsl.begin());
      CEBUG("th done\n");
    }

    // if changing something here, do not forget to change at the 100%
    //  trans place too
    AS_CUMADSLofstream << madsl.begin()->getWeight() << '\t'
		       << static_cast<int16>(direction) << '\t'
		       << flag_stronggood << '\t'
		       << flag_weakgood << '\t'
		       << flag_belowavgfreq << '\t'
		       << flag_norept << '\t'
		       << flag_rept << '\t';
    madsl.begin()->serialiseOut(AS_CUMADSLofstream);
    AS_CUMADSLofstream << '\n';
    if(AS_CUMADSLofstream.bad()){
      MIRANOTIFY(Notify::FATAL, "Could not write anymore to disk (at SWcomp). Disk full? Changed permissions?");
    }
    AS_numADSFacts_fromalignments++;
  }

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Checks an ADS for mismatches occurring at RepeatMarkerBase positions
 *
 * Function does not care about quality, those tags are set by the contig
 *  and are therefore not to be questioned
 *
 * Gives back -1 if a mismatch occurs at a SRMr or CRMr tag
 * Gives back 1 for WRMr/WRMr mismatches
 * Else gives back 0 for no problem
 *
 * also gives back 1 for SIOr against SRMr/WRMr mismatches if the
 *  strains are the same and we are not in "assumeSNP" mode
 *
 * Sets CRMr (Carbon-copy Repeat Marker) for reads that have no SRMr
 *  tag where the other read has on IF the strains are the same
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }

//#define CEBUG(bla)   {if(AS_readpool.getRead(id1).getName()=="368544_2432_0180" || AS_readpool.getRead(id2).getName()=="368544_2432_0180") cout << bla; cout.flush(); }

//#define CEBUG(bla)   {if(id1==9078 || id2==9078) cout << bla; cout.flush(); }

// due to new deferred tagging in checkADSForRepeatMismatches(), a "no problem"
//  result at the first call may not be correct if some CRMr tags were set
// In that case, we need a second call to be certain
int32 Assembly::checkADSForRepeatMismatches(AlignedDualSeq & ads)
{
  bool need2ndpass=false;
  int32 retval=checkADSForRepeatMismatches_wrapped(ads,need2ndpass);
  if(retval==0 && need2ndpass) retval=checkADSForRepeatMismatches_wrapped(ads,need2ndpass);
  return retval;
}

int32 Assembly::checkADSForRepeatMismatches_wrapped(AlignedDualSeq & ads, bool & need2ndpass)
{
  FUNCSTART("void Assembly::checkADSForRepeatMismatches(AlignedDualSeq & ads)");

  need2ndpass=false;

  int32 id1=ads.getID1();
  int32 id2=ads.getID2();

  CEBUG("ADSfrm: " << AS_readpool.getRead(id1).getName() << " (" << id1 << ") " << AS_readpool.getRead(id2).getName() << " (" << id2 << ")" << endl << endl);

  CEBUG(AS_readpool.getRead(id1) << endl << endl);
  CEBUG(AS_readpool.getRead(id2) << endl << endl);
  CEBUG(ads);

  int32 r1index, r2index, r1delta, r2delta;

  if(ads.getSequenceDirection(id1) > 0) {
    r1index=AS_readpool.getRead(id1).getLeftClipoff()-ads.getOffsetInAlignment(id1);
    r1delta=1;
  } else {
    CEBUG("-1: " << AS_readpool.getRead(id1).getRightClipoff() << " " << ads.getOffsetInAlignment(id1) << endl);
    r1index=AS_readpool.getRead(id1).getRightClipoff()-1+ads.getOffsetInAlignment(id1);
    r1delta=-1;
  }

  if(ads.getSequenceDirection(id2) > 0) {
    r2index=AS_readpool.getRead(id2).getLeftClipoff()-ads.getOffsetInAlignment(id2);
    r2delta=1;
  } else {
    CEBUG("-2: " << AS_readpool.getRead(id1).getRightClipoff() << " " << ads.getOffsetInAlignment(id1) << endl);
    r2index=AS_readpool.getRead(id2).getRightClipoff()-1+ads.getOffsetInAlignment(id2);
    r2delta=-1;
  }

  CEBUG("1: r1index: " << r1index << "\tr2index: " << r2index << endl);

  char const * cons=ads.getGapedConsensusSequence();
  char const * alseq1=ads.getAlignedSequence(id1);
  char const * alseq2=ads.getAlignedSequence(id2);

  bool hasWRMrWRMrmismatch=false;
  bool hasstrongRMBmismatch=false;
  bool hasweakRMBmismatch=false;
  bool penaliseSIOrmismatch=true;

  size_t numfrequencymismatch=0;

  // no SIO checks when we assume SNPs instead of RMBs
  // no SIO checks when strains are different
  if(AS_miraparams[0].getContigParams().con_assume_snp_insteadof_rmb
     || AS_readpool.getRead(id1).getStrainID() != AS_readpool.getRead(id2).getStrainID()) penaliseSIOrmismatch=false;

  CEBUG("Check SIOr? " << penaliseSIOrmismatch << endl);

  // advance the pointers and index to just check the overlapping part
  if(ads.getOffsetInAlignment(id1)==0){
      cons+=ads.getOffsetInAlignment(id2);
      alseq1+=ads.getOffsetInAlignment(id2);

      r1index+=ads.getOffsetInAlignment(id2)*r1delta;
      r2index+=ads.getOffsetInAlignment(id2)*r2delta;
  }else {
    cons+=ads.getOffsetInAlignment(id1);
    alseq2+=ads.getOffsetInAlignment(id1);

    r1index+=ads.getOffsetInAlignment(id1)*r1delta;
    r2index+=ads.getOffsetInAlignment(id1)*r2delta;
  }
  CEBUG("2: r1index: " << r1index << "\tr2index: " << r2index << endl);
  CEBUG("r1delta: " << r1delta << "\tr2delta: " << r2delta << endl);
  // Internal compiler error fur gcc2.95x! cout << "count: " << count << endl;

  // to enable this routine to set Carbon copy Marker tags which are wider than 1 base,
  //  need to use "deferred tagging" using these variables
  int32 r1newctagfrom=-1;
  int32 r1newctagto=-1;
  int32 r2newctagfrom=-1;
  int32 r2newctagto=-1;

  bool problemoccurred=false;

  // keep count in this loop forward, is used further downward!
  for(uint32 count=0; count<ads.getOverlapLen(); ++count, ++cons, r1index+=r1delta, r2index+=r2delta, ++alseq1, ++alseq2){

    if(*alseq1=='*') {
      r1index-=r1delta;
    }
    if(*alseq2=='*') {
      r2index-=r2delta;
    }

    if(r1index<0) {
      cout << "Megaproblem: r1index<0" << endl;
      problemoccurred=true;
    }
    if(r1index >= static_cast<int32>(AS_readpool.getRead(id1).getLenSeq())){
      cout << "Megaproblem: r1index >= AS_readpool.getRead(id1).getLenSeq()" << endl;
      problemoccurred=true;
    }
    if(r2index<0) {
      cout << "Megaproblem: r2index<0" << endl;
      problemoccurred=true;
    }
    if(r2index >= static_cast<int32>(AS_readpool.getRead(id2).getLenSeq())){
      cout << "Megaproblem: r2index >= AS_readpool.getRead(id2).getLenSeq()" << endl;
      problemoccurred=true;
    }
    if(problemoccurred){
      cout << "ADSfrm: (" << id1 << ") " << AS_readpool.getRead(id1).getName() << " (" << id2 << ") " << AS_readpool.getRead(id2).getName() << endl << endl;

      Read::setCoutType(Read::AS_CAF);
      cout << AS_readpool.getRead(id1) << endl << endl;
      cout << AS_readpool.getRead(id2) << endl << endl;
      cout << ads;
    }

    CEBUG(static_cast<char>(*cons) << " r1index: " << r1index);
    CEBUG("\tr2index: " << r2index);

    BUGIFTHROW(r1index < 0, "r1index < 0 ???");
    BUGIFTHROW(r1index >= static_cast<int32>(AS_readpool.getRead(id1).getLenSeq()), "r1index >= length of read ???");
    BUGIFTHROW(r2index < 0, "r2index < 0 ???");
    BUGIFTHROW(r2index >= static_cast<int32>(AS_readpool.getRead(id2).getLenSeq()), "r2index >= length of read ???");

    CEBUG("\t" << AS_readpool.getRead(id1).getBaseInSequence(r1index) << " " << AS_readpool.getRead(id2).getBaseInSequence(r2index));
    if(toupper(AS_readpool.getRead(id1).getBaseInSequence(r1index)) != toupper(AS_readpool.getRead(id2).getBaseInSequence(r2index))) {
      CEBUG(" !");
    }
    CEBUG("\t" << static_cast<char>(*alseq1) << " " << static_cast<char>(*alseq2));
    CEBUG("\t" << count);
    //if(contit) {
    //  //CEBUG("\tStarbase r1index & r2index not valid!" << endl);
    //  //continue;
    //}


    // now the real check whether the bases differ.
    // This gets a bit more difficult than initially thought,
    //  as the alignment objects convert gaps into "N" bases

    // FIXME: now alignments keep gaps, rethink this!
    //  So until this is fixed (never?), we need to check
    //  the complicated way

    // First, none of the bases may be a "X"

    if(!(toupper(*alseq1)=='X'
	 // TODO: check whether take N back into consideration
//	 || toupper(*alseq1)=='N'
//	 || toupper(*alseq2)=='N'
	 || toupper(*alseq2)=='X')){

      // if needed, carbon copy existing SRMr tags at this place

      bool tagr1p=AS_readpool.getRead(id1).hasTag(Read::REA_tagentry_idSRMr, r1index);
      bool tagr2p=AS_readpool.getRead(id2).hasTag(Read::REA_tagentry_idSRMr, r2index);
      bool tagr1pc=AS_readpool.getRead(id1).hasTag(Read::REA_tagentry_idCRMr, r1index);
      bool tagr2pc=AS_readpool.getRead(id2).hasTag(Read::REA_tagentry_idCRMr, r2index);

      // if both reads are from same strain
      //  AND only one read has SRMr at this position (tagr1p XOR tagr2p)
      // then try to create CRMr tag
      bool mustctagr1=false;
      bool mustctagr2=false;
      if((AS_readpool.getRead(id1).getStrainID()
	  == AS_readpool.getRead(id2).getStrainID())
	 && (tagr1p  || tagr2p)
	 && ((tagr1p  && !tagr2p && !tagr2pc)
	     || (tagr2p  && !tagr1p && !tagr1pc))){

	// ok, let's check whether the environmenat is clean
	// i.e., whether the surrounding bases have no mismatch

//	bool isclean=false;
//	uint32 leftright=4;
//	// we need to stay within the overlap bounds
//	if(count > leftright && count + leftright < ads.getOverlapLen()){
//	  CEBUG("check area");
//	  isclean=true;
//	  // ok, check left and right whether everything is the same
//	  for(int32 checki=-leftright; checki <= static_cast<int32>(leftright); checki++){
//	    CEBUG(static_cast<char>(toupper(*(alseq1+checki))) << static_cast<char>(toupper(*(alseq2+checki))) << " ");
//	    if(checki==0) continue;
//	    if(toupper(*(alseq1+checki)) != toupper(*(alseq2+checki))){
//	      isclean=false;
//	      break;
//	    }
//	  }
//	  CEBUG("isclean: " << isclean);
//	  if(isclean) {
//	    if(tagr1p && dptools::isValidBase(*alseq2)) {
//	      AS_readpool.getRead(id2).addTag(r2index,r2index,
//					      Read::REA_tagentry_idCRMr,
//					      Read::REA_tagentry_coEmpty);
//	    }else if(tagr2p && dptools::isValidBase(*alseq1)){
//	      AS_readpool.getRead(id1).addTag(r1index,r1index,
//					      Read::REA_tagentry_idCRMr,
//					      Read::REA_tagentry_coEmpty);
//	    }
//	  }
//	}


	bool isclean=true;
	int32 leftright=4;
	int32 runvar=static_cast<int32>(count);
	runvar-=leftright;  // effectively -4
	// we need to stay within the overlap bounds
	for(int32 checki=-leftright; checki <= leftright; ++checki, ++runvar){
	  // no check at base (duh, normal as they differ) but also not in immediate surrounding
	  // reason: the multibasetags (length 2 to 3) set when gaps are involved.
	  if(checki >=-1 && checki <=1) continue;

	  // are we in overlap? (==in the non-clipped area of the reads)
	  if(runvar>=0 && runvar < ads.getOverlapLen()){
	    if(toupper(*(alseq1+checki)) != toupper(*(alseq2+checki))){
	      isclean=false;
	      break;
	    }
	  }
	}

	CEBUG("isclean: " << isclean);
	if(isclean) {
	  if(tagr1p && dptools::isValidBase(*alseq2)) {
	    mustctagr2=true;
	  }else if(tagr2p && dptools::isValidBase(*alseq1)){
	    mustctagr1=true;
	  }
	}

      }

      // now see whether we must set something in read 1
      if(mustctagr1 || r1newctagfrom>=0){
	if(mustctagr1){
	  if(r1newctagfrom>=0){
	    // extend tag
	    r1newctagto=r1index;
	  }else{
	    // new tag
	    r1newctagfrom=r1index;
	    r1newctagto=r1index;
	  }
	}else{
	  // deferred setting of tags
	  if(r1newctagfrom>r1newctagto) std::swap(r1newctagfrom,r1newctagto);
	  AS_tmptag_CRMr.from=r1newctagfrom;
	  AS_tmptag_CRMr.to=r1newctagto;
	  AS_readpool.getRead(id1).addTagO(AS_tmptag_CRMr);
	  need2ndpass=true;
	  r1newctagfrom=-1;
	  r1newctagto=-1;
	}
      }
      // now see whether we must set something in read 2
      if(mustctagr2 || r2newctagfrom>=0){
	if(mustctagr2){
	  if(r2newctagfrom>=0){
	    // extend tag
	    r2newctagto=r2index;
	  }else{
	    // new tag
	    r2newctagfrom=r2index;
	    r2newctagto=r2index;
	  }
	}else{
	  // deferred setting of tags
	  if(r2newctagfrom>r2newctagto) std::swap(r2newctagfrom,r2newctagto);
	  AS_tmptag_CRMr.from=r2newctagfrom;
	  AS_tmptag_CRMr.to=r2newctagto;
	  AS_readpool.getRead(id2).addTagO(AS_tmptag_CRMr);
	  need2ndpass=true;
	  r2newctagfrom=-1;
	  r2newctagto=-1;
	}
      }

      // simpler, both bases are some IUPAC base != N and != *
      // do  not check for "equal" bases, but for "contained"

      if(!dptools::areBasesContained(*alseq1,*alseq2)){
	//// BaCh 23.03.2009
	//// NEW and test: check on hash frequency
	//// ... and do this only for normal bases
	//// temp taken out 29.03.2009
	//if(dptools::isValidBase(*alseq1) && dptools::isValidBase(*alseq2)){
	//  uint8 freq1, freq2;
	//  if(r1delta>0){
	//    freq1=AS_readpool.getRead(id1).getBPosHashStats(r1index).fwd.getFrequency();
	//  }else{
	//    freq1=AS_readpool.getRead(id1).getBPosHashStats(r1index).rev.getFrequency();
	//  }
	//  if(r2delta>0){
	//    freq2=AS_readpool.getRead(id2).getBPosHashStats(r2index).fwd.getFrequency();
	//  }else{
	//    freq2=AS_readpool.getRead(id2).getBPosHashStats(r2index).rev.getFrequency();
	//  }
	//
	//  CEBUG("\tfreq1: " << static_cast<uint16>(freq1) << "\tfreq2: " << static_cast<uint16>(freq2));
	//  if(freq1>=3 && freq2>=3){
	//    numfrequencymismatch++;
	//    CEBUG("\tnumfrequencymismatch: " << numfrequencymismatch);
	//  }
	//}


	CEBUG("\t tagr?p:" << tagr1p << " " << tagr2p);

	//tagr1p|=AS_readpool.getRead(id1).hasTag(Read::REA_tagentry_idCRMr, r1index);
	//tagr2p|=AS_readpool.getRead(id2).hasTag(Read::REA_tagentry_idCRMr, r2index);
	tagr1p|=tagr1pc;
	tagr2p|=tagr2pc;

	CEBUG("\t TRMr?p:" << tagr1p << " " << tagr2p);

	// weak mismatches set if only one of the reads has a SRMr / CRMr
	//   mismatch
	if (tagr1p || tagr2p) {
	  CEBUG("Weak RMB mismatch found\n");
	  hasweakRMBmismatch=true;
	}

	// sheeesh, mismatch at both RMB position strong
	if (tagr1p && tagr2p) {
	  CEBUG("Strong RMB mismatch found\n");
	  hasstrongRMBmismatch=true;
	}
	// if mismatch at RMB position strong and weak, also care
	bool tagr1w=AS_readpool.getRead(id1).hasTag(Read::REA_tagentry_idWRMr, r1index);
	if (tagr2p && tagr1w) {
	  CEBUG("Strong-weak 1 RMB mismatch found\n");
	  hasstrongRMBmismatch=true;
	}
	bool tagr2w=AS_readpool.getRead(id2).hasTag(Read::REA_tagentry_idWRMr, r2index);
	if (tagr1p && tagr2w) {
	  CEBUG("Strong-weak 2 RMB mismatch found\n");
	  hasstrongRMBmismatch=true;
	}
	if(tagr1w && tagr2w) hasWRMrWRMrmismatch=true;
	if(penaliseSIOrmismatch) {
	  bool tagr1s=AS_readpool.getRead(id1).hasTag(Read::REA_tagentry_idSIOr, r1index);
	  bool tagr2s=AS_readpool.getRead(id2).hasTag(Read::REA_tagentry_idSIOr, r2index);

	  if (tagr1s || tagr2s) {
	    CEBUG("Weak RMB mismatch found\n");
	    hasweakRMBmismatch=true;
	  }
	  if((tagr1s && (tagr2s || tagr2p || tagr2w))
	     || (tagr2s && (tagr1s || tagr1p || tagr1w))){
	    CEBUG("SIOrXXX mismatch found\n");
	    hasstrongRMBmismatch=true;
	  }
	}
      }
    }

    CEBUG(endl);

  }

  // deferred setting of tags
  // read1
  if(r1newctagfrom>=0){
    if(r1newctagfrom>r1newctagto) std::swap(r1newctagfrom,r1newctagto);
    AS_tmptag_CRMr.from=r1newctagfrom;
    AS_tmptag_CRMr.to=r1newctagto;
    AS_readpool.getRead(id1).addTagO(AS_tmptag_CRMr);
    need2ndpass=true;
  }
  // read2
  if(r2newctagfrom>=0){
    if(r2newctagfrom>r2newctagto) std::swap(r2newctagfrom,r2newctagto);
    AS_tmptag_CRMr.from=r2newctagfrom;
    AS_tmptag_CRMr.to=r2newctagto;
    AS_readpool.getRead(id2).addTagO(AS_tmptag_CRMr);
    need2ndpass=true;
  }

  if(hasstrongRMBmismatch || numfrequencymismatch>2){
    CEBUG("Returning strong RMB mismatch\n");
    return -1;
  } else if(hasweakRMBmismatch) {
    CEBUG("Returning weak RMB mismatch\n");
    return -1;
  } else if(hasWRMrWRMrmismatch) {
    CEBUG("Returning WRMrWRMrmismatch\n");
    return 1;
  }
  CEBUG("Nothing special found\n");
  FUNCEND();
  return 0;
}

#define CEBUG(bla)





/*************************************************************************
 *
 * Takes an ADS and fills the hit/miss vector of the reads involved
 * Later on, this vector will be used to clip ends of the reads
 *
 *************************************************************************/

void Assembly::transcribeHits(AlignedDualSeq & ads)
{
  FUNCSTART("void Assembly::transcribeHits(AlignedDualSeq & ads)");

  if(AS_readhmcovered.size()==0 || AS_readhitmiss.size() == 0) {
    FUNCEND();
    return;
  }

  //cout << ads << endl << endl;

  int32 id1=ads.getID1();
  int32 id2=ads.getID2();

  AS_count_rhm[id1]++;
  AS_count_rhm[id2]++;

  int32 r1index, r2index, r1delta, r2delta;
  char const * alseq1, * alseq2, * rs1, * rs2;
  char c1, c2;

  alseq1=ads.getAlignedSequence(id1);
  alseq2=ads.getAlignedSequence(id2);

  if(ads.getSequenceDirection(id1)*ads.getSequenceDirection(id2)>0) {
    r1index=-(ads.getOffsetInAlignment(id1));
    r2index=-(ads.getOffsetInAlignment(id2));

    rs1=AS_readpool.getRead(id1).getClippedSeqAsChar();
    rs2=AS_readpool.getRead(id2).getClippedSeqAsChar();

    r1delta=1;
    r2delta=1;

    //cout << "Ding!\n";
  } else {
    if(ads.getSequenceDirection(id1)<0) {
      r1index=static_cast<int32>(AS_readhitmiss[id1].size())+ads.getOffsetInAlignment(id1);
      r2index=-(ads.getOffsetInAlignment(id2));

      rs1=AS_readpool.getRead(id1).getClippedComplementSeqAsChar();
      rs2=AS_readpool.getRead(id2).getClippedSeqAsChar();

      r1delta=-1;
      r2delta=1;

      //cout << "Dong1!\n";
    } else {
      r1index=-(ads.getOffsetInAlignment(id1));
      // TODO: warum zum ... ist hier -1 notwendig???
      r2index=static_cast<int32>(AS_readhitmiss[id2].size())+ads.getOffsetInAlignment(id2)-1;

      rs1=AS_readpool.getRead(id1).getClippedSeqAsChar();
      rs2=AS_readpool.getRead(id2).getClippedComplementSeqAsChar();

      r1delta=1;
      r2delta=-1;

      //cout << "Dong2!\n";
    }
  }

  //cout << "r1index: " << r1index << endl;
  //cout << "r2index: " << r2index << endl;
  //cout << "rs1: " << rs1[0] << rs1[1] << rs1[2] << rs1[3] << endl;
  //cout << "rs2: " << rs2[0] << rs2[1] << rs2[2] << rs2[3] << endl;
  //cout << "alseq1: " << alseq1[0] << alseq1[1] << alseq1[2] << alseq1[3] << endl;
  //cout << "alseq2: " << alseq2[0] << alseq2[1] << alseq2[2] << alseq2[3] << endl;
  //cout << "AS_readhitmiss[id1].size(): " << AS_readhitmiss[id1].size() << endl;
  //cout << "AS_readhitmiss[id2].size(): " << AS_readhitmiss[id2].size() << endl;

  for(uint32 i=0; i<ads.getTotalLen(); i++) {
    //cout << r1index << "\t" << r2index << '\n';
    if(r1index<0 || r2index<0) {
      if(r1index>=0) {
	alseq1++;
	//AS_readhmcovered[id1][r1index]++;
      }
      if(r2index>=0) {
	alseq2++;
	//AS_readhmcovered[id2][r2index]++;
      }
      r1index+=r1delta;
      r2index+=r2delta;
      continue;
    }
    if (r1index >= static_cast<int32>(AS_readhitmiss[id1].size())
	|| r2index >= static_cast<int32>(AS_readhitmiss[id2].size())) {
      if(r1index < static_cast<int32>(AS_readhitmiss[id1].size())) {
	alseq1++;
	//AS_readhmcovered[id1][r1index]++;
      }
      if(r2index < static_cast<int32>(AS_readhitmiss[id2].size())) {
	alseq2++;
	//AS_readhmcovered[id2][r2index]++;
      }
      r1index+=r1delta;
      r2index+=r2delta;
      continue;
    }

    c1=static_cast<char>(toupper(*alseq1));
    alseq1++;
    c2=static_cast<char>(toupper(*alseq2));
    alseq2++;
    //cout << "\t" << c1 << "\t" << c2;
    if(c1==c2
       || (c1=='N' || c2=='N' || c1=='X' || c2=='X')) {
      //cout << '\n';
      AS_readhmcovered[id1][r1index]++;
      AS_readhmcovered[id2][r2index]++;
      r1index+=r1delta;
      r2index+=r2delta;
      continue;
    }
    //cout << "\txxx\n";
    if(c1=='*') {
      if(rs1[r1index]=='*') {
	AS_readhmcovered[id1][r1index]++;
	AS_readhitmiss[id1][r1index]++;
	r1index+=r1delta;
      }
    } else {
      AS_readhmcovered[id1][r1index]++;
      AS_readhitmiss[id1][r1index]++;
      r1index+=r1delta;
    }
    if(c2=='*') {
      if(rs2[r2index]=='*') {
	AS_readhitmiss[id2][r2index]++;
	AS_readhmcovered[id2][r2index]++;
	r2index+=r2delta;
      }
    } else {
      AS_readhitmiss[id2][r2index]++;
      AS_readhmcovered[id2][r2index]++;
      r2index+=r2delta;
    }
  }


  FUNCEND();

  return;
}


/*************************************************************************
 *
 * in skimhitforsave_t, the "percent_in_overlap" is only an approximation
 *  of the true SW score ration
 *
 * In mapping assemblies with many repetitive elements and a less than
 *  true reference sequence, this can lead to reads being mapped to a
 *  completely wrong copy of the repetitive element.
 *
 * E.g.:
 *  a read
 *           .......................*...................
 *  with a homopolymer indel error may map to
 *           ..XXX......................................
 *  which would be another copy, three bases difference at end, but also
 *   different number of bases in homopolymer
 *
 * Therefore, for matches pf a read against a rail, go through a SW align
 *  to calculate the real score ration (and save that as
 *  percent_in_overlap)
 *
 * Always work on the the full posfmatch files, immediately after skim
 *
 *************************************************************************/

void Assembly::recalcNonPerfectSkimMappingsBySW(int32 version)
{
  FUNCSTART("void Assembly::recalcNonPerfectSkimMappingsBySW(int32 version)");

  try{
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    cout << "\nBackbone SW checks forward:\n";
    std::string ourskimfilename=AS_posfmatch_full_filename;
    rnpskmbs_helper(ourskimfilename.c_str(),version,1);

    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    cout << "\nBackbone SW checks reverse:\n";
    ourskimfilename=AS_poscmatch_full_filename;
    rnpskmbs_helper(ourskimfilename.c_str(),version,-1);
  }
  catch(Notify n) {
    n.handleError(THISFUNC);
  }

  FUNCEND();
}

void Assembly::rnpskmbs_helper(const std::string & filename, const int32 version, const int8 direction)
{
  FUNCSTART("void Assembly::rnpskmbs_helper(const std::string & filename, const int32 version, const int8 direction)");

  std::ofstream fout;
  if(AS_logflag_swbbcheck){
    std::string lfilename=buildFileName(version, "", "",
				   "elog.swbbcheck",
				   ".lst");
    if(direction>0){
      fout.open(lfilename, std::ios::out);
    }else{
      fout.open(lfilename, std::ios::out|std::ios::app);
    }
  }

  std::vector<Align> chkalign;
  setupAlignCache(chkalign);
  std::list<AlignedDualSeq> madsl;

  // temporary skim container
  std::vector<skimhitforsave_t> tsc;
  tsc.reserve(500000);

  FILE * finfout;
  finfout = fopen(filename.c_str(),"r+");
  if(finfout == nullptr) {
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }
  myFSeek(finfout, 0, SEEK_END);
  auto finsize=myFTell(finfout);
  rewind(finfout);

  uint64 numrecalc=0;

  int64 skimsprocessed=0;
  int64 numskims=finsize/sizeof(skimhitforsave_t);
  ProgressIndicator<std::streamsize> P (0, numskims, 2000);

  long freadpos=0;
  long fwritepos=0;

  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();

  while(!feof(finfout)){
    tsc.resize(tsc.capacity());
    myFSeek(finfout, freadpos, SEEK_SET);
    size_t numread=myFRead(&tsc[0],sizeof(skimhitforsave_t),tsc.capacity(),finfout);

    if(numread==0) break;
    CEBUG("rnpskmbs_helper: read " << numread << endl;)

    freadpos=myFTell(finfout);
    CEBUG("new freadpos: " << freadpos << endl);

    if(numread<tsc.capacity()) tsc.resize(numread);

    for(auto readI=tsc.begin(); readI != tsc.end(); ++readI, ++skimsprocessed){
      P.progress(skimsprocessed);
      if(readI->percent_in_overlap != 100
	 && (AS_readpool[readI->rid1].isRail()
	     || AS_readpool[readI->rid2].isRail())){
	// don't check if both are rails
	if(AS_readpool[readI->rid1].isRail()
	   && AS_readpool[readI->rid2].isRail()) {
	  continue;
	}
	// don't check if both are rails
	if(AS_permanent_overlap_bans.checkIfBanned(readI->rid1,readI->rid2) > 0) {
	  continue;
	}
	// don't check if length of any sequence is 0 (should never be, but ...)
	if(AS_readpool[readI->rid1].getLenClippedSeq() == 0
	   || AS_readpool[readI->rid2].getLenClippedSeq() == 0) continue;

	// OK, looks like we need to recompute that one and update the estimations
	++numrecalc;
	computeSWAlign(madsl,readI->rid1,readI->rid2,readI->eoffset,direction,chkalign, -1);
	//{
	//  std::list<AlignedDualSeq>::const_iterator mI=madsl.begin();
	//  uint32 i=0;
	//  for(; mI!=madsl.end(); ++mI, ++i){
	//    cout << "ZZZZZZZZ " << i << endl << *mI << endl;
	//  }
	//}
	minimiseMADSL(madsl);
	if(!madsl.empty()){
	  //cout << "ZZZZZZZ chosen \n" << madsl.front() << endl;
	  if(AS_logflag_swbbcheck){
	    fout << static_cast<int32>(direction) << "\tOLD: " << *readI;
	  }
	  readI->percent_in_overlap=madsl.begin()->getScoreRatio();
	  // very crude recalc of "numhashes" (ahem)
	  uint32 numhashes=madsl.begin()->getOverlapLen();
	  numhashes-=skim_params.sk_basesperhash-1;
	  numhashes-=madsl.begin()->getNumMismatches();
	  numhashes-=madsl.begin()->getNumGaps();
	  readI->numhashes=numhashes;
	  if(AS_logflag_swbbcheck){
	    fout << static_cast<int32>(direction) << "\tNEW: " << *readI << endl;
	  }
	}
      }
    }
    if(!tsc.empty()){
      myFSeek(finfout, fwritepos, SEEK_SET);
      if(myFWrite(&tsc[0],
		  sizeof(skimhitforsave_t),
		  tsc.size(),
		  finfout) != tsc.size()){
	MIRANOTIFY(Notify::FATAL, "Could not overwrite part of file. Changed permissions?");
      }
      fwritepos=myFTell(finfout);
      CEBUG("new fwritepos: " << fwritepos << endl);
    }
  }

  P.finishAtOnce();

  fclose(finfout);

  cout << "\nHad to recalculate " << numrecalc << " skims (out of " << numskims << ") with Smith-Waterman.\n";

  FUNCEND();
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////        Obsolete         ///////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
