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

#include "assembly.H"

// BOOST
#include <boost/algorithm/string.hpp>

#include "errorhandling/errorhandling.H"
#include "util/progressindic.H"
#include "util/machineinfo.H"
#include "util/fileanddisk.H"

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


//#define FUNCSTART(bla)  static const char * THISFUNC = bla"  "; {cout << THISFUNC << "enter.\n"; cout.flush();}
//#define FUNCTRACE(bla) { cout << THISFUNC << bla; cout.flush();}
//#define FUNCEND() {cout << THISFUNC << "exit.\n"; cout.flush();}
//
//#define CEBUG(bla)
//#define CEBUGF(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// sort id1 low to high,
//  on same id1 by skimweight high to low,
//  on same skimweight by id2 low to high
//  on same id2 by eoffset low to high
bool skimedges_t::stdSortCmp(const skimedges_t & a, const skimedges_t & b)
{
  if(a.rid1 == b.rid1){
    if(a.skimweight == b.skimweight){
      if(a.linked_with == b.linked_with){
	return a.eoffset < b.eoffset;
      }else{
	return a.linked_with < b.linked_with;
      }
    }else{
      return a.skimweight > b.skimweight;
    }
  }
  return a.rid1 < b.rid1;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::reduceSkimHits4(int32 version, const std::string prefix, const std::string postfix, const std::string logname)
{
  FUNCSTART("void Assembly::reduceSkimHits4(int32 version, const std::string prefix, const std::string postfix, const std::string logname)");
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  //saveVector(AS_writtenskimhitsperid,"as_wshpi.bin");
  //saveVector(AS_overlapcritlevell,"as_ocll.bin");
  //saveVector(AS_overlapcritlevelr,"as_oclr.bin");

  if(!AS_chimeracutflag.empty()){
    cout << "Chimeras were searched for ... looking for hits to purge.\n";
    AS_writtenskimhitsperid.clear();
    AS_writtenskimhitsperid.resize(AS_readpool.size(),0);
    rsh4_purgeSkimsOfReadsCutByChimera(AS_posfmatch_filename);
    rsh4_purgeSkimsOfReadsCutByChimera(AS_poscmatch_filename);
  }

  size_t totalskimhits=0;
  // idblocks contains the read ids of block boundaries;
  std::list<int64> idblocks;
  int64 maxseenblock=0;

  {
    skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();
    // 512MiB as default max memory for the skim edges
    int64 memtouse=static_cast<int64>(skim_params.sk_memcaphitreduction)*1024*1024;

    // either reuse existing AS_skim_edges size or (if wished) calc from scratch if not present
    if(AS_skim_edges.capacity()){
      memtouse=AS_skim_edges.capacity()*sizeof(skimedges_t);
    }else if(as_fixparams.as_automemmanagement && AS_systemmemory > 0 ){
      int64 onegig=1024*1024*1024;
      int64 memusedbymira=MachineInfo::getVMSize();

      int64 mem2keepfree=0;

      if(as_fixparams.as_amm_keeppercentfree){
	mem2keepfree=(AS_systemmemory*as_fixparams.as_amm_keeppercentfree)/100;
      }
      if(as_fixparams.as_amm_maxprocesssize){
	mem2keepfree=std::max(mem2keepfree,AS_systemmemory-static_cast<int64>(onegig*as_fixparams.as_amm_maxprocesssize)) ;
      }

      int64 memavail=AS_systemmemory-mem2keepfree-memusedbymira;

      cout << "System memory: " << AS_systemmemory
	   << "\nMem2keepfree: " << mem2keepfree
	   << "\nUsed by MIRA: " << memusedbymira
	   << "\nMem avail: " << memavail << endl;
      if(memavail>memtouse) {
	memtouse=memavail;

	cout << "rsh increased memtouse to: " << memtouse << endl;
      }else{
	cout << "rsh not increased.\n";
      }
    }

#ifdef ENABLE64
    //// test: cap at 32 GiB
    //if(memtouse>34359738368L) memtouse=34359738368L;
    // test: cap at 16 GiB
    if(memtouse>17179869184LL) memtouse=17179869184LL;
#else
    // test: cap at 1 GiB
    if(memtouse>1073741824LL) memtouse=1073741824LL;
#endif

    // loop to allocate container for skim edges
    // first, try to allocate as much as possible respectively needed
    // if that's not possible, try to allocate less and less until even
    //  allocation of 1 million skim edges fails and in that case, give up
    uint64 maxnumskimedges=static_cast<uint64>(memtouse/sizeof(skimedges_t));
    //maxnumskimedges=1;

    // try circumventing the "would extend AS_skim_edges" error by faking less available space
    // really ugly, memory wasting hack
    if(AS_skim_edges.capacity()>0){
      if(maxnumskimedges>100000){
	maxnumskimedges-=5000;
      }else if(maxnumskimedges>10000){
	maxnumskimedges-=2000;
      }else{
	maxnumskimedges-=maxnumskimedges/10;
	if(maxnumskimedges<10) maxnumskimedges=10;
      }
    }


    while(AS_skim_edges.empty()){
      totalskimhits=0;
      idblocks.clear();
      maxseenblock=0;

      uint64 thisblock=0;
      for(size_t i=0; i<AS_readpool.size(); i++){
	CEBUG("DBG: AS_wshpid[" << i << "]\t" << AS_writtenskimhitsperid[i] << '\n');
	thisblock+=AS_writtenskimhitsperid[i];
	totalskimhits+=AS_writtenskimhitsperid[i];
	if(i!=0 && i+1!=AS_readpool.size()){
	  if(thisblock >= maxnumskimedges
	     || thisblock+AS_writtenskimhitsperid[i+1] > maxnumskimedges){
	    if(thisblock>maxseenblock) maxseenblock=thisblock;
	    idblocks.push_back(i+1);
	    thisblock=0;
	  }
	}
      }
      if(thisblock!=0) {
	if(thisblock>maxseenblock) maxseenblock=thisblock;
	idblocks.push_back(AS_readpool.size());
      }

      CEBUG("DBG: TSH " << totalskimhits << endl);
      CEBUG("Max seen block 1: " << maxseenblock << '\n');
      CEBUG("Max num skim edges: " << maxnumskimedges << '\n');

      try{
	// no space reserved yet? Good, first pass. Looks like that due to chimera clipping and some
	//  editing in the pass, the memory needed in the first pass is slightly lower than in the next.
	// Add 5% safety buffer ... if that is not enough, well then tough luck.
	if(AS_skim_edges.capacity()==0){
	  maxseenblock+=maxseenblock/20;
	}
	// again, ugly hack to prevent "would extend AS_skim_edges" error by reserving more space than actually computed
	AS_skim_edges.reserve(maxseenblock+5000);
	// on successful reserve, break the while loop
	break;
      }
      catch(const std::bad_alloc & e){
	// could not reserve this size ... try again with smaller size
	if(maxnumskimedges>=2000000) {
	  cout << "Need to reduce memory 1 ... retrying.\n";
	  maxnumskimedges-=1000000;
	}else{
	  throw e;
	}
      }
      catch(...){
	// some systems throw a different exception ... oh well, duplicate
	// and throw here if it fails
	if(maxnumskimedges>=2000000) {
	  cout << "Need to reduce memory 2 ... retrying.\n";
	  maxnumskimedges-=1000000;
	}else{
	  AS_skim_edges.reserve(maxseenblock);
	  // if the above does not throw ...
	  cout << "Ouch, something went very wrong." << endl;
	  abort();
	}
      }
    }
    cout << "Edge vector capacity: " << AS_skim_edges.capacity() << "\n";
    cout << "Can load up to " << maxseenblock << " skim edges at once.\n";
  }

  cout << "Partitioning into " << idblocks.size() << " blocks.\nBlocks: ";
  for(auto ibe : idblocks){
    cout << ibe << "\t";
  }
  cout << endl;

  CEBUG("Max seen block 2: " << maxseenblock << '\n');

  cout << "We have " << totalskimhits << " skims in file.\n";

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 10\n";
    dumpMemInfo();
#endif

  AS_skimstaken.clear();
  AS_skimstaken.resize(totalskimhits,false);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 20\n";
    dumpMemInfo();
#endif

  AS_readmaytakeskim.clear();
  AS_readmaytakeskim.resize(AS_readpool.size(),true);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 30\n";
    dumpMemInfo();
#endif

  AS_numskimoverlaps.clear();
  AS_numskimoverlaps.resize(AS_readpool.size(),0);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 40\n";
    dumpMemInfo();
#endif

  AS_numleftextendskims.clear();
  AS_numleftextendskims.resize(AS_readpool.size(),0);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 40\n";
    dumpMemInfo();
#endif

  AS_numrightextendskims.clear();
  AS_numrightextendskims.resize(AS_readpool.size(),0);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 50\n";
    dumpMemInfo();
#endif

  AS_skimleftextendratio.clear();
  AS_skimleftextendratio.resize(AS_readpool.size(),0);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 60\n";
    dumpMemInfo();
#endif

  AS_skimrightextendratio.clear();
  AS_skimrightextendratio.resize(AS_readpool.size(),0);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 70\n";
    dumpMemInfo();
#endif

  bool hasshortreads=AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA] | AS_seqtypespresent[ReadGroupLib::SEQTYPE_ABISOLID];

  CEBUG("Has short reads: " << hasshortreads << '\n');

  bool onlyagainstrails=false;
  // very first call will be with version=0 ... pre_assembly
  //  don't set to true there as this analysis is then needed
  //  for multicopy analysis
  if(AS_hasbackbones
     && version >= as_fixparams.as_startbackboneusage_inpass){
    onlyagainstrails=true;
  }

  std::vector<uint64> blockpos;
  std::vector<size_t> blocklen;

  std::string fskimname=buildFileName(version, prefix, postfix,
			     as_fixparams.as_tmpf_normalisedskim,
			     ".bin");

  if(as_fixparams.as_dateoutput) dateStamp(cout);
  rsh4_denormaliseSkimHits(fskimname, idblocks, blockpos, blocklen);
  if(as_fixparams.as_dateoutput) dateStamp(cout);


#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 80\n";
    dumpMemInfo();
#endif

  rsh4_flagMulticopyReads(fskimname, blockpos, blocklen);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 90\n";
    dumpMemInfo();
#endif

  if(as_fixparams.as_dateoutput) dateStamp(cout);

  CEBUG("Only against rails: " << onlyagainstrails << '\n');

  if(AS_hasbackbones) {
    cout << "Step rail\n";
    rsh4_takeRailHits(fskimname, blockpos, blocklen);
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    rsh4_countTotalSkimsTaken();

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 a0\n";
    dumpMemInfo();
#endif
  }

  //rsh2_takeNeedAllOverlaps(idblocks);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  if(!onlyagainstrails){

    // temporary vectors for counting "nbest" in some of rsh4_* routines below
    std::vector<uint32> nbestl;
    std::vector<uint32> nbestr;
    nbestl.reserve(AS_readpool.size());
    nbestr.reserve(AS_readpool.size());

    cout << "Step 0\n";
    rsh4_weedOutBadHits(fskimname, blockpos, blocklen);
    if(as_fixparams.as_dateoutput) dateStamp(cout);

    if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
      if(0){
	cout << "Step cheat, not for any master/development release!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	exit(999);
	rsh4_takeAll(fskimname, blockpos, blocklen);
      }else
      if(!hasshortreads){
	cout << "Only long reads\n";
	cout << "Step 10\n";
	rsh4_takeNonReptHitsThatExtend(3,nbestl,nbestr,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 20\n";
	rsh4_takeReptPEPEHitsThatExtend(20,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 30\n";
	rsh4_takeReptNPENPEHitsThatExtend(2,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 40\n";
	rsh4_takeReptPENPEHitsThatExtend(2,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 50\n";
	rsh4_take100PC100bpHitsThatExtend(3,nbestl,nbestr,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 53\n";
	rsh4_takeNonTakenReadsWithHitsThatExtend(5,nbestl,nbestr,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 55\n";
	rsh4_takeNonTakenSideExtends(5,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 57\n";
	rsh4_takeIdenticalReadHits(3,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 60\n";
	rsh4_takeNonTakenReadsWithHits(5,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
      }else if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_SANGER]
	       || AS_seqtypespresent[ReadGroupLib::SEQTYPE_454GS20]
	       || AS_seqtypespresent[ReadGroupLib::SEQTYPE_TEXT]
	       || AS_seqtypespresent[ReadGroupLib::SEQTYPE_PACBIOHQ]
	       || AS_seqtypespresent[ReadGroupLib::SEQTYPE_PACBIOLQ]
	){
	// mixed: long and short reads
	// TODO: optimise
	cout << "Mixed long/short reads\n";
	cout << "Step 01\n";
	rsh4_take100PCMappingHits(fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 10\n";
	rsh4_takeNonReptHitsThatExtend(3,nbestl,nbestr,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 20\n";
	rsh4_takeReptPEPEHitsThatExtend(20,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 30\n";
	rsh4_takeReptNPENPEHitsThatExtend(2,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 40\n";
	rsh4_takeReptPENPEHitsThatExtend(2,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 50\n";
	rsh4_take100PC100bpHitsThatExtend(3,nbestl,nbestr,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 53\n";
	rsh4_takeNonTakenReadsWithHitsThatExtend(5,nbestl,nbestr,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 55\n";
	rsh4_takeNonTakenSideExtends(5,0,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 57\n";
	rsh4_takeIdenticalReadHits(3,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
	cout << "Step 60\n";
	rsh4_takeNonTakenReadsWithHits(5,fskimname, blockpos, blocklen);
	rsh4_countTotalSkimsTaken();
      }else{
	// short only
	// TODO: optimise
	cout << "Only short reads\n";
	if(1){
	  cout << "Step 10\n";
	  rsh4_takeNonReptHitsThatExtend(3,nbestl,nbestr,100,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 15\n";
	  rsh4_takeReptNonReptBorders(98,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 20\n";
	  rsh4_takeReptPEPEHitsThatExtend(20,100,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 30\n";
	  rsh4_takeReptNPENPEHitsThatExtend(2,100,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 40\n";
	  rsh4_takeReptPENPEHitsThatExtend(2,100,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 50\n";
	  rsh4_take100PC100bpHitsThatExtend(3,nbestl,nbestr,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 53\n";
	  rsh4_takeNonTakenReadsWithHitsThatExtend(5,nbestl,nbestr,100,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 55\n";
	  rsh4_takeNonTakenSideExtends(5,100,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 56\n";
	  rsh4_takeNonTakenSideExtends(5,95,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 57\n";
	  rsh4_takeIdenticalReadHits(3,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	  cout << "Step 60\n";
	  rsh4_takeNonTakenReadsWithHits(5,fskimname, blockpos, blocklen);
	  rsh4_countTotalSkimsTaken();
	}else{
	  cout << "Step cheat\n";
	  rsh4_takeAll(fskimname, blockpos, blocklen);
	}
      }
    }else{
      // EST / RNASeq
      cout << "Step 05\n";
      rsh4_takeWellconnected(3,nbestl,nbestr,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 10\n";
      rsh4_take100PC100bpHitsThatExtend(3,nbestl,nbestr,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 20\n";
      rsh4_takeReptPEPEHitsThatExtend(12,0,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 30\n";
      rsh4_takeNonTakenReadsWithHitsThatExtend(5,nbestl,nbestr,100,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 40\n";
      rsh4_takeNonTakenReadsWithHitsThatExtend(5,nbestl,nbestr,80,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 50\n";
      rsh4_takeNonTakenReadsWithHitsThatExtend(5,nbestl,nbestr,0,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 60\n";
      rsh4_takeNonTakenSideExtends(5,0,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 70\n";
      rsh4_takeIdenticalReadHits(3,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
      cout << "Step 80\n";
      rsh4_takeNonTakenReadsWithHits(5,fskimname, blockpos, blocklen);
      rsh4_countTotalSkimsTaken();
    }

    cout << "Step solexa by critlevel\n";
    rsh4_takeSolexaByCritLevel(0,3,nbestl,nbestr,fskimname, blockpos, blocklen);
    rsh4_takeSolexaByCritLevel(1,3,nbestl,nbestr,fskimname, blockpos, blocklen);

    cout << "Step template overlaps\n";
    rsh4_takeTemplateOverlaps(fskimname, blockpos, blocklen);
  }

  cout << "Step NAO\n";
  rsh4_takeNeedAllOverlaps_weakgood(fskimname, blockpos, blocklen);

  {
    uint64 counter=0;
    for(uint64 stid=0; stid < AS_skimstaken.size(); stid++){
      if(AS_skimstaken[stid]) counter++;
    }
    cout << "Total skims taken: " << counter << '\n';

#if 0
    // not yet
    // newedges_t has a 32 bit adsfindex ... not sure I want to change
    if(counter>0xffffffffffff){
      MIRANOTIFY(Notify::FATAL, "More than 2^48 skim hits would be taken ... this was not foreseen yet. Please contact the author.");
#else
    if(counter>0xffffffff){
      MIRANOTIFY(Notify::FATAL, "More than 2^32-1 skim hits would be taken ... this was not foreseen yet. Please contact the author.");
#endif
    }
  }

#if TRACKMEMUSAGE
  cout << "\ndmi rsh2 b0\n";
  dumpMemInfo();
#endif

  cout << "\nFiltering forward skims." << endl;
  if(as_fixparams.as_dateoutput) dateStamp(cout);

  std::string newname=AS_posfmatch_filename+".reduced";
  uint64 skimindex=0;
  rsh4_filterSkimHits(AS_posfmatch_filename,
		      newname,
		      skimindex);
  AS_posfmatch_filename=newname;
  cout << "Done.\nFiltering complement skims." << endl;
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  newname=AS_poscmatch_filename+".reduced";
  rsh4_filterSkimHits(AS_poscmatch_filename,
		      newname,
		      skimindex);
  AS_poscmatch_filename=newname;
  cout << "Done all filtering." << endl;
  if(as_fixparams.as_dateoutput) dateStamp(cout);

#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 c0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_skimrightextendratio" << endl);
  nukeSTLContainer(AS_skimrightextendratio);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 d0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_skimleftextendratio" << endl);
  nukeSTLContainer(AS_skimleftextendratio);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 e0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_numrightextendskims" << endl);
  nukeSTLContainer(AS_numrightextendskims);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 f0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_numleftextendskims" << endl);
  nukeSTLContainer(AS_numleftextendskims);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 g0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_numskimoverlaps" << endl);
  nukeSTLContainer(AS_numskimoverlaps);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 h0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_readmaytakeskim" << endl);
  nukeSTLContainer(AS_readmaytakeskim);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 i0\n";
    dumpMemInfo();
#endif

  CEBUG("Nuking AS_skimstaken" << endl);
  nukeSTLContainer(AS_skimstaken);
#if TRACKMEMUSAGE
    cout << "\ndmi rsh2 j0\n";
    dumpMemInfo();
#endif

    // BaCh 05.11.2010
    // clear, don't nuke AS_skim_edges, will be re-used as is in the next pass.
    // cause: the memory allocator may or may not give back the memory to the
    //  OS memory pool. If not, then in the next pass the calculation will be
    //  completely wrong. Better to get a grip on that memory and hold on tight

    AS_skim_edges.clear();

    // BaCh 06.03.2011
    // but help the VM when swapping this out by filling everything with 0
    if(AS_skim_edges.capacity()>0){
      memset(&AS_skim_edges[0],0,AS_skim_edges.capacity()*sizeof(skimedges_t));
    }

//  CEBUG("Nuking AS_skim_edges" << endl);
//  nukeSTLContainer(AS_skim_edges);
//#if TRACKMEMUSAGE
//    cout << "\ndmi rsh2 k0\n";
//    dumpMemInfo();
//#endif

  FUNCEND();

}

#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::rsh4_countTotalSkimsTaken()
{
  size_t counter=0;
  for(size_t stid=0; stid < AS_skimstaken.size(); stid++){
    if(AS_skimstaken[stid]) counter++;
  }
  cout << "Total skims taken: " << counter << '\n';
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::rsh4_filterSkimHits(const std::string & oldfilename, const std::string & newfilename, uint64 & skimindex)
{
  FUNCSTART("void Assembly::rsh4_filterSkimHits(const std::string & oldfilename, const std::string & newfilename, uint64 & skimindex)");

  char buf[10000];

  FILE * fin;
  FILE * fout;

  fin = fopen(oldfilename.c_str(),"r");
  if(fin == nullptr) {
    MIRANOTIFY(Notify::FATAL, "File not found: " << oldfilename);
  }
  myFSeek(fin, 0, SEEK_END);
  auto finsize=myFTell(fin);
  rewind(fin);

  fout = fopen(newfilename.c_str(),"w");
  if(fout == nullptr) {
    fclose(fin);
    MIRANOTIFY(Notify::FATAL, "Could not open file: " << newfilename);
  }

  std::ofstream logfout;
  if(AS_miraparams[0].getSpecialParams().mi_extended_log){
    std::string path,justfilename;
    splitFullPathAndFileName(newfilename,path,justfilename);
    std::string logfilename=path+"/elog.rsh4.filter."+justfilename;
    logfout.open(logfilename, std::ios::out|std::ios::app|std::ios::ate);
  }


  skimhitforsave_t tmpshfs;

  cout << "Writing reduced skim file:\n";

  CEBUG("Starting at skimindex: " << skimindex << endl);

  ProgressIndicator<std::streamsize> P(0, finsize-1);

  while(!feof(fin)){
    auto numread=myFRead(&tmpshfs,sizeof(tmpshfs),1,fin);
    if(numread){
      if(AS_skimstaken[skimindex]){
	if(AS_miraparams[0].getSpecialParams().mi_extended_log){
	  logfout << "Selected:\t" << AS_readpool[tmpshfs.rid1].getName()
		  << '\t' << AS_readpool[tmpshfs.rid2].getName()
		  << '\t' << tmpshfs;
	}
	myFWrite(&tmpshfs, sizeof(tmpshfs), 1, fout);
	if(ferror(fout)){
	  fclose(fout);
	  fclose(fin);
	  MIRANOTIFY(Notify::FATAL, "Could not write anymore to normalised skim file. Disk full? Changed permissions?");
	}
      }else{
	logfout << "Dropped:\t" << AS_readpool[tmpshfs.rid1].getName()
		<< '\t' << AS_readpool[tmpshfs.rid2].getName()
		<< '\t' << tmpshfs;
      }
      ++skimindex;
    }
    if(skimindex%1000 == 0 ) P.progress(myFTell(fin));
  }

  P.finishAtOnce();

  fclose(fout);
  fclose(fin);

  cout << "\nDone.\n";

  FUNCEND();
  return;
}

/*************************************************************************
 *
 * return number of skims loaded / resp. in AS_skim_edges
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
size_t Assembly::rsh4_getNextNormalisedSkimBlock(std::list<int64> & idblocks, int64 & blockstartid, int64 & blockendid)
{

  // handling of special case: one block (fits entirely in memory)
  // on very call: load the data
  // Then: if called for first time in a loop (start and end are 0)
  //  simply return size of loaded data
  //  else return 0 (== no more data to "load")
  if(idblocks.size() == 1){
    CEBUG("Only one block ... ");
    if(blockstartid==0 && blockendid == 0){
      CEBUG("first loop ...");
      if(!AS_skim_edges.empty()){
	CEBUG("already loaded.\n");
	blockendid=static_cast<int32>(AS_skim_edges.size());
	return static_cast<uint32>(AS_skim_edges.size());
      }
      CEBUG("not loaded yet.\n");
    }else{
      CEBUG("not in first loop, return 0.\n");
      return 0;
    }
  }

  size_t loadedskims=0;
  AS_skim_edges.clear();

  blockstartid=blockendid;
  // find the next block end
  for(auto idbI=idblocks.cbegin(); idbI != idblocks.end() && blockendid<=blockstartid; ++idbI){
    blockendid=*idbI;
  }

  // if found a next valid block, load it
  if(blockendid>blockstartid){
    CEBUG("Loading skims in id range " << blockstartid << " to " << blockendid-1 << endl);
    size_t skimindex=rsh4_loadNormalisedSkimHitBlock(AS_posfmatch_filename,
						     0, blockstartid, blockendid,
						     1,1);
    rsh4_loadNormalisedSkimHitBlock(AS_poscmatch_filename,
				    skimindex, blockstartid, blockendid,
				    1,-1);

    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    CEBUG("Loaded " <<  AS_skim_edges.size() << " elements.\nSorting ... ");
    cout.flush();

    // Apply malus to overlaps we do not want to be taken early
    for(auto & see : AS_skim_edges){
      uint32 malus=getOverlapMalusDivider(see.rid1, see.linked_with);
      if(malus>1){
	see.skimweight/=malus;
      }
    }

    // that thing can be rather big :-)
    mstd::psort(AS_skim_edges, skimedges_t::stdSortCmp);
    CEBUG("done.\n")

//    {
//      std::vector<skimedges_t>::iterator seI=AS_skim_edges.begin();
//      for(; seI != AS_skim_edges.end(); seI++){
//	CEBUG(*seI);
//      }
//    }

    return static_cast<uint32>(AS_skim_edges.size());
  }else{
    // nothing more to load, we'll stop
  }

  return loadedskims;
}
#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

size_t Assembly::rsh4_loadNormalisedSkimHitBlock(const std::string & filename, uint64 skimindex, int64 blockstartid, int64 blockendid, int8 rid1dir, int8 rid2dir)
{
  FUNCSTART("size_t Assembly::rsh4_loadNormalisedSkimHitBlock(const std::string & filename, uint64 skimindex, int64 blockstartid, int64 blockendid, int8 rid1dir, int8 rid2dir)");


  FILE * fin;
  fin = fopen(filename.c_str(),"r");
  if(fin == nullptr) {
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  skimedges_t tmpsedge;
  tmpsedge.setRID1dir(rid1dir);
  tmpsedge.setRID2dir(rid2dir);

  uint64 lineno=0;
  uint32 bannedoverlapsfound=0;
  size_t totalhits=0;

  skimhitforsave_t tmpshfs;

  CEBUG("Starting at " << skimindex << endl);

  while(!feof(fin)){
    lineno++;

    auto numread=myFRead(&tmpshfs,sizeof(tmpshfs),1,fin);

    if(!feof(fin)) {
      if(numread==0) {
	MIRANOTIFY(Notify::FATAL,"In elemcount " << lineno << "of file " << filename << ": expected 1 structure, got 0?!");
      }
      tmpsedge.rid1=tmpshfs.rid1;
      tmpsedge.linked_with=tmpshfs.rid2;
      if((tmpsedge.rid1>=blockstartid && tmpsedge.rid1<blockendid)
	 || (tmpsedge.linked_with>=blockstartid && tmpsedge.linked_with<blockendid)){
	// insert only where the reads are not permanently banned
	//  from overlapping
	if(!AS_permanent_overlap_bans.checkIfBanned(tmpsedge.rid1,tmpsedge.linked_with)){
	  tmpsedge.eoffset=tmpshfs.eoffset;

	  tmpsedge.skimweight=tmpshfs.numhashes*tmpshfs.percent_in_overlap*tmpshfs.percent_in_overlap;
	  tmpsedge.scoreratio=tmpshfs.percent_in_overlap;

	  tmpsedge.ol_stronggood  = tmpshfs.ol_stronggood  ;
	  tmpsedge.ol_weakgood    = tmpshfs.ol_weakgood    ;
	  tmpsedge.ol_belowavgfreq= tmpshfs.ol_belowavgfreq;
	  tmpsedge.ol_norept      = tmpshfs.ol_norept      ;
	  tmpsedge.ol_rept        = tmpshfs.ol_rept        ;

	  tmpsedge.skimindex=skimindex;
	  //cout << "tmpse: " << tmpsedge;

	  if(AS_skim_edges.size()+1 >= AS_skim_edges.capacity()){
	    cout << "there's gonna be a problem ..." << endl;
	    cout << AS_writtenskimhitsperid.size() << endl;
	    cout << AS_skim_edges.size() << endl;

	    std::vector<uint32> tmpperid(AS_writtenskimhitsperid.size(),0);
	    for(const auto & see : AS_skim_edges){
	      ++tmpperid[std::min(see.rid1,see.linked_with)];
	    }
	    ++tmpperid[std::min(tmpsedge.rid1, tmpsedge.linked_with)];

//	    for(size_t tmpi=0; tmpi<tmpperid.size(); ++tmpi){
//	      if(tmpperid[tmpi]>AS_writtenskimhitsperid[tmpi]){
//		cout << "FAIL: "
//		     << tmpi << "\t"
//		     << AS_writtenskimhitsperid[tmpi]
//		     << '\t'
//		     << tmpperid[tmpi]
//		     << endl;
//	      }else{
//		cout << "GOOD: "
//		     << tmpi << "\t"
//		     << AS_writtenskimhitsperid[tmpi]
//		     << '\t'
//		     << tmpperid[tmpi]
//		     << endl;
//	      }
//	    }
	    BUGIFTHROW(true, "Would extend memory of AS_skim_edges? Shouldn't be. File: " << filename << "\tLine: " << lineno << '\n' << tmpsedge << endl);
	  }
	  // DO NOT insert ids from past or future blocks as rid1, this leads to wrong interpretation of the block
	  //  and wrong/superfluous/not-intended selection of overlaps in the rsh_* routines if
	  //  the skim table is partitioned in multiple blocks
	  if(tmpsedge.rid1>=blockstartid  && tmpsedge.rid1<blockendid) AS_skim_edges.push_back(tmpsedge);

	  std::swap(tmpsedge.rid1,tmpsedge.linked_with);
	  tmpsedge.swapRID12dirs();
	  tmpsedge.eoffset=-tmpsedge.eoffset;
	  // DO NOT insert ids from past or future blocks as rid1, this leads to wrong interpretation of the block
	  //  and wrong/superfluous/not-intended selection of overlaps in the rsh_* routines if
	  //  the skim table is partitioned in multiple blocks
	  if(tmpsedge.rid1>=blockstartid && tmpsedge.rid1<blockendid) AS_skim_edges.push_back(tmpsedge);

	  totalhits++;
	} else {
	  bannedoverlapsfound++;
	}
      }
      skimindex++;
    }
  }

  fclose(fin);

  CEBUG("Ending at " << skimindex << endl);

  FUNCEND();
  return skimindex;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


void Assembly::rsh4_takeThisSkim(const skimedges_t & se, ADSEstimator & adse, bool calcadse)
{
  if(AS_skimstaken[se.skimindex]==true) return;
  AS_skimstaken[se.skimindex]=true;
  ++AS_numskimoverlaps[se.rid1];
  ++AS_numskimoverlaps[se.linked_with];

  //if(se.scoreratio < AS_skimleftextendratio[se.rid1]
  //   || se.scoreratio < AS_skimrightextendratio[se.rid1]){
  //  cout << "OUCH! " <<  static_cast<uint16>(AS_skimleftextendratio[se.rid1])
  //	 << ' '  << static_cast<uint16>(AS_skimrightextendratio[se.rid1])
  //	 << ' ' << *seI;
  //}

  if(calcadse) {
    adse.calcNewEstimateFromSkim(
      se.eoffset,
      AS_readpool[se.rid1].getLenClippedSeq(),
      AS_readpool[se.linked_with].getLenClippedSeq(),
      se.rid1,
      se.linked_with,
      se.getRID1dir(),
      se.getRID2dir());
  }
  if(adse.getEstimatedLeftExpand(se.rid1)){
    ++AS_numleftextendskims[se.rid1];
  }
  if(adse.getEstimatedRightExpand(se.rid1)){
    ++AS_numrightextendskims[se.rid1];
  }
  if(adse.getEstimatedLeftExpand(se.linked_with)){
    ++AS_numleftextendskims[se.linked_with];
  }
  if(adse.getEstimatedRightExpand(se.linked_with)){
    ++AS_numrightextendskims[se.linked_with];
  }
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::rsh4_denormaliseSkimHits(const std::string & targetfile, std::list<int64> & idblocks, std::vector<uint64> & blockpos, std::vector<size_t> & blocklen)
{
  FUNCSTART("void Assembly::rsh4_denormaliseSkimHits(std::string & targetfile, std::list<int64> & idblocks, std::vector<uint64> & blockpos, std::vector<size_t> & blocklen)");

  cout << "De-normalising SKIM hits ... (this will take a while)" << endl;

  FILE * fout=fopen(targetfile.c_str(), "w");

  blockpos.clear();
  blocklen.clear();

  int64 blockstartid=0;
  int64 blockendid=0;
  while(rsh4_getNextNormalisedSkimBlock(idblocks, blockstartid, blockendid)>0) {

    blockpos.push_back(myFTell(fout));
    blocklen.push_back(AS_skim_edges.size());

    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    {
      std::ostringstream ostrstr;
      byteToHumanReadableSize(static_cast<double>(AS_skim_edges.size()*sizeof(skimedges_t)), ostrstr);
      cout << "Writing normalised skimblock " << blockstartid << " (" << std::setw(12) << ostrstr.str() << ") ... "; cout.flush();
    }

    if(myFWrite(&AS_skim_edges[0],
		sizeof(skimedges_t),
		AS_skim_edges.size(),
		fout) != AS_skim_edges.size()){
      MIRANOTIFY(Notify::FATAL, "Could not write anymore to normalised skim file. Disk full? Changed permissions?");
    }
    cout << "done." << endl;
  }

  fclose(fout);

  FUNCEND();
}
//#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_getNextSkimBlock(const std::string & dnsfile, uint32 blocki, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  FUNCSTART("void Assembly::rsh4_getNextSkimBlock(const std::string & dnsfile, uint64 blockpos, size_t blocklen)");

  CEBUG("Loading block " << blocki << '\n');
  if(blocki == 0 && blockpos.size()==1 && !AS_skim_edges.empty()){
    CEBUG("Only one block, already loaded.\n");
    return;
  }

  CEBUG("Loading " << blocklen[blocki] << " elements " << " at offset " << blockpos[blocki] << '\n');

  FILE * fin;
  fin=fopen(dnsfile.c_str(), "r");

  if(myFSeek(fin, blockpos[blocki],SEEK_SET)) {
    MIRANOTIFY(Notify::FATAL, "Could not seek " << blockpos[blocki] << " bytes in file " << dnsfile << ". Was the file deleted? Disk full?");
  }

  AS_skim_edges.resize(blocklen[blocki]);

  if(myFRead(&AS_skim_edges[0],sizeof(skimedges_t),blocklen[blocki],fin) != blocklen[blocki]) {
      MIRANOTIFY(Notify::FATAL, "Expected to read " << blocklen[blocki] << " elements in file " << dnsfile << " but read less. Was the file deleted? Disk full?");
  }
  fclose(fin);

  CEBUG("Loaded " << blocklen[blocki] << " elements\n");

  //std::vector<skimedges_t>::const_iterator seI=AS_skim_edges.begin();
  //for(; seI != AS_skim_edges.end(); ++seI){
  //  if(seI->rid1 >= AS_readpool.size()
  //     || seI->linked_with >= AS_readpool.size()){
  //    MIRANOTIFY(Notify::FATAL, "skimedge invalid???\n" << *seI << endl);
  //  }
  //}

  FUNCEND();
  return;
}
#define CEBUG(bla)



/*************************************************************************
 *
 * Go through written skim hits on disk. All skims involving reads
 *  cut by chimera search are removed (their ADS Estimator would be
 *  totally wrong)
 *
 * Careful: this rewrite the original file.
 * Length will awlays be <= initial length
 *
 * Also rewrite AS_writtenskimhitsperid[], which must be cleared outside
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_purgeSkimsOfReadsCutByChimera(std::string & filename)
{
  FUNCSTART("void Assembly::rsh4_purgeSkimsOfReadsCutByChimera(std::string & filename)");

  // temporary skim container
  std::vector<skimhitforsave_t> tsc;

  FILE * finfout;
  finfout = fopen(filename.c_str(),"r+");
  if(finfout == nullptr) {
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  myFSeek(finfout, 0, SEEK_END);
  auto finsize=myFTell(finfout);
  rewind(finfout);

  uint64 lineno=0;
  uint32 bannedoverlapsfound=0;
  size_t totalhits=0;

  long freadpos=0;
  long fwritepos=0;

  while(!feof(finfout)){
    tsc.resize(500000);
    myFSeek(finfout, freadpos, SEEK_SET);

    auto numread=myFRead(&tsc[0],sizeof(skimhitforsave_t),tsc.capacity(),finfout);

    if(numread==0) break;
    CEBUG("rsh4_pSOFRCBC read " << numread << endl;)
    lineno+=numread;

    freadpos=myFTell(finfout);
    CEBUG("new freadpos: " << freadpos << endl);

    if(numread<tsc.capacity()) tsc.resize(numread);

    auto writeI=tsc.begin();
    auto readI=tsc.cbegin();  // keep outside loop, needed down below
    for(; readI != tsc.cend(); ++readI){
      bool del=false;
      if(AS_chimeracutflag[readI->rid1] || AS_chimeracutflag[readI->rid2]){
	del=true;
      }
      if(readI != writeI){
	*writeI=*readI;
      }
      ++writeI;
      if(del){
	--writeI;
	CEBUG("Has chimera: " << *readI);
      }else{
	++AS_writtenskimhitsperid[readI->rid1];
	++AS_writtenskimhitsperid[readI->rid2];
      }
    }

    // the resize thing is really not optimal ... one could write to file only a
    //  subset. However, at the moment just keep it for 100% safety
    //cout << "Purge skim data. Old size: " << tsc.size() << endl;
    tsc.resize(tsc.size()-(readI-writeI));
    //cout << "New size: " << tsc.size() << endl;

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

  cout << "truncating " << filename << " from " << finsize << " to " << fwritepos << endl;
  if(truncate(filename.c_str(),fwritepos)){
    MIRANOTIFY(Notify::FATAL, "Could not truncate normalised skim file? Strange ...");
  }

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
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeRailHits(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  FUNCSTART("void Assembly::rsh4_takeRailHits(std::list<int64> & idblocks)");
  CEBUG("rsh4 Take rail hits." << endl);

  ADSEstimator adse;

  // first, look which reads are completely covered by one rail
  //  as for those, we'll disallow overlaps with partial rails
  //  ...
  std::vector<uint8> fullycovered(AS_readpool.size(),0);

  /* ...
   BUT: status fullycovered can be achieved only for Skims where the
    skimweight is at least 10% (or 50%?)  of the best available Skim or else
    some spurious hits (without real alignment) in a contig may turn on the
    "fully covered" status whereas there may be a perfect hit at one of the
    contig ends
  */
  // best skim weight encountered by a read
  //  should be searched normally in two-pass loop, but as reads in MIRA
  //  are sorted by rails first, then other reads, we can do this is one loop
  std::vector<uint32> bestskimweight(AS_readpool.size(),0);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    //uint32 bestweight=0;
    for(const auto & see : AS_skim_edges){
      CEBUG("eSEd " << AS_readpool[see.rid1].getName() << "\t" << AS_readpool[see.linked_with].getName() << "\t" << see);
      if(see.skimweight > bestskimweight[see.rid1]) bestskimweight[see.rid1]=see.skimweight;
      if(see.skimweight > bestskimweight[see.linked_with]) bestskimweight[see.linked_with]=see.skimweight;

      // we want first non rail
      if(AS_readpool.getRead(see.rid1).isRail()) continue;

      adse.calcNewEstimateFromSkim(
	see.eoffset,
	AS_readpool[see.rid1].getLenClippedSeq(),
	AS_readpool[see.linked_with].getLenClippedSeq(),
	see.rid1,
	see.linked_with,
	see.getRID1dir(),
	see.getRID2dir());

      // fully covered only for skim weights >= 95% of best weight
      if(see.skimweight >= static_cast<uint32>(static_cast<uint64>(bestskimweight[see.rid1])*95/100)
	 // getEstimated*Expand() should be evaluated >= 0, but as the estimated expand
	 //  may not be 100% accurate anyway, this builds in a (small) security buffer
	 && adse.getEstimatedLeftExpand(see.rid1) > 0
	 && adse.getEstimatedRightExpand(see.rid1) > 0){
	CEBUG("Set fully cover " << AS_readpool[see.rid1].getName() << endl);
	fullycovered[see.rid1]=1;
      }
    }
  }

  for(uint32 i=0; i<AS_readpool.size(); ++i){
    if(AS_readpool[i].isRail() || AS_readpool[i].isBackbone()) continue;
    if(fullycovered[i]){
      CEBUG("Fully covered: " << AS_readpool[i].getName() << endl);
    }else{
      CEBUG("Not fully covered: " << AS_readpool[i].getName() << endl);
    }
  }

  // reuse now unneeded bestskimweight for potentialmapcounter
  std::vector<uint32> & potentialmapcounter=bestskimweight;
  potentialmapcounter.clear();
  potentialmapcounter.resize(AS_readpool.size(),0);
  std::vector<uint32> hasmapcounter(AS_readpool.size(),0);
  uint32 counter=0;

  // then, count how many reads map to each rail

  cout << "Count mappings:\n";

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 thisrid1=-1;
    //uint32 bestweight=0;
    auto seS=AS_skim_edges.cbegin();
    decltype(seS) seE;

    ProgressIndicator<size_t> P(0,AS_skim_edges.size());
    size_t pcounter=0;
    // careful: seS is also changed in the loop!
    for(; seS != AS_skim_edges.end(); ++seS, ++pcounter){
      if(pcounter%1000) P.progress(pcounter);
      if(seS->rid1 != thisrid1){
	thisrid1=seS->rid1;
	//bestweight=seS->skimweight;
      }
      CEBUG("seS " << AS_readpool[seS->rid1].getName() << "\t" << AS_readpool[seS->linked_with].getName() << "\t" << *seS);
      // we want first non rail
      if(AS_readpool.getRead(seS->rid1).isRail()) continue;

      // if read has full coverage by a rail, make sure this rail is
      //  covering fully, too
      if(fullycovered[seS->rid1]){
	adse.calcNewEstimateFromSkim(
	  seS->eoffset,
	  AS_readpool[seS->rid1].getLenClippedSeq(),
	  AS_readpool[seS->linked_with].getLenClippedSeq(),
	  seS->rid1,
	  seS->linked_with,
	  seS->getRID1dir(),
	  seS->getRID2dir());
	// if containment-level == 0, then some sequence has an overhang
	//  and we don't want that
	if(adse.getContainmentLevel() == 0){
	  CEBUG("Cont because cont-level " << adse.getContainmentLevel() << endl);
	  continue;
	}
      }

      CEBUG("Passed fully covered check\n");

      // search best weight of link with a rail
      // (and at the same time set seS as endpoint for [seS-seS[ range)
      uint32 bestrailweight=0;
      seE=seS;
      for(; seE != AS_skim_edges.end() && seE->rid1==thisrid1; ++seE){
	if(AS_readpool.getRead(seE->linked_with).isRail()
	   && !AS_skimstaken[seE->skimindex]){
	  bestrailweight=std::max(bestrailweight,seE->skimweight);
	}
      }

      CEBUG("brw: " << bestrailweight << '\n');

      auto seI=seS;
      for(; seI != seE; ++seI){
	if(seI->skimweight == bestrailweight
	   && !AS_skimstaken[seI->skimindex]
	   && AS_readpool.getRead(seI->linked_with).isRail()){
	  potentialmapcounter[seI->linked_with]++;
	  CEBUG("inc pmc[" << seI->linked_with << "]: " << potentialmapcounter[seI->linked_with] << '\n');
	}
      }

      // we looked at a given "normal" read (rid1) and increased all map counters of rail reads (linked_with)
      // now advance seS to the next rid1
      // incidentally, seE is already there (as is seI)
      seS=seE;
      // but the for would go one further, so decrease by one
      --seS;
    }
    P.finishAtOnce();
  }

  for(uint32 i=0; i<potentialmapcounter.size(); i++){
    CEBUG("pmc[" << i << "]: " << potentialmapcounter[i] << '\n');
  }


  cout << "Now take maps;\n";

  // now map the reads: the the rails with the lowest number of mapped
  //  reads, map it to the lowest potential match

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 thisrid1=-1;
    //uint32 bestweight=0;
    auto seS=AS_skim_edges.cbegin();
    decltype(seS) seE;

    ProgressIndicator<size_t> P(0,AS_skim_edges.size());
    size_t pcounter=0;
    // careful: seS is also changed in the loop!
    for(; seS != AS_skim_edges.end(); ++seS, ++pcounter){
      if(seS->rid1 != thisrid1){
	thisrid1=seS->rid1;
	//bestweight=seS->skimweight;
      }
      // we want first non rail
      CEBUG("Dubidi: " << *seS);
      if(AS_readpool.getRead(seS->rid1).isRail()) continue;
      CEBUG("Dubido\n");

      // if read has full coverage by a rail, make sure this rail is
      //  covering fully, too
      if(fullycovered[seS->rid1]){
	adse.calcNewEstimateFromSkim(
	  seS->eoffset,
	  AS_readpool[seS->rid1].getLenClippedSeq(),
	  AS_readpool[seS->linked_with].getLenClippedSeq(),
	  seS->rid1,
	  seS->linked_with,
	  seS->getRID1dir(),
	  seS->getRID2dir());
	// if containment-level == 0, then some sequence has an overhang
	//  and we don't want that
	if(adse.getContainmentLevel() == 0){
	  CEBUG("Cont 1 because cont-level " << adse.getContainmentLevel() << endl);
	  continue;
	}
      }

      // search best weight of link with a rail
      // (and at the same time set seS as endpoint for [seS-seS[ range)
      uint32 bestrailweight=0;
      seE=seS;
      for(; seE != AS_skim_edges.end() && seE->rid1==thisrid1; seE++){
	if(AS_readpool.getRead(seE->linked_with).isRail()
	   && !AS_skimstaken[seE->skimindex]){

	  // if read has full coverage by a rail, make sure this rail is
	  //  covering fully, too
	  if(fullycovered[seE->rid1]){
	    adse.calcNewEstimateFromSkim(
	      seE->eoffset,
	      AS_readpool[seE->rid1].getLenClippedSeq(),
	      AS_readpool[seE->linked_with].getLenClippedSeq(),
	      seE->rid1,
	      seE->linked_with,
	      seE->getRID1dir(),
	      seE->getRID2dir());
	    // if containment-level == 0, then some sequence has an overhang
	    //  and we don't want that
	    if(adse.getContainmentLevel() == 0){
	      CEBUG("Cont 2 because cont-level " << adse.getContainmentLevel() << endl);
	      continue;
	    }
	  }

	  bestrailweight=std::max(bestrailweight,seE->skimweight);
	}
      }

      //

      uint32 lowestnummapped=0xffffffff;
      uint32 lowestpotentialmapped=0xffffffff;
      auto seI=seS;
      auto seC=AS_skim_edges.cend();
      for(; seI != seE; ++seI){
	if(seI->skimweight == bestrailweight
	   && AS_readmaytakeskim[seI->rid1]
	   && !AS_skimstaken[seI->skimindex]
	   && AS_readpool.getRead(seI->linked_with).isRail()){

	  // if read has full coverage by a rail, make sure this rail is
	  //  covering fully, too
	  if(fullycovered[seI->rid1]){
	    adse.calcNewEstimateFromSkim(
	      seI->eoffset,
	      AS_readpool[seI->rid1].getLenClippedSeq(),
	      AS_readpool[seI->linked_with].getLenClippedSeq(),
	      seI->rid1,
	      seI->linked_with,
	      seI->getRID1dir(),
	      seI->getRID2dir());
	    // if containment-level == 0, then some sequence has an overhang
	    //  and we don't want that
	    if(adse.getContainmentLevel() == 0){
	      CEBUG("Cont 3 because cont-level " << adse.getContainmentLevel() << endl);
	      continue;
	    }
	  }

	  CEBUG(*seI);
	  if(hasmapcounter[seI->linked_with]<lowestnummapped){
	    lowestnummapped=hasmapcounter[seI->linked_with];
	    lowestpotentialmapped=potentialmapcounter[seI->linked_with];
	    seC=seI;
	  }else if(hasmapcounter[seI->linked_with]==lowestnummapped
		   && potentialmapcounter[seI->linked_with] < lowestpotentialmapped) {
	    lowestpotentialmapped=potentialmapcounter[seI->linked_with];
	    seC=seI;
	  }
	}
      }

      if(seC!=AS_skim_edges.end()){
	CEBUG("!!!chosen: " << *seC);
	potentialmapcounter[seC->linked_with]--;
	hasmapcounter[seC->linked_with]++;
	rsh4_takeThisSkim(*seC,adse,true);
	AS_readmaytakeskim[seC->rid1]=false;
	++counter;
      }

      // advance seS to next rid1 ... incidentally seE is already there
      seS=seE;
      // but the for would go one further, so decrease by one
      --seS;
    }
    P.finishAtOnce();
  }

  CEBUG("Taken " << counter << " hits." << endl);
  FUNCEND();
}
//#define CEBUG(bla)





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::rsh4_flagMulticopyReads(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  FUNCSTART("void Assembly::rsh4_flagMulticopyReads(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)");

  AS_multicopies.clear();
  AS_multicopies.resize(AS_readpool.size(),0);

  CEBUG("rsh4 flag multicopy reads." << endl);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    for(const auto & see : AS_skim_edges){
      if(see.ol_rept){
	AS_multicopies[see.rid1] = 1;
	AS_multicopies[see.linked_with] = 1;
      }
    }
  }

  //cout << "Multicopies:\n";
  //for(uint32 i=0; i<AS_readpool.size(); i++){
  //  cout << AS_readpool[i].getName() << "\t" << static_cast<uint16>(AS_multicopies[i]) << '\n';
  //}
  //cout << "Multicopies end\n";

  FUNCEND();
}


/*************************************************************************
 *
 * The match with the best score which extends to either side is taken as
 *  reference (-5%) for the minimum score ratio other matches should have
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_weedOutBadHits(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_weedOutBadHits." << endl);
  ADSEstimator adse;

  std::vector<bool> hasleft80100extend(AS_readpool.size(),false);
  std::vector<bool> hasright80100extend(AS_readpool.size(),false);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      CEBUG("blr: " << bestlextendratio << " brr: " << bestrextendratio << '\n');
      CEBUG(*seI);
      CEBUG("AS_readmaytakeskim[see.rid1]: " << AS_readmaytakeskim[see.rid1]);
      CEBUG("\nAS_readmaytakeskim[see.linked_with]: " << AS_readmaytakeskim[see.linked_with]);
      CEBUG("\nAS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG('\n');

      adse.calcNewEstimateFromSkim(
	see.eoffset,
	AS_readpool[see.rid1].getLenClippedSeq(),
	AS_readpool[see.linked_with].getLenClippedSeq(),
	see.rid1,
	see.linked_with,
	see.getRID1dir(),
	see.getRID2dir());

      if(AS_skimleftextendratio[see.rid1]== 0
	 || AS_skimrightextendratio[see.rid1]== 0){
	CEBUG("check!\n");

	if(AS_skimleftextendratio[see.rid1]==0 && adse.getEstimatedLeftExpand(see.rid1)>0){
	  CEBUG(" lextend");
	  AS_skimleftextendratio[see.rid1]=see.scoreratio;
	  if(AS_skimleftextendratio[see.rid1]>=5){
	    AS_skimleftextendratio[see.rid1]-=5;
	  }
	}
	if(AS_skimrightextendratio[see.rid1]==0 && adse.getEstimatedRightExpand(see.rid1)>0){
	  CEBUG(" rextend");
	  AS_skimrightextendratio[see.rid1]=see.scoreratio;
	  if(AS_skimrightextendratio[see.rid1]>=5){
	    AS_skimrightextendratio[see.rid1]-=5;
	  }
	}
	CEBUG('\n');

	if(AS_readpool[see.rid1].getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA && see.scoreratio==100){
	  if(adse.getEstimatedOverlap() >= (AS_readpool[see.rid1].getLenClippedSeq()*8)/10){
	    if(adse.getEstimatedLeftExpand(see.rid1)>0){
	      hasleft80100extend[see.rid1]=true;
	    }
	    if(adse.getEstimatedRightExpand(see.rid1)>0){
	      hasright80100extend[see.rid1]=true;
	    }
	  }
	}
      }
    }
  }

  CEBUG("skim extend ratios needed:\n");
  for(uint32 i=0;i<AS_skimleftextendratio.size();i++){
    CEBUG(AS_readpool[i].getName()
	  << '\t' << static_cast<uint16>(AS_skimleftextendratio[i])
	  << '\t' << static_cast<uint16>(AS_skimrightextendratio[i])
	  << '\n');
  }
}
//#define CEBUG(bla)




/*************************************************************************
 *
 * takes overlaps of 100% mapping short reads (sxa) to longer reads
 *  (sanger, 454, PacBio)
 *
 * Works only as long SOLID or other sequencing type not implemented
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_take100PCMappingHits(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  FUNCSTART("void Assembly::rsh4_take100PCMappingHits(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)");
  CEBUG("rsh4 Take 100% mapping hits." << endl);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

//  if(as_fixparams.as_dateoutput) {
//    cout << "100pcm\n";
//    dateStamp(cout);
//  }

  std::vector<uint32> potentialmapcounter(AS_readpool.size(),0);
  std::vector<uint32> hasmapcounter(AS_readpool.size(),0);
  size_t counter=0;

  // first, count how many short reads map to each long read

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 thisrid1=-1;
    auto seS=AS_skim_edges.cbegin();
    decltype(seS) seE;

    for(; seS != AS_skim_edges.cend(); ++seS){
      if(seS->rid1 != thisrid1){
	thisrid1=seS->rid1;
      }
      CEBUG(*seS);
      // first non rail
      //if(AS_readpool.getRead(seS->rid1).isRail()) continue;

      // we want first short and linked with a long read,
      //  i.e. continue if not the case
      if(!AS_readpool[seS->rid1].isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	 || !AS_readmaytakeskim[seS->rid1]) continue;
      CEBUG("s1\n");

      // search best 100% map weight of this short read with a long read
      // (and at the same time set seS as endpoint for [seS-seS[ range)
      uint32 bestmapweight=0;
      seE=seS;
      for(; seE != AS_skim_edges.end() && seE->rid1==thisrid1; ++seE){
	if(!AS_readpool.getRead(seE->linked_with).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	   && seE->scoreratio == 100
	   && seE->ol_stronggood
	   && !AS_skimstaken[seE->skimindex]){
	  bestmapweight=std::max(bestmapweight,seE->skimweight);
	}
      }

      CEBUG("rid1: " << thisrid1 << "\tbmw: " << bestmapweight << '\n');

      // take an edge only if it's the best weight
      for(auto seI=seS; seI != seE; ++seI){
 	CEBUG((seI->skimweight == bestmapweight)
	      << ' ' << !AS_readpool.getRead(seI->linked_with).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	      << ' ' << (seI->scoreratio == 100)
	      << ' ' << seI->ol_stronggood
	      << ' ' << !AS_skimstaken[seI->skimindex] << endl);

 	if(seI->skimweight == bestmapweight
	   && !AS_readpool.getRead(seI->linked_with).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	   && seI->scoreratio == 100
	   && seI->ol_stronggood
	   && !AS_skimstaken[seI->skimindex]){
	  potentialmapcounter[seI->linked_with]++;
	  CEBUG("inc pmc[" << seI->linked_with << "]: " << potentialmapcounter[seI->linked_with] << '\n');
	}
      }
      seS=seE;
      --seS;
    }
  }

  for(uint32 i=0; i<potentialmapcounter.size(); i++){
    CEBUG("pmc[" << i << "]: " << potentialmapcounter[i] << '\n');
  }


//  if(as_fixparams.as_dateoutput) {
//    cout << "100pcm2\n";
//    dateStamp(cout);
//  }

  // now map the Solexa reads to the long reads with the lowest number of mapped
  //  reads, map it to the lowest potential match

  ADSEstimator adse;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 thisrid1=-1;
    auto seS=AS_skim_edges.cbegin();
    decltype(seS) seE;

    for(; seS != AS_skim_edges.cend(); ++seS){
      if(seS->rid1 != thisrid1){
	thisrid1=seS->rid1;
      }
      // we want first to be Solexa
      if(!AS_readpool[seS->rid1].isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	 || !AS_readmaytakeskim[seS->rid1]) continue;

      // search best 100% map weight of this short read with a long read
      // (and at the same time set seS as endpoint for [seS-seS[ range)
      uint32 bestmapweight=0;
      seE=seS;
      for(; seE != AS_skim_edges.end() && seE->rid1==thisrid1; seE++){
	if(!AS_readpool.getRead(seE->linked_with).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	   && seE->scoreratio == 100
	   && seE->ol_stronggood
	   && !AS_skimstaken[seE->skimindex]){
	  bestmapweight=std::max(bestmapweight,seE->skimweight);
	}
      }

      if(bestmapweight>0){
	uint32 lowestnummapped=0xffffffff;
	uint32 lowestpotentialmapped=0xffffffff;
	auto seI=seS;
	auto seC=AS_skim_edges.cend();
	for(; seI != seE; ++seI){
	  if(seI->skimweight == bestmapweight
	     && !AS_readpool.getRead(seI->linked_with).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	     && seI->scoreratio == 100
	     && seI->ol_stronggood
	     && !AS_skimstaken[seI->skimindex]){
	    CEBUG(*seI);
	    if(hasmapcounter[seI->linked_with]<lowestnummapped){
	      lowestnummapped=hasmapcounter[seI->linked_with];
	      lowestpotentialmapped=potentialmapcounter[seI->linked_with];
	      seC=seI;
	    }else if(hasmapcounter[seI->linked_with]==lowestnummapped
		     && potentialmapcounter[seI->linked_with] < lowestpotentialmapped) {
	      lowestpotentialmapped=potentialmapcounter[seI->linked_with];
	      seC=seI;
	    }
	  }
	}

	if(seC!=AS_skim_edges.end()){
	  CEBUG("!!!chosen: " << *seC);
	  --potentialmapcounter[seC->linked_with];
	  ++hasmapcounter[seC->linked_with];
	  rsh4_takeThisSkim(*seC,adse,true);
	  AS_readmaytakeskim[seC->rid1]=false;
	  ++counter;
	}
      }

      seS=seE;
      --seS;
    }
  }

  for(uint32 i=0; i<AS_readmaytakeskim.size(); i++){
    CEBUG("rmts[" << i << "]: " << AS_readmaytakeskim[i] << '\n');
  }

//  if(as_fixparams.as_dateoutput) {
//    cout << "100pcm3\n";
//    dateStamp(cout);
//  }

  CEBUG("Taken " << counter << " hits." << endl);

  FUNCEND();
}
//#define CEBUG(bla)





//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeNonReptHitsThatExtend(uint32 nbest, std::vector<uint32> & nbestl, std::vector<uint32> & nbestr, uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeAvgFreqHitsThatExtend." << endl);

  nbestl.clear();
  nbestl.resize(AS_readpool.size(),nbest);
  nbestr.clear();
  nbestr.resize(AS_readpool.size(),nbest);

  ADSEstimator adse;
  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    //uint32 needlextend=nbest;
    //uint32 needrextend=nbest;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	//needlextend=nbest;
	//needrextend=nbest;
      }
      //CEBUG("nle: " << needlextend << " nre: " << needrextend << '\n');
      CEBUG(*seI);
      CEBUG("AS_readmaytakeskim[see.rid1]: " << AS_readmaytakeskim[see.rid1]);
      CEBUG("\nAS_readmaytakeskim[see.linked_with]: " << AS_readmaytakeskim[see.linked_with]);
      CEBUG("\nAS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG('\n');
      if((!see.ol_rept || see.ol_weakgood)
	 && see.scoreratio >= minscoreratio
	 && (nbestl[see.rid1]>0 || nbestr[see.rid1]>0)
	 && AS_readmaytakeskim[see.rid1]
	 && AS_readmaytakeskim[see.linked_with]
	 && AS_wellconnected[see.rid1]
	 && AS_wellconnected[see.linked_with]
	 && AS_skimstaken[see.skimindex] == false){

	CEBUG("check!\n");
	adse.calcNewEstimateFromSkim(
	  see.eoffset,
	  AS_readpool[see.rid1].getLenClippedSeq(),
	  AS_readpool[see.linked_with].getLenClippedSeq(),
	  see.rid1,
	  see.linked_with,
	  see.getRID1dir(),
	  see.getRID2dir());

	takeskim=false;
	if(nbestl[see.rid1]>0
	   && adse.getEstimatedLeftExpand(see.rid1)>0
	   && see.scoreratio >= AS_skimleftextendratio[see.rid1]){
	  CEBUG(" lextend");
	  --nbestl[see.rid1];
	  takeskim=true;
	}
	if(nbestr[see.rid1]>0
	   && adse.getEstimatedRightExpand(see.rid1)>0
	   && see.scoreratio >= AS_skimrightextendratio[see.rid1]){
	  CEBUG(" rextend");
	  --nbestr[see.rid1];
	  takeskim=true;
	}
	if(takeskim){
	  CEBUG(" takeskim");
	  rsh4_takeThisSkim(see,adse,false);
	  counter++;
	}
	CEBUG('\n');
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * New 14.10.2014 to get rept/nonrept borders covered as best as possible
 * (for scaffolding)
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeReptNonReptBorders(uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeReptNonReptBorders()" << endl);
  uint32 counter=0;
  for(uint32 bi=0; bi<blockpos.size(); ++bi){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);
    for(const auto & see : AS_skim_edges){
      if(see.ol_stronggood
	 && see.ol_rept
	 && see.scoreratio >= minscoreratio
	 && AS_skimstaken[see.skimindex] == false){
	// a "silent" take skim:
	// does not appear in AS_numskimoverlaps and AS_numleft/rightextends
	// so that other thresholds are not disturbed by this
	CEBUG(" takeskim");
	AS_skimstaken[see.skimindex]=true;
	++counter;
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeReptPEPEHitsThatExtend(uint32 nbest, uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeReptPEPEHitsThatExtend." << endl);
  ADSEstimator adse;

  size_t counter=0;

  std::vector<uint32> lextendbins(200);
  std::vector<uint32> rextendbins(200);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;

	// these size are enough for binning reads up to 2000 bases
	//  (longer reads fall into the last bin)
	// change this in case PacBio sequences arrive
	mstd::fill(lextendbins,nbest);
	mstd::fill(rextendbins,nbest);
      }
      CEBUG(*seI);
      CEBUG("AS_readmaytakeskim[see.rid1]: " << AS_readmaytakeskim[see.rid1]);
      CEBUG("\nAS_readmaytakeskim[see.linked_with]: " << AS_readmaytakeskim[see.linked_with]);
      CEBUG("\nAS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG("\nrid1 tpartner: " << AS_readpool[see.rid1].getTemplatePartnerID());
      CEBUG("\nlwid tpartner: " << AS_readpool[see.linked_with].getTemplatePartnerID());
      CEBUG('\n');
      if(!see.ol_weakgood
	 && see.scoreratio >= minscoreratio
	 && AS_readpool[see.rid1].getTemplatePartnerID() != -1
	 && AS_readpool[see.linked_with].getTemplatePartnerID() != -1
	 && AS_readmaytakeskim[see.rid1]
	 && AS_readmaytakeskim[see.linked_with]
	 && AS_wellconnected[see.rid1]
	 && AS_wellconnected[see.linked_with]
	 && AS_skimstaken[see.skimindex] == false){

	CEBUG("check!\n");
	adse.calcNewEstimateFromSkim(
	  see.eoffset,
	  AS_readpool[see.rid1].getLenClippedSeq(),
	  AS_readpool[see.linked_with].getLenClippedSeq(),
	  see.rid1,
	  see.linked_with,
	  see.getRID1dir(),
	  see.getRID2dir());

	takeskim=false;

	int32 bin=adse.getEstimatedLeftExpand(see.rid1)/10;
	if(bin>=lextendbins.size()) bin=lextendbins.size()-1;
	if(lextendbins[bin]>0
	   && see.scoreratio >= AS_skimleftextendratio[see.rid1]){
	  CEBUG(" lextend");
	  lextendbins[bin]--;
	  takeskim=true;
	}
	bin=adse.getEstimatedRightExpand(see.rid1)/10;
	if(bin>=rextendbins.size()) bin=rextendbins.size()-1;
	if(rextendbins[bin]>0
	   && see.scoreratio >= AS_skimrightextendratio[see.rid1]){
	  CEBUG(" rextend");
	  rextendbins[bin]--;
	  takeskim=true;
	}
	if(takeskim){
	  CEBUG(" takeskim");
	  rsh4_takeThisSkim(see,adse,false);
	  counter++;
	}
	CEBUG('\n');
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)


//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeReptNPENPEHitsThatExtend(uint32 nbest, uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeReptNPENPEHitsThatExtend." << endl);
  ADSEstimator adse;

  size_t counter=0;

  std::vector<uint32> lextendbins(200);
  std::vector<uint32> rextendbins(200);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	mstd::fill(lextendbins,nbest);
	mstd::fill(rextendbins,nbest);
      }
      CEBUG("nle: " << needlextend << " nre: " << needrextend << '\n');
      CEBUG(*seI);
      CEBUG("AS_readmaytakeskim[see.rid1]: " << AS_readmaytakeskim[see.rid1]);
      CEBUG("\nAS_readmaytakeskim[see.linked_with]: " << AS_readmaytakeskim[see.linked_with]);
      CEBUG("\nAS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG('\n');
      if(!see.ol_weakgood
	 && see.scoreratio >= minscoreratio
	 && AS_readpool[see.rid1].getTemplatePartnerID() == -1
	 && AS_readpool[see.linked_with].getTemplatePartnerID() == -1
	 && AS_readmaytakeskim[see.rid1]
	 && AS_readmaytakeskim[see.linked_with]
	 && AS_wellconnected[see.rid1]
	 && AS_wellconnected[see.linked_with]
	 && AS_skimstaken[see.skimindex] == false){

	CEBUG("check!\n");
	adse.calcNewEstimateFromSkim(
	  see.eoffset,
	  AS_readpool[see.rid1].getLenClippedSeq(),
	  AS_readpool[see.linked_with].getLenClippedSeq(),
	  see.rid1,
	  see.linked_with,
	  see.getRID1dir(),
	  see.getRID2dir());

	takeskim=false;

	int32 bin=adse.getEstimatedLeftExpand(see.rid1)/10;
	if(bin>=lextendbins.size()) bin=lextendbins.size()-1;
	if(lextendbins[bin]>0
	   && see.scoreratio >= AS_skimleftextendratio[see.rid1]){
	  CEBUG(" lextend");
	  lextendbins[bin]--;
	  takeskim=true;
	}
	bin=adse.getEstimatedRightExpand(see.rid1)/10;
	if(bin>=rextendbins.size()) bin=rextendbins.size()-1;
	if(rextendbins[bin]>0
	   && see.scoreratio >= AS_skimrightextendratio[see.rid1]){
	  CEBUG(" rextend");
	  rextendbins[bin]--;
	  takeskim=true;
	}
	if(takeskim){
	  CEBUG(" takeskim");
	  rsh4_takeThisSkim(see,adse,false);
	  ++counter;
	}
	CEBUG('\n');
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)



//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeReptPENPEHitsThatExtend(uint32 nbest, uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeReptNPENPEHitsThatExtend." << endl);
  ADSEstimator adse;

  size_t counter=0;

  std::vector<uint32> lextendbins(200);
  std::vector<uint32> rextendbins(200);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	mstd::fill(lextendbins,nbest);
	mstd::fill(rextendbins,nbest);
      }
      CEBUG("nle: " << needlextend << " nre: " << needrextend << '\n');
      CEBUG(*seI);
      CEBUG("AS_readmaytakeskim[see.rid1]: " << AS_readmaytakeskim[see.rid1]);
      CEBUG("\nAS_readmaytakeskim[see.linked_with]: " << AS_readmaytakeskim[see.linked_with]);
      CEBUG("\nAS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG('\n');

      uint32 pecount=0;
      if(AS_readpool[see.rid1].getTemplatePartnerID() == -1) pecount++;
      if(AS_readpool[see.linked_with].getTemplatePartnerID() == -1) pecount++;

      if(!see.ol_weakgood
	 && see.scoreratio >= minscoreratio
	 && pecount==1
	 && AS_readmaytakeskim[see.rid1]
	 && AS_readmaytakeskim[see.linked_with]
	 && AS_wellconnected[see.rid1]
	 && AS_wellconnected[see.linked_with]
	 && AS_skimstaken[see.skimindex] == false){

	CEBUG("check!\n");
	adse.calcNewEstimateFromSkim(
	  see.eoffset,
	  AS_readpool[see.rid1].getLenClippedSeq(),
	  AS_readpool[see.linked_with].getLenClippedSeq(),
	  see.rid1,
	  see.linked_with,
	  see.getRID1dir(),
	  see.getRID2dir());

	takeskim=false;

	int32 bin=adse.getEstimatedLeftExpand(see.rid1)/10;
	if(bin>=lextendbins.size()) bin=lextendbins.size()-1;
	if(lextendbins[bin]>0
	   && see.scoreratio >= AS_skimleftextendratio[see.rid1]){
	  CEBUG(" lextend");
	  lextendbins[bin]--;
	  takeskim=true;
	}
	bin=adse.getEstimatedRightExpand(see.rid1)/10;
	if(bin>=rextendbins.size()) bin=rextendbins.size()-1;
	if(rextendbins[bin]>0
	   && see.scoreratio >= AS_skimrightextendratio[see.rid1]){
	  CEBUG(" rextend");
	  rextendbins[bin]--;
	  takeskim=true;
	}
	if(takeskim){
	  CEBUG(" takeskim");
	  rsh4_takeThisSkim(see,adse,false);
	  counter++;
	}
	CEBUG('\n');
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)




//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeNonTakenReadsWithHitsThatExtend(uint32 nbest, std::vector<uint32> & nbestl, std::vector<uint32> & nbestr, uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeNonTakenReadsWithHitsThatExtend." << endl);

  nbestl.clear();
  nbestl.resize(AS_readpool.size(),nbest);
  nbestr.clear();
  nbestr.resize(AS_readpool.size(),nbest);

  ADSEstimator adse;

  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    //uint32 needlextend=nbest;
    //uint32 needrextend=nbest;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	//needlextend=nbest;
	//needrextend=nbest;
      }
      if(!AS_skimstaken[see.skimindex]
	 && (nbestl[see.rid1]>0
	     || nbestr[see.rid1]>0)){
	CEBUG("nle: " << needlextend << " nre: " << needrextend << '\n');
	CEBUG(*seI);
	CEBUG("AS_readmaytakeskim[see.rid1]: " << AS_readmaytakeskim[see.rid1]);
	CEBUG("\nAS_readmaytakeskim[see.linked_with]: " << AS_readmaytakeskim[see.linked_with]);
	CEBUG("\nAS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
	CEBUG('\n');
	if((!see.ol_rept || see.ol_weakgood)
	   && see.scoreratio >= minscoreratio
	   && AS_readmaytakeskim[see.rid1]
	   && AS_readmaytakeskim[see.linked_with]
// Too dangerous for low cov projects (U13)
//	   && AS_wellconnected[see.rid1]
//	   && AS_wellconnected[see.linked_with]
	  ){

	  CEBUG("check!\n");
	  adse.calcNewEstimateFromSkim(
	    see.eoffset,
	    AS_readpool[see.rid1].getLenClippedSeq(),
	    AS_readpool[see.linked_with].getLenClippedSeq(),
	    see.rid1,
	    see.linked_with,
	    see.getRID1dir(),
	    see.getRID2dir());

	  takeskim=false;
	  if(nbestl[see.rid1]>0
	     && adse.getEstimatedLeftExpand(see.rid1)>0
	     && see.scoreratio >= AS_skimleftextendratio[see.rid1]){
	    CEBUG(" lextend");
	    --nbestl[see.rid1];
	    takeskim=true;
	  }
	  if(nbestr[see.rid1]>0
	     && adse.getEstimatedRightExpand(see.rid1)>0
	     && see.scoreratio >= AS_skimrightextendratio[see.rid1]){
	    CEBUG(" rextend");
	    --nbestr[see.rid1];
	    takeskim=true;
	  }
	  if(takeskim){
	    CEBUG(" takeskim");
	    rsh4_takeThisSkim(see,adse,false);
	    counter++;
	  }
	  CEBUG('\n');
	}
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)




//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_take100PC100bpHitsThatExtend(uint32 nbest, std::vector<uint32> & nbestl, std::vector<uint32> & nbestr, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  return;
  CEBUG("rsh4_take100PC100bpHitsThatExtend\n");

  nbestl.clear();
  nbestl.resize(AS_readpool.size(),nbest);
  nbestr.clear();
  nbestr.resize(AS_readpool.size(),nbest);

  ADSEstimator adse;

  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	//needlextend=nbest;
	//needrextend=nbest;
      }
      CEBUG(*seI);
      CEBUG("AS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG("\nnbl: " << nbestl[see.rid1] << "\tnbr: " << nbestr[see.rid1] << '\n');
      if(see.scoreratio == 100
	 && (nbestl[see.rid1]>0
	     || nbestr[see.rid1]>0)
	 && !AS_skimstaken[see.skimindex]){
	adse.calcNewEstimateFromSkim(
	  see.eoffset,
	  AS_readpool[see.rid1].getLenClippedSeq(),
	  AS_readpool[see.linked_with].getLenClippedSeq(),
	  see.rid1,
	  see.linked_with,
	  see.getRID1dir(),
	  see.getRID2dir());

	CEBUG(AS_readpool[see.rid1].getName() << "\t" <<AS_readpool[see.linked_with].getName() << "\n");

	CEBUG("ele: " << adse.getEstimatedLeftExpand(see.rid1) << "\tere: " << adse.getEstimatedRightExpand(see.rid1) << "\teov: " << adse.getEstimatedOverlap() << '\n');

	takeskim=false;
	if(adse.getEstimatedOverlap()>=95){
	  if(nbestl[see.rid1]>0
	     && adse.getEstimatedLeftExpand(see.rid1)>0){
	    CEBUG(" lextend");
	    --nbestl[see.rid1];
	    takeskim=true;
	  }
	  if(nbestr[see.rid1]>0
	     && adse.getEstimatedRightExpand(see.rid1)>0){
	    CEBUG(" rextend");
	    --nbestr[see.rid1];
	    takeskim=true;
	  }
	}
	if(takeskim){
	  CEBUG(" takeskim");
	  rsh4_takeThisSkim(see,adse,false);
	  counter++;
	}
	CEBUG('\n');
      }
    }
  }

  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)





//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeWellconnected(uint32 nbest, std::vector<uint32> & nbestl, std::vector<uint32> & nbestr, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_take100PC100bpHitsThatExtend\n");

  nbestl.clear();
  nbestl.resize(AS_readpool.size(),nbest);
  nbestr.clear();
  nbestr.resize(AS_readpool.size(),nbest);

  ADSEstimator adse;

  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    bool takeskim=false;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	//needlextend=nbest;
	//needrextend=nbest;
      }
      CEBUG(*seI);
      CEBUG("AS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG("\nnbl: " << nbestl[see.rid1] << "\tnbr: " << nbestr[see.rid1] << '\n');
      if((nbestl[see.rid1]>0
	  || nbestr[see.rid1]>0)
	 && AS_wellconnected[see.rid1]
	 && AS_wellconnected[see.linked_with]
	 && !AS_skimstaken[see.skimindex]){
	adse.calcNewEstimateFromSkim(
	  see.eoffset,
	  AS_readpool[see.rid1].getLenClippedSeq(),
	  AS_readpool[see.linked_with].getLenClippedSeq(),
	  see.rid1,
	  see.linked_with,
	  see.getRID1dir(),
	  see.getRID2dir());

	CEBUG(AS_readpool[see.rid1].getName() << "\t" <<AS_readpool[see.linked_with].getName() << "\n");

	CEBUG("ele: " << adse.getEstimatedLeftExpand(see.rid1) << "\tere: " << adse.getEstimatedRightExpand(see.rid1) << "\teov: " << adse.getEstimatedOverlap() << '\n');

	takeskim=false;
	if(nbestl[see.rid1]>0
	   && adse.getEstimatedLeftExpand(see.rid1)>0){
	  CEBUG(" lextend");
	  --nbestl[see.rid1];
	  takeskim=true;
	}
	if(nbestr[see.rid1]>0
	   && adse.getEstimatedRightExpand(see.rid1)>0){
	  CEBUG(" rextend");
	  --nbestr[see.rid1];
	  takeskim=true;
	}

	if(takeskim){
	  CEBUG(" takeskim");
	  rsh4_takeThisSkim(see,adse,false);
	  counter++;
	}
	CEBUG('\n');
      }
    }
  }

  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)




//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeNonTakenSideExtends(uint32 nbest, uint8 minscoreratio, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeNonTakenSideExtends." << endl);
  ADSEstimator adse;

  size_t counter=0;

  std::vector<std::vector<skimedges_t>::const_iterator> edl_stronggoods;
  decltype(edl_stronggoods) edl_weakgoods;
  decltype(edl_stronggoods) edl_others;
  decltype(edl_stronggoods) edr_stronggoods;
  decltype(edl_stronggoods) edr_weakgoods;
  decltype(edl_stronggoods) edr_others;

  edl_stronggoods.reserve(1000);
  edl_weakgoods.reserve(1000);
  edl_others.reserve(1000);
  edr_stronggoods.reserve(1000);
  edr_weakgoods.reserve(1000);
  edr_others.reserve(1000);

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    bool needlextend=false;
    bool needrextend=false;
    bool takeskim=false;
    for(auto seI=AS_skim_edges.cbegin(); seI != AS_skim_edges.cend(); ++seI){
      if(seI->rid1 != lastrid1){
	CEBUG("New rid1: " << seI->rid1 << '\n');
	if(edl_stronggoods.size()){
	  counter+=rsh4_tNTSEhelper(nbest,adse,edl_stronggoods);
	}else if(edl_weakgoods.size()){
	  counter+=rsh4_tNTSEhelper(nbest,adse,edl_weakgoods);
	}else if(edl_others.size()){
	  counter+=rsh4_tNTSEhelper(nbest,adse,edl_others);
	}
	if(edr_stronggoods.size()){
	  counter+=rsh4_tNTSEhelper(nbest,adse,edr_stronggoods);
	}else if(edr_weakgoods.size()){
	  counter+=rsh4_tNTSEhelper(nbest,adse,edr_weakgoods);
	}else if(edr_others.size()){
	  counter+=rsh4_tNTSEhelper(nbest,adse,edr_others);
	}

	lastrid1=seI->rid1;
	needlextend=false;
	needrextend=false;
	if(AS_numleftextendskims[seI->rid1]==0) needlextend=true;
	if(AS_numrightextendskims[seI->rid1]==0) needrextend=true;
      }
      if((needlextend || needrextend)
	 && !AS_skimstaken[seI->skimindex]){
	CEBUG("nle: " << needlextend << " nre: " << needrextend << '\n');
	CEBUG(*seI);
	CEBUG("AS_readmaytakeskim[seI->rid1]: " << AS_readmaytakeskim[seI->rid1]);
	CEBUG("\nAS_readmaytakeskim[seI->linked_with]: " << AS_readmaytakeskim[seI->linked_with]);
	CEBUG("\nAS_skimstaken[seI->skimindex]: " << AS_skimstaken[seI->skimindex]);
	CEBUG('\n');
	if(seI->scoreratio >= minscoreratio){
	  CEBUG("check!\n");
	  adse.calcNewEstimateFromSkim(
	    seI->eoffset,
	    AS_readpool[seI->rid1].getLenClippedSeq(),
	    AS_readpool[seI->linked_with].getLenClippedSeq(),
	    seI->rid1,
	    seI->linked_with,
	    seI->getRID1dir(),
	    seI->getRID2dir());

	  takeskim=false;
	  if(needlextend>0
	     && adse.getEstimatedLeftExpand(seI->rid1)>0){
	    takeskim=true;
	    if(seI->ol_stronggood){
	      edl_stronggoods.push_back(seI);
	    }else if(seI->ol_weakgood){
	      edl_weakgoods.push_back(seI);
	    }else{
	      edl_others.push_back(seI);
	    }
	  }
	  if(!takeskim
	     && needrextend>0
	     && adse.getEstimatedRightExpand(seI->rid1)>0){
	    if(seI->ol_stronggood){
	      edr_stronggoods.push_back(seI);
	    }else if(seI->ol_weakgood){
	      edr_weakgoods.push_back(seI);
	    }else{
	      edr_others.push_back(seI);
	    }
	  }

	  CEBUG('\n');
	}
      }
    }

    if(edl_stronggoods.size()){
      counter+=rsh4_tNTSEhelper(nbest,adse,edl_stronggoods);
    }else if(edl_weakgoods.size()){
      counter+=rsh4_tNTSEhelper(nbest,adse,edl_weakgoods);
    }else if(edl_others.size()){
      counter+=rsh4_tNTSEhelper(nbest,adse,edl_others);
    }
    if(edr_stronggoods.size()){
      counter+=rsh4_tNTSEhelper(nbest,adse,edr_stronggoods);
    }else if(edr_weakgoods.size()){
      counter+=rsh4_tNTSEhelper(nbest,adse,edr_weakgoods);
    }else if(edr_others.size()){
      counter+=rsh4_tNTSEhelper(nbest,adse,edr_others);
    }

  }

  CEBUG("Taken " << counter << " hits." << endl);
}



size_t Assembly::rsh4_tNTSEhelper(uint32 nbest, ADSEstimator & adse, std::vector<std::vector<skimedges_t>::const_iterator> & sev)
{
  if(sev.empty()) return 0;

  if(sev.size()>nbest) sev.resize(nbest);
  const size_t counter=sev.size();

  for(uint32 sevi=0; sevi<sev.size(); ++sevi){
    rsh4_takeThisSkim(*(sev[sevi]),adse,true);
  }

  sev.clear();
  return counter;
}





/*************************************************************************
 *
 * Sometimes, identical reads end up as singlets as they were not linked
 *  by previous "takes" (they do not extend each other). This happens for
 *  reads linking to a cluster of similar reads with very high coverage
 *  where other filters say "we have enough".
 * This helps to counterbalance.
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeIdenticalReadHits(uint32 nbest, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeIdenticalReadHits." << endl);

  ADSEstimator adse;

  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    uint32 ncounter=0;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
	CEBUG("New rid1: " << see.rid1 << '\n');
	lastrid1=see.rid1;
	ncounter=0;
      }
      CEBUG(*seI);
      CEBUG("AS_skimstaken[see.skimindex]: " << AS_skimstaken[see.skimindex]);
      CEBUG("\nnbl: " << nbestl[see.rid1] << "\tnbr: " << nbestr[see.rid1] << '\n');
      if(see.eoffset==0
	 && see.scoreratio==100
	 && ncounter<nbest
	 && !AS_skimstaken[see.skimindex]
	 && AS_readpool[see.rid1].getLenClippedSeq() == AS_readpool[see.linked_with].getLenClippedSeq()
	){
	CEBUG(" takeskim\n");
	rsh4_takeThisSkim(see,adse,true);
	++ncounter;
	++counter;
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeNonTakenReadsWithHits(uint32 nbest, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeNonTakenReadsWithHits." << endl);
  ADSEstimator adse;

  size_t counter=0;

  // see how which reads need to be taken care of
  std::vector<bool> hasoverlaps(AS_numskimoverlaps.size(),false);

  for(uint32 i=0;i<AS_numskimoverlaps.size(); i++){
    hasoverlaps[i]=(AS_numskimoverlaps[i]>0);
  }

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    for(auto seI=AS_skim_edges.cbegin(); seI != AS_skim_edges.cend(); ++seI){
      if(seI->rid1 != lastrid1){
	CEBUG("New rid1: " << seI->rid1 << '\n');
	lastrid1=seI->rid1;
      }

      // have a look only at reads who had no overlaps yet
      if(hasoverlaps[seI->rid1]) continue;

      hasoverlaps[seI->rid1]=true;

      // take the n best
      for(uint32 tcount=0; tcount<nbest && seI != AS_skim_edges.end() && seI->rid1==lastrid1; ++tcount, ++seI){
	rsh4_takeThisSkim(*seI,adse,true);
	++counter;
      }
      --seI;
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
//#define CEBUG(bla)


#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::rsh4_takeAll(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeAll." << endl);
  ADSEstimator adse;

  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    for(const auto & see : AS_skim_edges){
      rsh4_takeThisSkim(see,adse,true);
      ++counter;
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
#define CEBUG(bla)



#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::rsh4_takeNeedAllOverlaps_weakgood(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeNeedAllOverlaps." << endl);

  bool foundnao=false;
  for(auto & naoe : AS_needalloverlaps){
    if(naoe){
      foundnao=true;
      break;
    }
  }
  if(!foundnao){
    CEBUG("None needed.\n");
  }else{
    ADSEstimator adse;
    size_t counter=0;

    for(uint32 bi=0; bi<blockpos.size(); bi++){
      rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

      for(const auto & see : AS_skim_edges){
	// have a look only at reads who had no overlaps yet
	if(AS_needalloverlaps[see.rid1]
	   || AS_needalloverlaps[see.linked_with]){

	  if(see.ol_weakgood
	     && AS_skimstaken[see.skimindex] == false){
	    rsh4_takeThisSkim(see,adse,true);
	    ++counter;
	  }
	}
      }
    }
    CEBUG("Taken " << counter << " hits." << endl);
  }
}

#define CEBUG(bla)




#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::rsh4_takeSolexaByCritLevel(uint32 ocvi, uint32 nbest, std::vector<uint32> & nbestl, std::vector<uint32> & nbestr, const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeSolexaByCritLevel." << endl);

  nbestl.clear();
  nbestl.resize(AS_readpool.size(),nbest);
  nbestr.clear();
  nbestr.resize(AS_readpool.size(),nbest);

  ADSEstimator adse;
  size_t counter=0;

  uint8 ocll=255;
  uint8 oclr=255;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    int32 lastrid1=-1;
    //uint32 needlextend=nbest;
    //uint32 needrextend=nbest;
    for(const auto & see : AS_skim_edges){
      if(see.rid1 != lastrid1){
      	lastrid1=see.rid1;
      //	needlextend=nbest;
      //	needrextend=nbest;
      }
      if(!AS_skimstaken[see.skimindex]
	 && (nbestl[see.rid1]>0
	     || nbestr[see.rid1]>0)){
	if((AS_readpool.getRead(see.rid1).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	    && (AS_overlapcritlevelvl[ocvi][see.rid1] != 0
		|| AS_overlapcritlevelvr[ocvi][see.rid1] != 0))
	   || (AS_readpool.getRead(see.linked_with).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	       && (AS_overlapcritlevelvl[ocvi][see.linked_with] != 0
		   || AS_overlapcritlevelvr[ocvi][see.linked_with] != 0))){

	  adse.calcNewEstimateFromSkim(
	    see.eoffset,
	    AS_readpool[see.rid1].getLenClippedSeq(),
	    AS_readpool[see.linked_with].getLenClippedSeq(),
	    see.rid1,
	    see.linked_with,
	    see.getRID1dir(),
	    see.getRID2dir());

	  // templated call, but function does not need to be different anyway, so take the vhash64_t version
	  Skim<vhash64_t>::getOverlapCriterionLevel(see.rid1,
						    AS_readpool[see.rid1].getSequencingType(),
						    adse,see.scoreratio,
						    ocll,oclr);
	  bool takeit=false;
	  if(nbestl[see.rid1]>0
	     && ocll <= 29
	     && AS_overlapcritlevelvl[ocvi][see.rid1] != 0
	     && adse.getEstimatedLeftExpand(see.rid1)>0){
	    --nbestl[see.rid1];
	    takeit=true;
	  }else if(nbestr[see.rid1]>0
		   && oclr <= 29
		   && AS_overlapcritlevelvr[ocvi][see.rid1] != 0
		   && adse.getEstimatedRightExpand(see.rid1)>0){
	    --nbestr[see.rid1];;
	    takeit=true;
	  }

	  if(takeit){
	    rsh4_takeThisSkim(see,adse,true);
	    ++counter;
	  }
	}
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}

#define CEBUG(bla)


void Assembly::rsh4_takeTemplateOverlaps(const std::string & dnsfile, const std::vector<uint64> & blockpos, const std::vector<size_t> & blocklen)
{
  CEBUG("rsh4_takeTemplateOverlaps." << endl);

  ADSEstimator adse;
  size_t counter=0;

  for(uint32 bi=0; bi<blockpos.size(); bi++){
    rsh4_getNextSkimBlock(dnsfile, bi, blockpos, blocklen);

    for(const auto & see : AS_skim_edges){
      if(AS_readpool.getRead(see.rid1).getTemplateID()>=0
	 && AS_readpool.getRead(see.rid1).getTemplateID() == AS_readpool.getRead(see.linked_with).getTemplateID()
	 && AS_skimstaken[see.skimindex] == false){
	rsh4_takeThisSkim(see,adse,true);
	++counter;
      }
    }
  }
  CEBUG("Taken " << counter << " hits." << endl);
}
