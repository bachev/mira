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

#include "util/stlimprove.H"

#include "mira/hashstats.H"

#include "mira/skim.H"
#include "mira/seqtohash.H"
#include "mira/readgrouplib.H"
#include "mira/readpool_io.H"
#include "mira/vhash.H"

/* *sigh*
   This is such a f*cked up construct ... I'm probably doing things wrong, but I see
   no other way to get it linked on Apple gcc and still have moderate compile
   turnaround times on Linux.
   Problem:
    - HashStatistics needs explicit template instantiation for the class (see end of file)
    - HashStatistics became big, so the .C was splitted into two files
    - each .C now needs an explicit template instantiation
    - the linker on OSX then complains about duplicate symbols *argh*

    Solution:
     - make explicit template instantiation for every function -> not going to happen
     - make "interesting" construction which includes the other .C when on OSX gcc
       but keep it separate otherwise.
       NOT PRETTY!
       The #define for this is set in autoconf
 */

#ifdef BOTCHEDEXPLICITINSTANTIATIONLINKER
#define BOTCHEDEXPLICITINSTANTIATIONLINKER_HELPERDEF
#include "mira/hashstats_sdbg.C"
#endif



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


template<typename TVHASH_T>
size_t HashStatistics<TVHASH_T>::HS_numelementsperbuffer=0;
template<typename TVHASH_T>
uint32 HashStatistics<TVHASH_T>::HS_hsfilemagic=0x4D4C6873;  // magic: "MLhs" MiraLibHashStat


#ifdef HSVHM_var
template<typename TVHASH_T>
const TVHASH_T HashStatistics<TVHASH_T>::HS_MAXVHASHMASK(0xffffffUL);
#endif

template<typename TVHASH_T>
TVHASH_T HashStatistics<TVHASH_T>::HS_vhashmask(0);


/*************************************************************************
 *
 * static
 *
 *************************************************************************/
template<typename TVHASH_T>
uint32 HashStatistics<TVHASH_T>::byteSizeOfHash(uint32 hashlen)
{
  if(hashlen <= 32) return 8;
  if(hashlen <= 64) return 16;
  if(hashlen <= 128) return 32;
  if(hashlen <= 256) return 64;
  // should never arrive here
  return -1;
}

/*************************************************************************
 *
 * for MiraDiff
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::inANotB(HashStatistics & hsa, HashStatistics & hsb)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::inANotB(HashStatistics & hsa, HashStatistics & hsb)");

  BUGIFTHROW(hsa.HS_hs_basesperhash!=hsb.HS_hs_basesperhash,"hsa.HS_hs_basesperhash " << hsa.HS_hs_basesperhash << " != hsb.HS_hs_basesperhash " << hsb.HS_hs_basesperhash);
  HS_hs_basesperhash=hsa.HS_hs_basesperhash;

  if(hsa.HS_hsv_hsshortcuts.empty()) hsa.priv_makeHashStatArrayShortcuts();
  if(hsb.HS_hsv_hsshortcuts.empty()) hsb.priv_makeHashStatArrayShortcuts();

  HS_hsv_hashstats.clear();
  HS_hsv_hsshortcuts.clear();
  for(auto & hsae : hsa.HS_hsv_hashstats){
    if(hsb.findVHash(hsae)==nullptr){
      HS_hsv_hashstats.push_back(hsae);
    }
  }
  priv_makeHashStatArrayShortcuts();
}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::inAAndB(HashStatistics & hsa, HashStatistics & hsb)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::inAAndB(HashStatistics & hsa, HashStatistics & hsb)");

  BUGIFTHROW(hsa.HS_hs_basesperhash!=hsb.HS_hs_basesperhash,"hsa.HS_hs_basesperhash " << hsa.HS_hs_basesperhash << " != hsb.HS_hs_basesperhash " << hsb.HS_hs_basesperhash);
  HS_hs_basesperhash=hsa.HS_hs_basesperhash;

  if(hsa.HS_hsv_hsshortcuts.empty()) hsa.priv_makeHashStatArrayShortcuts();
  if(hsb.HS_hsv_hsshortcuts.empty()) hsb.priv_makeHashStatArrayShortcuts();

  HS_hsv_hashstats.clear();
  HS_hsv_hsshortcuts.clear();
  for(auto & hsae : hsa.HS_hsv_hashstats){
    if(hsb.findVHash(hsae)!=nullptr){
      HS_hsv_hashstats.push_back(hsae);
    }
  }
  priv_makeHashStatArrayShortcuts();
}



/*************************************************************************
 *
 * hashstats
 *
 *************************************************************************/

template<typename TVHASH_T>
HashStatistics<TVHASH_T>::HashStatistics()
{
  HS_logflag_hashcount=false;
  HS_abortall=false;
};

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::discard()
{
  HS_hsv_hashstats.clear();
  HS_hsv_hashstatnodes.clear();
  HS_hsv_dbgseqs.clear();
  HS_hsv_hsshortcuts.clear();
  HS_hs_basesperhash=0;
  HS_hs_sortstatus=HSSS_NOTSORTED;
  HS_avg_freq=avg_freq_t();
  digiNormReset();

  removeDirectory(HS_tmpdirectorytodelete,true,true);
  HS_tmpdirectorytodelete.clear();
}

//#define SORTCOUT(bla)   {cout << bla; cout.flush();}
#define SORTCOUT(bla)

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sorthelper(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr, uint8 finalstatus, const char * compname, bool (*comparator)(const hashstat_t & a, const hashstat_t & b))
{
  SORTCOUT("HSsort: " << compname << " " << &hashstats[0] << " " << static_cast<void *>(sortstatusptr) << " --> ");
  if(hashstats.empty()){
    SORTCOUT("not sorted, empty.\n");
  }else if(sortstatusptr!=nullptr && *sortstatusptr==finalstatus){
    SORTCOUT("not sorted, already correct final status.\n");
  }else{
    SORTCOUT("need sort.\n");
    mstd::psort(hashstats, comparator);
  }
  if(sortstatusptr!=nullptr){
    *sortstatusptr=finalstatus;
    HS_hsv_hsshortcuts.clear();     // TODO: really
  }
}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sortLow24Bit(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr)
{
  priv_sorthelper(hashstats,sortstatusptr,
		  HSSS_LOW24BIT,"sortLow24Bit",sortHashStatComparatorByLow24bit);
}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sortLexicographicallyUp(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr)
{
  priv_sorthelper(hashstats,sortstatusptr,
		  HSSS_LEXIUP,"sortLexicographicallyUp",sortHashStatComparatorLexicographicallyUp);
}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sortByCountUp(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr)
{
  priv_sorthelper(hashstats,sortstatusptr,
		  HSSS_BYCOUNTUP,"sortByCountUp",sortHashStatComparatorByCountUp);
}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sortByCountDown(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr)
{
  priv_sorthelper(hashstats,sortstatusptr,
		  HSSS_BYCOUNTDOWN,"sortByCountDown",sortHashStatComparatorByCountDown);
}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sortLexByCount(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr)
{
  priv_sorthelper(hashstats,sortstatusptr,
		  HSSS_LEXBYCOUNT,"sortLexByCount",sortHashStatComparatorLexByCount);
}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sortMaskUp(std::vector<hashstat_t> & hashstats, uint8 * sortstatusptr)
{
  // for a mask sort, we need to force that (as the mask could have changed)
  if(sortstatusptr != nullptr){
    *sortstatusptr=HSSS_NOTSORTED;
  }
  priv_sorthelper(hashstats,sortstatusptr,
		  HSSS_MASKUP,"sortMaskUp",sortHashStatComparatorByMaskUp);
}

#undef SORTCOUT

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::hash2string(TVHASH_T hash, uint32 basesperhash, std::string & str)
{
  const static TVHASH_T lasttwobits(3);
  const static char acgtc[4]={'A','C','G','T'};

  str.clear();
  str.resize(basesperhash,' ');
  auto srI=str.rbegin();
  for(auto ci=0; ci<basesperhash; ++ci, ++srI){
    *srI=acgtc[static_cast<uint64>(hash&lasttwobits)];
    hash>>=2;
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::setHashFrequencyRatios(double freqest_minnormal,
						      double freqest_maxnormal,
						      double freqest_repeat,
						      double freqest_heavyrepeat,
						      double freqest_crazyrepeat,
						      uint32 rarekmercount,
						      uint32 nastyrepeatratio,
						      uint32 nastyrepeatcoverage)
{
  HS_freqest_minnormal=freqest_minnormal;
  HS_freqest_maxnormal=freqest_maxnormal;
  HS_freqest_repeat=freqest_repeat;
  HS_freqest_heavyrepeat=freqest_heavyrepeat;
  HS_freqest_crazyrepeat=freqest_crazyrepeat;
  HS_rarekmercount=rarekmercount;
  HS_nastyrepeatratio=nastyrepeatratio;
  HS_nastyrepeatcoverage=nastyrepeatcoverage;
}




/*************************************************************************
 *
 * compute kmer statistics from reads in ReadPool
 *
 * Note: does not delete the final hash statistics file on disk (only the
 *  temporary files)
 *
 * Returns explicitly:
 *   nothing
 *
 * Returns implicitly:
 *  - the created hash statistics is in memory, ready to be used
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::computeHashStatistics(ReadPool & rp, int32 memtouse, bool checkusedinassembly, bool alsorails, bool fwdandrev, uint32 fwdrevmin, uint32 rarekmerearlykill, uint32 basesperhash, const std::string & hashstatfilename, const std::string & tmpdirectory)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::computeHashStatistics(ReadPool & rp, bool checkusedinassembly, bool alsorails, bool fwdandrev, uint32 fwdrevmin, uint32 rarekmerearlykill, uint32 basesperhash, uint32 millionhashesperbuffer, std::string & hashstatfilename, const std::string & directory)");

  priv_phsCommon(&rp,basesperhash,memtouse,tmpdirectory);
  HS_hashstatfilename=hashstatfilename;
  BUGIFTHROW(HS_hashstatfilename.empty(),"HS_hashstatfilename.empty() ???");

  //if(!HS_hashfilenames.empty()) return;

  cout << "Writing temporary hstat files:\n";
  priv_hashes2disk(rp,
		   checkusedinassembly,alsorails,
		   fwdandrev,
		   basesperhash);

  dateStamp(cout);

  if(HS_abortall){
    discard();
  }else{
    prepareStreamFinalise(fwdrevmin, rarekmerearlykill);
  }

  return;
}


/*************************************************************************
 *
 * WARNING: atm do NOT use if there's any readpool with valid reads/name
 *  you want to keep. This trashed the readNameContainer!!! Need to find
 *  solution if needed
 *
 * Also: hashstatfilename & directory have slightly different meaning atm
 *  from the readpool version.
 * hashstatfilename = path for the final result to be saved (not done if empty)
 * directory = directory in which to create tmp files
 *
 * Important: make sure the object desctructor is called or else temporary
 *  files and directories might remain on the disk (e.g. due to hard process kills)
 *
 * compute kmer statistics from reads in "seqfiles"
 *
 * Streaming version of computeHashStatistics() above
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
uint64 HashStatistics<TVHASH_T>::computeHashStatistics(const std::list<std::string> & seqfiles, int32 memtouse, bool fwdandrev, uint32 fwdrevmin, uint32 rarekmerearlykill, uint32 basesperhash, const std::string & hashstatfilename, const std::string & directory)
{
  FUNCSTART("uint64 HashStatistics<TVHASH_T>::computeHashStatistics(const std::list<std::string> & seqfiles, int32 memtouse, bool fwdandrev, uint32 fwdrevmin, uint32 rarekmerearlykill, uint32 basesperhash, const std::string & hashstatfilename, const std::string & directory)");

  uint64 numtotalreads=0;

  ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_SOLEXA);  // saves time, no adjustment vector used
  rgid.setReadNamingScheme(ReadGroupLib::SCHEME_EMPTY);  // saves time during read.setName()

  ReadPool baitrp;
  ReadPoolIO rpio(baitrp);
  rpio.setAttributeFASTAQualFileWanted(false); // in case we load FASTAs
  rpio.setAttributeFASTQAPreserveComment(false); // no need, faster
  rpio.setAttributeFASTQTransformName(false); // no name transform (+ loading faster)
  rpio.setAttributeProgressIndicator(true);

  HS_tmpdirectorytodelete=directory+"/mirakmc_";
  myTempFileName(HS_tmpdirectorytodelete);
  if(HS_tmpdirectorytodelete.empty()) MIRANOTIFY(Notify::FATAL,"Could not create name for temporary directory ???");

  if(ensureDirectory(HS_tmpdirectorytodelete,true,true,false)){
    MIRANOTIFY(Notify::FATAL,"Could not create temporary directory '"+HS_tmpdirectorytodelete+"'? Check write permissions.");
  }
  auto hashstattmpname=HS_tmpdirectorytodelete+"/mirakmc.mhs.gz";

  prepareStreamHashStatistics(
    basesperhash,
    memtouse,
    hashstattmpname,HS_tmpdirectorytodelete);

  for(auto & bfn : seqfiles){
    if(HS_abortall) break;
    uint8 ziptype=0;
    std::string ft;
    std::string dummyfromstem;
    std::string dummypathto;

    auto actfn=guessFileAndZipType(bfn,dummypathto,dummyfromstem,ft,ziptype);

    cout << "Loading baits from " << actfn << ":\n";

    rpio.registerFile(ft,actfn,"",rgid,false);
    while(rpio.loadNextSeqs(-1,-1,100000)){
      if(HS_abortall) break;
      numtotalreads+=baitrp.size();

      for(uint32 actreadid=0; actreadid<baitrp.size(); ++actreadid){
	Read & actread= baitrp.getRead(actreadid);
	prepareStreamAddNextSequence(
	  actread.getClippedSeqAsChar(),
	  actread.getLenClippedSeq(),
	  actread.getName().c_str(),
	  actread.getSequencingType(),
	  false
	  );
	if(fwdandrev){
	  prepareStreamAddNextSequence(
	    actread.getClippedComplementSeqAsChar(),
	    actread.getLenClippedSeq(),
	    actread.getName().c_str(),
	    actread.getSequencingType(),
	    true
	    );
	}
      }

      // TODO: ... ouch ouch ouch ... trashing that container is a big NONO in general!!!
      //  change that ASAP (but how???)
      Read::trashReadNameContainer();
      baitrp.discard();
    }
  }
  rpio.discard();
  cout << endl;

  if(!HS_abortall){
    prepareStreamFinalise(1,0);
    if(!HS_abortall){
      if(!hashstatfilename.empty()) fileRename(hashstattmpname,hashstatfilename);
    }
  }
  removeDirectory(HS_tmpdirectorytodelete,true,true);

  if(HS_abortall){
    discard();
    numtotalreads=0;
  }

  return numtotalreads;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::prepareStreamHashStatistics(uint32 basesperhash, int32 memtouse, const std::string & hashstatfilename, const std::string & tmpdirectory)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::prepareStreamHashStatistics(uint32 basesperhash, int32 memtouse, const std::string & hashstatfilename, const std::string & tmpdirectory)");

  priv_phsCommon(nullptr,basesperhash,memtouse,tmpdirectory);
  HS_hashstatfilename=hashstatfilename;
  BUGIFTHROW(HS_hashstatfilename.empty(),"HS_hashstatfilename.empty() ???");

  return;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::prepareStreamFinalise(uint32 fwdrevmin, uint32 rarekmerearlykill)
{
  cout << "Flushing buffers to disk:" << endl;;
  ProgressIndicator<int32> P(0, HS_hashfilebuffer.size());
  for(size_t hbi=0; hbi<HS_hashfilebuffer.size(); ++hbi){
    if(HS_abortall) break;
    P.progress(hbi);
    HS_elementsperfile[hbi]+=
      priv_writeCompressedHFB(hbi,
			      HS_hashfilebuffer[hbi],
			      0,        // all kmers, don't throw away yet!
			      HS_hashfiles[hbi],
			      true);
    nukeSTLContainer(HS_hashfilebuffer[hbi]);
    gzclose(HS_hashfiles[hbi]);
  }
  P.finishAtOnce();
  cout << "done.\n";
  HS_hashfiles.clear();

  //dateStamp(cout);
  //exit(100);

  if(!HS_abortall){
    cout << "\nAnalysing hstat files:\n";
    priv_createHashStatisticsFile(HS_hashstatfilename,
				  HS_hashfilenames,
				  HS_elementsperfile,
				  fwdrevmin,
				  rarekmerearlykill
      );
    cout << "\n";
    dateStamp(cout);
  }

  cout << "clean up temporary stat files..."; cout.flush();
  // clean up temporary stat files
  for(uint32 hfni=0; hfni<HS_hashfilenames.size();++hfni){
    fileRemove(HS_hashfilenames[hfni],true);
  }
  dateStamp(cout); cout.flush();

  if(!HS_abortall){
    priv_makeHashStatArrayShortcuts();
    dateStamp(cout);
  }

  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * common part of prepareHashStatistics()
 *
 * ReadPool * rpptr may be nullptr
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_phsCommon(ReadPool * rpptr, uint32 basesperhash, int32 megabytesforbuffer, const std::string & tmpdirectory)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_phsCommon(bool alsosavesinglehashes, uint32 basesperhash, int32 megabytesforbuffer, const std::string & tmpdirectory)");

  HS_abortall=false;

  HS_hs_basesperhash=basesperhash;

  HS_hsv_hashstats.clear();
  HS_hsv_hsshortcuts.clear();

  HS_hashfilenames.clear();
  HS_elementsperfile.clear();
  HS_hashfiles.clear();
  HS_hashfilebuffer.clear();

  HS_avg_freq.corrected=0;
  HS_avg_freq.raw=0;
  HS_avg_freq.taken=0;

  dateStamp(cout);

  const size_t upperbases=2;

  BUGIFTHROW(basesperhash==0,"basesperhash == 0 ???");
  BUGIFTHROW(upperbases>=basesperhash,"upperbases (" << upperbases << ") >=basesperhash " << basesperhash << ") ???");

  size_t numfiles=1<<(upperbases*2);
  HS_rightshift=(basesperhash-upperbases)*2;

  CEBUG("bph: " << basesperhash << ".\n");
  CEBUG("Must create " << numfiles << " files.\n");
  CEBUG("Rightshift:" << HS_rightshift << '\n');
  CEBUG("sizeof(TVHASH_T): " << sizeof(TVHASH_T) << '\n');

  if(1){
    cout << "hashstat_t: " << sizeof(hashstat_t) << endl;
    int64 freemb=MachineInfo::getMemAvail();
    cout << "freemem: " << freemb << endl;
    freemb/=1024*1024;
    if(megabytesforbuffer<0){
      // keep that free
      cout << "Keep free\n";
      megabytesforbuffer=freemb+megabytesforbuffer;
    }else if(megabytesforbuffer<=100){
      // that's percent of memory
      cout << "Percentage\n";
      megabytesforbuffer=freemb*megabytesforbuffer/100;
    }
    cout << "mbfb1: " << megabytesforbuffer << '\n';
    if(megabytesforbuffer<512) megabytesforbuffer=512;
    if(sizeof(void *)>4 && megabytesforbuffer<2048) megabytesforbuffer=2048;
    cout << "mbfb2: " << megabytesforbuffer << '\n';
    HS_numelementsperbuffer=static_cast<size_t>(megabytesforbuffer)*1024*1024/numfiles/sizeof(hashstat_t);

    // see if we can get away with a smaller buffer size
    if(rpptr!=nullptr){
      uint64 tnumhashes=0;
      for(size_t rpi=0; rpi<rpptr->size(); ++rpi){
	auto & actread = rpptr->getRead(rpi);
	if(actread.hasValidData()
	   && !actread.isBackbone()
	   && !actread.isRail()
	   && actread.getLenClippedSeq()>=basesperhash){
	  tnumhashes+=actread.getLenClippedSeq()-basesperhash+1;
	}
      }
      cout << "TNH: " << tnumhashes << endl;
      const double fillratio=1.5;
      double xmillionelem=static_cast<double>(tnumhashes)/(fillratio*1024*1024*sizeof(hashstat_t));
      cout << "XME 1: " << xmillionelem << endl;
      if(xmillionelem<0.1) xmillionelem=0.1;
      cout << "XME 2: " << xmillionelem << endl;
      if(xmillionelem*1024*1024 < HS_numelementsperbuffer) HS_numelementsperbuffer=xmillionelem*1024*1024;
    }
    cout << "HS_nepb: " << HS_numelementsperbuffer << endl;
  }


  HS_hashfilebuffer.resize(numfiles);
  for(size_t nfi=0; nfi<numfiles; ++nfi){
    HS_hashfilebuffer[nfi].reserve(HS_numelementsperbuffer);
  }
  for(size_t nfi=0; nfi<numfiles; ++nfi){
    std::string fname(tmpdirectory+"/stattmp"+str(boost::format("%x") % nfi )+".bin.gz");
    HS_hashfilenames.push_back(fname);
    HS_hashfiles.emplace_back(gzopen(fname.c_str(), "wb1"));
    if(HS_hashfiles.back()==nullptr){
      MIRANOTIFY(Notify::FATAL,"Could not open " << fname << " for temporary stat file output? Disk full? Wrong path? Access permissions?");
    }
  }

  HS_elementsperfile.clear();
  HS_elementsperfile.resize(numfiles,0);
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

/*

  No Bloom filter. They work exactly as advertised (keep single hashes out of
  the buffers), but are a major disappointment in terms if speed: influence is
  negligible (1 to max 2% faster) for Solexa data. The additional memory is
  not worth it.

 */

#if __GNUC__ >= 3
#define prefetchwrite(p)     __builtin_prefetch((p), 1, 0)
#else
#define prefetchwrite(p)
#endif

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_hashes2disk(ReadPool & rp, bool checkusedinassembly, bool alsorails, bool fwdandrev, uint32 basesperhash)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_hashes2disk(ReadPool & rp, bool checkusedinassembly, bool alsorails, bool fwdandrev, uint32 basesperhash)");

#ifdef CLOCKSTEPS
  gettimeofday(&HS_CHEAT_tvfill,nullptr);
#endif

  Read::setCoutType(Read::AS_TEXTSHORT);
  ProgressIndicator<int32> P(0, rp.size());

  for(uint32 actreadid=0; actreadid<rp.size(); ++actreadid){
    if(HS_abortall) break;
    P.progress(actreadid);

    //if(actreadid>100) return;

    Read & actread= rp.getRead(actreadid);

    if(!actread.getReadGroupID().wantStatisticsCalc()) continue;

    // Has been taken out as hash statistics now also used for mirabait
    // TODO: check whether this has big influence on "normal" assembly jobs
    //  !!! it has ... for mapping assemblies !!!

    if(!actread.hasValidData()
       || actread.isBackbone()
       || (!alsorails && actread.isRail())
       || (checkusedinassembly && !actread.isUsedInAssembly())) continue;

    CEBUG("hname: " << actread.getName() << endl);
    //CEBUG("h2d new read: " << actread << endl);

    prepareStreamAddNextSequence(actread.getClippedSeqAsChar(),
				 actread.getLenClippedSeq(),
				 actread.getName().c_str(),
				 actread.getSequencingType(),
				 false
      );
    if(fwdandrev){
      prepareStreamAddNextSequence(actread.getClippedComplementSeqAsChar(),
				   actread.getLenClippedSeq(),
				   actread.getName().c_str(),
				   actread.getSequencingType(),
				   true
	);
    }

  }

  P.finishAtOnce();
  cout << "done\n";

  TEBUG("\nTiming fill HFB: " << diffsuseconds(HS_CHEAT_tvfill) << endl);

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::prepareStreamAddNextSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::prepareStreamAddNextSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse");

  if(slen<HS_hs_basesperhash) return;

  // We will use prefetch in the loops below, therefore make sure we do not prefetch memory
  //  which we do not own by making sure the loops flush the buffer before reaching
  //  the capacity of the buffer
  const size_t capacityflush=HS_hashfilebuffer[0].capacity()-2;
  // both memory write prefetches save ~15 to 20% time (well, 1s for 4m Solexa reads at 100bp)

  const auto basesperhash=HS_hs_basesperhash;

  hashstat_t tmpdh;
  tmpdh.hsc.seqtype=seqtype;

  size_t hashfilesindex;

  const uint8 * seq=static_cast<const uint8 *>(seqvoid);

  if(!isreverse){
    CEBUG("PSANS fwd " << namestr << endl);
    tmpdh.hsc.fcount=1;
    tmpdh.hsc.rcount=0;
    SEQTOHASH_LOOPSTART(TVHASH_T);
    {
      tmpdh.vhash=acthash;
      tmpdh.hsc.setLowPos(seqi-(basesperhash-1));
      hashfilesindex=static_cast<uint64>(tmpdh.vhash>>HS_rightshift);
      CEBUG("Want to write fwd: " << hash2string(acthash,HS_hs_basesperhash) << " " << tmpdh << " to " << hashfilesindex << endl);
      BUGIFTHROW(hashfilesindex>=HS_hashfiles.size(),"hashfilesindex>=HS_hashfiles.size() ???");

      if(HS_hashfilebuffer[hashfilesindex].size()==capacityflush){
#ifdef CLOCKSTEPS
	timeval now;
	gettimeofday(&now,nullptr);
#endif
	HS_elementsperfile[hashfilesindex]+=
	  priv_writeCompressedHFB(hashfilesindex,
				  HS_hashfilebuffer[hashfilesindex],
				  0, // no rarekmerearlykill
				  HS_hashfiles[hashfilesindex],
				  false);
#ifdef CLOCKSTEPS
	timeval after,diff;
	gettimeofday(&after,nullptr);
	timersub(&after,&now,&diff);
	timeradd(&diff,&HS_CHEAT_tvfill,&now);
	HS_CHEAT_tvfill=now;
#endif
      }
      HS_hashfilebuffer[hashfilesindex].push_back(tmpdh);
#ifndef _GLIBCXX_DEBUG
      // _GLIBCXX_DEBUG will barf on the [size()+1], but in normal operation we are allowed to
      //   do that as prefetching on non-existent memory is silently ignored
      prefetchwrite(&(HS_hashfilebuffer[hashfilesindex][HS_hashfilebuffer[hashfilesindex].size()+1]));
#endif
    }
    SEQTOHASH_LOOPEND;
  }else{
    CEBUG("PSANS rev " << namestr << endl);
    tmpdh.hsc.fcount=0;
    tmpdh.hsc.rcount=1;

    SEQTOHASH_LOOPSTART(TVHASH_T);
    {
      tmpdh.vhash=acthash;
      tmpdh.hsc.setLowPos(slen-seqi+1);
      hashfilesindex=static_cast<uint64>(tmpdh.vhash>>HS_rightshift);
      CEBUG("Want to write rev: " << hash2string(acthash,HS_hs_basesperhash) << " " << tmpdh << " to " << hashfilesindex << endl);
      BUGIFTHROW(hashfilesindex>=HS_hashfiles.size(),"hashfilesindex>=HS_hashfiles.size() ???");

      if(HS_hashfilebuffer[hashfilesindex].size()==capacityflush){
#ifdef CLOCKSTEPS
	timeval now;
	gettimeofday(&now,nullptr);
#endif
	HS_elementsperfile[hashfilesindex]+=
	  priv_writeCompressedHFB(hashfilesindex,
				  HS_hashfilebuffer[hashfilesindex],
				  0, // no rarekmerearlykill
				  HS_hashfiles[hashfilesindex],
				  false);
#ifdef CLOCKSTEPS
	timeval after,diff;
	gettimeofday(&after,nullptr);
	timersub(&after,&now,&diff);
	timeradd(&diff,&HS_CHEAT_tvfill,&now);
	HS_CHEAT_tvfill=now;
#endif
      }
      HS_hashfilebuffer[hashfilesindex].push_back(tmpdh);
#ifndef _GLIBCXX_DEBUG
      // _GLIBCXX_DEBUG will barf on the [size()+1], but in normal operation we are allowed to
      //   do that as prefetching on non-existent memory is silently ignored
      prefetchwrite(&(HS_hashfilebuffer[hashfilesindex][HS_hashfilebuffer[hashfilesindex].size()+1]));
#endif
    }
    SEQTOHASH_LOOPEND;

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
size_t HashStatistics<TVHASH_T>::priv_writeCompressedHFB(size_t hfindex, std::vector<hashstat_t> & hfb, uint32 rarekmerearlykill, gzFile & gzf, bool force)
{
  FUNCSTART("size_t HashStatistics<TVHASH_T>::priv_writeCompressedHFB(size_t hfindex, std::vector<hashstat_t> & hfb, uint32 rarekmerearlykill, gzFile & gzf, bool force)");
  size_t retvalue=0;
  CEBUG("WCHFB " << hfindex << " " << force << " " << hfb.size() << " " << hfb.capacity() << endl);
  if(hfb.size()){
    priv_compressHashStatBufferInPlace(hfb, rarekmerearlykill);
    if(force || hfb.size()>=hfb.capacity()*2/3){
      CEBUG("WCHFB write buffer " << hfindex << " " << 100*hfb.size()/hfb.capacity() << endl);
      auto writtenbytes=myGZWrite(gzf,&(hfb[0]),sizeof(hashstat_t)*hfb.size());
      if(writtenbytes != sizeof(hashstat_t)*hfb.size()){
	MIRANOTIFY(Notify::FATAL, "Could not write anymore to hash file. Disk full? Changed permissions?");
      }
      retvalue=hfb.size();
      hfb.clear();
    }else{
      CEBUG("WCHFB no write buffer " << hfindex << " " << 100*hfb.size()/hfb.capacity() << endl);
    }
  }
  return retvalue;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_calcFwdRevMinThresholdOKFlag(std::vector<hashstat_t> & hsb, uint32 fwdrevmin)
{
  for(auto & hsbe : hsb){
    hsbe.hsc.hasfwdrevthresholdok=((hsbe.hsc.fcount>=fwdrevmin) & (hsbe.hsc.rcount>=fwdrevmin));
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_compressHashStatBufferInPlace(std::vector<hashstat_t> & hsb, uint32 rarekmerearlykill)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_compressHashStatBufferInPlace(std::vector<hashstat_t> & hsb, bool alsosavesinglehashes)");

  if(hsb.empty()) return;

  CEBUG("CHSBIP\n");

  CEBUG("Sorting " << hsb.size() << " elements ..."); cout.flush();

#ifdef CLOCKSTEPS
  timeval tv,tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  //sort(hsb.begin(), hsb.end(), sortHashStatComparatorLexByCount);
  priv_sortLexByCount(hsb,nullptr);

  CEBUG("done.\n");
  TEBUG("\nTiming sort HFB: " << diffsuseconds(tv) << endl);

#ifdef CLOCKSTEPS
  gettimeofday(&tv,nullptr);
#endif

  bool haskmerforkf=false;
  bool haskmerforkr=false;
  bool hasfrthresholdok=false;

  uint8 thisseqtype=0;
  uint32 thishashcounter=0;
  uint32 thishashfcounter=0;
  uint32 thishashrcounter=0;
  uint16 thislowpos=0;
  hashstat_t tmphs;
  auto srcI=hsb.cbegin();
  auto dstI=hsb.begin();
  TVHASH_T thishash=(srcI->vhash);
  // setting this leads the very first iteration of the main loop
  //  to set correct values
  --thishash;
  for(; srcI!=hsb.cend(); ++srcI){
    CEBUG("cl " << srcI-hsb.cbegin() << "\t" << srcI->vhash << " " << thishash << endl);
    if(srcI->vhash != thishash){
      // save only hashes that appeared at least 'rarekmerearlykill' time
      if(thishashcounter >0 && thishashcounter>=rarekmerearlykill){
	tmphs.vhash=thishash;
	tmphs.hsc.fcount=thishashfcounter;
	tmphs.hsc.rcount=thishashrcounter;
	tmphs.hsc.setLowPos(thislowpos);
	tmphs.hsc.seqtype=thisseqtype;
	tmphs.hsc.iskmerforkf=haskmerforkf;
	tmphs.hsc.iskmerforkr=haskmerforkr;
	tmphs.hsc.hasfwdrevthresholdok=hasfrthresholdok;
	CEBUG("Write mid to " << dstI-hsb.begin() << " from " << srcI-hsb.begin() << ": " << hash2string(tmphs.vhash,HS_hs_basesperhash) << "\t" << tmphs << '\n');

	*dstI=tmphs;
	++dstI;
      }
      thishashfcounter=0;
      thishashrcounter=0;
      hasfrthresholdok=false;
      haskmerforkf=false;
      haskmerforkr=false;
      thishash=srcI->vhash;
      thislowpos=srcI->hsc.getLowPos();
      thishashcounter=0;
      thishashfcounter=0;
      thishashrcounter=0;
      thisseqtype=srcI->hsc.seqtype;
      CEBUG("New vhash: " << hash2string(thishash,HS_hs_basesperhash) << "\t" << *srcI << endl);
    }else{
      CEBUG("Existing vhash: " << *srcI << endl);
    }
    thishashfcounter+=srcI->hsc.fcount;
    thishashrcounter+=srcI->hsc.rcount;
    thishashcounter+=srcI->hsc.getCount();
    if(srcI->hsc.getLowPos() < thislowpos) thislowpos=srcI->hsc.getLowPos();
    if(srcI->hsc.seqtype != thisseqtype) thisseqtype=MULTISEQTYPE;
    haskmerforkf|=srcI->hsc.iskmerforkf;
    haskmerforkr|=srcI->hsc.iskmerforkr;
    hasfrthresholdok|=srcI->hsc.hasfwdrevthresholdok;

    CEBUG("thc: " << thishashcounter << "\thf: " << thishashfcounter << "\thr: " << thishashrcounter << "\ttlp: " << thislowpos
	  << "\thfrto: " << hasfrthresholdok << "\t" << rarekmerearlykill << endl);
  }

  // we're out of the loop, write last elements if there were any
  if(thishashcounter>0 && thishashcounter>=rarekmerearlykill){
    tmphs.vhash=thishash;
    tmphs.hsc.fcount=thishashfcounter;
    tmphs.hsc.rcount=thishashrcounter;
    tmphs.hsc.setLowPos(thislowpos);
    tmphs.hsc.seqtype=thisseqtype;
    tmphs.hsc.iskmerforkf=haskmerforkf;
    tmphs.hsc.iskmerforkr=haskmerforkr;
    tmphs.hsc.hasfwdrevthresholdok=hasfrthresholdok;
    CEBUG("Write end to " << dstI-hsb.begin() << " from " << srcI-hsb.begin() << ": " << hash2string(tmphs.vhash,HS_hs_basesperhash) << "\t" << tmphs << '\n');
    *dstI=tmphs;
    ++dstI;
  }

  TEBUG("Timing compress HFB: " << diffsuseconds(tv) << endl);
  TEBUG("Timing compressHashStatBufferInPlace: " << diffsuseconds(tvtotal) << endl);

  hsb.resize(dstI-hsb.begin());

  CEBUG("New hsb size: " << hsb.size() << endl);

#ifndef PUBLICQUIET
  {
    uint64 numsingle=0;
    uint64 nummulti=0;
    auto eI=hsb.cend();
    for(auto iI=hsb.cbegin(); iI!=eI; ++iI){
      if(iI->hsc.getCount()==1){
	++numsingle;
      }else{
	++nummulti;
      }
    }
    cout << "\nnumsingle: " << numsingle << "\nnummulti:  " << nummulti << endl;
  }
#endif

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * 1) sorts every hashfile
 * 2) writes a first hash statistics file hashfile by hashfile
 * 3) loads the above created hash statistics file
 * 4) calculate statistics on that
 * 5) sorts hash statistics by low24 bits (directly usable by
 *    makeHashStatArrayShortcuts())
 * 6) saves final hash statistics file
 *
 * returns:
 *  - by value: number of elements in hash statistics file
 *  - name of the hash statistics file in the call by reference variable
 *  - the created hash statistics is in memory, ready to be used
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
size_t HashStatistics<TVHASH_T>::priv_createHashStatisticsFile(const std::string & hashstatfilename, std::vector<std::string> & hashfilenames, std::vector<size_t> & elementsperfile, uint32 fwdrevmin, uint32 rarekmerearlykill)
{
  FUNCSTART("size_t HashStatistics<TVHASH_T>::priv_createHashStatisticsFile(std::string & hashstatfilename, std::vector<std::string> & hashfilenames, std::vector<size_t> & elementsperfile, uint32 fwdrevmin, bool alsosavesinglehashes)");

  BUGIFTHROW(hashstatfilename.empty(),"hashstatfilename.empty() ???");

  size_t maxelementsperfile=0;

  for(size_t fi=0; fi< elementsperfile.size(); ++fi){
    maxelementsperfile=std::max(maxelementsperfile,elementsperfile[fi]);
  }
  CEBUG("Max elements per file: " << maxelementsperfile << '\n');

  auto gzf=gzopen(hashstatfilename.c_str(),"wb1");
  if(gzf==nullptr){
    MIRANOTIFY(Notify::FATAL,"Could not open " << hashstatfilename << " for writing, is the disk full? Are permissions set right?");
  }

  auto mhsh=priv_writeHashStatFileHeader(gzf,HS_hs_basesperhash,HSSS_NOTSORTED,-1);

  ProgressIndicator<int32> P(0, static_cast<int32>(elementsperfile.size()));
  uint64 totalelements=0;
  {
    std::vector<hashstat_t> hashpool;
    hashpool.reserve(maxelementsperfile+10);

    for(size_t fi=0; fi< elementsperfile.size(); fi++){
      P.increaseprogress();

      CEBUG("Loading " << hashfilenames[fi] << endl);
      CEBUG("elements in file: " << elementsperfile[fi] << endl);

      if(elementsperfile[fi]==0) continue;
      hashpool.clear();
      hashpool.resize(elementsperfile[fi]);

      auto gzhf=gzopen(hashfilenames[fi].c_str(), "rb");
      if(gzhf==nullptr){
	MIRANOTIFY(Notify::FATAL,"Could not open " << hashfilenames[fi] << " for reading? It was written just moments ago, something with your machine is broken I think.");
      }

      auto readbytes=myGZRead(gzhf,&hashpool[0],sizeof(hashstat_t)*elementsperfile[fi]);
      gzclose(gzhf);
      if(readbytes != sizeof(hashstat_t)*elementsperfile[fi]) {
	MIRANOTIFY(Notify::FATAL, "Expected to read " << sizeof(hashstat_t)*elementsperfile[fi] << " bytes in file " << hashfilenames[fi] << " but read " << readbytes << ". Was the file deleted? Disk full?");
      }

      //for(size_t i=0; i<hashpool.size(); ++i){
      //  CEBUG(hashpool[i] << '\n');
      //}

      priv_compressHashStatBufferInPlace(hashpool,rarekmerearlykill);
      priv_calcFwdRevMinThresholdOKFlag(hashpool,fwdrevmin);

      totalelements+=hashpool.size();
      CEBUG("after comp: " << hashpool.size() << endl);

      if(!hashpool.empty()){
	auto writtenbytes=myGZWrite(gzf,
				    reinterpret_cast<const char *>(&hashpool[0]),
				    sizeof(hashstat_t)*hashpool.size());
	if(static_cast<size_t>(writtenbytes) != sizeof(hashstat_t)*hashpool.size()){
	  gzclose(gzf);
	  MIRANOTIFY(Notify::FATAL, "Could not save anymore the hash statistics (1). Disk full? Changed permissions?");
	}
      }
    }
  }

  gzclose(gzf);
  P.finishAtOnce();

  mhsh.numelem=totalelements;

  // final read to rewrite header with correct number of elements
  // and to sort the hashstatistics to be directly usable for making shortcuts

  // cannot use loadHashStatistics(), wrong header, need to go by foot
  CEBUG("opening gz" << endl);
  gzf=gzopen(hashstatfilename.c_str(),"r");
  if(gzf==nullptr){
    MIRANOTIFY(Notify::FATAL,"Could not open " << hashstatfilename << " for reading although it was written just moments ago??? Somethings is broken on your machine.");
  }

  CEBUG("loading header" << endl);
  {
    auto dummy=loadHashStatisticsFileHeader(gzf); // yes, we're throwing away this header
  }
  CEBUG("loading main statistics" << endl);
  loadHashStatistics(mhsh,gzf);
  gzclose(gzf);

  // calc some statistics
  CEBUG("some statistics" << endl);
  HS_hs_sortstatus=HSSS_NOTSORTED; // just to be sure
  priv_calcAvgHashFreq();

  // Now sort and save
  HS_hs_sortstatus=HSSS_NOTSORTED; // just to be sure
  CEBUG("sort low24" << endl);
  priv_sortLow24Bit();
  CEBUG("save statistics" << endl);
  saveHashStatistics(hashstatfilename)

  FUNCEND();
  return totalelements;
}
#define CEBUG(bla)



/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::saveHashStatistics(const std::string & filename)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::saveHashStatistics(const std::string & filename");

  auto gzf=gzopen(filename.c_str(),"wb1");
  if(gzf==nullptr){
    MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is the disk full? Are permissions set right?");
  }
  try{
    saveHashStatistics(gzf);
  }
  catch(Notify n){
    gzclose(gzf);
    cout << "Error for file " << filename << endl;
    n.handleError(THISFUNC);
  }
  gzclose(gzf);
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::saveHashStatistics(gzFile & gzf)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::saveHashStatistics(gzFile & gzf)");

  priv_saveHashVStatistics(gzf);

  //if(HSN_hsum_hashstats.empty()){
  //  saveHashVStatistics(ostr);
  //}else{
  //  saveHashMStatistics(ostr);
  //}
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_saveHashVStatistics(gzFile & gzf)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_saveHashVStatistics(gzFile & gzf)");

  auto mhs=priv_writeHashStatFileHeader(gzf,HS_hs_basesperhash,HS_hs_sortstatus,HS_hsv_hashstats.size());
  if(!HS_hsv_hashstats.empty()){
    auto writtenbytes=myGZWrite(gzf,
			      reinterpret_cast<const char *>(&HS_hsv_hashstats[0]),
			      sizeof(hashstat_t)*HS_hsv_hashstats.size());
    if(static_cast<size_t>(writtenbytes) != sizeof(hashstat_t)*HS_hsv_hashstats.size()){
      gzclose(gzf);
      MIRANOTIFY(Notify::FATAL, "Could not save anymore the hash statistics (1). Disk full? Changed permissions?");
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
typename HashStatistics<TVHASH_T>::mhsheader_t HashStatistics<TVHASH_T>::priv_writeHashStatFileHeader(std::ostream & ostr, uint32 basesperhash, uint8 sortstatus, uint64 numelem)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_writeHashStatFileHeader(std::ostream & ostr, uint32 basesperhash, uint8 sortstatus, uint64 numelem)");

  cout << "This hash save should not be called anymore, aborting.\n";
  exit(999);

  mhsheader_t mhsh;
  mhsh.version=3;
  mhsh.sortstatus=sortstatus;
  mhsh.basesperhash=basesperhash;
  mhsh.sizeofhash=sizeof(TVHASH_T);
  mhsh.numelem=numelem;
  ostr.write(reinterpret_cast<const char *>(&HS_hsfilemagic),4);
  ostr.write(reinterpret_cast<const char *>(&mhsh),sizeof(mhsh));
  CEBUG("Written MHSh " << mhsh << endl);
  return mhsh;
}

template<typename TVHASH_T>
typename HashStatistics<TVHASH_T>::mhsheader_t HashStatistics<TVHASH_T>::priv_writeHashStatFileHeader(gzFile & gzf, uint32 basesperhash, uint8 sortstatus, uint64 numelem)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_writeHashStatFileHeader(gzFile & gzf, uint32 basesperhash, uint8 sortstatus, uint64 numelem)");

  mhsheader_t mhsh;
  mhsh.version=4;
  mhsh.sortstatus=sortstatus;
  mhsh.basesperhash=basesperhash;
  mhsh.sizeofhash=sizeof(TVHASH_T);
  mhsh.numelem=numelem;
  mhsh.freq=HS_avg_freq;
  auto writtenbytes=myGZWrite(gzf,reinterpret_cast<const char *>(&HS_hsfilemagic),4);
  writtenbytes=myGZWrite(gzf,reinterpret_cast<const char *>(&mhsh),sizeof(mhsh));
  if(writtenbytes != sizeof(mhsh)) {
    MIRANOTIFY(Notify::FATAL,"Could not write header information. Is the disk full or quota reached? Changed access permissions?\n");
  }
  CEBUG("Written MHSh " << mhsh << endl);
  return mhsh;
}



/*************************************************************************
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
const typename HashStatistics<TVHASH_T>::mhsheader_t HashStatistics<TVHASH_T>::loadHashStatisticsFileHeader(const std::string & filename)
{
  FUNCSTART("const typename HashStatistics<TVHASH_T>::mhsheader_t HashStatistics<TVHASH_T>::loadHashStatisticsFileHeader(const std::string & fn)");
  mhsheader_t mhs;

  auto gzf=gzopen(filename.c_str(),"rb");
  if(gzf==nullptr){
    MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is it present? Are permissions set right?");
  }
  try{
    // 128k larger buffer to speed up decompression
    // TODO: test whether 64k or 16k is enough for good speedup
    gzbuffer(gzf,128*1024);
    mhs=loadHashStatisticsFileHeader(gzf);
  }
  catch(Notify n){
    gzclose(gzf);
    cout << "Error while loading file " << filename << endl;
    n.handleError(THISFUNC);
  }
  gzclose(gzf);
  return mhs;
}



/*************************************************************************
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
const typename HashStatistics<TVHASH_T>::mhsheader_t HashStatistics<TVHASH_T>::loadHashStatisticsFileHeader(gzFile & gzf)
{
  FUNCSTART("bool HashStatistics<TVHASH_T>::loadHashStatisticsFileHeader(gzFile & gzf)");

  mhsheader_t ret;

  auto localmagic=HS_hsfilemagic;
  auto readbytes=myGZRead(gzf,reinterpret_cast<char *>(&localmagic),4);
  if(readbytes==0) return ret;
  if(readbytes != 4
     || localmagic!=HS_hsfilemagic) {
    MIRANOTIFY(Notify::FATAL,"No magic found or truncated?\n");
  }
  readbytes=myGZRead(gzf,reinterpret_cast<char *>(&ret),sizeof(ret));

  if(readbytes != sizeof(ret)) {
    MIRANOTIFY(Notify::FATAL,"Not enough bytes read for header information. File truncated?\n");
  }

  CEBUG("Loaded MHSh " << ret << endl);

  if(ret.version!=4) {
    MIRANOTIFY(Notify::FATAL,"The file looks to be a MIRA HashStatistics file, but version " << static_cast<uint16>(ret.version) << " and not 3?\n");
  }

  return ret;
}
//#define CEBUG(bla)

/*************************************************************************
 *
 * Note: not using boost::iostreams with the gzlib decompressor as that thing
 *  does not automatically detect uncompressed files and gives back
 *  nonsense. I don't want to write a wrapper that does this detection and
 *  then sets up the streams as needed ... boost::iostreams should do that for
 *  me. *sigh*reinterpret_cast<char *>(&localmagic),4);
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::loadHashStatistics(const std::string & filename)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::loadHashStatistics(const std::string & filename)");

  auto gzf=gzopen(filename.c_str(),"rb");
  if(gzf==nullptr){
    MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is it present? Are permissions set right?");
  }
  try{
    // 128k larger buffer to speed up decompression
    // TODO: test whether 64k or 16k is enough for good speedup
    gzbuffer(gzf,128*1024);
    loadHashStatistics(gzf);
  }
  catch(Notify n){
    gzclose(gzf);
    cout << "Error while loading file " << filename << endl;
    n.handleError(THISFUNC);
  }
  gzclose(gzf);
}


/*************************************************************************
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::loadHashStatistics(gzFile & gzf)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::loadHashStatistics(gzFile & gzf)");

  auto mhsh=loadHashStatisticsFileHeader(gzf);
  loadHashStatistics(mhsh,gzf);
  HS_avg_freq=mhsh.freq;

  return;
}


/*************************************************************************
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::loadHashStatistics(const mhsheader_t & mhsh, gzFile & gzf)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::loadHashStatistics(const mhsheader_t & mhsh, gzFile & gzf)");

  CEBUG("Got MHSh " << mhsh << endl);

  if(mhsh.sizeofhash > sizeof(TVHASH_T)){
    MIRANOTIFY(Notify::FATAL,"Hash size " << mhsh.sizeofhash << " does not fit into currently used hash " << sizeof(TVHASH_T) << " ???");
  }
  if(HS_hs_basesperhash != 0 && !HS_hsv_hashstats.empty() && mhsh.basesperhash != HS_hs_basesperhash){
    MIRANOTIFY(Notify::FATAL,"Current hashstat kmer size is " << HS_hs_basesperhash
	       << ", but kmer size in data to load is " << mhsh.basesperhash
	       << " ???\n");
   }else{
    HS_hs_basesperhash=mhsh.basesperhash;
  }

  if(!HS_hsv_hashstats.empty()){
//    HS_hs_sortstatus=HSSS_NOTSORTED;
//    HS_hs_needsconsolidation=true;
    MIRANOTIFY(Notify::FATAL,"Appending to existing hashstat not implemented yet\n");
  }else{
    HS_hs_sortstatus=mhsh.sortstatus;
  }

  HS_avg_freq=mhsh.freq;

  if(mhsh.numelem==0) return;

  if(mhsh.numelem){
//   HS_avg_freq.isvalid=false;
    HS_avg_freq.corrected=0;
    HS_avg_freq.raw=0;
    HS_avg_freq.taken=0;
    HS_hsv_hsshortcuts.clear();

    auto oldsize=HS_hsv_hashstats.size();
    CEBUG("Will resize to " << oldsize+mhsh.numelem << endl);
    HS_hsv_hashstats.resize(oldsize+mhsh.numelem);
    auto readbytes=myGZRead(gzf,reinterpret_cast<char *>(&HS_hsv_hashstats[oldsize]),mhsh.numelem*sizeof(hashstat_t));
    if(readbytes != mhsh.numelem*sizeof(hashstat_t)){
      MIRANOTIFY(Notify::FATAL,"Expected to read " << mhsh.numelem*sizeof(hashstat_t) << " bytes, but got " << readbytes << endl);
    }
  }
}


template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_calcAvgHashFreq(bool verbose)
{
  HS_avg_freq=avg_freq_t();
  priv_sortByCountUp();

  if(HS_hsv_hashstats.empty()) return;

  HS_avg_freq.raw=priv_calcMidHashStatIndex(HS_hsv_hashstats,0);
  HS_avg_freq.corrected=HS_avg_freq.raw;

  auto hsthreshold=HS_hsv_hashstats.size()-HS_hsv_hashstats.size()/10;

  if(verbose){
    cout << "Raw MHI: " << HS_avg_freq.raw << endl;
    cout << "Raw avg. freq. : " << HS_hsv_hashstats[HS_avg_freq.raw].hsc.fcount << " " << HS_hsv_hashstats[HS_avg_freq.raw].hsc.rcount << endl;
    cout << "HSS " << HS_hsv_hashstats.size() << "\tHSST: " << hsthreshold << endl;
  }


  //// if mh index is in last 10 % of the hashstats, we have a pretty skewed
  ////  distribution. In that case, recalc without last 10%
  //// TODO: check whether 40 or 50% wouldn't be better.
  if(HS_avg_freq.corrected >= hsthreshold){
    HS_avg_freq.corrected=priv_calcMidHashStatIndex(HS_hsv_hashstats,10);
    if(verbose){
      cout << "Corrected MHI: " << HS_avg_freq.corrected << endl;
      cout << "Corrected avg. freq. : " << HS_hsv_hashstats[HS_avg_freq.corrected].hsc.fcount << " " << HS_hsv_hashstats[HS_avg_freq.corrected].hsc.rcount << endl;
      cout << "HSS " << HS_hsv_hashstats.size() << "\tHSST: " << (HS_hsv_hashstats.size()-HS_hsv_hashstats.size()/10) << endl;
    }
  }

  HS_avg_freq.corrected=HS_hsv_hashstats[HS_avg_freq.corrected].hsc.getCount();
  HS_avg_freq.raw=HS_hsv_hashstats[HS_avg_freq.raw].hsc.getCount();

  HS_avg_freq.taken=HS_avg_freq.corrected;
  if(HS_avg_freq.taken < HS_avg_freq.min){
    HS_avg_freq.taken=HS_avg_freq.min;
    if(verbose){
      cout << "Forced avg. freq: " << HS_avg_freq.taken << endl;
    }
  }

  FUNCEND();
  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
size_t HashStatistics<TVHASH_T>::priv_calcMidHashStatIndex(const std::vector<hashstat_t> & hashstats, size_t dontcarepercent)
{
  FUNCSTART("size_t HashStatistics<TVHASH_T>::priv_calcMidHashStatIndex(const std::vector<hashstat_t> & hashstats, size_t dontcarepercent)");

  if(hashstats.empty()) return 0;

  size_t firsti=0;
  size_t lasti=hashstats.size();
  if(dontcarepercent){
    firsti=hashstats.size()*dontcarepercent/100;
    lasti-=hashstats.size()*dontcarepercent/100;
  }else{
    // 5% default
    firsti=hashstats.size()/20;
    lasti-=hashstats.size()/20;
  }

  size_t sumhashcounts=0;
  uint32 oldhashcount=hashstats[0].hsc.getCount()-1;
  size_t oldsumhashcounts=0;
  for(size_t i=firsti; i<lasti; i++){
    if(hashstats[i].hsc.getCount() != oldhashcount){
      BUGIFTHROW(oldhashcount>hashstats[i].hsc.getCount(),"haststat array not sorted by count???");
      oldhashcount=hashstats[i].hsc.getCount();
      CEBUG("count: " << oldhashcount << "\tsumhash: " << sumhashcounts << "\tdiff: " << sumhashcounts-oldsumhashcounts << endl);
      oldsumhashcounts=sumhashcounts;
    }
    if(hashstats[i].hsc.hasfwdrevthresholdok) sumhashcounts+=hashstats[i].hsc.getCount();
  }
  CEBUG("count: " << oldhashcount << "\tsumhash: " << sumhashcounts << endl);

  // Hmmm, pathological case. Maybe all reads were in the same direction.
  //  simply recalc without the "has fwd/rev" clause
  bool dontusefwdrev=false;
  if(sumhashcounts==0){
    dontusefwdrev=true;
    for(size_t i=firsti; i<lasti; i++){
      sumhashcounts+=hashstats[i].hsc.getCount();
    }
    CEBUG("recalc sumhash: " << sumhashcounts << endl);
  }

  size_t midhashstats=sumhashcounts/2;

  CEBUG("midhashstats: " << midhashstats << endl);

  sumhashcounts=0;
  for(size_t i=firsti; i<lasti; i++){
    if(dontusefwdrev || hashstats[i].hsc.hasfwdrevthresholdok) sumhashcounts+=hashstats[i].hsc.getCount();
    if(sumhashcounts>midhashstats) {
      return i;
    }
  }

  FUNCEND();

  return 0;
}
//#define CEBUG(bla)





/*************************************************************************
 *
 * Needs:
 *  - the hash statistics vector (sorted by count)
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::showHashStatisticsInfo()
{
  FUNCSTART("void HashStatistics<TVHASH_T>::showHashStatisticsInfo()");

  cout << "Kmer statistics:\n"
       << "=========================================================\n"
       << "Using kmers of size: " << HS_hs_basesperhash << endl
       << "Measured avg. raw frequency coverage: " << HS_avg_freq.raw << endl
       << "Corrected avg. raw frequency coverage: " << HS_avg_freq.corrected;

  if(HS_avg_freq.raw!=HS_avg_freq.corrected){
    cout << "\tSKEWED DISTRIBUTION!";
  }
  cout << '\n';

  if(HS_avg_freq.corrected<HS_avg_freq.min){
    cout << "Forced minimum average frequency: " << HS_avg_freq.min << endl;
  }

  cout << "\nFinal average frequency: " << HS_avg_freq.taken << endl;


  cout << "\nDeduced thresholds:\n"
       << "-------------------"
       << "\nRare freq: " << HS_freqest_minnormal*HS_avg_freq.taken
       << "\nMin normal freq: " << HS_freqest_minnormal*HS_avg_freq.taken
       << "\nMax normal freq " << HS_freqest_maxnormal*HS_avg_freq.taken
       << "\nRepeat freq: " << HS_freqest_repeat*HS_avg_freq.taken
       << "\nHeavy freq: " << HS_freqest_heavyrepeat*HS_avg_freq.taken
       << "\nCrazy freq: " << HS_freqest_crazyrepeat*HS_avg_freq.taken
       << "\nMask freq: " << HS_nastyrepeatratio*HS_avg_freq.taken
       << "\n\nRepeat ratio histogram:\n"
       << "-----------------------"
       << endl;

  std::vector<size_t> ratiocounts;
  ratiocounts.reserve(8192);
  for(size_t i=0; i<HS_hsv_hashstats.size(); i++){
    uint32 rci=static_cast<uint32>((static_cast<double>(HS_hsv_hashstats[i].hsc.getCount()) / HS_avg_freq.taken) + 0.5);
    if(rci>=ratiocounts.size()){
      ratiocounts.resize(rci+1,0);
    }
    ratiocounts[rci]++;
  }

  for(size_t i=0; i<ratiocounts.size(); i++){
    if(ratiocounts[i]) cout << i << '\t' << ratiocounts[i] << endl;
  }

  cout << "=========================================================\n\n";

  FUNCEND();

  return;
}





/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::trimHashStatsByFrequencyAND(uint32 minfwd, uint32 minrev, uint32 mintotal)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  priv_trimHashVStatsByFrequencyAND(minfwd,minrev,mintotal);

  //if(HSN_hsum_hashstats.empty()){
  //  trimHashVStatsByFrequency(minfwd,minrev,mintotal);
  //}else{
  //  trimHashMStatsByFrequency(minfwd,minrev,mintotal);
  //}

}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::trimHashStatsByFrequencyANDOR(uint32 minfwd, uint32 minrev, uint32 mintotal)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  priv_trimHashVStatsByFrequencyANDOR(minfwd,minrev,mintotal);

  //if(HSN_hsum_hashstats.empty()){
  //  trimHashVStatsByFrequency(minfwd,minrev,mintotal);
  //}else{
  //  trimHashMStatsByFrequency(minfwd,minrev,mintotal);
  //}

}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_trimHashVStatsByFrequencyAND(uint32 minfwd, uint32 minrev, uint32 mintotal)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_trimHashVStatsByFrequencyAND(int32 minfwd, int32 minrev, int32 mintotal)");

  auto srcI=HS_hsv_hashstats.begin();
  auto dstI=srcI;

  for(; srcI!=HS_hsv_hashstats.end(); ++srcI){
    bool ok=true;
    if(srcI->hsc.fcount<minfwd
       || srcI->hsc.rcount<minrev
       || srcI->hsc.fcount+srcI->hsc.rcount < mintotal){
      ok=false;
      CEBUG("rm\t");
    }else{
      CEBUG("keep\t");
    }
    CEBUG(srcI-HS_hsv_hashstats.begin() << "\t" << *srcI << endl);
    *dstI=*srcI;
    if(ok)++dstI;
  }
  HS_hsv_hashstats.resize(dstI-HS_hsv_hashstats.begin());
  HS_hsv_hsshortcuts.clear();
  //HS_hs_dist.clear();
}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_trimHashVStatsByFrequencyANDOR(uint32 minfwd, uint32 minrev, uint32 mintotal)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_trimHashVStatsByFrequencyANDOR(int32 minfwd, int32 minrev, int32 mintotal)");

  auto srcI=HS_hsv_hashstats.begin();
  auto dstI=srcI;

  for(; srcI!=HS_hsv_hashstats.end(); ++srcI){
    bool ok=false;
    if((srcI->hsc.fcount>=minfwd
	&& srcI->hsc.rcount>=minrev)
       || srcI->hsc.fcount+srcI->hsc.rcount >= mintotal){
      ok=true;
      CEBUG("keep\t");
    }else{
      CEBUG("rm\t");
    }
    CEBUG(srcI-HS_hsv_hashstats.begin() << "\t" << *srcI << endl);
    *dstI=*srcI;
    if(ok)++dstI;
  }
  HS_hsv_hashstats.resize(dstI-HS_hsv_hashstats.begin());
  HS_hsv_hsshortcuts.clear();
  //HS_hs_dist.clear();
}



/*************************************************************************
 *
 * This is a quick hack for kmer dust ends
 *
 * kmers with microrepeats have their counts set to 0, then removed
 *  from hashstats
 *
 * replen: length of microrepeat (0..4, 0==all)
 * Ratio()
 *   perc: total microrepeat length, as percentage of kmer
 * Fixed()
 *   numbases: total microrepeat length, as num bases in kmer
 *
 *
 * e.g: bph=17, replen=0, perc=67
 *  67% of 17 = 11.39 --> 12 bases
 *
 * kmers ending in
 *  - 12 same bases    (AAAAAAAAAAAAxxxxx & xxxxxAAAAAAAAAAAA)
 *  - 6 * 2 same bases (ACACACACACACxxxxx ...)
 *  - 4 * 3 same bases (ACTACTACTACTxxxxx ...)
 *  - 3 * 4 same bases (ACGTACGTACGTxxxxx ...)
 *
 * For 15 bp:
 *  - AAAAAAAAAAAAAAAxx ...
 *  - ACACACACACACACAxx & xxACACACACACACACA
 *  - ACTACTACTACTACTxx ...
 *  - ACGTACGTACGTACGxx ...
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::removeDustEndsRatio(uint8 replen, uint8 perc)
{
  uint32 numbases=HS_hs_basesperhash*perc/100;
  if((numbases==0 || HS_hs_basesperhash*perc)%100!=0) ++numbases;
  removeDustEndsFixed(replen,numbases);
}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::removeDustEndsFixed(uint8 replen, uint32 numbases)
{
  const static std::string acgt="ACGT";

  if(HS_hsv_hashstats.empty()) return;
  if(numbases>HS_hs_basesperhash) numbases=HS_hs_basesperhash;

  TVHASH_T kmerdust;
  TVHASH_T dustadd;

  std::vector<TVHASH_T> kdv1;
  std::vector<TVHASH_T> kdv2;

  TVHASH_T rightmask(1);
  if(numbases>=sizeof(TVHASH_T)*4){
    rightmask=0;
  }else{
    rightmask<<=numbases*2;
  }
  --rightmask;
  uint32 leftshift=(HS_hs_basesperhash-numbases)*2; // numbases always >= 1
  TVHASH_T leftmask(rightmask);
  leftmask<<=leftshift;

  if(replen==0 || replen==1){
    uint32 repeats=numbases;
    for(uint16 x0=0;x0<4;++x0){
      kmerdust=0;
      for(uint32 ci=0; ci<repeats; ++ci) {
	kmerdust<<=2;
	kmerdust|=seqtohash::hashaddmatrix[acgt[x0]]-1;
      }
      kmerdust&=rightmask; // depending on repeat me may have too much in our kmer!
      kdv1.push_back(kmerdust);
      kmerdust<<=leftshift;
      kdv2.push_back(kmerdust);
    }
    priv_rde_helper1_set2zero(kdv1,rightmask);
    priv_rde_helper1_set2zero(kdv2,leftmask);
  }

  if(replen==0 || replen==2){
    uint32 repeats=numbases/2;
    if(numbases%2) ++repeats;
    kdv1.clear();
    kdv2.clear();
    for(uint16 x0=0;x0<1;++x0){
      for(uint16 x1=0;x1<4;++x1){
	kmerdust=0;
	dustadd=seqtohash::hashaddmatrix[acgt[x0]]-1;
	dustadd<<=2;
	dustadd|=seqtohash::hashaddmatrix[acgt[x1]]-1;
	for(uint32 ci=0; ci<repeats; ++ci) {
	  kmerdust<<=4;
	  kmerdust|=dustadd;
	}
	kmerdust&=rightmask; // depending on repeat me may have too much in our kmer!
	kdv1.push_back(kmerdust);
	kmerdust<<=leftshift;
	kdv2.push_back(kmerdust);
      }
    }
    priv_rde_helper1_set2zero(kdv1,rightmask);
    priv_rde_helper1_set2zero(kdv2,leftmask);
  }

  if(replen==0 || replen==3){
    uint32 repeats=numbases/3;
    if(numbases%3) ++repeats;
    kdv1.clear();
    kdv2.clear();
    for(uint16 x0=0;x0<4;++x0){
      for(uint16 x1=0;x1<4;++x1){
	for(uint16 x2=0;x2<4;++x2){
	  kmerdust=0;
	  dustadd=seqtohash::hashaddmatrix[acgt[x0]]-1;
	  dustadd<<=2;
	  dustadd|=seqtohash::hashaddmatrix[acgt[x1]]-1;
	  dustadd<<=2;
	  dustadd|=seqtohash::hashaddmatrix[acgt[x2]]-1;
	  for(uint32 ci=0; ci<repeats; ++ci) {
	    kmerdust<<=6;
	    kmerdust|=dustadd;
	  }
	  kmerdust&=rightmask; // depending on repeat me may have too much in our kmer!
	  kdv1.push_back(kmerdust);
	  kmerdust<<=leftshift;
	  kdv2.push_back(kmerdust);
	}
      }
    }
    priv_rde_helper1_set2zero(kdv1,rightmask);
    priv_rde_helper1_set2zero(kdv2,leftmask);
  }

  // Yes, I've seen 4 bp repeats in SILVA data
  if(replen==0 || replen==4){
    uint32 repeats=numbases/4;
    if(numbases%4) ++repeats;
    kdv1.clear();
    kdv2.clear();
    for(uint16 x0=0;x0<4;++x0){
      for(uint16 x1=0;x1<4;++x1){
	for(uint16 x2=0;x2<4;++x2){
	  for(uint16 x3=0;x3<4;++x3){
	    kmerdust=0;
	    dustadd=seqtohash::hashaddmatrix[acgt[x0]]-1;
	    dustadd<<=2;
	    dustadd|=seqtohash::hashaddmatrix[acgt[x1]]-1;
	    dustadd<<=2;
	    dustadd|=seqtohash::hashaddmatrix[acgt[x2]]-1;
	    dustadd<<=2;
	    dustadd|=seqtohash::hashaddmatrix[acgt[x3]]-1;
	    for(uint32 ci=0; ci<repeats; ++ci) {
	      kmerdust<<=8;
	      kmerdust|=dustadd;
	    }
	    kmerdust&=rightmask; // depending on repeat me may have too much in our kmer!
	    kdv1.push_back(kmerdust);
	    kmerdust<<=leftshift;
	    kdv2.push_back(kmerdust);
	  }
	}
      }
    }
    priv_rde_helper1_set2zero(kdv1,rightmask);
    priv_rde_helper1_set2zero(kdv2,leftmask);
  }

  trimHashStatsByFrequencyAND(0,0,1);
}


// helper for above
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_rde_helper1_set2zero(const std::vector<TVHASH_T> & kdv, TVHASH_T & andmask)
{
/*
  std::string tmpstr;
  hash2string(andmask,HS_hs_basesperhash,tmpstr);
  cout << "andmask " << std::hex << andmask << "\t" << tmpstr << endl;
  for(auto & kdve : kdv){
    hash2string(kdve,HS_hs_basesperhash,tmpstr);
    cout << "kmerdustval " << kdve << "\t" << tmpstr << endl;
  }
  cout << std::dec << endl;
*/

  for(auto & hse : HS_hsv_hashstats) {
    auto andval=(hse.vhash & andmask);
    for(auto & kdve : kdv){
      if(unlikely(andval == kdve)){
	hse.hsc.fcount=0;
	hse.hsc.rcount=0;
	break;
      }
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

std::string laberbla;

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::calcKMerForks(uint32 mincount, bool needfwdrev)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::calcKMerForks(uint32 mincount, bool needfwdrev)");

  if(HS_hsv_hashstats.empty()) return;

  for(auto & hse : HS_hsv_hashstats) {
    hse.hsc.iskmerforkf=false;
    hse.hsc.iskmerforkr=false;
  }

  HS_vhashmask=1;
  // *grml* undefined behaviour of left shift for 64 shifts in a 64 bit type makes this cludge necessary
  // the same for 32 shift in 32 bit types etc.pp
  if(HS_hs_basesperhash>=sizeof(TVHASH_T)*4){
    HS_vhashmask=0;
  }else{
    auto rollbases=HS_hs_basesperhash-1;
    HS_vhashmask<<=(rollbases*2);
  }
  // vhashmask is now, e.g. for bph=31, 00010000000....
  --HS_vhashmask;
  // vhashmask is now, e.g. for bph=31, 000011111....

  // calc the status on ?..............
  priv_sortMaskUp();
  laberbla="rev ";
  CEBUG("KMERFORK REV " << hex << HS_vhashmask << dec << endl);
  if(needfwdrev){
    priv_ckmf_helper(HS_vhashmask,mincount);
  }else{
    priv_ckmf_relaxed_helper(HS_vhashmask,mincount,false);
  }

  // calc the status on ..............?
  HS_vhashmask<<=2;
  // vhashmask is now, e.g. for bph=31, 0011111....00
  priv_sortMaskUp();

  laberbla="fwd ";
  CEBUG("KMERFORK FWD " << hex << HS_vhashmask << dec << endl);
  if(needfwdrev){
    priv_ckmf_helper(HS_vhashmask,mincount);
  }else{
    priv_ckmf_relaxed_helper(HS_vhashmask,mincount,true);
  }

  priv_sortLow24Bit();

  return;
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_ckmf_helper(TVHASH_T hashmask, uint32 mincount)
{
  auto hsI=HS_hsv_hashstats.begin();
  auto hsJ=hsI+1;

  CEBUG("KMER DUMP:\n");

  for(; hsJ!=HS_hsv_hashstats.end(); ++hsI,++hsJ){
    CEBUG(" b " << (hsI->hsc.hasfwdrevthresholdok));
    CEBUG(" c " << (hsJ->hsc.hasfwdrevthresholdok));
    CEBUG(" d " << ((hsI->vhash&hashmask) == (hsJ->vhash&hashmask)));
    CEBUG(" e " << (hsI->vhash != hsJ->vhash));
    CEBUG(" f " << (hsI->hsc.fcount >= mincount));
    CEBUG(" g " << (hsI->hsc.rcount >= mincount));
    CEBUG(" h " << (hsJ->hsc.fcount >= mincount));
    CEBUG(" i " << (hsJ->hsc.rcount >= mincount));
    CEBUG(endl);
    if(hsI->hsc.hasfwdrevthresholdok
       && hsJ->hsc.hasfwdrevthresholdok
       && (hsI->vhash&hashmask) == (hsJ->vhash&hashmask)
       && hsI->vhash != hsJ->vhash
       && hsI->hsc.fcount >= mincount
       && hsI->hsc.rcount >= mincount
       && hsJ->hsc.fcount >= mincount
       && hsJ->hsc.rcount >= mincount){
      CEBUG("ISKMER!\n");
      hsI->hsc.iskmerforkf=true;
      hsI->hsc.iskmerforkr=true;
      hsJ->hsc.iskmerforkf=true;
      hsJ->hsc.iskmerforkr=true;
    }
    CEBUG(laberbla << hash2string(hsI->vhash,HS_hs_basesperhash) << ' ' << *hsI << '\n');
    CEBUG(laberbla << hash2string(hsJ->vhash,HS_hs_basesperhash) << ' ' << *hsJ << '\n');
  }
}


template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_ckmf_relaxed_helper(TVHASH_T hashmask, uint32 mincount, bool isfwd)
{
  auto hsI=HS_hsv_hashstats.begin();

  CEBUG("KMER RELAX:\n");

  std::vector<uint8> iskmer;
  while(hsI!=HS_hsv_hashstats.end()){
    // find range with same subhash
    auto hsE=hsI;
    for(;hsE!=HS_hsv_hashstats.end() && (hsE->vhash & hashmask) == (hsI->vhash & hashmask); ++hsE) {}
    if(hsE-hsI > 1){
      // check for strong kmers
      iskmer.clear(); iskmer.resize(hsE-hsI,0);
      uint32 count=0;
      uint32 runi=0;
      uint32 maxcount=0;
      for(auto rI=hsI; rI != hsE; ++rI, ++runi){
	maxcount=std::max(maxcount,rI->hsc.getCount());
	if(rI->hsc.fcount>2 && rI->hsc.rcount>2){
	  ++count;
	  iskmer[runi]=1;
	}
      }
      // if no strong, relax a bit: the ones with the highest count is defined as strong kmer
      if(count==0){
	runi=0;
	for(auto rI=hsI; rI != hsE; ++rI, ++runi){
	  if(rI->hsc.getCount() == maxcount) iskmer[runi]=1;
	}
      }

      // now, as this routine is specifically designed for proposed end clip,
      // set the iskmerfork flag for all non strong kmers, but do not set it for
      // the strong.
      // seems a bit illogic, but is really made to fit proposed end clipping
      //  who stops "saving" when detecting a kmerfork. So the strong, good
      //  will be saved while weak ones wont.

      CEBUG("KMER DECIDE " << count << " " << maxcount << endl);

      runi=0;
      for(auto rI=hsI; rI != hsE; ++rI, ++runi){
	if(!iskmer[runi]) {
	  if(isfwd) {
	    rI->hsc.iskmerforkf=true;
	  }else{
	    rI->hsc.iskmerforkr=true;
	  }
	}
	CEBUG(laberbla << hash2string(rI->vhash,HS_hs_basesperhash) << ' ' << *rI << '\n');
      }
    }

    hsI=hsE;
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * needs:
 *  - hashstats filled with entries (can be unsorted, will be re-sorted
 *    anyway)
 *
 * returns:
 *  - hashstats array sorted by low 24 bit (low to high), then by vhash
 *  - hsshortcuts_begin and ..._end pointing to start and end of each
 *    low 24 bit group of same value
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_makeHashStatArrayShortcuts()
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_makeHashStatArrayShortcuts()");

  CEBUG("makeHashStatArrayShortcuts: basesperhash: " << HS_hs_basesperhash << "\n");

  BUGIFTHROW(HS_hs_basesperhash==0, "HS_hs_basesperhash == 0 ???");

  for(size_t i=0; i<HS_hsv_hashstats.size(); i++){
    CEBUG(HS_hsv_hashstats[i] << '\n');
  }

  //cout << "Going on a sort. "; dateStamp(cout);

  //sort(HS_hsv_hashstats.begin(), HS_hsv_hashstats.end(), sortHashStatComparatorByLow24bit);
  priv_sortLow24Bit();

  //cout << "Sort ended. "; dateStamp(cout);

  HS_hsv_hsshortcuts.clear();
  {
    hsvbendit_t tmpb;
    tmpb.b=HS_hsv_hashstats.end();
    tmpb.e=HS_hsv_hashstats.end();
    HS_hsv_hsshortcuts.resize(
      1<<(std::min(static_cast<uint32>(12),HS_hs_basesperhash)*2),
      tmpb
      );
  }

  CEBUG("HS_hsv_hsshortcuts.size(): " << HS_hsv_hsshortcuts.size() << endl);

  auto hsI=HS_hsv_hashstats.begin();
  if(hsI==HS_hsv_hashstats.end()) return;


  TVHASH_T acthash= (hsI->vhash & HS_MAXVHASHMASK);
  while(hsI != HS_hsv_hashstats.end()){
    CEBUG("begin " << hex << acthash << dec << " is: " << *hsI << endl);
    HS_hsv_hsshortcuts[static_cast<uint64>(acthash)].b=hsI;
    for(;(hsI != HS_hsv_hashstats.end()) && ((hsI->vhash & HS_MAXVHASHMASK) == acthash); hsI++) {
      CEBUG("INC\n")
    }
    CEBUG("end " << hex << acthash << dec << " is: " << *hsI << endl);
    HS_hsv_hsshortcuts[static_cast<uint64>(acthash)].e=hsI;
    //cout << "vhash: " << hex << acthash << "\t" << dec << HS_hsv_hsshortcuts_end[acthash]-HS_hsv_hsshortcuts_begin[acthash] << '\n';
    if(hsI != HS_hsv_hashstats.end()) acthash= hsI->vhash & HS_MAXVHASHMASK;
  }

  //cout << "Done making shortcuts. " << endl; dateStamp(cout);

  FUNCEND();
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

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::assignReadBaseStatistics_MultiThread(ReadPool & rp, uint32 numthreads, bool masknastyrepeats, bool calckmerforks, uint32 mincountkmerforks, bool needfwdrev)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::assignReadBaseStatistics_MultiThread(ReadPool & rp, uint32 numthreads, bool masknastyrepeats, uint32 mincountkmerforks, bool needfwdrev)");

  if(calckmerforks){
    calcKMerForks(mincountkmerforks, needfwdrev);
    // calcKMerForks() uses the hashstat sort function which may have emptied shortcuts
    // redo them just to be sure
    priv_makeHashStatArrayShortcuts();
  }else if(HS_hsv_hsshortcuts.empty()){
    priv_makeHashStatArrayShortcuts();
  }

  if(HS_hsv_hsshortcuts.empty()){
    MIRANOTIFY(Notify::FATAL,"Tried to assign base statistics (1) though there are no statistics which can be assigned? Something's wrong with your data set.");
  }

  arbs_threadsharecontrol_t atsc;

  atsc.from=0;
  atsc.to=rp.size();
  atsc.todo=0;
  atsc.done=0;
  atsc.stepping=1000;

  // TODO: unneeded now as working on HS_* variables, reorganise
  // vvvvvvvvvvvvvvv
  atsc.rpptr=&rp;
  atsc.avghashcov=HS_avg_freq.taken;
  atsc.hashstatsptr=&HS_hsv_hashstats;
  atsc.basesperhash=HS_hs_basesperhash;
  atsc.hsscptr=&HS_hsv_hsshortcuts;
  // ^^^^^^^^^

  atsc.masknastyrepeats=masknastyrepeats;
  atsc.truekmerforks=!needfwdrev;

  CEBUG("minnormalhashcov: " << atsc.avghashcov << endl);


  //uint32 numthreads=8;
  boost::thread_group workerthreads;
  for(uint32 ti=0; ti<numthreads;++ti){
    workerthreads.create_thread(boost::bind(&HashStatistics<TVHASH_T>::priv_arb_thread, this, ti, &atsc));
  }

  ProgressIndicator<int64> pi(0,rp.size());
  while(atsc.done!=rp.size()){
    pi.progress(atsc.done);
    sleep(1);
  }
  pi.finishAtOnce(cout);

  // they normally should all have exited at this point, but be nice and play by the rules
  workerthreads.join_all();

}
//#define CEBUG(bla)

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_arb_thread(uint32 threadnum, arbs_threadsharecontrol_t * tscptr)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_arb_thread(uint32 threadnum, arbs_threadsharecontrol_t * tscptr)");

  try{
    int32 from;
    int32 to;
    while(true){
      {
	boost::mutex::scoped_lock lock(tscptr->accessmutex);
	if(tscptr->todo >= tscptr->to) break;
	from=tscptr->todo;
	tscptr->todo+=tscptr->stepping;
	if(tscptr->todo > tscptr->to) tscptr->todo = tscptr->to;
	to=tscptr->todo;
      }
      priv_arb_DoStuff(
	*(tscptr->rpptr),
	tscptr->avghashcov,
	*(tscptr->hashstatsptr),
	tscptr->basesperhash,
	*(tscptr->hsscptr),
	tscptr->masknastyrepeats,
	from,
	to,
	tscptr->truekmerforks
	);
      {
	boost::mutex::scoped_lock lock(tscptr->accessmutex);
	tscptr->done+=tscptr->stepping;
	if(tscptr->done > tscptr->to) tscptr->done=tscptr->to;
      }
    }
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
}


/*************************************************************************
 *
 * truekmerforks = false; traditional 'assembly' style fwd/rev kmers, for
 *             being cautious in assembly
 * truekmerforks = true; only real kmers in given direction, for pec clipping
 *
 * false
 * aaaaaaaaaaaaaaaaRRRRRRRRRRRRRRRRRRRRbbbbbbbbbbbbbbbbbb
 *               <KKKKKKKK>
 *                             <KKKKKKKK>
 *
 * true
 * aaaaaaaaaaaaaaaaRRRRRRRRRRRRRRRRRRRRbbbbbbbbbbbbbbbbbb
 *               <KKKKKKKK
 *                              KKKKKKKK>
 *
 * The 'true' needs more time as one needs to really go through the reverse
 *  sequence and analyse the kmers for real
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {if(docebug) {cout << bla; cout.flush();}}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_arb_DoStuff(ReadPool & rp, size_t avghashcov, std::vector<hashstat_t> & hashstats, const uint32 basesperhash, std::vector<hsvbendit_t> & hsshortcuts, bool masknastyrepeats, int32 fromid, int32 toid, bool truekmerforks)
{
  FUNCSTART("HashStatistics<TVHASH_T>::priv_arb_DoStuff(ReadPool & rp, size_t avghashcov, std::vector<hashstat_t> & hashstats, const uint32 basesperhash, std::vector<hsvbendit_t> & hsshortcuts, bool masknastyrepeats, int32 fromid, int32 toid, bool truekmerforks)");

  //bool docebug=false;

  BUGIFTHROW(hsshortcuts.empty(),"hsshortcuts.empty() ??? at " << &hsshortcuts);

  auto minnormalhashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_minnormal);
  auto maxnormalhashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_maxnormal);
  auto repeathashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_repeat);
  auto heavyrepthashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_heavyrepeat);
  auto crazyrepthashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_crazyrepeat);
  uint32 maskhashcov=0;

  if(masknastyrepeats){
    maskhashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_nastyrepeatratio);
    if(HS_nastyrepeatcoverage>0 && HS_nastyrepeatcoverage<maskhashcov){
      maskhashcov=HS_nastyrepeatcoverage;
    }
  }

  std::vector<vhrap_t> singlereadvhraparray;
  singlereadvhraparray.reserve(10000);

  // we will not use a mask, but
  //  we need to supply an empty one anyway
  std::vector<uint8> tagmaskvector;

  // stores in each read whether the given hash frequency was seen
  std::vector<uint8> hasfrequency(8);

  std::vector<uint8> mcmask;
  mcmask.reserve(10000);

  multitag_t tmpmt(Read::REA_defaulttag_MNRr);

  CEBUG("dostuff from " << fromid << " to " << toid << endl);
  Read::setCoutType(Read::AS_TEXT);

  for(int32 actreadid=fromid; actreadid<toid; ++actreadid){
    //if(actreadid>100) return;

    Read & actread= rp.getRead(actreadid);

    //if(actread.getName()=="mictbac_bg:1:2101:20987:13472"){
    //  docebug=true;
    //}else{
    //  docebug=false;
    //}

    CEBUG("dsloop " << fromid << " " << actread.getName() << " " << actread.getLenClippedSeq() << endl);

    // get rid of old values
    actread.clearAllBPosHashStats();
    actread.setHasFreqAvg(false);
    actread.setHasFreqRept(false);
    actread.setHasKMerFork(false);
    actread.deleteTag(tmpmt.identifier);

    // whatever happens: this read was looked upon by this routine, so technically we "have" base hashstats
    actread.setHasBaseHashStats(true);


    if(!actread.hasValidData()
      || !actread.isUsedInAssembly()) continue;

//#define CEBUG(bla)   {if(cebugok) cout << bla; cout.flush();}
//    bool cebugok=false;
//    if(actread.getName()=="E0K6C4E01CTNQI") cebugok=true;

    uint32 slen=actread.getLenClippedSeq();

    if(slen<basesperhash) continue;

    if(masknastyrepeats) {
      mcmask.clear();
      mcmask.resize(actread.getLenSeq(),0);
    }

    mstd::fill(hasfrequency,0);

    CEBUG("### Before ...\n" << actread << endl);

    singlereadvhraparray.resize(slen);
    // BaCh 23.05.2015: don't need that. Skim::transformSeqToVariableHash() tests for this
    // tagmaskvector.resize(slen,0);

    auto srvaI=singlereadvhraparray.begin();

    std::vector<Read::bposhashstat_t> & bposhashstats=const_cast<std::vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
    uint32 hashesmade;

    bool haskmerfork=false;

    // fwd
    if(true){
      {
	int32 bfpos=actread.calcClippedPos2RawPos(0);
	int32 bfposinc=1;

	hashesmade=Skim<TVHASH_T>::transformSeqToVariableHash(
	  actreadid,
	  actread,
	  actread.getClippedSeqAsChar(),
	  slen,
	  basesperhash,
	  srvaI,
	  false,
	  1,
	  tagmaskvector,
	  bposhashstats,
	  bfpos,
	  bfposinc
	  );
      }
      singlereadvhraparray.resize(hashesmade);

      CEBUG("hashesmade: " << hashesmade << endl);
      CEBUG("maskhashcov: " << maskhashcov << endl);

      typename std::vector<hashstat_t>::const_iterator lowerbound;

      typename std::vector<hashstat_t>::const_iterator hssearchI;
      srvaI=singlereadvhraparray.begin();

      int32 bfpos1,bfpos2;
      hashstat_t hstmp;
      bool foundit=false;
      for(; srvaI != singlereadvhraparray.end(); srvaI++){
	CEBUG(*srvaI << '\n');

	foundit=false;
	lowerbound=hsshortcuts[static_cast<uint64>(srvaI->vhash & HS_MAXVHASHMASK)].b;

	// "HS_empty_vector_hashstat_t.end()" is the "nullptr" replacement
	if(hashstats.end() != lowerbound){
	  if(basesperhash>12){
	    // with more than 12 bases in a hash, the array is subdivided
	    hstmp.vhash=srvaI->vhash;
	    hssearchI=lower_bound(lowerbound,
				  hsshortcuts[static_cast<uint64>(srvaI->vhash & HS_MAXVHASHMASK)].e,
				  hstmp,
				  HashStatistics::sortHashStatComparatorLexicographicallyUp);
	    if(hssearchI != hashstats.end()
	       && hssearchI->vhash == srvaI->vhash) foundit=true;
	  }else{
	    hssearchI=lowerbound;
	    foundit=true;
	  }
	}else{
	  CEBUG("---------- NO LB HIT??? -------\n");
	}

	if(foundit) {
	  CEBUG("VHRAP: " << *srvaI << '\n');
	  CEBUG("HashStat: " << *hssearchI << '\n');
	  CEBUG("srvaI->hashpos: " << srvaI->hashpos << '\n');

	  bfpos1=actread.calcClippedPos2RawPos(srvaI->hashpos-(basesperhash-1));
	  bfpos2=bfpos1+basesperhash-1;

	  CEBUG("b bfpos1: " << bfpos1 << '\t' << bposhashstats[bfpos1] << endl);
	  CEBUG("b bfpos2: " << bfpos2 << '\t' << bposhashstats[bfpos2] << endl);

	  bposhashstats[bfpos1].fwd.setValid();
	  bposhashstats[bfpos2].rev.setValid();

	  if(hssearchI->hsc.hasfwdrevthresholdok) {
	    //bhs|=Read::BFLAGS_CONFIRMED_FWDREV;
	    CEBUG("Set ConfFWDREV\n");
	    bposhashstats[bfpos1].fwd.setConfirmedFwdRev();
	    bposhashstats[bfpos2].rev.setConfirmedFwdRev();
	  }
	  if(hssearchI->hsc.iskmerforkf) {
	    haskmerfork=true;
	    bposhashstats[bfpos1].fwd.setKMerFork();
	    if(!truekmerforks) bposhashstats[bfpos2].rev.setKMerFork();
	  }
	  if(hssearchI->hsc.getLowPos()<=4){
	    //bhs|=Read::BFLAGS_SEENATLOWPOS;
	    CEBUG("Set SeenAtLowPos\n");
	    bposhashstats[bfpos1].fwd.setSeenAtLowPos();
	    bposhashstats[bfpos2].rev.setSeenAtLowPos();
	  }
	  if(hssearchI->hsc.seqtype==MULTISEQTYPE){
	    //bhs|=Read::BFLAGS_CONFIRMED_MULTIPLESEQTYPE;
	    CEBUG("Set ConfMultSeqType\n");
	    bposhashstats[bfpos1].fwd.setConfirmedMultipleSeqType();
	    bposhashstats[bfpos2].rev.setConfirmedMultipleSeqType();
	  }
	  uint8 frequency=2;
	  {
	    auto hgc=hssearchI->hsc.getCount();
	    if(hgc < minnormalhashcov) {
	      frequency=2;
	    }else if(hgc >= minnormalhashcov
		     && hgc <= maxnormalhashcov) {
	      frequency=3;
	      //}else if(hssearchI->count > minnormalhashcov*20){
	    }else if(hgc > crazyrepthashcov){
	      frequency=7;
	    }else if(hgc > heavyrepthashcov){
	      frequency=6;
	    }else if(hgc >= repeathashcov){
	      frequency=5;
	    }else{
	      frequency=4;
	    }
	    if(hgc == 1){
	      frequency=0;
	    }else if(hgc > 0
		     && hgc <= HS_rarekmercount ){
	      // maybe additional checks ... ?
	      frequency=1;
	    }
	    CEBUG("Set frequency: " << static_cast<uint16>(frequency) << endl);

	    if(masknastyrepeats && maskhashcov>0 && hgc >= maskhashcov){
	      CEBUG("mcmask " << bfpos1 << "\t" << bfpos1+basesperhash << '\n');
	      for(uint32 j=0; j<basesperhash; j++){
		mcmask[bfpos1+j]=1;
	      }
	    }
	  }

	  CEBUG("a1 bfpos1: " << bfpos1 << '\t' << bposhashstats[bfpos1] << endl);
	  CEBUG("a1 bfpos2: " << bfpos2 << '\t' << bposhashstats[bfpos2] << endl);

	  bposhashstats[bfpos1].fwd.setFrequency(frequency);
	  bposhashstats[bfpos2].rev.setFrequency(frequency);

	  CEBUG("a2 bfpos1: " << bfpos1 << '\t' << bposhashstats[bfpos1] << endl);
	  CEBUG("a2 bfpos2: " << bfpos2 << '\t' << bposhashstats[bfpos2] << endl);

	  hasfrequency[frequency]=1;

	  //cout.flush();
	  //actread.setBaseFlagsInClippedSequence(bhs,
	  //				      srvaI->hashpos-(basesperhash-1),
	  //				      basesperhash);
	  //actread.setHasBaseFlags(true);
	}
      }

      if(hasfrequency[3]){
	actread.setHasFreqAvg(true);
      }
      if(hasfrequency[5] || hasfrequency[6] || hasfrequency[7]){
	actread.setHasFreqRept(true);
      }

      //Read::setCoutType(Read::AS_TEXT);
      //CEBUG("### After ...\n" << actread << endl);


      // BaCh 07.04.2009 Bad Idea!!!
      // BaCh 12.07.2009 Why? Forgot ... :-(
      //// the fwd/rev of a read now looks like this (e.g.)
      //// (for better viewing dot == 0)
      ////
      //// f   ..........2222222233333....355555....................
      //// r   ................2222222....33333355555...............
      ////
      //// in dubio pro reo and to allow for potential matches,
      //// do this:
      ////
      //// f   ..........2222222233333....355555->..................
      //// r   ..............<-2222222....33333355555...............
      ////
      //// so that this
      ////
      //// f   ..........2222222233333....35555555555...............
      //// r   ..........2222222222222....33333355555...............
      ////
      //// is generated
      ////
      ////
      //
      //{
      //  uint32 bfposi=0;
      //  for(; bfposi<bposhashstats.size() && bposhashstats[bfposi].fwd.getFrequency()==0; bfposi++) {};
      //  uint32 bfpose=bfposi;
      //  for(; bfpose<bposhashstats.size() && bposhashstats[bfpose].rev.getFrequency()==0; bfpose++) {};
      //  if(bfposi<bposhashstats.size() && bfpose<bposhashstats.size()){
      //	for(uint32 i=bfposi; i<bfpose; i++){
      //	  bposhashstats[i].fwd=bposhashstats[bfpose].rev;
      //	}
      //  }
      //
      //  bfposi=bposhashstats.size()-1;
      //  for(; bfposi>0 && bposhashstats[bfposi].rev.getFrequency()==0; bfposi--) {};
      //  bfpose=bfposi;
      //  for(; bfpose>0 && bposhashstats[bfpose].fwd.getFrequency()==0; bfpose--) {};
      //  if(bfposi>0){
      //	for(uint32 i=bfposi; i>bfpose; i--){
      //	  bposhashstats[i].fwd=bposhashstats[bfpose].rev;
      //	}
      //  }
      //}


      // go through multicopy array and set MNRr tags for
      //  consecutive positions in read tagged as multicopy
      if(masknastyrepeats){
	bool inrun=false;
	uint32 runstart=0;
	uint32 pos=0;
	for(; pos<mcmask.size(); pos++){
	  CEBUG("pos: " << pos << '\t' << static_cast<uint16>(mcmask[pos]) << '\t' << inrun << '\n');
	  if(mcmask[pos]){
	    if(!inrun){
	      runstart=pos;
	      inrun=true;
	    }
	  }else{
	    if(inrun){
	      CEBUG("reprun " << actread.getName() << '\t' << runstart << '\t' << pos-1 << endl);
	      tmpmt.from=runstart;
	      tmpmt.to=pos-1;
	      actread.addTagO(tmpmt);
	      inrun=false;
	    }
	  }
	}
	if(inrun){
	  CEBUG("reprun " << actread.getName() << '\t' << runstart << '\t' << pos-1 << endl);
	  tmpmt.from=runstart;
	  tmpmt.to=pos-1;
	  actread.addTagO(tmpmt);
	}
      }
    }

    // this just sets kmer forks in reverse direction
    if(truekmerforks){
      CEBUG("REVkmerforks HWI-ST143:485:D0MR8ACXX:6:2204:4364:1904\n");
      {
	int32 bfpos=actread.calcClippedComplPos2RawPos(0);
	int32 bfposinc=-1;

	srvaI=singlereadvhraparray.begin();

	hashesmade=Skim<TVHASH_T>::transformSeqToVariableHash(
	  actreadid,
	  actread,
	  actread.getClippedComplementSeqAsChar(),
	  slen,
	  basesperhash,
	  srvaI,
	  false,
	  1,
	  tagmaskvector,
	  bposhashstats,
	  bfpos,
	  bfposinc
	  );
      }
      singlereadvhraparray.resize(hashesmade);

      CEBUG("hashesmade: " << hashesmade << endl);
      CEBUG("maskhashcov: " << maskhashcov << endl);

      typename std::vector<hashstat_t>::const_iterator lowerbound;

      typename std::vector<hashstat_t>::const_iterator hssearchI;
      srvaI=singlereadvhraparray.begin();

      hashstat_t hstmp;
      bool foundit=false;
      for(; srvaI != singlereadvhraparray.end(); srvaI++){
	CEBUG(*srvaI << '\n');

	foundit=false;
	lowerbound=hsshortcuts[static_cast<uint64>(srvaI->vhash & HS_MAXVHASHMASK)].b;

	// "HS_empty_vector_hashstat_t.end()" is the "nullptr" replacement
	if(hashstats.end() != lowerbound){
	  if(basesperhash>12){
	    // with more than 12 bases in a hash, the array is subdivided
	    hstmp.vhash=srvaI->vhash;
	    hssearchI=lower_bound(lowerbound,
				  hsshortcuts[static_cast<uint64>(srvaI->vhash & HS_MAXVHASHMASK)].e,
				  hstmp,
				  HashStatistics::sortHashStatComparatorLexicographicallyUp);
	    if(hssearchI != hashstats.end()
	       && hssearchI->vhash == srvaI->vhash) foundit=true;
	  }else{
	    hssearchI=lowerbound;
	    foundit=true;
	  }
	}else{
	  CEBUG("---------- NO LB HIT??? -------\n");
	}

	if(foundit) {
	  CEBUG("VHRAP: " << *srvaI << '\n');
	  CEBUG("HashStat: " << *hssearchI << '\n');
	  CEBUG("srvaI->hashpos: " << srvaI->hashpos << '\n');

	  if(hssearchI->hsc.iskmerforkr) {
	    haskmerfork=true;

	    int32 bfpos1=actread.calcClippedComplPos2RawPos(srvaI->hashpos-(basesperhash-1));
	    CEBUG("-------------------- b bfpos1: " << bfpos1 << '\t' << srvaI->hashpos << '\t' << bposhashstats[bfpos1] << endl);

	    bposhashstats[srvaI->hashpos].rev.setKMerFork();
	  }
	}
      }

    }

    actread.setHasKMerFork(haskmerfork);

    Read::setCoutType(Read::AS_TEXT);
    CEBUG("### After ...\n" << actread << endl);
  }

  FUNCEND();
  return;
}
//#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define prefetchrl(p)     __builtin_prefetch((p), 0, 3)

/*
template<typename TVHASH_T>
uint32 HashStatistics<TVHASH_T>::checkBaitHit(Read & actread, std::vector<vhrap_t> & baiting_singlereadvhraparray, std::vector<uint8> & baiting_tagmaskvector, bool changeseqcase)
{
  //, const uint32 basesperhash, std::vector<hashstat_t> & hashstats, std::vector<std::vector<hashstat_t>::const_iterator > & hsshortcuts_begin, std::vector<std::vector<hashstat_t>::const_iterator > & hsshortcuts_end)
  FUNCSTART("uint32 HashStatistics<TVHASH_T>::checkBaitHit(Read & actread, std::vector<vhrap_t> & baiting_singlereadvhraparray, std::vector<uint8> & baiting_tagmaskvector, bool changeseqcase)");

  if(changeseqcase) actread.upDownCase(127); // TODO: replace with a real tolower() asap

  if(HS_hsv_hashstats.empty()) return 0;
  if(!actread.hasValidData()) return 0;
  uint32 slen=actread.getLenClippedSeq();
  if(slen<HS_hs_basesperhash) return 0;

  if(unlikely(HS_hsv_hsshortcuts.empty())) priv_makeHashStatArrayShortcuts();

  // don't really need to clear out these re-used vectors
  //   - tagmask just needs to be empty
  //   - singlereadvhraparray needs to be big enough to be
  //     written into by transformSeqToVariableHash(), will be
  //     resized later on num hashes made
  baiting_tagmaskvector.clear();
  if(baiting_singlereadvhraparray.size() < slen){
    baiting_singlereadvhraparray.resize(slen);
  }

  auto srvaI=baiting_singlereadvhraparray.begin();

  std::vector<Read::bposhashstat_t> & bposhashstats=const_cast<std::vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
  uint32 hashesmade;

  {
    int32 bfpos=0;
    int32 bfposinc=1;

    uint32 actreadid=0;

    hashesmade=Skim<TVHASH_T>::transformSeqToVariableHash(
      actreadid,
      actread,
      actread.getClippedSeqAsChar(),
      slen,
      HS_hs_basesperhash,
      srvaI,
      false,
      1,
      baiting_tagmaskvector,
      bposhashstats,
      bfpos,
      bfposinc
      );
  }
  baiting_singlereadvhraparray.resize(hashesmade);

  CEBUG("hashesmade: " << hashesmade << endl);

  typename std::vector<hashstat_t>::const_iterator lowerbound;

  typename std::vector<hashstat_t>::const_iterator hssearchI;
  srvaI=baiting_singlereadvhraparray.begin();

  hashstat_t hstmp;
  bool foundit;
  uint32 numhits=0;
  for(; srvaI != baiting_singlereadvhraparray.end(); srvaI++){
    CEBUG(*srvaI << '\n');

    lowerbound=HS_hsv_hsshortcuts[static_cast<uint64>(srvaI->vhash & HS_MAXVHASHMASK)].b;

    foundit=false;

    if(HS_hsv_hashstats.end() != lowerbound){
      if(HS_hs_basesperhash>12){
	// with more than 12 bases in a hash, the array is subdivided
	hstmp.vhash=srvaI->vhash;
	hssearchI=lower_bound(lowerbound,
			      HS_hsv_hsshortcuts[static_cast<uint64>(srvaI->vhash & HS_MAXVHASHMASK)].e,
			      hstmp,
			      sortHashStatComparatorLexicographicallyUp);
	if(hssearchI != HS_hsv_hashstats.end()
	   && hssearchI->vhash == srvaI->vhash) foundit=true;
      }else{
	hssearchI=lowerbound;
	foundit=true;
      }
    }else{
      CEBUG("---------- NO LB HIT??? -------\n");
    }

    if(foundit) {
      ++numhits;

      if(changeseqcase){
	// TODO: quite inefficient, change ASAP
	for(uint32 pi=0;pi<HS_hs_basesperhash; ++pi){
	  auto cpos=pi+srvaI->hashpos-HS_hs_basesperhash+1;
	  actread.changeBaseInClippedSequence(toupper(actread.getBaseInClippedSequence(cpos)),255,cpos);
	}
      }
    }
  }

  //cout << "\nskim Needs redo!\n";
  //exit(0);

  FUNCEND();
  return numhits;
}
*/

template<typename TVHASH_T>
uint32 HashStatistics<TVHASH_T>::checkBaitHit(Read & actread, bool changeseqcase, char mask)
{
  FUNCSTART("baithit2");

  if(changeseqcase) actread.upDownCase(127); // TODO: replace with a real tolower() asap

  if(HS_hsv_hashstats.empty()) return 0;
  if(!actread.hasValidData()) return 0;
  uint64 slen=actread.getLenClippedSeq();
  if(slen<HS_hs_basesperhash) return 0;

  if(unlikely(HS_hsv_hsshortcuts.empty())) priv_makeHashStatArrayShortcuts();

  typename std::vector<hashstat_t>::const_iterator lowerbound;
  typename std::vector<hashstat_t>::const_iterator hssearchI;

  hashstat_t hstmp;
  bool foundit;
  uint32 numhits=0;

  const uint8 * seq = reinterpret_cast<const uint8 *>(actread.getClippedSeqAsChar());
  const char *  namestr=actread.getName().c_str();
  const uint32 basesperhash=HS_hs_basesperhash;

  SEQTOHASH_LOOPSTART(TVHASH_T);
  {
    lowerbound=HS_hsv_hsshortcuts[static_cast<uint64>(acthash & HS_MAXVHASHMASK)].b;

    foundit=false;

    if(HS_hsv_hashstats.end() != lowerbound){
      if(HS_hs_basesperhash>12){
	// with more than 12 bases in a hash, the array is subdivided
	hstmp.vhash=acthash;
	hssearchI=lower_bound(lowerbound,
			      HS_hsv_hsshortcuts[static_cast<uint64>(acthash & HS_MAXVHASHMASK)].e,
			      hstmp,
			      sortHashStatComparatorLexicographicallyUp);
	if(hssearchI != HS_hsv_hashstats.end()
	   && hssearchI->vhash == acthash) foundit=true;
      }else{
	hssearchI=lowerbound;
	foundit=true;
      }
    }else{
      CEBUG("---------- NO LB HIT??? -------\n");
    }

    if(foundit) {
      ++numhits;

      if(changeseqcase || mask){
	// TODO: quite inefficient, change ASAP
	for(uint32 pi=0;pi<HS_hs_basesperhash; ++pi){
	  auto cpos=pi+static_cast<uint32>(seqi)-HS_hs_basesperhash+1;
	  auto actbase=actread.getBaseInClippedSequence(cpos);
	  if(mask){
	    actbase=mask;
	  } else {
	    actbase=toupper(actbase);
	  }
	  actread.changeBaseInClippedSequence(actbase,255,cpos);
	}
      }
    }
  }SEQTOHASH_LOOPEND;

  //cout << "\nskim Needs redo!\n";
  //exit(0);

  FUNCEND();
  return numhits;
}





template<typename TVHASH_T>
const typename HashStatistics<TVHASH_T>::hashstat_t * HashStatistics<TVHASH_T>::findVHash(const hashstat_t & searchval)
{
  FUNCSTART("const typename HashStatistics<TVHASH_T>::hashstat_t * HashStatistics<TVHASH_T>::findVHash(const hashstat_t & searchval)");

  const hashstat_t * ret=nullptr;

  if(unlikely(HS_hsv_hsshortcuts.empty())){
    priv_makeHashStatArrayShortcuts();
    BUGIFTHROW(unlikely(HS_hsv_hsshortcuts.empty()),"no shortcuts made ... empty hashstats?");
  }

  auto hsindex=static_cast<uint64>(searchval.vhash & HS_MAXVHASHMASK);
  if(likely(!HS_hsv_hashstats.empty())
     && HS_hsv_hashstats.end() != HS_hsv_hsshortcuts[hsindex].b){
    auto hsI=HS_hsv_hsshortcuts[hsindex].b;
    // TODO: test with large & diverse data set effect of prefetch
    prefetchrl(&(*hsI));
    if(HS_hsv_hsshortcuts[hsindex].e-hsI > 1){
      // with more than 12 bases in a hash, the array is subdivided
	// For small arrays, linear search will still be faster (caches)
	//   from miranar: sweet spot is around 32
      if(HS_hsv_hsshortcuts[hsindex].e-hsI > 32){
	hsI=lower_bound(hsI,
			HS_hsv_hsshortcuts[hsindex].e, // upperbound
			searchval,
			sortHashStatComparatorLexicographicallyUp);
      }else{
	while(hsI!=HS_hsv_hsshortcuts[hsindex].e && hsI->vhash!=searchval.vhash){
	  ++hsI;
	}
      }
    }
    if(hsI != HS_hsv_hashstats.end()
       && hsI->vhash == searchval.vhash) ret=&(*hsI);
  }

  return ret;
}




/*************************************************************************
 *
 * test
 *
 * implicit return:
 *  - dn_vhashindexes with indexes to all valid vhashes in sequence
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {if(docebug) {cout << bla; cout.flush();}}
template<typename TVHASH_T>
bool HashStatistics<TVHASH_T>::priv_dn_TestSingleSeq(Read & actread, std::vector<uint8> & dn_allow, std::vector<size_t> & dn_vhashindexes)
{
  FUNCSTART("bool HashStatistics<TVHASH_T>::priv_dn_TestSingleSeq(Read & actread, std::vector<uint8> & dn_allow, std::vector<size_t> & dn_vhashindexes)");

  //bool docebug=false;
  BUGIFTHROW(HS_hsv_hsshortcuts.empty(),"no shortcuts made, not ready for searching?");

  const uint8 * seq = reinterpret_cast<const uint8 *>(actread.getClippedSeqAsChar());
  uint64 slen=actread.getLenClippedSeq();

  if(slen<HS_hs_basesperhash) return false;

  const char *  namestr=actread.getName().c_str();

  dn_vhashindexes.clear();
  dn_allow.clear();
  dn_allow.resize(slen,1);

  // TODO: option to have it only on MNRr or given HAFx stretches
  auto bhsI=actread.getBPosHashStats().begin();
  bhsI+=actread.getLeftClipoff();
  for(auto ri=0; ri<slen; ++ri){
    if(bhsI->fwd.getFrequency()<2 || !bhsI->fwd.hasConfirmedFwdRev()){
      dn_allow[ri]=0;
    }
  }


  hashstat_t searchval;
  bool takeread=false;

  auto basesperhash=HS_hs_basesperhash;

  SEQTOHASH_LOOPSTART(TVHASH_T){

    auto hssi=static_cast<uint64>(acthash & HS_MAXVHASHMASK);
    if(HS_hsv_hashstats.end() != HS_hsv_hsshortcuts[hssi].b){
      auto hsI=HS_hsv_hsshortcuts[hssi].b;
      // TODO: test with large & diverse data set effect of prefetch
      prefetchrl(&(*hsI));
      if(HS_hsv_hsshortcuts[hssi].e-hsI > 1){
	// with more than 12 bases in a hash, the array is subdivided
	// For small arrays, linear search will still be faster (caches)
	//   from miranar: sweet spot is around 32
	if(HS_hsv_hsshortcuts[hssi].e-hsI > 32){
	  searchval.vhash=acthash;
	  hsI=lower_bound(hsI,
			  HS_hsv_hsshortcuts[hssi].e, // upperbound
			  searchval,
			  sortHashStatComparatorLexicographicallyUp);
	}else{
	  while(hsI!=HS_hsv_hsshortcuts[hssi].e && hsI->vhash!=acthash){
	    ++hsI;
	  }
	}
      }

      if(hsI != HS_hsv_hashstats.end()
	 && hsI->vhash == acthash) {
	// hsI on valid valid hash
	size_t hsindex=hsI-HS_hsv_hashstats.begin();
	CEBUG("hashfound " << seqi << "\t" << hsindex << endl);
	dn_vhashindexes.push_back(hsindex);
	if(dn_allow[seqi] && HS_diginorm_count[hsindex]<10){
	  takeread=true;
	  // NOOOOO: do not break here! We need to test every hash to build the vhashindex with all of them!
	  // break;
	}
      }else{
	CEBUG("no hash? " << seqi << endl);
      }
    }

  }SEQTOHASH_LOOPEND;

  return takeread;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
bool HashStatistics<TVHASH_T>::digiNormTestRead(Read & actread, bool forcetake)
{
  FUNCSTART("bool HashStatistics<TVHASH_T>::digiNorm(Read & actread)");

  if(unlikely(HS_diginorm_count.empty())){
    HS_diginorm_count.resize(HS_hsv_hashstats.size(),0);
  }

  if(!actread.hasTag(Read::REA_defaulttag_MNRr.identifier)) return true;

  bool takeread=priv_dn_TestSingleSeq(actread,HS_diginorm_allow_s1,HS_diginorm_vhashindexes_s1);

  if(forcetake) {
    CEBUG("Forced take\n");
    takeread=true;
  }

  if(takeread){
    CEBUG("dntr take " << actread.getName() << ": " << HS_diginorm_vhashindexes_s1.size() << endl);
    for(auto hsi : HS_diginorm_vhashindexes_s1){
      ++HS_diginorm_count[hsi];
    }
  }

  return takeread;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * test
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<typename TVHASH_T>
uint32 HashStatistics<TVHASH_T>::estimDigiNormCov(Read & actread)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::estimDigiNormCov(Read & actread)");

  BUGIFTHROW(HS_hsv_hsshortcuts.empty(),"no shortcuts made, not ready for searching?");

  const uint8 * seq = reinterpret_cast<const uint8 *>(actread.getClippedSeqAsChar());
  uint64 slen=actread.getLenClippedSeq();

  if(slen<HS_hs_basesperhash) return 1;

  const char *  namestr=actread.getName().c_str();

  auto basesperhash=HS_hs_basesperhash;


  double dncmin=10000000.0;
  bool hasnewmin=false;
  //double dncmax=0.0;
  double totaladd=0.0;
  uint32 numtotal=0;


  hashstat_t searchval;

  SEQTOHASH_LOOPSTART(TVHASH_T){

    auto hssi=static_cast<uint64>(acthash & HS_MAXVHASHMASK);
    if(HS_hsv_hashstats.end() != HS_hsv_hsshortcuts[hssi].b){
      auto hsI=HS_hsv_hsshortcuts[hssi].b;
      // TODO: test with large & diverse data set effect of prefetch
      prefetchrl(&(*hsI));
      if(HS_hsv_hsshortcuts[hssi].e-hsI > 1){
	// with more than 12 bases in a hash, the array is subdivided
	// For small arrays, linear search will still be faster (caches)
	//   from miranar: sweet spot is around 32
	if(HS_hsv_hsshortcuts[hssi].e-hsI > 32){
	  searchval.vhash=acthash;
	  hsI=lower_bound(hsI,
			  HS_hsv_hsshortcuts[hssi].e, // upperbound
			  searchval,
			  sortHashStatComparatorLexicographicallyUp);
	}else{
	  while(hsI!=HS_hsv_hsshortcuts[hssi].e && hsI->vhash!=acthash){
	    ++hsI;
	  }
	}
      }

      if(hsI != HS_hsv_hashstats.end()
	 && hsI->vhash == acthash) {
	// hsI on valid valid hash

	auto hsindex=hsI-HS_hsv_hashstats.begin();

	// *sigh* this is very well possible
	// BUGIFTHROW(HS_diginorm_count[hsindex]==0,actread.getName() << " HS_diginorm_count[hsindex]==0 ???");
	// e.g. in first pass a DGNR tag gets set like so: ...TTTTT
	// and somewhere the read gets cut back like so: ...
	// then the tag still may exist.
	// therefore: if count is 0, then do as if no hash existed

	if(HS_diginorm_count[hsindex]!=0){
	  double actdnc=static_cast<double>(hsI->hsc.getCount())/HS_diginorm_count[hsindex];
	  ++numtotal;
	  totaladd+=actdnc;
	  //if(actdnc>dncmax) dncmax=actdnc;
	  if(actdnc<dncmin) {
	    dncmin=actdnc;
	    hasnewmin=true;
	  }
	  CEBUG(actread.getName() << "\t" << hsI->hsc.getCount() << "\t" << HS_diginorm_count[hsindex] << "\t" << actdnc << endl);
	}
      }else{
	//BUGIFTHROW(true,"Can't be?");
	// well ... can be, as the haststats does contain only "valid" hashes where we are sure they're present multiple times etc.,
	//  we might encounter a hash here which is not in hashstats.
	// therefore, this is a singlet event

	// but singlet events were, by default, not accounted for in digiNormTestRead(), therefore
	//  we should not return "1" but simply do nothing and continue calculating
	//return 1;
      }
    }

  }SEQTOHASH_LOOPEND;

  CEBUG("dncstats " << actread.getName() << ": " << dncmin);
  //CEBUG("\t" << dncmax(;
  CEBUG("\t" << totaladd/numtotal);
  CEBUG(endl);

  //if(dncmin<=6) return dncmin;

  if(!hasnewmin) {
    CEBUG("WTH??? " << actread.getName() << " has no dncmin ???\n");
    return 1;
  }

  // additional /2 because we are treating forward and reverse hashes separately
  return static_cast<uint32>((totaladd/numtotal/2)+.5);
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::dump(std::ostream & ostr)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::dump(std::ostream & ostr)");

  std::string tmpstr;
  for(auto & hse : HS_hsv_hashstats){
    hash2string(hse.vhash,HS_hs_basesperhash,tmpstr);
    ostr << tmpstr
	 << '\t' << hse.hsc.fcount
	 << '\t' << hse.hsc.rcount
	 << '\t' << hse.hsc.fcount + hse.hsc.rcount
	 << '\n';
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::dumpHSDebug(std::ostream & ostr)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::dumpHSDebug(std::ostream & ostr)");

  std::string tmpstr;
  for(auto & hse : HS_hsv_hashstats){
    hash2string(hse.vhash,HS_hs_basesperhash,tmpstr);
    ostr << tmpstr
	 << '\t' << hse
	 << '\n';
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::dumpAsFASTA(std::ostream & ostr)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::dumpAsFASTA(std::ostream & ostr)");

  std::string tmpstr;
  uint64 counter=1;
  for(auto & hse : HS_hsv_hashstats){
    hash2string(hse.vhash,HS_hs_basesperhash,tmpstr);
    ostr << ">" << counter++
	 << ' ' << hse.hsc.fcount
	 << ' ' << hse.hsc.rcount
	 << ' ' << hse.hsc.fcount + hse.hsc.rcount
	 << '\n' << tmpstr
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
