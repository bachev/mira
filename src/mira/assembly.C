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

#include "util/machineinfo.H"

#include "mira/align.H"

// BOOST
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>


#include <cstdlib>
#include <cstdio>



#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

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




void Assembly::test()
{
  if(AS_readpool[2].getLSClipoff()>0){
    Read::setCoutType(Read::AS_TEXT);
    cout << AS_readpool[2];
    exit(10);
  }

//  uint32 s=5000000;
//  AS_readpool.reserve(s);
//  for(uint32 i=0;i<s;i++) AS_readpool.addNewEmptyRead();
//  AS_permanent_overlap_bans.resize(s);
//  AS_istroublemaker.resize(s,0);
//  loadAlignmentsFromFile(-1,"miratmp","", ".ads_pass");
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

Assembly::Assembly(Manifest & manifest, std::vector<MIRAParameters> & params, bool resumeassembly): AS_dataprocessing(&params)
{
  FUNCSTART("Assembly::Assembly(MIRAParameters * params)");

  AS_manifest=manifest;
  AS_miraparams=params;
  AS_resumeasembly=resumeassembly;

  init();

  setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);

  // For resuming assemblies
  if(resumeassembly){
    if(!dirExists(AS_miraparams[0].getDirectoryParams().dir_top)){
      MIRANOTIFY(Notify::FATAL,"Could not not find directory " << AS_miraparams[0].getDirectoryParams().dir_top << " while trying to resume the assembly. If you intended to start a new assembly, do not use '-r' (resume) when starting MIRA.");
    }
    if(!dirExists(AS_miraparams[0].getDirectoryParams().dir_checkpoint)){
      MIRANOTIFY(Notify::FATAL,"Could not not find directory " << AS_miraparams[0].getDirectoryParams().dir_checkpoint << " while trying to resume the assembly. Resuming is impossible, sorry.");
    }
    // do purge the results directory!
    if(ensureDirectory(AS_miraparams[0].getDirectoryParams().dir_results, true)){
      MIRANOTIFY(Notify::FATAL,"Could not delete and recreate results directory while trying to resume an assembly?");
    }
  }

  // purge the the remaining directories (if we're not resuming)
  ensureStandardDirectories(!resumeassembly);

  // after ensureStandardDirectories() as usually located in info directory
  AS_warnings.setOutputPath(AS_miraparams[0].getDirectoryParams().dir_info+"/"+AS_miraparams[0].getAssemblyParams().as_outfile_stats_warnings);

  setContigBuiltCallback();

  AS_assemblyinfo.setLargeContigSize(500);

  AS_systemmemory=MachineInfo::getMemTotal();

  //makeTmpDir();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::init()
{
  AS_steps.resize(ASNUMOFSTEPS);

  for(int8 i=ASNUMOFSTEPS-1; i>=0;i--){
    AS_steps[i]=0;
  }

  AS_tmptag_CRMr=Read::REA_defaulttag_CRMr;

  // initialise some other variables
  zeroVars();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::zeroVars()
{
  AS_num_reads_valid=0;
  AS_num_reads_too_small=0;
  AS_numADSFacts_fromalignments=0;
  AS_numADSFacts_fromshreds=0;
  AS_seqtypespresent.clear();
  AS_hasbackbones=false;

  AS_hashstat_avghashfreq=0;

  AS_doneskimchimera=false;
  AS_resumeisok=false;

  AS_shouldrun_nfs_check=true;

  AS_coveragetotal=0;

  AS_everythingwentfine=false;

  //TODO: rest
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Assembly::~Assembly()
{
  FUNCSTART("Assembly::~Assembly()");

  dumpMemInfo();

  AS_warnings.dumpWarnings();

  cout << "Dynamic s allocs: " << Dynamic::DYN_alloccounts << endl;
  cout << "Dynamic m allocs: " << Dynamic::DYN_alloccountm << endl;
  cout << "Align allocs: " << Align::AL_alloccount << endl;

  discard();

  if(AS_everythingwentfine && AS_miraparams[0].getAssemblyParams().as_output_removetmpdir){
    removeDirectory(AS_miraparams[0].getDirectoryParams().dir_tmp,false,false);
  }

  // TODO: scandir on result and remove all if no results?

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO not complete

void Assembly::discard()
{
  FUNCSTART("Assembly::discard()");

  AS_readpool.discard();
  AS_contigs.clear();
  AS_bbcontigs.clear();
  //AS_ok_for_assembly.clear();

  nukeSTLContainer(AS_adsfacts);
  nukeSTLContainer(AS_confirmed_edges);

  nukeSTLContainer(AS_used_ids);
  nukeSTLContainer(AS_multicopies);
  nukeSTLContainer(AS_hasmcoverlaps);
  nukeSTLContainer(AS_maxcoveragereached);
  nukeSTLContainer(AS_steps);
  nukeSTLContainer(AS_istroublemaker);
  nukeSTLContainer(AS_allrmbsok);
  nukeSTLContainer(AS_probablermbsnotok);
  nukeSTLContainer(AS_weakrmbsnotok);

  AS_permanent_overlap_bans.nuke();

  nukeSTLContainer(AS_readhitmiss);
  nukeSTLContainer(AS_readhmcovered);
  nukeSTLContainer(AS_count_rhm);
  nukeSTLContainer(AS_clipleft);
  nukeSTLContainer(AS_clipright);

  init();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO not complete

void Assembly::dmi_dumpALine(std::ostream & ostr, const char * desc, size_t numelem, size_t bytes_size, size_t freecapacity, size_t lostbyalign)
{

  ostr << std::setw(30) << desc
       << std::setw(10) << numelem;

  {
    std::ostringstream ostrstr;
    byteToHumanReadableSize(static_cast<double>(bytes_size), ostrstr);
    ostr << std::setw(12) << ostrstr.str();
  }
  {
    std::ostringstream ostrstr;
    byteToHumanReadableSize(static_cast<double>(freecapacity), ostrstr);
    ostr << std::setw(12) << ostrstr.str();
  }
  {
    std::ostringstream ostrstr;
    byteToHumanReadableSize(static_cast<double>(lostbyalign), ostrstr);
    ostr << std::setw(12) << ostrstr.str();
  }

  ostr << endl;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpMemInfo()
{
  FUNCSTART("Assembly::dumpMemInfo()");

  cout << "\n\n========================== Memory self assessment ==============================\n";

  size_t bytes_size = 0;
  size_t tmp_bytes_size=0;

  size_t numelem=0;
  size_t tmp_numelem=0;

  size_t freecapacity=0;
  size_t tmp_freecapacity=0;

  size_t lostbyalign=0;
  size_t tmp_lostbyalign=0;

  // we currently do not use these
  (void) numelem;
  (void) freecapacity;
  (void) lostbyalign;


  if(sizeof(size_t) == sizeof(int32)){
    cout << "Running in 32 bit mode.\n\n";
  }else{
    cout << "Running in 64 bit mode.\n\n";
  }

  dumpFile("/proc/meminfo",cout);
  cout << '\n';
  dumpFile("/proc/self/status",cout);

  cout << "\nInformation on current assembly object:\n\n";

  cout << "AS_readpool: " << AS_readpool.size() << " reads.\n";
  cout << "AS_contigs: " << AS_contigs.size() << " contigs.\n";
  cout << "AS_bbcontigs: " << AS_bbcontigs.size() << " contigs.\n";

  bytes_size+=AS_readpool.estimateMemoryUsage();
  cout << "Mem used for reads: " << bytes_size  << " (";
  byteToHumanReadableSize(static_cast<double>(bytes_size), cout);
  cout << ")\n\nMemory used in assembly structures:\n"
       << std::setw(52) << "Eff. Size" << std::setw(12) << "Free cap." << std::setw(12) << "LostByAlign" << endl;


  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_writtenskimhitsperid,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_writtenskimhitsperid: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skim_edges,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skim_edges: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_adsfacts,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_adsfacts: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_confirmed_edges,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_confirmed_edges: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_permanent_overlap_bans,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_permanent_overlap_bans: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readhitmiss,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readhitmiss: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readhmcovered,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readhmcovered: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_count_rhm,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_count_rhm: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_clipleft,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_clipleft: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_clipright,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_clipright: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_used_ids,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_used_ids: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_multicopies,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_multicopies: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_hasmcoverlaps,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_hasmcoverlaps: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_maxcoveragereached,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_maxcoveragereached: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_coverageperseqtype,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_coverageperseqtype: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_istroublemaker,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_istroublemaker: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_isdebris,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_isdebris: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_needalloverlaps,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_needalloverlaps: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readsforrepeatresolve,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readsforrepeatresolve: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_allrmbsok,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_allrmbsok: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_probablermbsnotok,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_probablermbsnotok: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_weakrmbsnotok,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_weakrmbsnotok: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readmaytakeskim,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readmaytakeskim: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimstaken,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimstaken: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_numskimoverlaps,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_numskimoverlaps: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_numleftextendskims,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_numleftextendskims: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_numrightextendskims,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_rightextendskims: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimleftextendratio,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimleftextendratio: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimrightextendratio,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimrightextendratio: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimmegahubs,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimmegahubs: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_usedtmpfiles,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_usedtmpfiles: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;


  cout << "Total: " << bytes_size << " (";
  byteToHumanReadableSize(static_cast<double>(bytes_size), cout);
  cout << ")";
  cout << "\n\n================================================================================\n";

  FUNCEND();
}

//void Assembly::dmiAddBytesWithCapacity(size_t & bcap, size_t & bsize)
//{
//}

/*************************************************************************
 *
 * go through list of previous files and delete old with same base name
 *  but different file name
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

uint32 Assembly::cleanupOldFile(const std::string & basename, const std::string & filename)
{
  FUNCSTART("uint32 Assembly::cleanupOldFile(const std::string & basename, const std::string & filename)");

  CEBUG("\nCOF: ###" << basename << "### and ###"<<filename<<"###\n");

  boost::system::error_code ec;

  uint32 numdeleted=0;
  auto ulfI=AS_usedtmpfiles.begin();
  while(ulfI != AS_usedtmpfiles.end()){
    //cout << "cOF: " << ulfI->basename << '\t' << ulfI->filename << endl;
    if(ulfI->basename == basename
       && ulfI->filename != filename) {
      ++numdeleted;

      // the following two checks seem pointless as normally the list would
      //  arrive at an end sometime
      // however, I've seen one case where this loop became endless,
      //  hence this foolguard
      if(numdeleted>100) {
	cerr << "\n\nOUCH! something strange ... tried more than 100 deletes of " << basename << " ... list size is " << AS_usedtmpfiles.size() << '\n';
      }
      if(numdeleted>1200) {
	cerr << "\n\nOUCH! something weird ... tried more than 1200 deletes of " << basename << " ... list size is " << AS_usedtmpfiles.size() << '\n';
	cerr << "We'll stop that here.\n";
	return numdeleted;
      }

      fileRemove(ulfI->filename,true);

      if(boost::filesystem::exists(ulfI->filename,ec)){
	cerr << "WARNING: Could not delete old file " + ulfI->filename
	     << "\nThis can have a number of different reasons, none of them"
	     << "\nwarranting an abort, but this is strange anyway.\n\n";
	++ulfI;
      }else{
	ulfI=AS_usedtmpfiles.erase(ulfI);
      }
    }else{
      ++ulfI;
    }
  }

  FUNCEND();
  return numdeleted;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::buildFileName(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, const std::string & suffix, const std::string & dirname, bool removeold)
{
  FUNCSTART("std::string Assembly::buildFileName(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, const std::string & suffix, const std::string & dirname, bool removeold)");

  std::ostringstream ostr;

  if(version>=0){
    ostr << AS_miraparams[0].getDirectoryParams().dir_tmp << "/";
  } else if(!dirname.empty()){
    ostr << dirname << "/";
  }

  ostr << prefix << basename << postfix ;

  if(version>=0){
    ostr << "." << version;
  }
  ostr << suffix;

  std::string filename=ostr.str();
  std::string newbasename(basename+suffix);

  if(removeold && AS_miraparams[0].getAssemblyParams().as_output_removerollovertmps) {
    cleanupOldFile(newbasename,filename);
  }

  bool mustadd=true;
  for(const auto & utfe : AS_usedtmpfiles){
    if(utfe.basename == newbasename
       && utfe.filename == filename) {
      mustadd=false;
      break;
    }
  }

  if(mustadd){
    usedtmpfiles_t ulf;
    AS_usedtmpfiles.push_back(ulf);
    AS_usedtmpfiles.back().basename=newbasename;
    AS_usedtmpfiles.back().filename=filename;
  }

  FUNCEND();
  return filename;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::ensureStandardDirectories(bool purge){
  FUNCSTART("void Assembly::ensureStandardDirectories(bool purge)");

  auto & dparams = AS_miraparams[0].getNonConstDirectoryParams();

  std::string existingtmpdir;
  // special purge logic to handle eventually existing symlinks to -DI:trt
  if(!dparams.dir_tmp_symlink.empty()
     && boost::filesystem::exists(dparams.dir_tmp_symlink)
     && boost::filesystem::is_symlink(dparams.dir_tmp_symlink)){
    existingtmpdir=boost::filesystem::read_symlink(dparams.dir_tmp_symlink).string();
  }else if(boost::filesystem::exists(dparams.dir_tmp)
	   && boost::filesystem::is_symlink(dparams.dir_tmp)){
    existingtmpdir=boost::filesystem::read_symlink(dparams.dir_tmp).string();
  }
  if(purge
     && !existingtmpdir.empty()){
    boost::filesystem::remove_all(existingtmpdir);
    existingtmpdir.clear();
  }

  // make sure the main directories exist or are created
  if(ensureDirectory(dparams.dir_top, purge, true, false)
     || ensureDirectory(dparams.dir_results, purge, true, false)
     || ensureDirectory(dparams.dir_info, purge, true, false)
     || ensureDirectory(dparams.dir_checkpoint, purge, true, false)){

    MIRANOTIFY(Notify::FATAL, "Could not make sure that a needed directory exists (see log above for more info), aborting MIRA.");
  }

  // make sure the tmp directory exists or is created
  // either as directory or as symlinked directory (-DI:trt=...)
  if(dparams.dir_tmp_symlink.empty()){
    if(ensureDirectory(dparams.dir_tmp, purge, true, false)){
      MIRANOTIFY(Notify::FATAL, "Could not make sure that the MIRA tmp directory exists, aborting.");
    }
  }else{
    if(!existingtmpdir.empty() && !boost::filesystem::is_directory(existingtmpdir)){
      // existing ist not a directory??? Should be very, very ... very rare
      // but if yes, let's purge it
      boost::filesystem::remove_all(existingtmpdir);
      existingtmpdir.clear();
    }
    if(existingtmpdir.empty()){
      existingtmpdir=dparams.dir_tmp+"_XXXXXX";
      auto * ptr=mkdtemp(const_cast<char *>(existingtmpdir.c_str()));
      if(ptr==nullptr){
	perror(static_cast<std::string>("Could not create directory for temporary MIRA data "+existingtmpdir).c_str());
      }
      // contrary to what "man mkdtemp" may lead to believe, mkdtemp() does not return NULL (failure)
      //  if given a path to a non-existing directory ("/bla/bla_XXXXXX")
      // need to check ourselves
      if(!boost::filesystem::exists(existingtmpdir)){
	MIRANOTIFY(Notify::FATAL, "Could not create MIRA tmp directory \"" << existingtmpdir << "\": is some part of the path not existing or access protected?");
      }
    }

    dparams.dir_tmp=existingtmpdir;

    if(!boost::filesystem::exists(dparams.dir_tmp_symlink)){
      boost::filesystem::create_symlink(dparams.dir_tmp,dparams.dir_tmp_symlink);
    }

    cout << "Symlink " << dparams.dir_tmp_symlink << " now pointing to " << dparams.dir_tmp << endl;;
  }

  checkForNFSMountOnTmpDir();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::checkForNFSMountOnTmpDir()
{
  FUNCSTART("void Assembly::checkForNFSMountOnTmpDir()");

  if(!AS_shouldrun_nfs_check) return;

  auto res=checkForNFSMountOnDirectory(AS_miraparams[0].getDirectoryParams().dir_tmp, true);

  cout << '\n';
  if(res==0){
    cout << "Tmp directory is not on a NFS mount, good.\n\n";
  }else if(res==1){
    cout << "\nMake sure " << AS_miraparams[0].getDirectoryParams().dir_tmp
	 << " is *NOT* on a NFS mount or else MIRA will run *very* slowly.\n";
  }else{
    cout << "\n\n\n\n\nWARNING WARNING WARNING!\n\n"
      "It looks like the directory MIRA uses for temporary files\n    " << AS_miraparams[0].getDirectoryParams().dir_tmp <<
      "\nis on a NFS (Network File System) mount. This will slow down MIRA *considerably*\n"
      "... by about a factor of 10!\n\n"
      "If you don't want that, you have three possibilities:\n\n"
      "1) RECOMMENDED! Use -DI:trt to redirect the tmp directory somewhere else on a\n"
      "   local disk or even SSD.\n"
      "2) ALSO POSSIBLE: put the whole project somewhere else on your file system.\n"
      "3) ABSOLUTELY NOT RECOMMENDED AT ALL: use \"-NW:cnfs=warn\" to tell MIRA not\n"
      "   to stop when it finds the tmp directory on NFS.\n\n"
      "If you do not know what NFS is and which directory to use in \"-DI:trt\", ask\n"
      "your local system administrator to guide you.\n\n";

    if(AS_miraparams[0].getNagAndWarnParams().nw_check_nfs==NWSTOP){
      MIRANOTIFY(Notify::FATAL,"Tmp directory is on a NFS mount ... but we don't want that.");
    }
  }

  AS_shouldrun_nfs_check=false;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::preassembleTasks(bool usereadextension, bool clipvectorleftovers)
{
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  {
    std::string tmpfname;
    tmpfname=buildFileName(0,"","",
			  as_fixparams.as_tmpf_clippings,
			   ".txt","",0);

    // doing it twice catches a few outliers missed the first time
    std::string logprefix="proposed cutback 1a: ";
    uint64 numclipped=performPECCHIMSRERKM(tmpfname,logprefix);
    if(numclipped>0){
      logprefix="proposed cutback 1b: ";
      performPECCHIMSRERKM(tmpfname,logprefix);
    }else{
      cout << "No bases clipped in first pec round, skipping second round.\n";
    }
    dumpSomeStatistics();
  }

  //performSnapshot(0);

  //performHashEditing();

#if TRACKMEMUSAGE
  cout << "\ndmi pre 00\n";
  dumpMemInfo();
#endif

  //if(AS_454dosimpleedit) editSimple454Overcalls(0);

  if(clipvectorleftovers
	 || (usereadextension && as_fixparams.as_readextension_firstpassnum == 0)){
    cout << "Pre-assembly alignment search for read extension and / or vector clipping:\n";

    findPossibleOverlaps(0, "", "_preassembly");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 61\n";
    dumpMemInfo();
#endif

    // do not use the 100% trans rule for read extension and clipping!
    // needed vectors would not get filled
    makeAlignments(Assembly::ma_takeall, false, false, 0, "", "_preassembly1");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 62a\n";
    dumpMemInfo();
#endif

    priv_loadAlignmentsFromFile(0, "", "_preassembly1");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 62b\n";
    dumpMemInfo();
#endif

    if(usereadextension) {
      cout << "Pre-assembly read extension:\n";
      extendADS(0, "", "_preassembly1");
#if TRACKMEMUSAGE
      cout << "\ndmi pre 62c\n";
      dumpMemInfo();
#endif
    }
    if(clipvectorleftovers) {
      cout << "Pre-assembly vector clipping\n";
      performSeqVectorClippings();
#if TRACKMEMUSAGE
      cout << "\ndmi pre 62d\n";
      dumpMemInfo();
#endif
    }

    dumpSomeStatistics();

    // we need to throw away all permbans that might have appeared in
    //  this pre-assembly sequence massage

    AS_permanent_overlap_bans.nuke();
    AS_permanent_overlap_bans.resize(AS_readpool.size());

#if TRACKMEMUSAGE
    cout << "\ndmi pre 63\n";
    dumpMemInfo();
#endif

    priv_performHashAnalysis("",false,false, 0, "", "_preassembly2");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 64\n";
    dumpMemInfo();
#endif

    {
      std::string tmpfname;
      tmpfname=buildFileName(0,"","",
			    as_fixparams.as_tmpf_clippings,
			     ".txt","",false);
      std::string logprefix="proposed cutback preassembly: ";

      performPECCHIMSRERKM(tmpfname,logprefix);
      dumpSomeStatistics();
#if TRACKMEMUSAGE
      cout << "\ndmi pre 64b\n";
      dumpMemInfo();
#endif
    }

    //nukeSTLContainer(AS_adsfacts);
    //nukeSTLContainer(AS_confirmed_edges);
    AS_adsfacts.clear();
    AS_confirmed_edges.clear();

#if TRACKMEMUSAGE
    cout << "\ndmi pre 65a\n";
    dumpMemInfo();
#endif
//    findPossibleOverlaps(0, "", "_preassembly2");
//#if TRACKMEMUSAGE
//    cout << "\ndmi pre 65b\n";
//    dumpMemInfo();
//#endif
  }
}


/*************************************************************************
 *
 * also determines the pass in which diginorm is applied (if needed)
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::priv_setupAutoPasses()
{
  auto & ncaparams=AS_miraparams[0].getNonConstAssemblyParams();

  if(!ncaparams.as_bphseries.empty()){
    cout << "user defined kmer series.\n";
    ncaparams.as_numpasses=ncaparams.as_bphseries.size();
  }else{
    // no bph series setup by user, we need to calc that
    AS_fixed_rls_bytype=AS_current_rls_bytype;
    AS_fixed_rls_byrg=AS_current_rls_byrg;

    auto & sparams=AS_miraparams[0].getSkimParams();

    uint64 maxavg=0;
    for(uint32 rgi=1; rgi<AS_fixed_rls_byrg.size();++rgi){
      auto rgid=ReadGroupLib::getReadGroupID(rgi);
      if(rgid.isBackbone() || rgid.isRail()) continue;
      maxavg=std::max(maxavg,AS_fixed_rls_byrg[rgi].getAvgLenUsed());
    }

    // two modes: either user set numpasses, then calculate bph series by
    //  initial bph and increase per pass
    // or user wanted all auto (numpasses<0), then use the info from the different average
    //  readgroup read lengths to decide
    if(ncaparams.as_numpasses>0){
      cout << "user defined number of passes.\n";

      auto stepping=sparams.sk_bph_increasestep;
      CEBUG("initial stepping " << stepping << endl);
      if(stepping==0){
	// no user supplied stepping. calc it equidistant in the range of [kmer,maxbph2use]
	// (yes, rounding errors ... I don't care)
	if(ncaparams.as_numpasses>1){
	  auto maxbphtouse=maxavg*100/60;
	  if(maxbphtouse<sparams.sk_basesperhash) maxbphtouse=sparams.sk_basesperhash;
	  stepping=(maxbphtouse-sparams.sk_basesperhash)/(ncaparams.as_numpasses-1);
	}else{
	  stepping=1;
	}
      }
      CEBUG("used stepping " << stepping << endl);
      for(uint8 ap=0; ap<ncaparams.as_numpasses; ++ap){
	auto newbph=sparams.sk_basesperhash+(stepping*ap);
	if(sparams.sk_bph_max && newbph>sparams.sk_bph_max) newbph=sparams.sk_bph_max;
	if(newbph>256) newbph=256;
	ncaparams.as_bphseries.push_back(newbph);
      }
    }else{
      cout << "all-auto determination of passes and kmer series.\n";

      CEBUG("ncaparams.as_assemblyjob_accurate " << ncaparams.as_assemblyjob_accurate << endl);
      CEBUG("maxavg " << maxavg << endl);

      // TODO: get that cleaned up
      uint32 lastkmer=(maxavg*7)/10;
      if(lastkmer>255) lastkmer=255;
      if(lastkmer==0) lastkmer=31;   // should never happen
      bool addlast=true;
      if(ncaparams.as_assemblyjob_accurate){
	if(maxavg>383) {
	  ncaparams.as_bphseries={17,31,63,127};
	}else if(maxavg>290) {
	  ncaparams.as_bphseries={17,31,63,127};
	}else if(maxavg>240) {
	  ncaparams.as_bphseries={17,31,63,127};
	}else if(maxavg>190) {
	  ncaparams.as_bphseries={17,31,63,95};
	}else if(maxavg>140) {
	  ncaparams.as_bphseries={17,31,53,75};
	}else if(maxavg>90) {
	  ncaparams.as_bphseries={17,31,53};
	}else if(maxavg>70) {
	  ncaparams.as_bphseries={17,27,37};
	}else if(maxavg>45) {
	  ncaparams.as_bphseries={17,21,25,31};
	  addlast=false;
	}else{
	  ncaparams.as_bphseries={17,19,21,23};
	  addlast=false;
	}
      }else{
	if(maxavg>383) {
	  ncaparams.as_bphseries={17,31};
	}else if(maxavg>290) {
	  ncaparams.as_bphseries={17,31};
	}else if(maxavg>240) {
	  ncaparams.as_bphseries={17,31};
	}else if(maxavg>190) {
	  ncaparams.as_bphseries={17,31};
	}else if(maxavg>140) {
	  ncaparams.as_bphseries={17,31};
	}else if(maxavg>90) {
	  ncaparams.as_bphseries={17,31};
	}else if(maxavg>70) {
	  ncaparams.as_bphseries={17,27};
	}else if(maxavg>45) {
	  ncaparams.as_bphseries={17,21,31};
	  addlast=false;
	}else{
	  ncaparams.as_bphseries={17,19,23};
	  addlast=false;
	}
      }
      if(addlast){
	ncaparams.as_bphseries.push_back(lastkmer);
      }

      ncaparams.as_numpasses=ncaparams.as_bphseries.size();
    }
  }


  cout << "Number of passes used by MIRA: " << ncaparams.as_numpasses;
  cout << "\nkmer series: " << ncaparams.as_bphseries[0];
  for(uint8 i=1; i<ncaparams.as_bphseries.size(); ++i) cout << ", " << ncaparams.as_bphseries[i];
  cout << '\n';

  if(AS_miraparams[0].getHashStatisticsParams().hs_masknastyrepeats
     && AS_miraparams[0].getHashStatisticsParams().hs_apply_digitalnormalisation){
    AS_applydiginorminpass=0;
    // first kmer >=50 or second to last pass
    for(uint8 i=0; i<ncaparams.as_bphseries.size(); ++i){
      if(ncaparams.as_bphseries[i]>=50){
	AS_applydiginorminpass=i+1;
	break;
      }
    }
    if(AS_applydiginorminpass==0){
      if(ncaparams.as_bphseries.size()>1){
	AS_applydiginorminpass=ncaparams.as_bphseries.size()-1;
      }else{
	AS_applydiginorminpass=1;
      }
    }
    cout << "Digital normalisation scheduled for pass " << AS_applydiginorminpass << endl;
  }

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Assembly::checkTerminationRequest()
{
  struct stat st;
  std::string fname=AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/terminate";
  int rc=stat(fname.c_str(),&st);
  if(rc==0) {
    std::string command="mv "
      +AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/terminate"
      +" "
      +AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/terminate_acknowledged";
    int dummy=system(command.c_str());
    return true;
  }
  return false;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::assemble()
{
  FUNCSTART("Assembly::assemble()");

#ifdef MIRAMEMORC
  cout.flush();
  MemORC::statistics();
  MemORC::checkAllMemBlocks();
#endif

  basicDataChecks();

  uint32 startpass=1;
  AS_guessedtemplatevalues=false;

  ensureStandardDirectories(false);

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  //directory_parameters const & dir_params= AS_miraparams->getDirectoryParams();
  edit_parameters const & ed_params= AS_miraparams[0].getEditParams();

  {
    std::string dummyfn(MIRAParameters::getMHSLibDir()+"/filter_default_rrna.mhs.gz");
    if(as_fixparams.as_filter_rrna && !fileExists(dummyfn)){
      MIRANOTIFY(Notify::INSTALL, "You specified -CL:frrna=yes, but MIRA could not find the default rRNA hash statistics file which should be located at " + dummyfn + "\nThis is usually the case when you forgot to run the script to install the rRNA data from the MIRA package.\n\nWhat to do:\n- if you installed from source: use the 'make install' functionality.\n- if you installed from a precompiled binaries package: go to the 'dbinstall' directory in that package and type\n  ./mira-install-sls-rrna.sh rfam_rrna-21-12.sls.gz\n");
    }
  }


  AS_donequickdenovocoveragecheck=false;

  // fill quick lookup which sequencing types are present as well as whether backbones are there
  AS_seqtypespresent.clear();
  AS_seqtypespresent.resize(ReadGroupLib::SEQTYPE_END,false);
  for(uint32 rgi=0; rgi< ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(rgid.isBackbone()){
      AS_hasbackbones=true;
    }
    if(!rgid.isRail() && !rgid.isBackbone()){
      AS_seqtypespresent[rgid.getSequencingType()]=true;
    }
  }

  if (as_fixparams.as_numrmbbreakloops <1) {
    const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams()).as_numrmbbreakloops=1;
    cout << "Number of RMB break loops <1, setting to 1\n";
  }
  if (AS_hasbackbones
      && as_fixparams.as_startbackboneusage_inpass > static_cast<int32>(as_fixparams.as_numpasses)) {
    const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams()).as_startbackboneusage_inpass=as_fixparams.as_numpasses;
    cout << "Start of backbone usage > number of passes, correcting start to " << as_fixparams.as_numpasses << endl;
  }

  if(AS_systemmemory==0 && as_fixparams.as_automemmanagement){
    cout << "Can't find info about system or process memory, switching off automatic"
      "\nmemory management.\n";
    AS_miraparams[0].getNonConstAssemblyParams().as_automemmanagement=false;
  }

  AS_assemblyinfo.setLargeContigSize(AS_miraparams[0].getSpecialParams().mi_as_largecontigsize);
  AS_assemblyinfo.setLargeContigSizeForStats(AS_miraparams[0].getSpecialParams().mi_as_largecontigsize4stats);

  if((AS_seqtypespresent[ReadGroupLib::SEQTYPE_454GS20]
      || AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]
      || AS_seqtypespresent[ReadGroupLib::SEQTYPE_IONTORRENT])
     && (as_fixparams.as_output_gap4da
	 || as_fixparams.as_output_tmp_gap4da
	 || as_fixparams.as_output_exttmp_gap4da)){
    cout << "454, Solexa or IonTorrent data present, switching off GAP4DA type output results (you *DO NOT* want millions of files in a directory, really.)\n";
    MIRAParameters::parseQuickmode("-OUT:org=no:otg=no:oetg=no",
				   "", AS_miraparams);
  }

  // if we have a nag and warn in place for Phi X 174 sequence:
  //  prepare the PhiX174 hash statistics to compare contig consensus against PhiX174
  //  (those contigs will get renamed with the given prefix)
  if(!AS_miraparams[0].getNagAndWarnParams().nw_warn_x174prefix.empty()){
    ReadPool tmprp;
    DataProcessing::addPhiX174ToReadpool(tmprp);
    std::string dummyfn(AS_miraparams[0].getDirectoryParams().dir_tmp+"/phix174.mhs.gz");
    AS_phix174hashstatistics.computeHashStatistics(tmprp,512,false,false,true,1,0,31,dummyfn,AS_miraparams[0].getDirectoryParams().dir_tmp);
  }

  // look for template ids found
  AS_maxtemplateid=-1;
  for(uint32 rpi=0; rpi<AS_readpool.size(); ++rpi){
    if(AS_readpool[rpi].getTemplateID()>AS_maxtemplateid) AS_maxtemplateid=AS_readpool[rpi].getTemplateID();
  }
  cout << "PRED MAXTID " << AS_maxtemplateid << endl;

  // allocate or reserve memory that is quite static from the size,
  //  reducing memory fragmentations at least a bit
  {
#if TRACKMEMUSAGE
    cout << "\ndmi as_init 00\n";
    dumpMemInfo();
#endif
    if(AS_maxtemplateid>=0){
      AS_templateguesses.reserve(AS_maxtemplateid+1);
    }

    AS_estnochimerakill.clear();
    AS_multicopies.reserve(AS_readpool.size());  // do NOT init multicopies!
    AS_hasmcoverlaps.reserve(AS_readpool.size()); // does not need init

    AS_permanent_overlap_bans.nuke();
    AS_permanent_overlap_bans.resize(AS_readpool.size());

    AS_used_ids.resize(AS_readpool.size(),0);
    AS_clipleft.resize(AS_readpool.size(),0);
    AS_clipright.resize(AS_readpool.size(),0);
    AS_istroublemaker.resize(AS_readpool.size(),0);
    AS_debrisreason.resize(AS_readpool.size(),0);
    AS_incorchim.resize(AS_readpool.size(),0);
    AS_isdebris.resize(AS_readpool.size(),0);
    AS_needalloverlaps.resize(AS_readpool.size(),false);
    AS_maxcoveragereached.resize(AS_readpool.size(),0);

    // these 3 only temp filled, but block anyway to reduce fragmentation
    AS_allrmbsok.reserve(AS_readpool.size());
    AS_probablermbsnotok.reserve(AS_readpool.size());
    AS_weakrmbsnotok.reserve(AS_readpool.size());
#if TRACKMEMUSAGE
    cout << "\ndmi as_init 10\n";
    dumpMemInfo();
#endif
  }

  bool clipvectorleftovers=false;
  bool usereadextension=false;
  for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++){
    if(AS_seqtypespresent[i]){
      if(AS_miraparams[i].getAssemblyParams().as_clip_possible_vectors){
	if(i != ReadGroupLib::SEQTYPE_SANGER){
	  cout << "-CL:pvlc takes effect only for Sanger sequences.\n";
	}else{
	  clipvectorleftovers=true;
	}
      }
      if(AS_miraparams[i].getAssemblyParams().as_use_read_extension){
	if(i != ReadGroupLib::SEQTYPE_SANGER){
	  cout << "-DP:ure takes effect only for Sanger sequences.\n";
	}else{
	  usereadextension=true;
	}
      }
    }
  }

  AS_resumeisok=false;
  if(AS_resumeasembly){
    AS_resumeisok=true;
    try{
      loadSnapshotData(startpass);
    }
    catch(Notify n){
      cout << "Error while loading snapshot metadata, resuming assembly is not possible, sorry\n";
      n.handleError(THISFUNC);
    }
    // set a couple of values
    if(startpass>1) AS_doneskimchimera=true;
  }else{
    preassembleTasks(usereadextension,clipvectorleftovers);
    performSnapshot(1);
  }

  EDITParameters eparams;
#ifdef MIRA_HAS_EDIT
  //  eparams.setDoEval();
  eparams.setVerbose(0);
  eparams.setShowProgress(true);
#else
#endif

  priv_setupAutoPasses();

  uint32 actpass=startpass;
  if(as_fixparams.as_assemblyjob_preprocessonly){
    if (true) {
      // I think read statistics are not needed by anyone atm, not even me
      cout << "You have selected to have read preprocessing done, so that's it I guess.\n";
    }else{
      // only preprocess?
      // Just calculate some statistics and the hash analysis

      // currently not needed I think, just to be sure
      AS_resumeisok=false;

      cout << "You have selected to have 0 passes on the assembly, therefore"
	"\njust running a couple of hash statistics and info for read repeats."
	"\n(also performing rare kmer clips if wished)\n";
      dumpSomeStatistics();
      priv_performHashAnalysis("",true,false, 0, "", "_pass");

      dumpSomeStatistics();
      performSnapshot(1);
    }
  }else{
    for(; actpass<=as_fixparams.as_numpasses; actpass++){

      //if(actpass==2) {
      //	Contig::setMasterCEBUGFlag(true);
      //}else{
      //	Contig::setMasterCEBUGFlag(false);
      //}

      //dumpMemInfo();
      dumpFile("/proc/self/status",cout);

      cout << "\n\nPass: " << actpass << " / " << as_fixparams.as_numpasses << endl;
      if(as_fixparams.as_dateoutput) dateStamp(cout);

      AS_warnings.dumpWarnings();

#ifdef VALGRIND_LEAKCHECK
      cout << "\n==MEMTRACK1 debugging start\n";
      dumpMemInfo();
      VALGRIND_DO_LEAK_CHECK
	cout << "\n==MEMTRACK1 debugging end\n";

      if(actpass==20) {
	cout << "\n==MEMTRACK exiting\n";
	dumpMemInfo();
	exit(0);
      }
#endif
#ifdef MIRAMEMORC
      cout.flush();
      MemORC::statistics();
      MemORC::checkAllMemBlocks();
#endif

      AS_contigs.clear();

#ifdef MIRA_HAS_EDIT
      if(actpass<as_fixparams.as_numpasses) {
	// eventually strict editing in first passes
	//  eparams.setStrictEvaluation(false);
	eparams.setStrictEvaluation(true);
	eparams.setConfirmationThreshold(static_cast<float>(0.8));
      } else {
	// eventually lazy editing in last passes
	eparams.setStrictEvaluation(ed_params.ed_strict_editing_mode);
	eparams.setConfirmationThreshold(static_cast<float>(ed_params.ed_confirmation_threshold/100.0));
      }
#endif

      {
	auto newbph=priv_calcBasesPerHashOnPass(actpass);
	AS_miraparams[0].getNonConstSkimParams().sk_basesperhash=newbph;

        if(AS_hasbackbones){
	  uint32 bblen=0;
	  for(auto & tc : AS_bbcontigs) bblen+=tc.getContigLength();
	  if(bblen<16000000) {
	    // do nothing for < 16MB
	  }else if(bblen<20000000){
	    if(newbph<11) {
	      newbph=11;
	      cout << "reference sequence >=16 MB, forcing minimum bph of 11.\n";
	    }
	  }else if(bblen<2500000){
	    if(newbph<12) {
	      newbph=12;
	      cout << "reference sequence >=20 MB, forcing minimum bph of 12.\n";
	    }
	  }else{
	    if(newbph<13) {
	      newbph=13;
	      cout << "reference sequence >=25 MB, forcing minimum bph of 13.\n";
	    }
	  }
	  // the following is bad, does not take length of reads into account. Oh well.
	  if(AS_readpool.size()<10000000){
	  }else if(AS_readpool.size()<15000000){
	    if(newbph<11) {
	      newbph=11;
	      cout << "more than 10m reads, forcing minimum bph of 11.\n";
	    }
	  }else if(AS_readpool.size()<20000000){
	    if(newbph<12) {
	      newbph=12;
	      cout << "more than 15m reads, forcing minimum bph of 12.\n";
	    }
	  }else{
	    if(newbph<13) {
	      newbph=13;
	      cout << "more than 20m reads, forcing minimum bph of 13.\n";
	    }
	  }
	}

	cout << "Automatic -SK:bph set to " << newbph << endl;
      }

      dumpSomeStatistics();

#if TRACKMEMUSAGE
      cout << "\ndmi 20\n";
      dumpMemInfo();
#endif

      // preparing hash analysis

      // on the last pass, maybe the user wants a rare kmer final kill

      {
	bool rkfk=(actpass==as_fixparams.as_numpasses);
	if(rkfk){
	  // resolve auto mode: only for de-novo, set "hs_rare_kmer_final_kill", set it to 3
	  // all other cases: no change
	  if(AS_miraparams[0].getNonConstHashStatisticsParams().hs_rare_kmer_final_kill==0
	     && !as_fixparams.as_assemblyjob_mapping) {
	    AS_miraparams[0].getNonConstHashStatisticsParams().hs_rare_kmer_final_kill=3;
	  }
	}
	priv_performHashAnalysis("", true, rkfk, actpass, "", "_pass");
      }

      // if wished, perform continuous PEC
      if(as_fixparams.as_clip_pec_continuous){
	std::string tmpfname;
	tmpfname=buildFileName(0,"","",
			       as_fixparams.as_tmpf_clippings,
			       ".txt","",false);

	std::string logprefix("CPEC-pass"+boost::lexical_cast<std::string>(actpass));
	auto actbph=priv_calcBasesPerHashOnPass(actpass);

	// this PEC round is influenced by how priv_performHashAnalysis() was
	//  computed above: with or without rare kmer final kill
	//
	// in the usual setting (auto, see above) this means conservative for
	//  all mapping, and all de-novo passes except the last one
	auto numfound=performProposedEndClips(actbph,tmpfname,logprefix);
      }

      // on every pass, throw out dubious reads which might be chimeras
      // atm only for denovo
      if(!as_fixparams.as_assemblyjob_mapping
	 && as_fixparams.as_clip_kmer_junkdetection){
	std::string tmpfname;
	tmpfname=buildFileName(0,"","",
			       as_fixparams.as_tmpf_clippings,
			       ".txt","",false);
	AS_dataprocessing.startLogging(tmpfname,false);

	if(actpass<as_fixparams.as_numpasses){
	  // honestly, if continuous PEC is on, this should always have zero numfound
	  //  as PEC should clip things away
	  auto numfound=AS_dataprocessing.markReadsWithInvalidKMerEndsAsChimeras_Pool(
	    AS_readpool,
	    AS_miraparams[0].getSkimParams().sk_basesperhash,
	    as_fixparams.as_clip_kmer_junkkill,  //false, just mark reads; true, kill them
	    &AS_debrisreason,
	    AS_incorchim,
	    AS_estnochimerakill,
	    "pass"+boost::lexical_cast<std::string>(actpass));
	  cout << "Found " << numfound << " incorrectible read ends or possible chimeras.\n";
	}else{
	  auto numfound=AS_dataprocessing.markReadsWithRareKMersAsChimeras_Pool(
	    AS_readpool,
	    AS_miraparams[0].getSkimParams().sk_basesperhash,
	    as_fixparams.as_clip_kmer_junkkill,  //false, just mark reads; true, kill them
	    &AS_debrisreason,
	    AS_incorchim,
	    AS_estnochimerakill,
	    "pass"+boost::lexical_cast<std::string>(actpass));
	  cout << "Found " << numfound << " terminally incorrectible reads or possible chimeras.\n";
	}
	AS_dataprocessing.stopLogging();
      }

#if TRACKMEMUSAGE
      cout << "\ndmi 30\n";
      dumpMemInfo();
#endif

      findPossibleOverlaps(actpass, "", "_pass");

      // TODO: setting that here is messy
      AS_steps[ASVECTORSCLIPPED]=0;

#if TRACKMEMUSAGE
      cout << "\ndmi 40\n";
      dumpMemInfo();
#endif

      {
	std::string signalfile(buildFileName(actpass, "", "_pass",
					AS_miraparams[0].getAssemblyParams().as_tmpf_signal_mainalignments,
					".ok"));
	if(!AS_resumeasembly || !AS_resumeisok || !fileExists(signalfile)){
	  AS_resumeisok=false;
	  makeAlignments(Assembly::ma_takeall, false, true, actpass, "", "_pass");
	  saveResumeDataMA(actpass, "", "_pass");
	  std::ofstream fout(signalfile);  // create checkpoint signal file for main alignments
	}else{
	  cout << "Resume assembly: alignments already present, good.\n";
	  loadResumeDataMA(actpass, "", "_pass");
	}
      }

#if TRACKMEMUSAGE
      cout << "\ndmi 50\n";
      dumpMemInfo();
#endif

      priv_loadAlignmentsFromFile(actpass, "", "_pass");

#if TRACKMEMUSAGE
      cout << "\ndmi 60\n";
      dumpMemInfo();
#endif

      // count all SRMr tags in reads
      // idea: after building contigs, run the repeat resolver with
      //  all reads that have new tags

      std::vector<uint32> xrmrcount(AS_readpool.size(),0);
      for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	xrmrcount[rnr]=AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idSRMr);
      }

#if TRACKMEMUSAGE
      cout << "\ndmi 70\n";
      dumpMemInfo();
#endif

#ifdef VALGRIND_LEAKCHECK
      cout << "\n==MEMTRACK2 debugging start\n";
      dumpMemInfo();
      VALGRIND_DO_LEAK_CHECK
	cout << "\n==MEMTRACK2 debugging end\n";
#endif

      if(as_fixparams.as_dateoutput) dateStamp(cout);

      bool foundrepeats=buildFirstContigs(actpass,
					  eparams,
					  (actpass==as_fixparams.as_numpasses));

      if(!AS_bbcontigs.empty()){
	// mapping may have created new readgroups (CER reads)
	AS_current_rls_byrg.resize(ReadGroupLib::getNumReadGroups());
      }

#ifdef VALGRIND_LEAKCHECK
      cout << "\n==MEMTRACK3 debugging start\n";
      dumpMemInfo();
      VALGRIND_DO_LEAK_CHECK
	cout << "\n==MEMTRACK3 debugging end\n";
#endif

#if TRACKMEMUSAGE
      cout << "\ndmi 80\n";
      dumpMemInfo();
#endif

      // try to guess a good value for minimum coverage
      //  for large contigs
      AS_assemblyinfo.calcCurrentInfo();
      {
	auto avc=AS_assemblyinfo.ASI_avgcoverage[1];
	if(avc==0){
	  avc=AS_assemblyinfo.ASI_avgcoverage[0];
	}
	double multval=0.5;
	if(avc<40) multval=static_cast<double>(1.0)/3;
	AS_assemblyinfo.setLargeTotalCov(
	  static_cast<uint32>(static_cast<double>(.5)+avc*multval)
	  );
	for(uint8 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	  avc=AS_assemblyinfo.ASI_avgcoverage_perst[1][st];
	  if(avc==0) avc=AS_assemblyinfo.ASI_avgcoverage_perst[0][st];
	  AS_assemblyinfo.setLargeContigCovPerST(
	    static_cast<uint32>(static_cast<double>(.5)+avc*multval),
	    st
	    );
	}
      }
      saveAssemblyInfo();

      // save large contigs info only for genome de-novo,
      //  i.e., no backbones present and we're using genomic pathfinder (not EST)
      // ah, and hunt for smile shape
      if(AS_bbcontigs.empty() && AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	saveLargeContigsInfo();
	warnAtSmileCoverage();
      }

      // this also sets the info for RGSTInfo
      warnChimeraContent();

      if(actpass==as_fixparams.as_numpasses) {
	saveRGSTInfo();
      }else{
	saveRGSTInfo(actpass, "", "_pass");
      }


#if TRACKMEMUSAGE
      cout << "\ndmi 90\n";
      dumpMemInfo();
#endif

      // some things are worth doing only if it's not the last pass
      if(actpass!=as_fixparams.as_numpasses){

	if(foundrepeats){
	  cout << "Repeats found during contig building, adding additional alignment iteration\nfor quick repeat resolving.\n";
	  AS_steps[ASADSLISTOK]=0;
	  //makeAlignments(Assembly::ma_needSRMrOrTwoCRMr, true, actpass, "", "", "repeat_resolve");

	  // Step 1 of repeat resolving
	  // define which reads should be taken into the alignments phase
	  //  (== those who got additional SRMr tags)
	  AS_readsforrepeatresolve.clear();
	  AS_readsforrepeatresolve.resize(AS_readpool.size(),false);
	  for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	    if(AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idSRMr) != xrmrcount[rnr]) AS_readsforrepeatresolve[rnr]=true;
	  }

	  // prepare for second round of repeat resolve
	  for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	    xrmrcount[rnr]=AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idCRMr);
	  }

#if TRACKMEMUSAGE
	  cout << "\ndmi a0\n";
	  dumpMemInfo();
#endif
	  makeAlignments(Assembly::ma_needRRFlag, true, false, actpass, "", "", "repeat_resolve");
	  AS_steps[ASADSLISTOK]=0;

	  // Step 2 of repeat resolving
	  // define which reads should be taken into the alignments phase
	  //  (== those who got additional CRMr tags in alignment phase above)

	  AS_readsforrepeatresolve.clear();
	  AS_readsforrepeatresolve.resize(AS_readpool.size(),false);
	  for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	    if(AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idCRMr) != xrmrcount[rnr]) AS_readsforrepeatresolve[rnr]=true;
	  }

#if TRACKMEMUSAGE
	  cout << "\ndmi b0\n";
	  dumpMemInfo();
#endif
	  makeAlignments(Assembly::ma_needRRFlagAndBothCRMr, true, false, actpass, "", "", "repeat_resolve");

#if TRACKMEMUSAGE
	  cout << "\ndmi c0\n";
	  dumpMemInfo();
#endif
	  nukeSTLContainer(AS_readsforrepeatresolve);

#if TRACKMEMUSAGE
	  cout << "\ndmi d0\n";
	  dumpMemInfo();
#endif
	}

	if(usereadextension
	   && actpass >= as_fixparams.as_readextension_firstpassnum
	   && actpass <= as_fixparams.as_readextension_lastpassnum) {
	  extendADS(actpass, "", "_pass");
	}
	if(clipvectorleftovers) {
	  performSeqVectorClippings();
	}
      }

      if(checkTerminationRequest()){
	cout << "Seen termination request by user\n";
	if(actpass + 2 < as_fixparams.as_numpasses){
	  assembly_parameters & ap=const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams());
	  cout << "Changing number of passes from " << ap.as_numpasses << " to ";
	  ap.as_numpasses=actpass+2;
	  cout << ap.as_numpasses << endl;
	}else{
	  cout << "No further action necessary, will terminate anyway in 2 passes max.\n";
	}
      }

#if TRACKMEMUSAGE
      cout << "\ndmi e0\n";
      dumpMemInfo();
#endif

      performSnapshot(actpass+1);

    }

    if(AS_hasbackbones && AS_guessedtemplatevalues){
      priv_hackMergeTwoResultMAFs();
    }

  }

  AS_warnings.dumpWarnings();


#ifdef MIRAMEMORC
  cout.flush();
  MemORC::statistics();
  MemORC::checkAllMemBlocks();
#endif

  FUNCEND();
  return;
}


/*************************************************************************
 *
 * Given a pass number, return the bph (kmer) size which will be
 * used in Skim
 *
 *************************************************************************/

uint32 Assembly::priv_calcBasesPerHashOnPass(uint32 passnr) const
{
  if(passnr>0) --passnr;
  if(passnr>=AS_miraparams[0].getAssemblyParams().as_bphseries.size()){
    passnr=AS_miraparams[0].getAssemblyParams().as_bphseries.size()-1;
  }
  return AS_miraparams[0].getAssemblyParams().as_bphseries[passnr];
}

/*************************************************************************
 *
 * In mapping assemblies, the results for guessing template numbers
 *  (size, segment placement) are available only after the MAF has been written
 * Therefore, need to rewrite header of MAF (and CAF completely)
 * Easiest way out: MAF header from chkpoint, body of 'usual' result MAF
 *
 * delete CAF, have outer caller recreate it from new MAF
 *
 *************************************************************************/

void Assembly::priv_hackMergeTwoResultMAFs()
{
  std::string headermaf(buildDefaultCheckpointFileName(AS_miraparams[0].getFileParams().chkpt_readpool));
  if(!fileExists(headermaf)) return;

  std::string bodymaf(getMAFFilename());
  if(!fileExists(bodymaf)) return;

  std::string newmaf(bodymaf+"tmp");

  std::ifstream fin(headermaf, std::ios::in);
  if(!fin.is_open()){
    cout << "MAFmerge: Could not open headermaf " << headermaf << endl;
    return;
  }

  std::ofstream fout(newmaf,std::ios::out);
  if(!fout.is_open()){
    cout << "MAFmerge: Could not open newmaf " << newmaf << endl;
    return;
  }

  std::string actline;
  actline.reserve(1000);
  while(!fin.eof()){
    getline(fin,actline);
    if(!actline.empty()){
      if(actline[0]=='@'
	 || actline[0]=='#'){
	fout << actline << '\n';
      }else{
	break;
      }
    }
  }
  fin.close();
  fin.open(bodymaf, std::ios::in);
  if(!fin.is_open()){
    cout << "MAFmerge: Could not open bodymaf " << bodymaf << endl;
    return;
  }
  while(!fin.eof()){
    getline(fin,actline);
    if(actline.empty()) continue;
    if(actline[0]!='@'
       && actline[0]!='#') break;
  }
  fout << actline << '\n';
  fout << fin.rdbuf();

  if(fout.bad()){
    cout << "MAFmerge: Could not finish copying to " << newmaf << endl;
  }
  fin.close();
  fout.close();

  fileRename(newmaf,bodymaf);
  std::string caffile(getCAFFilename());
  if(fileExists(caffile)){
    fileRemove(caffile,true);
    fout.open(AS_miraparams[0].getDirectoryParams().dir_results+"/_tmprecreate");
    fout << caffile << endl;
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::performSnapshot(uint32 actpass)
{
  FUNCSTART("void Assembly::performSnapshot(uint32 actpass)");

  auto const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  auto const & ffp= AS_miraparams[0].getFileParams();

  std::list<std::string> filesanew;

  auto & dirchkpt=AS_miraparams[0].getDirectoryParams().dir_checkpoint;
  auto & dirtmpchkpt=AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp;

  cout << "Performing snapshot " << actpass << endl;
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  try{
    boost::filesystem::rename(dirchkpt,dirtmpchkpt);
  }
  catch(boost::filesystem::filesystem_error fse){
    MIRANOTIFY(Notify::FATAL,"Could not rename checkpoint directory in preparation of saving new checkpoint?\nError message is " << fse.what());
  }
  catch(...){
    MIRANOTIFY(Notify::FATAL,"Could not rename checkpoint directory in preparation of saving new checkpoint?\nUnspecified exception?");
  }
  if(ensureDirectory(dirchkpt, false)){
    // Oooops? Cannot create new? Well, roll back the rename and exit
    try{
      boost::filesystem::rename(dirtmpchkpt,dirchkpt);
    }
    catch(boost::filesystem::filesystem_error fse){
      MIRANOTIFY(Notify::FATAL,"Could not rollback initial checkpoint rename?\nError message is " << fse.what());
    }
    catch(...){
      MIRANOTIFY(Notify::FATAL,"Could not rollback initial checkpoint rename?\nUnspecified exception?");
    }
    MIRANOTIFY(Notify::FATAL,"Could not create new snapshot directory? Disk full? changed permissions?");
  }
  try{
    ssdReadPool(buildDefaultCheckpointFileName(ffp.chkpt_readpool));
    ssdPassInfo(buildDefaultCheckpointFileName(ffp.chkpt_passinfo),actpass);
    ssdMaxCovReached(buildDefaultCheckpointFileName(ffp.chkpt_maxcovreached));
    ssdBannedOverlaps(buildDefaultCheckpointFileName(ffp.chkpt_bannedoverlaps));

    filesanew.push_back(ffp.chkpt_readpool);
    filesanew.push_back(ffp.chkpt_passinfo);
    filesanew.push_back(ffp.chkpt_maxcovreached);
    filesanew.push_back(ffp.chkpt_bannedoverlaps);
  }
  catch(...){
    // Oooops? Error while writing new checkpoint data? Rollback ...
    try{
      boost::filesystem::remove_all(dirtmpchkpt);
    }
    catch(boost::filesystem::filesystem_error fse){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not delete failed checkpoint directory. Contact the author.\nError message is " << fse.what());
    }
    catch(...){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not delete failed checkpoint directory. Contact the author.\nnUnspecified exception?");
    }
    try{
      boost::filesystem::rename(dirtmpchkpt,dirchkpt);
    }
    catch(boost::filesystem::filesystem_error fse){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not rollback checkpoint rename after failed snapshot?\nError message is " << fse.what());
    }
    catch(...){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not rollback checkpoint rename after failed snapshot?\nUnspecified exception?");
    }
    MIRANOTIFY(Notify::FATAL,"Could not correctly write new snapshot data. Disk full? Changed permissions?");
  }

  // so we've written anew a couple of files, let's delete the old versions in the tmp checkpoint and move the remaining files
  try{
    for(auto & fn : filesanew){
      fileRemove(dirtmpchkpt+"/"+fn,true);
    }

    // now move the files remaining in the tmp chkpoint to the chkpoint
    // then remove checkpoint tmp dir

    for(boost::filesystem::directory_iterator fI(dirtmpchkpt);
        fI != boost::filesystem::directory_iterator(); ++fI){
      //cout << "renaming " << fI->path().std::string() << " to " << (dirchkpt+"/"+fI->path().filename().string()) << endl;
      fileRename(fI->path().string(),dirchkpt+"/"+fI->path().filename().string());
    }

    boost::filesystem::remove_all(dirtmpchkpt);
  }
  catch(...){
  }

  if(as_fixparams.as_dateoutput) dateStamp(cout);

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadSnapshotData(uint32 & actpass)
{
  FUNCSTART("void Assembly::loadSnapshotData(uint32 & actpass)");

  auto const & ffp= AS_miraparams[0].getFileParams();

  // readpool was loaded before, don't load here
  actpass=lsdPassInfo(buildDefaultCheckpointFileName(ffp.chkpt_passinfo));
  lsdMaxCovReached(buildDefaultCheckpointFileName(ffp.chkpt_maxcovreached));
  lsdBannedOverlaps(buildDefaultCheckpointFileName(ffp.chkpt_bannedoverlaps));

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::ssdReadPool(const std::string & filename)
{
  FUNCSTART("void Assembly::ssdReadPool(const std::string & filename)");
  std::ofstream fout(filename,std::ios::out|std::ios::trunc);
  Contig::dumpMAF_Head(fout);
  AS_readpool.dumpAs(fout,Read::AS_MAF,true);
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot readpool?");
  }
}

void Assembly::ssdPassInfo(const std::string & filename, uint32 actpass)
{
  FUNCSTART("void Assembly::ssdPassInfo(const std::string & filename, uint32 actpass)");
  std::ofstream fout(filename,std::ios::out|std::ios::trunc);
  fout << actpass << endl;
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot actpass?");
  }
}

void Assembly::ssdMaxCovReached(const std::string & filename)
{
  FUNCSTART("void Assembly::ssdMaxCovReached(const std::string & filename)");
  std::ofstream fout(filename,std::ios::out|std::ios::trunc);
  for(auto mc : AS_maxcoveragereached) fout << mc << endl;
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot maxcov?");
  }
}

void Assembly::ssdBannedOverlaps(const std::string & filename)
{
  FUNCSTART("void Assembly::ssdBannedOverlaps(const std::string & filename)");
  // 1st line: size of AS_permanent_overlap_bans
  // nth line: id1 id2 id2 id2 ...
  std::ofstream fout(filename,std::ios::out|std::ios::trunc);
  fout << AS_permanent_overlap_bans.size() << endl;
  for(size_t rid=0; rid < AS_permanent_overlap_bans.size(); ++rid){
    if(!AS_permanent_overlap_bans.bop[rid].empty()){
      fout << rid;
      for(auto id2 : AS_permanent_overlap_bans.bop[rid]){
	fout << '\t' << id2;
      }
      fout << '\n';
    }
  }
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot banned overlaps?");
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint32 Assembly::lsdPassInfo(const std::string & filename)
{
  FUNCSTART("uint32 Assembly::lsdPassInfo(const std::string & filename)");

  std::ifstream fin(filename, std::ios::in);
  if(!fin.good()) {
    MIRANOTIFY(Notify::FATAL,"Did not find " << filename);
  }
  std::string tmp;
  fin >> tmp;
  auto n=atoll(tmp.c_str());
  if(n<0){
    MIRANOTIFY(Notify::FATAL,"negative value in " << filename << " is not expected");
  }
  return static_cast<uint32>(n);
}

void Assembly::lsdMaxCovReached(const std::string & filename)
{
  FUNCSTART("uint32 Assembly::lsdMaxCovReached(const std::string & filename)");

  std::ifstream fin(filename, std::ios::in);
  if(!fin.good()) {
    MIRANOTIFY(Notify::FATAL,"Did not find " << filename);
  }
  AS_maxcoveragereached.clear();
  std::string tmpline;
  while(getline(fin,tmpline)){
    boost::trim(tmpline);
    if(tmpline.empty()){
      MIRANOTIFY(Notify::FATAL,"empty line in " << filename << " is not expected");
    }
    auto n=atoll(tmpline.c_str());
    AS_maxcoveragereached.push_back(static_cast<uint32>(n));
  }
}

void Assembly::lsdBannedOverlaps(const std::string & filename)
{
  FUNCSTART("void Assembly::lsdBannedOverlaps(const std::string & filename)");

  std::ifstream fin(filename, std::ios::in);
  if(!fin.good()) {
    MIRANOTIFY(Notify::FATAL,"Did not find " << filename);
  }

  size_t bopsize=0;
  std::string tmpline;
  std::vector<std::string> substrs;
  if(getline(fin,tmpline)){
    boost::trim(tmpline);
    if(tmpline.empty()){
      MIRANOTIFY(Notify::FATAL,"empty first line in " << filename << " ??");
    }
    boost::split(substrs, tmpline, boost::is_any_of(" \t"));
    if(substrs.size()>1){
      MIRANOTIFY(Notify::FATAL,"first line in " << filename << " should have one element only");
    }
    auto n=atoll(substrs[0].c_str());
    if(n<0){
      cout << "Faulty first line, has negative values:\n" << tmpline << endl;
      MIRANOTIFY(Notify::FATAL,"negative value in " << filename << " is not expected");
    }
    bopsize=static_cast<size_t>(n);
  }else{
    MIRANOTIFY(Notify::FATAL,"Error reading first line of " << filename);
  }
  if(bopsize==0){
    MIRANOTIFY(Notify::FATAL,"Ooooooops, first line should have positive number: " << filename);
  }
  if(bopsize!=AS_readpool.size()){
    MIRANOTIFY(Notify::FATAL,"Ooooooops, : size of read pool in " << filename  << " (" << bopsize << ") is not equal to size of current readpool (" << AS_readpool.size() << ") ???");
  }
  AS_permanent_overlap_bans.nuke();
  AS_permanent_overlap_bans.resize(bopsize);
  std::vector<uint32> tmpvec;
  while(getline(fin,tmpline)){
    substrs.clear();
    boost::trim(tmpline);
    if(tmpline.empty()){
      MIRANOTIFY(Notify::FATAL,"empty line in " << filename << " is not expected");
    }
    boost::split(substrs, tmpline, boost::is_any_of(" \t"));
    if(substrs.size()<2){
      cout << "Faulty line, need more elements:\n" << tmpline << endl;
      MIRANOTIFY(Notify::FATAL,"empty line in " << filename << " is not expected");
    }
    auto sI=substrs.begin();
    size_t rid=0;
    {
      auto n=atoll(sI->c_str());
      if(n<0){
	cout << "Faulty line, has negative first value:\n" << tmpline << endl;
	MIRANOTIFY(Notify::FATAL,"Ooooooops, first value should value >= 0: " << filename);
      }
      rid=static_cast<size_t>(n);
    }
    ++sI;
    for(; sI!=substrs.end(); ++sI){
      auto n=atoll(sI->c_str());
      if(n<0){
	cout << "Faulty line, has negative values:\n" << tmpline << endl;
	MIRANOTIFY(Notify::FATAL,"negative value in " << filename << " is not expected");
      }
      tmpvec.push_back(static_cast<uint32>(n));
    }
    AS_permanent_overlap_bans.bop[rid].swap(tmpvec);
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveResumeDataFPO(int32 version, const std::string & prefix, const std::string & postfix)
{
  FUNCSTART("void Assembly::saveResumeDataFPO()");

  std::string filename(buildFileName(version,
					 prefix,
					 postfix,
					 AS_miraparams[0].getAssemblyParams().as_tmpf_wellconnected,
					 ".bin"));

  bool allok=true;
  if(saveVector(AS_wellconnected,filename)){
    filename=buildFileName(version, "", "_pass",
			   AS_miraparams[0].getAssemblyParams().as_tmpf_debrisreason,
			   ".bin");
    if(saveVector(AS_debrisreason,filename)){
      filename=buildFileName(version, "", "_pass",
			     AS_miraparams[0].getAssemblyParams().as_tmpf_skimmegahubs,
			     ".bin");
      if(saveVector(AS_skimmegahubs,filename)){
	allok=true;
      }
    }
  }
  if(!allok){
    MIRANOTIFY(Notify::FATAL, "Error while writing file " << filename << ". Is the disk full? Are permissions right?");
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadResumeDataFPO(int32 version, const std::string & prefix, const std::string & postfix)
{
  FUNCSTART("void Assembly::loadResumeDataFPO()");

  std::string filename(buildFileName(version,
				prefix,
				postfix,
				AS_miraparams[0].getAssemblyParams().as_tmpf_wellconnected,
				".bin"));
  bool allok=false;
  AS_wellconnected.clear();
  AS_debrisreason.clear();
  AS_skimmegahubs.clear();
  if(loadVector(AS_wellconnected,filename,0)){
    filename=buildFileName(version, "", "_pass",
			   AS_miraparams[0].getAssemblyParams().as_tmpf_debrisreason,
			   ".bin");
    if(loadVector(AS_debrisreason,filename,0)){
      filename=buildFileName(version, "", "_pass",
			     AS_miraparams[0].getAssemblyParams().as_tmpf_skimmegahubs,
			     ".bin");
      if(loadVector(AS_skimmegahubs,filename,0)){
	allok=true;
      }
    }
  }
  if(allok){
    // checks?
  }else{
    MIRANOTIFY(Notify::FATAL, "Error while reading file " << filename << ". Is the file present and correct? Are permissions right?");
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveResumeDataMA(int32 version, const std::string & prefix, const std::string & postfix)
{
  FUNCSTART("void Assembly::saveResumeDataMA()");

  //  need to save current overlap bans to tmp as makeAlignments could have created new ones
  ssdBannedOverlaps(buildFileName(version, "", "_pass",
				  AS_miraparams[0].getAssemblyParams().as_tmpf_banned_overlaps,
				  ".txt"));

  std::string filename(buildFileName(version, "", "_pass",
				AS_miraparams[0].getAssemblyParams().as_tmpf_istroublemaker,
				".bin"));
  bool allok=false;
  if(saveVector(AS_istroublemaker,filename)){
    filename=buildFileName(version, "", "_pass",
			   AS_miraparams[0].getAssemblyParams().as_tmpf_needalloverlaps,
			   ".bin");
    if(saveVector(AS_needalloverlaps,filename)){
      filename=buildFileName(version, "", "_pass",
			     AS_miraparams[0].getAssemblyParams().as_tmpf_multicopies,
			     ".bin");
      if(saveVector(AS_multicopies,filename)){
	filename=buildFileName(version, "", "_pass",
			       AS_miraparams[0].getAssemblyParams().as_tmpf_hasmcoverlap,
			       ".bin");
	if(saveVector(AS_hasmcoverlaps,filename)){
	  allok=true;
	}
      }
    }
  }
  if(!allok){
    MIRANOTIFY(Notify::FATAL, "Error while writing file " << filename << ". Is the disk full? Are permissions right?");
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadResumeDataMA(int32 version, const std::string & prefix, const std::string & postfix)
{
  FUNCSTART("void Assembly::loadResumeDataMA()");
  //  need to load banned overlaps
  lsdBannedOverlaps(buildFileName(version, "", "_pass",
				  AS_miraparams[0].getAssemblyParams().as_tmpf_banned_overlaps,
				  ".txt"));

  std::string filename(buildFileName(version, "", "_pass",
				AS_miraparams[0].getAssemblyParams().as_tmpf_istroublemaker,
				".bin"));
  bool allok=false;
  AS_istroublemaker.clear();
  AS_needalloverlaps.clear();
  AS_multicopies.clear();
  AS_hasmcoverlaps.clear();
  if(loadVector(AS_istroublemaker,filename,0)){
    filename=buildFileName(version, "", "_pass",
			   AS_miraparams[0].getAssemblyParams().as_tmpf_needalloverlaps,
			   ".bin");
    if(loadVector(AS_needalloverlaps,filename,0)){
      filename=buildFileName(version, "", "_pass",
			     AS_miraparams[0].getAssemblyParams().as_tmpf_multicopies,
			     ".bin");
      if(loadVector(AS_multicopies,filename,0)){
	filename=buildFileName(version, "", "_pass",
			       AS_miraparams[0].getAssemblyParams().as_tmpf_hasmcoverlap,
			       ".bin");
	if(loadVector(AS_hasmcoverlaps,filename,0)){
	  allok=true;
	}
      }
    }
  }
  if(allok){
    // checks?
  }else{
    MIRANOTIFY(Notify::FATAL, "Error while reading file " << filename << ". Is the file present and correct? Are permissions right?");
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::makeNewReadPoolFromContigs()
{
  FUNCSTART("void Assembly::makeNewReadPoolFromContigs()");


  int32 ccounter=0;
  for(auto cI=AS_contigs.begin(); cI!=AS_contigs.end(); ++cI, ++ccounter){

    auto & cr = cI->getContigReads();

    for(auto crI = cr.begin(); crI!=cr.end(); ++crI){
      if(crI.getORPID() >= 0){
	if(crI->checkRead()){
	  cout << "Precheck failed: " << endl;
	  cout << *crI;
	  MIRANOTIFY(Notify::FATAL, crI->checkRead()) ;
	}

	AS_readpool[crI.getORPID()]=*crI;

	if(AS_readpool[crI.getORPID()].checkRead()){
	  cout << "Postcheck1 failed: " << endl;
	  cout << AS_readpool[crI.getORPID()];
	  MIRANOTIFY(Notify::FATAL, AS_readpool[crI.getORPID()].checkRead()) ;
	}

	const_cast<Read &>(*crI).discard();

	if(AS_readpool[crI.getORPID()].checkRead()){
	  cout << "Postcheck2 failed: " << endl;
	  cout << AS_readpool[crI.getORPID()];
	  MIRANOTIFY(Notify::FATAL, AS_readpool[crI.getORPID()].checkRead()) ;
	}
      }
    }
  }
  AS_contigs.clear();

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/


bool Assembly::markRepeats(Contig & con, std::vector<bool> & readsmarkedsrm, Contig::repeatmarker_stats_t & repstats)
{
  FUNCSTART("void Assembly::markRepeats(Contig & con)");

  repstats.init();

  cout << "Marking possibly misassembled repeats: ";
  cout.flush();

  con.newMarkPossibleRepeats(repstats, readsmarkedsrm);
  cout << "done step 1, starting step 2:";
  con.codonSingleBaseRepeatMarker(6,repstats, readsmarkedsrm);
  if(repstats.numSRMs>0 || repstats.numWRMs>0 || repstats.numSNPs>0){
    cout << "\nFound\n";
    cout << " - " << repstats.numSRMs << " new Strong RMB (SRMc)\n";
    cout << " - " << repstats.numWRMs << " new Weak RMB (WRMc)\n";
    cout << " - " << repstats.numSNPs << " SNP\npositions tagged.";
    cout.flush();
  }else{
    cout << "done. Found none." << endl;
  }

  FUNCEND();

  return repstats.numSRMs>0;
}





/*************************************************************************
 *
 *
 *
 *************************************************************************/

////#define CEBUG(bla)   {cout << bla; cout.flush(); }
//void Assembly::transferContigReadsToReadpool(const Contig & buildcon, std::list<Contig::pbdse_t> & pbdsev, int32 passnr)
//{
//  FUNCSTART("void Assembly::transferContigReadsToReadpool(const Contig & buildcon, std::vector<Contig::pbdse_t> & pbdsev)");
//
//  cout << "Transfering reads to readpool." << endl;
//
//  BUGIFTHROW(true,"need redo for PlacedContigReads");
//
//  const std::vector<Contig::contigread_t> & cr = buildcon.getContigReads();
//
//  // split one list pbdsev into sublists for each contigread if needed
//
//  CEBUG("tCRTR 1" << endl);
//  std::vector<std::list<Contig::pbdse_t> > vpbdsev;
//  if(!pbdsev.empty()){
//    vpbdsev.resize(AS_readpool.size());
//    for(; pbdsev.begin() != pbdsev.end(); ){
//      //cout << "ls: " << pbdsev.size() << "\tSplicing: " << pbdsev.front();
//      size_t trid=pbdsev.front().rid;
//      vpbdsev[trid].splice(vpbdsev[trid].begin(),pbdsev,pbdsev.begin());
//      //cout << "ts: " << vpbdsev[trid].size() << endl;
//    }
//
//    //for(uint32 vi=0;vi<vpbdsev.size(); ++vi){
//    //  cout << "vi: " << vi << endl;
//    //  std::list<Contig::pbdse_t>::const_iterator pI=vpbdsev[vi].begin();
//    //  for(; pI!=vpbdsev[vi].end(); ++pI){
//    //	cout << "\t" << *pI;
//    //  }
//    //}
//  }
//
//  CEBUG("tCRTR 2" << endl);
//  // create a temporary read with enough capacity to hold the
//  //  largest of the reads to transfer (to prevent re-allocation)
//  //
//  // used for PacBio where elastic N stretches need to
//  //  be corrected
//  // also used to remove the gaps from the reads
//  std::vector<Contig::contigread_t>::const_iterator crI = cr.begin();
//  Read tmpr;
//  {
//    uint32 reservelen=0;
//    for(;crI!=cr.end();crI++){
//      if(crI->orpid>=0){
//	reservelen=std::max(reservelen,crI->read.getLenSeq());
//      }
//    }
//    tmpr.reserve(reservelen);
//  }
//
//  CEBUG("tCRTR 3" << endl);
//  // now copy all contig reads to read pool, one by one, using the tmp read
//  //  and performing edits on this tmp read, leaving contig reads untouched
//  crI = cr.begin();
//  for(;crI!=cr.end();crI++){
//    if(crI->orpid>=0){
//      tmpr=crI->read;
//
//      // carefull, maybe there is nothing in the vpbdsev vector!
//      //  then one would get segfaults
//      if(tmpr.isSequencingType(ReadGroupLib::SEQTYPE_PACBIOHQ) && !vpbdsev.empty()){
//	std::list<Contig::pbdse_t>::const_iterator pRI=vpbdsev[crI->orpid].begin();
//	for(; pRI != vpbdsev[crI->orpid].end(); ++pRI){
//	  CEBUG("Apply: " << *pRI);
//	  BUGIFTHROW(pRI->rid >= AS_readpool.size(), "pRI->rid (" << pRI->rid << ") >= AS_readpool.size() ?");
//	  //Read::setCoutType(Read::AS_TEXT);
//	  //cout << CON_reads[pRI->cri].read;
//	  tmpr.correctNStretch(pRI->rposs,
//			       pRI->rpose,
//			       pRI->changeestim);
//	}
//
//	if(passnr==4){
//	  tmpr.transformGapsToNs();
//	}
//      }
//
//      tmpr.removeGapsFromRead();
//      AS_readpool[crI->orpid]=tmpr;
//    }
//  }
//
//  FUNCEND();
//}
////#define CEBUG(bla)


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::transferContigReadsToReadpool(const Contig & buildcon, std::list<Contig::pbdse_t> & pbdsev, int32 passnr)
{
  FUNCSTART("void Assembly::transferContigReadsToReadpool(const Contig & buildcon, std::vector<Contig::pbdse_t> & pbdsev)");

  cout << "Transfering reads to readpool." << endl;

  // create a temporary read with enough capacity to hold the
  //  largest of the reads to transfer (to prevent re-allocation)
  //
  // used for PacBio where elastic N stretches need to
  //  be corrected
  // also used to remove the gaps from the reads
  Read tmpr;
  {
    uint32 reservelen=0;
    for(auto pcrI=buildcon.getContigReads().begin(); pcrI!=buildcon.getContigReads().end(); ++pcrI){
      if(pcrI.getORPID() >= 0){
	reservelen=std::max(reservelen,pcrI->getLenSeq());
      }
    }
    tmpr.reserve(reservelen);
  }

  CEBUG("tCRTR 3" << endl);
  // now copy all contig reads to read pool, one by one, using the tmp read
  //  and performing edits on this tmp read, leaving contig reads untouched
  for(auto pcrI=buildcon.getContigReads().begin(); pcrI!=buildcon.getContigReads().end(); ++pcrI){
    if(pcrI.getORPID() >= 0){
      tmpr=*pcrI;
      tmpr.removeGapsFromRead();
      AS_readpool[pcrI.getORPID()]=tmpr;
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
void Assembly::transferContigReadTagsToReadpool(const Contig & con, const std::list<Contig>::const_iterator bbContigI)
{
  FUNCSTART("void Assembly::transferContigReadTagsToReadpool(const Contig & con, const std::list<Contig>::const_iterator bbContigI)");

  cout << "Transfering contig read tags to readpool reads." << endl;

  auto & cr = con.getContigReads();
  auto pcrI = cr.begin();
  uint32 tagnumber=0;

  multitag_t tmpmt;

  try{
    // Transfer the contigread RMB tags into the readpool only!
    // Go through all the contigreads, if they have SRMr or WRMr tags,
    //  check if they are at an edited place.
    // If not, transfer the tag to the readpool (if not already present
    //  there)
    for(;pcrI!=cr.end(); ++pcrI){
      if(pcrI.getORPID()==-1) continue;
      uint32 numtags=pcrI->getNumOfTags();

      //CEBUGF(pcrI->getName() << " has " << numtags << " tags.\n");

      for(tagnumber=0; tagnumber < numtags; tagnumber++) {
	const multitag_t & acttag=pcrI->getTag(tagnumber);
	if(acttag.identifier==Read::REA_tagentry_idSRMr
	   ||acttag.identifier==Read::REA_tagentry_idWRMr
	   ||acttag.identifier==Read::REA_tagentry_idSROr
	   ||acttag.identifier==Read::REA_tagentry_idSIOr
	   ||acttag.identifier==Read::REA_tagentry_idSAOr
	  ) {

	  tmpmt=acttag;

	  //CEBUGF("Tag " << tagnumber << " at " << acttag.from << " is " << acttag.identifier << endl);

	  bool foundedit=false;
	  for(uint32 i=0; i<numtags; i++) {
	    if(pcrI->getTag(i).from==acttag.from
	       && (pcrI->getTag(i).identifier==Read::REA_defaulttag_ED_C.identifier
		   || pcrI->getTag(i).identifier==Read::REA_defaulttag_ED_I.identifier
		   || pcrI->getTag(i).identifier==Read::REA_defaulttag_ED_D.identifier
		 )
	      ) {
	      //CEBUGF("Found " << pcrI->getTag(i).identifier << " at that position, skipping!");
	      foundedit=true;
	      break;
	    }
	  }
	  if(foundedit) continue;


	  // the next if is true for tags at gap bases;
	  //if(adj_tagpos < 0 ) {
	  //  cout << "Tag at gap base!!!!!\n";
	  //  continue;
	  //}

	  if(bbContigI!=AS_bbcontigs.end()){
	    if(pcrI->isBackbone() || pcrI->isRail()) {
	      //cout << "IsBackbone or read\n";

	      int32 adj_tagposl=pcrI->getAdjustmentPosOfReadPos(acttag.from);
	      int32 adj_tagposu=pcrI->getAdjustmentPosOfReadPos(acttag.to);
	      // don't do that for gaps in backbones or reads
	      if(adj_tagposl < 0 && adj_tagposu < 0) continue;
	      if(adj_tagposl < 0) {
		adj_tagposl=adj_tagposu-(acttag.to-acttag.from);
		if(adj_tagposl<0) adj_tagposl=0;
	      }else if(adj_tagposu < 0) {
		adj_tagposu=adj_tagposl+(acttag.to-acttag.from);
		if(adj_tagposu > static_cast<int32>(pcrI->getLenSeq()-1)) adj_tagposu=pcrI->getLenSeq()-1;
	      }

	      // This is an ugly, ugly hack: we set the P|WRMB tag in
	      //  the read of the backbonecontig.
	      // This should be a no no, but it's needed to work
	      //  with backbones

	      int32 bbcrtagposl=pcrI->getReadPosOfAdjustmentPos(adj_tagposl);
	      int32 bbcrtagposu=pcrI->getReadPosOfAdjustmentPos(adj_tagposu);
	      if(bbcrtagposl >= 0 && bbcrtagposu>=0) {
		try {
		  tmpmt.from=bbcrtagposl;
		  tmpmt.to=bbcrtagposu;
		  const_cast<Read &>(*pcrI).addTagO(tmpmt);
		}
		catch (Notify n) {
		  cout << "Tried to transfer tags to bbackbone contig read from:\n";
		  Read::setCoutType(Read::AS_TEXT);
		  cout << *pcrI;
		  cout << "Exiting.\n";
		  n.handleError(THISFUNC);
		}
	      }
	    }
	  }

	  Read & rpr=AS_readpool.getRead(pcrI.getORPID());

	  // TODO: FIXME
	  // The whole thing is botched for reads without adjustment
	  // One should perform a Smith-Waterman of contig read against readpool read
	  //  and calculate new adjustments
	  // Not sure whether worth the effort atm
	  //
	  // Simple foolguard atm: transfer only if lengths are equal

	  // TODO: bad fix
	  // in assemblies with higher coverage and different variant, that is
	  //  essentially always untrue, hence no tag transposition at all. *sigh*
	  // if(rpr.getLenClippedSeq()==pcrI->getLenClippedSeq()) {

	  if(true) {
	    // this is wrong for reads without adjustments.
	    //  -> no way to detect deletions in those reads!
	    int32 adj_tagposl=pcrI->getLowerNonGapAdjustmentPosOfReadPos(acttag.from);
	    int32 adj_tagposu=pcrI->getUpperNonGapAdjustmentPosOfReadPos(acttag.to);
	    int32 rprtagposl=rpr.getReadPosOfAdjustmentPos(adj_tagposl);
	    int32 rprtagposu=rpr.getReadPosOfAdjustmentPos(adj_tagposu);

	    try {
	      //cout << "Transfering tag for " << rpr.getName()
	      //	 << "\t" << rprtagposl << " " << rprtagposu
	      //	 << "\tadj " << adj_tagposl << " " << adj_tagposu << endl;
	      //CEBUGF("Transfering tag for " << rpr.getName() << "\t" << rprtagposl << " " << rprtagposu << endl);
	      //Read::setCoutType(Read::AS_TEXT);
	      //CEBUGF("Before:\n" << rpr << endl);

	      // add tag only when both positions are >=0 (i.e. not starting/stopping
	      //  on an inserted base, else addTag() would understandably barf
	      if(rprtagposl>=0 && rprtagposu>=0){
		tmpmt.from=rprtagposl;
		tmpmt.to=rprtagposu;
		rpr.addTagO(tmpmt);
	      }
	      //CEBUGF("After:\n" << rpr << endl);
	    }
	    catch (Notify n) {
	      // care about errors only if we have adjustments in the read
	      //  if not, the whole thing is more or less wrong anyway
	      if(rpr.usesAdjustments()){
		cout << "Tried to transfer tags to readpool read from:\n";
		Read::setCoutType(Read::AS_TEXT);
		cout << *pcrI;
		cout << "Exiting." << endl;
		n.handleError(THISFUNC);
	      }else{
		cout << "Dbg: wrong tag transfer for " << rpr.getName() << '\n';
	      }
	    }
	  }
	}
      }
    }
  }
  catch(Notify n) {
    cout << "General error while transfering tag " << tagnumber << " to readpool read from:\n";
    Read::setCoutType(Read::AS_TEXT);
    cout << *pcrI;
    cout << "Exiting.\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}

//#define CEBUG(bla)



/*************************************************************************
 *
 * New version to search for spoil sports
 * Works with one contig at a time
 *
 * Look through all reads. If end of read with (500) bases of end of contig
 *  and has permbans and no freq >4:
 *     throw out (later version perhaps cut back)
 *
 *
 *************************************************************************/
// idea: keep info whether read has been totally embedded (>1000 bases left, right)
// if true, then less likely to be an error

void Assembly::huntSpoilSports(Contig & chkcon)
{
  cout << "Hunting join spoiler" << endl;

  auto & cr=chkcon.getContigReads();

  const uint32 endrange=30;

  auto pcrI=cr.begin();
  for(; pcrI != cr.end() ; ++pcrI){
    const Read & actread=*pcrI;
    bool atfront=false;
    bool atback=false;
    if(pcrI.getReadStartOffset()<=endrange) atfront=true;
    if(pcrI.getReadStartOffset()+actread.getLenClippedSeq() > chkcon.getContigLength()-endrange) atback=true;
    if(atfront || atback){
      cout << "HSS At end: " << atfront << ' ' << atback << '\t' << actread.getName() << endl;

      if(pcrI.getORPID()>=0){
	// currently: always remove

	// if it has >2 permbans, remove
	//cout << "HSS Permbans: " << AS_permanent_overlap_bans[pcrI.getORPID()].size() << endl;
	//if(AS_permanent_overlap_bans[pcrI.getORPID()].size()>3) {
	//  if(!actread.hasTag(Read::REA_tagentry_idSRMr)
	//     && !actread.hasTag(Read::REA_tagentry_idCRMr)){
	//    //if(true){
	//    cout << "HSS remove " << actread.getName() << endl;
	//
	//    const std::vector<Read::bposhashstat_t> & bposhashstats=actread.getBPosHashStats();
	//    if(!bposhashstats.empty()){
	//      bool clipfront=false;
	//      bool clipback=false;
	//      if(atfront){
	//	if(pcrI->direction>0){
	//	  clipfront=true;
	//	}else{
	//	  clipback=true;
	//	}
	//      }
	//      if(atback){
	//	if(pcrI->direction>0){
	//	  clipback=true;
	//	}else{
	//	  clipfront=true;
	//	}
	//      }
	//      cout << "HSS clip front back " << clipfront << ' ' << clipback << endl;
	//      if(clipfront){
	//	uint32 bposfrom=actread.getLeftClipoff();
	//	uint32 maxcheck=std::max(actread.getLenSeq()/4,static_cast<uint32>(50));
	//	uint32 maxto=std::min(maxcheck,actread.getLenSeq()-1);
	//	bool foundinv=false;
	//	for(; bposfrom<maxto; bposfrom++){
	//	  if(!bposhashstats[bposfrom].fwd.isValid()){
	//	    foundinv=true;
	//	    break;
	//	  }else{
	//	    if(bposhashstats[bposfrom].fwd.getFrequency()>=4) break;
	//	  }
	//	}
	//	if(foundinv){
	//	  if(bposfrom-actread.getLeftClipoff()<50){
	//	    bposfrom=actread.getLeftClipoff()+50;
	//	    if(bposfrom>=actread.getLenSeq()) bposfrom=actread.getLenSeq()-1;
	//	  }
	//	  cout << "HSS moving left " <<  AS_readpool[pcrI.getORPID()].getLQClipoff();
	//	  AS_readpool[pcrI.getORPID()].setLQClipoff(bposfrom);
	//	  cout << " to " << AS_readpool[pcrI.getORPID()].getLQClipoff()
	//	       << "\t" << actread.getName() << endl;
	//	}
	//      }
	//      if(clipback){
	//	uint32 bposto=actread.getRightClipoff();
	//	uint32 maxcheck=std::max(actread.getLenSeq()/4,static_cast<uint32>(50));
	//	bool foundinv=false;
	//	for(; bposto>0 && maxcheck>0; --bposto, --maxcheck){
	//	  if(!bposhashstats[bposto].rev.isValid()){
	//	    foundinv=true;
	//	    break;
	//	  }
	//	  if(bposhashstats[bposto].rev.getFrequency()>=4) break;
	//	}
	//	if(foundinv){
	//	  if(actread.getRightClipoff()-bposto<50){
	//	    if(actread.getRightClipoff()<50){
	//	      bposto=0;
	//	    }else{
	//	      bposto=actread.getRightClipoff()-50;
	//	    }
	//	  }
	//	  cout << "HSS moving right " <<  AS_readpool[pcrI.getORPID()].getRQClipoff();
	//	  AS_readpool[pcrI.getORPID()].setRQClipoff(bposto);
	//	  cout << " to " << AS_readpool[pcrI.getORPID()].getRQClipoff()
	//	       << "\t" << actread.getName() << endl;
	//	}
	//      }
	//    }
	//  }
	//}
	AS_istroublemaker[pcrI.getORPID()]=true;
      }
    }
  }
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::findPossibleOverlaps(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)
{
  FUNCSTART("void Assembly::findPossibleOverlaps(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)");

  std::string signalfilename(buildFileName(version,
				  prefix,
				  postfix,
				  AS_miraparams[0].getAssemblyParams().as_tmpf_signal_findpossibleoverlaps,
				  ".ok"));

  if(!AS_resumeasembly || !AS_resumeisok || !fileExists(signalfilename)){
    cout << "AS_resumeasembly " << AS_resumeasembly << endl;
    cout << "AS_resumeisok " << AS_resumeisok << endl;
    cout << "fileExists(" << signalfilename << ") " << fileExists(signalfilename) << endl;

    AS_resumeisok=false;
    fpo_main(version,prefix,postfix,tmpfname);

    // initialiase well connected with overlap criterion levels
    // log overlap criterion levels if wanted

    AS_wellconnected.clear();
    AS_wellconnected.resize(AS_readpool.size(),false);

    {
      std::ofstream fout;
      //AS_logflag_oclevel=true;
      if(AS_logflag_oclevel){
	std::string filename=buildFileName(version, "", "",
				      "elog.oclevel_pass",
				      ".lst");
	fout.open(filename, std::ios::out);
      }

      for(uint32 i=0; i<AS_readpool.size();++i){
	if(AS_logflag_oclevel){
	  fout << AS_readpool[i].getName()
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvl[0][i])
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvr[0][i])
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvl[1][i])
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvr[1][i])
	       << '\n';
	}
	// for Solexa, critlevel 0 is too harsch (because of special rule to calc
	//  Solexa critlevels in Skim)
	for(uint32 ocvi=0; ocvi<2;++ocvi){
	  if(AS_readpool[i].isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	    if(AS_overlapcritlevelvl[ocvi][i] <= 2 && AS_overlapcritlevelvr[ocvi][i] <= 2 ){
	      AS_wellconnected[i]=true;
	    }
	  }else if(AS_overlapcritlevelvl[ocvi][i] == 0 && AS_overlapcritlevelvr[ocvi][i] == 0 ){
	    AS_wellconnected[i]=true;
	  }
	}
      }
    }
    // save some needed data
    saveResumeDataFPO(version,prefix,postfix);

    reduceSkimHits4(version, prefix, postfix, tmpfname);

    std::ofstream fout(signalfilename);  // create checkpoint signal file for findPossibleOverlaps
  }else{
    cout << "Resume assembly: skim and skim reduction already present, good.\n";
    // set AS_pos?match_filename*
    fpo_buildFileNames(version,prefix,postfix,tmpfname);
    AS_posfmatch_full_filename=AS_posfmatch_filename+".reduced";
    AS_poscmatch_full_filename=AS_poscmatch_filename+".reduced";

    loadResumeDataFPO(version,prefix,postfix);
  }
}

void Assembly::fpo_buildFileNames(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)
{
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  if(tmpfname.size()){
    AS_posfmatch_filename=buildFileName(version,
					prefix,
					postfix,
					tmpfname+"f",
					".bin");
    AS_poscmatch_filename=buildFileName(version,
					prefix,
					postfix,
					tmpfname+"c",
					".bin");
    // the following two only to also automatically remove
    //  older versions of the reduced skim files
    buildFileName(version,
		  prefix,
		  postfix,
		  tmpfname+"f",
		  ".bin.reduced");
    buildFileName(version,
		  prefix,
		  postfix,
		  tmpfname+"c",
		  ".bin.reduced");
  }else{
    AS_posfmatch_filename=buildFileName(version,
					prefix,
					postfix,
					as_fixparams.as_tmpf_posmatch+"f",
					".bin");
    AS_poscmatch_filename=buildFileName(version,
					prefix,
					postfix,
					as_fixparams.as_tmpf_posmatch+"c",
					".bin");
    // the following two only to also automatically remove
    //  older versions of the reduced skim files
    buildFileName(version,
		  prefix,
		  postfix,
		  as_fixparams.as_tmpf_posmatch+"f",
		  ".bin.reduced");
    buildFileName(version,
		  prefix,
		  postfix,
		  as_fixparams.as_tmpf_posmatch+"c",
		  ".bin.reduced");
  }
}

void Assembly::fpo_main(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)
{
  FUNCSTART("void Assembly::fpo_main(int32 version, const std::string prefix, const std::string postfix, const std::string tmpfname)");

  //directory_parameters const & dir_params= AS_miraparams->getDirectoryParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();

  if(as_fixparams.as_dateoutput) dateStamp(cout);

  cout << "\n\nSearching for possible overlaps";


#if TRACKMEMUSAGE
    cout << "\ndmi fpo 00\n";
    dumpMemInfo();
#endif

  // save memory for this step, those structures will have to be recomputed anyway
  //nukeSTLContainer(AS_adsfacts);
  //nukeSTLContainer(AS_confirmed_edges);
  AS_adsfacts.clear();
  AS_confirmed_edges.clear();

  nukeSTLContainer(AS_readhitmiss);
  nukeSTLContainer(AS_readhmcovered);
  nukeSTLContainer(AS_count_rhm);

#if TRACKMEMUSAGE
    cout << "\ndmi fpo 05\n";
    dumpMemInfo();
#endif

  {
    std::vector<uint32> overlapcounter(AS_readpool.size(),0);
    bool onlyagainstrails=false;

    // very first call will be with version=0 ... pre_assembly
    //  don't set to true there as this analysis is then needed
    //  for multicopy analysis
    if(AS_hasbackbones
       && version >= as_fixparams.as_startbackboneusage_inpass){
      onlyagainstrails=true;
    }

    if(onlyagainstrails) cout << " (only against backbone, the progress bar will be skewed)";
    cout << ":\n";

    //std::string rawhitsfilename;
    std::string megahubtmpfname;
    {
      fpo_buildFileNames(version,prefix,postfix,tmpfname);
      AS_posfmatch_full_filename=AS_posfmatch_filename;
      AS_poscmatch_full_filename=AS_poscmatch_filename;

      if(tmpfname.size()){
	megahubtmpfname=buildFileName(version,
				      prefix,
				      postfix,
				      tmpfname+"_megahubs",
				      ".lst");
      }else{
	megahubtmpfname=buildFileName(version,
				      prefix,
				      postfix,
				      as_fixparams.as_tmpf_posmatch+"_megahubs",
				      ".lst");
      }

      //cout << "Only against rails? " << onlyagainstrails << endl;

      if(!onlyagainstrails) {
	if(skim_params.sk_basesperhash <=12
	   || (skim_params.sk_basesperhash <=14
	       && skim_params.sk_hashsavestepping <3)){
	  cout << "\n\nWARNING!!!!!!\nYou are not performing a 'mapping only' assembly and the parameters"
	    "\n -SK:bph=" << skim_params.sk_basesperhash << " and -SK:hss="
	       << static_cast<uint16>(skim_params.sk_hashsavestepping)
	       << "\nare quite low. If SKIM takes ages, stop this assembly and restart while"
	    "\nincreasing these parameters.\n\n";
	}
      }

      std::vector<int32> overlaplenrequired;
      std::vector<int32> prrequired;
      for(uint32 i=0;i<ReadGroupLib::getNumSequencingTypes(); i++){
	overlaplenrequired.push_back(AS_miraparams[i].getAlignParams().al_min_overlap);
	prrequired.push_back(AS_miraparams[i].getSkimParams().sk_percentrequired);
      }

      std::vector<int32> chuntleftcut;
      std::vector<int32> chuntrightcut;

      // trigger chimera search in pre-assembly skim;
      if(!AS_doneskimchimera &&
	 (as_fixparams.as_clip_skimchimeradetection || as_fixparams.as_clip_skimjunkdetection)){
	chuntleftcut.resize(1);
	chuntrightcut.resize(1);
	AS_doneskimchimera=true;
      }

      AS_chimeracutflag.clear();

      std::vector<uint8> * megahubsptr=nullptr;
      if(skim_params.sk_filtermegahubs){
	megahubsptr=&AS_skimmegahubs;
      }

      // force taking strong good overlaps only for denovo genome
      // TODO: currently not active in skim.
      bool forcetakestronggood=!as_fixparams.as_assemblyjob_mapping & AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms;


      uint32 nummegahubs=0;
      if(skim_params.sk_basesperhash <= 32){
	Skim<vhash64_t> s2;
	s2.setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);
	nummegahubs=s2.skimGo(AS_readpool,
			      AS_posfmatch_filename,
			      AS_poscmatch_filename,
			      AS_permanent_overlap_bans,
			      overlapcounter,
			      AS_writtenskimhitsperid,
			      chuntleftcut,
			      chuntrightcut,
			      AS_overlapcritlevelvl,
			      AS_overlapcritlevelvr,
			      megahubsptr,
			      skim_params.sk_numthreads,
			      skim_params.sk_maxhashesinmem,
			      onlyagainstrails,
			      skim_params.sk_alsoskimrevcomp,
			      skim_params.sk_basesperhash,
			      skim_params.sk_hashsavestepping,
			      prrequired,
			      overlaplenrequired,
			      skim_params.sk_maxhitsperread,
			      skim_params.sk_megahubcap,
			      forcetakestronggood
	  );
      }else if(skim_params.sk_basesperhash <= 64){
	Skim<vhash128_t> s2;
	s2.setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);
	nummegahubs=s2.skimGo(AS_readpool,
			      AS_posfmatch_filename,
			      AS_poscmatch_filename,
			      AS_permanent_overlap_bans,
			      overlapcounter,
			      AS_writtenskimhitsperid,
			      chuntleftcut,
			      chuntrightcut,
			      AS_overlapcritlevelvl,
			      AS_overlapcritlevelvr,
			      megahubsptr,
			      skim_params.sk_numthreads,
			      skim_params.sk_maxhashesinmem,
			      onlyagainstrails,
			      skim_params.sk_alsoskimrevcomp,
			      skim_params.sk_basesperhash,
			      skim_params.sk_hashsavestepping,
			      prrequired,
			      overlaplenrequired,
			      skim_params.sk_maxhitsperread,
			      skim_params.sk_megahubcap,
			      forcetakestronggood
	  );
      }else if(skim_params.sk_basesperhash <= 128){
	Skim<vhash256_t> s2;
	s2.setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);
	nummegahubs=s2.skimGo(AS_readpool,
			      AS_posfmatch_filename,
			      AS_poscmatch_filename,
			      AS_permanent_overlap_bans,
			      overlapcounter,
			      AS_writtenskimhitsperid,
			      chuntleftcut,
			      chuntrightcut,
			      AS_overlapcritlevelvl,
			      AS_overlapcritlevelvr,
			      megahubsptr,
			      skim_params.sk_numthreads,
			      skim_params.sk_maxhashesinmem,
			      onlyagainstrails,
			      skim_params.sk_alsoskimrevcomp,
			      skim_params.sk_basesperhash,
			      skim_params.sk_hashsavestepping,
			      prrequired,
			      overlaplenrequired,
			      skim_params.sk_maxhitsperread,
			      skim_params.sk_megahubcap,
			      forcetakestronggood
	  );
      }else if(skim_params.sk_basesperhash <= 256){
	Skim<vhash512_t> s2;
	s2.setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);
	nummegahubs=s2.skimGo(AS_readpool,
			      AS_posfmatch_filename,
			      AS_poscmatch_filename,
			      AS_permanent_overlap_bans,
			      overlapcounter,
			      AS_writtenskimhitsperid,
			      chuntleftcut,
			      chuntrightcut,
			      AS_overlapcritlevelvl,
			      AS_overlapcritlevelvr,
			      megahubsptr,
			      skim_params.sk_numthreads,
			      skim_params.sk_maxhashesinmem,
			      onlyagainstrails,
			      skim_params.sk_alsoskimrevcomp,
			      skim_params.sk_basesperhash,
			      skim_params.sk_hashsavestepping,
			      prrequired,
			      overlaplenrequired,
			      skim_params.sk_maxhitsperread,
			      skim_params.sk_megahubcap,
			      forcetakestronggood
	  );
      }else{
	MIRANOTIFY(Notify::FATAL,"Cannot perform a skim analysis with -SK:bph > 256 (you used " << skim_params.sk_basesperhash << "), though MIRA should've failed earlier, I admit.\n");
      }

      cout << "Total megahubs: " << nummegahubs << endl;

#if TRACKMEMUSAGE
    cout << "\ndmi fpo 10\n";
    dumpMemInfo();
#endif

      if(nummegahubs){

	{
	  std::ofstream mout(megahubtmpfname, std::ios::out| std::ios::trunc);
	  for(uint32 mhi=0; mhi<AS_skimmegahubs.size(); mhi++){
	    if(AS_skimmegahubs[mhi]>0) {
	      mout << AS_readpool[mhi].getName() << '\n';
	    }
	  }
	}

	cout << "\n\nMIRA has detected megahubs in your data."
	  "This may not be a problem, but most probably is, especially for eukaryotes.\n\n";
	if(100.0/AS_num_reads_valid*nummegahubs > skim_params.sk_maxmegahubratio){
	  cout << "\n\nYou have more than " << skim_params.sk_maxmegahubratio << "% of your reads found to be megahubs."
	    "\n\nYou should check the following:\n\n"
	    "\t1) for Sanger sequences: are all the sequencing vectors masked / clipped?\n"
	    "\t2) for 454 sequences: are all the adaptors masked / clipped?\n\n";
	  if(AS_miraparams[0].getHashStatisticsParams().hs_repeatlevel_in_infofile){
	    cout << "You will find in the info directory a file called\n"
	      "    '*_info_readrepeats.lst',\n"
	      "consult the MIRA manual on how to extract repeat information from there.\n\n";
	  }else{
	    cout << "To learn more on the types of repeats you have and how MIRA\n"
	      " can help you to find them, please consult the manual on the\n"
	      " usage of -HS:rliif and the tmp files they create.\n";
	  }
	  cout << "*ONLY* when you are sure that no (or only a very negligible number) of sequencing"
	    "\nvector / adaptor sequence is remaining, try this:\n\n"
	    "\t3) for organisms with complex repeats (eukaryots & some bacteria):\n";
	  if(!AS_miraparams[0].getHashStatisticsParams().hs_masknastyrepeats) cout << "\t\t- use -HS:mnr=yes\n";
	  cout << "\t\t- reduce the -HS:nrr parameter (divide by 2)\n"
	    // huh?
	    //"4) for EST projects, -SK:nrr will not really work, use -SK:nrr (start at\n"
	    //"   10 and increase in steps of of 10)\n"
	    "\n"
	    "*ONLY* if the above fails, try increasing the -SK:mmhr parameter\n"
	    "Note that the number of present megahubs will increase computation time in\n"
	    "an exponential way, so be careful when changing -SK:mmhr.\n";
	}

	if(100.0/AS_readpool.size()*nummegahubs >= skim_params.sk_maxmegahubratio){
	  cout << "\n\nYou have " << 100.0/AS_readpool.size()*nummegahubs
	       << "% of your reads as megahubs.\n"
	       << "You have set a maximum allowed ratio of: " << skim_params.sk_maxmegahubratio
	       << "\nTo change this ratio, read the text above (it's kind of important) and also see -SK:mmhr at multiple location in the docs\n"
	       << "\n\nEnding the assembly because the maximum ratio has been reached/surpassed.\n";
	  exit(10);
	}
      }

      if(chuntleftcut.size()){
	std::string cltmpfname;
	cltmpfname=buildFileName(0,"","",
				as_fixparams.as_tmpf_clippings,
				 ".txt",
				 "",
				 false);
	std::string logprefix="skim detect: ";
	AS_chimeracutflag.resize(1);
	cutBackPossibleChimeras(cltmpfname, logprefix, chuntleftcut,chuntrightcut,AS_chimeracutflag);
      }

    }

    // in mapping assemblies, correct the matches not being 100%
    if(AS_miraparams[0].getSkimParams().sk_swcheckonbackbones
       && AS_hasbackbones){
      recalcNonPerfectSkimMappingsBySW(version);
    }


    //// log the raw hash hits
    //{
    //  std::ofstream ofs;
    //  ofs.open(rawhitsfilename, std::ios::out| std::ios::trunc);
    //  for(size_t rhhci=0; rhhci < rawhashhitcounter.size(); rhhci++){
    //	ofs << rhhci << '\t' << rawhashhitcounter[rhhci] << '\n';
    //  }
    //  ofs.close();
    //}

    // Do this only once
    // TODO: check if ok to do more, i.e. if skim can be adapted to
    //       still count banned overlaps.
    if(AS_multicopies.size()==0){
      std::string filename;
      if(tmpfname.size()){
	filename=buildFileName(version, prefix, postfix, tmpfname+"_multicopystat", ".txt");
      }else{
	filename=buildFileName(version, prefix, postfix,
			       AS_miraparams[0].getAssemblyParams().as_tmpf_posmatch+"_multicopystat",
			       ".txt");
      }
      std::ofstream fout(filename, std::ios::out);
      fout.close();
      //flagMulticopyReads(overlapcounter, filename);
    }

    //if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    //cout << '\n';
    //exit(0);
  }


// TODO: adapt to data in files instead of the old posXmatch multimaps
/*
#if 0
  {
    std::vector<int32> cluster;
    cluster.resize(AS_readpool.size(),-1);
    uint32 clustercount=0;

    auto I=AS_posfmatch.cbegin();
    while(I!=AS_posfmatch.end()){
      int32 cnum1=cluster[I->first];
      int32 cnum2=cluster[I->second.otherid];
      if(cnum1==-1 && cnum2==-1) {
	cluster[I->first]=clustercount;
	cluster[I->second.otherid]=clustercount;
	clustercount++;
      } else if(cnum1==-1) {
	cluster[I->first]=cluster[I->second.otherid];
      } else if(cnum2==-1) {
	cluster[I->second.otherid]=cluster[I->first];
      } else {
	if (cnum1!=cnum2) {
	  // uh oh ... we have to merge both these clusters
	  // simply change all cnum1 into cnum2 in cluster vector
	  for(uint32 j=0; j<AS_readpool.size(); j++) {
	    if(cluster[j]==cnum1) cluster[j]=cnum2;
	  }
	}
      }

      I++;
    }

    I=AS_poscmatch.begin();
    while(I!=AS_poscmatch.end()){
      int32 cnum1=cluster[I->first];
      int32 cnum2=cluster[I->second.otherid];
      if(cnum1==-1 && cnum2==-1) {
	cluster[I->first]=clustercount;
	cluster[I->second.otherid]=clustercount;
	clustercount++;
      } else if(cnum1==-1) {
	cluster[I->first]=cluster[I->second.otherid];
      } else if(cnum2==-1) {
	cluster[I->second.otherid]=cluster[I->first];
      } else {
	if (cnum1!=cnum2) {
	  // uh oh ... we have to merge both these clusters
	  // simply change all cnum1 into cnum2 in cluster vector
	  for(uint32 j=0; j<AS_readpool.size(); j++) {
	    if(cluster[j]==cnum1) cluster[j]=cnum2;
	  }
	}
      }

      I++;
    }

    std::string filename;
    if(tmpfname.size()){
      filename=buildFileName(version, prefix, postfix, tmpfname+"_pcluster", ".lst");
    }else{
      filename=buildFileName(version, prefix, postfix,
			     AS_miraparams->getAssemblyParams().as_tmpf_posmatch+"_pcluster",
			     ".lst");
    }

    std::ofstream fout;
    fout.open(filename, std::ios::out);
    uint32 outputcount=0;
    for(uint32 i=0; i<clustercount; i++) {
      bool found=false;
      for(uint32 j=0; j<AS_readpool.size(); j++) {
	if(cluster[j]==static_cast<int32>(i)) {
	  found=true;
	  fout << outputcount << " " << AS_readpool.getRead(j).getName() << endl;
	}
      }
      if(found) outputcount++;
    }
    for(uint32 j=0; j<AS_readpool.size(); j++) {
      if(cluster[j]==-1) {
	fout << "-1 " << AS_readpool.getRead(j).getName() << endl;
      }
    }
  }
#endif
*/

#if TRACKMEMUSAGE
    cout << "\ndmi fpo 20\n";
    dumpMemInfo();
#endif

  AS_steps[ASADSLISTOK]=0;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

/*
void Assembly::flagMulticopyReads(const std::vector<uint32> & overlapcounter, const std::string & tmpfilename)
{
  AS_multicopies.clear();
  AS_multicopies.resize(AS_readpool.size(),0);

  std::ofstream fout;
  fout.open(tmpfilename, std::ios::out);

  //cout << "Searching for multicopy reads:\n";

  for(uint8 actseqtype=0; actseqtype<ReadGroupLib::SEQTYPE_END; actseqtype++){
    std::vector<uint32> sortedoverlapcounter=overlapcounter;

    // set overlapcounter of backbones and railreads and read that are
    //  not of currently analysed sequencing type to 0 so as not to
    //  count them
    uint32 numreadsinactseqtype=0;
    for(uint32 i=0; i<overlapcounter.size(); i++){
      if(AS_readpool[i].isRail()
	 || AS_readpool[i].isBackbone()
	 || AS_readpool[i].getSequencingType()!=actseqtype) {
	sortedoverlapcounter[i]=0;
      }else{
	numreadsinactseqtype++;
      }
    }

    if(numreadsinactseqtype==0) continue;

    //cout << numreadsinactseqtype << " in sequencing type " << static_cast<uint16>(actseqtype) << endl;
    sort(sortedoverlapcounter.begin(), sortedoverlapcounter.end());

    // 5% quantil
    uint32 quantilnumber=5*numreadsinactseqtype/100;
    uint32 ifrom=0;
    // well, start the 5% quantil only at reads which have at least
    //  1 overlap
    while(ifrom<sortedoverlapcounter.size()
	  && sortedoverlapcounter[ifrom]==0) ifrom++;

    ifrom+=quantilnumber;
    uint32 ito=static_cast<uint32>(sortedoverlapcounter.size())-quantilnumber;

    if(ito<=ifrom || ifrom>=sortedoverlapcounter.size()){
      ifrom=0;
      ito=static_cast<uint32>(sortedoverlapcounter.size());
    }

    uint32 nonsinglets=0;
    uint32 totaloverlaps=0;
    for(uint32 i=ifrom; i<ito; i++){
      totaloverlaps+=sortedoverlapcounter[i];
      nonsinglets++;
    }


    if(nonsinglets==0) nonsinglets=1;
    uint32 avgoverlaps=static_cast<uint32>(.5+static_cast<double>(totaloverlaps)/static_cast<double>(nonsinglets));
    // strictly speaking, this is not median. But close enough.
    uint32 medianoverlaps=sortedoverlapcounter[ifrom+((ito-ifrom)/2)];
    //uint32 multicopythreshold=avgoverlaps*2;
    uint32 multicopythreshold=medianoverlaps*2;

    fout << "Hitstatistics (" << ReadGroupLib::getShortNameOfSequencingType(actseqtype) << "): nonsinglets:" << nonsinglets << "\ttotaloverlaps: " << totaloverlaps << "\tavgoverlaps: " << avgoverlaps << "\tmedianoverlaps: " << medianoverlaps << endl;

    // now set multicopy flags for affected reads of this sequencing type
    for(uint32 i=0; i<overlapcounter.size(); i++){
      if(AS_readpool[i].isRail()
	 || AS_readpool[i].isBackbone()
	 || AS_readpool[i].getSequencingType()!=actseqtype) continue;
      if(overlapcounter[i]>multicopythreshold) {
	AS_multicopies[i]=static_cast<uint8>(actseqtype+1);
	if(overlapcounter[i]>multicopythreshold*10) {
	  AS_multicopies[i]=static_cast<uint8>(actseqtype+1+100);
	}
      }
    }
  }

  for(uint32 i=0; i<overlapcounter.size(); i++){
    if(AS_readpool[i].isRail()
       || AS_readpool[i].isBackbone()) continue;
    fout << i << "\t" << AS_readpool[i].getName() << "\t" << overlapcounter[i];
    if(AS_multicopies[i]>100) {
      fout << "\tst: " << static_cast<uint16>(AS_multicopies[i]-101);
      fout << "\tmulticopy / insane (forgotten clone vector?)";
    }else if(AS_multicopies[i]>0) {
      fout << "\tst: " << static_cast<uint16>(AS_multicopies[i]-1);
      fout << "\tmulticopy";
    }
    fout << '\n';
  }
}
*/





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const Read & Assembly::getRead(uint32 index)
{
  FUNCSTART("const Read & Assembly::getRead(uint32 index)");

  if(index>=AS_readpool.size()){
    MIRANOTIFY(Notify::INTERNAL,"index: " << index << " greater than AS_readpool.size():" << AS_readpool.size() << "  (out of bounds)");
  }

  FUNCEND();

  return AS_readpool.getRead(static_cast<int32>(index));
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::refreshContigAndReadpoolValuesAfterLoading(ReadPool & rp, std::list<Contig> & contigs, uint8 nagwarn_templateproblems)
{
  FUNCSTART("void Assembly::refreshContigAndReadpoolValuesAfterLoading(ReadPool & rp, std::list<Contig> & contigs, uint8 nagwarn_templateproblems)");

  try{
    rp.makeTemplateIDs(nagwarn_templateproblems);

    // Not needed anymore for readgroup approach: RG keeps track of this
    //rp.makeStrainIDs();

    for(auto & cle : contigs) cle.recalcTemplateIDsAndStrainPresent();
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
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


/*************************************************************************
 *
 * dead code, not used anymore?
 *
 *
 *************************************************************************/

//void Assembly::banReadPairGroups(const std::vector<int32> & g1, const std::vector<int32> & g2)
//{
//  FUNCSTART("void Assembly::banReadPairGroups(const std::vector<int32> & g1, const std::vector<int32> & g2)");
//
//  int32 id1, id2;
//  newedges_t tmpegde;
//  for(uint32 i1=0; i1 < g1.size(); i1++) {
//    id1=g1[i1];
//    for(uint32 i2=0; i2 < g2.size(); i2++) {
//      id2=g2[i2];
//
//      // skip it if there's alreaddy a permban on it
//      if(AS_permanent_overlap_bans.checkIfBanned(id1,id2)) continue;
//
//      if(AS_readpool.getRead(id1).getTemplatePartnerID()==id2) {
//#ifndef PUBLICQUIET
//	cout << "Do not ban template partners "  << id1 << ": " << AS_readpool.getRead(id1).getName() << "\t" << id2 << ": " << AS_readpool.getRead(id2).getName() << endl;
//#endif
//	  continue;
//      }
//
//#ifndef PUBLICQUIET
//      cout << "Banning: " << id1 << ": " << AS_readpool.getRead(id1).getName() << "\t" << id2 << ": " << AS_readpool.getRead(id2).getName() << endl;
//#endif
//
//      // put both ids in permanent overlap banlist
//      AS_permanent_overlap_bans.insertBan(id1,id2);
//
//      // now remove overlap edges between these reads
//      // first in one direction ...
//      //overlap_edges_t::iterator Irun=AS_edges_forward.lower_bound(id1);
//      tmpegde.rid1=id1;
//      std::vector<newedges_t>::iterator Irun=lower_bound(AS_confirmed_edges.begin(),
//						    AS_confirmed_edges.begin(),
//						    tmpegde,
//						    Assembly__compareNewEdges_t_);
//      while(Irun != AS_confirmed_edges.end()
//	    && Irun->rid1 == id1) {
//	if(Irun->linked_with==id2) {
//	  // erase() doesn't give back an iterator as advertised?
//	  AS_confirmed_edges.erase(Irun);
//	  //Irun=AS_edges_forward.lower_bound(id1);
//	  Irun=lower_bound(AS_confirmed_edges.begin(),
//			   AS_confirmed_edges.begin(),
//			   tmpegde,
//			   Assembly__compareNewEdges_t_);
//	  // Irun points to element after (or end) now, no need to increment
//	  continue;
//	}
//	Irun++;
//      }
//
//      // .. then in the other one
//      //Irun=AS_edges_forward.lower_bound(id2);
//      tmpegde.rid1=id2;
//      Irun=lower_bound(AS_confirmed_edges.begin(),
//		       AS_confirmed_edges.begin(),
//		       tmpegde,
//		       Assembly__compareNewEdges_t_);
//      while(Irun != AS_confirmed_edges.end()
//	    && Irun->rid1 == id2) {
//	if(Irun->linked_with==id1) {
//	  AS_confirmed_edges.erase(Irun);
//	  //Irun=AS_edges_forward.lower_bound(id2);
//	  Irun=lower_bound(AS_confirmed_edges.begin(),
//			   AS_confirmed_edges.begin(),
//			   tmpegde,
//			   Assembly__compareNewEdges_t_);
//	  // Irun points to element after (or end) now, no need to increment
//	  continue;
//	}
//	Irun++;
//      }
//    }
//  }
//
//  FUNCEND();
//}



/////////////////////////////////////////////////////////////////////////
////////////////////        Dead since a while         //////////////////
/////////////////////////////////////////////////////////////////////////
