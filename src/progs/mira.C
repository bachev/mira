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


#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <sys/resource.h>     // defines rlimit for getMoreStack()

#include <boost/algorithm/string.hpp>

#include "errorhandling/errorhandling.H"
#include "util/machineinfo.H"
#include "util/fileanddisk.H"
#include "progs/quirks.H"

#include "modules/mod_mira.H"
#include "modules/mod_bait.H"
#include "modules/mod_sqt.H"
#include "modules/mod_convert.H"
#include "modules/mod_diff.H"
#include "modules/mod_memestim.H"
#include "modules/mod_mer.H"
#include "modules/mod_scaffold.H"
#include "modules/mod_sigconex.H"
#include "modules/mod_tagsnp.H"

#include "version.H"

#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

using std::cout;
using std::cerr;
using std::endl;


extern const char * compileinfo;

void doAbort()
{
#ifndef PUBLICQUIET
  Read::dumpStringContainerStats(cout);
#endif

  cout << "\n\nVCODE: " << miraversion << endl;

  cout << "\n\nFor general help, you will probably get a quicker response on the\n"
    "    MIRA talk mailing list\n"
    "than if you mailed the author directly.\n"
    "\nTo report issues or ask for features, please use the GitHub issue\nsystem at:\n"
    "\thttps://github.com/bachev/mira/issues\n"
    "This ensures that requests do not get lost.\n";
  exit(100);
}



// usual linux stack size of 8Mb will lead to segfaults in very long
//  alignments (>15-18k) in align.C
// therefore, get some more stack, atm 64 Mb

void getMoreStack()
{
  struct rlimit rl;
  const rlim_t wantstacksize=64L*1024L*1024L;

  auto result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0){
    //cout << "Have cur " << rl.rlim_cur << endl;
    //cout << "Have max " << rl.rlim_max << endl;
    if(rl.rlim_cur<wantstacksize){
      rl.rlim_cur=wantstacksize;
      result = setrlimit(RLIMIT_STACK, &rl);
    }
  }else{
    cout << "could not query stack size?\n";
  }
}

int main(int argc, char ** argv)
{
#ifdef MIRAMEMORC
  MemORC::setChecking(true);
#endif

  fixQuirks();
  getMoreStack();

#ifdef TIMERESTRICTED

  cout << "Compiled on " << __DATE__ << ". Will run until "<< TR_OUT_MAXDAY << "." << TR_OUT_MAXMONTH << "." << TR_OUT_MAXYEAR << " (dd.mm.yyyy)\n";

  {
    struct stat st;
    int rc=stat(dir_params.dir_log.c_str(),&st);
    struct tm *mytm;
    mytm=localtime(&st.st_mtime);

    if(mytm->tm_year > TR_MAXYEAR
       || (mytm->tm_year == TR_MAXYEAR
	   && mytm->tm_mon >= TR_OUT_MAXMONTH)){
      cerr << "\n\nThis version of MIRA is definitively old, please get a newer version of the assembler.\n";
      exit(0);
    }
  }
#endif

  //cpu_set_t mask;
  //cout << "####" << sizeof(cpu_set_t) << endl;
  //sched_getaffinity(0,sizeof(cpu_set_t),&mask);
  //
  //cout << "Affinity: " << (uint8 *) &mask[0] << endl;

#ifdef HAVE_OPENMP
  omp_set_num_threads(MachineInfo::getCoresTotal());
#endif


  std::string path;
  std::string miraprog;

  splitFullPathAndFileName(argv[0],path,miraprog);
  boost::to_lower(miraprog);

  int retvalue=0;
  try{
    if(miraprog=="mira"
      || miraprog=="mira4"
      || miraprog=="mira5"){
      miraMain(argc,argv);
//    } else if(miraprog=="mirasearchestsnps"){
//      //miraEST(argc,argv);
//    } else if(miraprog=="miraclip"){
//      //miraClip(argc,argv);
//    } else if(miraprog=="mirapre"){
//      //miraPre(argc,argv);
//    } else if(miraprog=="tagsnp"){
//      tagsnp t;
//      t.mainTagSNP(argc, argv);
    } else if(miraprog=="mirapbcorrect"){
      mainPBCorrect(argc,argv);
    } else if(miraprog=="miramem"
	      ||miraprog=="mira4mem"){
      miraMemEstimate(argc, argv);
    } else if(miraprog=="dbgreplay"){
      dbgReplayMain(argc, argv);
    }else if(miraprog=="mirascaff"
	     || miraprog=="mira4scaff"){
      MiraScaffold m;
      m.mainMiraScaffold(argc, argv);
    }else if(miraprog=="mirabait"
	     || miraprog=="mira4bait"){
      MiraBait m;
      retvalue=m.mainMiraBait(argc, argv);
    }else if(miraprog=="mirasqt"
	     || miraprog=="mira4sqt"){
      MiraSQT m;
      retvalue=m.mainMiraSQT(argc, argv);
    }else if(miraprog=="miradiff"
	     || miraprog=="mira4diff"){
      MiraDiff m;
      m.mainMiraDiff(argc, argv);
    }else if(miraprog=="mirasce"
	     || miraprog=="mira4dsce"){
      MiraSCE m;
      m.mainSignificantContigExtraction(argc, argv);
    }else if(miraprog=="miramer"
	     || miraprog=="mira4mer"){
      MiraMer m;
      m.mainMiraMer(argc, argv);
    }else if(miraprog=="miraconvert"
	     || miraprog=="mira4convert"
	     || miraprog=="convert_project"
	     || miraprog=="convert_projectd"){
      if(miraprog.front()=='c'){
	cout << "convert_project is a deprecated name for miraconvert, please switch to new name.\n";
      }
      ConvPro cp;
      cp.mainConvPro(argc, argv);
    } else {
      cout << miraprog << " is a non-recognised program name of MIRA.\n"
	"The programs SHOULD be named either\n"
	"\"mira\", \"miraconvert\", \"miramem\", \"mirabait\""
	"\nAssuming 'mira'\n" << endl;

      miraMain(argc,argv);
    }
  }
  catch(Notify n){
    cout.flush(); cerr.flush();
    n.handleError("main");
  }
  catch(Flow f){
    cout.flush(); cerr.flush();
    cout << "INTERNAL ERROR: Unexpected exception: Flow()\n";
    doAbort();
  }
  catch(const std::bad_alloc & e){
    cout.flush(); cerr.flush();
    cout << "Out of memory detected, exception message is: ";
    cout << e.what() << endl;

    if(sizeof(size_t) == sizeof(int32)){
      cout << "\nYou are running a 32 bit executable. Please note that the maximum"
	"\ntheoretical memory a 32 bit programm can use (be it in Linux, Windows or"
	"\nother) is 4 GiB, in practice less: between 2.7 and 3.3 GiB. This is valid"
	"\neven if your machine has hundreds of GiB."
	"\nShould your machine have more that 4 GiB, use a 64 bit OS and a 64 bit"
	"\nversion of MIRA.";
    }

    cout << "\n\nIf you have questions on why this happened, please send the last 1000"
      "\nlines of the output log (or better: the complete file) to the author"
      "\ntogether with a short summary of your assembly project.\n\n";

    doAbort();
  }
  catch(const std::ios_base::failure & e){
    cout.flush(); cerr.flush();
    cout << "Failure in IO stream detected, exception message is: "
	 << e.what() << endl
	 << "\nWe perhaps ran out of disk space or hit a disk quota?\n";
    doAbort();
  }
  catch (std::exception& e) {
    cout.flush(); cerr.flush();
    cout << "A 'standard' exception occurred (that's NOT normal):\n" << e.what() << "\n\nIf the cause is not immediately obvious, please contact: bach@chevreux.org\n\n";
    doAbort();
  }
  catch(...){
    cout.flush(); cerr.flush();
    cout << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    doAbort();
  }

#ifndef PUBLICQUIET
  Read::dumpStringContainerStats(cout);
#endif

  return retvalue;
}
