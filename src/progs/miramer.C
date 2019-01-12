/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2013 and later by Bastien Chevreux
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


#include <getopt.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>


#include "io/generalio.H"
#include "util/fileanddisk.H"
#include "util/machineinfo.H"

#include "mira/readpool_io.H"
#include "mira/assembly.H"
#include "mira/hashstats.H"
#include "mira/manifest.H"
#include "mira/parameters.H"
#include "mira/seqtohash.H"

#include "version.H"



using namespace std;



// make the "tcmalloc: large alloc" messages from TCMallom disappear
// by setting the reporting environment variable to a very large value
// see: http://groups.google.com/group/google-perftools/browse_thread/thread/24a003fc35f3d470?pli=1
void quietenTCMALLOC()
{
  setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD",
	 "1099511627776", // 1 TiB, should be enough to quieten almost everything
	 0); // "0" means "do not overwrite if already existing"
}

/*
  On some systems, Boost::filesystem (at least up to 1.50) throws a standard
  exception when using some functions: locale::facet::_S_create_c_locale name not valid
  Only countermeasure possible: setenv
  Using setlocale or std::locale::global doe *NOT* work as workaround (tried and tested)
 */
void fixLocaleQuirk()
{
  try{
    // this must work
    boost::filesystem::path fp = boost::filesystem::current_path();
  }
  catch(...){
    // if not, we're on a system with quirks
    // maybe using LC_CTYPE???
    setenv("LC_ALL",
	   "C",
	   1);

    cout << "Your system seems to be older or have some quirks with locale settings."
      "\nUsing the LC_ALL=C workaround."
      "\nIf you don ot want that, upgrade or fix your system ;-)\n";
  }
}

void doAbort()
{
#ifndef PUBLICQUIET
  Read::dumpStringContainerStats(cout);
#endif

  cout << "\n\nFor general help, you will probably get a quicker response on the\n"
    "    MIRA talk mailing list\n"
    "than if you mailed the author directly.\n"
    "\nTo report bugs or ask for features, please use the new ticketing system at:\n"
    "\thttp://sourceforge.net/apps/trac/mira-assembler/\n"
    "This ensures that requests don't get lost.\n";
  abort();
}

namespace miramer {
  uint8 MER_basesperhash=31;
  uint16 MER_numlearnsteps=2;
  string MER_job("create");

  void main(int argc, char ** argv);
  void merCreateHashStats(int argc, char ** argv);
  void merSortHashStats(int argc, char ** argv);
  void merDumpHashStats(int argc, char ** argv);
  void merDumpHashDistrib(int argc, char ** argv);
};

void miramer::main(int argc, char ** argv)
{
  FUNCSTART("void miramer::main(int argc, char ** argv)");

  // that loop is straight from the GNU getopt_long example
  // http://www.gnu.org/s/hello/manual/libc/Getopt-Long-Option-Example.html
  while (1){
    static struct option mlong_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"job", required_argument,         0, 'j'},
	{"kmersize", required_argument,         0, 'k'},
	{"version", no_argument,         0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hj:k:v",
		     mlong_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 'h':
      cout << "mira\t\tMIRALIB version " << MIRAVERSION << "\n"
	"Author:\t\tBastien Chevreux (bach@chevreux.org)\n"
	"Purpose:\thandle k-mer statistics of a data set\n\n";

      cout << "Usage:\n"
	"miradiff ...\n";
      cout << "\nOptions:\n";
      cout <<
	"  -h / --help\t\t\t\tPrint short help and exit\n"
	"  -v / --version\t\t\tPrint version and exit\n"
	;
      exit(0);
    case 'j': {
      MER_job=optarg;
      boost::to_lower(MER_job);
      break;
    }
    case 'k': {
      uint64 bla=atoi(optarg);
      if(bla>32) bla=32;
      MER_basesperhash=bla;
      break;
    }
    case 'v':
      cout << MIRAVERSION << endl;
      exit(0);
    default:
      abort();
    }
  }

  if(MER_job=="create"){
    merCreateHashStats(argc,argv);
  }else if(MER_job=="sort"){
    merSortHashStats(argc,argv);
  }else if(MER_job=="dumpcounts"){
    merDumpHashStats(argc,argv);
  }else if(MER_job=="dumpdistrib"){
    merDumpHashDistrib(argc,argv);
  }else{
    cout << argv[0] << ": unknown job '" << MER_job << "'???" << endl;
    exit(1);
  }
}

void miramer::merCreateHashStats(int argc, char ** argv)
{
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  string loadfn(argv[optind++]);

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  auto rgid = ReadGroupLib::getReadGroupID(0);

  NHashStatistics nhs;
  nhs.setupNewAnalysis(32,4,MER_basesperhash,MER_numlearnsteps);
  {
    uint8 ziptype=0;
    string ft,pathto,stem;
    guessFileAndZipType(loadfn,pathto,stem,ft,ziptype);

    ReadPool rp1;
    {
      ReadPoolIO rpio(rp1);
      rpio.registerFile(
	"fastq",
	loadfn,
	"",
	rgid,
	false);
      rpio.loadNextSeqs(-1);
    }

    nhs.analyseReadPool(rp1);
    nhs.dumpHealth(cout);
    nhs.deleteBloomFilter();
    nhs.saveHashStatistics(stem+".mhs",true);
  }

}

void miramer::merSortHashStats(int argc, char ** argv)
{
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  string loadfn(argv[optind++]);
  NHashStatistics nhs;
  nhs.loadHashStatistics(loadfn);
  nhs.sortLow24Bit();
  nhs.saveHashStatistics(loadfn+".sorted",true);
}

void miramer::merDumpHashStats(int argc, char ** argv)
{
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  string loadfn(argv[optind++]);
  NHashStatistics nhs;
  nhs.loadHashStatistics(loadfn);
  nhs.sortLexicographically();
  nhs.dumpHashCount(cout);
}

void miramer::merDumpHashDistrib(int argc, char ** argv)
{
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  string loadfn(argv[optind++]);
  NHashStatistics nhs;
  nhs.loadHashStatistics(loadfn);
  nhs.dumpHashDistrib(cout);
}





namespace miradiff {
  void main(int argc, char ** argv);
  void prepHashStats(NHashStatistics & nhs);
};


void miradiff::prepHashStats(NHashStatistics & nhs)
{
  cout << "trim by 4 4 -1\n";
  nhs.trimHashStatsByFrequency(2,2,-1);
  cout << "hash distrib:\n";
  nhs.dumpHashDistrib(cout);
  // force nhs to calc offset table
  NHashStatistics::nhashstat_t tmp;
  nhs.findVHash(tmp);
}

void miradiff::main(int argc, char ** argv)
{
  FUNCSTART("void miradiff::main(int argc, char ** argv)");

  // that loop is straight from the GNU getopt_long example
  // http://www.gnu.org/s/hello/manual/libc/Getopt-Long-Option-Example.html
  while (1){
    static struct option long_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"version", no_argument,         0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hv",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 'h':
      cout << "mira\t\tMIRALIB version " << MIRAVERSION << "\n"
	"Author:\t\tBastien Chevreux (bach@chevreux.org)\n"
	"Purpose:\tdiff two data sets\n\n";

      cout << "Usage:\n"
	"miradiff ...\n";
      cout << "\nOptions:\n";
      cout <<
	"  -h / --help\t\t\t\tPrint short help and exit\n"
	"  -v / --version\t\t\tPrint version and exit\n"
	;
      exit(0);
    case 'v':
      cout << MIRAVERSION << endl;
      exit(0);
    default:
      abort();
    }
  }


  string refdata("pa1077_4mr");
  string querydata("pa1077_cilv2_4_4mr");

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  auto rgid = ReadGroupLib::newReadGroup();
  rgid.setSequencingType("solexa");
  rgid.fillInSensibleDefaults();

  NHashStatistics nhs1;
  nhs1.loadHashStatistics(refdata+".mhs");
  prepHashStats(nhs1);

  NHashStatistics nhs2;
  nhs2.loadHashStatistics(querydata+".mhs");
  prepHashStats(nhs2);

  {
    ReadPool rp;
    {
      ReadPoolIO rpio(rp);
      rpio.registerFile(
	"fastq",
	querydata+".fastq",
	"",
	rgid,
	false);
      rpio.loadNextSeqs(-1);
    }
    NHashStatistics::nhashstat_t tmphs;

    ofstream fout("bingo.txt");

    dateStamp(cout);

    vector<uint64> yyseqi;
    vector<uint64> ynseqi;

    ProgressIndicator<int64> pi(0,rp.size());
    for(uint32 rpi=0; rpi<rp.size(); ++rpi){
      pi.progress(rpi);

      bool bingo=false;
      uint32 yesyes=0;
      uint32 yesno=0;

      auto basesperhash=31;
      uint64 slen=rp[rpi].getLenClippedSeq();
      auto namestr=rp[rpi].getName().c_str();
      const uint8 * seq=reinterpret_cast<const uint8 *>(rp[rpi].getClippedSeqAsChar());
      yyseqi.clear();
      ynseqi.clear();
      SEQTOHASH_LOOPSTART(vhash_t){
	tmphs.vhash=acthash;
	auto hptr=nhs2.findVHash(tmphs);
	if(hptr!=nullptr){
	  hptr=nhs1.findVHash(tmphs);
	  if(hptr==nullptr){
	    ynseqi.push_back(seqi);
	    ++yesno;
	  }else{
	    yyseqi.push_back(seqi);
	    ++yesyes;
	  }
	}
      }SEQTOHASH_LOOPEND;

      if(yesno){
	fout << rp[rpi].getName()
	     << '\t' << yesyes
	     << '\t' << yesno
	     << "\t(" << slen;
	if(ynseqi.size()>=21) fout << ",bingo" << ynseqi.size();
	fout << "):";
	for(auto & si : ynseqi) fout << '\t' << si;
	fout << endl;
      }
    }
  }

  dateStamp(cout);



  exit(0);

  Manifest manifest;
//  MIRAParameters::generateProjectNames(Pv,manifest.getProjectName());
//
//  string mparams(manifest.getFullMIRAParameterString());
//  cout << "Seen parameters in manifest: " << mparams << endl;
//  MIRAParameters::parse(mparams, Pv);
//
//
//  MIRAParameters::postParsingChanges(Pv);
//  MIRAParameters::dumpAllParams(Pv, cout);

//...

  return;
}


int main(int argc, char ** argv)
{
  fixLocaleQuirk();
  quietenTCMALLOC();

  string path;
  string miraprog;

  splitFullPathAndFileName(argv[0],path,miraprog);
  boost::to_lower(miraprog);

  try{
    if(miraprog=="miramer"){
      miramer::main(argc,argv);
    } else if(miraprog=="miradiff"){
      miradiff::main(argc,argv);
    } else {
      cout << miraprog << " is a non-recognised program name of MIRAmer.\n"
	"The programs SHOULD be named either\n"
	"\"miramer\", or \"miradiff\""
	;
      exit(100);
    }
  }
  catch(Notify n){
    n.handleError("main");
  }
  catch(Flow f){
    cout << "INTERNAL ERROR: Unexpected exception: Flow()\n";
    doAbort();
  }
  catch(const std::bad_alloc & e){
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
  catch(const ios_base::failure & e){
    cout << "Failure in IO stream detected, exception message is: "
	 << e.what() << endl
	 << "\nWe perhaps ran out of disk space or hit a disk quota?\n";
    doAbort();
  }
  catch (exception& e)
  {
    cout << "A 'standard' exception occurred (that's NOT normal):\n" << e.what() << "\n\nIf the cause is not immediatly obvious, please contact: bach@chevreux.org\n\n";
    doAbort();
  }
  catch(...){
    cout << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    doAbort();
  }

  return 0;
}
