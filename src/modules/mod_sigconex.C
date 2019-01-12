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

#include "modules/mod_sigconex.H"
#include "util/fmttext.H"
#include "util/fileanddisk.H"

#include <getopt.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <regex>

#include "version.H"


using namespace std;




void MiraSCE::usage()
{
//hdiIpPrvb:f:t:o:a:k:L:n:
  cout << "miramer\t(MIRALIB version " << miraversion << ")\n";
  cout << "Author: Bastien Chevreux\t(bach@chevreux.org)\n\n";
  cout << "...\n";
}


struct scinfo_t {
  uint32 rid=0;
  uint32 len=0;
  double cov=0.0;
  bool learned=false;

  inline static bool sortComparatorByLen(const scinfo_t & a, const scinfo_t & b){
    if(a.len==b.len){
      return a.cov<b.cov;
    }
    return a.len<b.len;
  }
};

void MiraSCE::excon(int argc, char ** argv)
{
  FUNCSTART("void MiraSCE::excon(int argc, char ** argv)");

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  list<string> loadfn;
  for(;optind<argc;++optind){
    loadfn.push_back(argv[optind]);
  }

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  auto rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);

  cout << "Loading data into memory ...";
  ReadPool loadrp;
  ReadPoolIO rpio(loadrp);
  rpio.setAttributeFASTAQualFileWanted(false); // in case we load FASTAs
  rpio.setAttributeFASTQAPreserveComment(true); // want to store comments in COMM read tags

  for(auto & dfn : loadfn){
    uint8 ziptype=0;
    string ft;
    string dummyfromstem;
    string dummypathto;
    guessFileAndZipType(dfn,dummypathto,dummyfromstem,ft,ziptype);

    rpio.registerFile(ft,dfn,"",rgid,false);
    rpio.loadNextSeqs(-1,-1);
  }

  bool allcovfound=true;
  std::regex rgx("cov=([\\d\\.]+)");
  std::smatch match;
  vector<scinfo_t> scinfo(loadrp.size());
  for(uint32 rpi=0; rpi<loadrp.size(); ++rpi){
    Read & actread=loadrp[rpi];
    scinfo[rpi].rid=rpi;
    scinfo[rpi].len=actread.getLenSeq();
    string rcomm;
    for(auto & rt : actread.getTags()){
      if(rt.identifier==Read::REA_tagentry_idCOMM){
	rcomm=rt.getCommentStr();
	break;
      }
    }
    if(!rcomm.empty()){
      if (std::regex_search(rcomm, match, rgx)) {
	scinfo[rpi].cov=std::stod(match[1]);
      }else{
	cout << actread.getName() << " is missing a coverage number in the form of 'cov=xx.xxx' in its comments.\n";
	allcovfound=false;
      }
    }
  }

  if(!allcovfound){
    cout << "Sorry, fix that.\n";
    exit(100);
  }

  double mincov=34.00;
  for(auto & sce : scinfo){
    if(sce.cov>=mincov){
      // learn seq
      sce.learned=true;
    }
  }

  sort(scinfo.begin(),scinfo.end(),scinfo_t::sortComparatorByLen);

  for(uint32 downlen=32; downlen>=8; downlen/=2){
    for(uint32 downcov=32; downcov>=8; downcov/=2){
      auto actmincov=static_cast<double>(downcov);
      for(auto & sce : scinfo){
	if(!sce.learned && sce.cov>=actmincov && sce.len>=downlen){
	  // learn seq
	  sce.learned=true;
	}
      }
    }
  }

  for(uint32 downlen=6; downlen>=2; downlen-=2){
    for(uint32 downcov=32; downcov>=8; downcov/=2){
      auto actmincov=static_cast<double>(downcov);
      for(auto & sce : scinfo){
	if(!sce.learned && sce.cov>=actmincov && sce.len>=downlen){
	  // learn seq
	  sce.learned=true;
	}
      }
    }
  }

  // check & select remaining seqs

}




int MiraSCE::mainSignificantContigExtraction(int argc, char ** argv)
{
  FUNCSTART("int MiraSCE::mainSignificantContigExtraction(int argc, char ** argv)");

  // that loop is straight from the GNU getopt_long example
  // http://www.gnu.org/s/hello/manual/libc/Getopt-Long-Option-Example.html
  while (1){
    static struct option mlong_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"kmersize", required_argument,         0, 'k'},
	{"version", no_argument,         0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hc:j:k:v",
		     mlong_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 'h':
      cout << "mira\t\tMIRALIB version " << miraversion << "\n"
	"Author:\t\tBastien Chevreux (bach@chevreux.org)\n"
	"Purpose:\t...\n\n";

      cout << "Usage:\n"
	",,,, ...\n";
      cout << "\nOptions:\n";
      cout <<
	"  -h / --help\t\t\t\tThis help text\n"
	"  -v / --version\t\t\tPrint version and exit\n"
	;
      exit(0);
//    case 'j': {
//      SCE_job=optarg;
//      boost::to_lower(SCE_job);
//      break;
//    }
    case 'k': {
      uint64 bla=atoi(optarg);
      //if(bla>32) bla=32;
      SCE_basesperhash=bla;
      break;
    }
    case 'v':
      cout << miraversion << endl;
      exit(0);
    default:
      abort();
    }
  }

  if(SCE_basesperhash==0) SCE_basesperhash=32;
  if(SCE_basesperhash>32){
    cout << "Sorry, -k for kmer size must be <= 32 for the time being.\n";
    exit(100);
  }

  try {
    excon(argc,argv);
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
  catch(Flow f){
    cerr << "Unexpected exception: Flow()\n";
  }
  catch(...){
    cerr << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    abort();
  }

  FUNCEND();
  return 0;
}
