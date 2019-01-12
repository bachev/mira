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

#include <getopt.h>

#include "modules/mod_diff.H"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "util/fmttext.H"
#include "mira/readpool_io.H"

#include "version.H"


using namespace std;




MiraDiff::MiraDiff()
{
}

void MiraDiff::usage()
{
//hdiIpPrvb:f:t:o:a:k:L:n:
  cout << "miramer\t(MIRALIB version " << miraversion << ")\n";
  cout << "Author: Bastien Chevreux\t(bach@chevreux.org)\n\n";
  cout << "...\n";
}


void MiraDiff::priv_clipsAfterLoad(ReadPool & myrp)
{
  MIRAParameters::setupStdMIRAParameters(MD_miraparams);
  MIRAParameters::parse("--job=genome,denovo,accurate,solexa",MD_miraparams,false);
  MD_miraparams[0].getNonConstDirectoryParams().dir_tmp=".";  // for DataProcessing
  //MIRAParameters::dumpAllParams(MD_miraparams,cout);
  //cout.flush();

  //auto numthreads=MD_miraparams[0].getAssemblyParams().as_numthreads;
  auto numthreads=4;

  vector<uint8> debrisreason;
  debrisreason.clear();
  debrisreason.resize(myrp.size(),Assembly::DEBRIS_NOTDEBRIS);

  vector<unique_ptr<DataProcessing>> dpv(numthreads);
  DataProcessing dpcollector(&MD_miraparams);

  for(uint32 ti=0; ti<dpv.size(); ++ti){
    dpv[ti]=unique_ptr<DataProcessing>(new DataProcessing(&MD_miraparams));
    string logname("mdclip_t" + boost::lexical_cast<string>(ti) + ".txt");
    cout << logname << endl;
    dpv[ti]->startLogging(logname,false);
  }

  string logprefix("loadclip: ");
  cout << "Post-load clips:\n";

  DataProcessing::stdTreatmentPool_MultiThread(MD_miraparams,dpcollector,dpv,myrp,&debrisreason,logprefix,true);

  cout << endl;
  if(dpcollector.DP_stats.cphix174){
    cout << "SEARCH MSG: PhiX 174 found: " << dpcollector.DP_stats.cphix174 << endl;
  }
  if(dpcollector.DP_stats.cadapright){
    cout << "CLIP MSG: Adaptor right found: " << dpcollector.DP_stats.cadapright << endl;
  }
  if(dpcollector.DP_stats.cadaprightpartial){
    cout << "CLIP MSG: Partial adaptor right found: " << dpcollector.DP_stats.cadaprightpartial << endl;
  }
}


void MiraDiff::priv_simpleDiff(int argc, char ** argv)
{
  list<string> loadfs1;

  uint32 numfiles=1;
  uint8 basesperhash=31;

  int32 trimfr=4;
  int32 trimtot=20;

  string nameseta("A");
  string namesetb("B");

  for(auto fc=0;optind<argc && fc<numfiles;++optind,++fc){
    loadfs1.push_back(argv[optind]);
  }

  list<string> loadfs2;
  for(auto fc=0;optind<argc && fc<numfiles;++optind,++fc){
    loadfs2.push_back(argv[optind]);
  }

  auto rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_SOLEXA);

  cout << "Loading data set " << nameseta << " into memory:\n";
  ReadPool myrp;
  ReadPoolIO rpio(myrp);
  rpio.setAttributeFASTAQualFileWanted(false); // in case we load FASTAs
  rpio.setAttributeFASTQQualOffset(33); // in case we load FASTQs

  for(auto & dfn : loadfs1){
    cout << dfn << endl;
    uint8 ziptype=0;
    string ft;
    string dummyfromstem;
    string dummypathto;
    guessFileAndZipType(dfn,dummypathto,dummyfromstem,ft,ziptype);

    rpio.registerFile(ft,dfn,"",rgid,false);
    rpio.loadNextSeqs(-1,-1);
  }
  cout << "\nDone\n";
  cout << "Clipping pool " << nameseta << ":\n";
  priv_clipsAfterLoad(myrp);
  cout << "\nDone\n";
  cout << "Calculating hash statistics for pool " << nameseta << " :\n";

  string dummyfn;
  bool MB_fwdandrev=true;
  HashStatistics<vhash64_t> hs1;
  hs1.computeHashStatistics(myrp,
			    75,
			    false,false,MB_fwdandrev,1,0,basesperhash,
			    dummyfn,".");
  cout << "\nDone\nSaving ...";cout.flush();
  hs1.saveHashStatistics(nameseta+".mhs.gz");
  cout << "done." << endl;

  hs1.trimHashStatsByFrequencyANDOR(trimfr,trimfr,trimtot);

  cout << "hash distrib " << nameseta << " :\n";
  //hs1.dumpHashDistrib(cout);
  hs1.sortByCountDown();
  //hs1.dump(cout);

//  hs1.sortLexicographically();
//  cout << "########## HS1\n";
//  hs1.dump(cout);


  cout << "Loading data set " << namesetb << " into memory:\n";

  myrp.discard();
  for(auto & dfn : loadfs2){
    cout << dfn << endl;
    uint8 ziptype=0;
    string ft;
    string dummyfromstem;
    string dummypathto;
    guessFileAndZipType(dfn,dummypathto,dummyfromstem,ft,ziptype);

    rpio.registerFile(ft,dfn,"",rgid,false);
    rpio.loadNextSeqs(-1,-1);
  }
  cout << "\nDone\n";
  cout << "Clipping pool " << namesetb << " :\n";
  priv_clipsAfterLoad(myrp);
  cout << "\nDone\n";

  cout << "Calculating hash statistics for pool " << namesetb << " :\n";
  HashStatistics<vhash64_t> hs2;
  hs2.computeHashStatistics(myrp,
			    75,
			    false,false,MB_fwdandrev,1,0,basesperhash,
			    dummyfn,".");
  cout << "\nDone\nSaving ...";cout.flush();
  hs2.saveHashStatistics(namesetb+".mhs.gz");
  cout << "done." << endl;

  myrp.discard();

  hs2.trimHashStatsByFrequencyANDOR(trimfr,trimfr,trimtot);
  cout << "hash distrib " << namesetb << " :\n";

//  cout << "########## HS2\n";
  hs2.sortByCountDown();
//  hs2.dump(cout);

//  hs2.sortLexicographically();
//  cout << "########## HS2\n";
//  hs2.dump(cout);

  cout << "creating subhs" << endl;
  HashStatistics<vhash64_t> in_a_not_b;
  in_a_not_b.inANotB(hs1,hs2);
  in_a_not_b.sortByCountDown();
  cout << "\n######## " << nameseta << " not " << namesetb << "\n";
  in_a_not_b.dump(cout);

  HashStatistics<vhash64_t> in_b_not_a;
  in_b_not_a.inANotB(hs2,hs1);
  in_b_not_a.sortByCountDown();
  cout << "\n######## " << namesetb << " not " << nameseta << "\n";
  in_b_not_a.dump(cout);

//  HashStatistics<vhash64_t> in_a_and_b;
//  in_a_and_b.inAAndB(hs1,hs2);
//  in_a_and_b.sortByCountDown();
//  cout << "\n######## A and B\n";
//  in_a_and_b.dump(cout);

  cout << "\nDone\n";

}

int MiraDiff::mainMiraDiff(int argc, char ** argv)
{
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
      cout << "mira\t\tMIRALIB version " << miraversion << "\n"
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
      //MER_job=optarg;
      //boost::to_lower(MER_job);
      break;
    }
    case 'k': {
      uint64 bla=atoi(optarg);
      if(bla>32) bla=32;
      //MER_basesperhash=bla;
      break;
    }
    case 'v':
      cout << miraversion << endl;
      exit(0);
    default:
      abort();
    }
  }

  try {
    priv_simpleDiff(argc,argv);
  }
  catch(Notify n){
    n.handleError("main");
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
