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

#include "modules/mod_mer.H"
#include "util/fmttext.H"
#include "util/fileanddisk.H"
#include "util/machineinfo.H"

#include <getopt.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "version.H"


using std::cout;
using std::cin;
using std::cerr;
using std::endl;



void MiraMer::usage()
{
//hdiIpPrvb:f:t:o:a:k:L:n:
  cout << "miramer\t(MIRALIB version " << miraversion << ")\n";
  cout << "Author: Bastien Chevreux\t(bach@chevreux.org)\n\n";
  cout << "...\n";
}


/*
void MiraMer::merCreateHashStats(int argc, char ** argv)
{
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  std::string loadfn(argv[optind++]);

  std::vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  auto rgid = ReadGroupLib::getReadGroupID(0);

  NHashStatistics nhs;
  nhs.setupNewAnalysis(32,4,MER_basesperhash,MER_numlearnsteps);
  {
    uint8 ziptype=0;
    std::string ft,pathto,stem;
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
*/

void MiraMer::merCreateHashStats(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merCreateHashStats(int argc, char ** argv)");

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  std::list<std::string> loadfn;
  for(;optind<argc;++optind){
    loadfn.push_back(argv[optind]);
  }

  bool fwdandrev=true;
  int32 memtouse=75;
  cout << "MER_basesperhash " << MER_basesperhash << endl;
  auto bytes=HashStatistics<vhash64_t>::byteSizeOfHash(MER_basesperhash);
  if(bytes==8){
    MER_hs64.computeHashStatistics(loadfn,memtouse,fwdandrev,1,MER_rarekmerearlykill,MER_basesperhash,
				   MER_outmhs,".");
  }else if(bytes==16){
    MER_hs128.computeHashStatistics(loadfn,memtouse,fwdandrev,1,MER_rarekmerearlykill,MER_basesperhash,
				    MER_outmhs,".");
  }else if(bytes==32){
    MER_hs256.computeHashStatistics(loadfn,memtouse,fwdandrev,1,MER_rarekmerearlykill,MER_basesperhash,
				    MER_outmhs,".");
  }else if(bytes==64){
    MER_hs512.computeHashStatistics(loadfn,memtouse,fwdandrev,1,MER_rarekmerearlykill,MER_basesperhash,
				    MER_outmhs,".");
  }else{
    MIRANOTIFY(true,"Kmer size " << MER_basesperhash << " with " << bytes << " bytes are not expected here.\n");
  }
}


void MiraMer::merInfoHashStats(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merInfoHashStats(int argc, char ** argv)");


  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  std::list<std::string> loadfn;
  for(;optind<argc;++optind){
    loadfn.push_back(argv[optind]);

    std::string hfn(argv[optind]);
    cout << "File " << hfn << ":";
    try {
      auto mhs=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(hfn);
      cout << "\n  File format version:\t" << static_cast<uint16>(mhs.version)
	   << "\n  Kmer length:\t" << mhs.basesperhash
	   << "\n  Kmer bytes:\t" << mhs.sizeofhash
	   << "\n  Num. kmers:\t" << mhs.numelem
	   << "\n  Sort status:\t" << static_cast<uint16>(mhs.sortstatus)
	   << endl;
    }
    catch(Notify n){
      cout << " not readable or not a mhs file.\n";
    }
  }

}



template<typename TVHASH_T>
void MiraMer::mer_fhs_helper1(int argc, char ** argv,HashStatistics<TVHASH_T> & hs)
{
  int32 trimf=0;
  int32 trimr=0;
  int32 trimtot=0;

  auto numopt=argc-optind;
  if(numopt == 2) {
    trimtot=atoi(argv[optind++]);
  }else if(numopt == 3) {
    trimf=atoi(argv[optind++]);
    trimr=trimf;
    trimtot=atoi(argv[optind++]);
  }else if(numopt == 4) {
    trimf=atoi(argv[optind++]);
    trimr=atoi(argv[optind++]);
    trimtot=atoi(argv[optind++]);
  }

  std::string loadfn(argv[optind++]);
  if(loadfn==MER_outmhs){
    cerr << "Outfile cannot be the same as infile.\n";
    exit(99);
  }
  cout << "Loading " << loadfn << endl;
  hs.loadHashStatistics(loadfn);
  cout << "Trimming" << endl;
  hs.trimHashStatsByFrequencyAND(trimf,trimr,trimtot);
  cout << "Recalc frequencies" << endl;
  hs.recalcFreqStats();
  cout << "24 bit sorting" << endl;
  hs.sortLow24Bit();
  cout << "Saving " << MER_outmhs << endl;
  hs.saveHashStatistics(MER_outmhs);
}

void MiraMer::merFilterHashStats(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merFilterHashStats(int argc, char ** argv)");

  auto numopt=argc-optind;
  if(numopt < 2 || numopt > 4) {
    cerr << argv[0] << ": " << "Usage: filter [-o out] [[minfwd] minrev] mintotal in\n";
    exit(1);
  }

  auto bytes=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(argv[argc-1]).sizeofhash;
  if(bytes==8){
    mer_fhs_helper1(argc,argv,MER_hs64);
  }else if(bytes==16){
    mer_fhs_helper1(argc,argv,MER_hs128);
  }else if(bytes==32){
    mer_fhs_helper1(argc,argv,MER_hs256);
  }else if(bytes==64){
    mer_fhs_helper1(argc,argv,MER_hs512);
  }
}



template<typename TVHASH_T>
void MiraMer::mer_bdbg_helper1(int argc, char ** argv,HashStatistics<TVHASH_T> & hs)
{
  std::string loadfn(argv[optind++]);

  if(loadfn==MER_outmhs){
    cerr << "Outfile cannot be the same as infile.\n";
    exit(99);
  }
  cout << "Loading kmer statistics." << endl;
  hs.loadHashStatistics(loadfn);
  cout << "Calculatig forks." << endl;
  hs.calcKMerForks(1,false);
  hs.buildSDBGraphs();

  {
    std::ofstream fout(MER_outmhs+".raw");
    hs.dumpDBGSeqs(fout);
  }
  {
    std::ofstream fout(MER_outmhs+".bylen");
    hs.sortDBGSeqsByLenDown();
    hs.dumpDBGSeqs(fout);
  }
  {
    std::ofstream fout(MER_outmhs+".bymhc");
    hs.sortDBGSeqsByMHCDown();
    hs.dumpDBGSeqs(fout);
  }
}

void MiraMer::merBuildDBGHashStats(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merBuildDBGHashStats(int argc, char ** argv)");

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Usage: dbg [-o out] in\n";
    exit(1);
  }

  auto bytes=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(argv[argc-1]).sizeofhash;
  if(bytes==8){
    mer_bdbg_helper1(argc,argv,MER_hs64);
  }else if(bytes==16){
    mer_bdbg_helper1(argc,argv,MER_hs128);
  }else if(bytes==32){
    mer_bdbg_helper1(argc,argv,MER_hs256);
  }else if(bytes==64){
    mer_bdbg_helper1(argc,argv,MER_hs512);
  }
}




void MiraMer::merSortHashStats(int argc, char ** argv)
{
  if(argc-optind < 2) {
    cerr << argv[0] << ": " << "Missing name of at least one file.\n";
    exit(1);
  }

//  std::string loadfn(argv[optind]);
//  NHashStatistics nhs;
//  nhs.loadHashStatistics(loadfn);
//  nhs.sortLow24Bit();
//  nhs.saveHashStatistics(loadfn+".sorted",true);
}




template<typename TVHASH_T>
void MiraMer::mer_diff_helper1(int argc, char ** argv)
{
  std::string fn1(argv[optind++]);
  std::string fn2(argv[optind++]);

  boost::filesystem::path fp(fn1);
  std::string nameseta(fp.stem().string());
  fp=fn2;
  std::string namesetb(fp.stem().string());

  cout << nameseta << " " << namesetb << endl;

  int32 trimfr=4;
  int32 trimtot=10;

  HashStatistics<TVHASH_T> hs1;
  cout << "load " << fn1 << endl;
  dateStamp(cout);
  hs1.loadHashStatistics(fn1);

  HashStatistics<TVHASH_T> hs2;
  cout << "load " << fn2 << endl;
  dateStamp(cout);
  hs2.loadHashStatistics(fn2);
  dateStamp(cout);

  cout << "creating subhs" << endl;
  {
    HashStatistics<TVHASH_T> in_a_not_b;
    in_a_not_b.inANotB(hs1,hs2);
    in_a_not_b.sortByCountDown();
    dateStamp(cout);
    {
      std::string outname("in_"+nameseta+"_notin_"+namesetb+".fasta");
      cout << "Saving hashes to FASTA file " << outname << endl;
      std::ofstream fout(outname);
      in_a_not_b.dumpAsFASTA(fout);
    }
    {
      std::string outname("in_"+nameseta+"_notin_"+namesetb+".txt");
      cout << "Saving hashes to text file " << outname << endl;
      std::ofstream fout(outname);
      in_a_not_b.dump(fout);
    }
  }
  dateStamp(cout);
  {
    HashStatistics<TVHASH_T> in_b_not_a;
    in_b_not_a.inANotB(hs2,hs1);
    in_b_not_a.sortByCountDown();
    {
      std::string outname("in_"+namesetb+"_notin_"+nameseta+".fasta");
      cout << "Saving hashes to FASTA file " << outname << endl;
      std::ofstream fout(outname);
      in_b_not_a.dumpAsFASTA(fout);
    }
    {
      std::string outname("in_"+namesetb+"_notin_"+nameseta+".txt");
      cout << "Saving hashes to text file " << outname << endl;
      std::ofstream fout(outname);
      in_b_not_a.dump(fout);
    }
    dateStamp(cout);
  }
}


void MiraMer::merDiffHashStats(int argc, char ** argv)
{
  if(argc-optind < 2) {
    cerr << argv[0] << ": " << "Missing name two files to diff.\n";
    exit(1);
  }

  std::string fn1(argv[optind]);
  std::string fn2(argv[optind+1]);

  auto bytes=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(fn1).sizeofhash;
  auto bytes2=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(fn2).sizeofhash;

  if(bytes!=bytes2){
    cerr << "Unequal kmer sizes. "
	 << fn1 << ": " << ", "
	 << fn2 << ": " << ".\n";
    exit(10);
  }

  if(bytes==8){
    mer_diff_helper1<vhash64_t>(argc,argv);
  }else if(bytes==16){
    mer_diff_helper1<vhash128_t>(argc,argv);
  }else if(bytes==32){
    mer_diff_helper1<vhash256_t>(argc,argv);
  }else if(bytes==64){
    mer_diff_helper1<vhash512_t>(argc,argv);
  }

  return;
}


void MiraMer::merDumpHashStats(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merDumpHashStats(int argc, char ** argv)");

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  std::string loadfn(argv[optind]);
  auto bytes=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(loadfn).sizeofhash;
  if(bytes==8){
    MER_hs64.loadHashStatistics(loadfn);
    MER_hs64.dump(cout);
  }else if(bytes==16){
    MER_hs128.loadHashStatistics(loadfn);
    MER_hs128.dump(cout);
  }else if(bytes==32){
    MER_hs256.loadHashStatistics(loadfn);
    MER_hs256.dump(cout);
  }else if(bytes==64){
    MER_hs512.loadHashStatistics(loadfn);
    MER_hs512.dump(cout);
  }else{
    MIRANOTIFY(true,"Kmer size " << MER_basesperhash << " with " << bytes << " bytes are not expected here.\n");
  }
}

void MiraMer::merDumpDebug(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merDumpDebug(int argc, char ** argv)");

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  std::string loadfn(argv[optind]);
  auto bytes=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(loadfn).sizeofhash;
  if(bytes==8){
    MER_hs64.loadHashStatistics(loadfn);
    MER_hs64.dumpHSDebug(cout);
  }else if(bytes==16){
    MER_hs128.loadHashStatistics(loadfn);
    MER_hs128.dumpHSDebug(cout);
  }else if(bytes==32){
    MER_hs256.loadHashStatistics(loadfn);
    MER_hs256.dumpHSDebug(cout);
  }else if(bytes==64){
    MER_hs512.loadHashStatistics(loadfn);
    MER_hs512.dumpHSDebug(cout);
  }else{
    MIRANOTIFY(true,"Kmer size " << MER_basesperhash << " with " << bytes << " bytes are not expected here.\n");
  }
}

void MiraMer::merDumpHashDistrib(int argc, char ** argv)
{
  FUNCSTART("void MiraMer::merDumpHashDistrib(int argc, char ** argv)");
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  std::string loadfn(argv[optind]);
  auto bytes=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(loadfn).sizeofhash;
  if(bytes==8){
    MER_hs64.loadHashStatistics(loadfn);
    MER_hs64.showHashStatisticsInfo();
  }else if(bytes==16){
    MER_hs128.loadHashStatistics(loadfn);
    MER_hs128.showHashStatisticsInfo();
  }else if(bytes==32){
    MER_hs256.loadHashStatistics(loadfn);
    MER_hs256.showHashStatisticsInfo();
  }else if(bytes==64){
    MER_hs512.loadHashStatistics(loadfn);
    MER_hs512.showHashStatisticsInfo();
  }else{
    MIRANOTIFY(true,"Kmer size " << MER_basesperhash << " with " << bytes << " bytes are not expected here.\n");
  }
}

void MiraMer::merDeltaTest(int argc, char ** argv)
{
  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing name of input.\n";
    exit(1);
  }

  cerr << argv[0] << ": " << "merDeltaTest() Currently de-activated.\n";
  exit(1);
/*
  TODO
  De-activated until KMer gets a - operator

  std::string loadfn(argv[optind]);
  HashStatistics<vhash512_t> hs;
  hs.loadHashStatistics(loadfn);
  hs.sortLexicographically();
  auto & hsd=hs.getHashStats();
  if(!hsd.empty()){
    uint64 oldvh=0;
    for(auto & hsde : hsd){
      //cout << oldvh << " " << hsde.vhash << "\t->\t";
      uint64 newvh=hsde.vhash-oldvh;
      oldvh=hsde.vhash;
      hsde.vhash=newvh;
      //cout << oldvh << " " << newvh << endl;
    }
    hs.sortLexicographically();
    hs.saveHashStatistics(loadfn+".delta.sorted",true);
  }
*/
}


int MiraMer::mainMiraMer(int argc, char ** argv)
{
  // that loop is straight from the GNU getopt_long example
  // http://www.gnu.org/s/hello/manual/libc/Getopt-Long-Option-Example.html
  while (1){
    static struct option mlong_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"job", required_argument,         0, 'j'},
	{"kmersize", required_argument,         0, 'k'},
	{"out", required_argument,         0, 'o'},
	{"version", no_argument,         0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hc:j:k:o:v",
		     mlong_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 'h':
      cout << "mira\t\tMIRALIB version " << miraversion << "\n"
	"Author:\t\tBastien Chevreux (bach@chevreux.org)\n"
	"Purpose:\thandle k-mer statistics of a data set\n\n";

      cout << "Usage:\n"
	"miramer [-hv] [-j job] [-o outfile] ...\n";
      cout << "\nOptions:\n";
      cout <<
	"  -h / --help\t\t\t\tPrint short help and exit\n"
	"  -v / --version\t\t\tPrint version and exit\n"
	"  -j / --job\t\t\t\tJob type. Currently:\n"
	"            \t\t\t\tcreate (default)\n"
	"            \t\t\t\tfilter\n"
	"            \t\t\t\tinfo\n"
	"            \t\t\t\tsort\n"
	"            \t\t\t\tdiff\n"
	"            \t\t\t\tdumpcounts\n"
	"            \t\t\t\tdumpdistrib\n"
	"            \t\t\t\tdebug\n"
	"            \t\t\t\tdtest\n"
	"  -o / --out\t\t\t\tOutfile (MHS)\n"
	;
      exit(0);
    case 'j': {
      MER_job=optarg;
      boost::to_lower(MER_job);
      break;
    }
    case 'o': {
      MER_outmhs=optarg;
      break;
    }
    case 'k': {
      uint64 bla=atoi(optarg);
      //if(bla>32) bla=32;
      MER_basesperhash=bla;
      break;
    }
    case 'c': {
      uint64 bla=atoi(optarg);
      //if(bla>32) bla=32;
      MER_rarekmerearlykill=bla;
      break;
    }
    case 'v':
      cout << miraversion << endl;
      exit(0);
    default:
      abort();
    }
  }

  if(MER_optthreads==0){
    MER_optthreads=MachineInfo::getCoresTotal();
  }
#ifdef HAVE_OPENMP
  omp_set_num_threads(MER_optthreads);
#endif


  if(MER_basesperhash==0) MER_basesperhash=32;
  if(MER_basesperhash>256){
    cout << "Sorry, -k for kmer size must be <= 256 for the time being.\n";
    exit(100);
  }

  try {
    if(MER_job=="create"){
      merCreateHashStats(argc,argv);
    }else if(MER_job=="dbg"){
      merBuildDBGHashStats(argc,argv);
    }else if(MER_job=="diff"){
      merDiffHashStats(argc,argv);
    }else if(MER_job=="dumpcounts"){
      merDumpHashStats(argc,argv);
    }else if(MER_job=="dumpdistrib"){
      merDumpHashDistrib(argc,argv);
    }else if(MER_job=="filter"){
      merFilterHashStats(argc,argv);
    }else if(MER_job=="info"){
      merInfoHashStats(argc,argv);
    }else if(MER_job=="sort"){
      merSortHashStats(argc,argv);
    }else if(MER_job=="debug"){
      merDumpDebug(argc,argv);
    }else if(MER_job=="dtest"){
      merDeltaTest(argc,argv);
    }else{
      cout << argv[0] << ": unknown job '" << MER_job << "'???" << endl;
      exit(1);
    }
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
