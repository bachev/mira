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

#include <signal.h>

#include <getopt.h>

#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "modules/mod_sqt.H"
#include "util/fmttext.H"
#include "util/machineinfo.H"

#include "version.H"


using std::cout;
using std::cin;
using std::cerr;
using std::endl;


struct sigaction msqt_sigIntHandler;

std::vector<MIRAParameters> MiraSQT::MS_Pv;
std::vector<std::unique_ptr<DataProcessing>> MiraSQT::MS_dpv;
DataProcessing MiraSQT::MS_dataprocessing(&MiraSQT::MS_Pv);

uint8 MiraSQT::MS_tortype;

std::list<std::string> MiraSQT::MS_infiles;
std::string MiraSQT::MS_hitpath;
std::string MiraSQT::MS_nameprefix;

std::string MiraSQT::MS_tmpdirname{"."};
std::string MiraSQT::MS_hashstatfname;

bool   MiraSQT::MS_dustfilter=false;
bool   MiraSQT::MS_mustdeletetargetfiles=true;
bool   MiraSQT::MS_fwdandrev=true;
bool   MiraSQT::MS_changeseqcase=true;

bool   MiraSQT::MS_testflag=false;

int16  MiraSQT::MS_dustperc=67;

uint32 MiraSQT::MS_optthreads=0;
int32  MiraSQT::MS_optmbtouse=75;

std::list<char> MiraSQT::MS_filepairinfo;    // p = 2 files, P = 1 file interleave
uint32 MiraSQT::MS_basesperhash=0;
uint64 MiraSQT::MS_numreadsread=0;

std::list<MiraSQT::wqueueunit_t> MiraSQT::MS_workqueue;
MiraSQT::files_t MiraSQT::MS_files;

std::unordered_map<std::string,uint8> MiraSQT::MS_fromtypemap = {
  {"fastq",Read::AS_FASTQ},
  {"fasta",Read::AS_FASTA},
};

std::unordered_map<std::string,uint8> MiraSQT::MS_totypemap = {
  {"fastq",Read::AS_FASTQ},
  {"fasta",Read::AS_FASTA},
};

boost::mutex msqt_mutex_ctrlc;

//MIRATimer dbg_mt;
//auto dbg_baittime=dbg_mt.diff();

HashStatistics<vhash64_t> msqt_hs64;
HashStatistics<vhash128_t> msqt_hs128;
HashStatistics<vhash256_t> msqt_hs256;
HashStatistics<vhash512_t> msqt_hs512;

bool MiraSQT::MS_signal_ctrlc=false;


MiraSQT::~MiraSQT()
{
}

// Do not use exit()!
// Do not throw an exeption! It sometimes works and the exeption is caught in the main program,
//  but sometimes crashes to command line with
//     terminate called after throwing an instance of Notify
// So, HashStatistics has own abort functionality and we'll exit the regular way
void MiraSQT::ctrlCHandler(int s)
{
  FUNCSTART("void MiraSQT::ctrlCHandler(int s)");
  cout << "Control-C was caught. Cleaning up, this may take a few seconds." << endl;
  MS_signal_ctrlc=true;
  msqt_hs64.abortAll();
  msqt_hs128.abortAll();
  msqt_hs256.abortAll();
  msqt_hs512.abortAll();
}


void MiraSQT::usage()
{
  cout << "mirasqt\t(MIRALIB version " << miraversion << ")\n";
  cout << "Author: Bastien Chevreux\t(bach@chevreux.org)\n\n";

  cout << FmtText::wordWrap("MIRA Streaming Quality\n");
  cout << "\nUsage:\n";
  cout << "mirasqt [options]"
    "\n\t\t{-b baitfile [-b ...] | -B file}"
    "\n\t\t{-p file_1 file_2 | -P file3}* [file4 ...]\n\n";
  cout << "Main options:\n";
  cout <<
    "\t-B file\t\tLoad kmer statistics file\n"
    ;

  cout << "\n"
    "\t-p file1 file2\tLoad paired sequences to search from file1 and file2\n"
    "\t\t\t Files must contain same number of sequences, sequence \n"
    "\t\t\t names must be in same order.\n"
    "\t\t\t Multiple -p allowed, but must come before non-paired\n"
    "\t\t\t files.\n"
    "\t-P file\t\tLoad paired sequences from file\n"
    "\t\t\t File must be interleaved: pairs must follow each other,\n"
    "\t\t\t non-pairs are not allowed.\n"
    "\t\t\t Multiple -p allowed, but must come before non-paired\n"
    "\t\t\t files.\n";

//  cout << "\n"
//    "\t-k int\t\tkmer length of bait in bases (<=256, default=31)\n"
//    "\t-n int\t\tIf >0: minimum number of k-mer baits needed (default=1)\n"
//    "\t\t\tIf <=0: allowed number of missed kmers over sequence\n"
//    "\t\t\t        length\n";

  cout << "\n"
    "\t-d\t\tDo not use kmers with microrepeats (DUST-like, see also -D)\n"
    "\t-D int\t\tSet length of microrepeats in kmers to discard from bait.\n"
    "\t\t\t int > 0 microrepeat len in percentage of kmer length.\n"
    "\t\t\t       E.g.: -k 17 -D 67 --> 11.39 bases --> 12 bases.\n"
    "\t\t\t int < 0 microrepeat len in bases.\n"
    "\t\t\t int != 0 implies -d, int=0 turns DUST filter off.\n"
    "\t-t\t\tNumber of threads to use (default=0 -> up to 4 CPU cores)\n";

  cout << "\nOptions for output definition:\n";
  cout << FmtText::wordWrap("Bla bla");
  cout << "\t-N name\t\tChange the prefix 'msqt' to <name>\n"
    "\t\t\t Has no effect if -o/-O is used and targets are not\n"
    "\t\t\t directories\n";
  cout << "\t-o <path>\tSave sequences to path\n"
    "\t\t\t If path is a directory, write separate files into this\n"
    "\t\t\t directory. If not, combine all matching sequences from\n"
    "\t\t\t the input file(s) into a single file specified by the\n"
    "\t\t\t path.\n";

  cout << "\nOther options:\n";
  cout << "\t-T dir\t\tUse 'dir' as directory for temporary files instead of\n"
    "\t\t\t current working directory.\n";

  cout << "\nDefining files types to load/save:\n";
  cout << FmtText::wordWrap("Normally mirabait recognises the file types according to the file extension (even when packed). In cases you need to force a certain file type because the file extension is non-standard, use the EMBOSS notation to force a type: <filetype>::<name_of_file>. E.g., to tell that \"somefile.dat\" is FASTQ, use: fastq::somefile.dat\nRecognised types are: caf, fasta, fastq, gbf, gbk, gbff, maf and phd.\n\nMIRABAIT will write files in the same file type as the corresponding input files.");

  cout << "\nExamples:\n"
    "  mirabait -b b.fasta file.fastq"
    "\n  mirabait -I -j rrna -p file_1.fastq file_2.fastq"
    "\n";
}


uint8 MiraSQT::checkFromType(std::string & fromtype)
{
  auto mI=MS_fromtypemap.find(fromtype);
  if(mI!=MS_fromtypemap.end()) return mI->second;
  return Read::AS_ILLEGAL;
}

uint8 MiraSQT::checkToType(std::string & totype)
{
  auto mI=MS_totypemap.find(totype);
  if(mI!=MS_totypemap.end()) return mI->second;
  return Read::AS_ILLEGAL;
}



uint8 MiraSQT::setupRPIO(std::string filename, ReadGroupLib::ReadGroupID rgid, ReadPoolIO & rpio, uint8 & ziptype)
{
  ziptype=0;
  std::string filetype;
  std::string dummytostem;
  std::string dummypathto;

  filename=guessFileAndZipType(filename,dummypathto,dummytostem,filetype,ziptype);
  boost::to_lower(filetype);
  uint8 rtype=checkFromType(filetype);

  cout << "Loading data from " << filename << " (" << filetype << ")";

  std::string fn2;
  if(filetype=="fasta"){
    fn2=filename+".qual";
  }

  std::string loadtype(filetype);
  if(loadtype=="fasta"){
    loadtype="fastanoqual";
  }

  rpio.registerFile(loadtype, filename, fn2, rgid, false);
  return rtype;
}

void MiraSQT::setupOutfiles(const std::string & fname, uint8 rtype, uint8 ziptype, std::ofstream & tfout)
{
  if(tfout.is_open()) tfout.close();

  std::string pname(MS_nameprefix);
  if(pname.empty()) pname="msqt_";

  boost::filesystem::path fp(fname);
  // if we work on zipped files, the ouput must get away the zip extension
  if(ziptype){
    fp=fp.parent_path() / fp.stem();  // boost::filesystem append operation
  }

  std::string hitname{MS_hitpath};
  if(MS_hitpath.empty()){
    hitname=pname+fp.filename().string();
  }

  auto omode=std::ios::out;
  tfout.open(hitname,omode);
  if(!tfout.is_open()){
    cout.flush();
    cerr << "\n\nCould not open " << hitname << ".\nDoes the path exist, is it writable? Is the disk full?\n";
    exit(10);
  }
}


void MiraSQT::setFakeFlags(ReadPool & rp)
{
  for(uint32 rpi=0; rpi<rp.size(); ++rpi){
    rp[rpi].setUsedInAssembly(true);
  }
}

void MiraSQT::parallelSimpleClip(ReadPool & rp)
{
  static std::string dummy;
  DataProcessing::stdTreatmentPool_MultiThread(MS_Pv,MS_dataprocessing,MS_dpv,rp,NULL,
					       dummy, // no log prefix
					       false  // no progress
    );
}


template<typename TVHASH_T>
void MiraSQT::pecSingleRead(HashStatistics<TVHASH_T> & mbhs, Read & actread)
{
}

// Need to do this here "by hand" as I do not want to use the Read frequency structures / facilities
//  of MIRA (probably not suited / too slow)
template<typename TVHASH_T>
void MiraSQT::parallelProposedEndClip(HashStatistics<TVHASH_T> & mbhs, ReadPool & rp)
{
  cout << "dubidu" << endl;
  mbhs.assignReadBaseStatistics_MultiThread(rp, MS_optthreads,
//  mbhs.assignReadBaseStatistics_MultiThread(rp, 1,
					    false, // mask nasty repeats?
					    false,  // calc kmer forks?
					    0,     // kmerforks: minkmer=0, so take everything! TODO: cmd-line argument??
					    false  // kmerforks: we don't want fwd/rev, relaxed TODO: cmd-line argument??
      );

  cout << "starting PEC" << endl;

  MS_dataprocessing.proposedEndClipping_Pool(rp, MS_basesperhash);

//#pragma omp parallel
//  {
//#pragma omp for
//    for(uint32 rpi=0; rpi<rp.size(); ++rpi){
//      pecSingleRead(mbhs,rp[rpi]);
//    }
//  }
}


void MiraSQT::saveWQueueElement(wqueueunit_t & wqu)
{
  FUNCSTART("void MiraSQT::saveWQueueElement(wqueueunit_t & wqu)");

  Read::setCoutType(MS_files.writetype);

  cout << "\nSaving..." << endl;

  for(uint32 rpi=0; rpi<wqu.rp1.size(); ++rpi){
    wqu.rp1[rpi].performHardTrim();
    MS_files.tfout1 << wqu.rp1[rpi];

    if(wqu.pairstatus==PS_2FILES){
      wqu.rp2[rpi].performHardTrim();
      MS_files.tfout2 << wqu.rp2[rpi];
    }
  }

  MS_numreadsread+=wqu.rp1.size()+wqu.rp2.size();
}



template<typename TVHASH_T>
void MiraSQT::doTreatWithHS(HashStatistics<TVHASH_T> & mbhs)
{
  FUNCSTART("void MiraSQT::doTreatWithHS(HashStatistics<TVHASH_T> & mbhs)");

  MS_dpv.resize(MS_optthreads);
  for(uint32 ti=0; ti<MS_dpv.size(); ++ti){
    MS_dpv[ti]=std::unique_ptr<DataProcessing>(new DataProcessing(&MS_Pv));
  }

//  // setup the CTRL-C handler only NOW
//  msqt_sigIntHandler.sa_handler = MiraSQT::ctrlCHandler;
//  sigemptyset(&msqt_sigIntHandler.sa_mask);
//  msqt_sigIntHandler.sa_flags = 0;
//  sigaction(SIGINT, &msqt_sigIntHandler, NULL);

  if(!MS_hashstatfname.empty()){
    cout << "Loading from existing hashstat file ... "; cout.flush();
    mbhs.loadHashStatistics(MS_hashstatfname);
    cout << "done.\n";
    MS_basesperhash=mbhs.getBasesPerHash();

    mbhs.showHashStatisticsInfo();
    if(MS_dustfilter){
      if(MS_dustperc>=0){
	mbhs.removeDustEndsRatio(0,MS_dustperc);
      }else{
	mbhs.removeDustEndsFixed(0,-MS_dustperc);
      }
      mbhs.showHashStatisticsInfo();
    }
  }

  if(MS_signal_ctrlc) return;

  //dbg_baittime=0;

  MS_workqueue.resize(1);
  auto qI=MS_workqueue.begin(); // fixed atm

  // calling this on an empty readpool won't do anything to the readpool
  //  but internally, HashStats will assign kmer fork flags to the kmers
  // Only needs to be done once
  if(mbhs.getNumHashEntries()){
    cout << "Creating kmer forks ..."; cout.flush();
    mbhs.assignReadBaseStatistics_MultiThread(qI->rp1, MS_optthreads,
					      false, // mask nasty repeats?
					      true,  // calc kmer forks?
					      0,     // kmerforks: minkmer=0, so take everything! TODO: cmd-line argument??
					      false  // kmerforks: we don't want fwd/rev, relaxed TODO: cmd-line argument??
      );
    cout << " done" << endl;
  }

  ReadPoolIO rpio1(qI->rp1);
  ReadPoolIO rpio2(qI->rp2);
  rpio1.setAttributeFASTAQualFileWanted(false); // in case we load FASTAs
  rpio2.setAttributeFASTAQualFileWanted(false); // in case we load FASTAs
  rpio1.setAttributeFASTQQualOffset(33); // in case we load FASTQs
  rpio2.setAttributeFASTQQualOffset(33); // in case we load FASTQs
  rpio1.setAttributeFASTQTransformName(false); // no name transform (+ loading faster)
  rpio2.setAttributeFASTQTransformName(false); // no name transform (+ loading faster)
  rpio1.setAttributeFASTQAPreserveComment(true); // in case we load FASTA / Qs
  rpio2.setAttributeFASTQAPreserveComment(true); // in case we load FASTA / Qs
  rpio1.setAttributeFASTQACheckQuals(false); // no checks, make loading faster
  rpio2.setAttributeFASTQACheckQuals(false); // no checks, make loading faster
  ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_SOLEXA);
  rgid.setReadNamingScheme(ReadGroupLib::SCHEME_EMPTY);  // saves time during read.setName()

  auto ifI=MS_infiles.begin();
  while(ifI!=MS_infiles.end()){
    if(MS_signal_ctrlc) return;
    uint8 ziptype=0;
    rpio1.setAttributeProgressIndicator(true);
    MS_files.infilename1.clear();
    MS_files.infilename2.clear();

    MS_files.intype1=setupRPIO(*ifI,rgid,rpio1,ziptype);
    MS_files.writetype=MS_files.intype1;
    setupOutfiles(*ifI,MS_files.intype1,ziptype,MS_files.tfout1);
    MS_files.infilename1=*ifI;

    qI->pairstatus=PS_NOPAIR;
    if(!MS_filepairinfo.empty()){
      if(MS_filepairinfo.front()=='P'){
	qI->pairstatus=PS_INTERLEAVE;
      }else if(MS_filepairinfo.front()=='p'){
	++ifI;
	if(ifI==MS_infiles.end()){
	  MIRANOTIFY(Notify::FATAL,"Something's wrong here: -p says to expect one further file, but " << *(--ifI) << " is the last file seen on the command line?");
	}
	rpio1.setAttributeProgressIndicator(false);
	rpio2.setAttributeProgressIndicator(true);
	MS_files.intype2=setupRPIO(*ifI,rgid,rpio2,ziptype);
	setupOutfiles(*ifI,MS_files.intype2,ziptype,MS_files.tfout2);
	MS_files.infilename2=*ifI;
	qI->pairstatus=PS_2FILES;
      }
      MS_filepairinfo.pop_front();
    }

    ++ifI;

    uint32 poolsize=100000;

    while(rpio1.loadNextSeqs(poolsize)){
      if(MS_signal_ctrlc) return;
      if(qI->pairstatus==PS_2FILES){
	rpio2.loadNextSeqs(poolsize);
	if(qI->rp1.size() != qI->rp2.size()){
	  MIRANOTIFY(Notify::FATAL,"Something's wrong here: -p says you have two files with reads paired across both files. But file " << MS_files.infilename1 << " does not have the same number of reads as file " << MS_files.infilename2 << " ???");
	}
      }

      if(qI->rp1.size()==0) break;

      // name or template checks
      {
	if(qI->pairstatus==PS_2FILES){
	  for(uint32 rpi=0; rpi<qI->rp1.size(); ++rpi){
	    if(qI->rp1[rpi].getName()!=qI->rp2[rpi].getName()
	       && qI->rp1[rpi].getTemplate()!=qI->rp2[rpi].getTemplate()){
	      MIRANOTIFY(Notify::FATAL,"Paired end files not synchronised: read name " << qI->rp1[rpi].getName() << " not equal to " << qI->rp2[rpi].getName() << " and templates also do not match: " << qI->rp1[rpi].getTemplate() << " vs " << qI->rp2[rpi].getTemplate());
	    }
	  }
	}else if(qI->pairstatus==PS_INTERLEAVE){
	  if(qI->rp1.size()%2){
	    MIRANOTIFY(Notify::FATAL,"Interleaved paired end file apparently not cleanly interleaved: last read " << qI->rp1[qI->rp1.size()-1].getName() << " does not have a partner.");
	  }
	  for(uint32 rpi=0; rpi<qI->rp1.size(); rpi+=2){
	    if(qI->rp1[rpi].getTemplate()!=qI->rp1[rpi+1].getTemplate()){
	      MIRANOTIFY(Notify::FATAL,"Interleaved paired end file apparently not cleanly interleaved: read template " << qI->rp1[rpi].getTemplate() << " not equal to " << qI->rp1[rpi+1].getTemplate());
	    }
	  }
	}
      }

      setFakeFlags(qI->rp1);
      parallelSimpleClip(qI->rp1);
      if(mbhs.getNumHashEntries()) {
	parallelProposedEndClip(mbhs,qI->rp1);
      }
      if(qI->pairstatus==PS_2FILES){
	parallelSimpleClip(qI->rp2);
	if(mbhs.getNumHashEntries()) {
	  parallelProposedEndClip(mbhs,qI->rp2);
	}
      }
      saveWQueueElement(*qI);

      Read::trashReadNameContainer();
      qI->rp1.discard();
      qI->rp2.discard();
    }
    if(qI->pairstatus==PS_2FILES){
      rpio2.loadNextSeqs(1);
      if(qI->rp2.size()){
	MIRANOTIFY(Notify::FATAL,"File (bla2) more reads than file (bla)???");
      }
    }
  }

  cout << endl;

  //cout << "DBG_baittime: " << dbg_baittime << endl;
}


int MiraSQT::mainMiraSQT(int argc, char ** argv)
{
  //CALLGRIND_STOP_INSTRUMENTATION;

  FUNCSTART("int mainMiraSQT(int argc, char ** argv)");

  int c;
  extern char *optarg;
  extern int optind;


  std::string path;
  std::string convertprog;
  splitFullPathAndFileName(argv[0],path,convertprog);

  std::string miraparams;

  MS_basesperhash=0;

  Read::setCoutLen(0);

  MS_filepairinfo.clear();    // p = 2 files, P = 1 file interleave

  while (true){
    static struct option long_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"version", no_argument,         0, 'v'},
	{"throw", no_argument,         0, '{'},
	{"loadhsf", no_argument,    0, 'L'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "{hcdpPrva:B:D:k:l:m:N:t:T:",
			 long_options, &option_index);
    if(c == -1) break;

    switch (c) {
    case '{': {
      Notify::setBangOnThrow(true);
      break;
    }
    case '1': {
      MS_testflag=true;
      break;
    }
    case 'a': {
      miraparams=optarg;
      break;
    }
    case 'B': {
      MS_hashstatfname=optarg;
      break;
    }
    case 'c': {
      MS_changeseqcase=false;
      break;
    }
    case 'd': {
      MS_dustfilter=true;
      break;
    }
    case 'D': {
      MS_dustfilter=true;
      int64 bla=atoi(optarg);
      if(bla>100
	 || bla<-1000){
	usage();
	cerr << "-D must be >-1000 and <= 100 in the posi\n";
	exit(2);
      }
      if(bla==0) {
	MS_dustfilter=false;
      }else{
	MS_dustperc=bla;
      }
      break;
    }
    case 'k': {
      uint32 bla=atoi(optarg);
      MS_basesperhash=bla;
      break;
    }
    case 'l': {
      Read::setCoutLen(atoi(optarg));
      break;
    }
    case 'm': {
      MS_optmbtouse=atoi(optarg);
      break;
    }
    case 'N': {
      MS_nameprefix=optarg;
      break;
    }
    case 'p' :
    case 'P' : {
      MS_filepairinfo.push_back(c);
      break;
    }
    case 'r': {
      MS_fwdandrev=false;
      break;
    }
    case 't': {
      MS_optthreads=atoi(optarg);
      break;
    }
    case 'T': {
      MS_tmpdirname=optarg;
      break;
    }
    case 'h':
    case '?': {
      usage();
      exit(0);
    }
    case 'v':
      cout << miraversion << endl;
      exit(0);
    default : {}
    }
  }

  if(MS_optthreads==0){
    MS_optthreads=MachineInfo::getCoresTotal();
    if(MS_optthreads>6) MS_optthreads--;
  }
#ifdef HAVE_OPENMP
  omp_set_num_threads(MS_optthreads);
#endif

  if(MS_hashstatfname.empty()){
    cout << argv[0] << ": " << "No -B given (rename that option!)!\nNo PEC clipping.\n";
  }

  if(argc-optind < 1){
    usage();
    cout << endl;
    cerr << argv[0] << ": " << "Missing files to do stream quality treatment to!\n";
    exit(1);
  }

  for(;optind<argc;++optind){
    MS_infiles.push_back(argv[optind]);
  }

  for(auto & fname : MS_infiles){
    uint8 ziptype=0;
    std::string ft;
    std::string dummyfromstem;
    std::string dummypathto;
    getCanonicalFileAndZipType(fname,dummypathto,dummyfromstem,ft,ziptype);
    if(ft.empty()) {
      cerr << "Unknown or illegal file extension '" << ft << "' in file name " << fname << "\n";
      exit(1);
    }
  }

  MS_tortype=Read::AS_ILLEGAL;
//  if(!MS_totype.empty()){
//    MS_tortype=checkToType(MS_totype);
//    if(MS_tortype==Read::AS_ILLEGAL){
//      cerr << "Unknown or illegal format '" << MS_totype << "' defined as to-type\n";
//      exit(1);
//    }
//  }

//  if(!MS_hitpath.empty()){
//    boost::filesystem::path finaldest(walkSymLinks(MS_hitpath));
//    bool setmerge=true;
//    if(boost::filesystem::exists(finaldest)){
//      if(boost::filesystem::is_directory(finaldest)) {
//	setmerge=false;
//	if(MS_hitpath.back()!='/') MS_hitpath+='/';
//      }
//    }
//    MS_mergeoutput=setmerge;
//    if(setmerge) fileRemove(finaldest.string(),false);
//    if(boost::filesystem::exists(finaldest)
//       && !boost::filesystem::is_directory(finaldest)) {
//      cout.flush();
//      cerr << "\n\nCould not remove file " << finaldest << "\nIs it writable?";
//      exit(10);
//     }
//  }

  MIRAParameters::setupStdMIRAParameters(MS_Pv);
  if(!miraparams.empty()){
    cout << "Parsing special MIRA parameters: " << miraparams << endl;
    MIRAParameters::parse(miraparams,MS_Pv,false);
    cout << "Ok.\n";
  }else{
    MIRAParameters::parse("--job=Solexa",MS_Pv,false);
  }

  MS_Pv[0].getNonConstDirectoryParams().dir_tmp=".";
  MIRAParameters::postParsingChanges(MS_Pv);

  cout << "Pv size: " << MS_Pv.size() << endl;

  MIRAParameters::dumpAllParams(MS_Pv, cout);

  try{
    // find out which size of hash we are going to work with
    uint32 sizeofhash=0;
    if(MS_hashstatfname.empty()){
      if(MS_basesperhash==0) MS_basesperhash=31;
      sizeofhash=HashStatistics<vhash64_t>::byteSizeOfHash(MS_basesperhash);
    }else{
      auto mhs=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(MS_hashstatfname);
      sizeofhash=mhs.sizeofhash;
      if(MS_basesperhash>0 && MS_basesperhash != mhs.basesperhash){
	MIRANOTIFY(Notify::FATAL,"-k specified on the command line (" << MS_basesperhash << ") is different than the kmer size saved (" << mhs.basesperhash << ") in the hash statistics file. This is treated as error, bailing out.");
      }
      MS_basesperhash=mhs.basesperhash;
      cout << "Size of kmers in file " << MS_hashstatfname << ": " << MS_basesperhash << endl;
    }

    if(MS_basesperhash>256){
      cout << "Sorry, the max. kmer size supported atm is 256.\n";
      exit(10);
    }

    // saving time (a lot) when adding comments while loading reads
    // just remember to clear out the stringcontainers from time to time
    multitag_t::MT_fastnew=true;

    if(sizeofhash==8){
      doTreatWithHS(msqt_hs64);
    }else if(sizeofhash==16){
      doTreatWithHS(msqt_hs128);
    }else if(sizeofhash==32){
      doTreatWithHS(msqt_hs256);
    }else if(sizeofhash==64){
      doTreatWithHS(msqt_hs512);
    }else{
      BUGIFTHROW(true,"sizeofhash == " << sizeofhash << " is rather unexpected.");
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
  }

  int retvalue=0;
  if(MS_signal_ctrlc){
    retvalue=1;
  }else{
    if(!MS_infiles.empty()){
      cout << "\nProcess finished.\n\n";

      cout << "Total number of sequences read: " << MS_numreadsread << endl;
    }
  }

  FUNCEND();
  return retvalue;
}
