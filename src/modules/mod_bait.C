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

#include "modules/mod_bait.H"
#include "util/fmttext.H"
#include "util/machineinfo.H"

#include "version.H"


using std::cout;
using std::cin;
using std::cerr;
using std::endl;


struct sigaction sigIntHandler;

std::vector<MIRAParameters> MiraBait::MB_Pv;

//std::string MiraBait::MB_fromtype;
//std::string MiraBait::MB_totype;
//std::string MiraBait::MB_baitfromtype;

uint8 MiraBait::MB_tortype;

std::list<std::string> MiraBait::MB_infiles;
std::list<std::string> MiraBait::MB_baitfiles;
std::string MiraBait::MB_hitpath;
std::string MiraBait::MB_misspath;
std::string MiraBait::MB_nameprefix;

std::string MiraBait::MB_tmpdirname{"."};
std::string MiraBait::MB_hashstatfname;
std::string MiraBait::MB_savehashstatfname;

bool   MiraBait::MB_dustfilter=false;
bool   MiraBait::MB_mustdeletetargetfiles=true;
bool   MiraBait::MB_wantbaithits=true;
bool   MiraBait::MB_wantbaitmiss=false;
bool   MiraBait::MB_fwdandrev=true;
bool   MiraBait::MB_mergeoutput=false;
bool   MiraBait::MB_changeseqcase=true;

char   MiraBait::MB_maskchar=0;

bool   MiraBait::MB_testflag=false;

int16  MiraBait::MB_dustperc=67;

int32  MiraBait::MB_numbaithits=1;
uint32 MiraBait::MB_optthreads=0;
int32  MiraBait::MB_optmbtouse=75;

std::list<char> MiraBait::MB_filepairinfo;    // p = 2 files, P = 1 file interleave
uint32 MiraBait::MB_basesperhash=0;
uint64 MiraBait::MB_baitpoolsize=0;
uint64 MiraBait::MB_numreadsread=0;
uint64 MiraBait::MB_numpairsbaited=0;
uint64 MiraBait::MB_numpairsmissed=0;
uint64 MiraBait::MB_numunpairedbaited=0;
uint64 MiraBait::MB_numunpairedmissed=0;

uint64 MiraBait::MB_numclippedreadsinload; // number of reads which have clips already when loaded (CAF/MAF)

std::list<MiraBait::wqueueunit_t> MiraBait::MB_workqueue;
MiraBait::files_t MiraBait::MB_files;

std::unordered_map<std::string,uint8> MiraBait::MB_fromtypemap = {
  {"fastq",Read::AS_FASTQ},
  {"fasta",Read::AS_FASTA},
  {"gb",Read::AS_GBF},
  {"gbf",Read::AS_GBF},
  {"gbk",Read::AS_GBF},
  {"gbff",Read::AS_GBF},
  {"caf",Read::AS_CAF},
  {"maf",Read::AS_MAF},
  {"exp",Read::AS_GAP4DA},
  {"phd",Read::AS_PHD},
};

std::unordered_map<std::string,uint8> MiraBait::MB_totypemap = {
  {"fastq",Read::AS_FASTQ},
  {"fasta",Read::AS_FASTA},
  {"caf",Read::AS_CAF},
  {"scaf",Read::AS_CAF},
  {"maf",Read::AS_MAF},
  {"txt",Read::AS_READNAME},
};

boost::mutex mutex_ctrlc;

//MIRATimer dbg_mt;
//auto dbg_baittime=dbg_mt.diff();

HashStatistics<vhash64_t> hs64;
HashStatistics<vhash128_t> hs128;
HashStatistics<vhash256_t> hs256;
HashStatistics<vhash512_t> hs512;

bool MiraBait::MB_signal_ctrlc=false;


MiraBait::~MiraBait()
{
}

// Do not use exit()!
// Do not throw an exeption! It sometimes works and the exeption is caught in the main program,
//  but sometimes crashes to command line with
//     terminate called after throwing an instance of Notify
// So, HashStatistics has own abort functionality and we'll exit the regular way
void MiraBait::ctrlCHandler(int s)
{
  FUNCSTART("void MiraBait::ctrlCHandler(int s)");
  cout << "Control-C was caught. Cleaning up, this may take a few seconds." << endl;
  MB_signal_ctrlc=true;
  hs64.abortAll();
  hs128.abortAll();
  hs256.abortAll();
  hs512.abortAll();
}


void MiraBait::usage()
{
  cout << "\nUsage: ";
  cout << "mirabait [options]"
    " {-b baitfile [-b ...] | -B file | -j joblibrary}"
    " {-p file_1 file_2 | -P file3}* [file4 ...]\n\n";

  cout << FmtText::wordWrap("MIRAbait: a 'grep' like tool for kmers up to 256 bp\n\nmirabait selects reads from a read collection which are partly similar or equal to sequences defined as target baits. Similarity is defined by finding a user-adjustable number of common k-mers (sequences of k consecutive bases) which are the same in the bait sequences and the screened sequences to be selected, either in forward or forward/reverse complement direction. Adding a DUST-like repeat filter for repeats up 4 bases is optional. \nWhen used on paired files, selects sequences where at least one mate matches.\n\n");
  cout << "Options:\n";
  cout << "\t-b file\t\tLoad bait sequences from file\n"
    "\t\t\t (multiple -b allowed)\n"
    "\t-B file\t\tLoad baits from kmer statistics file, not from sequence files.\n"
    "\t\t\t Only one -B allowed, cannot be combined with -b.\n"
    "\t\t\t (see -K for creating such a file)\n"
    "\t-j job\t\tSet options for predefined job from supplied MIRA library\n"
    "\t\t\t Currently available jobs:\n"
    "\t\t\t   rrna\tBait rRNA sequences";

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

  cout << "\n"
    "\t-k int\t\tkmer length of bait in bases (<=256, default=31)\n"
    "\t-n int\t\tIf >0: minimum number of k-mer baits needed (default=1)\n"
    "\t\t\tIf <=0: allowed number of missed kmers over sequence\n"
    "\t\t\t        length\n";

  cout << "\n"
    "\t-d\t\tDo not use kmers with microrepeats (DUST-like, see also -D)\n"
    "\t-D int\t\tSet length of microrepeats in kmers to discard from bait.\n"
    "\t\t\t int > 0 microrepeat len in percentage of kmer length.\n"
    "\t\t\t       E.g.: -k 17 -D 67 --> 11.39 bases --> 12 bases.\n"
    "\t\t\t int < 0 microrepeat len in bases.\n"
    "\t\t\t int != 0 implies -d, int=0 turns DUST filter off.\n"
    "\t-i\t\tSelects sequences that do not hit bait\n"
    "\t-I\t\tSelects sequences that hit and do not hit bait (to\n"
    "\t\t\t different files)\n"
    "\t-r\t\tNo checking of reverse complement direction\n"
    "\t-t\t\tNumber of threads to use (default=0 -> up to 4 CPU cores)\n";

  cout << "\nOptions for output definition:\n";
  cout << FmtText::wordWrap("Normally mirabait writes separate result files (named 'bait_match_*' and 'bait_miss_*') for each input to the current directory. For changing this behaviour and other relating to output, use these options:\n");
  cout << "\t-c char\t\tNormally, mirabait lowercases bases a kmer hit.\n"
    "\t\t\t Using this option, one can instead mask those bases with the\n"
    "\t\t\t given character.\n"
    "\t\t\t Use a blank to neither mask nor lowercase hits.\n";
  cout << "\t-l int\t\tlength of a line (FASTA only, default 0=unlimited)\n";
  cout << "\t-K file\t\tSave kmer statistics to 'file' (see also -B)\n";
  cout << "\t-N name\t\tChange the prefix 'bait' to <name>\n"
    "\t\t\t Has no effect if -o/-O is used and targets are not\n"
    "\t\t\t directories\n";
  cout << "\t-o <path>\tSave sequences matching bait to path\n"
    "\t\t\t If path is a directory, write separate files into this\n"
    "\t\t\t directory. If not, combine all matching sequences from\n"
    "\t\t\t the input file(s) into a single file specified by the\n"
    "\t\t\t path.\n";
  cout << "\t-O <path>\tLike -o, but for sequences not matching\n";

  cout << "\nOther options:\n";
  cout << "\t-T dir\t\tUse 'dir' as directory for temporary files instead of\n"
    "\t\t\t current working directory.\n";
  cout << "\t-m integer\tMemory to use for computing kmer statistics\n"
    "\t\t\t 0..100 = use percentage of free system memory\n"
    "\t\t\t >100 = amount of MiB to use (e.g. 16384 for 16 GiB)\n"
    "\t\t\t Default 75 (75% of free system memory).\n";


  cout << "\nDefining files types to load/save:\n";
  cout << FmtText::wordWrap("Normally mirabait recognises the file types according to the file extension (even when packed). In cases you need to force a certain file type because the file extension is non-standard, use the EMBOSS notation to force a type: <filetype>::<name_of_file>. E.g., to tell that \"somefile.dat\" is FASTQ, use: fastq::somefile.dat\nRecognised types are: caf, fasta, fastq, gbf, gbk, gbff, maf and phd.\n\nMIRABAIT will write files in the same file type as the corresponding input files.");

  cout << "\n\nExamples:\n"
    "  mirabait -b b.fasta file.fastq"
    "\n\n  mirabait -I -j rrna -p file_1.fastq file_2.fastq"
    "\n\n  mirabait -b b1.fasta -b b2.gbk file.fastq"
    "\n\n  mirabait -b fasta::baits.dat -p fastq::file_1.dat fastq::file_2.dat"
    "\n\n  mirabait -b b.fasta -p file_1.fastq file_2.fastq -P file3.fasta file4.caf"
    "\n\n  mirabait -I -b b.fasta -p file_1.fastq file_2.fastq -P file3.fasta file4.caf"
    "\n\n  mirabait -k 27 -n 10 -b b.fasta file.fastq"
    "\n\n  mirabait -b fasta::b.dat fastq::file.dat"
    "\n\n  mirabait -o /dev/shm/ -b b.fasta -p file_1.fastq file_2.fastq"
    "\n\n  mirabait -o /dev/shm/match -b b.fasta -p file_1.fastq file_2.fastq"
    "\n\n  mirabait -b human_genome.fasta -K HG_kmerstats.mhs.gz -p file1.fastq file2.fastq"
    "\n\n  mirabait -B HG_kmerstats.mhs.gz -p file1.fastq file2.fastq"
    "\n\n  mirabait -d -B HG_kmerstats.mhs.gz -p file1.fastq file2.fastq"
    "\n";
}


uint8 MiraBait::checkFromType(std::string & fromtype)
{
  auto mI=MB_fromtypemap.find(fromtype);
  if(mI!=MB_fromtypemap.end()) return mI->second;
  return Read::AS_ILLEGAL;
}

uint8 MiraBait::checkToType(std::string & totype)
{
  auto mI=MB_totypemap.find(totype);
  if(mI!=MB_totypemap.end()) return mI->second;
  return Read::AS_ILLEGAL;
}



uint8 MiraBait::setupRPIO(std::string filename, ReadGroupLib::ReadGroupID rgid, ReadPoolIO & rpio, uint8 & ziptype)
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

void MiraBait::setupOutfiles(const std::string & fname, uint8 rtype, uint8 ziptype, std::ofstream & hitfout, std::ofstream & missfout)
{
  if(hitfout.is_open()) hitfout.close();
  if(missfout.is_open()) missfout.close();

  std::string pname(MB_nameprefix);
  if(pname.empty()) pname="bait";
  std::string matchpre(pname+"_match_");
  std::string misspre(pname+"_miss_");

  boost::filesystem::path fp(fname);
  // if we work on zipped files, the ouput must get away the zip extension
  if(ziptype){
    fp=fp.parent_path() / fp.stem();  // boost::filesystem append operation
  }

  std::string hitname;
  if(MB_mergeoutput && !MB_hitpath.empty()){
    hitname=MB_hitpath;
  }else{
    hitname=MB_hitpath;
    if(!hitname.empty()) hitname+='/';
    hitname+=matchpre+fp.filename().string();
  }
  std::string missname;
  if(MB_mergeoutput && !MB_misspath.empty()){
    missname=MB_misspath;
  }else{
    missname=MB_misspath;
    if(!missname.empty()) missname+='/';
    missname+=misspre+fp.filename().string();
  }

  auto omode=std::ios::out;
  if(MB_mergeoutput) omode|=std::ios::app;
  if(MB_wantbaithits){
    hitfout.open(hitname,omode);
    if(!hitfout.is_open()){
      cout.flush();
      cerr << "\n\nCould not open " << hitname << ".\nDoes the path exist, is it writable? Is the disk full?";
      exit(10);
    }
  }
  if(MB_wantbaitmiss){
    missfout.open(missname,omode);
    if(!missfout.is_open()){
      cout.flush();
      cerr << "\n\nCould not open " << missname << ".\nDoes the path exist, is it writable? Is the disk full?";
      exit(10);
    }
  }

  if(MB_wantbaithits && MB_wantbaitmiss){
    cout << ", sorting" << endl;
  }else{
    cout << ", filtering" << endl;
  }
  if(MB_wantbaithits){
    cout << "+++ matches to " << hitname << endl;
  }
  if(MB_wantbaitmiss){
    cout << "--- non-matches to " << missname << endl;
  }
}


template<typename TVHASH_T>
void MiraBait::baitReads(HashStatistics<TVHASH_T> & hs, const ReadPool & rp, std::vector<uint8> & take)
{
  //dbg_mt.reset();

  take.clear();
  take.resize(rp.size(),0);
  for(uint32 rpi=0; rpi<take.size(); ++rpi){
    if(rp[rpi].getLenClippedSeq() != rp[rpi].getLenSeq()){
      ++MB_numclippedreadsinload;
    }
    uint32 neededhashes=static_cast<uint32>(MB_numbaithits);
    if(MB_numbaithits<=0){
      auto tmpnh=
	static_cast<int32>(rp[rpi].getLenClippedSeq())-static_cast<int32>(MB_basesperhash-1)+MB_numbaithits;
      if(tmpnh<0) tmpnh=0;
      neededhashes=static_cast<uint32>(tmpnh);
    }
    if(hs.checkBaitHit(rp[rpi],MB_changeseqcase,MB_maskchar) >= neededhashes){
      take[rpi]=1;
    }
  }

  //dbg_baittime+=dbg_mt.diff();
}


template<typename TVHASH_T>
void MiraBait::parallelBaitReads(HashStatistics<TVHASH_T> & hs, const ReadPool & rp, std::vector<uint8> & take)
{
  take.clear();
  take.resize(rp.size(),0);
  for(uint32 rpi=0; rpi<take.size(); ++rpi){
    if(rp[rpi].getLenClippedSeq() != rp[rpi].getLenSeq()){
      ++MB_numclippedreadsinload;
    }
  }

#pragma omp parallel
  {
#pragma omp for
    for(uint32 rpi=0; rpi<take.size(); ++rpi){
      uint32 neededhashes=static_cast<uint32>(MB_numbaithits);
      if(MB_numbaithits<=0){
	auto tmpnh=
	  static_cast<int32>(rp[rpi].getLenClippedSeq())-static_cast<int32>(MB_basesperhash-1)+MB_numbaithits;
	if(tmpnh<0) tmpnh=0;
	neededhashes=static_cast<uint32>(tmpnh);
      }
      if(hs.checkBaitHit(rp[rpi],MB_changeseqcase,MB_maskchar) >= neededhashes){
	take[rpi]=1;
      }
    }
  }
}


void MiraBait::saveWQueueElement(wqueueunit_t & wqu)
{
  FUNCSTART("void MiraBait::saveWQueueElement(wqueueunit_t & wqu)");

  Read::setCoutType(MB_files.writetype);

  uint64 numbaited=0;
  uint64 numnonbaited=0;

  auto t1I=wqu.take1.begin();
  for(uint32 rpi=0; rpi<wqu.rp1.size(); ++rpi, ++t1I){
    if(*t1I && MB_wantbaithits) {
      MB_files.hitfout1 << wqu.rp1[rpi];
      ++numbaited;
    }else if(!*t1I && MB_wantbaitmiss){
      MB_files.missfout1 << wqu.rp1[rpi];
      ++numnonbaited;
    }
    if(wqu.pairstatus==PS_2FILES){
      if(*t1I && MB_wantbaithits) {
	if(MB_mergeoutput){
	  MB_files.hitfout1 << wqu.rp2[rpi];
	}else{
	  MB_files.hitfout2 << wqu.rp2[rpi];
	}
      }else if(!*t1I && MB_wantbaitmiss){
	if(MB_mergeoutput){
	  MB_files.missfout1 << wqu.rp2[rpi];
	}else{
	  MB_files.missfout2 << wqu.rp2[rpi];
	}
      }
    }
  }

  MB_numreadsread+=wqu.rp1.size()+wqu.rp2.size();

  if(wqu.pairstatus==PS_NOPAIR){
    MB_numunpairedbaited+=numbaited;
    MB_numunpairedmissed+=numnonbaited;
  }else if(wqu.pairstatus==PS_INTERLEAVE){
    MB_numpairsbaited+=numbaited/2;
    MB_numpairsmissed+=numnonbaited/2;
  }else if(wqu.pairstatus==PS_2FILES){
    MB_numpairsbaited+=numbaited;
    MB_numpairsmissed+=numnonbaited;
  }else{
    MIRANOTIFY(Notify::INTERNAL,"Ummm ... unknown pairstatus " << static_cast<uint16>(wqu.pairstatus));
  }

//  MB_numreadsread+=rp.size();
//  MB_numreadswritten+=rp.size();
}


template<typename TVHASH_T>
void MiraBait::doBaitWithHS(HashStatistics<TVHASH_T> & mbhs)
{
  FUNCSTART("void MiraBait::doBaitWithHS(HashStatistics<TVHASH_T> & mbhs)");
  {
    // setup the CTRL-C handler only NOW
    sigIntHandler.sa_handler = MiraBait::ctrlCHandler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    if(!MB_hashstatfname.empty()){
      cout << "Loading from existing hashstat file ... "; cout.flush();
      mbhs.loadHashStatistics(MB_hashstatfname);
      cout << "done.\n";
      if(MB_basesperhash!=0){
	if(mbhs.getBasesPerHash()!=MB_basesperhash){
	  cout << "Error: kmer size set for mirabait (" << MB_basesperhash
	       << ") is not equal to the kmer size loaded from file ("
	       << mbhs.getBasesPerHash() << ")!\nDid you know you can leave out -k when using -B?\nAborting!\n";
	  exit(10);
	}
      }else{
	cout << "No -k given, using k from the loaded file: "
	     << mbhs.getBasesPerHash() << endl;
	MB_basesperhash=mbhs.getBasesPerHash();
      }
    }else{
      mbhs.computeHashStatistics(
	MB_baitfiles,
	MB_optmbtouse, // 75% of free memory as buffers
	true, //fwdandrev
	1, //fwdrevmin
	0, // rare kmer early kill
	MB_basesperhash,
	MB_savehashstatfname,MB_tmpdirname);

      if(MB_signal_ctrlc) return;
      if(mbhs.getNumHashEntries()==0){
	cout << FmtText::makeTextSign("WARNING: not a single kmer bait could be generated. This is due to the sequences you are using to bait are all either too short or contain too many closely located IUPAC codes.\nThis may be right, but most probably is not. If not: either check your bait sequences in the input files or choose a lower kmer size.") << endl;
      }
    }
  }

  if(MB_signal_ctrlc) return;

  mbhs.showHashStatisticsInfo();
  if(MB_dustfilter){
    if(MB_dustperc>=0){
      mbhs.removeDustEndsRatio(0,MB_dustperc);
    }else{
      mbhs.removeDustEndsFixed(0,-MB_dustperc);
    }
    mbhs.showHashStatisticsInfo();
    if(!MB_savehashstatfname.empty()){
      cout << "KMER statistics dusted and wants save ... resaving now.\n";
      mbhs.saveHashStatistics(MB_savehashstatfname);
    }
  }

  //dbg_baittime=0;

  // in this part I've never seen faster execution than with 4, on the contrary
  // to be monitored
#ifdef HAVE_OPENMP
  if(MB_optthreads>4) omp_set_num_threads(4);
#endif

  MB_workqueue.resize(1);
  auto qI=MB_workqueue.begin(); // fixed atm

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

  auto ifI=MB_infiles.begin();
  while(ifI!=MB_infiles.end()){
    if(MB_signal_ctrlc) return;
    uint8 ziptype=0;
    rpio1.setAttributeProgressIndicator(true);
    MB_files.infilename1.clear();
    MB_files.infilename2.clear();

    MB_files.intype1=setupRPIO(*ifI,rgid,rpio1,ziptype);
    MB_files.writetype=MB_files.intype1;
    if(MB_mergeoutput){
      if(ifI==MB_infiles.begin() && MB_tortype==Read::AS_ILLEGAL){
	MB_tortype=MB_files.intype1;
      }
    }
    if(MB_tortype!=Read::AS_ILLEGAL){
      MB_files.writetype=MB_tortype;
    }
    setupOutfiles(*ifI,MB_files.intype1,ziptype,MB_files.hitfout1,MB_files.missfout1);
    MB_files.infilename1=*ifI;

    qI->pairstatus=PS_NOPAIR;
    if(!MB_filepairinfo.empty()){
      if(MB_filepairinfo.front()=='P'){
	qI->pairstatus=PS_INTERLEAVE;
      }else if(MB_filepairinfo.front()=='p'){
	++ifI;
	if(ifI==MB_infiles.end()){
	  MIRANOTIFY(Notify::FATAL,"Something's wrong here: -p says to expect one further file, but " << *(--ifI) << " is the last file seen on the command line?");
	}
	rpio1.setAttributeProgressIndicator(false);
	rpio2.setAttributeProgressIndicator(true);
	MB_files.intype2=setupRPIO(*ifI,rgid,rpio2,ziptype);
	setupOutfiles(*ifI,MB_files.intype2,ziptype,MB_files.hitfout2,MB_files.missfout2);
	MB_files.infilename2=*ifI;
	qI->pairstatus=PS_2FILES;
      }
      MB_filepairinfo.pop_front();
    }

    ++ifI;

    while(rpio1.loadNextSeqs(500)){
      if(MB_signal_ctrlc) return;
      if(qI->pairstatus==PS_2FILES){
	rpio2.loadNextSeqs(500);
	if(qI->rp1.size() != qI->rp2.size()){
	  MIRANOTIFY(Notify::FATAL,"Something's wrong here: -p says you have two files with reads paired across both files. But file " << MB_files.infilename1 << " does not have the same number of reads as file " << MB_files.infilename2 << " ???");
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

      parallelBaitReads(mbhs,qI->rp1,qI->take1);
      if(qI->pairstatus==PS_2FILES){
	parallelBaitReads(mbhs,qI->rp2,qI->take2);

	// take both reads in dual load
	auto t1I=qI->take1.begin();
	auto t2I=qI->take2.begin();
	for(; t1I != qI->take1.end(); ++t1I, ++t2I){
	  if(*t2I) *t1I=1;
	  if(*t1I) *t2I=1;
	}
      }

      if(qI->pairstatus==PS_INTERLEAVE){
	// take both reads in interleaved
	auto sI=qI->take1.begin();
	auto eI=sI+1;
	for(; sI!=qI->take1.end(); sI+=2, eI+=2){
	  if(*sI | *eI){
	    *sI=1;
	    *eI=1;
	  }
	}
      }

      saveWQueueElement(*qI);

      Read::trashReadNameContainer();
      multitag_t::trashContainers();  // saving memory (and time)
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


int MiraBait::mainMiraBait(int argc, char ** argv)
{
  //CALLGRIND_STOP_INSTRUMENTATION;

  FUNCSTART("int mainMiraBait(int argc, char ** argv)");

  int c;
  extern char *optarg;
  extern int optind;


  std::string path;
  std::string convertprog;
  splitFullPathAndFileName(argv[0],path,convertprog);

  std::string miraparams;

  MB_basesperhash=0;
  MB_baitpoolsize=0;

  Read::setCoutLen(0);

  MB_filepairinfo.clear();    // p = 2 files, P = 1 file interleave

  while (true){
    static struct option long_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"version", no_argument,         0, 'V'},
	{"throw", no_argument,         0, '{'},
	{"loadhsf", no_argument,    0, 'L'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "{hdiIpPrva:b:c:B:D:j:k:K:l:m:n:N:o:O:t:T:",
			 long_options, &option_index);
    if(c == -1) break;

    switch (c) {
    case '{': {
      Notify::setBangOnThrow(true);
      break;
    }
    case '1': {
      MB_testflag=true;
      break;
    }
    case 'a': {
      miraparams=optarg;
      break;
    }
    case 'b': {
      MB_baitfiles.push_back(optarg);
      break;
    }
    case 'B': {
      MB_hashstatfname=optarg;
      break;
    }
    case 'c': {
      std::string tmpc=optarg;
      if(tmpc.size()!=1){
	usage();
	cout << endl;
	cerr << "ERROR: -c must be a single character\n";
	exit(1);
      }
      MB_changeseqcase=false;
      if(tmpc[0]!=' '){
	MB_maskchar=tmpc[0];
      }
      break;
    }
    case 'd': {
      MB_dustfilter=true;
      break;
    }
    case 'D': {
      MB_dustfilter=true;
      int64 bla=atoi(optarg);
      if(bla>100
	 || bla<-1000){
	usage();
	cerr << "-D must be >-1000 and <= 100 in the posi\n";
	exit(2);
      }
      if(bla==0) {
	MB_dustfilter=false;
      }else{
	MB_dustperc=bla;
      }
      break;
    }
    case 'i': {
      MB_wantbaithits=false;
      MB_wantbaitmiss=true;
      break;
    }
    case 'I': {
      MB_wantbaithits=true;
      MB_wantbaitmiss=true;
      break;
    }
    case 'j': {
      std::string job(optarg);
      boost::to_lower(job);
      if(job=="rrna"){
	MB_hashstatfname=MIRAParameters::getMHSLibDir()+"/filter_default_rrna.mhs.gz";
	MB_numbaithits=20;
	MB_dustfilter=true;
	MB_nameprefix="rRNA";
	cout << "Found -j rrna. Equivalent settings:\n"
	     << "-B " << MB_hashstatfname
	     << " -n " << MB_numbaithits
	     << " -d"
	     << " -N " << MB_nameprefix << "\n\n";
	if(!fileExists(MB_hashstatfname)){
	  MIRANOTIFY(Notify::INSTALL, "You specified -j rrna, but MIRA could not find the default rRNA hash statistics file which should be located at " + MB_hashstatfname + "\nThis is usually the case when you forgot to run the script to install the rRNA data from the MIRA package.\n\nWhat to do:\n- if you installed from source: use the 'make install' functionality.\n- if you installed from a precompiled binaries package: go to the 'dbinstall' directory in that package and type\n  ./mira-install-sls-rrna.sh rfam_rrna-21-12.sls.gz\n");
	}
      }else{
	usage();
	cerr << "\nUnknown job -j " << job << "\nCurrently known jobs: 'rrna'" << endl;
	exit(2);
      }
      break;
    }
    case 'k': {
      uint32 bla=atoi(optarg);
      MB_basesperhash=bla;
      break;
    }
    case 'K': {
      MB_savehashstatfname=optarg;
      break;
    }
    case 'l': {
      Read::setCoutLen(atoi(optarg));
      break;
    }
    case 'm': {
      MB_optmbtouse=atoi(optarg);
      break;
    }
    case 'n': {
      MB_numbaithits=atoi(optarg);
      break;
    }
    case 'N': {
      MB_nameprefix=optarg;
      break;
    }
    case 'o': {
      MB_wantbaithits=true;
      MB_hitpath=optarg;
      break;
    }
    case 'O': {
      MB_wantbaitmiss=true;
      MB_misspath=optarg;
      break;
    }
    case 'p' :
    case 'P' : {
      MB_filepairinfo.push_back(c);
      break;
    }
    case 'r': {
      MB_fwdandrev=false;
      break;
    }
    case 't': {
      MB_optthreads=atoi(optarg);
      break;
    }
    case 'T': {
      MB_tmpdirname=optarg;
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
    case 'V':
      cout << "MIRABAIT " << miraversion
	   << "\n\nWritten by Bastien Chevreux\n";
      exit(0);
    default : {}
    }
  }

  if(MB_optthreads==0){
    MB_optthreads=MachineInfo::getCoresTotal();
  }
#ifdef HAVE_OPENMP
  omp_set_num_threads(MB_optthreads);
#endif

  if(MB_baitfiles.empty() && MB_hashstatfname.empty()){
    usage();
    cout << endl;
    cerr << argv[0] << ": " << "No bait files defined via -b and no -B given!\nDid you use the command line for the old mirabait (<= 4.0.2)?\n";
    exit(1);
  }

  if(argc-optind < 1){
    if(!MB_savehashstatfname.empty()) {
      cout << "-K given and no sequence files, will just create a kmer statistics file.\n";
    }else{
      usage();
      cout << endl;
      cerr << argv[0] << ": " << "Missing files to bait sequences from!\n";
      exit(1);
    }
  }

  for(;optind<argc;++optind){
    MB_infiles.push_back(argv[optind]);
  }

  for(auto & fname : MB_infiles){
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

  MB_tortype=Read::AS_ILLEGAL;
//  if(!MB_totype.empty()){
//    MB_tortype=checkToType(MB_totype);
//    if(MB_tortype==Read::AS_ILLEGAL){
//      cerr << "Unknown or illegal format '" << MB_totype << "' defined as to-type\n";
//      exit(1);
//    }
//  }

  if(!MB_hitpath.empty()){
    boost::filesystem::path finaldest(walkSymLinks(MB_hitpath));
    bool setmerge=true;
    if(boost::filesystem::exists(finaldest)){
      if(boost::filesystem::is_directory(finaldest)) {
	setmerge=false;
	if(MB_hitpath.back()!='/') MB_hitpath+='/';
      }
    }
    MB_mergeoutput=setmerge;
    if(setmerge) fileRemove(finaldest.string(),false);
    if(boost::filesystem::exists(finaldest)
       && !boost::filesystem::is_directory(finaldest)) {
      cout.flush();
      cerr << "\n\nCould not remove file " << finaldest << "\nIs it writable?";
      exit(10);
     }
  }
  if(!MB_misspath.empty()){
    boost::filesystem::path finaldest(walkSymLinks(MB_misspath));
    bool setmerge=true;
    if(boost::filesystem::exists(finaldest)){
      if(boost::filesystem::is_directory(finaldest)) {
	setmerge=false;
	if(MB_misspath.back()!='/') MB_misspath+='/';
      }
    }
    MB_mergeoutput=setmerge;
    if(setmerge) fileRemove(finaldest.string(),false);
    if(boost::filesystem::exists(finaldest)){
      cout.flush();
      cerr << "\n\nCould not remove file " << finaldest << "\nIs it writable?";
      exit(10);
     }
  }

  MIRAParameters::setupStdMIRAParameters(MB_Pv);
  if(!miraparams.empty()){
    cout << "Parsing special MIRA parameters: " << miraparams << endl;
    MIRAParameters::parse(miraparams,MB_Pv,false);
    cout << "Ok.\n";
  }

  if(MB_numbaithits>0){
    cout << "Baiting sequences with at least " << MB_numbaithits << " exact kmer matches.\n";
  }else{
    cout << "Baiting sequences allowing for " << -MB_numbaithits << " missed kmer matches over the sequence length.\n";
  }

  try{
    // find out which size of hash we are going to work with
    uint32 sizeofhash=0;
    if(MB_hashstatfname.empty()){
      if(MB_basesperhash==0) MB_basesperhash=31;
      sizeofhash=HashStatistics<vhash64_t>::byteSizeOfHash(MB_basesperhash);
    }else{
      auto mhs=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(MB_hashstatfname);
      sizeofhash=mhs.sizeofhash;
      if(MB_basesperhash>0 && MB_basesperhash != mhs.basesperhash){
	MIRANOTIFY(Notify::FATAL,"-k specified on the command line (" << MB_basesperhash << ") is different than the kmer size saved (" << mhs.basesperhash << ") in the hash statistics file. This is treated as error, bailing out.");
      }
      MB_basesperhash=mhs.basesperhash;
      cout << "Size of kmers in file " << MB_hashstatfname << ": " << MB_basesperhash << endl;
    }

    if(MB_basesperhash>256){
      cout << "Sorry, the max. kmer size supported atm is 256.\n";
      exit(10);
    }

    // saving time (a lot) when adding comments while loading reads
    // just remember to clear out the stringcontainers from time to time
    multitag_t::MT_fastnew=true;

    if(sizeofhash==8){
      doBaitWithHS(hs64);
    }else if(sizeofhash==16){
      doBaitWithHS(hs128);
    }else if(sizeofhash==32){
      doBaitWithHS(hs256);
    }else if(sizeofhash==64){
      doBaitWithHS(hs512);
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
  catch (const std::exception& e) { // reference to the base of a polymorphic object
    std::cerr << "Std exception caught. Message: " << e.what() << std::endl; // information from length_error printed
    exit(10);
  }
  catch(...){
    cerr << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    exit(10);
  }

  int retvalue=0;
  if(MB_signal_ctrlc){
    retvalue=1;
  }else{
    if(!MB_infiles.empty()){
      if(MB_numclippedreadsinload){
	cout << FmtText::wordWrap("\nNOTICE! You baited sequences which had clipping information (CAF or MAF). Mirabait will have baited *only* in the unclipped parts of the sequences (which are thought to represent 'good, viable' sequence).\n");
      }

      cout << "\nBaiting process finished.\n\n";
      cout << "Number of bait sequences:   " << MB_baitpoolsize << endl;
      cout << "Total number of sequences read: " << MB_numreadsread << endl;
      if (MB_numreadsread) {
	cout << "Pairs baited: " << MB_numpairsbaited << " (" << std::fixed << std::setprecision(2) << 100.0f/MB_numreadsread*MB_numpairsbaited*2 << "%)\n";
	cout << "Pairs missed: " << MB_numpairsmissed << " (" << std::fixed << std::setprecision(2) << 100.0f/MB_numreadsread*MB_numpairsmissed*2 << "%)\n";
	cout << "Unpaired baited: " << MB_numunpairedbaited << " (" << std::fixed << std::setprecision(2) << 100.0f/MB_numreadsread*MB_numunpairedbaited << "%)\n";
	cout << "Unpaired missed: " << MB_numunpairedmissed << " (" << std::fixed << std::setprecision(2) << 100.0f/MB_numreadsread*MB_numunpairedmissed << "%)\n";
      }
    }

    if(!MB_savehashstatfname.empty()){
      cout << "Kmer statistics file saved to " << MB_savehashstatfname << " and is ready to reuse via -B.\n";
    }
  }

  FUNCEND();
  return retvalue;
}
