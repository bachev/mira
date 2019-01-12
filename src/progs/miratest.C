#include <typeinfo>
#include <iostream>
#include <string>

#include <sys/times.h>
#include <limits.h>
#include <unistd.h>

#include "stdinc/types.H"
#include "stdinc/defines.H"

#include "io/fastq-mira.H"
#include "io/fastq-lh.H"

#include "util/misc.H"
#include "util/machineinfo.H"

#include "errorhandling/errorhandling.H"


#include <boost/filesystem.hpp>

KSEQ_INIT(gzFile, gzread)

using namespace std;


#include "mira/assembly.H"
#include "mira/readpool_io.H"
#include "mira/skim.H"
#include "mira/seqtohash.H"
#include "util/dptools.H"
#include "mira/hashstats.H"



class A {
  class B {
    ~B();
  };
};

A::B::~B() {
  // ...
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string revComp(string & s)
{
  string ret;
  for(auto rI=s.rbegin(); rI!=s.rend(); ++rI){
    if(*rI=='A'){
      ret.push_back('T');
    }else if(*rI=='C'){
      ret.push_back('G');
    }else if(*rI=='G'){
      ret.push_back('C');
    }else{
      ret.push_back('A');
    }
  }
  return ret;
}

uint16 str2hash(string & s)
{
  uint8 ret=0;
  for(auto rI=s.begin(); rI!=s.end(); ++rI){
    ret<<=2;
    if(*rI=='A'){
    }else if(*rI=='C'){
      ret+=1;
    }else if(*rI=='G'){
      ret+=2;
    }else{
      ret+=3;
    }
  }
  return ret;
}

/**************************************************************************************************
 **************************************************************************************************/

template<class TVHASH_T>
void ddbgCorrectReadsInReadpool(ReadPool & rp, HashStatistics<TVHASH_T> & hs)
{
  vector<typename HashStatistics<TVHASH_T>::dbgedits_t> edits;
  for(uint32 rpi=0; rpi<rp.size(); ++rpi){
    Read & actread=rp[rpi];
    hs.proposeSDBGEditsForSequence(actread.getClippedSeqAsChar(),
				   actread.getLenClippedSeq(),
				   actread.getName().c_str(),
				   edits);
    if(!edits.empty()){
      cout << "Have " << edits.size() << " edits for " << actread.getName() << endl;
      for(auto erI=edits.rbegin(); erI!=edits.rend(); ++erI){
	string olds(static_cast<string>(actread.getClippedSeqAsChar()).substr(erI->pos,erI->len));
	cout << "Could correct " << erI-edits.rbegin() << " " << actread.getName() << " " << erI->pos << " " << erI->pos+erI->len << " " << erI->len;
	if(olds.size() != erI->replacement.size()) cout << " LENDIFF";
	cout << "\n" << olds
	     << "\n" << erI->replacement << endl;
	//Read::setCoutType(Read::AS_TEXT);
	//cout << "OLD:\n" << actread << endl;
	actread.smoothSequenceReplace(erI->pos+actread.getLeftClipoff(),erI->len,erI->replacement);
	//cout << "NEW:\n" << actread << endl;
      }
    }
  }
}


template<class TVHASH_T>
void doCorrect(uint32 trimfreq, const string & merfile, HashStatistics<TVHASH_T> & hs)
{
  hs.loadHashStatistics(merfile);

  hs.trimHashStatsByFrequencyAND(trimfreq,trimfreq,0);
  hs.calcKMerForks(trimfreq,true);
  hs.buildSDBGraphs();

  ReadPool rp1;
  ReadPoolIO rpio1(rp1);
  rpio1.setAttributeFASTQQualOffset(33); // in case we load FASTQs
  ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_SOLEXA);

  string fn("reads2correct.fastq");
  //string fn("NG-5413_03_SP1.fastq");

  ofstream fout("corrected.fastq");

  rpio1.setAttributeProgressIndicator(true);
  rpio1.registerFile("fastq", fn, "", rgid, false);

  while(rpio1.loadNextSeqs(1)){
    cout << rp1[0].getName() << endl;
    Assembly::findDegeneratePolymerases(rp1);
    //rp1[0].performHardTrim();
    ddbgCorrectReadsInReadpool(rp1,hs);
    rp1.dumpAs(fout,Read::AS_FASTQ,true);
    rp1.discard();
  }
}


/**************************************************************************************************
 **************************************************************************************************/

template<class TVHASH_T>
void hcHuntInPool(ReadPool & rp, HashStatistics<TVHASH_T> & hs)
{
  for(uint32 rpi=0; rpi<rp.size(); ++rpi){
    Read & actread=rp[rpi];
    auto res=hs.checkSequenceForSDBGChimeras(actread.getClippedSeqAsChar(),
					     actread.getLenClippedSeq(),
					     actread.getName().c_str());
    if(res){
      cout << "CHIMERA! " << actread.getName() << endl;
    }
  }
}


template<class TVHASH_T>
void huntChimeras(uint32 trimfreq, const string & merfile, HashStatistics<TVHASH_T> & hs)
{
  hs.loadHashStatistics(merfile);

  hs.trimHashStatsByFrequencyAND(trimfreq,trimfreq,0);
  hs.calcKMerForks(trimfreq,true);
  hs.buildSDBGraphs();

  ReadPool rp1;
  ReadPoolIO rpio1(rp1);
  rpio1.setAttributeFASTQQualOffset(33); // in case we load FASTQs
  ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_SOLEXA);

  string fn("reads2correct.fastq");
  //string fn("NG-5413_03_SP1.fastq");

  ofstream fout("corrected.fastq");

  rpio1.setAttributeProgressIndicator(true);
  rpio1.registerFile("fastq", fn, "", rgid, false);

  while(rpio1.loadNextSeqs(1)){
    cout << rp1[0].getName() << endl;
    hcHuntInPool(rp1,hs);
    rp1.dumpAs(fout,Read::AS_FASTQ,true);
    rp1.discard();
  }
}


/**************************************************************************************************
 **************************************************************************************************/

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>

#include "util/timer.H"

// the function f() does some time-consuming work
void fH()
{
  volatile double d = 0;

  HRTimer ttime;
  HRTimer mtime;
  auto mdiff=mtime.diff();
  for(int n=0; n<10000000; ++n){
    mtime.reset();
    for(int m=0; m<10; ++m){
      d += d*n*m;
    }
    mdiff+=mtime.diff();
  }
  auto tdiff=ttime.diff();

  cout << "HR Total mtime " << HRTimer::toMicro(mdiff) << endl;
  cout << "HR Total ttime " << HRTimer::toMicro(tdiff) << endl;
  cout << "HR Diffs " << HRTimer::toMicro(tdiff)-HRTimer::toMicro(mdiff) << endl;
}

void fL()
{
  volatile double d = 0;

  LRTimer ttime;
  LRTimer mtime;
  auto mdiff=mtime.diff();
  for(int n=0; n<10000000; ++n){
    mtime.reset();
    for(int m=0; m<10; ++m){
      d += d*n*m;
    }
    mdiff+=mtime.diff();
  }
  auto tdiff=ttime.diff();

  cout << "LR Total mtime " << LRTimer::toMicro(mdiff) << endl;
  cout << "LR Total ttime " << LRTimer::toMicro(tdiff) << endl;
  cout << "LR Diffs " << LRTimer::toMicro(tdiff)-LRTimer::toMicro(mdiff) << endl;
}

void fC()
{
  volatile double d = 0;

  CLRTimer ttime;
  CLRTimer mtime;
  auto mdiff=mtime.diff();
  for(int n=0; n<10000000; ++n){
    mtime.reset();
    for(int m=0; m<10; ++m){
      d += d*n*m;
    }
    mdiff+=mtime.diff();
  }
  auto tdiff=ttime.diff();

  cout << "CLR Total mtime " << CLRTimer::toMicro(mdiff) << endl;
  cout << "CLR Total ttime " << CLRTimer::toMicro(tdiff) << endl;
  cout << "CLR Diffs " << CLRTimer::toMicro(tdiff)-CLRTimer::toMicro(mdiff) << endl;
}

void fM()
{
  volatile double d = 0;

  MIRATimer ttime;
  MIRATimer mtime;
  auto mdiff=mtime.diff();
  for(int n=0; n<10000000; ++n){
    mtime.reset();
    for(int m=0; m<10; ++m){
      d += d*n*m;
    }
    mdiff+=mtime.diff();
  }
  auto tdiff=ttime.diff();

  cout << "CM Total mtime " << MIRATimer::toMicro(mdiff) << endl;
  cout << "CM Total ttime " << MIRATimer::toMicro(tdiff) << endl;
  cout << "CM Diffs " << MIRATimer::toMicro(tdiff)-MIRATimer::toMicro(mdiff) << endl;
}

void fS()
{
  volatile double d = 0;

  timeval mstart;
  timeval tstart;

  uint64 mdiff=0;

  gettimeofday(&tstart,nullptr);
  for(int n=0; n<10000000; ++n){
    gettimeofday(&mstart,nullptr);
    for(int m=0; m<10; ++m){
      d += d*n*m;
    }
    mdiff+=diffsuseconds(mstart);
  }
  auto tdiff=diffsuseconds(tstart);

  cout << "SR Total mtime " << mdiff << endl;
  cout << "SR Total ttime " << tdiff << endl;
  cout << "SR Diffs " << tdiff-mdiff << endl;
}

namespace xstd {
  template <class Container>
  void sort(Container & c) {
    std::sort(begin(c),end(c));
  }
  template <class Container, class Compare>
  void sort(Container & c, Compare comp) {
    std::sort(begin(c),end(c),comp);
  }
}

int main(int argc, char ** argv)
{
  FUNCSTART("int main(int argc, char ** argv)");

  vector<uint32> x(10);
  xstd::sort(x);                       // works
  xstd::sort(x,std::greater<uint32>());   // does not compile
  cout << x[0] << endl;

  auto c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  fH();
  fL();
  fS();
  fC();
  fM();
  auto c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
	    << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
	    << "Wall clock time passed: "
	    << std::chrono::duration<double, std::milli>(t_end-t_start).count()
	    << " ms\n";
  exit(0);

  //vhash64_t vh(0xf0c06dbcfb02800);
  vhash64_t vh(18);
  uint32 bph=4;

  auto rvh(nsvhash::reverseComplement(vh,bph));
  cout << vh
       << "\t" << HashStatistics<vhash64_t>::hash2string(vh,bph)
       << "\t" << rvh
       << "\t" << HashStatistics<vhash64_t>::hash2string(rvh,bph)
       << endl;

  vh<<=2;
  cout << vh
       << "\t" << HashStatistics<vhash64_t>::hash2string(vh,bph)
       << endl;

//  string acgt="ACGT";
//  string kmf="    ";
//  string kmr;
//  uint32 counter=0;
//  for(uint8 x0=0;x0<4;++x0){
//    for(uint8 x1=0;x1<4;++x1){
//      for(uint8 x2=0;x2<4;++x2){
//	for(uint8 x3=0;x3<4;++x3){
//	  kmf[0]=acgt[x0];
//	  kmf[1]=acgt[x1];
//	  kmf[2]=acgt[x2];
//	  kmf[3]=acgt[x3];
//	  kmr=revComp(kmf);
////	  cout << kmf << "\t" << hex <<  str2hash(kmf) << "\t"
////	       << kmr << "\t" << str2hash(kmr)
////	       << endl;
//	  cout << "0x" << hex << str2hash(kmr) << ", ";
//	  if(++counter%16 == 0 ){
//	    cout << endl;
//	  }
//	}
//      }
//    }
//  }

  try{
    if(argc!=3) {
      cout << "Need exactly 2 parameters: trimfreq merfile\n";
      exit(10);
    }

    uint32 trimfreq=atoi(argv[1]);
    string merfile(argv[2]);

    if(0){
      auto mhs=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(merfile);
      auto bytes=mhs.sizeofhash;
      if(bytes==8){
	HashStatistics<vhash64_t> hs;
	doCorrect(trimfreq,merfile,hs);
      }else if(bytes==16){
	HashStatistics<vhash128_t> hs;
	doCorrect(trimfreq,merfile,hs);
      }else if(bytes==32){
	HashStatistics<vhash256_t> hs;
	doCorrect(trimfreq,merfile,hs);
      }else if(bytes==64){
	HashStatistics<vhash512_t> hs;
	doCorrect(trimfreq,merfile,hs);
      }else{
	MIRANOTIFY(true,"Kmer size " << mhs.basesperhash << " with " << bytes << " bytes are not expected here.\n");
      }
    }

    if(1){
      auto mhs=HashStatistics<vhash64_t>::loadHashStatisticsFileHeader(merfile);
      auto bytes=mhs.sizeofhash;
      if(bytes==8){
	HashStatistics<vhash64_t> hs;
	huntChimeras(trimfreq,merfile,hs);
      }else if(bytes==16){
	HashStatistics<vhash128_t> hs;
	huntChimeras(trimfreq,merfile,hs);
      }else if(bytes==32){
	HashStatistics<vhash256_t> hs;
	huntChimeras(trimfreq,merfile,hs);
      }else if(bytes==64){
	HashStatistics<vhash512_t> hs;
	huntChimeras(trimfreq,merfile,hs);
      }else{
	MIRANOTIFY(true,"Kmer size " << mhs.basesperhash << " with " << bytes << " bytes are not expected here.\n");
      }
    }

  }
  catch(Notify n){
    n.handleError("main");
  }

  exit(0);

  vector<MIRAParameters> Pv;
  ReadPool rp1;
  ReadPool rp2;
  ReadPoolIO rpio(rp1);
  ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
  rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);

  try{
    rpio.registerFile("fastq",
		      "bla.fastq",
		      "",
		      rgid,
		      false);
    rpio.loadNextSeqs(100);
    cout << rp1.size() << " " << rp2.size() << endl;;
    rpio.loadNextSeqs(100);
    cout << rp1.size() << " " << rp2.size()  << endl;

    rpio.setAttributeReadPool(rp2);
    rpio.loadNextSeqs(1000);
    cout << rp1.size() << " " << rp2.size()  << endl;
    rpio.loadNextSeqs(1000);


    exit(0);
  }
  catch(Notify n){
    n.handleError("main");
  }

  exit(0);

  cout << "Have " << MachineInfo::getCoresTotal() << " cores" << endl;
  cout << "Have " << MachineInfo::getMemTotal() << " mem total" << endl;
  cout << "Have " << MachineInfo::getMemAvail() << " mem avail" << endl;

  srand(1234567);

  string bla("bla");
  boost::system::error_code ec;
  try{
    boost::filesystem::remove_all(bla,ec);
  }
  catch(boost::filesystem::filesystem_error fse){
    cout << fse.what() << endl;
  }

  exit(0);

  try {
    timeval tv;
    suseconds_t sus=0;

//    gzFile fp;
//    kseq_t *seq;
//    fp = gzopen("bla.fastq", "r");
//    seq = kseq_init(fp);
//    gettimeofday(&tv,nullptr);
//    int l;
//    while ((l = kseq_read(seq)) >= 0);
//    sus=diffsuseconds(tv);
//    cout << "timing fastq-lh: " << sus << endl;


    FastQ fq;

    gettimeofday(&tv,nullptr);
    fq.openFile("bla.fastq");
    while(fq.loadNext()>=0){
      //cout << fq.getLineCount() << endl;
      //if(fq.getLineCount()%10000==0) cout << fq.getLineCount() << endl;
    }
    cout << fq.getLineCount() << endl;
    sus=diffsuseconds(tv);
    cout << "timing fastq-mira: " << sus << endl;



  }
  catch(Notify n){
    n.handleError("main");
  }

  return 0;
}
