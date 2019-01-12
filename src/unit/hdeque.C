#include <iostream>
#include <string>

#include <sys/times.h>
#include <limits.h>
#include <unistd.h>

#include "stdinc/types.H"

#include "util/machineinfo.H"


#include "mira/hdeque.H"


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
template <class TT>
void testHDequeItTiming()
{
  TT hd;
  size_t numelem=40*1000*1000;

  timeval tv;
  suseconds_t sus=0;

  int dummy = 0;

  gettimeofday(&tv,nullptr);
  hd.resize(numelem);
  sus=diffsuseconds(tv);
  cout << "timing resize: " << sus << endl;

  auto aI=hd.begin();
  auto bI=aI;
  gettimeofday(&tv,nullptr);
  for(uint32 n=0; n<numelem-1; ++n){
    aI=bI;
    advance(aI,n);
    if(aI==hd.end()){
      cout << "End for " << n << " ????" << endl;
      exit(1);
    }
  }
  sus=diffsuseconds(tv);
  cout << "timing advance to every pos: " << sus << endl;
  dummy+=*reinterpret_cast<int *>(&(*aI));

  aI=hd.begin();
  gettimeofday(&tv,nullptr);
  for(; aI!=hd.end(); ++aI){
    ++aI;
  }
  sus=diffsuseconds(tv);
  cout << "timing iterate: " << sus << endl;
  --aI;
  dummy+=*reinterpret_cast<int *>(&(*aI));

  cout << dummy << endl;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDequeResize()
{
  //HDeque<Contig::consensus_counts_t> hd;
  HDeque<int> hd;
  hd.setBinSize(4);

  Contig::consensus_counts_t ccx;

  hd.resize(14);
  hd.debugDump(false);
  hd.resize(10);
  hd.debugDump(false);
  hd.resize(27);
  hd.debugDump(false);

  hd.clear();

  hd.setBinSize(8192);
  hd.resize(12);
  hd.resize(15);
  hd.debugDumpDeep();


  for(;;);
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDequeErase()
{
  HDeque<int32> hd;
  hd.setBinSize(4);

  for(int32 i=0; i<18; ++i){
    hd.push_back(i);
  }
  hd.debugDumpDeep();

  auto fromI=hd.begin()+5;
  auto toI=fromI+7;
  hd.erase(fromI,toI);

  hd.debugDumpDeep();

  fromI=hd.begin()+4;
  toI=fromI+1;
  hd.erase(fromI,toI);

  hd.debugDumpDeep();

  hd.clear();
  for(int32 i=0; i<18; ++i){
    hd.push_back(i);
  }
  fromI=hd.begin()+5;
  for(int32 i=0; i<6; ++i){
    fromI=hd.erase(fromI);
  }
  hd.debugDumpDeep();

  hd.clear();
  for(int32 i=0; i<18; ++i){
    hd.push_back(i);
  }

  auto hd2=hd;
  fromI=hd2.end();
  for(int32 i=0; i<18; ++i){
    --fromI;
    cout << *fromI << endl;
  }

  for(;;);
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDequeInsert()
{
  HDeque<int32> hd;
  hd.setBinSize(4);

  hd.insert(hd.end(),7,7);
  hd.debugDumpDeep();
  hd.insert(hd.begin()+1,7,6);
  hd.debugDumpDeep();
  hd.insert(hd.begin(),7,5);
  hd.debugDumpDeep();
  hd.insert(hd.begin(),7,4);
  hd.debugDumpDeep();

  cout << "dideldei" << endl;
  hd.resize(10);
  cout << "dideldu" << endl;
  hd.debugDumpDeep();

  auto hdI=hd.begin();
  for(int32 i=0; i>-10; --i){
    hdI=hd.insert(hdI,i);
  }
  cout << "eh 1" << endl;
  hd.debugDumpDeep();
  cout << "eh 2" << endl;

  hdI=hd.end();
  for(int32 i=0; i<10; ++i){
    hdI=hd.insert(hdI,i);
    cout << "inserted " << i << endl;
    hd.debugDumpDeep();
    cout << "got " << hdI << endl;
  }
  cout << "eh 3" << endl;
  hd.debugDumpDeep();
  cout << "eh 4" << endl;

  for(;;);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDeque5()
{
  Contig::consensus_counts_t ccx;
  HDeque<Contig::consensus_counts_t> hd;

  hd.resize(30000,ccx);

  hd.debugDumpDeep();

  for(;;);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDeque4()
{
  size_t numelem=4*1000*1000;
  uint32 reps=100;
  size_t advfact=2;

  Contig::consensus_counts_t ccx;

  timeval tv;
  suseconds_t sus=0;

  deque<Contig::consensus_counts_t> dd;
  dd.resize(numelem,ccx);
  gettimeofday(&tv,nullptr);

  auto dI=dd.begin()+dd.size()/advfact+4096;
  for(uint32 i=0;i<reps; ++i){
    dI=dd.insert(dI,ccx);
  }
  sus=diffsuseconds(tv);
  cout << "timing insert deque: " << sus << endl;

  gettimeofday(&tv,nullptr);

  dI=dd.begin()+dd.size()/advfact+4096;
  for(uint32 i=0;i<reps; ++i){
    dI=dd.erase(dI);
  }
  sus=diffsuseconds(tv);
  cout << "timing erase deque: " << sus << endl;

  dd.clear();



  HDeque<Contig::consensus_counts_t> hd;

  hd.resize(numelem,ccx);

  gettimeofday(&tv,nullptr);

  auto iI=hd.begin()+hd.size()/advfact+4096;
  for(uint32 i=0;i<reps; ++i){
    iI=hd.insert(iI,ccx);
  }
  sus=diffsuseconds(tv);
  cout << "timing insert hdeque: " << sus << endl;


  gettimeofday(&tv,nullptr);

  iI=hd.begin()+hd.size()/advfact+4096;
  for(uint32 i=0;i<reps; ++i){
    iI=hd.erase(iI);
  }
  sus=diffsuseconds(tv);
  cout << "timing erase hdeque: " << sus << endl;


  for(;;);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDeque3()
{
  HDeque<int32> hd;
  hd.setBinSize(4);

  for(int32 i=-1; i>=-20; --i){
    hd.push_front(i);
  }
  hd.debugDumpDeep();
  for(int32 i=0; i<6; ++i){
    hd.pop_back();
  }
  hd.debugDumpDeep();
  for(int32 i=0; i<16; ++i){
    hd.push_back(i);
  }
  hd.debugDumpDeep();
  hd.resize(40);
  hd.debugDumpDeep();

  for(;;);
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDeque2()
{
  size_t numelem=20*1000*1000;

  cout << "Testing " << numelem << " elements" << endl;

  HDeque<Contig::consensus_counts_t> hd;

  Contig::consensus_counts_t ccx;

  timeval tv;
  suseconds_t sus=0;

  gettimeofday(&tv,nullptr);

  for(size_t i=0; i<numelem; ++i){
    //myrand=rand();
    hd.push_back(ccx);
    if(i%1000000==0) cout << i << endl;
  }
  sus=diffsuseconds(tv);
  cout << "timing push_back hdeque consensus_counts_t: " << sus << endl;
  for(;;);
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void testHDeque()
{
  size_t numelem=20*1000*1000;
  //size_t numelem=100000;

  cout << "Testing " << numelem << " elements" << endl;

  HDeque<int32> hd;
  hd.setBinSize(8192);
  hd.clear();

  {
    auto hd2(hd);
  }

  auto P=hd.begin();
  cout << "b: " << P << endl;
  P=hd.end();
  cout << "e: " << P << endl;

  auto Q=hd.cbegin();
  cout << "dist P-Q " << Q-P << endl;

  for(auto I=hd.begin(); I!=hd.end(); ++I){
    if(*I<0) cout << *I << endl;
  }
  cout << "bla" << endl;
  for(auto x : hd){
    if(x<0) cout << x << endl;
  }

  srand(1234567);
  int32 myrand=1;

  timeval tv;
  suseconds_t sus=0;

  gettimeofday(&tv,nullptr);

  {
    deque<int32> d;
    for(size_t i=0; i<numelem; ++i){
      d.push_front(myrand);
    }
    sus=diffsuseconds(tv);
    cout << "timing push_front deque: " << sus << endl;
    gettimeofday(&tv,nullptr);

    for(auto x : d){
      if(x<0) cout << x << endl;
    }
    sus=diffsuseconds(tv);
    cout << "timing traverse deque: " << sus << endl;
    gettimeofday(&tv,nullptr);
  }

  {
    deque<int32> d;
    for(size_t i=0; i<numelem; ++i){
      d.push_back(myrand);
    }
    sus=diffsuseconds(tv);
    cout << "timing push_back deque: " << sus << endl;
    gettimeofday(&tv,nullptr);

    for(auto x : d){
      if(x<0) cout << x << endl;
    }
    sus=diffsuseconds(tv);
    cout << "timing traverse deque: " << sus << endl;
    gettimeofday(&tv,nullptr);
  }

  for(size_t i=0; i<numelem; ++i){
    //myrand=rand();
    hd.push_back(myrand);
  }
  sus=diffsuseconds(tv);
  cout << "timing push_back: " << sus << endl;
  gettimeofday(&tv,nullptr);

  for(auto x : hd){
    if(x<0) cout << x << endl;
  }
  sus=diffsuseconds(tv);
  cout << "timing traverse: " << sus << endl;

  gettimeofday(&tv,nullptr);
  {
    size_t sum=0;
    auto tI=hd.begin();
    while(tI!=hd.end()){
      sum+=*tI;
      tI+=1;
    }
    cout << "Sum " << sum << endl;
  }
  sus=diffsuseconds(tv);
  cout << "timing traverse by foot +: " << sus << endl;

  gettimeofday(&tv,nullptr);
  {
    size_t sum=0;
    auto tI=hd.end();
    while(tI!=hd.begin()){
      tI-=1;
      sum+=*tI;
    }
    cout << "Sum " << sum << endl;
  }
  sus=diffsuseconds(tv);
  cout << "timing traverse by foot -: " << sus << endl;


  cout << "size computed: " << hd.end()-hd.begin() << endl;


  gettimeofday(&tv,nullptr);
  while(hd.size()){
    hd.pop_back();
  }
  sus=diffsuseconds(tv);
  cout << "timing pop_back: " << sus << endl;

  for(size_t i=0; i<numelem; ++i){
    //myrand=rand();
    hd.push_front(myrand);
  }
  sus=diffsuseconds(tv);
  cout << "timing push_front: " << sus << endl;
  gettimeofday(&tv,nullptr);

  for(auto x : hd){
    if(x<0) cout << x << endl;
  }
  sus=diffsuseconds(tv);
  cout << "timing traverse: " << sus << endl;
  gettimeofday(&tv,nullptr);

  while(hd.size()){
    hd.pop_front();
  }
  sus=diffsuseconds(tv);
  cout << "timing pop_front: " << sus << endl;
  gettimeofday(&tv,nullptr);

  hd.resize(numelem);
  sus=diffsuseconds(tv);
  cout << "timing resize to " << numelem << ": " << sus << endl;
  gettimeofday(&tv,nullptr);

  hd.resize(2000);
  sus=diffsuseconds(tv);
  cout << "timing resize to 2000: " << sus << endl;
  gettimeofday(&tv,nullptr);


//  hd.debugDump(false);

  for(;;);
}


class MyTest {
private:
  PlacedContigReads mt_pcr;
public:
  MyTest(ReadPool & rp) : mt_pcr(rp) {};

  PlacedContigReads & getPCR() {
    return mt_pcr;
  };
};

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
int main(int argc, char ** argv)
{
  FUNCSTART("int main(int argc, char ** argv)");

  cout << "Have " << MachineInfo::getCoresTotal() << " cores" << endl;
  cout << "Have " << MachineInfo::getMemTotal() << " mem total" << endl;
  cout << "Have " << MachineInfo::getMemAvail() << " mem avail" << endl;

  try {
    cout << "deque:\n------------------\n";
    testHDequeItTiming<deque<Contig::consensus_counts_t> >();
    cout << "\n\n\nHDeque:\n------------------\n";
    testHDequeItTiming<HDeque<Contig::consensus_counts_t> >();
    exit(0);
    //testHDequeResize();
    //testHDequeInsert();
    //testHDequeErase();
    //testHDeque();
    //testHDeque2();
    //testHDeque3();
    //testHDeque4();
    //testHDeque5();

    exit(0);
  }
  catch(Notify n){
    n.handleError("main");
  }

  return 0;
}
