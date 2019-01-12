#include "memorc/memorc.H"

#include <cstdlib>      // new(), abort(), free() on Cygwin
#include <cstdio>       // printf() on Cygwin
#include <cstring>
#include <iostream>


using std::cout;


bool MemORC::MOC_readytouse=false;
bool MemORC::MOC_newallocsgetthisnocheckflag=true;
bool MemORC::MOC_exitset=false;

const int32 MemORC::MOC_fencesize=10;

uint8 MemORC::MOC_clowfencemagic  [4] = {0xab, 0xad, 0xca, 0xfe};
uint8 MemORC::MOC_chighfencemagic [4] = {0xde, 0xad, 0xbe, 0xef};

int32 MemORC::MOC_lowfencemagic= *((int32 *) MemORC::MOC_clowfencemagic);
int32 MemORC::MOC_highfencemagic= *((int32 *) MemORC::MOC_chighfencemagic);

uint8 MemORC::MOC_maskmagic=0xaa;

uint32 MemORC::MOC_innewcount=0;

bool MemORC::MOC_mostfatalerroroccurred=false;
bool MemORC::MOC_allhot=false;

bool MemORC::MOC_fillondelete=false;

uint64 MemORC::MOC_alloccounter=0;
uint64 MemORC::MOC_checksperformed=0;

size_t MemORC::MOC_reqmemalloc=0;
size_t MemORC::MOC_totalmemalloc=0;



//#define LEBUG(bla) {bla;}
#define LEBUG(bla)

#if 1
/////////////////////////////////////////////////////////////////////

void * operator new[](size_t n)
{
  LEBUG(printf("new [] called\n"));

  auto * newmem=MemORC::newMemBlock(n);
  if(newmem==nullptr) throw std::bad_alloc();

  LEBUG(printf("giving back: %p\n", newmem));
  return newmem;
}

void operator delete[](void *mb)
{
  LEBUG(printf("delete [] called with adress: %p\n", mb));
  MemORC::deleteMemBlock(mb);
}

void * operator new(size_t n)
{
  LEBUG(printf("new called\n"));

  auto * newmem=MemORC::newMemBlock(n);
  if(newmem==nullptr) throw std::bad_alloc();

  LEBUG(printf("giving back: %p\n", newmem));
  return newmem;
}

void operator delete(void *mb)
{
  LEBUG(printf("delete called with adress: %p\n", mb));
  MemORC::deleteMemBlock(mb);
}

/////////////////////////////////////////////////////////////////////
#else
/////////////////////////////////////////////////////////////////////

void * operator new[](size_t n)
{
  return malloc(n);
}

void operator delete[](void *mb)
{
  free(mb);
}

void * operator new(size_t n)
{
  return malloc(n);
}

void operator delete(void *mb)
{
  free(mb);
}
/////////////////////////////////////////////////////////////////////
#endif


MemORC::MemORC()
{
  LEBUG(printf("MemORC is ready to use, fencesize %d\n",MOC_fencesize));
  MOC_readytouse=true;
}

void MemORC::atExit()
{
  LEBUG(printf("exiting, stopping MemORC\n"));
  statistics();
}

MemORC::~MemORC()
{
  MOC_readytouse=false;
}

void MemORC::setChecking(bool b)
{
  MOC_newallocsgetthisnocheckflag=!b;
  if(!MOC_exitset){
    atexit(&MemORC::atExit);
    MOC_exitset=true;
  }
  if(b){
    printf("memory allocation now checks fences, this will cause some slow down\n");
  }else{
    printf("switching off check of fences, checks now apply only to existing memory\n");
  }
  LEBUG(printf("MOC_innewcount %d\n", MOC_innewcount));
}

void MemORC::setAllHot(bool b)
{
  MOC_allhot=b;
  if(b){
    printf("all hot memory checks, expect major slowdown on large programs\n");
  }else{
    printf("all hot now off, things should go faster\n");
  }
}

void MemORC::statistics()
{
  printf("MemORC statistics:\n");
  printf("Current alloc counter  : %llu\n", MOC_alloccounter);
  printf("Total checks performed : %llu\n", MOC_checksperformed);
  printf("Num. memory blocks     : %llu\n", MOC_memblocks.size());
  printf("Num. active hot blocks : %llu\n", MOC_hotblocks.size());
  printf("Num. pending hot blocks: %llu\n", MOC_hotaidsrequested.size());
  printf("Total mem reqested     : %llu\n", MOC_reqmemalloc);
  printf("Total mem allocated    : %llu (including pads and fences)\n", MOC_totalmemalloc);
}

void MemORC::myexit(int32 n)
{
  LEBUG(printf("shouldn't check anything anymore\n"));
  MemORC::MOC_readytouse=false;
  statistics();
  abort();
  exit(n);
}


void * MemORC::internalMAlloc(size_t n)
{
  if(MOC_mostfatalerroroccurred){
    abort();
  }
  void * newmem=malloc(n);
  // Most significant bit set? Come on ... not on x86_64! Must be a magic number
  //  for invalid mem
  if(reinterpret_cast<uint64>(newmem)&0x8000000000000000ull){
    MOC_mostfatalerroroccurred=true;
    printf("Houston, we have a malloc pointer problem: %p\n",newmem);
    statistics();
    MemORC::checkAllMemBlocks();
    printf("No error found yet, exiting\n");
    myexit(1000);
  }
  return newmem;
}


void * MemORC::newMemBlock(size_t n)
{
  LEBUG(printf("MOC_innewcount %d\n", MOC_innewcount));
  void * newmem;
  if(!MOC_readytouse || MOC_innewcount){
    if(!MOC_readytouse){
      LEBUG(printf("not ready malloc\n"));
    }else {
      LEBUG(printf("innew malloc\n"));
    }
    newmem=internalMAlloc(n);
  }else{
    ++MOC_innewcount;
    LEBUG(printf("my malloc, want %llu\n",n));

    if(MOC_allhot){
      checkAllMemBlocks();
    } else if(!MOC_hotblocks.empty()){
      checkAllHotBlocks();
    }

    size_t totaln= n +  2*(sizeof(MOC_lowfencemagic) * MOC_fencesize);

    if(totaln%4){
      totaln+=4-(totaln%4);
    }

    //if(totaln<256)totaln=256;
    LEBUG(printf("my malloc, getting %llu\n",totaln));
    newmem=internalMAlloc(totaln);
    LEBUG(printf("my malloc, addr %p to %p\n",newmem,newmem+totaln));

    if(newmem==nullptr){
      printf("Used malloc()\nTried to allocate %llu byte.\nSorry, memory allocation failed. Bailing out.\n", totaln);
      myexit(1);
    }

    MOC_reqmemalloc+=n;
    MOC_totalmemalloc+=totaln;

    int32 * magicptr=static_cast<int32*>(newmem);
    for(int32 i=0; i<MOC_fencesize; ++i) *magicptr++=MOC_lowfencemagic;

    newmem=magicptr;
    magicptr+=(n/sizeof(int32));
    *magicptr=0xbbbbbbbb;

    uint8 * maskptr=(uint8 *) magicptr;

    switch(n%4){
    case 3:{
      *maskptr++=0;
      *maskptr++=0;
      *maskptr++=0;
      *maskptr++=MOC_maskmagic;
      break;
    }
    case 2:{
      *maskptr++=0;
      *maskptr++=0;
      *maskptr++=MOC_maskmagic;
      *maskptr++=MOC_maskmagic;
      break;
    }
    case 1:{
      *maskptr++=0;
      *maskptr++=MOC_maskmagic;
      *maskptr++=MOC_maskmagic;
      *maskptr++=MOC_maskmagic;
      break;
    }
    case 0:{
      break;
    }
    }

    magicptr=(int32 *) maskptr;

    for(int32 i=0; i<MOC_fencesize; i++) *magicptr++=MOC_highfencemagic;

    LEBUG(printf("done init\n"));
    bool hashot=false;
    ++MOC_alloccounter;
    for(auto hvI=MOC_hotaidsrequested.begin(); hvI!=MOC_hotaidsrequested.end(); ++hvI){
      if(*hvI==MOC_alloccounter){
	hashot=true;
	MOC_hotaidsrequested.erase(hvI);
	break;
      }
    }
    MOC_memblocks.insert(mbmapelem_t(newmem,memblockinfo_t(n,hashot,MOC_newallocsgetthisnocheckflag,MOC_alloccounter)));
    if(hashot){
      MOC_hotblocks.insert(mbmapelem_t(newmem,memblockinfo_t(n,hashot,MOC_newallocsgetthisnocheckflag,MOC_alloccounter)));
    }

    --MOC_innewcount;

  }
  return newmem;
}

void MemORC::deleteMemBlock(void * usraddr)
{
  LEBUG(printf("dmb\n"));
  LEBUG(printf("MOC_innewcount %d\n", MOC_innewcount));
  if(usraddr==nullptr) return;
  if(!MOC_readytouse){
    LEBUG(printf("not ready free %p\n", usraddr));
    free(usraddr);
    return;
  }
  if(MOC_memblocks.find(usraddr)==MOC_memblocks.end()){
    if(MOC_innewcount==0 && !MOC_memblocks.empty()){
      ++MOC_innewcount;
      printf("Houston, trying to free something we haven't seen: %p\n",usraddr);
      auto mbI=MOC_memblocks.lower_bound(usraddr);
      if(mbI==MOC_memblocks.end()){
	printf("No nearest memblock above\n");
      }else{
	printf("Nearest memblock above (%d):\t",reinterpret_cast<uint64>(mbI->first)-reinterpret_cast<uint64>(usraddr));
	printMemBlockInfo(mbI);
      }
      if(mbI==MOC_memblocks.begin()){
	printf("No nearest memblock below\n");
      }else{
	--mbI;
	printf("Nearest memblock below (-%d):\t",reinterpret_cast<uint64>(usraddr)-reinterpret_cast<uint64>(mbI->first));
	printMemBlockInfo(mbI);
      }
      statistics();
      MemORC::checkAllMemBlocks();
      printf("No error found yet, exiting\n");
      myexit(1000);
    }
    ++MOC_innewcount;
    if(reinterpret_cast<uint64>(usraddr)&0x8000000000000000ull){
      printf("Houston, we have a pointer problem for free: %p\n",usraddr);
      statistics();
      MemORC::checkAllMemBlocks();
      printf("No error found yet, exiting\n");
      myexit(1000);
    }
    LEBUG(printf("oooops orig free %p\n", usraddr));
    free(usraddr);
    --MOC_innewcount;
    return;
  }

  ++MOC_innewcount;
  LEBUG(printf("my free\n"));

  auto mbI=MOC_memblocks.end();  // .end() because I'm lazy, want to use "auto"
  auto haserror=checkMemBlock(usraddr,mbI);
  if(haserror){
    printf("Overruns detected, exiting!\n");
    myexit(100);
  }

  if(MOC_allhot){
    checkAllMemBlocks();
  } else if(!MOC_hotblocks.empty()){
    checkAllHotBlocks();
  }

  size_t size=mbI->second.size;
  MOC_reqmemalloc-=size;

  uint8 * ptr=static_cast<uint8 *>(mbI->first);
  // align size to be filled
  if(size%4){
    size+=4-(size%4);
  }
  uint8 * ptrfh=ptr+size;
  uint8 * ptrfl=ptr-MOC_fencesize*sizeof(MOC_lowfencemagic);

  // add fence size, magic sizes and size of a size_t to it
  MOC_totalmemalloc-=(size+2*MOC_fencesize*sizeof(MOC_lowfencemagic));


  if(MOC_fillondelete){
    memset(ptrfl,0xaa,MOC_fencesize*sizeof(MOC_lowfencemagic));
    memset(ptr,0xbb,size);
    memset(ptrfh,0xcc,MOC_fencesize*sizeof(MOC_highfencemagic));
  }

  if(mbI->second.hashot){
    auto tmpI=MOC_hotblocks.erase(usraddr);
  }
  MOC_memblocks.erase(usraddr);

  LEBUG(printf("my free %p\n", ptrfl));
  free(ptrfl);

  --MOC_innewcount;

  return;
}


bool MemORC::checkMemBlock(void *mb, mbmap_t::iterator & retI)
{
  if(!MOC_readytouse) return false;

  LEBUG(printf("check raw user addr %p\n",mb));

  ++MOC_innewcount;
  retI=MOC_memblocks.find(mb);
  if(retI==MOC_memblocks.end()){
    printf("checking memory block that was not allocated by MemORC ???: %p\n", (uint8 *) mb);
    printf("Bailing out\n");
    myexit(1);
    return false;
  }

  auto ret=checkMemBlock(retI);
  --MOC_innewcount;
  return ret;
}

bool MemORC::checkMemBlock(mbmap_t::iterator mbI)
{
  if(!MOC_readytouse) return false;
  if(mbI->second.nocheck) return false;

  LEBUG(printf("checkMemBlock mbI %p\n", mbI->first));

  ++MOC_checksperformed;

  ++MOC_innewcount;

  bool haserror=false;

  int32 * fenceptr=static_cast<int32 *>(mbI->first);
  fenceptr-=MOC_fencesize;
  for(int32 i=-MOC_fencesize; i<0; ++i, ++fenceptr){
    if(*fenceptr!=MOC_lowfencemagic){
      uint8 * fencecptr=static_cast<uint8 *>(static_cast<void *>(fenceptr));
      if(*fencecptr!=MOC_clowfencemagic[0]){
	printf("check: lower fence .0 destroyed at position: %d\tExpected: %x\tGot: %x\n", i*sizeof(int32),MOC_clowfencemagic[0], *fencecptr);
      }
      if(*(++fencecptr)!=MOC_clowfencemagic[1]){
	printf("check: lower fence .1 destroyed at position: %d\tExpected: %x\tGot: %x\n", i*sizeof(int32)+1,MOC_clowfencemagic[1], *fencecptr);
      }
      if(*(++fencecptr)!=MOC_clowfencemagic[2]){
	printf("check: lower fence .2 destroyed at position: %d\tExpected: %x\tGot: %x\n", i*sizeof(int32)+2,MOC_clowfencemagic[2], *fencecptr);
      }
      if(*(++fencecptr)!=MOC_clowfencemagic[3]){
	printf("check: lower fence .3 destroyed at position: %d\tExpected: %x\tGot: %x\n", i*sizeof(int32)+3,MOC_clowfencemagic[3], *fencecptr);
      }
      haserror=true;
    }
  }

  uint32 ufpos=1;
  uint8 * maskptr=static_cast<uint8 *>(mbI->first);
  maskptr+=mbI->second.size;
  switch(mbI->second.size%4){
  case 0: break;
  case 1:{
    if(*(maskptr) != MOC_maskmagic){
      printf("check: upper fence (mm) destroyed at position: %d\tExpected: %x\tGot: %x\n",ufpos,MOC_maskmagic,*maskptr);
      haserror=true;
    }
    ++ufpos;
    ++maskptr;
  }
  case 2:{
    if(*(maskptr) != MOC_maskmagic){
      printf("check: upper fence (mm) destroyed at position: %d\tExpected: %x\tGot: %x\n",ufpos,MOC_maskmagic,*maskptr);
      haserror=true;
    }
    ++ufpos;
    ++maskptr;
  }
  case 3:{
    if(*(maskptr) != MOC_maskmagic){
      printf("check: upper fence (mm) destroyed at position: %d\tExpected: %x\tGot: %x\n",ufpos,MOC_maskmagic,*maskptr);
      haserror=true;
    }
    ++ufpos;
    ++maskptr;
  }
  }
  fenceptr=static_cast<int32 *>(static_cast<void*>(maskptr));

  for(int32 i=0; i<MOC_fencesize; ++i, ++fenceptr, ufpos+=4){
    if(*fenceptr!=MOC_highfencemagic){
      uint8 * fencecptr=static_cast<uint8 *>(static_cast<void *>(fenceptr));
      if(*fencecptr!=MOC_chighfencemagic[0]){
	printf("check: upper fence .0 destroyed at position: %d\tExpected: %x\tGot: %x\n", ufpos,MOC_chighfencemagic[0], *fencecptr);
      }
      if(*(++fencecptr)!=MOC_chighfencemagic[1]){
	printf("check: upper fence .1 destroyed at position: %d\tExpected: %x\tGot: %x\n", ufpos+1,MOC_chighfencemagic[1], *fencecptr);
      }
      if(*(++fencecptr)!=MOC_chighfencemagic[2]){
	printf("check: upper fence .2 destroyed at position: %d\tExpected: %x\tGot: %x\n", ufpos+2,MOC_chighfencemagic[2], *fencecptr);
      }
      if(*(++fencecptr)!=MOC_chighfencemagic[3]){
	printf("check: upper fence .3 destroyed at position: %d\tExpected: %x\tGot: %x\n", ufpos+3,MOC_chighfencemagic[3], *fencecptr);
      }
      haserror=true;
    }
  }

  if(haserror){
    cout.flush();
    printf("\nError while checking a memoryblock:\n");
    printf("Addr: %p\n", mbI->first);
    printf("Size: %llu\n", mbI->second.size);
    printf("AllocID: %llu\n", mbI->second.allocid);
  }

  --MOC_innewcount;
  return haserror;
}

void MemORC::checkAllMemBlocks()
{
  if(!MOC_readytouse) return;
  LEBUG(printf("check all\n"));

  ++MOC_innewcount;

  bool haserror=false;
  for(auto mbI=MOC_memblocks.begin(); mbI!=MOC_memblocks.end(); ++mbI){
    haserror|=checkMemBlock(mbI);
  }
  if(haserror){
    printf("Overruns detected, exiting!\n");
    myexit(100);
  }
  --MOC_innewcount;
}

void MemORC::checkAllHotBlocks()
{
  if(!MOC_readytouse) return;
  LEBUG(printf("check hb\n"));

  ++MOC_innewcount;

  bool haserror=false;
  for(auto mbI=MOC_hotblocks.begin(); mbI!=MOC_hotblocks.end(); ++mbI){
    haserror|=checkMemBlock(mbI);
  }
  if(haserror){
    printf("Overruns in hot blocks detected, exiting!\n");
    myexit(100);
  }
  --MOC_innewcount;
}


void MemORC::requestHotAlloicID(uint64 aid)
{
  ++MOC_innewcount;
  MOC_hotaidsrequested.push_back(aid);
  --MOC_innewcount;
};


void MemORC::printMemBlockInfo(mbmap_t::iterator mbI)
{
  if(mbI==MOC_memblocks.end()){
    printf("no memblock, end of list");
  }
  printf("%p\t%llu\t",mbI->first,mbI->second.size);
  if(mbI->second.hashot) {
    printf("hot");
  }else{
    printf("normal");
  }
  printf("\t%llu\n",mbI->second.allocid);
}


