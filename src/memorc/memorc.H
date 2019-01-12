#ifndef _bas_memorc_h_
#define _bas_memorc_h_

#include "stdinc/defines.H"

#include <vector>
#include <map>

class MemORC
{
private:
  static MemORC MOC_semaphore;
  static bool MOC_readytouse;
  static bool MOC_trackingon;
  static bool MOC_exitset;

  static const int32 MOC_fencesize;

  static uint32 MOC_innewcount;

  static bool MOC_mostfatalerroroccurred;
  static bool MOC_newallocsgetthisnocheckflag;
  static bool MOC_allhot;
  static bool MOC_fillondelete; // unused atm

  static uint8 MOC_clowfencemagic[4];
  static uint8 MOC_chighfencemagic[4];

  static int32 MOC_lowfencemagic;
  static int32 MOC_highfencemagic;

  static uint8 MOC_maskmagic;

  static uint64 MOC_alloccounter;
  static uint64 MOC_checksperformed;

  static size_t MOC_reqmemalloc;
  static size_t MOC_totalmemalloc;

  struct memblockinfo_t {
    size_t size;
    uint64 allocid;

    bool   hashot:1;
    bool   nocheck:1;

    memblockinfo_t(size_t s, bool hot, bool nc, uint64 aid) : size(s), allocid(aid), hashot(hot), nocheck(nc) {};
  };
  typedef std::map<void *,memblockinfo_t> mbmap_t;
  typedef std::pair<void *,memblockinfo_t> mbmapelem_t;

  // Note: these 3 are not defined in memorc.C!
  // to prevent the typical C++ init fiasco, they were integrated in preventinitfiasco.C
  static mbmap_t MOC_memblocks;
  static mbmap_t MOC_hotblocks;
  static std::vector<uint64> MOC_hotaidsrequested;

private:
  static void myexit(int32 n);

  static bool staticInit();

  static bool checkMemBlock(void * usraddr, mbmap_t::iterator & retI);
  static bool checkMemBlock(mbmap_t::iterator mbI);

  static void checkAllHotBlocks();

  static void * internalMAlloc(size_t n);

  static void printMemBlockInfo(mbmap_t::iterator mbI);

  // No error: contructor and destructor are private, no one else can instantiate
  MemORC();
  ~MemORC();
  static void atExit();

public:
  static void * newMemBlock(size_t n);
  static void deleteMemBlock(void * usraddr);

  static void setTracking(bool b);
  static void statistics();
  static void checkAllMemBlocks();
  static void requestHotAlloicID(uint64 aid);
  static void setAllHot(bool b);
  static void setChecking(bool b);
};



#endif
