// to be included by assembly, pathfinder ...

#ifndef _mira_overlapedges_h_
#define _mira_overlapedges_h_

#include <iostream>
#include "stdinc/defines.H"
#include "mira/adsfacts.H"

struct newedges_t{
  int32 rid1;                   // number of read in seqpool
  int32 linked_with;                   // number of read in seqpool
  uint32 best_weight;                   // score ratio^2*score. Make sure it is >0 !!! (for pathfinder)
  uint32 adsfindex;              // index of adsfacts elements this refers to
  int16 direction;

  bool pf_banned:1;              // temp use by pathfinder: banned overlap

  bool ol_stronggood:1;  // frequency: 2*bph-1 pos at 3, thereof bph-1 contiguous
  bool ol_weakgood:1;   // frequency: bph-1 positions contiguous at 3
  bool ol_belowavgfreq:1; // frequency: bph-1 positions contiguous at <=3
  bool ol_norept:1;      // nothing >3 (but can contain 1 (single hashes == errors)
  bool ol_rept:1;      // bph-1 positions >=5

  friend std::ostream & operator<<(std::ostream &ostr, const newedges_t & e){
    ostr << "NE:\t" << e.rid1
	 << '\t' << e.linked_with
	 << '\t' << e.best_weight
	 << '\t' << e.adsfindex
	 << "\tdir " << e.direction
	 << "\tban " << e.pf_banned
	 << "\tsg  " << e.ol_stronggood
	 << "\twg  " << e.ol_weakgood
	 << "\tbaf " << e.ol_belowavgfreq
	 << "\tnrp " << e.ol_norept
	 << "\trep " << e.ol_rept
	 << '\n';
    return ostr;
  }

  // sort rid1 from low to high
  inline static bool sortComparatorByRIDUp(const newedges_t & a, const newedges_t & b){
    return a.rid1 < b.rid1;
  };

  // sort rid1 from low to high
  // on equality, sort on bestweight from high to low
  inline static bool sortComparatorByRIDUpByWeightDown(const newedges_t & a, const newedges_t & b) {
    if(a.rid1 == b.rid1){
      return a.best_weight > b.best_weight;
    }
    return a.rid1 < b.rid1;
  }
};


struct skimedges_t{
  int32  rid1;                   // number of read in seqpool
  int32  linked_with;                   // number of read in seqpool
  int32  eoffset;

  uint32 skimweight;                   // score ratio^2*score

  uint8  scoreratio;

  bool ol_stronggood:1;  // frequency: 2*bph-1 pos at 3, thereof bph-1 contiguous
  bool ol_weakgood:1;   // frequency: bph-1 positions contiguous at 3
  bool ol_belowavgfreq:1; // frequency: bph-1 positions contiguous at <=3
  bool ol_norept:1;      // nothing >3 (but can contain 1 (single hashes == errors)
  bool ol_rept:1;      // bph-1 positions >= 5 freq

private:
  bool   rid1dirb:1;  // new for reduceSkimHits2
  bool   rid2dirb:1;

public:
  uint64 skimindex:48;    // index of line in original skim file, plenty of space, should never been maxed out


public:
  inline int8 getRID1dir() const { return rid1dirb ? 1 : -1;}
  inline int8 getRID2dir() const { return rid2dirb ? 1 : -1;}
  inline void setRID1dir(int8 b) {
    if(b>0) {
      rid1dirb=true;
    }else if(likely(b<0)) {
      rid1dirb=false;
    }else{
      std::cout << "internal error: setRID1dir() b==0???" << std::endl;
      exit(100);
    }
  }
  inline void setRID2dir(int8 b) {
    if(b>0) {
      rid2dirb=true;
    }else if(likely(b<0)) {
      rid2dirb=false;
    }else{
      std::cout << "internal error: setRID2dir() b==0???" << std::endl;
      exit(100);
    }
  }
  void swapRID12dirs() { bool tmp=rid1dirb; rid1dirb=rid2dirb; rid2dirb=tmp;}


  friend std::ostream & operator<<(std::ostream &ostr, const skimedges_t & e){
    ostr << "SE:\t" << e.rid1
	 << '\t' << e.linked_with
	 << '\t' << e.eoffset
	 << '\t' << e.skimweight
	 << '\t' << e.skimindex
	 << "\td1 " << static_cast<int16>(e.getRID1dir())
	 << "\td2 " << static_cast<int16>(e.getRID2dir())
	 << "\tsr " << static_cast<uint16>(e.scoreratio)
	 << "\tsg  " << e.ol_stronggood
	 << "\twg  " << e.ol_weakgood
	 << "\tbaf " << e.ol_belowavgfreq
	 << "\tnrp " << e.ol_norept
	 << "\trep " << e.ol_rept
	 << '\n';
    return ostr;
  }

  static bool stdSortCmp(const skimedges_t & a, const skimedges_t & b);
};


// not sure whether worth the effort
// Eco 1655 2*3m Solexa 100bp: 37 GiB instead of 39 GiB
// On the other hand: that was on a table size of ~10 GiB!
//
// If activated again, need to rework some functions which
//  currently use [] to use iterators (when possible)
#if 0
typedef std::deque<AlignedDualSeqFacts> adsfcontainer_t;
typedef std::deque<newedges_t> necontainer_t;
#else
#define MIRA_ADSSTORE_CAPACITY
typedef std::vector<AlignedDualSeqFacts> adsfcontainer_t;
typedef std::vector<newedges_t> necontainer_t;
#endif

#endif
