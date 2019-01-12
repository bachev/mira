/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2008 and later by Bastien Chevreux
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


#ifndef _mira_offsetcontainer_h_
#define _mira_offsetcontainer_h_


#include <vector>
#include <deque>
#include <list>


#if 0
#define IDCEBUG(bla)   {cout << bla; cout.flush();}
#define IDBUGSTAT(statement) { statement;}
#else
#define IDCEBUG(bla)
#define IDBUGSTAT(statement)
#endif

template<class KK, class VV>
class OffsetContainer
{
public:
  typedef size_t size_type;

  void dumpStats(){
    cout << "Binsize: " << OC_binsize << endl;
    cout << "OC_size: " << OC_size << endl;
    cout << "OC_outerc.size(): " << OC_outerc.size() << endl;
    for(uint32 mi=0; mi<OC_map.size();++mi){
      cout << mi << " " << OC_map[mi];
    }
  }


private:

  struct pair_t {
    KK   offset;   // key
    VV   value;

    pair_t(const KK & o, const VV & v) : offset(o), value(v) {};

    // Need this for containers using .resize() without default element
    pair_t() {};

//    inline bool operator<(const pair_t & other){
//      return offset < other.offset;
//    }
    inline bool operator<(const KK & otheroffset){
      return offset < otheroffset;
    }
  };


  // Using a deque here uses less memory, but is slower
  //typedef deque<pair_t> innerc_t;
  typedef vector<pair_t> innerc_t;

  typedef list<innerc_t> outerc_t;

  struct mapinfo_t {
    KK         first_offsetbase;
    KK         last_offsetbase;
    innerc_t * icptr;

    mapinfo_t(const KK & fo, const KK & lo, innerc_t * i) : first_offsetbase(fo), last_offsetbase(lo), icptr(i) {};
    mapinfo_t(const KK & o, innerc_t * i) : first_offsetbase(o), last_offsetbase(o), icptr(i) {};
    mapinfo_t(const KK & o) : first_offsetbase(o), last_offsetbase(o), icptr(nullptr) {};
    friend ostream & operator<<(ostream &ostr, const mapinfo_t & i){
      ostr << "fob: " << i.first_offsetbase << "\tlob: " << i.last_offsetbase ;
      ostr.flush();
      ostr << "\ts: " << (i.icptr)->size()
	   << endl;
      return ostr;
    }

    //inline static bool lobcomp(const mapinfo_t & a, const mapinfo_t & b) {
    //  return a.last_offsetbase < b.last_offsetbase;
    //}
    inline static bool lobcompKK(const mapinfo_t & a, const KK & b) {
      return a.last_offsetbase < b;
    }
    //inline bool operator<(const mapinfo_t & other){
    //  return offsetbase<other.offsetbase;
    //}
  };

  size_t     OC_binsize;  // Note: may never(!!!) be < 2  ----> some functions assume this, and a segfault would ensue

  size_t     OC_size;
  outerc_t   OC_outerc;

  vector<mapinfo_t> OC_map;



  //inline static bool lmicomp(const mapinfo_t & a, const mapinfo_t & b) {
  //  return a.to <= b.to;
  //}
  inline typename vector<mapinfo_t>::iterator locateMapInfo(KK offset){
    return lower_bound(OC_map.begin(),
		       OC_map.end(),
		       offset,
		       mapinfo_t::lobcompKK);
  }




public:
  void dumpDebug(){
    for(uint32 mi=0; mi<OC_map.size();++mi){
      cout << "Map " << mi << " " << OC_map[mi];
      for(auto & ie: *(OC_map[mi].icptr)){
	cout << "Offset: " << ie.offset << "\tVal: " << ie.value << endl;
      }
    }
  }


  /**********************************************************************************************************
   *
   *
   *  Iterator
   *
   *
   **********************************************************************************************************/

  class iterator
    : public boost::iterator_facade<
    iterator                          // Derived
    , pair_t              // Value
    , std::random_access_iterator_tag  //boost::forward_traversal_tag    // CategoryOrTraversal
    > {

  private:
    //typename vector<T>::iterator cI;
    OffsetContainer * OCI_origin;
    size_t OCI_elemcounter;
    size_t OCI_mapindex;
    typename innerc_t::iterator OCI_icI;
    typename innerc_t::iterator OCI_icIb;
    typename innerc_t::iterator OCI_icIe;

  public:
    inline iterator(OffsetContainer * origin, size_t elemcounter, size_t mapindex, typename innerc_t::iterator icI, typename innerc_t::iterator icIb, typename innerc_t::iterator icIe) : OCI_origin(origin), OCI_elemcounter(elemcounter), OCI_mapindex(mapindex), OCI_icI(icI), OCI_icIb(icIb), OCI_icIe(icIe) { };

    inline KK getOffset() const {
      return OCI_origin->OC_map[OCI_mapindex].first_offsetbase+OCI_icI->offset;
    }
    inline VV & getValue() const {
      return OCI_icI->value;
    }

    void dumpDebug() {
      cout << "OCI_origin " << OCI_origin
	   << "\nOCI_elemcounter: " << OCI_elemcounter
	   << "\nOCI_mapindex: " << OCI_mapindex
	   << endl;
    }

//  private:
//    inline iterator() {};


  private:
    // Everything Boost's iterator facade needs
    friend class boost::iterator_core_access;

//    iterator & dereference() const
//    {
//      return this;
//    }

    bool equal(iterator const & other) const
    {
      return OCI_elemcounter == other.OCI_elemcounter && OCI_origin == other.OCI_origin;
    }

    void increment()
    {
      ++OCI_elemcounter;
      ++OCI_icI;
      if(OCI_icI == OCI_icIe) {
	++OCI_mapindex;
	if(OCI_mapindex < OCI_origin->OC_map.size()){
	  OCI_icI=OCI_origin->OC_map[OCI_mapindex].icptr->begin();
	  OCI_icIb=OCI_icI;
	  OCI_icIe=OCI_origin->OC_map[OCI_mapindex].icptr->end();
	}
      }
    }

    void decrement()
    {
      if(OCI_icI == OCI_icIb) {
	--OCI_mapindex;
	OCI_icI=OCI_origin->OC_map[OCI_mapindex].icptr->end();
	OCI_icIe=OCI_icI;
	OCI_icIb=OCI_origin->OC_map[OCI_mapindex].icptr->begin();
      }
      --OCI_icI;
      --OCI_elemcounter;
    }
//
//    void advance(int64 n)
//    {
//      cI+=n;
//    }
//
//    int64 distance_to(iterator const & other) const
//    {
//      return other.cI - cI;
//    }
  };


  /**********************************************************************************************************
   *
   *
   *  OffsetContainer
   *
   *
   **********************************************************************************************************/

public:
  iterator begin() {
    return iterator(this, 0, 0, OC_map[0].icptr->begin(), OC_map[0].icptr->begin(), OC_map[0].icptr->end());
  }
  iterator end() {
    return iterator(this, OC_size, OC_map.size()-1, OC_map[OC_map.size()-1].icptr->end(), OC_map[OC_map.size()-1].icptr->begin(), OC_map[OC_map.size()-1].icptr->end());
  }

  OffsetContainer() : OC_binsize(1024) { clear(); };

  ~OffsetContainer() {
    //dumpStats();
  }

  void setBinSize(size_t s) {
    if(s<2){
      __throw_out_of_range("internal error OffsetContainer::setBinSize: size may not be < 2");
    }
    OC_binsize=s;
  }

  inline bool empty() const { return OC_size==0; }
  inline size_t size() const { return OC_size; }

  inline void swap(OffsetContainer & other) {
    swap(OC_size,other.OC_size);
    swap(OC_binsize,other.OC_binsize);
    swap(OC_outerc,other.OC_outerc);
    swap(OC_map,other.OC_map);
  }

  inline void clear() {
    IDBUGSTAT(cout << "clear enter\n"; dumpStats(););

    OC_size=0;
    OC_outerc.clear();
    OC_outerc.push_back(innerc_t());
    OC_map.clear();
    OC_map.push_back(mapinfo_t(0,0,&const_cast<innerc_t &>(OC_outerc.back())));

    IDBUGSTAT(cout << "clear at end\n"; dumpStats(); cout << "clear end" << endl;);
  }

  void insert(KK offset, const VV & val){
    auto mI=locateMapInfo(offset);
    if(mI==OC_map.end()){
      IDBUGSTAT(cout << "Map end!\n");
      --mI;
    }

    IDBUGSTAT(cout << "I found " << *mI);

    typename innerc_t::iterator iI;

    if(mI->icptr->size() >= OC_binsize){
      IDBUGSTAT(cout << "full, "; cout.flush());
      innerc_t tmp;
      OC_outerc.push_back(tmp);
      if(mI+1==OC_map.end()
	 && offset >= mI->last_offsetbase){
	IDBUGSTAT(cout << "push_back new"; cout.flush());
	OC_map.push_back(mapinfo_t(offset,offset,&const_cast<innerc_t &>(OC_outerc.back())));
	mI=OC_map.end();
	--mI;
	iI=mI->icptr->end();
      }else if(mI==OC_map.begin()
	 && offset <= mI->first_offsetbase){
	IDBUGSTAT(cout << "push_front new"; cout.flush());
	OC_map.insert(OC_map.begin(),mapinfo_t(offset,offset,&const_cast<innerc_t &>(OC_outerc.back())));
	mI=OC_map.begin();
	iI=mI->icptr->end();
      }else{
	IDBUGSTAT(cout << "need a split"; cout.flush());

	auto blaI=mI->icptr->begin();
	advance(blaI,mI->icptr->size()/2);

	KK suboffset=blaI->offset;
	mapinfo_t newinfo(mI->first_offsetbase+suboffset,mI->last_offsetbase,&const_cast<innerc_t &>(OC_outerc.back()));

	mI->last_offsetbase=mI->first_offsetbase+(blaI-1)->offset;

	// No! Do by hand, need to adjust offset anyway
	// OC_outerc.back().insert(OC_outerc.back().end(),
	//			  blaI,
	//			  mI->icptr->end());

	for(; blaI != mI->icptr->end(); ++blaI){
	  blaI->offset-=suboffset;
	  OC_outerc.back().push_back(*blaI);
	}
	mI->icptr->resize(mI->icptr->size()/2);

	++mI;
	OC_map.insert(mI,newinfo);

	IDBUGSTAT(cout << "Done with split, I now have:\n"; dumpDebug());

	// simple and efficient: one level of recursion to really insert the element

	insert(offset,val);
	return;
      }
    }else{
      iI=lower_bound(mI->icptr->begin(),
		     mI->icptr->end(),
		     offset-mI->first_offsetbase);
    }

    pair_t tmpp(offset-mI->first_offsetbase,val);

    if(mI->icptr->size()==0){
      mI->first_offsetbase=offset;
      mI->last_offsetbase=offset;
      tmpp.offset=0;
    }else{
      if(offset<mI->first_offsetbase){
	KK tmp=mI->first_offsetbase-offset;
	for(auto & ao : *(mI->icptr)){
	  ao.offset+=tmp;
	}
	mI->first_offsetbase=offset;
	tmpp.offset=0;
      }
      mI->last_offsetbase=max(mI->last_offsetbase,offset);
    }
    mI->icptr->insert(iI,tmpp);

    ++OC_size;
  }


  void pushOffsets(KK threshold, KK distance){
    IDBUGSTAT(cout << "### " << threshold << "\t" << distance << endl);

    auto mI=locateMapInfo(threshold);

    IDBUGSTAT(cout << "po: found index " << mI-OC_map.begin() << '\n');

    for(;mI!=OC_map.end() && mI->first_offsetbase < threshold; ++mI){
      IDBUGSTAT(cout << "po: working on index " << mI-OC_map.begin());
      IDBUGSTAT(cout << "\nBefore " << *mI << '\n');

      // TODO: not sure this lower bound search actually saves time, need to check
      auto pI=lower_bound(mI->icptr->begin(),mI->icptr->begin(),threshold-mI->first_offsetbase);
      //auto pI=mI->icptr->begin();
      for(; pI!=mI->icptr->end(); ++pI){
	IDBUGSTAT(cout << "ieob: " << ie.offset);
	if(mI->first_offsetbase + pI->offset >= threshold){
	  pI->offset+=distance;
	  IDBUGSTAT(cout << "\tieoa: " << pI->offset);
	}
	IDBUGSTAT(cout << '\n');
      }
      if(mI->last_offsetbase >= threshold) mI->last_offsetbase += distance;
      IDBUGSTAT(cout << "After " << *mI << '\n');
    }
    for(;mI!=OC_map.end(); ++mI){
      IDBUGSTAT(cout << "Pushb: " << *mI);
      mI->first_offsetbase += distance;
      mI->last_offsetbase += distance;
      IDBUGSTAT(cout << "Pusha: " << *mI);
    }
  }

  iterator lowerbound(KK offset){
    auto mI=locateMapInfo(offset);
    if(mI==OC_map.end()) return end();

    auto pI=lower_bound(mI->icptr->begin(),mI->icptr->end(),offset-mI->first_offsetbase);

    // need to find out the n-th element we are
    size_t n=pI-mI->icptr->begin();
    for(auto ocmI=OC_map.begin(); ocmI!=mI; ++ocmI){
      n+=ocmI->icptr->size();
    }

    return iterator(this,
		    n,
		    mI-OC_map.begin(),
		    pI,
		    mI->icptr->begin(),
		    mI->icptr->end());
  }

};


#endif
