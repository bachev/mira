/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2012 and later by Bastien Chevreux
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


#ifndef _mira_hdeque_h_
#define _mira_hdeque_h_

#include <deque>
#include <list>
#include <vector>

#include <boost/iterator_adaptors.hpp>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#ifndef likely
#if __GNUC__ >= 4
#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)
#else
#define likely(x) (x)
#define unlikely(x) (x)
#endif
#endif

template <class TT>
class HDeque
{
private:

  typedef std::deque<TT> innercontainer_t ;

  struct mapinfo_t{
    int64_t  from;
    typename std::list<innercontainer_t>::iterator icI;       // iterator back to rposbin_t in list

    inline mapinfo_t(size_t afrom, typename std::list<innercontainer_t>::iterator I) : from(afrom),icI(I) {};
    static inline bool lt_indexfrom(const mapinfo_t & a, const mapinfo_t & b) { return a.from<b.from;}
  };

/*************************************************************************
 *
 * bins and map are always the same size and never(!) empty
 *
 *************************************************************************/

  std::list<innercontainer_t>  HD_bins;
  std::deque<mapinfo_t>        HD_map;

  size_t HD_bo_binsize;
  size_t HD_size;


  /*************************************************************************************
   *
   *
   *  Iterators
   *
   *
   *************************************************************************************/
public:
  template <class HDIValue>
  class hditer
    : public boost::iterator_facade<
      hditer<HDIValue>                          // Derived
    , HDIValue             // Value
    , std::random_access_iterator_tag  //bidirectional_iterator_tag   //random_access_iterator_tag  //boost::forward_traversal_tag    // CategoryOrTraversal
    > {

  private:
    HDeque * hdptr;
    innercontainer_t * icptr;
    typename innercontainer_t::iterator icI;
    size_t mapindex;

  public:
    inline hditer(HDeque * ahdptr, innercontainer_t * aicptr, typename innercontainer_t::iterator aicI, size_t amapindex) : hdptr(ahdptr),icptr(aicptr),icI(aicI),mapindex(amapindex) { };

    template <class OType>
    hditer(hditer<OType> const & other) : hdptr(other.hdptr), icptr(other.icptr), icI(other.icI), mapindex(other.mapindex) {}

    friend std::ostream & operator<<(std::ostream &ostr, hditer & mI) {
      ostr << "hdptr: " << mI.hdptr
	   << "\ticptr: " << mI.icptr
	   << "\ticsize: " << mI.icptr->size()
	   << "\ticIpos: " << mI.icI-mI.icptr->begin()
	   << "\tmapindex: " << mI.mapindex
	   << std::endl;
      if(mI==mI.hdptr->end()){
	ostr << "This is end()" << std::endl;
      }else{
	//ostr << *this;
      }
      return ostr;
    }


  private:
    // Everything Boost's iterator facade needs
    friend class boost::iterator_core_access;

    // friend with itself means that "iterator" and "const_iterator" can access each other
    template <class> friend class hditer;

    // friend with parent HDeque class so that some quick checks can be performed
    template <class> friend class HDeque;

    HDIValue & dereference() const
    {
      // TODO: eventually remove the BUGIFTHROW if prooves to be to slow (check that
      //  on projects). Current tests indicate it's really minimal, e.g., in the order of
      //  a few seconds for a mapping project
      FUNCSTART("inline TT & HDeque::hditer::dereference() const");
      BUGIFTHROW(*this==hdptr->end(),"Trying to dereference an hditer pointing to end()???");
      return *icI;
    }

    template <class OType>
    inline bool equal(hditer<OType> const & other) const
    {
      //std::cout << "equal? " << mapindex << "\t" << other.mapindex << std::endl;
      return mapindex == other.mapindex && icI == other.icI;
    }

    inline void increment(){
      //if(*this!=hdptr->end()){       // this is 10% slower than the version below
      if(likely(mapindex<hdptr->HD_map.size())){    // likely() brings ~ 2%
	++icI;
	if(unlikely(icI==icptr->end())){
	  ++mapindex;
	  if(unlikely(mapindex==hdptr->HD_map.size())){
	    // special treatment: we're at end
	    *this=hdptr->end();
	  }else{
	    icptr=&(*(hdptr->HD_map[mapindex].icI));
	    icI=icptr->begin();
	  }
	}
      }
    }

    inline void decrement(){
      // are we at end? special treatment
      if(likely(icI != icptr->begin() && mapindex<hdptr->HD_map.size())){
	// simple case: not at start of an inner container and not at end()
	--icI;
      }else{
	if(likely(mapindex>0)){
	  // normal case
	  --mapindex;
	  icptr=&(*(hdptr->HD_map[mapindex].icI));
	  icI=icptr->end();
	  --icI;
	}else{
	  // hmmm, we were at beginning of container anyway
	  // let's do nothing, that should break the least things
	}
      }
    }

    int64 distance_to(hditer const & other) const {
      return (other.hdptr->priv_getPosOfIterator(other) - hdptr->priv_getPosOfIterator(*this));
    }

    // currently used by advance() for larger jumps
    void setToPos(size_t pos){
      if(unlikely(pos>=hdptr->HD_size)) {
	*this=hdptr->end();
      }else{
	int64 ipos=static_cast<int64>(pos);
	ipos+=hdptr->HD_map.front().from;
	auto mI=lower_bound(hdptr->HD_map.begin(),
			    hdptr->HD_map.end(),
			    mapinfo_t(ipos, hdptr->HD_bins.begin()),  // comparator object, only ipos matters
			    mapinfo_t::lt_indexfrom);
	// move the returned map iterator to the right position
	if(likely(mI != hdptr->HD_map.begin()
		  && (mI == hdptr->HD_map.end()
		      || mI->from > ipos))){
	     --mI;
	}
	// and now initialise this iterator to now values
	icptr=&(*(mI->icI));
	icI=mI->icI->begin()+(ipos-mI->from);
	mapindex=mI-hdptr->HD_map.begin();
      }
    }

  public:
    void advance(int64 dist){
      if(abs(dist)>hdptr->HD_bo_binsize){
	int64 thispos=hdptr->priv_getPosOfIterator(*this);
	thispos+=dist;
	setToPos(static_cast<size_t>(thispos));
      }else if(dist>0){
	// TODO: perhaps remove the [] operators
	increment();
	--dist;
	if(dist){
	  bool icivalid=true;
	  if(dist >= icptr->end()-icI){
	    dist-=icptr->end()-icI;
	    ++mapindex;
	    if(mapindex>=hdptr->HD_map.size()){
	      icI=hdptr->HD_bins.front().end();
	      return;
	    }
	    icptr=&(*(hdptr->HD_map[mapindex].icI));
	    icivalid=false;
	  }
	  while(dist>=static_cast<int64>(icptr->size())){
	    dist-=icptr->size();
	    ++mapindex;
	    if(mapindex>=hdptr->HD_map.size()){
	      icI=hdptr->HD_bins.front().end();
	      return;
	    }
	    icptr=&(*(hdptr->HD_map[mapindex].icI));
	    icivalid=false;
	  }
	  if(!icivalid) icI=icptr->begin();
	  icI+=dist;
	}
      }else if(dist<0){
	// TODO: perhaps remove the [] operators
	// hack to get end() play nicely
	decrement();
	dist=-1-dist;
	// rest
	if(dist){
	  bool icivalid=true;
	  if(dist>icI-icptr->begin()){
	    dist-=icI-icptr->begin()+1;
	    if(mapindex==0){
	      icI=icptr->begin();
	      return;
	    }
	    --mapindex;
	    icptr=&(*(hdptr->HD_map[mapindex].icI));
	    icivalid=false;
	  }
	  while(dist>=static_cast<int64>(icptr->size())){
	    if(mapindex==0){
	      icI=icptr->begin();
	      return;
	    }
	    dist-=icptr->size();
	    --mapindex;
	    icptr=&(*(hdptr->HD_map[mapindex].icI));
	    icivalid=false;
	  }
	  if(!icivalid){
	    icI=icptr->end();
	    --icI;
	  }
	  icI-=dist;
	}
      }
    }


  };


  // pure magic to construct "iterator" and "const_iterator" :-)
  typedef hditer<TT> iterator;
  typedef hditer<TT const> const_iterator;


private:


//  vector<offsettile_t>::iterator searchOffsetTileForPlacement(int32 position);
//  //std::list<PlacedContigReads::rposbin_t>::iterator findReadpositionIterator(int32 position);
//
//  void splitBin(uint32 binindex);
//  std::list<rposbin_t>::iterator placeRead_helper(int32 rpid,
//					     int32 position,
//					     int8 direction,
//					     int32 urdid,
//					     iterator::raoindex_t & araoindex);

  inline void priv_checksize(){
    size_t csize=0;
    for(auto & b : HD_bins){
      csize+=b.size();
    }
    if(csize!=HD_size){
      std::cout << "csize " << csize << " != HD_size " << HD_size << std::endl;
      exit(1);
    }
  }

  template <class OType>
  size_t priv_getPosOfIterator(hditer<OType> const & hdI){
    size_t ipos=HD_size;
    if(hdI!=end()){
      ipos=static_cast<size_t>(HD_map[hdI.mapindex].from + (hdI.icI-hdI.icptr->begin()) - HD_map.front().from);
    }
    return ipos;
  }

  inline void priv_resize_shrink(size_t newsize){
    if(newsize<HD_size){
      if(newsize==0) {
	clear();
      }else{
	// reduce container
	// note that bins can have different sizes
	//
	// first, try to get rid of complete bins
	auto diff=HD_size-newsize;
	while(diff>=HD_bins.back().size()){
	  diff-=HD_bins.back().size();
	  HD_size-=HD_bins.back().size();
	  HD_bins.pop_back();
	  HD_map.pop_back();
	}
	// now reduce last bin
	// last bin size will never drop to zero, else it would have
	//  been handled above
	auto fI=HD_bins.back().end();
	HD_size-=diff;
	for(;diff; --fI, --diff) {};
	HD_bins.back().erase(fI,HD_bins.back().end());
      }
    }
  }
  inline void priv_resize_grow_back(size_t newsize, const TT & x) {
    int64 diff=static_cast<int64>(newsize-HD_size);
    while(diff>0){
      if(HD_bins.back().size()>=HD_bo_binsize){
	auto lastsize=HD_bins.back().size();
	HD_bins.push_back(innercontainer_t());
	HD_map.push_back(mapinfo_t(HD_map.back().from+lastsize,--(HD_bins.end())));
      }
      auto numins=HD_bo_binsize-HD_bins.back().size();
      if(numins>diff) numins=diff;
      diff-=numins,
      HD_bins.back().resize(HD_bins.back().size()+numins,x);
    }
    HD_size=newsize;
  }
  inline void priv_resize_grow_front(size_t newsize, const TT & x) {
    auto diff=newsize-HD_size;
    // fill up first bin first
    if(HD_bins.front().size()<HD_bo_binsize){
      auto numins=HD_bo_binsize-HD_bins.front().size();
      if(numins>diff) numins=diff;
      HD_size+=numins;
      diff-=numins;
      HD_map.front().from-=numins;
      HD_bins.front().insert(HD_bins.front().begin(),numins,x);
    }
    // fill up rest of bins using resize()
    while(diff){
      HD_bins.push_front(innercontainer_t());
      HD_map.push_front(mapinfo_t(HD_map.front().from,HD_bins.begin()));
      auto numins=HD_bo_binsize;
      if(diff<HD_bo_binsize) numins=diff;
      diff-=numins;
      HD_size+=numins;
      HD_map.front().from-=numins;
      HD_bins.front().resize(numins,x);
    }
  }

  // for MIRA and the intended use in the CON_counts container, we cheat big time on insert:
  //  insertion at begin() / end() is done via resize_grow (which is OK), but there is
  //  currently no "re-distribution" of elements to bins for insertions in the middle of a container
  // It's probably not necessary, because how many bases will be inserted amidst a contig? Even
  //  for IonTorrent at very high coverages this will perhaps be 1 base every two or three bases,
  //  therefore the inner bins should not really grow much more than double the bin size set
  //  for the HDeque
  iterator priv_insert(iterator where, size_t num, const TT & x){
    FUNCSTART("HDeque::insert(iterator where, size_t num, const TT & x)");
    BUGIFTHROW(where.hdptr!=this,"called with an iterator not belonging to this container?");
    auto wpos=priv_getPosOfIterator(where);
    //std::cout << "where " << where << "\tnum: " << num << "\twpos: " << wpos << std::endl;
    if(wpos>=HD_size){
      priv_resize_grow_back(HD_size+num,x);
    }else if(where==where.hdptr->begin()){
      priv_resize_grow_front(HD_size+num,x);
    }else{
      where.icptr->insert(where.icI,num,x);
      // just need to adjust the map.from values
      auto mI=HD_map.begin()+(where.mapindex+1);
      for(; mI!=HD_map.end(); ++mI){
	mI->from+=num;
      }
      HD_size+=num;
    }
    return (begin()+wpos);
  }


  iterator priv_erase(iterator from, iterator to){
    FUNCSTART("iterator priv_erase(iterator from, iterator to)");
    BUGIFTHROW(from.hdptr!=this,"where called with an iterator not belonging to this container?");
    BUGIFTHROW(to.hdptr!=this,"to called with an iterator not belonging to this container?");
    auto frompos=priv_getPosOfIterator(from);
    int64 diff=static_cast<int64>(priv_getPosOfIterator(to)-frompos);
    BUGIFTHROW(diff<0,"diff < 0 ... 'from' was > 'to' ???");
    BUGIFTHROW(diff>static_cast<int64>(HD_size),"diff>HD_size");
    HD_size-=diff;
    //std::cout << "diff: " << diff << "\tHD_size " << HD_size << std::endl;
    if(HD_size==0){
      clear();
    }else{
      // erase in first container
      // tempmapindex needed to reach the correct "first middle container"
      //  in case the first container also gets completely erased
      auto tempmapindex=from.mapindex;
      if(diff){
	int64 maxerase=from.icptr->end()-from.icI;
	size_t thiserase=std::min(diff,maxerase);
	auto ictoI=from.icI+thiserase;
	from.icptr->erase(from.icI,ictoI);
	diff-=thiserase;
	if(from.icptr->empty()){
	  HD_bins.erase((HD_map.begin()+tempmapindex)->icI);
	  HD_map.erase(HD_map.begin()+tempmapindex);
	  --tempmapindex;
	}
      }
      // erase in middle and end containers
      if(diff){
	++tempmapindex;
	auto mI=HD_map.begin()+tempmapindex;
	while(diff){
	  size_t thiserase=std::min(static_cast<int64>(diff),static_cast<int64>(mI->icI->size()));
	  auto ictoI=mI->icI->begin()+thiserase;
	  mI->icI->erase(mI->icI->begin(),ictoI);
	  diff-=thiserase;
	  if(mI->icI->empty()){
	    HD_bins.erase(mI->icI);
	    mI=HD_map.erase(mI);
	  }else{
	    ++mI;
	  }
	}
      }
      // shift down all "from" values in the map which need adjustment
      auto mI=HD_map.begin();
      int64 lastfrom=0;
      int64 lastsize=0;
      if(from.mapindex){
	advance(mI,from.mapindex-1);
	lastfrom=mI->from;
	lastsize=mI->icI->size();
	++mI;
      }
      for(; mI != HD_map.end(); ++mI){
	mI->from=lastfrom+lastsize;
	lastfrom=mI->from;
	lastsize=mI->icI->size();
      }
    }
    from.setToPos(frompos);
    return from;
  }


public:
  //HDeque() : HD_bo_binsize(2048), HD_size(0) { HD_bins.resize(1); HD_map.push_back(mapinfo_t(0,HD_bins.begin()));};
  HDeque() : HD_bo_binsize(8192), HD_size(0) { HD_bins.resize(1); HD_map.push_back(mapinfo_t(0,HD_bins.begin()));};
  ~HDeque() {} ;

  inline size_t size() const {return HD_size;};
  inline bool empty() const {return HD_size==0;};

  // Copy constructor & copy operator
  HDeque(const HDeque & other){
    *this=other;
  }
  HDeque const & operator=(HDeque const & other) {
    HD_bins=other.HD_bins;
    HD_bo_binsize=other.HD_bo_binsize;
    HD_size=other.HD_size;
    HD_map=other.HD_map;
    // of course, need to adapt iterators in map
    auto iclI=HD_bins.begin();
    auto mI=HD_map.begin();
    for(; iclI!=HD_bins.end(); ++iclI, ++mI){
      mI->icI=iclI;
    }
    return *this;
  }
  inline void clear() {
    HD_bins.resize(1);
    HD_bins.begin()->clear();
    HD_map.clear();
    HD_map.push_back(mapinfo_t(0,HD_bins.begin()));
    HD_size=0;
  }
  inline void swap(HDeque & other) {
    HD_bins.swap(other.HD_bins);
    HD_map.swap(other.HD_map);
    std::swap(HD_bo_binsize,other.HD_bo_binsize);
    std::swap(HD_size,other.HD_size);
  }

  inline void resize(size_t newsize, const TT & x) {
    priv_checksize();
    if(newsize<=HD_size){
      priv_resize_shrink(newsize);
    }else{
      priv_resize_grow_back(newsize,x);
    }
  }
  inline void resize(size_t newsize) {
    priv_checksize();
    if(newsize<=HD_size){
      priv_resize_shrink(newsize);
    }else{
      priv_resize_grow_back(newsize,TT());
    }
  }


  void setBinSize(uint32 bs) {
    if(bs<2) throw std::out_of_range("internal error HDeque::setBinSize: size may not be < 2");
    HD_bo_binsize=bs;
  };

  void debugDump(bool shortdbg) {
    std::cout << "HD_bins.size(): " << HD_bins.size()
	 << "\tHD_map.size(): " << HD_map.size()
	 << "\tHD_bo_binsize: " << HD_bo_binsize
	 << "\tHD_size: " << HD_size << std::endl;
    if(!shortdbg){
      for(size_t imi=0; imi<HD_map.size(); ++imi){
	std::cout << "Bin " << imi << ": from " << HD_map[imi].from << "\tsize " << HD_map[imi].icI->size() << std::endl;
      }
    }
  }
  void debugDumpDeep() {
    debugDump(false);
    for(size_t imi=0; imi<HD_map.size(); ++imi){
      std::cout << "Bin " << imi << ": from " << HD_map[imi].from << "\tsize " << HD_map[imi].icI->size() << std::endl;
      uint32 l=0;
      for(auto & x : *HD_map[imi].icI){
	std::cout << "B " << imi << " " << l << " ";
	++l;
	std::cout << x << std::endl;
      }
    }
  }

  inline iterator begin() {
    return iterator(this,
		    &(*(HD_bins.begin())),
		    HD_bins.begin()->begin(),
		    HD_size==0); // is == 1 for empty container (which allows for begin==end()), else == 0
  }
  inline iterator end() {
    return iterator(this,&(*(HD_bins.begin())),HD_bins.begin()->end(),HD_map.size());
  }
  inline const_iterator cbegin() {
    return const_iterator(this,
		    &(*(HD_bins.begin())),
		    HD_bins.begin()->begin(),
		    HD_size==0); // is == 1 for empty container (which allows for begin==end()), else == 0
  }
  inline const_iterator cend() {
    return const_iterator(this,&(*(HD_bins.begin())),HD_bins.begin()->end(),HD_map.size());
  }


  void push_back(const TT & x){
    if(HD_map.back().icI->size()>=HD_bo_binsize){
      HD_bins.push_back(innercontainer_t());
      HD_map.push_back(mapinfo_t(HD_map.back().from+HD_bo_binsize,--(HD_bins.end())));
    }
    HD_map.back().icI->push_back(x);
    ++HD_size;
  }
  void pop_back(){
    if(HD_map.back().icI->size()>1 || HD_size==1){
      HD_map.back().icI->pop_back();
    }else{
      HD_map.pop_back();
      HD_bins.pop_back();
    }
    --HD_size;
  }

  void push_front(const TT & x){
    if(HD_map.front().icI->size()>=HD_bo_binsize){
      HD_bins.push_front(innercontainer_t());
      HD_map.push_front(mapinfo_t(HD_map.front().from,HD_bins.begin()));
    }
    HD_map.front().icI->push_front(x);
    --HD_map.front().from;
    ++HD_size;
  }
  void pop_front(){
    if(HD_map.front().icI->size()>1 || HD_size==1){
      HD_map.front().icI->pop_front();
      ++HD_map.front().from;
    }else{
      HD_map.pop_front();
      HD_bins.pop_front();
    }
    --HD_size;
  }

  // front() and back() of course crash when used on empty container, but so does
  //  vectory, deque etc.pp
  TT & back(){
    return HD_bins.back().back();
  }
  const TT & back() const{
    return HD_bins.back().back();
  }
  TT & front(){
    return HD_bins.front().front();
  }
  const TT & front() const{
    return HD_bins.front().front();
  }

  void insert(iterator where, size_t num, const TT & x){
    priv_insert(where, num, x);
    return;
  }
  iterator insert(iterator where, const TT & x){
    return priv_insert(where,1,x);
  }

  iterator erase(iterator from, iterator to){
    return priv_erase(from,to);
  }

  iterator erase(iterator elemI){
    FUNCSTART("HDeque::erase(iterator elemI)");
    BUGIFTHROW(elemI==end(),"elemi==end()?");
    return erase(elemI,elemI+1);
  }

};


#endif
