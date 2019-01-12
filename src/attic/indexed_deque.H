/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
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


#ifndef _mira_IndexedDeque_h_
#define _mira_IndexedDeque_h_




#include <vector>
#include <algorithm>
#include <iterator>
#include <signal.h>


// This is a drop-in replacement for the STL deque
//
// It replicates most of the functionality of deque,
//  but currently not implemented are:
//   - function variants with InputIterators (constructor, assign(), insert())
//   - operator==()
//   - get_allocator()
//   - probably some iterator comparisons like <, >, etc. (and other ???)
//
// Major difference to deque:
//   - random access: O(1) (albeit using one indirection) instead of O(log N)
//   - insertion/deletion which are not at front or back: like vector O(N/2) on average,
//     but 5 to 6 times faster than deque. Could be improved to O(N/4).
//
// Class is meant to store *vast* amounts of elements (1 million, 10 million or more)
//  but still have constant random access and insertion/deletion which is quicker
//  than deque
//
// Iterator validity: unlike deque, iterators remain valid after insertion/deletion
//  into the IndexedDeque. They may point to deleted items though.
// However, to permit exchange of IndexedDeque with deque, one should not rely on
//  these attributes of IndexedDeque iterators
//
// What a letdown. While good in synthetic tests, real life with MIRA on a 45x 454
//  FLX project used 40%-50% more memory. Memory fragmentation even increased?
// Even with TC malloc there's ~25% more memory used
// One week-end worth of work down the drain. Oh well, at least I found a small
//  bug in contig.C which had gone unnoticed with deque.


//#define IDCEBUG(bla)   {cout << bla; cout.flush();}
//#define IDBUGSTAT(statement) { statement;}
#define IDCEBUG(bla)
#define IDBUGSTAT(statement)

template <typename T>
class IndexedDeque
{
public:

  typedef size_t size_type;

  // Modelled after
  //  http://www.dreamincode.net/forums/topic/58468-making-your-own-iterators/
  // and stl_deque.h from STL
  //
  // STL magic ...
  // using
  //   typedef std::random_access_iterator_tag iterator_category;
  // instead of
  //   typedef std::forward_iterator_tag iterator_category;
  //
  // lets STL algorithms like advance() use iterator operator+=() instead of
  //  repeated operator++()

  class iterator
  {
  public:
    typedef iterator self_type;
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_type difference_type;
    inline iterator(vector<T *> * vbase, size_type offset) : vbase_(vbase), index_(offset) { }
    inline iterator() : vbase_(nullptr), index_(0) { }
    inline self_type operator++() { self_type i = *this; ++index_; return i; }
    inline self_type operator++(int) { ++index_; return *this; }
    inline self_type operator--() { self_type i = *this; --index_; return i; }
    inline self_type operator--(int) { --index_; return *this; }
    inline self_type & operator+=(difference_type n) { index_+=n; return *this; }
    inline self_type operator+(difference_type n) const { self_type tmp=*this; return tmp+=n; }
    inline size_type operator+(iterator rhs) const { return index_+rhs.index_; }
    inline self_type & operator-=(difference_type n) { index_-=n; return *this; }
    inline self_type operator-(difference_type n) const { self_type tmp=*this; return tmp-=n; }
    inline size_type operator-(iterator rhs) const { return index_-rhs.index_; }
    inline reference operator*() {
      IDBUGSTAT(cout << "os:\tvb " << vbase_ << "\tvbs: " << vbase_->size() << "\ti: " << index_ << endl;);
      IDBUGSTAT(if(index_ >= vbase_->size()) { cout << "WARNING! OUT OF BOUNDS." << endl; }
	);
      return *((*vbase_)[index_]);
    }
    inline pointer operator->() {
      IDBUGSTAT(cout << "op:\tvb " << vbase_ << "\tvbs: " << vbase_->size() << "\ti: " << index_ << endl);
      IDBUGSTAT(if(index_ >= vbase_->size()) {cout << "WARNING! OUT OF BOUNDS." << endl; }
	);
      return (*vbase_)[index_];
    }
    inline bool operator==(const self_type& rhs) { return index_ == rhs.index_ && vbase_==rhs.vbase_ ;}
    inline bool operator!=(const self_type& rhs) { return index_ != rhs.index_ || vbase_!=rhs.vbase_ ;}

    inline vector<T *> * getVBase() const {return vbase_;}
    inline size_type getIndex() const {return index_;}

  private:
    size_type index_;
    vector<T *> * vbase_;
  };

  class const_iterator
  {
  public:
    typedef const_iterator self_type;
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef int difference_type;
    inline const_iterator(vector<T *> * vbase, size_type offset) : vbase_(vbase), index_(offset) { }
    inline const_iterator() : vbase_(nullptr), index_(0) { }
    inline const_iterator(iterator it) : vbase_(it.getVBase()), index_(it.getIndex()) { }
    inline self_type operator++() { self_type i = *this; ++index_; return i; }
    inline self_type operator++(int junk) { ++index_; return *this; }
    inline self_type operator--() { self_type i = *this; --index_; return i; }
    inline self_type operator--(int junk) { --index_; return *this; }
    inline self_type & operator+=(difference_type n) { index_+=n; return *this; }
    inline self_type operator+(difference_type n) const { self_type tmp=*this; return tmp+=n; }
    inline size_type operator+(iterator rhs) const { return index_+rhs.index_; }
    inline self_type & operator-=(difference_type n) { index_-=n; return *this; }
    inline self_type operator-(difference_type n) const { self_type tmp=*this; return tmp-=n; }
    inline size_type operator-(iterator rhs) const { return index_-rhs.index_; }

    inline const reference operator*() {
      IDBUGSTAT(cout << "cos:\tvb " << vbase_ << "\tvbs: " << vbase_->size()
	      << "\ti: " << index_ << endl;
	      if(index_ >= vbase_->size()) cout << "WARNING! OUT OF BOUNDS." << endl;);
      return *((*vbase_)[index_]);
    }
    inline pointer operator->() {
      IDBUGSTAT(cout << "cop:\tvb " << vbase_ << "\tvbs: " << vbase_->size()
	      << "\ti: " << index_ << endl;
	      if(index_ >= vbase_->size()) cout << "WARNING! OUT OF BOUNDS." << endl;);
      return (*vbase_)[index_];
    }

    inline bool operator==(const self_type& rhs) { return index_ == rhs.index_ && vbase_==rhs.vbase_ ;}
    inline bool operator!=(const self_type& rhs) { return index_ != rhs.index_ || vbase_!=rhs.vbase_ ;}
    inline size_type getIndex() const {return index_;}
    inline vector<T *> * getVBase() const {return vbase_;}

  private:
    size_type index_;
    vector<T *> * vbase_;
  };



  class reverse_iterator
  {
  public:
    typedef reverse_iterator self_type;
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_type difference_type;
    inline reverse_iterator(vector<T *> * vbase, size_type offset) : vbase_(vbase), index_(offset) { }
    inline reverse_iterator() : vbase_(nullptr), index_(0) { }
    inline self_type operator++() { self_type i = *this; --index_; return i; }
    inline self_type operator++(int junk) { --index_; return *this; }
    inline self_type operator--() { self_type i = *this; ++index_; return i; }
    inline self_type operator--(int junk) { ++index_; return *this; }
    inline self_type & operator+=(difference_type n) { index_-=n; return *this; }
    inline self_type operator+(difference_type n) const { self_type tmp=*this; return tmp-=n; }
    inline size_type operator+(iterator rhs) const { return index_+rhs.index_; }
    inline self_type & operator-=(difference_type n) { index_+=n; return *this; }
    inline self_type operator-(difference_type n) const { self_type tmp=*this; return tmp+=n; }
    inline size_type operator-(iterator rhs) const { return index_-rhs.index_; }
    inline reference operator*() {return *((*vbase_)[index_]); }
    inline pointer operator->() { return (*vbase_)[index_]; }
    inline bool operator==(const self_type& rhs) {
      return index_ == rhs.index_ && vbase_==rhs.vbase_ ;
    }
    inline bool operator!=(const self_type& rhs) {
      return index_ != rhs.index_ || vbase_!=rhs.vbase_ ;
    }
    inline size_type getIndex() const {return index_;}
    inline vector<T *> * getVBase() const {return vbase_;}

  private:
    size_type index_;
    vector<T *> * vbase_;
  };

  class const_reverse_iterator
  {
  public:
    typedef const_reverse_iterator self_type;
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef int difference_type;
    inline const_reverse_iterator(vector<T *> * vbase, size_type offset) : vbase_(vbase), index_(offset) { }
    inline const_reverse_iterator() : vbase_(nullptr), index_(0) { }
    inline const_reverse_iterator(reverse_iterator it) : vbase_(it.getVBase()), index_(it.getIndex()) { }
    inline self_type operator++() { self_type i = *this; --index_; return i; }
    inline self_type operator++(int junk) { --index_; return *this; }
    inline self_type operator--() { self_type i = *this; ++index_; return i; }
    inline self_type operator--(int junk) { ++index_; return *this; }
    inline self_type & operator+=(difference_type n) { index_-=n; return *this; }
    inline self_type operator+(difference_type n) const { self_type tmp=*this; return tmp-=n; }
    inline size_type operator+(iterator rhs) const { return index_+rhs.index_; }
    inline self_type & operator-=(difference_type n) { index_+=n; return *this; }
    inline self_type operator-(difference_type n) const { self_type tmp=*this; return tmp+=n; }
    inline size_type operator-(iterator rhs) const { return index_-rhs.index_; }
    inline const reference operator*() { return *((*vbase_)[index_]); }
    inline pointer operator->() { return (*vbase_)[index_]; }
    inline bool operator==(const self_type& rhs) {
      return index_ == rhs.index_ && vbase_==rhs.vbase_ ;
    }
    inline bool operator!=(const self_type& rhs) {
      return index_ != rhs.index_ || vbase_!=rhs.vbase_ ;
    }
    inline size_type getIndex() const {return index_;}
    inline vector<T *> * getVBase() const {return vbase_;}

  private:
    size_type index_;
    vector<T *> * vbase_;
  };
  reverse_iterator rbegin() {
    return reverse_iterator(&ID_indexvector,ID_indexvector.size()-1);
  }

  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(&ID_indexvector,ID_indexvector.size()-1);
  }

  reverse_iterator rend() {
    return reverse_iterator(&ID_indexvector,ID_iv_freeelementsatstart-1);
  }

  const_reverse_iterator rend() const {
    return reverse_const_iterator(&ID_indexvector,ID_iv_freeelementsatstart-1);
  }



  T& operator[](size_type index) {
    return *(ID_indexvector[index+ID_iv_freeelementsatstart]);
  }

  const T& operator[](size_type index) const {
    return *(ID_indexvector[index+ID_iv_freeelementsatstart]);
  }

  void assign (size_type n, const T& u ){
    IDBUGSTAT(cout << "assign enter\n"; dumpStats(););
    clear();
    if(n>ID_indexvector.capacity()){
      ID_indexvector.reserve(n);
    }
    ID_iv_freeelementsatstart=(ID_indexvector.capacity()-n)/2;
    ID_indexvector.resize(ID_iv_freeelementsatstart);
    for(size_type i=0; i<n; ++i){
      push_back(u);
    }
    IDBUGSTAT(cout << "assign at end\n"; dumpStats(); cout << "assign end" << endl;);
  }

  T& at(size_type index) {
    rangeCheck(index);
    return *(ID_indexvector[index+ID_iv_freeelementsatstart]);
  }

  const T& at(size_type index) const {
    rangeCheck(index);
    return *(ID_indexvector[index+ID_iv_freeelementsatstart]);
  }

  iterator begin() {
    return iterator(&ID_indexvector,ID_iv_freeelementsatstart);
  }

  const_iterator begin() const {
    return const_iterator(&ID_indexvector,ID_iv_freeelementsatstart);
  }

  iterator end() {
    return iterator(&ID_indexvector,ID_indexvector.size());
  }

  const_iterator end() const {
    return const_iterator(&ID_indexvector,ID_indexvector.size());
  }

  T & front() {
    return *(ID_indexvector[0]);
  }

  const T & front() const {
    return *(ID_indexvector[0]);
  }

  T & back() {
    return *(ID_indexvector[ID_indexvector.size()-1]);
  }

  const T & back() const {
    return *(ID_indexvector[ID_indexvector.size()-1]);
  }

  T const & operator=(T const & other){
    if(this != &other){
      IDBUGSTAT(cout << "op= enter\n"; dumpStats(););
      clear();
      ID_indexvector.reserve(other.ID_indexvector.size());
      ID_iv_freeelementsatstart=0;
      typename vector<T*>::iterator vI=other.ID_indexvector.begin();
      for(; vI != other.ID_indexvector.end(); vI++){
	push_back(*vI);
      }
      IDBUGSTAT(cout << "op= at end\n"; dumpStats(); cout << "op= end" << endl;);
    }
  }



  IndexedDeque() : ID_iv_freeelementsatstart(0) {}

  ~IndexedDeque() {
    IDBUGSTAT(cout << "~ enter\n"; dumpStats(););
    clear();
    IDBUGSTAT(cout << "~ at end\n"; dumpStats(); cout << "~ end" << endl;);
  }

  inline void clear() {
    IDBUGSTAT(cout << "clear enter\n"; dumpStats(););
//    typename vector<T*>::reverse_iterator ivRI = ID_indexvector.rbegin();
//    size_type i=0;
//    for(; ivRI!=ID_indexvector.rend(); ivRI++){
    typename vector<T*>::iterator ivRI = ID_indexvector.begin();
    size_type i=0;
    for(; ivRI!=ID_indexvector.end(); ivRI++){
      delete *ivRI;
      ++i;
      if(i==ID_indexvector.size()-ID_iv_freeelementsatstart){
	break;
      }
    }
    ID_iv_freeelementsatstart=ID_indexvector.size()/2;
    ID_indexvector.resize(ID_iv_freeelementsatstart);
    IDBUGSTAT(cout << "clear at end\n"; dumpStats(); cout << "clear end" << endl;);
  }

  inline bool empty() {
    return (ID_indexvector.size()-ID_iv_freeelementsatstart == 0);
  }

  iterator erase(const iterator & it){
    IDBUGSTAT(cout << "erase it enter\n"; dumpStats(););
    delete ID_indexvector[it.getIndex()];

    typename vector<T *>::iterator tmpvI=ID_indexvector.begin()+it.getIndex();
    ID_indexvector.erase(tmpvI);
    IDBUGSTAT(cout << "erase it at end\n"; dumpStats(); cout << "erase it end" << endl;);
    return it;
  }

  iterator erase(const iterator & itf,const iterator & itl){
    IDBUGSTAT(cout << "erase it it enter\n"; dumpStats(););
    for(iterator ti=itf; ti!=itl; ++ti){
      delete ID_indexvector[ti.getIndex()];
    }
    typename vector<T *>::iterator tmpvfI=ID_indexvector.begin()+itf.getIndex();
    typename vector<T *>::iterator tmpvlI=ID_indexvector.begin()+itl.getIndex();
    ID_indexvector.erase(tmpvfI,tmpvlI);
    IDBUGSTAT(cout << "erase it it at end\n"; dumpStats(); cout << "erase it it end" << endl;);
    return itf;
  }

  inline iterator insert(const iterator & it, const T & elem){
    if(it.getIndex()-ID_iv_freeelementsatstart==0){
      push_front(elem);
    }else if(it.getIndex()==ID_indexvector.size()){
      push_back(elem);
    }else{
      IDBUGSTAT(cout << "insert middle enter\n"; dumpStats(););
      T * ne=new T(elem);

      typename vector<T *>::iterator tmpvI=ID_indexvector.begin()+it.getIndex();
      ID_indexvector.insert(tmpvI,ne);
      IDBUGSTAT(cout << "insert middle at end\n"; dumpStats(); cout << "insert middle end" << endl;);
    }
    return it;
  }

  inline iterator insert(const iterator & it, size_type n, const T & elem){
    if(it.getIndex()-ID_iv_freeelementsatstart==0){
      for(size_type i=0; i<n; ++i){
	push_front(elem);
      }
    }else if(it.getIndex()==ID_indexvector.size()){
      for(size_type i=0; i<n; ++i){
	push_back(elem);
      }
    }else{
      IDBUGSTAT(cout << "insert middle n enter\n"; dumpStats(););
      typename vector<T *>::iterator tmpvI=ID_indexvector.begin()+it.getIndex();
      ID_indexvector.insert(tmpvI,n,nullptr);
      tmpvI=ID_indexvector.begin()+it.getIndex();
      for(size_type i=0; i<n; ++i, ++tmpvI){
	*tmpvI=new T(elem);
      }
      IDBUGSTAT(cout << "insert middle n at end\n"; dumpStats(); cout << "insert middle n end" << endl;);
    }
    return it;
  }

  size_type max_size () const {return static_cast<size_type>(-1);}

  inline void push_back(const T & elem) {
    IDBUGSTAT(cout << "push_back enter\n"; dumpStats(););
    T * ne=new T(elem);

    if(ID_indexvector.capacity()==ID_indexvector.size()){
      ID_indexvector.reserve(calcAdditionalCapacity());
    }
    ID_indexvector.push_back(ne);
    IDBUGSTAT(cout << "push_back at end\n"; dumpStats(); cout << "push_back end" << endl;);
  }

  void pop_back() {
    IDBUGSTAT(cout << "pop_back enter\n"; dumpStats(););
    delete ID_indexvector.back();
    ID_indexvector.pop_back();
    IDBUGSTAT(cout << "pop_back at end\n"; dumpStats(); cout << "pop_back end" << endl;);
  }

  inline void push_front(const T & elem) {
    IDBUGSTAT(cout << "push_front enter\n"; dumpStats(););

    T * ne=new T(elem);

    if(ID_iv_freeelementsatstart==0){
      ID_iv_freeelementsatstart=calcAdditionalCapacity();
      ID_indexvector.insert(ID_indexvector.begin(),ID_iv_freeelementsatstart,nullptr);
    }
    --ID_iv_freeelementsatstart;
    ID_indexvector[ID_iv_freeelementsatstart]=ne;
    IDBUGSTAT(cout << "push_front at end\n"; dumpStats(); cout << "push_front end" << endl;);
  }

  void pop_front() {
    IDBUGSTAT(cout << "pop_front enter\n"; dumpStats(););
    delete ID_indexvector[ID_iv_freeelementsatstart++];
    IDBUGSTAT(cout << "pop_front at end\n"; dumpStats(); cout << "pop_front end" << endl;);
  }

  inline void resize(size_type sz, T c = T()) {
    IDBUGSTAT(cout << "resize enter\n"; dumpStats(););
    if(ID_indexvector.size()-ID_iv_freeelementsatstart < sz){
      for(; ID_indexvector.size()-ID_iv_freeelementsatstart < sz;){
	push_back(c);
      }
    }else{
      for(; ID_indexvector.size()-ID_iv_freeelementsatstart > sz;){
	pop_back();
      }
    }
    IDBUGSTAT(cout << "resize at end\n"; dumpStats(); cout << "resize end" << endl;);
  }

  size_type size() const { return ID_indexvector.size()-ID_iv_freeelementsatstart; }

  void swap(IndexedDeque<T> & other){
    if(this != &other){
      IDBUGSTAT(cout << "swap enter\n"; dumpStats(););
      std::swap(ID_iv_freeelementsatstart, other.ID_iv_freeelementsatstart);
      ID_indexvector.swap(other.ID_indexvector);
      IDBUGSTAT(cout << "swap at end\n"; dumpStats(); cout << "swap end" << endl;);
    }
  }



private:
  vector<T*> ID_indexvector;
  size_type ID_iv_freeelementsatstart;

  inline size_type calcAdditionalCapacity() {
    return max(static_cast<size_t>(1),(ID_indexvector.size()-ID_iv_freeelementsatstart)/2);
  }

  void rangeCheck(size_type n){
    if(n>=size()){
      __throw_out_of_range("IndexedDeque::rangeCheck");
    }
  }

  void dumpStats(){
    cout << "ds: s " << ID_indexvector.size()
	 << "\tID_iv " << ID_iv_freeelementsatstart
	 << "\tc " << ID_indexvector.capacity()
	 << endl;
  }
};


#endif
