#include <typeinfo>
#include <string>

#include "stdinc/types.H"
#include "stdinc/stlincludes.H"

#include <list>
#include <map>
#include <set>

#ifndef _memusage_h
#define _memusage_h

#define MMEBUG(bla)



/*************************************************************************
 *
 * Calculation routines
 *
 *************************************************************************/

// for strings
// this is a real cludge
inline void emuoc_statsOnType(const std::string & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer);

// for simple types (int, structs, unions, objects etc)
template<class myType>
inline void emuoc_statsOnType(const myType & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer);

// for containers without capacity
template<class myType>
inline void emuoc_statsOnType_sizeonly(const myType & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer, size_t overhead_elem);

// for containers with capacity (at the moment only vector?)
template<class myType>
inline void emuoc_statsOnType_sizeAndCapacity(const myType & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer);


/*************************************************************************
 *
 * Specialised hooks
 *
 *************************************************************************/

template<class myType>
inline void emuoc_statsOnType(const std::vector<myType> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{

  emuoc_statsOnType_sizeAndCapacity(
    c,total_size,total_freecap, total_lostba, num_size,
    bytes_size, overhead, lostbyalign, iscontainer);
  // overhead of each elem is 0
}

template<class myType>
inline void emuoc_statsOnType(const std::vector<std::string> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{

  emuoc_statsOnType_sizeAndCapacity(
    c,total_size,total_freecap, total_lostba, num_size,
    bytes_size, overhead, lostbyalign, iscontainer);
  // overhead of each elem is 0
}

template<class myType>
inline void emuoc_statsOnType(const std::list<myType> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{

  emuoc_statsOnType_sizeonly(c,total_size,total_freecap, total_lostba, num_size,
		bytes_size, overhead, lostbyalign, iscontainer,
		2*sizeof(void *));
  // overhead of each elem is a forward and a reverse pointer
}


//template<class myType>
//inline void emuoc_statsOnType(const slist<myType> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
//{
//
//  emuoc_statsOnType_sizeonly(c,total_size,total_freecap, total_lostba, num_size,
//		bytes_size, overhead, lostbyalign, iscontainer,
//		sizeof(void *));
//  // overhead of each elem is a forward pointer
//}

template<class myType>
inline void emuoc_statsOnType(const std::set<myType> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{

  emuoc_statsOnType_sizeonly(c,total_size,total_freecap, total_lostba, num_size,
		bytes_size, overhead, lostbyalign, iscontainer,
		4*sizeof(void *));
  // overhead of each elem is
  //  - the tree color (padded to fill up void *)
  //  - a parent, a left and a right child
}

template<class myType>
inline void emuoc_statsOnType(const std::multiset<myType> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{

  emuoc_statsOnType_sizeonly(c,total_size,total_freecap, total_lostba, num_size,
		bytes_size, overhead, lostbyalign, iscontainer,
		4*sizeof(void *));
  // overhead of each elem is
  //  - the tree color (padded to fill up void *)
  //  - a parent, a left and a right child
}


template<class myType1,class myType2>
inline void emuoc_statsOnType(const std::map<myType1,myType2> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{

  emuoc_statsOnType_sizeonly(c,total_size,total_freecap, total_lostba, num_size,
		bytes_size, overhead, lostbyalign, iscontainer,
		4*sizeof(void *));
  // overhead of each elem is
  //  - the tree color (padded to fill up void *)
  //  - a parent, a left and a right child
}

template<class myType1,class myType2>
inline void emuoc_statsOnType(const std::multimap<myType1,myType2> & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{
  emuoc_statsOnType_sizeonly(c,total_size,total_freecap, total_lostba, num_size,
		bytes_size, overhead, lostbyalign, iscontainer,
		4*sizeof(void *));
  // overhead of each elem is
  //  - the tree color (padded to fill up void *)
  //  - a parent, a left and a right child
}




template<class myType>
size_t estimateMemoryUsageOfContainer(const myType & c, bool withoverhead, size_t & numelem, size_t & bytes_size, size_t & free_capacity, size_t & lostbyalign)
{
  MMEBUG("dido called with: " << typeid(myType).name() << endl);

  numelem=0;
  bytes_size=0;
  free_capacity=0;
  lostbyalign=0;

  {
    size_t stats_s, stats_o, stats_lba;
    bool stats_isc;

    // specifying <myType> in the function call helps the compiler
    //  to resolve to the correct objetc type and only switch to
    //  generic if everything else fails
    emuoc_statsOnType(c,numelem,free_capacity,lostbyalign,
		stats_s,bytes_size,stats_o,stats_lba,stats_isc);

    if(!withoverhead) {
      bytes_size-=sizeof(myType);
    }

    MMEBUG("Stats total size: " << numelem << endl);
    MMEBUG("Stats byte size: " << bytes_size << endl);
    MMEBUG("Stats total free cap: " << free_capacity << endl);
    MMEBUG("Stats total lba: " << lostbyalign << endl);

    MMEBUG("Stats size: " << stats_s << endl);
    MMEBUG("Stats overhead: " << stats_o << endl);
    MMEBUG("Stats lost b align: " << stats_lba << endl);
    MMEBUG("Stats is container: " << stats_isc << endl);
  }

  return bytes_size;
}


template<class myType>
size_t estimateMemoryUsageOfContainer(const myType & c, bool withoverhead)
{
  MMEBUG("dido called with: " << typeid(myType).name() << endl);

  size_t numelem, bytes_size, free_capacity, lostbyalign;

  numelem=0;
  bytes_size=0;
  free_capacity=0;
  lostbyalign=0;

  {
    size_t stats_s, stats_o, stats_lba;
    bool stats_isc;

    // specifying <myType> in the function call helps the compiler
    //  to resolve to the correct objetc type and only switch to
    //  generic if everything else fails
    emuoc_statsOnType(c,numelem,free_capacity,lostbyalign,
		stats_s,bytes_size,stats_o,stats_lba,stats_isc);

    if(!withoverhead) {
      bytes_size-=sizeof(myType);
    }

    MMEBUG("Stats total size: " << numelem << endl);
    MMEBUG("Stats byte size: " << bytes_size << endl);
    MMEBUG("Stats total free cap: " << free_capacity << endl);
    MMEBUG("Stats total lba: " << lostbyalign << endl);

    MMEBUG("Stats size: " << stats_s << endl);
    MMEBUG("Stats overhead: " << stats_o << endl);
    MMEBUG("Stats lost b align: " << stats_lba << endl);
    MMEBUG("Stats is container: " << stats_isc << endl);
  }

  return bytes_size;
}

//#define MMEBUG(bla)





/*************************************************************************
 *
 * Calculation routines
 *
 *************************************************************************/

// for strings
// this is a real cludge
inline void emuoc_statsOnType(const std::string & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{
  MMEBUG("emuoc_statsOnType(string): " << typeid(c).name() << endl);

  (void) c;
  (void) total_freecap;
  (void) total_lostba;

  total_size++;
  num_size=1;
  // at least in gcc we seem to have:
  //  one pointer
  overhead=sizeof(c);
  //  and if the string has capacity, 3 class variables
  //  (round to 4 for alignment)
  if(c.capacity()>0) overhead+=4 * sizeof(size_t);
  bytes_size=overhead+c.capacity();
  iscontainer=true;

  lostbyalign=0;
  if(bytes_size%sizeof(void *)){
    lostbyalign=sizeof(void *)-bytes_size%sizeof(void *);
  }
  MMEBUG("Lost by align: " <<lostbyalign << endl);

  return;
}



// for simple types (int, structs, unions, objects etc)
template<class myType>
inline void emuoc_statsOnType(const myType & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{
  MMEBUG("emuoc_statsOnType(generic): " << typeid(myType).name() << endl);

  (void) c;
  (void) total_freecap;
  (void) total_lostba;

  total_size++;
  num_size=1;
  bytes_size=sizeof(c);
  overhead=0;
  iscontainer=false;

  lostbyalign=0;
  if(bytes_size%sizeof(void *)){
    lostbyalign=sizeof(void *)-bytes_size%sizeof(void *);
  }
  MMEBUG("Lost by align: " <<lostbyalign << endl);

  return;
}


// for containers without capacity
template<class myType>
inline void emuoc_statsOnType_sizeonly(const myType & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer, size_t overhead_elem)
{
  MMEBUG("emuoc_statsOnType(simple): " << typeid(std::list<myType>).name() << endl);

  (void) c;
  (void) total_freecap;
  num_size=0;
  overhead=sizeof(myType);
  bytes_size=overhead;
  lostbyalign=0;
  iscontainer=true;

  typename myType::const_iterator cI=c.begin();
  for(; cI != c.end(); cI++){
    size_t stats_s, stats_bs, stats_o, stats_lba;
    bool stats_isc;
    MMEBUG("loop\n");
    emuoc_statsOnType(*cI,total_size,total_freecap,total_lostba,
		stats_s,stats_bs,stats_o,stats_lba,stats_isc);
    num_size+=stats_s;
    bytes_size+=stats_bs+stats_lba+overhead_elem;
  }
  return;
}

// for containers with capacity (at the moment only vector?)
template<class myType>
inline void emuoc_statsOnType_sizeAndCapacity(const myType & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{
  MMEBUG("emuoc_statsOnType(vector): " << typeid(std::vector<myType>).name() << endl);

  (void) c;
  num_size=0;
  overhead=sizeof(myType);
  bytes_size=overhead;
  lostbyalign=0;
  iscontainer=true;

  MMEBUG("overhead: " << overhead
	<< "\nbytes_size: " << bytes_size
	<< "\ntotal_freecap: " << total_freecap
	<< endl);

  typename myType::const_iterator cI=c.begin();

  bool subelem_is_container=true;
  for(; cI != c.end(); cI++){
    size_t stats_s, stats_bs, stats_o, stats_lba;
    MMEBUG("loop ");
    emuoc_statsOnType(
      *cI,total_size,total_freecap,total_lostba,
      stats_s,stats_bs,stats_o,stats_lba,
      subelem_is_container);
    // no lost by align and no overhead of elements in vector
    bytes_size+=stats_bs;
    MMEBUG(stats_bs << " " << bytes_size << endl);
  }
  num_size+=c.size();
  MMEBUG("size: " << num_size << endl);
  MMEBUG("capacity: " << c.capacity() << endl);
  if(c.capacity() != num_size){
    size_t vtype=sizeof(typename myType::value_type);
    MMEBUG("empty element: " << vtype << endl);
    total_freecap+=vtype*(c.capacity()-num_size);
    bytes_size+=vtype*(c.capacity()-num_size);
  }
  if(bytes_size%sizeof(void *)){
    lostbyalign=sizeof(void *)-bytes_size%sizeof(void *);
    bytes_size+=lostbyalign;
  }

  total_lostba+=lostbyalign;

  MMEBUG("at end:"
	<< "\n\tbytes_size: " << bytes_size
	<< "\n\ttotal_freecap: " << total_freecap
	<< "\n\tlost by align: " << lostbyalign
	<< "\n\ttotal lba: " << total_lostba);
  MMEBUG("\nv_out\n");

  return;
}


#endif
