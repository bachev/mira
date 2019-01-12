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


#ifndef _mira_StringContainer_h_
#define _mira_StringContainer_h_

#include <iostream>

#include "stdinc/stlincludes.H"

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include <boost/bind.hpp>

//#define SCCEBUG(bla) {std::cout << bla; std::cout.flush();}
#define SCCEBUG(bla)

template <class T>
class StringContainer
{
public:

  class Entry {
  private:
    T stringid;
  public:
    Entry() : stringid(0) {};
    inline bool operator==(const Entry &other) const {return (stringid==other.stringid);};
    inline bool operator!=(const Entry &other) const {return (stringid!=other.stringid);};

    inline T getSCID() const {return stringid;};
    inline void setSCID(T id) {stringid=id;};
  };

private:
  std::string SC_name;                /* name of the container, used for
				    error messages and debugging */
  std::vector<std::string> SC_thestrings;  // all the strings of that container
  std::vector<T> SC_stringsorder;     // order of lexically sorted strings in "thestrings"
  T SC_maxnumentries;            /* max num of entries this container can have
				    (governed by size of template type T) */

  bool SC_issorted;

// functions

public:
  void dump(std::ostream & ostr);
  void status(std::ostream & ostr);

private:
  void init(const char * name);

  typename std::vector<T>::iterator searchInsertIndex(const std::string & s);

  inline bool sortcontainer_funccmp(const T & i, const T & j) const {
    return SC_thestrings[i]<SC_thestrings[j];
  }
  inline void sortContainer() {
    if(SC_issorted || SC_stringsorder.size()<=1) return;

    SCCEBUG("Must sort\n");
    if(SC_thestrings.size()==0){
      std::cout << "Oooooops? SC_thestrings.size()==0 and still want to sort???" << std::endl;
    }

    //dump(std::cout);
    sort(SC_stringsorder.begin(),SC_stringsorder.end(),
	 bind(&StringContainer<T>::sortcontainer_funccmp, this, _1, _2));
    SC_issorted=true;
    //dump(std::cout);
  }




public:
  StringContainer(const char * name)  {init(name);}
  StringContainer(const std::string & name) {init(name.c_str());}

  inline size_t size() const {return SC_thestrings.size();};
  inline void reserve(size_t s) {SC_thestrings.reserve(s);SC_stringsorder.reserve(s); return;};
  inline void trash() {SC_thestrings.clear();SC_stringsorder.clear(); init(SC_name.c_str());return;};

  Entry addEntry(const std::string & s);
  Entry addEntryNoDoubleCheck(const std::string & s);
  Entry hasEntry(const std::string & s);

  inline const std::string & getEntry(Entry e) const {
    FUNCSTART("const std::string & StringContainer::getEntry(Entry e)");
    BUGIFTHROW(e.getSCID() >= SC_thestrings.size(), "StringContainer: " << SC_name << "\nTried to get non-existent entry: " << static_cast<uint64>(e.getSCID()) << std::endl);
    FUNCEND();
    return SC_thestrings[e.getSCID()];
  }

  inline Entry emptyEntry() const {
    Entry e;
    e.setSCID(0);
    return e;
  }

};


template <class T>
void StringContainer<T>::init(const char * name)
{
  FUNCSTART("template <class T> void StringContainer<T>::init(const char * name)");

  SCCEBUG("SC init: " << name << std::endl);
  SC_name=name;
  SC_thestrings.resize(1);
  SC_stringsorder.push_back(0);
  SC_maxnumentries=0;
  SC_maxnumentries--;
  SC_issorted=true;

  //status(std::cout);
  SCCEBUG("SC end init: " << name << std::endl);

  FUNCEND();
}

template <class T>
typename std::vector<T>::iterator StringContainer<T>::searchInsertIndex(const std::string & s)
{
  FUNCSTART("template <class T> typename std::vector<T>::iterator StringContainer<T>::searchInsertIndex(const string & s)");

//  std::cout << "StringContainer<T>::searchInsertIndex '" << s << "'" << std::endl;

  if(!SC_issorted) sortContainer();

  typename std::vector<T>::iterator first, last, it;

  first=SC_stringsorder.begin();
  last=SC_stringsorder.end();

  ptrdiff_t count, step;
  count = distance(first,last);
  while (count>0){
    it = first; step=count/2; advance (it,step);
    if (SC_thestrings[*it]<s) {
      first=++it; count-=step+1;
    }else{
      count=step;
    }
  }

  FUNCEND();
  return first;
}

template <class T>
typename StringContainer<T>::Entry StringContainer<T>::hasEntry(const std::string & s)
{
  FUNCSTART("template <class T> StringContainer<T>::Entry<T> StringContainer<T>::hasEntry(const string & s)");

  // empty entries get quick treatment
  if(s.empty()){
    FUNCEND();
    return emptyEntry();
  }

  Entry e;

  typename std::vector<T>::iterator inserti=searchInsertIndex(s);
  if(inserti!=SC_stringsorder.end()
     && SC_thestrings[*inserti]==s) {
    e.setSCID(*inserti);
    FUNCEND();
    return e;
  }

  FUNCEND();
  return emptyEntry();
}

template <class T>
typename StringContainer<T>::Entry StringContainer<T>::addEntry(const std::string & s)
{
  FUNCSTART("template <class T> StringContainer<T>::Entry<T> StringContainer<T>::addEntry(const std::string & s)");

  // empty entries get quick treatment
  if(s.empty()){
    FUNCEND();
    return emptyEntry();
  }

  SCCEBUG("SC " << SC_name << ": adding non-empty entry '" << s << "'" << std::endl);

  if(SC_thestrings.size()==SC_maxnumentries){
    std::cout << "Oooops? Going to throw in addEntry because of this: " << SC_thestrings.size() << " " << SC_maxnumentries << " for " << s << std::endl;
    dump(std::cout);
    MIRANOTIFY(Notify::INTERNAL, "Tried to add '" << s << "', but max number of entries (" << SC_thestrings.size() << ") reached.\n");
  }

  Entry e;

  typename std::vector<T>::iterator inserti=searchInsertIndex(s);
  if(inserti!=SC_stringsorder.end()
     && SC_thestrings[*inserti]==s) {
    e.setSCID(*inserti);
    FUNCEND();
    return e;
  }

  //if(inserti==SC_stringsorder.end()){
  //  std::cout << "insert at end\n";
  //}else{
  //  std::cout << "insert at pos" << inserti-SC_stringsorder.begin() << std::endl;;
  //}


  e.setSCID(static_cast<T>(SC_thestrings.size()));
  SC_stringsorder.insert(inserti,static_cast<T>(SC_thestrings.size()));
  SC_thestrings.push_back(s);

  FUNCEND();
  return e;
}



template <class T>
typename StringContainer<T>::Entry StringContainer<T>::addEntryNoDoubleCheck(const std::string & s)
{
  FUNCSTART("template <class T> StringContainer<T>::Entry<T> StringContainer<T>::addEntryNoDoubleCheck(const std::string & s)");

  // empty entries get quick treatment
  if(s.empty()){
    FUNCEND();
    return emptyEntry();
  }

  SCCEBUG("SC " << SC_name << ": adding non-empty entry '" << s << "'" << std::endl);

  if(SC_thestrings.size()==SC_maxnumentries){
    std::cout << "Oooops? Going to throw in addEntryNoDoubleCheck because of this:\n";
    dump(std::cout);
    MIRANOTIFY(Notify::INTERNAL, "Tried to add '" << s << "', but max number of entries (" << SC_thestrings.size() << ") reached.\n");
  }

  Entry e;
  e.setSCID(static_cast<T>(SC_thestrings.size()));
  SC_stringsorder.push_back(static_cast<T>(SC_thestrings.size()));
  SC_thestrings.push_back(s);

  SC_issorted=false;

  FUNCEND();
  return e;
}


template <class T>
void StringContainer<T>::status(std::ostream & ostr)
{
  FUNCSTART("template <class T> void StringContainer<T>::status(std::ostream & ostr)");

  ostr << "SC " << SC_name
       << " issorted (" << SC_issorted << ") capacity "
       << static_cast<uint64>(SC_maxnumentries) << "(" << sizeof(SC_maxnumentries) << ") size "
       << SC_thestrings.size() << std::endl;

  FUNCEND();
}

template <class T>
void StringContainer<T>::dump(std::ostream & ostr)
{
  FUNCSTART("template <class T> void StringContainer<T>::dump(std::ostream & ostr)");

  status(ostr);
  for(size_t i=0; i<SC_thestrings.size(); i++){
    ostr << i << "\traw:" << SC_thestrings[i]
	 << "\to:" << static_cast<size_t>(SC_stringsorder[i])
	 << "\tsrt:" << SC_thestrings[SC_stringsorder[i]]
	 << '\n';
  }
  ostr.flush();
  FUNCEND();
}



#endif
