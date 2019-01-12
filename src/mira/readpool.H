/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 *   and Thomas Pfisterer
 * Copyright (C) 2000 and later by Bastien Chevreux
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


#ifndef _readpool_h_
#define _readpool_h_

#include <string>

#include <unordered_map>

#include <boost/lambda/bind.hpp>

#include "stdinc/defines.H"

#include "mira/parameters.H"
#include "mira/read.H"


#undef yyFlexLexer
#define yyFlexLexer FAFlexLexer
#include <FlexLexer.h>


class ReadPool
{
public:
  class ReadContainer {
  private:
    std::deque<Read> RC_thepool;
    std::vector<Read *> RC_poolrptr;
    std::vector<uint32> RC_releasedidx; // index of free elements in the std::deque
  private:
    // sort criterion for standard MIRA readpool order
    // rails first (need that for skim!)
    // then reads with template partners by template id (sorted by segment number)
    // then rest
    inline bool sortcompMIRAstd(const uint32 i1, const uint32 i2){
      //std::cout << "scm " << i1 << ' ' << i2 << std::endl;
      if(RC_poolrptr[i1]->isRail()){
	if(RC_poolrptr[i2]->isRail()) return i1<i2;
	return true;
      }else if(RC_poolrptr[i2]->isRail()){
	return false;
      }
      if(RC_poolrptr[i1]->getTemplatePartnerID() >= 0
	 || RC_poolrptr[i2]->getTemplatePartnerID() >=0){
	if(RC_poolrptr[i1]->getTemplatePartnerID() >= 0
	   && RC_poolrptr[i2]->getTemplatePartnerID() >=0){
	  if(RC_poolrptr[i1]->getTemplateID() == RC_poolrptr[i2]->getTemplateID()){
	    return RC_poolrptr[i1]->getTemplateSegment() < RC_poolrptr[i2]->getTemplateSegment();
	  }else{
	    return RC_poolrptr[i1]->getTemplateID() < RC_poolrptr[i2]->getTemplateID();
	  }
	}else if(RC_poolrptr[i1]->getTemplatePartnerID() >=0) {
	  return true;
	}
      }
      return false;
    }
  public:
    ReadContainer() = default;
    // Copy operator
    ReadContainer(const ReadContainer&) = delete;
    ReadContainer const & operator=(ReadContainer const & other);

    inline size_t size() const { return RC_poolrptr.size();}
    inline size_t getNumActiveReads() const { return RC_poolrptr.size() - RC_releasedidx.size();}
    inline void clear() { RC_thepool.clear(); RC_poolrptr.clear(); RC_releasedidx.clear();}
    inline size_t provideEmptyRead() {
      size_t readidx=-1;
      if(RC_releasedidx.size()){
	// let's try without
	// RC_freerptr.back()->~Read();            // destructor, followed by ...
	// retread=new(RC_freerptr.back()) Read(); // ... a "placement new"
	readidx=RC_releasedidx.back();

	RC_releasedidx.pop_back();
      }else{
	readidx=size();
	RC_thepool.resize(RC_thepool.size()+1);
	RC_poolrptr.push_back(&(RC_thepool[readidx]));
      }
      return readidx;
    }

    // Read not physically deleted (also no destructor called),
    //  an empty, invalid read then exists at this index which will be re-used
    //  as soon as possible
    void releaseRead(size_t index) {
      FUNCSTART("Read & ReadContainer::realeaseRead(size_t index)");
      BUGIFTHROW(index >= RC_poolrptr.size(),"index (" << index << ") >= RC_poolrptr.size() (" << RC_poolrptr.size() << ") ?");
      RC_poolrptr[index]->discard();
      RC_releasedidx.push_back(static_cast<uint32>(index));
    }
    // returns read at the index in pool
    inline Read & getRead(size_t index) const {
      FUNCSTART("Read & ReadContainer::getRead(uint32 index)");
      BUGIFTHROW(index >= RC_poolrptr.size(),"error: index (" << index << ") >= RC_poolrptr.size() (" << RC_poolrptr.size() << ") ?");
      return *(RC_poolrptr[index]);
    }
    inline Read & operator [] (size_t index) const {return getRead(index);};
    inline Read & back() const {return *(RC_poolrptr.back());};

    void sortContainerToMIRAStandard(std::vector<uint32> & sortorder) {
      //std::cout << "sTTMS 1 " << RC_poolrptr.size() << std::endl;
      sortorder.resize(RC_poolrptr.size());
      //std::cout << "sTTMS 2" << std::endl;
      for(uint32 i=0; i<sortorder.size(); ++i) sortorder[i]=i;
      //std::cout << "sTTMS 3" << std::endl;
      sort(sortorder.begin(), sortorder.end(), boost::lambda::bind(&ReadContainer::sortcompMIRAstd, this, boost::lambda::_1, boost::lambda::_2));
      //std::cout << "sTTMS 4" << std::endl;
      std::vector<Read *> newrcpoolptr;
      newrcpoolptr.reserve(RC_poolrptr.size());
      //std::cout << "sTTMS 5" << std::endl;
      for(uint32 i=0; i<sortorder.size(); ++i) newrcpoolptr.push_back(RC_poolrptr[sortorder[i]]);
      //std::cout << "sTTMS 6" << std::endl;
      RC_poolrptr.swap(newrcpoolptr);
      //std::cout << "sTTMS 7" << std::endl;
    }

    void sortContainerToGivenOrder(const std::vector<uint32> & sortorder) {
      FUNCSTART("void ReadContainer::sortContainerToMIRAStandard(const std::vector<uint32> & sortorder)");
      BUGIFTHROW(4294967295 == size(), "uh oh ... 4294967295 != size()");

      std::vector<Read *> tmp;
      tmp.reserve(RC_poolrptr.size());
      for(uint32 i=0; i<sortorder.size(); ++i) {
	if(sortorder[i]<size()){
	  if(RC_poolrptr[sortorder[i]]==nullptr){
	    MIRANOTIFY(Notify::FATAL,"RC_poolrptr[sortorder[i]]==nullptr ??? Did you specify a read twice in your input for sorting?");
	  }
	  tmp.push_back(RC_poolrptr[sortorder[i]]);
	  RC_poolrptr[sortorder[i]]=nullptr;
	}
      }
      for(uint32 i=0; i<RC_poolrptr.size(); ++i) {
	if(RC_poolrptr[i]!=nullptr) tmp.push_back(RC_poolrptr[i]);
      }
      RC_poolrptr.swap(tmp);

      FUNCEND();
    }
    void dumpDebug(){
      std::cout << "RC_thepool: " << RC_thepool.size() << std::endl;
      std::cout << "RC_poolrptr: " << RC_poolrptr.size() << std::endl;
      std::cout << "RC_releasedidx: " << RC_releasedidx.size() << std::endl;
    }
  };

  // Variables
private:
  static std::string RP_missingfastaqual_resolvemsg;
  const static base_quality_t RP_sxa2phredmap[256];

  ReadContainer REP_thepool3;

  bool REP_allownameindex;
  std::unordered_map<std::string, uint32> REP_nameindex;

  // Functions
private:
  void foolCompiler();

  static void mergeXMLTraceInfo_callback(void * callerobjptr, std::list<std::string> & elements, std::list<std::string> & elements_cdata);

  inline void rpDateStamp() {
    dateStamp(std::cout);
  }


public:
  //ReadPool(std::vector<MIRAParameters> * params);
  ReadPool();
  ~ReadPool();

  // currently do not allow copying of readpool
  ReadPool & operator=(const ReadPool&) = delete;

  void discard();
  size_t estimateMemoryUsage() const;

  inline size_t provideEmptyRead() { REP_nameindex.clear(); return REP_thepool3.provideEmptyRead();}
  inline void releaseRead(size_t index) {REP_thepool3.releaseRead(index);}
  inline void dumpRCDebug() {REP_thepool3.dumpDebug(); }

  inline size_t size() const {return REP_thepool3.size();};
  inline size_t getNumActiveReads() const {return REP_thepool3.getNumActiveReads();};

  // returns read at the index in pool
  inline Read & getRead(uint32 index) const {
    return REP_thepool3.getRead(index);
  }
  inline Read & operator [] (uint32 index) const {return REP_thepool3.getRead(index);};

  void allowNameIndex(bool b) {REP_allownameindex=b; if(!b) REP_nameindex.clear();}
  void refreshNameIndex();

  // returns (first) read with name i.p.
  Read & getRead(const std::string & name);

  // returns id of (first) read with name
  int32  getReadIndex(const std::string & name);

  void adaptFASTQQualValues(uint32 rpstart,
			    uint32 rpend,
			    base_quality_t fastqoffset,
			    bool verbose,
			    uint32 longestread=0,
			    base_quality_t seenminqual=255,
			    base_quality_t seenmaxqual=0);

  void deleteReadsByName(const std::string & nfile,
			 bool invertselection);

  void loadQualitiesFromSCF(uint32 force,
			    bool loadfailnoerror,
			    const std::string & log1name,
			    const std::string & log2name);
  void loadQualitiesFromSCF(uint32 force, bool loadfailnoerror);

  void mergeXMLTraceInfo(const std::string & xmlfile, uint8 nagwarn_templateproblems);
  void mergeSSAHA2SMALTVecScreenData(const std::string & ssahafile,
				     bool issmalt,
				     std::vector<MIRAParameters> & miraparams,
				     const std::string & logname,
				     const std::string & logprefix);

  void dumpAs(std::ostream & ostr, uint8 astype, bool alsoinvalids) const;
  void saveAsMAF(std::ostream & ostr, bool alsoinvalids) const;
  void dumpAsEXPs(std::string & dirname) const;
  void dumpPoolInfo(std::ostream & ostr) const;


  bool makeTemplateIDs(uint8 nagwarn_templateproblems, bool verbose=true);
  void checkTemplateIDs(const std::string & errmsg);

  void sortPoolToMIRAStandard(std::vector<uint32> & sortorder) {REP_thepool3.sortContainerToMIRAStandard(sortorder);}
  void sortPoolToGivenOrder(const std::vector<uint32> & sortorder) {REP_thepool3.sortContainerToGivenOrder(sortorder);}

  bool checkForDuplicateReadNames(uint8 nagwarn_templateproblems);

  void adjustIllegalQualities(base_quality_t bq);
};

#endif
