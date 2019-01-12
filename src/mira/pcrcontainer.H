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


#ifndef _mira_pcrcont_h_
#define _mira_pcrcont_h_

#include "stdinc/types.H"
#include "mira/readpool.H"


class PlacedContigReads
{
private:
  ReadPool * PCR_originalrp;  // original readpool

  // additional offset per read
  struct addoff_t {
    int32  addoffset;        // additional offset
    int32  urdid;            // read id in the unsorted readdump, >=0 OK, <0 should never happen

    inline addoff_t(int32 ao, int32 uid) : addoffset(ao),urdid(uid) {};
    static inline bool lt(const addoff_t & a, const addoff_t & b) { return a.addoffset<b.addoffset;}
  };

  // for list, sorted by position
  struct rposbin_t {
    uint32           offsetmapindex;   // index to offsetbin
    std::vector<addoff_t> readao;

    inline rposbin_t(uint32 omi, size_t reservesize) : offsetmapindex(omi) {readao.reserve(reservesize);};
  };


  // unsorted vector, linked to ReadPool::ReadContainer (same index)
  struct ancillaryinfo_t {
    int32 orpid;              // read id in the original readpool
    int8  direction;

    std::list<rposbin_t>::iterator rpbI;

    inline ancillaryinfo_t(int32 rpid, int8 dir, std::list<rposbin_t>::iterator I): orpid(rpid),direction(dir),rpbI(I) {};
  };

  // sorted by baseoffset (vector)
  struct offsettile_t{
    int32  from;
    std::list<rposbin_t>::iterator rpbI;       // iterator back to rposbin_t in list

    inline offsettile_t(int32 bofr, std::list<rposbin_t>::iterator I) : from(bofr),rpbI(I) {};
    static inline bool lt_offsetfrom(const offsettile_t & a, const offsettile_t & b) { return a.from<b.from;}
  };

/*************************************************************************
 *
 * Main logic for the containers and placement
 *
 *  - readdump & ancillaryinfo are linked: always the same size (except
 *    for a short time in placeRead() and placeRead_helper())
 *  - readdump and baseoffset link back via iterators to readpositioning
 *  - readpositioning links back via readao.urdid (index) to readdump/ancillaryinfo
 *                    links back via offsetmapindex to baseoffset
 *  - after having removed a read, ancillary info may contain entries paired
 *    to empty / invalid reads.
 *
 * readpositioning: the first element of each new bin has *always* an
 *  addoffset==0!
 *
 *************************************************************************/
  ReadPool::ReadContainer   PCR_readdump;
  std::vector<ancillaryinfo_t>   PCR_ancillaryinfo;

  std::list<rposbin_t>           PCR_readposbins;
  std::vector<offsettile_t>      PCR_offsetmap;

  uint32                    PCR_bo_binsize;

/*************************************************************************
 *
 * Main logic for the containers and placement
 *
 *************************************************************************/
  std::vector<std::list<rposbin_t>::iterator>    PCR_maprpids_to_rpb_v;      // map RP id to readposbin_t
  /* same as above, but with a map behaviour: for less than (8k?)  reads, use
     map, else will switch on the fly to the vector. Rationale: a vector is bad
     for projects with millions of reads contigs with low number of reads,
     a map is bad for large contigs with many reads */
  std::map<int32,std::list<rposbin_t>::iterator> PCR_maprpids_to_rpb_m;

  size_t PCR_numreads; /* number of placed reads. Cannot use PCR_ancillaryinfo because
			  removeRead() and the ReadContainer PCR_readdump will
			  not physically remove, but rather flag for re-use */

  suseconds_t PCR_time_sr_lb1;
  suseconds_t PCR_time_sr_lb2;
  suseconds_t PCR_time_sr_aoadj;
  suseconds_t PCR_time_sr_omadj;

  suseconds_t PCR_time_sb_total;
  suseconds_t PCR_time_sb_c2h;

  suseconds_t PCR_time_prh_pf;
  suseconds_t PCR_time_prh_a2b1;
  suseconds_t PCR_time_prh_a2b2;
  suseconds_t PCR_time_prh_a2b3;
  suseconds_t PCR_time_prh_a2b;


  /*************************************************************************************
   *
   *
   *  Iterators
   *
   *
   *************************************************************************************/
public:
  class const_iterator
    : public boost::iterator_facade<
    const_iterator                          // Derived
    , Read const             // Value
    , std::bidirectional_iterator_tag   //random_access_iterator_tag  //boost::forward_traversal_tag    // CategoryOrTraversal
    > {

  friend class PlacedContigReads;

  public:
    typedef uint32 raoindex_t;

  private:
    PlacedContigReads const * pcr;
    std::list<rposbin_t>::const_iterator rpbI;
    raoindex_t raoindex;

  public:
    inline const_iterator(const PlacedContigReads * apcr, std::list<rposbin_t>::const_iterator arpbI, uint32 araoindex) : pcr(apcr),rpbI(arpbI),raoindex(araoindex) { };

    inline uint32 getReadStartOffset() const {
      return pcr->PCR_offsetmap[rpbI->offsetmapindex].from + rpbI->readao[raoindex].addoffset;
    }
    inline int32 getOriginalReadPoolID() const {
      return pcr->PCR_ancillaryinfo[rpbI->readao[raoindex].urdid].orpid;
    }
    inline int32 getORPID() const {
      return getOriginalReadPoolID();
    }
    inline void setORPID(int32 neworpid) {
      // some const gymnastics to be able to assign the value in an otherwise const object
      *(const_cast<int32 *>(&(pcr->PCR_ancillaryinfo[rpbI->readao[raoindex].urdid].orpid)))=neworpid;
    }
    inline int8 getReadDirection() const {
      return pcr->PCR_ancillaryinfo[rpbI->readao[raoindex].urdid].direction;
    }
    inline int32 getURDID() const {
      return rpbI->readao[raoindex].urdid;
    }

    friend std::ostream & operator<<(std::ostream &ostr, const const_iterator & ci) {
      ostr << "pcr: " << ci.pcr << "\trpbI: " << &(*(ci.rpbI)) << "\traoi " << ci.raoindex << "\t";
      if(ci==ci.pcr->end()){
	ostr << "This is end()" << std::endl;
      }else{
	ci.pcr->dumpPCRIElement(ostr,ci);
      }
      return ostr;
    }


   /*******************************************************************
    *
    *  unclippedReadPosToContigPos
    *  We have a position in an UNCLIPPED read and we want to know the corresponding
    *  position in the contig
    *
    ******************************************************************/
    inline int32 unclippedReadPos2ContigPos(const int32 readpos) const
      {
	if (getReadDirection() > 0) {
	  return readpos - (*this)->getLeftClipoff() + getReadStartOffset() ;
	} else {
	  return (*this)->getRightClipoff()- readpos + getReadStartOffset() -1 ;
	}
      }

    /******************************************************************
     *
     * contigPosToUnclippedReadPos
     * We have a contig position and we want to know which position in the
     * UNCLIPPED FORWARD read corresponds to the contig position
     *
     ******************************************************************/

    inline int32 contigPos2UnclippedReadPos(const int32 contigpos) const
      {
	if (getReadDirection() > 0) {
	  return contigpos - getReadStartOffset() + (*this)->getLeftClipoff();
	} else {
	  // FALSCH return ric.offset + ric.read.getRightClipoff() -1 - contigpos;
	  return (*this)->getLenSeq() - (*this)->getRightClipoff() + contigpos - getReadStartOffset();
	}
      };
    inline int32 contigPos2RawReadPos(const int32 contigpos) const
      {
	if (getReadDirection() > 0) {
	  return contigPos2UnclippedReadPos(contigpos);
	} else {
	  return (*this)->calcComplPos(contigPos2UnclippedReadPos(contigpos));
	}
      };


    /******************************************************************
     *
     *
     *
     ******************************************************************/
    inline char getBase(const int32 contigpos) const
      {
	if(getReadDirection()>0){
	  return (*this)->getBaseInSequence(contigPos2UnclippedReadPos(contigpos));
	}else{
	  return (*this)->getBaseInComplementSequence(contigPos2UnclippedReadPos(contigpos));
	}
      }
    inline base_quality_t getQuality(const int32 contigpos) const
      {
	if(getReadDirection()>0){
	  return (*this)->getQualityInSequence(contigPos2UnclippedReadPos(contigpos));
	}else{
	  return (*this)->getQualityInComplementSequence(contigPos2UnclippedReadPos(contigpos));
	}
      }

    // this is a bad idea ... incredibly slow compared to currently used methods outside pcrI
    //   (calc read pos, then access base and qual directly on read)
    //
    //inline void getBaseAndQuality(const int32 contigpos, char & base, base_quality_t & qual) const
    //  {
    //	auto urp=contigPos2UnclippedReadPos(contigpos);
    //	if(getReadDirection()>0){
    //	  base=(*this)->getBaseInSequence(urp);
    //	  qual=(*this)->getQualityInSequence(urp);
    //	}else{
    //	  base=(*this)->getBaseInComplementSequence(urp);
    //	  qual=(*this)->getQualityInComplementSequence(urp);
    //
    //  }


  private:
    // Everything Boost's iterator facade needs
    friend class boost::iterator_core_access;

    // behaviour of dereference might seem a bit strange: it returns const Read &
    // i.e., when dereferenced one cannot access position or read direction of read
    //
    // Reason: the PlacedContigReads::const_iterator is behaving a bit like a SQL VIEW,
    //  i.e., it aggregates information which is not present together as "struct"
    //  I could copy all these things into a temporary structure, but why should this be
    //  done? For many use cases (iterating through a contig etc.) this would just waste
    //  time
    //
    // In the end this means: if one wants to get info on the placement, use the iterator
    //  functions above (getReadPosition(), getReadDirection(), getOriginalReadPoolID() etc.pp)
    // Also the other functions allowing to translate contig positions to read positions,
    //  get a base of a read in a contig etc., use the iterator functions.
    inline Read const & dereference() const
    {
      // TODO: eventually remove the BUGIFTHROW if prooves to be to slow (check that
      //  on projects)
      FUNCSTART("inline Read const & PlacedContigReads::const_iterator::dereference() const");
      BUGIFTHROW(*this==pcr->end(),"Trying to dereference an iterator pointing to end()???");
      return pcr->PCR_readdump[rpbI->readao[raoindex].urdid];
    }

    inline bool equal(const_iterator const & other) const
    {
      return rpbI == other.rpbI && raoindex == other.raoindex;
    }

    inline void increment()
    {
      ++raoindex;
      if(raoindex==rpbI->readao.size()){
	++rpbI;
	raoindex=0;
      }
    }

    inline void decrement()
    {
      if(raoindex>0){
	--raoindex;
      }else{
	--rpbI;
	raoindex=static_cast<raoindex_t>(rpbI->readao.size()-1);
      }
    }
  };


private:
  void addORPID2Map(int32 rpid, std::list<rposbin_t>::iterator rpbI);
  void delORPIDFromMap(int32 rpid);
  inline void updateMapBinOfORPID(int32 rpid, std::list<rposbin_t>::iterator rpbI){
    // careful with rpid: may be -1 for CER reads
    if(rpid>=0){
      if(!PCR_maprpids_to_rpb_v.empty()){
	if(rpid>=PCR_maprpids_to_rpb_v.size()){
	  std::cout << "updateMapBinOfORPID(): should never arrive here. rpid " << rpid << " >= PCR_maprpids_to_rpb_v.size() " << PCR_maprpids_to_rpb_v.size() << std::endl;
	  exit(1000);
	}
	PCR_maprpids_to_rpb_v[rpid]=rpbI;
      }else{
	auto mI=PCR_maprpids_to_rpb_m.find(rpid);
	if(mI==PCR_maprpids_to_rpb_m.end()){
	  std::cout << "updateMapBinOfRPID(): should never arrive here, did not find rpid " << rpid << std::endl;
	  exit(1000);
	}
	mI->second=rpbI;
      }
    }
  };


  std::vector<offsettile_t>::iterator searchOffsetTileForPlacement(int32 position);
  //std::list<PlacedContigReads::rposbin_t>::iterator findReadpositionIterator(int32 position);

  void splitBin(uint32 binindex);
  std::list<rposbin_t>::iterator placeRead_helper(int32 rpid,
						  int32 position,
						  int8 direction,
						  int32 urdid,
						  const_iterator::raoindex_t & araoindex);

public:
  PlacedContigReads(ReadPool & rp) : PCR_originalrp(&rp), PCR_bo_binsize(2048), PCR_numreads(0),
				     PCR_time_sr_lb1(0),
				     PCR_time_sr_lb2(0),
				     PCR_time_sr_aoadj(0),
				     PCR_time_sr_omadj(0),
				     PCR_time_sb_c2h(0),
				     PCR_time_sb_total(0),
				     PCR_time_prh_pf(0),
				     PCR_time_prh_a2b1(0),
				     PCR_time_prh_a2b2(0),
				     PCR_time_prh_a2b3(0),
				     PCR_time_prh_a2b(0)
    {};
  ~PlacedContigReads();

  // Copy constructor & copy operator
  PlacedContigReads(const PlacedContigReads&) = delete;
  PlacedContigReads const & operator=(PlacedContigReads const & other);


  inline const_iterator begin() const {
    return const_iterator(this,PCR_readposbins.begin(),0);
  }
  inline const_iterator end() const {
    return const_iterator(this,PCR_readposbins.end(),0);
  }
//  inline reverse_iterator rbegin() {
//    return reverse_iterator(this,PCR_readposbins.rbegin(),PCR_readposbins.rbegin()->readao.size()-1);
//  }
//  inline reverse_iterator rend() {
//    return reverse_iterator(this,PCR_readposbins.rend(),0);
//  }

  inline size_t size() const {return PCR_numreads;};
  inline bool empty() const {return PCR_numreads==0;};

  std::ostream & dumpPCRIElement(std::ostream &ostr, const const_iterator & ci) const;

  const_iterator getIteratorOfReadpoolID(int32 rpid);
  inline const Read & getReadAtURDID(int32 urdid) const {
    return PCR_readdump[urdid];
  }
  inline int8 getReadDirectionAtURDID(int32 urdid) const {
    FUNCSTART("inline int8 getReadDirectionAtURDID(int32 urdid) const");
    BUGIFTHROW(urdid >= static_cast<int32>(PCR_ancillaryinfo.size()),"urdid >= PCR_ancillaryinfo.size() ???");
    return PCR_ancillaryinfo[urdid].direction;
  }
  inline int32 getReadORPIDAtURDID(int32 urdid) const {
    FUNCSTART("inline int8 getReadORPIDAtURDID(int32 urdid) const");
    BUGIFTHROW(urdid >= static_cast<int32>(PCR_ancillaryinfo.size()),"urdid >= PCR_ancillaryinfo.size() ???");
    return PCR_ancillaryinfo[urdid].orpid;
  }

  void clear() {
    PCR_readdump.clear();
    PCR_ancillaryinfo.clear();
    PCR_readposbins.clear();
    PCR_offsetmap.clear();
    PCR_maprpids_to_rpb_v.clear();
    PCR_maprpids_to_rpb_m.clear();
    PCR_numreads=0;
  }

  void setBinSize(uint32 bs) {PCR_bo_binsize=bs;};
  const_iterator placeRead(const Read & theread, int32 rpid, int32 position, int8 direction);
  const_iterator removeRead(const_iterator pcrI);

  void shiftReads(int32 position, int32 offsetdiff);
  void shiftReadsBounceZero(int32 position, int32 offsetdiff);

  const_iterator getPCRIForReadsStartingAtPos(int32 position) const;


  void debugDump(bool shortdbg);

};


#endif
