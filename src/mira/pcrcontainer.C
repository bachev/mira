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



#include "mira/pcrcontainer.H"

#include "util/stlimprove.H"


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
PlacedContigReads::~PlacedContigReads()
{
//    cout << "~PlacedContigReads()"
//	 << "\nPCR_readdump: " << PCR_readdump.size()
//	 << "\nPCR_offsetmap: " << PCR_offsetmap.size()
//	 << "\nPCR_readposbins: " << PCR_readposbins.size()
//	 << "\nPCR_ancillaryinfo: " << PCR_ancillaryinfo.size()
//	 << "\nsr_lb1: " << PCR_time_sr_lb1
//	 << "\nsr_lb2: " << PCR_time_sr_lb2
//	 << "\nsr_aoadj: " << PCR_time_sr_aoadj
//	 << "\nsr_omadj: " << PCR_time_sr_omadj
//	 << "\nsb_c2h: " << PCR_time_sb_c2h
//	 << "\nsb_total: " << PCR_time_sb_total
//	 << "\nprh_pf: " << PCR_time_prh_pf
//	 << "\nprh_a2b1: " << PCR_time_prh_a2b1
//	 << "\nprh_a2b2: " << PCR_time_prh_a2b2
//	 << "\nprh_a2b3: " << PCR_time_prh_a2b3
//	 << "\nprh_a2b: " << PCR_time_prh_a2b
//	 << endl;
};



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

PlacedContigReads const & PlacedContigReads::operator=(PlacedContigReads const & other)
{
  if(this != &other){
    PlacedContigReads(*(other.PCR_originalrp)); // sets readpool, zeroes all timing counters
    PCR_readdump.clear();
    PCR_ancillaryinfo.clear();
    PCR_readposbins.clear();
    PCR_offsetmap.clear();
    PCR_bo_binsize=other.PCR_bo_binsize;
    for(auto opcrI=other.begin(); opcrI != other.end(); ++opcrI){
      placeRead(*opcrI,opcrI.getORPID(),opcrI.getReadStartOffset(),opcrI.getReadDirection());
    }
  }

  return *this;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void PlacedContigReads::debugDump(bool shortdbg)
{
  FUNCSTART("void PlacedContigReads::debugDump(bool shortdbg)");
  cout << "debugDump PlacedContigReads\nbinsize: " << PCR_bo_binsize << endl;
  cout << "rd size: " << PCR_readdump.size() << endl;
  cout << "rd active: " << PCR_readdump.getNumActiveReads() << endl;
  cout << "anc size: " << PCR_ancillaryinfo.size() << endl;
  cout << "size(): " << size() << endl;
  cout << "readposbins size: " << PCR_readposbins.size() << endl;
  if(shortdbg) return;
  for(uint32 i=0; i<PCR_offsetmap.size(); ++i){
    cout << "om " << i
	 << "\tf: " << PCR_offsetmap[i].from
	 << "\tnum: " << PCR_offsetmap[i].rpbI->readao.size()
	 << "\trpbI: " << &(*(PCR_offsetmap[i].rpbI))
	 << endl;
  }
  cout << endl;

  uint32 elemnum=0;
  for(auto & rpbe : PCR_readposbins){
    cout << "rpbe " << elemnum << "\tomi: " << rpbe.offsetmapindex << " (from: " << PCR_offsetmap[rpbe.offsetmapindex].from << ")" << endl;
    for(uint32 i=0; i<rpbe.readao.size(); ++i){
      cout << "aoi " << i << "\tao: " << rpbe.readao[i].addoffset
	   << "\turdid: " << rpbe.readao[i].urdid << endl;
    }
    ++elemnum;
  }

  for(uint32 i=0; i<PCR_readdump.size(); ++i){
    cout << "rd rn " << i << "\tname: ";
    //cout.flush();
    cout << PCR_readdump[i].getName();
    //cout.flush();
    BUGIFTHROW(i>PCR_ancillaryinfo.size(),"i>PCR_ancillaryinfo.size() ???");
    if(i==PCR_ancillaryinfo.size()){
      cout << "ancillary info not available yet (OK while debugging PCR, not OK else!";
    }else{
      cout << "\torpid: " << PCR_ancillaryinfo[i].orpid << "\tdir: " << static_cast<int16>(PCR_ancillaryinfo[i].direction);
      if(PCR_ancillaryinfo[i].rpbI == PCR_readposbins.end()){
	cout << "\trpbI.end()";
      }else{
	cout << "\trpbI->omi: " <<PCR_ancillaryinfo[i].rpbI->offsetmapindex
	     << " (from: " << PCR_offsetmap[PCR_ancillaryinfo[i].rpbI->offsetmapindex].from << ")";
      }
    }
    cout << endl;
  }
  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void PlacedContigReads::addORPID2Map(int32 rpid, std::list<rposbin_t>::iterator rpbI)
{
  if(rpid>=0){
    if(!PCR_maprpids_to_rpb_v.empty()){
      if(rpid>=static_cast<int32>(PCR_maprpids_to_rpb_v.size())){
	PCR_maprpids_to_rpb_v.resize(static_cast<uint64>(rpid)*2,PCR_readposbins.end());
      }
      PCR_maprpids_to_rpb_v[rpid]=rpbI;
    }else{
      if(PCR_maprpids_to_rpb_m.size()<8192){
	PCR_maprpids_to_rpb_m.insert(std::pair<int32,std::list<rposbin_t>::iterator>(rpid,rpbI));
      }else{
	// switch from map to vector
	PCR_maprpids_to_rpb_v.resize(PCR_originalrp->size(),PCR_readposbins.end());
	for(auto & m : PCR_maprpids_to_rpb_m){
	  PCR_maprpids_to_rpb_v[m.first]=m.second;
	}
	PCR_maprpids_to_rpb_m.clear();
	PCR_maprpids_to_rpb_v[rpid]=rpbI;
      }
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void PlacedContigReads::delORPIDFromMap(int32 rpid)
{
  FUNCSTART("void PlacedContigReads::delRPIDFromMap(int32 rpid)");
  if(!PCR_maprpids_to_rpb_v.empty()){
    PCR_maprpids_to_rpb_v[rpid]=PCR_readposbins.end();
  }else{
    auto tmp=PCR_maprpids_to_rpb_m.erase(rpid);
    BUGIFTHROW(tmp!=1,"Erased " << tmp << " instances of rpid " << rpid << " from map???");
  }
  FUNCEND();
}


/*************************************************************************
 *
 * Returns end() if rpid not present
 *
 *
 *************************************************************************/

PlacedContigReads::const_iterator PlacedContigReads::getIteratorOfReadpoolID(int32 rpid)
{
  FUNCSTART("PlacedContigReads::iterator PlacedContigReads::getIteratorOfReadpoolID(int32 rpid)");

  const_iterator retI(end());

  if(rpid>=0){
    auto rpbI(PCR_readposbins.end());
    if(!PCR_maprpids_to_rpb_v.empty()){
      if(static_cast<size_t>(rpid)<PCR_maprpids_to_rpb_v.size()) rpbI=PCR_maprpids_to_rpb_v[rpid];
    }else{
      auto mI=PCR_maprpids_to_rpb_m.find(rpid);
      if(mI!=PCR_maprpids_to_rpb_m.end()) rpbI=mI->second;
    }

    if(rpbI!=PCR_readposbins.end()){
      auto raoI=rpbI->readao.begin();
      for(; raoI!=rpbI->readao.end(); ++raoI){
	if(PCR_ancillaryinfo[raoI->urdid].orpid==rpid) break;
      }
      BUGIFTHROW(raoI==rpbI->readao.end(), "Should never happen, did not find rpid " << rpid);
      retI=const_iterator(this,rpbI,
			  static_cast<const_iterator::raoindex_t>(raoI-rpbI->readao.begin()));
    }
  }

  FUNCEND();
  return retI;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla; cout.flush();}
#ifdef CEBUG
#define GIMMECEBUGBACK
#undef CEBUG
#endif
std::vector<PlacedContigReads::offsettile_t>::iterator PlacedContigReads::searchOffsetTileForPlacement(int32 position)
{
  auto fI=mstd::upper_bound(PCR_offsetmap,
			    offsettile_t(position,PCR_readposbins.end()),       // comparator object, only position matters
			    offsettile_t::lt_offsetfrom);

#ifndef CEBUG
  while(fI!=PCR_offsetmap.begin() && (fI-1)->from + (fI-1)->rpbI->readao.back().addoffset >= position) --fI;
#else
  while(fI!=PCR_offsetmap.begin()){
    CEBUG("PCR sotfp w1 ls " << fI-PCR_offsetmap.begin() << " s " << PCR_offsetmap.size() << endl);
    CEBUG("PCR sotfp w1 -1 from " << (fI-1)->from << endl);
    CEBUG("PCR sotfp w1 -1 omi " << (fI-1)->rpbI->offsetmapindex << endl);
    CEBUG("PCR sotfp w1 -1 aos " << (fI-1)->rpbI->readao.size() << endl);
    CEBUG("rpbI: " << &(*((fI-1)->rpbI)) << endl);

    if(!((fI-1)->from + (fI-1)->rpbI->readao.back().addoffset >= position)) break;
    --fI;
  }
  CEBUG("PCR sotfp 3" << endl);
#endif
  while(fI!=PCR_offsetmap.end() && fI->from + fI->rpbI->readao.back().addoffset < position) ++fI;
  if(fI == PCR_offsetmap.end()
    && fI!=PCR_offsetmap.begin()){
    --fI;
    if(fI->rpbI->readao.size() == fI->rpbI->readao.capacity()) ++fI;
  }
  return fI;
};
#ifdef GIMMECEBUGBACK
#undef GIMMECEBUGBACK
#define CEBUG(bla)
#endif



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla; cout.flush();}
void PlacedContigReads::splitBin(uint32 binindex)
{
  FUNCSTART("void PlacedContigReads::splitBin(uint32 binindex)");

  timeval tvp;
  timeval tvt;

  gettimeofday(&tvt,nullptr);

  CEBUG("want to split offsetbin " << binindex << endl);
  BUGIFTHROW(binindex>=PCR_offsetmap.size(), "binindex>=PCR_offsetmap.size() ???");

  // num elements of first partial bin after split
  auto numfirst=PCR_offsetmap[binindex].rpbI->readao.size()/2;

  // make a duplication of the offsettile
  {
    auto boiI=PCR_offsetmap.begin();
    std::advance(boiI,binindex);
    auto tmp=PCR_offsetmap[binindex];
    PCR_offsetmap.insert(boiI,tmp);

    // set "from" in offsetmap to correct values
    PCR_offsetmap[binindex+1].from=PCR_offsetmap[binindex].from+PCR_offsetmap[binindex].rpbI->readao[numfirst].addoffset;

  }

  CEBUG("after DUP \n"; debugDump());

  // insert new rposbin_t after current one
  auto rpbI=PCR_offsetmap[binindex].rpbI;
  {
    ++rpbI;
    rpbI=PCR_readposbins.insert(rpbI,rposbin_t(binindex+1,PCR_bo_binsize));
    PCR_offsetmap[binindex+1].rpbI=rpbI;

    // and copy 2nd half of readao of first to second readao
    //  already adjusting the additional offsets
    //  & rpbI of ancillaryinfo_t
    //  & the rpid2rpbI map
    int32 offsetdiff=PCR_offsetmap[binindex+1].from-PCR_offsetmap[binindex].from;
    auto raoI=PCR_offsetmap[binindex].rpbI->readao.begin();
    std::advance(raoI,numfirst);
    gettimeofday(&tvp,nullptr);
    // a simple unrolling of this is pointless
    // furthermore, total runtime for splitBin on 10m inserts is 0.25s
    // keep as is
    for(; raoI != PCR_offsetmap[binindex].rpbI->readao.end(); ++raoI){
      rpbI->readao.push_back(addoff_t(raoI->addoffset-offsetdiff,raoI->urdid));
      PCR_ancillaryinfo[raoI->urdid].rpbI=rpbI;
      updateMapBinOfORPID(PCR_ancillaryinfo[raoI->urdid].orpid,rpbI);
    }
    PCR_time_sb_c2h+=diffsuseconds(tvp);

    // adjust size of first readao (it's an adjustment down, but need to pass default object anyway
    PCR_offsetmap[binindex].rpbI->readao.resize(numfirst,addoff_t(0,0));
  }

  // adjust offsetmap indexes
  rpbI=PCR_offsetmap[binindex+1].rpbI;
  ++rpbI;
  for(; rpbI!=PCR_readposbins.end(); ++rpbI){
    rpbI->offsetmapindex+=1;
  }

  PCR_time_sb_total+=diffsuseconds(tvt);

  CEBUG("after split\n"; debugDump());

  FUNCEND();

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla; cout.flush();}
PlacedContigReads::const_iterator PlacedContigReads::placeRead(const Read & theread, int32 rpid, int32 position, int8 dir)
{
  FUNCSTART("void PlacedContigReads::insertRead(Read theread, int32 rpid, int32 position, int8 dir)");

//#ifdef PARANOIABUGTRACKFLAG
// really just for paranoia?
//  BUGIFTHROW(getURDIDOfReadpoolID(rpid)>=0,"rpid " << rpid << " already present ???");
//#endif

  BUGIFTHROW(dir!=-1 && dir!=1,"dir == " << static_cast<int16>(dir) << " ???");
  BUGIFTHROW(position<0,"position " << position << " < 0 ???");

  //CEBUG("PCR PR BEFORE\n"; debugDump());
  //if(PCR_readdump.size() > 8200) {
  //  cout << "pcr before" << endl;
  //  debugDump(false);
  //}
  int32 urdid=static_cast<int32>(PCR_readdump.provideEmptyRead());
  CEBUG("got urdid " << urdid << endl);
  //if(urdid==8216) debugDump(false);
  //cout << "dbg end" << endl;
  PCR_readdump[urdid]=theread;
  //CEBUG("read assigned" << endl);
  //if(urdid==8216) debugDump(false);

  // PCR_readposbins.end() below is a placeholder, is overwritten a couple of lines later
  //  but splitBin() in placeRead_helper()
  //  expects readdump and ancillaryinfo to exist and be consistent
  if(urdid>=static_cast<int32>(PCR_ancillaryinfo.size())){
    PCR_ancillaryinfo.push_back(ancillaryinfo_t(rpid,dir,PCR_readposbins.end()));
  }else{
    PCR_ancillaryinfo[urdid]=ancillaryinfo_t(rpid,dir,PCR_readposbins.end());
  }
  CEBUG("ancillary created" << endl);
  //if(urdid==8216) debugDump(false);

  const_iterator::raoindex_t araoindex;
  auto rpbI=placeRead_helper(rpid,position,dir,urdid,araoindex);
  CEBUG("helper done" << endl);

  PCR_ancillaryinfo[urdid].rpbI=rpbI;

  CEBUG("ancillary finished\n");

  addORPID2Map(rpid,rpbI);

  CEBUG("ORPID added 2 map\n");

  ++PCR_numreads;

  //CEBUG("PCR PR INSERTED\n"; debugDump(););

  CEBUG("PCR PR end" << endl);

  FUNCEND();
  return const_iterator(this,rpbI,araoindex);
}
//#define CEBUG(bla)

/*************************************************************************
 *
 * Note:
 *  Using the "const_iterator::raoindex_t & araoindex" to give back
 *   the raoindex of the newly placed read is a bit of a hack
 *  Normally I'd let _helper return a fully valid const_iterator,
 *   but the caller needs a non-const std::list<PlacedContigReads::rposbin_t>::iterator
 *   and there is absolutely no way to cast a std::list<>::const_iterator (which
 *   is stored in the PlacedContigReads::const_iterator) into a
 *   std::list<>::iterator. Hence this hack.
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla; cout.flush();}
std::list<PlacedContigReads::rposbin_t>::iterator PlacedContigReads::placeRead_helper(int32 rpid, int32 position, int8 dir, int32 urdid, const_iterator::raoindex_t & araoindex)
{
  FUNCSTART("void PlacedContigReads::insertRead(Read r, int32 rpid, int32 position, int8 dir, int32 urdid)");

  decltype(PCR_readposbins.end()) rpbI;
  araoindex=0;

  auto otI=searchOffsetTileForPlacement(position);
  CEBUG("got otI " << otI-PCR_offsetmap.begin() << " (size " << PCR_offsetmap.size() << ") " << endl);

  if(otI == PCR_offsetmap.end()){
    CEBUG("push back");
    // easy case: this just needs a new bin

    // the following line is a push_back, but it does not take more time and we get back the iterator
    //  to the last element gratis

    rpbI= PCR_readposbins.insert(PCR_readposbins.end(),rposbin_t(PCR_offsetmap.size(),PCR_bo_binsize));
    PCR_offsetmap.push_back(offsettile_t(position,rpbI));

    rpbI->readao.push_back(addoff_t(0,urdid));

    CEBUG("new rp bin\n");
  }else if((otI==PCR_offsetmap.begin() && position <= otI->from) && otI->rpbI->readao.size() == PCR_bo_binsize){
    CEBUG("push front");
    // easy case: this just needs a new bin

    timeval tvp;
    gettimeofday(&tvp,nullptr);

    // the following line is a push_front
    rpbI= PCR_readposbins.insert(PCR_readposbins.begin(),rposbin_t(0,PCR_bo_binsize));
    rpbI->readao.push_back(addoff_t(0,urdid));
    PCR_offsetmap.insert(PCR_offsetmap.begin(),offsettile_t(position,rpbI));
    ++rpbI;
    for(; rpbI!=PCR_readposbins.end(); ++rpbI) rpbI->offsetmapindex+=1;
    rpbI=PCR_readposbins.begin();
    PCR_time_prh_pf+=diffsuseconds(tvp);

    CEBUG("new rp bin\n");
  }else{
    CEBUG("must add to bin\n");

    CEBUG("have otI " << otI-PCR_offsetmap.begin() << endl);

    if(otI->rpbI->readao.size() < otI->rpbI->readao.capacity()){
      timeval tvt;
      gettimeofday(&tvt,nullptr);

      timeval tvp;
      gettimeofday(&tvp,nullptr);

      int32 positionadditionaloffset=position-(otI->from);

      rpbI=otI->rpbI;
      auto raoI=lower_bound(rpbI->readao.begin(),rpbI->readao.end(),
			    addoff_t(positionadditionaloffset,0),
			    addoff_t::lt);

      araoindex=static_cast<const_iterator::raoindex_t>(raoI-rpbI->readao.begin());

      PCR_time_prh_a2b1+=diffsuseconds(tvp);
      gettimeofday(&tvp,nullptr);

      // if we're adding to the end of the container, things are fast and easy ...
      if(raoI==rpbI->readao.end()){
	rpbI->readao.push_back(addoff_t(positionadditionaloffset,urdid));
      }else{
	// ... else we need to copy a bit around

	// The push_back/memmove is a 33% faster replacement for
	//    raoI=rpbI->readao.insert(raoI,addoff_t(positionadditionaloffset,urdid));

	rpbI->readao.push_back(addoff_t(0,0)); // just dummy, will be overwritten either by memmove or afterwards

	// we should recalc raoI. Even if not really needed in this case
	//  ... but it's fast enough to not be of any concern.
	raoI=rpbI->readao.begin();
	std::advance(raoI,araoindex);

	CEBUG("Before memmove. Have araoindex " << araoindex << " and raoI at " << raoI-rpbI->readao.begin() << " with readao size " << rpbI->readao.size() << "\n");
	CEBUG("memadr  : " << &(*(raoI)) << endl);
	CEBUG("memadr+1: " << &(*(raoI+1)) << endl);
	memmove(&(*(raoI+1)),
		&(*raoI),
		sizeof(addoff_t)*(rpbI->readao.end()-raoI-1));
	CEBUG("After memmove\n");
	raoI->addoffset=positionadditionaloffset;
	raoI->urdid=urdid;
	CEBUG("After assign\n");
      }

      PCR_time_prh_a2b2+=diffsuseconds(tvp);
      gettimeofday(&tvp,nullptr);

      // have we been inserted at the very front of the bin, with an offset lower than the previous one?
      if(positionadditionaloffset<0){
	CEBUG("ADJUST B1" << endl);// debugDump(false));
	// if yes: adjust additional offset of reads for this bin

	// unroll: 2x faster
	//__builtin_prefetch(&(*(omI+32)), 1, 3);
	for(uint32 count=static_cast<uint32>((rpbI->readao.end()-raoI)/4); count; --count){
	  raoI->addoffset-=positionadditionaloffset;
	  ++raoI;
	  raoI->addoffset-=positionadditionaloffset;
	  ++raoI;
	  raoI->addoffset-=positionadditionaloffset;
	  ++raoI;
	  raoI->addoffset-=positionadditionaloffset;
	  ++raoI;
	}
	CEBUG("ADJUST B2" << endl);// debugDump(false));
	for(;raoI!=rpbI->readao.end(); ++raoI) {
	  raoI->addoffset-=positionadditionaloffset;
	}
	CEBUG("ADJUST B3" << endl);// debugDump(false));
	otI->from+=positionadditionaloffset;
	PCR_time_prh_a2b3+=diffsuseconds(tvp);
      }
      PCR_time_prh_a2b+=diffsuseconds(tvt);

    }else{
      CEBUG("Must split");
      splitBin(otI-PCR_offsetmap.begin());
      // recursion: simplest way to rerun placeRead_helper() from start,
      //  will recurse only once anyway
      rpbI=placeRead_helper(rpid,position,dir,urdid,araoindex);
    }
  }
  CEBUG("returning" << endl);
  return rpbI;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * shift all reads at "position" or higher by the given offset
 *
 * This does *NOT* sort the container anew if some reads are suddenly
 *  wrongly placed. E.g., shift(5,-2) on
 *   1,2,3,4,5,6,7,8,9
 *  gives
 *   1,2,3,4,3,4,5,6,7
 * It's up to the caller to make sure the logic is OK.
 *
 * Having negative positions after the shift is possible. E.g.: shift(0,5) on
 *   0,1,2,3,4,5,6,7,8
 *  gives
 *   -5,-4,-3,-2,-1,0,1,2,3
 *
 *************************************************************************/

void PlacedContigReads::shiftReads(int32 position, int32 offsetdiff)
{
  FUNCSTART("void PlacedContigReads::shiftReads(int32 position, int32 offsetdiff)");

  if(offsetdiff==0) return;

  timeval tvp;
  gettimeofday(&tvp,nullptr);


  auto omI=lower_bound(PCR_offsetmap.begin(),PCR_offsetmap.end(),
		       offsettile_t(position,PCR_readposbins.end()),
		       offsettile_t::lt_offsetfrom);
  for(;omI!=PCR_offsetmap.begin(); --omI){
    if((omI-1)->from + (omI-1)->rpbI->readao.back().addoffset < position) break;
  }

  PCR_time_sr_lb1+=diffsuseconds(tvp);

  //if(omI == PCR_offsetmap.end()){
  //  cout << "I'm at PCR_offsetmap.end()\n";
  //}
  if(omI != PCR_offsetmap.end()){
    int targetao=position-omI->from;
    //cout << "My from is " << omI->from << " and targetao is " << targetao << endl;
    auto rpbI=omI->rpbI;

    // Time of this simple loop is the same as for lower_bound(),
    //  but lower_bound should be more cache friendly.
    //auto raoI=rpbI->readao.begin();
    //while(raoI!=rpbI->readao.end() && raoI->addoffset < targetao) ++raoI;

    gettimeofday(&tvp,nullptr);
    auto raoI=lower_bound(rpbI->readao.begin(),rpbI->readao.end(),
			  addoff_t(targetao,0),
			  addoff_t::lt);

    PCR_time_sr_lb2+=diffsuseconds(tvp);

    // even on large containers (>1m reads), this loop takes almost no time
    //  ...
    if(raoI!=rpbI->readao.begin()){
      gettimeofday(&tvp,nullptr);
      for(; raoI!=rpbI->readao.end(); ++raoI){
	raoI->addoffset+=offsetdiff;
      }
      PCR_time_sr_aoadj+=diffsuseconds(tvp);
      ++omI;
    }
    PCR_time_sr_aoadj+=diffsuseconds(tvp);

    // ... when compared to this loop
    //
    // unrolling saves ~5% (from & to, from only 2.5%)
    // unrolling with prefetch saves ~19% (from & to, from only 17 %)
    gettimeofday(&tvp,nullptr);
    for(uint32 count=(PCR_offsetmap.end()-omI)/4; count; --count){
      __builtin_prefetch(&(*(omI+32)), 1, 3);
      omI->from+=offsetdiff;
      ++omI;
      omI->from+=offsetdiff;
      ++omI;
      omI->from+=offsetdiff;
      ++omI;
      omI->from+=offsetdiff;
      ++omI;
    }
    for(; omI!=PCR_offsetmap.end(); ++omI){
      omI->from+=offsetdiff;
    }
    PCR_time_sr_omadj+=diffsuseconds(tvp);

  }

  FUNCEND();
}



/*************************************************************************
 *
 * shift all reads at "position" or higher by the given offset, but bounce
 *  negative offsets back to 0
 *
 * E.g.: shiftReadsBounceZero(0,5) on
 *   0,1,2,3,4,5,6,7,8
 * becomes an intermediate (from shiftReads(0,-5))
 *   -5,-4,-3,-2,-1,0,1,2,3
 * and finally a
 *   0,0,0,0,0,0,1,2,3
 *
 * This is a hacky thing helping to implement cutting back overhanging
 *  front reads in mapping assemblies (Contig::trimMapOverhang())
 *
 * Like for shiftReads(), there is no reordering of the container, so doing
 *  a shiftReadsBounceZero(2,-5) will wreak havoc on
 *   0,1,2,3,4,5,6,7,8
 * producing
 *   0,1,-3,-2,-1,0,1,2,3
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void PlacedContigReads::shiftReadsBounceZero(int32 position, int32 offsetdiff)
{
  FUNCSTART("void PlacedContigReads::shiftReads(int32 position, int32 offsetdiff)");

  //CEBUG("BEFORE1\n");
  //debugDump(false);

  shiftReads(position,offsetdiff);

  //CEBUG("BEFORE2\n");
  //debugDump(false);

  // need to bounce?
  if(PCR_offsetmap.front().from<=0){
    // yes
    CEBUG("BOUNCE!\n");
    auto omI=PCR_offsetmap.begin();
    for(; omI!=PCR_offsetmap.end(); ++omI){     // TODO: eventually more intelligent, earlier stop, but don't care ATM
      for(auto & raoe : omI->rpbI->readao){
	CEBUG("test bounce " << PCR_readdump[raoe.urdid].getName());
	CEBUG("\tof: " << omI->from);
	CEBUG("\tao: " << raoe.addoffset);
	CEBUG(endl);
	if(omI->from+raoe.addoffset < 0){
	  CEBUG("bounceit!\n");
	  raoe.addoffset=-omI->from;
	}
      }
    }
  }

  //CEBUG("AFTER\n");
  //debugDump(false);

  FUNCEND();
}
//#define CEBUG(bla)

/*************************************************************************
 *
 * remove a placed read, return iterator to element after
 *
 *************************************************************************/

PlacedContigReads::const_iterator PlacedContigReads::removeRead(PlacedContigReads::const_iterator pcrI)
{
  FUNCSTART("void removeRead(const_iterator pcrI)");

  BUGIFTHROW(pcrI==end(),"pcrI=end() ??");

  auto retrpbI=pcrI.rpbI;
  auto retaoi=pcrI.raoindex;

  auto urdid=pcrI.rpbI->readao[pcrI.raoindex].urdid;

  // delete orpid from PCR_map*
  delORPIDFromMap(PCR_ancillaryinfo[urdid].orpid);

  // invalidate ancillary info
  PCR_ancillaryinfo[urdid]=ancillaryinfo_t(-1,0,PCR_readposbins.end());

  if(pcrI.rpbI->readao.size()==1){
    auto omI=PCR_offsetmap.begin();
    std::advance(omI,pcrI.rpbI->offsetmapindex);
    PCR_offsetmap.erase(omI);

    // terrible hack: cannot assign a std::list<>::const_iterator to a std::list<>::iterator
    //  though we need that for erase()
    // I don't want to change PlacedContigReads::const_iterator.rpbI to a std::list<>::iterator
    // Therefore:
    // recreate artificially a std::list<>::iterator by doing a list traversal *sigh*
    //  this is so sick ... !
    auto rpbI=PCR_readposbins.begin();
    for(; rpbI!=pcrI.rpbI; ++rpbI) {};

    rpbI=PCR_readposbins.erase(rpbI);
    retrpbI=rpbI;
    retaoi=0;
    for(; rpbI!=PCR_readposbins.end(); ++rpbI){
      rpbI->offsetmapindex-=1;
    }
  }else{
    // more sickness:
    // pcrI.rpbI->readao is constified, this is the easiest way I found to cast the const away
    //  as const_cast<std::vector<addoff_t>>(pcrI.rpbI->readao) does not work
    std::vector<addoff_t> * raovptr=const_cast<std::vector<addoff_t> *>(&pcrI.rpbI->readao);

    if(pcrI.raoindex == 0){
      auto offsetdiff=pcrI.rpbI->readao[1].addoffset;
      PCR_offsetmap[pcrI.rpbI->offsetmapindex].from+=offsetdiff;
      for(auto tmpI=raovptr->begin(); tmpI!=raovptr->end(); ++tmpI){
	tmpI->addoffset-=offsetdiff;
      }
    }

    auto eI=raovptr->begin();
    std::advance(eI,pcrI.raoindex);
    raovptr->erase(eI);

    if(retaoi>=raovptr->size()){
      ++retrpbI;
      retaoi=0;
    }
  }

  // release read in read dump
  PCR_readdump.releaseRead(urdid);

  // if the whole thing gets empty, clear() the PCR to get things clean
  //  (e.g. the ReadContainer PCR_readdump gets faster when re-filled)
  if(--PCR_numreads==0){
    clear();
  }

  FUNCEND();

  return const_iterator(pcrI.pcr,retrpbI,retaoi);
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

PlacedContigReads::const_iterator PlacedContigReads::getPCRIForReadsStartingAtPos(int32 position) const
{
  FUNCSTART("PlacedContigReads::const_iterator PlacedContigReads::getPCRIForReadsStartingAtPos(int32 position) const");

  // TODO: this is very much like the beginning of shiftRead()
  // consolidate?


  // const gymnastics to get this function to compile as "const"
  auto ncthis=const_cast<PlacedContigReads *>(this);

  auto rpbI=ncthis->PCR_readposbins.end();
  const_iterator::raoindex_t araoindex=0;

  auto otI=mstd::lower_bound(ncthis->PCR_offsetmap,
			     offsettile_t(position,ncthis->PCR_readposbins.end()),       // comparator object, only position matters
			     offsettile_t::lt_offsetfrom);
  for(;otI!=PCR_offsetmap.begin(); --otI){
    if((otI-1)->from + (otI-1)->rpbI->readao.back().addoffset < position) break;
  }

  if(otI != PCR_offsetmap.end()){
    int targetao=position-otI->from;
    rpbI=otI->rpbI;
    auto raoI=mstd::lower_bound(rpbI->readao,
				addoff_t(targetao,0),           // comparator object, only position matters
				addoff_t::lt);
    araoindex=static_cast<const_iterator::raoindex_t>(raoI - rpbI->readao.begin());
  }

  return const_iterator(this,rpbI,araoindex);

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

std::ostream & PlacedContigReads::dumpPCRIElement(std::ostream &ostr, const const_iterator & pcrI) const
{
  ostr << pcrI->getLenClippedSeq() << " "
       << pcrI.getReadStartOffset()
       << "-" << pcrI.getReadStartOffset() + pcrI->getLenClippedSeq()-1;
  if(pcrI.getReadDirection() < 0){
    ostr << "\t-";
  }else{
    ostr << "\t+";
  }
  ostr << pcrI->getName();
  ostr << '\t' << pcrI.getORPID();

  return ostr;
}
