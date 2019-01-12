/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 and later by Bastien Chevreux
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 *
 */

#ifdef BOTCHEDEXPLICITINSTANTIATIONLINKER
// the following define is set by hashstats.C when including this file on OSX gcc configuration
// so, the helperdef define is NOT set when this file is compiled by itself
#ifndef BOTCHEDEXPLICITINSTANTIATIONLINKER_HELPERDEF
#define NOCOMPILE
#endif
#endif

#ifdef NOCOMPILE
// do absolutely nothing
#else

#include "errorhandling/errorhandling.H"

#include "util/progressindic.H"

#include "util/dptools.H"
#include "util/stlimprove.H"

#include "mira/seqtohash.H"
#include "mira/hashstats.H"
#include "mira/vhash.H"


using std::cout;
using std::cerr;
using std::endl;


//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif

#ifndef PUBLICQUIET
#define CLOCKSTEPS
#endif

//#define CEBUG(bla)   {cout << bla; cout.flush();}



/*************************************************************************
 *
 * Builds HS_hsv_hashstatnodes and populates .next and .prev, but not
 *  .graphid
 *
 * TODO:
 * kind of dumb to use makeHashStatArrayShortcuts() and low24Bit sorting
 * lexicographical sort would make much more sense (faster inner-4 loop!)
 * but atm: reuse as much as possible
 *
 * Careful here: due to subgraph joining, graphs still may contain
 *  loops at end of this function
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sdbg_populateHashStatNodes()
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_sdbg_populateHashStatNodes()");

  CEBUG("priv_populateHashStatNodes\n");

  HS_hsv_hashstatnodes.clear();
  if(HS_hsv_hashstats.empty()) return;

  for(auto & hse : HS_hsv_hashstats){
    CEBUG("HSE0: " << hse << endl);
  }

  priv_makeHashStatArrayShortcuts();

  for(auto & hse : HS_hsv_hashstats){
    CEBUG("HSEx: " << hse << " " << hash2string(hse.vhash,HS_hs_basesperhash) << endl);
  }

  CEBUG("bph: " << HS_hs_basesperhash << endl);
  HS_hsv_hashstatnodes.resize(HS_hsv_hashstats.size());

  {
    TVHASH_T maxmask(1);
    // *grml* undefined behaviour of left shift for 64 shifts in a 64 bit type makes this cludge necessary
    // the same for 32 shift in 32 bit types etc.pp
    if(HS_hs_basesperhash>=sizeof(TVHASH_T)*4){
      maxmask=0;
    }else{
      auto rollbases=HS_hs_basesperhash;
      maxmask<<=(rollbases*2);
    }
    --maxmask;
    CEBUG("maxmask: " << maxmask << "\t" << hash2string(maxmask, HS_hs_basesperhash) << endl);

    uint64 hsi=0;
    hashstat_t hstmp;

    cout << "Populating HSN:" << endl;
    ProgressIndicator<int64> P(0, HS_hsv_hashstats.size());
    for(auto hsI=HS_hsv_hashstats.begin(); hsI!=HS_hsv_hashstats.end(); ++hsI, ++hsi){
      P.progress(hsi);
      CEBUG("hsI: " << *hsI << endl);
      if(HS_hsv_hashstats[hsi].hsc.iskmerforkf) continue;
      if(HS_hsv_hashstats[hsi].hsc.iskmerforkr) continue;

      // found current start
      auto currenti=hsi;
      while(true){
	  BUGIFTHROW(currenti<0 || currenti>= HS_hsv_hashstats.size(),"Ooops? " << currenti << " currenti<0 || currenti>= HS_hsv_hashstats.size() ???");
	if(HS_hsv_hashstatnodes[currenti].next >= 0) {
	  CEBUG("BREAK: NEXT nonnull\n");
	  break;
	}
	auto vhash=HS_hsv_hashstats[currenti].vhash;
	CEBUG("vh: " << hex << vhash << "\tas: " << hash2string(vhash,HS_hs_basesperhash));
	vhash<<=2;
	vhash&=maxmask;
	CEBUG("\tvh2: " << vhash << "\tas: " << hash2string(vhash,HS_hs_basesperhash) << endl);
	auto numnext=0;
	int64 nextfoundi=-1;
	for(uint8 trials=0;trials<4;++trials,++vhash){
	  CEBUG("Look for: " << vhash << "\tas: " << hash2string(vhash,HS_hs_basesperhash) << endl);
	  hstmp.vhash=vhash;

	  auto hsptr=findVHash(hstmp);
	  if(hsptr!=nullptr) {
	    CEBUG("foundit\n");
	    if(!hsptr->hsc.iskmerforkf
	       && !hsptr->hsc.iskmerforkr){
	      ++numnext;
	      nextfoundi=hsptr-&(HS_hsv_hashstats[0]);
	    }else{
	      CEBUG("fork\n");
	    }
	  }
	}
	CEBUG(dec);

	CEBUG("nn: " << numnext << endl);
	// found exactly 1 successor
	if(numnext==1){
	  BUGIFTHROW(nextfoundi<0 || nextfoundi>= HS_hsv_hashstats.size(),"Ooops? " << nextfoundi << " nextfoundi<0 || nextfoundi>= HS_hsv_hashstats.size() ???");
	  if(HS_hsv_hashstatnodes[nextfoundi].prev >= 0){
	    // oops, potential next already already linked to a previous?
	    // must stop this list
	    numnext=0;
	    CEBUG("non-null " << nextfoundi << ": " << HS_hsv_hashstatnodes[nextfoundi].prev << "\t" << HS_hsv_hashstatnodes[nextfoundi].next << endl);
	  }else{
	    // we don't care whether the .next of the next is already set ... this joins subgraphs :-)
	    //  (and can create loops, but we'll deal with that somewhere else)
	    CEBUG("Stored " << currenti << " <-> " << nextfoundi << endl);
	    HS_hsv_hashstatnodes[nextfoundi].prev=currenti;
	    HS_hsv_hashstatnodes[currenti].next=nextfoundi;
	    currenti=nextfoundi;
	  }
	}
	if(numnext!=1) {
	  CEBUG("BREAKOUT\n");
	  break; // out of the while(true) loop
	}
      }
    }
    P.finishAtOnce(); cout << endl;
  }

  CEBUG("### nodes\n" << dec);
  for(auto hsnI=HS_hsv_hashstatnodes.cbegin(); hsnI!=HS_hsv_hashstatnodes.cend(); ++hsnI){
    CEBUG(hsnI-HS_hsv_hashstatnodes.cbegin() << "\t");
    CEBUG(hsnI->prev << "\t" << hsnI->next);
    CEBUG(endl);
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Needs priv_sdbg_populateHashStatNodes() to be run before
 *
 * Result:
 * - populates .graphid of HS_hsv_hashstatnodes
 * - builds HS_hsv_dbgseqs: .seq, .medianhashcov,
 *                          .hsni_first, .hsni_last and .graphid
 *                  BUT NOT .expectedcopynumber
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sdbg_assignGraphIDsMakeSeqsAndCollectStats()
{
  FUNCSTART("void HashStatistics<TVHASH_T>::priv_sdbg_assignGraphIDsMakeSeqsAndCollectStats");

  CEBUG("start ts" << endl);
  std::vector<uint8> taken(HS_hsv_hashstatnodes.size(),0);
  std::vector<uint8> runvisit(HS_hsv_hashstatnodes.size(),0);    // for breaking loops

  HS_hsv_dbgseqs.clear();
  HS_hsv_seqcontainer.clear();

  hashstat_t hstmp;
  std::vector<uint32> count4median;
  count4median.reserve(512000);

  uint32 graphid=0;

  cout << "Collecting DBG stats:" << endl;
  ProgressIndicator<int64> P(0, HS_hsv_hashstatnodes.size());
  for(int64 hsi=0; hsi<HS_hsv_hashstatnodes.size(); ++hsi){
    P.progress(hsi);
//    for(auto v : runvisit){
//      BUGIFTHROW(v,"runvisit not completely empty??? " << hsi);
//    }
    hstmp.vhash=nsvhash::reverseComplement(HS_hsv_hashstats[hsi].vhash,HS_hs_basesperhash);
    size_t revhsi=-1;
    auto revhsptr=findVHash(hstmp);
    BUGIFTHROW(revhsptr==nullptr,"1) revhsptr==nullptr ??? for " << hash2string(HS_hsv_hashstats[hsi].vhash,HS_hs_basesperhash) << " versus " << hash2string(hstmp.vhash,HS_hs_basesperhash));
    revhsi=(revhsptr-(&HS_hsv_hashstats[0]));

    CEBUG("Check1: #" << hsi << "#" << revhsi << "#" << endl);
    if(!taken[hsi] && !taken[revhsi]){
      auto acti=hsi;
      CEBUG("New start: " << hsi << "\t" << HS_hsv_hashstats[acti]);
      while(HS_hsv_hashstatnodes[acti].prev >= 0
	    && !taken[HS_hsv_hashstatnodes[acti].prev]
	    && !runvisit[HS_hsv_hashstatnodes[acti].prev]) {
	runvisit[acti]=1;
	acti=HS_hsv_hashstatnodes[acti].prev;
	CEBUG("\tmoved to " << acti << "\t" << hash2string(HS_hsv_hashstats[acti].vhash, HS_hs_basesperhash) << "\t" << HS_hsv_hashstats[acti] << "\n");
      }
      CEBUG("\n");
      HS_hsv_seqcontainer.resize(HS_hsv_seqcontainer.size()+1);
      HS_hsv_dbgseqs.resize(HS_hsv_dbgseqs.size()+1);
      HS_hsv_dbgseqs.back().graphid=graphid;
      graphid+=2;  // increase graphid for next loop
      HS_hsv_dbgseqs.back().seqstrptr=&(HS_hsv_seqcontainer.back());
      //auto & actseq=*(HS_hsv_dbgseqs.back().seqstrptr);
      auto assptr=HS_hsv_dbgseqs.back().seqstrptr;
      CEBUG("assptr " << HS_hsv_dbgseqs.back().graphid << " : " << assptr << endl);
      HS_hsv_dbgseqs.back().hsni_first=acti;
      count4median.clear();
      do{
	hstmp.vhash=nsvhash::reverseComplement(HS_hsv_hashstats[acti].vhash,HS_hs_basesperhash);
	revhsptr=findVHash(hstmp);
	BUGIFTHROW(revhsptr==nullptr,"2) revhsptr==nullptr ???");
	revhsi=(revhsptr-(&HS_hsv_hashstats[0]));

	CEBUG("Check2: #" << acti << "#" << revhsi << "#" << endl);
	BUGIFTHROW(taken[acti]!=taken[revhsi],"taken[acti]!=taken[revhsi] ???");

	if(taken[acti]){
	  // this may indeed be legitimate if a reverse complement of
	  //  an already taken kmer is in this graph. Tricky.
	  // in that case ... we end reconstruction of this sequence here
	  // simply break the for loop
	  break;
	}
	taken[acti]=1;
	taken[revhsi]=1;

	HS_hsv_hashstatnodes[acti].graphid=graphid;
	HS_hsv_hashstatnodes[revhsi].graphid=graphid+1;

	HS_hsv_dbgseqs.back().hsni_last=acti;

	CEBUG("gtaken #" << acti << "# (next: " << HS_hsv_hashstatnodes[acti].next << ")\t" << graphid << "\t" << hash2string(HS_hsv_hashstats[acti].vhash,HS_hs_basesperhash) << "\t" << HS_hsv_hashstats[acti] << endl);
	CEBUG("revtaken #" << revhsi << "#" << endl);

	if(unlikely(assptr->empty())){
	  *assptr=hash2string(HS_hsv_hashstats[acti].vhash, HS_hs_basesperhash);
	}else{
	  switch(static_cast<uint64>(HS_hsv_hashstats[acti].vhash) & 3){
	  case 0 : (*assptr)+='A'; break;
	  case 1 : (*assptr)+='C'; break;
	  case 2 : (*assptr)+='G'; break;
	  case 3 : (*assptr)+='T'; break;
	  }
	}
	count4median.push_back(HS_hsv_hashstats[acti].hsc.getCount());

	acti=HS_hsv_hashstatnodes[acti].next;
      }while(acti >= 0);
      CEBUG("slen: " << assptr->size() << endl);
      BUGIFTHROW(assptr->empty(),"assptr->empty() ???");
      // use nth_element for partial sorting to get the median
      std::nth_element(count4median.begin(), count4median.begin()+count4median.size()/2,count4median.end() );
      HS_hsv_dbgseqs.back().medianhashcov=count4median[count4median.size()/2];

      // and reset runvisit to all 0 (can't do that in loop above because of break statement
      acti=hsi;
      while(HS_hsv_hashstatnodes[acti].prev >= 0
	    && runvisit[acti]) {
	runvisit[acti]=0;
	acti=HS_hsv_hashstatnodes[acti].prev;
      }

    }
  }
  for(auto v : runvisit){
    BUGIFTHROW(v,"end check runvisit not completely empty???");
  }

  P.finishAtOnce(); cout << endl;
}
//#define CEBUG(bla)

/*************************************************************************
 *
 * Needs priv_sdbg_populateHashStatNodes() and
 *   priv_sdbg_assignGraphIDsMakeSeqsAndCollectStats() to be run before
 *
 * Result:
 * - populates .expectedcopynumber of HS_hsv_dbgseqs
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sdbg_calcExpectedCopyNumbers()
{
  std::vector<uint32> count4median;
  count4median.reserve(HS_hsv_dbgseqs.size());

  std::vector<uint32> trials{5000,2000,1000,500,200,100};
  trials.push_back(HS_hs_basesperhash);
  uint32 medianlengthused;
  for(auto tn : trials){
    medianlengthused=tn;
    for(auto & dse : HS_hsv_dbgseqs){
      CEBUG("dse.graphid " << dse.graphid << " " << dse.seqstrptr << endl);
      if(dse.seqstrptr->size()>=tn){
	count4median.push_back(dse.medianhashcov);
      }
    }
    if(!count4median.empty()) break;
  }

  if(count4median.empty()) {
    for(auto & dse : HS_hsv_dbgseqs){
      dse.expectedcopynumber=1;
    }
  }else{
    std::nth_element(count4median.begin(), count4median.begin()+count4median.size()/2,count4median.end() );
    auto mhc=count4median[count4median.size()/2];
    cout << "Median >= " << medianlengthused << ": " << mhc << endl;

    for(auto & dse : HS_hsv_dbgseqs){
      dse.expectedcopynumber=static_cast<uint32>((static_cast<double>(dse.medianhashcov)/mhc)+static_cast<double>(.5));
      if(dse.expectedcopynumber<1) dse.expectedcopynumber=1;
      if(dse.expectedcopynumber==1
	 && (HS_hsv_hashstats[dse.hsni_first].hsc.iskmerforkf
	     || HS_hsv_hashstats[dse.hsni_first].hsc.iskmerforkr
	     || HS_hsv_hashstats[dse.hsni_last].hsc.iskmerforkf
	     || HS_hsv_hashstats[dse.hsni_last].hsc.iskmerforkr)){
	dse.expectedcopynumber=2;
      }
    }
  }
}
//#define CEBUG(bla)

/*************************************************************************
 *
 * Needs:
 * - hashstatistics to be loaded
 * - hs must be already trimmed to needs of user
 * - hs must have kmer forks already marked
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::buildSDBGraphs()
{
  FUNCSTART("void HashStatistics<TVHASH_T>::buildSDBGraphs()");

  cout << "Starting priv_sdbg_populateHashStatNodes" << endl;
  priv_sdbg_populateHashStatNodes();
  cout << "Starting priv_sdbg_assignGraphIDsMakeSeqsAndCollectStats" << endl;
  priv_sdbg_assignGraphIDsMakeSeqsAndCollectStats();

  priv_sdbg_calcExpectedCopyNumbers();

#if 0
  cout << "###\n";
  for(auto & dse : HS_hsv_dbgseqs) {
    cout << ">con_" << dse.graphid << " " << dse.seqstrptr->size() << "\t" << dse.medianhashcov
	 << "\t"
	 << dse.expectedcopynumber
	 << "\t"
	 << HS_hsv_hashstats[dse.hsni_first].hsc.iskmerforkr
	 << HS_hsv_hashstats[dse.hsni_first].hsc.iskmerforkf
	 << " "
	 << HS_hsv_hashstats[dse.hsni_last].hsc.iskmerforkf
	 << HS_hsv_hashstats[dse.hsni_last].hsc.iskmerforkr
	 << endl;
    cout << dse.seqstrptr->size() << endl;
  }
#endif
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Builds sequence as string between hashstatnodes hsi1 and hsi2
 * Returns empty sequence if hsi1 and hsi2 are not on the same graph
 *
 * TODO: quite slow. See if can be sped up
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::priv_sdbg_buildSequenceFromPath(int64 hsi1, int64 hsi2, std::string & retseq)
{
  retseq.clear();
  if(HS_hsv_hashstatnodes[hsi1].graphid != HS_hsv_hashstatnodes[hsi2].graphid) return;

  int64 acthsi=hsi1;
  for(;acthsi>=0;){
    CEBUG("HSN: " << HS_hsv_hashstatnodes[acthsi] << endl);
    if(unlikely(retseq.empty())){
      retseq=hash2string(HS_hsv_hashstats[acthsi].vhash,HS_hs_basesperhash);
    }else{
      switch(static_cast<uint64>(HS_hsv_hashstats[acthsi].vhash) & 3){
      case 0 : retseq+='A'; break;
      case 1 : retseq+='C'; break;
      case 2 : retseq+='G'; break;
      case 3 : retseq+='T'; break;
      }
    }
    CEBUG(hsi1 << " " << acthsi << " " << hsi2 << " " << retseq << endl);
    if(acthsi==hsi2) break;
    acthsi=HS_hsv_hashstatnodes[acthsi].next;
  }
  if(acthsi!=hsi2){
    // Ooooooops ... not found
    retseq.clear();
  }

  return;
}
//#define CEBUG(bla)   {cout << bla; cout.flush();}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::dumpDBGSeqs(std::ostream & ostr)
{
  for(auto & dse : HS_hsv_dbgseqs){
    ostr << ">" << dse.graphid
	 << " mhc=" << dse.medianhashcov
	 << " ecn=" << dse.expectedcopynumber
	 << " len=" << dse.seqstrptr->size()
	 << '\n' << *(dse.seqstrptr) << '\n';
  }
}


template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::sortDBGSeqsByLenDown()
{
  mstd::psort(HS_hsv_dbgseqs, sortDBGComparatorLenDown);
}

template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::sortDBGSeqsByMHCDown()
{
  mstd::psort(HS_hsv_dbgseqs, sortDBGComparatorMHCDown);
}


/*************************************************************************
 *
 * Propose edits for a sequence, given back in vector 'edits'
 * Edits are sorted high to low in vector, all positions refer to unedited seq
 *
 * Needs:
 * - buildSDBGraphs() must be run before to fill needed sdbg data structures
 *   (HS_hsv_hashstatnodes)
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void HashStatistics<TVHASH_T>::proposeSDBGEditsForSequence(const void * seqvoid, uint64 slen, const char * namestr, std::vector<dbgedits_t> & edits)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::proposeSDBGEdits()");

  edits.clear();

  // while rare and this should not happen for normal data sets, HS_hsv_hashstatnodes may be empty
  // E.g.: tiny mapping projects
  if(slen<HS_hs_basesperhash
     ||HS_hsv_hashstatnodes.empty()) return;

  CEBUG("testStuff2\n");

  const auto basesperhash=HS_hs_basesperhash;
  const uint8 * seq=static_cast<const uint8 *>(seqvoid);

  int32 lastgoodpos=-1;
  int64 lastgoodhsi=-1;

  std::string repstring;
  repstring.reserve(200);

  hashstat_t hstmp;
  bool insearch=false;
  {
    SEQTOHASH_LOOPSTART(TVHASH_T);
    {
      hstmp.vhash=acthash;
      auto hsptr=findVHash(hstmp);

      {
	int32 actpos=seqi-(basesperhash-1);
	CEBUG("si: " << actpos << "\t" << hash2string(acthash,HS_hs_basesperhash) << "\t");
	if(hsptr==nullptr){
	  CEBUG("not found\n");
	}else{
	  int64 acthsi=hsptr-&(HS_hsv_hashstats[0]);
	  auto actgid=HS_hsv_hashstatnodes[acthsi].graphid;
	  CEBUG("ahsi: " << acthsi << "\tagid: " << actgid << "\t" << HS_hsv_hashstats[acthsi] << endl);
	}
      }

      if(insearch){
	if(hsptr==nullptr){
	  // do nothing?
	}else{
	  // new good found ... go on
	  int64 acthsi=hsptr-&(HS_hsv_hashstats[0]);
	  int32 actpos=seqi-(basesperhash-1);

	  priv_sdbg_buildSequenceFromPath(lastgoodhsi,acthsi,repstring);
	  if(!repstring.empty()){
	    // we could correct!
	    edits.push_back(dbgedits_t(lastgoodpos,actpos-lastgoodpos+basesperhash,repstring));
	  }

	  lastgoodhsi=acthsi;
	  lastgoodpos=actpos;
	  insearch=false;
	}
      }else{
	if(hsptr==nullptr){
	  if(lastgoodpos>=0){
	    insearch=true;
	  }
	}else{
	  lastgoodhsi=hsptr-&(HS_hsv_hashstats[0]);
	  int32 actpos=seqi-(basesperhash-1);
	  lastgoodpos=actpos;
	}
      }
    }
    SEQTOHASH_LOOPEND;
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Checks whether a sequence is a chimera according to SDBG
 * Returns true if yes
 *
 * Needs:
 * - buildSDBGraphs() must be run before to fill needed sdbg data structures
 *   (HS_hsv_hashstatnodes)
 *
 * TODO:
 * Terribly slow because priv_sdbg_checkPathDistance()
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
bool HashStatistics<TVHASH_T>::checkSequenceForSDBGChimeras(const void * seqvoid, uint64 slen, const char * namestr)
{
  FUNCSTART("void HashStatistics<TVHASH_T>::checkSequenceForSDBGChimeras()");

  // while rare and this should not happen for normal data sets, HS_hsv_hashstatnodes may be empty
  // E.g.: tiny mapping projects
  if(slen<HS_hs_basesperhash
     ||HS_hsv_hashstatnodes.empty()) return false;

  CEBUG("chimcheck on " << namestr << endl);

  bool retvalue=false;

  const auto basesperhash=HS_hs_basesperhash;
  const uint8 * seq=static_cast<const uint8 *>(seqvoid);

  int32 lastgoodpos=-1;
  int64 lastgoodhsi=-1;

  std::string repstring;
  repstring.reserve(200);

  hashstat_t hstmp;
  bool insearch=false;
  {
    SEQTOHASH_LOOPSTART(TVHASH_T);
    {
      hstmp.vhash=acthash;
      auto hsptr=findVHash(hstmp);

      {
	int32 actpos=seqi-(basesperhash-1);
	CEBUG("si: " << actpos << "\t" << hash2string(acthash,HS_hs_basesperhash) << "\t");
	if(hsptr==nullptr){
	  CEBUG("not found\n");
	}else{
	  int64 acthsi=hsptr-&(HS_hsv_hashstats[0]);
	  auto actgid=HS_hsv_hashstatnodes[acthsi].graphid;
	  CEBUG("ahsi: " << acthsi << "\tagid: " << actgid << "\t" << HS_hsv_hashstats[acthsi] << endl);
	}
      }

      if(insearch){
	if(hsptr==nullptr){
	  // do nothing?
	}else{
	  // new good found ... go on
	  int64 acthsi=hsptr-&(HS_hsv_hashstats[0]);
	  int32 actpos=seqi-(basesperhash-1);

	  if(!priv_sdbg_checkPathDistance(lastgoodhsi,acthsi,actpos-lastgoodpos+10)){
	    // we have a chimera, let's stop here
	    retvalue=true;
	    break;
	  }

	  lastgoodhsi=acthsi;
	  lastgoodpos=actpos;
	  insearch=false;
	}
      }else{
	if(hsptr==nullptr){
	  if(lastgoodpos>=0){
	    insearch=true;
	  }
	}else{
	  lastgoodhsi=hsptr-&(HS_hsv_hashstats[0]);
	  int32 actpos=seqi-(basesperhash-1);
	  lastgoodpos=actpos;
	}
      }
    }
    SEQTOHASH_LOOPEND;
  }
  return retvalue;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Checks whether a hashstatsnodes hsi1 and hsi2 are at most maxdist away
 * Returns true if yes or if path ends before maxdist was searched
 *
 * Needs:
 * - buildSDBGraphs() must be run before to fill needed sdbg data structures
 *   (HS_hsv_hashstatnodes)
 *
 * TODO:
 *  Terribly slow because needs to walk paths EVERY.SINGLE.TIME!
 *  Solution would be to add a position element to hashstatnode_t
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
bool HashStatistics<TVHASH_T>::priv_sdbg_checkPathDistance(int64 hsi1, int64 hsi2, int64 maxdist)
{
  bool retval=false;

  // check fwd
  int64 acthsi=hsi1;
  uint32 searchdist=0;
  for(; searchdist < maxdist && acthsi>=0; ++searchdist){
    CEBUG("fwalk: " << searchdist << "\t" << hash2string(HS_hsv_hashstats[acthsi].vhash,HS_hs_basesperhash) << "\t" << HS_hsv_hashstats[acthsi] << endl);
    if(acthsi==hsi2) break;
    acthsi=HS_hsv_hashstatnodes[acthsi].next;
  }
  if(acthsi<0
     || (searchdist<maxdist && acthsi==hsi2)) retval=true;
  CEBUG(acthsi << " " << maxdist << " RETVAL " << retval << endl);

  return retval;
}
//#define CEBUG(bla)


// explicit template instantiations needed for the linker to create these for library files
// but not with botched linker
#ifndef BOTCHEDEXPLICITINSTANTIATIONLINKER
template class HashStatistics<vhash64_t>;
#ifndef KMER_INTERNALTYPE
template class HashStatistics<vhash128_t>;
template class HashStatistics<vhash256_t>;
template class HashStatistics<vhash512_t>;
#endif
#endif

#endif // NOCOMPILE
