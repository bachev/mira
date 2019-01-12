/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
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

// functions to process reads
// currently in namespace and object assembly


#include "mira/dataprocessing.H"

#include "mira/assembly.H"
#include "mira/hashstats.H"
#include "mira/parameters.H"
#include "mira/readpool.H"
#include "mira/skim.H"

#include "util/dptools.H"
#include "util/progressindic.H"

#include "util/stlimprove.H"

#include <regex>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using std::cout;
using std::cerr;
using std::endl;


//#define CEBUG(bla)   {if(CEBUGFLAG) {cout << bla; cout.flush();}}
#define CEBUG(bla)



std::string DataProcessing::DP_ggcstring("ggc");

std::vector<DataProcessing::poolskim_t> DataProcessing::DP_adapskims;
boost::mutex DataProcessing::DP_ps_changemutex;

HashStatistics<vhash64_t> DataProcessing::DP_phix174hashstatistics;
HashStatistics<vhash64_t> DataProcessing::DP_rrnahashstatistics;    // maybe dangerous if ppl change the kmer size to >32 for that

bool DataProcessing::DP_px174hs_init=false;
boost::mutex DataProcessing::DP_px174hs_changemutex;         // exclusive mutex for write access to hashstatistics, we need that as it's a static variable

bool DataProcessing::DP_rrnahs_init=false;
boost::mutex DataProcessing::DP_rrnahs_changemutex;         // exclusive mutex for write access to hashstatistics, we need that as it's a static variable


DataProcessing::poolskim_t::~poolskim_t(){
  if(poolptr) delete poolptr;
  if(skimptr) delete skimptr;
}


DataProcessing::DataProcessing(std::vector<MIRAParameters> * params) : DP_miraparams_ptr(params), DP_tmpmtpolyAT(Read::REA_defaulttag_SOFApolyA_sequence)
{
  DP_threadid=-1;
  DP_tmpvu8.reserve(16300); // bit less than 16kb

  // vector with enough capacity so that it does not get reallocated
  // -> multiple threads won't get their data removed under them during the run
  DP_adapskims.reserve(1024);

};



DataProcessing::~DataProcessing()
{
  stopLogging();

//  for(auto & ps : DP_adapskims){
//    if(ps.poolptr!=nullptr) {
//      delete ps.poolptr;
//      delete ps.skimptr;
//    }
//  }
}


/*************************************************************************
 *
 * Needs frequencies to be set in reads
 * returns number of reads marked as chimeras
 * if killreads==true: reads killed (length set to 0, isUsedInAssembly()==false)
 *   else only marked in ischimera
 * no kill if nochimerakill[readindex] is set
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
uint32 DataProcessing::markReadsWithInvalidKMerEndsAsChimeras_Pool(ReadPool & rp, uint32 bph, bool killreads, std::vector<uint8> * debrisreasonptr, std::vector<uint8> & ischimera, std::vector<uint8> & nochimerakill, const std::string & logprefix)
{
  uint32 retvalue=0;
  for(readid_t rpi=0; rpi<rp.size(); ++rpi){
    Read & actread=rp[rpi];

    CEBUG("CHIMP 1:"
	  << " " << actread.hasValidData()
	  << " " << actread.isUsedInAssembly()
	  << "\t" << actread.getName() << endl);
    CEBUG(actread);

    if(actread.hasValidData()
       && actread.isUsedInAssembly()
       && !(actread.isBackbone()
	    || actread.isRail())){
      CEBUG("CHIMP 2 " << actread.getName() << endl);
      if(checkReadForInvalidKMerEndsAsChimera(actread,bph)){
	CEBUG("CHIMP 3 " << actread.getName() << endl);
	++retvalue;

	bool killthis=killreads && (nochimerakill.empty() || !nochimerakill[rpi]);
	auto oldclip=actread.getRQClipoff();
	if(killthis){
	  actread.setRQClipoff(0);
	  actread.setUsedInAssembly(false);
	}

	if(!ischimera[rpi]){
	  DP_logfout << logprefix
		     << " incorrectible end or chimera. ";
	  if(killthis){
	    DP_logfout << "Kill "
		     << killthis << " " << killreads << " "
		     << actread.getName() << '\t'
		     << oldclip << " -> " << actread.getRightClipoff() << '\n';
	  }else{
	    DP_logfout << "Mark "
		       << actread.getName() << '\n';
	  }
	}
	ischimera[rpi]=true;

	if(debrisreasonptr!=nullptr) (*debrisreasonptr)[rpi]=Assembly::DEBRIS_CLIP_INCORRECTIBLEENDORCHIMERA;
      }
    }
  }

  return retvalue;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Needs frequencies to be set in read
 * returns true if considered chimera
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool DataProcessing::checkReadForInvalidKMerEndsAsChimera(Read & actread, uint32 bph)
{
  bool retvalue=false;
  CEBUG("CHIMR 1:"
	<< " " << actread.hasValidData()
	<< " " << actread.getLenClippedSeq()
	<< " " << actread.getRightClipoff()
	<< " " << actread.hasBaseHashStats()
	<< "\t" << actread.getName() << endl);

  if(!actread.hasValidData()
     || actread.getLenClippedSeq()==0
     || actread.getRightClipoff()==0
     || !actread.hasBaseHashStats()) return false;

  // no chimera if read < kmer size
  // or if kmer does not cover at least 50%
  if(actread.getLenClippedSeq()<bph
     || 100.0/actread.getLenClippedSeq()*bph < 50.0) return false;

  bool hasvalid=false;
  auto bphI=actread.getBPosHashStats().begin();
  bphI+=actread.getLeftClipoff();
  auto bphE=bphI+actread.getLenClippedSeq()-bph+1;
  for(; bphI!=bphE; ++bphI){
    if(bphI->fwd.getFrequency()>1){
      hasvalid=true;
      break;
    }
  }

  CEBUG("CHIMR 2:"
	<< " " << hasvalid
	<< "\t" << actread.getName() << endl);

  bool maybechimera=false;
  if(hasvalid){
    CEBUG("CHIMR 3:"
	  << " " << static_cast<uint16>(actread.getBPosHashStats()[actread.getLeftClipoff()].fwd.getFrequency())
	  << " " << static_cast<uint16>(actread.getBPosHashStats()[actread.getRightClipoff()-1].rev.getFrequency())
	  << "\t" << actread.getName() << endl);

    if(actread.getBPosHashStats()[actread.getLeftClipoff()].fwd.getFrequency()<2
       || actread.getBPosHashStats()[actread.getRightClipoff()-1].rev.getFrequency()<2){
      maybechimera=true;
    }else{
      // quick way did not find chimera, do it the tedious way: invalid on fwd/rev at same pos
      // or maybe not atm ...
    }
  }

  if(!hasvalid || maybechimera){
    retvalue=true;
  }

  return retvalue;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * Meaner than the version looking just at ends. Should probably be called
 *  only on the last pass
 *
 * Needs frequencies to be set in reads
 * returns number of reads clipped and marked as chimeras
 * no length check of reads if bph==0
 * no kill if nochimerakill[readindex] is set
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
uint32 DataProcessing::markReadsWithRareKMersAsChimeras_Pool(ReadPool & rp, uint32 bph, bool killreads, std::vector<uint8> * debrisreasonptr, std::vector<uint8> & ischimera, std::vector<uint8> & nochimerakill, const std::string & logprefix)
{
  uint32 retvalue=0;
  for(readid_t rpi=0; rpi<rp.size(); ++rpi){
    Read & actread=rp[rpi];

    CEBUG("RKACp 1:"
	  << " " << actread.hasValidData()
	  << " " << actread.isUsedInAssembly()
	  << "\t" << actread.getName() << endl);
    CEBUG(""; Read::setCoutType(Read::AS_TEXT); cout << actread);

    if(actread.hasValidData()
       && actread.isUsedInAssembly()
       && !(actread.isBackbone()
	    || actread.isRail())){
      CEBUG("RKACp 2 " << actread.getName() << endl);
      if(checkReadForRareKMersAsChimera(actread,bph)){
	CEBUG("RKACp 3 " << actread.getName() << endl);
	++retvalue;

	bool killthis=killreads && (nochimerakill.empty() || !nochimerakill[rpi]);
	auto oldclip=actread.getRQClipoff();
	if(killthis){
	  actread.setRQClipoff(0);
	  actread.setUsedInAssembly(false);
	}

	if(!ischimera[rpi]){
	  DP_logfout << logprefix
		     << " terminally incorrectible end or chimera. ";
	  if(killthis){
	    DP_logfout << "Kill "
		     << killthis << " " << killreads << " "
		     << actread.getName() << '\t'
		     << oldclip << " -> " << actread.getRightClipoff() << '\n';
	  }else{
	    DP_logfout << "Mark "
		       << actread.getName() << '\n';
	  }
	}
	ischimera[rpi]=true;

	if(debrisreasonptr!=nullptr) (*debrisreasonptr)[rpi]=Assembly::DEBRIS_CLIP_TERMINALLYINCORRECTIBLEORCHIMERA;
      }
    }
  }

  return retvalue;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * Needs frequencies to be set in read
 * returns true if chimera
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool DataProcessing::checkReadForRareKMersAsChimera(Read & actread, uint32 bph)
{
  CEBUG("RKACr 1:"
	<< " " << actread.hasValidData()
	<< " " << actread.getLenClippedSeq()
	<< " " << actread.getRightClipoff()
	<< " " << actread.hasBaseHashStats()
	<< "\t" << actread.getName() << endl);

  if(!actread.hasValidData()
     || actread.getLenClippedSeq()==0
     || actread.getRightClipoff()==0
     || !actread.hasBaseHashStats()) return false;

  // no chimera if read < kmer size
  if(actread.getLenClippedSeq()<bph) return false;

  bool hasrare=false;
  auto bphI=actread.getBPosHashStats().begin();
  bphI+=actread.getLeftClipoff();
  auto bphE=bphI+actread.getLenClippedSeq()-bph+1;
  for(; bphI!=bphE; ++bphI){
    if(!bphI->fwd.isValid() || bphI->fwd.getFrequency()<2){
      hasrare=true;
      break;
    }
  }

  return hasrare;
}
//#define CEBUG(bla)





/*************************************************************************
 *
 * Note: hashstatistics will be modified by this if trimfreq>0
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
uint32 DataProcessing::performSDBGChimeraSearch_Pool(ReadPool & rp, HashStatistics<TVHASH_T> & hsd, uint32 trimfreq, std::vector<uint8> * debrisreasonptr, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::performSDBGChimeraSearch_Pool(ReadPool & rp, HashStatistics<TVHASH_T> & hsd, uint32 trimfreq, std::vector<uint8> * debrisreasonptr, const std::string & logprefix)");

  uint32 retvalue=0;

  if(debrisreasonptr != nullptr){
    BUGIFTHROW(debrisreasonptr->size() != rp.size(),"debrisreasonptr->size() != rp.size()");
  }

  if(trimfreq>0){
    hsd.trimHashStatsByFrequencyAND(trimfreq,trimfreq,0);
    hsd.calcKMerForks(trimfreq,true);
  }

  hsd.buildSDBGraphs();

  ProgressIndicator<int32>  pi(0,rp.size());
  for(readid_t rpi=0; rpi<rp.size(); ++rpi){
    pi.increaseprogress();
    Read & actread=rp[rpi];

    if(actread.hasValidData()
//	   && AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_clip_proposeendclips
       && !(actread.isBackbone()
	    || actread.isRail())){
      if(performSDBGChimeraSearch_Read(actread,hsd,logprefix)){
	++retvalue;
	if(debrisreasonptr!=nullptr) (*debrisreasonptr)[rpi]=Assembly::DEBRIS_CLIP_CHIMERA;
	actread.setUsedInAssembly(false);
      }
    }
  }

  pi.finishAtOnce();
  return retvalue;
}


template<typename TVHASH_T>
bool DataProcessing::performSDBGChimeraSearch_Read(Read & actread, HashStatistics<TVHASH_T> & hsd, const std::string & logprefix)
{
  auto res=hsd.checkSequenceForSDBGChimeras(actread.getClippedSeqAsChar(),
					    actread.getLenClippedSeq(),
					    actread.getName().c_str());
  if(res){
    DP_logfout << logprefix
	       << " SDBG chimera kill "
	       << actread.getName() << '\t'
	       << actread.getRightClipoff() << " -> 0\n";
    actread.setRQClipoff(0);
    // TODO: maybe set tags?
  }
  return res;
}



/*************************************************************************
 *
 * Note: hashstatistics will be modified by this if trimfreq>0
 *
 * TODO: slow, should be parallelised. In the outer for loop.
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
uint32 DataProcessing::performSDBGEdits_Pool(ReadPool & rp, HashStatistics<TVHASH_T> & hsd, uint32 trimfreq)
{
  FUNCSTART("uint32 DataProcessing::performSDBGEdits_Pool(ReadPool & rp, HashStatistics<TVHASH_T> & hsd, uint32 trimfreq)");

  uint32 retvalue=0;
  uint32 totaledits=0;
  uint32 numfishyedits=0;

  if(trimfreq>0){
    hsd.trimHashStatsByFrequencyAND(trimfreq,trimfreq,3); // but at least 3 kmers
    hsd.calcKMerForks(trimfreq,true);
  }

  hsd.buildSDBGraphs();

  std::vector<typename HashStatistics<TVHASH_T>::dbgedits_t> edits;
  ProgressIndicator<int32>  pi(0,rp.size());
  for(readid_t rpi=0; rpi<rp.size(); ++rpi){
    pi.increaseprogress();
    Read & actread=rp[rpi];

    if(actread.hasValidData()
//	   && AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_clip_proposeendclips
       && !(actread.isBackbone()
	    || actread.isRail())){

      hsd.proposeSDBGEditsForSequence(actread.getClippedSeqAsChar(),
				      actread.getLenClippedSeq(),
				      actread.getName().c_str(),
				      edits);
      bool hasedits=false;
      if(!edits.empty()){
	//cout << "Have " << edits.size() << " edits for " << actread.getName() << endl;
	for(auto erI=edits.rbegin(); erI!=edits.rend(); ++erI){
	  auto lendiff=static_cast<int32>(erI->len)-static_cast<int32>(erI->replacement.size());
	  if(abs(lendiff)>2){
	    // TODO: really stay silent?
	    //cout << "Fishy edit: length " << lendiff << " " << actread.getName() << " " << erI->pos << " " << erI->pos+erI->len << " " << erI->len;
	    //std::string olds(static_cast<std::string>(actread.getClippedSeqAsChar()).substr(erI->pos,erI->len));
	    //cout << "\t" << olds
	    //	 << "\t" << erI->replacement << endl;
	    ++numfishyedits;
	  }else{
	    //cout << "Could correct " << erI-edits.rbegin() << " " << actread.getName() << " " << erI->pos << " " << erI->pos+erI->len << " " << erI->len;
	    //std::string olds(static_cast<std::string>(actread.getClippedSeqAsChar()).substr(erI->pos,erI->len));
	    //cout << "\n" << olds
	    //     << "\n" << erI->replacement << endl;
	    //Read::setCoutType(Read::AS_TEXT);
	    //cout << "OLD:\n" << actread << endl;
	    actread.clearTags();
	    actread.smoothSequenceReplace(erI->pos+actread.getLeftClipoff(),erI->len,erI->replacement);
	    //cout << "NEW:\n" << actread << endl;
	    ++totaledits;
	    hasedits=true;
	  }
	}
	if(hasedits) ++retvalue;
      }
    }
  }

  pi.finishAtOnce();

  cout << "Performed " << totaledits << " edits.\n";
  cout << "Rejected " << numfishyedits << " fishy edits.\n";

  return retvalue;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
template<typename TVHASH_T>
void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics<TVHASH_T> & hsd, std::vector<uint8> * debrisreasonptr)
{
  FUNCSTART("void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics & hsd)");

  uint32 numtaken=0;
  uint32 numnormout=0;

  std::vector<bool> goodread(rp.size(),false);  // reads which pass test: all freq >= 2, fwd, rev and no N
  for(uint32 rpi=0; rpi<rp.size(); ++rpi){
    auto & actread=rp[rpi];
    if(!actread.hasValidData()
       || !actread.isUsedInAssembly()
       || actread.isRail()
       || actread.isBackbone()) continue;
    bool isgood=true;
    auto bhsI=actread.getBPosHashStats().begin();
    bhsI+=actread.getLeftClipoff();
    for(auto ri=0; ri<actread.getLenClippedSeq(); ++ri){
      if(bhsI->fwd.getFrequency()<2 || !bhsI->fwd.hasConfirmedFwdRev()) {
	isgood=false;
	break;
      }
    }
    // getClippedSeqAsChar may throw if the MIRA clipping set the left cutoff to the length of the sequence
    // too lazy to get things otherwise;
    if(isgood && actread.getLenClippedSeq()>0){
      auto sptr=actread.getClippedSeqAsChar();
      auto eptr=sptr+actread.getLenClippedSeq();
      for(; sptr!=eptr; ++sptr){
	if(toupper(*sptr)=='N') {
	  isgood=false;
	  break;
	}
      }
    }
    goodread[rpi]=isgood;
  }

  std::vector<bool> normdone(rp.size(),false);
  std::vector<bool> normout(rp.size(),false);
  std::vector<bool> normthisrg(rp.size(),false);


  // do the normalisation for every readgroup so that we independently get reads from every rg
  for(auto rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    // auto rgid=ReadGroupLib::getReadGroupID(rgi);

    hsd.digiNormReset();
    mstd::fill(normthisrg,false);
    cout << "\nReadgroup " << rgi << ":\n";
    ProgressIndicator<int32>  pi(0,rp.size()*3); // three steps, see below

    // step 0 : pairs where both reads are good
    // step 1 : reads (unpaired or paired) which are good
    // step 2 : remaining reads
    for(uint32 step=0; step<3; ++step){
      for(uint32 rpi=0; rpi<rp.size(); ++rpi){
	pi.increaseprogress();
	if(normdone[rpi]) continue; // checking "|| normout[rpi]" not needed because isUsedInAssembly() will be false below anyway
	auto & actread=rp[rpi];
	if(!actread.hasValidData()
	   || !actread.isUsedInAssembly()
	   || actread.isRail()
	   || actread.isBackbone()) continue;
	//Read::setCoutType(Read::AS_TEXT);
	//cout << "### bla\n";
	//cout << actread << endl;
	bool lookatread=goodread[rpi];
	if(step==0){
	  lookatread=false;
	  if(actread.getTemplatePartnerID()>=0){
	    lookatread=goodread[actread.getTemplatePartnerID()];
	  }
	}else if(step==1){
	  // further tests? I do not think so
	}else if(step==2){
	  lookatread=true;
	}else{
	  BUGIFTHROW(true,"Oooops, step " << step << " not foreseen?");
	}
	if(lookatread){
	  bool taken1=hsd.digiNormTestRead(actread,false);
	  bool taken2=false;
	  if(step==0 && actread.getTemplatePartnerID()>=0){
	    auto actmate=rp[actread.getTemplatePartnerID()];
	    taken2=hsd.digiNormTestRead(actmate,false);
	  }
	  if(taken1 || taken2){
	    if(taken1){
	      ++numtaken;
	      CEBUG("Kept " << actread.getName() << endl);
	      normthisrg[rpi]=true;
	      normdone[rpi]=true;
	      for(auto & te : const_cast<std::vector<multitag_t> &>(actread.getTags())){
		if(te.identifier == Read::REA_tagentry_idMNRr) te.identifier = Read::REA_tagentry_idDGNr;
	      }
	    }
	    if(taken2){
	      ++numtaken;
	      auto mrpi=actread.getTemplatePartnerID();
	      auto actmate=rp[mrpi];
	      CEBUG("Kept " << actmate.getName() << endl);
	      normthisrg[mrpi]=true;
	      normdone[mrpi]=true;
	      for(auto & te : const_cast<std::vector<multitag_t> &>(actmate.getTags())){
		if(te.identifier == Read::REA_tagentry_idMNRr) te.identifier = Read::REA_tagentry_idDGNr;
	      }
	    }
	  }else{
	    ++numnormout;
	    CEBUG("NormOut " << actread.getName() << endl);
	    normout[rpi]=true;
	    actread.setRQClipoff(0);
	    actread.setUsedInAssembly(false);
	  }
	}
      }
    }
    pi.finishAtOnce();

    cout << "Calculating replacement coverage";
    uint32 chkall=0;
    for(uint32 rpi=0; rpi<rp.size(); ++rpi){
      if(normthisrg[rpi]){
	auto & actread=rp[rpi];
	auto & tv = const_cast<std::vector<multitag_t> &>(actread.getTags());
	for(auto & te : tv){
	  if(te.identifier == Read::REA_tagentry_idDGNr){
	    CEBUG(actread.getName() << "\tDGNr: " << te.to-te.from+1 << "\t" << actread.getLenClippedSeq() << endl);
	    CEBUG(actread.getName() << "\tLIML: " << actread.getLeftClipoff() << '\t' << te.from << endl;);
	    CEBUG(actread.getName() << "\tLIMR: " << actread.getRightClipoff() << '\t' << te.to << endl);
	    //Read::setCoutType(Read::AS_TEXT);
	    //CEBUG(actread);

	    auto perccovered=static_cast<uint8>(100.0f/actread.getLenClippedSeq()*(te.to-te.from+1));
	    if(perccovered>=80){
	      ++chkall;
	      CEBUG("Next read:\n");
	      auto repcov=hsd.estimDigiNormCov(actread);
	      CEBUG("repcov: " << repcov << endl);
	      if(repcov>1){
		auto newtag=te;
		// WARNING: with this we'll probably break the for(auto & te ...) functionality
		//  (depending on the container it is in)
		// we MUST get out of the loop afterwards with a "break"!
		actread.deleteTag(Read::REA_tagentry_idDGNr);
		std::string comnum(boost::lexical_cast<std::string>(repcov));
		newtag.setCommentStr(comnum);
		newtag.commentisgff3=false;
		actread.addTagO(newtag);
		CEBUG("DGN repcov " << actread.getName() << ": " << repcov << endl; Read::setCoutType(Read::AS_TEXT););
		CEBUG(actread);
		break;
	      }
	    }
	  }
	}
      }
    }
    CEBUG("Chkall: " << chkall << endl);
  }

  cout << "\nDigital normalisation: removed " << numnormout << " reads.\n";

  if(debrisreasonptr != nullptr){
    auto & db=*debrisreasonptr;
    BUGIFTHROW(db.size()!=rp.size(),"db.size()!=rp.size() ???");
    for(auto rpi=0; rpi<rp.size(); ++rpi){
      if(db[rpi]==0 && normout[rpi]){
	db[rpi]=Assembly::DEBRIS_DIGITAL_NORMALISATION;
      }
    }
  }

}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsureAdapRegexes(ReadGroupLib::ReadGroupID rgid)
{
  if(DP_adapres.size()>rgid.getLibId()
     && DP_adapres[rgid.getLibId()].areinit) return;

  if(DP_adapres.size()<=rgid.getLibId()) DP_adapres.resize(rgid.getLibId()+1);

  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){
    static const char regexfile[] = {
#include "adaptorsregex.solexa.xxd.H"
      ,0
    };
    addAdapRegexes(rgid,regexfile);
  }else if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){
    static const char regexfile[] = {
#include "adaptorsregex.iontor.xxd.H"
      ,0
    };
    addAdapRegexes(rgid,regexfile);
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::addAdapRegexes(ReadGroupLib::ReadGroupID rgid, const char * regexfile)
{
  FUNCSTART("void priv_constructorAdapRes(uint8 seqtype, const char * regexfile)");

  BUGIFTHROW(rgid.getLibId()>=ReadGroupLib::getNumReadGroups(),"Oooops, readgroupid " << static_cast<uint16>(rgid.getLibId()) << " is unknown?");

  if(DP_adapres.size()<=rgid.getLibId()) DP_adapres.resize(rgid.getLibId()+1);

  CEBUG("prepping regexp for " << rgid.getLibId() << endl);

  masterslavere_t tmpmsre;
  std::istringstream tmpis(regexfile);
  std::string line;
  while(true){
    getline(tmpis,line);
    if(tmpis.eof()) break;
    if(line[0]=='>'){
      DP_adapres[rgid.getLibId()].adapres.push_back(tmpmsre);
      line.erase(0,1);         // get away the ">"
      boost::trim(line);
      if(!line.empty()){
	boost::to_upper(line);
	DP_adapres[rgid.getLibId()].adapres.back().masterre=std::regex(line);
	DP_adapres[rgid.getLibId()].adapres.back().hasmaster=true;
      }
    }else{
      BUGIFTHROW(DP_adapres[rgid.getLibId()].adapres.empty(),"Oooops, no master expression found?");
      boost::to_upper(line);
      DP_adapres[rgid.getLibId()].adapres.back().slaveres.push_back(std::regex(line));
    }
  }
  DP_adapres[rgid.getLibId()].areinit=true;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsureAdapSkims(ReadGroupLib::ReadGroupID rgid)
{
  if(DP_adapskims.size()>rgid.getLibId()
     && DP_adapskims[rgid.getLibId()].skimptr) return;

  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){
    static const char adapfile[] = {
#include "adaptorsforclip.solexa.xxd.H"
      ,0
    };
    priv_constructorSkimPool(rgid,DP_adapskims,7,adapfile);
  }else if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){
    static const char adapfile[] = {
#include "adaptorsforclip.iontor.xxd.H"
      ,0
    };
    priv_constructorSkimPool(rgid,DP_adapskims,7,adapfile);
  }else if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_454GS20){
    static const char adapfile[] = {
#include "adaptorsforclip.454.xxd.H"
      ,0
    };
    priv_constructorSkimPool(rgid,DP_adapskims,7,adapfile);
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::priv_constructorSkimPool(ReadGroupLib::ReadGroupID rgid, std::vector<poolskim_t> & skimpool, const uint32 basesperhash, const char * adapfile)
{
  FUNCSTART("void DataProcessing::priv_constructorAdapPool(uint8 seqtype, const char * adapfile)");
  BUGIFTHROW(rgid.getLibId()>=ReadGroupLib::getNumReadGroups(),"Oooops, readgroupid " << static_cast<uint16>(rgid.getLibId()) << " is unknown, have " << ReadGroupLib::getNumReadGroups() << " read groups ?");

  boost::mutex::scoped_lock lock(DP_ps_changemutex);
  if(skimpool.size()<=rgid.getLibId()) skimpool.resize(rgid.getLibId()+1);

  if(skimpool[rgid.getLibId()].poolptr==nullptr){
    CEBUG("prepping adappool " << rgid.getLibId() << endl);
    skimpool[rgid.getLibId()].poolptr= new ReadPool;

    std::istringstream tmpis(adapfile);
    std::string nline,sline;
    bool addedsomething=false;
    while(true){
      getline(tmpis,nline);
      if(tmpis.eof()) break;
      boost::trim(nline);
      nline.erase(0,1);         // get away the ">"
      if(!nline.empty()){
	getline(tmpis,sline);
	boost::trim(sline);
	boost::to_upper(sline);
	size_t ereadidx=0;
	for( ; ereadidx < skimpool[rgid.getLibId()].poolptr->size(); ++ereadidx){
	  //CEBUG("COMP: " << skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getName()<< "\t" << skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getSeqAsChar() << "\t" << sline << endl);
	  if(skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getLenSeq() == sline.size()
	     && strncmp(skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getSeqAsChar(),sline.c_str(),sline.size())==0) {
	    //CEBUG("BINGO!\n");
	    break;
	  }
	}
	if(ereadidx==skimpool[rgid.getLibId()].poolptr->size()){
	  CEBUG("New for " << nline<<endl);
	  ereadidx=skimpool[rgid.getLibId()].poolptr->provideEmptyRead();
	  Read & actread=skimpool[rgid.getLibId()].poolptr->getRead(ereadidx);
	  actread.disallowAdjustments();
	  // anonymous read, don't give it a name
	  // This was for mod_sqt being able to trash the read name string container
	  // if the read had a name, sigfault or similar ensues.
	  //actread.setName(nline);
	  if(tmpis.eof()) break;
	  actread.setSequenceFromString(sline);
	  addedsomething=true;
	  actread.setUsedInAssembly(false);
	}else{
	  Read & actread=skimpool[rgid.getLibId()].poolptr->getRead(ereadidx);

	  CEBUG("Extend for " << actread.getName() << ": " << nline<<endl);

	  // anonymous read, don't give it a name
	  // This was for mod_sqt being able to trash the read name string container
	  // if the read had a name, sigfault or similar ensues.
	  //actread.setName(actread.getName()+"_/_"+nline);
	  //actread.setUsedInAssembly(true);
	}
      }
    }

    //for(auto idx=0; idx < skimpool[rgid.getLibId()].poolptr->size(); ++idx){
    //  Read & actread=skimpool[rgid.getLibId()].poolptr->getRead(idx);
    //  if(actread.isUsedInAssembly()){
    //	cout << actread.getName() << endl;
    //  }
    //}
    //exit(0);

    if(addedsomething){
      CEBUG("prepping skim " << rgid.getLibId() << endl);
      if(skimpool[rgid.getLibId()].skimptr!=nullptr) delete skimpool[rgid.getLibId()].skimptr;

      // IMPORTANT: keep assignment of skimpool[rgid.getLibId()].skimptr
      //  as last thing to be done in this branch, i.e., when the Skim is
      //  completely initialised
      // Reason: the checks in priv_ensureAdapSkims() do not use any mutexes, just
      //  the vector size AND the skimptr being nullptr. It could be that
      //  the constructing thread therefore did not finish constructing
      //  the Skim but that another thread would jump forward using it because
      //  everything points to that it's ready
      auto skimptr = new Skim<vhash64_t>();
      skimptr->skimStreamPrepare(*skimpool[rgid.getLibId()].poolptr,basesperhash,1);
      skimptr->prepareForMultithreadFarc((*DP_miraparams_ptr)[0].getSkimParams().sk_numthreads);
      skimpool[rgid.getLibId()].skimptr = skimptr;

      CEBUG("Done\n");
    }
  }
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsurePhiX174Statistics()
{
  if(DP_px174hs_init) return;

  boost::mutex::scoped_lock lock(DP_px174hs_changemutex);
  if(DP_px174hs_init) return;

  ReadPool baitrp;

  if(addPhiX174ToReadpool(baitrp)){
    std::string dummyfn((*DP_miraparams_ptr)[0].getDirectoryParams().dir_tmp+"/phix174.mhs.gz");
    DP_phix174hashstatistics.computeHashStatistics(baitrp,512,false,false,true,1,0,31,dummyfn,(*DP_miraparams_ptr)[0].getDirectoryParams().dir_tmp);
    DP_px174hs_init=true;
    CEBUG("Done\n");
  }

}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool DataProcessing::addPhiX174ToReadpool(ReadPool & baitrp)
{
  static const char adapfile[] = {
#include "seqforfilter_phix174.solexa.xxd.H"
    ,0
  };

  std::istringstream tmpis(adapfile);
  std::string line;
  bool addedsomething=false;
  while(true){
    getline(tmpis,line);
    if(tmpis.eof()) break;
    line.erase(0,1);         // get away the ">"
    if(!line.empty()){
      size_t ereadidx=baitrp.provideEmptyRead();
      Read & actread=baitrp[ereadidx];
      actread.disallowAdjustments();
      // anonymous read, don't give it a name
      // This was for mod_sqt being able to trash the read name string container
      // if the read had a name, sigfault or similar ensues.
      //actread.setName(line);
      getline(tmpis,line);
      //CEBUG("For " << actread.getName() << ": " << line<<endl);
      if(tmpis.eof()) break;
      actread.setSequenceFromString(line);
      addedsomething=true;
    }
  }

  // honestly, this should always be true;
  return addedsomething;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsureRRNAStatistics()
{
  if(DP_rrnahs_init) return;

  boost::mutex::scoped_lock lock(DP_rrnahs_changemutex);
  if(DP_rrnahs_init) return;

  std::string dummyfn(MIRAParameters::getMHSLibDir()+"/filter_default_rrna.mhs.gz");
  DP_rrnahashstatistics.loadHashStatistics(dummyfn);
  DP_rrnahs_init=true;
  CEBUG("Done\n");
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::startLogging(const std::string & filename, bool newfile)
{
  FUNCSTART("void DataProcessing::startLogging(const std::string & filename, bool newfile)");
  stopLogging();
  if(!filename.empty()){
    DP_logname=filename;
    if(newfile){
      DP_logfout.open(filename, std::ios::out|std::ios::trunc);
    }else{
      DP_logfout.open(filename, std::ios::out|std::ios::app);
    }
    if(!DP_logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open " << filename << " for logging.");
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::stopLogging()
{
  if(DP_logfout.is_open()){
    DP_logfout.close();
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::performRareKMERMasking_Pool(ReadPool & rpool, uint32 basesperhash, const std::string & logprefix)
{
  auto & miraparams = *DP_miraparams_ptr;

  bool needtomask=false;
  for(auto & mp : miraparams){
    if(mp.getAssemblyParams().as_clipmask_rarekmers){
      needtomask=true;
    }
  }

  if(!needtomask) return;

  cout << "Rare kmer masking ... ";cout.flush();
  for(uint32 rpi=0; rpi<rpool.size(); ++rpi){
    Read & actread= rpool[rpi];
    if(!actread.hasValidData()
       || !actread.isUsedInAssembly()
       || actread.isBackbone()
       || actread.isRail()
       || !miraparams[actread.getSequencingType()].getAssemblyParams().as_clipmask_rarekmers ) continue;
    performRareKMERMasking_Read(actread,basesperhash, logprefix);

  }

  cout << "done\n";
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::performRareKMERMasking_Read(Read & actread, uint32 basesperhash, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::performRareKMERMasking_Read(Read & actread, uint32 basesperhash, const std::string & logprefix)");

  if(actread.getLenSeq()==0) return;

  if(!actread.hasBaseHashStats()){
    Read::setCoutType(Read::AS_TEXT);
    cout << actread << endl;
    MIRANOTIFY(Notify::FATAL,"!actread.hasBaseHashStats() ??? ");
  }

  DP_tmpvu8.clear();

  Read::setCoutType(Read::AS_TEXT);
  CEBUG("rare kmer, bph " << static_cast<uint16>(basesperhash) << " looking at:" << actread << endl);
  bool mustmask=false;

  if(actread.getLenClippedSeq()<basesperhash){
    CEBUG("complete kill " << actread.getName() << endl);
    actread.setRQClipoff(actread.getLQClipoff());
  }else{
    DP_tmpvu8.resize(actread.getLenSeq(),0);

    auto bhsI=actread.getBPosHashStats().cbegin()+actread.getLeftClipoff();
    auto bhsE=actread.getBPosHashStats().cend();
    auto tfI=DP_tmpvu8.begin()+actread.getLeftClipoff();
    auto tfE=DP_tmpvu8.end();

    for(uint32 readpos=actread.getLeftClipoff(); readpos<actread.getRightClipoff()-(basesperhash-1); ++bhsI, ++tfI, ++readpos){
      if(bhsI->fwd.getFrequency()==1
	 || !bhsI->fwd.hasConfirmedFwdRev()){
	auto ttfI=tfI;
	uint32 tmpreadpos=readpos;
	CEBUG("new mask: " << tmpreadpos << "..");
	for(uint32 i=0; i<basesperhash && ttfI!=tfE; ++i, ++ttfI, ++tmpreadpos){
	  if(tmpreadpos>=actread.getLeftClipoff() && tmpreadpos<actread.getRightClipoff()){
	    mustmask=true;
	    *ttfI=1;
	  }
	}
	CEBUG(tmpreadpos-1 << endl);
      }
    }

    if(mustmask){
      // Gnaw away what was masked too much
      // also turns
      // .............xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx...........
      // into
      // .............xxxxx.............................................xxxxx...........
      if(1){
	auto bhsI=actread.getBPosHashStats().cbegin()+actread.getLeftClipoff();
	auto bhsE=actread.getBPosHashStats().cend();
	auto tfI=DP_tmpvu8.begin()+actread.getLeftClipoff();
	auto tfE=DP_tmpvu8.end();

	for(uint32 readpos=actread.getLeftClipoff(); readpos<actread.getRightClipoff()-(basesperhash-1); ++bhsI, ++tfI, ++readpos){
	  if(!(bhsI->fwd.getFrequency()==1
	       || !bhsI->fwd.hasConfirmedFwdRev())){
	    auto ttfI=tfI+1;
	    uint32 tmpreadpos=readpos+1;
	    for(uint32 i=0; i<basesperhash-2 && ttfI!=tfE; ++i, ++ttfI, ++tmpreadpos){
	      if(tmpreadpos>=actread.getLeftClipoff() && tmpreadpos<actread.getRightClipoff()){
		CEBUG("cleaner " << tmpreadpos << endl);
		*ttfI=0;
	      }
	    }
	  }
	}
      }

      // and log changes
      bool mustoutputprefix=true;
      for(uint32 readpos=0; readpos<DP_tmpvu8.size(); ++readpos){
	if(DP_tmpvu8[readpos]){
	  auto start=readpos;
	  while(readpos<DP_tmpvu8.size() && DP_tmpvu8[readpos]) ++readpos;
	  if(readpos>start) --readpos;
	  if(mustoutputprefix){
	    DP_logfout << logprefix
		       << " rare kmer "
		       << actread.getName();
	    mustoutputprefix=false;
	  }
	  DP_logfout << "\t["
		     << start
		     << ".."
		     << readpos
		     << ']';
	}
      }
      if(!mustoutputprefix) DP_logfout << '\n';

      // masking with 'x'

      for(uint32 readpos=0; readpos<DP_tmpvu8.size(); ++readpos){
	if(DP_tmpvu8[readpos]) actread.changeBaseInSequence('x',0,readpos);
      }

      auto gaplen=basesperhash;
      maskClips_Read(actread, logprefix,gaplen,gaplen,gaplen);

      // TODO: first mask with X so that maskClips_Read() can do it's work, then use 'n' is not
      //   really elegant. Alternative would be yet another parameter to maskClips_Read()
      for(uint32 readpos=0; readpos<DP_tmpvu8.size(); ++readpos){
	if(DP_tmpvu8[readpos]) actread.changeBaseInSequence('n',0,readpos);
      }

      Read::setCoutType(Read::AS_TEXT);
      CEBUG("rare kmer result:" << actread << endl);
    }else{
      CEBUG("nomask " << actread.getName() << endl);
    }
  }
}
//#define CEBUG(bla)






/*************************************************************************
 *
 * expects reads to have baseflags set
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

struct tmpbhentry_t{
  uint32 from;
  uint32 to;
  uint8 freq;
};

void DataProcessing::buntifyReadsByHashFreq_Pool(ReadPool & rp, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::buntifyReadsByHashFreq()");

  cout << "Buntifying reads";
  if(rp.size()>500000) cout << " (this may take a while)";
  cout << " ... "; cout.flush();

  for(uint32 actid=0; actid<rp.size(); actid++){
    Read & actread=rp[actid];

    buntifyReadsByHashFreq_Read(rp[actid],basesperhash);
  }

  cout << "done." << endl;

  FUNCEND();

}
//#define CEBUG(bla)


void DataProcessing::buntifyReadsByHashFreq_Read(Read & actread, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::buntifyReadsByHashFreq_Read(Read & actread, uint32 basesperhash)");

  //Read::setCoutType(Read::AS_TEXT);
  //cout << actread;

  // remove old hash frequence tags
  //  but only if we will be able to set new HAF tags (read is long enough)
  if(actread.getLenClippedSeq()>=basesperhash){
    for(uint32 i=0; i<Read::REA_allhaftags.size(); ++i){
      actread.deleteTag(Read::REA_allhaftags[i]);
    }
  }

  if(actread.hasValidData()
     && actread.hasBaseHashStats()
     && actread.getLenClippedSeq()>=basesperhash){

    static multitag_t tmpmt("","","MIRA");

    DP_tmpvu8.clear();
    DP_tmpvu8.resize(actread.getLenSeq(),0);

    auto bhsI=actread.getBPosHashStats().cbegin();
    auto bhsE=actread.getBPosHashStats().cend();
    auto tfI=DP_tmpvu8.begin();
    auto tfE=DP_tmpvu8.end();

    priv_buntifyHelper(2, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(3, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(4, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(5, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(6, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(7, basesperhash, bhsI, bhsE, tfI, tfE);

    std::vector<tmpbhentry_t> telist;
    telist.reserve(20);

    //{
    //	cout << "bfr: " << actread.getName() << endl;
    //	for(uint32 i=0;i<DP_tmpvu8.size(); i++){
    //	  cout << "i: " << i << '\t' << static_cast<uint16>(DP_tmpvu8[i]) << endl;
    //	}
    //}

    uint32 from=0;
    uint32 to=0;
    for(; from<actread.getLenSeq(); from=to+1){
      to=from;
      uint8 actfreq=DP_tmpvu8[to];
      for(; to<actread.getLenSeq() && DP_tmpvu8[to]==actfreq; to++) {} ;
      to--;
      if(actfreq>0){
	telist.resize(telist.size()+1);
	telist.back().from=from;
	telist.back().to=to;
	telist.back().freq=actfreq;
      }
    }

    // for first or last entry, do not put tags for frequencies
    //  >=2 if their length is < basesperhash
    // BaCh 03.06.2011: hmmm, why not. OK, makes CAF/MAF bigger, but else?
    for(uint32 ti=0; ti<telist.size(); ti++){
      bool settag=true;
//	if(telist[ti].freq>=2 &&
//	   (ti==0 || ti==telist.size()-1)){
//	  if(telist[ti].to - telist[ti].from < basesperhash-1){
//	    settag=false;
//	  }
//	}
      if(settag) {
	tmpmt.identifier=Read::REA_allhaftags[telist[ti].freq];
	tmpmt.from=telist[ti].from;
	tmpmt.to=telist[ti].to;
	actread.addTagO(tmpmt);
      }
    }
  }

}
//#define CEBUG(bla)


void DataProcessing::priv_buntifyHelper(uint8 allowedfreq, uint32 basesperhash, std::vector<Read::bposhashstat_t>::const_iterator bhsI, std::vector<Read::bposhashstat_t>::const_iterator bhsE, std::vector<uint8>::iterator tfI, std::vector<uint8>::iterator tfE)
{
  for(; bhsI!= bhsE; bhsI++, tfI++){
    uint8 actfreq=bhsI->fwd.getFrequency();
    if(allowedfreq==actfreq){
      if(actfreq>0){
	auto ttfI=tfI;
	for(uint32 i=0; i<basesperhash && ttfI!=tfE; ++i, ++ttfI){
	  *ttfI=actfreq;
	}
      }
    }
  }
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/


void DataProcessing::addKMerForkTags_Pool(ReadPool & rp, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::addKMerForkTags_Pool(ReadPool & rp, uint32 basesperhash)");

  cout << "Adding fork tags";
  if(rp.size()>500000) cout << " (this may take a while)";
  cout << " ... "; cout.flush();

  static multitag_t tmpmt("","","MIRA");
  tmpmt.identifier=Read::REA_tagentry_idKMRF;

  for(uint32 actid=0; actid<rp.size(); actid++){
    Read & actread=rp[actid];

    //Read::setCoutType(Read::AS_TEXT);
    //cout << actread;

    // remove old KMRF tags
    actread.deleteTag(tmpmt.identifier);

    if(actread.hasValidData()
       && actread.hasBaseHashStats()){

      DP_tmpvu8.clear();
      DP_tmpvu8.resize(actread.getLenSeq(),0);

      auto bhsI=actread.getBPosHashStats().cbegin();
      auto bhsE=actread.getBPosHashStats().cend();
      auto tfI=DP_tmpvu8.begin();
      auto tfE=DP_tmpvu8.end();

      for(; bhsI!= bhsE; bhsI++, tfI++){
	if(bhsI->fwd.isKMerFork()){
	  auto ttfI=tfI;
	  for(uint32 i=0; i<basesperhash && ttfI!=tfE; ++i, ++ttfI){
	    *ttfI=1;
	  }
	}
	if(bhsI->rev.isKMerFork()){
	  auto ttfI=tfI;
	  for(uint32 i=0; i<basesperhash; ++i, --ttfI){
	    *ttfI=1;
	    if(ttfI!=DP_tmpvu8.begin()) break;
	  }
	}
      }

      uint32 from=0;
      uint32 to=0;
      for(; from<actread.getLenSeq(); from=to+1){
	to=from;
	if(DP_tmpvu8[to]){
	  for(; to<actread.getLenSeq() && DP_tmpvu8[to]; ++to) {} ;
	  to--;
	  tmpmt.from=from;
	  tmpmt.to=to;
	  actread.addTagO(tmpmt);
	}
      }
    }

    //if(actid==273250 || actid==273252){
    //  Read::setCoutType(Read::AS_TEXT);
    //  cout << actread;
    //}

  }

  cout << "done." << endl;

  FUNCEND();

}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/


void DataProcessing::performKMERRepeatTagging_Pool(ReadPool & rp, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::buntifyReadsByHashFreq()");

  cout << "Adding RMB tags by fork";
  if(rp.size()>500000) cout << " (this may take a while)";
  cout << " ... "; cout.flush();

  for(uint32 actid=0; actid<rp.size(); actid++){
    if(1){
      performKMERRepeatTagging_Read(rp[actid],basesperhash);
    }
  }

  cout << "done." << endl;

  FUNCEND();

}
//#define CEBUG(bla)

void DataProcessing::performKMERRepeatTagging_Read(Read & actread, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::addBla_Read(Read & actread, uint32 basesperhash)");
  if(actread.hasValidData()
     && actread.hasBaseHashStats()
     && actread.getLenSeq() >= 2*basesperhash){

    static multitag_t tmpmt("","addBla","MIRA");
    tmpmt.identifier=Read::REA_tagentry_idCRMr;

    DP_tmpvu8.clear();
    DP_tmpvu8.resize(actread.getLenSeq(),0);

    auto bhsE=actread.getBPosHashStats().cend();
    auto bhsIf=actread.getBPosHashStats().cbegin();
    uint32 tagpos=basesperhash-1;
    auto bhsIr=bhsIf+2*tagpos;
    auto tfI=DP_tmpvu8.begin()+tagpos;


    for(; bhsIr!= bhsE; ++bhsIf, ++bhsIr, ++tfI, ++tagpos){
      if(bhsIf->fwd.isKMerFork()
	 && bhsIr->rev.isKMerFork()){
	*tfI=1;
      }
    }

    uint32 runcount=0;
    for(tfI=DP_tmpvu8.begin(); tfI!=DP_tmpvu8.end(); ++tfI){
      if(*tfI){
	++runcount;
      }else{
	if(runcount==1){
	  tmpmt.from=tfI-1-DP_tmpvu8.begin();
	  tmpmt.to=tmpmt.from;
	  actread.addTagO(tmpmt);
	}
	runcount=0;
      }
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::proposedEndClipping_Pool(ReadPool & rp, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::proposedEndClipping_Pool(ReadPool & rp, uint32 basesperhash)");

#pragma omp parallel
  {
#pragma omp for
    for(uint32 rpi=0;rpi<rp.size();++rpi){
      Read & actread=rp[rpi];
      if(actread.hasValidData()
	 && actread.getLenClippedSeq()>=basesperhash // if read smaller it won't be cut anyway (in dubio pro reo) and having the check here is faster
	 && (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_proposeendclips
	 && actread.hasBaseHashStats()
	 && !(actread.isBackbone()
	      || actread.isRail())){
	proposedEndClipping_Read(actread, basesperhash);
      }
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::proposedEndClipping_Read(Read & actread, uint32 basesperhash)
{
  FUNCSTART("void DataProcessing::proposedEndClipping_Read(Read & actread, uint32 basesperhash)");

  pechelper_fwd(actread, basesperhash);
  pechelper_rev(actread, basesperhash);
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::pechelper_fwd(Read & actread, uint32 basesperhash)
{
  FUNCSTART("int32 DataProcessing::pechelper_fwd(Read & actread, uint32 basesperhash)")

  int32 lpos=actread.getLeftClipoff();
  auto bhsI=actread.getBPosHashStats().begin();
  advance(bhsI,lpos);
  for(int32 lend=static_cast<int32>(actread.getLenClippedSeq()); lpos<lend; ++lpos, ++bhsI) {
    CEBUG("lpos tst " << lpos << "\t" << *bhsI << endl);
    CEBUG("1: " << (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq
	  << "\t2: " << static_cast<uint16>(bhsI->fwd.getFrequency())
	  << "\t3: " << static_cast<uint16>(bhsI->rev.getFrequency()) << endl);
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq >0
       && (bhsI->fwd.getFrequency() > (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq
	   || bhsI->rev.getFrequency() > (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq)) {
      CEBUG("ffreq stop at " << lpos << "\n");
      break;
    }
    CEBUG("4: " << (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_ffr
	  << "\t5: " << static_cast<uint16>(bhsI->fwd.hasConfirmedFwdRev())
	  << "\t6: " << static_cast<uint16>(bhsI->rev.hasConfirmedFwdRev()) << endl);
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_ffr
       && ( bhsI->fwd.hasConfirmedFwdRev()
	    || bhsI->rev.hasConfirmedFwdRev())) {
      CEBUG("ffore stop at " << lpos << "\n");
      break;
    }
    CEBUG("7: " << (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_fcmst
	  << "\t8: " << static_cast<uint16>(bhsI->fwd.hasConfirmedMultipleSeqType())
	  << "\t9: " << static_cast<uint16>(bhsI->rev.hasConfirmedMultipleSeqType()) << endl);
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_fcmst
       && ( bhsI->fwd.hasConfirmedMultipleSeqType()
	    || bhsI->rev.hasConfirmedMultipleSeqType())) {
      CEBUG("fcmst stop at " << lpos << "\n");
      break;
    }
    CEBUG("A: " << (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_fsalp
	  << "\tB: " << static_cast<uint16>(bhsI->fwd.hasSeenAtLowPos())
	  << "\tC: " << static_cast<uint16>(bhsI->rev.hasSeenAtLowPos()) << endl);
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_fsalp
       && ( bhsI->fwd.hasSeenAtLowPos()
	    || bhsI->rev.hasSeenAtLowPos())) {
      CEBUG("fsalp stop at " << lpos << "\n");
      break;
    }
  }

  // trying to save things
  // on high GC, some contigs end due to lousy good coverage in both directions,
  // bit there would be enough coverage in one direction
  // So, we try to go back the clipping as long as there's
  //  - no singlet
  //  - no kmer fork
  // in that direction, back to the original lclip at most

  CEBUG("l+ " << lpos << "\trco " << actread.getRightClipoff() << endl);
  if(lpos != actread.getLeftClipoff()){

    CEBUG("trying to save\n");
    int32 xpos=lpos+basesperhash-1;
    if(xpos>=actread.getRightClipoff()){
      xpos=actread.getRightClipoff()-1;
      if(xpos<0) xpos=0;
      lpos=xpos-(basesperhash-1);
      if(lpos<0) lpos=0;
    }
    bhsI=actread.getBPosHashStats().begin();
    advance(bhsI,xpos);

    CEBUG("xpos: " << xpos << "\tlpos now: " << lpos << endl);
    for(;lpos > actread.getLeftClipoff() && bhsI->rev.isValid() && bhsI->rev.getFrequency() < 2 && !bhsI->rev.isKMerFork(); --lpos, --bhsI){}
    CEBUG("lpos after jump: " << lpos << endl);

    bool advanced=false;
    for(;lpos > actread.getLeftClipoff(); --lpos, --bhsI){
      CEBUG((bhsI->rev.getFrequency() < 2) << " " << bhsI->rev.isKMerFork());
      if(bhsI->rev.getFrequency() < 2
	 || bhsI->rev.isKMerFork()){
	if(advanced) ++lpos; // go back one if we went too far
	CEBUG(" nonfork lext stop at " << lpos << "\n");
	break;
      }
      CEBUG(" back\n");
      advanced=true;
    }
    if(lpos<actread.getLeftClipoff()) lpos=actread.getLeftClipoff();
    CEBUG("saved to " << lpos << endl);
  }

  if(lpos != actread.getLeftClipoff()){
    if(lpos>0 && lpos>actread.getLenSeq()) lpos=actread.getLenSeq();
    CEBUG("pcb l: " << actread.getName() << " " << actread.getLeftClipoff()
	  << " " << lpos << endl);
    if(lpos==actread.getLenSeq()){
      actread.setRQClipoff(actread.getLeftClipoff());
    }else{
      actread.setLQClipoff(lpos);
    }
  }

  return;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::pechelper_rev(Read & actread, uint32 basesperhash)
{
  FUNCSTART("int32 DataProcessing::pechelper_rev(Read & actread, uint32 basesperhash)");

  int32 rpos=actread.getRightClipoff();
  auto bhsI=actread.getBPosHashStats().cbegin();
  advance(bhsI,rpos);

  for(int32 rend=actread.getLeftClipoff(); rpos >rend; --rpos){
    --bhsI;
    CEBUG("rpos tst " << rpos << "\t" << *bhsI << endl);
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_bfreq
       && (bhsI->fwd.getFrequency() > (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_bfreq
	   || bhsI->rev.getFrequency() > (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_bfreq)) {
      CEBUG("bfreq stop at " << rpos << "\n");
      break;
    }
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_bfr
       && ( bhsI->fwd.hasConfirmedFwdRev()
	    || bhsI->rev.hasConfirmedFwdRev())) {
      CEBUG("bfore stop at " << rpos << "\n");
      break;
    }
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_bcmst
       && ( bhsI->fwd.hasConfirmedMultipleSeqType()
	    || bhsI->rev.hasConfirmedMultipleSeqType())) {
      CEBUG("bcmst stop at " << rpos << "\n");
      break;
    }
    if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_pec_bsalp
       && ( bhsI->fwd.hasSeenAtLowPos()
	    || bhsI->rev.hasSeenAtLowPos())) {
      CEBUG("fsalp stop at " << rpos << "\n");
      break;
    }
  }

  // trying to save things
  // on high GC, some contigs end due to lousy good coverage in both directions,
  // bit there would be enough coverage in one direction
  // So, we try to go back the clipping as long as there's
  //  - no singlet
  //  - no kmer fork
  // in that direction, back to the original rclip at most
  if(rpos != actread.getRightClipoff()){
    CEBUG("trying to save\n");
    int32 xpos=rpos-(basesperhash-1);
    if(xpos<actread.getLeftClipoff()){
      xpos=actread.getLeftClipoff();
      rpos=xpos+(basesperhash-1);
      if(rpos>actread.getRightClipoff()) rpos=actread.getRightClipoff();
    }

    bhsI=actread.getBPosHashStats().begin();
    advance(bhsI,xpos);

    CEBUG("rpos now: " << rpos << endl);
    for(;rpos < actread.getRightClipoff() && bhsI->fwd.isValid() && bhsI->fwd.getFrequency() < 2 && !bhsI->fwd.isKMerFork(); ++rpos, ++bhsI){}
    CEBUG("after jump: " << rpos << endl);

    for(;rpos < actread.getRightClipoff(); ++rpos, ++bhsI){
      CEBUG("rext " << rpos << "\n");
      if(bhsI->fwd.getFrequency() < 2
	 || bhsI->fwd.isKMerFork()){
	CEBUG("nonfork rext stop at " << rpos << "\n");
	break;
      }
    }
    if(rpos>actread.getRightClipoff()) rpos=actread.getRightClipoff();
    CEBUG("saved to " << rpos << endl);
  }

  if(rpos != actread.getRightClipoff()){
    CEBUG("pcb r: " << actread.getName() << " " << actread.getRightClipoff()
	  << " " << rpos << endl);

    BUGIFTHROW(rpos>actread.getRightClipoff(),"rpos>actread.getRightClipoff() ???");

    actread.setRQClipoff(rpos);

    // special handling of Solexa GGC.G error
    // from point of right clip, 15 bases backwards:
    //  search for first ggc.g and clip there
    if(actread.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA
       && (*DP_miraparams_ptr)[0].getAssemblyParams().as_clip_pec_sxaggcxg
       && actread.getLenClippedSeq() >=15){
      //Read::setCoutType(Read::AS_TEXTSHORT);
      //cout << r;
      std::string searchstr=actread.getSeqAsChar();
      boost::to_lower(searchstr);
      int64 searchstart=actread.getRightClipoff()-15;
      if(searchstart<0) searchstart=0;
      size_t found;
      do{
	found=searchstr.find(DP_ggcstring, searchstart);
	if (found!=std::string::npos){
	  searchstart=found+1;
	  if(found < actread.getRightClipoff()
	     && found+4<actread.getRightClipoff()
	     && searchstr[found+4]=='g'){
	    actread.setRQClipoff(static_cast<int32>(found+4));
	    found=std::string::npos; // stop the loop
	  }
	}
      }while(found!=std::string::npos);
    }
  }

  return;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::clipBadSolexaEnds_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipBadSolexaEnds_Pool(ReadPool & rp, const std::string & logprefix)");

  for(uint32 i=0;i<rp.size();++i){
    Read & actread=rp[i];
    if(actread.hasValidData()
       && actread.isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
       && !(actread.isBackbone()
	    || actread.isRail())){
      clipBadSolexaEnds_Read(actread,logprefix);
    }
  }
}


/*************************************************************************
 *
 * TODO: not really good, rethink that for eukaryotes
 *
 *************************************************************************/

void DataProcessing::clipBadSolexaEnds_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipBadSolexaEnds_Read(Read & actread, const std::string & logprefix)");

  // invalidate all Solexa reads that have a stretch of 20 A (or T)
  //  or if has stretch >= 12 and non-A (or T) bases < 20%
  // N in-between do not reset the counter
  // invalidate by setting left seq vec to length of read

  int32 runindex=actread.getLeftClipoff();
  char actbase=' ';
  uint32 bcount=0;

  uint32 arun=0;
  uint32 maxarun=0;
  uint32 nona=0;

  uint32 trun=0;
  uint32 maxtrun=0;
  uint32 nont=0;

  for(; runindex<actread.getRightClipoff(); runindex++) {
    actbase=static_cast<char>(toupper(actread.getBaseInSequence(runindex)));
    if(actbase!='N'){
      bcount++;
      if(actbase=='A'){
	arun++;
	if(arun>maxarun) maxarun=arun;
	nont++;
	trun=0;
      }else if(actbase=='T'){
	trun++;
	if(trun>maxtrun) maxtrun=trun;
	nona++;
	arun=0;
      }else{
	nona++;
	nont++;
	arun=0;
	trun=0;
      }
    }
  }
  if(maxarun>=20){
    actread.setLSClipoff(actread.getLenSeq());
    DP_logfout << logprefix << " bad solexa end: A hard "
	       << actread.getName()
	       << '\n';
  }else if(maxarun>=12){
    uint32 ratio= static_cast<uint32>((static_cast<double>(100.0)/bcount)*nona);
    if(ratio<20) {
      actread.setLSClipoff(actread.getLenSeq());
      DP_logfout << logprefix << " bad solexa end: A soft "
		 << actread.getName()
		 << '\n';
    }
  }

  if(maxtrun>=20){
    actread.setLSClipoff(actread.getLenSeq());
    DP_logfout << logprefix << " bad solexa end: T (hard) "
	       << actread.getName()
	       << '\n';
  }else if(maxtrun>=12){
    uint32 ratio= static_cast<uint32>((static_cast<double>(100.0)/bcount)*nont);
    if(ratio<20) {
      actread.setLSClipoff(actread.getLenSeq());
      DP_logfout << logprefix << " bad solexa end: T (soft) "
		 << actread.getName()
		 << '\n';
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 * clip all lowercase at the end of reads
 *
 *************************************************************************/

void DataProcessing::lowerCaseClipping_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::lowerCaseClipping_Pool(ReadPool & rp, const std::string & logprefix)");

  uint64 totallen=0;
  uint64 lowercaselen=0;
  for(uint32 i=0;i<rp.size();i++){
    Read & actread=rp.getRead(i);
    if(actread.hasValidData()
       && ((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_front
	    || (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_back)
       && !(actread.isBackbone()
	    || actread.isRail())){
      totallen+=actread.getLenClippedSeq();
      int32 runindex=actread.getLeftClipoff();
      for(; runindex<actread.getRightClipoff(); ++runindex){
	if(islower(actread.getBaseInSequence(runindex))) lowercaselen++;
      }
    }
  }

  if(totallen==lowercaselen) {
    cout << "Lowercase clip: all sequences to be clipped are lowercase?! Failsafe: no clipping performed.\n";
    return;
  }

  for(uint32 i=0;i<rp.size();i++){
    Read & actread=rp.getRead(i);
    if(actread.hasValidData()
       && !(actread.isBackbone()
	    || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_front){
	lowerCaseClippingFront_Read(actread,logprefix);
      }
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_back){
	lowerCaseClippingBack_Read(actread,logprefix);
      }
    }
  }

  FUNCEND();
}

/*************************************************************************
 *
 * clip all lowercase at the end of reads
 *
 *************************************************************************/

void DataProcessing::lowerCaseClippingFront_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::lowerCaseClippingFront_Read(Read & actread, const std::string & logprefix)");

  // TODO: for streaming, implement check counter by rgid for all sequence lowercase

  int32 runindex=actread.getLeftClipoff();
  for(; runindex<actread.getRightClipoff(); ++runindex){
    char ab=actread.getBaseInSequence(runindex);
    if(!islower(ab)
       && ab != 'N'
       && ab != 'X') break;
  }
  // TODO: 01.01.2013 check this, changed != to >
  if(runindex>actread.getLeftClipoff()) {
    actread.setLSClipoff(runindex);
    DP_logfout << logprefix << " changed left (lowercase) "
	       << actread.getName() << " to " << actread.getLeftClipoff() << '\n';
  }
}


void DataProcessing::lowerCaseClippingBack_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::lowerCaseClippingBack_Read(Read & actread, const std::string & logprefix)");

  // TODO: for streaming, implement check counter by rgid for all sequence lowercase

  int32 runindex=actread.getRightClipoff()-1;
  for(; runindex>=actread.getLeftClipoff() && islower(actread.getBaseInSequence(runindex)); --runindex) ;
  // TODO: implement jumping over N,X like for front
  // TODO: 01.01.2013 really check this, changed != to <
  if(runindex<actread.getRightClipoff()-1) {
    actread.setRSClipoff(runindex+1);
    DP_logfout << logprefix << " changed right (lowercase) "
	       << actread.getName() <<  " to " << actread.getRightClipoff() << '\n';
  }
  //cout << actread;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::qualClips_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::qualClips_Pool(ReadPool & rp, const std::string & logprefix)");

  cout << "Starting qual clips: ";

  for(uint32 i=0;i<rp.size();i++){
    Read & r=rp.getRead(i);
    if(r.hasValidData()
       && !(r.isBackbone()
	    || r.isRail())){
      if((*DP_miraparams_ptr)[r.getSequencingType()].getAssemblyParams().as_clip_quality) {
	qualClips_Read(r,logprefix);
      }
    }
  }
}


void DataProcessing::qualClips_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::qualClips_Read(Read & actread, const std::string & logprefix)");

  int32 oldlq=actread.getLQClipoff();
  int32 oldrq=actread.getRQClipoff();
  actread.performQualityClip(
    (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_minqual,
    (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_winlen);
  int changed=2;
  if(oldlq>actread.getLQClipoff()){
    // old LQ was more conservative, change back
    actread.setLQClipoff(oldlq);
    --changed;
  }
  if(oldrq<actread.getRQClipoff()){
    // old RQ was more conservative, change back
    actread.setRQClipoff(oldrq);
    --changed;
  }
  if(changed){
    DP_logfout << logprefix
	       << " changed qual. "
	       << actread.getName()
	       << "\tfrom: " << oldlq << ' ' << oldrq << "\tto: "
	       << actread.getLQClipoff()
	       << ' '
	       << actread.getRQClipoff() << '\n';
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::maskClips_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::maskClips_Pool(ReadPool & rp, const std::string & logprefix)");

  cout << "Starting qual clips: ";

  for(uint32 i=0;i<rp.size();i++){
    Read & r=rp.getRead(i);
    if(r.hasValidData()
       && !(r.isBackbone()
	    || r.isRail())){
      if((*DP_miraparams_ptr)[r.getSequencingType()].getAssemblyParams().as_clip_maskedbases) {
	maskClips_Read(r,logprefix);
      }
    }
  }
}


void DataProcessing::maskClips_Read(Read & actread, const std::string & logprefix, int32 gapsize, int32 maxfrontgap, int32 maxendgap)
{
  FUNCSTART("void DataProcessing::maskClips_Read(Read & actread, const std::string & logprefix, int32 gapsize, int32 maxfrontgap, int32 maxendgap)");

  if(gapsize<0) gapsize=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_maskedbase_gapsize;
  if(maxfrontgap<0) maxfrontgap=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_maskedbase_maxfrontgap;
  if(maxendgap<0) maxendgap=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_maskedbase_maxendgap;

  auto oldlm=actread.getLMClipoff();
  auto oldrm=actread.getRMClipoff();
  auto oldlc=actread.getLeftClipoff();
  auto oldrc=actread.getRightClipoff();
  actread.setClipoffsToMaskedChars(
    gapsize,
    maxfrontgap,
    maxendgap,
    false);
  actread.setClipoffsToMaskedChars(
    1,
    1,
    1,
    true);
  int changed=2;
  if(oldlm>=actread.getLMClipoff()){
    // old LM was more conservative, change back
    if(oldlm>actread.getLMClipoff()) actread.setLMClipoff(oldlm);
    --changed;
  }
  if(oldrm<=actread.getRMClipoff()){
    // old RM was more conservative, change back
    if(oldrm<actread.getRMClipoff()) actread.setRMClipoff(oldrm);
    --changed;
  }
  if(changed){
    DP_logfout << logprefix
	       << " changed mask. "
	       << actread.getName()
	       << "\tfrom: " << oldlm << ' ' << oldrm << " (" << oldlc << " " << oldrc << ")\tto: "
	       << actread.getLMClipoff()
	       << ' '
	       << actread.getRMClipoff() << " (" << actread.getLeftClipoff() << " " << actread.getRightClipoff() << ")\n";
  }

  FUNCEND();
  return;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::minimumQualityThreshold_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::minimumQualityThreshold_Pool(ReadPool & rp, const std::string & logprefix)");

  cout << "Starting minimum quality threshold clip ... "; cout.flush();

  uint32 numkilled=0;

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if(!minimumQualityThreshold_Read(actread,logprefix)){
	++numkilled;
      }
    }
  }
  cout << "done. Killed " << numkilled << " reads.\n";

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool DataProcessing::minimumQualityThreshold_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("bool DataProcessing::minimumQualityThreshold_Read(Read & actread, const std::string & logprefix)");

  base_quality_t minqual=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_minthreshold;
  uint32 minnum=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_numminthreshold;
  bool mustkill=true;
  for(const auto & qe : actread.getQualities()){
    if(qe>=minqual
       && --minnum==0) {
      mustkill=false;
      break;
    }
  }

  if(mustkill){
    actread.setLQClipoff(actread.getLenSeq());
    actread.setRQClipoff(actread.getLenSeq());
    DP_logfout << logprefix << ' ' << actread.getName() << ": min qual threshold not met, killed\n";
  }

  FUNCEND();
  return !mustkill;
}






/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::minimumLeftClip_Pool(ReadPool & rp, const std::string & logprefix, bool qual, bool seqvec, bool mask)
{
  for(uint32 ri=0;ri<rp.size();++ri){
    Read & actread=rp[ri];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_ensureminimumleftclipoff){
	minimumLeftClip_Read(actread,logprefix,qual,seqvec,mask);
      }
    }
  }
}

void DataProcessing::minimumLeftClip_Read(Read & actread, const std::string & logprefix, bool qual, bool seqvec, bool mask)
{
  auto oldlc=actread.getLeftClipoff();
  if(oldlc < (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minslrequired){
    if(qual) actread.setLQClipoff((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minqlsetto);
    if(seqvec) actread.setLSClipoff((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minqlsetto);
    if(mask) actread.setLMClipoff((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minqlsetto);

    DP_logfout << logprefix
	       << " changed minleft. "
	       << actread.getName()
	       << "\tLeft: "
	       << oldlc
	       << "\t -> "
	       << actread.getLeftClipoff()
	       << '\n';
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::minimumRightClip_Pool(ReadPool & rp, const std::string & logprefix, bool qual, bool seqvec, bool mask)
{
  for(uint32 ri=0;ri<rp.size();++ri){
    Read & actread=rp[ri];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_ensureminimumrightclipoff){
	minimumRightClip_Read(actread,logprefix,qual,seqvec,mask);
      }
    }
  }
}

void DataProcessing::minimumRightClip_Read(Read & actread, const std::string & logprefix, bool qual, bool seqvec, bool mask)
{
  auto oldrc=actread.getRightClipoff();
  if(actread.getLenSeq()-oldrc < (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minsrrequired){
    int32 newr=actread.getLenSeq()-(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minsrrequired;
    if(qual) actread.setRQClipoff(newr);
    if(seqvec) actread.setRSClipoff(newr);
    if(mask) actread.setRMClipoff(newr);

    DP_logfout << logprefix
	       << " changed minRight. "
	       << actread.getName()
	       << "\tRightt: "
	       << oldrc
	       << "\t -> "
	       << actread.getRightClipoff()
	       << '\n';
  }
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::badSequenceSearch_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::badSequenceSearch_Pool(ReadPool & rp, const std::string & logprefix)");

  cout << "Performing search for bad sequence quality ... "; cout.flush();

  for(uint32 ri=0;ri<rp.size();++ri){
    Read & actread=rp[ri];
    if(actread.hasValidData()
       && actread.hasQuality()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_badstretchquality){
      }
    }
  }
}

void DataProcessing::badSequenceSearch_Read(Read & actread, const std::string & logprefix)
{
  uint32 winlen=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_badstretchquality_winlen;
  base_quality_t minqual=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_badstretchquality_minqual;

  const std::vector<base_quality_t> & bquals=actread.getQualities();
  int32 runi=actread.getLeftClipoff();
  int32 endi=actread.getRightClipoff();

  uint32 qualsbelow=0;
  bool foundbad=false;
  for(; runi < endi; runi++){
    if(bquals[runi] < minqual){
      ++qualsbelow;
      if(qualsbelow >= winlen){
	foundbad=true;
	break;
      }
    }else{
      qualsbelow=0;
    }
  }

  if(foundbad) {
    int32 newrclip=runi-qualsbelow+1;
    int32 shortened=actread.getRightClipoff()-newrclip;
    //cout << actread.getName() << " has bad stretch, shortening by " << actread.getRightClipoff()-newrclip << '\n';
    if(newrclip < actread.getLQClipoff()) newrclip=actread.getLQClipoff();
    actread.setRQClipoff(newrclip);
    DP_logfout << logprefix << " bad seq. "
	       << actread.getName()
	       << "\tShortened by " << shortened
	       << "\tNew right: "
	       << actread.getRQClipoff()
	       << '\n';
  }

  FUNCEND();
  return;
}





/*************************************************************************
 *
 * clip poly-A in forward and poly-T in reverse direction
 * or: clip only after the poly-stretches, and tag the stretches with Fpas
 *
 * If poly stretches are kept, they are also "cleaned", i.e., bases are
 *  forced to be A (or T) in the detected stretch
 *
 *************************************************************************/

void DataProcessing::clipPolyATAtEnds_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipPolyATAtEnds_Pool(ReadPool & rp, const std::string & logprefix)");

  cout << "Clipping or tagging poly A/T stretches at ends of reads ... ";
  cout.flush();

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_polyat){
	clipPolyATAtEnds_Read(actread,logprefix);
      }
    }
  }
}

void DataProcessing::clipPolyATAtEnds_Read(Read & actread, const std::string & logprefix)
{
  CEBUG(actread.getName() << endl);
  CEBUG(actread << endl);

  auto & as_params=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams();
  uint32 mincount=as_params.as_clip_polyat_len;
  uint32 maxbad=as_params.as_clip_polyat_maxerrors;
  int32 grace=static_cast<int32>(as_params.as_clip_polyat_maxgap);
  bool keepstretch=as_params.as_clip_polyat_keeppolystretch;

  // search poly-a in forwarddirection
  {
    int32 lpolystart=-1;
    int32 rpolyend=-1;
    if(searchPolyBaseFrom5Prime(actread,'a',lpolystart,rpolyend,mincount,maxbad,grace)){
      if(keepstretch){
	actread.setRMClipoff(rpolyend+1);
	CEBUG("setting rm " << rpolyend+1 << endl);
	for(auto i=lpolystart; i<=rpolyend; ++i){
	  if(toupper(actread.getBaseInSequence(i))!='A'){
	    actread.changeBaseInSequence('a',0,i);
	  }
	}
      }else{
	actread.setRMClipoff(lpolystart);
	CEBUG("setting rm " << lpolystart << endl);
      }
      CEBUG("taggingl " << lpolystart << " " << rpolyend << endl);

      DP_logfout << logprefix << " poly-A fwd. "
		 << actread.getName()
		 << "\tMask right: "
		 << actread.getRMClipoff()
		 << '\n';

      DP_tmpmtpolyAT.from=lpolystart;
      DP_tmpmtpolyAT.to=rpolyend;
      actread.addTagO(DP_tmpmtpolyAT);
    }
  }

  // search poly-t in reverse direction
  {
    int32 lpolystart=-1;
    int32 rpolyend=-1;

    if(searchPolyBaseFrom3Prime(actread,'t',lpolystart,rpolyend,mincount,maxbad,grace)){
      if(keepstretch){
	actread.setLMClipoff(lpolystart);
	CEBUG("setting lm " << lpolystart << endl);
	for(auto i=lpolystart; i<=rpolyend; ++i){
	  if(toupper(actread.getBaseInSequence(i))!='T'){
	    actread.changeBaseInSequence('t',0,i);
	  }
	}
      }else{
	actread.setLMClipoff(rpolyend+1);
	CEBUG("setting lm " << rpolyend+1 << endl);
      }
      CEBUG("taggingl " << lpolystart << " " << rpolyend << endl);

      DP_logfout << logprefix << " poly-T rev. "
		 << actread.getName()
		 << "\tMask left: "
		 << actread.getLMClipoff()
		 << '\n';

      DP_tmpmtpolyAT.from=lpolystart;
      DP_tmpmtpolyAT.to=rpolyend;
      actread.addTagO(DP_tmpmtpolyAT);
    }
  }
}




/*************************************************************************
 *
 * Search poly-base (mincount length and maximum maxbad other bases) from left
 *  side of read (with 'grace' length grace if not encountered),
 *
 *  return:
 *    - true if found and return left and right coordinates in lpolystart and
 *      rpolyend
 *    - false if not found (lpolystart and rpolyend undefined)
 *
 *************************************************************************/

bool DataProcessing::searchPolyBaseFrom5Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)
{
  FUNCSTART("bool DataProcessing::searchPolyBaseFrom5Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)");


  BUGIFTHROW(!dptools::isValidACGTBase(polybase),"Ummm ... " << polybase << " is not ACGT?");
  BUGIFTHROW(grace<0,"grace (" << grace << ") < 0 ?");
  BUGIFTHROW(maxbad>=mincount,"maxbad (" << maxbad << ") >= mincount (" << mincount << ") ?");

  CEBUG(actread.getName() << endl);
  CEBUG(actread << endl);

  lpolystart=-1;
  rpolyend=-1;

  int32 runindex=actread.getLeftClipoff();
  int32 lastgoodrunindex=runindex;
  char actbase=' ';
  bool found=false;

  for(; grace >=0 && runindex<actread.getRightClipoff(); ++runindex, --grace) {
    actbase=actread.getBaseInSequence(runindex);
    CEBUG("###1 : " << grace << " " << runindex << "\t" << actbase << endl);
    if(dptools::areBasesContained(polybase,actbase)){
      lpolystart=runindex;
      lastgoodrunindex=runindex;
      uint32 acount=0;
      uint32 othercount=0;
      char cbase;
      for(; lastgoodrunindex<actread.getRightClipoff(); lastgoodrunindex++){
	cbase=actread.getBaseInSequence(lastgoodrunindex);
	if(dptools::areBasesContained(polybase,cbase)){
	  acount++;
	}else if(tolower(cbase)!='n'){
	  othercount++;
	  if(othercount>maxbad) break;
	}
      }
      if(acount>=mincount) {
	found=true;
	// get off non-poly characters as far as possible
	if(lastgoodrunindex==actread.getRightClipoff()) lastgoodrunindex--;
	while(lastgoodrunindex>runindex && !dptools::areBasesContained(polybase,actread.getBaseInSequence(lastgoodrunindex))) lastgoodrunindex--;
	rpolyend=lastgoodrunindex;
	break;
      }
      lpolystart=-1;
    }else{
      lpolystart=-1;
    }
  }

  if(rpolyend >=0 && lpolystart != -1) {
    FUNCEND();
    return true;
  }

  FUNCEND();
  return false;
}



/*************************************************************************
 *
 * Search poly-base (mincount length and maximum maxbad other bases) from left
 *  side of read (with 'grace' length grace if not encountered),
 *
 *  return:
 *    - true if found and return left and right coordinates in lpolystart and
 *      rpolyend
 *    - false if not found (lpolystart and rpolyend undefined)
 *
 * TODO: how dumb to have an own function instead of a templated function
 *       shared with 5p and working on iterators
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool DataProcessing::searchPolyBaseFrom3Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)
{
  FUNCSTART("bool DataProcessing::searchPolyBaseFrom3Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)");


  BUGIFTHROW(!dptools::isValidACGTBase(polybase),"Ummm ... " << polybase << " is not ACGT?");
  BUGIFTHROW(grace<0,"grace (" << grace << ") < 0 ?");
  BUGIFTHROW(maxbad>=mincount,"maxbad (" << maxbad << ") >= mincount (" << mincount << ") ?");

  CEBUG(actread.getName() << endl);
  CEBUG(actread << endl);

  lpolystart=-1;
  rpolyend=-1;


  int32 runindex=actread.getRightClipoff()-1;
  int32 lastgoodrunindex=runindex;
  char actbase=' ';
  bool found=false;

  CEBUG("Reverse " << actread.getName() << '\n');

  for(; grace >=0 && runindex>=actread.getLeftClipoff(); --runindex, --grace) {
    actbase=static_cast<char>(tolower(actread.getBaseInSequence(runindex)));
    CEBUG("###1 : " << grace << " " << runindex << "\t" << actbase << endl);
    if(dptools::hasNucleicAcidInIUPAC(polybase,actbase)){
      rpolyend=runindex;
      lastgoodrunindex=runindex;
      uint32 tcount=0;
      uint32 othercount=0;
      uint32 runcount=0;
      char cbase;
      char dbase;
      for(; lastgoodrunindex>=actread.getLeftClipoff(); --lastgoodrunindex, ++runcount){
	cbase=actread.getBaseInSequence(lastgoodrunindex);
	CEBUG("###2 : " << runcount << " " << lastgoodrunindex << "\t" << cbase << " " << tcount << " " << othercount << endl);
	if(dptools::areBasesContained(polybase,cbase)){
	  tcount++;
	  if(othercount>0 && runcount>=mincount){
	    dbase=actread.getBaseInSequence(lastgoodrunindex+mincount);
	    if(dptools::areBasesContained(polybase,dbase)){
	      --othercount;
	    }
	  }
	}else if(tolower(cbase)!='n'){
	  othercount++;
	  if(othercount>maxbad) break;
	}
      }
      if(tcount>=mincount) {
	CEBUG("Found tcount\n");
	found=true;
	// get off non-t characters as far as possible
	if(lastgoodrunindex<actread.getLeftClipoff()) lastgoodrunindex++;
	while(lastgoodrunindex<rpolyend && !dptools::areBasesContained(polybase,actread.getBaseInSequence(lastgoodrunindex))) lastgoodrunindex++;
	lpolystart=lastgoodrunindex;
	break;
      }
      rpolyend=-1;
    }else{
      rpolyend=-1;
    }
  }

  CEBUG("LPOLYSTART: " << lpolystart << "\tRPOLYEND: " << rpolyend << endl);

  if(lpolystart >=0 && rpolyend != -1) {
    FUNCEND();
    return true;
  }

  FUNCEND();
  return false;
}
//#define CEBUG(bla)








/*************************************************************************
 *
 * clip poly-base at right end of read
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::clipPolyBaseAtEnd_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipPolyBaseAtEnd_Pool(ReadPool & rpool, const std::string & logprefix)");

  cout << "Clipping dubious poly-base stretches at end of reads ... ";
  cout.flush();

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_3ppolybase_len){
	clipPolyBaseAtEnd_Read(actread,logprefix);
      }
    }
  }
}


void DataProcessing::clipPolyBaseAtEnd_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipPolyBaseAtEnd_Read(Read & actread, const std::string & logprefix)");

  assembly_parameters const & as_params= (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams();
  CEBUG(actread.getName() << endl);

  uint32 mincount=as_params.as_clip_3ppolybase_len;
  if(mincount==0){
    MIRANOTIFY(Notify::FATAL, "-AS:c3ppmsl may not be 0");
  }
  if(actread.getLenClippedSeq() < mincount) return;

  uint32 maxbad=as_params.as_clip_3ppolybase_maxerrors;
  int32 grace=static_cast<int32>(as_params.as_clip_3ppolybase_maxgap);

  // first guess which base might be a polybase
  //
  // count occurrences of bases in last (mincount+grace or mincount?) positions of read
  // the largest count >=30% of real bases (no 'n') wins

  Read::setCoutType(Read::AS_FASTA);
  CEBUG(actread << endl);
  Read::setCoutType(Read::AS_TEXTSHORT);
  CEBUG(actread << endl);

  DP_tmpau32_128['a']=0;
  DP_tmpau32_128['c']=0;
  DP_tmpau32_128['g']=0;
  DP_tmpau32_128['n']=0;
  DP_tmpau32_128['t']=0;

  int32 runindex=actread.getRightClipoff()-1;
  for(uint32 ri=0; ri<mincount && runindex>=actread.getLeftClipoff(); --runindex, ++ri) {
    char actbase=static_cast<char>(tolower(actread.getBaseInSequence(runindex)));
    ++DP_tmpau32_128[actbase];
  }

  CEBUG("CV: " << DP_tmpau32_128['a'] << " " << DP_tmpau32_128['c'] << " " << DP_tmpau32_128['g'] << " " << DP_tmpau32_128['t'] << endl);

  uint32 realbases=DP_tmpau32_128['a']+DP_tmpau32_128['c']+DP_tmpau32_128['g']+DP_tmpau32_128['t'];
  uint32 maxreal=std::max(DP_tmpau32_128['a'],std::max(DP_tmpau32_128['c'],std::max(DP_tmpau32_128['g'],DP_tmpau32_128['t'])));

  CEBUG("RB: " << realbases << "\tMR: " << maxreal << endl);

  char tentativepolybase='?';
  if(realbases>0 && 100*maxreal/realbases >= 30){
    CEBUG("MRThresh\n");
    for(uint32 testi=0; testi<4; ++testi){
      if(DP_tmpau32_128["acgt"[testi]]==maxreal){
	CEBUG("MRThreshHit\n");
	tentativepolybase="acgt"[testi];
	break;
      }
    }
  }

  // so, if a tentative polybase was found, try to find some clips and clip if found

  if(tentativepolybase!='?') {
    int32 lpolystart=-1;
    int32 rpolyend=-1;
    CEBUG("looking...\n");
    if(searchPolyBaseFrom3Prime(actread,tentativepolybase,lpolystart,rpolyend,mincount,maxbad,grace)){
      actread.setRMClipoff(lpolystart);
      CEBUG("setting rm " << lpolystart << endl);

      DP_logfout << logprefix << " poly-base " << tentativepolybase << " at end "
		 << actread.getName()
		 << "\tMask right: "
		 << actread.getRMClipoff()
		 << '\n';
    }
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::adaptorRightClip_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::adaptorRightClip_Pool(ReadPool & rpool, const std::string & logprefix)");

  cout << "Searching for sequencing adaptors.\n";
  cout.flush();

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      adaptorRightClip_Read(actread,logprefix);
    }
  }
}

void DataProcessing::adaptorRightClip_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::adaptorRightClip_Read(Read & actread, const std::string & logprefix)");

  priv_EnsureAdapRegexes(actread.getReadGroupID());
  BUGIFTHROW(actread.getReadGroupID().getLibId()>=DP_adapres.size(),"Huh? no re lib " << actread.getReadGroupID().getLibId());

  priv_EnsureAdapSkims(actread.getReadGroupID());

  //assembly_parameters const & as_params= (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams();
  auto oldrsclip=actread.getRSClipoff();
  auto newclip=-1;
  if(DP_adapskims[actread.getReadGroupID().getLibId()].skimptr==nullptr){
    cout << "Bah? nullptr for " << actread.getReadGroupID().getLibId() << "?" << endl;
    return;
  }
  readid_t ridadapfound=-1;
  if(DP_adapskims[actread.getReadGroupID().getLibId()].skimptr){
    newclip=DP_adapskims[actread.getReadGroupID().getLibId()].skimptr->findAdaptorRightClip(actread,9,ridadapfound,DP_threadid);
  }
  if(newclip>=0){
    ++DP_stats.cadapright;
    actread.setRSClipoff(newclip);
    DP_logfout << logprefix << " "
	       << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
	       << " adaptor: ";
    if(ridadapfound>=0) {
      DP_logfout << DP_adapskims[actread.getReadGroupID().getLibId()].poolptr->getRead(ridadapfound).getName() << " in ";
    }
    DP_logfout << actread.getName()
	       << " changed right clip from " << oldrsclip << " to " << newclip << "\n";
  }else{
    std::string seq(actread.getSeqAsChar());
    boost::to_upper(seq);

    std::match_results<std::string::const_iterator> what;
    auto flags = std::regex_constants::match_default;
    decltype(seq.cbegin()) start, end;

    for(auto & msre : DP_adapres[actread.getReadGroupID().getLibId()].adapres){
      bool dosearch=true;
      if(msre.hasmaster){
	if(!regex_search(start, end, what, msre.masterre, flags)) {
	  dosearch=false;
	}
      }
      bool breakit=false;
      if(dosearch){
	for(auto & thisre : msre.slaveres){
	  start = seq.begin();
	  end = seq.end();
	  if(regex_search(start, end, what, thisre, flags)) {
	    if(what.position()< oldrsclip){
	      ++DP_stats.cadaprightpartial;
	      actread.setRSClipoff(what.position());
	      DP_logfout << logprefix << " "
			 << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
			 << " partial end adaptor: " << actread.getName()
			 << " changed right clip from " << oldrsclip << " to " << what.position() << "\n";
	      breakit=true;
	      break;
	    }
	  }
	}
      }
      if(breakit) break;
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::seqMatchPhiX174_Read(Read & actread, const std::string & logprefix, bool filter)
{
  FUNCSTART("void DataProcessing::seqMatchPhiX174_Read(Read & actread, const std::string & logprefix)");

  BUGIFTHROW(!DP_px174hs_init,"Phi X 174 search structure not initialised.");

  //auto numbaithits=DP_phix174hashstatistics.checkBaitHit(actread,DP_baiting_singlereadvhraparray,DP_baiting_tagmaskvector,false);
  auto numbaithits=DP_phix174hashstatistics.checkBaitHit(actread,false,0);
  if(numbaithits>10){
    ++DP_stats.cphix174;
    DP_logfout << logprefix << " "
	       << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
	       << " phix174 in "
	       << actread.getName();
    if(filter) {
      actread.setRSClipoff(0);
      DP_logfout << " ... killed read\n";
    }else{
      DP_logfout << " ... matched\n";
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::seqMatchRRNA_Read(Read & actread, const std::string & logprefix, uint32 numkmers)
{
  FUNCSTART("void DataProcessing::seqMatchRRNA_Read(Read & actread, const std::string & logprefix)");

  BUGIFTHROW(!DP_rrnahs_init,"rRNA search structure not initialised.");

  uint32 numbaithits=DP_rrnahashstatistics.checkBaitHit(actread,false,0);
  if(numbaithits>0 && numbaithits>=numkmers){
    ++DP_stats.crrna;
    DP_logfout << logprefix << " "
	       << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
	       << " rRNA in "
	       << actread.getName()
	       << " ... killed read\n";
    actread.setRSClipoff(0);
  }
}


/*************************************************************************
 *
 * careful, rrna pair not done here!!!
 * then again: maybe retire singlethread callability from outside.
 *
 *************************************************************************/

void DataProcessing::stdTreatmentPool_SingleThread(std::vector<MIRAParameters> & mp, DataProcessing & dp, ReadPool & rpool, std::vector<uint8> * debrisreasonptr, std::string & logprefix, bool progress, int32 fromid, int32 toid)
{
  FUNCSTART("void DataProcessing::stdTreatmentPool_SingleThread(std::vector<MIRAParameters> & mp, DataProcessing & dp, ReadPool & rpool, std::string & logprefix, bool progress, int32 fromid, int32 toid)");

  if(fromid<0) fromid=0;
  if(toid<0) toid=rpool.size();
  BUGIFTHROW(fromid>toid,"fromid>toid ?");
  BUGIFTHROW(toid>rpool.size(),"toid>rpool.size()?");

  auto & agp = mp[0].getAssemblyParams();

  std::unique_ptr<ProgressIndicator<int64> > pi;
  if(progress) pi=std::unique_ptr<ProgressIndicator<int64>>(new ProgressIndicator<int64>(fromid,toid));
  for(uint32 actid=fromid; actid < toid; ++actid){
    if(progress) pi->increaseprogress();
    Read & actread = rpool[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_SHORTONLOAD);

      auto & asp = mp[actread.getSequencingType()].getAssemblyParams();

      if(asp.as_search_phix174){
	dp.seqMatchPhiX174_Read(actread,logprefix,asp.as_filter_phix174);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_PHIX174);
      }
      if(agp.as_filter_rrna){
	dp.seqMatchRRNA_Read(actread,logprefix,agp.as_filter_rrna_numkmers);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_RRNA);
      }
      if(asp.as_clip_knownadaptorsright){
	dp.adaptorRightClip_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_KNOWNADAPTORRIGHT);
      }
      // do the clip for bad solexa ends after adaptor clip
      // reason: in quite a number of illumina reads, one has
      //     xxxxxxxx - ADAPTOR - AAAAAAAAAAA...AAAAAAAAAAA - junk
      // removing the adaptor first helps to keep the "xxxxxxxx" part
      //  as the bad solexa ends kill the reads completely :-(
      if(asp.as_clip_badsolexaends && actread.isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	dp.clipBadSolexaEnds_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_BADSOLEXAEND);
      }
      if(asp.as_clip_lowercase_front){
	dp.lowerCaseClippingFront_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_LOWERCASEFRONT);
      }
      if(asp.as_clip_lowercase_back){
	dp.lowerCaseClippingBack_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_LOWERCASEBACK);
      }
      if(asp.as_clip_quality_minthreshold){
	dp.minimumQualityThreshold_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_QUALMINTHRESHOLD);
      }
      if(asp.as_clip_quality){
	dp.qualClips_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_QUALCLIPS);
      }
      if(asp.as_clip_maskedbases){
	dp.maskClips_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MASKEDBASES);
      }
      bool mlc=asp.as_clip_ensureminimumleftclipoff;
      if(asp.as_clip_badstretchquality){
	if(mlc){
	  dp.maskClips_Read(actread,logprefix);
	  mlc=false;
	  priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MASKEDBASES);
	}
	dp.badSequenceSearch_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_BADSEQUENCESERACH);
      }
      if(asp.as_clip_3ppolybase){
	dp.clipPolyBaseAtEnd_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_POLYBASEATEND);
      }
      if(asp.as_clip_polyat){
	dp.clipPolyATAtEnds_Read(actread,logprefix);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_POLYAT);
      }
      if(mlc){
	dp.minimumLeftClip_Read(actread,logprefix,true,false,false);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MINLEFTCLIP);
      }
      if(asp.as_clip_ensureminimumrightclipoff){
	dp.minimumRightClip_Read(actread,logprefix,false,true,false);
	priv_stp_helperDebris(mp,rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MINRIGHTCLIP);
      }
    }
  }
  if(progress) pi->finishAtOnce();
}

void DataProcessing::priv_stp_helperDebris(std::vector<MIRAParameters> & mp, ReadPool & rpool, Read & actread, int32 rid, std::vector<uint8> * debrisreasonptr, uint8 reason)
{
  FUNCSTART("void DataProcessing::priv_stp_helperDebris(ReadPool & rpool, Read & actread, int32 rid, std::vector<uint8> * debrisreasonptr, uint8 reason)");

  if(debrisreasonptr == nullptr) return;
  if(debrisreasonptr->empty()) return;
  BUGIFTHROW(rid>=debrisreasonptr->size(),"rid " << rid << " >= debrisreasonptr->size() " << debrisreasonptr->size());
  if(actread.getLenClippedSeq() < mp[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength){
    if((*debrisreasonptr)[rid]==0) (*debrisreasonptr)[rid]=reason;
  }
}


/*************************************************************************
 *
 * rRNA pair in addition to single thread, so both are different atm!
 *
 *************************************************************************/

void DataProcessing::stdTreatmentPool_MultiThread(std::vector<MIRAParameters> & mp, DataProcessing & dpcollector, std::vector<std::unique_ptr<DataProcessing>> & dpv, ReadPool & rpool, std::vector<uint8> * debrisreasonptr, std::string & logprefix, bool progress, int32 fromid, int32 toid)
{
  FUNCSTART("void DataProcessing::stdTreatmentPool_MultiThread(std::vector<DataProcessing> & dpv, ReadPool & rpool, std::string & logprefix, bool progress)");

  if(fromid<0) fromid=0;
  if(toid<0) toid=rpool.size();
  BUGIFTHROW(dpv.empty(),"dpv.empty() ?");
  BUGIFTHROW(fromid>toid,"fromid>toid ?");
  BUGIFTHROW(toid>rpool.size(),"toid>rpool.size()?");

  dpv[0]->priv_EnsurePhiX174Statistics();
  if(mp[0].getAssemblyParams().as_filter_rrna){
    dpv[0]->priv_EnsureRRNAStatistics();
  }

  // make sure data for adaptor clipping is ready.
  for(auto rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(mp[rgid.getSequencingType()].getAssemblyParams().as_clip_knownadaptorsright){
      dpv[0]->priv_EnsureAdapRegexes(rgid);
      dpv[0]->priv_EnsureAdapSkims(rgid);
    }
  }

  threadsharecontrol_t tsc;

  tsc.from=fromid;
  tsc.to=toid;
  tsc.todo=fromid;
  tsc.done=fromid;
  tsc.stepping=1000;

  uint32 numthreads=dpv.size();
  boost::thread_group workerthreads;
  for(uint32 ti=0; ti<numthreads;++ti){
    dpv[ti]->setThreadID(ti);
    workerthreads.create_thread(boost::bind(&DataProcessing::priv_stdTreatmentThread, ti, &tsc, &mp, &(*dpv[ti]), &rpool, debrisreasonptr, &logprefix));
  }

  ProgressIndicator<int64> pi(fromid,toid);
  while(tsc.done!=toid){
    if(progress) pi.progress(tsc.done);
    sleep(1);
  }
  if(progress) pi.finishAtOnce(cout);

  // they normally should all have exited at this point, but be nice and play by the rules
  workerthreads.join_all();

  // collect all stats
  for(auto & dpvp : dpv){
    dpcollector.DP_stats.cphix174+=dpvp->DP_stats.cphix174;
    dpcollector.DP_stats.crrna+=dpvp->DP_stats.crrna;
    dpcollector.DP_stats.cadapright+=dpvp->DP_stats.cadapright;
    dpcollector.DP_stats.cadaprightpartial+=dpvp->DP_stats.cadaprightpartial;
  }

  // do this here.
  // maybe push to single thread, but that might lead to conflicts
  auto & as_fixparams= mp[0].getAssemblyParams();
  if(as_fixparams.as_filter_rrna_pairs
     && debrisreasonptr != nullptr
     && !debrisreasonptr->empty()) {
    for(uint32 rid=fromid; rid < toid; ++rid){
      auto tprid=rpool[rid].getTemplatePartnerID();
      if(tprid>=0
	 && (*debrisreasonptr)[rid] == Assembly::DEBRIS_CLIP_RRNA
	 && (*debrisreasonptr)[tprid] == 0) {
	++dpcollector.DP_stats.crrna;
	Read & actread = rpool[tprid];
	dpv[0]->DP_logfout << logprefix << " "
			   << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
			   << " rRNA by pair for "
			   << actread.getName()
			   << " ... killed read\n";
	actread.setRSClipoff(0);
	priv_stp_helperDebris(mp,rpool,actread,tprid,debrisreasonptr,Assembly::DEBRIS_CLIP_RRNA_PAIR);
      }
    }
  }
}


void DataProcessing::priv_stdTreatmentThread(uint32 threadnum, threadsharecontrol_t * tscptr, std::vector<MIRAParameters> * mpptr, DataProcessing * dpptr, ReadPool * rpoolptr, std::vector<uint8> * debrisreasonptr, std::string * logprefixptr)
{
  FUNCSTART("void DataProcessing::priv_stdTreatmentThread(uint32 threadnum, threadsharecontrol_t * tscptr, DataProcessing * dpptr, ReadPool * rpoolptr, std::string * logprefixptr)");

  try{
    int32 from;
    int32 to;
    while(true){
      {
	boost::mutex::scoped_lock lock(tscptr->accessmutex);
	if(tscptr->todo >= tscptr->to) break;
	from=tscptr->todo;
	tscptr->todo+=tscptr->stepping;
	if(tscptr->todo > tscptr->to) tscptr->todo = tscptr->to;
	to=tscptr->todo;
      }
      stdTreatmentPool_SingleThread(*mpptr, *dpptr,*rpoolptr,debrisreasonptr,*logprefixptr,false,from,to);
      {
	boost::mutex::scoped_lock lock(tscptr->accessmutex);
	tscptr->done+=tscptr->stepping;
	if(tscptr->done > tscptr->to) tscptr->done=tscptr->to;
      }
    }
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
}


// explicit template instantiations needed for the linker to create these for library files
template void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics<vhash64_t> & hsd, std::vector<uint8> * debrisreasonptr);
template uint32 DataProcessing::performSDBGChimeraSearch_Pool(ReadPool & rp, HashStatistics<vhash64_t> & hsd, uint32 trimfreq, std::vector<uint8> * debrisreasonptr, const std::string & logprefix);
template bool DataProcessing::performSDBGChimeraSearch_Read(Read & actread, HashStatistics<vhash64_t> & hsd, const std::string & logprefix);
template uint32 DataProcessing::performSDBGEdits_Pool(ReadPool & rp, HashStatistics<vhash64_t> & hsd, uint32 trimfreq);
#ifndef KMER_INTERNALTYPE
template void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics<vhash128_t> & hsd, std::vector<uint8> * debrisreasonptr);
template void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics<vhash256_t> & hsd, std::vector<uint8> * debrisreasonptr);
template void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics<vhash512_t> & hsd, std::vector<uint8> * debrisreasonptr);

template uint32 DataProcessing::performSDBGChimeraSearch_Pool(ReadPool & rp, HashStatistics<vhash128_t> & hsd, uint32 trimfreq, std::vector<uint8> * debrisreasonptr, const std::string & logprefix);
template uint32 DataProcessing::performSDBGChimeraSearch_Pool(ReadPool & rp, HashStatistics<vhash256_t> & hsd, uint32 trimfreq, std::vector<uint8> * debrisreasonptr, const std::string & logprefix);
template uint32 DataProcessing::performSDBGChimeraSearch_Pool(ReadPool & rp, HashStatistics<vhash512_t> & hsd, uint32 trimfreq, std::vector<uint8> * debrisreasonptr, const std::string & logprefix);

template bool DataProcessing::performSDBGChimeraSearch_Read(Read & actread, HashStatistics<vhash128_t> & hsd, const std::string & logprefix);
template bool DataProcessing::performSDBGChimeraSearch_Read(Read & actread, HashStatistics<vhash256_t> & hsd, const std::string & logprefix);
template bool DataProcessing::performSDBGChimeraSearch_Read(Read & actread, HashStatistics<vhash512_t> & hsd, const std::string & logprefix);

template uint32 DataProcessing::performSDBGEdits_Pool(ReadPool & rp, HashStatistics<vhash128_t> & hsd, uint32 trimfreq);
template uint32 DataProcessing::performSDBGEdits_Pool(ReadPool & rp, HashStatistics<vhash256_t> & hsd, uint32 trimfreq);
template uint32 DataProcessing::performSDBGEdits_Pool(ReadPool & rp, HashStatistics<vhash512_t> & hsd, uint32 trimfreq);
#endif
