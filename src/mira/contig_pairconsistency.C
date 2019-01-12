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


#include "contig.H"
#include "simple_2Dsignalprocessing.H"

using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)   {cout << bla; cout.flush();}
std::pair<int32,int32> Contig::findBestPairConsistencyRange()
{
  FUNCSTART("void Contig::checkPairConsistency()");

#define MAXINSIZETOCHECK 80000

  std::pair<int32,int32> retval(0,getContigLength());

  //return retval;

  if(getContigLength()==0) return retval;

  std::vector<ReadGroupLib::ReadGroupID> rgtocheck;
  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(rgid.hasTemplateInfo()
       && rgid.getInsizeFrom()>=0
       && rgid.getInsizeTo()>=0){
      rgtocheck.push_back(rgid);
    }
  }

  if(rgtocheck.size()){
    CEBUG("Need to check readgroups:\n");
    std::vector<uint8> breakmarker(getContigLength(),0);
    std::vector<std::vector<uint64> > spanner(rgtocheck.size());
    std::vector<uint64> tmpforstats;
    std::vector<std::vector<int32> > valleydepths(rgtocheck.size());
    std::vector<std::vector<int32> > peakheights(rgtocheck.size());

    bool hascuts=false;

    // create vectors with coverage data
    auto spanner12=spanner;
    auto nonspanner=spanner;
    for(auto rgi=0; rgi<spanner.size(); ++rgi){
      priv_cprForRGID(rgtocheck[rgi],spanner[rgi],spanner12[rgi],nonspanner[rgi]);
      //std::string desc=getContigName();
      //desc+='_';
      //desc+=rgtocheck[rgi].getGroupName();
      //dbgContainerToWiggle(spanner[rgi],getContigName(),desc+"_sp");
      //dbgContainerToWiggle(spanner12[rgi],getContigName(),desc+"_sp12");
      //dbgContainerToWiggle(nonspanner[rgi],getContigName(),desc+"_nsp");
    }

    // calculate peaks
    for(auto rgi=0; rgi<spanner.size(); ++rgi){
      // do this only for libraries with small template size
      if(rgtocheck[rgi].getInsizeTo()>=0 && rgtocheck[rgi].getInsizeTo()<=MAXINSIZETOCHECK){
	// We won't search here, just need some numbers for the
	// correct threshold in non-spanner down below
	tmpforstats=spanner12[rgi];

	coverageinfo_t tci;
	calcStatsOnContainer(tci,tmpforstats);
	cout << "spanner12 - normal:\n" << tci << endl;
	calcSecondOrderStatsOnContainer(tci,tmpforstats);
	cout << "spanner12 - 2n:\n" << tci << endl;

	uint64 nsthreshold=tci.mean/4;

	tmpforstats=nonspanner[rgi];
	calcStatsOnContainer(tci,tmpforstats);
	cout << "nonspanner - normal:\n" << tci << endl;
	calcSecondOrderStatsOnContainer(tci,tmpforstats);
	cout << "nonspanner - 2n:\n" << tci << endl;

	bool inpeak=false;
	auto vsI=nonspanner[rgi].cbegin();
	uint64 pos=0;
	uint64 threshold=nsthreshold;
	if(tci.mean*10>threshold) threshold=tci.mean*10;

	uint64 firsttpos=1;
	uint64 lasttpos=getContigLength()-1;

	cout << "threshold=" << threshold << endl;
	cout << "ftp: " << firsttpos << endl;
	cout << "ltp: " << lasttpos << endl;
	for(auto sI=nonspanner[rgi].cbegin(); sI!=nonspanner[rgi].cend(); ++sI, ++pos){
	  bool peakstop=false;
	  if(pos<firsttpos) continue;
	  if(pos>=lasttpos){
	    if(inpeak) {
	      peakstop=true;
	    }else{
	      break;
	    }
	  }
	  if(inpeak){
	    if(*sI<threshold){
	      peakstop=true;
	      inpeak=false;
	    }
	  }else{
	    if(*sI>threshold){
	      inpeak=true;
	      vsI=sI;
	    }
	  }
	  if(peakstop){
	    auto veI=sI;
	    cout << "Potential peak: " << vsI-nonspanner[rgi].cbegin() << "\t" << veI-nonspanner[rgi].cbegin() << "\t" << veI-vsI << endl;
	    if(veI-nonspanner[rgi].cbegin() <= rgtocheck[rgi].getInsizeTo()){
	      cout << "Near front, skip.\n";
	    }else if(nonspanner[rgi].cend()-vsI <= rgtocheck[rgi].getInsizeTo()){
	      cout << "Near end, skip.\n";
	    }else{
	      auto maxval=*vsI;
	      for(auto tI=vsI; tI<veI; ++tI){
		if(*tI>maxval) maxval=*tI;
	      }
	      for(auto tI=vsI; tI<veI; ++tI){
		if(*tI==maxval) {
		  peakheights[rgi].push_back(tI-nonspanner[rgi].cbegin());
		  cout << "New peak max " << maxval << " at " << tI-nonspanner[rgi].cbegin() << endl;
		  hascuts=true;
		}
	      }
	    }
	  }
	}
      }
    }

    // calculate valley info
    for(auto rgi=0; rgi<spanner.size(); ++rgi){
      coverageinfo_t tci;
      tmpforstats=spanner[rgi];
      calcStatsOnContainer(tci,tmpforstats);
      cout << "spanner - normal:\n" << tci << endl;
      calcSecondOrderStatsOnContainer(tci,tmpforstats);
      cout << "spanner - 2n:\n" << tci << endl;

      // spanner are searched for valleys in coverage
      // valley is <= mean cov/3, size >=50
      if(tci.mean>=9){
	bool invalley=false;
	auto vsI=spanner[rgi].cbegin();
	uint64 pos=0;
	uint64 threshold=tci.mean/3;

	uint64 firsttpos=getContigLength();
	uint64 lasttpos=0;
	{
	  auto dist=distToFirstThreshold(spanner[rgi].begin(),spanner[rgi].end(),tci.mean/2);
	  if(dist>=0) firsttpos=dist;
	  dist=distToFirstThreshold(spanner[rgi].rbegin(),spanner[rgi].rend(),tci.mean/2);
	  if(dist>=0) lasttpos=getContigLength()-dist;
	}

	cout << "threshold=" << threshold << endl;
	cout << "ftp: " << firsttpos << endl;
	cout << "ltp: " << lasttpos << endl;
	for(auto sI=spanner[rgi].cbegin(); sI!=spanner[rgi].cend(); ++sI, ++pos){
	  bool valleystop=false;
	  if(pos<firsttpos) continue;
	  if(pos>=lasttpos){
	    if(invalley) {
	      valleystop=true;
	    }else{
	      break;
	    }
	  }
	  if(invalley){
	    if(*sI>threshold){
	      valleystop=true;
	      invalley=false;
	    }
	  }else{
	    if(*sI<threshold){
	      invalley=true;
	      vsI=sI;
	    }
	  }
	  if(valleystop){
	    auto veI=sI;
	    int32 vstart=vsI-spanner[rgi].cbegin();
	    int32 vstop=veI-spanner[rgi].cbegin();
	    cout << "Potential valley: " << vstart << "\t" << vstop << "\t" << vstop-vstart << endl;
	    uint32 confirmedbypeak=0;
	    int32 peakvalleydist=500;
	    for(auto & phe : peakheights[rgi]){
	      if((vstart-phe >= 0 && vstart-phe < peakvalleydist)
		 || (phe >= vstart && phe <= vstop)
		 || (phe-vstop>=0 && phe-vstop < peakvalleydist)){
		++confirmedbypeak;
		cout << "confirmedbypeak " << phe << endl;
		phe=-1;
	      }
	    }
	    if(confirmedbypeak){
	      // temove all -1 entries in this area as the valley will give better cutting places
	      // yeah, remove_if with a lambda, but let's wait a bit more so that newer C++
	      // compilers get a broader base
	      for(auto phI=peakheights[rgi].begin(); phI!=peakheights[rgi].end();){
		if(*phI==-1) {
		  phI=peakheights[rgi].erase(phI);
		}else{
		  ++phI;
		}
	      }

	      auto minval=*vsI;
	      for(auto tI=vsI; tI<veI; ++tI){
		if(*tI<minval) minval=*tI;
	      }
	      for(auto tI=vsI; tI<veI; ++tI){
		if(*tI==minval) {
		  valleydepths[rgi].push_back(tI-spanner[rgi].cbegin());
		  hascuts=true;
		  cout << "New valley min " << minval << " at " << tI-spanner[rgi].cbegin() << endl;
		}
	      }
	    }
	  }
	}
      }
    }

    if(hascuts){
      // now create the break ranges +/-50 from the valleys
      for(auto rge : valleydepths){
	for(auto pos : rge){
	  auto from=pos-50;
	  if(from<0) from=0;
	  auto to=pos+50;
	  if(to>getContigLength()) to=getContigLength();
	  for(; from<to; ++from) breakmarker[from]=1;
	}
      }

      // add in the peak markers
      for(auto rge : peakheights){
	for(auto pos : rge){
	  breakmarker[pos]=1;
	}
      }

      // last step: find longest range without problems
      retval.first=0;
      retval.second=1;
      auto tpair=retval;
      bool inrange=!breakmarker[0];
      if(inrange) tpair.first=0;
      for(int32 bmi=0; bmi<breakmarker.size(); ++bmi){
	//cout << bmi << "\t" << inrange << " " << static_cast<int16>(breakmarker[bmi]);
	if(inrange){
	  if(breakmarker[bmi]) {
	    //cout << " stop range";
	    tpair.second=bmi;
	    if(tpair.second-tpair.first > retval.second-retval.first){
	      retval=tpair;
	      //cout << " best";
	    }
	    inrange=false;
	  }
	}else{
	  if(!breakmarker[bmi]) {
	    //cout << " starting range";
	    inrange=true;
	    tpair.first=bmi;
	  }
	}
	//cout << "\n";
      }
      if(inrange){
	tpair.second=getContigLength();
	if(tpair.second-tpair.first > retval.second-retval.first){
	  retval=tpair;
	}
      }
    }
  }

  CEBUG("returning ret: " << retval.first << "\t" << retval.second << endl);
  return retval;
}
#define CEBUG(bla)

#define CEBUG(bla)   {cout << bla; cout.flush();}
void Contig::priv_cprForRGID(ReadGroupLib::ReadGroupID rgid, std::vector<uint64> & spanner, std::vector<uint64> & spanner12, std::vector<uint64> & nonspanner)
{
  FUNCSTART("void Contig::checkPairConsistency()");

  spanner.clear();
  spanner.resize(getContigLength(),0);
  spanner12=spanner;
  nonspanner=spanner;

  CEBUG("cprForRGID on readgroup:\n" << rgid);
  auto pcrI=CON_reads.begin();
  auto opcrI=pcrI;
  auto crE=CON_reads.end();
  for(; pcrI != crE; ++pcrI){
    if(pcrI->getReadGroupID()!=rgid
       || !pcrI->hasTemplateInfo()
       || pcrI->isBackbone()
       || pcrI->isRail()
       || pcrI->isCoverageEquivalentRead()) continue;

    auto tpid=pcrI->getTemplatePartnerID();
    opcrI=CON_reads.getIteratorOfReadpoolID(tpid);
    bool tpartnerincontig=opcrI!=crE;
    //CEBUG("Looking at " << pcrI.getORPID() << "\t" << pcrI->getName() << endl);
    //CEBUG(pcrI->isBackbone() << '\t' << pcrI->isRail() << '\t' << pcrI->isCoverageEquivalentRead() << '\t' << pcrI->hasTemplateInfo()  << '\t' << pcrI->getTemplatePartnerID() << '\t' << tpartnerincontig << endl);

    if(tpartnerincontig){
      auto start1=pcrI.getReadStartOffset();
      auto end1=start1+pcrI->getLenClippedSeq();
      auto start2=opcrI.getReadStartOffset();
      auto end2=start2+opcrI->getLenClippedSeq();

      if(start1<=start2){
	int64 addval=2;
	if(start1==start2) addval=1;

	auto start=std::min(start1,start2);
	auto end=std::max(end1,end2);

	for(auto * ptr=&spanner[start]; ptr!=&spanner[end]; ++ptr) *ptr+=addval;
	for(auto * ptr=&spanner12[start1]; ptr!=&spanner12[end1]; ++ptr) *ptr+=addval;
	for(auto * ptr=&spanner12[start2]; ptr!=&spanner12[end2]; ++ptr) *ptr+=addval;
      }
    }else{
      // template partner not in contig
      //
      // but maybe the template partner was completely clipped away?
      // then this read would not be a non-spanner, just a paired read
      //  which became a single read. How to treat this?
      // TODO: cleanest solution: after loading / clipping, move
      //  widowed reads into own readgroup
      //
      // for the time being, just make a small workaround: if clipped
      //  length of partner <25, do not treat it as non-spanner

      bool isnonspanner=true;

      if(tpid>=0){
	if(!CON_readpool->getRead(tpid).hasValidData()
	   || CON_readpool->getRead(tpid).getLenClippedSeq()<=25){
	  isnonspanner=false;
	}
      }

      if(isnonspanner){
	auto start=pcrI.getReadStartOffset();
	auto end=start+pcrI->getLenClippedSeq();
	for(auto * ptr=&nonspanner[start]; ptr!=&nonspanner[end]; ++ptr) *ptr+=1;
      }
    }

  }

  // hacky hack
  // templates having both reads in a contig should be counted only once and the
  //  "if(start1<=start2)" takes care of most cases ... but not when templates start
  //  at the same position.
  // therefore, in the loop above the standard add is +2 for templates not starting
  //  at same position and +1 at same position (but these will run twice)
  // aftererwards (well, here), we divide by 2 and get back what we want
  for(auto & se : spanner) se/=2;
  for(auto & se : spanner12) se/=2;
}
#define CEBUG(bla)
