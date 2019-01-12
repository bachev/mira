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
 *
 */


#include <iostream>
#include <sstream>

#include "mira/sam_collect.H"
#include "mira/maf_parse.H"

#include "errorhandling/errorhandling.H"
#include "util/progressindic.H"

#include "util/stlimprove.H"

// boost::split
#include <boost/algorithm/string.hpp>

#define CEBUG(bla)


using std::cout;
using std::cerr;
using std::endl;


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void SAMCollect::processMAF(const std::string & mafname)
{
  FUNCSTART("void SAMCollect::processMAF(const std::string & mafname)");

  collectInfoFromMAF(mafname);
  processSAMRIs();

  FUNCEND();
}


void SAMCollect::processSAMRIs()
{
  FUNCSTART("void SAMCollect::processSAMRIs(const std::string & mafname)");

  cout << "Sorting read info ... "; cout.flush();
  if(SAMC_samris.size()<=65535){
    mstd::ssort(SAMC_samris,SAMCollect::samrinfo_t::lt_templateid);
  }else{
    mstd::psort(SAMC_samris,SAMCollect::samrinfo_t::lt_templateid);
  }
  cout << "done" << endl;

  // this has changed the ordering in SAMC_samris, therefore the samri_ids in
  //  SAMC_rname2samriid are not valid anymore. Generate the new ones and
  //  put them into SAMC_rname2samriid
  {
    std::vector<std::pair<size_t,size_t>> neworder(SAMC_samris.size());
    auto noI=neworder.begin();
    size_t counter=0;
    for(auto sI=SAMC_samris.begin(); sI!=SAMC_samris.end(); ++sI, ++noI, ++counter){
      *noI=std::pair<size_t,size_t>(sI->initial_samri_id,counter);
    }
    if(neworder.size()<=65535){
      mstd::ssort(neworder);
    }else{
      mstd::psort(neworder);
    }
    for(auto & rnme : SAMC_rname2samriid){
      rnme.second=neworder[rnme.second].second;
    }
  }

  // now go through the SAMC_samris template by template to set flags
  auto sI=SAMC_samris.begin();
  auto sE=sI;
  for(; sI!=SAMC_samris.end(); sI=sE){
    BUGIFTHROW(sI->rgid.isDefaultNonValidReadGroupID(),"Template " << SAMC_templatenames[sI->templateid] << " has no readgroup set? Should not be.");
    for(; sE!=SAMC_samris.end() && sE->templateid==sI->templateid; ++sE) {};

    CEBUG("LAT: " << SAMC_templatenames[sI->templateid] << endl);
    CEBUG(sI->rgid << endl);
    if(sE-sI > 1){
      uint32 samflags=1;
      bool samecontig=true;
      uint32 numnonzerosegments=0;
      uint32 numzerosegments=0;
      int32 tlen=0;
      int32 seg1left=true;
      bool eachsegmentproperlyaligned=true;
      bool seenfirstsegment=false;
      bool seenlastsegment=false;
      for(auto iI=sI; iI!=sE; ++iI) {
	if(iI->contigid!=sI->contigid){
	  samecontig=false;
	}
	if(iI->template_segment){
	  ++numnonzerosegments;
	  if(iI->template_segment==1){
	    seenfirstsegment=true;
	  }else if(iI->template_segment==255){
	    seenlastsegment=true;
	  }
	}else{
	  ++numzerosegments;
	}
      };
      if(numzerosegments && numnonzerosegments){
	if(sI->templateid>0) {
	  cout << "WARNING: template " << SAMC_templatenames[sI->templateid] << " (" << sI->templateid << ") has " << sE-sI << " segments, but some of them have no segment order???\n";
	  eachsegmentproperlyaligned=false;
	}
      }else if(!samecontig || numzerosegments || !seenfirstsegment || !seenlastsegment){
	eachsegmentproperlyaligned=false;
	CEBUG("BAH 1" << "\tsc " << samecontig << "\tzs " << numzerosegments << "\tfs " << seenfirstsegment << "\tls " << seenlastsegment << "\tpc " << static_cast<int16>(sI->rgid.getSegmentPlacementCode()) << "\n");
      }

      tlen=0;
      if(eachsegmentproperlyaligned){
	BUGIFTHROW(sI->template_segment==0,"sI->template_segment==0 ???");
	switch(sI->rgid.getSegmentPlacementCode()){
	case ReadGroupLib::SPLACE_UNKNOWN : {
	  // intentionally do nothing
	  // we don't know the placement code, so assembled reads from a
	  //  template should be seen as correct
	  break;
	}
	case ReadGroupLib::SPLACE_SF : {
	  auto initialstartpos=sI->clippedstartpos;
	  auto dirtaken=sI->dir;
	  for(auto iI=sI; iI!=sE; ++iI){
	    if(iI->dir != dirtaken
	       || (dirtaken>0 && iI->clippedstartpos < initialstartpos)
	       || (dirtaken<0 && iI->clippedstartpos > initialstartpos)){
	      eachsegmentproperlyaligned=false;
	      CEBUG("BAH 2sf\n");
	      break;
	    }
	  }
	  break;
	}
	case ReadGroupLib::SPLACE_SB : {
	  auto initialstartpos=sI->clippedstartpos;
	  auto dirtaken=sI->dir;
	  for(auto iI=sI; iI!=sE; ++iI){
	    if(iI->dir != dirtaken
	       || (dirtaken>0 && iI->clippedstartpos > initialstartpos)
	       || (dirtaken<0 && iI->clippedstartpos < initialstartpos)){
	      eachsegmentproperlyaligned=false;
	      CEBUG("BAH 2sb\n");
	      break;
	    }
	  }
	  break;
	}
	case ReadGroupLib::SPLACE_SU : {
	  auto initialstartpos=sI->clippedstartpos;
	  auto dirtaken=sI->dir;
	  for(auto iI=sI; iI!=sE; ++iI){
	    if(iI->dir != dirtaken){
	      eachsegmentproperlyaligned=false;
	      CEBUG("BAH 2su\n");
	      break;
	    }
	  }
	  break;
	}
	case ReadGroupLib::SPLACE_RF :
	case ReadGroupLib::SPLACE_FR : {
	  if(sE-sI>2){
	    eachsegmentproperlyaligned=false;
	    cout << "WARNING: template " << SAMC_templatenames[sI->templateid] << " has " << sE-sI << " segments, but should have only 2???\n";
	  }else{
	    // for overlapping reads, it makes not much sense to check RF/FR orientation
	    bool isoverlapping=true;
	    if(sI->clippedstartpos+sI->clippedlen < (sI+1)->clippedstartpos
	       || (sI+1)->clippedstartpos+(sI+1)->clippedlen < sI->clippedstartpos){
	      isoverlapping=false;
	    }
	    if(!isoverlapping){
	      // check orientation
	      bool isok=true;
	      if(sI->dir>0 && sI->clippedstartpos > (sI+1)->clippedstartpos) isok=false;
	      if(sI->dir<0 && sI->clippedstartpos < (sI+1)->clippedstartpos) isok=false;
	      if(sI->rgid.getSegmentPlacementCode()==ReadGroupLib::SPLACE_RF) isok=!isok;
	      eachsegmentproperlyaligned=isok;
	      if(!eachsegmentproperlyaligned) {CEBUG("BAH 3\n");}
	    }
	  }
	  break;
	}
	default : {
	  MIRANOTIFY(Notify::FATAL,"Placement code " << static_cast<int16>(sI->rgid.getSegmentPlacementCode()) << " is unknown???");
	}
	}

	if(eachsegmentproperlyaligned){
	  // check distance
	  int32 mins=std::min(sI->clippedstartpos,(sE-1)->clippedstartpos);
	  int32 maxs=std::max(sI->clippedstartpos+sI->clippedlen,(sE-1)->clippedstartpos+(sE-1)->clippedlen);
	  tlen=maxs-mins;
	  if((sE-1)->clippedstartpos < sI->clippedstartpos) seg1left=false;
	  int32 grace=(sI->rgid.getInsizeTo()-sI->rgid.getInsizeFrom())/20; // 5% grace (for gap columns etc.)
	  if((sI->rgid.getInsizeFrom() >= 0 && tlen<sI->rgid.getInsizeFrom()-grace)
	     || (sI->rgid.getInsizeTo() >=0 && tlen>sI->rgid.getInsizeTo()+grace)){
	    eachsegmentproperlyaligned=false;
	    CEBUG("BAH 4\n");
	  }
	}

      }

      auto inextI=sI;
      ++inextI;
      for(auto iactI=sI; true; ++iactI, ++inextI){
	bool nextsegmentmapped=false;
	if(inextI!=sE && iactI->template_segment!=0){
	  if(iactI->template_segment >0
	     && (inextI->template_segment == 255
		 || (inextI->template_segment == iactI->template_segment+1))){
	    nextsegmentmapped=true;
	  }
	}
	auto actsamflags=samflags;
	if(eachsegmentproperlyaligned) actsamflags|=0x2;
	if(inextI->dir<0) actsamflags|=0x20;
	if(iactI->dir<0) actsamflags|=0x10;
	if(iactI->template_segment==1) actsamflags|=0x40;  // first segment in template
	if(iactI->template_segment==255) actsamflags|=0x80;  // last segment in template
	if(iactI->template_segment>1 && iactI->template_segment<255) actsamflags|=0x40 | 0x80;  // part of linear segment, set 0x40 & 0x80
	if(!nextsegmentmapped && iactI->template_segment!=255) actsamflags|=0x8;
	iactI->samflags=actsamflags;

	// rnext, pnext
	if(nextsegmentmapped){
	  iactI->rnext_conid=inextI->contigid;
	  iactI->pnext=inextI->clippedstartpos;
	}else if(iactI->template_segment==255 && sI->template_segment==1){
	  iactI->rnext_conid=sI->contigid;
	  iactI->pnext=sI->clippedstartpos;
	}

	// tlen
	if(eachsegmentproperlyaligned){
	  if(iactI->template_segment==1){
	    if(seg1left){
	      iactI->tlen=tlen;
	    }else{
	      iactI->tlen=-tlen;
	    }
	  }else if(iactI->template_segment==255){
	    if(seg1left){
	      iactI->tlen=-tlen;
	    }else{
	      iactI->tlen=tlen;
	    }
	  }
	}

	if(inextI==sE) break;
      }
    }

/*
    cout << sE-sI << " reads for " << SAMC_templatenames[sI->templateid] << endl;
    for(auto xI=sI; xI!=sE; ++xI) {
      cout << "sI tid: " << xI->templateid << "\tsegid: " << static_cast<uint32>(xI->template_segment)
	   << "\tpos: " << xI->clippedstartpos
	   << "\tdir: " << static_cast<int16>(xI->dir)
	   << "\tflag: " << hex << xI->samflags << dec
	   << endl;
    };
*/

  }

  FUNCEND();
}
#define CEBUG(bla)


// having define instead of inline function makes BOOST header disappear from .H
#define SPLITMAFLINE(numexpected) {boost::split(mafsplit, mafline, boost::is_any_of("\t")); if(mafsplit.size()!=numexpected){cout << "Oooops, expected " << numexpected << " elements but found " << mafsplit.size() << "???\n";errorMsgMAFFormat(mafname,linenumber,mafline,"wrong number of elements in line");}};

void SAMCollect::collectInfoFromMAF(const std::string & mafname)
{
  FUNCSTART("void SAMCollect::collectInfoFromMAF(const std::string & mafname)");

  static const std::string cpsHReadGroupShort("@R");
  static const std::string cpsHReadGroup("@ReadGroup");

  static const std::string cpsCO("CO");
  static const std::string cpsCS("CS");
  static const std::string cpsEC("EC");

  static const std::string cpsRD("RD");
  static const std::string cpsRG("RG");
  static const std::string cpsTN("TN");
  static const std::string cpsTS("TS");
  static const std::string cpsER("ER");
  static const std::string cpsAT("AT");

  std::ifstream mafin(mafname, std::ios::in|std::ios::ate);
  if(!mafin) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << mafname << " not found for loading.");
  }
  if(!mafin.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << mafname << " is empty.");
  }
  mafin.seekg(0, std::ios::beg);

  // template id 0 for no template
  if(SAMC_templatenames.empty()) SAMC_templatenames.push_back("");
  if(SAMC_tname2tid.empty()) SAMC_tname2tid.insert(std::pair<std::string,size_t>("",0));

  bool sawreadrg=false;
  size_t actcontigid=0;
  size_t acttemplateid=0;
  int32  acttemplatesegment=-1;
  std::string actcontigname;
  std::string actreadname;
  ReadGroupLib::ReadGroupID actrgid;

  std::string mafline;
  std::string maftoken;
  std::vector<std::string> mafsplit;
  uint64 linenumber=0;

  std::pair<std::unordered_map<std::string, size_t>::iterator,bool> lastsamri;
  lastsamri.second=false;

  while(true){
    ++linenumber;
    getline(mafin,mafline);
    if(mafin.eof()) break;

    CEBUG("l: " << linenumber << "\tt: ###" << mafline << "###" << endl);

    if(mafline.size()<2) continue;
    maftoken=mafline.substr(0,2);
    if(maftoken==cpsRD){
      if(actreadname.size()>0){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found RD token while already in read, did not close previous read?");
      }
      SPLITMAFLINE(2);
      actreadname=std::move(mafsplit[1]);
      if(actreadname.empty()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found RD token without a read name?");
      }
    }else if(maftoken==cpsER){
      if(actreadname.empty()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found ER token while not in read (missing RD token?)");
      }
      if(actrgid.isDefaultNonValidReadGroupID()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"read has no read group (missing RG token?)");
      }
      if(acttemplatesegment==-1){
	acttemplatesegment=0;
      }
      BUGIFTHROW(actcontigid==0,"ER token, actcontigid==0 ??");
      SAMC_samris.push_back(samrinfo_t(SAMC_samris.size(),actcontigid-1,acttemplateid,acttemplatesegment,actrgid,-1,-1,0,0));
      lastsamri=SAMC_rname2samriid.insert(std::pair<std::string,size_t>(actreadname,SAMC_samris.size()-1));
      BUGIFTHROW(!lastsamri.second,"!lastsamri.second ???")
      actreadname.clear();
      acttemplateid=0;
      acttemplatesegment=-1;
      actrgid.resetLibId();
    }else if(maftoken==cpsRG){
      sawreadrg=true;
      if(actreadname.empty()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found RG token while not in read (missing RD token?)");
      }
      if(!actrgid.isDefaultNonValidReadGroupID()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"read already in read group (double RG token?)");
      }
      SPLITMAFLINE(2);
      // convoluted way to assign readgroupid, but hey
      actrgid=ReadGroupLib::ReadGroupID(atoi(mafsplit[1].c_str()));
    }else if(maftoken==cpsTN){
      SPLITMAFLINE(2);
      auto tnI=SAMC_tname2tid.find(mafsplit[1]);
      if(tnI==SAMC_tname2tid.end()){
	acttemplateid=SAMC_templatenames.size();
	SAMC_templatenames.push_back(mafsplit[1]);
	SAMC_tname2tid.insert(std::pair<std::string,size_t>(mafsplit[1],acttemplateid));
      }else{
	acttemplateid=tnI->second;
      }
    }else if(maftoken==cpsTS){
      if(acttemplatesegment!=-1){
	errorMsgMAFFormat(mafname,linenumber,mafline,"TS already set for this read?");
      }
      SPLITMAFLINE(2);
      acttemplatesegment=atoi(mafsplit[1].c_str());
    }else if(maftoken==cpsAT){
      if(!actreadname.empty()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found AT token while in read (missing ER token?)");
      }
      BUGIFTHROW(actcontigid==0,"ER token, actcontigid==0 ??");
      if(SAMC_samris.empty()){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found AT token without previous read ?)");
      }
      if(SAMC_samris.back().clippedstartpos>=0){
	errorMsgMAFFormat(mafname,linenumber,mafline,"last read seems already to have positions set ... double AT line ?)");
      }
      SPLITMAFLINE(5);
      int32 first=atoi(mafsplit[1].c_str());
      int32 second=atoi(mafsplit[2].c_str());
      SAMC_samris.back().clippedstartpos=std::min(first,second);
      SAMC_samris.back().dir=1;
      if(first>second) SAMC_samris.back().dir=-1;
      first=atoi(mafsplit[3].c_str());
      second=atoi(mafsplit[4].c_str());
      SAMC_samris.back().clippedlen=second-first;
    }else if(maftoken==cpsCO){
      if(actcontigid!=0){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found CO token while already in contig? Missed an EC token.");
      }
      SPLITMAFLINE(2);
      actcontigname=std::move(mafsplit[1]);
      auto oldsize=SAMC_namecheck_contig.size();
      SAMC_namecheck_contig.insert(actcontigname);
      if(SAMC_namecheck_contig.size() == oldsize){
	cout << "Duplicate: " << actcontigname << endl;
	errorMsgMAFFormat(mafname,linenumber,mafline,"duplicate contig name?");
      }
      SAMC_contignames.push_back(actcontigname);
      actcontigid=SAMC_contignames.size();
    }else if(maftoken==cpsCS){
      if(actcontigid==0){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found CS token while not in contig?");
      }
      // I'm lazy
      mafline[0]=' ';
      mafline[1]=' ';
      boost::trim(mafline);
      SAMC_contiglengths.push_back(mafline.size());
    }else if(maftoken==cpsEC){
      if(actcontigid==0){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found EC token without preceding CO token?");
      }
      actcontigid=0;
    }else if(maftoken==cpsHReadGroupShort
	     && mafline.size()>=10
	     && mafline.substr(0,10)==cpsHReadGroup
      ){
      if(sawreadrg){
	errorMsgMAFFormat(mafname,linenumber,mafline,"found @Readgroup in MAF while other reads have already been defined. Cannot do that yet.\n\nUse 'grep ^@ old.maf >new.maf; grep -v ^@ old.maf >>new.maf; ' to create a MAF which can be converted to SAM.");
      }
      std::vector<ReadGroupLib::ReadGroupID> dummy_externalidmapper;
      auto dummy_rgid=ReadGroupLib::newReadGroup();
      MAFParse::parseReadGroup(mafin, dummy_rgid, dummy_externalidmapper, linenumber);
    }
  }

  if(SAMC_contiglengths.size() != SAMC_contignames.size()){
    MIRANOTIFY(Notify::FATAL,"different number of contigs and contig lengths? were there contigs without sequences?");
  }

  FUNCEND();
}


void SAMCollect::errorMsgMAFFormat(const std::string & filename, size_t linenumber, const std::string & line, const char * msg)
{
  cout << "MAF file " << filename << " at line " << linenumber << ":\n" << line << "\n\n" << msg << endl;
  cout << "Error in MAF format, aborting." << endl;
  exit(100);
}


void SAMCollect::createSAMHeader()
{
  std::stringstream ostr;
  // readgroup 0 should never be present anyway
  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    ReadGroupLib::dumpReadGroupAsSAM(rgi,ostr);
  }

  for(size_t ci=0; ci<SAMC_contignames.size(); ++ci){
    ostr << "@SQ\tSN:" << SAMC_contignames[ci]
	 << "\tLN:" << SAMC_contiglengths[ci]
	 << '\n';
  }
  SAMC_headerstring=ostr.str();
}


/*************************************************************************
 *
 * returns false if not found
 * true if found and the samrinfo_t structure in samri
 *
 *************************************************************************/

bool SAMCollect::getSAMRInfo(const std::string & readname, samrinfo_t & samri) const
{
  auto srI=SAMC_rname2samriid.find(readname);
  if(srI==SAMC_rname2samriid.end()) return false;
  samri=SAMC_samris[srI->second];
  return true;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

const std::string & SAMCollect::getContigName(samrinfo_t & samri) const
{
  FUNCSTART("const std::string & SAMCollect::getContigName(samrinfo_t & samri) const");
  BUGIFTHROW(samri.contigid >= SAMC_contignames.size(),"samri.contigid >= SAMC_contignames.size() ???");
  return SAMC_contignames[samri.contigid];
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

const std::string & SAMCollect::getRNextEntry(samrinfo_t & samri) const
{
  static const std::string sstar("*");
  static const std::string sequal("=");
  if(samri.rnext_conid>=0) {
    if(samri.contigid==samri.rnext_conid) return sequal;
    return SAMC_contignames[samri.rnext_conid];
  }
  return sstar;
}
