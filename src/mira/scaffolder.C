/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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


#include "mira/scaffolder.H"

#include "mira/contig.H"


using std::cout;
using std::cerr;
using std::endl;


// Plain vanilla constructor
Scaffolder::Scaffolder()
{
  FUNCSTART("Scaffolder::Scaffolder()");

  zeroVars();
  init();

  FUNCEND();
}

void Scaffolder::zeroVars()
{
  FUNCSTART("void Scaffolder::zeroVars()");
  FUNCEND();
}

void Scaffolder::init()
{
  FUNCSTART("void Scaffolder::init()");

  SCA_tname2tid.insert(std::pair<std::string,uint32>("",0));

  FUNCEND();
}



Scaffolder::~Scaffolder()
{
  FUNCSTART("Scaffolder::~Scaffolder()");

  discard();

  FUNCEND();
}


void Scaffolder::discard()
{
  FUNCSTART("Scaffolder::discard()");

  zeroVars();

  FUNCEND();
}


#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Scaffolder::storeInfoFreshContig(Contig & con)
{
  FUNCSTART("void Scaffolder::storeInfoFreshContig(contig & con)");

  readlinkinfo_t tmplink;
  tmplink.contigid=SCA_contiginfo.size();

  SCA_contiginfo.resize(SCA_contiginfo.size()+1);
  auto & actcinfo=SCA_contiginfo.back();
  actcinfo.name=con.getContigName();

  auto & creads=con.getContigReads();
  {
    auto crE=creads.end();
    for(auto crI=creads.begin(); crI!=crE; ++crI){
      auto & actread=*crI;
      int64 startpos=crI.getReadStartOffset();
      int64 endpos=startpos+actread.getLenClippedSeq();
      auto rgid=actread.getReadGroupID();
      auto wanteddist=std::max(rgid.getInsizeTo(),rgid.getInsizeFrom()); // if wanteddist <= 0: just go by directions
      CEBUG("Wanted: " << wanteddist << endl);
      if(actread.getTemplateSegment() != 0         // do we have a valid template segment?
	 && rgid.getSegmentPlacementCode() != ReadGroupLib::SPLACE_UNKNOWN
	 && (wanteddist <= 0
	     || startpos <= wanteddist
	     || endpos >= static_cast<int64>(con.getContigLength())-wanteddist)){
	// OK ... we are at the right distance from one of the contig ends
	// find out which is the correct one:
	//  from the direction of the read in the contig its segment placement scheme
	tmplink.wantdir=0;
	bool wantleft=false;    // partner read should be on left
	bool wantright=false;   // partner read should be on right
	auto rdir=crI.getReadDirection();
	// this is the logic for rdir > 0
	switch(rgid.getSegmentPlacementCode()){
	case ReadGroupLib::SPLACE_RF : {
	  tmplink.wantdir=-1;
	  wantleft=true;
	  break;
	}
	case ReadGroupLib::SPLACE_FR : {
	  tmplink.wantdir=-1;
	  wantright=true;
	  break;
	}
	case ReadGroupLib::SPLACE_SF : {
	  tmplink.wantdir=1;
	  if(actread.getTemplateSegment()==1){
	    wantright=true;
	  }else if(actread.getTemplateSegment()==255){
	    wantleft=true;
	  }
	  break;
	}
	case ReadGroupLib::SPLACE_SB : {
	  tmplink.wantdir=1;
	  if(actread.getTemplateSegment()==1){
	    wantleft=true;
	  }else if(actread.getTemplateSegment()==255){
	    wantright=true;
	  }
	  break;
	}
	case ReadGroupLib::SPLACE_SU : {
	  tmplink.wantdir=1;
	  wantleft=true;
	  wantright=true;
	  break;
	}
	default : {
	  BUGIFTHROW(true,"rgid.getSegmentPlacement() " << static_cast<int16>(rgid.getSegmentPlacementCode()) << " (" << rgid.getSegmentPlacement() << ") is unexpected here.");
	}
	}
	// this is the logic for rdir < 0
	if(rdir<0){
	  tmplink.wantdir=-tmplink.wantdir;
	  if(rgid.getSegmentPlacementCode() != ReadGroupLib::SPLACE_SU){
	    wantleft=!wantleft;
	    wantright=!wantright;
	  }
	}

	// Good, we know where we want the partner to be, let's see whether we are at the correct distances
	// If not, forget again about what we want
	if(wanteddist>0){
	  if(wantleft && endpos > wanteddist) wantleft=false;
	  if(wantright && startpos < static_cast<int64>(con.getContigLength())-wanteddist) wantright=false;
	}

	// Yay, now store what we want ... and we may want both
	tmplink.dir=rdir;
	// tmplink.templateid=...
	{
	  auto tnI=SCA_tname2tid.find(actread.getTemplate());
	  if(tnI==SCA_tname2tid.end()){
	    tmplink.templateid=SCA_tname2tid.size();
	    SCA_tname2tid.insert(std::pair<std::string,uint32>(actread.getTemplate(),tmplink.templateid));
	  }else{
	    tmplink.templateid=tnI->second;
	  }
	}

	if(wantleft){
	  tmplink.distance=startpos+endpos;
	  SCA_readlinks.push_back(tmplink);
	  CEBUG("Storing " << con.getContigName() << " left " << actread.getName() << " " << tmplink << endl);
	}
	if(wantright){
	  tmplink.distance=con.getContigLength()-startpos;
	  SCA_readlinks.push_back(tmplink);
	  CEBUG("Storing " << con.getContigName() << " right " << actread.getName() << " " << tmplink << endl);
	}
      }
    }
  }

  FUNCEND();
}
#define CEBUG(bla)


void Scaffolder::dumpDebug()
{
  for(const auto & rlie : SCA_readlinks){
    cout << rlie << endl;
  }
}
