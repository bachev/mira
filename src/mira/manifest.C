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

#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <regex>

#include "util/fileanddisk.H"
#include "mira/enums.H"
#include "mira/manifest.H"


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)


// Plain vanilla constructor
Manifest::Manifest()
{
  FUNCSTART("Manifest::Manifest()");

  init();

  FUNCEND();
}

void Manifest::zeroVars()
{
  FUNCSTART("void Manifest::zeroVars()");
  FUNCEND();
}

void Manifest::init()
{
  FUNCSTART("void Manifest::init()");

  MAN_seentechnology.clear();
  MAN_seentechnology.resize(ReadGroupLib::getNumSequencingTypes(),false);

  FUNCEND();
}



Manifest::~Manifest()
{
  FUNCSTART("Manifest::~Manifest()");

  FUNCEND();
}


void Manifest::discard()
{
  FUNCSTART("Manifest::discard()");

  zeroVars();

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

/*
readgroup = Bla0
file fastq = x1.fastq x2.fastq ...
file xml = x1.xml ...
technology = sanger 454 illumina iontorrent pacbio lq|hq
templatesize = 340 780
templateordering = >< <> >> strictordering
templatenaming = solexa fr
strainname = xxxyyyzzz

include = manifestfile

// TODO: fastq33 fastq...

*/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Manifest::loadManifestFile(const std::string & mfilename, bool resume)
{
  FUNCSTART("void Manifest::loadManifestFile(std::string & mfilename)");

  std::stringstream buffer;
  {
    std::ifstream mfin(mfilename);
    if(!mfin.good()) {
      cout << "Did not find manifest file " << mfilename << endl;
      exit(10);
      return;
    }
    buffer << mfin.rdbuf();
  }

  slurpInManifest(buffer,mfilename,resume);

  //cout << "Readgroups after loading manifest:\n";
  //ReadGroupLib::debugDumpReadGroupInfo(cout);

  if(!resume){
    cout << "Looking for files named in data ..."; cout.flush();
    if(provideFileNames()){
      MIRANOTIFY(Notify::FATAL,"Some 'data' entries named in the manifest file could not be verified, see the log above.\nMaybe some files are missing, not readable or there is a typo in the manifest file?");
    }
  }

  FUNCEND();
}

void Manifest::slurpInManifest(std::stringstream & mfin, const std::string & origsource, bool resume)
{
  FUNCSTART("void Manifest::slurpInManifest(std::stringstream & mfin, const std::string & origsource, bool resume))");

  boost::char_separator<char> separator(" \t");

  auto valuestart=std::string::npos;

  std::string lastkeyok;
  std::string key;
  std::string value;

  std::string tmpline;
  std::string actline;
  bool foundreadgroup=false;
  uint32 numbackslashes=0;
  while(getline(mfin,tmpline)){
    boost::trim(tmpline);
    if(tmpline.empty()
       || tmpline[0]=='#') continue;

    actline=tmpline;

    while(!actline.empty() && actline[actline.size()-1]=='\\'){
      ++numbackslashes;
      actline.resize(actline.size()-1);
      tmpline.clear();
      if(getline(mfin,tmpline)){
	boost::trim(tmpline);
	actline+=' ';
	actline+=tmpline;
      }
    }

    valuestart=actline.find_first_of('=',0);
    //if(valuestart==std::string::npos){
    //  MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nhas no equal sign?");
    //}

    lastkeyok.swap(key);
    key=actline.substr(0, valuestart);
    boost::trim(key);
    if(valuestart!=std::string::npos){
      value=actline.substr(valuestart+1);
      boost::trim(value);
    }else{
      value.clear();
    }

    CEBUG("k: " << key << endl);
    CEBUG("v: " << value << endl);

    boost::to_lower(key);

    if(key=="include"){
      loadManifestFile(value,resume);
    }else if(key=="projectname"
	     || key=="project_name"
	     || key=="project"){
      if(!MAN_projectname.empty()){
	MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nFound keyword 'projectname', but there already has been a job defined previously ('" << MAN_projectname << "')");
      }
      MAN_projectname=value;
    }else if(key=="job"){
      if(!MAN_job.empty()){
	MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nFound keyword 'job', but there already has been a job defined previously ('" << MAN_job << "')");
      }
      MAN_job=value;
    }else if(key=="parameters"){
      if(!MAN_parameters.empty()) MAN_parameters+=' ';
      MAN_parameters+=value;
    }else if(key=="readgroup"
	     || key=="read_group"){
      if(foundreadgroup){
	// some defaults and checks for previous readgroup
	MAN_manifestdata2load.back().rgid.fillInSensibleDefaults();
	MAN_manifestdata2load.back().rgid.checkValidity();
      }
      foundreadgroup=true;
      MAN_manifestdata2load.resize(MAN_manifestdata2load.size()+1);
      MAN_manifestdata2load.back().rgid=ReadGroupLib::newReadGroup();
      MAN_manifestdata2load.back().rgid.setGroupName(value);
      MAN_manifestdata2load.back().loadasbackbone=false;
    }else{
      if(!foundreadgroup){
	if(lastkeyok=="parameters"){
	  if(numbackslashes==0){
	    MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nhas keyword '" << key << "' which is not recognised.\nThe last recognised keyword was 'parameters', no continuation line was seen. Is this a continuation line and did you eventually forget the backslash in the line above?");
	  }else{
	    MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nhas keyword '" << key << "' which is not recognised.\nThe last recognised keyword was 'parameters' and there were continuation lines with backslashes. Did you eventually erroneously put a backslash in the last continuation line above?");
	  }
	}else{
	  MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nExpected keyword 'readgroup' not found. Maybe you mistyped '" << key << "' and meant something else?\nAccepted keywords at this stage are 'include', 'project', 'job', 'parameters' and 'readgroup'.");
	}
      }

      if(key=="as_reference"
	 || key=="is_reference"){
	MAN_manifestdata2load.back().loadasbackbone=true;
	MAN_manifestdata2load.back().rgid.setBackbone(true);
      }else if(key=="autopairing"
	       || key=="autopair"
	       || key=="autotemplate"){
	MAN_manifestdata2load.back().rgid.setInsizeFrom(-1);
	MAN_manifestdata2load.back().rgid.setInsizeTo(-1);
	MAN_manifestdata2load.back().rgid.setWantTemplateSizeEstimate(true);
	MAN_manifestdata2load.back().rgid.setWantSegmentPlacementEstimate(true);
	MAN_manifestdata2load.back().rgid.setExpectsReadPairs(true);
      }else if(key=="nostats"
	       || key=="nostatistics"){
	MAN_manifestdata2load.back().rgid.setStatisticsCalc(false);
      }else if(key=="technology"){
	try{
	  MAN_manifestdata2load.back().rgid.setSequencingType(value);
	}
	catch(Notify n){
	  MIRANOTIFY(Notify::FATAL,"In file " << origsource << ": line\n" << actline << "\nUnknown technology '" << value << "' found. Did you make a typo?");
	}
	if(MAN_manifestdata2load.back().rgid.getSequencingType()==ReadGroupLib::SEQTYPE_PACBIOLQ){
	  //MIRANOTIFY(Notify::FATAL,"Sorry, PacBio low quality data currently not supported.\n");
	}
      }else if(key=="defaultqual"
	       || key=="default_qual"){
	int64 dq=atoll(value.c_str());
	if(dq<0 || dq>100){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nKeyword '" << key << "' values 0 <= x <= 100, found " << dq);
	}
	MAN_manifestdata2load.back().rgid.setDefaultQual(static_cast<base_quality_t>(dq));
      }else if(key=="templatesize"
	       || key=="template_size"){
	boost::tokenizer<boost::char_separator<char> > tok(value,separator);
	auto tI=tok.begin();
	if(tI==tok.end()){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nKeyword '" << key << "' expects two or three values, found none.");
	}
	int32 dummy=atoi(tI->c_str());
	MAN_manifestdata2load.back().rgid.setInsizeFrom(dummy);
	MAN_manifestdata2load.back().rgid.setWantTemplateSizeEstimate(false);
	MAN_manifestdata2load.back().rgid.setExpectsReadPairs(true);
	++tI;
	if(tI==tok.end()){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nKeyword '" << key << "' expects exactly two or three values, found only one.");
	}
	dummy=atoi(tI->c_str());
	MAN_manifestdata2load.back().rgid.setInsizeTo(dummy);
	bool foundautorefine=false;
	while(++tI!=tok.end()){
	  // oops, 3 or 4 values
	  std::string tmpiostring(*tI);
	  boost::to_lower(tmpiostring);
	  if(tmpiostring=="infoonly" || tmpiostring=="exclusion_criterion"
	    || tmpiostring=="autorefine"){
	    MAN_manifestdata2load.back().ts_infoonlygiven=true;
	    if(tmpiostring=="autorefine"){
	      foundautorefine=true;
	    }else if(tmpiostring=="infoonly"){
	      MAN_manifestdata2load.back().rgid.setTSInfoOnly(true);
	    }else{
	      MAN_manifestdata2load.back().rgid.setTSInfoOnly(false);
	    }
	  }else{
	    MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nKeyword '" << key << "' expects exactly 'infoonly', 'exclusion_criterion' or 'autorefine' as third or fourth value, found '" << tmpiostring << "'");
	  }
	}
	if(foundautorefine){
	  MAN_manifestdata2load.back().rgid.setWantTemplateSizeEstimate(true);
	}
      }else if(key=="prefix_rename"
	       || key=="rename_prefix"){
	boost::tokenizer<boost::char_separator<char> > tok(value,separator);
	auto tI=tok.begin();
	if(tI==tok.end()){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nKeyword '" << key << "' expects two values, found none.");
	}
	auto prefix=*tI;
	++tI;
	if(tI==tok.end()){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nKeyword '" << key << "' expects exactly two values, found only one.");
	}
	MAN_manifestdata2load.back().rgid.addReadRenamePrefix(prefix,*tI);
      }else if(key=="segmentplacement"
	       || key=="segment_placement"){
	// handle "infoonly"/"exclusion_criterion" by parsing it out
	static std::regex iostr( "infoonly" ) ;
	static std::regex ecstr( "exclusion_criterion" ) ;
	static std::regex arstr( "autorefine" ) ;
	static std::string iostrrep; // empty
	auto newval = std::regex_replace(value, iostr, iostrrep);
	if(newval!=value){
	  MAN_manifestdata2load.back().rgid.setSPInfoOnly(true);
	  MAN_manifestdata2load.back().sp_infoonlygiven=true;
	  swap(value,newval); // we just want value=newval
	}else{
	  newval = std::regex_replace(value, ecstr, iostrrep);
	  if(newval!=value){
	    MAN_manifestdata2load.back().rgid.setSPInfoOnly(false);
	    MAN_manifestdata2load.back().sp_infoonlygiven=true;
	    swap(value,newval); // we just want value=newval
	  }
	}

	MAN_manifestdata2load.back().rgid.setWantSegmentPlacementEstimate(false);
	MAN_manifestdata2load.back().rgid.setExpectsReadPairs(true);
	newval = std::regex_replace(value, arstr, iostrrep);
	if(newval!=value){
	  MAN_manifestdata2load.back().rgid.setWantSegmentPlacementEstimate(true);
	  swap(value,newval); // we just want value=newval
	}

	if(!MAN_manifestdata2load.back().rgid.setSegmentPlacement(value)){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nUnrecognised value '" << value << "' for key '" << key << "'");
	}
	if(MAN_manifestdata2load.back().rgid.getSegmentPlacementCode()!=ReadGroupLib::SPLACE_UNKNOWN
	   && MAN_manifestdata2load.back().rgid.getSegmentPlacementCode()!=ReadGroupLib::SPLACE_SU
	   && MAN_manifestdata2load.back().rgid.wantSegmentPlacementEstimate()){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nUsing 'autorefine' for segment_placement other than 'unknown' or 'samedir' does not make much sense. Please fix your manifest file.");
	}
      }else if(key=="segmentnaming"
	       || key=="segment_naming"){
	static std::regex csstr( "rollcomment" ) ;
	static std::string repstr; // empty
	auto newval = std::regex_replace(value, csstr, repstr);
	if(newval!=value){
	  MAN_manifestdata2load.back().rgid.setUseReadNameFromComment(true);
	  swap(value,newval); // we just want value=newval
	  boost::trim(value);
	}
	if(!MAN_manifestdata2load.back().rgid.setSegmentNaming(value)){
	  MIRANOTIFY(Notify::FATAL,  "In file " << origsource << ": line\n" << actline << "\nUnrecognised value '" << value << "' for key '" << key << "'");
	}
      }else if(key=="strainname"
	       || key=="strain_name"
	       || key=="strain"){
	boost::trim(value);
	MAN_manifestdata2load.back().rgid.setStrainName(value);
      }else if(key=="datadir_scf"){
	boost::trim(value);
	MAN_manifestdata2load.back().rgid.setDataDir(value);
      }else if(key.substr(0,4)=="data"){
	boost::tokenizer<boost::char_separator<char> > tok(value,separator);
	for(const auto & te : tok){
	  MAN_manifestdata2load.back().datanames.push_back(te);
	}
      }else{
	MIRANOTIFY(Notify::FATAL, "In file " << origsource << ": line\n" << actline << "\nhas keyword '" << key << "' which is not recognised.");
      }
    }
  }

  if(foundreadgroup){
    // check last defined read group for validity
    MAN_manifestdata2load.back().rgid.fillInSensibleDefaults();
    MAN_manifestdata2load.back().rgid.checkValidity();
  }

  // quick consistency checks
  bool seenbb=false;
  bool seennotbb=false;
  for(auto & me : MAN_manifestdata2load){
    // check for presence of filenames
    if(me.datanames.empty()){
      MIRANOTIFY(Notify::FATAL,"Oooops, the readgroup '" << me.rgid.getGroupName() << "' has no data defined for loading? Check your manifest.");
    }

    // check for setting of sequencing tech
    if(me.rgid.getSequencingType()==ReadGroupLib::SEQTYPE_END){
      if(me.loadasbackbone) {
	// backbones do not need technolgies (but user can choose so)
	// by default, they're TEXT
	me.rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);
      }else{
	cout << "Oooops, the readgroup '" << me.rgid.getGroupName() << "' has no sequencing technology defined.\n"
	     << "Files in this readgroup:";
	for(const auto & fn : me.datanames){
	  cout << " " << fn;
	}
	cout << endl;
	MIRANOTIFY(Notify::FATAL,"Error in defining readgroups: missing technology setting. See output above for more info.");
      }
    }

    // put together technology string
    if(!me.loadasbackbone){
      seennotbb=true;
      uint8 st=me.rgid.getSequencingType();
      if(!MAN_seentechnology[st]){
	MAN_seentechnology[st]=true;
	if(!MAN_technologystring.empty()) MAN_technologystring+=',';
	MAN_technologystring+=ReadGroupLib::getNameOfSequencingType(st);
      }
    }else{
      seenbb=true;
    }
  }

  if(!seennotbb){
    if(seenbb){
      MIRANOTIFY(Notify::FATAL,"Error in defining readgroups: there is no readgroup which is not defined as reference/backbone. Basically this says that while you defined a reference, you did not define what data should map to it.");
    }else{
      MIRANOTIFY(Notify::FATAL,"Error in defining readgroups: there is no readgroup defined containing data to load? What should MIRA assemble???");
    }
  }

  // handle missing "infoonly" for non-backbone readgroups in mapping assemblies:
  // set default "infoonly"
  if(seenbb){
    for(auto & me : MAN_manifestdata2load){
      if(!me.loadasbackbone){
	if((me.rgid.getInsizeFrom()>=0 || me.rgid.getInsizeTo()>=0)
	   && !me.ts_infoonlygiven) {
	  cout << "For mapping assembly: readgroup " << me.rgid.getLibId() << " named '" << me.rgid.getGroupName() << "' has no 'infoonly' or 'exclusion_criterion' set for 'template_segment',\nassuming 'infoonly'.\n";
	  me.rgid.setTSInfoOnly(true);
	}
	if(me.rgid.getSegmentPlacementCode()!=ReadGroupLib::SPLACE_UNKNOWN
	   && !me.sp_infoonlygiven) {
	  cout << "For mapping assembly: readgroup " << me.rgid.getLibId() << " named '" << me.rgid.getGroupName() << "' has no 'infoonly' or 'exclusion_criterion' set for 'segment_placement',\nassuming 'infoonly'.\n";
	  me.rgid.setSPInfoOnly(true);
	}
      }
    }
  }

  FUNCEND();
}


std::string Manifest::getFullMIRAParameterString()
{
  FUNCSTART("std::string Manifest::getFullMIRAParameterString()");

  std::string ret("--job=");

  uint32 numchars=0;
  for(auto & je : MAN_job){
    if(!isspace(je)) {
      ret.push_back(je);
      ++numchars;
    }
  }

  if(!numchars){
    MIRANOTIFY(Notify::FATAL,"Found no job description in the manifest loaded.\nYou need to provide MIRA with information what it should do!\nE.g.: job=denovo,genome,accurate");
  }

  if(MAN_technologystring.empty()){
    MIRANOTIFY(Notify::FATAL,"Found no technology descriptions in the manifest loaded.\nYou need to provide MIRA with information what type of data it is loading. Although somehow MIRA should have caught this earlier ...");
  }

  ret+=',';
  ret+=MAN_technologystring;
  ret+=' ';

  ret+=MAN_parameters;

  FUNCEND();

  return ret;
}



std::ostream & operator<<(std::ostream &ostr, Manifest const &m)
{
  FUNCSTART("friend std::ostream & Manifest::operator<<(std::ostream &ostr, const  &m)");

  cout << "Manifest:\nprojectname: " << m.MAN_projectname
       << "\njob: " << m.MAN_job
       << "\nparameters: " << m.MAN_parameters << endl;

  cout << "Manifest load entries: " << m.MAN_manifestdata2load.size() << endl;
  uint32 mlei=0;
  for(const auto & mle : m.MAN_manifestdata2load){
    ++mlei;
    cout << "MLE " << mlei << ":\n";
    cout << mle.rgid;
    for(const auto & fn :  mle.datanames){
      cout << fn << " ";
    }
  }
  cout << endl;

  FUNCEND();
  return ostr;
}



bool Manifest::provideFileNames()
{
  FUNCSTART("bool Manifest::provideFileNames()");

  // for delayed transformation reads to contigs
  // needed to simplify life when loading .fna and .gff3 with annotations
  // this gives the Readpool loader the opportunity to load the sequence,
  //  then the annotation, map the annotation to the sequence
  // then only make the contigs
  //
  // if this wasn't done, one would have to annotate both the readpool read
  //   AND search contigs for this read and annotate that one too
  // TODO: might be an idea for "re-annotation"

  std::list<fnft_t> tmpgff3;
  bool hassomeerror=false;
  for(auto & mle : MAN_manifestdata2load){
    mle.mainfilesfoundfordata.clear();
    mle.ancillaryfilesfoundfordata.clear();

    std::list<fnft_t> fnftl;
    for(const auto & dn : mle.datanames){
      hassomeerror|=globWalkPath(dn,fnftl);
    }

    // now fill in maindata and ancillary, leaving out GFF3 file types
    //  then append the remaining files (well, GFF3) to the main data
    //  this effectively splits the files into maindata and ancillary AND
    //  shuffles GFF3 to the end of maindata
    //  ... and keeps the order of files (not really needed, but good for
    //  user to see his files get loaded in the order he named them, except
    //  GFF3 that is)

    auto ntI=fnftl.begin();
    while(ntI!=fnftl.end()){
      auto tmpI=ntI;
      ++ntI;
      if(tmpI->ft!="gff3"){
	if(tmpI->ft=="xml"
	   || tmpI->ft=="ssaha2"
	   || tmpI->ft=="smalt"){
	  mle.ancillaryfilesfoundfordata.splice(mle.ancillaryfilesfoundfordata.end(),fnftl,tmpI);
	}else{
	  mle.mainfilesfoundfordata.splice(mle.mainfilesfoundfordata.end(),fnftl,tmpI);
	}
      }
    }
    // what's left over is GFF3
    mle.mainfilesfoundfordata.splice(mle.mainfilesfoundfordata.end(),fnftl);
  }

  return hassomeerror;
}
