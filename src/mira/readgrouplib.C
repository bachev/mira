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


#include "readgrouplib.H"

#include <unordered_map>

#include <boost/algorithm/string.hpp>


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)

bool ReadGroupLib::staticInitialiser()
{
  // not really nice to do the following, but what the heck
  // --
  std::vector<std::string> & lcnames=const_cast<std::vector<std::string>&>(RG_lcnamesofseqtypes);
  lcnames=RG_namesofseqtypes;
  for(size_t i=0; i<lcnames.size(); ++i){
    boost::to_lower(lcnames[i]);
  }
  //--

  RG_static_infolib.reserve(255);

  zeroVars();

  return true;
}



// Plain vanilla constructor
ReadGroupLib::ReadGroupLib()
{
  FUNCSTART("ReadGroupLib::ReadGroupLib()");

  zeroVars();

  FUNCEND();
}

void ReadGroupLib::zeroVars()
{
  FUNCSTART("void ReadGroupLib::zeroVars()");

  RG_strainids_clean=false;
  RG_numstrains=-1;

  RG_static_infolib.clear();
  newReadGroup();
  completeDefaultsForReadGroups();

  FUNCEND();
}

void ReadGroupLib::init()
{
  FUNCSTART("void ReadGroupLib::init()");

  FUNCEND();
}



ReadGroupLib::~ReadGroupLib()
{
  FUNCSTART("ReadGroupLib::~ReadGroupLib()");

  discard();

  FUNCEND();
}


void ReadGroupLib::discard()
{
  FUNCSTART("ReadGroupLib::discard()");

  zeroVars();

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//ReadGroupLib::ReadGroupLib(ReadGroupLib const &other)
//{
//  FUNCSTART("ReadGroupLib::ReadGroupLib(ReadGroupLib const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//ReadGroupLib const & ReadGroupLib::operator=(ReadGroupLib const & other)
//{
//  FUNCSTART("ReadGroupLib const & ReadGroupLib::operator=(ReadGroupLib const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//std::ostream & operator<<(std::ostream &ostr, ReadGroupLib const &???)
//{
//  FUNCSTART("friend std::ostream & ReadGroupLib::operator<<(std::ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}



ReadGroupLib::ReadGroupID ReadGroupLib::newReadGroup()
{
  FUNCSTART("ReadGroupLib::ReadGroupID ReadGroupLib::newReadGroup()");
  BUGIFTHROW(RG_static_infolib.size()==255,"Maximum number of readgroups reached? What kind of data do you have there?!");
  RG_static_infolib.resize(RG_static_infolib.size()+1);
  RG_strainids_clean=false;
  return ReadGroupID(RG_static_infolib.size()-1);
}



void ReadGroupLib::checkLibValidity(rgid_t libid)
{
  FUNCSTART("static void checkLibValidity(rgid_t libid)");
  if(libid==0) abort();
  BUGIFTHROW(RG_static_infolib[libid].seqtype>=SEQTYPE_END, "Readgroup " << static_cast<uint16>(libid) << " (named: '" << RG_static_infolib[libid].groupname << "') has no valid sequencing technology set.");

  BUGIFTHROW(RG_static_infolib[libid].defaultqual>=100, "Readgroup " << static_cast<uint16>(libid) << " (named: '" << RG_static_infolib[libid].groupname << "') has a default quality of " << static_cast<uint16>(RG_static_infolib[libid].defaultqual) << " which is >= 100 ... invalid.");
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadGroupLib::fillInSensibleDefaults(rgid_t libid)
{
  FUNCSTART("void ReadGroupLib::fillInSensibleDefaults(rgid_t libid)");

  if(libid==0) {
    // no "sensible" defaults for libid 0, try to set values which will let MIRA puke
    RG_static_infolib[0].seqtype=SEQTYPE_END;
    RG_static_infolib[0].namingscheme=SCHEME_UNKNOWN;
    RG_static_infolib[0].defaultqual=212;
    RG_static_infolib[0].groupname="UnassigneReads";
    RG_static_infolib[0].machine_type="ShouldNeverBeSeen";
    return;
  }

  // check for sequencing type
  if(RG_static_infolib[libid].seqtype==SEQTYPE_END){
    if(RG_static_infolib[libid].is_backbone) {
      // backbones do not need technolgies (but user can choose so)
      // by default, they're TEXT
      RG_static_infolib[libid].seqtype=SEQTYPE_TEXT;
    }else{
      MIRANOTIFY(Notify::FATAL,"Oooops, the readgroup '" << RG_static_infolib[libid].groupname << "' has no sequencing technology defined, nor is it defined as reference (which would excuse the missing technology definition).");
    }
  }
  // check for strain name
  if(RG_static_infolib[libid].strainname.empty()){
    if(RG_static_infolib[libid].is_backbone) {
      RG_static_infolib[libid].strainname="ReferenceStrain";
      if(RG_static_infolib[libid].defaultqual>100) RG_static_infolib[libid].defaultqual=30;
    }else{
      RG_static_infolib[libid].strainname="StrainX";
    }
  }
  // check for default qual
  RG_static_infolib[libid].has_userdefaultqual=true;
  if(RG_static_infolib[libid].defaultqual>=100){
    RG_static_infolib[libid].has_userdefaultqual=false;
    base_quality_t def=0;
    if(RG_static_infolib[libid].is_backbone) {
      def=30;
    }else{
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      switch(RG_static_infolib[libid].seqtype){
      case ReadGroupLib::SEQTYPE_TEXT : {
	def=10;
	break;
      }
      case ReadGroupLib::SEQTYPE_SANGER :
      case ReadGroupLib::SEQTYPE_454GS20 :
      case ReadGroupLib::SEQTYPE_IONTORRENT :
      case ReadGroupLib::SEQTYPE_PACBIOHQ : {
	def=20;
	break;
      }
      case ReadGroupLib::SEQTYPE_SOLEXA : {
	def=30;
	break;
      }
      default : {
	def=5;
      }
      }
    }
    RG_static_infolib[libid].defaultqual=def;
  }

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  // read naming scheme
  if(RG_static_infolib[libid].namingscheme==SCHEME_UNKNOWN){
    if(RG_static_infolib[libid].is_backbone) {
      RG_static_infolib[libid].namingscheme=SCHEME_NONE;
    }else{
      switch(RG_static_infolib[libid].seqtype){
      case ReadGroupLib::SEQTYPE_TEXT : {
	RG_static_infolib[libid].namingscheme=SCHEME_NONE;
	break;
      }
      case ReadGroupLib::SEQTYPE_SANGER : {
	RG_static_infolib[libid].namingscheme=SCHEME_SANGER;
	break;
      }
      case ReadGroupLib::SEQTYPE_454GS20 : {
	RG_static_infolib[libid].namingscheme=SCHEME_SOLEXA;
	break;
      }
      case ReadGroupLib::SEQTYPE_IONTORRENT : {
	RG_static_infolib[libid].namingscheme=SCHEME_SOLEXA;
	break;
      }
      case ReadGroupLib::SEQTYPE_PACBIOHQ :  {
	RG_static_infolib[libid].namingscheme=SCHEME_SOLEXA;
	break;
      }
      case ReadGroupLib::SEQTYPE_SOLEXA : {
	RG_static_infolib[libid].namingscheme=SCHEME_SOLEXA;
	break;
      }
      default : {
	MIRANOTIFY(Notify::FATAL,"Unknown seqtype " << static_cast<uint32>(RG_static_infolib[libid].seqtype));
      }
      }
    }
  }

  // templateinfo
  if(RG_static_infolib[libid].insize_from>=0
     || RG_static_infolib[libid].insize_to>=0
     || RG_static_infolib[libid].segmentplacementcode!=SPLACE_UNKNOWN){
    RG_static_infolib[libid].has_templateinfo=true;
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadGroupLib::completeDefaultsForReadGroups()
{
  FUNCSTART("void ReadGroupLib::fillInDefaultParams()");
  for(size_t rglibid=0; rglibid<RG_static_infolib.size(); ++rglibid){
    fillInSensibleDefaults(rglibid);
  }
}




ReadGroupLib::rgid_t ReadGroupLib::searchExactRGMatch(std::string & rgname, uint8 seqtype, int32 isfrom, int32 isto, int8 tpcode, std::string & strainname, bool isbb, bool israil, bool iscer, std::string & seqvecname, std::string & machine_type, std::string & basecaller)
//ReadGroupLib::rgid_t ReadGroupLib::searchExactRGMatch(std::string & rgname, uint8 seqtype, int32 isfrom, int32 isto, std::string tpcode, std::string & strainname, bool isbb, bool israil, bool iscer, std::string & seqvecname, std::string & machine_type, std::string & basecaller)
{
  for(rgid_t libid=1; libid < RG_static_infolib.size(); ++libid){
    /*
    cout << RG_static_infolib.size() << '\t' << static_cast<int64>(libid) << endl;
    cout << "RG_static_infolib[libid].seqtype                 " << (RG_static_infolib[libid].seqtype == seqtype);
    cout << "\nRG_static_infolib[libid].groupname	      " << (RG_static_infolib[libid].groupname == rgname);
    cout << "\nRG_static_infolib[libid].strainname	      " << (RG_static_infolib[libid].strainname == strainname);
    cout << "\nRG_static_infolib[libid].insize_from	      " << (RG_static_infolib[libid].insize_from == isfrom);
    cout << "\nRG_static_infolib[libid].insize_to	      " << (RG_static_infolib[libid].insize_to == isto);
    cout << "\nRG_static_infolib[libid].segmentplacementcode  " << (RG_static_infolib[libid].segmentplacementcode == tpcode);
    cout << "\nRG_static_infolib[libid].is_backbone	      " << (RG_static_infolib[libid].is_backbone == isbb);
    cout << "\nRG_static_infolib[libid].is_rail		      " << (RG_static_infolib[libid].is_rail == israil);
    cout << "\nRG_static_infolib[libid].is_coverageequivalent " << (RG_static_infolib[libid].is_coverageequivalent == iscer);
    cout << "\nRG_static_infolib[libid].seqvecname	      " << (RG_static_infolib[libid].seqvecname == seqvecname);
    cout << "\nRG_static_infolib[libid].machine_type	      " << (RG_static_infolib[libid].machine_type == machine_type);
    cout << "\nRG_static_infolib[libid].basecaller            " << (RG_static_infolib[libid].basecaller == basecaller);
    cout << '\n';
    */

    if(RG_static_infolib[libid].seqtype == seqtype
       && RG_static_infolib[libid].groupname == rgname
       && RG_static_infolib[libid].strainname == strainname
       && RG_static_infolib[libid].insize_from == isfrom
       && RG_static_infolib[libid].insize_to == isto
       && RG_static_infolib[libid].segmentplacementcode == tpcode
       && RG_static_infolib[libid].is_backbone == isbb
       && RG_static_infolib[libid].is_rail == israil
       && RG_static_infolib[libid].is_coverageequivalent == iscer
       && RG_static_infolib[libid].seqvecname == seqvecname
       && RG_static_infolib[libid].machine_type == machine_type
       && RG_static_infolib[libid].basecaller == basecaller
      ){
      return libid;
    }
  }
  return 0;
}



/*************************************************************************
 *
 *
 *
 *
 *
 *************************************************************************/

const std::string & ReadGroupLib::getNameOfSequencingType(uint32 st)
{
  FUNCSTART("const std::string & ReadGroupLib::getNameOfSequencingType(uint32 st)");
  BUGIFTHROW(st>=RG_namesofseqtypes.size(),"Asking for name of unknown sequencing type " << st << " ?");
  FUNCEND();
  return RG_namesofseqtypes[st];
}

const std::string & ReadGroupLib::getShortNameOfSequencingType(uint32 st)
{
  FUNCSTART("const std::string & ReadGroupLib::getShortNameOfSequencingType(uint32 st)");
  BUGIFTHROW(st>=RG_namesofseqtypes.size(),"Asking for name of unknown sequencing type " << st << " ?");
  FUNCEND();
  return RG_shortnamesofseqtypes[st];
}

const std::string & ReadGroupLib::getSAMNameOfSequencingType(uint32 st)
{
  FUNCSTART("const std::string & ReadGroupLib::getSAMNameOfSequencingType(uint32 st)");
  BUGIFTHROW(st>=RG_namesofseqtypes.size(),"Asking for name of unknown sequencing type " << st << " ?");
  FUNCEND();
  return RG_samnamesofseqtypes[st];
}


/*************************************************************************
 *
 * Returns code for sequencing type given in value string
 * or ReadGroupLib::SEQTYPE_END if string was not recognised
 *
 *************************************************************************/

uint8 ReadGroupLib::stringToSeqType(const std::string & value)
{
  if(!value.empty()){
    for(uint8 i=0; i<RG_namesofseqtypes.size(); i++){
      if(value==RG_namesofseqtypes[i]) return i;
    }

    // *sigh* retain compatibility with mira < 3rc3
    if(value=="454GS") return SEQTYPE_454GS20;
    if(value=="454gs") return SEQTYPE_454GS20;

    std::string lcv(value);
    boost::to_lower(lcv);
    for(uint8 i=0; i<RG_lcnamesofseqtypes.size(); i++){
      if(value==RG_lcnamesofseqtypes[i]) return i;
    }
  }

  return getNumSequencingTypes();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void ReadGroupLib::makeStrainIDs()
{
  FUNCSTART("void ReadGroupLib::makeStrainIDs()");

  //std::vector<int32> tidcounter;
  //tidcounter.resize(size(), 0);

  CEBUG("ReadGroupLib::makeStrainIDs dd\n"; debugDumpReadGroupInfo(cout));

  typedef std::unordered_map<std::string, int8> strintmap;
  strintmap sidmap;
  decltype(sidmap.end()) sidI;

  // go through readpool in two rounds, first just looking at Solexa reads,
  //  then at all remaining reads
  // reason: let Solexa have the low strain IDs 0-7, so that they can
  //  be mapped with merge option in the contig

  RG_numstrains=0;
  for(uint32 round=0; round < 2;++round){
    for(uint32 rgid=1; rgid<RG_static_infolib.size();++rgid){
      if(round==0 && RG_static_infolib[rgid].seqtype!=ReadGroupLib::SEQTYPE_SOLEXA) continue;
      if(round>0 && RG_static_infolib[rgid].seqtype==ReadGroupLib::SEQTYPE_SOLEXA) continue;
      sidI=sidmap.find(RG_static_infolib[rgid].strainname);
      if(sidI!=sidmap.end()){
	RG_static_infolib[rgid].strainid=sidI->second;
      }else{
	if(RG_numstrains==127){
	  MIRANOTIFY(Notify::FATAL, "More than 127 strains encountered? Sorry, not possible. Strain " << RG_static_infolib[rgid].strainname << " in readgroup " << RG_static_infolib[rgid].groupname << "\n");
	}
	RG_static_infolib[rgid].strainid=RG_numstrains;
	sidmap[RG_static_infolib[rgid].strainname]=RG_numstrains;
	++RG_numstrains;
      }
    }
  }

  RG_strainids_clean=true;

  CEBUG("Recalc of strain ids yielded " << static_cast<size_t>(RG_numstrains) << " strains." << endl);

  return;
}
//#define CEBUG(bla)


bool ReadGroupLib::getStrainIDOfStrain(const std::string & strainname, int32 & sid)
{
  for(sid=1; sid<static_cast<int32>(RG_static_infolib.size()); ++sid){
    if(RG_static_infolib[sid].strainname==strainname) return true;
  }
  sid=-1;
  return false;
}

const std::string & ReadGroupLib::getStrainOfStrainID(int32 sid)
{
  FUNCSTART("const std::string & ReadGroupLib::getStrainOfStrainID(int32 & sid)");
  for(size_t rgid=1; rgid<RG_static_infolib.size(); ++rgid){
    if(RG_static_infolib[rgid].strainid==sid) return RG_static_infolib[rgid].strainname;
  }

  MIRANOTIFY(Notify::INTERNAL, "Did not find strain id " << sid);
  FUNCEND();
  return RG_emptystring;
}





void ReadGroupLib::setSequencingType(rgid_t libid, uint8 st)
{
  FUNCSTART("void ReadGroupLib::setSequencingType(rgid_t libid, uint8 st)");
  checkLibExistence(libid);
  RG_static_infolib[libid].seqtype=st;

  if(RG_static_infolib[libid].namingscheme==SCHEME_UNKNOWN){
    if(st == ReadGroupLib::SEQTYPE_SANGER){
      setReadNamingScheme(libid,SCHEME_SANGER);
    }else if( st == ReadGroupLib::SEQTYPE_IONTORRENT
	      || st == ReadGroupLib::SEQTYPE_454GS20){
      setReadNamingScheme(libid,SCHEME_SOLEXA);
    }else if(st == ReadGroupLib::SEQTYPE_PACBIOHQ
	     || st == ReadGroupLib::SEQTYPE_PACBIOLQ) {
      setReadNamingScheme(libid,SCHEME_NONE);
    }else if (st==ReadGroupLib::SEQTYPE_SOLEXA){
      setReadNamingScheme(libid,SCHEME_SOLEXA);
    }else if (st==ReadGroupLib::SEQTYPE_TEXT){
      setReadNamingScheme(libid,SCHEME_NONE);
    }else if(st>= ReadGroupLib::SEQTYPE_END)
      MIRANOTIFY(Notify::INTERNAL, "Sequencing technology " << static_cast<int16>(st) << " is unknown to MIRA");
  }else{
    // Scheme already known ... nothing to do?
  }
}



void ReadGroupLib::saveAllReadGroupsAsMAF(std::ostream & ostr)
{
  for(uint32 rgi=1; rgi<RG_static_infolib.size(); ++rgi){
    saveReadGroupAsMAF(rgi,ostr);
  }
}


void ReadGroupLib::saveReadGroupAsMAF(uint32 rgi, std::ostream & ostr)
{
  FUNCSTART("void ReadGroupLib::saveReadGroupAsMAF(uint32 rgi, std::ostream & ostr)");
  BUGIFTHROW(rgi>=RG_static_infolib.size(),"rgi (" << rgi << ") >= RG_static_infolib.size() " << RG_static_infolib.size() << ") ???");
  if(!RG_static_infolib[rgi].wassaved){
    RG_static_infolib[rgi].wassaved=true;
    dumpReadGroupAsMAF(rgi,ostr);
  }

}


void ReadGroupLib::dumpAllReadGroupsAsMAF(std::ostream & ostr)
{
  for(uint32 rgi=1; rgi<RG_static_infolib.size(); ++rgi){
    dumpReadGroupAsMAF(rgi,ostr);
  }
}



void ReadGroupLib::dumpReadGroupAsMAF(uint32 rgi, std::ostream & ostr)
{
  FUNCSTART("void ReadGroupLib::dumpReadGroupAsMAF(uint32 rgi, std::ostream & ostr)");
  BUGIFTHROW(rgi>=RG_static_infolib.size(),"rgi (" << rgi << ") >= RG_static_infolib.size() " << RG_static_infolib.size() << ") ???");

  ostr << "@ReadGroup\n";
  if(RG_static_infolib[rgi].groupname.size()){
    ostr << "@RG\tname\t" << RG_static_infolib[rgi].groupname << '\n';
  }
  ostr << "@RG\tID\t" << rgi << '\n';
  ostr << "@RG\ttechnology\t" << getNameOfSequencingType(RG_static_infolib[rgi].seqtype) << '\n';
  if(RG_static_infolib[rgi].strainname.size()){
    ostr << "@RG\tstrainname\t" << RG_static_infolib[rgi].strainname << '\n';
  }
  if(RG_static_infolib[rgi].insize_from>=0 || RG_static_infolib[rgi].insize_to>=0 ){
    ostr << "@RG\ttemplatesize\t" << RG_static_infolib[rgi].insize_from
	 << '\t' << RG_static_infolib[rgi].insize_to << '\n';
  }
  if(RG_static_infolib[rgi].segmentplacementcode != SPLACE_UNKNOWN){
    ostr << "@RG\tsegmentplacement\t" << getNameOfSegmentplacement(RG_static_infolib[rgi].segmentplacementcode) << '\n';
  }
  ostr << "@RG\tsegmentnaming\t" << getNameOfNamingScheme(RG_static_infolib[rgi].namingscheme) << '\n';

  if(!RG_static_infolib[rgi].machine_type.empty()){
    ostr << "@RG\tmachinetype\t" << RG_static_infolib[rgi].machine_type << '\n';
  }
  if(!RG_static_infolib[rgi].basecaller.empty()){
    ostr << "@RG\tbasecaller\t" << RG_static_infolib[rgi].basecaller << '\n';
  }
  if(!RG_static_infolib[rgi].dye.empty()){
    ostr << "@RG\tdye\t" << RG_static_infolib[rgi].dye << '\n';
  }
  if(!RG_static_infolib[rgi].primer.empty()){
    ostr << "@RG\tprimer\t" << RG_static_infolib[rgi].primer << '\n';
  }
  if(!RG_static_infolib[rgi].clonevecname.empty()){
    ostr << "@RG\tclonevecname\t" << RG_static_infolib[rgi].clonevecname << '\n';
  }
  if(!RG_static_infolib[rgi].seqvecname.empty()){
    ostr << "@RG\tseqvecname\t" << RG_static_infolib[rgi].seqvecname << '\n';
  }
  if(!RG_static_infolib[rgi].adaptorleft.empty()){
    ostr << "@RG\tadaptorleft\t" << RG_static_infolib[rgi].adaptorleft << '\n';
  }
  if(!RG_static_infolib[rgi].adaptorright.empty()){
    ostr << "@RG\tadaptorright\t" << RG_static_infolib[rgi].adaptorright << '\n';
  }
  if(!RG_static_infolib[rgi].adaptorsplit.empty()){
    ostr << "@RG\tadaptorsplit\t" << RG_static_infolib[rgi].adaptorsplit << '\n';
  }
  if(!RG_static_infolib[rgi].datadir.empty()){
    ostr << "@RG\tdatadir\t" << RG_static_infolib[rgi].datadir << '\n';
  }
  if(!RG_static_infolib[rgi].datafile.empty()){
    ostr << "@RG\tdatafile\t" << RG_static_infolib[rgi].datafile << '\n';
  }

  if(RG_static_infolib[rgi].is_backbone){
    ostr << "@RG\tisbackbone\n";
  }
  if(RG_static_infolib[rgi].is_rail){
    ostr << "@RG\tisrail\n";
  }
  if(RG_static_infolib[rgi].is_coverageequivalent){
    ostr << "@RG\tiscoverageequivalent\n";
  }
  ostr << "@EndReadGroup\n";
}

void ReadGroupLib::dumpReadGroupAsSAM(uint32 rgi, std::ostream & ostr)
{
  FUNCSTART("void ReadGroupLib::dumpReadGroupAsSAM(uint32 rgi, std::ostream & ostr)");
  BUGIFTHROW(rgi>=RG_static_infolib.size(),"rgi (" << rgi << ") >= RG_static_infolib.size() " << RG_static_infolib.size() << ") ???");

  ostr << "@RG\tID:" << rgi
       << "\tPL:" << getSAMNameOfSequencingType(RG_static_infolib[rgi].seqtype);
  if(RG_static_infolib[rgi].groupname.size()){
    ostr << "\tLB:" << RG_static_infolib[rgi].groupname;
  }
  if(RG_static_infolib[rgi].strainname.size()){
    ostr << "\tSM:" << RG_static_infolib[rgi].strainname;
  }
  if(RG_static_infolib[rgi].insize_from>=0 || RG_static_infolib[rgi].insize_to>=0 ){
    ostr << "\tPI:" << RG_static_infolib[rgi].insize_to-RG_static_infolib[rgi].insize_from;
  }
  ostr << endl;
}


void ReadGroupLib::ReadGroupID::setSequencingType(std::string t)
{
  boost::to_lower(t);
  ReadGroupLib::setSequencingType(rgid_id,ReadGroupLib::stringToSeqType(t));
}

void ReadGroupLib::debugDumpReadGroupInfo(std::ostream & ostr)
{
  for(uint32 rgi=0; rgi<RG_static_infolib.size(); ++rgi){
    ostr << "\n\nRGI: " << rgi << "\t" << ReadGroupLib::getReadGroupID(rgi) << endl;
  }
}


/*************************************************************************
 *
 *  Accepts: FR / RF / ...
 *   or graphical representation like:
 *          ---------> <----------
 *          ---------> ---------->
 *          1---------> 2---------->
 *          2---------> 1---------->
 *  returns:
 *    true / false to signal whether code recognise (empty string is also OK)
 *    placementcode
 *
 *************************************************************************/


bool ReadGroupLib::parseSegmentPlacement(const std::string & sps, int8 & placementcode)
{
  FUNCSTART("bool ReadGroupLib::parseSegmentPlacement(const std::string & sps, int8 & placementcode)");

  placementcode=SPLACE_UNKNOWN;

  if(!sps.empty()){
    std::string shortsps;
    for(auto cs : sps){
      if(!isspace(cs)){
	if(cs=='>'){
	  shortsps+="FORWARD";
	}else if(cs=='<'){
	  shortsps+="REVERSE";
	}else if(cs=='='){
	  shortsps+="SAMEDIR";
	}else if(cs=='?'){
	  shortsps+="UNKNOWN";
	}else if(cs=='-'){
	  // do nothing, munch the character to allow for "--->" etc.
	}else{
	  shortsps+=toupper(cs);
	}
      }
    }
    if(shortsps=="SF"
       || shortsps=="SAMEDIRFORWARD"
       || shortsps=="FORWARDFORWARD"
       || shortsps=="1FORWARD2FORWARD"
       || shortsps=="SAMEDIRECTIONFORWARD"
       || shortsps=="LEFTIES"
       || shortsps=="LEFTIE"
       || shortsps=="LEFTY") {
      placementcode=SPLACE_SF;
    }else if(shortsps=="SB"
	     || shortsps=="2FORWARD1FORWARD"
	     || shortsps=="SAMEDIRREVERSE"
	     || shortsps=="SAMEDIRECTIONREVERSE"
	     || shortsps=="SAMEDIRBACKWARD"
	     || shortsps=="SAMEDIRECTIONBACKWARD"
	     || shortsps=="RIGHTIES"
	     || shortsps=="RIGHTIE"
	     || shortsps=="RIGHTY") {
      placementcode=SPLACE_SB;
    }else if(shortsps=="SAMEDIR"
	     || shortsps=="SU"
	     || shortsps=="SAMEDIRECTION"
	     || shortsps=="SAMEDIRUNKNOWN"
	     || shortsps=="SAMEDIRECTIONUNKNOWN") {
      placementcode=SPLACE_SU;
    }else if(shortsps=="FORWARDREVERSE"
	     || shortsps=="FR"
	     || shortsps=="INNIES"
	     || shortsps=="INNIE"
	     || shortsps=="INNY") {
      placementcode=SPLACE_FR;
    }else if(shortsps=="REVERSEFORWARD"
	     || shortsps=="RF"
	     || shortsps=="OUTIES"
	     || shortsps=="OUTIE"
	     || shortsps=="OUTY") {
      placementcode=SPLACE_RF;
    }else if(shortsps=="UNKNOWN"){
      // intentionally do nothing, SPLACE_UNKNOWN is default anyway
    }else{
      return false;
    }
  }
  return true;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool ReadGroupLib::parseSegmentNaming(const std::string & sns, uint8 & namingscheme)
{
  FUNCSTART("bool ReadGroupLib::parseSegmentNaming(const std::string & sns, uint8 & namingscheme)");

  namingscheme=SCHEME_UNKNOWN;

  std::string tmpsns(sns);
  boost::to_lower(tmpsns);

  if(tmpsns=="solexa"){
    namingscheme=SCHEME_SOLEXA;
  }else if(tmpsns=="sanger"){
    namingscheme=SCHEME_SANGER;
  }else if(tmpsns=="stlouis"){
    namingscheme=SCHEME_STLOUIS;
  }else if(tmpsns=="tigr"){
    namingscheme=SCHEME_TIGR;
  }else if(tmpsns=="fr"){
    namingscheme=SCHEME_FR;
  }else if(tmpsns=="sra"){
    namingscheme=SCHEME_SRARAW;
  }else if(tmpsns=="unknown"){
    namingscheme=SCHEME_UNKNOWN;
  }else if(tmpsns=="none"){
    namingscheme=SCHEME_NONE;
  }else{
    return false;
  }

  FUNCEND();
  return true;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool ReadGroupLib::hasLibWithSeqType(uint8 seqtype)
{
  for(uint32 rgi=1; rgi<RG_static_infolib.size(); ++rgi){
    if(RG_static_infolib[rgi].seqtype==seqtype) return true;
  }
  return false;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

const std::string & ReadGroupLib::getNameOfNamingScheme(uint8 nscode)
{
  FUNCSTART("const std::string & ReadGroupLib::getNameOfNamingScheme(uint8 nscode)");
  BUGIFTHROW(nscode >= RG_namesofnamingschemes.size(),"illegal nscode " << static_cast<uint16>(nscode));
  return RG_namesofnamingschemes[nscode];
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

const std::string & ReadGroupLib::getNameOfSegmentplacement(int8 spcode)
{
  FUNCSTART("const std::string & ReadGroupLib::getNameOfSequencingType(int8 spcode)");
  auto idx=static_cast<size_t>(spcode-1-SPLACE_UNUSED_LOW);
  BUGIFTHROW(idx >= RG_namesofsegmentplacements.size(),"illegal spcode " << static_cast<int16>(spcode) << " " << idx);
  return RG_namesofsegmentplacements[idx];
}
