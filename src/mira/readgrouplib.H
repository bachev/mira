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

#ifndef _bas_readgrouplib_h_
#define _bas_readgrouplib_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

class ReadGroupLib
{
public:
#include "mira/readseqtypes.H"

  typedef uint8 rgid_t;

  static const std::vector<std::string> RG_namesofnamingschemes;
  static const std::vector<std::string> RG_namesofsegmentplacements;
  static const std::vector<std::string> RG_namesofseqtypes;
  static const std::vector<std::string> RG_lcnamesofseqtypes;
  static const std::vector<std::string> RG_shortnamesofseqtypes;
  static const std::vector<std::string> RG_samnamesofseqtypes;

  static const std::string RG_emptystring;

  class ReadGroupID
  {
  private:
    rgid_t rgid_id;

  public:
    ReadGroupID() : rgid_id(0) {};
    ReadGroupID(size_t rgid) : rgid_id(static_cast<rgid_t>(rgid)) {};

    friend std::ostream & operator<<(std::ostream &ostr, const ReadGroupID & r){
      ostr << "RGID: " << static_cast<uint16>(r.rgid_id) << "\n";
      r.dumpDebug(ostr);
      ostr << std::endl;
      return ostr;
    }

    inline bool operator==(const ReadGroupID & other) const {return rgid_id==other.rgid_id;}
    inline bool operator!=(const ReadGroupID & other) const {return rgid_id!=other.rgid_id;}
    inline bool operator<(const ReadGroupID & other) const {return rgid_id<other.rgid_id;}

    inline bool isDefaultNonValidReadGroupID() const {return rgid_id==0;}
    inline void resetLibId() {rgid_id=0;};
    inline uint32 getLibId() const {return rgid_id;};
    void checkValidity() {ReadGroupLib::checkLibValidity(rgid_id);}
    void fillInSensibleDefaults() {ReadGroupLib::fillInSensibleDefaults(rgid_id);}

    inline const std::string & getBaseCaller() const  {return ReadGroupLib::getBaseCaller(rgid_id);}
    inline int8  getBuildDirection() const {return ReadGroupLib::getBuildDirection(rgid_id);}
    inline const std::string & getCloneVecName() const  {return ReadGroupLib::getCloneVecName(rgid_id);}
    inline uint8 getDefaultQual() const {return ReadGroupLib::getDefaultQual(rgid_id);}
    inline bool hasUserDefaultQual() const {return ReadGroupLib::hasUserDefaultQual(rgid_id);}
    inline const std::string & getDataDir() const  {return ReadGroupLib::getDataDir(rgid_id);}
    inline const std::string & getDataFile() const  {return ReadGroupLib::getDataFile(rgid_id);}
    inline const std::string & getDye() const  {return ReadGroupLib::getDye(rgid_id);}
    inline const std::string & getGroupName() const  {return ReadGroupLib::getGroupName(rgid_id);}
    inline int32 getInsizeFrom() const { return ReadGroupLib::getInsizeFrom(rgid_id);}
    inline int32 getInsizeTo() const { return ReadGroupLib::getInsizeTo(rgid_id);}
    inline const std::string & getMachineType() const  {return ReadGroupLib::getMachineType(rgid_id);}
    inline const std::string & getPrimer() const  {return ReadGroupLib::getPrimer(rgid_id);}
    inline uint8 getReadNamingScheme() const {return ReadGroupLib::getReadNamingScheme(rgid_id);}
    inline uint8 getSequencingType() const {return ReadGroupLib::getSequencingType(rgid_id);}
    inline const std::string & getSeqVecName() const  {return ReadGroupLib::getSeqVecName(rgid_id);}
    inline int8  getStrainID() const {return ReadGroupLib::getStrainID(rgid_id);}
    inline const std::string & getStrainName() const  {return ReadGroupLib::getStrainName(rgid_id);}
    inline const std::string & getSegmentPlacement() const  {return ReadGroupLib::getSegmentPlacement(rgid_id);}
    inline int8 getSegmentPlacementCode() const  {return ReadGroupLib::getSegmentPlacementCode(rgid_id);}

    inline bool  hasTemplateInfo() const {return ReadGroupLib::hasTemplateInfo(rgid_id);}
    inline bool  wantTemplateInfoEstimate() const {return ReadGroupLib::wantTemplateInfoEstimate(rgid_id);}
    inline bool  wantSegmentPlacementEstimate() const {return ReadGroupLib::wantSegmentPlacementEstimate(rgid_id);}
    inline bool  wantTemplateSizeEstimate() const {return ReadGroupLib::wantTemplateSizeEstimate(rgid_id);}
    inline bool  expectsReadPairs() const {return ReadGroupLib::expectsReadPairs(rgid_id);}
    inline bool  getSPInfoOnly() const {return ReadGroupLib::getSPInfoOnly(rgid_id);}
    inline bool  getTSInfoOnly() const {return ReadGroupLib::getTSInfoOnly(rgid_id);}

    inline bool  wantStatisticsCalc() const {return ReadGroupLib::wantStatisticsCalc(rgid_id);}
    inline void  setStatisticsCalc(bool b) {ReadGroupLib::setStatisticsCalc(rgid_id,b);}

    inline bool  wantUseReadNameFromComment() const {return ReadGroupLib::wantUseReadNameFromComment(rgid_id);}
    inline void  setUseReadNameFromComment(bool b) {ReadGroupLib::setUseReadNameFromComment(rgid_id,b);}

    inline bool  isBackbone() const {return ReadGroupLib::isBackbone(rgid_id);}
    inline bool  isRail() const {return ReadGroupLib::isRail(rgid_id);}
    inline bool  isCoverageEquivalentRead() const {return ReadGroupLib::isCoverageEquivalentRead(rgid_id);}

    inline void setWantSegmentPlacementEstimate(bool b) {ReadGroupLib::setWantSegmentPlacementEstimate(rgid_id,b);}
    inline void setWantTemplateSizeEstimate(bool b) {ReadGroupLib::setWantTemplateSizeEstimate(rgid_id,b);}
    inline void setExpectsReadPairs(bool b) {ReadGroupLib::setExpectsReadPairs(rgid_id,b);}

    inline void setBackbone(bool b) {ReadGroupLib::setBackbone(rgid_id,b);}
    inline void setBaseCaller(const std::string & s) {ReadGroupLib::setBaseCaller(rgid_id,s);}
    inline void setCoverageEquivalentRead(bool b) {ReadGroupLib::setCoverageEquivalentRead(rgid_id,b);}
    inline void setCloneVecName(const std::string & s) {ReadGroupLib::setCloneVecName(rgid_id,s);}
    inline void setDataDir(const std::string & s) {ReadGroupLib::setDataDir(rgid_id,s);}
    inline void setDataFile(const std::string & s) {ReadGroupLib::setDataFile(rgid_id,s);}
    inline void setDefaultQual(base_quality_t i) {ReadGroupLib::setDefaultQual(rgid_id,i);}
    inline void setDye(const std::string & s) {ReadGroupLib::setDye(rgid_id,s);}
    inline void setGroupName(const std::string & s){ReadGroupLib::setGroupName(rgid_id,s);}
    inline void setInsizeFrom(int32 i) {ReadGroupLib::setInsizeFrom(rgid_id,i);}
    inline void setInsizeTo(int32 i) {ReadGroupLib::setInsizeTo(rgid_id,i);}
    inline void setMachineType(const std::string & s) {ReadGroupLib::setMachineType(rgid_id,s);}
    inline void setPrimer(const std::string & s) {ReadGroupLib::setPrimer(rgid_id,s);}
    inline void setRail(bool b) {ReadGroupLib::setRail(rgid_id,b);}
    inline void setReadNamingScheme(uint8 t) {ReadGroupLib::setReadNamingScheme(rgid_id,t);}
    inline void setSequencingType(uint8 t) {ReadGroupLib::setSequencingType(rgid_id,t);}
    void setSequencingType(std::string t);
    inline void setSeqVecName(const std::string & s) {ReadGroupLib::setSeqVecName(rgid_id,s);}
    inline void setStrainID(int8 i) {ReadGroupLib::setStrainID(rgid_id,i);}
    inline void setStrainName(const std::string & s) {ReadGroupLib::setStrainName(rgid_id,s);}
    inline bool setSegmentPlacement(const std::string & s) {return ReadGroupLib::setSegmentPlacement(rgid_id,s);}
    inline void setSegmentPlacementCode(const int8 code) {ReadGroupLib::setSegmentPlacementCode(rgid_id,code);}
    inline bool setSegmentNaming(const std::string & s) {return ReadGroupLib::setSegmentNaming(rgid_id,s);}
    inline void setSPInfoOnly(bool b) {ReadGroupLib::setSPInfoOnly(rgid_id,b);}
    inline void setTSInfoOnly(bool b) {ReadGroupLib::setTSInfoOnly(rgid_id,b);}

    inline const std::string & getNameOfSequencingType() const {return ReadGroupLib::getNameOfSequencingType(RG_static_infolib[rgid_id].seqtype);};
    inline const std::string & getShortNameOfSequencingType() const {return ReadGroupLib::getShortNameOfSequencingType(RG_static_infolib[rgid_id].seqtype);};

    inline const std::vector<std::string> & getReadRenamePrefixList() const {return ReadGroupLib::getReadRenamePrefixList(rgid_id);};
    inline const void addReadRenamePrefix(const std::string & pre, const std::string & ren) {ReadGroupLib::addReadRenamePrefix(rgid_id,pre,ren);};

    std::ostream & dumpDebug(std::ostream & ostr) const {return ostr << ReadGroupLib::RG_static_infolib[rgid_id];}

  };


  //Variables
private:
  static const bool RG_initialisedstatics;

  struct rginfo_t {
    std::string groupname;
    std::string strainname;
    std::string segmentplacement;
    int32 insize_from;
    int32 insize_to;
    uint8 seqtype;
    uint8 namingscheme;

    base_quality_t defaultqual;

    std::string machine_type;
    std::string basecaller;
    std::string dye;
    std::string primer;
    std::string clonevecname;
    std::string seqvecname;

    std::string adaptorleft;
    std::string adaptorright;
    std::string adaptorsplit;

    std::string datadir;
    std::string datafile;

    std::vector<std::string> readrenameprefixlist;

    bool is_backbone;
    bool is_rail;
    bool is_coverageequivalent;

    //*****************************************************************
    // Maybe create new class "rgparams" which would also hold
    //  all clipping etc. parameters which are specific to readgroups?
    // This is again a long term TODO though
    bool expectsreadpairs;
    bool sp_wantestimate;  // segment placement: wants MIRA to estimate value
    bool ts_wantestimate;  // template size: wants MIRA to estimate value
    bool sp_infoonly;  // segment placement: info only, do not use in assembly
    bool ts_infoonly;  // template size: info only, do not use in assembly

    bool hack_calcstats;  // MIRA can use Rg to calc statistics
    //*****************************************************************

    bool sn_rollcommenttoname;  // roll first comment word to be name of the read

    // computed, then query by outside
    bool has_templateinfo;
    bool has_userdefaultqual;
    int8 segmentplacementcode;

    int8  strainid;

    bool wassaved;   // used by readgrouplib to show which Rg was already saved

  public:
    inline rginfo_t() : insize_from(-1),
			insize_to(-1),
			seqtype(SEQTYPE_END),
			namingscheme(SCHEME_UNKNOWN),
			defaultqual(101),
			strainid(0),
			wassaved(false),
			is_backbone(false),
			is_rail(false),
			is_coverageequivalent(false),
			has_templateinfo(false),
			has_userdefaultqual(false),
			expectsreadpairs(false),
			sp_wantestimate(true),
			ts_wantestimate(true),
			sp_infoonly(false),
			ts_infoonly(false),
			hack_calcstats(true),
			sn_rollcommenttoname(false),
			segmentplacementcode(SPLACE_UNKNOWN) {};

    friend std::ostream & operator<<(std::ostream &ostr, const rginfo_t & rgi){
      ostr << "RGN: " << rgi.groupname << "\tSN: " << rgi.strainname
	   << "\nSP: " << rgi.segmentplacement << "\tSPio: " << rgi.sp_infoonly
	   << "\tSPC: " << static_cast<int32>(rgi.segmentplacementcode)
	   << "\tIF: " << rgi.insize_from << "\tIT: " << rgi.insize_to << "\tTSio: " << rgi.ts_infoonly
	   << "\nST: " << static_cast<uint16>(rgi.seqtype);
      ostr.flush();
      ostr << " (" << getNameOfSequencingType(rgi.seqtype)
	   << ")\tnamschem: " << static_cast<uint16>(rgi.namingscheme) << "\tSID: " << static_cast<uint16>(rgi.strainid)
	   << "\nDQ: " << static_cast<uint16>(rgi.defaultqual)
	   << "\nBB: " << rgi.is_backbone << "\tRail: " << rgi.is_rail << "\tCER: " << rgi.is_coverageequivalent
	   << std::endl;
      return ostr;
    }
  };
  // sooner or later, clipping options should move here

  static std::vector<rginfo_t> RG_static_infolib;

  static bool RG_strainids_clean;
  static int8 RG_numstrains;

  //Functions
private:
  void foolCompiler();
  static bool staticInitialiser();
  void init();
  static void zeroVars();

  static void checkLibValidity(rgid_t libid);
  static void fillInSensibleDefaults(rgid_t libid);

  inline static void checkLibExistence(rgid_t libid) {
    FUNCSTART("void checkLibExistence(gid_t libid)");
    BUGIFTHROW(libid > RG_static_infolib.size(), "libid " << static_cast<uint32>(libid) << " > RG_static_infolib.size() " << RG_static_infolib.size() << " ???");
    FUNCEND();
  }

  inline static const std::string & getBaseCaller(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].basecaller;}
  inline static int8 getBuildDirection(rgid_t libid){
    checkLibExistence(libid);
    if(RG_static_infolib[libid].segmentplacementcode>=0) return 1;
    return -1;
  }
  inline static const std::string & getCloneVecName(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].clonevecname;}
  inline static const std::string & getDataDir(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].datadir;}
  inline static const std::string & getDataFile(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].datafile;}
  inline static base_quality_t getDefaultQual(rgid_t libid) {checkLibExistence(libid); return RG_static_infolib[libid].defaultqual;}
  inline static bool hasUserDefaultQual(rgid_t libid) {checkLibExistence(libid); return RG_static_infolib[libid].has_userdefaultqual;}
  inline static const std::string & getDye(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].dye;}
  inline static const std::string & getGroupName(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].groupname;}
  inline static int32 getInsizeFrom(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].insize_from;}
  inline static int32 getInsizeTo(rgid_t libid) {checkLibExistence(libid); return RG_static_infolib[libid].insize_to;}
  inline static const std::string & getMachineType(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].machine_type;}
  inline static const std::string & getPrimer(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].primer;}
  inline static uint8 getReadNamingScheme(rgid_t libid) {checkLibExistence(libid); return RG_static_infolib[libid].namingscheme;}
  inline static std::vector<std::string> & getReadRenamePrefixList(rgid_t libid) {checkLibExistence(libid); return RG_static_infolib[libid].readrenameprefixlist;}
  inline static const std::string & getStrainName(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].strainname;}
  inline static uint8 getSequencingType(rgid_t libid) {checkLibExistence(libid); return RG_static_infolib[libid].seqtype;}
  inline static const std::string & getSeqVecName(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].seqvecname;}
  inline static int8 getStrainID(rgid_t libid){
    checkLibExistence(libid);
    //std::cout << "RG_strainids_clean: " << RG_strainids_clean << endl;
    if(!RG_strainids_clean) makeStrainIDs();
    //std::cout << "RG_static_infolib[" << static_cast<int32>(libid) << "] = " << static_cast<int32>(RG_static_infolib[libid].strainid) << endl;
    return RG_static_infolib[libid].strainid;
  }
  inline static const std::string & getSegmentPlacement(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].segmentplacement;}
  inline static int8 getSegmentPlacementCode(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].segmentplacementcode;}

  inline static bool hasTemplateInfo(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].has_templateinfo;}
  inline static bool wantTemplateInfoEstimate(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].sp_wantestimate | RG_static_infolib[libid].ts_wantestimate;}
  inline static bool wantSegmentPlacementEstimate(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].sp_wantestimate;}
  inline static bool wantTemplateSizeEstimate(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].ts_wantestimate;}
  inline static bool expectsReadPairs(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].expectsreadpairs;}
  inline static bool getSPInfoOnly(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].sp_infoonly;}
  inline static bool getTSInfoOnly(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].ts_infoonly;}

  inline static bool isBackbone(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].is_backbone;}
  inline static bool isRail(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].is_rail;}
  inline static bool isCoverageEquivalentRead(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].is_coverageequivalent;}

  inline static bool  wantUseReadNameFromComment(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].sn_rollcommenttoname;}
  inline static void  setUseReadNameFromComment(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].sn_rollcommenttoname=v;}

  inline static bool  wantStatisticsCalc(rgid_t libid){checkLibExistence(libid); return RG_static_infolib[libid].hack_calcstats;}
  inline static void  setStatisticsCalc(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].hack_calcstats=v;}

  inline static void setWantSegmentPlacementEstimate(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].sp_wantestimate=v;}
  inline static void setWantTemplateSizeEstimate(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].ts_wantestimate=v;}
  inline static void setExpectsReadPairs(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].expectsreadpairs=v;}
  inline static void addReadRenamePrefix(rgid_t libid, const std::string & pre, const std::string & ren) {
    checkLibExistence(libid);
    RG_static_infolib[libid].readrenameprefixlist.push_back(pre);
    RG_static_infolib[libid].readrenameprefixlist.push_back(ren);
  }
  inline static void setBackbone(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].is_backbone=v;}
  inline static void setBaseCaller(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].basecaller=s;}
  inline static void setCloneVecName(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].clonevecname=s;}
  inline static void setCoverageEquivalentRead(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].is_coverageequivalent=v;}
  inline static void setDataDir(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].datadir=s;}
  inline static void setDataFile(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].datafile=s;}
  inline static void setDefaultQual(rgid_t libid, int8 i) {checkLibExistence(libid); RG_static_infolib[libid].defaultqual=i;}
  inline static void setDye(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].dye=s;}
  inline static void setGroupName(rgid_t libid, const std::string & s){checkLibExistence(libid); RG_static_infolib[libid].groupname=s;}
  inline static void setInsizeFrom(rgid_t libid, int32 i) {
    checkLibExistence(libid); RG_static_infolib[libid].insize_from=i;
    if(i>=0){
      RG_static_infolib[libid].has_templateinfo=true;
    }else if(RG_static_infolib[libid].insize_to<0
	     && RG_static_infolib[libid].segmentplacementcode==SPLACE_UNKNOWN){
      RG_static_infolib[libid].has_templateinfo=false;
    }
  }
  inline static void setInsizeTo(rgid_t libid, int32 i) {
    checkLibExistence(libid); RG_static_infolib[libid].insize_to=i;
    if(i>=0){
      RG_static_infolib[libid].has_templateinfo=true;
    }else if(RG_static_infolib[libid].insize_from<0
	     && RG_static_infolib[libid].segmentplacementcode==SPLACE_UNKNOWN){
      RG_static_infolib[libid].has_templateinfo=false;
    }
  }
  inline static void setMachineType(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].machine_type=s;}
  inline static void setPrimer(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].primer=s;}
  inline static void setRail(rgid_t libid, bool v) {checkLibExistence(libid); RG_static_infolib[libid].is_rail=v;}
  inline static void setReadNamingScheme(rgid_t libid, uint8 t) {checkLibExistence(libid); RG_static_infolib[libid].namingscheme=t;}
  static void setSequencingType(rgid_t libid, uint8 st);
  inline static void setSeqVecName(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].seqvecname=s;}
  inline static void setStrainID(rgid_t libid, int8 i) {checkLibExistence(libid); RG_static_infolib[libid].strainid=i;}
  inline static void setStrainName(rgid_t libid, const std::string & s) {checkLibExistence(libid); RG_static_infolib[libid].strainname=s; RG_strainids_clean=false;}

  static void setSegmentPlacementCode(rgid_t libid, const int8 code){
    FUNCSTART("static bool ReadGroupLib::setSegmentPlacement(rgid_t libid, const int8 code)");
    BUGIFTHROW(code<=SPLACE_UNUSED_LOW || code>=SPLACE_UNUSED_HIGH,"code " << static_cast<int16>(code) << " is invalid.");
    RG_static_infolib[libid].segmentplacementcode=code;
    if(RG_static_infolib[libid].segmentplacementcode!=SPLACE_UNKNOWN){
      RG_static_infolib[libid].has_templateinfo=true;
    }else if(RG_static_infolib[libid].insize_from<0
	     && RG_static_infolib[libid].insize_to<0){
      RG_static_infolib[libid].has_templateinfo=false;
    }
  }
  static bool setSegmentPlacement(rgid_t libid, const std::string & sps){
    RG_static_infolib[libid].segmentplacement=sps;
    auto x=parseSegmentPlacement(sps,RG_static_infolib[libid].segmentplacementcode);
    if(RG_static_infolib[libid].segmentplacementcode!=SPLACE_UNKNOWN){
      RG_static_infolib[libid].has_templateinfo=true;
    }else if(RG_static_infolib[libid].insize_from<0
	     && RG_static_infolib[libid].insize_to<0){
      RG_static_infolib[libid].has_templateinfo=false;
    }
    return x;
  }

  static bool setSegmentNaming(rgid_t libid, const std::string & sns){
    return parseSegmentNaming(sns,RG_static_infolib[libid].namingscheme);
  }

  inline static void setSPInfoOnly(rgid_t libid, bool b) {checkLibExistence(libid); RG_static_infolib[libid].sp_infoonly=b;}
  inline static void setTSInfoOnly(rgid_t libid, bool b) {checkLibExistence(libid); RG_static_infolib[libid].ts_infoonly=b;}


  static void makeStrainIDs();

public:
  ReadGroupLib();
  ReadGroupLib(ReadGroupLib const &other);
  ~ReadGroupLib();

  ReadGroupLib const & operator=(ReadGroupLib const & other);
  friend std::ostream & operator<<(std::ostream &ostr, ReadGroupLib const &rgl);

  static void discard();

  static rgid_t setReadGroupInfo2(rgid_t libid,
				  std::string & rgname,
				  uint8 seqtype,
				  int32 isfrom,
				  int32 isto,
				  int8  dir);

  static ReadGroupID newReadGroup();
  static inline size_t getNumReadGroups() { return RG_static_infolib.size();}
  static inline ReadGroupID getReadGroupID(rgid_t libid) {
    checkLibExistence(libid);
    return ReadGroupID(libid);
  }

  static inline uint8 getNumSequencingTypes() {return SEQTYPE_END;}
  static const std::string & getNameOfSequencingType(uint32 st);
  static const std::string & getShortNameOfSequencingType(uint32 st);
  static const std::string & getSAMNameOfSequencingType(uint32 st);
  static uint8 stringToSeqType(const std::string & value);

  static bool parseSegmentPlacement(const std::string & sps, int8 & tpcode);
  static bool parseSegmentNaming(const std::string & sns, uint8 & schemecode);

  static const std::string & getNameOfNamingScheme(uint8 ns);
  static const std::string & getNameOfSegmentplacement(int8 sp);

  static inline int8 getNumOfStrains() {
    if(!RG_strainids_clean) makeStrainIDs();
    return RG_numstrains;
  }
  static const std::string & getStrainOfStrainID(int32 sid);
  static bool getStrainIDOfStrain(const std::string & strainname, int32 & sid);

  static void dumpStrainIDSummary() {std::cout << "TODO: Like Readpool: strain x has y reads\n";}

  static bool hasLibWithSeqType(uint8 seqtype);

  static rgid_t searchExactRGMatch(
    std::string & rgname,
    uint8 seqtype,
    int32 isfrom,
    int32 isto,
    int8 segmentplacementcode,
    std::string & strainname,
    bool isbb,
    bool israil,
    bool iscer,
    std::string & seqvecname,
    std::string & machine_type,
    std::string & basecaller);

  static void completeDefaultsForReadGroups();

  static void dumpAllReadGroupsAsMAF(std::ostream & ostr);
  static void dumpReadGroupAsMAF(uint32 rgi, std::ostream & ostr);
  static void saveAllReadGroupsAsMAF(std::ostream & ostr);
  static void saveReadGroupAsMAF(uint32 rgi, std::ostream & ostr);
  static void resetSaveStatus() {for(auto & rg : RG_static_infolib) {rg.wassaved=false;} }

  static void dumpReadGroupAsSAM(uint32 rgi, std::ostream & ostr);

  static void debugDumpReadGroupInfo(std::ostream & ostr);
};


#endif
