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

#ifndef _mira_parameters_H
#define _mira_parameters_H

#include "stdinc/defines.H"

#include <iostream>
#include <fstream>

#include <sstream>
#include <algorithm>                     // STL min max template

// for strlen()
#include <cstring>



#include "errorhandling/errorhandling.H"

#include "mira/structs.H"

#undef yyFlexLexer
#define yyFlexLexer MPFlexLexer
#include <FlexLexer.h>

class MIRAParameters
{
private:

  enum qsj_jobareas_t {JA_QUALITY, JA_TYPE, JA_METHOD, JA_TECH};

  static bool MP_errorinparams;
  static std::string MP_currentparametersection;
  static std::vector<std::string>        MP_loadfilename;

  static std::string MP_binpath;
  static std::string MP_bindir;
  static std::string MP_homedir;
  static std::string MP_sharedir;
  static std::string MP_mhslibdir;

  special_parameters    mp_special_params;
  nagandwarn_parameters mp_nagandwarn_params;
  skim_parameters       mp_skim_params;
  hashstatistics_parameters  mp_hashstatistics_params;
  align_parameters      mp_align_params;
  finalmapping_parameters  mp_finalmap_params;
  dynamic_parameters    mp_dynamic_params;
  assembly_parameters   mp_assembly_params;
  pathfinder_parameters mp_pathfinder_params;
  contig_parameters     mp_contig_params;
  directory_parameters  mp_directory_params;
  file_parameters       mp_file_params;
  edit_parameters       mp_edit_params;

  bool MP_parsedsomesettings;
  bool MP_jobtechsettings;

  //Functions
private:
  void foolCompiler();

  template <class T>
  void checkMin(const char * parameter, T & var, const T & minval, const T & defaultval, bool verbose)
    {
      if(var < minval){
	if(verbose) std::cout << parameter << " is " << var << ", but must be at least " << minval << ". Setting to " << defaultval << std::endl;
	var=defaultval;
      }
    }
  template <class T>
  void checkMax(const char * parameter, T & var, const T & maxval, const T & defaultval, bool verbose)
    {
      if(var > maxval){
	if(verbose) std::cout << parameter << " is " << var << ", but must be no more than " << maxval << ". Setting to " << defaultval << std::endl;
	var=defaultval;
      }
    }

#include "mira/readseqtypes.H"

  template <class T>
  static void multiParamPrint(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr, T & varInPv0, const char * indent, const char * desc, int32 fieldlength){
    multiParamPrint_(Pv, indexesInPv, ostr, varInPv0,
		     indent, desc, fieldlength,
		     0);
  }
  template <class T>
  static void multiParamPrintNumericChar(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr, T & varInPv0, const char * indent, const char * desc, int32 fieldlength){
    multiParamPrint_(Pv, indexesInPv, ostr, varInPv0,
		     indent, desc, fieldlength,
		     1);
  }
  template <class T>
  static void multiParamPrintBool(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr, T & varInPv0, const char * indent, const char * desc, int32 fieldlength){
    multiParamPrint_(Pv, indexesInPv, ostr, varInPv0,
		     indent, desc, fieldlength,
		     2);
  }
  template <class T>
  static void multiParamPrintNagWarn(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr, T & varInPv0, const char * indent, const char * desc, int32 fieldlength){
    multiParamPrint_(Pv, indexesInPv, ostr, varInPv0,
		     indent, desc, fieldlength,
		     3);
  }

  template <class T>
  static void multiParamPrint_(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr, T & varInPv0, const char * indent, const char * desc, int32 fieldlength, uint8 printhow){
//bool asboolean, bool numericcast){

    int32 offsetOfVarInPv0=static_cast<int32>(reinterpret_cast<char *>(&varInPv0) - reinterpret_cast<char *>(&(Pv[0])));

    ostr << indent << desc;

    for(int32 firstpad=fieldlength-static_cast<int32>(strlen(desc));firstpad>1;firstpad--) ostr << ' ';

    ostr << ':';
    if(indexesInPv.size()==1){
      char * addrPv0= reinterpret_cast<char *>(&Pv[indexesInPv[0]]);
      char * valInPv0= reinterpret_cast<char *>(&(addrPv0[offsetOfVarInPv0]));
      if(printhow==3) {
	if(*valInPv0==0){
	  ostr << " no\n";
	}else if(*valInPv0==1){
	  ostr << " stop\n";
	}else{
	  ostr << " warn\n";
	}
      }else if(printhow==2) {
	if(*(reinterpret_cast<bool *>(valInPv0))) {
	  ostr << " yes\n";
	}else{
	  ostr << " no\n";
	}
      }else if(printhow==1) {
	ostr << ' ';
	charAsNum(ostr, valInPv0);
	ostr << '\n';
	//ostr << ' ' << showpos << varInPv0 << '\n';
      }else{
	ostr << ' ' << *(reinterpret_cast<T*>(valInPv0)) << '\n';
      }
    }else{
      //ostr << '\n';
      for(uint32 ipvindex=0; ipvindex<indexesInPv.size(); ipvindex++){
	if(indexesInPv[ipvindex] == SEQTYPE_END) continue;
	if(ipvindex>0) {
	  ostr << indent;
	  for(int32 pad=fieldlength;pad>0;pad--) ostr << ' ';
	}
	if(indexesInPv[ipvindex] >= static_cast<int>(Pv.size())) {
	  std::cerr << "BARF! .H MIRAParameters broken, wanted to read elem " << indexesInPv[ipvindex] << " for vector of size " << Pv.size() << " STE: " << SEQTYPE_END << "\n";
	  abort();
	}
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
	switch(indexesInPv[ipvindex]){
	case SEQTYPE_SANGER:{
	  ostr << "  [san] ";
	  break;
	}
	case SEQTYPE_454GS20:{
	  ostr << "  [454] ";
	  break;
	}
	case SEQTYPE_IONTORRENT:{
	  ostr << "  [ion] ";
	  break;
	}
	case SEQTYPE_PACBIOHQ:{
	  ostr << "  [pbh] ";
	  break;
	}
	case SEQTYPE_PACBIOLQ:{
	  ostr << "  [pbl] ";
	  break;
	}
	case SEQTYPE_TEXT:{
	  ostr << "  [txt] ";
	  break;
	}
	case SEQTYPE_SOLEXA:{
	  ostr << "  [sxa] ";
	  break;
	}
	case SEQTYPE_ABISOLID:{
	  ostr << "  [sid] ";
	  break;
	}
	default : {
	  std::cerr << "BARF! MIRAParameters broken, wanted readtype " << indexesInPv[ipvindex] << "\n";
	  abort();
	}
	}

	char * valInPvX;
	{
	  char * addrPvX= reinterpret_cast<char *>(&Pv[indexesInPv[ipvindex]]);
	  valInPvX = reinterpret_cast<char *>(&(addrPvX[offsetOfVarInPv0]));
	}

	ostr << ' ';
	if(printhow==3) {
	  if(*valInPvX==0){
	    ostr << "no\n";
	  }else if(*valInPvX==1){
	    ostr << "stop\n";
	  }else{
	    ostr << "warn\n";
	  }
	}else if(printhow==2) {
	  if(*(reinterpret_cast<bool*>(valInPvX))==true) {
	    ostr << "yes\n";
	  }else{
	    ostr << "no\n";
	  }
	}else if(printhow==1) {
	  charAsNum(ostr, valInPvX);
	  ostr << '\n';
	}else{
	  ostr << *(reinterpret_cast<T*>(valInPvX)) << '\n';
	}
      }
    }
  }

  static void charAsNum(std::ostream & ostr, char * c){
    char x=*c;
    ostr << static_cast<int>(x);
  }

  static int32 gimmeAnInt(FlexLexer * lexer, std::stringstream & errstream);
  static double gimmeADouble(FlexLexer * lexer, std::stringstream & errstream);
  static int32 getFixedStringMode(FlexLexer * lexer, std::stringstream & errstream);
  static void checkCOMMON(const std::string & currentst,
			  FlexLexer * lexer,
			  std::stringstream & errstream);
  static void checkNONCOMMON(const std::string & currentst,
			     FlexLexer * lexer,
			     std::stringstream & errstream);
  static std::string createAllTechString(const std::string & settings);


  static void interpretJobDefs(std::vector<MIRAParameters> & Pv,
			       std::vector<std::vector<uint32> > & jobdefs,
			       std::stringstream & errstream);

  static void addModifiersForMapping(std::string & modestring,
				     std::vector<std::vector<uint32> > & jobdefs,
				     std::vector<MIRAParameters> & Pv,
				     bool hasSHORTREADS);

  static void saveParsedSettingsValues(std::vector<MIRAParameters> & Pv,
				       std::vector<bool> & saved);
  static void restoreParsedSettingsValues(std::vector<MIRAParameters> & Pv,
					  std::vector<bool> & saved);

  static void initStdLocations();

public:
  MIRAParameters();
  ~MIRAParameters();

  friend std::ostream & operator<<(std::ostream &ostr, MIRAParameters const &mp);

  static void setupStdMIRAParameters(std::vector<MIRAParameters> & Pv,
				     bool verbose=false);
  static void postParsingChanges(std::vector<MIRAParameters> & Pv);

  static void dumpAllParams(std::vector<MIRAParameters> & Pv,
			    std::ostream & ostr);
  static void dumpAllParams(std::vector<MIRAParameters> & Pv,
			    const std::vector<int> & indexesInPvm,
			    std::ostream & ostr);
  static void dumpAlignParams(std::vector<MIRAParameters> & Pv,
			      const std::vector<int> & indexesInPvm,
			      std::ostream & ostr);
  static void dumpFinalMappingParams(std::vector<MIRAParameters> & Pv,
				     const std::vector<int> & indexesInPvm,
				     std::ostream & ostr);

  static void dumpGeneralParams(std::vector<MIRAParameters> & Pv,
				const std::vector<int> & indexesInPvm,
				std::ostream & ostr);
  static void dumpMiscParams(std::vector<MIRAParameters> & Pv,
			     const std::vector<int> & indexesInPvm,
			     std::ostream & ostr);
  static void dumpNagAndWarnParams(std::vector<MIRAParameters> & Pv,
			     const std::vector<int> & indexesInPvm,
			     std::ostream & ostr);
  static void dumpLoadParams(std::vector<MIRAParameters> & Pv,
			     const std::vector<int> & indexesInPvm,
			     std::ostream & ostr);
  //static void dumpShortReadParams(std::vector<MIRAParameters> & Pv,
  //			   const std::vector<int> & indexesInPvm,
  //			   std::ostream & ostr);
  static void dumpAssemblyParams(std::vector<MIRAParameters> & Pv,
				 const std::vector<int> & indexesInPvm,
				 std::ostream & ostr);
  static void dumpStrainBackboneParams(std::vector<MIRAParameters> & Pv,
				       const std::vector<int> & indexesInPvm,
				       std::ostream & ostr);
  static void dumpDataProcessingParams(std::vector<MIRAParameters> & Pv,
				       const std::vector<int> & indexesInPvm,
				       std::ostream & ostr);
  static void dumpClippingParams(std::vector<MIRAParameters> & Pv,
				 const std::vector<int> & indexesInPvm,
				 std::ostream & ostr);
  static void dumpSkimParams(std::vector<MIRAParameters> & Pv,
			     const std::vector<int> & indexesInPvm,
			     std::ostream & ostr);
  static void dumpHashStatisticsParams(std::vector<MIRAParameters> & Pv,
				       const std::vector<int> & indexesInPvm,
				       std::ostream & ostr);
  static void dumpPathfinderParams(std::vector<MIRAParameters> & Pv,
				   const std::vector<int> & indexesInPvm,
				   std::ostream & ostr);

  static void dumpContigParams(std::vector<MIRAParameters> & Pv,
			       const std::vector<int> & indexesInPvm,
			       std::ostream & ostr);
  static void dumpEditParams(std::vector<MIRAParameters> & Pv,
			     const std::vector<int> & indexesInPvm,
			     std::ostream & ostr);
  static void dumpDirectoryParams(std::vector<MIRAParameters> & Pv,
				  const std::vector<int> & indexesInPvm,
				  std::ostream & ostr);
  static void dumpFileInParams(std::vector<MIRAParameters> & Pv,
			       const std::vector<int> & indexesInPvm,
			       std::ostream & ostr);
  static void dumpFileOutParams(std::vector<MIRAParameters> & Pv,
				const std::vector<int> & indexesInPvm,
				std::ostream & ostr);
  static void dumpFileTempParams(std::vector<MIRAParameters> & Pv,
				 const std::vector<int> & indexesInPvm,
				 std::ostream & ostr);
  static void dumpOutputCustomisationParams(std::vector<MIRAParameters> & Pv,
					    const std::vector<int> & indexesInPvm,
					    std::ostream & ostr);
  static void dumpFileDirectoryOutNamesParams(std::vector<MIRAParameters> & Pv,
					      const std::vector<int> & indexesInPvm,
					      std::ostream & ostr);


  void dumpGeneralParams(std::ostream &ostr) const;
  void dump454Params(std::ostream &ostr) const;
  void dumpShortReadParams(std::ostream &ostr) const;
  void dumpAssemblyParams(std::ostream &ostr) const;
  void dumpStrainBackboneParams(std::ostream &ostr) const;
  void dumpDataProcessingParams(std::ostream &ostr) const;
  void dumpClippingParams(std::ostream &ostr) const;
  void dumpSkimParams(std::ostream &ostr) const;
  void dumpPathfinderParams(std::ostream &ostr) const;
  void dumpAlignParams(std::ostream &ostr) const;
  void dumpContigParams(std::ostream &ostr) const;
  void dumpEditParams(std::ostream &ostr) const;
  void dumpDirectoryParams(std::ostream &ostr) const;
  void dumpFileInParams(std::ostream &ostr) const;
  void dumpFileOutParams(std::ostream &ostr) const;
  void dumpFileTempParams(std::ostream &ostr) const;
  void dumpOutputCustomisationParams(std::ostream &ostr) const;
  void dumpFileDirectoryOutNamesParams(std::ostream &ostr) const;

  template <class T>
  static void multiParamPrint(std::vector<MIRAParameters> & Pv, const T & var, const std::string & desc);

  void consistencyCheck(bool verbose);

  inline special_parameters const & getSpecialParams() const
    { return mp_special_params; }
  inline nagandwarn_parameters const & getNagAndWarnParams() const
    { return mp_nagandwarn_params; }
  inline skim_parameters const & getSkimParams() const
    { return mp_skim_params; }
  inline hashstatistics_parameters const & getHashStatisticsParams() const
    { return mp_hashstatistics_params; }
  inline align_parameters const & getAlignParams() const
    { return mp_align_params; }
  inline finalmapping_parameters const & getFinalMappingParameters() const
    { return mp_finalmap_params; }
  inline dynamic_parameters const & getDynamicParams() const
    { return mp_dynamic_params; }
  inline assembly_parameters const & getAssemblyParams() const
    { return mp_assembly_params; }
  inline pathfinder_parameters const & getPathfinderParams() const
    { return mp_pathfinder_params; }
  inline contig_parameters const & getContigParams() const
    { return mp_contig_params; }
  inline directory_parameters const & getDirectoryParams() const
    { return mp_directory_params; }
  inline file_parameters const & getFileParams() const
    { return mp_file_params; }
  inline edit_parameters const & getEditParams() const
    { return mp_edit_params; }

  inline special_parameters & getNonConstSpecialParams()
    { return mp_special_params; }
  inline nagandwarn_parameters & getNonConstNagAndWarnParams()
    { return mp_nagandwarn_params;}
  inline skim_parameters & getNonConstSkimParams()
    { return mp_skim_params; }
  inline hashstatistics_parameters & getNonConstHashStatisticsParams()
    { return mp_hashstatistics_params; }
  inline align_parameters & getNonConstAlignParams()
    { return mp_align_params; }
  inline dynamic_parameters & getNonConstDynamicParams()
    { return mp_dynamic_params; }
  inline finalmapping_parameters & getNonConstFinalMappingParameters()
    { return mp_finalmap_params; }
  inline assembly_parameters & getNonConstAssemblyParams()
    { return mp_assembly_params; }
  inline pathfinder_parameters & getNonConstPathfinderParams()
    { return mp_pathfinder_params; }
  inline contig_parameters & getNonConstContigParams()
    { return mp_contig_params; }
  inline directory_parameters & getNonConstDirectoryParams()
    { return mp_directory_params; }
  inline edit_parameters & getNonConstEditParams()
    { return mp_edit_params; }

  static void loadParams(const std::string & pfile, std::vector<MIRAParameters> & Pv);
  static void parse(int argc, char **argv, std::vector<MIRAParameters> & Pv, bool verbose=true);
  static void parse(const char * params, std::vector<MIRAParameters> & Pv, bool verbose=true);
  static void parse(const std::string & params, std::vector<MIRAParameters> & Pv, bool verbose=true) {parse(params.c_str(),Pv,verbose);};
  static void parse(std::istream & is,
		    std::vector<MIRAParameters> & Pv,
		    bool verbose=true);
  static void parseQuickmode(const char * params,
			     const char * qm,
			     std::vector<MIRAParameters> & Pv,
			     bool verbose=true);
  static void parseQuickmodeNoTechSettingsChange(
    const char * params,
    const char * qm,
    std::vector<MIRAParameters> & Pv,
    bool verbose=true);

  static void generateProjectNames(std::vector<MIRAParameters> & Pv, std::string name="");
  static void generateProjectInNames(std::vector<MIRAParameters> & Pv, std::string name="");
  static void generateProjectOutNames(std::vector<MIRAParameters> & Pv, std::string name="");
  static void correctTmpDirectory(std::vector<MIRAParameters> & Pv);

  void setPathfinderMaxContigTime(uint32 t);
  void setPathfinderNextReadMaxTimeToESS(uint32 t);
  bool setAlignGapPenaltyLevel(const std::string & gplevels, std::stringstream * errstream=nullptr);
  void setContigForceNonIUPAC(bool perseq, bool amongseq);

  static const std::string & getBinPath();
  static const std::string & getBinDir();
  static const std::string & getHomeDir();
  static const std::string & getShareDir();
  static const std::string & getMHSLibDir();
};

#endif
