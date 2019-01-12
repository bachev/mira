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

#ifndef _bas_annomap_h_
#define _bas_annomap_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include <unordered_set>
#include <unordered_map>

#include <boost/bimap.hpp>


class AnnotationMappings
{
  //Variables
private:
  static const bool AM_staticinit;

  static const std::string AM_emptystring;

  typedef boost::bimap< std::string, std::string > strstr_bimap;
  typedef strstr_bimap::value_type strstr_bimap_valuetype;

  static strstr_bimap AM_map_soidANDsoterm;

  static strstr_bimap AM_map_xgap4ANDsoterm;  // eXtended gap4 identifiers

  // cannot use bimaps for these as some multiple gap4/gbf terms map to one SO term
  // e.g.: gbf: prim_transcript & precursor_RNA -> SO: primary_transcript
  // *sigh*

  static std::unordered_map<std::string,std::string> AM_map_gap4TOsoterm;
  static std::unordered_map<std::string,std::string> AM_map_sotermTOgap4;

  static std::unordered_map<std::string,std::string> AM_map_gbfTOsoterm;
  static std::unordered_map<std::string,std::string> AM_map_sotermTOgbf;

  static std::unordered_map<std::string,std::string> AM_map_oldsotermTOsoterm;

  static std::unordered_set<std::string> AM_set_miratags;

  //static strstr_bimap AM_map_gap4ANDsoterm;
  //static strstr_bimap AM_map_gbfANDsoterm;



public:


  //Functions
private:
  void init();
  void zeroVars();
  static bool staticInitialiser();

  static const std::string & translateFeat2Other(const strstr_bimap::right_map & map, const std::string & text) {
    // first use of C++11 feature in MIRA!
    // and that's because I initially could not figure out the return type of
    //   AM_map_gap4ANDsoterm.right.find(feature);
    // from the BOOST::bimap docs
    //
    // oh well, let's start using C++11 right now

    auto I=map.find(text);
    if(I==map.end()) return AM_emptystring;
    return I->second;
  }
  static const std::string & translateOther2Feat(const strstr_bimap::left_map & map, const std::string & text) {
    auto I=map.find(text);
    if(I==map.end()) return AM_emptystring;
    return I->second;
  }

  static const std::string & translateStrToStr(const std::unordered_map<std::string,std::string> & map, const std::string & text) {
    auto I=map.find(text);
    if(I==map.end()) return AM_emptystring;
    return I->second;
  }

public:
  AnnotationMappings();
  AnnotationMappings(AnnotationMappings const &other);
  ~AnnotationMappings();

  AnnotationMappings const & operator=(AnnotationMappings const & other);
  friend std::ostream & operator<<(std::ostream &ostr, AnnotationMappings const &am);

  void discard();

  inline static const std::string & translateSOfeat2GAP4feat(const std::string & feature){
    return translateStrToStr(AM_map_sotermTOgap4,feature);
  }
  inline static const std::string & translateGAP4feat2SOfeat(const std::string & feature) {
    return translateStrToStr(AM_map_gap4TOsoterm,feature);
  };


  inline static const std::string & translateSOfeat2GenBankfeat(const std::string & feature){
    return translateStrToStr(AM_map_sotermTOgbf,feature);
  }
  inline static const std::string & translateGenBankfeat2SOfeat(const std::string & feature){
    return translateStrToStr(AM_map_gbfTOsoterm,feature);
  };


  inline static const std::string & translateSOfeat2XGAP4feat(const std::string & feature){
    return translateFeat2Other(AM_map_xgap4ANDsoterm.right,feature);
  }
  inline static const std::string & translateXGAP4feat2SOfeat(const std::string & feature) {
    return translateOther2Feat(AM_map_xgap4ANDsoterm.left,feature);
  };



  inline static const std::string & translateSOfeat2SOID(const std::string & feature){
    return translateFeat2Other(AM_map_soidANDsoterm.right,feature);
  }
  inline static const std::string & translateSOID2SOfeat(const std::string & feature){
    return translateOther2Feat(AM_map_soidANDsoterm.left,feature);
  };

  inline static const std::string & translateOldSOfeat2SOfeat(const std::string & feature){
    return translateStrToStr(AM_map_oldsotermTOsoterm,feature);
  }


  static bool isValidGFF3SOEntry(const std::string & text){
    if(!translateSOfeat2SOID(text).empty()) return true;
    if(!translateSOID2SOfeat(text).empty()) return true;
    return false;
  }

  static bool isValidGAP4Entry(const std::string & text){
    if(translateXGAP4feat2SOfeat(text).empty()) return false;
    return true;
  }

  static bool isMIRAEntry(const std::string & text){
    return AM_set_miratags.find(text)!=AM_set_miratags.end();
  }
};


#endif
