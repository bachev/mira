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

#include "annotationmappings.H"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

using std::cout;
using std::cerr;
using std::endl;

const std::string AnnotationMappings::AM_emptystring;

AnnotationMappings::strstr_bimap AnnotationMappings::AM_map_xgap4ANDsoterm;
AnnotationMappings::strstr_bimap AnnotationMappings::AM_map_soidANDsoterm;
//AnnotationMappings::strstr_bimap AnnotationMappings::AM_map_gap4ANDsoterm;
//AnnotationMappings::strstr_bimap AnnotationMappings::AM_map_gbfANDsoterm;

std::unordered_map<std::string,std::string> AnnotationMappings::AM_map_gap4TOsoterm;
std::unordered_map<std::string,std::string> AnnotationMappings::AM_map_sotermTOgap4;

std::unordered_map<std::string,std::string> AnnotationMappings::AM_map_gbfTOsoterm;
std::unordered_map<std::string,std::string> AnnotationMappings::AM_map_sotermTOgbf;

std::unordered_map<std::string,std::string> AnnotationMappings::AM_map_oldsotermTOsoterm;

std::unordered_set<std::string> AnnotationMappings::AM_set_miratags;

const bool AnnotationMappings::AM_staticinit=AnnotationMappings::staticInitialiser();

// Plain vanilla constructor
AnnotationMappings::AnnotationMappings()
{
  FUNCSTART("AnnotationMappings::AnnotationMappings()");

  zeroVars();
  init();

  ERROR("Not implemented yet.");
  FUNCEND();
}

void AnnotationMappings::zeroVars()
{
  FUNCSTART("void AnnotationMappings::zeroVars()");
  FUNCEND();
}

void AnnotationMappings::init()
{
  FUNCSTART("void AnnotationMappings::init()");
  FUNCEND();
}



AnnotationMappings::~AnnotationMappings()
{
  FUNCSTART("AnnotationMappings::~AnnotationMappings()");

  discard();

  ERROR("Not implemented yet.");
  FUNCEND();
}


void AnnotationMappings::discard()
{
  FUNCSTART("AnnotationMappings::discard()");

  zeroVars();

  FUNCEND();
}

// remember: this is run before main(), so cannot really throw something here
//  as nothing would catch it, hence no error message printed
bool AnnotationMappings::staticInitialiser()
{
  FUNCSTART("bool AnnotationMappings::staticInitialiser()");

  {
    std::istringstream tmpis;
    static const char mappingfile[] = {
#include "mira_ft_set.xxd.H"
      ,0
    };
    tmpis.str(mappingfile);

    std::string tmpline;

    while(getline(tmpis,tmpline)){
      boost::trim(tmpline);
      if(tmpline.empty()
	 || tmpline[0]=='#') continue;

      AM_set_miratags.insert(tmpline);
    }
  }

  {
    std::istringstream tmpis;
    static const char mappingfile[] = {
#include "gap4_ft_so_map.xxd.H"
      ,0
    };
    tmpis.str(mappingfile);

    std::string tmpline;
    std::vector<std::string> tmpline_splitted;

    while(getline(tmpis,tmpline)){
      //boost::trim(tmpline);
      if(tmpline.empty()
	 || tmpline[0]=='#') continue;

      tmpline_splitted.clear();
      boost::split(tmpline_splitted, tmpline, boost::is_any_of("\t"));
      if(tmpline_splitted.size()!=6){
	std::string tmptmp(tmpline);
	boost::trim(tmptmp);
	if(tmptmp.empty()) continue;
	cerr << "gap4_ft_so_map.xxd.H: line contains " << tmpline_splitted.size() << " entries instead of 6. Line:\n"<<tmpline<<"\nOuch.\n";
	exit(1000);
      }
      // SO term empty? nothing to map
      if(tmpline_splitted[2].empty()) continue;
      if(!tmpline_splitted[0].empty()){
	//AM_map_gap4ANDsoterm.insert(strstr_bimap_valuetype(tmpline_splitted[0],tmpline_splitted[2]));
	AM_map_gap4TOsoterm[tmpline_splitted[0]]=tmpline_splitted[2];
	AM_map_sotermTOgap4[tmpline_splitted[2]]=tmpline_splitted[0];
      }
      if(!tmpline_splitted[1].empty()){
	//AM_map_gbfANDsoterm.insert(strstr_bimap_valuetype(tmpline_splitted[1],tmpline_splitted[2]));
	AM_map_gbfTOsoterm[tmpline_splitted[1]]=tmpline_splitted[2];
	AM_map_sotermTOgbf[tmpline_splitted[2]]=tmpline_splitted[1];
      }
    }
  }


  {
    std::istringstream tmpis;
    static const char mappingfile[] = {
#include "so.obo.xxd.H"
      ,0
    };
    tmpis.str(mappingfile);

    std::string soterm;
    std::string soid;
    std::string tmpline;
    std::vector<std::string> tmpline_splitted;

    while(getline(tmpis,tmpline)){
      //cout << "Line: " << tmpline << endl;
      boost::trim(tmpline);
      if(tmpline.empty()
	 || tmpline[0]=='#') continue;
      if(tmpline[0]=='['){
	if(!soterm.empty() && !soid.empty()){
	  //cout << "Adding ###" << soterm << "###" << soid << "###\n";
	  AM_map_soidANDsoterm.insert(strstr_bimap_valuetype(soid,soterm));
	}
	soterm.clear();
	soid.clear();
      }else{
	tmpline_splitted.clear();
	boost::split(tmpline_splitted, tmpline, boost::is_any_of("\t "));
	if(tmpline_splitted.size()!=2){
	  // Do nothing
	  // apparently sometimes the OBO contains errors where spaces have not been replaced by underscores
	  // In this case, silently drop this term

	  // cerr << "so.obo.xxd.H: line contains " << tmpline_splitted.size() << " entries instead of 2. Line:\n"<<tmpline<<"\nOuch.\n";
	  //exit(1000);

	}else{
	  if(tmpline_splitted[0]=="id:"){
	    soid=tmpline_splitted[1];
	  }else if(tmpline_splitted[0]=="name:"){
	    soterm=tmpline_splitted[1];
	  }else{
	    cerr << "so.obo.xxd.H: unexpected key in line:\n"<<tmpline<<"\nOuch.\n";
	  }
	}
      }
    }
  }

  {
    std::istringstream tmpis;
    static const char mappingfile[] = {
#include "oldso2so.xxd.H"
      ,0
    };
    tmpis.str(mappingfile);

    std::string oldsoterm;
    std::string soterm;
    std::string tmpline;
    std::vector<std::string> tmpline_splitted;

    while(getline(tmpis,tmpline)){
      //cout << "Line: " << tmpline << endl;
      boost::trim(tmpline);
      if(tmpline.empty()
	 || tmpline[0]=='#') continue;
      tmpline_splitted.clear();
      boost::split(tmpline_splitted, tmpline, boost::is_any_of("\t"));
      if(tmpline_splitted.size()!=2){
	std::string tmptmp(tmpline);
	boost::trim(tmptmp);
	if(tmptmp.empty()) continue;
	cerr << "oldso.xxd.H: line contains " << tmpline_splitted.size() << " entries instead of 2. Line:\n"<<tmpline<<"\nOuch.\n";
	exit(1000);
      }
      if(translateSOfeat2SOID(tmpline_splitted[1]).empty()){
	cerr << "oldso2so.xxd.H: new SO term " << tmpline_splitted[1] << " does not exist in SO?\n";
	exit(1000);
      }
      AM_map_oldsotermTOsoterm[tmpline_splitted[0]]=tmpline_splitted[1];
    }
  }


  {
    std::istringstream tmpis;
    static const char mappingfile[] = {
#include "so2xgap4map.xxd.H"
      ,0
    };
    tmpis.str(mappingfile);

    std::string tmpline;
    std::vector<std::string> tmpline_splitted;

    while(getline(tmpis,tmpline)){
      //cout << "Line: " << tmpline << endl;
      boost::trim(tmpline);
      if(tmpline.empty()
	 || tmpline[0]=='#') continue;
      tmpline_splitted.clear();
      boost::split(tmpline_splitted, tmpline, boost::is_any_of("\t "));
      if(tmpline_splitted.size()!=2){
	cerr << "so2xgap4map.xxd.H: line contains " << tmpline_splitted.size() << " entries instead of 2. Line:\n"<<tmpline<<"\nOuch.\n";
	exit(1000);
      }
      if(tmpline_splitted[0].empty()){
	cerr << "so2xgap4map.xxd.H: line has empty key field? Line:\n"<<tmpline<<"\nOuch.\n";
	exit(1000);
      }
      if(tmpline_splitted[1].empty()){
	cerr << "so2xgap4map.xxd.H: line has empty value field? Line:\n"<<tmpline<<"\nOuch.\n";
	exit(1000);
      }
      AM_map_xgap4ANDsoterm.insert(strstr_bimap_valuetype(tmpline_splitted[1],tmpline_splitted[0]));
    }
    //cout << "xgap 2 so map: " << AM_map_xgap4ANDsoterm.size() << endl;
  }


  {
    char buffer[24];
    uint32 xnum=4095;
    for(;true; --xnum){
      sprintf(buffer, "%03x", xnum);
      std::string xgap4="X"+static_cast<std::string>(&buffer[0]);
      if(!translateXGAP4feat2SOfeat(xgap4).empty()) break;
    }
    ++xnum;

    //cout << "highest xnum:" << xnum << endl;

    for(auto sotI=AM_map_soidANDsoterm.right.begin(); sotI!=AM_map_soidANDsoterm.right.end(); ++sotI){
      //cout << sotI->first << '\t'  << translateSOfeat2GAP4feat(sotI->first) << endl;
      if(translateSOfeat2GAP4feat(sotI->first).empty() && translateSOfeat2XGAP4feat(sotI->first).empty()){
	sprintf(buffer, "%03x", xnum);
	std::string xgap4="X"+static_cast<std::string>(&buffer[0]);
	if(translateXGAP4feat2SOfeat(xgap4).empty()){
	  AM_map_xgap4ANDsoterm.insert(strstr_bimap_valuetype(xgap4,sotI->first));
	}else{
	  cerr << "Tried to add " << sotI->first << " as new xgap4 term, but xgap4 id " << xgap4 << " already exists?\nOuch.\n";
	  exit(1000);
	}
	cout << "New gap4 xdb: " << sotI->first << "\t" << xgap4 << endl;
	++xnum;
      }
    }
  }

  return true;
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//AnnotationMappings::AnnotationMappings(AnnotationMappings const &other)
//{
//  FUNCSTART("AnnotationMappings::AnnotationMappings(AnnotationMappings const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//AnnotationMappings const & AnnotationMappings::operator=(AnnotationMappings const & other)
//{
//  FUNCSTART("AnnotationMappings const & AnnotationMappings::operator=(AnnotationMappings const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, AnnotationMappings const &???)
//{
//  FUNCSTART("friend ostream & AnnotationMappings::operator<<(ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}
