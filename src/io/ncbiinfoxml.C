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

/*
 * Loads a XML TraceInfo file (in NCBI standard)
 *  Calls a callback for every complete trace it loaded with the
 *  elements and their cdata values in two lists of strings
 *
 * Note: this was a quick'n dirty hack, not well thought out but works
 *       should be cleaned up if more functionality is to be entered
 *
 */


#include "ncbiinfoxml.H"

#include <cstring>

using std::cerr;
using std::endl;


// Plain vanilla constructor
NCBIInfoXML::NCBIInfoXML()
{
  FUNCSTART("NCBIInfoXML::NCBIInfoXML()");

  zeroVars();
  init();

  NIX_valid=1;
  FUNCEND();
}

void NCBIInfoXML::zeroVars()
{
  FUNCSTART("void NCBIInfoXML::zeroVars()");

  NIX_parsedata.cdata.clear();

  FUNCEND();
}

void NCBIInfoXML::init()
{
  FUNCSTART("void NCBIInfoXML::init()");
  FUNCEND();
}



NCBIInfoXML::~NCBIInfoXML()
{
  FUNCSTART("NCBIInfoXML::~NCBIInfoXML()");

  discard();

  FUNCEND();
}


void NCBIInfoXML::discard()
{
  FUNCSTART("NCBIInfoXML::discard()");

  if(NIX_valid==0){
    ERROR("Not valid yet???");
  }

  zeroVars();

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//NCBIInfoXML::NCBIInfoXML(NCBIInfoXML const &other)
//{
//  FUNCSTART("NCBIInfoXML::NCBIInfoXML(NCBIInfoXML const &other)");
//
//  NIX_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//NCBIInfoXML const & NCBIInfoXML::operator=(NCBIInfoXML const & other)
//{
//  FUNCSTART("NCBIInfoXML const & NCBIInfoXML::operator=(NCBIInfoXML const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, NCBIInfoXML const &nix)
//{
//  FUNCSTART("friend ostream & NCBIInfoXML::operator<<(ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}



void NCBIInfoXML::cdata(void *userData, const char *data, int len)
{
  XML_parsedata *pd = static_cast<XML_parsedata *>(userData);
  int i=0;
  for(;i<len;i++) {
    pd->cdata+=static_cast<char>(data[i]);
  }
}

void NCBIInfoXML::startElement(void *userData, const char *name, const char **atts)
{
  // turn off that warning about unused atts
  (void) atts;

  XML_parsedata *pd = static_cast<XML_parsedata *>(userData);

  pd->cdata.clear();
  if(strcmp(name,"trace")==0) {
    pd->acttrace.elements.clear();
    pd->acttrace.elements_cdata.clear();
  }
}

void NCBIInfoXML::endElement(void *userData, const char *name)
{
  XML_parsedata *pd = static_cast<XML_parsedata *>(userData);

  {
    char * i=const_cast<char *>(name);
    for(;*i;i++) *i=static_cast<char>(tolower(*i));
  }

  if(strcmp(name,"trace")==0) {
    pd->targetlist->push_back(pd->acttrace);
  } else {
    uint32 element_code=0;
    if(strcmp(name,"trace_name") == 0) {
      element_code=NCBIXML_TRACE_NAME;
    } else
    if(strcmp(name,"trace_file") == 0) {
      element_code=NCBIXML_TRACE_FILE;
    } else
    if(strcmp(name,"trace_type_code") == 0) {
      element_code=NCBIXML_TRACE_TYPE_CODE;
    } else
    if(strcmp(name,"clip_quality_left") == 0) {
      element_code=NCBIXML_CLIP_QUALITY_LEFT;
    } else
    if(strcmp(name,"clip_quality_right") == 0) {
      element_code=NCBIXML_CLIP_QUALITY_RIGHT;
    } else
    if(strcmp(name,"clip_vector_left") == 0) {
      element_code=NCBIXML_CLIP_VECTOR_LEFT;
    } else
    if(strcmp(name,"clip_vector_right") == 0) {
      element_code=NCBIXML_CLIP_VECTOR_RIGHT;
    } else
    if(strcmp(name,"insert_size") == 0) {
      element_code=NCBIXML_INSERT_SIZE;
    } else
    if(strcmp(name,"insert_stdev") == 0) {
      element_code=NCBIXML_INSERT_STDEV;
    } else
//    if(strcmp(name,"insert_size_min") == 0) {
//      element_code=NCBIXML_INSERT_SIZE_MIN;
//    } else
//    if(strcmp(name,"insert_size_max") == 0) {
//      element_code=NCBIXML_INSERT_SIZE_MAX;
//    } else
    if(strcmp(name,"template_id") == 0) {
      element_code=NCBIXML_TEMPLATE_ID;
    } else
    if(strcmp(name,"trace_end") == 0) {
      element_code=NCBIXML_TRACE_END;
    } else
    if(strcmp(name,"mate_ti") == 0) {
      element_code=NCBIXML_MATE_TI;
    } else
    if(strcmp(name,"mate_name") == 0) {
      element_code=NCBIXML_MATE_NAME;
    } else
    if(strcmp(name,"ti") == 0) {
      element_code=NCBIXML_TI;
    } else
    if(strcmp(name,"base_file") == 0) {
      element_code=NCBIXML_BASE_FILE;
    } else
    if(strcmp(name,"qual_file") == 0) {
      element_code=NCBIXML_QUAL_FILE;
    } else
    if(strcmp(name,"program_id") == 0) {
      element_code=NCBIXML_PROGRAM_ID;
    } else
    if(strcmp(name,"machine_type") == 0) {
      element_code=NCBIXML_MACHINE_TYPE;
    } else
    if(strcmp(name,"strain") == 0) {
      element_code=NCBIXML_STRAIN;
    }

    if(element_code>0) {
      pd->acttrace.elements.push_back(element_code);
      pd->acttrace.elements_cdata.push_back(pd->cdata);
    }
  }
  pd->cdata.clear();
//http://www.cplusplus.com/doc/tutorial/tut3-3.html
}




void NCBIInfoXML::readXMLFile(const std::string & filename, std::list<ncbitraceelements_t> & traces)
{
  FUNCSTART("void NCBIInfoXML::readXMLFile(std::string filename)");

  traces.clear();
  NIX_parsedata.targetlist=&traces;

  std::ifstream fin(filename, std::ios::in|std::ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::FATAL, "TraceInfo XML file not found for loading: " << filename);
  }

  auto xmlfilelen=fin.tellg();
  if(xmlfilelen==0){
    MIRANOTIFY(Notify::FATAL, "Zero length TraceInfo XML file: " << filename);
  }
  fin.seekg(0, std::ios::beg);
  fin.clear();


  char buf[10000];
  XML_Parser parser = XML_ParserCreate(nullptr);
  bool done;
  //int depth = 0;
  XML_SetUserData(parser, &NIX_parsedata);
  XML_SetElementHandler(parser,
			NCBIInfoXML::startElement,
			NCBIInfoXML::endElement);
  XML_SetCharacterDataHandler(parser,NCBIInfoXML::cdata);
  do {
    // got this trick from usenet: use unformatted stream i/o
    auto len = fin.rdbuf()->sgetn(buf, 10000);
    //cout << len << " " << buf;
    done = len < static_cast<int>(sizeof(buf));
    if (XML_Parse(parser, buf, static_cast<int>(len), done) == XML_STATUS_ERROR) {
      cerr << "Error: " << XML_ErrorString(XML_GetErrorCode(parser)) << " at line " << XML_GetCurrentLineNumber(parser) << endl;
      MIRANOTIFY(Notify::FATAL, "The XML Parser had a severe semantical problem, see error message and line number above: " << filename);
    }
  } while (!done);
  XML_ParserFree(parser);

  FUNCEND();
  return;
}
