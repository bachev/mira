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


#ifndef _bas_nix_h_
#define _bas_nix_h_

#include <fstream>

#include <list>

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"


#if __GNUC__ >= 3
#include <sstream>
#include <istream>

#else

#include <strstream>
#endif

#include <string>

#include "errorhandling/errorhandling.H"


#include <expat.h>



class NCBIInfoXML
{

public:

  // INSERT_SIZE_MIN and INSERT_SIZE_MAX are not NCBI!
  //  own expansion

  enum ncbixml_elements {
    NCBIXML_TRACE_NAME=1,
    NCBIXML_TRACE_FILE,
    NCBIXML_TRACE_TYPE_CODE,
    NCBIXML_CLIP_QUALITY_LEFT ,
    NCBIXML_CLIP_QUALITY_RIGHT ,
    NCBIXML_CLIP_VECTOR_LEFT ,
    NCBIXML_CLIP_VECTOR_RIGHT ,
    NCBIXML_INSERT_SIZE ,
    NCBIXML_INSERT_STDEV ,
    NCBIXML_INSERT_SIZE_MIN ,
    NCBIXML_INSERT_SIZE_MAX ,
    NCBIXML_TEMPLATE_ID ,
    NCBIXML_TRACE_END ,
    NCBIXML_TI ,
    NCBIXML_MATE_TI ,
    NCBIXML_MATE_NAME ,
    NCBIXML_BASE_FILE ,
    NCBIXML_QUAL_FILE,
    NCBIXML_PROGRAM_ID,
    NCBIXML_MACHINE_TYPE,
    NCBIXML_STRAIN
  };

  struct ncbitraceelements_t {
    std::list<uint32> elements;
    std::list<std::string> elements_cdata;
  };

  //Variables
private:
  uint32 NIX_valid;

  struct XML_parsedata {
    std::list<ncbitraceelements_t> * targetlist;

    std::string cdata;
    ncbitraceelements_t acttrace;
  } NIX_parsedata;


  //Functions
private:
  void init();
  void zeroVars();

  static void cdata(void *userData, const char *data, int len);
  static void startElement(void *userData, const char *name, const char **atts);
  static void endElement(void *userData, const char *name);


public:
  NCBIInfoXML();
  NCBIInfoXML(NCBIInfoXML const &other);
  ~NCBIInfoXML();

  NCBIInfoXML const & operator=(NCBIInfoXML const & other);
  friend std::ostream & operator<<(std::ostream &ostr, NCBIInfoXML const &nix);

  void discard();

  void readXMLFile(const std::string & filename,
		   std::list<ncbitraceelements_t> & traceelems
    );
};


#endif
