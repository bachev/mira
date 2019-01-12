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

#ifndef _bas_gffparse_h_
#define _bas_gffparse_h_

#include <iostream>
#include <string>

#include <unordered_map>

#include <regex>

#include "stdinc/defines.H"

class ReadPool;
class multitag_t;

typedef std::vector<std::string> gff3values_t;
struct gff3tagvalue_t {
  std::string       tag;
  gff3values_t values;
  friend std::ostream & operator<<(std::ostream & ostr, gff3tagvalue_t const & gf3v) {
    std::cout << "tag: " << gf3v.tag << "\n";
    for(auto & x : gf3v.values){
      std::cout << "  value: " << x << "\n";
    }
    return ostr;
  }
};
typedef std::list<gff3tagvalue_t> gff3attributes_t;

class GFFParse
{
  //Variables
private:
  static const std::string GFFP_emptystring;

  static const bool   GFFP_staticfeaturesinit;

  static const std::vector<std::string> GFFP_gff3scankeys;

  // storage for regex already used on GFF strings to extract different key/tags (Name=, Alias=, ...)
  static std::unordered_map<std::string,std::regex> GFFP_regex_extractkeys;

  std::vector<std::string> GFFP_seqnames;
  std::vector<std::string> GFFP_sequences;
  std::vector<std::vector<multitag_t> > GFFP_sequencetags;

  // map to sequence names (GFFP_seqnames)
  typedef std::unordered_map<std::string, size_t> strintmap;
  strintmap GFFP_snmap;

  // in case we deal with GFFs without sequence data, try to find sequences in this pool
  // (if given)
  ReadPool * GFFP_readpool;

  int32 GFFP_errorstatus;

public:


  //Functions
private:
  void foolCompiler();
  static bool staticInitialiser();
  void init();
  void zeroVars();

  void parseNormalGFFLine(const std::string & line,
			  const uint64 lineno,
			  std::vector<std::string> & substrs);
  void checkTagsOnceLoaded();

public:
  GFFParse(ReadPool * rp=nullptr);
  ~GFFParse();

  GFFParse(GFFParse const &other) = delete;
  GFFParse const & operator=(GFFParse const & other) = delete;

  friend std::ostream & operator<<(std::ostream &ostr, GFFParse const &gffp);

  void discard();

  void loadFile(const std::string & filename);
  size_t getNumSequences() const { return GFFP_sequences.size();};
  const std::string & getSequenceName(uint32 i) const;
  const std::string & getSequence(uint32 i) const;
  const std::vector<multitag_t> & getTags(uint32 i) const;

//  static const std::string & translateGFFfeat2GAP4feat(const std::string & feature);
//  static const std::string & translateGAP4feat2GFFfeat(const std::string & feature);

  static bool checkCommentForGFF3(const std::string & comment);
  static void extractMIRAGFF3InfoFromGFF3Attributes(const std::string & src,
						    std::string & dst,
						    std::string & source,
						    char & strand,
						    uint8 & phase,
						    std::string & miraitag);
  static void parseGFF3Attributes(const std::string & attributestr, gff3attributes_t & parsed);
  static void createGFF3AttributeString(const gff3attributes_t & attributes, std::string & attributestr);

  static std::string extractCommonName(const std::string & attributestr, bool extended=false);
  static std::string extractKeytag(const std::string & key, const std::string & attributestr);

};


#endif
