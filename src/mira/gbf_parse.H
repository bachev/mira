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

#ifndef _bas_gbf_h_
#define _bas_gbf_h_


#include <iostream>
#include <fstream>
#include <string>

#include <unordered_map>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

class multitag_t;

//// apparently, some systems may lack isblank() although it's part
////  of C99
//// patch kindly provided by
////  Alfredo Jose Hernandez Alvarez <ajhernan@lcg.unam.mx>

#if HAVE_ISBLANK == 0
#define isblank(c) ((c) == ' ' || (c) == '\t')
#endif


class GBF
{
  //Variables
private:
  static const bool GBF_staticfeaturesinit;

  uint32 GBF_valid;

  std::vector<std::string>          GBF_sequencenames;
  std::vector<std::string>          GBF_sequences;
  std::vector<std::vector<multitag_t> >  GBF_tags;

  static const char * GBF_miragbfscankeys[];        /* keys used by MIRA to check if comments were stored
							 in "old GenBank" manner */

  static const char * GBF_oldmiragbf2gff3translations[];

  // maps to convert GFF features to "gap4" features"
  typedef std::unordered_map<std::string, std::string> strstrmap;
  static strstrmap GBF_mapoldmiragbf2gff3;



public:


  //Functions
private:
  void foolCompiler();
  void init();
  void zeroVars();

  void loadTheFile(std::ifstream & gbfin, uint64 lentoreserve);
  void getNextToken(const std::string & line,
		    std::string::size_type & linepos,
		    std::string & token,
		    std::string::size_type & postokenstart) const;

  //const char * getFeatureIdentifier(const std::string & feature) const;

//  void addKeyVal2TagComment(const std::string & actkey,
//			    const std::string & actval,
//			    tag_t & thetag) const;
//  void storeAccumulatedTags(tag_t & acttag, std::vector<int32> & fromto);

  void addKeyVal2TagComment(const std::string & actkey,
			    const std::string & actval,
			    std::string & comment) const;
  void storeAccumulatedTags(const std::string & guessedname,
			    const std::string & identifier,
			    const std::string & comment,
			    std::vector<int32> & fromto);

  void parseGBFLocation(const std::string & location,
			std::vector<int32> & fromto,
			int32 direction) const;
  std::string::size_type getCorrespondingRightBracket(const std::string & chars,
						 std::string::size_type start) const;
  void correctForTagPositionErrors(const std::string & gbfname);

  static void fillFeatureTranslationVectors();

  static bool staticInitialiser();

public:
  GBF();
  GBF(GBF const &other);
  ~GBF();

  GBF const & operator=(GBF const & other);
  friend std::ostream & operator<<(std::ostream &ostr, GBF const & gbf);

  void discard();

  void load(const std::string & gbfname);
  size_t getNumSequences() const {return GBF_sequences.size();}
  const std::string & getSequenceName(uint32 i) const;
  const std::string & getSequence(uint32 i) const;
  const std::vector<multitag_t> & getTags(uint32 i) const;

  void transferGeneInfoToCDSInfo();

  static bool checkIfCommentInOldMIRAGBFstyle(const std::string & comment);

  static void liftOldMIRAGBFCommentToGFF3(const std::string & src, std::string & dst);
};


#endif
