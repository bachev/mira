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

#ifndef _mira_multitag_h_
#define _mira_multitag_h_


#include <string>

#include "stdinc/stlincludes.H"

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "io/exp.H"
#include "io/annotationmappings.H"
#include "mira/stringcontainer.H"



class multitag_t {
public:
  typedef StringContainer<uint8>::Entry mte_src_t;
  typedef StringContainer<uint16>::Entry mte_id_t;
  typedef StringContainer<uint32>::Entry mte_co_t;

private:
  // string container for multitag tagsource, identifier and comments
  // tagsource = either MIRA, GenBank (when read from gbk) or read from GFF
  // identifier = id of the tag (gap4-IDs with MIRA additions)
  // comment = free string
  static StringContainer<uint8> MT_sc_mttagsource;
  static StringContainer<uint16> MT_sc_mtidentifier;
  static StringContainer<uint32> MT_sc_mtcomment;

  static const std::vector<uint32> MT_zerosizeuint32v;


  static std::vector<int8> MT_cache_identifier_isgff3so_entry;



public:
  static const mte_src_t MT_tagsrcentry_idEmpty;
  static const mte_src_t MT_tagsrcentry_idMIRA;
  static const mte_src_t MT_tagsrcentry_idGenBank;
  static const mte_src_t MT_tagsrcentry_idGFF3;

  /* if true, newComment will use a non checking version for the stringcontainer
     i.e., no every string gets added as own copy
     much faster, mostly made for throwaway reads (miraconvert/mirabait)
     only available via setting value "by hand" */

  static bool MT_fastnew;

public:
  uint32    from;                  // including
  uint32    to;                    // including
  mte_co_t  comment;               // comment to that tag (32 bit)
  mte_id_t  identifier;            // gap4 identifier (extended by mira) (16bit)
  mte_src_t source;                // where does the tag come from? (8bit)

  uint8 phase:2;    // 0,1,2 and 3 is "don't care"
  uint8 strandc:2;    // 0== =, 1== +, 2== -, 3== ?, use get/setStrand() to get/set a char

  bool  commentisgff3:1;

private:
  void dumpCommentAsAttributes(std::ostream & ostr) const;

public:

  multitag_t() : from(0),
		 to(0),
		 identifier(newIdentifier("")),
		 comment(newComment("")),
		 source(newSource("")),
		 phase(3),
		 strandc(0),
		 commentisgff3(false){
  };

  multitag_t(const std::string & pident, const std::string & pcom,
	     const std::string & psource) : from(0),
					    to(0),
					    comment(newComment(pcom)),
					    identifier(newIdentifier(pident)),
					    source(newSource(psource)),
					    phase(3),
					    strandc(0),
					    commentisgff3(false){
  };

  multitag_t(uint32 pfrom, uint32 pto, const std::string & pident, const std::string & pcom,
	     const std::string & psource, char pstrand) : from(pfrom),
							  to(pto),
							  comment(newComment(pcom)),
							  identifier(newIdentifier(pident)),
							  source(newSource(psource)),
							  phase(3),
							  commentisgff3(false){
    setStrand(pstrand);
  };

  friend std::ostream & operator<<(std::ostream &ostr, const multitag_t &t){
    t.dumpDebug(ostr);
    return ostr;
  }

  void dumpDebug(std::ostream & ostr) const;
  void dumpAsCAF(std::ostream & ostr) const;
  void dumpAsMAF(std::ostream & ostr, const char * type) const;
  void dumpAsSAM(std::ostream & ostr, int32 rlen=0) const;  // for read tags
  void dumpAsSAM(std::ostream & ostr, const std::string & contigname) const; // contig tags

  void dumpAsGFF3(std::ostream & ostr, const char * seqid) const;
  void dumpAsGFF3(std::ostream & ostr, const std::string & seqid) const {dumpAsGFF3(ostr,seqid.c_str());};


  inline const std::string & getIdentifierStr() const {return MT_sc_mtidentifier.getEntry(identifier);}
  inline static const std::string & getIdentifierStr(mte_id_t e) {return MT_sc_mtidentifier.getEntry(e);}
  inline const std::string & getCommentStr() const {return MT_sc_mtcomment.getEntry(comment);}
  inline static const std::string & getCommentStr(mte_co_t e) {return MT_sc_mtcomment.getEntry(e);}
  inline const std::string & getSourceStr() const {return MT_sc_mttagsource.getEntry(source);}
  inline static const std::string & getSourceStr(mte_src_t e) {return MT_sc_mttagsource.getEntry(e);}

  inline void setIdentifierStr(const std::string & s) {identifier=newIdentifier(s);}
  inline void setCommentStr(const std::string & s) {comment=newComment(s);}
  inline void setSourceStr(const std::string & s) {source=newSource(s);}

  static void trashContainers();

  inline bool isSourceMIRA() const {return source==MT_tagsrcentry_idMIRA;}

  void setStrand(char s) {
    FUNCSTART("void setStrand(char s)");
    switch(s) {
    case '.' :
    case '=' : {
      strandc=0;
      break;
    }
    case '+' : {
      strandc=1;
      break;
    }
    case '-' : {
      strandc=2;
      break;
    }
    case '?' : {
      strandc=3;
      break;
    }
    default : {
      MIRANOTIFY(Notify::FATAL,"Trying to set illegal strand character " << s << " (" << static_cast<uint16>(s) << ")");
      break;
    }
    }
    FUNCEND();
  };

  inline char getStrand() const {
    switch(strandc){
    case 0: {
      return '=';
    }
    case 1: {
      return '+';
    }
    case 2: {
      return '-';
    }
    default : {
      return '?';
    }
    }
  };

  inline int8 getStrandDirection() const {
    if(strandc==1) return 1;
    if(strandc==2) return -1;
    return 0;
  }

  // get that from tag_t
  inline void serialiseComment(std::string & result) const {
    tag_t::serialiseComment(getCommentStr(),result);
  }

  inline bool extractGenBankKeyValueFromComment(const std::string & whatkey,
					 std::string & result) const {
    return tag_t::extractGenBankKeyValueFromComment(getCommentStr(), whatkey, result);
  };

  inline static mte_id_t newIdentifier(const char * s) {
    //cout << "nI*: " << s << std::endl;
    return MT_sc_mtidentifier.addEntry(s);
  };
  inline static mte_id_t newIdentifier(const std::string & s) {
    //cout << "nIs: " << s << endl;
    return newIdentifier(s.c_str());
  };
  inline static mte_id_t hasIdentifier(const char * s) {
    //cout << "nI*: " << s << endl;
    return MT_sc_mtidentifier.hasEntry(s);
  };
  inline static mte_id_t hasIdentifier(const std::string & s) {
    //cout << "nIs: " << s << endl;
    return hasIdentifier(s.c_str());
  };
  inline static mte_co_t newComment(const char * s) {
    //cout << "nC*: " << s << endl;
    mte_co_t retval;
    if(MT_fastnew){
      retval=MT_sc_mtcomment.addEntryNoDoubleCheck(s);
    }else{
      retval=MT_sc_mtcomment.addEntry(s);
    }
    return retval;
  };
  inline static mte_co_t newComment(const std::string & s) {
    //cout << "nCs: " << s << endl;
    return newComment(s.c_str());
  };
  inline static mte_src_t newSource(const char * s) {
    //cout << "nC*: " << s << endl;
    return MT_sc_mttagsource.addEntry(s);
  };
  inline static mte_src_t newSource(const std::string & s) {
    //cout << "nCs: " << s << endl;
    return newSource(s.c_str());
  };

  inline bool identifierIsValidGFF3SOEntry() const {
    if(MT_cache_identifier_isgff3so_entry.size() < MT_sc_mtidentifier.size()){
      MT_cache_identifier_isgff3so_entry.resize(MT_sc_mtidentifier.size(),0);
    }
    if(MT_cache_identifier_isgff3so_entry[identifier.getSCID()]){
      if(MT_cache_identifier_isgff3so_entry[identifier.getSCID()]>0) return true;
    }else{
      if(AnnotationMappings::isValidGFF3SOEntry(getIdentifierStr())){
	MT_cache_identifier_isgff3so_entry[identifier.getSCID()]=1;
	return true;
      }
      MT_cache_identifier_isgff3so_entry[identifier.getSCID()]=-1;
    }

    return false;
  }

};

#endif
