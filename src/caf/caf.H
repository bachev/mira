/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
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
 *
 * Reading CAF-Files
 *
 * Version 1.4   19.11.99     Thomas Pfisterer
 */


// 	$Id$



#ifndef _caf_caf_h
#define _caf_caf_h


#include "stdinc/stlincludes.H"
#include "io/scf.H"

#include "util/progressindic.H"
#include "util/fileanddisk.H"

#include "mira/read.H"
#include "mira/readpool.H"
#include "mira/contig.H"
#include "errorhandling/errorhandling.H"


#include <string>         // added by BaCh  25.01.2002
#include <unordered_map>  // added by BaCh  31.10.2010

#include <vector>

#undef  yyFlexLexer
#define yyFlexLexer CAFFlexLexer

#include <FlexLexer.h>



#undef yyFlexLexer
#define yyFlexLexer CAFFlexLexer


#define DYE_UNDEFINED -1
#define DYE_PRIMER     0
#define DYE_TERMINATOR 1


enum caf_objecttype {
  caftype_unknown,
  caftype_read,
  caftype_contig,
  caftype_group,
  caftype_assembly
};


struct caf_intervall {
  int32 leftBorder;
  int32 rightBorder;
};

typedef uint8 caf_token;

class CAF
{
  // Variables
private:
  /*
    Added by BaCh, 22.07.05
    Functionality to store read readnames in a hash
     to allow for fast loading of CAF files
    Problem before: searching a readpool for a readname gets
     VERY slow for larger CAF files (>10000 reads). The readpool
     itself cannot use hashes and provide the functionality without
     larger tweaks, so this is a workaround: the CAF object itself
     keeps track of the readnames and the assigned read IDs

     31.10.2010: changed to unordered_map
  */

  typedef std::unordered_map<std::string, int32> strmap;
  strmap CAF_hashedreadnames;

  std::list<std::string> CAF_listofreadnames;

  /*******************************************************************
   *
   *
   *
   *
   *******************************************************************/

  bool  verbose;
  bool  strict_seq_vector;
  // if strict_seq_vector is true all seqVec tags are treated as valid
  // even if they are in the middle of the read.

  FlexLexer*        CAF_lexer;           // The flexer
  caf_token         CAF_token;           // the last token read

  std::vector<MIRAParameters> *   CAF_miraparams;   // parameters for the contigs

  // for CAF streaming
  std::ifstream CAF_fin;
  std::string CAF_filename;
  ProgressIndicator<int64> * CAF_piptr;
  uint64 CAF_numseqsloaded;
  uint64 CAF_lenseqsloaded;
  uint64 CAF_numconsloaded;


  // Result Variables
  ReadPool          * CAF_readpool;         // where to store the reads
  std::list<Contig> * CAF_clist;      //
  void (*CAF_ccallbackfunc)(std::list<Contig> &, ReadPool &);
  void (*CAF_rcallbackfunc)(ReadPool &);
  bool              CAF_recalcconsensus;

  std::vector<uint32> CAF_lrperseqtype;  // longest read per seqtype
  uint8 CAF_defaultseqtype;

  int32 reReadToken;

  std::string CAF_state;
  std::string CAF_readname;
  std::string CAF_asped_date;
  std::string CAF_scf_file;
  std::string CAF_primer;
  std::string CAF_template;
  std::string CAF_ligation;
  std::string CAF_clone_vector_text;
  std::string CAF_base_caller;
  std::string CAF_stolen;
  std::string CAF_clone;
  std::string CAF_sequencing_vector;

  //  char  *CAF_seq_vector_text;
  //char  *CAF_clipping_text;


  int32  CAF_staden_id;

  std::vector<char>  CAF_DNA;
  std::vector<uint8> CAF_quality;
  std::vector<int32> CAF_align_scf;
  std::vector<char>  CAF_unprocessed;
  std::vector<multitag_t> CAF_taglist;
  std::list<Contig::contig_init_read_t> CAF_assembledFrom;
  std::list<caf_intervall> CAF_no_seq_intervall;

  caf_objecttype CAF_type;
  bool  CAF_isPadded;
  int32 CAF_dye;
  int32 CAF_length;
  int32 CAF_insert_size_min;
  int32 CAF_insert_size_max;
  int32 CAF_strand_given;
  int32 CAF_qual_clip_min;
  int32 CAF_qual_clip_max;
  int32 CAF_clone_min;
  int32 CAF_clone_max;
  int32 CAF_seq_min;
  int32 CAF_seq_max;

  // Methods
  caf_token getNextToken();  // advance in rule
  caf_token readToken();     // consume next; evtl. reread a token
  char* getIdentifier();
  std::string getText();
  std::string getReadname(std::string old_name);
  int32 get2Numbers(int &x1, int &x2);
  void  addUnprocessed(void);
  void  showVector(std::vector<char> v);

  int32 getCafReadname();
  int32 getCafProcessStatus();
  int32 getCafAspedDate();
  int32 getCafDye();
  int32 getCafStadenId();
  int32 getCafBaseCaller();
  int32 getCafStolen();
  int32 getCafClone();
  int32 getCafSCFFile();
  int32 getCafPrimer();
  int32 getCafTemplate();
  int32 getCafInsertSize();
  int32 getCafLigation();
  int32 getCafSeqVector();
  int32 getCafCloneVector();
  int32 getCafClipping();
  int32 getCafSequence();
  int32 getCafDNA();
  int32 getCafSequencingVector();
  int32 getCafQuality();
  int32 getCafAlignSCF();
  int32 getCafTag();
  int32 getCafAssembledFrom();
private:
  void cleanup();
  void cleanupContig();
  void createCafRead();
  void createCafContig();
  void createCafObject();
  void discard();

  void deescapeString(std::string & s);

  size_t countReadsBeforeLoad(const std::string & fileName);

public:
  CAF(ReadPool * aPool,
      std::list<Contig>  * clist,
      std::vector<MIRAParameters> * mp);
  ~CAF();

  size_t load(const std::string & fileName,
	      const uint8 seqtype,
	      const uint8 loadaction,
	      std::vector<uint32> & lrperseqtype,
	      bool recalcconsensus = true,
	      void (*ccallback)(std::list<Contig> &, ReadPool &) = nullptr, // callback after each contig loaded
	      void (*rcallback)(ReadPool &) = nullptr, // callback after each read loaded
	      bool isVerbose = false);

  void registerFile(const std::string & fileName);
  uint64 loadNextSeqs(uint64 numseqstoload, uint64 numconsstoload, uint64 lenseqstoload);
  void setNewContainers(ReadPool * rpool,
			std::list<Contig>  * clist,
			std::vector<MIRAParameters> * mp);
  bool checkIfEOF() {return CAF_token==0;}
  void setProgressIndicator(bool b);


  char* getVersion();

  void setStrictSeqVector(bool strict);

  friend std::ostream &operator<<(std::ostream &ostr, CAF const &i);

};





#endif
