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

#ifndef _mira_read_h_
#define _mira_read_h_


#include <string>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "mira/enums.H"
#include "mira/readgrouplib.H"
#include "mira/multitag.H"

#include "util/memusage.H"

template <class T> class StringContainer;
class multitag_t;

#define BQ_UNDEFINED 255

class Read
{

public:
  enum read_output_type {AS_ILLEGAL=0, AS_TEXT, AS_TEXTSHORT, AS_TEXTCLIPS, AS_CAF, AS_BAF, AS_MAF, AS_FASTA, AS_FASTQ, AS_CLIPPEDFASTA, AS_MASKEDMASKFASTA, AS_SEQVECMASKEDFASTA, AS_FASTAQUAL, AS_CLIPPEDFASTAQUAL, AS_MASKEDMASKFASTAQUAL, AS_SEQVECMASKEDFASTAQUAL, AS_ACE, AS_ACE_COMPLEMENT, AS_GAP4DA, AS_READNAME, AS_GBF, AS_PHD};

  enum bhashstat_bits {BHSBITS_NONE=0,
		       BHSBITS_DEFAULT=0,
		       BHSBITS_FREQUENCY=7,
		       BHSBITS_SEENATLOWPOS=8,
		       BHSBITS_CONFIRMED_FWDREV=0x10,
		       BHSBITS_CONFIRMED_MULTIPLESEQTYPE=0x20,
		       BHSBITS_KMERFORK=0x40,
		       BHSBITS_BITSAREVALID=0x80,
		       BHSBITS_ALL=0xff};

  /* the three bits for frequency above are used like this:
     0 no info
     1 single hash or actively masked out k-mer
     2 below nmc
     3 non-mc
     4 above nmc
     5 mc
     6 crazy


     the routines below rely on bit numbering 76543210
  */

  /*
    for each position in read, the bposhashstat_t holds a status
    for the hash made of sequence starting at that position and
    ending k bases *downstream!* of hash calculation

    I.e.: pos 20 and 10-mer
    The fwd holds status of hash for seqence from 20-29
    The rev holds status of hash for rev compl. seq from 20 to 11
  */

  struct bhashstat_t {
    uint8 flags;

    bhashstat_t(): flags(BHSBITS_DEFAULT) {};

    inline void clear() { flags=BHSBITS_DEFAULT; }
    inline bool isValid() const {
      return flags & BHSBITS_BITSAREVALID;
    }
    inline void setValid() {
      flags=flags|Read::BHSBITS_BITSAREVALID;
    }
    inline uint8 getFrequency() const {
      return static_cast<uint8>(flags & BHSBITS_FREQUENCY);
    }
    inline void setFrequency(uint8 frequency) {
      if(frequency>7) frequency=7;
      flags=static_cast<uint8>((flags&(~BHSBITS_FREQUENCY))|frequency);
    }

    inline void setSeenAtLowPos() {
      flags=flags|Read::BHSBITS_SEENATLOWPOS;
    }
    inline bool hasSeenAtLowPos() const {
      return flags&Read::BHSBITS_SEENATLOWPOS;
    }
    inline void setConfirmedFwdRev() {
      flags=flags|Read::BHSBITS_CONFIRMED_FWDREV;
    }
    inline bool hasConfirmedFwdRev() const {
      return flags&Read::BHSBITS_CONFIRMED_FWDREV;
    }
    inline void setConfirmedMultipleSeqType() {
      flags=flags|Read::BHSBITS_CONFIRMED_MULTIPLESEQTYPE;
    }
    inline bool hasConfirmedMultipleSeqType() const {
      return flags&Read::BHSBITS_CONFIRMED_MULTIPLESEQTYPE;
    }
    inline void setKMerFork() {
      flags=flags|Read::BHSBITS_KMERFORK;
    }
    inline bool isKMerFork() const {
      return flags&Read::BHSBITS_KMERFORK;
    }

    friend std::ostream & operator<<(std::ostream &ostr, const bhashstat_t & bhs){
      //ostr << hex << static_cast<uint16>(bhs.flags) << dec << ' ';
      if(bhs.isValid()){
	ostr << "valid ";
      }else{
	ostr << "invld ";
      }
      if(bhs.hasSeenAtLowPos()) ostr << "slp ";
      if(bhs.hasConfirmedFwdRev()) ostr << "fr ";
      if(bhs.hasConfirmedMultipleSeqType()) ostr << "mst ";
      if(bhs.isKMerFork()) ostr << "kmrf ";
      ostr << '#' << static_cast<uint16>(bhs.getFrequency()) << ' ';

      return ostr;
    }

  };

  struct bposhashstat_t {
    bhashstat_t fwd;
    bhashstat_t rev;

    inline bhashstat_t getBHashStat(int8 direction) const {
      if(direction>0) return fwd;
      return rev;
    }

    inline bool isValid(int8 direction) const {
      if(direction>0) return fwd.isValid();
      return rev.isValid();
    }

    inline uint8 getFrequency(int8 direction) const {
      if(direction>0) return fwd.getFrequency();
      return rev.getFrequency();
    }
    inline void setFrequency(int8 direction, uint8 frequency) {
      if(frequency>7) frequency=7;
      if(direction>0) {
	fwd.setFrequency(frequency);
      }else{
	rev.setFrequency(frequency);
      }
    }

    inline void setSeenAtLowPos(int8 direction) {
      if(direction>0) fwd.setSeenAtLowPos();
      rev.setSeenAtLowPos();
    }
    inline bool hasSeenAtLowPos(int8 direction) const {
      if(direction>0) return fwd.hasSeenAtLowPos();
      return rev.hasSeenAtLowPos();
    }
    inline void setConfirmedFwdRev(int8 direction) {
      if(direction>0) fwd.setConfirmedFwdRev();
      rev.setConfirmedFwdRev();
    }
    inline bool hasConfirmedFwdRev(int8 direction) const {
      if(direction>0) return fwd.hasConfirmedFwdRev();
      return rev.hasConfirmedFwdRev();
    }
    inline void setConfirmedMultipleSeqType(int8 direction) {
      if(direction>0) fwd.setConfirmedMultipleSeqType();
      rev.setConfirmedMultipleSeqType();
    }
    inline bool hasConfirmedMultipleSeqType(int8 direction) const {
      if(direction>0) return fwd.hasConfirmedMultipleSeqType();
      return rev.hasConfirmedMultipleSeqType();
    }
    inline void setKMerFork(int8 direction) {
      if(direction>0) fwd.setKMerFork();
      rev.setKMerFork();
    }
    inline bool isKMerFork(int8 direction) const {
      if(direction>0) return fwd.isKMerFork();
      return rev.isKMerFork();
    }

    friend std::ostream & operator<<(std::ostream &ostr, const bposhashstat_t & bhs){
      ostr << "f: " << bhs.fwd << ' ';
      ostr << "r; " << bhs.rev << ' ';
      return ostr;
    }

    bposhashstat_t() {};
  };


  // Static variables
public:
  static const bposhashstat_t REA_bposhashstat_default;

  // Note: due to possible static initialization fiasco,
  //       the REA_tagentry_* variables are defined/initialised
  //       in "mira/preventinitfiasco.C" and not in "mira/read.C"
  static const multitag_t::mte_id_t REA_tagentry_idEmpty;

  static const multitag_t::mte_id_t REA_tagentry_idMINF;
  static const multitag_t::mte_id_t REA_tagentry_idMIT2;

  static const multitag_t::mte_id_t REA_tagentry_idCOMM;

  static const multitag_t::mte_id_t REA_tagentry_idESDN;

  static const multitag_t::mte_id_t REA_tagentry_idSRMr;
  static const multitag_t::mte_id_t REA_tagentry_idCRMr;
  static const multitag_t::mte_id_t REA_tagentry_idWRMr;
  static const multitag_t::mte_id_t REA_tagentry_idSAOr;
  static const multitag_t::mte_id_t REA_tagentry_idSROr;
  static const multitag_t::mte_id_t REA_tagentry_idSIOr;
  /* m version: "Marker"
    these get written to GFF3 to indicate the position of 'c' version in the consensus tags so
    that resulting GFF files can be used directly for subsequent mapping assemblies
    without having the 'c' version in the then new backbone
  */
  static const multitag_t::mte_id_t REA_tagentry_idSAOm;
  static const multitag_t::mte_id_t REA_tagentry_idSROm;
  static const multitag_t::mte_id_t REA_tagentry_idSIOm;
  static const multitag_t::mte_id_t REA_tagentry_idMCVm;
  static const multitag_t::mte_id_t REA_tagentry_idSRMm;
  static const multitag_t::mte_id_t REA_tagentry_idWRMm;

  static const multitag_t::mte_id_t REA_tagentry_idUNSr;
  static const multitag_t::mte_id_t REA_tagentry_idMNRr;

  static const multitag_t::mte_id_t REA_tagentry_idHAF0;
  static const multitag_t::mte_id_t REA_tagentry_idHAF1;
  static const multitag_t::mte_id_t REA_tagentry_idHAF2;
  static const multitag_t::mte_id_t REA_tagentry_idHAF3;
  static const multitag_t::mte_id_t REA_tagentry_idHAF4;
  static const multitag_t::mte_id_t REA_tagentry_idHAF5;
  static const multitag_t::mte_id_t REA_tagentry_idHAF6;
  static const multitag_t::mte_id_t REA_tagentry_idHAF7;


  static const multitag_t::mte_id_t REA_tagentry_idRLE1;
  static const multitag_t::mte_id_t REA_tagentry_idRLE2;
  static const multitag_t::mte_id_t REA_tagentry_idRLE3;
  static const multitag_t::mte_id_t REA_tagentry_idRLE4;
  static const multitag_t::mte_id_t REA_tagentry_idRLE5;
  static const multitag_t::mte_id_t REA_tagentry_idRLE6;
  static const multitag_t::mte_id_t REA_tagentry_idRLE7;
  static const multitag_t::mte_id_t REA_tagentry_idRLE8;

  static const multitag_t::mte_id_t REA_tagentry_idKMRF;

  static const multitag_t::mte_id_t REA_tagentry_idDGNr;

  static const multitag_t::mte_id_t REA_tagentry_idMFSM;


  static const multitag_t::mte_id_t REA_tagentry_idALUS;
  static const multitag_t::mte_id_t REA_tagentry_idREPT;

  static const multitag_t::mte_id_t REA_tagentry_idSVEC;

  static const multitag_t::mte_id_t REA_tagentry_idSOFAdatabank_entry;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAcontig;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAgene;
  static const multitag_t::mte_id_t REA_tagentry_idSOFACDS;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAexon;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAintron;
  static const multitag_t::mte_id_t REA_tagentry_idSOFApolyA_sequence;

  static const multitag_t::mte_id_t REA_tagentry_idSOFAmRNA;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAtranscript;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAprimary_transcript;
  static const multitag_t::mte_id_t REA_tagentry_idSOFArRNA;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAscRNA;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAsnRNA;
  static const multitag_t::mte_id_t REA_tagentry_idSOFAtRNA;

  static const multitag_t::mte_id_t REA_tagentry_idSOFAmiraitagstore;

  static const multitag_t::mte_co_t REA_tagentry_coEmpty;
  static const multitag_t::mte_co_t REA_tagentry_coUnknown;
  static const multitag_t::mte_co_t REA_tagentry_coSRMr;
  static const multitag_t::mte_co_t REA_tagentry_coCRMr;
  static const multitag_t::mte_co_t REA_tagentry_coWRMr;
  static const multitag_t::mte_co_t REA_tagentry_coSAOr;
  static const multitag_t::mte_co_t REA_tagentry_coSROr;
  static const multitag_t::mte_co_t REA_tagentry_coSIOr;
  static const multitag_t::mte_co_t REA_tagentry_coPSHP;
  static const multitag_t::mte_co_t REA_tagentry_coUNSr;

  static const multitag_t REA_defaulttag_SRMr;
  static const multitag_t REA_defaulttag_CRMr;
  static const multitag_t REA_defaulttag_WRMr;
  static const multitag_t REA_defaulttag_SAOr;
  static const multitag_t REA_defaulttag_SROr;
  static const multitag_t REA_defaulttag_SIOr;
  static const multitag_t REA_defaulttag_PSHP;
  static const multitag_t REA_defaulttag_CJSP;

  static const multitag_t REA_defaulttag_UNSr;
  static const multitag_t REA_defaulttag_MNRr;

  static const multitag_t REA_defaulttag_ED_C;
  static const multitag_t REA_defaulttag_ED_D;
  static const multitag_t REA_defaulttag_ED_I;

  static const multitag_t REA_defaulttag_SOFApolyA_sequence;

  static const std::vector<multitag_t::mte_id_t> REA_allhaftags;

private:

  // a dummy variable, used in read.C to call the static initialiser
  // before entering main;
  static const bool REA_initialisedstatics;
  static std::vector<double> REA_bqual2errorrate;


  static StringContainer<uint32> REA_sc_readname;

  // StringContainer for additional info
  static StringContainer<uint8> REA_sc_processstatus;
  static StringContainer<uint32> REA_sc_asped;

  // controls the behaviour of std::cout << Read
  // this is not thread safe! calling functions must use mutexes
  static uint64 REA_outlen;   // this atm only for FASTA output
  static uint8  REA_outtype;

  static const char  REA_zerostring;

private:
  std::string REA_template;

  // The padded sequences
  //  changes allowed afterwards
  // NOTE: IMPORTANT
  //  insertions/deletions made to one sequence are not made directly in
  //   the other. Instead, the dirty flag of the other will be set and the
  //   other sequence will be recomputed when needed.
  mutable std::vector<char> REA_padded_sequence;
  mutable std::vector<char> REA_padded_complementsequence;


  // The qualities of the bases called, for each trace one byte
  //  (as suggested for the staden package)
  // This vector is always padded
  std::vector<base_quality_t> REA_qualities;

  // The adjustment vector contains the relations between the bases
  //  in this class and the ones in the SCF file
  // Inserted bases have -1 as relation.
  std::vector<int32>          REA_adjustments;

  // NULL if unused
  std::vector<uint8> * REA_rlevptr;

  // new in 2.9.41x4
  // each base can have flags set to it, at the moment 8 suffice
  std::vector<bposhashstat_t>    REA_bposhashstats;

  // the tags set to the sequence
  std::vector<multitag_t> REA_tags;


  // left cutoffs: 10 meaning bases[0..9] cut , take [10..]

  int32 REA_ql;   // REA_qual_left_clipoff;
  int32 REA_sl;	  // REA_svec_left_clipoff;
  int32 REA_cl;	  // REA_clon_left_clipoff;
  int32 REA_ml;	  // REA_maskbase_left_clipoff;

  // right cutoffs: 100 meaning take bases [..99], cut [100..]
  int32 REA_qr;	  // REA_qual_right_clipoff;
  int32 REA_sr;	  // REA_svec_right_clipoff;
  int32 REA_cr;   // REA_clon_right_clipoff;
  int32 REA_mr;   // REA_maskbase_right_clipoff;

  /* ergo: xL and xR make normal c delimiters like in
     for(int i=xL, i<xR, i++) ... */


  int32  REA_templateid;               // ID of template
  int32  REA_templatepartnerid;        // ID of read in same template


  StringContainer<uint32>::Entry REA_nameentry;
  StringContainer<uint32>::Entry REA_asped;
//  StringContainer<uint16>::Entry REA_clonevec_name;
//  StringContainer<uint16>::Entry REA_seqvec_name;

//  StringContainer<uint8>::Entry REA_machine_type;
//  StringContainer<uint8>::Entry REA_primer;
//  StringContainer<uint8>::Entry REA_strain;
//  StringContainer<uint8>::Entry REA_basecaller;
//  StringContainer<uint8>::Entry REA_dye;
  StringContainer<uint8>::Entry REA_processstatus;

  StringContainer<uint8>::Entry REA_scf_pathname;
  StringContainer<uint8>::Entry REA_exp_pathname;

  ReadGroupLib::ReadGroupID REA_rgid;

//  char REA_template_end;    /* equivalent TRACE_END from TRACEINFO file
//			       or to the ...F or ...R from TIGR naming scheme
//			       or to p1cx and q1cx from Sanger scheme
//
//			       stores whether the read is from Forward,
//			       or Reverse end of template
//			       or unknowN
//			       'F', 'R', 'N' (N is default) */

  /* like old REA_template_end, but more versatile
     0 : unknown / undefined (REA_template_end = 'N')
     1 : 1st segment  (REA_template_end = 'F')
     2..254: nth segment
     255: last segment (REA_template_end = 'R') */
  uint8 REA_template_segment;

  // the dirty flags for the padded and padded complement sequence
  mutable bool REA_ps_dirty:1;
  mutable bool REA_pcs_dirty:1;

  bool REA_has_quality:1;
  bool REA_has_basehashstats:1;
  bool REA_has_freqavg:1;
  bool REA_has_freqrept:1;

  bool REA_has_kmerfork:1;

  // has the SCF been found and do the data seem valid?
  bool REA_scf_available:1;
  bool REA_scf_loadattempted:1;

//  bool REA_name_scheme_valid:1;       /* was the name scheme valid, i.e.,
//				       is it paired-end? */

  // Last not least, has this read valid data?
  bool REA_has_valid_data:1;

  // is it used?
  bool REA_used_in_assembly:1;

//  bool REA_isbackbone:1;
//  bool REA_israil:1;
//  bool REA_iscoverageequivalentread:1;

  // space saving
  //
  // true by default
  // is set to false by disallowAdjustments() and from then on
  //  the read does not use the REA_adjustments vector anymore
  //  until re-initialised by discard()

  bool REA_uses_adjustments:1;


// Functions

private:
  static bool staticInitialiser();

  void foolCompiler();
  void init();
  void zeroVars();

  const char * sanityCheck() const;

  void checkQualities();

  void postLoadEXPFASTA();

  //void construct_names(const char * givenname);
  void makeComplement(const std::vector<char> & source, std::vector<char> & destination) const;
  inline void updateClipoffs() const {};
  void updateTagBaseInserted(uint32 position);
  void updateTagBaseDeleted(uint32 position);
  inline void refreshPaddedSequence() const {if(REA_ps_dirty) helper_refreshPaddedSequence();}
  inline void refreshPaddedComplementSequence() const {if(REA_pcs_dirty) helper_refreshPaddedComplementSequence();}
  void helper_refreshPaddedSequence() const;
  void helper_refreshPaddedComplementSequence() const;
  const char * processEXPTag(const char * ptr);

  void  getInternalTemplateName_Sanger(std::string & tname, uint8 & segment);
  void  getInternalTemplateName_TIGR(std::string & tname, uint8 & segment);
  void  getInternalTemplateName_FR(std::string & tname, uint8 & segment);
  void  getInternalTemplateName_Solexa(std::string & tname, uint8 & segment);
  void  getInternalTemplateName_StLouis(std::string & tname, uint8 & segment);
  void  getInternalTemplateName_SRARaw(std::string & tname, uint8 & segment);

  // sort on "from" ASCENDING and on "to" DESCENDING and on "identifier" ASCENDING
  //
  static bool tag_t_comparator(const multitag_t & t1, const multitag_t & t2)
  {
    if(t1.from < t2.from) return true;
    if(t1.from > t2.from) return false;
    if(t1.to < t2.to) return true;
    if(t1.to > t2.to) return false;
    return (t1.getIdentifierStr() < t2.getIdentifierStr());
  };

  // sort on "from" ASCENDING and on "some special entries" FIXED ORDER and on "to" DESCENDING and on "identifier" ASCENDING
  //
  static bool tag_t_comparatorGFF3(const multitag_t & t1, const multitag_t & t2)
  {
    if(t1.from < t2.from) return true;
    if(t1.from > t2.from) return false;
    if(t1.identifier!=t2.identifier){
      if(t1.identifier==Read::REA_tagentry_idSOFAdatabank_entry
	 || t2.identifier==Read::REA_tagentry_idSOFAdatabank_entry) return t1.identifier==Read::REA_tagentry_idSOFAdatabank_entry;
      if(t1.identifier==Read::REA_tagentry_idSOFAcontig
	 || t2.identifier==Read::REA_tagentry_idSOFAcontig) return t1.identifier==Read::REA_tagentry_idSOFAcontig;
      if(t1.identifier==Read::REA_tagentry_idSOFAgene
	 || t2.identifier==Read::REA_tagentry_idSOFAgene) return t1.identifier==Read::REA_tagentry_idSOFAgene;
      if(t1.identifier==Read::REA_tagentry_idSOFAmRNA
	 || t2.identifier==Read::REA_tagentry_idSOFAmRNA) return t1.identifier==Read::REA_tagentry_idSOFAmRNA;
      if(t1.identifier==Read::REA_tagentry_idSOFAexon
	 || t2.identifier==Read::REA_tagentry_idSOFAexon) return t1.identifier==Read::REA_tagentry_idSOFAexon;
      if(t1.identifier==Read::REA_tagentry_idSOFACDS
	 || t2.identifier==Read::REA_tagentry_idSOFACDS) return t1.identifier==Read::REA_tagentry_idSOFACDS;
    }
    if(t1.to < t2.to) return true;
    if(t1.to > t2.to) return false;
    return (t1.getIdentifierStr() < t2.getIdentifierStr());
  };

  static double qualityToErrorRate_compute(base_quality_t qual);

  //void transferMINFTagsToReadInfo();

  void moderateContainerGrowth();

public:
  Read();
  Read(Read const &other);
  ~Read();

  Read const & operator=(Read const & other);
  friend std::ostream & operator<<(std::ostream &ostr, Read const &read);

  size_t estimateMemoryUsage() const;

  // reserve works like STL function for vector etc.
  // reserves memory in data structures so that a read
  //  can have up to 'len' bases without need for internal re-allocation
  void reserve(uint32 len);
  size_t capacity() const;
  size_t getWastage() const;

  void integrityCheck() const;

  static void dumpStringContainerStats(std::ostream & ostr);
  static void setCoutType(uint8 type);
  static void setCoutLen(uint64 len) {REA_outlen=len;}
  void dumpAsGAP4DA(std::ostream & ostr, std::string & APline, bool outputTags=true);
  void dumpTagsAsGFF3(std::ostream & ostr);
  void dumpAsReadname(std::ostream & ostr);
  void dumpAsBAF(std::ostream & ostr);
  void dumpAsCAF(std::ostream & ostr);
  void dumpAsMAF(std::ostream & ostr);
  void dumpAsFASTQ(std::ostream & ostr, bool clippedonly, bool maskedseqvec, bool maskedmask);
  void dumpAsFASTA(std::ostream & ostr, bool clippedonly, bool maskedseqvec, bool maskedmask);
  void dumpAsFASTAQual(std::ostream & ostr, bool clippedonly, bool maskedseqvec, bool maskedmask);
  void dumpAsACE(std::ostream & ostr, int32 direction);

  void getCIGARString(std::string & s) {s.clear();};

  void discard();
  const char * checkRead() const;

  void setSequenceFromString(const char * sequence, int32 len);
  inline void setSequenceFromString(const char * sequence) { setSequenceFromString(sequence,strlen(sequence));}
  inline void setSequenceFromString(const std::string & sequence) {
    setSequenceFromString(sequence.c_str(),sequence.size());
  };

  const std::vector<char> & getActualSequence() const;
  const std::vector<char> & getActualComplementSequence() const;

  inline void disallowAdjustments(){
    REA_uses_adjustments=false;
    if(REA_adjustments.capacity()>0){
      nukeSTLContainer(REA_adjustments);
    }
  }
  inline bool usesAdjustments() const {return REA_uses_adjustments;};
  inline const std::vector<int32> & getAdjustments() const
    {return REA_adjustments;}

  inline const std::vector<uint8> * getRLEValues() const
    {return REA_rlevptr;}

  // NOT! 0 terminated. getLenClippedSeq() to get its length;
  const char * getClippedSeqAsChar() const;
  const char * getClippedComplementSeqAsChar() const;
  const char * getSeqAsChar() const;
  const char * getComplementSeqAsChar() const;

  void getSeqAsString(std::string & result);

  std::vector<char>::const_iterator getSeqIteratorBegin() const;
  std::vector<char>::const_iterator getSeqIteratorEnd() const;
  std::vector<char>::const_iterator getComplementSeqIteratorBegin() const;
  std::vector<char>::const_iterator getComplementSeqIteratorEnd() const;

  std::vector<char>::const_iterator getClippedSeqIterator() const;
  std::vector<char>::const_iterator getClippedComplementSeqIterator() const;

  char getBaseInSequence(uint32 pos) const ;
  char getBaseInComplementSequence(uint32 pos) const ;
  base_quality_t getQualityInSequence(uint32 pos) const ;
  base_quality_t getQualityInComplementSequence(uint32 pos) const ;
  char getBaseInClippedSequence(uint32 pos) const;
  char getBaseInClippedComplementSequence(uint32 pos);

  inline char nocheckGetBaseInSequence(uint32 pos) const {
    refreshPaddedSequence();
    return REA_padded_sequence[pos];
  }
  inline char nocheckGetBaseInComplementSequence(uint32 pos) const {
    refreshPaddedComplementSequence();
    return REA_padded_complementsequence[pos];
  }
  inline base_quality_t nocheckGetQualityInSequence(uint32 pos) const {return REA_qualities[pos];}
  inline base_quality_t nocheckGetQualityInComplementSequence(uint32 pos) const {return REA_qualities[REA_qualities.size()-1-pos];}


  // Every change in the sequence also changes the complement (and
  //  vice versa), but the changes are not made immediately
  //  (caching mechanism in effect), only when needed.
  // inserts a base in front of the given position
  // deletes a base at the given position
  void insertBaseInSequence(char base,
			    base_quality_t quality,
			    uint32 position,
			    bool extends_clipped_area);
  void deleteBaseFromSequence(uint32 position);
  void insertBaseInComplementSequence(char base,
				      base_quality_t quality,
				      uint32 position,
				      bool extends_clipped_area);
  void deleteBaseFromComplementSequence(uint32 uposition);

  // Every change here is made directly in both vectors (costs almost
  //  no time), _IF_ they're valid.
  // changes a base at the given position
  void changeBaseInSequence(char base,
			    base_quality_t quality,
			    uint32 position);
  void changeBaseInComplementSequence(char base,
				      base_quality_t quality,
				      uint32 position);

  void changeAdjustment(uint32 position, int32 newadjustment);

  // The 'Clipped' functions work basically like the ones above,
  //  except the positions has to be given relative to the
  //  clipped sequence
  void insertBaseInClippedSequence(char base,
				   base_quality_t quality,
				   uint32 position,
				   bool extends_clipped_area);
  void deleteBaseFromClippedSequence(uint32 position);
  void insertBaseInClippedComplementSequence(char base,
					     base_quality_t quality,
					     uint32 position,
					     bool extends_clipped_area);
  void deleteBaseFromClippedComplementSequence(uint32 position);

  // Every change here is made directly in both vectors (costs almost
  //  no time), _IF_ they're valid.
  // changes a base at the given position
  void changeBaseInClippedSequence(char base,
				   base_quality_t quality,
				   uint32 position);
  void changeBaseInClippedComplementSequence(char base,
					     base_quality_t quality,
					     uint32 position);

  // This one is new, needed by disjoint de-Bruijn-Graph corrector
  void smoothSequenceReplace(uint32 position, uint32 len, std::string & replacement);

  inline const std::vector<base_quality_t> & getQualities() const
    {return REA_qualities;}
  void setQualities(const std::vector<base_quality_t> & quals);
  void setQualities(base_quality_t qual);

  inline const std::vector<bposhashstat_t> & getBPosHashStats() const
    {return REA_bposhashstats;}
  bposhashstat_t getBPosHashStats(uint32 pos) const;
  inline bposhashstat_t getBPosHashStatsInClippedSequence(uint32 pos) const {
    return getBPosHashStats(pos+getLeftClipoff());
  }
  void setBPosHashStats(bposhashstat_t bf, uint32 from, uint32 len);
  inline void setBPosHashStatsInClippedSequence(bposhashstat_t bf, uint32 from, uint32 len) {
    setBPosHashStats(bf,from+getLeftClipoff(),len);
  }

  bool hasTheseBaseHashStatsSet(bposhashstat_t bf, uint32 pos) const;
  inline void clearAllBPosHashStats() {
    REA_bposhashstats.clear();
    REA_bposhashstats.resize(getLenSeq(),REA_bposhashstat_default);
    REA_has_basehashstats=false;
  }
  inline void setHasBaseHashStats(bool b) { REA_has_basehashstats=b;}
  inline bool hasBaseHashStats() const {return REA_has_basehashstats;}

  inline void setHasFreqAvg(bool b) { REA_has_freqavg=b;}
  inline bool hasFreqAvg() const {return REA_has_freqavg;}
  inline void setHasFreqRept(bool b) { REA_has_freqrept=b;}
  inline bool hasFreqRept() const {return REA_has_freqrept;}

  inline void setHasKMerFork(bool b) { REA_has_kmerfork=b;}
  inline bool hasKMerFork() const {return REA_has_kmerfork;}

  void setTags(const std::vector<multitag_t> & tags);
  void addTagO(const multitag_t & tag);          // TODO: rename back to addTag()
  const multitag_t & getTag(uint32 tagnumber) const;
  inline const std::vector<multitag_t> & getTags() const {return REA_tags;}
  inline uint32 getNumOfTags() const { return static_cast<uint32>(REA_tags.size());}; // TODO: rename


  //bool hasTag(const string & identifier, int32 pos=-1) const;
  bool hasTag(const multitag_t::mte_id_t identifier, int32 pos=-1) const;
  bool hasTag(const std::string & identifier, int32 pos) const;
  static multitag_t::mte_id_t getTagID(const std::string & identifier);


  //uint32 countTags(const std::string & identifier, int32 pos=-1) const;
  uint32 countTags(const multitag_t::mte_id_t identifier, int32 pos=-1) const;
  void sortTags() {sort(REA_tags.begin(),REA_tags.end(),tag_t_comparator);};
  void sortTagsForGFF3() {sort(REA_tags.begin(),REA_tags.end(),tag_t_comparatorGFF3);};

  static bool extractMINFTagInfo(std::vector<multitag_t> & tags,
				 const std::string & readname,
				 std::string & strainname,
				 std::string & seqtypename,
				 std::string & machinetype,
				 int8 & tplacementcode,
				 bool & isbb,
				 bool & israil,
				 bool & isCER);

  uint32 deleteTag(const multitag_t::mte_id_t identifier);

  void clearTags() { REA_tags.clear(); };

  void setName(const std::string & name);
  inline const std::string & getName() const { return REA_sc_readname.getEntry(REA_nameentry);};

  inline const std::string & getSCFName() const { std::cout << "Read::getSCFName() killed\n"; exit(0);};

  void getSCFFullPathName(std::string & path) const;
  void getEXPFullPathName(std::string & path) const;

  void loadDataFromEXP(std::string filename, std::string path);
  void loadDataFromSCF(const std::string & scfname);
  void checkSCFAndLoadQual(bool justcheck, bool loadfailnoerror=false);

  // "sequence" has to be in fasta format (valid IUPAC codes,* is gap)
  // Concession to Staden: - is treated as N
  void initialiseRead(bool preserve_originals,  // if true, vectors are
		                                // copied, else swaped (destr.)
		      bool iscomplement,  // true, iv seq, qual, adj are compl.
		      bool ispadded,
		      ReadGroupLib::ReadGroupID rgid,
		      std::vector<char>           & sequence,
		      std::vector<base_quality_t> & qualities,
		      std::vector<int32>          & adjustments,
		      std::vector<multitag_t>     & tags,
		      const std::string & name,
		      const std::string & SCFname,
		      int32 ql, int32 qr,        // quality clipping left/right
		      int32 sl, int32 sr,        // sequence vector clipping
		      int32 cl, int32 cr         // clone vector clipping
		      );

  void removeGapsFromRead();
  void fixZeroGapQuals();
  void transformGapsToNs();
  void correctNStretch(const int32 from, const int32 to, int32 changeestim);

  inline void setTemplate(const std::string & s) {REA_template=s;};
  inline void setTemplatePartnerID(int32 id) {REA_templatepartnerid=id;}
  inline void setTemplateID(int32 id) {REA_templateid=id;};

  inline void setAsped(const std::string & s) {REA_asped=REA_sc_asped.addEntry(s);};
  inline void setProcessStatus(const std::string & s) {REA_processstatus=REA_sc_processstatus.addEntry(s);};

  inline int32 getInsizeFrom() const {return REA_rgid.getInsizeFrom();};
  inline int32 getInsizeTo() const {return REA_rgid.getInsizeTo();};

  inline bool hasTemplateInfo() const {return REA_rgid.hasTemplateInfo();};
  inline bool wantTemplateInfoEstimate() const {return REA_rgid.wantTemplateInfoEstimate();};
  inline bool hasUserDefaultQuality() const {return REA_rgid.hasUserDefaultQual();};

  inline int32 getTemplatePartnerID() const {return REA_templatepartnerid;};
  inline const std::string & getTemplate() const { return REA_template;};
  inline const std::string & getClonevecName() const { return REA_rgid.getCloneVecName();};
  inline const std::string & getSeqvecName() const { return REA_rgid.getSeqVecName();};
  inline const std::string & getStrainName() const { return REA_rgid.getStrainName();};
  inline const std::string & getBaseCaller() const { return REA_rgid.getBaseCaller();};
  inline const std::string & getAsped() const { return REA_sc_asped.getEntry(REA_asped);};
  inline const std::string & getDye() const { return REA_rgid.getDye();};
  inline const std::string & getProcessStatus() const { return REA_sc_processstatus.getEntry(REA_processstatus);};
  inline const std::string & getMachineType() const { return REA_rgid.getMachineType();};
  inline const std::string & getPrimer() const { return REA_rgid.getPrimer();};
  inline int8 getStrainID() const {return REA_rgid.getStrainID();};

  inline void setReadGroupID(ReadGroupLib::ReadGroupID rgid) {REA_rgid=rgid;}
  inline ReadGroupLib::ReadGroupID getReadGroupID() const {return REA_rgid;}

  inline uint8 getReadNamingScheme() const {return REA_rgid.getReadNamingScheme();};
  //std::string getInternalTemplateName();
  void calcTemplateInfo();
  inline int32 getTemplateID() const {return REA_templateid;};

  //inline char getTemplateEnd() const {return REA_template_end;}
  //void setTemplateEnd(char e);

  inline uint8 getTemplateSegment() const {return REA_template_segment;}
  inline void setTemplateSegment(uint8 s) {REA_template_segment=s;}

  inline int8 getTemplateBuildDirection() const {return REA_rgid.getBuildDirection();}

  inline int32 getLQClipoff() const { return REA_ql;}
  inline int32 getRQClipoff() const { return REA_qr;}
  inline int32 getLSClipoff() const { return REA_sl;}
  inline int32 getRSClipoff() const { return REA_sr;}
  inline int32 getLMClipoff() const { return REA_ml;}
  inline int32 getRMClipoff() const { return REA_mr;}

  inline int32 getLeftClipoff()  const { return std::max(REA_ql, REA_sl);}
  inline int32 getRightClipoff() const { return std::min(REA_qr, REA_sr);}
  inline uint32 getLenClippedSeq() const {
    if(unlikely(getLeftClipoff()>getRightClipoff())) return 0;
    return getRightClipoff()-getLeftClipoff();
  };

  inline uint32 getLenSeq() const {
    if(REA_ps_dirty){
      return static_cast<uint32>(REA_padded_complementsequence.size());
    }else{
      return static_cast<uint32>(REA_padded_sequence.size());
    }
  }


  int32 getLeftExtend() const;
  int32 getRightExtend() const;

  void setClipoffs(int32 lclip, int32 rclip, bool force);
  void setMinimumLeftClipoff(int32, int32, bool, bool);
  void setMinimumRightClipoff(int32, int32);
  void setLMClipoff(int32 l);
  void setRMClipoff(int32 r);
  void setLSClipoff(int32 l);
  void setRSClipoff(int32 r);
  void setLQClipoff(int32 l);
  void setRQClipoff(int32 r);

  int32 getAdjustmentPosOfReadPos(const uint32 position) const;
  int32 getLowerNonGapAdjustmentPosOfReadPos(uint32 position) const;
  int32 getUpperNonGapAdjustmentPosOfReadPos(uint32 position) const;
  int32 getReadPosOfAdjustmentPos(const int32 position) const;

  int32 getLowerNonGapPosOfReadPos(uint32 position) const;
  int32 getUpperNonGapPosOfReadPos(uint32 position) const;

  uint32 getLowerBoundPosOfBaseRun(uint32 pos,
				   char base,
				   const bool alsotakegap) const;
  uint32 getUpperBoundPosOfBaseRun(uint32 pos,
				   char base,
				   const bool alsotakegap) const;
  uint32 getLenOfGapRun(uint32 pos) const;

  bool hasSCFData(bool loadfailnoerror=false);
  //inline bool hasValidNameScheme() const {return REA_name_scheme_valid;}
  inline bool hasQuality() const {return REA_has_quality;}
  inline void setQualityFlag(bool v) {REA_has_quality=v;};
  inline void setValidData(bool v) {REA_has_valid_data=v;}
  inline bool hasValidData() const {return REA_has_valid_data;}
  inline bool isUsedInAssembly() const {return REA_used_in_assembly;}
  inline void setUsedInAssembly(bool used) {REA_used_in_assembly=used;}
  inline bool isBackbone() const {return REA_rgid.isBackbone();}
  inline bool isRail() const {return REA_rgid.isRail();}
  inline bool isCoverageEquivalentRead() const {return REA_rgid.isCoverageEquivalentRead();}

  static void trashReadNameContainer() {REA_sc_readname.trash();}

  inline uint8 getSequencingType() const {return REA_rgid.getSequencingType();}
  inline bool isSequencingType(uint8 t) const {return REA_rgid.getSequencingType()==t;}


  base_quality_t queryAverageQualInSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const;
  base_quality_t queryAverageQualInComplementSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const;
  base_quality_t queryAverageQualInClippedSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const;
  base_quality_t queryAverageQualInClippedComplementSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const;

  // inlining these makes things slower
  int32 calcComplPos(int32 pos) const;
  int32 calcClippedPos2RawPos(int32 pos) const;
  int32 calcRawPos2ClippedPos(int32 pos) const;
  int32 calcClippedComplPos2RawPos(int32 pos) const;
  //

  inline static double qualityToErrorRate(base_quality_t qual) {
    if(qual>100) qual=100;
    return REA_bqual2errorrate[qual];
  };

  static base_quality_t errorRateToQuality(double errProb);

  void performQualityClip(uint32 minqual, uint32 winlen);
  void performHardTrim();
  void setClipoffsToMaskedChars(int32 gapsize, int32 frontgs=-1, int32 endgs=-1, bool allowN=false);
  void deleteWeakestBaseInRun(const char base, const uint32 position, const bool insertgap);
  bool getPosOfWeakestBaseInRun(char base, uint32 & position);

  void transferSVTagsToClip(int32 tolerancelen, int32 clipsearchlen) ;

  void exchangeNsWithGaps();
  void blindSeqData(char newbase);
  void upDownCase(base_quality_t threshold);
  void upDownCaseClips();

  int32 getDigiNormMultiplier() const;

  static void upgradeOldTagToMultitagWithGFF3(const multitag_t & oldtag, multitag_t & mt);
  void rleRead();
  void deRLERead();

  std::string renameReadPrefix(const std::string & actname);
  static bool checkStringForSRANamingScheme(const std::string & rname);
};


// for Read
inline void emuoc_statsOnType(const Read & c, size_t & total_size, size_t & total_freecap, size_t & total_lostba, size_t & num_size, size_t & bytes_size, size_t & overhead, size_t & lostbyalign, bool & iscontainer)
{
  MMEBUG("emuoc_statsOnType(Read): " << typeid(c).name() << std::endl);

  (void) c;
  (void) total_freecap;
  (void) total_lostba;

  total_size++;
  num_size=1;
  overhead=sizeof(c);
  bytes_size=overhead;//continue
  bytes_size+=c.estimateMemoryUsage();

  iscontainer=false;

  lostbyalign=0;
  if(bytes_size%sizeof(void *)){
    lostbyalign=sizeof(void *)-bytes_size%sizeof(void *);
  }
  MMEBUG("Lost by align: " <<lostbyalign << std::endl);

  return;
}


#endif
