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


#include "mira/read.H"

// for boost::trim, split
#include <boost/algorithm/string.hpp>

#include "util/misc.H"
#include "util/dptools.H"

#include "io/exp.H"
#include "io/scf.H"

#include "mira/readgrouplib.H"
#include "mira/gff_parse.H"
#include "mira/gbf_parse.H"

#include "util/stlimprove.H"

using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)

//#define PARANOIABUGTRACKFLAG
#ifdef PARANOIABUGTRACKFLAG
#define paranoiaBUGIF(ifcond, statement) { if(ifcond) {statement;}}
#else
#define paranoiaBUGIF(ifcond, statement)
#endif


const Read::bposhashstat_t Read::REA_bposhashstat_default;


// initialise the filename logic
// the read expects the name of the EXP file from which it gets
//  the SCF and CAF names

const char Read::REA_zerostring=0;

uint8 Read::REA_outtype=AS_TEXT;
uint64 Read::REA_outlen=80;



/*************************************************************************
 *
 * a static class variables initialised by static class initialiser
 *
 *
 *************************************************************************/

std::vector<double> Read::REA_bqual2errorrate;
// keep this last!
const bool Read::REA_initialisedstatics=Read::staticInitialiser();


/*************************************************************************
 *
 *  Given a position, a length and a replacement sequence,
 *   replaces the sequence in the read defined by pos and pos+len
 *   by the replacement sequence.
 *  Read length can change.
 *  Replaced sequence gets lowercase and quality values set to zero
 *   except the stretches and ends which match the replaced sequence.
 *
 * E.g.:    aaaaaaaatttatcgatgcaaaaaaa
 *          qqqqqqqqqqqqqqqqqqqqqqqqqq
 * Replace "atgcatgc" with "atgcttatgc" yields
 *          aaaaaaaatttatcgttatgcaaaaaaa
 *          qqqqqqqqqqqqqqq00qqqqqqqqqqq
 *
 * Bugs: simplistic uppercasing / quality replacement can lead to
 *       wrongl assigned qualities
 * E.g.:    aaaaaaaattaaaaaaa
 *          qqqqqqqqqqqqqqqqq
 * Replace "tt" with "tttt" yields
 * E.g.:    aaaaaaaattttaaaaaaa
 *          qqqqqqqqqqqqqqqqqqq
 * where the qualities of the new Ts are copied from existing Ts
 * Same for microrepeats (E.g.: "acgacg" with "acgacgacg")
 * ATM not really a big problem though.
 *
 * Has lots of side effects atm
 * - no 'complement' or 'clipped' version
 * - does not work on rle sequences
 * - all tags deleted
 * - base hashstats deleted
 *
 *************************************************************************/

void Read::smoothSequenceReplace(uint32 position, uint32 len, std::string & replacement)
{
  FUNCSTART("void Read::smoothSequenceReplace(uint32 position, uint32 len, std::std::string & replacement)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(REA_rlevptr!=nullptr,"Currently doesn't work on RLE");
  BUGIFTHROW(!REA_tags.empty(),"Tags not empty, not implemented atm");

  refreshPaddedSequence();

  BUGIFTHROW(position>=REA_padded_sequence.size(), getName() << ": position (" << position << ") >= size of read (" << REA_padded_sequence.size() << ")?");

  REA_pcs_dirty=true;

  int32 lendiff=static_cast<int32>(replacement.size())-static_cast<int32>(len);
  int32 newlen=static_cast<int32>(REA_padded_sequence.size())+lendiff;
  std::vector<char> newseq;
  std::vector<base_quality_t> newqual;
  std::vector<int32> newadjustments;
  newseq.reserve(newlen);
  newqual.reserve(newlen);
  if(!REA_adjustments.empty()) newadjustments.reserve(newlen);

  // create unchanged head
  uint32 pi=0;
  for(; pi<position; ++pi){
    newseq.push_back(REA_padded_sequence[pi]);
    newqual.push_back(REA_qualities[pi]);
    if(!REA_adjustments.empty()) newadjustments.push_back(REA_adjustments[pi]);
  }

  // append replaced mid section
  bool copyqualsfwd=true;
  for(char c : replacement){
    char upc=toupper(c);
    char upps=toupper(REA_padded_sequence[pi]);
    if(REA_padded_sequence[pi]==upps){
      c=upc;
    }else{
      c=tolower(c);
    }
    if(pi<REA_padded_sequence.size() && upc != upps) {
      copyqualsfwd=false;
    }
    newseq.push_back(c);
    if(copyqualsfwd){
      newqual.push_back(REA_qualities[pi]);
    }else{
      newqual.push_back(0);
    }
    if(!REA_adjustments.empty()) newadjustments.push_back(-1);
    ++pi;
  }

  // ok, transfer the quality values and sequence case at the back of replaced area
  {
    pi=position+len-1;
    auto srI=newseq.rbegin();
    auto qrI=newqual.rbegin();
    for(uint32 ti=0; ti<len; ++ti, --pi, ++srI, ++qrI){
      if(toupper(REA_padded_sequence[pi]) == toupper(*srI)){
	*qrI=REA_qualities[pi];
	//*srI=toupper(*srI);
      }else{
	break;
      }
    }
  }

  // append tail
  for(pi=position+len; pi<REA_padded_sequence.size(); ++pi){
    newseq.push_back(REA_padded_sequence[pi]);
    newqual.push_back(REA_qualities[pi]);
    if(!REA_adjustments.empty()) newadjustments.push_back(REA_adjustments[pi]);
  }

  REA_padded_sequence.swap(newseq);
  REA_qualities.swap(newqual);
  if(!newadjustments.empty()) REA_adjustments.swap(newadjustments);

  if(lendiff){
    uint32 threshold=position+len;
    if(REA_ql>=threshold) REA_ql+=lendiff;
    if(REA_sl>=threshold) REA_sl+=lendiff;
    if(REA_cl>=threshold) REA_cl+=lendiff;
    if(REA_ml>=threshold) REA_ml+=lendiff;
    if(REA_qr>=threshold) REA_qr+=lendiff;
    if(REA_sr>=threshold) REA_sr+=lendiff;
    if(REA_cr>=threshold) REA_cr+=lendiff;
    if(REA_mr>=threshold) REA_mr+=lendiff;
  }

  REA_tags.clear();
  clearAllBPosHashStats();
}


/*************************************************************************
 *
 * a static class initialiser
 *
 * initialises a few static variables / arrays _before_ main is called
 *
 *
 *************************************************************************/

bool Read::staticInitialiser()
{
  REA_bqual2errorrate.clear();
  REA_bqual2errorrate.reserve(101);
  REA_bqual2errorrate.push_back(1);
  for(base_quality_t i=1; i<=100; i++){
    REA_bqual2errorrate.push_back(qualityToErrorRate_compute(i));
  }

  return true;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// Plain vanilla constructor
Read::Read()
{
  FUNCSTART("Read::Read()");

  init();
  // Hmmm, do what?
  zeroVars();

  FUNCEND();
}

void Read::init()
{
  REA_rlevptr=nullptr;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::zeroVars()
{
  REA_rgid.resetLibId();

  nukeSTLContainer(REA_padded_sequence);
  nukeSTLContainer(REA_padded_complementsequence);

  REA_ps_dirty=false;
  REA_pcs_dirty=false;

  REA_has_quality=false;
  REA_has_basehashstats=false;
  REA_has_freqavg=false;
  REA_has_freqrept=false;

  REA_has_kmerfork=false;

  REA_scf_available=false;
  REA_scf_loadattempted=false;

  REA_has_valid_data=false;

  REA_used_in_assembly=false;

  REA_uses_adjustments=true;
  //REA_uses_adjustments=false;

  REA_template_segment=0;

  nukeSTLContainer(REA_qualities);
  nukeSTLContainer(REA_adjustments);
  nukeSTLContainer(REA_bposhashstats);
  nukeSTLContainer(REA_tags);
  if(REA_rlevptr!=nullptr){
    delete REA_rlevptr;
    REA_rlevptr=nullptr;
  }

  REA_ql=0;
  REA_sl=0;
  REA_cl=0;
  REA_ml=0;

  REA_qr=0;
  REA_sr=0;
  REA_cr=0;
  REA_mr=0;

  //REA_leftclip=0;
  //REA_rightclip=0;
  //REA_len_clipped=0;

  //REA_stadenid=-1;
  REA_templateid=-1;
  REA_templatepartnerid=-1;

  nukeSTLContainer(REA_template);

  REA_nameentry=REA_sc_readname.emptyEntry();
  REA_asped=REA_sc_asped.emptyEntry();

  REA_processstatus=REA_sc_processstatus.emptyEntry();

  //REA_template_seqtry=" ";
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Read::~Read()
{
  FUNCSTART("Read::~Read()");

  discard();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// Copy constructor
//  no discard needed as this object will be freshly created when
//  called through this constructor
Read::Read(Read const &other)
{
  FUNCSTART("Read::Read(Read const &other)");

  init();
  zeroVars();
  *this=other;                               // call the copy operator

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpStringContainerStats(std::ostream & ostr)
{
  REA_sc_readname.status(ostr);
  REA_sc_processstatus.status(ostr);
  REA_sc_asped.status(ostr);
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
const char * Read::sanityCheck() const
{
  FUNCSTART("const char * Read::sanityCheck() const");

  // no name? no sanity check
  if(getName().empty()) return nullptr;

  BUGIFTHROW(!REA_rgid.isDefaultNonValidReadGroupID() && REA_rgid.getSequencingType() >= ReadGroupLib::getNumSequencingTypes(), getName() << ": undefined technology " << static_cast<uint16>(REA_rgid.getSequencingType()));

  BUGIFTHROW(REA_has_valid_data==false, getName() << ": read has no valid data?");

  BUGIFTHROW(REA_ps_dirty==true && REA_pcs_dirty==true, getName() << "REA_ps_dirty and REA_pcs_dirty both true?");
  if(REA_ps_dirty==false && REA_pcs_dirty==false){
    BUGIFTHROW(REA_padded_sequence.size()!=REA_padded_complementsequence.size(), getName() << "Sizes of forward " << REA_padded_sequence.size() << " and complement padded " << REA_padded_complementsequence.size() << " differ.");
  }
  BUGIFTHROW(REA_ql<0, getName() << ": REA_ql " << REA_ql << " <0 ?");
  BUGIFTHROW(REA_cl<0, getName() << ": REA_cl " << REA_cl << " <0 ?");
  BUGIFTHROW(REA_sl<0, getName() << ": REA_sl " << REA_sl << " <0 ?");
  BUGIFTHROW(REA_ml<0, getName() << ": REA_ml " << REA_ml << " <0 ?");

  int32 actlen=static_cast<int32>(getLenSeq());
  BUGIFTHROW(REA_cr>actlen+1,getName() << ": REA_cr " << REA_cr << " > actlen+1 " << actlen+1 << " ?");
  BUGIFTHROW(REA_qr>actlen+1,getName() << ": REA_qr " << REA_qr << " > actlen+1 " << actlen+1 << " ?");
  BUGIFTHROW(REA_sr>actlen+1,getName() << ": REA_sr " << REA_sr << " > actlen+1 " << actlen+1 << " ?");
  BUGIFTHROW(REA_mr>actlen+1,getName() << ": REA_mr " << REA_mr << " > actlen+1 " << actlen+1 << " ?");

  //if(REA_leftclip != std::max(REA_ql, REA_sl)) return "REA_leftclip!=std::max(REA_ql, REA_sl) ?";
  //if(REA_rightclip != std::min(REA_qr, REA_sr)) return "REA_rightclip!=std::min(REA_qr, REA_sr) ?";
  //if(REA_len_clipped != REA_rightclip-REA_leftclip){
  //  if(REA_rightclip-REA_leftclip>0) return "REA_len_clipped != REA_rightclip-REA_leftclip ?";
  //}

  BUGIFTHROW(REA_uses_adjustments && static_cast<int32>(REA_adjustments.size()) != actlen, getName() << ": REA_adjustments " << REA_adjustments.size() << " != expected size " << actlen << " ?");
  BUGIFTHROW(static_cast<int32>(REA_qualities.size()) != actlen, getName() << ": REA_qualities " << REA_qualities.size() << " != expected size " << actlen << " ?");
  BUGIFTHROW(static_cast<int32>(REA_bposhashstats.size()) != actlen, getName() << ": REA_bposhashstats " << REA_bposhashstats.size() << " != " << " expected size " << actlen << " ?");
  BUGIFTHROW(REA_rlevptr!=nullptr && static_cast<int32>(REA_rlevptr->size()) != actlen, getName() << ": REA_rlevptr->size() " << REA_rlevptr->size() << " != actlen " << actlen << " ? ");

  // TODO: configure this only when bughunting
#if 0
  {
    for(const auto & te : REA_tags){
      if(te.from > actlen) {
	te.dump();
	return "From > actlen?";
      }
      if(te.to > actlen) {
	te.dump();
	return "To > actlen?";
      }
    }
  }
#endif

  return nullptr;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
const char * Read::checkRead() const
{
  //  return nullptr;


  if(sanityCheck()) {
    cout << "Sanity check of read '" << getName() << "' failed" << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    return sanityCheck();
  }

#ifdef PARANOIABUGTRACKFLAG
  {
    for(auto I=REA_padded_sequence.cbegin(); I!=REA_padded_sequence.cend(); ++I){
      if(!dptools::isValidStarBase(*I)){
	cout << "Invalid base: " << *I << "\t(" << static_cast<uint16>(*I) << ")" << endl;
	cout << "At pos: " << I-REA_padded_sequence.begin() << endl;
	return "Invalid base in padded sequence!";
      }
    }
    for(auto I=REA_padded_complementsequence.cbegin(); I!=REA_padded_complementsequence.cend(); ++I){
      if(!dptools::isValidStarBase(*I)){
	cout << "Invalid base: " << *I << "\t(" << static_cast<uint16>(*I) << ")" << endl;
	cout << "At pos: " << I-REA_padded_complementsequence.begin() << endl;
	return "Invalid base in padded complement sequence!";
      }
    }
  }

  {
    for(auto I=REA_qualities.cbegin(); I!=REA_qualities.cend(); ++I){
      if(*I>100){
	cout << "Invalid quality: " << static_cast<uint16>(*I) << endl;
	cout << "At pos: " << I-REA_qualities.begin() << endl;
	return "Invalid quality value!";
      }
    }
  }

  if(REA_uses_adjustments){
    for(auto I=REA_adjustments.cbegin(); I!=REA_adjustments.cend(); ++I){
      if(*I<-1 || *I >= static_cast<int32>(REA_qualities.size())){
	cout << "Invalid adjustment: " << *I << endl;
	cout << "At pos: " << I-REA_adjustments.begin() << endl;
	return "Invalid adjustment value!";
      }
    }
  }else{
    cout << "Does not use adjustments.\n";
  }
#endif

  return nullptr;
}


void Read::integrityCheck() const
{
  FUNCSTART("void Read::integrityCheck() const ");

  if(REA_has_valid_data==false){
    MIRANOTIFY(Notify::INTERNAL, "Read " << getName() << " has no valid data?");
  }

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::reserve(uint32 lentoreserve)
{
  REA_padded_sequence.reserve(lentoreserve);
  REA_padded_complementsequence.reserve(lentoreserve);
  REA_qualities.reserve(lentoreserve);
  REA_bposhashstats.reserve(lentoreserve);
  if(REA_uses_adjustments) REA_adjustments.reserve(lentoreserve);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
size_t Read::getWastage() const
{
  size_t w=REA_padded_sequence.capacity()-REA_padded_sequence.size();
  w+=REA_padded_complementsequence.capacity()-REA_padded_complementsequence.size();
  w+=REA_qualities.capacity()-REA_qualities.size();
  w+=sizeof(bposhashstat_t)*(REA_bposhashstats.capacity()-REA_bposhashstats.size());
  if(REA_adjustments.capacity()) {
    w+=4*(REA_adjustments.capacity()-REA_adjustments.size());
  }
  if(REA_rlevptr!=nullptr){
    w+=4*(REA_rlevptr->capacity()-REA_rlevptr->size());
  }
  return w;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
size_t Read::capacity() const
{
  return REA_bposhashstats.capacity();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::moderateContainerGrowth()
{
  static const size_t divval=10; // == 10%
  static const size_t minim=5;

  if(REA_padded_sequence.capacity()==REA_padded_sequence.size()){
    REA_padded_sequence.reserve(
      std::max(minim,
	  REA_padded_sequence.size()+REA_padded_sequence.size()/divval)
      );
  }
  if(REA_padded_complementsequence.capacity()==REA_padded_complementsequence.size()){
    REA_padded_complementsequence.reserve(
      std::max(minim,
	  REA_padded_complementsequence.size()+REA_padded_complementsequence.size()/divval)
      );
  }
  if(REA_qualities.capacity()==REA_qualities.size()){
    REA_qualities.reserve(
      std::max(minim,
	  REA_qualities.size()+REA_qualities.size()/divval)
      );
  }
  if(REA_uses_adjustments){
    if(REA_adjustments.capacity()==REA_adjustments.size()){
      REA_adjustments.reserve(
	std::max(minim,
	    REA_adjustments.size()+REA_adjustments.size()/divval)
	);
    }
  }

  if(REA_bposhashstats.capacity()==REA_bposhashstats.size()){
    REA_bposhashstats.reserve(
      std::max(minim,
	  REA_bposhashstats.size()+REA_bposhashstats.size()/divval)
      );
  }

  if(REA_rlevptr!=nullptr
     && REA_rlevptr->capacity()==REA_rlevptr->size()){
    REA_rlevptr->reserve(
      std::max(minim,
	  REA_rlevptr->size()+REA_rlevptr->size()/divval)
      );
  }

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// Copy operator, needed by copy-constructor
Read const & Read::operator=(Read const & other)
{
  FUNCSTART("Read const & Read::operator=(Read const & other)");

  if(this != &other){
    discard();

    REA_uses_adjustments=other.REA_uses_adjustments;

    //uint32 lentoreserve=0;
    //if(REA_ps_dirty==false){
    //  lentoreserve=other.REA_padded_sequence.size();
    //} else if(REA_pcs_dirty==false){
    //  lentoreserve=other.REA_padded_complementsequence.size();
    //}
    //reserve(lentoreserve+lentoreserve/10);

    //REA_name=        other.REA_name;
    REA_nameentry=        other.REA_nameentry;

    //REA_exp_filename=other.REA_exp_filename;
    //REA_caf_filename=other.REA_caf_filename;

    REA_scf_pathname=    other.REA_scf_pathname;
    REA_exp_pathname=	   other.REA_exp_pathname;

    REA_ps_dirty=other.REA_ps_dirty;
    if(REA_ps_dirty==false){
      REA_padded_sequence=other.REA_padded_sequence;
    }
    REA_pcs_dirty=other.REA_pcs_dirty;
    if(REA_pcs_dirty==false){
      REA_padded_complementsequence=other.REA_padded_complementsequence;
    }

    REA_qualities=other.REA_qualities;

    REA_adjustments=other.REA_adjustments;
    REA_bposhashstats=other.REA_bposhashstats;

    if(other.REA_rlevptr!=nullptr){
      REA_rlevptr= new std::vector<uint8>;
      *REA_rlevptr= *(other.REA_rlevptr);
    }

    REA_ql=other.REA_ql;
    REA_qr=other.REA_qr;
    REA_sl=other.REA_sl;
    REA_sr=other.REA_sr;
    REA_cl=other.REA_cl;
    REA_cr=other.REA_cr;
    REA_ml=other.REA_ml;
    REA_mr=other.REA_mr;

    //REA_leftclip=other.REA_leftclip;
    //REA_rightclip=other.REA_rightclip;
    //REA_len_clipped=other.REA_len_clipped;

    REA_tags=other.REA_tags;

    REA_has_quality=other.REA_has_quality;
    REA_has_basehashstats=other.REA_has_basehashstats;
    REA_has_freqavg=other.REA_has_freqavg;
    REA_has_freqrept=other.REA_has_freqrept;

    REA_has_kmerfork=other.REA_has_kmerfork;

    REA_scf_available=other.REA_scf_available;
    REA_scf_loadattempted=other.REA_scf_loadattempted;

    REA_has_valid_data=other.REA_has_valid_data;
    REA_used_in_assembly=other.REA_used_in_assembly;

    REA_rgid=other.REA_rgid;
    //REA_stadenid=other.REA_stadenid;
    REA_templateid=other.REA_templateid;
    REA_templatepartnerid=other.REA_templatepartnerid;
    REA_template_segment=other.REA_template_segment;

    REA_template=other.REA_template;
    REA_asped=other.REA_asped;
    REA_processstatus=other.REA_processstatus;
  }

  FUNCEND();
  return *this;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// Copy operator, needed by copy-constructor
size_t Read::estimateMemoryUsage() const
{
  FUNCSTART("size_t Read::estimateMemoryUsage()");

  size_t components=0;

  size_t cnum,cbytes,freecap,clba;

  components+=estimateMemoryUsageOfContainer(REA_padded_sequence,false,cnum,cbytes,freecap,clba);
  components+=estimateMemoryUsageOfContainer(REA_padded_complementsequence,false,cnum,cbytes,freecap,clba);
  components+=estimateMemoryUsageOfContainer(REA_qualities,false,cnum,cbytes,freecap,clba);
  components+=estimateMemoryUsageOfContainer(REA_adjustments,false,cnum,cbytes,freecap,clba);
  components+=estimateMemoryUsageOfContainer(REA_bposhashstats,false,cnum,cbytes,freecap,clba);
  components+=estimateMemoryUsageOfContainer(REA_tags,false,cnum,cbytes,freecap,clba);

  if(REA_rlevptr!=nullptr){
    components+=estimateMemoryUsageOfContainer(*REA_rlevptr,false,cnum,cbytes,freecap,clba);
  }

  FUNCEND();
  return components;
}



/*************************************************************************
 *
 * only ACGT are RLEd, IUPACs and Ns are left untouched
 *
 * Big problem: the adjustment vector
 * Initial version tried to preserve the adjustments of the original reads,
 *  but the PacBio corrector routines depend on continuous adjustment values
 *  (0,1,2,3,...) to find out where during correction insertions or deletions
 *  occurred.
 * As I am not sure what the original adjustments could be used for and
 *  I do not want to keep two vectors, I'm switching to continuous for
 *  the moment.
 * In case "near original" would be needed in future, would need to find
 *  other solution
 *
 *************************************************************************/

void Read::rleRead()
{
  FUNCSTART("void Read::rleRead()");

  BUGIFTHROW(REA_rlevptr!=nullptr,"Read " << getName() << " already RLEd?");

  REA_rlevptr=new std::vector<uint8>;
  REA_rlevptr->resize(getLenSeq());

  refreshPaddedSequence();
  REA_padded_complementsequence.clear();
  REA_pcs_dirty=true;

  uint32 rlectr=0;
  uint32 run=0;
  auto aI=REA_adjustments.cbegin();
  auto qI=REA_qualities.cbegin();
  auto sI=REA_padded_sequence.cbegin();
  char runchar=*sI;
  uint32 sumqual=0;
  int32 startadj=*aI;
  for(uint32 ri=0; ri<REA_padded_sequence.size(); ++sI, ++qI, ++aI, ++ri){
    if(unlikely(ri==0) || (*sI==runchar && dptools::isValidACGTBase(*sI))){
      ++run;
      sumqual+=*qI;
    }else{
      REA_padded_sequence[rlectr]=runchar;
      runchar=*sI;
      REA_qualities[rlectr]=sumqual/run;
      sumqual=*qI;
      (*REA_rlevptr)[rlectr]=run;
      run=1;
      // REA_adjustments[rlectr]=startadj;   // "near original adjustments"
      REA_adjustments[rlectr]=rlectr;       // continuous RLE adjustments
      startadj=*aI;

      ++rlectr;
    }
    if(REA_ql==ri) REA_ql=rlectr;
    if(REA_sl==ri) REA_sl=rlectr;
    if(REA_cl==ri) REA_cl=rlectr;
    if(REA_ml==ri) REA_ml=rlectr;
    if(REA_qr==ri) REA_qr=rlectr;
    if(REA_sr==ri) REA_sr=rlectr;
    if(REA_cr==ri) REA_cr=rlectr;
    if(REA_mr==ri) REA_mr=rlectr;

    for(auto & te : REA_tags){
      if(te.from==ri) te.from=rlectr;
      if(te.to==ri) te.to=rlectr;
    }
  }
  REA_padded_sequence[rlectr]=runchar;
  (*REA_rlevptr)[rlectr]=run;
  REA_qualities[rlectr]=sumqual/run;
  //REA_adjustments[rlectr]=startadj;
  REA_adjustments[rlectr]=rlectr;

  ++rlectr;
  REA_padded_sequence.resize(rlectr);
  REA_rlevptr->resize(rlectr);
  REA_qualities.resize(rlectr);
  REA_adjustments.resize(rlectr);
  REA_bposhashstats.clear();
  REA_bposhashstats.resize(rlectr);

  if(REA_ql>=rlectr) REA_ql=rlectr;
  if(REA_sl>=rlectr) REA_sl=rlectr;
  if(REA_cl>=rlectr) REA_cl=rlectr;
  if(REA_ml>=rlectr) REA_ml=rlectr;
  if(REA_qr>=rlectr) REA_qr=rlectr;
  if(REA_sr>=rlectr) REA_sr=rlectr;
  if(REA_cr>=rlectr) REA_cr=rlectr;
  if(REA_mr>=rlectr) REA_mr=rlectr;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

int32 drr_helper(std::vector<int32> & v, int32 p)
{
  if(v.empty() || p<0) return -1;
  if(p>=v.size()) return v.back();
  return v[p];
}

void Read::deRLERead()
{
  FUNCSTART("void Read::deRLERead()");

//  if(REA_rlevptr==nullptr){
//    setCoutType(AS_TEXT);
//    cout << "Read " << getName() << " not RLEd?\n" << *this;
//  }
  BUGIFTHROW(REA_rlevptr==nullptr,"Read " << getName() << " not RLEd?");

  uint32 totalbases=0;
  for(auto rvI=REA_rlevptr->begin(); rvI!=REA_rlevptr->end(); ++rvI){
    BUGIFTHROW(*rvI==0,"0 RLE value at pos " << rvI-REA_rlevptr->begin() << " in " << getName());
    totalbases+=*rvI;
  }

  std::vector<char> newseq;
  newseq.reserve(totalbases);
  std::vector<base_quality_t> newqual;
  newqual.reserve(totalbases);
  std::vector<int32> posmap;
  posmap.reserve(totalbases);

  refreshPaddedSequence();
  REA_pcs_dirty=true;
  {
    auto sI=REA_padded_sequence.begin();
    auto qI=REA_qualities.begin();
    for(auto rvI=REA_rlevptr->begin(); rvI!=REA_rlevptr->end(); ++rvI, ++sI, ++qI){
      posmap.push_back(newseq.size());
      for(uint8 rlei=0; rlei<*rvI; ++rlei){
	newseq.push_back(*sI);
	newqual.push_back(*qI);
      }
    }
    posmap.push_back(newseq.size());
  }

  REA_bposhashstats.clear();
  REA_bposhashstats.resize(newseq.size());

  REA_ql=drr_helper(posmap,REA_ql);
  REA_sl=drr_helper(posmap,REA_sl);
  REA_cl=drr_helper(posmap,REA_cl);
  REA_ml=drr_helper(posmap,REA_ml);
  REA_qr=drr_helper(posmap,REA_qr);
  REA_sr=drr_helper(posmap,REA_sr);
  REA_cr=drr_helper(posmap,REA_cr);
  REA_mr=drr_helper(posmap,REA_mr);

  for(auto & te : REA_tags){
    te.from=drr_helper(posmap,te.from);
    te.to=drr_helper(posmap,te.to);
  }

  REA_padded_sequence.swap(newseq);
  REA_qualities.swap(newqual);
  REA_adjustments.clear();
  REA_adjustments.resize(REA_padded_sequence.size());
  for(int32 i=0; i< REA_adjustments.size(); ++i){
    REA_adjustments[i]=i;
  }

  delete REA_rlevptr;
  REA_rlevptr=nullptr;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Read::setCoutType(uint8 type)
{
  FUNCSTART("void Read::setCoutType(uint8 type)");
  switch(type){
  case AS_TEXT:
  case AS_TEXTSHORT:
  case AS_TEXTCLIPS:
  case AS_CAF:
  case AS_MAF:
  case AS_ACE:
  case AS_ACE_COMPLEMENT:
  case AS_FASTA:
  case AS_FASTQ:
  case AS_GAP4DA:
  case AS_CLIPPEDFASTA:
  case AS_SEQVECMASKEDFASTA:
  case AS_MASKEDMASKFASTA:
  case AS_FASTAQUAL:
  case AS_CLIPPEDFASTAQUAL:
  case AS_SEQVECMASKEDFASTAQUAL:
  case AS_MASKEDMASKFASTAQUAL:
  case AS_READNAME: {
    REA_outtype=type;
    break;
  }
  case AS_ILLEGAL:
  default:{
    MIRANOTIFY(Notify::INTERNAL, "Wrong type " << static_cast<uint16>(type) << " is not known.");
  }
  }

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
std::ostream & operator<<(std::ostream &ostr, Read const &read)
{
  FUNCSTART("friend std::ostream & Read::operator<<(std::ostream &ostr, const  &read)");

//  BUGSTAT(read.checkRead());

//  if(read.REA_outtype==Read::AS_CAF){
//    Read & nonconstread = const_cast<Read &>(read);
//    nonconstread.dumpAsCAF(ostr);
//    FUNCEND();
//    return ostr;
//  } else if(read.REA_outtype==Read::AS_MASKEDFASTA){
//    Read & nonconstread = const_cast<Read &>(read);
//    nonconstread.dumpAsFASTA(ostr,false,true);
//    FUNCEND();
//    return ostr;
//  } else if(read.REA_outtype==Read::AS_CLIPPEDFASTA){
//    Read & nonconstread = const_cast<Read &>(read);
//    nonconstread.dumpAsFASTA(ostr,true, false);
//    FUNCEND();
//    return ostr;
//  } else if(read.REA_outtype==Read::AS_FASTA){
//    Read & nonconstread = const_cast<Read &>(read);
//    nonconstread.dumpAsFASTA(ostr,false, false);
//    FUNCEND();
//    return ostr;
//  } else if(read.REA_outtype==Read::AS_ACE){
//    Read & nonconstread = const_cast<Read &>(read);
//    nonconstread.dumpAsACE(ostr, 1);
//    FUNCEND();
//    return ostr;
//  } else if(read.REA_outtype==Read::AS_ACE_COMPLEMENT){
//    Read & nonconstread = const_cast<Read &>(read);
//    nonconstread.dumpAsACE(ostr, -1);
//    FUNCEND();
//    return ostr;
//  }

  bool stopit=true;
  switch(read.REA_outtype) {
  case Read::AS_CAF : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsCAF(ostr);
    break;
  }
  case Read::AS_MAF : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsMAF(ostr);
    break;
  }
  case Read::AS_BAF : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsBAF(ostr);
    break;
  }
  case Read::AS_FASTA : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTA(ostr,false, false, false);
    break;
  }
  case Read::AS_FASTQ : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTQ(ostr,false, false, false);
    break;
  }
  case Read::AS_SEQVECMASKEDFASTA : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTA(ostr,false,true, false);
    break;
  }
  case Read::AS_MASKEDMASKFASTA : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTA(ostr,false,false, true);
    break;
  }
  case Read::AS_CLIPPEDFASTA : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTA(ostr,true, false, false);
    break;
  }
  case Read::AS_FASTAQUAL : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTAQual(ostr,false, false, false);
    break;
  }
  case Read::AS_SEQVECMASKEDFASTAQUAL : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTAQual(ostr,false,true, false);
    break;
  }
  case Read::AS_MASKEDMASKFASTAQUAL : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTAQual(ostr,false,false, true);
    break;
  }
  case Read::AS_CLIPPEDFASTAQUAL : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsFASTAQual(ostr,true, false, false);
    break;
  }
  case Read::AS_GAP4DA : {
    Read & nonconstread = const_cast<Read &>(read);
    std::string gna;
    nonconstread.dumpAsGAP4DA(ostr, gna, true);
    break;
  }
  case Read::AS_ACE : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsACE(ostr, 1);
    break;
  }
  case Read::AS_ACE_COMPLEMENT : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsACE(ostr, -1);
    break;
  }
  case Read::AS_READNAME : {
    Read & nonconstread = const_cast<Read &>(read);
    nonconstread.dumpAsReadname(ostr);
    break;
  }
  default: {
    stopit=false;
  }
  }

  if(stopit) {
    FUNCEND();
    return ostr;
  }

  ostr << std::dec;

  ostr << "\nName: " << read.getName();
  ostr << "\nql: " << read.REA_ql << "\tqr: " << read.REA_qr;
  ostr << "\nsl: " << read.REA_sl << "\tsr: " << read.REA_sr;
  ostr << "\ncl: " << read.REA_cl << "\tcr: " << read.REA_cr;
  ostr << "\nml: " << read.REA_ml << "\tmr: " << read.REA_mr;
  ostr << "\nLeftclip: " << read.getLeftClipoff() << "\tRightclip: " << read.getRightClipoff();
  ostr << "\tLen: " << read.getLenClippedSeq() << " (" << (read.getRightClipoff()-read.getLeftClipoff()) << ")\tLenSeq: " << read.getLenSeq() << '\n';
  ostr << "\nLeftextend: " << read.getLeftExtend() << "\tRightextend: " << read.getRightExtend();
  ostr << "\nTemplate: ";
  if(!read.REA_template.empty()){
    ostr << read.REA_template;
  }
  ostr << "\tTemplate segment: " << static_cast<uint32>(read.REA_template_segment);
  //ostr << "\tInt. Tname: " << const_cast<Read &>(read).getInternalTemplateName();
  ostr << "\nT-ID: " << read.REA_templateid;
  ostr << "\tTPartner-ID: " << read.REA_templatepartnerid << endl;
  ostr << "RG-Info\n" << read.REA_rgid << endl;
  ostr << "RLE: " << (read.REA_rlevptr!=nullptr) << endl;

  ostr << "\nTags:\n";
  {
    Read & nonconstread=const_cast<Read &>(read);
    nonconstread.sortTags();
  }
  for(uint32 i=0; i < read.getNumOfTags(); i++) {
    const auto & acttag=read.getTag(i);
    ostr << "Tag " << i << ":\t" << acttag.getIdentifierStr() << ' ' << acttag.from << ' ' << acttag.to;
    ostr << " \"" << acttag.getCommentStr() << "\"\n";
  }

  if(read.REA_outtype!=Read::AS_TEXTCLIPS
    && read.REA_outtype!=Read::AS_TEXTSHORT){

    if(read.REA_ps_dirty==false){
      ostr << "\n\nRead size padded: " << read.REA_padded_sequence.size();
      ostr << "\nRead padded sequence:\n";
      {
	auto & seq=read.getActualSequence();
	for(const auto & se : seq) ostr << se;
	ostr << endl;
      }
    }else{
      ostr << "Forward padded dirty.\n";
    }


    if(read.REA_pcs_dirty==false){
      ostr << "\n\nRead complement size padded: " << read.REA_padded_complementsequence.size();
      ostr << "\nRead padded complement sequence:\n";

      auto & seq=read.getActualComplementSequence();
      for(const auto & se : seq) ostr << se;
      ostr << endl;
    }else{
      ostr << "Complement padded dirty.\n";
    }

    if(read.REA_outtype==Read::AS_TEXT) {
      ostr << "\nRead padded sequence, adjustment, quality, baseflags:\n";
      {
	uint8 * rleptr=nullptr;
	if(read.REA_rlevptr!=nullptr){
	  rleptr=&(read.REA_rlevptr->front());
	}

	auto sI=read.getActualSequence().cbegin();
	auto aI=read.getAdjustments().cbegin();
	auto qI=read.getQualities().cbegin();
	auto fI=read.getBPosHashStats().cbegin();

	uint32 actpos=0;
	for(;sI!=read.REA_padded_sequence.cend(); ++sI,++qI,++fI,++actpos,++rleptr){
	  ostr << actpos << ":\t" << *sI;
	  // RLE?
	  if(read.REA_rlevptr!=nullptr){
	    ostr << "(" << static_cast<uint16>(*rleptr) << ')';
	  }

	  if(!dptools::isValidStarBase(*sI)) ostr << " inv!";

	  if(read.REA_uses_adjustments) {
	    ostr << '\t' << *aI;
	    // manually counting up aI. Not in for loop as adjustments
	    //  may be empty
	    ++aI;
	  }else{
	    ostr << "\tNoAdj";
	  }
	  ostr << '\t' << static_cast<uint16>(*qI)
	       << '\t' << *fI << '\n';
	}
      }
    }
  }

  FUNCEND();
  return ostr;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpAsReadname(std::ostream & ostr)
{
  FUNCSTART("void Read::dumpAsReadname(std::ostream & ostr)");

  if(checkRead()!=nullptr || getName().empty()){
    // This is some kind of invalid read ... try to save something minimal

    if(getName().empty()){
      char buffer[24] ;
      sprintf(buffer, "%x", this) ;

      ostr << "readwithoutname_" << buffer << '\n';
    }else{
      ostr << getName() << '\n';
    }

  }else{
    ostr << getName() << '\n';
  }

  FUNCEND();
  return;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpAsFASTA(std::ostream & ostr, bool clippedonly, bool maskedseqvec, bool maskedmask)
{
  FUNCSTART("void Read::dumpAsFASTA(std::ostream & ostr, bool clippedonly, bool maskedseqvec)");

  if(checkRead()!=nullptr || getName().empty()){
    // This is some kind of invalid read ... try to save something minimal

    if(getName().empty()){
      char buffer[24] ;
      sprintf(buffer, "%x", this) ;

      ostr << ">ivld_" << buffer << "\nn\n";
    }else{
      ostr << ">ivld_" << getName() << "\nn\n";
    }

    FUNCEND();
    return;
  }

  //cout << "f:\t" << clippedonly << maskedseqvec << maskedmask << endl;

  if(getLenSeq()>0) {
    if(!clippedonly || (clippedonly && getLenClippedSeq() > 0)){

      // TODO: u.U. die ComplSeq nutzen
      refreshPaddedSequence();
      {
	ostr << ">" << getName();
	for(auto & ce : REA_tags){
	  if(ce.identifier == REA_tagentry_idCOMM) {
	    ostr << " " << ce.getCommentStr();
	  }
	}
	ostr << endl;

	auto I=REA_padded_sequence.cbegin();
	auto J=REA_padded_sequence.cend();
	int32 seqcharcount=0;
	uint64 cpl=0;
	bool dooutput;
	while(I!=J){
	  dooutput=true;
	  if(clippedonly){
	    if(seqcharcount<getLeftClipoff() || seqcharcount >= getRightClipoff()) dooutput=false;
	  }
	  if(dooutput){
	    if((maskedseqvec || maskedmask)
	       && (seqcharcount<REA_sl || seqcharcount >= REA_sr
		 || seqcharcount<REA_ml || seqcharcount >=REA_mr)) {
	      if(*I=='x'){
		ostr << "x";
	      } else {
		ostr << "X";
	      }
	    } else {
	      ostr << *I;
	    }
	    if(++cpl==REA_outlen){
	      cpl=0;
	      ostr << '\n';
	    }
	  }
	  ++I;
	  ++seqcharcount;
	}
	if(cpl!=0) ostr << '\n';
      }
    }
  }

  FUNCEND();

  return;
}


/*************************************************************************
 *
 * Output as FASTQ
 *
 * If length of (clipped/unclipped) sequence to output is 0,
 * writes a read consisting of a single N
 * TODO: adapt FASTA output to same behaviour
 *
 *************************************************************************/

void Read::dumpAsFASTQ(std::ostream & ostr, bool clippedonly, bool maskedseqvec, bool maskedmask)
{
  FUNCSTART("void Read::dumpAsFASTQ(std::ostream & ostr, bool clippedonly, bool maskedseqvec)");

  if(checkRead()!=nullptr || getName().empty()){
    // This is some kind of invalid read ... try to save something minimal

    if(getName().empty()){
      char buffer[24] ;
      sprintf(buffer, "%x", this) ;

      ostr << "@ivld_" << buffer << "\nn\n";
    }else{
      ostr << "@ivld_" << getName() << "\nn\n";
    }

    FUNCEND();
    return;
  }

  //cout << "f:\t" << clippedonly << maskedseqvec << maskedmask << endl;

  auto dolen=getLenSeq();
  if(clippedonly) dolen=getLenClippedSeq();

  ostr << '@' << getName();
  for(auto & ce : REA_tags){
    if(ce.identifier == REA_tagentry_idCOMM) {
      ostr << ' ' << ce.getCommentStr();
    }
  }
  ostr << '\n';

  if(dolen==0){
    ostr << "N\n+\nA\n";
  }else{
    // TODO: u.U. die ComplSeq nutzen
    refreshPaddedSequence();
    {
      auto I=REA_padded_sequence.cbegin();
      auto J=REA_padded_sequence.cend();

      int32 seqcharcount=0;
      bool dooutput;
      while(I!=J){
	dooutput=true;
	if(clippedonly){
	  if(seqcharcount<getLeftClipoff() || seqcharcount >= getRightClipoff()) dooutput=false;
	}
	if(dooutput){
	  if((maskedseqvec || maskedmask)
	     && (seqcharcount<REA_sl || seqcharcount >= REA_sr
		 || seqcharcount<REA_ml || seqcharcount >=REA_mr)) {
	    if(*I=='x'){
	      ostr << 'x';
	    } else {
	      ostr << 'X';
	    }
	  } else {
	    ostr << *I;
	  }
	}
	++I;
	++seqcharcount;
      }
    }
    {
      ostr << "\n+\n";
      auto I=REA_qualities.cbegin();
      auto J=REA_qualities.cend();

      int32 seqcharcount=0;
      bool dooutput;
      while(I!=J){
	dooutput=true;
	if(clippedonly){
	  if(seqcharcount<getLeftClipoff() || seqcharcount >= getRightClipoff()) dooutput=false;
	}
	if(dooutput){
	  ostr << static_cast<char>((*I)+33);
	}
	++I;
	++seqcharcount;
      }
    }
    ostr << '\n';
  }

  FUNCEND();

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpAsFASTAQual(std::ostream & ostr, bool clippedonly, bool maskedseqvec, bool maskedmask)
{
  FUNCSTART("void Read::dumpAsFASTAQual(std::ostream & ostr, bool clippedonly, bool maskedseqvec)");

  ostr << std::dec;

  if(checkRead()!=nullptr || getName().empty()){
    // This is some kind of invalid read ... try to save something minimal

    if(getName().empty()){
      char buffer[24] ;
      sprintf(buffer, "%x", this) ;

      ostr << ">ivld_" << buffer << "\n0\n";
    }else{
      ostr << ">ivld_" << getName() << "\n0\n";
    }

    FUNCEND();
    return;
  }

  //cout << "fq:\t" << clippedonly << maskedseqvec << maskedmask << endl;

  if(getLenSeq()>0) {
    if(!clippedonly || (clippedonly && getLenClippedSeq() > 0)){

      {
	ostr << ">" << getName() << endl;
	auto I=REA_qualities.cbegin();
	auto J=REA_qualities.cend();

	int32 seqcharcount=0;
	uint64 cpl=0;
	bool dooutput;
	while(I!=J){
	  dooutput=true;
	  if(clippedonly){
	    if(seqcharcount<getLeftClipoff() || seqcharcount >= getRightClipoff()) dooutput=false;
	  }
	  if(dooutput){
	    if((maskedseqvec || maskedmask)
	       && (seqcharcount<REA_sl || seqcharcount >= REA_sr
		 || seqcharcount<REA_ml || seqcharcount >=REA_mr)) {
	      ostr << "00 ";
	    } else {
	      ostr << static_cast<uint16>(*I) << ' ';
	    }
	    if(++cpl==REA_outlen){
	      cpl=0;
	      ostr << '\n';
	    }
	  }
	  ++I;
	  ++seqcharcount;
	}
	if(cpl!=0) ostr << '\n';
      }
    }
  }

  FUNCEND();

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpAsCAF(std::ostream & ostr)
{
  FUNCSTART("void Read::dumpAsCAF(std::ostream & ostr)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  if(checkRead()!=nullptr || getName().empty()){
    if(getName().empty()){
      char buffer[24] ;
      sprintf(buffer, "%x", this) ;
      ostr << "Sequence : ivld_" << buffer;
    }else{
      ostr << "Sequence : ivld_" << getName();
    }
    ostr << "Is_read\n\n";
    FUNCEND();
    return;
  }

  // TODO: u.U. die ComplSeq nutzen
  refreshPaddedSequence();
  {
    ostr << "\nDNA : " << getName() << '\n';
    auto I=REA_padded_sequence.cbegin();
    uint32 cpl=0;
    while(I!=REA_padded_sequence.cend()){
      if(*I=='*'){
	ostr << '-';
      }else{
	ostr << *I;
      }
      if(cpl++==59){
	cpl=0;
	ostr << '\n';
      }
      ++I;
    }
    ostr << '\n';
  }

#ifdef PARANOIABUGTRACKFLAG
  if(checkRead()!=nullptr){
    cout << "Ouch, failed second check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }
#endif

  {
    ostr << "\nBaseQuality : " << getName() << '\n';

    auto I=REA_qualities.cbegin();
    uint32 cpl=0;
    while(I!=REA_qualities.cend()){
      ostr << static_cast<uint16>(*I);
      if(cpl++==30){
	ostr << '\n';
	cpl=0;
      }else{
	ostr << ' ';
      }
      ++I;
    }
    ostr << "\n\n";
  }


  ostr << "Sequence : " << getName() << '\n';
  ostr << "Is_read\nPadded\n";

  if(!REA_rgid.getSeqVecName().empty()){
    ostr << "Sequencing_vector \"" << REA_rgid.getSeqVecName() << "\"" << '\n';
  }
  if(!REA_template.empty()){
    ostr << "Template \"" << REA_template <<  "\"" << '\n';
  }
  if(REA_template_segment != 0){
    if(REA_rgid.getSequencingType() != ReadGroupLib::SEQTYPE_SOLEXA
       || REA_templatepartnerid!=-1){
      if(REA_template_segment == 1) {
	ostr << "Strand Forward\n";
      }else if(REA_template_segment == 255) {
	ostr << "Strand Reverse\n";
      }
    }
  }

  if(REA_rgid.getInsizeFrom() != -1
     && REA_rgid.getInsizeTo() != -1){
    ostr << "Insert_size " << REA_rgid.getInsizeFrom()
	 << ' ' << REA_rgid.getInsizeTo() << '\n';
  }
  // TODO: readgroups. SCF filename?
  //if(!REA_scf_filename.empty()){
  //  ostr << "SCF_File \"" << REA_scf_filename <<  "\"" << '\n';
  //}else if(getSequencingType()==ReadGroupLib::SEQTYPE_454GS20){
  //  ostr << "SCF_File \"" << getName() <<  "\"" << '\n';
  //}

  if(!REA_rgid.getBaseCaller().empty()){
    ostr << "Base_caller \"" << REA_rgid.getBaseCaller() <<  "\"" << '\n';
  }

  if(!REA_sc_asped.getEntry(REA_asped).empty()){
    ostr << "Asped \"" << REA_sc_asped.getEntry(REA_asped) <<  "\"" << '\n';
  }
  if(!REA_sc_processstatus.getEntry(REA_processstatus).empty()){
    ostr << "ProcessStatus \"" << REA_sc_processstatus.getEntry(REA_processstatus) <<  "\"" << '\n';
  }
  if(!REA_rgid.getDye().empty()){
    ostr << "Dye \"" << REA_rgid.getDye() <<  "\"" << '\n';
  }
  if(!REA_rgid.getPrimer().empty()){
    ostr << "Primer " << REA_rgid.getPrimer() << '\n';
  }


  if(REA_sl>0){
    ostr << "Seq_vec SVEC 1 " << REA_sl;
    if(!REA_rgid.getSeqVecName().empty()) ostr << " \"" << REA_rgid.getSeqVecName() << "\"";
    ostr << '\n';;
  }
  if(REA_sr < static_cast<int32>(getLenSeq())){
    ostr << "Seq_vec SVEC "<< REA_sr+1 <<' ' << getLenSeq();
    if(!REA_rgid.getSeqVecName().empty()) ostr << " \"" << REA_rgid.getSeqVecName() << "\"";
    ostr << '\n';
  }

  ostr << "Clipping QUAL " << REA_ql+1 << ' ' << REA_qr << '\n'; //<<"NameOfSeq\n";

  if(REA_cl>0){
    ostr << "Clone_vec CVEC 1 " << REA_cl;
    if(!REA_rgid.getCloneVecName().empty()) ostr << " \"" << REA_rgid.getCloneVecName() <<  "\"";
    ostr << '\n';
  }
  if(REA_cr < static_cast<int32>(getLenSeq())){
    ostr << "Clone_vec CVEC "<< REA_cr+1 <<' ' << getLenSeq();
    if(!REA_rgid.getCloneVecName().empty()) ostr << " \"" << REA_rgid.getCloneVecName() <<  "\"";
    ostr << '\n';
  }

//  if(REA_cl > 0){
//    ostr << "Clone_vec CVEC "<< REA_cl+1 <<' ' << REA_cr <<  ' ' << tmpexp.getCF() << '\n';;
//  }

  //  ostr << "// Need the align_to_scf\n";

  if(REA_uses_adjustments){
    uint32 i;
    uint32 imin=0;
    for(i=1; i< getLenSeq(); i++){
      if(REA_adjustments[i-1]!=-1){
	if(REA_adjustments[i]!=REA_adjustments[i-1]+1){
	  ostr << "Align_to_SCF " << imin+1 << ' ' << i << ' ' << REA_adjustments[imin]+1 << ' ' << REA_adjustments[i-1]+1 << '\n';
	  imin=i;
	}
      }else{
	imin=i;
      }
    }
    if(REA_adjustments[i-1]!=-1){
      ostr << "Align_to_SCF " << imin+1 << ' ' << i << ' ' << REA_adjustments[imin]+1 << ' ' << REA_adjustments[i-1]+1 << '\n';
    }
  }else{
    ostr << "Align_to_SCF " << 1 << ' ' << getLenSeq() << ' ' << 1 << ' ' << getLenSeq() << '\n';
  }

  // Tags
  for(auto & rte : REA_tags) rte.dumpAsCAF(ostr);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  // Store read type and strain name in a special tag
  {
    ostr << "Tag MIT2 1 1 \"st=";
    switch(getSequencingType()) {
    case ReadGroupLib::SEQTYPE_SANGER :
    case ReadGroupLib::SEQTYPE_454GS20 :
    case ReadGroupLib::SEQTYPE_IONTORRENT :
    case ReadGroupLib::SEQTYPE_PACBIOLQ :
    case ReadGroupLib::SEQTYPE_PACBIOHQ :
    case ReadGroupLib::SEQTYPE_TEXT :
    case ReadGroupLib::SEQTYPE_SOLEXA : {
      ostr << REA_rgid.getNameOfSequencingType();
      break;
    }
    case ReadGroupLib::SEQTYPE_ABISOLID : {
      MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 7.");
      break;
    }
    default : {
      cerr << "For read " << getName() << ": unknown seqtype code " << static_cast<uint16>(getSequencingType()) << " ?\n";
      cerr << "Allowed are: Sanger, 454, IonTor, PacBioLQ, PacBioHQ, Text, Solexa, SOLiD.\n";
      MIRANOTIFY(Notify::FATAL, "Illegal readtype?");
    }
    }
    //if(REA_strainid > 0){
    if(!REA_rgid.getStrainName().empty()){
      ostr << ";sn=" << REA_rgid.getStrainName();
    }
    if(!REA_rgid.getMachineType().empty()){
      ostr << ";mt=" << REA_rgid.getMachineType();
    }
    if(isBackbone()){
      ostr << ";bb=1";
    }
    if(isRail()){
      ostr << ";rr=1";
    }
    if(isCoverageEquivalentRead()){
      ostr << ";cer=1";
    }
    if(REA_rgid.getSegmentPlacementCode()!=ReadGroupLib::SPLACE_UNKNOWN){
      if(REA_rgid.getSegmentPlacementCode()==ReadGroupLib::SPLACE_SF){
	ostr << ";pc=sf";
      }else if(REA_rgid.getSegmentPlacementCode()==ReadGroupLib::SPLACE_SB){
	ostr << ";pc=sb";
      }else if(REA_rgid.getSegmentPlacementCode()==ReadGroupLib::SPLACE_SU){
	ostr << ";pc=su";
      }else if(REA_rgid.getSegmentPlacementCode()==ReadGroupLib::SPLACE_FR){
	ostr << ";pc=fr";
      }else if(REA_rgid.getSegmentPlacementCode()==ReadGroupLib::SPLACE_RF){
	ostr << ";pc=rf";
      }else{
	BUGIFTHROW(true,"Oooops? Illegal teplate placement code " << static_cast<int32>(REA_rgid.getSegmentPlacementCode()) << " in read " << getName());
      }
    }
    ostr << "\"\n";
  }

  ostr << "\n\n";


  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpAsBAF(std::ostream & ostr)
{
  FUNCSTART("void Read::dumpAsBAF(std::ostream & ostr)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  if(checkRead()!=nullptr){
    cout << "Ouch, failed first check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }

  /*
    To be written outside:
    AP=
    DR=
    SI= & SS=
    SO=

    Not written by class yet:
    AL=


    Problems:
    - need hard clips (SL/SR)
    - SO as definition and used in reads: bad
    - RD/CO in annotations: bad, double RD from definition)
    - RD/CO in annotations: bad, consensus annotations? AT= (attach to)
    - direction of annotation
    - multiline annotations in TX?

    Perhaps include???
    optional: LN= (for preassigned sizes)


    Propose:
    DR=+/-/=
   */


  // TODO: u.U. die ComplSeq nutzen
  refreshPaddedSequence();

  ostr << "RD=" << getName()
       << "\nSQ=";

  for(const auto & pse : REA_padded_sequence){
    ostr << pse;
  }

#ifdef PARANOIABUGTRACKFLAG
  if(checkRead()!=nullptr){
    cout << "Ouch, failed second check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }
#endif

  ostr << "\nFQ=";
  for(const auto & qe : REA_qualities){
    ostr << qe+33;
  }
  ostr << '\n';

  // TODO: readgroups. SCF filename?
  //if(!REA_scf_filename.empty()){
  //  ostr << "TR=" << REA_scf_filename << '\n';
  //}else if(getSequencingType()==ReadGroupLib::SEQTYPE_454GS20){
  //  ostr << "TR=" << getName() << '\n';
  //}

  // for the time being, we can't use REA_ql an qr as MIRA
  //  has other clips that may meddle
  if(getLeftClipoff()>0) ostr << "QL=" << getLeftClipoff() << '\n';
  if(getRightClipoff()<getLenSeq()) ostr << "QR=" << getRightClipoff()+1 << '\n';

  if(!REA_template.empty()){
    ostr << "TN=" << REA_template << '\n';
  }


  if(!REA_tags.empty()){
    for(const auto & te : REA_tags){
      ostr << "AN=" << te.getIdentifierStr()
	   << "\nLO=" << (te.from)+1;
      if((te.to)-(te.from) > 0) ostr << "\nLL=" << (te.to)-(te.from)+1 << "\n" ;
      if(te.getStrandDirection()==-1) {
	ostr << "DR=-1\n";
      }

      // use serialiseComment()?
      const std::string & tmpcs=te.getCommentStr();
      if(!tmpcs.empty()){
	ostr << "TX=";
	for(uint32 i=0; i< tmpcs.size(); i++){
	  if(tmpcs[i]!='\n'){
	    ostr << tmpcs[i];
	  }else{
	    ostr << " :: ";
	  }
	}
	ostr << "\n";
      }
    }
  }

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  // Store read type and strain name in a special tag
  {
    ostr << "AN=MIT2\nTX=st=";
    switch(getSequencingType()) {
    case ReadGroupLib::SEQTYPE_SANGER :
    case ReadGroupLib::SEQTYPE_454GS20 :
    case ReadGroupLib::SEQTYPE_IONTORRENT :
    case ReadGroupLib::SEQTYPE_PACBIOHQ :
    case ReadGroupLib::SEQTYPE_PACBIOLQ :
    case ReadGroupLib::SEQTYPE_TEXT :
    case ReadGroupLib::SEQTYPE_SOLEXA : {
      ostr << REA_rgid.getNameOfSequencingType() << '\n';
      break;
    }
    case ReadGroupLib::SEQTYPE_ABISOLID : {
      MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 7.");
      break;
    }
    default : {
      cerr << "For read " << getName() << ": unknown readtype code " << static_cast<uint16>(getSequencingType()) << " ?\n";
      cerr << "Allowed are: Sanger, 454, IonTor, Solexa, PacBioLQ, PacBioHQ, Text, SOLiD.\n";
      MIRANOTIFY(Notify::FATAL, "Illegal readtype?");
    }
    }
    //if(REA_strainid > 0){
    if(!REA_rgid.getStrainName().empty()){
      ostr << ";sn=" << REA_rgid.getStrainName();
    }
    if(!REA_rgid.getMachineType().empty()){
      ostr << ";mt=" << REA_rgid.getMachineType();
    }
    if(isBackbone()){
      ostr << ";bb=1";
    }
    if(isRail()){
      ostr << ";rr=1";
    }
    if(isCoverageEquivalentRead()){
      ostr << ";cer=1";
    }
    ostr << "\n";
  }

  ostr << "\n\n";


  FUNCEND();
}



/*************************************************************************
 *
 * MAF is mixture of CAF & EXP, but single line per entry and adapted to MIRA
 *
 *
 *************************************************************************/

void Read::dumpAsMAF(std::ostream & ostr)
{
  FUNCSTART("void Read::dumpAsMAF(std::ostream & ostr)");

  if(checkRead()!=nullptr || getName().empty()){
    // This is some kind of invalid read ... try to save something minimal

    if(getName().empty()){
      char buffer[24] ;
      sprintf(buffer, "%x", this) ;

      ostr << "RD\tivld_" << buffer << "\nER\n";
    }else{
      ostr << "RD\tivld_" << getName() << "\nER\n";
    }

    FUNCEND();
    return;
  }

  // TODO: u.U. die ComplSeq nutzen
  refreshPaddedSequence();

  // RD = Read name
  // RG = Read Group ID
  // LR = Read length (length of unclipped sequence) (optional, only written on len >=2000)
  // RS = sequence (in one line, gaps as '*')
  // RQ = base quality (in one line, FASTQ-33 format)
  // SV = sequencing vector (name)
  // TN = template (name)
  // TS = template segment

  // SF = Sequencing File (file name with raw data like SCF file etc.pp)

  // SL, SR = sequencing vector left / right cutoffs
  // QL, QR = quality left / right cutoffs
  // CL, CR = (other) clipping left / right cutoffs

  // AO = Align to Original (is align_to_scf from CAF: from to from to)

  // RT = Read Tag (identifier from to comment).

  // ER = End Read (marker for MAF parsing, mandatory)
  {
    ostr << "RD\t" << getName() << '\n';
    ostr << "RG\t" << getReadGroupID().getLibId() << '\n';

    if(getLenSeq()>=2000) ostr << "LR\t" << getLenSeq() << '\n';
    ostr << "RS\t";
    for(auto & cv : REA_padded_sequence) ostr << cv;
    ostr << '\n';
  }

#ifdef PARANOIABUGTRACKFLAG
  if(checkRead()!=nullptr){
    cout << "Ouch, failed second check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }
#endif

  {
    ostr << "RQ\t";
    for(auto & qv: REA_qualities) ostr << static_cast<char>(qv+33);
    ostr << '\n';
  }

  if(!REA_rgid.getSeqVecName().empty()){
    ostr << "SV\t" << REA_rgid.getSeqVecName() << '\n';
  }
  if(!REA_template.empty()){
    ostr << "TN\t" << REA_template << '\n';
  }
  if(REA_template_segment != 0){
    ostr << "TS\t" << static_cast<uint32>(REA_template_segment) << '\n';
  }

  if(REA_sl>0){
    ostr << "SL\t" << REA_sl+1 << '\n';
  }
  if(REA_sr < static_cast<int32>(getLenSeq())){
    ostr << "SR\t"<< REA_sr  << '\n';
  }

  if(REA_ql>0){
    ostr << "QL\t" << REA_ql+1 << '\n';
  }
  if(REA_qr < static_cast<int32>(getLenSeq())){
    ostr << "QR\t"<< REA_qr  << '\n';
  }

  if(REA_cl>0){
    ostr << "CL\t" << REA_cl+1 << '\n';
  }
  if(REA_cr < static_cast<int32>(getLenSeq())){
    ostr << "CR\t"<< REA_cr  << '\n';
  }

  if(REA_uses_adjustments){
    uint32 i;
    uint32 imin=0;
    bool wroteAO=false;
    for(i=1; i< getLenSeq(); i++){
      if(REA_adjustments[i-1]!=-1){
	if(REA_adjustments[i]!=REA_adjustments[i-1]+1){
	  ostr << "AO\t" << imin+1 << '\t' << i << '\t' << REA_adjustments[imin]+1 << '\t' << REA_adjustments[i-1]+1 << '\n';
	  wroteAO=true;
	  imin=i;
	}
      }else{
	imin=i;
      }
    }
    // we're at the end of the sequence
    // if we previously wrote an AO line, there were edits and we need to write the last
    //  AO range too
    // but if AO was not written until now, do not write any
    if(REA_adjustments[i-1]!=-1
       && wroteAO){
      ostr << "AO\t" << imin+1 << '\t' << i << '\t' << REA_adjustments[imin]+1 << '\t' << REA_adjustments[i-1]+1 << '\n';
    }
  }

  for(auto & rt : REA_tags) rt.dumpAsMAF(ostr,"RT");

  ostr << "ER\n";

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::dumpAsACE(std::ostream & ostr, int32 direction)
{
  FUNCSTART("void Read::dumpAsACE(std::ostream & ostr, int32 direction)");

  if(checkRead()!=nullptr){
    cout << "Ouch, failed first check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }

  {
    decltype(REA_padded_sequence.cend()) I;
    decltype(REA_padded_sequence.cend()) IEnd;
    if (direction > 0) {
      refreshPaddedSequence();
      ostr << "\nRD " << getName() << ' ' << REA_padded_sequence.size() << " 0 0" << endl;
      I=REA_padded_sequence.cbegin();
      IEnd=REA_padded_sequence.cend();
    } else {
      refreshPaddedComplementSequence();
      ostr << "\nRD " << getName() << ' ' << REA_padded_complementsequence.size() << " 0 0" << endl;
      I=REA_padded_complementsequence.cbegin();
      IEnd=REA_padded_complementsequence.cend();
    }

    uint32 cpl=0;
    for(;I!=IEnd; ++I){
      ostr << *I;
      if(cpl++==59){
	cpl=0;
	ostr << '\n';
      }
    }
    ostr << '\n';
  }

  if(checkRead()!=nullptr){
    cout << "Ouch, failed second check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }

  // ACE has no base qualities *sigh*

  if (direction > 0) {
    ostr << "\nQA " << getLeftClipoff()+1 << ' ' << getRightClipoff();
    ostr << ' ' << getLeftClipoff()+1 << ' ' << getRightClipoff() << endl;
  } else {
    ostr << "\nQA " << REA_padded_complementsequence.size()-getRightClipoff()+1 << ' ' << REA_padded_complementsequence.size()-getLeftClipoff();
    ostr << ' ' << REA_padded_complementsequence.size()-getRightClipoff()+1 << ' ' << REA_padded_complementsequence.size()-getLeftClipoff() << endl;
  }

  ostr << "\nDS ";
  if(!REA_template.empty()){
    ostr << "TEMPLATE: " << REA_template << ' ';
  }
  // TODO: readgroups. SCF filename?
  //if(!REA_scf_filename.empty()){
  //  ostr << "SCF_FILE: " << REA_scf_filename <<  ' ' ;
  //}
  ostr << "CHROMAT_FILE: " << getName() << ' ';
  ostr << "PHD_FILE: " << getName() << ".phd.1 ";
  ostr << "TIME: Sat Jan  1 11:11:11 MEST 2000";
  ostr << '\n';

  // TODO: consed does not support tags in reverse direction
  //  MIRA sets tags always in forward direction or on both
  //  strands ('='), but tags read in from GenBank files
  //  NEED to allow for reverse (genes can be in reverse direction
  //  on the genome after all).
  // So, this is something that consed must fix.
  if(!REA_tags.empty()){
    std::string serialc;
    for(const auto & te : REA_tags){
      ostr << "RT{\n" << getName() << ' ' << te.getIdentifierStr();
      if(te.getStrandDirection()==-1){
	ostr << " MIRA " << (te.to)+1 << ' ' << (te.from)+1;
      } else{
	ostr << " MIRA " << (te.from)+1 << ' ' << (te.to)+1;
      }
      // TODO: make it a real date
      ostr << " 020202:121212\n";

      const std::string & tmpcs=te.getCommentStr();
      if(!tmpcs.empty()){
	te.serialiseComment(serialc);
	if(!serialc.empty()) {
	  ostr << "COMMENT{\n" << serialc << "\nC}\n";
	}
      }
      ostr << "}\n\n";
    }
  }else{
    ostr << '\n';
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::dumpAsGAP4DA(std::ostream & ostr, std::string & APline, bool outputTags)
{
  FUNCSTART("void Read::dumpAsGAP4DA(std::ostream & ostr, std::string & APline, bool outputTags)");

  if(checkRead()!=nullptr){
    cout << "Ouch, failed first check." << endl;
    setCoutType(AS_TEXT);
    cout << *this;
    MIRANOTIFY(Notify::FATAL, checkRead());
  }

  ostr << "ID   " << getName() << endl;
  ostr << "EN   " << getName() << endl;
  ostr << "DR   +\n";
  if(REA_ql>0) {
    ostr << "QL   " << REA_ql << endl;
  }
  if(REA_qr<static_cast<int32>(REA_qualities.size())) {
    ostr << "QR   " << REA_qr+1 << endl;
  }
  if(REA_sl>0) {
    ostr << "SL   " << REA_sl << endl;
  }
  if(REA_sr<static_cast<int32>(REA_qualities.size())) {
    ostr << "SR   " << REA_sr+1 << endl;
  }
  if(!REA_rgid.getSeqVecName().empty()){
    ostr << "SV   " << REA_rgid.getSeqVecName() << endl;
  }
  if(!REA_rgid.getCloneVecName().empty()){
    ostr << "CV   " << REA_rgid.getCloneVecName() << endl;
  }

  if(!REA_template.empty()){
    ostr << "TN   " << REA_template << endl;
  }
  if(REA_rgid.getInsizeFrom() != -1
     && REA_rgid.getInsizeTo() != -1){
    ostr << "SI   " << REA_rgid.getInsizeFrom()
	 << ".." << REA_rgid.getInsizeTo() << endl;
  }
  // TODO: readgroups. SCF filename?
  //if(!REA_scf_filename.empty()){
  //  ostr << "LN   " << REA_scf_filename;
  //  ostr << "\nLT   SCF" << endl;
  //}
  if(!REA_rgid.getBaseCaller().empty()){
    ostr << "BC   " << REA_rgid.getBaseCaller() << endl;
  }
  if(!REA_sc_processstatus.getEntry(REA_processstatus).empty()){
    ostr << "PS   " << REA_sc_processstatus.getEntry(REA_processstatus) << endl;
  }
  if(!REA_rgid.getPrimer().empty()){
    ostr << "PN   " << REA_rgid.getPrimer() << endl;
  }

  // ON
  // Basically, that's the align to scf from CAF
  // except it is not written if the read has no adjustments
  if(REA_uses_adjustments){
    ostr << "ON   ";
    bool hasrange=false;
    uint32 i;
    int32 rmin=-1;

    // James said there is a line limit for EXP
    //  so I need this to make an output of at most
    //  "8" (may be a few more) values per line
    bool mustwriteON=false;
    uint32 valuesinline=0;
    uint32 maxvaluesinline=8;

    for(i=0; i< getLenSeq(); ++i){
      if(REA_adjustments[i]==-1){
	if(mustwriteON){
	  ostr << "\nON        ";
	  valuesinline=0;
	  mustwriteON=false;
	}
	if(hasrange) {
	  ostr << rmin+1 << ".." << REA_adjustments[i-1]+1 << ' ';
	  valuesinline++;
	  hasrange=false;
	}
	ostr << "0 ";
	++valuesinline;
      } else {
	if(hasrange) {
	  if(REA_adjustments[i]!=REA_adjustments[i-1]+1) {
	    if(mustwriteON){
	      ostr << "\nON        ";
	      valuesinline=0;
	      mustwriteON=false;
	    }
	    ostr << rmin+1 << ".." << REA_adjustments[i-1]+1 << ' ';
	    valuesinline++;
	    rmin=REA_adjustments[i];
	  }
	} else {
	  rmin=REA_adjustments[i];
	  hasrange=true;
	}
      }
      if(valuesinline>=maxvaluesinline) {
	mustwriteON=true;
      }
    }
    if(hasrange) {
      ostr << rmin+1 << ".." << REA_adjustments[i-1]+1 << ' ';
    }
    ostr << endl;
  }

  // AV
  if(!REA_qualities.empty()){
    ostr << "AV   ";
    uint32 cpl=21;
    for(const auto & qe : REA_qualities){
      ostr << static_cast<uint16>(qe);
      if(--cpl){
	ostr << ' ';
      }else{
	ostr << "\nAV        ";
	cpl=19;
      }
    }
    if(cpl) ostr << '\n';
  }

  if(APline.size()) ostr << "AP   " << APline << endl;

  // SQ
  {
    refreshPaddedSequence();
    ostr << "SQ\n";
    uint32 cpl=0;
    for(const auto & pse : REA_padded_sequence){
      if(cpl==0) ostr << "    ";
      if(cpl%10==0) ostr << ' ';
      if(pse=='N' || pse=='n') {
	ostr << "-";
      }else{
	ostr << pse;
      }
      if(cpl++==59){
	cpl=0;
	ostr << '\n';
      }
    }
    if(cpl) ostr << '\n';
    ostr << "//\n";
  }

  // TG
  if(outputTags && !REA_tags.empty()){
    for(const auto & te : REA_tags){
      ostr << "TG   " << te.getIdentifierStr();
      ostr << ' ' << te.getStrand() << ' ' << (te.from)+1 << ".." << (te.to)+1;
      const char * cs=te.getCommentStr().c_str();
      if (*cs!=0) ostr << "\nTG        ";
      while(*cs) {
	ostr << *cs;
	if(*cs == '\n') ostr << "TG        ";
	++cs;
      }
      ostr << endl;
    }
  }

  FUNCEND();
  return;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::dumpTagsAsGFF3(std::ostream & ostr)
{
  FUNCSTART("void Read::dumpTagsAsGFF3(std::ostream & ostr)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  for(const auto & te : REA_tags){
    te.dumpAsGFF3(ostr,getName());
  }

  // create MIT2 info for that sequence
  ostr << getName() << "\tMIRA\texperimental_feature\t1\t1\t.\t.\t.\tst=" << REA_rgid.getNameOfSequencingType();
  if(!REA_rgid.getStrainName().empty()){
    ostr << ";sn=" << REA_rgid.getStrainName();
  }
  if(!REA_rgid.getMachineType().empty()){
    ostr << ";mt=" << REA_rgid.getMachineType();
  }
  if(isBackbone()){
    ostr << ";bb=1";
  }
  if(isRail()){
    ostr << ";rr=1";
  }
  if(isCoverageEquivalentRead()){
    ostr << ";cer=1";
  }
  ostr << ";miraitag=MIT2\n";

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// frees everything possible, making the class ready for reuse
void Read::discard()
{
  FUNCSTART("Read::discard()");

  nukeSTLContainer(REA_template);

  nukeSTLContainer(REA_padded_sequence);
  nukeSTLContainer(REA_padded_complementsequence);
  nukeSTLContainer(REA_qualities);
  nukeSTLContainer(REA_adjustments);
  nukeSTLContainer(REA_bposhashstats);
  nukeSTLContainer(REA_tags);

  if(REA_rlevptr!=nullptr){
    delete REA_rlevptr;
    REA_rlevptr=nullptr;
  }

  zeroVars();

  FUNCEND();
}





/*************************************************************************
 *
 * starting from the actual quality clips, clip more from the ends while
 *  the characters occurring are X, N, *. That is, results will differ
 *  whether or not a previous quality clipping might have happened or not
 *
 * stretches with other bases will be merged while being <= 'gapsize',
 *  the very first stretch may be 'frontgs' away and the last
 *  'endgs' away from the ends
 *
 *************************************************************************/
void Read::setClipoffsToMaskedChars(int32 gapsize, int32 frontgs, int32 endgs, bool allowN)
{
  FUNCSTART("void Read::setClipoffstoMaskedChars(int32 gapsize, int32 frontgs, int32 endgs, bool allowN)");

  BUGIFTHROW(gapsize < 0, getName() << ": gapsize (" << gapsize << ") < 0 ?");

  refreshPaddedSequence();

  if(frontgs < 0) frontgs=gapsize;
  if(endgs < 0) endgs=gapsize;

  int32 grace;

  // calc the left clip
  int32 icount;
  int32 ql;

  // 1st loop: calc from ends of read
  // 2nd loop: calc from current left clipoff
  for(uint32 loop=0; loop<2; loop++){
    grace=frontgs+1;
    if(loop) {
      ql=getLeftClipoff();
    } else {
      ql=0;
    }
    icount=ql+1;
    auto psI=REA_padded_sequence.cbegin();
    advance(psI,ql);
    bool inXrun=false;
    bool foundX=false;
    char actbase;
    while(grace>0 && psI!=REA_padded_sequence.end()) {
      actbase=toupper(*psI);
      if(allowN && actbase=='N') actbase='X';
      switch(actbase) {
      case '*' :
      case 'N' : {
	if(inXrun) {
	  ql=icount;
	  grace=gapsize+1;
	}
	break;
      }
      case 'X' : {
	inXrun=true;
	foundX=true;
        ql=icount;
	grace=gapsize+1;
	break;
      }
      default: {
	inXrun=false;
	--grace;
      }
      }
      ++psI;
      ++icount;
    }
    REA_ql=std::max(ql,REA_ql);
    if(foundX) REA_ml=ql;
  }

  // getLenSeq needs sane clipoffs
  updateClipoffs();


  // calc the right clip

  int32 qr;
  // 1st loop: calc from ends of read
  // 2nd loop: calc from current right clipoff
  for(uint32 loop=0; loop<2; loop++){
    grace=endgs+1;
    if(loop) {
      qr=getRightClipoff();
    } else {
      qr=getLenSeq();
    }
    icount=qr-1;
    auto psI=REA_padded_sequence.cbegin();
    if(qr > 0) advance(psI,qr-1);
    bool inXrun=false;
    bool foundX=false;
    char actbase;
    while(grace>0 && psI!=REA_padded_sequence.cbegin()) {
      actbase=toupper(*psI);
      if(allowN && actbase=='N') actbase='X';
      switch(actbase) {
      case 'N' :
      case '*' : {
	if(inXrun) {
	  qr=icount;
	  grace=gapsize+1;
	}
	break;
      }
      case 'X' : {
	inXrun=true;
	foundX=true;
        qr=icount;
	grace=gapsize+1;
	break;
      }
      default: {
	inXrun=false;
	--grace;
      }
      }
      --psI;
      --icount;
    }
    if(qr==-1) REA_qr=0;       // to check: necessary?

    REA_qr=std::min(qr,REA_qr);
    if(foundX) REA_mr=qr;
  }

  if(REA_qr<REA_ql) REA_qr=REA_ql;
  updateClipoffs();

  //if(REA_qr<REA_mr && REA_qr < REA_sr ) cout << getName() << " rextend possible " << REA_qr << ' ' <<REA_mr << ' ' <<REA_sr << endl;

  FUNCEND();

  return;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::setSequenceFromString(const char * sequence, int32 slen)
{
  FUNCSTART("void Read::setSequenceFromString(const char * sequence)");

  REA_has_valid_data=false;

  REA_padded_sequence.clear();

  // initialise the sequence vector
  {
    REA_padded_sequence.reserve(slen+3);
    REA_padded_sequence.resize(slen,0);
    auto dI=REA_padded_sequence.begin();
    const char * sptr=sequence;
    for(auto eptr= sptr+slen; sptr!=eptr; ++dI, ++sptr){
      if(likely(dptools::isValidIUPACStarBase(*sptr))){
	*dI=*sptr;
      }else{
	switch(*sptr){
	case 'u':
	case 'U': {   // read / set RNA sequence
	  *dI='T';
	  break;
	}
	case '@':      // formerly used by MIRA for "no coverage"
	case '-': {
	  *dI='N';
	  break;
	}
	// be lenient: jump over spaces, tabs, newlines, ^M from DOS, etc.pp
	case ' ' :
	case '\t' :
	case '\n' :
	case '\r' :{
	  break;
	}
	case 0 :
	default: {
	  MIRANOTIFY(Notify::FATAL,"tried to set a base '" << *sptr << "' (ASCII: " << static_cast<uint16>(*sptr) << "), which is not a valid IUPAC base nor N, X, - or @.");
	}
	}
      }
    }
  }

  REA_ql=0;
  REA_qr=static_cast<int32>(REA_padded_sequence.size());
  REA_sl=0;
  REA_sr=static_cast<int32>(REA_padded_sequence.size());
  REA_cl=0;
  REA_cr=static_cast<int32>(REA_padded_sequence.size());
  REA_ml=0;
  REA_mr=static_cast<int32>(REA_padded_sequence.size());

  postLoadEXPFASTA();

  REA_has_valid_data=true;

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::checkQualities()
{
  FUNCSTART("void Read::checkQualities()");

  for(auto & qv : REA_qualities){
    if(qv>100){
      MIRANOTIFY(Notify::FATAL,"Read " << getName() << ": have quality <0 or > 100???\n");
    }
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::setQualities(const std::vector<base_quality_t> & quals)
{
  FUNCSTART("void Read::setQuality(std::vector<base_quality_t> & quals)");

  if(quals.size() != REA_qualities.size()) {
    MIRANOTIFY(Notify::FATAL,"Read " << getName() << ": tried to set " << quals.size() << " qualities although the read has " << REA_qualities.size() << " bases.\n");
  }

  REA_qualities=quals;

  setQualityFlag(true);

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::setQualities(base_quality_t qual)
{
  FUNCSTART("void Read::setQuality(base_quality_t qual)");

  mstd::fill(REA_qualities,qual);
  setQualityFlag(true);

  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::loadDataFromSCF(const std::string & scfname)
{
  FUNCSTART("void Read::loadDataFromSCF(const std::string & scfname)");

  discard();

  if(scfname.empty()) {
    MIRANOTIFY(Notify::FATAL, "No scfname given?");
  }
  if(!hasSCFData()){
    MIRANOTIFY(Notify::FATAL, "No scf data available for " << scfname);
  }

  //REA_scf_filename=scfname;
  //setFileNamesFromSCFFileName(scfname);

  SCF scf;
  std::string emsg;
  if(scf.load(scfname.c_str(),emsg)>0){
    MIRANOTIFY(Notify::FATAL,"Unexpected error where it should not happen: " << emsg);
  }

  // initialise the sequence vector
  {
    REA_padded_sequence.clear();
    REA_padded_sequence.reserve(scf.getNumBases()+3);

    for(uint32 i=0; i<scf.getNumBases(); i++) {
      char tmpc=scf.getBase(i);
      if(dptools::isValidIUPACBase(tmpc)){
	REA_padded_sequence.push_back(tmpc);
      }else{
	switch(toupper(tmpc)){
	case 'N':
	case 'X':{
	  REA_padded_sequence.push_back(tmpc);
	  break;
	}
	case '*':
	case '-': {
	  REA_padded_sequence.push_back('N');
	  break;
	}
	default: {
	  cerr << "At base position " << i << ": base " << tmpc << " unknown IUPAC code?\n";
	  MIRANOTIFY(Notify::FATAL, "Illegal base found: " << scfname);
	}
	}
      }
    }
  }

  REA_ql=0;
  REA_qr=static_cast<int32>(REA_padded_sequence.size());
  REA_sl=0;
  REA_sr=static_cast<int32>(REA_padded_sequence.size());
  REA_cl=0;
  REA_cr=static_cast<int32>(REA_padded_sequence.size());
  REA_ml=0;
  REA_mr=static_cast<int32>(REA_padded_sequence.size());

  updateClipoffs();

  postLoadEXPFASTA();

  for(uint32 i=0; i<scf.getNumBases(); i++) {
    REA_qualities[i]=scf.getCalledBaseProb(i);
  }

  REA_has_valid_data=true;

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
  return;

}



/*************************************************************************
 *
 * void postLoadEXPFASTA()
 * build complements, build padded vectors, set qualities to defaults
 *  set adjustments to defaults, set bposhashstats to default
 *  set template if necessary and possible
 *
 *************************************************************************/
void Read::postLoadEXPFASTA()
{
  FUNCSTART("void postLoadEXPFASTA();");

  // build the complement padded sequence
  // NO, don't do that now. Defer until really needed.
  //makeComplement(REA_padded_sequence, REA_padded_complementsequence);

  REA_ps_dirty=false;
  REA_pcs_dirty=true;

  // When loading from an EXP-file, there might be no qualities given
  // So, fill the quality-vector with 0, which will be overwritten if
  //  quals are available
  REA_qualities.clear();
  REA_qualities.reserve(REA_padded_sequence.capacity());
  REA_qualities.resize(REA_padded_sequence.size(),0);

  // base flags
  REA_bposhashstats.clear();
  REA_bposhashstats.reserve(REA_padded_sequence.capacity());
  clearAllBPosHashStats();

  // When loading an EXP-file, the bases correspond to the bases in the
  //  SCF file.
  if(REA_uses_adjustments){
    REA_adjustments.reserve(REA_padded_sequence.capacity());
    REA_adjustments.resize(REA_padded_sequence.size());
    auto aI=REA_adjustments.begin();
    auto sI=REA_padded_sequence.cbegin();
    int32 actadjust=0;
    for(; aI<REA_adjustments.end(); aI++, sI++){
      if(*sI!='*') {
	*aI=actadjust;
	++actadjust;
      }else{
	*aI=-1;
      }
    }
  }

  // and don't forget that this read was possibly used before and had RLE ...
  if(REA_rlevptr!=nullptr){
    delete REA_rlevptr;
    REA_rlevptr=nullptr;
  }

  FUNCEND();
}


/*************************************************************************
 *
 * TODO; check 'frommaskedchar'
 *
 *
 *************************************************************************/

void Read::setMinimumLeftClipoff(int32 minreq, int32 setto, bool fromseqvec, bool frommaskedchar)
{
  FUNCSTART("void Read::setMinimumLeftClipoff(int32 minreq, int32 setto, bool fromseqvec, bool frommaskedchar)");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  if( REA_ql < setto ) {
    if(REA_ql < minreq) {
      REA_ql=setto;
    }
    if(fromseqvec && REA_sl < minreq) {
      REA_ql=setto;
    }
    if(REA_ml>0 && REA_ml < minreq) {
      REA_ql=setto;
    }
    if(REA_ql>=static_cast<int32>(REA_qualities.size())) REA_ql=static_cast<int32>(REA_qualities.size())-1;
    if(REA_ql<0) REA_ql=0;
    if(REA_qr<REA_ql) REA_qr=REA_ql;
    updateClipoffs();
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::setMinimumRightClipoff(int32 minreq, int32 setto)
{
  FUNCSTART("void Read::setMinimumRightClipoff(int32 minreq, int32 setto)");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  if( REA_qr > static_cast<int32>(REA_qualities.size())-minreq ) {
    REA_qr=static_cast<int32>(REA_qualities.size())-setto;
    if(REA_qr<0) REA_qr=0;
    if(REA_qr<REA_ql) REA_qr=REA_ql;
    updateClipoffs();
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Read::hasSCFData(bool loadfailnoerror)
{
  FUNCSTART("bool Read::hasSCFData(bool loadfailnoerror=false)");

  if(!REA_has_valid_data
     || getSequencingType() != ReadGroupLib::SEQTYPE_SANGER) {
    FUNCEND();
    return false;
  }

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  checkSCFAndLoadQual(hasQuality(),loadfailnoerror);

  FUNCEND();
  return REA_scf_available;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::checkSCFAndLoadQual(bool justcheck, bool loadfailnoerror)
{
  FUNCSTART("void Read::checkSCFAndLoadQual(bool justcheck)");

  if(hasValidData()==false){
    MIRANOTIFY(Notify::INTERNAL, "Trying to load SCF from unitialised Read.");
  }
  if(REA_scf_loadattempted
     || isRail()
     || isCoverageEquivalentRead()
     || getSequencingType()!=ReadGroupLib::SEQTYPE_SANGER){
    FUNCEND();
    return;
  }

  REA_scf_loadattempted=true;

  SCF scf;

  std::string dummy;
  getSCFFullPathName(dummy);
  if(dummy.empty()) return;

  std::string emsg;
  auto loadres=scf.load(dummy.c_str(),emsg);
  if(loadres>0){
    if(!loadfailnoerror){
      cout << "Warning: " << emsg << endl;
    }
    return;
  }

  scf.transposeAmbiguityCodes();

  if(justcheck==false){
    if(scf.getNumBases()!=getLenSeq()){
      cout << "Warning: Number of bases in SCF file (" << dummy << ") does not correspond to the number of bases expected in the read (read from fasta, exp, phd or caf file).\nRead will _not_ be used in assembly!" << endl;
      return;
    }
    for(uint32 i=0; i<REA_padded_sequence.size(); ++i){
      if(toupper(scf.getBase(i))!=toupper(REA_padded_sequence[i])){
	if(!(scf.getBase(i)=='-' || toupper(REA_padded_sequence[i])=='N' || toupper(REA_padded_sequence[i])=='X')){
	  cout << "Warning: A base (" << scf.getBase(i) << ") in the SCF (" << dummy << ") does not correspond to the the one (" << REA_padded_sequence[i] << ") in the read read from fasta, exp, phd or caf (position: " << i << ").\nRead will _not_ be used in assembly!" << endl;
	  return;
	}
      }
      if(justcheck==false){
	REA_qualities[i]=scf.getCalledBaseProb(i);
      }
    }
    setQualityFlag(true);
  }
  REA_scf_available=true;

  FUNCEND();
}


/*************************************************************************
 *
 * The "const"-ness of this function is a lie
 * Destination is usually a member variable, which, however, is mutable :-)
 *  (for non-mutable variables the compiler rightly complains)
 *
 * Not the most performant solution, but I don't care atm
 *
 *************************************************************************/
void Read::makeComplement(const std::vector<char> & source, std::vector<char> & destination) const
{
  FUNCSTART("void Read::makeComplement(std::vector<char> & source, std::vector<char> & destination);");

  destination.reserve(source.capacity());
  destination.resize(source.size());
  auto dI=destination.begin();
  auto sI=source.crbegin();
  for(;sI!=source.crend();++dI,++sI){
    char cbase=dptools::getComplementIUPACBase(*sI);
    if(likely(cbase!=0)){
      *dI=cbase;
    }else{
      cout << "Argh! Found illegal base " << *sI << " and there's no complement for it!\n";
      MIRANOTIFY(Notify::INTERNAL, "Illegal base in a phase where there should be none.");
    }
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

int32 Read::calcComplPos(int32 pos) const
{
  FUNCSTART("int32 Read::calcComplPos(int32 pos) const");
  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  FUNCEND();
  return static_cast<int32>(getLenSeq())-1-pos;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

int32 Read::calcClippedPos2RawPos(int32 pos) const
{
  FUNCSTART("int32 Read::calcComplPos2RawPos(int32 pos) const");
  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  FUNCEND();
  return (pos+getLeftClipoff());
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

int32 Read::calcRawPos2ClippedPos(int32 pos) const
{
  FUNCSTART("int32 Read::calcComplPos2RawPos(int32 pos) const");
  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  FUNCEND();
  return (pos-getLeftClipoff());
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

int32 Read::calcClippedComplPos2RawPos(int32 pos) const
{
  FUNCSTART("int32 Read::calcClippedComplPos2RawPos(int32 pos) const");
  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  FUNCEND();
  return (getRightClipoff()-1-pos);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::setName(const std::string & name)
{
  FUNCSTART("void Read::setName(const std::string & name)");

  const char * c=name.c_str();
  const char * emsg=nullptr;
  const char * wmsg=nullptr;
  char wchar=' ';
  for(; *c !=0; ++c){
    if(*c <32 || *c==127){
      emsg="This is a control code, names should not contain that!";
    }
    switch(*c){
    case 32 : {
      emsg="This is the space character, names should not contain that!";
      break;
    }
    case 34 : //"
    case 38 : //&
    case 39 : //'
    case 59 : //;
    case 63 : //?
    case 96 : //`
    {
      wchar=*c;
      wmsg="This character may pose problems in downstream processing by other programs, it is suggested you choose to change that name.";
      break;
    }
    default : {
    }
    }
    if(emsg != nullptr) break;
  }

  if(emsg!=nullptr){
    MIRANOTIFY(Notify::INTERNAL,"While trying to set the name of read \n" << name << "\nEncountered character with ASCII code " << static_cast<uint16>(*c) << ". " << *emsg << "\nIt is probably due to your input data, but normally, MIRA should have caught that earlier!");
  }
  if(wmsg!=nullptr){
    cout << "Warning while setting the name of read '"
	 << name
	 << "':\nEncountered character with ASCII code "
	 << static_cast<uint16>(wchar) << " (" << wchar << ").\n"
	 << *wmsg << "\n";
  }

  //std::string finalname(renameReadPrefix(name));
  if(REA_rgid.getReadRenamePrefixList().empty()){
    REA_nameentry=REA_sc_readname.addEntryNoDoubleCheck(name);
  }else{
    REA_nameentry=REA_sc_readname.addEntryNoDoubleCheck(renameReadPrefix(name));
  }

  if(!REA_rgid.isDefaultNonValidReadGroupID()) calcTemplateInfo();

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool Read::checkStringForSRANamingScheme(const std::string & rname)
{
  if(rname.size()>10
     && rname[0] >= 'A' && rname[0] <= 'Z'
     && rname[1] >= 'A' && rname[1] <= 'Z'
     && rname[2] >= 'A' && rname[2] <= 'Z'
     && isdigit(rname[3])
     && isdigit(rname[4])
     && isdigit(rname[5])
     && isdigit(rname[6])
     && isdigit(rname[7])
     && isdigit(rname[8])
     && rname[9] == '.'
    ){
    return true;
  }
  return false;
}


/*************************************************************************
 *
 * renamelist: (prefix2search replacewith) (... ...)
 *
 *
 *************************************************************************/

std::string Read::renameReadPrefix(const std::string & actname)
{
  FUNCSTART("bool Read::renameReadPrefix(std::string & actname)");

  auto & renamelist=REA_rgid.getReadRenamePrefixList();
  if(renamelist.empty()) return actname;
  BUGIFTHROW(renamelist.size()%2,"renamelist % 2 != 0???");

  for(auto rlI=renamelist.begin(); rlI!=renamelist.end(); ++rlI){
    auto rmI = std::mismatch(rlI->begin(), rlI->end(), actname.begin());
    if(rmI.first == rlI->end()){
      // found a prefix, now rename read;
      std::string newname=*(++rlI);  // new prefix
      newname+=actname.substr((--rlI)->size(),std::string::npos); // old postfix
      return newname;
    }
    ++rlI;
  }

  return actname;
}



void Read::getSCFFullPathName(std::string & path) const
{
  path.clear();

  if(!REA_rgid.getDataDir().empty()){
    path+=REA_rgid.getDataDir();
    if(path.back()!='/') path+='/';
  }
  if(!REA_rgid.getDataFile().empty()){
    path+=REA_rgid.getDataFile();
  }else{
    path+=getName()+".scf";
  }
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::initialiseRead(bool preserve_originals,
			  bool iscomplement,
			  bool ispadded,
			  ReadGroupLib::ReadGroupID rgid,
			  std::vector<char>           & sequence,
			  std::vector<base_quality_t> & qualities,
			  std::vector<int32>          & adjustments,
			  std::vector<multitag_t>     & tags,
			  const std::string & name,
			  const std::string & SCFname,
			  int32 ql, int32 qr,    // quality clipping left/right
			  int32 sl, int32 sr,   // sequence vector clipping
			  int32 cl, int32 cr   // clone vector clipping
			  )
{
  FUNCSTART("void Read::initialiseRead(bool preserve_originals, ..., ...<base_quality_t>, ...)");

  discard();

  // handle the names
  BUGIFTHROW(name.empty(), "no name given");

  REA_rgid=rgid;  // do this before calling setName()!
  setName(name);

  BUGIFTHROW(qualities.size() != sequence.size(),"qualities.size() " << qualities.size() << " != sequence.size() " << sequence.size() << " for: " << name);
  BUGIFTHROW(!adjustments.empty() && adjustments.size() != sequence.size(),"adjustments.size() " << adjustments.size() << " != sequence.size() " << sequence.size() << " for: " << name);

  // force not using adjustements if those are empty
  if(!sequence.empty() && adjustments.empty()){
    disallowAdjustments();
  }

  // handle the sequence
  if(iscomplement==true){
    if(ispadded==true){
      if(preserve_originals==true){
	REA_padded_complementsequence=sequence;
      }else{
	REA_padded_complementsequence.swap(sequence);
      }
      makeComplement(REA_padded_complementsequence, REA_padded_sequence);
    }else{
      if(preserve_originals==true){
	REA_padded_complementsequence=sequence;
      }else{
	REA_padded_complementsequence.swap(sequence);
      }
      makeComplement(REA_padded_complementsequence, REA_padded_sequence);
    }
  }else{
    if(ispadded==true){
      if(preserve_originals==true){
	REA_padded_sequence=sequence;
      }else{
	REA_padded_sequence.swap(sequence);
      }
      makeComplement(REA_padded_sequence, REA_padded_complementsequence);
    }else{
      if(preserve_originals==true){
	REA_padded_sequence=sequence;
      }else{
	REA_padded_sequence.swap(sequence);
      }
      makeComplement(REA_padded_sequence, REA_padded_complementsequence);
    }
  }

  // handle the qualities
  if(iscomplement==true){
    REA_qualities.reserve(qualities.size());
    auto I=qualities.rbegin();
    while(I!=qualities.rend()){
      REA_qualities.push_back(*I++);
    }
  }else{
    if(preserve_originals==true){
      REA_qualities=qualities;
    }else{
      REA_qualities.swap(qualities);
    }
  }
  setQualityFlag(true);

  // handle the adjustments
  if(REA_uses_adjustments){
    if(iscomplement==true){
      REA_adjustments.reserve(adjustments.size());
      auto I=adjustments.rbegin();
      while(I!=adjustments.rend()){
	REA_adjustments.push_back(*I++);
      }
    }else{
      if(preserve_originals==true){
	REA_adjustments=adjustments;
      }else{
	REA_adjustments.swap(adjustments);
      }
    }
  }

  // handle the tags
  setTags(tags);

  // handle the baseflags (not given, so set everything to default)
  clearAllBPosHashStats();


  // handle the normal clippings

  REA_qr=qr;
  REA_ql=ql;
  REA_sr=sr;
  REA_sl=sl;
  REA_cr=cr;
  REA_cl=cl;
  updateClipoffs();

  // getLenSeq needs a valid read!
  REA_has_valid_data=true;
  REA_ml=0;
  REA_mr=getLenSeq();

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::removeGapsFromRead()
{
  FUNCSTART("void Read::removeGapsFromRead()");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  refreshPaddedSequence();

  // TODO: check what happens with gaps at front of reads?
  //bool del=false;
  for(int32 i=0; i<static_cast<int32>(REA_padded_sequence.size());i++){
    if(REA_padded_sequence[i]=='*'){
      //del=true;
      deleteBaseFromSequence(i);
      i--;
    }
  }

  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::fixZeroGapQuals()
{
  FUNCSTART("void Read::fixZeroGapQuals()");

  // cannot do anything in these cases
  if(REA_qualities.empty() || getLenSeq()<2) return;

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  refreshPaddedSequence();

  for(int32 i=0; i<static_cast<int32>(REA_padded_sequence.size());i++){
    if(REA_padded_sequence[i]=='*' && REA_qualities[i]==0){
      int32 lowpos=i;
      for(; lowpos>=0 && REA_padded_sequence[lowpos]=='*'; --lowpos);
      if(lowpos<0) lowpos=0;
      int32 highpos=i;
      for(; highpos<REA_padded_sequence.size() && REA_padded_sequence[highpos]=='*'; ++highpos);
      if(highpos>=REA_padded_sequence.size()) highpos=REA_padded_sequence.size()-1;
      REA_qualities[i]=(REA_qualities[lowpos]+REA_qualities[highpos])/2;
    }
  }

  FUNCEND();
  return;
}



/*************************************************************************
 *
 * for strobed PacBio reads which have stretches like ...x*****NNNNNNNNNx...
 *
 * in stretch given, insert or delete the number of Ns given
 *
 * from, to: including! from can be > to
 * changeestim <0: delete Ns; >0 insert Ns
 *
 *************************************************************************/

void Read::correctNStretch(int32 from, int32 to, int32 changeestim)
{
  FUNCSTART("void Read::correctNStretch(const int32 from, const int32 to, int32 changeestim)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  if(from>to) std::swap(from,to);

  refreshPaddedSequence();

  BUGIFTHROW(from<0,getName() << ": from (" << from << ") <0 ?");
  BUGIFTHROW(from>=REA_padded_sequence.size(),getName() << ": from (" << from << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");
  BUGIFTHROW(to<0,getName() << ": to (" << to << ") <0 ?");
  BUGIFTHROW(to>=REA_padded_sequence.size(),getName() << ": to (" << to << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  if(changeestim==0) return;

  bool ok=true;
  for(int acti=from; acti <= to; acti++){
    if(REA_padded_sequence[acti] != '*'
       && REA_padded_sequence[acti] != 'N'
       && REA_padded_sequence[acti] != 'n'){
      ok=false;
      break;
    }
  }

  BUGIFTHROW(!ok, getName() << " from: " << from << " to: " << to << " not exclusively N or *?");

  if(changeestim>0){
    // insert Ns
    // first, try to convert gaps to N
    // TODO: slow, inefficient
    for(int acti=from; acti <= to && changeestim!=0; acti++){
      if(REA_padded_sequence[acti] == '*'){
	changeBaseInSequence('N', 1, static_cast<uint32>(acti));
	--changeestim;
      }
    }
    // if still need to insert, insert them
    // TODO: slow, inefficient
    if(changeestim!=0){
      for(; changeestim != 0; --changeestim){
	insertBaseInSequence('N', 1, static_cast<uint32>(from),false);
      }
    }
  }else{
    // delete Ns
    for(int acti=from; acti <= to && changeestim!=0; acti++){
      if(REA_padded_sequence[acti] == 'N'
	 || REA_padded_sequence[acti] == 'n'){
	deleteBaseFromSequence(static_cast<uint32>(acti));
	--to;
	++changeestim;
      }
    }
  }

  FUNCEND();
}



/*************************************************************************
 *
 * for strobed PacBio reads which have stretches like ...x*****NNNNNNNNNx...
 *
 * if minimum of (10) N in a gap/N stretch, transform all gaps of that
 *  stretch to N
 *
 *************************************************************************/

void Read::transformGapsToNs()
{
  FUNCSTART("void Read::transformGapsToNs()");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  refreshPaddedSequence();

  const uint32 minns=10;

  auto cS=REA_padded_sequence.begin();
  auto qS=REA_qualities.begin();
  for(; cS != REA_padded_sequence.cend(); ++cS, ++qS){
    if(*cS == '*'
       || *cS == 'N'
       || *cS == 'n'){
      auto cE=cS;
      uint32 ncount=0;
      for(; cE != REA_padded_sequence.cend(); ++cE){
	if(*cE == 'N' || *cE == 'n') {
	  ncount++;
	}else if(*cE != '*'){
	  break;
	}
      }
      uint32 ccount=0;
      for(; cS != cE; ++cS, ++qS){
	if(ncount>=minns
	   && *cS == '*'
	   && ccount <5 ) {
	  *cS='N';
	  *qS=1;
	  --ccount;
	}
      }
      // we might have arrived at the end of the read in the inner loop
      //  decrease by one so that the outer for will arrive at the end
      if(cS==REA_padded_sequence.cend()){
	--cS;
	--qS;
      }
    }
  }

  FUNCEND();
  return;
}



/*************************************************************************
 *
 * return the adjustment position of a read position
 *
 *************************************************************************/

int32 Read::getAdjustmentPosOfReadPos(const uint32 position) const
{
  FUNCSTART("int32 getAdjustmentPosOfReadPos(int32 position)");

  if(REA_uses_adjustments){
    BUGIFTHROW(position>=REA_adjustments.size(), getName() << ": readpos (" << position << ") >= size of read (" << REA_adjustments.size() << ")?");

    FUNCEND();
    return REA_adjustments[position];
  }

  refreshPaddedSequence();

  if(REA_padded_sequence[position]=='*') return -1;
  return position;
}


/*************************************************************************
 *
 * return the nearest adjustment position that is not a -1 of a read position
 * (searching toward lower)
 *
 *************************************************************************/

int32 Read::getLowerNonGapAdjustmentPosOfReadPos(uint32 position) const
{
  FUNCSTART("int32 Read::getLowerNonGapAdjustmentPosOfReadPos(const uint32 position) const");

  if(REA_uses_adjustments){
    BUGIFTHROW(position>=REA_adjustments.size(), getName() << ": readpos (" << position << ") >= size of read (" << REA_adjustments.size() << ")?");

    while(position>0 && REA_adjustments[position] == -1) --position;

    FUNCEND();
    return REA_adjustments[position];
  }

  refreshPaddedSequence();
  while(position>0 && REA_padded_sequence[position] == '*') --position;

  int32 adjpos=0;
  auto cI=REA_padded_sequence.cbegin();
  for(int32 i=0; i<position; ++i, ++cI){
    if(*cI!='*') ++adjpos;
  }

  FUNCEND();
  return adjpos;
}


/*************************************************************************
 *
 * return the nearest adjustment position that is not a -1 of a read position
 * (searching toward upper)
 *
 *************************************************************************/

int32 Read::getUpperNonGapAdjustmentPosOfReadPos(uint32 position) const
{
  FUNCSTART("int32 Read::getUpperNonGapAdjustmentPosOfReadPos(const uint32 position) const");

  if(REA_uses_adjustments){
    BUGIFTHROW(position>=REA_adjustments.size(), getName() << ": readpos (" << position << ") >= size of read (" << REA_adjustments.size() << ")?");

    while(position<REA_adjustments.size()-1 && REA_adjustments[position] == -1) ++position;

    FUNCEND();
    return REA_adjustments[position];
  }

  refreshPaddedSequence();
  while(position<REA_padded_sequence.size()-1 && REA_padded_sequence[position] == '*') ++position;

  int32 adjpos=0;
  auto cI=REA_padded_sequence.cbegin();
  for(int32 i=0; i<position; ++i, ++cI){
    if(*cI!='*') ++adjpos;
  }

  return adjpos;
}




/*************************************************************************
 *
 * return the nearest position that is not a gap
 * (searching toward lower)
 *
 *************************************************************************/

int32 Read::getLowerNonGapPosOfReadPos(uint32 position) const
{
  FUNCSTART("int32 Read::getLowerNonGapPosOfReadPos(const uint32 position) const");

  refreshPaddedSequence();

  BUGIFTHROW(position>=REA_padded_sequence.size(),getName() << ": position (" << position << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  while(position>0 && REA_padded_sequence[position] == '*') --position;

  FUNCEND();
  return position;
}


/*************************************************************************
 *
 * return the nearest position that is not a gap
 * (searching toward upper)
 *
 *************************************************************************/

int32 Read::getUpperNonGapPosOfReadPos(uint32 position) const
{
  FUNCSTART("int32 Read::getUpperNonGapPosOfReadPos(const uint32 position) const");

  refreshPaddedSequence();
  BUGIFTHROW(position>=REA_padded_sequence.size(),getName() << ": position (" << position << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  while(position<REA_padded_sequence.size()-1 && REA_padded_sequence[position] == '*') ++position;

  return position;
}







/*************************************************************************
 *
 * given an adjustment position, return a read position (if it exists)
 * else -1
 *
 *************************************************************************/

int32 Read::getReadPosOfAdjustmentPos(const int32 position) const
{
  FUNCSTART("int32 getReadPosOfAdjustmentPos(int32 position)");

  if(REA_uses_adjustments){
    for(auto adjI=REA_adjustments.begin(); adjI!=REA_adjustments.end(); ++adjI) {
      if(*adjI==position) {
	FUNCEND();
	return static_cast<int32>(adjI-REA_adjustments.begin());
      }
    }

    FUNCEND();
    return -1;
  }

  return position;
}


/*************************************************************************
 *
 * Not multithread capable because of the static variables,
 *  but multithread really not needed here and that's faster
 *
 * If need to be, remove "static"
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}
void Read::upgradeOldTagToMultitagWithGFF3(const multitag_t & oldtag, multitag_t & mt)
{
  std::string emptystr;
  std::string dummystr;

  std::string setcomment;
  std::string setidentifier;

  std::string miraitag;
  std::string g3source;
  char   g3strand;
  uint8  g3phase;

  bool isgff3=GFFParse::checkCommentForGFF3(oldtag.getCommentStr());

  //cout << "###Checked oldmiragbf: " << oldmiraGBFstyle << "\tisgff3: " << isgff3 << "\n" << oldtag;

  if(!isgff3 && oldtag.identifier!=REA_tagentry_idMIT2){
    CEBUG("non-gff3" << endl);

    setidentifier=oldtag.getIdentifierStr();

    bool oldmiraGBFstyle=GBF::checkIfCommentInOldMIRAGBFstyle(oldtag.getCommentStr());

    // retain compatibility to old MINF tags by doing nothing with the comment itself
    //  and allow read finction transferMINFTagsToReadInfo() to get what it wants
    if(oldmiraGBFstyle){
      GBF::liftOldMIRAGBFCommentToGFF3(oldtag.getCommentStr(),setcomment);
    }else if(!oldtag.getCommentStr().empty()){
      setcomment="Note="+oldtag.getCommentStr();
    }

    // there will be no GFF3 fields in those ... just called with empty string
    //  to make sure the g3* variables get default values
    GFFParse::extractMIRAGFF3InfoFromGFF3Attributes(emptystr,
						    dummystr,
						    g3source,
						    g3strand,
						    g3phase,
						    miraitag);
    g3strand=oldtag.getStrand();
    if(oldtag.from==oldtag.to) g3strand='.';
    g3source="MIRA";
  }else if(oldtag.getCommentStr().empty()){
    //cout << "isempty\n";
    // there will be no GFF3 fields in those ... just called with empty string
    //  to make sure the g3* variables get default values
    GFFParse::extractMIRAGFF3InfoFromGFF3Attributes(emptystr,
						    dummystr,
						    g3source,
						    g3strand,
						    g3phase,
						    miraitag);
    g3source="MIRA";
  }else{
    CEBUG("is-gff3" << endl);
    GFFParse::extractMIRAGFF3InfoFromGFF3Attributes(oldtag.getCommentStr(),
						    setcomment,
						    g3source,
						    g3strand,
						    g3phase,
						    miraitag);
    if(g3strand=='*'){
      g3strand=oldtag.getStrand();
      if(oldtag.from==oldtag.to) g3strand='.';
    }
    if(!miraitag.empty()) {
      //cout << "seen itag\n";
      setidentifier=miraitag;
    }else{
      //cout << "no itag\n";
      setidentifier=oldtag.getIdentifierStr();
    }
  }
  mt.from=oldtag.from;
  mt.to=oldtag.to;
  mt.comment=multitag_t::newComment(setcomment);
  mt.identifier=multitag_t::newIdentifier(setidentifier);
  mt.source=multitag_t::newSource(g3source);
  mt.setStrand(g3strand);
  mt.phase=g3phase;
  mt.commentisgff3=true;

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


void Read::setTags(const std::vector<multitag_t> & tags)
{
  FUNCSTART("void Read::setTags(const std::vector<multitag_t> & tags)");

  for(const auto & te : tags){
    if(te.identifier==REA_tagentry_idMINF
       || te.identifier==REA_tagentry_idMIT2){
      BUGIFTHROW(true,"Read " << getName() << ": found MINF / MIT2 tag while setting tags. With read groups, should not be anymore, must be parsed out before.");
    }
    try{
      addTagO(te);
    }
    catch (Notify n){
      cout << te;
      n.handleError(THISFUNC);
    }
  }
  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::addTagO(const multitag_t & tag)
{
  FUNCSTART("void Read::addTag(multitag_t & tag)");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  BUGIFTHROW(tag.from>=getLenSeq(),"read " << getName() << " : tag.from " << tag.from << " >= getLenSeq() " << getLenSeq() <<" ???");
  BUGIFTHROW(tag.to>=getLenSeq(),"read " << getName() << " : tag.to " << tag.to << " >= getLenSeq() " << getLenSeq() << " ???");

  // simple check for tag duplicates
  // TODO: MIRA started to use multiple base tags at this time,
  //  so this should be perhaps reworked.
  // New: if tag is already there, then backpack information (comment
  //  at the time) is copied into the location.
  //  Cannot use string.swap() as tag might be present multiple times
  //   (should not happen, but might be due to loading non-MIRA tags)

  bool mustaddtag=true;
  for(auto & te : REA_tags) {
    if(te.from == tag.from
       && te.to == tag.to
       && te.identifier == tag.identifier
       && te.getStrand() == tag.getStrand()){
      mustaddtag=false;
      te.comment=tag.comment;
    }
  }

  if(mustaddtag) {
    REA_tags.push_back(tag);
    if(tag.to<tag.from) std::swap(REA_tags.back().to,REA_tags.back().from);
  }

  BUGIFTHROW(tag.to>=getLenSeq(), "Read " << getName() << " for tag: " << tag << "\nto (" << tag.to << ") >= len of sequence (" << getLenSeq() << ")?");

  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

multitag_t::mte_id_t Read::getTagID(const std::string & identifier)
{
  return multitag_t::hasIdentifier(identifier);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Read::hasTag(const std::string & identifier, int32 pos) const
{
  multitag_t::mte_id_t id=multitag_t::hasIdentifier(identifier);
  if(!multitag_t::getIdentifierStr(id).empty()){
    return hasTag(id,pos);
  }
  return false;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Read::hasTag(const multitag_t::mte_id_t identifier, int32 pos) const
{
  FUNCSTART("Read::hasTag(const multitag_t::mte_id_t identifier, int32 pos) const");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  BUGIFTHROW(pos>=0 && pos>=static_cast<int32>(REA_qualities.size()),
	     getName() << ": pos(" << pos << ") > read size " << REA_qualities.size() << ")?");

  auto tI=REA_tags.cbegin();
  // handle case where we don't care about tag position
  if(pos<0) {
    for(;tI!=REA_tags.cend();++tI) {
      if(tI->identifier==identifier){
	// we found this id in the read, we don't need to examine further
	FUNCEND();
	return true;
      }
    }
  }else{
    // here we need to care about tag position
    for(;tI!=REA_tags.cend();++tI) {
      if(static_cast<int32>(tI->to) < pos
	 || static_cast<int32>(tI->from) > pos) continue;
      // tag at  specified position, check identifier
      if(tI->identifier==identifier){
	// we found this id in the read at that pos,
	//  we don't need to examine further
	FUNCEND();
	return true;
      }
    }
  }

  FUNCEND();
  return false;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint32 Read::countTags(const multitag_t::mte_id_t identifier, int32 pos) const
{
  FUNCSTART("uint32 Read::countTags(const multitag_t::mte_id_t identifier, int32 pos) const");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  BUGIFTHROW(pos>=0 && pos>=static_cast<int32>(REA_qualities.size()),
	     getName() << ": pos(" << pos << ") > read size " << REA_qualities.size() << ")?");

  // handle case where we don't care about tag position
  uint32 retcount=0;
  auto tI=REA_tags.cbegin();
  if(pos<0) {
    for(;tI!=REA_tags.end();++tI) {
      if(tI->identifier==identifier){
	// we found this id in the read, we don't need to examine further
	++retcount;
      }
    }
  }else{
    // here we need to care about tag position
    for(;tI!=REA_tags.end();++tI) {
      if(static_cast<int32>(tI->to) < pos
	 || static_cast<int32>(tI->from) > pos) continue;
      // tag at  specified position, check identifier
      if(tI->identifier==identifier){
	++retcount;
      }
    }
  }

  FUNCEND();
  return retcount;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const multitag_t & Read::getTag(uint32 tagnumber) const
{
  FUNCSTART("const multitag_t & Read::getTag(uint32 tagnumber) const");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  BUGIFTHROW(tagnumber >= REA_tags.size(),
	     getName() << ": tagnumber (" << tagnumber << ") >= REA_tags.size (" << REA_tags.size() << ") ?");

  FUNCEND();

  return REA_tags[tagnumber];
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::helper_refreshPaddedSequence() const
{
  FUNCSTART("void Read::helper_refreshPaddedSequence()");

  if(REA_ps_dirty==true){
    BUGIFTHROW(REA_pcs_dirty==true, "Both seq and compl.seq. are tagged dirty.");
    BUGIFTHROW(checkRead()!=nullptr, checkRead());

    makeComplement(REA_padded_complementsequence, REA_padded_sequence);
    CEBUG("Refreshed.\n");
    REA_ps_dirty=false;

    if(sanityCheck()!=nullptr){
      setCoutType(AS_TEXT);
      cout << *this;
      MIRANOTIFY(Notify::FATAL, sanityCheck());
    }
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Read::helper_refreshPaddedComplementSequence() const
{
  FUNCSTART("void Read::helper_refreshPaddedComplementSequence()");

  if(REA_pcs_dirty){
    BUGIFTHROW(REA_ps_dirty, "Both seq and compl.seq. are tagged dirty.");
    BUGIFTHROW(checkRead()!=nullptr, checkRead());

    CEBUG("Refreshed.\n");
    makeComplement(REA_padded_sequence, REA_padded_complementsequence);
    REA_pcs_dirty=false;
    if(sanityCheck()!=nullptr){
      setCoutType(AS_TEXT);
      cout << *this;
      MIRANOTIFY(Notify::FATAL, sanityCheck());
    }
  }

  FUNCEND();
}




// Tested: Insert/Delete/Change with caching methods + clips

/*************************************************************************
 *
 * Inserts a base at the given position with the given base qualities.
 *    E.g.:   ATGC   inserting N at pos 2 gives    ATNGC
 *
 * If extends_clipped_area is true and the position given is on the border
 * of a clip, the clip will be moved in that way to extend the clipped
 * area.
 *    E.g.:   ATGC with CL=2   (meaning AT GC)
 *            inserting N at pos 2 and extends_clipped_area=true gives
 *            ATNGC with CL=2   (meaning AT NGC)
 *    E.g.:   inserting N at pos 2 and extends_clipped_area=false gives
 *            ATNGC with CL=3   (meaning ATN GC)
 *    E.g.:   ATGC with CR=2   (meaning AT GC)    Note the R for Right!!!
 *            inserting N at pos 2 and extends_clipped_area=true gives
 *            ATNGC with CR=3   (meaning ATN GC)
 *
 * Tags will NEVER! be extended!
 *
 *tested
 *************************************************************************/

void Read::insertBaseInSequence(char base, base_quality_t quality, uint32 position, bool extends_clipped_area)
{
  FUNCSTART("void Read::insertBaseInSequence(char base, base_quality_t quality, int32 position, bool extends_clipped_area)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  CEBUG("Position: " << position << endl);

  if(REA_ps_dirty==true && REA_pcs_dirty==false){
    insertBaseInComplementSequence(dptools::getComplementIUPACBase(base),
				   quality,
				   static_cast<uint32>(REA_padded_complementsequence.size())-position,
				   extends_clipped_area);
  }else{
    refreshPaddedSequence();
    BUGIFTHROW(position>REA_padded_sequence.size(),getName() << ": position (" << position << ") > REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

    moderateContainerGrowth();

    {
      auto I=REA_padded_sequence.begin();
      advance(I, position);
      REA_padded_sequence.insert(I, base);
    }

    {
      BUGIFTHROW(position>REA_qualities.size(),getName() << ": position (" << position << ") > REA_qualities.size (" << REA_qualities.size() << ") ?");
      auto J=REA_qualities.begin();
      advance(J,position);
      REA_qualities.insert(J, quality);
    }

    // insert -1 into adjustment vector, as this base is not in the
    //  original sequence
    if(REA_uses_adjustments){
      BUGIFTHROW(position > REA_adjustments.size(),getName() << ": position (" << position << " > REA_adjustments.size (" << REA_adjustments.size() << ") ?");
      auto K=REA_adjustments.begin();
      advance(K, position);
      REA_adjustments.insert(K, -1);
    }

    // insert 1 into RLE value vector, as this base is not in the
    //  original sequence
    // strictly speaking, this may now not be an RLE vector now ... but it's how MIRA needs it
    if(REA_rlevptr!=nullptr){
      BUGIFTHROW(position > REA_rlevptr->size(),getName() << ": position (" << position << " > REA_rlevptr->size (" << REA_rlevptr->size() << ") ?");
      auto K=REA_rlevptr->begin();
      advance(K, position);
      REA_rlevptr->insert(K, 1);
    }

    // insert baseflags into baseflags_t vector (with defaults)
    {
      BUGIFTHROW(position > REA_bposhashstats.size(),getName() << ": position (" << position << " > REA_bposhashstats.size (" << REA_bposhashstats.size() << ") ?");
      auto B=REA_bposhashstats.begin();
      advance(B, position);
      REA_bposhashstats.insert(B, REA_bposhashstat_default);
    }


    REA_pcs_dirty=true;

    if(extends_clipped_area==true){
      if(static_cast<int32>(position) <  REA_ql) REA_ql++;
      if(static_cast<int32>(position) <  REA_sl) REA_sl++;
      if(static_cast<int32>(position) <  REA_ml) REA_ml++;
      if(static_cast<int32>(position) <= REA_qr) REA_qr++;
      if(static_cast<int32>(position) <= REA_sr) REA_sr++;
      if(static_cast<int32>(position) <= REA_mr) REA_mr++;
      if(REA_cl>=0 &&
	 static_cast<int32>(position) <  REA_cl) REA_cl++;
      if(REA_cr>=0 &&
	 static_cast<int32>(position) <= REA_cr) REA_cr++;
    }else{
      if(static_cast<int32>(position) <= REA_ql) REA_ql++;
      if(static_cast<int32>(position) <= REA_sl) REA_sl++;
      if(static_cast<int32>(position) <= REA_ml) REA_ml++;
      if(static_cast<int32>(position) <  REA_qr) REA_qr++;
      if(static_cast<int32>(position) <  REA_sr) REA_sr++;
      if(static_cast<int32>(position) <  REA_mr) REA_mr++;
      if(REA_cl>=0 &&
	 static_cast<int32>(position) <= REA_cl) REA_cl++;
      if(REA_cr>=0 &&
	 static_cast<int32>(position) <  REA_cr) REA_cr++;
    }

    updateClipoffs();
    updateTagBaseInserted(position);

  }

  BUGIFTHROW(sanityCheck()!=nullptr, sanityCheck());

  FUNCEND();
}


// TODO: check and test tag handling
/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::updateTagBaseInserted(uint32 position)
{
  FUNCSTART("void Read::updateTagInsert(uint32 position)");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  auto tI=REA_tags.begin();
  for(;tI!=REA_tags.end();++tI){
    if(tI->from >= position) tI->from+=1;
    if(tI->to >= position) tI->to+=1;
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::updateTagBaseDeleted(uint32 position)
{
  FUNCSTART("void Read::updateTagBaseDeleted(int32 position)");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  auto tI=REA_tags.begin();
  while(tI!=REA_tags.end()){
    if(tI->from == tI->to &&
       tI->to == position){
      tI=REA_tags.erase(tI);
    }else{
      if(position < tI->from) tI->from--;
      if(position <= tI->to) tI->to--;
      ++tI;
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *tested
 *************************************************************************/
void Read::deleteBaseFromSequence(uint32 position)
{
  FUNCSTART("void Read::deleteBaseFromSequence(int32 position)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  CEBUG("Position: " << position << endl);

  if(REA_ps_dirty==true && REA_pcs_dirty==false){
    deleteBaseFromComplementSequence(static_cast<uint32>(REA_padded_complementsequence.size())-position-1);
  }else{
    refreshPaddedSequence();
    BUGIFTHROW(position>=REA_padded_sequence.size(),getName() << ": position (" << position << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

    {
      auto I=REA_padded_sequence.begin();
      advance(I, position);
      REA_padded_sequence.erase(I);
    }

    {
      BUGIFTHROW(position>=REA_qualities.size(), getName() << ": position (" << position << ") >= REA_qualities.size (" << REA_qualities.size() << ") ?");
      auto J=REA_qualities.begin();
      advance(J,position);
      REA_qualities.erase(J);
    }

    if(REA_uses_adjustments){
      BUGIFTHROW(position >= REA_adjustments.size(),getName() << ": position (" << position << " >= REA_adjustments.size (" << REA_adjustments.size() << ") ?");
      auto K=REA_adjustments.begin();
      advance(K, position);
      REA_adjustments.erase(K);
    }

    // for RLE reads, this does not only delete a base, but eventually a whole
    //  stretch
    // That's OK for now, MIRA needs it that way.
    if(REA_rlevptr!=nullptr){
      BUGIFTHROW(position >= REA_rlevptr->size(),getName() << ": position (" << position << " >= REA_rlevptr->size (" << REA_rlevptr->size() << ") ?");
      auto K=REA_rlevptr->begin();
      advance(K, position);
      REA_rlevptr->erase(K);
    }

    {
      BUGIFTHROW(position >= REA_bposhashstats.size(),getName() << ": position (" << position << " >= REA_bposhashstats.size (" << REA_bposhashstats.size() << ") ?");
      auto B=REA_bposhashstats.begin();
      advance(B, position);
      REA_bposhashstats.erase(B);
    }

    REA_pcs_dirty=true;

    // FALSCH! was <, but I think that was wrong. <= should be right
    // TODO: was <, but I think that was wrong. <= should be right
    if(static_cast<int32>(position)<REA_ql && REA_ql > 0) REA_ql--;
    if(static_cast<int32>(position)<REA_qr && REA_qr > 0) REA_qr--;
    if(static_cast<int32>(position)<REA_ml && REA_ml > 0) REA_ml--;
    if(static_cast<int32>(position)<REA_mr && REA_mr > 0) REA_mr--;
    if(static_cast<int32>(position)<REA_sl && REA_sl > 0) REA_sl--;
    if(static_cast<int32>(position)<REA_sr && REA_sr > 0) REA_sr--;
    if(REA_cl>=0 &&
       static_cast<int32>(position)<REA_cl && REA_cl > 0) REA_cl--;
    if(REA_cr>=0 &&
       static_cast<int32>(position)<REA_cr && REA_cr > 0) REA_cr--;

    updateClipoffs();
    updateTagBaseDeleted(position);

  }

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *tested
 *************************************************************************/
void Read::insertBaseInComplementSequence(char base, base_quality_t quality, uint32 position, bool extends_clipped_area)
{
  FUNCSTART("void Read::insertBaseInComplementSequence(char base, base_quality_t quality, int32 position, bool extends_clipped_area)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  CEBUG("Position: " << position << endl);

  if(REA_pcs_dirty==true && REA_ps_dirty==false){
    insertBaseInSequence(dptools::getComplementIUPACBase(base),
			 quality,
			 static_cast<uint32>(REA_padded_sequence.size())-position,
			 extends_clipped_area);
  }else{
    refreshPaddedComplementSequence();

    BUGIFTHROW(position>=REA_padded_complementsequence.size(),getName() << ": position (" << position << ") >= REA_padded_complementsequence.size (" << REA_padded_sequence.size() << ") ?");

    moderateContainerGrowth();

    uint32 complement_position=static_cast<uint32>(REA_padded_complementsequence.size())-position;

    CEBUG("compl_pos: " << complement_position << endl);

    {
      BUGIFTHROW(position>REA_padded_complementsequence.size(),getName() << ": position (" << position << ") > REA_padded_complementsequence.size (" << REA_padded_complementsequence.size() << ") ?");
      auto I=REA_padded_complementsequence.begin();
      advance(I, position);
      REA_padded_complementsequence.insert(I, base);
    }

    {
      BUGIFTHROW(complement_position>REA_qualities.size(),getName() << ": complement_position (" << complement_position << " > REA_qualities.size (" << REA_qualities.size() << ") ?");
      auto J=REA_qualities.begin();
      advance(J,complement_position);
      REA_qualities.insert(J, quality);
    }

    // insert -1 into adjustment vector, as this base is not in the
    //  original sequence
    if(REA_uses_adjustments){
      BUGIFTHROW(complement_position > REA_adjustments.size(),getName() << ": complement_position (" << complement_position << ") > REA_adjustments.size (" << REA_adjustments.size() << ")");
      auto K=REA_adjustments.begin();
      advance(K, complement_position);
      REA_adjustments.insert(K, -1);
    }

    // insert 1 into RLE value vector, as this base is not in the
    //  original sequence
    // strictly speaking, this may now not be an RLE vector now ... but it's how MIRA needs it
    if(REA_rlevptr!=nullptr){
      BUGIFTHROW(complement_position > REA_rlevptr->size(),getName() << ": complement_position (" << complement_position << " > REA_rlevptr->size (" << REA_rlevptr->size() << ") ?");
      auto K=REA_rlevptr->begin();
      advance(K, complement_position);
      REA_rlevptr->insert(K, 1);
    }

    // insert baseflags into baseflags_t vector (with defaults)
    {
      BUGIFTHROW(complement_position > REA_bposhashstats.size(),getName() << ": complement_position (" << position << " > REA_bposhashstats.size (" << REA_bposhashstats.size() << ") ?");
      auto B=REA_bposhashstats.begin();
      advance(B, complement_position);
      REA_bposhashstats.insert(B, REA_bposhashstat_default);
    }

    REA_ps_dirty=true;

    if(extends_clipped_area==true){
      if(static_cast<int32>(complement_position) <  REA_ql) REA_ql++;
      if(static_cast<int32>(complement_position) <  REA_sl) REA_sl++;
      if(static_cast<int32>(complement_position) <  REA_ml) REA_ml++;
      if(static_cast<int32>(complement_position) <= REA_qr) REA_qr++;
      if(static_cast<int32>(complement_position) <= REA_sr) REA_sr++;
      if(static_cast<int32>(complement_position) <= REA_mr) REA_mr++;
      if(REA_cl>=0 &&
	 static_cast<int32>(complement_position) <  REA_cl) REA_cl++;
      if(REA_cr>=0 &&
	 static_cast<int32>(complement_position) <= REA_cr) REA_cr++;
    }else{
      if(static_cast<int32>(complement_position) <= REA_ql) REA_ql++;
      if(static_cast<int32>(complement_position) <= REA_sl) REA_sl++;
      if(static_cast<int32>(complement_position) <= REA_ml) REA_ml++;
      if(static_cast<int32>(complement_position) <  REA_qr) REA_qr++;
      if(static_cast<int32>(complement_position) <  REA_sr) REA_sr++;
      if(static_cast<int32>(complement_position) <  REA_mr) REA_mr++;
      if(REA_cl>=0 &&
	 static_cast<int32>(complement_position) <= REA_cl) REA_cl++;
      if(REA_cr>=0 &&
	 static_cast<int32>(complement_position) <  REA_cr) REA_cr++;
    }

    updateClipoffs();
    updateTagBaseInserted(complement_position);

  }

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *tested
 *************************************************************************/
void Read::deleteBaseFromComplementSequence(uint32 position)
{
  FUNCSTART("void Read::deleteBaseFromComplementSequence(int32 uposition)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  CEBUG("Position: " << position << endl);

  if(REA_pcs_dirty==true && REA_ps_dirty==false){
    deleteBaseFromSequence(static_cast<uint32>(REA_padded_sequence.size())-position-1);
  }else{
    refreshPaddedComplementSequence();

    BUGIFTHROW(position>=REA_padded_complementsequence.size(),getName() << ": position (" << position << ") >= REA_padded_complementsequence.size (" << REA_padded_sequence.size() << ") ?");

    uint32 complement_position=static_cast<uint32>(REA_padded_complementsequence.size())-position-1;

    {
      BUGIFTHROW(position>REA_padded_complementsequence.size(),getName() << ": position (" << position << ") > REA_padded_complementsequence.size (" << REA_padded_complementsequence.size() << ") ?");
      auto I=REA_padded_complementsequence.begin();
      advance(I, position);
      REA_padded_complementsequence.erase(I);
    }

    {
      BUGIFTHROW(complement_position>=REA_qualities.size(),getName() << ": complement_position (" << complement_position << ") >= REA_qualities.size (" << REA_qualities.size() << ") ?");
      auto J=REA_qualities.begin();
      advance(J,complement_position);
      REA_qualities.erase(J);
    }

    if(REA_uses_adjustments){
      BUGIFTHROW(complement_position>=REA_adjustments.size(),getName() << " complement_position (" << complement_position << ") >= REA_adjustments.size (" << REA_adjustments.size() << ") ?");
      auto K=REA_adjustments.begin();
      advance(K, complement_position);
      REA_adjustments.erase(K);
    }

    // for RLE reads, this does not only delete a base, but eventually a whole
    //  stretch
    // That's OK for now, MIRA needs it that way.
    if(REA_rlevptr!=nullptr){
      BUGIFTHROW(complement_position >= REA_rlevptr->size(),getName() << ": position (" << complement_position << " >= REA_rlevptr->size (" << REA_rlevptr->size() << ") ?");
      auto K=REA_rlevptr->begin();
      advance(K, complement_position);
      REA_rlevptr->erase(K);
    }

    {
      BUGIFTHROW(complement_position>=REA_bposhashstats.size(),getName() << " complement_position (" << complement_position << ") >= REA_bposhashstats.size (" << REA_bposhashstats.size() << ") ?");
      auto B=REA_bposhashstats.begin();
      advance(B, complement_position);
      REA_bposhashstats.erase(B);
    }

    REA_ps_dirty=true;

    // TODO: was <, but I think that was wrong. <= should be right
    // update: <= is wrong, due to new complpos calc, < should be right
    if(static_cast<int32>(complement_position)<REA_ql && REA_ql > 0) REA_ql--;
    if(static_cast<int32>(complement_position)<REA_sl && REA_sl > 0) REA_sl--;
    if(static_cast<int32>(complement_position)<REA_ml && REA_ml > 0) REA_ml--;
    if(static_cast<int32>(complement_position)<REA_qr && REA_qr > 0) REA_qr--;
    if(static_cast<int32>(complement_position)<REA_sr && REA_sr > 0) REA_sr--;
    if(static_cast<int32>(complement_position)<REA_mr && REA_mr > 0) REA_mr--;
    if(REA_cl>=0 &&
       static_cast<int32>(complement_position)<REA_cl && REA_cl> 0) REA_cl--;
    if(REA_cr>=0 &&
       static_cast<int32>(complement_position)<REA_cr && REA_cr > 0) REA_cr--;

    updateClipoffs();
    updateTagBaseDeleted(complement_position);

  }

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}






/*************************************************************************
 *
 * qual==255 -> no qual change
 *
 *
 *tested
 *************************************************************************/

void Read::changeBaseInSequence(char base, base_quality_t quality, uint32 position)
{
  FUNCSTART("void changeBaseInSequence(char base, base_quality_t quality, int32 position)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  CEBUG("Position: " << position << endl);

  if(REA_ps_dirty==false){
    CEBUG("ps clean\n");
    BUGIFTHROW(position>=REA_padded_sequence.size(), getName() << ": position (" << position << ") >= size of read (" << REA_padded_sequence.size() << ")?");

    REA_padded_sequence[position]=base;
  }

  if(REA_pcs_dirty==false){
    CEBUG("pcs clean\n");
    BUGIFTHROW(position>=REA_padded_complementsequence.size(),getName() << ": position (" << position << ") >= REA_padded_complementsequence.size (" << REA_padded_complementsequence.size() << ") ?");

    uint32 cposition=static_cast<uint32>(REA_padded_complementsequence.size())-position-1;

    BUGIFTHROW(cposition>=REA_padded_complementsequence.size(),getName() << ": cposition (" << " ) >= REA_padded_complementsequence.size (" << REA_padded_complementsequence.size() << ") ?");

    REA_padded_complementsequence[cposition]=dptools::getComplementIUPACBase(base);
  }

  if(quality!=255){
    BUGIFTHROW(position>=REA_qualities.size(),getName() << ": position (" << position << ") >= REA_qualities.size (" << REA_qualities.size() << ") ?");
    REA_qualities[position]=quality;
  }
  REA_bposhashstats[position]=REA_bposhashstat_default;

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

void Read::changeBaseInComplementSequence(char base, base_quality_t quality, uint32 position)
{
  FUNCSTART("void changeBaseInComplementSequence(char base, base_quality_t quality, int32 position)");

  CEBUG("Position: " << position << endl);

  changeBaseInSequence(dptools::getComplementIUPACBase(base), quality, getLenSeq()-position-1);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *untested
 *************************************************************************/

void Read::changeAdjustment(uint32 position, int32 newadjustment)
{
  FUNCSTART("void Read::changeAdjustment(uint32 position, int32 newadjustment)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  BUGIFTHROW(REA_adjustments.empty(), getName() << ": trying to change adjustment in read which does not have any?");
  BUGIFTHROW(position>=REA_adjustments.size(), getName() << ": trying to change adjustment at position " << position << " but size of adjustment is only " << REA_adjustments.size());

  REA_adjustments[position]=newadjustment;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

void Read::insertBaseInClippedSequence(char base, base_quality_t quality, uint32 position, bool extends_clipped_area)
{
  FUNCSTART("void Read::insertBaseInClippedSequence(char base, base_quality_t quality, uint32 position, bool extends_clipped_area)");

  CEBUG("Position: " << position << endl);

  BUGIFTHROW(static_cast<int32>(position)>getLenClippedSeq(), "Position > len of clipped read ?");

  insertBaseInSequence(base, quality, position+getLeftClipoff(),
		       extends_clipped_area);

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

void Read::deleteBaseFromClippedSequence(uint32 position)
{
  FUNCSTART("void Read::deleteBaseFromClippedSequence(uint32 position)");

  CEBUG("Position: " << position << endl);

  BUGIFTHROW(static_cast<int32>(position)>=getLenClippedSeq(), "Position >= len of clipped read ?");

  deleteBaseFromSequence(position+getLeftClipoff());

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

void Read::insertBaseInClippedComplementSequence(char base, base_quality_t quality, uint32 position, bool extends_clipped_area)
{
  FUNCSTART("void Read::insertBaseInClippedComplementSequence(char base, base_quality_t quality, int32 position, bool extends_clipped_area)");

  CEBUG("Position: " << position << endl);

  BUGIFTHROW(static_cast<int32>(position)>getLenClippedSeq(), "Position > len of clipped read ?");

  insertBaseInComplementSequence(base, quality,
				 getLenSeq()-getRightClipoff()+position,
				 extends_clipped_area);

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

void Read::deleteBaseFromClippedComplementSequence(uint32 position)
{
  FUNCSTART("void Read::deleteBaseFromClippedComplementSequence(int32 position)");

  CEBUG("Position: " << position << endl);

  BUGIFTHROW(static_cast<int32>(position)>=getLenClippedSeq(), "Position >= len of clipped read ?");

  deleteBaseFromComplementSequence(getLenSeq()-getRightClipoff()+position);

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::changeBaseInClippedSequence(char base, base_quality_t quality, uint32 position)
{
  FUNCSTART("void Read::changeBaseInClippedSequence(char base, base_quality_t quality, uint32 position)");

  CEBUG("Position: " << position << endl);

  BUGIFTHROW(static_cast<int32>(position)>=getLenClippedSeq(), "Position >= len of clipped read ?");

  changeBaseInSequence(base, quality, getLeftClipoff()+position);

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::changeBaseInClippedComplementSequence(char base, base_quality_t quality, uint32 position)
{
  FUNCSTART("void Read::changeBaseInClippedComplementSequence(char base, base_quality_t quality, uint32 position)");

  CEBUG("Position: " << position << endl);

  BUGIFTHROW(static_cast<int32>(position)>=getLenClippedSeq(), "Position >= len of clipped read ?");

  changeBaseInComplementSequence(base, quality, getLenSeq()-getRightClipoff()+position);

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

// This should be
//    const std::vector<char> & Read::getActualSequence()  const
// too, but SGI CC 7.1 doesn't know the keyword mutable
const std::vector<char> & Read::getActualSequence() const
{
  FUNCSTART("const std::vector<char> & Read::getActualSequence() const");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  refreshPaddedSequence();

  FUNCEND();
  return REA_padded_sequence;
}


// This should be
//    const std::vector<char> & Read::getActualComplementSequence  const
// too, but SGI CC 7.1 doesn't know the keyword mutable
const std::vector<char> & Read::getActualComplementSequence() const
{
  FUNCSTART("const std::vector<char> & Read::getActualComplementSequence() const");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));
  refreshPaddedComplementSequence();

  FUNCEND();
  return REA_padded_complementsequence;
}



/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/
// This should be const
// too, but SGI CC 7.1 doesn't know the keyword mutable
const char * Read::getClippedSeqAsChar() const
{
  FUNCSTART("const char * Read::getClippedSeqAsChar()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedSequence();

  if(REA_padded_sequence.empty()||getLeftClipoff()==REA_padded_sequence.size()) return &REA_zerostring;

  BOUNDCHECK(getLeftClipoff(), 0, static_cast<int32>(REA_padded_sequence.size()));

  auto cI=REA_padded_sequence.begin();
  advance(cI, getLeftClipoff());
  return &(*cI);
}



/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/
const char * Read::getSeqAsChar() const
{
  FUNCSTART("const char * Read::getSeqAsChar()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedSequence();

  if(REA_padded_sequence.empty()) return &REA_zerostring;

  auto cI=REA_padded_sequence.cbegin();
  return &(*cI);
}


/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/
const char * Read::getComplementSeqAsChar() const
{
  FUNCSTART("const char * Read::getComplementSeqAsChar()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedComplementSequence();

  if(REA_padded_complementsequence.empty()) return &REA_zerostring;

  auto cI=REA_padded_complementsequence.cbegin();
  return &(*cI);
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Read::getSeqAsString(std::string & result)
{
  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedSequence();

  result.resize(getLenSeq());
  auto cI=REA_padded_sequence.cbegin();
  auto sI=result.begin();
  for(; cI != REA_padded_sequence.cend(); ++cI, ++sI){
    *sI=*cI;
  }
  return;
}




/*************************************************************************
 *
 *
 *
 *tested
 *************************************************************************/

// This should be const
// too, but SGI CC 7.1 doesn't know the keyword mutable
const char * Read::getClippedComplementSeqAsChar() const
{
  FUNCSTART("const char * Read::getClippedComplementSeqAsChar()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedComplementSequence();
  if(REA_padded_complementsequence.empty()) return &REA_zerostring;

  BOUNDCHECK(static_cast<int32>(REA_padded_complementsequence.size())-getRightClipoff(), 0, static_cast<int32>(REA_padded_complementsequence.size())+1);

  auto cI=REA_padded_complementsequence.cbegin();
  advance(cI, REA_padded_complementsequence.size()-getRightClipoff());

  return &(*cI);
}





/*************************************************************************
 *
 * Returns an iterator to the sequencing vector and quality clipped read
 *
 *tested
 *************************************************************************/

std::vector<char>::const_iterator Read::getClippedSeqIterator() const
{
  FUNCSTART("const std::vector<char>::iterator Read::getClippedSeqIterator()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedSequence();

  auto I=REA_padded_sequence.cbegin();

  if(getLeftClipoff() < 0 ||
     getLeftClipoff() >= static_cast<int32>(REA_padded_sequence.size())){
    setCoutType(AS_TEXT);
    cout << '\n' << *this << '\n';
    cout.flush();
  }

  BOUNDCHECK(getLeftClipoff(), 0, static_cast<int32>(REA_padded_sequence.size()));
  advance(I, getLeftClipoff());

  FUNCEND();

  return I;
}



/*************************************************************************
 *
 * Returns an iterator to the read
 *
 *tested
 *************************************************************************/

std::vector<char>::const_iterator Read::getSeqIteratorBegin() const
{
  FUNCSTART("const std::vector<char>::iterator Read::getSeqIteratorBegin()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedSequence();

  FUNCEND();

  return REA_padded_sequence.cbegin();
}

std::vector<char>::const_iterator Read::getSeqIteratorEnd() const
{
  FUNCSTART("const std::vector<char>::iterator Read::getSeqIteratorEnd()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedSequence();

  FUNCEND();

  return REA_padded_sequence.cend();
}



/*************************************************************************
 *
 * Returns an iterator to the
 *  complement read
 *
 *tested
 *************************************************************************/

std::vector<char>::const_iterator Read::getComplementSeqIteratorBegin() const
{
  FUNCSTART("const std::vector<char>::iterator Read::getComplementSeqIteratorBegin()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedComplementSequence();

  FUNCEND();

  return REA_padded_complementsequence.cbegin();
}

std::vector<char>::const_iterator Read::getComplementSeqIteratorEnd() const
{
  FUNCSTART("const std::vector<char>::iterator Read::getComplementSeqIteratorEnd()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedComplementSequence();

  FUNCEND();

  return REA_padded_complementsequence.cend();
}



/*************************************************************************
 *
 * Returns an iterator to the sequencing vector and quality clipped
 *  complement read
 *
 *tested
 *************************************************************************/

std::vector<char>::const_iterator Read::getClippedComplementSeqIterator() const
{
  FUNCSTART("const std::vector<char>::iterator Read::getClippedComplementSeqIterator()");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  refreshPaddedComplementSequence();
  BOUNDCHECK(static_cast<int32>(REA_padded_complementsequence.size())-getRightClipoff(), 0, static_cast<int32>(REA_padded_complementsequence.size()+1));

  auto I=REA_padded_complementsequence.cbegin();
  advance(I, REA_padded_complementsequence.size()-getRightClipoff());

  FUNCEND();

  return I;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

char Read::getBaseInSequence(uint32 pos) const
{
  FUNCSTART("char Read::getBaseInSequence(uint32 pos)");

  refreshPaddedSequence();

  BUGIFTHROW(pos >= REA_padded_sequence.size(), getName() << ": pos (" << pos << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  FUNCEND();

  return REA_padded_sequence[pos];
}

char Read::getBaseInClippedSequence(uint32 pos) const
{
  FUNCSTART("char Read::getBaseInClippedSequence(uint32 pos)");

  refreshPaddedSequence();

  BUGIFTHROW(pos >= REA_padded_sequence.size(), getName() << ": pos (" << pos << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  FUNCEND();

  return REA_padded_sequence[pos+getLeftClipoff()];
}

char Read::getBaseInComplementSequence(uint32 pos) const
{
  FUNCSTART("char Read::getBaseInComplementSequence(int32 pos)");

  refreshPaddedComplementSequence();

  FUNCEND();

  return REA_padded_complementsequence[pos];
}

base_quality_t Read::getQualityInSequence(uint32 pos) const
{
  FUNCSTART("char Read::getQualityInSequence(int32 pos)");

  //refreshPaddedSequence();

  BUGIFTHROW(pos >= REA_qualities.size(), getName() << ": pos (" << pos << ") >= REA_qualities.size (" << REA_qualities.size()  << ") ?");

  FUNCEND();

  return REA_qualities[pos];
}

base_quality_t Read::getQualityInComplementSequence(uint32 pos) const
{
  FUNCSTART("char Read::getQualityInComplementSequence(int32 pos)");

  //refreshPaddedComplementSequence();

  BUGIFTHROW(pos >= static_cast<uint32>(REA_qualities.size()), getName() << ": pos (" << pos << ") >= REA_qualities.size (" << REA_qualities.size()  << ") ?");

  FUNCEND();

  return REA_qualities[REA_qualities.size()-1-pos];
}


/*************************************************************************
 *
 * pos in unclipped sequence
 *
 *
 *************************************************************************/

Read::bposhashstat_t Read::getBPosHashStats(uint32 pos) const
{
  FUNCSTART("baseflags_t Read::getBaseFlags(uint32 pos) const");

  BUGIFTHROW(pos >= REA_bposhashstats.size(), getName() << ": pos (" << pos << ") >= REA_bposhashstats.size (" << REA_bposhashstats.size() << ") ?");

  FUNCEND();

  return REA_bposhashstats[pos];
}


/*************************************************************************
 *
 * returns the position in unclipped sequence of first base of a
 *  consecutive run of bases (including)
 *
 *
 *
 *************************************************************************/

uint32 Read::getLowerBoundPosOfBaseRun(uint32 pos, char base, const bool alsotakegap) const
{
  FUNCSTART("uint32 Read::getLowerBoundPosOfBaseRun(uint32 pos, char base, const bool alsotakegap) const");

  refreshPaddedSequence();

  BUGIFTHROW(pos >= REA_padded_sequence.size(), getName() << ": pos (" << pos << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  if(pos==0) return 0;
  if(!alsotakegap && REA_padded_sequence[pos]=='*') return pos;

  base=static_cast<char>(toupper(base));
  for(; pos>0; pos--){
    if(!alsotakegap && REA_padded_sequence[pos-1]=='*') return pos;
    if(REA_padded_sequence[pos-1]!='*'
       && toupper(REA_padded_sequence[pos-1]) !=base) return pos;
  }

  FUNCEND();
  return 0;
}

/*************************************************************************
 *
 * returns the position in unclipped sequence of last base of a
 *  consecutive run of bases (including)
 *
 *
 *************************************************************************/

uint32 Read::getUpperBoundPosOfBaseRun(uint32 pos, char base, const bool alsotakegap) const
{
  FUNCSTART("uint32 Read::getUpperBoundPosOfBaseRun(uint32 pos, char base, const bool alsotakegap) const");

  refreshPaddedSequence();

  BUGIFTHROW(pos >= REA_padded_sequence.size(), getName() << ": pos (" << pos << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  if(pos==REA_padded_sequence.size()-1) return pos;
  if(!alsotakegap && REA_padded_sequence[pos]=='*') return pos;

  base=static_cast<char>(toupper(base));
  for(; pos<REA_padded_sequence.size()-1; pos++){
    if(!alsotakegap && REA_padded_sequence[pos+1]=='*') return pos;
    if(REA_padded_sequence[pos+1]!='*'
       && toupper(REA_padded_sequence[pos+1])!=base) return pos;
  }

  FUNCEND();
  return static_cast<uint32>(REA_padded_sequence.size())-1;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

uint32 Read::getLenOfGapRun(uint32 pos) const
{
  FUNCSTART("uint32 Read::getLenOfGapRun(uint32 pos) const");

  refreshPaddedSequence();

  BUGIFTHROW(pos >= REA_padded_sequence.size(), getName() << ": pos (" << pos << ") >= REA_padded_sequence.size (" << REA_padded_sequence.size() << ") ?");

  if(REA_padded_sequence[pos]!='*') {
    FUNCEND();
    return 0;
  }

  while(pos>0 && REA_padded_sequence[pos-1]=='*') pos--;
  uint32 spos=pos;
  while(pos<REA_padded_sequence.size()-1 && REA_padded_sequence[pos+1]=='*') pos++;

  FUNCEND();
  return pos-spos+1;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Read::calcTemplateInfo()
{
  FUNCSTART("std::string Read::calcTemplateInfo()");

  BUGIFTHROW(!REA_template.empty(),getName() << ": template already set to " << REA_template << " ?");
  BUGIFTHROW(REA_template_segment!=0,getName() << ": template segment already set to " << static_cast<uint32>(REA_template_segment) << " ?");
  //BUGIFTHROW(REA_rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA && REA_rgid.getReadNamingScheme() == ReadGroupLib::SCHEME_NONE,"Oooops, Solexa type but no naming scheme?");

  switch(REA_rgid.getReadNamingScheme()) {
  case ReadGroupLib::SCHEME_SANGER: {
    getInternalTemplateName_Sanger(REA_template,REA_template_segment);
    break;
  }
  case ReadGroupLib::SCHEME_TIGR: {
    getInternalTemplateName_TIGR(REA_template,REA_template_segment);
    break;
  }
  case ReadGroupLib::SCHEME_SOLEXA: {
    getInternalTemplateName_Solexa(REA_template,REA_template_segment);
    break;
  }
  case ReadGroupLib::SCHEME_FR: {
    getInternalTemplateName_FR(REA_template,REA_template_segment);
    break;
  }
  case ReadGroupLib::SCHEME_STLOUIS: {
    getInternalTemplateName_StLouis(REA_template,REA_template_segment);
    break;
  }
  case ReadGroupLib::SCHEME_SRARAW: {
    getInternalTemplateName_SRARaw(REA_template,REA_template_segment);
    break;
  }
  case ReadGroupLib::SCHEME_NONE: {
    REA_template_segment=0;
    REA_template=getName();
    break;
  }
  case ReadGroupLib::SCHEME_EMPTY: {
    REA_template_segment=0;
    REA_template.clear();
    break;
  }
  default: {
    MIRANOTIFY(Notify::INTERNAL, getName() << ": unknown read naming scheme " << static_cast<uint32>(REA_rgid.getReadNamingScheme()) << " set??? ReadGroupLib::SCHEME_NONE " << ReadGroupLib::SCHEME_NONE);
  }
  }

  //cout << "### " << getName() << "\t" << REA_template << "\t" << static_cast<uint32>(REA_template_segment) << endl;
  //cout << *this;

  FUNCEND();
  return;
}



/*************************************************************************
 *
 * Deduce template name from sanger centre naming scheme
 *  also set template_segment
 *
 *************************************************************************/

void Read::getInternalTemplateName_Sanger(std::string & tname, uint8 & segment)
{
  std::string templbase;
  std::string sangerext;

  segment=0;
  tname.clear();

  auto bpos = getName().rfind(".");

  if (bpos != std::string::npos) {
    templbase=getName().substr(0,bpos);
    sangerext=getName().substr(bpos,getName().size());
    if(!sangerext.empty()){
      if(toupper(sangerext[1])=='P'
	 || toupper(sangerext[1])=='F'){
	segment=1;
      }else if(toupper(sangerext[1])=='Q'
	       || toupper(sangerext[1])=='R'){
	segment=255;
      }
    }
    tname=templbase;
    if(sangerext.size()>2) tname+='_'+sangerext.substr(2,sangerext.size());
  }

  if(tname.empty()) {
    tname=getName();
  }

  //cout << "getName(): " << getName()<< endl;
  //cout << "templbase: " << templbase << endl;
  //cout << "sangerext: " << sangerext << endl;
  //cout << "sangerext.substr(2,sangerext.size()): " << sangerext.substr(2,sangerext.size()) << endl;
  //cout << "return; " << (templbase+sangerext.substr(2,sangerext.size())) << endl;

  return;
}



/*************************************************************************
 *
 *  Deduce from TIGR centre naming scheme
 *  also set template_segment
 *
 *************************************************************************/

void Read::getInternalTemplateName_TIGR(std::string & tname, uint8 & segment)
{
  std::string templbase;
  std::string otherext;

  //cout << "Scheme TIGR";

  segment=0;
  tname.clear();

  auto bpos = getName().find("TF");

  if (bpos != std::string::npos) {
    templbase=getName().substr(0,bpos);
    segment=1;
  } else {
    bpos = getName().find("TR");
    if (bpos != std::string::npos) {
      templbase=getName().substr(0,bpos);
      segment=255;
    }
  }
  if(bpos!=std::string::npos){
    otherext=getName().substr(bpos,getName().size());
    tname=templbase;
    if(otherext.size()>2) tname+='_'+otherext.substr(2,otherext.size());
  }
  if(tname.empty()){
    tname=getName();
  }

  //cout << "otherext.substr(2,otherext.size()): " << otherext.substr(2,otherext.size()) << endl;
  //cout << "return; " << (templbase+otherext.substr(2,otherext.size())) << endl;

  return;
}



/*************************************************************************
 *
 *  Deduce from FR naming scheme (.f* / .r*)
 *  also set template_segment
 *
 *************************************************************************/

void Read::getInternalTemplateName_FR(std::string & tname, uint8 & segment)
{
  std::string templbase;
  std::string otherext;

  segment=0;
  tname.clear();

  auto bpos = getName().rfind(".");

  if (bpos != std::string::npos) {
    templbase=getName().substr(0,bpos);
    otherext=getName().substr(bpos,getName().size());

    if(otherext.size()>=2){
      if(toupper(otherext[1])=='F'){
	segment=1;
      }else if(toupper(otherext[1])=='R'){
	segment=255;
      }
      tname=templbase;
      if(otherext.size()>2) tname+='_'+otherext.substr(2,otherext.size());
    }
  }
  if(tname.empty()){
    tname=getName();
  }

  return;
}


/*************************************************************************
 *
 *  Deduce from .1/.2 naming scheme (.1* / .2*)
 *  also set template_segment
 *
 *************************************************************************/

void Read::getInternalTemplateName_SRARaw(std::string & tname, uint8 & segment)
{
  std::string templbase;
  std::string otherext;

  segment=0;
  tname.clear();

  auto bpos = getName().rfind(".");

  if (bpos != std::string::npos) {
    templbase=getName().substr(0,bpos);
    otherext=getName().substr(bpos,getName().size());

    if(otherext.size()>=2){
      if(otherext[1]=='1'){
	segment=1;
      }else if(otherext[1]=='2'){
	segment=255;
      }
      tname=templbase;
      if(otherext.size()>2) tname+='_'+otherext.substr(2,otherext.size());
    }
  }
  if(tname.empty()){
    tname=getName();
  }

  return;
}


/*************************************************************************
 *
 *  Deduce from Solexa naming scheme ( /1 or /2)
 *  If unsuccessful, try FR scheme from 454 and Ion reads
 *  also set template_segment
 *
 *************************************************************************/

void Read::getInternalTemplateName_Solexa(std::string & tname, uint8 & segment)
{
  std::string templbase;
  std::string otherext;

  segment=0;
  tname.clear();

  auto bpos = getName().rfind("/");

  if (bpos != std::string::npos) {
    templbase=getName().substr(0,bpos);
    otherext=getName().substr(bpos,getName().size());

    if(otherext.size()>=2){
      if(toupper(otherext[1])=='1'){
	segment=1;
      }else if(toupper(otherext[1])=='2'){
	segment=255;
      }
      tname=templbase;
      if(otherext.size()>2) tname+='_'+otherext.substr(2,otherext.size());
    }
  }else{
    // Hmmm ... no "/" found in the read name.
    // To retain compatibility to earlier MIRA versions and other programs:
    //  if sequencing tech is 454 or Ion, try the FR scheme
    if(REA_rgid.getSequencingType() == ReadGroupLib::SEQTYPE_454GS20
       || REA_rgid.getSequencingType() == ReadGroupLib::SEQTYPE_IONTORRENT){
      getInternalTemplateName_FR(tname,segment);
    }
  }

  if(tname.empty()){
    tname=getName();
  }

  //cout << "getName(): " << getName()<< endl;
  //cout << "templbase: " << templbase << endl;
  //cout << "otherext: " << otherext << endl;
  //cout << "segment: " << static_cast<uint32>(segment) << endl;
  //cout << "tname: " << tname << endl;

  return;
}



/*************************************************************************
 *
 *  Deduce from stlouis naming scheme
 *  also set template_segment
 *
 *************************************************************************/

void Read::getInternalTemplateName_StLouis(std::string & tname, uint8 & segment)
{
  std::string templbase;
  std::string otherext;

  segment=0;
  tname.clear();

  auto bpos = getName().rfind(".");

  if (bpos != std::string::npos) {
    templbase=getName().substr(0,bpos);
    otherext=getName().substr(bpos,getName().size());
    if(!otherext.empty()){
      switch(tolower(otherext[1])){
      case 'r' :
      case 'y' :
      case 'g' : {
	segment=255;
	break;
      }
      case 's' :
      case 'f' :
      case 'x' :
      case 'z' :
      case 'i' :
      case 'b' : {
	segment=1;
      }
      default: {
	// everything else treat as unknown
      }
      }
    }
    tname=templbase;
    if(otherext.size()>2) tname+='_'+otherext.substr(2,otherext.size());
  }
  if(tname.empty()){
    tname=getName();
  }

  //cout << "getName(): " << getName()<< endl;
  //cout << "templbase: " << templbase << endl;
  //cout << "otherext: " << otherext << endl;
  //cout << "otherext.substr(2,otherext.size()): " << otherext.substr(2,otherext.size()) << endl;
  //cout << "return; " << (templbase+otherext.substr(2,otherext.size())) << endl;

  return;
}




/*************************************************************************
 *
 * looks at the tags of the read: if SVEC is present, looks if the read
 *  clips left and right have to be adjusted to cover those tags
 *
 * first: go from left(right) clip and search within 'tolerance' distance
 *  if there's some SVEC
 * if nothing found there, search within first(last) 'clipsearchlen' for
 *  SVEC tags
 *************************************************************************/

void Read::transferSVTagsToClip(int32 tolerancelen, int32 clipsearchlen) {
  FUNCSTART("void Read::transferSVTagsToClip(int32 tolerancelen)");

  BUGIFTHROW(tolerancelen<0, "tolerancelen <0");
  BUGIFTHROW(clipsearchlen<0, "clipsearchlen <0");

  if(REA_has_valid_data==false) return;

  refreshPaddedSequence();

  if(REA_tags.empty()) {
    FUNCEND();
    return;
  }

  //std::string tagSVEC="SVEC";

  bool found=false;
  for(uint32 i=0; i<REA_tags.size(); i++){
    if(REA_tags[i].identifier==REA_tagentry_idSVEC) found=true;
  }
  if(found){
    std::vector<char> tagfield(REA_padded_sequence.size(),0);
    for(uint32 i=0; i<REA_tags.size(); i++){
      if(REA_tags[i].identifier==REA_tagentry_idSVEC) {
	uint32 from=REA_tags[i].from;
	uint32 to=REA_tags[i].to;
	if(to<from) std::swap(from,to);
	if(from>=REA_padded_sequence.size()) continue;
	if(to>=REA_padded_sequence.size()) to=static_cast<uint32>(REA_padded_sequence.size());
	for(uint32 index=from; index<to; index++) {
	  tagfield[index]=1;
	}

	if(getLeftClipoff()<static_cast<int32>(REA_padded_sequence.size())){
	  // first, search within left clip + tolerance
	  auto cI=tagfield.cbegin();
	  cI+=getLeftClipoff();
	  found=false;
	  for(int32 c=0; c<tolerancelen && cI!=tagfield.cend(); ++c, ++cI) {
	    if(*cI==1) {
	      found=true;
	      break;
	    }
	  }
	  // if no SV found, search from begining + some length
	  if(!found) {
	    cI=tagfield.begin();
	    for(int32 c=0; c<clipsearchlen && cI!=tagfield.cend(); ++c, ++cI) {
	      if(*cI==1) {
		found=true;
		break;
	      }
	    }
	  }
	  if(found) {
	    while(found) {
	      if(cI==tagfield.cend()) break;
	      if(*cI==0) break;
	      ++cI;
	    }
	    //REA_sl=cI-tagfield.begin();
	    // +1 experimentell ermittelt :)
	    REA_sl=static_cast<int32>(cI-tagfield.begin())+1;
	    if(REA_sl>REA_sr) REA_sl=REA_sr;
	    updateClipoffs();
	  }
	}

	// right clip
	if(getRightClipoff()>0){
	  // first, search within right clip - tolerance
	  auto cI=tagfield.cbegin();
	  cI+=getRightClipoff();
	  found=false;
	  for(int32 c=0; c<tolerancelen && cI!=tagfield.begin(); ++c) {
	    if(*(--cI)==1) {
	      found=true;
	      break;
	    }
	  }
	  if(!found) {
	    cI=tagfield.begin();
	    for(int32 c=0; c<clipsearchlen && cI!=tagfield.begin(); ++c) {
	      if(*(--cI)==1) {
		found=true;
		break;
	      }
	    }
	  }
	  if(found) {
	    while(found) {
	      if(cI==tagfield.begin()) break;
	      if(*cI==0) {
		++cI;
		break;
	      }
	      --cI;
	    }
	    REA_sr=static_cast<int32>(cI-tagfield.begin());
	    if(REA_sr<REA_sl) REA_sr=REA_sl;
	    updateClipoffs();
	  }
	}
      }
    }
  }



  FUNCEND();
}




/*************************************************************************
 *
 * Sets the (quality) clipoffs to the given values. If force is true,
 *  then the sequence vector clipoffs are also set if needed.
 *
 *************************************************************************/

void Read::setClipoffs(int32 lclip, int32 rclip, bool force)
{
  FUNCSTART("void Read::setClipoffs(uint32 lclip, uint32 rclip, bool force)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(lclip<0, "lclip < 0?");
  BUGIFTHROW(rclip<0, "rclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(rclip>lenseq) rclip=lenseq;
  if(lclip>lenseq) lclip=lenseq;

  REA_ql=lclip;
  REA_qr=rclip;
  if(force==true){
    if(REA_sl>lclip) REA_sl=lclip;
    if(REA_sr<rclip) REA_sr=rclip;
  }

  updateClipoffs();

  FUNCEND();
}



/*************************************************************************
 *
 * Sets the lq quality clipoff to the given values
 *
 *************************************************************************/

void Read::setLQClipoff(int32 lclip)
{
  FUNCSTART("void Read::setLQClipoff(int32 lclip)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(lclip<0, "lclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(lclip>lenseq) lclip=lenseq;

  REA_ql=lclip;
  updateClipoffs();

  FUNCEND();
}

/*************************************************************************
 *
 * Sets the rq quality clipoff to the given values
 *
 *************************************************************************/

void Read::setRQClipoff(int32 rclip)
{
  FUNCSTART("void Read::setRQClipoff(int32 rclip)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(rclip<0, "rclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(rclip>lenseq) rclip=lenseq;

  REA_qr=rclip;
  updateClipoffs();

  FUNCEND();
}


/*************************************************************************
 *
 * Sets the l seq clipoff to the given values
 *
 *************************************************************************/

void Read::setLSClipoff(int32 lclip)
{
  FUNCSTART("void Read::setLSClipoffs(int32 lclip)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(lclip<0, "lclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(lclip>lenseq) lclip=lenseq;

  REA_sl=lclip;
  updateClipoffs();

  FUNCEND();
}

/*************************************************************************
 *
 * Sets the r seq clipoff to the given values
 *
 *************************************************************************/

void Read::setRSClipoff(int32 rclip)
{
  FUNCSTART("void Read::setRSClipoffs(int32 rclip)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(rclip<0, "rclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(rclip>lenseq) rclip=lenseq;

  REA_sr=rclip;
  updateClipoffs();

  FUNCEND();
}



/*************************************************************************
 *
 * Sets the masked (and quality) clipoffs to the given values.
 *
 *************************************************************************/

void Read::setLMClipoff(int32 lclip)
{
  FUNCSTART("int32 Read::setLMClipoff(int32 lclip)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(lclip<0, "lclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(lclip>lenseq) lclip=lenseq;

  REA_ml=lclip;
  REA_ql=lclip;
  updateClipoffs();
}

void Read::setRMClipoff(int32 rclip)
{
  FUNCSTART("int32 Read::setRMClipoff(int32 rclip)");

  BUGIFTHROW(checkRead()!=nullptr, checkRead());
  BUGIFTHROW(rclip<0, "rclip < 0?");

  auto lenseq=static_cast<int32>(getLenSeq());
  if(rclip>lenseq) rclip=lenseq;

  REA_mr=rclip;
  REA_qr=rclip;
  updateClipoffs();
}


/*************************************************************************
 *
 * Returns the number of bases the read could be extended to the left
 *  before the sequence vector or masked bases clipped part begins. 0 if none.
 *
 *************************************************************************/

int32 Read::getLeftExtend() const
{
  FUNCSTART("int32 Read::getLeftExtend() const");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  FUNCEND();

  int32 hardclip=std::max(REA_sl,REA_ml);
  if(REA_ql > hardclip) return REA_ql-hardclip;
  return 0;
}

/*************************************************************************
 *
 * Returns the number of bases the read could be extended to the right
 *  before the sequence vector or masked bases clipped part begins. 0 if none.
 *
 *************************************************************************/

int32 Read::getRightExtend() const
{
  FUNCSTART("int32 Read::getRightExtend() const");
  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  FUNCEND();

  int32 hardclip=std::min(REA_sr,REA_mr);
  if(REA_qr < hardclip) return hardclip-REA_qr;
  return 0;
}



// Routines for phred style quality<->errProbability conversion.

/*************************************************************************
 *
 *
 *
 *************************************************************************/

double Read::qualityToErrorRate_compute(base_quality_t qual)
{
  return pow(10.0,static_cast<double>(qual)/10.0);
}


// Qual is caped at 100
base_quality_t Read::errorRateToQuality(double errProb)
{
  //FUNCSTART("static const double Read::qualityToErrorProb(base_quality_t qual)");
  if(errProb>1000000000.0) errProb=1000000000.0;
  base_quality_t tmp=static_cast<base_quality_t>(log10(errProb)*10.0+.5);
  return tmp;

  //FUNCEND();
}


// Alert! These routines are not prepared for phred-style qualities!
//  (but work quite well anyway)

/*************************************************************************
 *
 *
 *
 *************************************************************************/

base_quality_t Read::queryAverageQualInSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const
{
  FUNCSTART("base_quality_t Read::queryAverageQualInSequence(int32 posl, int32 posr, bool skipNs, bool skipStars)");

  CEBUG("posl: " << posl << endl);
  CEBUG("posr: " << posr << endl);
  CEBUG("skipNs: " << skipNs << endl);
  CEBUG("skipStars: " << skipStars << endl);

  uint32 avgqual=0;
  uint32 countavg=0;
  base_quality_t retval=0;

  if(REA_ps_dirty==true && REA_pcs_dirty==false){
    retval=queryAverageQualInComplementSequence(
      static_cast<int32>(REA_padded_complementsequence.size())-posr-1,
      static_cast<int32>(REA_padded_complementsequence.size())-posl-1,
      skipNs,
      skipStars);
  }else{
    refreshPaddedSequence();

    if(posl<0) posl=0;
    if(posr >= static_cast<int32>(REA_padded_sequence.size())) posr=static_cast<int32>(REA_padded_sequence.size())-1;

    BOUNDCHECK(posl, 0, static_cast<int32>(REA_qualities.size()));
    auto qI=REA_qualities.cbegin();
    advance(qI, posl);

    BOUNDCHECK(posl, 0, static_cast<int32>(REA_padded_sequence.size()));
    auto sI=REA_padded_sequence.cbegin();
    advance(sI, posl);

    while((((*sI=='N'
	       || *sI=='n'
	       || *sI=='X'
	       || *sI=='x'
	       || *sI=='-') && skipNs)
	     || (*sI=='*' && skipStars))
	    && sI!=REA_padded_sequence.cbegin()){
	--sI; --qI;
    }

    auto tsI=REA_padded_sequence.cbegin();
    BOUNDCHECK(posr, 0, static_cast<int32>(REA_padded_sequence.size()));
    advance(tsI, posr);

    while((((*tsI=='N'
	       || *tsI=='n'
	       || *tsI=='X'
	       || *tsI=='x'
	       || *tsI=='-') && skipNs)
	     || (*tsI=='*' && skipStars))
	    && (tsI+1)!=REA_padded_sequence.end()){
	tsI++;
    }

    for(;sI<=tsI; ++sI, ++qI){
      CEBUG("char: " << *sI<< endl);
      if(dptools::isValidIUPACBase(*sI)){
	if(*qI!=BQ_UNDEFINED){
	  avgqual+=static_cast<uint32>(*qI);
	  countavg++;
	}
      }else{
	switch(toupper(*sI)){
	case '-':
	case 'N':
	case 'X':{
	  if(!skipNs && *qI!=BQ_UNDEFINED){
	    avgqual+=static_cast<uint32>(*qI);
	    countavg++;
	  }
	  break;
	}
	case '*':{
	  if(!skipStars && *qI!=BQ_UNDEFINED){
	    avgqual+=static_cast<uint32>(*qI);
	    countavg++;
	  }
	  break;
	}
	default:{
	  cout << "Illegal base: " << *sI << "(" << std::hex << static_cast<uint16>(*sI) << std::dec << ")" << endl;
	  MIRANOTIFY(Notify::FATAL, "Illegal base found: " << getName());
	}
	}
      }
    }
  }

  if(countavg>0) retval=static_cast<base_quality_t>(avgqual/countavg);

  CEBUG("Counted " << countavg << " non-skipped bases.\n");
  CEBUG("Returning: " << (uint16) retval << endl);

  FUNCEND();

  return retval;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

base_quality_t Read::queryAverageQualInClippedSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const
{
  FUNCSTART("base_quality_t Read::queryAverageQualInClippedSequence(int32 posl, int32 posr, bool skipNs, bool skipStars)");

  CEBUG("posl: " << posl << endl);
  CEBUG("posr: " << posr << endl);
  CEBUG("skipNs: " << skipNs << endl);
  CEBUG("skipStars: " << skipStars << endl);

  posl+=getLeftClipoff();
  posr+=getLeftClipoff();
  if(posl<getLeftClipoff()) posl=getLeftClipoff();
  if(posr>=getRightClipoff()) posr=getRightClipoff();

  FUNCEND();

  return queryAverageQualInSequence(posl, posr, skipNs, skipStars);
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

base_quality_t Read::queryAverageQualInComplementSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const
{
  FUNCSTART("base_quality_t Read::queryAverageQualInComplementSequence(int32 posl, int32 posr, bool skipNs, bool skipStars)");

  CEBUG("posl: " << posl << endl);
  CEBUG("posr: " << posr << endl);
  CEBUG("skipNs: " << skipNs << endl);
  CEBUG("skipStars: " << skipStars << endl);

  uint32 avgqual=0;
  uint32 countavg=0;
  base_quality_t retval=0;

  if(REA_pcs_dirty==true && REA_ps_dirty==false){
    retval=queryAverageQualInSequence(
      static_cast<int32>(REA_padded_complementsequence.size())-posr-1,
      static_cast<int32>(REA_padded_complementsequence.size())-posl-1,
      skipNs,
      skipStars);
  }else{
    refreshPaddedComplementSequence();

    if(posl<0) posl=0;
    if(posr >= static_cast<int32>(REA_padded_complementsequence.size())) posr=static_cast<int32>(REA_padded_complementsequence.size())-1;

    uint32 complement_posl=static_cast<uint32>(REA_padded_complementsequence.size())-posl-1;

    CEBUG("complement_posl: " << complement_posl << endl);

    BUGIFTHROW(complement_posl>=REA_qualities.size(),getName() << ": complement_posl (" << complement_posl << ") >= REA_qualities.size(" << REA_qualities.size() << ") ?");
    auto qI=REA_qualities.cbegin();
    advance(qI, complement_posl);

    BOUNDCHECK(posl, 0, static_cast<int32>(REA_padded_complementsequence.size()));
    auto sI=REA_padded_complementsequence.cbegin();
    advance(sI, posl);

    while((((*sI=='N'
	       || *sI=='n'
	       || *sI=='X'
	       || *sI=='x'
	       || *sI=='-') && skipNs)
	     || (*sI=='*' && skipStars))
	    && sI!=REA_padded_complementsequence.cbegin()){
	--sI; ++qI;
    }

    BOUNDCHECK(posr, 0, static_cast<int32>(REA_padded_complementsequence.size()));
    auto tsI=REA_padded_complementsequence.begin();
    advance(tsI, posr);

    while((((*tsI=='N'
	       || *tsI=='n'
	       || *tsI=='X'
	       || *tsI=='x'
	       || *tsI=='-') && skipNs)
	     || (*tsI=='*' && skipStars))
	    && (tsI+1)!=REA_padded_complementsequence.end()){
	++tsI;
    }

    for(;sI<=tsI; ++sI, --qI){
      CEBUG("char: " << *sI<< endl);
      if(dptools::isValidIUPACBase(*sI)){
	if(*qI!=BQ_UNDEFINED){
	  avgqual+=static_cast<int32>(*qI);
	  countavg++;
	}
      }else{
	switch(toupper(*sI)){
	case '-':
	case 'N':
	case 'X':{
	  if(!skipNs && *qI!=BQ_UNDEFINED){
	    avgqual+=static_cast<int32>(*qI);
	    countavg++;
	  }
	  break;
	}
	case '*':{
	  if(!skipStars && *qI!=BQ_UNDEFINED){
	    avgqual+=static_cast<int32>(*qI);
	    countavg++;
	  }
	  break;
	}
	default:{
	  cout << "Illegal base: " << *sI << "(" << std::hex << static_cast<uint16>(*sI) << std::dec << ")" << endl;
	  MIRANOTIFY(Notify::FATAL, "Illegal base found: " << getName());
	}
	}
      }
    }
  }

  if(countavg>0) retval=static_cast<base_quality_t>(avgqual/countavg);

  CEBUG("Counted " << countavg << " non-skipped bases.\n");
  CEBUG("Returning: " << (uint16) retval << endl);

  FUNCEND();

  return retval;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

base_quality_t Read::queryAverageQualInClippedComplementSequence(int32 posl, int32 posr, bool skipNs, bool skipStars) const
{
  FUNCSTART("base_quality_t Read::queryAverageQualInClippedComplementSequence(int32 posl, int32 posr, bool skipNs, bool skipStars)");

  CEBUG("posl: " << posl << endl);
  CEBUG("posr: " << posr << endl);
  CEBUG("skipNs: " << skipNs << endl);
  CEBUG("skipStars: " << skipStars << endl);

  FUNCEND();

  return queryAverageQualInComplementSequence(getLenSeq()-getRightClipoff()+posl,
						getLenSeq()-getRightClipoff()+posr,
						skipNs,
						skipStars);

}



/*************************************************************************
 *
 * Clip off bad qual at start and end of sequence
 *
 *
 * forward:
 *  -------F2-----------------------F1-----------    winlen, minqual (avg)
 *         >>                       >>
 *
 *  then
 *  -------F2-----RF1-----RF2-------F1-----------    winlen/2, minqual-5 (avg)
 *                <<<     <<<
 *  then
 *  -------F2-----RF1-X---RF2-------F1-----------
 *                >>>
 *
 * backward:
 *  -------F1-----------------------F2-----------
 *         <<                       <<
 *
 *  then
 *  -------F1-----RF2-----RF1-------F2-----------
 *                >>>     >>>
 *  then
 *  -------F1-----RF2--X--RF1-------F2-----------
 *                        <<<
 *
 *************************************************************************/

void Read::performQualityClip(uint32 minqual, uint32 winlen)
{
  FUNCSTART("void Read::performQualityClip(uint32 avgqual, uint32 winlen)");

  if(!REA_has_valid_data) {
    FUNCEND();
    return;
  }


  if(!REA_has_quality){
    FUNCEND();
    return;
  }

  if(minqual<5) minqual=5;
  if(winlen<10) winlen=10;

  double minquald=minqual;
  double winlend=winlen;

  uint32 winlen_2=winlen/2;
  uint32 minqualm5=minqual-5;
  double minqualm5d=minqualm5;
  double winlen_2d=winlen_2;

  bool exitprematurely=false;

  if(REA_qualities.size() > winlen) {
    uint32 qualtotal=0;
    auto f1I=REA_qualities.cbegin();
    for(uint32 i=0; i<winlen; i++, f1I++) {
      qualtotal+=*f1I;
    }
    auto f2I=REA_qualities.cbegin();
    double avgquald=static_cast<double>(qualtotal)/winlend;
    while(f1I != REA_qualities.end() && avgquald < minquald){
      qualtotal+=*f1I;
      qualtotal-=*f2I;
      avgquald=static_cast<double>(qualtotal)/winlend;
      ++f1I;
      ++f2I;
    }
    if(avgquald < minquald) {
      // Sheesh, we went through that Read without finding good quality,
      //  set the REA_ql and REA_qr accordingly
      REA_ql=static_cast<int32>(REA_qualities.size());
      REA_qr=0;
      exitprematurely=true;
    } else {
      // some good qual start found ... iterate a bit backward to find
      //  an acceptable start
      // halve winlen and take minqual -5
      if(f1I==REA_qualities.end()) f1I--;
      auto rf1I=f1I;
      qualtotal=0;
      for(uint32 i=0; i<winlen_2; i++, rf1I--) {
	qualtotal+=*rf1I;
      }
      auto rf2I=f1I;
      avgquald=static_cast<double>(qualtotal)/winlen_2d;
      while(rf1I!=REA_qualities.begin() && rf1I >= f2I && avgquald < minqualm5d){
	qualtotal+=*rf1I;
	qualtotal-=*rf2I;
	avgquald=static_cast<double>(qualtotal)/winlend;
	--rf1I;
	--rf2I;
      }
      // ok, average qual on halved window is now below minqual-5
      // reiterate forward until first base with qual >= minqual-5
      while(rf1I!=f1I && *rf1I<minqualm5) rf1I++;

      // rf1I now points to the first base of good region
      REA_ql=static_cast<int32>(rf1I-REA_qualities.begin());
    }

    // ok, we have start of good sequence at the front, now we'll look at the end
    if(!exitprematurely){
      qualtotal=0;
      f1I=REA_qualities.end()-1;
      for(uint32 i=0; i<winlen; i++, f1I--) {
	qualtotal+=*f1I;
      }
      f2I=REA_qualities.end()-1;
      avgquald=static_cast<double>(qualtotal)/winlend;
      while(f1I >= REA_qualities.begin() && avgquald < minquald){
	qualtotal+=*f1I;
	qualtotal-=*f2I;
	avgquald=static_cast<double>(qualtotal)/winlend;
	--f1I;
	--f2I;
      }
      if(avgquald < minquald) {
	// Sheesh, we went through that Read without finding good quality,
	//  this ... eh ... cannot be! (remember, we found forward!)
	updateClipoffs();
	refreshPaddedSequence();
	cout << *this;
	MIRANOTIFY(Notify::INTERNAL, "Quality clipping error: no backward cutoff, but forward found?");
      } else {
	if(f1I<REA_qualities.begin()) f1I=REA_qualities.begin();
	auto rf1I=f1I;
	qualtotal=0;
	for(uint32 i=0; i<winlen_2; i++, rf1I++) {
	  qualtotal+=*rf1I;
	}
	auto rf2I=f1I;
	avgquald=static_cast<double>(qualtotal)/winlen_2d;
	while(rf1I <= f2I && avgquald < minqualm5d){
	  qualtotal+=*rf1I;
	  qualtotal-=*rf2I;
	  avgquald=static_cast<double>(qualtotal)/winlend;
	  ++rf1I;
	  ++rf2I;
	}
	if(rf1I==REA_qualities.end()) rf1I--;
	while(rf1I!=f1I && *rf1I<minqualm5) rf1I--;

	// rf1I now points to the first base of good region
	REA_qr=static_cast<int32>(rf1I-REA_qualities.begin());
      }
    }
  }

  // it might happen that no good qualities were found,
  //  in this case, REA_ql=sizeofseq and REA_qr=0
  // we still want a read that won't throw errors if we ask it
  //  to give us the sequence and other things
  // so: set ql=qr=0
  if(REA_ql>REA_qr) REA_ql=REA_qr;

  updateClipoffs();

  FUNCEND();
  return;
}




/*************************************************************************
 *
 * physically trims the read to the left and right clips
 * tags completely outside new bounds are erased, others get positions adjusted
 *
 * reads of length 0 get to read of length 1, containing only 1 "N"
 *
 *************************************************************************/

void Read::performHardTrim()
{
  FUNCSTART("void Read::performHardTrim()");

  refreshPaddedSequence();

  //setCoutType(AS_TEXT);
  //cout << *this; cout.flush();

  // see whether to completely delete
  if(getLenClippedSeq()==0){
    static const char * nseq = "N";

    REA_tags.clear();
    setSequenceFromString(nseq);
  }else{
    // see whether to trim on right
    if(static_cast<int32>(getLenSeq()) - getRightClipoff() >0){
      BUGIFTHROW(REA_padded_sequence.size()<getRightClipoff()-1,"1 REA_padded_sequence.size()<getRightClipoff()-1 ? " << REA_padded_sequence.size() << " " << getRightClipoff());
      REA_padded_sequence.resize(getRightClipoff()-1);
      REA_padded_complementsequence.clear();
      REA_pcs_dirty=true;
      if(!REA_qualities.empty()){
	BUGIFTHROW(REA_qualities.size()<getRightClipoff()-1,"REA_qualities.size()<getRightClipoff()-1 ?");
	REA_qualities.resize(getRightClipoff()-1);
      }
      if(!REA_adjustments.empty()){
	BUGIFTHROW(REA_adjustments.size()<getRightClipoff()-1,"REA_adjustments.size()<getRightClipoff()-1 ?");
	REA_adjustments.resize(getRightClipoff()-1);
      }
      if(REA_rlevptr!=nullptr){
	BUGIFTHROW(REA_rlevptr->size()<getRightClipoff()-1,"REA_rlevptr->size()<getRightClipoff()-1 ?");
	REA_rlevptr->resize(getRightClipoff()-1);
      }

      if(!REA_bposhashstats.empty()){
	BUGIFTHROW(REA_bposhashstats.size()<getRightClipoff()-1,"REA_bposhashstats.size()<getRightClipoff()-1 ?");
	REA_bposhashstats.resize(getRightClipoff()-1);
      }

      REA_qr=REA_padded_sequence.size();
      REA_sr=REA_padded_sequence.size();
      REA_cr=REA_padded_sequence.size();
      REA_mr=REA_padded_sequence.size();

      auto tI=REA_tags.begin();
      for(uint32 tagpos=0; tI != REA_tags.end(); ++tagpos){
	if(tI->from >= REA_padded_sequence.size()) {
	  tI=REA_tags.erase(tI);
	  --tagpos;
	} else {
	  if(tI->to > REA_padded_sequence.size()) tI->to=REA_padded_sequence.size();
	  ++tI;
	}
      }
    }


    // see whether to trim on left
    if(getLeftClipoff()>0){
      BUGIFTHROW(REA_padded_sequence.size()<getRightClipoff()-1,"2 REA_padded_sequence.size()<getRightClipoff()-1 ?");
      {
	auto bI=REA_padded_sequence.begin();
	auto eI=bI;
	advance(eI,getLeftClipoff());
	REA_padded_sequence.erase(bI,eI);
      }
      REA_padded_complementsequence.clear();
      REA_pcs_dirty=true;
      if(!REA_qualities.empty()){
	BUGIFTHROW(REA_qualities.size()<getRightClipoff()-1,"REA_qualities.size()<getRightClipoff()-1 ?");
	auto bI=REA_qualities.begin();
	auto eI=bI;
	advance(eI,getLeftClipoff());
	REA_qualities.erase(bI,eI);
      }
      if(!REA_adjustments.empty()){
	BUGIFTHROW(REA_adjustments.size()<getRightClipoff()-1,"REA_adjustments.size()<getRightClipoff()-1 ?");
	auto bI=REA_adjustments.begin();
	auto eI=bI;
	advance(eI,getLeftClipoff());
	REA_adjustments.erase(bI,eI);

	for(auto & ae : REA_adjustments){
	  if(ae >=0 ) ae-=getLeftClipoff();
	}
      }
      if(REA_rlevptr!=nullptr){
	BUGIFTHROW(REA_rlevptr->size()<getRightClipoff()-1,"REA_rlevptr->size()<getRightClipoff()-1 ?");
	auto bI=REA_rlevptr->begin();
	auto eI=bI;
	advance(eI,getLeftClipoff());
	REA_rlevptr->erase(bI,eI);
      }
      if(!REA_bposhashstats.empty()){
	BUGIFTHROW(REA_bposhashstats.size()<getRightClipoff()-1,"REA_bposhashstats.size()<getRightClipoff()-1 ?");
	auto bI=REA_bposhashstats.begin();
	auto eI=bI;
	advance(eI,getLeftClipoff());
	REA_bposhashstats.erase(bI,eI);
      }

      auto tI=REA_tags.begin();
      for(uint32 tagpos=0; tI != REA_tags.end(); ++tagpos){
	if(tI->to < getLeftClipoff()) {
	  tI=REA_tags.erase(tI);
	  --tagpos;
	} else {
	  tI->from-=getLeftClipoff();
	  tI->to-=getLeftClipoff();
	  ++tI;
	}
      }

      REA_qr-=getLeftClipoff();
      REA_sr-=getLeftClipoff();
      REA_cr-=getLeftClipoff();
      REA_mr-=getLeftClipoff();

      REA_ql=0;
      REA_sl=0;
      REA_cl=0;
      REA_ml=0;
    }

  }

  FUNCEND();
}



/*************************************************************************
 *
 * deletes the weakest base of a base run
 * can also insert a gap with quality 0 at the *original* position
 *  given in the run (helps to delete a weak base also in aligned reads
 *  without breaking the alignment)
 *
 *************************************************************************/

void Read::deleteWeakestBaseInRun(const char base, const uint32 position, const bool insertgap)
{
  FUNCSTART("");

  uint32 p1=position;
  bool found=getPosOfWeakestBaseInRun(base,p1);
  if(found){
    deleteBaseFromSequence(p1);
    if(insertgap) {
      insertBaseInSequence('*',0,position,true);
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 * Why give a char base to the function if it simply could read out
 * the base from its sequence?
 *
 * Easy: imagine ..AAA*TTT... and we gave the position of the gap character
 *  as only parameter? (this can happen during assembly in later stages)
 *
 *
 *************************************************************************/

bool Read::getPosOfWeakestBaseInRun(char base, uint32 & position)
{
  FUNCSTART("");

  if(!REA_has_valid_data) {
    FUNCEND();
    return false;
  }

  // force a check of quality status
  //paranoiaBUGIF(setQualityFlag(),);

  BUGIFTHROW(position>=getLenSeq(), getName() << ": position (" << position << ") >= size of read (" << getLenSeq() << ") ?");

  base=static_cast<char>(toupper(base));

  refreshPaddedSequence();

  base_quality_t minqual=101;
  uint32 minqualpos=position;
  bool foundmin=false;

  // check stretch of char base to left
  // if several bases have the same loq qual, takes farthest right
  {
    int32 rrpi=position-1;

    for(; rrpi>=0; rrpi--){
      if(REA_padded_sequence[rrpi] != '*') {
	if(toupper(REA_padded_sequence[rrpi]) == base) {
	  if(REA_qualities[rrpi] < minqual) {
	    minqual=REA_qualities[rrpi];
	    minqualpos=rrpi;
	    foundmin=true;
	  }
	} else {
	  break;
	}
      }
    }
  }

  // check stretch of char base to right
  // if several bases have the same loq qual, takes farthest right
  {
    int32 rrpi=position;

    for(; rrpi < static_cast<int32>(REA_padded_sequence.size()); rrpi++){
      if(REA_padded_sequence[rrpi] != '*'){
	if(toupper(REA_padded_sequence[rrpi]) == base) {
	  if(REA_qualities[rrpi] <= minqual) {
	    minqual=REA_qualities[rrpi];
	    minqualpos=rrpi;
	    foundmin=true;
	  }
	}else{
	  break;
	}
      }
    }
  }

  if(foundmin){
    position=minqualpos;
    FUNCEND();
    return true;
  }

  FUNCEND();
  return false;
}




/*************************************************************************
 *
 * extract read type and strain names from MINF or MIT2 tags
 *
 * deletes the tags after processing
 *
 * currently, israil is always false (as not written to CAF anyway)
 *
 * returns false if no MINF/MIT2 tags were encountered; else true
 *
 *************************************************************************/

bool Read::extractMINFTagInfo(std::vector<multitag_t> & tags, const std::string & readname, std::string & dummy_strainname, std::string & dummy_seqtypename, std::string & dummy_machinetype, int8 & dummy_tplacementcode, bool & dummy_isbb, bool & dummy_israil, bool & dummy_isCER)
{
  FUNCSTART("void Read::extractMINFTagInfo(std::vector<multitag_t> & tags, const std::string & readname, std::string & dummy_strainname, std::string & dummy_seqtypename, std::string & dummy_machinetype, int8 & dummy_tplacementcode, bool & dummy_isbb, bool & dummy_isCER)");

  bool retval=false;
  dummy_strainname.clear();
  dummy_seqtypename.clear();
  dummy_machinetype.clear();
  dummy_tplacementcode=ReadGroupLib::SPLACE_UNKNOWN;
  dummy_isbb=false;
  dummy_israil=false;
  dummy_isCER=false;

  auto tI=tags.begin();
  while(tI!=tags.end()){
    if(tI->identifier==REA_tagentry_idMINF){
      retval=true;
      dummy_strainname.clear();
      dummy_seqtypename.clear();
      dummy_machinetype.clear();
      dummy_isbb=false;
      dummy_israil=false;
      dummy_isCER=false;

      // extract values from tag comment
      // they're not GenBank comments, but stored the same way
      //  let's misuse the function then :-)
      std::string value;
      tI->extractGenBankKeyValueFromComment("ST",dummy_seqtypename);

      tI->extractGenBankKeyValueFromComment("SN",dummy_strainname);

      tI->extractGenBankKeyValueFromComment("MT",dummy_machinetype);

      tI->extractGenBankKeyValueFromComment("BB",value);
      if(!value.empty() && value!="0"){
	dummy_isbb=true;
      }
      tI->extractGenBankKeyValueFromComment("RR",value);
      if(!value.empty() && value!="0"){
	dummy_israil=true;
      }

      tI->extractGenBankKeyValueFromComment("CER",value);
      //cout << tI-> comment << "\nmsrextract: " << value << endl;
      if(!value.empty() && value!="0"){
	dummy_isCER=true;
      }

      // Ok, extracted everything. delete this tag
      tI=tags.erase(tI);
    }else if(tI->identifier==REA_tagentry_idMIT2){
      retval=true;
      dummy_strainname.clear();
      dummy_seqtypename.clear();
      dummy_machinetype.clear();
      dummy_isbb=false;
      dummy_israil=false;
      dummy_isCER=false;

      std::string tmpdecode;
      std::vector<std::string> attributes;

      // TODO: adapt to readgroups

      attributes.reserve(6);
      boost::split(attributes, tI->getCommentStr(), boost::is_any_of(";"),boost::token_compress_on);
      std::vector<std::string> keyvalue;
      keyvalue.reserve(2);
      for(const auto & ae : attributes) {
	keyvalue.clear();
	boost::split(keyvalue, ae, boost::is_any_of("="));
	if(keyvalue.size()==2){
	  if(keyvalue[0] == "st"){
	    dummy_seqtypename=keyvalue[1];
	  }else if(keyvalue[0] == "sn"){
	    gff3Decode(keyvalue[1],tmpdecode);
	    dummy_strainname=tmpdecode;
	  }else if(keyvalue[0] == "mt"){
	    gff3Decode(keyvalue[1],tmpdecode);
	    dummy_machinetype=tmpdecode;
	  }else if(keyvalue[0] == "bb"){
	    if(!keyvalue[1].empty() && keyvalue[1]!="0"){
	      dummy_isbb=true;
	    }
	  }else if(keyvalue[0] == "rr"){
	    if(!keyvalue[1].empty() && keyvalue[1]!="0"){
	      dummy_israil=true;
	    }
	  }else if(keyvalue[0] == "pc"){
	    if(!ReadGroupLib::parseSegmentPlacement(keyvalue[1],dummy_tplacementcode)){
	      MIRANOTIFY(Notify::FATAL, "Error in " << readname << ": key type 'pc' has unknown value in MIT2 tag: " << tI->getCommentStr() << "\n");
	    }
	  }else if(keyvalue[0] == "cer"){
	    if(!keyvalue[1].empty() && keyvalue[1]!="0"){
	      dummy_isCER=true;
	    }
	  }else{
	    MIRANOTIFY(Notify::FATAL, "Error in " << readname << " with unknown key type '" << keyvalue[0] << "' in MIT2 tag: " << tI->getCommentStr() << "\n");
	  }
	}
      }

      // Ok, extracted everything. delete this tag
      tI=tags.erase(tI);
    }else{
      ++tI;
    }
  }

  FUNCEND();
  return retval;
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint32 Read::deleteTag(const multitag_t::mte_id_t identifier)
{
  FUNCSTART("uint32 Read::deleteTag(const multitag_t::mte_id_t identifier)");

  paranoiaBUGIF(checkRead()!=nullptr, MIRANOTIFY(Notify::FATAL, checkRead()));

  uint32 deleted=0;
  auto tI=REA_tags.begin();
  while(tI!=REA_tags.end()){
    if(identifier==tI->identifier){
      tI=REA_tags.erase(tI);
      ++deleted;
    }else{
      ++tI;
    }
  }

  FUNCEND();
  return deleted;
}


/*************************************************************************
 *
 * Exchanges Ns in a read with gaps, takling care of adjustments and
 *  qualities
 *
 * Convenience function for Contig object
 *
 *************************************************************************/

void Read::exchangeNsWithGaps()
{
  FUNCSTART("void Read::exchangeNsWithGaps()");
  if(!REA_has_valid_data) {
    FUNCEND();
    return;
  }

  // TODO: perhaps also work with reverse? but not needed now
  refreshPaddedSequence();

  auto cI=REA_padded_sequence.begin();
  auto aI=REA_adjustments.begin();
  auto qI=REA_qualities.begin();
  for(; cI!=REA_padded_sequence.end(); ++cI, ++qI){
    if(*cI=='N'){
      REA_pcs_dirty=true;
      *cI='*';
      if(REA_uses_adjustments) *aI=-1;
      switch(getSequencingType()){
      case ReadGroupLib::SEQTYPE_SANGER :
      case ReadGroupLib::SEQTYPE_SOLEXA :
      case ReadGroupLib::SEQTYPE_ABISOLID : {
	base_quality_t newqual=0;
	base_quality_t poslooked=0;
	if(cI!=REA_padded_sequence.begin()){
	  newqual=*(qI-1);
	  poslooked++;
	}
	if((cI+1)!=REA_padded_sequence.end()){
	  // gcc 4.3.2 warns as it think the values could grow
	  //  higher than 255 ... they can't, base quals go to max
	  //  of 100.
	  //newqual+=*(qI+1);
	  // writing instead the version with cast to get gcc quiet
	  newqual=static_cast<base_quality_t>(newqual+(*(qI+1)));
	  poslooked++;
	}
	if(poslooked) {
	  *qI=static_cast<base_quality_t>(newqual/poslooked);
	}else{
	  // should not happen
	  *qI=0;
	}
	break;
      }
      case ReadGroupLib::SEQTYPE_454GS20 : {
	*qI=1;
      }
      default: {
	*qI=0;
      }
      }
    }
    // manually increasing aI. not in for loop as adjustments may be empty
    if(REA_uses_adjustments) aI++;
  }

  FUNCEND();
}







/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Read::blindSeqData(char base)
{
  FUNCSTART("void Read::blindSeqData(char base)")
  if(!REA_has_valid_data) {
    FUNCEND();
    return;
  }

  if(!REA_ps_dirty){
    for(auto & pse : REA_padded_sequence){
      if(pse != '*' && pse != 'n' && pse != 'N') pse=base;
    }
  }

  if(!REA_pcs_dirty){
    base=dptools::getComplementIUPACBase(base);
    for(auto & pcse : REA_padded_complementsequence){
      if(pcse != '*' && pcse != 'n' && pcse != 'N') pcse=base;
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Read::upDownCase(base_quality_t threshold)
{
  FUNCSTART("void Read::upDownCase(base_quality_t threshold)");
  if(!REA_has_valid_data) {
    FUNCEND();
    return;
  }

  // TODO: perhaps also work with reverse? but not needed now
  refreshPaddedSequence();

  auto cI=REA_padded_sequence.begin();
  auto qI=REA_qualities.cbegin();
  for(; cI!=REA_padded_sequence.cend(); ++cI, ++qI){
    *cI=static_cast<char>(toupper(*cI));
    if(*qI<threshold) *cI=static_cast<char>(tolower(*cI));
  }

  REA_pcs_dirty=true;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Read::upDownCaseClips()
{
  FUNCSTART("void Read::upDownCaseClips()");
  if(!REA_has_valid_data) {
    FUNCEND();
    return;
  }

  // TODO: perhaps also work with reverse? but not needed now
  refreshPaddedSequence();

  auto lc=getLeftClipoff();
  auto rc=getRightClipoff();
  int32 si=0;
  for(auto & pse : REA_padded_sequence){
    if(si>=lc && si<rc){
      pse=static_cast<char>(toupper(pse));
    }else{
      pse=static_cast<char>(tolower(pse));
    }
    ++si;
  }

  REA_pcs_dirty=true;

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Read::loadDataFromEXP(std::string filename, std::string path)
{
  FUNCSTART("void Read::loadDataFromEXP()");

  REA_has_valid_data=false;

  EXP tmpexp;

  boost::trim(path);
  boost::trim(filename);
  std::string fullpath(path);
  if(!fullpath.empty()){
    fullpath+="/";
  }
  fullpath+=filename;

  tmpexp.load(fullpath.c_str());

  // ok, get the name and SCFname from the exp file
  {
    // First, the name of the experiment
    //REA_name.clear();
    //REA_name=tmpexp.getID();
    setName(tmpexp.getID());
    if(getName().empty()){
      MIRANOTIFY(Notify::FATAL, "The experiment has no name (no ID field found): " << filename);
    }

    // Now the SCF filename of the experiment
    // TODO: readgroups. SCF filename
    //REA_scf_filename.clear();
    //REA_scf_filename=tmpexp.getLN();

    //if(REA_scf_filename.empty()){
    //  // FIXME: how to treat exp files without SCF data?
    //  //throw Notify(Notify::FATAL, THISFUNC, REA_exp_filename, ": The experiment has no SCF file (no LN field found)");
    //}
  }


  REA_ql=tmpexp.getQL();
  REA_qr=tmpexp.getQR();
  REA_sl=tmpexp.getSL();
  REA_sr=tmpexp.getSR();
  REA_cl=tmpexp.getCL();
  REA_cr=tmpexp.getCR();

  // TODO: adapt to readgroups
  //REA_machine_type=REA_sc_machine_type.addEntry(tmpexp.getMA());
  //setInsize(tmpexp.getSIfrom(),tmpexp.getSIto());

  // TODO: reactivate/emulate
  //tmpexp.swapTags(REA_tags);      // the EXP will be destructed, so save time

  // subtract 1 from the positions read in exp file (which are basis 1)
  for(uint32 i=0; i<REA_tags.size(); i++ ){
    REA_tags[i].from--;
    REA_tags[i].to--;
  }

  if(REA_ql>REA_qr){
#ifndef PUBLICQUIET
    WARNING("WARNING: " << getName() << "   left quality clip > right quality clip (QL > QR), adjusted.");
#endif
    REA_qr=REA_ql;
  }
  if(REA_sl>REA_sr){
#ifndef PUBLICQUIET
    WARNING("WARNING: " << getName() << "   sequence vector left clip > sequence vector right clip (SL > SR), adjusted.");
#endif
    REA_sr=REA_sl;
  }

  updateClipoffs();

  if(getRightClipoff()-getLeftClipoff()<0){
#ifndef PUBLICQUIET
    WARNING("WARNING: " << getName() << "   quality and sequence vector clippings form a negative length (clip right - clip left <0), adjusted");
#endif
    REA_qr=REA_sl;
    updateClipoffs();
  }

  // initialise the sequence vector
  {
    const char * exp_sequence=tmpexp.getSequence().c_str();
    if(exp_sequence==nullptr){
      MIRANOTIFY(Notify::WARNING, "No bases were found: " << getName());
    }

    REA_padded_sequence.clear();
    REA_padded_sequence.reserve(strlen(exp_sequence)+3);
    do{
      char tmpc=*exp_sequence++;
      if(dptools::isValidIUPACBase(tmpc)){
	REA_padded_sequence.push_back(tmpc);
      }else{
	switch(toupper(tmpc)){
	case 'N':
	case 'X':
	case '*':{
	  REA_padded_sequence.push_back(tmpc);
	  break;
	}
	case '-': {
	  REA_padded_sequence.push_back('N');
	  break;
	}
	case 0 : {
	  // FIXME: how to treat emtpy reads????
	  //throw Notify(Notify::WARNING, THISFUNC, REA_exp_filename, ": No bases were found?");
	  // cannot throw, isn't gentle enough for big projects :(
	  WARNING("WARNING: " << tmpexp.getID() << "   has no bases in the experiment file?");
	  // FIXME: pfusch!
	  REA_padded_sequence.push_back('N');
	  REA_ql=0;
	  REA_qr=1;
	  REA_sl=0;
	  REA_sr=1;
	  updateClipoffs();
	  exp_sequence--;
	  break;
	}
	default: {
	  cout << "Illegal base: " << tmpc << "(" << std::hex << static_cast<uint16>(tmpc) << std::dec << ")" << endl;
	  MIRANOTIFY(Notify::FATAL, "Illegal base found: " << tmpexp.getID());
	}
	}
      }

    }while(*exp_sequence!=0);
  }

  postLoadEXPFASTA();


  REA_template=tmpexp.getTN();

  // TODO: adapt to readgroups
//  if(tmpexp.getSV().empty()){
//    REA_seqvec_name=REA_sc_seqvec_name.addEntry(tmpexp.getSF());
//  }else{
//    REA_seqvec_name=REA_sc_seqvec_name.addEntry(tmpexp.getSV());
//  }
//  REA_basecaller=REA_sc_basecaller.addEntry(tmpexp.getBC());
//  REA_asped=REA_sc_asped.addEntry(tmpexp.getDT());
//  //  REA_dye;
//  REA_processstatus=REA_sc_processstatus.addEntry(tmpexp.getPS());
//  REA_primer=REA_sc_primer.addEntry(tmpexp.getPN());

  //  read AV qualities if any
  {
    const std::vector<base_quality_t> & equal=tmpexp.getAV();
    uint32 expavsize=static_cast<uint32>(equal.size());
    if(expavsize>0 && expavsize!=REA_padded_sequence.size()){
      cout << "Sequence length: " << REA_padded_sequence.size() << endl;
      cout << "Number of quality values: " << expavsize << endl;
      MIRANOTIFY(Notify::FATAL, "The experiment has an unequal number of bases and quality values (SQ vs AV fields): " << filename);
    }
    if(expavsize>0){
      REA_has_quality=true;
      for(uint32 i=0; i<REA_qualities.size(); i++){
	REA_qualities[i]=equal[i];
      }
    }
  }

  //  read ON adjustments
  if(REA_uses_adjustments){
    const std::vector<int32> onvals= tmpexp.getON();
    if(onvals.size()){
      if(onvals.size()%2) {
	MIRANOTIFY(Notify::INTERNAL, "there should be a even number of ON adjustments: " << filename);
      }
      //std::vector<int32> newadjustments[REA_adjustments.size()];
      std::vector<int32> newadjustments;
      newadjustments.resize(REA_adjustments.size());
      uint32 nai=0;
      uint32 ovi=0;
      while(ovi<onvals.size()){
	//cout << onvals[ovi] << "  " << onvals[ovi+1] << endl;
	if(onvals[ovi]==0){
	  newadjustments[nai]=-1;
	  nai++;
	} else {
	  uint32 onlow=onvals[ovi];
	  uint32 onhigh=onvals[ovi+1];
	  for(uint32 i=onlow; i<=onhigh; i++, nai++){
	    if(nai>=newadjustments.size()) {
	      MIRANOTIFY(Notify::FATAL, "the adjustments given in the ON tag exceed the size of the sequence: " << filename);
	    }
	    newadjustments[nai]=i-1;
	  }
	}
	ovi+=2;
      }
      if(nai!=newadjustments.size()){
	cout << "nai: " << nai << "  newadjustments.size(): " << newadjustments.size() << endl;
	MIRANOTIFY(Notify::FATAL, "the adjustments given in the ON tag do not cover the whole range of sequence: " << filename);
      }
      REA_adjustments.swap(newadjustments);
    }
  }

  //// transfer Mira INFormation in MINF tags to read variables
  //transferMINFTagsToReadInfo();
  {
    std::string dummy_strainname,dummy_seqtype,dummy_machinetype;
    int8 dummy_tplacementcode;
    bool dummy_issbb, dummy_israil, dummy_isCER;
    extractMINFTagInfo(REA_tags,getName(),
		       dummy_strainname,dummy_seqtype,dummy_machinetype,
		       dummy_tplacementcode,
		       dummy_issbb, dummy_israil, dummy_isCER);
  }

  REA_has_valid_data=true;

  REA_ml=0;
  REA_mr=getLenSeq();

  if(REA_qr<static_cast<int32>(REA_qualities.size())) REA_qr--;
  if(REA_sr<static_cast<int32>(REA_qualities.size())) REA_sr--;
  if(REA_cr<static_cast<int32>(REA_qualities.size())) REA_cr--;
  updateClipoffs();

  //  CEBUG(*this);
  BUGIFTHROW(checkRead()!=nullptr, checkRead());

  //tmpexp.dump();
  //cout << *this;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

int32 Read::getDigiNormMultiplier() const
{
  FUNCSTART("int32 Read::getDigiNormMultiplier() const");

  static std::string notestr("Note");

  int32 retvalue=1;

  std::string dnrstr;
  for(auto & te : REA_tags){
    // Note: must also check whether DGNr tag covers the complete current sequence
    // because sometimes reads get cut back in parts (e.g. those covered by DGNr)
    //  and this then gives a completely wrong impression.
    if(te.identifier==REA_tagentry_idDGNr
       && !te.getCommentStr().empty()
       && getLeftClipoff() >= te.from
       && getRightClipoff()-1 <= te.to){
      if(te.commentisgff3){
	retvalue=atoi(GFFParse::extractKeytag(notestr,te.getCommentStr()).c_str());
      }else{
	retvalue=atoi(te.getCommentStr().c_str());
      }
      break;
    }
  }

  BUGIFTHROW(retvalue==0,getName() << ": Found diginorm multiplier of 0?");

  return retvalue;
}
