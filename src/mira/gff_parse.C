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

// for boost::trim, split
#include <boost/algorithm/string.hpp>

#include <regex>

#include "mira/gff_parse.H"

#include "io/annotationmappings.H"
#include "errorhandling/errorhandling.H"

#include "util/progressindic.H"
#include "mira/readpool.H"
#include "mira/multitag.H"

using std::cout;
using std::cerr;
using std::endl;

#define CEBUG(bla)


// Plain vanilla constructor
GFFParse::GFFParse(ReadPool * rp)
{
  FUNCSTART("GFFParse::GFFParse(ReadPool * rp)");

  zeroVars();
  init();

  GFFP_readpool=rp;

  FUNCEND();
}

void GFFParse::zeroVars()
{
  FUNCSTART("void GFFParse::zeroVars()");

  GFFP_readpool=nullptr;

  FUNCEND();
}

void GFFParse::init()
{
  FUNCSTART("void GFFParse::init()");
  FUNCEND();
}



GFFParse::~GFFParse()
{
  FUNCSTART("GFFParse::~GFFParse()");

  discard();

  FUNCEND();
}


void GFFParse::discard()
{
  FUNCSTART("GFFParse::discard()");

  zeroVars();

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//GFFParse::GFFParse(GFFParse const &other)
//{
//  FUNCSTART("GFFParse::GFFParse(GFFParse const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//GFFParse const & GFFParse::operator=(GFFParse const & other)
//{
//  FUNCSTART("GFFParse const & GFFParse::operator=(GFFParse const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, GFFParse const &???)
//{
//  FUNCSTART("friend ostream & GFFParse::operator<<(ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::string & GFFParse::getSequenceName(uint32 i) const
{
  FUNCSTART("const std::string & GFFParse::getSequenceName(uint32 i) const");

  if(i>=GFFP_seqnames.size()){
    MIRANOTIFY(Notify::WARNING, ": Tried to get out of range sequence name.");
  }

  FUNCEND();
  return GFFP_seqnames[i];
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::string & GFFParse::getSequence(uint32 i) const
{
  FUNCSTART("const std::string & GFFParse::getSequence(uint32 i) const");

  if(i>=GFFP_sequences.size()){
    MIRANOTIFY(Notify::WARNING, ": Tried to get out of range sequence.");
  }

  FUNCEND();
  return GFFP_sequences[i];
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::vector<multitag_t> & GFFParse::getTags(uint32 i) const
{
  FUNCSTART("const std::vector<tag_t> & GFFParse::getTags(uint32 i) const");

  if(i>=GFFP_sequencetags.size()){
    MIRANOTIFY(Notify::WARNING, ": Tried to get out of range tags.");
  }

  FUNCEND();
  return GFFP_sequencetags[i];
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/



// substr vector passed by reference: avoid constructing vector n times
void GFFParse::parseNormalGFFLine(const std::string & line, const uint64 lineno, std::vector<std::string> & substrs)
{
  substrs.clear();
  boost::split(substrs, line, boost::is_any_of("\t"));

  if(substrs.size() != 9) {
    cout << "Line " << lineno << ": expected 9 elements, found " << substrs.size() << "\nBad line: " << line << "\n";
    return;
  }

  auto snI=GFFP_snmap.find(substrs[0]);
  size_t snmindex=0;
  if(snI != GFFP_snmap.end()){
    snmindex=snI->second;
  }else{
    GFFP_snmap[substrs[0]]=GFFP_seqnames.size();
    snmindex=GFFP_seqnames.size();
    GFFP_seqnames.push_back(substrs[0]);
    GFFP_sequences.resize(GFFP_sequences.size()+1);
    GFFP_sequencetags.resize(GFFP_sequences.size()+1);
  }

  GFFP_sequencetags[snmindex].resize(GFFP_sequencetags[snmindex].size()+1);

  // add tag defined by this line to sequence just found
  multitag_t & newtag=GFFP_sequencetags[snmindex].back();

  newtag.source=multitag_t::MT_tagsrcentry_idGFF3;
  if(substrs[1]!="."){
    newtag.setSourceStr(substrs[1]);
  }

  newtag.from=atoi(substrs[3].c_str());
  newtag.to=atoi(substrs[4].c_str());

  if(newtag.from==0){
    cout << "Line " << lineno
	 << ": position 'from' (field 4," << newtag.from << ") is 0? Coordinates in GFF files should have 1 as lowest value.\n";
    if(GFFP_errorstatus<2) GFFP_errorstatus=2;
  }else{
    --newtag.from;
  }
  if(newtag.to==0){
    cout << "Line " << lineno
	 << ": position 'to' (field 5," << newtag.to << ") is 0? Coordinates in GFF files should have 1 as lowest value.\n";
    if(GFFP_errorstatus<2) GFFP_errorstatus=2;
  }else{
    --newtag.to;
  }

  if(!substrs[2].empty()){
    newtag.identifier=multitag_t::newIdentifier(substrs[2]);
    if(!AnnotationMappings::isValidGFF3SOEntry(substrs[2])){
      if(AnnotationMappings::translateOldSOfeat2SOfeat(substrs[2]).empty()){
	cout << "Line " << lineno
	     << ": MIRA does not know type '" << substrs[2] << "' in column 3 of the GFF3 fle\n";
	if(GFFP_errorstatus<1) GFFP_errorstatus=1;
      }else{
	newtag.identifier=multitag_t::newIdentifier(AnnotationMappings::translateOldSOfeat2SOfeat(substrs[2]));
      }
    }
  }

  if(substrs[6].empty()){
    cout << "Line " << lineno
	 << ": field 7 may only be '+', '-', '.' or '?', but found empty string\n";
    newtag.setStrand('=');
    if(GFFP_errorstatus<2) GFFP_errorstatus=2;
  }else if(substrs[6]!="+" && substrs[6]!="-"){
    if(substrs[6]!="." && substrs[6]!="?"){
      cout << "Line " << lineno
	   << ": direction in field 7 may only be '+', '-', '.' or '?', but found '"
	   << substrs[6] << "'\n";
      if(GFFP_errorstatus<2) GFFP_errorstatus=2;
    }
    newtag.setStrand('=');
  }else{
    newtag.setStrand(substrs[6][0]);
  }

  if(substrs[7].size() != 1){
    cout << "Line " << lineno
	 << ": for CDS, field 8 may only be '0', '1' or '2', but found empty string\n";
    if(GFFP_errorstatus<2) GFFP_errorstatus=2;
  }else{
    switch(substrs[7][0]){
    case '0' : { newtag.phase=0; break; }
    case '1' : { newtag.phase=1; break; }
    case '2' : { newtag.phase=2; break; }
    case '.' : { newtag.phase=3; break; }
    default : {
      cout << "Line " << lineno
	   << ": field 8 may only be '0', '1' , '2' or '.', but found string '"
	   << substrs[7] << "'\n";
      if(GFFP_errorstatus<2) GFFP_errorstatus=2;
    }
    }
  }

  std::string setcomment;
  std::string tmpcomment;

  if(!substrs[8].empty() || !substrs[5].empty()){
    if(substrs[5].empty() || substrs[5] == "."){
      tmpcomment=substrs[8];
    }else{
      if(substrs[8].empty()){
	tmpcomment="gff3sco="+substrs[5];
      }else{
	tmpcomment=substrs[8]+";gff3sco="+substrs[5];
      }
    }
  }

  std::string miraitag;
  {
    std::string g3source;
    char   g3strand;
    uint8  g3phase;
    extractMIRAGFF3InfoFromGFF3Attributes(tmpcomment,
					  setcomment,
					  g3source,
					  g3strand,
					  g3phase,
					  miraitag);
  }

  if(!miraitag.empty()) {
    newtag.identifier=multitag_t::newIdentifier(miraitag);
  }

  newtag.comment=multitag_t::newComment(setcomment);


  newtag.commentisgff3=true;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void GFFParse::loadFile(const std::string & filename)
{
  FUNCSTART("void GFFParse::loadFile(const std::string & filename)");

  //size_t numseqsloaded=0;
  uint64 lineno=0;

  GFFP_errorstatus=0;

  std::ifstream gffin(filename, std::ios::in|std::ios::ate);
  if(!gffin) {
    MIRANOTIFY(Notify::FATAL, "GFF file not found for loading:" << filename);
  }
  if(!gffin.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "GFF file empty? " << filename);
  }

  ProgressIndicator<std::streamsize> P(0, gffin.tellg(),1000);

  gffin.seekg(0, std::ios::beg);

  std::string actline;
  actline.reserve(10000);
  std::vector<std::string> substrs;
  substrs.reserve(9);

  bool fastamode=false;
  while(!gffin.eof()){
    if(gffin.eof()) break;
    getline(gffin,actline);
    ++lineno;

    if(!actline.empty()){
      if(actline[0]=='#') {
	if(actline == "##FASTA") {
	  fastamode=true;
	  break;
	}
      }else{
	parseNormalGFFLine(actline,lineno,substrs);
	if(GFFP_errorstatus>0) {
	  MIRANOTIFY(Notify::FATAL,"Ooops?! Please have a look at that line:\n" << actline << "\n");
	}
      }
    }
    if(P.delaytrigger()) P.progress(gffin.tellg());
  }
  if(fastamode){
    std::string blanks=" \t\n";
    std::string tmpseq;
    std::string tmpsname;
    bool saveseq=false;
    while(!gffin.eof()){
      if(gffin.eof()) break;
      getline(gffin,actline);
      ++lineno;

      if(!actline.empty()){
	if(actline[0]=='>') {
	  // check name of sequence in sequence name map
	  // if present, save current sequence
	  // if not, do not save as this is either a protein or a bogus file
	  if(saveseq){
	    auto snI=GFFP_snmap.find(tmpsname);
	    BUGIFTHROW(snI == GFFP_snmap.end(),"snI == GFFP_snmap.end() ???");
	    //cout << "Must save " << tmpsname << "\t" << snI->second << endl;
	    GFFP_sequences[snI->second].swap(tmpseq);
	  }
	  tmpseq.clear();
	  auto tokenend=actline.find_first_of(blanks,0);
	  if(tokenend==std::string::npos) tokenend=actline.size();
	  tmpsname=actline.substr(1, tokenend-1);
	  CEBUG("tmpsname: " << tmpsname << endl);

	  saveseq=false;
	  // see whether we will need to save this sequence
	  auto snI=GFFP_snmap.find(tmpsname);
	  if(snI != GFFP_snmap.end()){
	    saveseq=true;
	  }else{
	    // sequence name not found ... yet
	    // it may be a NCBI format line though, so check for that
	    if(tmpsname.size() > 3
	       && tmpsname[0]=='g'
	       && tmpsname[1]=='i'
	       && tmpsname[2]=='|'){
	      // OK, this *might* be a GI line
	      std::vector<std::string> subnames;
	      boost::split(subnames, tmpsname, boost::is_any_of("|"));
	      if(subnames.size()==5){
		// bloody well looks like an NCBI style name line. Extract the real sequence name
		swap(tmpsname,subnames[3]);
		saveseq=true;
	      }
	    }
	  }
	  CEBUG("Seeing new seqname ###"<<tmpsname<<"###\n");
	}else{
	  // only spend time if it's a sequence we will keep
	  // (and not a protein or something)
	  if(saveseq){
	    // append trimmed actline to tmpseq
	    boost::trim(actline);
	    tmpseq+=actline;
	  }
	}
      }
      if(P.delaytrigger()) P.progress(gffin.tellg());
    }

    // there might be some unsaved sequences still
    if(saveseq){
      auto snI=GFFP_snmap.find(tmpsname);
      cout << "Must save " << tmpsname << "\t" << snI->second << endl;
      GFFP_sequences[snI->second].swap(tmpseq);
    }
  }

  P.finishAtOnce();
  cout << '\n';

  gffin.close();

  if(GFFP_errorstatus==0) {
    cout << "Loading finished OK, checking tags\n";
    checkTagsOnceLoaded();
  }else{
    cout << "Loading finished not OK: " << GFFP_errorstatus << endl;
  }

  if(GFFP_errorstatus > 0){
    if(GFFP_errorstatus == 1){
      cout << "GFF file '" << filename << "' had errors (see output above), but they seem minor.\n";
    }else{
      MIRANOTIFY(Notify::FATAL,"GFF file '" << filename << "' had unrecoverable errors (see output above). Fix your file!\n");
    }
  }

  FUNCEND();
}

/*************************************************************************
 *
 * Side effect: if readpool given and the GFF had annotations for a read
 *  in that readpool, then a name index is created in the read pool
 *
 *************************************************************************/

void GFFParse::checkTagsOnceLoaded()
{
  FUNCSTART("void GFFParse::checkTagsOnceLoaded()");

  bool errorsfound=false;
  for(size_t snmindex=0; snmindex<GFFP_seqnames.size(); ++snmindex){
    if(GFFP_sequences[snmindex].empty() && !GFFP_sequencetags.empty()){
      if(GFFP_readpool!=nullptr){
	GFFP_readpool->allowNameIndex(true);
	auto rid=GFFP_readpool->getReadIndex(GFFP_seqnames[snmindex]);
	if(rid>=0) continue;
	cout << "Sequence " << GFFP_seqnames[snmindex] << " has elements defined on sequence, but no sequence and none of the previously loaded sequences matches. Please make sure that files with sequence data belonging to " << GFFP_seqnames[snmindex] << " were loaded before the annotation in this GFF file.\n";
      }else{
	cout << "Sequence " << GFFP_seqnames[snmindex] << " has elements defined on sequence, but no sequence?\n";
      }
      errorsfound=true;
      if(GFFP_errorstatus<2) GFFP_errorstatus=2;
    }
  }

  if(errorsfound) {
    if(GFFP_readpool!=nullptr){
      GFFP_readpool->allowNameIndex(false);
    }
    return;
  }

  for(size_t snmindex=0; snmindex<GFFP_seqnames.size(); ++snmindex){
    std::vector<multitag_t>::iterator mtI=GFFP_sequencetags[snmindex].begin();
    //cout << "checking " << GFFP_seqnames[snmindex] <<endl;
    size_t seqsize=GFFP_sequences[snmindex].size();
    if(seqsize==0){
      auto rid=GFFP_readpool->getReadIndex(GFFP_seqnames[snmindex]);
      BUGIFTHROW(rid < 0,"rid < 0 not expected here at this stage");
      seqsize=GFFP_readpool->getRead(rid).getLenSeq();
    }
    for(; mtI != GFFP_sequencetags[snmindex].end(); ++mtI){
      if(mtI->from>=seqsize){
	errorsfound=true;
	cout << "Sequence " << GFFP_seqnames[snmindex]
	     << ": position 'from' (field 4, value " << mtI->from+1 << ") is larger than the sequence size (" << seqsize << "," << GFFP_seqnames[snmindex] << "): " << mtI->getCommentStr() << '\n';
	if(GFFP_errorstatus<2) GFFP_errorstatus=2;
      }

      if(mtI->to>=seqsize){
	errorsfound=true;
	cout << "Sequence " << GFFP_seqnames[snmindex]
	     << ": position 'to' (field 5, value " << mtI->to+1 << ") is larger than the sequence size (" << seqsize << "," << GFFP_seqnames[snmindex] << "): " << mtI->getCommentStr() << '\n';
	if(GFFP_errorstatus<2) GFFP_errorstatus=2;
      }
    }
  }

  FUNCEND();
}




/*************************************************************************
 *
 * TODO: better / faster with regex!
 *
 *
 *************************************************************************/

bool GFFParse::checkCommentForGFF3(const std::string & comment)
{
  FUNCSTART("void GFFParse::checkCommentForGFF3(const std::string & comment)");

  if(comment.empty()) return false;

  for(auto & s : GFFP_gff3scankeys) {
    auto fpos=comment.find(s);
    if(fpos!=std::string::npos) return true;
  }
  return false;
}




/*************************************************************************
 *
 * In:
 *   src = attributes string
 * Out:
 *   dst = attributes string minus the encodes attributes
 *   source, strand, phase
 *
 * Note: attribute score (gff3sco) remains in string!
 * special: if strand info was not found, return value of strand == '*'
 *
 *************************************************************************/

void GFFParse::extractMIRAGFF3InfoFromGFF3Attributes(const std::string & src, std::string & dst, std::string & source, char & strand, uint8 & phase, std::string & miraitag)
{
  FUNCSTART("void GFFParse::extractMIRAGFF3InfoFromGFF3Attributes(const std::string & src, std::string & dst, std::string & source, char & strand, uint8 & phase)");

  dst.clear();
  source.clear();
  miraitag.clear();
  strand='*';
  phase=3;

  if(src.empty()) return;

  CEBUG("Working on " << src << endl);
  std::vector<std::string> attributes;
  attributes.reserve(20);
  boost::split(attributes, src, boost::is_any_of(";"),boost::token_compress_on);

  std::vector<std::string> keyvalue;
  keyvalue.reserve(2);
  for(const auto & ae : attributes){
    CEBUG("Doing " << ae << endl);
    keyvalue.clear();
    boost::split(keyvalue, ae, boost::is_any_of("="));

    CEBUG("kv:");
    for(const auto & kve : keyvalue){
      CEBUG(" " << *kvI);
    }
    CEBUG(endl);

    BUGIFTHROW(keyvalue.size()>2, "Found " << keyvalue.size()-1 << " '=' signs while trying to parse " << ae << " in attributes of comment " << src << "\nthis should not be.");
    if(keyvalue.size()==2){
      if(keyvalue[0] == "gff3str"){
	boost::trim(keyvalue[1]);
	if(!keyvalue[1].empty()) strand=keyvalue[1][0];
      }else if(keyvalue[0] == "gff3pha"){
	boost::trim(keyvalue[1]);
	if(!keyvalue[1].empty()) {
	  phase=keyvalue[1][0];
	  switch(phase){
	  case '1':{
	    phase=0;
	    break;
	  }
	  case '2':{
	    phase=1;
	    break;
	  }
	  case '3':{
	    phase=2;
	    break;
	  }
	  default :{
	    phase=3;
	  }
	  }
	}
      }else if(keyvalue[0] == "gff3src"){
	source.swap(keyvalue[1]);
      }else if(keyvalue[0] == "miraitag"){
	miraitag.swap(keyvalue[1]);
      }else{
	if(!dst.empty()){
	  dst+=";";
	}
	dst+=ae;
      }
    }
  }

  FUNCEND();
  return;
}




/*************************************************************************
 *
 * In:
 *   attributestr = attributes string
 * Out:
 *   parsed = attributes parsed into gff3attributes_t, the value strings
 *            decoded
 *
 *
 *************************************************************************/

//#define CEBUG(bla)  {cout << bla; cout.flush();}
void GFFParse::parseGFF3Attributes(const std::string & attributestr, gff3attributes_t & parsed)
{
  parsed.clear();
  if(attributestr.empty()) return;

  std::vector<std::string> splitatts;
  splitatts.reserve(20);
  boost::split(splitatts, attributestr, boost::is_any_of(";"),boost::token_compress_on);
  if(splitatts.empty()) return;

  gff3tagvalue_t tmptv;
  gff3values_t   tmpgff3vals;

  std::vector<std::string> keyvalue;
  keyvalue.reserve(2);
  for(const auto & sae : splitatts){
    CEBUG("Doing " << sae << endl);
    keyvalue.clear();
    boost::split(keyvalue, sae, boost::is_any_of("="));
    // there should always be exactly one '='
    // if not, it's an illegal tag/value in the GFF3 attributes
    if(keyvalue.size()!=2) {
      if(keyvalue.size()==0) {
	cout << "In attributes '" << attributestr << "': empty tag/value? Very strange ...\n";
      }else if(keyvalue.size()==1) {
	cout << "In attributes '" << attributestr << "': encountered tag '" << keyvalue[0] << "' but no value? Not allowed by GFF3 standard, skipping.\n";
      }else{
	cout << "In attributes '" << attributestr << "': encountered tag '" << keyvalue[0] << "' and the value field behind contains " << keyvalue.size()-1 << " '='-signs. Not allowed by GFF3 standard, skipping.\n";
      }
      continue;
    }

    parsed.push_back(tmptv);
    parsed.back().tag=keyvalue[0];
    tmpgff3vals.clear();
    CEBUG("Whatcha do on: " << keyvalue[1] << endl);
    boost::split(tmpgff3vals, keyvalue[1], boost::is_any_of(","));
    CEBUG("tmpgff3vals.size(): " << tmpgff3vals.size() << endl);
    for(const auto & g3ve : tmpgff3vals){
      parsed.back().values.resize(parsed.back().values.size()+1);
      gff3Decode(g3ve,parsed.back().values.back());
      CEBUG("gff3decoded: " << g3ve << "\ninto       : " << parsed.back().values.back());
    }
  }

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * In:
 *   parsed = attributes as gff3attributes_t
 * Out:
 *   attributestr = attributes string, properly encoded
 *
 *
 *************************************************************************/
void GFFParse::createGFF3AttributeString(const gff3attributes_t & attributes, std::string & attributestr)
{
  attributestr.clear();
  std::string tmpcode;
  tmpcode.reserve(128);

  auto aI=attributes.begin();
  for(uint32 acounter=0; aI != attributes.end(); ++aI, ++acounter){
    if(acounter>0) attributestr+=';';
    attributestr+=aI->tag;
    attributestr+='=';
    auto vI=aI->values.cbegin();
    for(uint32 vcounter=0; vI != aI->values.cend(); ++vI, ++vcounter){
      if(vcounter>0) attributestr+=',';
      gff3Code(*vI,tmpcode);
      attributestr+=tmpcode;
    }
  }
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

std::string GFFParse::extractCommonName(const std::string & attributestr, bool extended)
{
  static const std::vector<std::string> esearches = {
    "name", "alias"
  };
  static const std::vector<std::string> ssearches = {
    "Name", "Alias", "gene", "locus_tag"
  };
  static const std::regex ipmatch("[\\.:]");

  std::string retvalue;
  //for(auto sI=searches.begin(); sI!=searches.end(); ++sI){

  for(auto & s : ssearches){
    retvalue=extractKeytag(s,attributestr);
    if(!retvalue.empty()) break;
  }
  if(retvalue.empty()){
    for(auto & s : esearches){
      retvalue=extractKeytag(s,attributestr);
      if(!retvalue.empty()) break;
    }
    if(!retvalue.empty() && !extended){
      if(std::regex_search(retvalue, ipmatch)){
	retvalue.clear();
      }
    }
  }

  return retvalue;
}


/*************************************************************************
 *
 * e.g.: given GFF3 string with ...;Name=XXXXX;...
 *       returns XXXXX (decoded) when queried with key=="Name"
 *
 *************************************************************************/

std::string GFFParse::extractKeytag(const std::string & key, const std::string & attributestr)
{
  auto reK=GFFP_regex_extractkeys.find(key);
  if(reK==GFFP_regex_extractkeys.end()){
    std::regex tmp;
    GFFP_regex_extractkeys[key]=tmp;
    reK=GFFP_regex_extractkeys.find(key);
    //reK->second="[^|;]"+key+"=(.*?)[;|$]";

    // Non-capture group
    // The "?:" in the second parenthesis means that the parenthesis just overrides the
    //  precedence (the "|" is weak), but does NOT capture
    reK->second="(?:^|;)"+key+"=(.*?)(?:;|$)";
  }

  std::string retvalue;
  std::match_results<std::string::const_iterator> matches;
  if(std::regex_search(attributestr, matches, reK->second)){
    gff3Decode(matches[1],retvalue);
  }

  return retvalue;
}
