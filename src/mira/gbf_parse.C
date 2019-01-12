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

#include <boost/algorithm/string.hpp>

#include "mira/gbf_parse.H"

#include "io/annotationmappings.H"
#include "util/progressindic.H"
#include "util/misc.H"
#include "mira/multitag.H"


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)

const char * GBF::GBF_miragbfscankeys[]= {
  "/gene",
  "/locus_tag",
  "/product",
  "/note",
  "/function",
  "/EC_number",
  "/operon",
  "/protein_id",
  "/codon_start",
  "/transl_table",
  "/translation",
  "/standard_name",
  ""
};

// add pair here if you want to have it written out to GFF3 format later
const char * GBF::GBF_oldmiragbf2gff3translations[]= {
  "/gene", "gene=",
  "/locus_tag", "locus_tag=",
  "/product", "product=",
  "/note", "Note=",
  "/function", "function=",
  "/EC_number", "eC_number=",
  "/operon", "operon=",
  "/protein_id", "protein_id=",
  "/codon_start", "codon_start=",
  "/transl_table", "transl_table=",
  "/translation", "translation=",
  "/standard_name", "standard_name=",

  // BaCh 23.08.2017
  "/inference","inference=",
  "/regulatory_class","regulatory_class=",
  "/gene_synonym","gene_synonym=",

  // BaCh 23.08.2017: These would need to be translated to, e.g, Ontology_term="GO:0046703" ... but I don't feel like it atm
  "/GO_process","gO_process=",
  "/GO_function","gO_function=",
  "/GO_component","gO_component=",
  ""
};

GBF::strstrmap GBF::GBF_mapoldmiragbf2gff3;

// keep this last
const bool GBF::GBF_staticfeaturesinit=GBF::staticInitialiser();


bool GBF::staticInitialiser()
{
  FUNCSTART("bool GBF::staticInitialiser()");

  for(uint32 i=0; GBF_oldmiragbf2gff3translations[i][0] != 0; i+=2) {
    GBF_mapoldmiragbf2gff3[GBF_oldmiragbf2gff3translations[i]]=GBF_oldmiragbf2gff3translations[i+1];
  }

  FUNCEND();
  return true;
}


// Plain vanilla constructor
GBF::GBF()
{
  FUNCSTART("GBF::GBF()");

  zeroVars();
  init();

  FUNCEND();
}

void GBF::zeroVars()
{
  FUNCSTART("void GBF::zeroVars()");

  GBF_sequencenames.clear();
  GBF_sequences.clear();
  GBF_tags.clear();

  FUNCEND();
}

void GBF::init()
{
  FUNCSTART("void GBF::init()");
  FUNCEND();
}



GBF::~GBF()
{
  FUNCSTART("GBF::~GBF()");

  discard();

  FUNCEND();
}


void GBF::discard()
{
  FUNCSTART("GBF::discard()");

  zeroVars();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::string & GBF::getSequenceName(uint32 i) const
{
  FUNCSTART("const std::string & GBF::getSequenceName(uint32 i) const");

  if(i>=GBF_sequencenames.size()){
    MIRANOTIFY(Notify::WARNING, ": Tried to get out of range sequence name.");
  }

  FUNCEND();
  return GBF_sequencenames[i];
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::string & GBF::getSequence(uint32 i) const
{
  FUNCSTART("const std::string & GBF::getSequence(uint32 i) const");

  if(i>=GBF_sequencenames.size()){
    MIRANOTIFY(Notify::WARNING, ": Tried to get out of range sequence.");
  }

  FUNCEND();
  return GBF_sequences[i];
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const std::vector<multitag_t> & GBF::getTags(uint32 i) const
{
  FUNCSTART("const std::vector<tag_t> & GBF::getTags(uint32 i) const");

  if(i>=GBF_sequencenames.size()){
    MIRANOTIFY(Notify::WARNING, ": Tried to get out of range tags.");
  }

  FUNCEND();
  return GBF_tags[i];
}

/*************************************************************************
 *
 * Copies keys in "thingstotransfer" from /gene features
 *  to all other features
 *
 * Expects that /gene features precedes the /otherfeature in the GBF
 *
 *************************************************************************/

void GBF::transferGeneInfoToCDSInfo()
{
  FUNCSTART("void GBF::transferGeneInfoToCDSInfo()");

//  BUGIFTHROW(GBF_sequences.size()!=GBF_sequencenames.size()
//	     || GBF_sequencenames.size() != GBF_tags.size(),
//	     "GBF data read is fishy, unequal number of sequences, names and tag vectors?");
//
//  std::vector<std::string> thingstotransfer;
//  thingstotransfer.push_back("/gene");
//  thingstotransfer.push_back("/locus_tag");
//
//  for(uint32 i=0; i<GBF_tags.size(); i++){
//    auto J=GBF_tags[i].begin();
//    for(uint32 j=1; j< GBF_tags[i].size(); j++){
//      if(GBF_tags[i][j].identifier!="Fgen") {
//	if(GBF_tags[i][j].from==GBF_tags[i][j-1].from
//	   && GBF_tags[i][j].to==GBF_tags[i][j-1].to) {
//	  for(uint32 t=0; t<thingstotransfer.size(); t++){
//	    auto fpos=GBF_tags[i][j].comment.find(thingstotransfer[t],0);
//	    if(fpos==std::string::npos){
//	      fpos=GBF_tags[i][j-1].comment.find(thingstotransfer[t],0);
//	      if(fpos!=std::string::npos){
//		auto lpos=GBF_tags[i][j-1].comment.find('\n',0);
//		if(lpos==std::string::npos) lpos=GBF_tags[i][j-1].comment.size();
//		std::string newcomment=GBF_tags[i][j-1].comment.substr(fpos, lpos-fpos);
//		if(GBF_tags[i][j].comment.size()){
//		  if(newcomment[newcomment.size()-1]!='\n') newcomment+='\n';
//		  newcomment+=GBF_tags[i][j].comment;
//		}
//		GBF_tags[i][j].comment=newcomment;
//	      }
//	    }
//	  }
//	}
//      }
//    }
//  }

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


void GBF::load(const std::string & gbfname)
{
  FUNCSTART("GBF::load(const std::string & gbfname)");

  discard();

  std::ifstream fin(gbfname, std::ios::in|std::ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::WARNING, "File not found: " << gbfname);
  }

  auto lenfile=fin.tellg();
  // static cast needed for gcc 3.3.3
  if(lenfile==0){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << gbfname);
  }
  fin.seekg(0, std::ios::beg);

  loadTheFile(fin, lenfile);

  correctForTagPositionErrors(gbfname);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void GBF::correctForTagPositionErrors(const std::string & gbfname)
{
  FUNCSTART("void GBF::correctForTagPositionErrors()");

  BUGIFTHROW(GBF_sequences.size()!=GBF_sequencenames.size()
	     || GBF_sequencenames.size() != GBF_tags.size(),
	     "GBF data read is fishy, unequal number of sequences, names and tag vectors?");

  bool unrecoverable=false;
  for(uint32 i=0; i<GBF_tags.size(); i++){
    size_t seqlen=GBF_sequences[i].size();
    if(seqlen==0){
      MIRANOTIFY(Notify::FATAL,"Entry " << GBF_sequencenames[i] << " in file " << gbfname << " has no DNA sequence??? Did you download an incomplete GenBank file without 'ORIGIN' entry (aka the DNA sequence)?");
    }
    for(auto & thistag : GBF_tags[i]){
      if(thistag.from > seqlen){
	cout << "Fishy tag which starts completely outside the sequence:\n" << thistag;
	unrecoverable=true;
      }
      if(thistag.to > seqlen) {
	cout << "Fishy tag (to is greater than sequence length of " << seqlen << "):\n" << thistag;
	thistag.to=static_cast<uint32>(seqlen);
      }
    }
  }

  if(unrecoverable){
    MIRANOTIFY(Notify::FATAL,"Some data in the loaded GenBank file had unrecoverable errors (see log above), please fix the file.");
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

// silly gcc 4.3: if lentoreserve is "streampos" type, there's a
// conversion error

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void GBF::loadTheFile(std::ifstream & fin, uint64 lentoreserve)
{
  FUNCSTART("void GBF::loadTheFile(std::ifstream & fin, streampos lentoreserve)");

  ProgressIndicator<uint64> P(0, lentoreserve,1000);

  std::string::size_type linepos;
  std::string::size_type postokenstart;

  uint32 numloci=0;

  //tag_t emptytag;
  //emptytag.from=0;
  //emptytag.to=0;
  //emptytag.strand='=';
  //emptytag.identifier="";
  //emptytag.comment="";
  //
  //tag_t acttag=emptytag;

  std::string guessedname;

  std::string tagidentifier;
  std::string tagcomment;

  std::string actkey;
  std::string actval;
  std::vector<int32> fromtopos;

  std::string line;
  std::string token="";
  bool haveLocus=false;
  bool haveFeature=false;
  bool haveOrigin=false;
  uint64 linecount=0;
  while (getline(fin,line,'\n')){
    ++linecount;
    if(P.delaytrigger()) P.progress(fin.tellg());

    // get rid of '\r' from DOS
    while(!line.empty() && line[line.size()-1]=='\r') line.resize(line.size()-1);

    linepos=0;
    getNextToken(line, linepos, token, postokenstart);
    CEBUG("line: " << line);
    CEBUG("token\t" << token << endl);
    CEBUG("linepos: " << linepos << endl);
    if(token.size()==0) continue;
    if(token=="//") {
      haveLocus=false;
      haveFeature=false;
      haveOrigin=false;

      guessedname.clear();
      tagidentifier.clear();
      tagcomment.clear();

      actkey.clear();
      actval.clear();
      fromtopos.clear();
    } else if(haveLocus && haveOrigin) {
      // in ORIGIN
      // act token contains the linenumbers ... not interesting to us

      CEBUG("In ORIGIN: line: " << line<<endl);
      std::string & actseq=GBF_sequences.back();
      for(;linepos<line.size(); linepos++){
	if(!isblank(line[linepos])) {
	  actseq+=line[linepos];
	  //CEBUG("Tadaaa: " << line[linepos] << endl);
	}
      }
    } else if(haveLocus && haveFeature) {
      // in feature
      if(postokenstart==0 && token=="ORIGIN") {
	haveOrigin=true;
	addKeyVal2TagComment(actkey,actval,tagcomment);
	storeAccumulatedTags(guessedname, tagidentifier, tagcomment, fromtopos);
      }else if(postokenstart==0 && token=="BASE") {
	// DDBJ file: BASE COUNT ... just ignore that line
      } else if(postokenstart==0){
	cout << "Fishy line (" << linecount << "): " << line << "\nViolating GenBank standard guidelines: we are in features, token is not ORIGIN, but starts at column 0.\n";
      } else if(postokenstart==5 && linepos<=20){
	// new feature
	CEBUG("New feature " << postokenstart << " " << linepos << endl);
	addKeyVal2TagComment(actkey,actval,tagcomment);
	storeAccumulatedTags(guessedname, tagidentifier, tagcomment, fromtopos);

	guessedname.clear();
	tagidentifier.clear();
	tagcomment.clear();

	actkey.clear();
	actval.clear();

	std::string location;
	getNextToken(line, linepos, location, postokenstart);
	// handle non-compliance of CloneManager which has a blank between complement and bracket: "complement (x..y)"
	// *sigh*
	// in those cases, location only contains "complement" and not "complement(...)"
	// make this a general resolution (in case other elements are also badly writen by CloneManager)
	if(location.back()!=')'){
	  auto newloc=location;
	  getNextToken(line, linepos, location, postokenstart);
	  location=newloc+location;
	}
	fromtopos.clear();
	parseGBFLocation(location,fromtopos,1);

	tagidentifier=token;
      } else {
	// feature continuation line
	if(token[0]=='/') {
	  // new key
	  CEBUG("New key : " << token << endl);

	  addKeyVal2TagComment(actkey,actval,tagcomment);

	  actkey.clear();
	  actval.clear();

	  auto equpos=token.find('=');
	  if(equpos!=std::string::npos && equpos != token.size()-1){
	    actkey=token.substr(0,equpos);
	    actval=token.substr(equpos+1,10000000)+line.substr(linepos,1000000);;

	    // *sigh* workaround for Vector NTI "/label"
	    auto tI=GBF_mapoldmiragbf2gff3.find(actkey);
	    if(tI==GBF_mapoldmiragbf2gff3.end() && actkey!="/label"){
	      actkey.clear();
	      actval.clear();
	    }else{
	      if(actkey!="/label"){
		actkey=tI->second;
	      }else{
		// transform the "\" from Vector NTI to regular spaces *sigh*
		actkey="label=";
		for(auto & c : actval){
		  if(c=='\\') c=' ';
		}
	      }
	      // BaCh 05.04.2013: Hmmmm ... this is double coding (addKeyVal2TagComment() already does it)
	      // Why was it in?
	      //
              //std::string tmp;
	      //gff3Code(actval,tmp);
	      //tmp.swap(actval);
	      if(actkey=="gene="){
		boost::trim_if(actval,boost::is_any_of("\""));
		guessedname=actval;
		actkey.clear();
		actval.clear();
	      }else{
		if(guessedname.empty()
		   && (actkey=="locus_tag="
		       || actkey=="protein_id")){
		  boost::trim_if(actval,boost::is_any_of("\""));
		  guessedname=actval;
		}
	      }
	    }
	  }
	}else{
	  // key continuation line
	  // additionally insert a space to anything but /translation
	  CEBUG("cont " << token << endl);
	  if(!actkey.empty()){
	    if(actkey!="/translation") actval+=' ';
	    actval+=token;
	    if(linepos!=std::string::npos) actval+=line.substr(linepos,100000);
	  }
	}
      }
    } else if(haveLocus) {
      // after locus
      if(token=="FEATURES") {
	haveFeature=true;
      }else if(token=="ORIGIN") {
	haveOrigin=true;
      }
    } else {
      // nothing here yet
      if(token=="LOCUS") {
	haveLocus=true;
	numloci++;
	std::string seqname;
	getNextToken(line, linepos,seqname, postokenstart);
	GBF_sequencenames.push_back(seqname);
	GBF_sequences.push_back("");
	GBF_sequences.back().reserve(lentoreserve);
	GBF_tags.resize(numloci);
	GBF_tags.back().reserve(10000);
      }else{
	MIRANOTIFY(Notify::FATAL, "Missing LOCUS token as first entry in file ... are you sure that this is a Genbank file?");
      }
    }
  }

  P.finishAtOnce();
  cout << endl;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void GBF::addKeyVal2TagComment(const std::string & actkey, const std::string & actval, std::string & comment) const
{
  if(!actkey.empty()) {
    if(!comment.empty()) comment+=';';
    comment+=actkey;
    std::string tmp(boost::trim_copy_if(actval,boost::is_any_of("\"")));

    // also remove "-characters in the string.
    // While technically allowed in GBK files using double quotes (""label""),
    //  they are a pain in the ass for CAF format. caf2gap e.g. has problems.
    // Replace "" with '. Quite rare, so in-place replace is OK.
    boost::replace_all(tmp,"\"\"","'");
    // and single double quotes. Albeit not allowed, they seem to exist here and there.
    boost::replace_all(tmp,"\"","'");

    std::string tmp2;
    gff3Code(tmp,tmp2);
    comment+=tmp2;
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


void GBF::storeAccumulatedTags(const std::string & guessedname, const std::string & identifier, const std::string & comment, std::vector<int32> & fromto)
{
  FUNCSTART("void GBF::storeAccumulatedTags(tag_t & acttag, std::vector<int32> & fromto)");

  if(!identifier.empty()){
    std::string namecomment;
    if(guessedname.empty()){
      namecomment=comment;
    }else{
      namecomment="Name="+guessedname+";"+comment;
    }
    if(fromto.size()%3 != 0) {
      MIRANOTIFY(Notify::INTERNAL, "Could not parse feature location in GBF file. This also might be an FATAL error due to buggy GBF!");
    }

    std::string soidentifier(AnnotationMappings::translateGenBankfeat2SOfeat(identifier));

    auto I=fromto.cbegin();
    for(;I!=fromto.cend(); I+=3) {
      multitag_t acttag;

      acttag.setSourceStr("GenBank");
      if(soidentifier.empty()){
	acttag.setIdentifierStr("located_sequence_feature");
      }else{
	acttag.setIdentifierStr(soidentifier);
      }
      acttag.setCommentStr(namecomment);
      acttag.commentisgff3=true;

      if(*I>0) {
	acttag.setStrand('+');
      }else{
	acttag.setStrand('-');
      }
      // tags have positions beginning at 0 in MIRA
      acttag.from=*(I+1)-1;
      acttag.to=*(I+2)-1;

      GBF_tags.back().push_back(acttag);

      CEBUG("Stored: "<< acttag << endl);
    }
  }

#if 0
  if(!acttag.identifier.empty()){
    //cout << "Write tag(s):"<<endl;
    //cout << "identifier: " << acttag.identifier << endl;
    //cout << "comment: " << acttag.comment << endl;
    //cout << "for ranges: ";
    if(fromto.size()%3 != 0) {
      MIRANOTIFY(Notify::INTERNAL, "Could not parse feature location in GBF file. This also might be an FATAL error due to buggy GBF!");
    }
    auto ftI=fromto.cbegin();
    for(;ftI!=fromto.cend(); ftI+=3) {
      if(*ftI>0) {
	acttag.strand='+';
      }else{
	acttag.strand='-';
      }
      // tags have positions beginning at 0 in MIRA
      acttag.from=*(ftI+1)-1;
      acttag.to=*(ftI+2)-1;
      //cout << "Storing: ";
      //acttag.dump();

      //cout << "strand: " << acttag.strand << endl;
      //cout << "from: " << acttag.from << "\tto: " << acttag.to << endl;
    }
  }

#endif

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void GBF::getNextToken(const std::string & line, std::string::size_type & linepos, std::string & token, std::string::size_type & postokenstart) const
{
  token="";
  if(linepos>=line.size()) {
    linepos=std::string::npos;
    return;
  }

  std::string blanks=" \t";

  postokenstart=std::string::npos;
  auto tokenend=std::string::npos;

  postokenstart=line.find_first_not_of(blanks,linepos);
  if(postokenstart==std::string::npos){
    linepos=std::string::npos;
    return;
  }
  tokenend=line.find_first_of(blanks,postokenstart);
  if(tokenend==std::string::npos) tokenend=line.size();

  token=line.substr(postokenstart, tokenend-postokenstart);
  linepos=tokenend;
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void GBF::parseGBFLocation(const std::string & location, std::vector<int32> & fromto, int32 direction) const
{
  if(location.empty()) return;

  //cout << "GBF need to parse loc: #" << location << "#\n";

  if(location[0]=='c'
     || location[0]=='j'
     || location[0]=='o') {
    if(location[0]=='c') direction*=-1;
    auto leftbracket=location.find('(');
    if(leftbracket==std::string::npos) return;
    auto rightbracket=getCorrespondingRightBracket(location, leftbracket);

    parseGBFLocation(location.substr(leftbracket+1,rightbracket-leftbracket-1),
		     fromto,
		     direction);
    rightbracket++;
    if(rightbracket<location.size())
      parseGBFLocation(location.substr(rightbracket+1,1000000),
		       fromto,
		       1);
  } else {
    std::string::size_type locpos=0;
    long firstnum=0;
    long secondnum=0;
    if(location[0]=='<') locpos++;
    if(location[0]=='>') locpos++; //??
    if(location[0]==',') locpos++;
    if(location[0]=='('){
      //cout << "Analyse bracket val 1st: " << location << endl;
      auto rightbracket=getCorrespondingRightBracket(location, 0);
      parseGBFLocation(location.substr(1,rightbracket-1),
		       fromto,
		       direction);
//      werte poppen und kleinsten nehmen
      int32 choice2=fromto.back();
      fromto.pop_back();
      int32 choice1=fromto.back();
      fromto.pop_back();
      fromto.pop_back();
      if(choice1<choice2) {
	firstnum=choice1;
	secondnum=choice2;
      } else {
	firstnum=choice2;
	secondnum=choice1;
      }

      locpos+=rightbracket+1;
      //cout << "Rest after bracket: " << (char *) &location[locpos] << endl;
    }else{
      //cout << "Analyse firstnum: " << location << endl;

      firstnum=strtol(&location[locpos],nullptr,10);
      //cout << "Firstnum is " << firstnum << endl;
      while(locpos<location.size() && isdigit(location[locpos])) locpos++;
      //cout << "Eaten until " << location[locpos] << endl;
      //std::string numbers="0123456789";
      //auto numend=location.find_first_not_of(numbers,locpos);
      //std::string firstnum=location.substr(locpos,numend);
      secondnum=firstnum;
    }
    if(locpos<location.size()){
      //cout << "Analyse secondpart: " << (char *) &location[locpos] << endl;
      while(locpos<location.size()
	    && (location[locpos] == '.'
		|| location[locpos]=='^'
		|| location[locpos]=='<'
		|| location[locpos]=='>')) {
	//cout << "Eaten " << location[locpos] << endl;
	locpos++;
      }
      if(location[locpos]=='('){
	//cout << "Analyse bracket val 2nd: " << location[locpos] << endl;
	auto rightbracket=getCorrespondingRightBracket(location, locpos);
	//cout << location.substr(locpos+1,rightbracket-locpos-1) << " rightbrack: " << location[rightbracket] << endl;
	parseGBFLocation(location.substr(locpos+1,rightbracket-locpos-1),
			 fromto,
			 direction);
//      werte poppen und größten nehmen
	int32 choice2=fromto.back();
	fromto.pop_back();
	int32 choice1=fromto.back();
	fromto.pop_back();
	// pop direction
	fromto.pop_back();
	//cout << "popped: " << choice1 << " " << choice2 << "\n";
	if(choice1>choice2) {
	  secondnum=choice1;
	} else {
	  secondnum=choice2;
	}
      }else{
	secondnum=strtol(&location[locpos],nullptr,10);
	while(locpos<location.size() && isdigit(location[locpos])) locpos++;
	//cout << "Eaten until " << location[locpos] << endl;
      }
    }
    //cout << "Secondnum is " << secondnum << endl;
    fromto.push_back(direction);
    fromto.push_back(static_cast<uint32>(firstnum));
    fromto.push_back(static_cast<uint32>(secondnum));
    CEBUG("Push direction: " << direction << endl);
    CEBUG("Push firstnum: " << firstnum << endl);
    CEBUG("Push secondnum: " << secondnum << endl);
    if(locpos<location.size() && location[locpos]==','){
      parseGBFLocation(location.substr(locpos+1,100000000),
		       fromto,
		       direction);
    }
  }
}



/*************************************************************************
 *
 * start must be on first ( bracket
 * gives back position of last character _before_ right bracket
 *
 *************************************************************************/

std::string::size_type GBF::getCorrespondingRightBracket(const std::string & chars, std::string::size_type start) const
{
  FUNCSTART("std::string::size_type GBF::getCorrespondingRightBracket(const std::string & chars, const size_type start) const");

  if(start<chars.size()) {
    if(chars[start]!='('){
      cerr << chars << "At start " << start << " in " << chars << " is " << chars[start] << " and not a bracket?\n";
      MIRANOTIFY(Notify::FATAL, "expected a open bracket (. This also might be an INTERNAL error in the parsing routine!");
    }
  } else {
    MIRANOTIFY(Notify::INTERNAL, ": tried to parse after string???");
  }
  start++;
  int32 leftopen=1;
  while(start<chars.size()){
    //cout << (char) chars[start];
    if(leftopen==1 && chars[start]==')') break;
    if(chars[start]=='(') leftopen++;
    if(chars[start]==')') leftopen--;
    start++;
  }

  //cout << endl;

  //cout << "start " << start << "\tchars.size() " << chars.size() << endl;

  // this line is for "saving" forgotten close brackets at end of lines
  if(start==chars.size() && chars[start-1]==')')  start--;

  FUNCEND();
  return start;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool GBF::checkIfCommentInOldMIRAGBFstyle(const std::string & comment)
{
  if(comment.empty()) return false;

  for(uint32 i=0; GBF_miragbfscankeys[i][0] != 0; ++i) {
    auto fpos=comment.find(GBF_miragbfscankeys[i],0);
    if(fpos!=std::string::npos) return true;
  }
  return false;
}




//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//GBF::GBF(GBF const &other)
//{
//  FUNCSTART("GBF::GBF(GBF const &other)");
//
//  GBF_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//GBF const & GBF::operator=(GBF const & other)
//{
//  FUNCSTART("GBF const & GBF::operator=(GBF const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, GBF const &gbf)
//{
//  FUNCSTART("friend ostream & GBF::operator<<(ostream &ostr, const  &gbf)");
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

void GBF::liftOldMIRAGBFCommentToGFF3(const std::string & src, std::string & dst)
{
  dst.clear();
  if(src.empty()) return;

  std::string tmpstr1;
  std::string tmpstr2;

  auto mI=GBF_mapoldmiragbf2gff3.begin();
  for(; mI!=GBF_mapoldmiragbf2gff3.end(); ++mI){
    if(tag_t::extractGenBankKeyValueFromComment(src,mI->first,tmpstr1)){
      gff3Code(tmpstr1,tmpstr2);
      if(!dst.empty()) dst+=';';
      dst+=mI->second;
      dst+=tmpstr2;
    }
  }

  if(dst.empty()){
    gff3Code(src,tmpstr1);
    dst+="Note=";
    dst+=tmpstr1;
  }
}
