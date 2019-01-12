/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2009 and later by Bastien Chevreux
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
 */


#include "mira/maf_parse.H"

#include <boost/algorithm/string.hpp>


#include "errorhandling/errorhandling.H"

#include "io/annotationmappings.H"
#include "util/fileanddisk.H"


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)



MAFParse::MAFParse(ReadPool * rpool, std::list<Contig> * clist, std::vector<MIRAParameters> * mp)
{
  FUNCSTART("MAFParse::MAFParse(ReadPool * rpool, std::list<Contig> * clist, std::vector<MIRAParameters> * mp)");

  setNewContainers(rpool,clist,mp);
  MAF_piptr= nullptr;

  reset();
}


MAFParse::~MAFParse() {
  //discard();
  if(MAF_piptr!=nullptr) delete MAF_piptr;
  if(MAF_fin.is_open()) MAF_fin.close();
}

void MAFParse::setNewContainers(ReadPool * rpool, std::list<Contig> * clist, std::vector<MIRAParameters> * mp)
{
  FUNCSTART("void MAFParse::setNewContainers(ReadPool * rpool, std::list<Contig> * clist, std::vector<MIRAParameters> * mp)");

  BUGIFTHROW(rpool==nullptr,"rpool==nullptr???");
  BUGIFTHROW(clist!=nullptr && mp==nullptr,"clist!=nullptr && mp==nullptr???");

  MAF_readpool   = rpool;
  MAF_contiglist = clist;
  MAF_miraparams = mp;
}

void MAFParse::reset()
{
  MAF_isinread=false;
  MAF_isincontig=false;
  MAF_isinreadgroup=false;

  MAF_haserror=false;

  cleanupHeaderData();
  cleanupReadData();
  cleanupContigData();
}


void MAFParse::countElements(const std::string & fileName, size_t & numreads, size_t & numcontigs)
{
  FUNCSTART("void MAFParse::countReadsBeforeLoad()");

  numreads=0;
  numcontigs=0;

  std::ifstream mafin(fileName, std::ios::in|std::ios::ate);
  if(!mafin) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << fileName << " not found for loading.");
  }
  if(!mafin.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << fileName << " is empty.");
  }

  ProgressIndicator<std::streamsize> P(0, mafin.tellg(),5000);

  mafin.seekg(0, std::ios::beg);

  std::string actline;
  std::string acttoken;

  actline.reserve(10000);

  while(!mafin.eof()){
    if(mafin.eof()) break;
    getline(mafin,actline);

    if(actline.size()>=2){
      if(actline[0]=='R'
	 && actline[1]=='D') {
	++numreads;
      }else if(actline[0]=='C'
	 && actline[1]=='T') {
	++numcontigs;
      }
    }
    if(P.delaytrigger()) P.progress(mafin.tellg());
  }
  P.finishAtOnce();

  FUNCEND();

  return;
}



void MAFParse::cleanupHeaderData()
{
  MAF_vmajor=-1;
  MAF_vminor=-1;
}

void MAFParse::cleanupContigData()
{
  MAF_contig_assembledfrom.clear();
  MAF_contig_sequence.clear();
  MAF_contig_qualities.clear();
  MAF_contig_taglist.clear();

  MAF_contig_name.clear();

  MAF_contig_numreads=0;
  MAF_contig_len=0;
}

void MAFParse::cleanupReadData()
{
  MAF_read_sequence.clear();
  MAF_read_qualities.clear();
  MAF_read_align_origin.clear();
  MAF_read_taglist.clear();

  MAF_read_name.clear();
  MAF_read_scf_file.clear();
  MAF_read_template.clear();
  MAF_read_base_caller.clear();
  MAF_read_sequencing_vector.clear();
  MAF_read_strain.clear();
  MAF_read_machinetype.clear();

  MAF_read_len=-1;
  MAF_read_insert_size_min=-1;
  MAF_read_insert_size_max=-1;

  MAF_read_ql=-1;
  MAF_read_qr=-1;
  MAF_read_cl=-1;
  MAF_read_cr=-1;
  MAF_read_sl=-1;
  MAF_read_sr=-1;

  MAF_read_strand_given='N';
  MAF_read_tsegment_given=0;
  MAF_read_seqtype=ReadGroupLib::SEQTYPE_SANGER;

  MAF_read_isbackbone=false;
  MAF_read_israil=false;
  MAF_read_isCER=false;

  MAF_read_seenATline=false;

  MAF_readpoolid=-1;
}


//void MAFParse::open(const std::string & fileName)
//{
//}


void MAFParse::registerFile(const std::string & fileName)
{
  FUNCSTART("void MAFParse::registerFile(const std::string & fileName)");

  reset();

  MAF_ccallbackfunc=nullptr;
  MAF_rcallbackfunc=nullptr;
  MAF_filename=fileName;
  MAF_recalcconsensus=false; // TODO!

  MAF_linenumber=0;

  if(MAF_fin.is_open()) MAF_fin.close();

  // Setting a stream buffer of 20MB helps to radically improve read performance
  //  when several concurrent programs read from different files
  //  --> less disk trashing
  // Most useful for miraconvert when reading MAF
  {
    uint32 bsize=20*1024*1024UL;
    MAF_streambuffer.resize(bsize);
    MAF_fin.rdbuf()->pubsetbuf(&(MAF_streambuffer[0]), bsize);
  }

  MAF_fin.open(MAF_filename, std::ios::in);
  if(!MAF_fin) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << MAF_filename << " not found for loading.");
  }
  if(getFileSize(fileName)==0) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << MAF_filename << " is empty.");
  }
}

void MAFParse::setProgressIndicator(bool b)
{
  if(b && MAF_piptr==nullptr){
    uint64 ms=1;
    if(MAF_fin.is_open()) ms=getFileSize(MAF_filename);
    MAF_piptr= new ProgressIndicator<int64>(0, ms,5000);
  }else if(!b && MAF_piptr!=nullptr){
    delete MAF_piptr;
  }
}

void MAFParse::checkCorrectFileEnd()
{
  FUNCSTART("void MAFParse::checkCorrectFileEnd()");

  try{
    if(MAF_isinread){
      MIRANOTIFY(Notify::FATAL, "MAF file ends without closing an open read: " << MAF_read_name << " .... file truncated?");
    }
    if(MAF_isincontig){
      MIRANOTIFY(Notify::FATAL, "MAF file ends without closing an open contig: " << MAF_contig_name << " ... file truncated?");
    }
    if(MAF_isinreadgroup){
      MIRANOTIFY(Notify::FATAL, "MAF file ends without closing an open readgroup. File truncated?");
    }
  }
  catch(Notify n){
    MAF_haserror=true;
    cout << "\nError at end of file " << MAF_filename << endl;
    if(!MAF_contig_name.empty()) cout << "Last contig name read: " << MAF_contig_name << endl;
    if(!MAF_read_name.empty()) cout << "Last read name read: " << MAF_read_name << endl;
    n.handleError(THISFUNC);
  }
}

/*
  seqtype = default seqtype of sequences if not encoded in the MAF
  loadaction:
    //  0 = count only
    //  1 = load
  lrperseqtype = longest read per seqtype

  returns:
    1) number of sequences loaded
    2) when loading: size of longest read per seqtype in lrperseqtype
 */

/*
  BaCh: 15.06.2014
  Uhhh ... seqtype & lrperseqtype not used anymore? isverbose also?
  loadaction also not really useful anymore

  TODO: default seqtype not used??? Make it like for CAF
 */

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
size_t MAFParse::load(const std::string & fileName, const uint8 seqtype, const uint8 loadaction, std::vector<uint32> & lrperseqtype, bool recalcconsensus, void (*ccallback)(std::list<Contig> &, ReadPool &), void (*rcallback)(ReadPool &), bool isVerbose)
{
  FUNCSTART("void MAFParse::load()");

  BUGIFTHROW(loadaction>1,"loadaction>1??");

  if(loadaction==0){
    cout << "Counting reads:\n";
    size_t numreads=0;
    size_t numcontigs=0;
    countElements(MAF_filename,numreads,numcontigs);
    return numreads;
  }

  registerFile(fileName);

  MAF_ccallbackfunc=ccallback;
  MAF_rcallbackfunc=rcallback;
  MAF_recalcconsensus=recalcconsensus;

  cout << "Loading MAF " << MAF_filename << " :\n";

  MAF_lrperseqtype.clear();
  MAF_lrperseqtype.resize(ReadGroupLib::SEQTYPE_END,0);

  if(MAF_piptr==nullptr) MAF_piptr= new ProgressIndicator<int64>(0, 1,5000);
  MAF_piptr->reset(0,getFileSize(MAF_filename));

  MAF_actline.reserve(10000);

  size_t numseqsloaded=loadNextSeqs(-1,-1,-1);

  checkCorrectFileEnd();

  lrperseqtype=MAF_lrperseqtype;

  return numseqsloaded;
}

uint64 MAFParse::loadNextSeqs(uint64 numseqstoload,uint64 numconstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 MAFParse::loadNextSeqs(uint64 numseqsloaded)");

  static const std::string cpsHVersion("@Version");
  static const std::string cpsHProgram("@Program");
  static const std::string cpsHReadGroup("@ReadGroup");
  static const std::string cpsHEndReadGroup("@EndReadGroup");
  static const std::string cpsHRG("@RG");

  static const std::string cpsRS("RS");
  static const std::string cpsRG("RG");
  static const std::string cpsRQ("RQ");
  static const std::string cpsRD("RD");
  static const std::string cpsLR("LR");
  static const std::string cpsSV("SV");
  static const std::string cpsTN("TN");
  static const std::string cpsDI("DI");
  static const std::string cpsTF("TF");
  static const std::string cpsTT("TT");
  static const std::string cpsTS("TS");
  static const std::string cpsSF("SF");
  static const std::string cpsBC("BC");
  static const std::string cpsSL("SL");
  static const std::string cpsSR("SR");
  static const std::string cpsQL("QL");
  static const std::string cpsQR("QR");
  static const std::string cpsCL("CL");
  static const std::string cpsCR("CR");
  static const std::string cpsAO("AO");
  static const std::string cpsRT("RT");
  static const std::string cpsST("ST");
  static const std::string cpsSN("SN");
  static const std::string cpsMT("MT");
  static const std::string cpsIB("IB");
  static const std::string cpsIC("IC");
  static const std::string cpsIR("IR");
  static const std::string cpsER("ER");
  static const std::string cpsCS("CS");
  static const std::string cpsCQ("CQ");
  static const std::string cpsCO("CO");
  static const std::string cpsNR("NR");
  static const std::string cpsLC("LC");
  static const std::string cpsCT("CT");
  static const std::string cpsSLSL("//");
  static const std::string cpsAT("AT");
  static const std::string cpsEC("EC");
  static const std::string cpsBSBS("\\\\");

  uint64 numseqsloaded=0;
  uint64 lenseqsloaded=0;
  uint64 numconsloaded=0;
  bool hascontig=false;
  try {
    while((numconstoload==0 || numconsloaded<numconstoload) && (numseqstoload==0 || numseqsloaded<numseqstoload) && lenseqsloaded<lenseqstoload){
      ++MAF_linenumber;
      MAF_fin >> MAF_acttoken;
      if(MAF_fin.eof()) break;

      CEBUG("l: " << MAF_linenumber << "\tt: ###" << MAF_acttoken << "###" << endl);

      if(MAF_acttoken.empty()) continue;


/* here for read*/

      if(MAF_acttoken==cpsRD){
	// read name
	parseLineRD(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsER){
	// End Read
	parseLineER(MAF_acttoken,MAF_actline);
	++numseqsloaded;
	if(!MAF_isincontig) lenseqsloaded+=MAF_read_sequence.size();
      }else if(MAF_acttoken==cpsRG){
	// Read Group
	parseLineRG(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsRS){
	// Read Sequence
	parseLineRS(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsRQ){
	// Read Qualities
	parseLineRQ(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsLR){
	// length read
	parseLineLR(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsSV){
	// sequencing vector
	parseLineSV(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsTN){
	// template name
	parseLineTN(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsDI){
	// Direction (strand)
	parseLineDI(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsTF){
	// template insize from
	parseLineTF(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsTT){
	// template insize to
	parseLineTT(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsTS){
	// template segment
	parseLineTS(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsSF){
	// Sequencing File
	parseLineSF(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsBC){
	// base caller
	parseLineBC(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsSL){
	//
	parseLineSL(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsSR){
	//
	parseLineSR(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsQL){
	//
	parseLineQL(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsQR){
	//
	parseLineQR(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsCL){
	//
	parseLineCL(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsCR){
	//
	parseLineCR(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsAO){
	//
	parseLineAO(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsRT){
	//
	parseLineRT(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsST){
	//
	parseLineST(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsSN){
	//
	parseLineSN(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsMT){
	//
	parseLineMT(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsIB){
	//
	parseLineIB(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsIC){
	//
	parseLineIC(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsIR){
	//
	parseLineIR(MAF_acttoken,MAF_actline);

/* here for contig*/
      }else if(MAF_acttoken==cpsCO){
	// COntig name
	parseLineCO(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsEC){
	// end contig
	parseLineEC(MAF_acttoken,MAF_actline);
	++numconsloaded;
	lenseqsloaded+=MAF_contig_sequence.size();
      }else if(MAF_acttoken==cpsCS){
	// Consensus Sequence
	parseLineCS(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsCQ){
	// Consensus Qualities
	parseLineCQ(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsNR){
	// Num Reads
	parseLineNR(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsLC){
	// Length Contig
	parseLineLC(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsCT){
	// Contig Tag
	parseLineCT(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsSLSL){
	// start of contig reads

      }else if(MAF_acttoken==cpsBSBS){
	// end of contig reads

      }else if(MAF_acttoken==cpsAT){
	// Assembled From
	parseLineAT(MAF_acttoken,MAF_actline);


/* here for header*/

      }else if(MAF_acttoken==cpsHVersion){
	// file version
	parseLineHeaderVersion(MAF_acttoken,MAF_actline);
      }else if(MAF_acttoken==cpsHProgram){
	// simply read the rest of the line and do nothing
	getline(MAF_fin,MAF_acttoken);
      }else if(MAF_acttoken==cpsHReadGroup){
	// file version
	parseLineHeaderReadGroup(MAF_acttoken,MAF_linenumber);
      }else if(MAF_acttoken==cpsHEndReadGroup
	       || MAF_acttoken==cpsHRG){
	cout << "File " << MAF_filename << ": around line " << MAF_linenumber
	     << "\n" << MAF_acttoken << " occurred without being in an open read group?\n";
	MIRANOTIFY(Notify::FATAL, "Error while reading MAF file.");
      }else{
	cout << "File " << MAF_filename << ": around line " << MAF_linenumber
	     << "\ndid not recognize token " << MAF_acttoken << '\n';
	MIRANOTIFY(Notify::FATAL, "Error while reading MAF file.");
      }

      if(MAF_piptr!=nullptr){
	if(MAF_piptr->delaytrigger()) MAF_piptr->progress(MAF_fin.tellg());
      }
    }
  }
  catch(Notify n){
    MAF_haserror=true;
    cout << "\nError around line " << MAF_linenumber << " of file " << MAF_filename << endl;
    if(!MAF_contig_name.empty()) cout << "Last contig name read: " << MAF_contig_name << endl;
    if(!MAF_read_name.empty()) cout << "Last read name read: " << MAF_read_name << endl;
    n.handleError(THISFUNC);
  }

  if(checkIfEOF() && MAF_piptr!=nullptr) {
    MAF_piptr->finishAtOnce();
    cout << endl;
  }

  FUNCEND();

  return numseqsloaded;
}
//#define CEBUG(bla)


void MAFParse::checkParseIsInReadGroup(std::string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsInReadGroup(std::string & acttoken)");
  if(!MAF_isinreadgroup) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while not in readgroup (@ReadGroup line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsNotInReadGroup(std::string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsNotInReadGroup(std::string & acttoken)");
  if(MAF_isinreadgroup) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while in readgroup (@EndReadGroup line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsInRead(std::string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsInRead(std::string & acttoken)");
  if(!MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while not in read (RD line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsNotInRead(std::string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsNotInRead(std::string & acttoken)");
  if(MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while bein in read (ER line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsInContig(std::string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsInContig(std::string & acttoken)");
  if(!MAF_isincontig) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while not in contig (CO line missing?)");
  }
  if(MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while being in read (RD line not closed by ER?)");
  }
  FUNCEND();
}




void MAFParse::parseLineRD(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineRD(std::string & acttoken, std::string & actline)");
  if(MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered new " << acttoken << " line when the previous read " << MAF_read_name << " was not closed with 'ER'");
  }
  checkParseIsNotInReadGroup(acttoken);

  cleanupReadData();
  MAF_fin >> MAF_read_name;
  MAF_isinread=true;

  FUNCEND();
}

void MAFParse::parseLineRG(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineRG(std::string & acttoken, std::string & actline)");
  checkParseIsInRead(acttoken);

  MAF_fin >> MAF_tmp_str;
  int32 dummy=atoi(MAF_tmp_str.c_str());
  BUGIFTHROW(dummy<0 || dummy >65535,"Line RG: id must be >=0 and <= 65535, but " << dummy << " was given.");
  BUGIFTHROW(dummy>=MAF_readgroup_externalidmapper.size()+1,"Line RG: id of " << dummy << " was given, but not readgroup with this id was defined (@RG ID)");
  BUGIFTHROW(MAF_readgroup_externalidmapper[dummy].isDefaultNonValidReadGroupID(),"Line RG: id of " << dummy << " was given, but not readgroup with this id was defined (@RG ID)");
  MAF_read_rgid=MAF_readgroup_externalidmapper[dummy];

  FUNCEND();
}

void MAFParse::parseLineRS(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineRS(std::string & acttoken, std::string & actline)");

  checkParseIsInRead(acttoken);

  if(!MAF_read_sequence.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered RS line when there already was one for read " << MAF_read_name);
  }

  MAF_fin >> actline;

  MAF_read_sequence.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_read_sequence.push_back(*seq);
  }

  FUNCEND();
}

void MAFParse::parseLineRQ(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineRQ(std::string & acttoken, std::string & actline)");

  checkParseIsInRead(acttoken);

  if(!MAF_read_qualities.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered RQ line when there already was one for read " << MAF_read_name);
  }

  MAF_fin >> actline;

  MAF_read_qualities.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_read_qualities.push_back(*seq-33);
  }

  FUNCEND();
}

void MAFParse::parseLineLR(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_len;
  if(MAF_read_len>actline.capacity()) actline.reserve(MAF_read_len+10);
}

void MAFParse::parseLineSV(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_sequencing_vector;
}

void MAFParse::parseLineTN(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_template;
}

void MAFParse::parseLineDI(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_strand_given;
}

void MAFParse::parseLineTF(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_insert_size_min;
}

void MAFParse::parseLineTT(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_insert_size_max;
}

void MAFParse::parseLineTS(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  uint32 dummy;     // sigh ... detour to read a number and NOT a 'char';
  MAF_fin >> dummy;
  MAF_read_tsegment_given=static_cast<uint8>(dummy);
}

void MAFParse::parseLineSF(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_scf_file;
}

void MAFParse::parseLineBC(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_base_caller;
}

void MAFParse::parseLineSL(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_sl;
  MAF_read_sl--;
}

void MAFParse::parseLineSR(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_sr;
}

void MAFParse::parseLineQL(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_ql;
  MAF_read_ql--;
}

void MAFParse::parseLineQR(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_qr;
}

void MAFParse::parseLineCL(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_cl;
  MAF_read_cl--;
}

void MAFParse::parseLineCR(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_cr;
}

void MAFParse::parseLineAO(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineAO(std::string & acttoken, std::string & actline)");

  checkParseIsInRead(acttoken);

  if(MAF_read_sequence.empty()){
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ": sequence (SQ line) must be defined before AO line");
  }
  if(MAF_read_align_origin.empty()){
    MAF_read_align_origin.resize(MAF_read_sequence.size(),-1);
  }

  int32 seqfrom, seqto, origfrom, origto;
  MAF_fin >> seqfrom;
  MAF_fin >> seqto;
  MAF_fin >> origfrom;
  MAF_fin >> origto;

  //cout << "xxx " << seqfrom << " " << seqto << " " << origfrom << " " << origto << "\n";

  if(seqfrom<1
     || seqto<1
     || origfrom<1
     || origto<1){
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ":  values may not be <1");
  }

  int32 seqinc=1;
  int32 originc=1;
  if(seqto<seqfrom) seqinc=-1;
  if(origto<origfrom) originc=-1;


  if (abs(seqto - seqfrom) != abs(origto - origfrom)) {
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ":  distance between seqfrom/to (" << seqfrom << " " << seqto << ") is unequal to originalfrom/to (" << origfrom << " " << origto << ")");
  }
  if (std::max(seqfrom, seqto) > static_cast<int32>(MAF_read_align_origin.size())) {
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ":  seqfrom/to (" << seqfrom << " " << seqto << ") is larger than size of read (" << MAF_read_align_origin.size() << ")");
  }

  int32 seqi=seqfrom;
  int32 origi=origfrom;
  for(int32 loopi=0; loopi < abs(seqto-seqfrom)+1; loopi++){
    if(seqi-1 <0 || seqi-1 >= MAF_read_align_origin.size()) {
      MIRANOTIFY(Notify::FATAL,"While reading AO line for read: " << MAF_read_name << " with AO values " << seqfrom << " " << seqto << " " <<origfrom << " " << origto << "\nThis makes for an illegal alignment of the sequence to the original, wrong values in this line?\n");
    }

    MAF_read_align_origin[seqi-1] = origi - 1;
    seqi += seqinc;
    origi += originc;
  }

  FUNCEND();
}

void MAFParse::parseLineRT(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);

  MAF_read_taglist.resize(MAF_read_taglist.size()+1);
  if(MAF_vmajor<2){
    parseTagData(acttoken,MAF_read_taglist.back());
  }else{
    parseTagDataV2(acttoken,MAF_read_taglist.back());
  }

  //cout << "Stored Rtag: " << MAF_read_taglist.back();
  //MAF_read_taglist.back().dumpDebug(cout);
}

void MAFParse::parseTagData(std::string & acttoken, multitag_t & targettag)
{
  FUNCSTART("void MAFParse::parseTagData(std::string & acttoken, multitag_t & tag)");

  multitag_t tmptag;

  MAF_fin >> MAF_tmp_str;

  if(!AnnotationMappings::isValidGFF3SOEntry(MAF_tmp_str)){
    std::string soident(AnnotationMappings::translateGAP4feat2SOfeat(MAF_tmp_str));
    if(soident.empty()){
      soident=AnnotationMappings::translateXGAP4feat2SOfeat(MAF_tmp_str);
      if(soident.empty()){
	tmptag.setIdentifierStr(MAF_tmp_str);
      }
    }
    if(!soident.empty()){
      tmptag.setIdentifierStr(soident);
    }
  }else{
    tmptag.setIdentifierStr(MAF_tmp_str);
  }

  MAF_fin >> tmptag.from;
  MAF_fin >> tmptag.to;

  if(tmptag.from<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in tmptag line " << acttoken << ": (" << tmptag.from << " " << tmptag.to << ") -> " << tmptag.from << " is <1, not allowed.");
  }
  if(tmptag.to<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in tmptag line " << acttoken << ": (" << tmptag.from << " " << tmptag.to << ") -> " << tmptag.to << " is <1, not allowed.");
  }

  tmptag.from-=1;
  tmptag.to-=1;


  if(tmptag.from<=tmptag.to){
    tmptag.setStrand('+');
  }else{
    tmptag.setStrand('-');
    std::swap(tmptag.from, tmptag.to);
  }

  // comment may be present or not
  char nextchar;
  MAF_fin.get(nextchar);
  if(nextchar=='\n') {
    targettag=tmptag;
    return;
  }
  if(nextchar=='\r') {
    // also eat \n
    MAF_fin.get(nextchar);
    targettag=tmptag;
    return;
  }
  getline(MAF_fin,MAF_tmp_str);
  tmptag.setCommentStr(MAF_tmp_str);

//  cout << "BEFORE: ";
//  tmptag.dumpDebug(cout);
  Read::upgradeOldTagToMultitagWithGFF3(tmptag,targettag);
//  cout << "AFTER: ";
//  targettag.dumpDebug(cout);

  if(targettag.getSourceStr().empty()){
    if(AnnotationMappings::translateSOfeat2GAP4feat(targettag.getIdentifierStr()).empty()
       && AnnotationMappings::translateSOfeat2XGAP4feat(targettag.getIdentifierStr()).empty()){
      targettag.source = multitag_t::MT_tagsrcentry_idMIRA;
    }else{
      targettag.source = multitag_t::MT_tagsrcentry_idGFF3;
    }
    targettag.source = multitag_t::MT_tagsrcentry_idGFF3;
  }else if(targettag.source == multitag_t::MT_tagsrcentry_idMIRA){
    if(!AnnotationMappings::translateSOfeat2SOID(targettag.getIdentifierStr()).empty()){
      targettag.source = multitag_t::MT_tagsrcentry_idGFF3;
    }
  }

  FUNCEND();
}



void MAFParse::parseTagDataV2(std::string & acttoken, multitag_t & targettag)
{
  // Raw speed to parse tab delimited line
  // both boost::split and boost::tokenizer are 40% to 50% slower.
  //
  // Cheating big time: the line std::string gets rewritten, replacing tabs by 0 and
  //  writing substring pointer (char *) to a fixed array
  //
  // When throwing with MAF_tmp_str in the message, we need to re-tabify beforehand!

  FUNCSTART("void MAFParse::parseTagData(std::string & acttoken, multitag_t & tag)");

  targettag.commentisgff3=true;

  {
    char nextchar;
    MAF_fin.get(nextchar); // eat away the following \t
  }
  getline(MAF_fin,MAF_tmp_str);

  static char * sarr[8];

  //////////// Overwrite \t with 0, populate the array pf char * to subparts of line
  char * sptr=const_cast<char *>(MAF_tmp_str.c_str());
  char * rptr=sptr;

  uint32 numtabs=0;
  while(*rptr){
    if(*rptr=='\t'){
      *rptr=0;
      if(numtabs>8){
	for(auto & x : MAF_tmp_str) if(x==0) x='\t';
	MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found more");
      }
      sarr[numtabs]=sptr;
      ++numtabs;
      sptr=rptr+1;
    }
    ++rptr;
  }
  if(rptr!=MAF_tmp_str.c_str()){
    sarr[numtabs]=sptr;
    ++numtabs;
  }

  if(numtabs<4 || numtabs>7){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found " << numtabs);
  }

  for(;numtabs<8;++numtabs) sarr[numtabs]=nullptr;
  ///// done

  static std::string acts;
  acts=sarr[0];

  if(AnnotationMappings::isMIRAEntry(acts)){
    targettag.setIdentifierStr(acts);
  }else{
    const std::string * soident=&AnnotationMappings::translateGAP4feat2SOfeat(acts);
    if(soident->empty()){
      soident=&AnnotationMappings::translateXGAP4feat2SOfeat(acts);
    }
    if(soident->empty()){
      targettag.setIdentifierStr(acts);
    }else{
      targettag.setIdentifierStr(*soident);
    }
  }

  targettag.from=atoi(sarr[1]);
  targettag.to=atoi(sarr[2]);

  if(targettag.from<1){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.from << " is <1, not allowed.");
  }
  if(targettag.to<1){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.to << " is <1, not allowed.");
  }

  targettag.from-=1;
  targettag.to-=1;

  if(*sarr[3]==0){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": the entry for the strand is empty, not allowed.");
  }

  targettag.setStrand(*sarr[3]);

  if(sarr[4]!=nullptr){
    targettag.setSourceStr(sarr[4]);

    if(sarr[5]!=nullptr){
      // hmmmmmm .... something was very wrong with the previous implementation
      switch(*sarr[5]){
      case '.':
      case '1':{
	targettag.phase=0;
	break;
      }
      case '2':{
	targettag.phase=1;
	break;
      }
      case '3':{
	targettag.phase=2;
	break;
      }
      default : {
	for(auto & x : MAF_tmp_str) if(x==0) x='\t';
	MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << "\t" << MAF_tmp_str << "\nthe entry for the strand (" << *sarr[5] << ") is not 0, 1, 2, 3 or .");
      }
      }

      if(sarr[6]!=nullptr){
	targettag.setCommentStr(sarr[6]);
      }
    }
  }

  FUNCEND();
}


void MAFParse::parseLineST(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineST(std::string & acttoken, std::string & actline)");
  checkParseIsInRead(acttoken);
  MAF_fin >> actline;

  MAF_read_seqtype=ReadGroupLib::stringToSeqType(actline);

  if(MAF_read_seqtype==ReadGroupLib::SEQTYPE_END){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in tag " << acttoken << ": unknown sequencing type '" << actline << "'?");
  }
  FUNCEND();
}

void MAFParse::parseLineSN(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_strain;
}

void MAFParse::parseLineMT(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_machinetype;
}

void MAFParse::parseLineIB(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_isbackbone;
}

void MAFParse::parseLineIC(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_isCER;
}

void MAFParse::parseLineIR(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  MAF_fin >> MAF_read_israil;
}

void MAFParse::parseLineER(std::string & acttoken, std::string & actline)
{
  checkParseIsInRead(acttoken);
  addReadToReadPool();
  MAF_isinread=false;
  if(!MAF_isincontig && MAF_rcallbackfunc!=nullptr) {
    (*MAF_rcallbackfunc)(*MAF_readpool);
  }
  MAF_read_rgid.resetLibId();
}


void MAFParse::parseLineCO(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineCO(std::string & acttoken, std::string & actline)");

  if(MAF_isincontig){
    MIRANOTIFY(Notify::FATAL, "Seen new CO line while previous CO was not closed by EC");
  }
  checkParseIsNotInReadGroup(acttoken);

  cleanupContigData();
  MAF_fin >> MAF_contig_name;
  MAF_isincontig=true;
  FUNCEND();
}

void MAFParse::parseLineCS(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineSQ(std::string & acttoken, std::string & actline)");

  checkParseIsInContig(acttoken);

  if(!MAF_contig_sequence.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered CS line when there already was one for contig " << MAF_contig_name);
  }

  MAF_fin >> actline;

  MAF_contig_sequence.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_contig_sequence.push_back(*seq);
  }

  FUNCEND();
}

void MAFParse::parseLineCQ(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineCQ(std::string & acttoken, std::string & actline)");

  checkParseIsInContig(acttoken);

  if(!MAF_contig_qualities.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered CQ line when there already was one for contig " << MAF_contig_name);
  }

  MAF_fin >> actline;

  MAF_contig_qualities.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_contig_qualities.push_back(*seq-33);
  }

  FUNCEND();
}

void MAFParse::parseLineNR(std::string & acttoken, std::string & actline)
{
  checkParseIsInContig(acttoken);
  MAF_fin >> MAF_contig_numreads;
}

void MAFParse::parseLineLC(std::string & acttoken, std::string & actline)
{
  checkParseIsInContig(acttoken);
  MAF_fin >> MAF_contig_len;

  if(MAF_contig_len>actline.capacity()) actline.reserve(MAF_contig_len+10);
}

void MAFParse::parseLineCT(std::string & acttoken, std::string & actline)
{
  checkParseIsInContig(acttoken);

  MAF_contig_taglist.resize(MAF_contig_taglist.size()+1);
  if(MAF_vmajor<2){
    parseTagData(acttoken,MAF_contig_taglist.back());
    // for consensus tags, change strand to '='
    // BaCh 21.01.2012: Why?
    //MAF_contig_taglist.back().setStrand('=');
  }else{
    parseTagDataV2(acttoken,MAF_contig_taglist.back());
  }

  //cout << "Stored Ctag: ";
  //MAF_contig_taglist.back().dumpDebug(cout);
}

void MAFParse::parseLineAT(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineAT(std::string & acttoken, std::string & actline)");

  checkParseIsInContig(acttoken);
  checkParseIsNotInRead(acttoken);

  if(MAF_readpoolid<0){
    MIRANOTIFY(Notify::FATAL, "Seen AT line but no read in contig defined before? (RD/ER block in a CO block)");
  }
  if(MAF_read_seenATline){
    MIRANOTIFY(Notify::FATAL, "Seen AT line, but either no read before or multiple AT lines.");
  }
  MAF_read_seenATline=true;

  int32 cfrom,cto,rfrom,rto;
  Contig::contig_init_read_t tmpcr;
  int8 direction=1;

  MAF_fin >> cfrom;
  MAF_fin >> cto;
  MAF_fin >> rfrom;
  MAF_fin >> rto;

  if(cfrom > cto){
    direction=-1;
    tmpcr.offset_in_contig=cto-1;
  }else{
    tmpcr.offset_in_contig=cfrom-1;
  }

  if(rfrom>rto){
    tmpcr.read_rclip=rfrom;
    tmpcr.read_lclip=rto-1;
    tmpcr.direction= -direction;
  }else{
    tmpcr.read_rclip=rto;
    tmpcr.read_lclip=rfrom-1;
    tmpcr.direction= direction;
  }

  tmpcr.id=MAF_readpoolid;
  MAF_contig_assembledfrom.push_back(tmpcr);

  FUNCEND();
}

void MAFParse::parseLineEC(std::string & acttoken, std::string & actline)
{
  FUNCSTART("void MAFParse::parseLineEC(std::string & acttoken, std::string & actline)");

  checkParseIsInContig(acttoken);
  MAF_isincontig=false;

  if(MAF_contiglist!=nullptr){
    {
      Contig dummy(MAF_miraparams, *MAF_readpool);
      MAF_contiglist->push_back(dummy);
    }

    try{
      if(MAF_recalcconsensus){
	std::string dummy1;
	std::vector<base_quality_t> dummy2;
	MAF_contiglist->back().initialiseContig(MAF_contig_assembledfrom,
						MAF_contig_taglist,
						MAF_contig_name,
						dummy1,dummy2);
      }else{
	std::string dummy1;
	dummy1.reserve(MAF_contig_sequence.size()+2);
	for(auto dnaI=MAF_contig_sequence.cbegin(); dnaI != MAF_contig_sequence.end(); dnaI++) dummy1+=*dnaI;
	MAF_contiglist->back().initialiseContig(MAF_contig_assembledfrom,
						MAF_contig_taglist,
						MAF_contig_name,
						dummy1,
						MAF_contig_qualities);
      }
    }
    catch(Notify n){
      cout << "Error for contig " << MAF_contig_name << endl;
      n.handleError(THISFUNC);
    }

    try{
      if(MAF_ccallbackfunc!=nullptr) {
	(*MAF_ccallbackfunc)(*MAF_contiglist, *MAF_readpool);
      }
    }
    catch(Notify n){
      cout << "Error while calling callback!\n";
      cout << "Error for contig " << MAF_contig_name << endl;
      n.handleError(THISFUNC);
    }
  }

  FUNCEND();
}



void MAFParse::checkReadData()
{
  FUNCSTART("void MAFParse::checkReadData()");

  BUGIFTHROW(MAF_read_seqtype>=ReadGroupLib::SEQTYPE_END, "Illegal seqtype in checkReadData()???");

  if(MAF_read_len>=0
     && MAF_read_sequence.size() != MAF_read_len){
    MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << ": size of sequence (" << MAF_read_sequence.size() << ") is not equal to size given in LR line (" << MAF_read_len << ")");
  }
  if(!MAF_read_qualities.empty()){
    if(MAF_read_sequence.size() != MAF_read_qualities.size()){
      MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << ": size of sequence (" << MAF_read_sequence.size() << ") is not equal to size of qualities (" << MAF_read_qualities.size() << ")");
    }
  }else{
    // sequence but no qualities ... then fake some according to the sequencing type
    // give them the standard qual for this sequencing type
    MAF_read_qualities.resize(MAF_read_sequence.size(),MAF_read_rgid.getDefaultQual());
  }

  if(!MAF_read_align_origin.empty() && MAF_read_align_origin.size() != MAF_read_sequence.size()){
    MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << ": the align to origin (AO) data led to a larger or smaller array that the length of the sequence?");
  }

  if(MAF_read_ql < 0) MAF_read_ql = 0;
  if(MAF_read_sl < 0) MAF_read_sl = 0;
  if(MAF_read_cl < 0) MAF_read_cl = 0;

  if(MAF_read_qr < 0) MAF_read_qr = MAF_read_sequence.size();
  if(MAF_read_sr < 0) MAF_read_sr = MAF_read_sequence.size();
  if(MAF_read_cr < 0) MAF_read_cr = MAF_read_sequence.size();

  FUNCEND();
}


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void MAFParse::addReadToReadPool()
{
  FUNCSTART("void MAFParse::addReadToReadPool()");

  checkReadData();

  MAF_readpoolid=MAF_readpool->size();
  Read & newread = MAF_readpool->getRead(MAF_readpool->provideEmptyRead());

  ReadGroupLib::ReadGroupID rgid;

  if(MAF_vmajor>=2 && MAF_read_rgid.isDefaultNonValidReadGroupID()){
    MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << " has no RG line to define the read group. This is needed for MAF version 2");
  }
  if(MAF_read_rgid.isDefaultNonValidReadGroupID()){
    std::string rgname;
    rgid=ReadGroupLib::searchExactRGMatch(
      rgname,
      MAF_read_seqtype,
      MAF_read_insert_size_min,
      MAF_read_insert_size_max,
      ReadGroupLib::SPLACE_UNKNOWN,
      MAF_read_strain,
      MAF_read_isbackbone,
      MAF_read_israil,
      MAF_read_isCER,
      MAF_read_sequencing_vector,
      MAF_read_machinetype,
      MAF_read_base_caller);

    if(rgid.isDefaultNonValidReadGroupID()){
      rgid=ReadGroupLib::newReadGroup();
      rgid.setGroupName(rgname);
      rgid.setSequencingType(MAF_read_seqtype);
      rgid.setInsizeFrom(MAF_read_insert_size_min);
      rgid.setInsizeTo(MAF_read_insert_size_max);
      rgid.setSegmentPlacement("?");
      rgid.setStrainName(MAF_read_strain);
      rgid.setBackbone(MAF_read_isbackbone);
      rgid.setRail(MAF_read_israil);
      rgid.setCoverageEquivalentRead(MAF_read_isCER);
      rgid.setSeqVecName(MAF_read_sequencing_vector);
      rgid.setMachineType(MAF_read_machinetype);
      rgid.setBaseCaller(MAF_read_base_caller);
    }
  }else{
    rgid=MAF_read_rgid;
  }

  //cout << "Setting rgid to:\n" << rgid << endl;

  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA
     || rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){
    newread.disallowAdjustments();
  }

  if(newread.usesAdjustments() && MAF_read_align_origin.empty() && !MAF_read_sequence.empty()){
    // no AO line given, i.e., no insertions/deletions
    // create vector which represents that
    MAF_read_align_origin.resize(MAF_read_sequence.size());
    uint32 num=0;
    for(auto & x : MAF_read_align_origin){
      x=num++;
    }
  }

  newread.initialiseRead(false,
			 false,
			 true,     // always padded
			 rgid,
			 MAF_read_sequence,
			 MAF_read_qualities,
			 MAF_read_align_origin,
			 MAF_read_taglist,
			 MAF_read_name,
			 MAF_read_scf_file,
			 MAF_read_ql,
			 MAF_read_qr,
			 MAF_read_sl,
			 MAF_read_sr,
			 MAF_read_cl,
			 MAF_read_cr);


  if(MAF_read_tsegment_given!=0){
    newread.setTemplateSegment(MAF_read_tsegment_given);
  }else if(MAF_read_strand_given!='N'){
    if(MAF_read_strand_given=='F'){
      newread.setTemplateSegment(1);
    }else{
      newread.setTemplateSegment(255);
    }
  }
  if (!MAF_read_template.empty()) {
    newread.setTemplate(MAF_read_template);
  }

  CEBUG("ER read is:\n" << newread << endl);
}
#define CEBUG(bla)



void MAFParse::parseLineHeaderVersion(std::string & acttoken, std::string & actline)
{
  MAF_fin >> MAF_tmp_str;
  MAF_vmajor=atoi(MAF_tmp_str.c_str());
  MAF_fin >> MAF_tmp_str;
  MAF_vminor=atoi(MAF_tmp_str.c_str());
}

void MAFParse::parseLineHeaderReadGroup(std::string & acttoken, uint64 & linenumber)
{
  FUNCSTART("void MAFParse::parseLineHeaderReadGroup(std::string & acttoken, uint64 & linenumber)");

  checkParseIsNotInReadGroup(acttoken);

  MAF_readgroup_rgid=ReadGroupLib::newReadGroup();
  parseReadGroup(MAF_fin,MAF_readgroup_rgid,MAF_readgroup_externalidmapper,linenumber);
  MAF_readgroup_rgid.fillInSensibleDefaults();
  MAF_readgroup_rgid.resetLibId();

  FUNCEND();
}


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void MAFParse::parseReadGroup(std::ifstream & mafin, ReadGroupLib::ReadGroupID & rgid, std::vector<ReadGroupLib::ReadGroupID> & externalidmapper, uint64 & linenumber)
{
  FUNCSTART("void MAFParse::parseReadGroup(std::ifstream & mafin, ReadGroupLib::ReadGroupID & rgid, std::vector<ReadGroupLib::ReadGroupID> & externalidmapper, uint64 & linenumber)");

  std::string mafline;
  std::vector<std::string> mafsplit;

  while(true){
    ++linenumber;
    getline(mafin,mafline);
    CEBUG("MLRG: " << mafline << endl);
    if(mafin.eof()) break;
    if(mafline.empty()) continue;
    boost::split(mafsplit, mafline, boost::is_any_of("\t"));
    if(mafsplit.empty()) continue;
    if(mafsplit.size()==1){
      if(mafsplit[0]=="@EndReadGroup") break;
      cout << "\nOuch, erroneous line: " << mafline << endl;
      MIRANOTIFY(Notify::FATAL,"Did not find a tab character in line " << linenumber << " and keyword is not @EndReadGroup? Something is wrong.");
    }

    std::string & rgtoken=mafsplit[1];
    CEBUG("read rgtoken #" << rgtoken << "#\n");

    if(rgtoken=="isbackbone"){
      rgid.setBackbone(true);
      continue;
    }else if(rgtoken=="israil"){
      rgid.setRail(true);
      continue;
    }else if(rgtoken=="iscoverageequivalent"){
      rgid.setCoverageEquivalentRead(true);
      continue;
    }

    if(mafsplit.size()<3){
      MIRANOTIFY(Notify::FATAL,"Line " << mafline << "\nexpected at least 3 elements, found " << mafsplit.size() << endl);
    }

    std::string & rgval1=mafsplit[2];
    CEBUG("rgval1 #" << rgval1 << "#\n");

    if(rgtoken=="name"){
      rgid.setGroupName(rgval1);
    }else if(rgtoken=="segmentnaming"
	     || rgtoken=="templatenaming"){
      //TODO:  implement templatenaming
    }else if(rgtoken=="ID"){
      int32 dummy=atoi(rgval1.c_str());
      if(dummy<0 || dummy >65535){
	MIRANOTIFY(Notify::FATAL,"Line @RG ID: id must be >=0 and <= 65535, but " << dummy << " was given.");
      }
      if(dummy>=externalidmapper.size()){
	externalidmapper.resize(dummy+1);
      }
      externalidmapper[dummy]=rgid;
    }else if(rgtoken=="technology"){
      rgid.setSequencingType(rgval1);
    }else if(rgtoken=="strainname"){
      rgid.setStrainName(rgval1);
    }else if(rgtoken=="segmentplacement"
	     || rgtoken=="templateplacement"){
      if(!rgid.setSegmentPlacement(rgval1)){
	MIRANOTIFY(Notify::FATAL,"Line @RG segmentplacement: did not recognise '" << rgval1 << "' as valid placement code.");
      }
    }else if(rgtoken=="templatesize"){
      int32 dummy=atoi(rgval1.c_str());
      rgid.setInsizeFrom(dummy);
      dummy=atoi(mafsplit[3].c_str());
      rgid.setInsizeTo(dummy);
    }else if(rgtoken=="machinetype"){
      rgid.setMachineType(rgval1);
    }else if(rgtoken=="basecaller"){
      rgid.setBaseCaller(rgval1);
    }else if(rgtoken=="dye"){
      rgid.setDye(rgval1);
    }else if(rgtoken=="primer"){
      rgid.setPrimer(rgval1);
    }else if(rgtoken=="clonevecname"){
      rgid.setCloneVecName(rgval1);
    }else if(rgtoken=="seqvecname"){
      rgid.setSeqVecName(rgval1);
    }else if(rgtoken=="adaptorleft"){
//    rgid.set(rgval1);
    }else if(rgtoken=="adaptorright"){
//    rgid.setSeqVecName(rgval1);
    }else if(rgtoken=="adaptorsplit"){
//    rgid.setSeqVecName(rgval1);
    }else if(rgtoken=="datadir"){
      rgid.setDataDir(rgval1);
    }else if(rgtoken=="datafile"){
      rgid.setDataFile(rgval1);
    }else{
      MIRANOTIFY(Notify::FATAL,"For line @RG: did not recognize token " << rgtoken);
    }
  }

  FUNCEND();
}

