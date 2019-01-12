/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 *   and Thomas Pfisterer
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


#include "mira/readpool.H"

#include "errorhandling/errorhandling.H"

#include "io/generalio.H"
#include "io/fasta.H"
#include "io/fastq-lh.H"
#include "io/phd.H"
#include "io/ncbiinfoxml.H"

#include "util/fileanddisk.H"
#include "util/progressindic.H"

#include "caf/caf.H"

#include "mira/gbf_parse.H"
#include "mira/gff_parse.H"
#include "mira/maf_parse.H"

#include <unordered_map>
#include <unordered_set>

KSEQ_INIT(gzFile, myGZRead)


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)

std::string ReadPool::RP_missingfastaqual_resolvemsg("name your files with a .fna postfix");

const base_quality_t ReadPool::RP_sxa2phredmap[256]= {
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 2, 2, 3,
  3, 4, 4, 5, 5, 6, 7, 8,
  9, 10, 10, 11, 12, 13, 14, 15,
  16, 17, 18, 19, 20, 21, 22, 23,
  24, 25, 26, 27, 28, 29, 30, 31,
  32, 33, 34, 35, 36, 37, 38, 39,
  40, 41, 42, 43, 44, 45, 46, 47,
  48, 49, 50, 51, 52, 53, 54, 55,
  56, 57, 58, 59, 60, 61, 62, 63,
  64, 65, 66, 67, 68, 69, 70, 71,
  72, 73, 74, 75, 76, 77, 78, 79,
  80, 81, 82, 83, 84, 85, 86, 87,
  88, 89, 90, 91, 92, 93, 94, 95,
  96, 97, 98, 99, 100, 101, 102, 103,
  104, 105, 106, 107, 108, 109, 110, 111,
  112, 113, 114, 115, 116, 117, 118, 119,
  120, 121, 122, 123, 124, 125, 126, 127,
  128, 129, 130, 131, 132, 133, 134, 135,
  136, 137, 138, 139, 140, 141, 142, 143,
  144, 145, 146, 147, 148, 149, 150, 151,
  152, 153, 154, 155, 156, 157, 158, 159,
  160, 161, 162, 163, 164, 165, 166, 167,
  168, 169, 170, 171, 172, 173, 174, 175,
  176, 177, 178, 179, 180, 181, 182, 183,
  184, 185, 186, 187, 188, 189, 190, 191
};



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
ReadPool::ReadContainer const & ReadPool::ReadContainer::operator=(ReadPool::ReadContainer const & other)
{
  // 1 to 1 copy, including eventually free elements so that
  // operator[] has the same result on source and copy object
  if(this != &other){
    RC_thepool.clear();
    for(auto rptr : other.RC_poolrptr){
      RC_thepool.push_back(*rptr);
      RC_poolrptr.push_back(&(RC_thepool.back()));
    }
    RC_releasedidx=other.RC_releasedidx;
  }
  //cout << "Copied ReadContainer size " << size() << " and with " << getNumActiveReads() << " active reads\n";
  return *this;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
ReadPool::ReadPool()
{
  FUNCSTART("ReadPool::ReadPool()");

  REP_thepool3.clear();

  REP_allownameindex=false;
  REP_nameindex.clear();

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void ReadPool::discard()
{
  FUNCSTART("ReadPool::discard()");

  REP_thepool3.clear();
  REP_nameindex.clear();

  //delete [] REP_filenames;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
ReadPool::~ReadPool()
{
  FUNCSTART("ReadPool::~ReadPool()");

  discard();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

size_t ReadPool::estimateMemoryUsage() const
{
  FUNCSTART("size_t ReadPool::estimateMemoryUsage()");

  size_t ret=sizeof(ReadPool);
  // TODO: estimateMemoryUsage for deque
  //ret+=estimateMemoryUsageOfContainer(REP_thepool2,false);

  FUNCEND();
  return ret;
}


/*************************************************************************
 *
 * makes template IDs for the reads in the pool
 *
 * returns
 *  true if there are usable templates (less templates than valid reads)
 *  false otherwise
 *
 *************************************************************************/

bool ReadPool::makeTemplateIDs(uint8 nagwarn_templateproblems, bool verbose)
{
  FUNCSTART("void ReadPool::makeTemplateIDs()");

// TODO: bla auf backbone und rails eingehen

  if(verbose){
    rpDateStamp();
    cout << endl;
  }

  std::vector<int32> tidcounter;
  std::vector<int32> tid_firstpartner;
  tidcounter.resize(size(), 0);
  tid_firstpartner.resize(size(), -1);

  // will we need to check template ends?
  // rationale: if no info about template ends was available in the
  //  ancillary data, the no_te_check will ensure that we'll then still use
  //  the template information
  bool notecheck=true;
  for(size_t i=0; i<REP_thepool3.size(); ++i){
    if(getRead(i).getTemplateSegment()!=0) {
      notecheck=false;
    }
  }


  typedef std::unordered_map<std::string, int32> strintmap;
  strintmap tnmap;
  decltype(tnmap.end()) tnI;

  int32 acttid=0;
  int32 validreads=0;

  uint32 outlines_fatal=0;
  uint32 outlines_warn=0;
  uint32 ol_emptytn=0;
  uint32 ol_gt2reads=0;
  uint32 ol_unknownsegment=0;
  uint32 ol_unknownsegmentfixed=0;
  uint32 ol_samesegment=0;
  for(uint32 rid=0; rid<size();++rid){
      //cout << "acttid: " << acttid << "\t";
    if(getRead(rid).hasValidData()==false) continue;
    ++validreads;
    if(getRead(rid).getTemplate().empty()) {
      getRead(rid).setTemplateID(acttid);
      tidcounter[acttid]++;
      ++acttid;
      if(getRead(rid).getReadGroupID().hasTemplateInfo()){
	// Ooooops ... empty template but readgroupinfo?
	if(ol_emptytn<200){
	  ++outlines_fatal;
	  cout << "Read " << getRead(rid).getName() << " has template info given, but the internal template name is empty? This is fishy!" << endl;
	  if(++ol_emptytn==200){
	    cout << "More than 200 cases like the above, will not report more.\n";
	  }
	}
      }
    }else{
      tnI=tnmap.find(getRead(rid).getTemplate());
      if(tnI!=tnmap.end()){
	getRead(rid).setTemplateID(tnI->second);
	int32 firstpartner=tid_firstpartner[tnI->second];
	if(++tidcounter[tnI->second]==2){
	  getRead(rid).setTemplatePartnerID(firstpartner);
	  getRead(firstpartner).setTemplatePartnerID(rid);

	  // look out for unknown segments ("0"), try to fix if possible, else dump error message
	  uint8 numunknownseg=(getRead(rid).getTemplateSegment()==0)+(getRead(firstpartner).getTemplateSegment()==0);
	  if(numunknownseg){
	    bool couldfix=false;
	    if(numunknownseg==1){
	      Read * tofix=&getRead(rid);
	      Read * other=&getRead(firstpartner);
	      uint8 segmenttoset=0;
	      if(tofix->getTemplateSegment()>0) std::swap(tofix,other);
	      if(other->getTemplateSegment()==1){
		segmenttoset=255;
	      }else if(other->getTemplateSegment()==255){
		segmenttoset=1;
	      }
	      if(segmenttoset){
		couldfix=true;
		tofix->setTemplateSegment(segmenttoset);
		ol_unknownsegmentfixed++;
		if(ol_unknownsegmentfixed<200){
		  cout << "Read " << tofix->getName() << ": fixed unrecognised template segment.\n";
		}
		if(++ol_unknownsegmentfixed==200){
		  cout << "More than 200 cases like the above, will not report more.\n";
		}
	      }
	    }
	    if(!couldfix){
	      if(ol_unknownsegment<200){
		++outlines_fatal;
		if(getRead(rid).getTemplateSegment()==0){
		  cout << "DNA template " << tnI->first << ", read " << getRead(rid).getName() << ": template segment not recognised.\n";
		}
		if(getRead(firstpartner).getTemplateSegment()==0){
		  cout << "DNA template " << tnI->first << ", read " << getRead(firstpartner).getName() << ": template segment not recognised.\n";
		}
		if(++ol_unknownsegment==200){
		  cout << "More than 200 cases like the above, will not report more.\n";
		}
	      }
	    }
	  }

	  // the folowing is almost impossible when reading FASTA/FASTQ and other simple files,
	  //  but it may happen in CAF/MAF (and maybe SAM?)
	  if(getRead(rid).getTemplateSegment()==getRead(firstpartner).getTemplateSegment()){
	    if(ol_samesegment<200){
	      ++outlines_fatal;
	      cout << "Reads " << getRead(rid).getName() << " and " << getRead(firstpartner).getName() << " have the same template segments.\n";
	      if(++ol_samesegment==200){
		cout << "More than 200 cases like the above, will not report more.\n";
	      }
	    }
	  }

	}else{
	  // Ooooops ... not good: more than two reads for this template
	  getRead(rid).setTemplatePartnerID(-1);
	  getRead(firstpartner).setTemplatePartnerID(-1);
	  if(ol_gt2reads<200){
	    ++outlines_fatal;
	    cout << tidcounter[tnI->second] << " ";
	    cout << "DNA template " << tnI->first << " has more than two reads, template info not used. Read at fault: " << getRead(rid).getName() << endl;
	    if(++ol_gt2reads==200){
	      cout << "More than 200 cases like the above, will not report more.\n";
	    }
	  }
	}
      }else{
	getRead(rid).setTemplateID(acttid);

	tid_firstpartner[acttid]=rid;

	tidcounter[acttid]++;
	tnmap[getRead(rid).getTemplate()]=acttid;
	++acttid;
      }
    }
  }

  if(outlines_fatal>0 && nagwarn_templateproblems){
    static std::string emsg="Problems found with the template data, see output log for more info. This points some serious problem either with read naming (like unrecognised read naming scheme) or broken template information, please fix your input files!\nOr switch off this warning with -NW:ctp=no (but you'll do this at own risk!)";
    if(nagwarn_templateproblems==NWSTOP) {
      MIRANOTIFY(Notify::FATAL,emsg);
    }else{
      cout << "WARNING!\n" << emsg << endl;
    }
  }

  if(verbose){
    cout << "Generated " << acttid << " unique DNA template ids for " << validreads << " valid reads.\n";
  }

  FUNCEND();
  return acttid!=validreads;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadPool::checkTemplateIDs(const std::string & errmsg)
{
  FUNCSTART("void ReadPool::checkTemplateIDs(std::string & errmsg)");
  for(uint32 ri=0;ri< REP_thepool3.size(); ++ri){
    if(REP_thepool3.getRead(ri).getTemplatePartnerID()>=0
       && REP_thepool3.getRead(ri).getTemplateID() != REP_thepool3.getRead(REP_thepool3.getRead(ri).getTemplatePartnerID()).getTemplateID()){
      cout << "Ouch, template problem for read " << ri << " " << REP_thepool3.getRead(ri).getName() << ", dumping readpool for debug\n";
      dumpAs(cout,Read::AS_TEXT,true);
      BUGIFTHROW(true, errmsg);
    }
  }
}

/*************************************************************************
 *
 * return true if no duplicate was found
 *
 *************************************************************************/

bool ReadPool::checkForDuplicateReadNames(uint8 nagwarn_duplicatenames)
{
  FUNCSTART("bool ReadPool::checkForDuplicateReadNames(bool verbose)");

  typedef std::unordered_set<std::string> strset;
  strset rnset;
  decltype(rnset.end()) rnI;

  bool allok=true;
  uint32 outlines=0;
  for(uint32 i=0; i<size();i++){
    if(getRead(i).hasValidData()==false) continue;

    rnI=rnset.find(getRead(i).getName());
    if(rnI!=rnset.end()){
      allok=false;
      if(outlines<2000){
	cout << "Read " << *rnI << " is present more than once in the data set. Did you load a file twice in the manifest? Is a read present more than once in your file(s)?\n";
	if(++outlines==2000){
	  cout << "More than 2000 cases like the above, will not report more. Fix your input!\n";
	}
      }
    }else{
      rnset.insert(getRead(i).getName());
    }
  }

  if(outlines>0 && nagwarn_duplicatenames){
    static std::string emsg="Some read names were found more than once (see log above). This usually hints to a serious problem with your input and should really, really be fixed. You can choose to ignore this error with '-NW:cdrn=no', but this will almost certainly lead to problems with result files (ACE and CAF for sure, maybe also SAM) and probably to other unexpected effects.";
    if(nagwarn_duplicatenames==NWSTOP){
      MIRANOTIFY(Notify::FATAL,emsg);
    }else{
      cout << "WARNING!\n" << emsg << endl;
    }
  }

  FUNCEND();
  return allok;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadPool::adaptFASTQQualValues(uint32 rpstart, uint32 rpend, base_quality_t fastqoffset, bool verbose, uint32 longestread, base_quality_t seenminqual, base_quality_t seenmaxqual)
{
  if(verbose) cout << "Looking at FASTQ type ... ";cout.flush();

  base_quality_t qualcorrector=0;
  bool needoldsxamapping=false;

  bool checkqual=false;
  if(seenminqual==255 && seenmaxqual==0) checkqual=true;

  if(checkqual || longestread==0){
    for(size_t rpi=rpstart; rpi<rpend; ++rpi){
      longestread=std::max(longestread,getRead(rpi).getLenSeq());
      if(checkqual){
	for(auto & qv : const_cast<std::vector<base_quality_t> &>(getRead(rpi).getQualities())){
	  seenminqual=std::min(seenminqual,qv);
	  seenmaxqual=std::max(seenmaxqual,qv);
	}
      }
    }
  }

  if(fastqoffset<33){
    // longestread criterion catches everything non-Solexa (including 454, or contigs) as Illumina
    //  switched to Sanger style before attaining 200bp
    if(seenminqual<59 || seenmaxqual <= 73 || seenminqual==seenmaxqual || longestread>=200){
      if(verbose) cout << "guessing FASTQ-33 (Sanger)\n";
      qualcorrector=33;
    }else if(seenminqual>=59 && seenminqual<64){
      if(verbose) cout << "guessing FASTQ-59 (old Solexa)\n";
      needoldsxamapping=false;
    }else{
      if(verbose) cout << "guessing FASTQ-64 (Illumina)\n";
      qualcorrector=64;
    }
  }else{
    if(verbose) cout << "told it to be FASTQ-" << static_cast<uint16>(fastqoffset) << '\n';
    qualcorrector=fastqoffset;
  }

  if(verbose) cout << "Running quality values adaptation ... "; cout.flush();
  if(needoldsxamapping){
    for(size_t rpi=rpstart; rpi<rpend; ++rpi){
      for(auto & qv : const_cast<std::vector<base_quality_t> &>(getRead(rpi).getQualities())){
	qv=RP_sxa2phredmap[qv];
      }
    }
  }else{
    for(size_t rpi=rpstart; rpi<rpend; ++rpi){
      for(auto & qv : const_cast<std::vector<base_quality_t> &>(getRead(rpi).getQualities())){
	qv-=qualcorrector;
      }
    }
  }
  if(verbose) cout << "done." << endl;
}




/*************************************************************************
 *
 * loads names of reads from external file and deletes them from pool
 * if invertselection true, then delets those not in the file
 *
 * BEWARE SIDE EFFECT: deletes all reads from pool that have invalid data
 *
 * the data must be in a key (a value in one line)
 * key can be: read name, or filename of exp read or file name of caf read
 *   (may not contain spaces, sorry)
 *
 * line with # as first nonwhitespace character are comments and read over
 *
 *************************************************************************/

void ReadPool::deleteReadsByName(const std::string & nfile, bool invertselection)
{
  FUNCSTART("void ReadPool::InvalidateReadsByName(const std::string & nfile)");

  BUGIFTHROW(true,"this needs to be adapted to indirection");

//  //stringhash_t M;
//  typedef std::unordered_map<std::string, uint32> strintmap;
//  strintmap rnmap;
//  decltype(rnmap.end()) rnI;
//
//  for(uint32 i=0; i<size();i++){
//    if(!REP_thepool[i].getName().empty()) {
//      rnmap[REP_thepool[i].getName()]=i;
//    }
//  }
//
//  ifstream fin;
//  fin.open(nfile, std::ios::in|std::ios::ate);
//  if(!fin){
//    MIRANOTIFY(Notify::FATAL, "File not found: " << nfile);
//  }
//  fin.seekg(0, std::ios::beg);
//
//  if(invertselection){
//    for(uint32 i=0; i<size();i++){
//      REP_thepool[i].setValidData(false);
//    }
//  }
//
//  string readname, dummy;
//  while(GeneralIO::readKeyValue(fin, readname,dummy)){
//    rnI=rnmap.find(readname);
//    if(rnI!=rnmap.end()) {
//      if(invertselection){
//	REP_thepool[rnI->second].setValidData(true);
//      }else{
//	REP_thepool[rnI->second].setValidData(false);
//      }
//    }
//  }
//  fin.close();
//
//  if(size()){
//    for(int32 i=size()-1; i >= 0; i--){
//      if(!REP_thepool[i].hasValidData()){
//	auto I=REP_thepool.begin();
//	advance(I,i);
//	REP_thepool.erase(I);
//      }
//    }
//  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void ReadPool::mergeXMLTraceInfo(const std::string & xmlfile, uint8 nagwarn_templateproblems)
{
  FUNCSTART("void ReadPool::mergeXMLTraceInfo(const std::string & filename)");

  rpDateStamp();

  cout << "Merging data from XML trace info file " << xmlfile << " ...";
  cout.flush();

  //make this called!

  std::string id454="454";
  uint32 numfound=0;

  NCBIInfoXML nix;

  std::list<NCBIInfoXML::ncbitraceelements_t> traces;

  try{
    nix.readXMLFile(xmlfile, traces);
  }
  catch(Notify n) {
    cout << "\n\n\nMIRA tried to load a XML TRACEINFO file containing ancillary data, but failed.\n"
      "Loading ancillary data when using FASTA files as input is\n"
      "really,\n"
      "        really,\n"
      "                REALLY encouraged, and therefore MIRA sets this as default.\n"
      "\nHowever, if you are really sure that you do not want to load ancillary data\n"
      "in TRACEINFO files, you can switch it off.\n"
      "Either use '<technology>_SETTINGS -LR:mxti=no' (e.g. SANGER_SETTING -LR:mxti=no),\n"
      "or use the '-notraceinfo' quickswitch to kill loading traceinfo files for all\n"
      "types of sequencing technologies. (place it after -fasta and -job quickswitches)\n\n\n";
    n.handleError(THISFUNC);
  }


  cout << "Num reads: " << traces.size() << endl;

  typedef std::unordered_map<std::string, int32> strintmap;
  strintmap rnmap;
  decltype(rnmap.end()) rnI;

  cout << "Building hash table ... "; cout.flush();

  for(uint32 i=0; i<size();++i){
    if(!getRead(i).getName().empty()) {
      rnmap[getRead(i).getName()]=i;
    }
  }
  cout << "done." << endl;

  std::string acttracename;

  // loop through std::list<NCBIInfoXML::ncbitraceelements_t>
  for(const auto & tre : traces){
    rnI=rnmap.end();

    auto E=tre.elements.cbegin();
    auto ECD=tre.elements_cdata.cbegin();
    bool found=false;
    for(;!found && E!=tre.elements.end(); ++E, ++ECD) {
      if((*E == NCBIInfoXML::NCBIXML_TRACE_NAME
	  || *E == NCBIInfoXML::NCBIXML_TI)
	 && !ECD->empty()){
	if(*E == NCBIInfoXML::NCBIXML_TRACE_NAME){
	  acttracename=*ECD;
	} else if(*E == NCBIInfoXML::NCBIXML_TI){
	  acttracename="gnl|ti|"+*ECD;
	}
	rnI=rnmap.find(acttracename);

	if(rnI!=rnmap.end()){
	  numfound++;
	  found=true;
	}
      }
    }

    if(found){
      int32 idoffound=rnI->second;
      // cout << "Found " << REP_thepool[idoffound].getName() << endl;
      // Read::setCoutType(Read::AS_TEXTCLIPS);
      // cout << REP_thepool[idoffound];

      int32 insertsize=-1;
      int32 insertstdev=-1;
      //int32 inssizemin=-1;
      //int32 inssizemax=-1;

      uint8 seqtype=ReadGroupLib::SEQTYPE_SANGER;

      E=tre.elements.begin();
      ECD=tre.elements_cdata.begin();
      for(;E!=tre.elements.end(); ++E, ++ECD) {
	switch(*E) {
	case NCBIInfoXML::NCBIXML_TRACE_NAME : {
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_QUALITY_LEFT  : {
	  getRead(idoffound).setLQClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_QUALITY_RIGHT  : {
	  getRead(idoffound).setRQClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_VECTOR_LEFT  : {
	  getRead(idoffound).setLSClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_VECTOR_RIGHT  : {
	  getRead(idoffound).setRSClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_TEMPLATE_ID  : {
	  getRead(idoffound).setTemplate(ECD->c_str());
	  //cout<< *ECD << endl;
	  break;
	}
	case NCBIInfoXML::NCBIXML_TRACE_END  : {
	  if(strlen(ECD->c_str())>0){
	    switch(toupper(ECD->c_str()[0])){
	    case 'F': {
	      getRead(idoffound).setTemplateSegment(1);
	      break;
	    }
	    case 'R': {
	      getRead(idoffound).setTemplateSegment(255);
	      break;
	    }
	    case 'U' : // fall through
	    case 'N' : {
	      getRead(idoffound).setTemplateSegment(0);
	      break;
	    }
	    default : {
	      MIRANOTIFY(Notify::FATAL, "Illegal trace_end, it's not one of F/R/N(U) (or empty): " << ECD->c_str(););
	    }
	    }
	  }
	  break;
	}
	default : {
	  // Ooooops?
	}
	}
      }
      // cout << "After:\n";
      // cout << getRead(idoffound);

    }
  }

  makeTemplateIDs(nagwarn_templateproblems);

  cout << "Done merging XML data, matched " << numfound << " reads." << endl;

  rpDateStamp();

  FUNCEND();
  return;
}


/*************************************************************************
 *
 * ugly and slow, but works and is fast enough
 *
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void ReadPool::mergeSSAHA2SMALTVecScreenData(const std::string & ssahafile, bool issmalt, std::vector<MIRAParameters> & miraparams, const std::string & logname, const std::string & logprefix)
{
  FUNCSTART("void ReadPool::mergeSSAHA2VecScreenData(const std::string & ssahafile, bool issmalt, std::vector<MIRAParameters> & miraparams, const std::string & logname, const std::string & logprefix)");

  cout << "Merging vector screen data from ";
  if(issmalt){
    cout << "SMALT";
  }else{
    cout << "SSAHA2";
  }
  cout << " results file " << ssahafile << ":\n";

  CEBUG("Building hash table ... "); cout.flush();

  typedef std::unordered_map<std::string, int32> strmap;
  strmap rnmap;
  decltype(rnmap.end()) rnI;

  for(uint32 rpi=0; rpi<size();++rpi){
    if(!getRead(rpi).getName().empty()) {
      rnmap[getRead(rpi).getName()]=rpi;
    }
  }
  CEBUG("done." << endl);

  std::ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname, std::ios::out|std::ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  std::ifstream ssahafin(ssahafile, std::ios::in|std::ios::ate);
  if(!ssahafin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << ssahafile);
  }

  ProgressIndicator<std::streamsize> P (0, ssahafin.tellg(),1000);
  ssahafin.seekg(0, std::ios::beg);

  uint32 sd_score;
  std::string sd_readname;
  std::string sd_vecname;
  uint32 sd_rfrom;
  uint32 sd_rto;
  uint32 sd_vfrom;
  uint32 sd_vto;
  std::string sd_dir;
  uint32 sd_totalmatchsize;
  float  sd_percentmatch;
  uint32 sd_rlen;

  std::string token;
  std::string alstring;
  if(issmalt){
    alstring="alignment:";
  }else{
    alstring="ALIGNMENT:";
  }

  bool haserrors=false;

  while(!ssahafin.eof()){
    ssahafin >> token;
    if(ssahafin.eof()) break;
    if(P.delaytrigger()) P.progress(ssahafin.tellg());
    if(token.compare(0,alstring.size(),alstring) != 0) {
      getline(ssahafin,token);
      continue;
    }
    ssahafin >> sd_score >> sd_readname;

    if(ssahafin.eof()) break;

    // *sigh* allow for empty names
    sd_vecname.clear();
    {
      bool loopit=true;
      char tmp;

      ssahafin.get(tmp);
      loopit=true;
      do{
	ssahafin.get(tmp);
	if(ssahafin.eof()) break;
	if(tmp==' ' || tmp=='\t'){
	  loopit=false;
	}else{
	  sd_vecname.push_back(tmp);
	}
      }while(loopit);
    }

    if(ssahafin.eof()) break;

    ssahafin >> sd_rfrom
	     >> sd_rto
	     >> sd_vfrom
	     >> sd_vto
	     >> sd_dir
	     >> sd_totalmatchsize
	     >> sd_percentmatch
	     >> sd_rlen;

    if(ssahafin.eof()) break;

    CEBUG(sd_readname << '\t' << sd_rfrom << '\t' << sd_rto << '\n');

    bool foundname=false;
    rnI=rnmap.find(sd_readname);
    if(rnI==rnmap.end()) {
      CEBUG("Not found: " << sd_readname << endl);
      continue;
    }
    uint32 foundreadid=rnI->second;
    if(!getRead(foundreadid).hasValidData()) continue;

    Read actread(getRead(foundreadid));
    assembly_parameters const & as_params= miraparams[actread.getSequencingType()].getAssemblyParams();

    if(actread.getLenSeq() != sd_rlen){
      if(actread.isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA) && actread.getLenSeq() != sd_rlen+1) {
	cout << "\nError! The length of read " << actread.getName()
	     << " (" << actread.getLenSeq()
	     << ") does not match the length given in the SSAHA2/SMALT file ("
	     << sd_rlen << ")\nSSAHA2 line:"
	     << ' ' << token
	     << ' ' << sd_score
	     << ' ' << sd_readname
	     << ' ' << sd_vecname
	     << ' ' << sd_rfrom
	     << ' ' << sd_rto
	     << ' ' << sd_vfrom
	     << ' ' << sd_vto
	     << ' ' << sd_dir
	     << ' ' << sd_totalmatchsize
	     << ' ' << sd_percentmatch
	     << ' ' << sd_rlen << endl;
	haserrors=true;
      }
    }

    CEBUG("SSAHA2/SMALT line:"
	  << ' ' << token
	  << ' ' << sd_score
	  << " r: " << sd_readname
	  << " v: " << sd_vecname
	  << " # " << sd_rfrom
	  << ' ' << sd_rto
	  << ' ' << sd_vfrom
	  << ' ' << sd_vto
	  << ' ' << sd_dir
	  << ' ' << sd_totalmatchsize
	  << ' ' << sd_percentmatch
	  << ' ' << sd_rlen << endl);

    //Read::setCoutType(Read::AS_FASTA);
    //CEBUG(actread);
    //Read::setCoutType(Read::AS_CLIPPEDFASTA);
    //CEBUG(actread);

    // in SSAHA2 output, from rfrom may be > rto for reverse matches
    // swap in these cases
    if(sd_rfrom > sd_rto) std::swap(sd_rfrom,sd_rto);

    for(uint32 i=sd_rfrom-1; i<sd_rto; i++){
      bool domask=false;
      if(as_params.as_clip_ssahamerge_strictfrontclip >0
	 || as_params.as_clip_ssahamerge_strictendclip >0){
	if(as_params.as_clip_ssahamerge_strictfrontclip >0
	   && static_cast<int32>(i)<as_params.as_clip_ssahamerge_strictfrontclip) domask=true;
	if(as_params.as_clip_ssahamerge_strictendclip>0
	   && i>=actread.getLenSeq()-as_params.as_clip_ssahamerge_strictendclip) domask=true;
      }else{
	domask=true;
      }
      if(domask) actread.changeBaseInSequence('X',0,i);
    }
    //Read::setCoutType(Read::AS_FASTA);
    //CEBUG(actread);
    //Read::setCoutType(Read::AS_CLIPPEDFASTA);
    //CEBUG(actread);

    actread.setClipoffsToMaskedChars(as_params.as_clip_ssahamerge_gapsize,
				     as_params.as_clip_ssahamerge_maxfrontgap,
				     as_params.as_clip_ssahamerge_maxendgap,
				     false);
    //Read::setCoutType(Read::AS_CLIPPEDFASTA);
    //CEBUG(actread);

    if(actread.getLMClipoff() > getRead(foundreadid).getLSClipoff()){
      getRead(foundreadid).setLSClipoff(actread.getLMClipoff());
      CEBUG("clippyl\n");
      if(!logname.empty()){
	logfout << logprefix << " SSAHA2/SMALT clip left "
		<< actread.getName()
		<< " to: "
		<< getRead(foundreadid).getLSClipoff() << '\n';
      }
    }else{
      if(!logname.empty()){
	logfout << logprefix << "unchanged SSAHA2/SMALT clip left "
		<< actread.getName()
		<< " stays: "
		<< getRead(foundreadid).getLSClipoff() << '\n';
      }
    }
    if(actread.getRMClipoff() < getRead(foundreadid).getRSClipoff()){
      getRead(foundreadid).setRSClipoff(actread.getRMClipoff());
      CEBUG("clippyr\n");
      if(!logname.empty()){
	logfout << logprefix << " SSAHA2/SMALT clip right "
		<< actread.getName()
		<< " to: "
		<< getRead(foundreadid).getRSClipoff() << '\n';
      }
    }else{
      if(!logname.empty()){
	logfout << logprefix << "unchanged SSAHA2/SMALT clip right "
		<< actread.getName()
		<< " stays: "
		<< getRead(foundreadid).getRSClipoff() << '\n';
      }
    }

    //Read::setCoutType(Read::AS_TEXTSHORT);
    //CEBUG(getRead(foundreadid));
  }
  P.finishAtOnce();

  ssahafin.close();

  if(!logname.empty()){
    logfout.close();
  }

  cout << "\nDone merging SSAHA2 vector screen data." << endl;

  if(haserrors){
    MIRANOTIFY(Notify::FATAL,"There were errors in the SSAHA2 data, most probably the sequences used to screen are different from\nthe ones loaded now (see log above). Sorry, MIRA has to abort, please check your data.");
  }

  FUNCEND();
  return;
}
//#define CEBUG(bla)





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
Read & ReadPool::getRead(const std::string & name)
{
  return getRead(getReadIndex(name));
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void ReadPool::dumpAs(std::ostream & ostr, uint8 astype, bool alsoinvalids) const
{
  FUNCSTART("void ReadPool::dumpAs(std::ostream & ostr, uint8 astype, bool alsoinvalids) const)");

  if(astype==Read::AS_MAF){
    for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
      // use dumpReadGroupAsMAF() instead saveReadGroupAsMAF!
      ReadGroupLib::dumpReadGroupAsMAF(rgi,ostr);
    }
  }
  Read::setCoutType(astype);

  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    try{
      if(getRead(i).hasValidData() || alsoinvalids) ostr << getRead(i);
    }
    catch(Notify n){
      cout << "Ouch, ouch, ouch, ouch ... ouch. Error while output of read " << i << ", " << getRead(i).getName() << endl;
      n.handleError(THISFUNC);
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 * like above, but using save()
 * This is for miraconvert
 *
 *************************************************************************/
void ReadPool::saveAsMAF(std::ostream & ostr, bool alsoinvalids) const
{
  FUNCSTART("void ReadPool::dumpAs(std::ostream & ostr, uint8 astype, bool alsoinvalids) const)");

  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    ReadGroupLib::saveReadGroupAsMAF(rgi,ostr);
  }

  Read::setCoutType(Read::AS_MAF);
  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    try{
      if(getRead(i).hasValidData() || alsoinvalids) ostr << getRead(i);
    }
    catch(Notify n){
      cout << "Ouch, ouch, ouch, ouch ... ouch. Error while output of read " << i << ", " << getRead(i).getName() << endl;
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
void ReadPool::dumpPoolInfo(std::ostream & ostr) const
{
  FUNCSTART("void ReadPool::dumpPoolInfo(std::ostream & ostr)");

  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    if(getRead(i).hasValidData()) {
      ostr << i << '\t' << getRead(i).getName() << '\n';
    }else{
      ostr << i << "\tinvalid\n";
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
void ReadPool::dumpAsEXPs(std::string & dirname) const
{
  FUNCSTART("void ReadPool::dumpAsEXPs(std::string & dirname) const");

  if(ensureDirectory(dirname,true)){
    MIRANOTIFY(Notify::FATAL, "Could not make sure that directory '" << dirname << "' exists, aborting MIRA.");
  }

  std::ofstream fofnout((dirname+"/fofn"), std::ios::out | std::ios::trunc);

  std::string dummyAP="";
  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    if(getRead(i).hasValidData()) {
      {
	std::ofstream expout((dirname+"/"+getRead(i).getName()+".exp"), std::ios::out | std::ios::trunc);
	(const_cast<Read &>(getRead(i))).dumpAsGAP4DA(expout, dummyAP);
      }
      fofnout << getRead(i).getName() << ".exp" << endl;
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

void ReadPool::refreshNameIndex()
{
  FUNCSTART("void ReadPool::refreshNameIndex()");
  BUGIFTHROW(!REP_allownameindex,"Pool not allowed to make a name index?");

  REP_nameindex.clear();
  for(uint32 ri=0; ri<size(); ++ri){
    REP_nameindex[getRead(ri).getName()]=ri;
  }

  // check if no double names exist
  if(REP_nameindex.size() != size()){
    // Ooooooops, there is at least one double name
    // Rebuild and this time check to find out which reads
    // Then bail out.
    REP_nameindex.clear();
    for(uint32 ri=0; ri<size(); ++ri){
      auto mI=REP_nameindex.find(getRead(ri).getName());
      if(mI!=REP_nameindex.end()){
	MIRANOTIFY(Notify::FATAL,"Read '" << mI->first << "' present at least twice in read pool? Not good, you may not have more than one read with the same name.");
      }
      REP_nameindex[getRead(ri).getName()]=ri;
    }
  }

  FUNCEND();
}

//  // returns (first) read with name i.p.
//  Read & getRead(const std::string & name);
//

/*************************************************************************
 *
 * returns id of (first) read with name or -1
 *
 *
 *************************************************************************/

int32 ReadPool::getReadIndex(const std::string & name)
{
  if(REP_nameindex.empty()) refreshNameIndex();

  auto mI=REP_nameindex.find(name);
  if(mI!=REP_nameindex.end()){
    return mI->second;
  }
  return -1;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadPool::adjustIllegalQualities(base_quality_t bq)
{
  for(size_t i=0; i<REP_thepool3.size(); ++i){
    std::vector<base_quality_t> & bqv =const_cast<std::vector<base_quality_t>&>(getRead(i).getQualities());
    bool mustadjust=true;
    for(auto tbq: bqv){
      if(tbq<=100) {
	mustadjust=false;
	break;
      }
    }
    if(mustadjust){
      auto s=bqv.size();
      bqv.clear();
      bqv.resize(s,bq);
    }
  }

}
