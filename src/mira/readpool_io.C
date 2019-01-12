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

#include "mira/readpool_io.H"

#include "util/fileanddisk.H"
#include "caf/caf.H"
#include "mira/maf_parse.H"


using std::cout;
using std::cerr;
using std::endl;


// Plain vanilla constructor
ReadPoolIO::ReadPoolIO(ReadPool & rp)
{
  FUNCSTART("ReadPoolIO::ReadPoolIO(ReadPool & rp)");

  priv_zeroVars();
  priv_init();

  RPIO_rpptr=&rp;
  RPIO_clistptr=nullptr;
  RPIO_mpptr=nullptr;
  RPIO_loadstatus=LS_UNDEFINED;
  RPIO_loadtype=LT_UNDEFINED;

  RPIO_maf_parse=nullptr;
  RPIO_caf_parse=nullptr;

  RPIO_fastq_transformname=false;
  RPIO_fastqa_preservecomments=false;
  RPIO_fastqa_checkquals=true;

  RPIO_progressindic=nullptr;

  RPIO_missingfastaqual_resolvemsg="name your files with a .fna postfix";

  FUNCEND();
}

void ReadPoolIO::setNewReadPool(ReadPool & rp)
{
  RPIO_rpptr=&rp;
  if(RPIO_maf_parse!=nullptr) RPIO_maf_parse->setNewContainers(RPIO_rpptr,RPIO_clistptr,RPIO_mpptr);
  if(RPIO_caf_parse!=nullptr) RPIO_caf_parse->setNewContainers(RPIO_rpptr,RPIO_clistptr,RPIO_mpptr);
}

void ReadPoolIO::priv_zeroVars()
{
  FUNCSTART("void ReadPoolIO::zeroVars()");
  FUNCEND();
}

void ReadPoolIO::priv_init()
{
  FUNCSTART("void ReadPoolIO::init()");
  FUNCEND();
}



ReadPoolIO::~ReadPoolIO()
{
  FUNCSTART("ReadPoolIO::~ReadPoolIO()");

  discard();

  FUNCEND();
}


void ReadPoolIO::discard()
{
  FUNCSTART("ReadPoolIO::discard()");

  priv_zeroVars();

  priv_closeFiles();

  FUNCEND();
}


void ReadPoolIO::registerFile(const std::string & filename1, const std::string & optfilename2, const ReadGroupLib::ReadGroupID rgid, bool countonly)
{
  uint8 ziptype=0;
  std::string filetype;
  std::string dummytostem;
  std::string dummypathto;
  std::string fn2(optfilename2);

  guessFileAndZipType(filename1,dummypathto,dummytostem,filetype,ziptype);
  boost::to_lower(filetype);

  if(fn2.empty() && filetype=="fasta"){
    fn2=filename1+".qual";
  }
  registerFile(filetype,filename1, fn2, rgid, countonly);
}

void ReadPoolIO::registerFile(const std::string & filetype, const std::string & filename1, const std::string & optfilename2, const ReadGroupLib::ReadGroupID rgid, bool countonly)
{
  FUNCSTART("void ReadPoolIO::registerFile(const std::string & filetype, const std::string & filename1, const std::string & optfilename2, const ReadGroupLib::ReadGroupID rgid, bool countonly)");

  if(RPIO_loadstatus==LS_OPEN) {
    priv_closeFiles();
  }
  RPIO_loadstatus=LS_NOTOPEN;
  RPIO_filetype=filetype;
  RPIO_filename1=filename1;
  RPIO_optfilename2=optfilename2;
  RPIO_rgid=rgid;
  RPIO_countonly=countonly;

  if(!fileExists(RPIO_filename1)){
    MIRANOTIFY(Notify::FATAL,"Could not open " << RPIO_filetype << " file '" << RPIO_filename1 << "'. Is it present? Is it readable?");
  }
  if(getFileSize(RPIO_filename1)==0){
    MIRANOTIFY(Notify::FATAL,"Zero length  " << RPIO_filetype << " file '" << RPIO_filename1 << "'. This is fishy.");
  }

  priv_openFiles();
}


void ReadPoolIO::setAttributesForContigs(std::list<Contig>  * clist, std::vector<MIRAParameters> * mp)
{
  RPIO_clistptr=clist;
  RPIO_mpptr=mp;
  if(RPIO_maf_parse!=nullptr) RPIO_maf_parse->setNewContainers(RPIO_rpptr,RPIO_clistptr,RPIO_mpptr);
  if(RPIO_caf_parse!=nullptr) RPIO_caf_parse->setNewContainers(RPIO_rpptr,RPIO_clistptr,RPIO_mpptr);
}

void ReadPoolIO::priv_openFiles()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles()");
  if(RPIO_loadstatus!=LS_NOTOPEN) return;

  if(RPIO_filetype=="fastq"
     || RPIO_filetype=="fq"){
    priv_openFiles_fastq();
    RPIO_loadtype=LT_FASTQ;
  }else if(RPIO_filetype=="fasta"
	   || RPIO_filetype=="fna"
	   || RPIO_filetype=="fastanoqual"
	   || RPIO_filetype=="fa"){
    priv_openFiles_fasta();
    RPIO_loadtype=LT_FASTA;
  }else if(RPIO_filetype=="maf"){
    priv_openFiles_maf();
    RPIO_loadtype=LT_MAF;
  }else if(RPIO_filetype=="caf"){
    priv_openFiles_caf();
    RPIO_loadtype=LT_CAF;
  }else if(RPIO_filetype=="gbf"
	   || RPIO_filetype=="gbk"
	   || RPIO_filetype=="gb"
	   || RPIO_filetype=="gbff"
    ){
    priv_openFiles_gbf();
    RPIO_loadtype=LT_GBF;
  }else if(RPIO_filetype=="gff3"){
    priv_openFiles_gff3();
    RPIO_loadtype=LT_GFF3;
  }else if(RPIO_filetype=="fofnexp"){
    priv_openFiles_fofnexp();
    RPIO_loadtype=LT_FOFNEXP;
  }else if(RPIO_filetype=="exp"){
    priv_openFiles_exp();
    RPIO_loadtype=LT_EXP;
  }else{
    MIRANOTIFY(Notify::FATAL, "Unknown file type '" << RPIO_filetype << "'");
  }

  RPIO_totalreadsloaded=0;
  RPIO_loadstatus=LS_OPEN;
  RPIO_fsize=getFileSize(RPIO_filename1);
}

void ReadPoolIO::priv_closeFiles()
{
  FUNCSTART("void ReadPoolIO::priv_closeFiles()");
  if(RPIO_loadstatus!=LS_OPEN) return;

  if(RPIO_progressindic!=nullptr) RPIO_progressindic->finishAtOnce();

  switch(RPIO_loadtype){
  case LT_FASTQ : {
    kseq_destroy(RPIO_fastq_seq);
    gzclose(RPIO_fastq_fp);
    break;
  }
  case LT_FASTA : {
    RPIO_fasta_fin.close();
    if(RPIO_fasta_hasqualfile){
      RPIO_fasta_qin.close();
    }
    break;
  }
  case LT_GBF : {
    RPIO_gbf_ioobj.discard();
    break;
  }
  case LT_GFF3 : {
    RPIO_gff_ioobj.discard();
    break;
  }
  case LT_EXP : {
    // really nothing to do
    break;
  }
  case LT_FOFNEXP : {
    RPIO_fofnexp_names.clear();
    break;
  }
  case LT_MAF : {
    if(RPIO_maf_parse!=nullptr) delete RPIO_maf_parse;
    break;
  }
  case LT_CAF : {
    if(RPIO_caf_parse!=nullptr) delete RPIO_caf_parse;
    break;
  }
  default:{
    BUGIFTHROW(true,"close not implemented for " << RPIO_filetype);
  }
  }
  RPIO_loadstatus=LS_CLOSED;
}

void ReadPoolIO::priv_openFiles_fastq()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_fastq()");
  if(!RPIO_optfilename2.empty()){
    MIRANOTIFY(Notify::INTERNAL,"quality in filename2 not supported anymore");
    RPIO_fastq_qualoffset=static_cast<base_quality_t>(atoi(RPIO_optfilename2.c_str()));
  }

  //auto fsize=getFileSize(RPIO_filename1);

  RPIO_fastq_fp = gzopen(RPIO_filename1.c_str(), "r");
  if(RPIO_fastq_fp==Z_NULL) {
    MIRANOTIFY(Notify::FATAL,"Could not open FASTQ file '" << RPIO_filename1 << "' though it was possible just moments ago? Was it deleted?");
  }

  RPIO_fastq_seq = kseq_init(RPIO_fastq_fp);
}


void ReadPoolIO::priv_openFiles_fasta()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_fasta()");

  RPIO_fasta_fin.open(RPIO_filename1, std::ios::in);
  if(!RPIO_fasta_fin){
    MIRANOTIFY(Notify::FATAL, "Could not open: " << RPIO_filename1);
  }

  RPIO_fasta_hasqualfile=false;
  if(!RPIO_optfilename2.empty()){
    if(!fileExists(RPIO_optfilename2)){
      //MIRANOTIFY(Notify::FATAL,"Could not open " << RPIO_filetype << " quality file '" << RPIO_optfilename2 << "'. Is it present? Is it readable?");
      cout << "Could not find FASTA quality file " << RPIO_optfilename2;
      if(RPIO_fasta_wantsqualfiletoexist){
	cout << ", aborting. If you want to work without qualities, " << RPIO_missingfastaqual_resolvemsg;
	MIRANOTIFY(Notify::FATAL, "File not found: " << RPIO_optfilename2);
      }else{
	cout << ", using default values for these reads.\n";
      }
    }else{
      RPIO_fasta_hasqualfile=true;
      if(RPIO_fasta_wantsqualfiletoexist && getFileSize(RPIO_optfilename2)==0){
	MIRANOTIFY(Notify::FATAL, "FASTA quality file " << RPIO_optfilename2 << " has zero length? Seems fishy.");
      }
      RPIO_fasta_qin.open(RPIO_optfilename2, std::ios::in);
    }
  }else{
    if(RPIO_fasta_wantsqualfiletoexist){
      MIRANOTIFY(Notify::FATAL, "FASTA quality file expected to exist, but no quality filename given??? Refer to\nhttp://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html#sect_ref_manifest_readgroups\non how to load FASTA data without corresponding quality files.");
    }
  }
}

void ReadPoolIO::priv_openFiles_maf()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_maf()");

  if(RPIO_maf_parse!=nullptr) delete RPIO_maf_parse;
  RPIO_maf_parse=new MAFParse(RPIO_rpptr,RPIO_clistptr,RPIO_mpptr);
  RPIO_maf_parse->registerFile(RPIO_filename1);
  if(RPIO_progressindic!=nullptr) setAttributeProgressIndicator(true);
}

void ReadPoolIO::priv_openFiles_caf()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_caf()");
  if(RPIO_caf_parse!=nullptr) delete RPIO_caf_parse;
  RPIO_caf_parse=new CAF(RPIO_rpptr,RPIO_clistptr,RPIO_mpptr);
  RPIO_caf_parse->registerFile(RPIO_filename1);
  if(RPIO_progressindic!=nullptr) setAttributeProgressIndicator(true);
}

void ReadPoolIO::priv_openFiles_gbf()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_gbf()");

  RPIO_gbf_gbfloaded=false;
  RPIO_gbf_numtransferred=0;
}

void ReadPoolIO::priv_openFiles_gff3()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_gff3()");

  RPIO_gff_gffloaded=false;
  RPIO_gff_numtransferred=0;
}

void ReadPoolIO::priv_openFiles_fofnexp()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_fofnexp()");

  std::string justfilenamedummy;
  splitFullPathAndFileName(RPIO_filename1,RPIO_fofnexp_justpath,justfilenamedummy);

   // Load the file of filenames
  {
    std::ifstream fin(RPIO_filename1, std::ios::in|std::ios::ate);
    if(!fin){
      MIRANOTIFY(Notify::FATAL, "File not found: " << RPIO_filename1);
    }
    if(fin.rdstate()){
      MIRANOTIFY(Notify::FATAL, "Failed to open file: " << RPIO_filename1);
    }
    if(fin.tellg()==0){
      MIRANOTIFY(Notify::FATAL, "Zero length file: " << RPIO_filename1);
    }
    fin.seekg(0, std::ios::beg);

    RPIO_fofnexp_names.clear();
    std::string filename, dummy;
    while(GeneralIO::readKeyValue(fin, filename, dummy)){
      RPIO_fofnexp_names.push_back(filename);
    }
  }

  bool stopprocessing=false;
  {
    typedef std::unordered_set<std::string> strset;
    strset namemap;
    decltype(namemap.end()) nI;

    for(uint32 i=0; i< RPIO_fofnexp_names.size(); i++){
      nI=namemap.find(RPIO_fofnexp_names[i]);
      if(nI!=namemap.end()){
	cout << "WARNING: file " << RPIO_fofnexp_names[i] << " is present more than once in your file of filenames." << endl;
	stopprocessing=true;
      }else{
	namemap.insert(RPIO_fofnexp_names[i]);
      }
    }
  }

  if(stopprocessing){
    MIRANOTIFY(Notify::FATAL, "Some entries in your file of filenames lead to unrecoverable error: duplicate names. Aborting, see log above for further information.");
  }

  RPIO_fofnexp_nameiloaded=0;
}

void ReadPoolIO::priv_openFiles_exp()
{
  FUNCSTART("void ReadPoolIO::priv_openFiles_exp()");
  // really nothing to do?
}

uint64 ReadPoolIO::loadNextSeqs(uint64 numseqs, uint64 numcons, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::loadNextSeqs(uint64 numseqs, uint64 numcons)");

  uint64 retvalue=0;
  if(RPIO_loadstatus==LS_NOTOPEN){
    priv_openFiles();
  }
  if(RPIO_loadstatus==LS_OPEN){
    // load things
    if(RPIO_loadtype==LT_FASTQ){
      auto oldrpsize=RPIO_rpptr->size();
      retvalue=priv_loadNextSeqs_fastq(numseqs, lenseqstoload);
      if(RPIO_fastq_qualoffset) RPIO_rpptr->adaptFASTQQualValues(oldrpsize,RPIO_rpptr->size(),RPIO_fastq_qualoffset,false);
    }else if(RPIO_loadtype==LT_FASTA){
      retvalue=priv_loadNextSeqs_fasta(numseqs, lenseqstoload);
    }else if(RPIO_loadtype==LT_MAF){
      retvalue=priv_loadNextSeqs_maf(numseqs,numcons,lenseqstoload);
    }else if(RPIO_loadtype==LT_CAF){
      retvalue=priv_loadNextSeqs_caf(numseqs,numcons,lenseqstoload);
    }else if(RPIO_loadtype==LT_GBF){
      retvalue=priv_loadNextSeqs_gbf(numseqs, lenseqstoload);
    }else if(RPIO_loadtype==LT_GFF3){
      retvalue=priv_loadNextSeqs_gff3(numseqs, lenseqstoload);
    }else if(RPIO_loadtype==LT_FOFNEXP){
      retvalue=priv_loadNextSeqs_fofnexp(numseqs, lenseqstoload);
    }else if(RPIO_loadtype==LT_EXP){
      retvalue=priv_loadNextSeqs_exp();
    }else{
      BUGIFTHROW(true,"not implemented for " << RPIO_filetype);
    }
  }
  return retvalue;
}


void ReadPoolIO::setAttributeProgressIndicator(bool b)
{
  FUNCSTART("void ReadPoolIO::setAttributeProgressIndicator(bool b)");

  if(b && RPIO_progressindic==nullptr){
    RPIO_progressindic= new ProgressIndicator<int64>(0, 1,10000);
    if(RPIO_maf_parse!=nullptr) RPIO_maf_parse->setProgressIndicator(true);
    if(RPIO_caf_parse!=nullptr) RPIO_caf_parse->setProgressIndicator(true);
  }else if(!b && RPIO_progressindic!=nullptr){
    delete RPIO_progressindic;
    RPIO_progressindic=nullptr;
    if(RPIO_maf_parse!=nullptr) RPIO_maf_parse->setProgressIndicator(false);
    if(RPIO_caf_parse!=nullptr) RPIO_caf_parse->setProgressIndicator(false);
  }
  return;
}


uint64 ReadPoolIO::priv_loadNextSeqs_fastq(uint64 numseqstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_fastq(uint64 numseqstoload)");

  static multitag_t tmpcomm("COMM","","");
  tmpcomm.from=0;
  tmpcomm.to=0;

  bool fatalloaderror=false;
  bool qualerror=false;
  uint64 numseqsloaded=0;
  uint64 lenseqsloaded=0;
  int kseqretvalue=-1;

  std::vector<base_quality_t> fastq_bq;
  fastq_bq.reserve(1000);

  std::string fastq_tmpname;
  fastq_tmpname.reserve(100);
  std::string fastq_tmpcomment;
  fastq_tmpcomment.reserve(100);

  if(RPIO_totalreadsloaded==0 && RPIO_progressindic!=nullptr) RPIO_progressindic->reset(0,RPIO_fsize-1);

  while (numseqsloaded<numseqstoload && lenseqsloaded<lenseqstoload && (kseqretvalue = kseq_read(RPIO_fastq_seq)) >= 0) {
#ifdef HAVE_GZOFFSET
    if(RPIO_progressindic!=nullptr && RPIO_progressindic->delaytrigger()) RPIO_progressindic->progress(gzoffset(RPIO_fastq_fp));
#endif
    ++numseqsloaded;
    ++RPIO_totalreadsloaded;
    lenseqsloaded+=RPIO_fastq_seq->seq.l;
    if(RPIO_countonly) continue;

    Read & actread = RPIO_rpptr->getRead(RPIO_rpptr->provideEmptyRead());
    actread.setReadGroupID(RPIO_rgid);

    fastq_tmpcomment.clear();
    fastq_tmpname.clear();
    if(RPIO_fastq_seq->name.s!=nullptr) fastq_tmpname=RPIO_fastq_seq->name.s;
    if(RPIO_fastq_seq->comment.s!=nullptr) fastq_tmpcomment=RPIO_fastq_seq->comment.s;

    if(RPIO_rgid.wantUseReadNameFromComment()
       && !fastq_tmpcomment.empty()){
      auto bpos = fastq_tmpcomment.find_first_of(" \t");
      fastq_tmpname=fastq_tmpcomment.substr(0,bpos);
      if (bpos != std::string::npos) {
	fastq_tmpcomment=fastq_tmpcomment.substr(bpos+1,std::string::npos);
	boost::trim(fastq_tmpcomment);
      }
    }

    bool maybesolexa=(RPIO_rgid.getSequencingType() == ReadGroupLib::SEQTYPE_SOLEXA);

    if(RPIO_rgid.getSequencingType() == ReadGroupLib::SEQTYPE_TEXT
       && RPIO_fastq_seq->comment.l>0){
      uint32 numcolons=0;
      char * sptr=RPIO_fastq_seq->comment.s;
      for(; *sptr; ++sptr){
	if(*sptr==':') ++numcolons;
      }
      if(numcolons==3){
	sptr=RPIO_fastq_seq->comment.s;
	if(sptr[1]==':'
	   && (sptr[2]=='Y' || sptr[2]=='N')
	   && sptr[3]==':') {
	  maybesolexa=true;
	}
      }
    }

    if(RPIO_fastq_transformname && maybesolexa){
      auto bpos = fastq_tmpname.rfind("/");
      //cout << "tmpname: " << fastq_tmpname << endl;
      if (bpos == std::string::npos && RPIO_fastq_seq->comment.l>0) {
	//cout << "No / for " << fastq_tmpname << " ... need to make one:" << RPIO_fastq_seq->comment.s << endl;
	char * colonptr=RPIO_fastq_seq->comment.s;
	for(; *colonptr!=0; ++colonptr){
	  if(*colonptr==':') break;
	}
	if(*colonptr){
	  fastq_tmpname+='/';
	  colonptr=RPIO_fastq_seq->comment.s;
	  while(*colonptr!=':') {
	    fastq_tmpname+=*colonptr;
	    ++colonptr;
	  }
	}
      }
    }

    if(maybesolexa){
      actread.disallowAdjustments();
    }

    actread.setName(fastq_tmpname);

    if(actread.getName().empty()){
      cout << "Ouch, there's a read without a name? This is illegal. The sequence\n  "
	   << RPIO_fastq_seq->seq.s
	   << "\nmust have a name!\n";
      fatalloaderror=true;
    }

    if(RPIO_fastq_seq->seq.l==0){
      actread.setValidData(false);
    }else{
      actread.setSequenceFromString(RPIO_fastq_seq->seq.s,RPIO_fastq_seq->seq.l);

      if(RPIO_fastqa_preservecomments && !fastq_tmpcomment.empty()){
	tmpcomm.setCommentStr(fastq_tmpcomment);
	actread.addTagO(tmpcomm);
      }

      if(RPIO_fastq_seq->qual.l){
	if(RPIO_fastq_seq->qual.l != RPIO_fastq_seq->seq.l){
	  cout << actread.getName()
	       << ": different number of quality values than bases?\n";
	  qualerror=true;
	}else{
	  if(RPIO_fastqa_checkquals){
	    base_quality_t * qi = reinterpret_cast<base_quality_t *>(RPIO_fastq_seq->qual.s);
	    bool qualok=true;
	    for(;*qi; ++qi) {
	      if(unlikely(*qi<33 || *qi>164)){
		cout << "Read " << actread.getName() << ": invalid quality " << static_cast<uint16>(*qi) << '\n';
		qualok=false;
		break;
	      }
	    }
	    if(qualok) {
	      // This is such a hack ... all in the name of speed for mirabait loading
	      qi = const_cast<base_quality_t *>(&(actread.getQualities().front()));
	      memcpy(qi,RPIO_fastq_seq->qual.s,RPIO_fastq_seq->qual.l);
	      actread.setQualityFlag(true);
	    }else{
	      qualerror=true;
	    }
	  }else{
	    // This is such a hack ... all in the name of speed for mirabait loading
	    base_quality_t * qptr = const_cast<base_quality_t *>(&(actread.getQualities().front()));
	    memcpy(qptr,RPIO_fastq_seq->qual.s,RPIO_fastq_seq->qual.l);
	    actread.setQualityFlag(true);
	  }
	}
      }else{
	fastq_bq.clear();
	if(RPIO_fastq_qualoffset<33){
	  // ooops, trying to guess automatically ... not good if there's no sequence
	  // most probable nowadays: Sanger style FASTQ
	  fastq_bq.resize(RPIO_fastq_seq->seq.l,RPIO_rgid.getDefaultQual()+33);
	}else{
	  fastq_bq.resize(RPIO_fastq_seq->seq.l,RPIO_rgid.getDefaultQual()+RPIO_fastq_qualoffset);
	}
	actread.setQualities(fastq_bq);
	actread.setQualityFlag(RPIO_rgid.getDefaultQual()>0);
      }
    }
  }

  if(kseqretvalue<-1){
    cout << "Whoooops, something seems fishy with the last sequence loaded, the FASTQ parser returned " << kseqretvalue << " instead of the expected >= -1.\n";
    cout << "\nThis could be read: " << RPIO_fastq_seq->name.s << endl;
    cout << "Sequence string length: " << RPIO_fastq_seq->seq.l << endl;
    cout << "Quality string length: " << RPIO_fastq_seq->qual.l << endl;
    if(numseqsloaded>0 && !RPIO_countonly && RPIO_rpptr->size()>0){
      cout << "\nLast read which seemed OK: " << RPIO_rpptr->getRead(RPIO_rpptr->size()-1).getName() << endl;
    }
    if(RPIO_fastq_seq->seq.l != RPIO_fastq_seq->qual.l){
      MIRANOTIFY(Notify::FATAL,"FASTQ seems broken, there are reads where length of sequence does not match length of quality string. See log above.\n");
    }
    MIRANOTIFY(Notify::FATAL,"FASTQ seems broken, see log above.\n");
  }

  if(kseqretvalue==-1){
    priv_closeFiles();
  }

  if(qualerror){
    MIRANOTIFY(Notify::FATAL,"Unrecoverable error while loading data from FASTQ (see output above) Fix your input please.");
  }

  if(fatalloaderror) {
    MIRANOTIFY(Notify::FATAL, "Fatal error encountered during load of data (see log), aborting.\n") ;
  }

  return numseqsloaded;
}

uint64 ReadPoolIO::priv_loadNextSeqs_fasta(uint64 numseqstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_fasta(uint64 numseqstoload)");

  static multitag_t tmpcomm("COMM","","");
  tmpcomm.from=0;
  tmpcomm.to=0;

  uint64 numseqsloaded=0;
  uint64 lenseqsloaded=0;

  if(RPIO_totalreadsloaded==0 && RPIO_progressindic!=nullptr) RPIO_progressindic->reset(0,RPIO_fsize-1);

  auto formerpoolsize=RPIO_rpptr->size();

  while(numseqsloaded<numseqstoload && lenseqsloaded<lenseqstoload && !RPIO_fasta_fin.eof()) {
    RPIO_fasta_ioobj.loadNextSeq(RPIO_fasta_fin);
    if(RPIO_progressindic!=nullptr && RPIO_progressindic->delaytrigger()) RPIO_progressindic->progress(RPIO_fasta_fin.tellg());
    ++numseqsloaded;
    ++RPIO_totalreadsloaded;
    lenseqsloaded+=RPIO_fasta_ioobj.getSequence().size();
    if(RPIO_countonly) continue;
    Read & actread = RPIO_rpptr->getRead(RPIO_rpptr->provideEmptyRead());
    actread.setReadGroupID(RPIO_rgid);
    actread.setName(RPIO_fasta_ioobj.getSeqName());

    if(RPIO_rgid.getSequencingType() == ReadGroupLib::SEQTYPE_SOLEXA){
      actread.disallowAdjustments();
    }

    if(RPIO_fasta_ioobj.getSequence().empty()){
      actread.setValidData(false);
      cout << "\nWarning: read '" << actread.getName() << "' has no bases?! This usually points at some error in the processing of data before it arrives to MIRA.\n";
    }else{
      actread.setSequenceFromString(RPIO_fasta_ioobj.getSequence());
    }
    if(RPIO_fastqa_preservecomments && !RPIO_fasta_ioobj.getComment().empty()){
      tmpcomm.setCommentStr(RPIO_fasta_ioobj.getComment());
      actread.addTagO(tmpcomm);
    }

    if(actread.hasValidData()){
      actread.setQualities(RPIO_rgid.getDefaultQual());
      actread.setQualityFlag(false);
    }
  }

  bool endoffile=RPIO_fasta_ioobj.testIfEmpty();

  if(RPIO_fasta_hasqualfile && !RPIO_countonly){
    for(auto rpi=formerpoolsize; rpi<RPIO_rpptr->size(); ++rpi){
      RPIO_fasta_ioobj.loadNextINTSeq(RPIO_fasta_qin,255);
      //if(P.delaytrigger()) P.progress(fin.tellg());
      Read & actread=RPIO_rpptr->getRead(rpi);
      if(RPIO_fasta_ioobj.testIfEmpty()) {
	MIRANOTIFY(Notify::FATAL,"Premature end of FASTA quality file: did not find qualities for " << actread.getName());
      }
      if(actread.getName()!=RPIO_fasta_ioobj.getQualName()){
	MIRANOTIFY(Notify::FATAL,"FASTA quality file '" << RPIO_optfilename2 << "' is not in the same order as the FASTA file '" << RPIO_filename1 << "' itself: expected qualities for read '" << actread.getName() << "' but found qualities for '" << RPIO_fasta_ioobj.getQualName() << "'.\nPlease fix your input files.");
      }
      for(auto & iv : RPIO_fasta_ioobj.getINTValues()){
	if(iv<0) {
	  MIRANOTIFY(Notify::FATAL,"Whooops ... FASTA quality values <0? That's the very old Solexa scoring scheme and not supported anymore, sorry.\n");
	}
      }
      actread.setQualities(RPIO_fasta_ioobj.getQualities());
    }
  }

  if(endoffile){
    priv_closeFiles();
  }

  return numseqsloaded;
}


/*************************************************************************
 *
 * MAF is one of the formats which may contain contigs and/or reads
 * When loaded via the readpool io, all reads are loaded, but no contigs
 *  created.
 * E.g.: fle contains 1 contig (2 reads) and one 1 read without contig
 *  -> 3 reads are added to readpool
 *
 *************************************************************************/

uint64 ReadPoolIO::priv_loadNextSeqs_maf(uint64 numseqstoload, uint64 numconstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_maf(uint64 numseqstoload)");

  auto retvalue=RPIO_maf_parse->loadNextSeqs(numseqstoload,numconstoload,lenseqstoload);
  if(RPIO_maf_parse->checkIfEOF()){
    priv_closeFiles();
  }

  return retvalue;
}


/*************************************************************************
 *
 * CAF is one of the formats which may contain contigs and/or reads
 * When loaded via the readpool io, all reads are loaded, but no contigs
 *  created.
 * E.g.: fle contains 1 contig (2 reads) and one 1 read without contig
 *  -> 3 reads are added to readpool
 *
 *************************************************************************/

uint64 ReadPoolIO::priv_loadNextSeqs_caf(uint64 numseqstoload, uint64 numconstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_caf(uint64 numseqstoload, uint64 numconstoload)");
  auto retvalue=RPIO_caf_parse->loadNextSeqs(numseqstoload,numconstoload,lenseqstoload);
  if(RPIO_caf_parse->checkIfEOF()){
    priv_closeFiles();
  }
  return retvalue;
}


/*************************************************************************
 *
 * At the moment GBF does not provide a mechanism for loading sequences
 *  one by one. So the ReadPoolIO will have a GBF which contains a fully
 *  loaded object and release it as desired to the ReadPool.
 *
 * Probably less important as GBFs will not be huge files.
 *
 *************************************************************************/

uint64 ReadPoolIO::priv_loadNextSeqs_gbf(uint64 numseqstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_gbf(uint64 numseqstoload)");

  if(!RPIO_gbf_gbfloaded){
    RPIO_gbf_ioobj.load(RPIO_filename1);
    RPIO_gbf_ioobj.transferGeneInfoToCDSInfo();
    RPIO_gbf_gbfloaded=true;
  }

  uint64 numseqsloaded=0;
  uint64 lenseqsloaded=0;
  for(; RPIO_gbf_numtransferred<RPIO_gbf_ioobj.getNumSequences() && numseqsloaded<numseqstoload && lenseqsloaded<lenseqstoload; ++numseqsloaded,++RPIO_gbf_numtransferred){
    Read & actread=RPIO_rpptr->getRead(RPIO_rpptr->provideEmptyRead());
    actread.setReadGroupID(RPIO_rgid);

    //cout << "Read GBF " << i << endl;
    //cout << "\tName   : " << RPIO_gbf_ioobj.getSequenceName(i) << endl;
    //cout << "\tLenseq : " << RPIO_gbf_ioobj.getSequence(i).size() << endl;
    //cout << "\tNumtags: " << RPIO_gbf_ioobj.getTags(i).size() << endl;

    actread.setName(RPIO_gbf_ioobj.getSequenceName(RPIO_gbf_numtransferred));
    actread.setSequenceFromString(RPIO_gbf_ioobj.getSequence(RPIO_gbf_numtransferred));
    actread.setTags(RPIO_gbf_ioobj.getTags(RPIO_gbf_numtransferred));
    actread.setQualities(RPIO_rgid.getDefaultQual());
    lenseqsloaded+=actread.getLenSeq();
  }

  if(RPIO_gbf_numtransferred==RPIO_gbf_ioobj.getNumSequences()){
    priv_closeFiles();
  }

  return numseqsloaded;
}


/*************************************************************************
 *
 * Like GBF, GFF3 does not provide a mechanism for loading sequences
 *  one by one. So the ReadPoolIO will have a GFF which contains a fully
 *  loaded object and release it as desired to the ReadPool.
 *
 * Probably less important as GBFs will not be huge files.
 *
 *************************************************************************/
uint64 ReadPoolIO::priv_loadNextSeqs_gff3(uint64 numseqstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_gff3(uint64 numseqstoload)");

  //GFFParse thegff(this);
  //thegff.loadFile(filename);

  if(!RPIO_gff_gffloaded){
    RPIO_gff_ioobj.loadFile(RPIO_filename1);
    RPIO_gff_gffloaded=true;
  }


  uint64 numseqsloaded=0;
  uint64 lenseqsloaded=0;
  for(; RPIO_gff_numtransferred<RPIO_gff_ioobj.getNumSequences() && numseqsloaded<numseqstoload && lenseqsloaded<lenseqstoload; ++numseqsloaded,++RPIO_gff_numtransferred){
    if(RPIO_gff_ioobj.getSequence(RPIO_gff_numtransferred).size()==0 && !RPIO_gff_ioobj.getTags(RPIO_gff_numtransferred).empty()){
      // special case: empty sequence but tags means: the GFF had descriptions
      //  but no sequence for it. One of the already loaded reads however has a name
      //  matching, so simply add the tags to that read
      auto rid=RPIO_rpptr->getReadIndex(RPIO_gff_ioobj.getSequenceName(RPIO_gff_numtransferred));
      BUGIFTHROW(rid<0,"rid<0 ? should not happen here");
      Read & actread=RPIO_rpptr->getRead(rid);
      for(auto & ttag : RPIO_gff_ioobj.getTags(RPIO_gff_numtransferred)){
	actread.addTagO(ttag);
      }
    }else{
      Read & actread=RPIO_rpptr->getRead(RPIO_rpptr->provideEmptyRead());
      actread.setReadGroupID(RPIO_rgid);

      std::string minft_strainname;
      std::string minft_seqtypename;
      std::string minft_machinetype;
      int8 minft_tplacementcode;
      bool minft_isbb;
      bool minft_israil;
      bool minft_isCER;

      // ugly, but we know what we do and do not want to copy around the tags
      // unnecessarily
      std::vector<multitag_t> & gff3tags=const_cast<std::vector<multitag_t> & >(RPIO_gff_ioobj.getTags(RPIO_gff_numtransferred));

      if(Read::extractMINFTagInfo(gff3tags,
				  RPIO_gff_ioobj.getSequenceName(RPIO_gff_numtransferred),
				  minft_strainname,
				  minft_seqtypename,
				  minft_machinetype,
				  minft_tplacementcode,
				  minft_isbb,
				  minft_israil,
				  minft_isCER)){

	std::string dummy_empty;

	uint8 st=ReadGroupLib::stringToSeqType(minft_seqtypename);
	if(st==ReadGroupLib::SEQTYPE_END && !minft_seqtypename.empty()){
	  MIRANOTIFY(Notify::FATAL,"File " << RPIO_filename1 << " has invalid sequencing type '" << minft_seqtypename << "' in MIRA MIT2 tag.\n");
	}

	ReadGroupLib::ReadGroupID rgidt=ReadGroupLib::searchExactRGMatch(
	  dummy_empty,
	  st,
	  -1,
	  -1,
	  minft_tplacementcode,
	  minft_strainname,
	  minft_isbb,
	  minft_israil,
	  minft_isCER,
	  dummy_empty,
	  dummy_empty,
	  dummy_empty);

	if(rgidt.isDefaultNonValidReadGroupID()){
	  rgidt=ReadGroupLib::newReadGroup();
	  rgidt.setGroupName(dummy_empty);
	  rgidt.setSequencingType(st);
	  rgidt.setStrainName(minft_strainname);
	  rgidt.setBackbone(minft_isbb);
	  rgidt.setRail(minft_israil);
	  rgidt.setCoverageEquivalentRead(minft_isCER);
	  rgidt.setMachineType(minft_machinetype);
	}
      }

      actread.setName(RPIO_gff_ioobj.getSequenceName(RPIO_gff_numtransferred));
      actread.setSequenceFromString(RPIO_gff_ioobj.getSequence(RPIO_gff_numtransferred));
      actread.setTags(RPIO_gff_ioobj.getTags(RPIO_gff_numtransferred));
      actread.setQualities(RPIO_rgid.getDefaultQual());
      lenseqsloaded+=actread.getLenSeq();
    }
  }

  RPIO_rpptr->allowNameIndex(false);

  FUNCEND();
  return numseqsloaded;
}


/*************************************************************************
 *
 * In theory, EXP *could* contain multiple reads, but never encountered
 * in the wild and it's not the time to start that.
 *
 *************************************************************************/

uint64 ReadPoolIO::priv_loadNextSeqs_fofnexp(uint64 numseqstoload, uint64 lenseqstoload)
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_fofnexp(uint64 numseqstoload)");

  uint64 numseqsloaded=0;

  for(;RPIO_fofnexp_nameiloaded<RPIO_fofnexp_names.size() && numseqsloaded<numseqstoload; ++RPIO_fofnexp_nameiloaded){
    //P.progress(RPIO_fofnexp_nameiloaded);
    if(RPIO_fofnexp_names[RPIO_fofnexp_nameiloaded].empty()) continue;
    try{
      Read & actread = RPIO_rpptr->getRead(RPIO_rpptr->provideEmptyRead());
      actread.setReadGroupID(RPIO_rgid);

      actread.loadDataFromEXP(RPIO_fofnexp_names[RPIO_fofnexp_nameiloaded],RPIO_fofnexp_justpath);
      actread.transferSVTagsToClip(20,60);

      ++numseqsloaded;
    }
    catch(Flow){
    }
    catch(Notify n){
      n.handleError(THISFUNC);
    }
  }

  if(RPIO_fofnexp_nameiloaded>=RPIO_fofnexp_names.size()){
    //P.finishAtOnce();
    priv_closeFiles();
  }

  return numseqsloaded;
}


/*************************************************************************
 *
 * In theory, EXP *could* contain multiple reads, but never encountered
 * in the wild and it's not the time to start that.
 *
 *************************************************************************/
uint64 ReadPoolIO::priv_loadNextSeqs_exp()
{
  FUNCSTART("uint64 ReadPoolIO::priv_loadNextSeqs_exp()");

  uint64 retvalue=0;
  try{
    Read & actread = RPIO_rpptr->getRead(RPIO_rpptr->provideEmptyRead());
    actread.setReadGroupID(RPIO_rgid);

    actread.loadDataFromEXP(RPIO_filename1,"");
    actread.transferSVTagsToClip(20,60);

    retvalue=1;
  }
  catch(Flow){
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
  return retvalue;
}


/*
    loadDataFromGFF3_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="fofnexp"){
    loadDataFromFOFNEXP_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="exp"){
    loadDataFromEXP_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="scf"){
    // scf is actually just a data directory
    cout << "Setting to " << filename1 << endl;
    const_cast<ReadGroupLib::ReadGroupID &>(rgid).setDataDir(filename1);
  }else if(filetype.empty()){
    MIRANOTIFY(Notify::FATAL,"While trying to load data for readgroup, an empty string was given as file type (should have been fasta, fastq, etc.pp).\nReadgroup:" << rgid << endl);
  }else{
    MIRANOTIFY(Notify::FATAL,"While trying to load data for readgroup, type " << filetype << " not known.");
  }

  FUNCEND();

  //TODO: should we return the number of leaded reads?
  return 0;
}
*/



//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//ReadPoolIO::ReadPoolIO(ReadPoolIO const &other)
//{
//  FUNCSTART("ReadPoolIO::ReadPoolIO(ReadPoolIO const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//ReadPoolIO const & ReadPoolIO::operator=(ReadPoolIO const & other)
//{
//  FUNCSTART("ReadPoolIO const & ReadPoolIO::operator=(ReadPoolIO const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//std::ostream & operator<<(std::ostream &ostr, ReadPoolIO const &???)
//{
//  FUNCSTART("friend std::ostream & ReadPoolIO::operator<<(std::ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}
