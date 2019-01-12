/*
 * Copyright (c) Bastien Chevreux 2000
 * All rights reserved.
 *
 *
 * Written by Bastien Chevreux (BaCh)
 */


#include "io/scf.H"
#include "io/phd.H"
#include "io/fasta.H"
#include "util/dptools.H"
#include "util/progressindic.H"
#include "mira/read.H"

#include "version.H"

uint32 CLASS_trace;




/*
  Functionality to store read qualities in memory and hash
  readnames to allow for fast lookup of qualities
*/

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) == 0;
    }
};

typedef hash_multimap<const char*, int, hash<const char *>, eqstr> stringhash_t;
typedef stringhash_t::value_type stringhash_entry_t;

stringhash_t GLO_Fhashedreadnames;
vector<string> GLO_Flistofreadnames;
vector< vector<base_quality_t> > GLO_Fqualities;





void usage()
{
  cerr << "ac2caf V1.1 (MIRALIB version "MIRALIBVERSION")\n";
  cerr << "Author: Bastien Chevreux\t(bach@chevreux.org)\n\n";

  cerr << "Converts phrap-generated assembly (ACE) files to a CAF file written to stdout.\n\
Qualities for the reads are fetched either from PHD, fastaqual or SCF files. \n\
Tags in the ace file are currently not supported.\n\
\n\
Usage:\n\
      ace2caf [-d defaultqual] \n\
              [-f fastaqualpostfix] \n\
              [-F fastafile] \n\
              [-p phdpostfix] \n\
              -s scfpostfix <ace_file >caf_file\n\
\n\
The 'postfix' strings - which can also be an empty string (\"\") - will be \n\
appended to the name of the reads gathered from the ACE file to construct the\n\
corresponding PHD|fastaqual|SCF-filename.\n\
\n\
\t-d\tnumber\tdefault base quality if no quality file could be loaded\n\
\t\t\t must be >=0 and <= 100\n\
\t-f\tstring\tload base qualities from qual files (fasta format,\n\
\t\t\t containing qualities only) with given postfix\n\
\t-F\tstring\tload all base qualities from given qual file (fasta format,\n\
\t\t\t containing qualities only)\n\
\t-p\tstring\tload base qualities from PHD files with given postfix\n\
\t-s\tstring\tpostfix for SCF files. Additionally load qualities from \n\
\t\t\t SCF if -p or -q are not set or the PHD or FASTA files\n\
\t\t\t contained errors\n\
\n\
The order of -f and -p defines the order of the load attempts.\n\
\n\
E.g.:\n\
\tace2caf -s .scf <bla.ace >bla.caf\n\
will convert the ace file 'bla.ace' into the caf file 'bla.caf'. Qualities \n\
will be loaded the hypothetical SCF file K606.scf for the hypothetical read\n\
K606 contained in the ACE file bla.ace.\n\
\n\
\tace2caf -d 15 -F bla.fasta.qual -f .qual.1 -s .SCF <bla.ace >bla.caf\n\
will first try to see whether quality information was found in the (big)\n\
fasta quality file named 'bla.fasta.qual'. If not found there, it\n\
will try to load the qualities from the hypothetical qual file K606.qual.1 for\n\
the hypothetical read K606 contained in the ACE file bla.ace. The\n\
corresponding SCF file (which will only be loaded if the fastaqual was not\n\
found, but set in the CAF file) is K606.SCF. If no quality could be loaded from\n\
anywhere, the base quality values are all set to 15.\n";
  exit(2);
}


string GLO_phdpostfix;
string GLO_qualpostfix;
string GLO_scfpostfix;
string GLO_bigfastaqual;

enum {LOAD_PHDa2c=0, LOAD_QUALa2c, LOAD_SCFa2c, MAKE_DEFAULTQUALa2c};
vector<uint32> GLO_loadfrom;
uint8 GLO_defaultqual=10;

struct AF_t {
  string name;
  char   orientation;
  int32  startpos;
};
struct contig_t {
  string name;
  uint32 len;
  uint32 numreads;
  uint32 numbasesegments;
  char   orientation;
  string sequence;
  vector<uint8> qualities;
  vector<AF_t> afs;
  int32 minpos;

  bool ok;
};

struct read_t{
  string name;
  uint32 len;
  uint32 numinfoitems;
  uint32 numreadtags;

  string sequence;
  int32 qa_qualclipstart;
  int32 qa_qualclipend;
  int32 qa_alignclipstart;
  int32 qa_alignclipend;

  string desc;
  string scfname;
};

vector<read_t> GLO_reads;

contig_t GLO_contig;
uint32 GLO_stadenid=1;



void dumpContig()
{
  FUNCSTART("dumpContig()");

  if(!GLO_contig.ok) return;
//    throw Notify(Notify::INTERNAL, THISFUNC, ": wanted to dump a contig, but contig data is not available.");
//  }

  cout << "\nSequence : " << GLO_contig.name << "\nIs_contig\nPadded\n";

  // TODO: this is dog-slow for contigs with many reads (454 type data)
  //  rework!
  for(uint32 i=0; i<GLO_contig.afs.size(); i++) {
    for(uint32 j=0;j<GLO_reads.size(); j++){
      if(GLO_contig.afs[i].name.compare(GLO_reads[j].name)==0) {
	if(GLO_reads[j].qa_alignclipstart == -1 && GLO_reads[j].qa_alignclipend== -1) continue;
	cout << "Assembled_from " << GLO_contig.afs[i].name;
	if(GLO_contig.afs[i].orientation=='c') {
	  int32 cliplen=GLO_reads[j].qa_alignclipend-GLO_reads[j].qa_alignclipstart;
//	    int32 s1=GLO_contig.afs[i].startpos+GLO_reads[j].len-1;
//	    int32 s2=s1-cliplen;
//	    cout << " " << s1 << " " << s2;
//	    cout << " " << GLO_reads[j].qa_alignclipstart-GLO_reads[j].qa_alignclipstart+1;
//	    cout << " " << GLO_reads[j].qa_alignclipend-GLO_reads[j].qa_alignclipstart+1;
	  int32 r1=GLO_reads[j].len-GLO_reads[j].qa_alignclipend+1;
	  int32 r2=r1+cliplen;
	  int32 s1=GLO_contig.afs[i].startpos+GLO_reads[j].qa_alignclipend-1;
	  int32 s2=s1-cliplen;
	  cout << " " << s1 << " " << s2;
	  cout << " " << r1 << " " << r2;
	} else {
	  int32 s1=GLO_contig.afs[i].startpos+GLO_reads[j].qa_alignclipstart-1;
	  int32 s2=GLO_reads[j].qa_alignclipend-GLO_reads[j].qa_alignclipstart+s1;

	  cout << " " << s1 << " " << s2;
	  cout << " " << GLO_reads[j].qa_alignclipstart;
	  cout << " " << GLO_reads[j].qa_alignclipend;
	}
	cout << "\n";
	break;
      }
    }
  }

  cout << "\nDNA : " << GLO_contig.name;
  for(uint32 i=0; i<GLO_contig.sequence.size();i++){
    if(i%60 == 0) cout << endl;
    if(GLO_contig.sequence[i]=='*'){
      cout << '-';
    } else {
      cout << GLO_contig.sequence[i];
    }
  }

  cout << endl;

  GLO_contig.ok=false;

  FUNCEND();
}

string getContigFromFile()
{
  FUNCSTART("void getContigFromFile()");

  if(GLO_contig.ok) dumpContig();

  GLO_reads.clear();

  GLO_contig.name.erase();
  GLO_contig.len=0;
  GLO_contig.numreads=0;
  GLO_contig.numbasesegments=0;
  GLO_contig.orientation=' ';
  GLO_contig.sequence.erase();
  GLO_contig.qualities.clear();
  GLO_contig.afs.clear();
  GLO_contig.ok=false;

  cin >> GLO_contig.name;
  cin >> GLO_contig.len;
  cin >> GLO_contig.numreads;
  cin >> GLO_contig.numbasesegments;
  cin >> GLO_contig.orientation;

  GLO_contig.orientation=tolower(GLO_contig.orientation);

  //cout << "### Contig: " << GLO_contig.name << " " << GLO_contig.len << " " << GLO_contig.numreads << " " << GLO_contig.numbasesegments << " " << GLO_contig.orientation << endl;

  string tmpstr;

  while(GLO_contig.sequence.size()!=GLO_contig.len){
    cin >> tmpstr;
    if(cin.eof()){
      cerr <<  "Read " << GLO_contig.sequence.size() << " bases, but expected " << GLO_contig.len << endl;
      cerr <<  "Unexpected end of file while reading sequence for contig " << GLO_contig.name << endl;
      throw Notify(Notify::FATAL, THISFUNC);
    }
    GLO_contig.sequence+=tmpstr;
    for(uint32 i=GLO_contig.sequence.size()-tmpstr.size(); i<GLO_contig.sequence.size(); i++){
      if(!dptools::isValidIUPACStarBase(GLO_contig.sequence[i])) {
	cerr << "Base " << GLO_contig.sequence[i] << "(" << hex << static_cast<uint16>(GLO_contig.sequence[i]) << dec << ") at position " << i << " in contig " << GLO_contig.name << " is not a valid IUPAC base or a gap (*).\n";
	  throw Notify(Notify::FATAL, THISFUNC);
      }
    }
    if(GLO_contig.sequence.size() > GLO_contig.len){
      cerr << GLO_contig.sequence << endl;
      cerr <<  "Read " << GLO_contig.sequence.size() << " bases, but expected only " << GLO_contig.len << endl;
      cerr << "Unexpected sequence bases for contig " << GLO_contig.name << endl;
      cerr << "If the last bases of the sequence above are 'funny', then it's probable that\nthe sequence contained _not enough_ bases instead of too many.\n";
      throw Notify(Notify::FATAL, THISFUNC);
    }
  }

  //cout << "### Sequence: " << GLO_contig.sequence << endl;

  cin >> tmpstr;
  if(tmpstr.compare("BQ") != 0){
    cerr << "Missing BQ after sequence for contig " << GLO_contig.name << endl;
    throw Notify(Notify::FATAL, THISFUNC);
  }

  uint16 tmpp=0;
  for(uint32 i=0; i<GLO_contig.len; i++) {
    if(GLO_contig.sequence[i]=='*') {
      tmpp=0;
      //      cout << " *";
    } else {
      cin >> tmpp;
      //cout << tmpp << " ";
    }
    GLO_contig.qualities.push_back(static_cast<uint8>(tmpp));
    if(cin.eof()){
      cerr <<  "Unexpected end of file while reading base qualities for contig " << GLO_contig.name << endl;
      throw Notify(Notify::FATAL, THISFUNC);
    }
  }

  AF_t tmpaf;
  while(!cin.eof()){
    cin >> tmpstr;
    //    cout << "!!!!!!!!!!!!!!!!!! " << tmpstr;
    if(tmpstr.compare("AF") != 0) break;
    cin >> tmpaf.name;
    cin >> tmpaf.orientation;
    cin >> tmpaf.startpos;
    tmpaf.orientation=tolower(tmpaf.orientation);

    GLO_contig.afs.push_back(tmpaf);
    //cout << "AF " << tmpaf.name << " " << tmpaf.orientation << " " << tmpaf.startpos << endl;
  }

  while(!cin.eof()){
    if(tmpstr.compare("BS")!=0) break;
    //cout << "BS";
    cin >> tmpstr;
    //cout << " " << tmpstr;
    cin >> tmpstr;
    //cout << " " << tmpstr;
    cin >> tmpstr;
    //cout << " " << tmpstr << endl;
    cin >> tmpstr;
  }

  FUNCEND();

  GLO_contig.ok=true;

  return tmpstr;
}



void dumpQualValues(const read_t & read, const vector<uint8> quals)
{
  FUNCSTART("void dumpQualValues(const read_t & read, const vector<uint8> quals)");

  uint32 qualindex=0;
  for(uint32 i=0; i<read.sequence.size(); i++, qualindex++){
    if(read.sequence[i]=='*') {
      cout << "1 ";
      qualindex--;
    }else if (tolower(read.sequence[i])=='n'){
      cout << "0 ";
    } else {
      cout << (uint16) quals[qualindex] << " ";
    }
    if(i>0 && i%20==0) cout << "\n";
  }

  FUNCEND();
}


string getReadFromFile()
{
  FUNCSTART("string getReadFromFile()");

  GLO_reads.resize(GLO_reads.size()+1);
  read_t & read=GLO_reads.back();

  cin >> read.name;
  cin >> read.len;
  cin >> read.numinfoitems;
  cin >> read.numreadtags;

  read.qa_qualclipstart=1;
  read.qa_qualclipend=read.len;
  read.qa_alignclipstart=1;
  read.qa_alignclipend=read.len;

  string tmpstr;
  while(read.sequence.size()!=read.len){
    cin >> tmpstr;
    if(cin.eof()){
      cerr <<  "Read " << read.sequence.length() << " bases, but expected " << read.len << endl;
      cerr << "Unexpected end of file while reading sequence for read " << read.name << endl;
      throw Notify(Notify::FATAL, THISFUNC);
    }
    read.sequence+=tmpstr;
    for(uint32 i=read.sequence.size()-tmpstr.size(); i<read.sequence.size(); i++){
      if(!dptools::isValidIUPACStarBase(read.sequence[i])) {
	cerr << "Base " << read.sequence[i] << "(" << hex << static_cast<uint16>(read.sequence[i]) << dec << ") at position " << i << " in read " << read.name << " is not a valid IUPAC base or a gap (*).\n";
	  throw Notify(Notify::FATAL, THISFUNC);
      }
    }
    if(read.sequence.size()>read.len){
      cerr << read.sequence << endl;
      cerr <<  "Read " << read.sequence.length() << " bases, but expected only " << read.len << endl;
      cerr << "Unexpected sequence bases for read " << read.name << endl;
      cerr << "If the last bases of the sequence above are 'funny', then it's probable that\nthe sequence contained _not enough_ bases instead of too many.\n";
      throw Notify(Notify::FATAL, THISFUNC);
    }
  }

  char orientation='u';
  {
    for(uint32 i=0; i<GLO_contig.afs.size(); i++) {
      if(GLO_contig.afs[i].name.compare(read.name)==0) {
	orientation=GLO_contig.afs[i].orientation;
	break;
      }
    }

    //cerr << "Orient: " << orientation<< endl;

    if(orientation=='c') {
      //cerr << "Dooooooooooooooooo\n";
      tmpstr=read.sequence;
      uint32 rsi=read.sequence.size()-1;
      for(uint32 i=0;i<read.sequence.size(); i++, rsi--) {
	read.sequence[rsi]=dptools::getComplementIUPACBase(tmpstr[i]);
      }
    }
  }

  tmpstr.erase();

  while(!cin.eof()){
    cin >> tmpstr;
    if(tmpstr.length()==0) break;
    if(tmpstr.compare("CO") == 0) break;
    if(tmpstr.compare("RD") == 0) break;
    if(tmpstr.compare("CT{") == 0) break;
    if(tmpstr.compare("RT{") == 0) break;
    if(tmpstr.compare("WA{") == 0) break;
    if(tmpstr.compare("WR{") == 0) break;    // What's that one?
    if(tmpstr.compare("QA") == 0){
      int32 cs,ce;
      cin >> cs;
      cin >> ce;
      if(cs!=-1 && ce!=-1){
	read.qa_qualclipstart=cs;
	read.qa_qualclipend=ce;
      } else {
	if(orientation=='c') {
	  read.qa_qualclipstart=read.len-ce;
	  read.qa_qualclipend=read.len-cs;
	}
      }
      cin >> cs;
      cin >> ce;
      //if(cs!=-1 && ce!=-1){
	read.qa_alignclipstart=cs;
	read.qa_alignclipend=ce;
	//}

      // 5_002_e01.p1ca

//	if(orientation=='c') {
//	  read.qa_alignclipstart=read.len-ce;
//	  read.qa_alignclipend=read.len-cs;
//	} else {
//	  read.qa_alignclipstart=cs;
//	  read.qa_alignclipend=ce;
//	}
      continue;
    }
    if(tmpstr.compare("DS") == 0){
      getline(cin,read.desc);
      continue;
    }
    cerr << "Label \"" << tmpstr << "\" after read " << read.name << " not recognised, aborting.\nCheck also if sequence length of the read is correct.\n";
    throw Notify(Notify::FATAL, THISFUNC);
  }

  cout << "\nSequence : " << read.name;
  cout << "\nIs_read\nPadded\n";
  cout << "Staden_id " << GLO_stadenid++;
  cout << "\nClipping QUAL " << read.qa_qualclipstart << " " << read.qa_qualclipend << endl;
  //Align_to_SCF
  uint32 s1=1;
  uint32 s2=1;
  uint32 r1=1;
  uint32 r2=1;
  char oldbase=' ';
  uint32 numstars=0;
  for(uint32 i=1; i<read.sequence.size()+1; i++){
    if(read.sequence[i-1]=='*') {
      numstars++;
      if(oldbase!='*') {
	// need this because of 454 funny assemblies where
	//  aligned reads also start with a '*'?
	if(s2-1>0 && r2-1>0) {
	  cout << "Align_to_SCF " << s1 << " " << s2-1 << " " << r1 << " " << r2-1 << endl;
	}
	r1=i+1-numstars;
      }
      s1=i+1;
      s2++;
    } else {
      s2++;
      r2++;
    }
    oldbase=read.sequence[i-1];
  }
  cout << "Align_to_SCF " << s1 << " " << s2-1 << " " << r1 << " " << r2-1 << endl;

  string scfname=read.name+GLO_scfpostfix;
  cout << "SCF_File " << scfname << endl;
//  cout << "//Dye Look in DS!!!!\n";
//  cout << "//Template Look in DS!!!!\n";
//  cout << "//Tag *sigh*\n";

  cout << "\nDNA : " << read.name;


  for(uint32 i=0; i<read.sequence.size(); i++) {
    if(i%60 == 0) cout << endl;
    if(read.sequence[i]=='*'){
      cout << '-';
    } else {
      cout << read.sequence[i];
    }
  }

  cout << "\n\nBaseQuality : " << read.name << "\n";


  bool valuesdumped=false;

  // first try to lookup qualities that could have been loaded
  pair<stringhash_t::const_iterator, stringhash_t::const_iterator> p;
  p=GLO_Fhashedreadnames.equal_range(read.name.c_str());
  if(p.first!=GLO_Fhashedreadnames.end()){
    stringhash_t::const_iterator I = p.first;
    bool found=false;

    for (I = p.first; I != p.second; ++I) {
      if(GLO_Flistofreadnames[I->second]==read.name) {
	//cout << "found";
	found=true;
	break;
	//}else{
	//  cout << "Looked at " << thePool->getRead(I->second).getName() << endl;
      }
    }

    //cout << endl;

    dumpQualValues(read, GLO_Fqualities[I->second]);
    valuesdumped=true;
  }


  for(uint32 loadattempt=0; loadattempt<GLO_loadfrom.size() && !valuesdumped; loadattempt++){
    switch(GLO_loadfrom[loadattempt]) {
    case LOAD_SCFa2c: {
      SCF bla;

      try{
	bla.load(scfname.c_str());
      }
      catch (...){
	// error while loading
	break;
      }

      if(read.len-numstars != bla.getNumBases()) {
	cerr << "Read " << read.name << " in .ace file has " << read.len << " bases (containing " << numstars << " gaps).\nBut the SCF file " << scfname << " has " << bla.getNumBases() << " bases (without gaps).\nCannot pad because files are probably different reads or the assembly was\nmade with clipped sequences.\n";
	break;
      }

      uint32 baseindex=0;
      for(uint32 i=0; i<read.sequence.size(); i++, baseindex++){
	switch(tolower(read.sequence[i])){
	case 'a':{
	  //cout << "0 ";
	  cout << (uint16) bla.getAProb(baseindex)<< " ";
	  break;
	}
	case 'c':{
	  //cout << "0 ";
	  cout << (uint16) bla.getCProb(baseindex)<< " ";
	  break;
	}
	case 'g':{
	  //cout << "0 ";
	  cout << (uint16) bla.getGProb(baseindex)<< " ";
	  break;
	}
	case 't':{
	  //cout << "0 ";
	  cout << (uint16) bla.getTProb(baseindex)<< " ";
	  break;
	}
	case 'x':{
	  cout << (uint16) max(bla.getAProb(baseindex), max(bla.getCProb(baseindex), max(bla.getGProb(baseindex), bla.getTProb(baseindex)))) << " ";
	  break;
	}
	case 'n':{
	  cout << "0 ";
	  //      cout << "10 ";
	  break;
	}
	case '*':{
	  cout << "1 ";
	  baseindex--;
	  break;
	}
	default:
	  //cerr << "Unknown base???" << endl;
	  //throw Notify(Notify::FATAL, THISFUNC);
	  cout << (uint16) max(bla.getAProb(baseindex), max(bla.getCProb(baseindex), max(bla.getGProb(baseindex), bla.getTProb(baseindex)))) << " ";
	}
      }
      valuesdumped=true;
      break;
    }
    case LOAD_PHDa2c : {
      string phdfilename=read.name+GLO_phdpostfix;

      PHD thephd;
      try{
	thephd.load(phdfilename.c_str());
      }
      catch (...){
	// error while loading
	break;
      }

      if(read.len-numstars != thephd.getQualities().size()) {
	cerr << "Read " << read.name << " in .ace file has " << read.len << " bases (containing " << numstars << " gaps).\nBut the PHD file " << phdfilename << " has " << thephd.getQualities().size() << " bases (without gaps).\nCannot pad because files are probably different reads or the assembly was\nmade with clipped sequences.\n";
	break;
      }
      dumpQualValues(read, thephd.getQualities());
      valuesdumped=true;
      break;
    }
    case LOAD_QUALa2c :{
      string qualname=read.name+GLO_qualpostfix;

      FASTA thefasta;
      try{
	thefasta.loadQual(qualname.c_str());
      }
      catch (...){
	// error while loading
	break;
      }

      if(read.len-numstars != thefasta.getQualities().size()) {
	cerr << "Read " << read.name << " in .ace file has " << read.len << " bases (containing " << numstars << " gaps).\nBut the read in file " << qualname << " has " << thefasta.getQualities().size() << " bases (without gaps).\nCannot pad because files are probably different reads or the assembly was\nmade with clipped sequences.\n";
	break;
      }
      dumpQualValues(read, thefasta.getQualities());
      valuesdumped=true;

      break;
    }
    //case MAKE_DEFAULTQUALa2c :{
    //  vector<uint8> quals;
    //  quals.resize(read.len-numstars, GLO_defaultqual);
    //  dumpQualValues(read, quals);
    //  valuesdumped=true;
    //  break;
    //}
    }

    if(valuesdumped) break;
  }

  if(!valuesdumped){
    vector<uint8> quals;
    quals.resize(read.len-numstars, GLO_defaultqual);
    dumpQualValues(read, quals);
    cerr << "No quality values found for read " << read.name << ", using default.\n";
    //throw Notify(Notify::FATAL, THISFUNC);
  }

  cout << "\n\n";

  FUNCEND();

  return tmpstr;
}


void translateACE2CAF()
{
  FUNCSTART("void translateACE2CAF()");


  GLO_contig.ok=false;

  string tmpstr;

  bool found=false;
  while(!cin.eof()){
    cin >> tmpstr;
    if(tmpstr.compare("AS") == 0){
      found=true;
      break;
    }
  }
  if(!found) {
    throw Notify(Notify::FATAL, THISFUNC, "Inputfile does not seem to be an ACE file.");
  }

  uint32 numcontigs;
  uint32 numreads;
  cin >> numcontigs;
  cin >> numreads;

  cin >> tmpstr;
  while(!cin.eof() && tmpstr.length()!=0){
    if(tmpstr.compare("CO") == 0){
      tmpstr=getContigFromFile();
      continue;
    }
    if(tmpstr.compare("RD") == 0){
      tmpstr=getReadFromFile();
      continue;
    }
    if(tmpstr.compare("CT{") == 0){
      while(!cin.eof()){
	getline(cin,tmpstr);
	if(tmpstr.compare("}") == 0) break;
	if(cin.eof()){
	  throw Notify(Notify::FATAL, THISFUNC, "Unexpected end of file while waiting for \"}\" after a CT tag.");
	}
      }
      cin >> tmpstr;
      continue;
    }
    if(tmpstr.compare("RT{") == 0){
      while(!cin.eof()){
	getline(cin,tmpstr);
	if(tmpstr.compare("}") == 0) break;
	if(cin.eof()){
	  throw Notify(Notify::FATAL, THISFUNC, "Unexpected end of file while waiting for \"}\" after a RT tag.");
	}
      }
      cin >> tmpstr;
      continue;
    }
    if(tmpstr.compare("WA{") == 0){
      while(!cin.eof()){
	getline(cin,tmpstr);
	if(tmpstr.compare("}") == 0) break;
	if(cin.eof()){
	  throw Notify(Notify::FATAL, THISFUNC, "Unexpected end of file while waiting for \"}\" after a WA tag.");
	}
      }
      cin >> tmpstr;
      continue;
    }
    // wth is WR ???
    if(tmpstr.compare("WR{") == 0){
      while(!cin.eof()){
	getline(cin,tmpstr);
	if(tmpstr.compare("}") == 0) break;
	if(cin.eof()){
	  throw Notify(Notify::FATAL, THISFUNC, "Unexpected end of file while waiting for \"}\" after a WR tag.");
	}
      }
      cin >> tmpstr;
      continue;
    }
    cerr << "Label " << tmpstr << " not recognised, aborting.\n";
    exit(100);
  }

  dumpContig();


  FUNCEND();
}


void loadBigFastaQual (const string & filename)
{
  FUNCSTART("void loadBigFastaQual (const string & filename)");

  ifstream fin;
  fin.open(filename.c_str(), ios::in|ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  if(!fin.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << filename);
  }

  streamsize fsize=fin.tellg();

  {
    cerr << "Counting sequences in FASTA quality file:\n";

    uint32 numfasta=0;
    fin.seekg(0, ios::beg);
    ProgressIndicator<streamsize> P(0, fsize-1);
    FASTA thefasta;
    for(;true;numfasta++){
      thefasta.loadNextINTSeq(fin,255);
      P.progress(fin.tellg(),cerr);
      if(thefasta.testIfEmpty()) break;
    }

    P.finishAtOnce();
    cerr << "\n";

    // safety margin
    numfasta+=20;
    GLO_Flistofreadnames.reserve(numfasta);
    GLO_Fqualities.reserve(numfasta);
  }


  {
    cerr << "Loading quality data from FASTA quality file:\n";

    fin.clear();
    fin.seekg(0, ios::beg);

    ProgressIndicator<streamsize> P(0, fsize-1);
    FASTA thefasta;

    uint32 numquals=0;
    pair<stringhash_t::const_iterator, stringhash_t::const_iterator> p;
    for(;true;numquals++){
      thefasta.loadNextINTSeq(fin,255);
      P.progress(fin.tellg(),cerr);
      if(thefasta.testIfEmpty())  break;
      GLO_Fqualities.push_back(thefasta.getQualities());

      p=GLO_Fhashedreadnames.equal_range(thefasta.getQualName().c_str());
      if(p.first!=GLO_Fhashedreadnames.end()){
	cerr << "Read " << thefasta.getQualName() << " already loaded before!" << endl;
	MIRANOTIFY(Notify::FATAL, "Duplicate readname in FASTA quality file: " << thefasta.getQualName());
      }else{
	GLO_Flistofreadnames.push_back(thefasta.getQualName());
	GLO_Fhashedreadnames.insert(stringhash_entry_t(GLO_Flistofreadnames.back().c_str(), numquals));
      }

    }
    P.finishAtOnce(cerr);
    cerr << "\n";
  }
}


int main(int argc, char **argv)
{
  int c;
  extern char *optarg;
  extern int optind;

  bool dflag=false;
  bool pflag=false;
  bool fflag=false;
  bool sflag=false;

  bool flast=false;
  bool plast=false;

  while (1){
    c = getopt(argc, argv, "d:f:F:p:s:");
    if(c == -1) break;

    switch (c) {
    case 'd': {
      uint32 qual=atoi(optarg);
      if(qual>100) {
	cerr << argv[0] << ": " << "-d defaultquality must be >= 0 and <=100\n";
	usage();
      }
      GLO_defaultqual=atoi(optarg);
      dflag=true;
      break;
    }
    case 'f': {
      GLO_qualpostfix=optarg;
      fflag=true;
      flast=true;
      plast=false;
      break;
    }
    case 'F': {
      GLO_bigfastaqual=optarg;
      break;
    }
    case 'p': {
      GLO_phdpostfix=optarg;
      pflag=true;
      flast=false;
      plast=true;
      break;
    }
    case 's': {
      GLO_scfpostfix=optarg;
      sflag=true;
      break;
    }
    case '?':
      usage();
    }
  }
  if(!sflag) {
    cerr << argv[0] << ": " << "missing -s flag\n";
    usage();
  }
  if(fflag && pflag) {
    if(plast) {
      GLO_loadfrom.push_back(LOAD_QUALa2c);
      GLO_loadfrom.push_back(LOAD_PHDa2c);
    } else {
      GLO_loadfrom.push_back(LOAD_PHDa2c);
      GLO_loadfrom.push_back(LOAD_QUALa2c);
    }
  } else if(fflag) {
    GLO_loadfrom.push_back(LOAD_QUALa2c);
  } else if(pflag){
    GLO_loadfrom.push_back(LOAD_PHDa2c);
  }
  GLO_loadfrom.push_back(LOAD_SCFa2c);

  if(dflag){
    GLO_loadfrom.push_back(MAKE_DEFAULTQUALa2c);
  }

  try{
    if(GLO_bigfastaqual.size()!=0) loadBigFastaQual(GLO_bigfastaqual);
    translateACE2CAF();
  }
  catch(Notify n){
    n.handleError("main");
    exit(1);
  }


  return 0;
}
