/*
 * Written by Bastien Chevreux (BaCh)
 *
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

// 	$Id$

#include "io/fasta.H"

// for boost::trim, split
#include <boost/algorithm/string.hpp>

using std::cerr;
using std::endl;


// Plain vanilla constructor
FASTA::FASTA()
{
  FUNCSTART("FASTA::FASTA()");

  zeroVars();
  init();

  FUNCEND();
}

void FASTA::zeroVars()
{
  FUNCSTART("void FASTA::zeroVars()");

  FA_fastaseqname.clear();
  FA_qualseqname.clear();
  FA_comment.clear();
  FA_sequence.clear();
  FA_intvalues.clear();

  FUNCEND();
}

void FASTA::init()
{
  FUNCSTART("void FASTA::init()");
  FUNCEND();
}



FASTA::~FASTA()
{
  FUNCSTART("FASTA::~FASTA()");

  discard();

  FUNCEND();
}


void FASTA::discard()
{
  FUNCSTART("FASTA::discard()");

  zeroVars();

  FUNCEND();
}



/*************************************************************************
 *
 * data downloaded from NCBI has one problem: the name of the sequences
 *  in fasta file is the gnl|ti number ... the real read name is some-
 *  what later, preceded by the string " name:"
 * this function searches for reads with that characteristics and sets
 *  the name of the string given to the "real" one
 *
 * extracts name
 * also splits line into name and comment if comment!=nullptr
 *
 *************************************************************************/

void FASTA::adjustNameNCBIHack(std::string & name, std::string * comment)
{
  FUNCSTART("void FASTA::adjustNameNCBIHack(std::string & name)");

  static const std::string blanks=" \t";
  if(name.size()>7
     && name[0] == 'g'
     && name[1] == 'n'
     && name[2] == 'l'
     && name[3] == '|'
     && name[4] == 't'
     && name[5] == 'i'
     && name[6] == '|'){

    auto tokenstart=std::string::npos;

    tokenstart=name.find(" name:",0);

    if(tokenstart!=std::string::npos){
      // great, found a name: description
      auto tokenend=std::string::npos;
      tokenstart+=6;
      tokenend=name.find_first_of(blanks,tokenstart);
      if(tokenend==std::string::npos) tokenend=name.size();

      auto token=name.substr(tokenstart, tokenend-tokenstart);
      swap(name,token);
      return;
    }
  }

  // ok, no gnl|ti|
  // but maybe a NCBI GI "name"?
  if(name.size() > 3
     && name[0]=='g'
     && name[1]=='i'
     && name[2]=='|'){
    // OK, this *might* be a GI line
    std::vector<std::string> subnames;
    boost::split(subnames, name, boost::is_any_of("|"));
    if(subnames.size()==5){
      // bloody well looks like an NCBI style name line. Extract the real sequence name
      swap(name,subnames[3]);
      return;
    }
  }

  // nothing of the above, get only the first "token"
  //  as the name, the rest as comment
  if(!name.empty()){
    auto tokenend=std::string::npos;
    tokenend=name.find_first_of(blanks,0);
    if(tokenend!=std::string::npos){
      if(comment!=nullptr && tokenend<name.size()-1){
	*comment=name.substr(tokenend+1, std::string::npos);
	boost::trim(*comment);
      }
      auto token=name.substr(0, tokenend);
      name=token;
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

void FASTA::loadNextSeq(std::ifstream & fin)
{
  FUNCSTART("void FASTA::loadNextSeq(std::ifstream & fin)");

  FA_fastaseqname.clear();
  FA_sequence.clear();
  FA_comment.clear();

  bool read_read=false;

  bool nl=true;
  bool inseq=false;
  char inchar;

  std::vector<char> badchars;
  std::vector<std::streampos> badpos;

  //cout << "------------\n";

  while(!(fin.get(inchar)).eof() && !read_read) {
    //cout << "++" << inchar << "++ "; cout.flush();
    switch(inchar) {
    case ' ' : continue;
    case '\t' : continue;
    case '\r' : continue;      // that's from MS-DOS Files, where \n is \r\n
    case '\n' : {
      nl=true;
      continue;
    }
    case '>' : {
      if(inseq==true) {
	read_read=true;
	fin.unget();
	// FIXME: why? where does this char get eaten?
	fin.unget();

	//fin.putback('>');
	//// FIXME: why? where does this char get eaten?
	//fin.putback(' ');
	break;
      }
      if(nl==true) {
	while(!(fin.get(inchar)).eof()) {
	  //cout << "--" << inchar << "-- "; cout.flush();
	  if(inchar!='\n' && inchar!='\r'){
	    FA_fastaseqname+=inchar;
	  } else {
	    fin.unget();
	    //fin.putback(inchar);
	    break;
	  }
	}
	if(FA_fastaseqname.size() == 0) {
	  MIRANOTIFY(Notify::FATAL,"Missing name of fasta sequence at file byte position " << fin.tellg());
	}
	//cout << "RawName: " << FA_fastaseqname << endl;
	adjustNameNCBIHack(FA_fastaseqname, &FA_comment);
	nl=true;
	inseq=true;
      } else {
	MIRANOTIFY(Notify::FATAL,"Illegal character (" << inchar << ": " << std::hex << static_cast<uint16>(inchar) << std::dec << ") in fasta sequence name at file byte position " << fin.tellg() << endl);
      }
      break;
    }
    default : {
      if(inseq==false || nl==false) {
	MIRANOTIFY(Notify::FATAL,"Illegal character (" << inchar << ": " << std::hex << static_cast<uint16>(inchar) << std::dec << ") at begin of fasta sequence at file byte position " << fin.tellg());
      }
      fin.unget();
      //fin.putback(inchar);
      badchars.clear();
      badpos.clear();
      while(!(fin.get(inchar)).eof()) {
	//cout << "IC: " << inchar<<endl;
	if(inchar== '>' && nl==true) break;
	if(inchar== ' ') continue;
	if(inchar== '\t') continue;
	if(inchar== '\r' ) continue;           // from MS-DOS ... bla bla *sigh*
	nl=false;
	if(inchar== '\n') {
	  nl=true;
	  //cout << "dingNL:" << nl;
	  continue;
	}
	//cout <<"...";
	switch(toupper(inchar)){
	case 'A' :
	case 'C' :
	case 'G' :
	case 'T' :
	case 'N' :
	case 'X' :
	case 'R' :
	case 'Y' :
	case 'M' :
	case 'S' :
	case 'K' :
	case 'W' :
	case 'H' :
	case 'B' :
	case 'V' :
	case 'D' :
	case '*' : break;
	case '-' :
	  inchar='*';
	  break;
	case 'U' :
	  inchar='T';
	  break;
	default : {
	  if(badchars.size()<=100){
	    badchars.push_back(inchar);
	    badpos.push_back(fin.tellg());
	  }
	}
	}
	FA_sequence+=inchar;
      }
      inseq=false;
      if(!fin.eof()) {
	fin.unget();
	// FIXME: why? where does this character get eaten?
	fin.unget();

	// that doesn't work with gcc3.2
	//fin.putback(inchar);
	//// FIXME: why? where does this character get eaten?
	//fin.putback(' ');
      }
      //cout << FA_sequence << endl;
      if(badchars.size()>0){
	cerr << '\n';
	for(uint32 i=0; i<badchars.size(); i++){
	  cerr << "-- 2 Illegal character (" << badchars[i] << ": " << std::hex << static_cast<uint16>(badchars[i]) << std::dec << ") in fasta sequence at file byte position " << badpos[i] << endl;
	}
	if(badchars.size()==100){
	  cerr << "\nThere may be more errors like the above, but stopping reporting here.\n";
	}
	cerr << "This happened in sequence: " << FA_fastaseqname << "\nPlease fix your file.\n";

	MIRANOTIFY(Notify::FATAL,"The sequence " << FA_sequence << " in file " << FA_fastaseqname << " showed unrecoverable errors while trying to load it (see also log above). Is it a valid FASTA sequence? Please double check ... and fix your file if necessary.");
      }

      read_read=true;
    }
    }
    //    cout << inchar;
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::loadQual(const char * qualin)
{
  loadINT(qualin,255);
}

void FASTA::loadINT(const char * qualin, int32 maxvalue)
{
  FUNCSTART("void FASTA::loadQual(const char * qualin)");

  std::ifstream qualfin(qualin, std::ios::in|std::ios::ate);
  if(!qualfin){
    MIRANOTIFY(Notify::WARNING, "File not found: " << qualin);
  }

  if(!qualfin.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << qualin);
  }
  qualfin.seekg(0, std::ios::beg);

  loadNextINTSeq(qualfin,maxvalue);

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::loadNextINTSeq(std::ifstream & fin, int32 maxvalue)
{
  FUNCSTART("void FASTA::loadNextINTSeq(std::ifstream & fin, int32 maxvalue)");

  FA_qualseqname.clear();
  FA_intvalues.clear();

  bool read_read=false;

  bool nl=true;
  bool inseq=false;
  bool isnegative=false;
  char inchar=' ';

  while(!(fin.get(inchar)).eof() && !read_read) {
    switch(inchar) {
    case ' ' :
    case '\t' :
    case '\r' : {      // that's from MS-DOS Files, where \n is \r\n
      isnegative=false;
      break;
    }
    case '\n' : {
      nl=true;
      isnegative=false;
      continue;
    }
    case '-' : {
      isnegative=true;
      break;
    }
    case '>' : {
      if(inseq==true) {
	read_read=true;
	fin.unget();
	// FIXME: why? where does this char get eaten?
	fin.unget();
	break;
      }
      if(nl==true) {
	while(!(fin.get(inchar)).eof()) {
	  if(inchar!='\n' && inchar != '\r'){
	    FA_qualseqname+=inchar;
	  } else {
	    fin.unget();
	    break;
	  }
	}
	if(FA_qualseqname.size() == 0) {
	  MIRANOTIFY(Notify::FATAL,"Missing name of fasta sequence in quality file at byte position " << fin.tellg());
	}
      	//cout << "Name: " << FA_qualseqname << endl;

	adjustNameNCBIHack(FA_qualseqname, nullptr);

	nl=true;
	inseq=true;
      } else {
	MIRANOTIFY(Notify::FATAL,"Illegal character (" << inchar << ": " << std::hex << static_cast<uint16>(inchar) << std::dec << ") in fasta sequence name in integer value file at byte position " << fin.tellg());
      }
      break;
    }
    default : {
      if(inseq==false || nl==false) {
	MIRANOTIFY(Notify::FATAL,"Illegal character (" << inchar << ": " << std::hex << static_cast<uint16>(inchar) << std::dec << ") at begin of fasta integer value sequence in file at byte position " << fin.tellg());
      }
      fin.unget();
      std::string tmp;
      while(!fin.eof()) {
	tmp.clear();
	fin >> tmp;
	if(tmp.size()==0) break;
	if(tmp[0]=='>') break;
	if(tmp[0]=='-') {
	  isnegative=true;
	  tmp=tmp.substr(1);
	}
	if(tmp[0]<'0' || tmp[0]>'9') {
	  MIRANOTIFY(Notify::FATAL,"Illegal character (" <<  tmp[0] << ": " << std::hex << static_cast<uint16>(tmp[0]) << std::dec << ") in integer value sequence in file at byte position " << (static_cast<uint32>(fin.tellg()))-tmp.size());
      	}
	int32 thequal=static_cast<int32>(atoi(tmp.c_str()));
	if(thequal>maxvalue){
	  MIRANOTIFY(Notify::FATAL,"Illegal value " << thequal << " (>" << maxvalue << ") in sequence " << FA_qualseqname << " at byte position in file " << fin.tellg() << "\nPlease fix your file.");
      	}
	if(isnegative) thequal=-thequal;
	isnegative=false;
	FA_intvalues.push_back(thequal);
      }
      if(!fin.eof() && tmp.size()) {
	//cout << "Putting back:\n";
	for(size_t tmpi=tmp.size()-1; tmpi!=0; tmpi--){
	  fin.unget();
	  //cout << tmp[tmpi];
	}
	fin.unget();
	// FIXME: why? where does this char get eaten?
	fin.unget();
      }
      read_read=true;
    }
    }
    //    cout << inchar;
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::testIfSeqAndQualMatch()
{
  FUNCSTART("void FASTA::testIfSeqAndQualMatch()");
  if(FA_fastaseqname!=FA_qualseqname) {
    MIRANOTIFY(Notify::FATAL,"Name of read in fasta file (" << FA_fastaseqname << ") and in quality file (" << FA_qualseqname << ") do not match.");
  }

  if(FA_sequence.size()!=FA_intvalues.size()){
    MIRANOTIFY(Notify::FATAL,"Read " << FA_fastaseqname << " has " << FA_sequence.size() << " bases in fasta file, but " << FA_intvalues.size() << " quality values. Cannot be.");
  }
  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool FASTA::testIfEmpty()
{
  if(FA_fastaseqname.size()) return false;
  if(FA_qualseqname.size()) return false;
  if(FA_sequence.size()) return false;
  if(FA_intvalues.size()) return false;
  return true;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::loadNext(std::ifstream & fastafin, std::ifstream & qualin)
{
  FUNCSTART("void FASTA::loadNext(std::ifstream & fastafin, std::ifstream & qualin)");

  try {
    loadNextSeq(fastafin);
  }
  catch(Notify n){
    loadNextINTSeq(qualin,255);
    throw n;
  }
  loadNextINTSeq(qualin,255);

  testIfSeqAndQualMatch();

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::load(const char * fastain)
{
  FUNCSTART("void FASTA::load(const char * fastain)");

  std::ifstream fin(fastain, std::ios::in|std::ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::WARNING, "File not found: " << fastain);
  }

  if(!fin.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << fastain);
  }
  fin.seekg(0, std::ios::beg);

  loadNextSeq(fin);

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::load(const char * fastain, const char * qualin)
{
  FUNCSTART("void FASTA::load(const char * fastain, const char * qualin)");

  std::ifstream fin1(fastain, std::ios::in|std::ios::ate);
  if(!fin1){
    MIRANOTIFY(Notify::WARNING, "File not found: " << fastain);
  }
  if(!fin1.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << fastain);
  }
  fin1.seekg(0, std::ios::beg);

  std::ifstream fin2(qualin, std::ios::in|std::ios::ate);
  if(!fin2){
    MIRANOTIFY(Notify::WARNING, "File not found: " << qualin);
  }
  if(!fin2.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << qualin);
  }
  fin2.seekg(0, std::ios::beg);

  loadNext(fin1,fin2);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::dumpSequence(std::ostream & fout)
{
  FUNCSTART("void FASTA::dumpSequence(ofstream & fout)");
  if(FA_fastaseqname.size()){
    fout << ">" << FA_fastaseqname;
    if(!FA_comment.empty()){
      fout << FA_comment;
    }
    for(uint32 i=0; i<FA_sequence.size(); i++){
      if(i%60==0) fout << "\n";
      fout << FA_sequence[i];
    }
    fout << endl;
  }
  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void FASTA::dumpQuality(std::ostream & fout)
{
  FUNCSTART("void FASTA::dumpQuality(ofstream & fout)");
  if(FA_qualseqname.size()){
    fout << ">" << FA_qualseqname;
    for(uint32 i=0; i<FA_intvalues.size(); i++){
      if(i%25==0) fout << "\n";
      fout << FA_intvalues[i] << " ";
    }
    fout << endl;
  }
  FUNCEND();
}

//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//FASTA::FASTA(FASTA const &other)
//{
//  FUNCSTART("FASTA::FASTA(FASTA const &other)");
//
//  FA_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//FASTA const & FASTA::operator=(FASTA const & other)
//{
//  FUNCSTART("FASTA const & FASTA::operator=(FASTA const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//std::ostream & operator<<(std::ostream &ostr, FASTA const &fas)
//{
//  FUNCSTART("friend std::ostream & FASTA::operator<<(std::ostream &ostr, const  &fas)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}
