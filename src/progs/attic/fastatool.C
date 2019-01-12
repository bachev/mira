/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2005 and later by Bastien Chevreux
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


#include <unistd.h>

#include <fstream>
#include <iostream>
#include "io/fasta.H"

#include "version.H"


using namespace std;

#if __GNUC__ >= 3
#define NO_NOCREATE
#endif


void usage_header()
{
  cout << "fastatool V1.0 (MIRALIB version " MIRALIBVERSION ")\n";
  cout << "Written by Bastien Chevreux (bach@chevreux.org)\n\n";
}

void usageTool()
{
  usage_header();
  cout << "Provides a set of tools useful when working with FASTA files.\n\n";
  cout << "Usage:\n";
  cout << "  fastatool <toolname> <tool parameters>\n\n";
  cout << "Available tools:\n";
  cout << "\tclip\tClip left and right parts of a FASTA sequence.\n";
  cout << "\tsanitize\tDeletes invalid sequences from FASTA files.\n";
  cout << "\nTo get help on a specific tool, type 'fastatool <toolname>'. E.g.: fastatool cut\n\n";
}

void usageSanitize()
{
  usage_header();
  cout << "fastatool sanitize\n\n\
Deletes sequences with non-IUPAC bases or empty sequences from a FASTA file,\n\
writes the 'surviving' sequences to a new file.\n\
If a quality file is given, a cleaned version of that is also written.\n\
The sequences in the quality files (if given) MUST be in the same order than\n\
in the fasta file.\n\
\n\
Usage:\n\
\tsanitizeFASTA fastainfile fastaoutfile [fastaqualin fastaqualoutfile]\n";
}


void usageClip()
{
  usage_header();
  cout << "fastatool clip\n\n\
Clips bases on the left and right of a FASTA sequence, writes result to STDOUT.\n\
\n\
Usage:\n\
\tclipFASTA -l leftclip -r rightclip fastainfile\n";
}


int mainSanitize (int argc, char **argv)
{
  if(argc != 4 && argc != 6) {
    if(argc>2) cout << "ERROR: Wrong number of arguments for fastatool sanitize!\n";
    usageSanitize();
    exit(1);
  }

  try{
    ifstream fin;
    ofstream fout;

#ifdef NO_NOCREATE
    fin.open(argv[2], ios::in|ios::ate);
#else
    fin.open(argv[2], ios::in|ios::ate|ios::nocreate);
#endif
    if(!fin){
      throw Notify(Notify::FATAL, "mainSanitize", (static_cast<std::string>("File not found: ")+argv[2]).c_str());
    }
    if(!fin.tellg()){
      throw Notify(Notify::FATAL, "mainSanitize", (static_cast<std::string>("Zero length file: ")+argv[2]).c_str());
    }
    fin.seekg(0, ios::beg);
    fout.open(argv[3], ios::out);

    ifstream qin;
    ofstream qout;
    if(argc == 6) {
#ifdef NO_NOCREATE
      qin.open(argv[4], ios::in|ios::ate);
#else
      qin.open(argv[4], ios::in|ios::ate|ios::nocreate);
#endif
      if(!qin){
	throw Notify(Notify::FATAL, "mainSanitize", (static_cast<std::string>("File not found: ")+argv[4]).c_str());
      }
      if(!qin.tellg()){
	throw Notify(Notify::FATAL, "mainSanitize", (static_cast<std::string>("Zero length file: ")+argv[4]).c_str());
      }
      qin.seekg(0, ios::beg);

      qout.open(argv[5], ios::out);
    }

    FASTA thefasta;
    while(1){
      try{
	if(argc==4){
	  thefasta.loadNextSeq(fin);
	}else{
	  thefasta.loadNext(fin,qin);
	}
	if(thefasta.testIfEmpty()) {
	  // no more sequences.
	  break;
	}
	if(thefasta.getSequence().size()>0) {
	  //      REP_thepool.back().setFileNames(thefasta.getSeqName().c_str());
	  const char * seq= thefasta.getSequence().c_str();
	  const vector<uint8> quals=thefasta.getQualities();
	  uint32 slen=thefasta.getSequence().size();
	  fout << ">" << thefasta.getSeqName().c_str();
	  for(uint32 i=0; i<slen; i++, seq++){
	    if(i%60==0){
	      fout << endl;
	    }
	    fout << *seq;
	  }
	  fout << "\n";
	  if(argc==6){
	    qout << ">" << thefasta.getSeqName().c_str();
	    for(uint32 i=0; i<slen; i++, seq++){
	      if(i%20==0){
		qout << endl;
	      }
	      qout << (uint16) quals[i] << " ";
	    }
	    qout << "\n";
	  }
	} else{
	  cerr << "Thrown out " << thefasta.getSeqName().c_str() << ": no bases" << endl;
	}

      }
      catch(Notify n){
	n.setGravity(Notify::WARNING);
	n.handleError("mainSanitize");
	cerr << "Read " << thefasta.getSeqName().c_str() << " thrown out." << endl;
      }
    }

    fin.close();
  }
  catch(Notify n) {
    n.handleError("mainSanitize");
    exit(1);
  }

  return 0;
}


int mainClip(int argc, char **argv)
{
  if(argc==2) {
    usageClip();
    exit(0);
  }

  int c;
  extern char *optarg;
  extern int optind;

  bool lflag=false;
  bool rflag=false;

  int32 GLO_cl=0;
  int32 GLO_cr=0;

  while (1){
    c = getopt(argc, argv, "l:r:");
    if(c == -1) break;

    switch (c) {
    case 'l': {
      int32 cl=atoi(optarg);
      if(cl < 0) {
	cerr << "ERROR: -l clip left must be >=0 in fastatool clip\n";
	usageClip();
	exit(1);
      }
      GLO_cl=cl;
      lflag=true;
      break;
    }
    case 'r': {
      int32 cr=atoi(optarg);
      if(cr<0) {
	cerr << "ERROR: -r clip right must be >=0 in fastatool clip\n";
	usageClip();
	exit(1);
      }
      GLO_cr=cr;
      rflag=true;
      break;
    }
    case '?':
      usageClip();
      exit(0);
    }
  }

  if(!lflag) {
    cerr << "ERROR: missing -l flag in fastatool clip\n";
    usageClip();
    exit(1);
  }
  if(!rflag) {
    cerr << "ERROR: missing -r flag in fastatool clip\n";
    usageClip();
    exit(1);
  }


  optind++;

  try{
    ifstream fin;
    ofstream fout;

#ifdef NO_NOCREATE
    fin.open(argv[optind], ios::in|ios::ate);
#else
    fin.open(argv[optind], ios::in|ios::ate|ios::nocreate);
#endif
    if(!fin){
      throw Notify(Notify::FATAL, "main", (static_cast<std::string>("File not found: ")+argv[optind]).c_str());
    }
    if(!fin.tellg()){
      throw Notify(Notify::FATAL, "main", (static_cast<std::string>("Zero length file: ")+argv[optind]).c_str());
    }
    fin.seekg(0, ios::beg);

    FASTA thefasta;
    while(1){
      try{
	thefasta.loadNextSeq(fin);

	if(thefasta.testIfEmpty()) {
	  // no more sequences.
	  break;
	}
	cout << ">" << thefasta.getSeqName().c_str();
	if(thefasta.getSequence().size()>0) {

	  const char * seq= thefasta.getSequence().c_str();
	  uint32 slen=thefasta.getSequence().size();

	  uint32 charsout=0;
	  for(uint32 i=0; i<slen; i++, seq++){
	    if(i<GLO_cl) continue;
	    if(i>=GLO_cr) break;
	    if(charsout%60==0){
	      cout << "\n";
	    }
	    cout << *seq;
	    charsout++;
	  }
	}
	cout << "\n";

      }
      catch(Notify n){
	n.setGravity(Notify::WARNING);
	n.handleError("main");
	cerr << "Read " << thefasta.getSeqName().c_str() << " thrown out." << endl;
	exit(10);
      }
    }

    fin.close();
  }
  catch(Notify n) {
    n.handleError("main");
    exit(1);
  }

  return 0;
}

int main(int argc, char **argv)
{
  if(argc<2){
    usageTool();
    exit(1);
  }

  string toolchosen=argv[1];
  std::transform(toolchosen.begin(),
		 toolchosen.end(),
		 toolchosen.begin(),
		 (int(*)(int))std::tolower); // now, that's what I call ugly

  if(toolchosen=="sanitize"){
    mainSanitize(argc,argv);
  }else if(toolchosen=="clip"){
    mainClip(argc,argv);
    //}else if(toolchosen=="bla"){
  }else{
    cout << "ERROR: unknown tool '" << toolchosen << "'\n\n";
    usageTool();
    exit(1);
  }


  return 0;
}
