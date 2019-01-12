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

//TODO
// in SCF an Kommentare Hinweis auf remix setzen


#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <string>

#include "io/generalio.H"
#include "io/fasta.H"
#include "io/phd.H"
#include "io/scf.H"

#include "util/misc.H"
#include "mira/read.H"


#include "version.H"


using namespace std;


void usage_header()
{
  cout << "scftool V1.0 (MIRALIB version " MIRALIBVERSION ")\n";
  cout << "Written by Bastien Chevreux (bach@chevreux.org)\n\n";
}

void usageTool()
{
  usage_header();
  cout << "Provides a set of tools useful when working with SCF trace files.\n\n";
  cout << "Usage:\n";
  cout << "  scftool <toolname> <tool parameters>\n\n";
  cout << "Available tools:\n";
  cout << "\tconvert\tConverts SCF file(s) to other formats.\n";
  cout << "\tcut\tCuts a given range of a SCF file into a new SCF file.\n";
  cout << "\tremix\tCombines trace information of a SCF file with new bases, qualities and\n\t\tpeak values.\n";
  cout << "\nTo get help on a specific tool, type 'scftool <toolname>'. E.g.: scftool remix\n\n";
}

void usageRemix()
{
  usage_header();
  cout << "scftool remix\n\n";
  cout << "Combines trace information of a SCF file with new bases, qualities and peak\nvalues (either in FASTA or PHD format) into a new SCF file.\n\n";
  cout << "Usage (when using PHD files as input for bases, quals and peaks):\n";
  cout << "  scftool remix scf_infile phd_infile scf_outfile\n";
  cout << "Usage (when using FASTA files as input for bases, quals and peaks):\n";
  cout << "  scftool remix scf_infile bases_infile quals_infile peaks_infile scf_outfile\n";
  cout << "\n";
}

void usageCut()
{
  usage_header();
  cout << "scftool cut\n\n";
  cout << "Cuts the given range out of an SCF file and makes a new SCF out of it.\n\n";
  cout << "Usage:\n";
  cout << "  scftool cut infile lower_base_bound upper_base_bound outfile\n";
  cout << "\n\tbounds: array index is 0! lower bound including, upper bound excluding\n";
}

void usageConvert()
{
  usage_header();
  cout << "scftool convert\n\n\
Extracts the information from SCFs given in a file of filenames and converts \n\
them to other formats.\n\
\n\
Usage:\tscftool convert [-a|-c|-e|-f|-p|-s|-t] [-i] in out\n\
\n\
\t-i\tswitch\ttreat <in> as single file name instead\n\
\t\t\t of a file of filenames\n\
\t-e\tswitch\tconvert to EXP format \n\
\t\t\t EXP files will be written to 'out' directory\n\
\t-f\tswitch\tconvert to FASTA format (default)\n\
\t\t\t Sequence will be written to 'out', quality\n\
\t\t\t values to 'out'.qual, peak index values to\n\
\t\t\t 'out'.peak\n\
\t-a\tswitch\tconvert to ACE format\n\
\t-c\tswitch\tconvert to CAF format\n\
\t-p\tswitch\tconvert to PHD format (not implemented yet)\n\
\t-s\tswitch\tconvert to short readable text format\n\
\t-t\tswitch\tconvert to readable text format\n\
\n\
The file of filenames may contain comment lines (starting with a # sign).\n\
The filenames in the file of filenames may not contain spaces.\n\
";
  exit(2);
}

int mainRemix(int argc, char **argv)
{
  FUNCSTART("mainRemix(int argc, char **argv)");

  string emsg;
  if(argc==5) {
    // remix with PHD

    try{
      SCF thescf;
      PHD thephd;

      if(thescf.load(argv[2],emsg)>0){
	MIRANOTIFY(Notify::FATAL,emsg);
      }

      thephd.load(argv[3]);

      const string & seq=thephd.getSequence();
      thescf.setBases(seq, thephd.getQualities(), thephd.getPeakIndex() );
      const SCF_Comments * comments =thescf.getComments();
      if(comments != nullptr) {
	uint32 comlen=strlen(comments);
	if(comments[comlen-1] != '\n') thescf.addComment("\n");
      }
      thescf.addComment("TPSW=scftool remix\n");

      thescf.save(argv[4]);
    }
    catch(Notify n){
      n.handleError("mainremix PHD");
      exit(1);
    }
    return 0;

  } else if(argc!=7) {
    if(argc>2) cout << "ERROR: wrong number of parameters for tool 'remix'!\n\n";
    usageRemix();
    exit(1);
  }

  // remix with FASTA

  try{
    SCF thescf;

    if(thescf.load(argv[2],emsg)>0){
      MIRANOTIFY(Notify::FATAL,emsg);
    }

    FASTA thefasta;
    thefasta.load(argv[3], argv[4]);

    string bases=thefasta.getSequence();
    vector<uint8> quals=thefasta.getQualities();
    if(thefasta.getSeqName() != thefasta.getQualName()){
      cerr << "Name of sequences in fasta files for bases and for qualities do not match!\nAborting." << endl;
      exit(1);
    }

    thefasta.loadINT(argv[5],1000000);
    vector<uint32> peaks;
    {
      const vector<int32> & lp=thefasta.getINTValues();
      peaks.reserve(lp.size());
      for(uint32 i=0; i<lp.size(); i++){
	peaks.push_back(lp[i]);
      }
    }
    if(thefasta.getSeqName() != thefasta.getQualName()){
      cerr << "Name of sequences in fasta files for bases and for peaks do not match!\nAborting." << endl;
      exit(1);
    }

    thescf.setBases(bases, quals, peaks );

    const SCF_Comments * comments =thescf.getComments();
    if(comments != nullptr) {
      uint32 comlen=strlen(comments);
      if(comments[comlen-1] != '\n') thescf.addComment("\n");
    }
    thescf.addComment("TPSW=scftool remix\n");
    thescf.save(argv[6]);
  }
  catch(Notify n){
    n.handleError("mainRemix FASTA");
    exit(1);
  }

  return 0;
}

int mainCut(int argc, char **argv)
{
  FUNCSTART("int mainCut(int argc, char **argv)");
  if(argc!=6) {
    if(argc>2) cout << "ERROR: wrong number of parameters for tool 'cut'!\n\n";
    usageCut();
    return 1;
  }

  if(atol(argv[3])>= atol(argv[4])){
    cout << "ERROR: lowerbound >= upperbound in parameters of tool 'cut'!\n\n";
    usageCut();
    return 1;
  }

  uint32 low=atol(argv[3]);
  uint32 high=atol(argv[4]);

  try{
    SCF thescf;
    string emsg;

    if(thescf.load(argv[2],emsg)>0){
      MIRANOTIFY(Notify::FATAL,emsg);
    }
    thescf.cutBases(low,high);
    thescf.save(argv[5]);
  }
  catch(Notify n){
    n.handleError("mainCut");
    exit(1);
  }


  return 0;
}





int mainConvert(int argc, char **argv)
{
  FUNCSTART("int mainConvert(int argc, char **argv)");

  int c;
  extern char *optarg;
  extern int optind;

  char opt_targetformat='f';
  bool opt_iswitch=false;

  while (1){
    c = getopt(argc, argv, "acefptsih");
    if(c == -1) break;

    switch (c) {
    case 'i': {
      opt_iswitch=true;
      break;
    }
    case 'p': {
      cout << argv[0] << " " << argv[1] << ": " << "p switch not implemented!\n";
      usageConvert();
    }
    case 'a':
    case 'c':
    case 'e':
    case 'f':
    case 's':
    case 't': {
      opt_targetformat=c;
      break;
    }
    case 'h': {
      usageConvert();
      exit(0);
    }
    case '?': {
      usageConvert();
      exit(1);
    }
    }
  }

  if(argc-optind < 2) {
    usageConvert();
  }

  if(argc-optind < 3) {
    cout << argv[0] << " " << argv[1]  << ": " << "Missing either <in> or <out>  file as arguments!\n";
    usageConvert();
  }

  if(argc-optind > 3) {
    cout << argv[0] << " " << argv[1] << ": " << "Whoops, found more than <in> and <out> as arguments left on the command line!\nPerhaps you forgot some \"\" around strings that contain spaces or other special\ncharacters?\n";
    cout << "Unparsed command line: ";
    for(optind+=1;optind<argc;optind++) cout <<argv[optind] << " ";
    cout << endl;
    usageConvert();
  }

  optind++;
  string fofninfile=argv[optind++];
  string outfile=argv[optind];

  //cout << "fofn: " << fofninfile << endl;
  //cout << "out; " << outfile << endl;
  //exit(0);

  try{

    vector<string> scfnames;
    if(opt_iswitch) {
      scfnames.push_back(fofninfile);
    }else{
      // Load the file of filenames
      ifstream fin;
      fin.open(fofninfile.c_str(), ios::in|ios::ate);
      if(!fin){
	MIRANOTIFY(Notify::FATAL, "File not found: " << fofninfile);
      }

      uint32 len_fofn=fin.tellg();
      if(len_fofn==1){
	MIRANOTIFY(Notify::FATAL, "Zero length file: "<<fofninfile);
      }
      fin.seekg(0, ios::beg);

      string filename, dummy;
      while(GeneralIO::readKeyValue(fin, filename, dummy)){
	scfnames.push_back(filename);
      }
      fin.close();
    }

    ofstream fout1;
    ofstream fout2;
    ofstream fout3;

    if(opt_targetformat=='e') {

      string system_rmdir = (string) "rm -rf "+outfile;
      string system_mkdir = (string) "mkdir "+outfile;

      if(system(system_rmdir.c_str())) {
	MIRANOTIFY(Notify::FATAL, "Could not delete old directory: " << outfile);
      }
      if(system(system_mkdir.c_str())) {
	MIRANOTIFY(Notify::FATAL, "Could not make new directory:"  << outfile);
      }
    } else {
      fout1.open(outfile.c_str(), ios::out);
      if(!fout1){
	MIRANOTIFY(Notify::FATAL, "Could not open file for saving: " << outfile);
      }
      if(opt_targetformat=='f') {
	string of=outfile+".qual";
	fout2.open(of.c_str(), ios::out);
	if(!fout2){
	  MIRANOTIFY(Notify::FATAL, "Could not open file for saving: " << of);
	}

	of=outfile+".peak";
	fout3.open(of.c_str(), ios::out);
	if(!fout3){
	  MIRANOTIFY(Notify::FATAL, "Could not open file for saving: " << of);
	}
      }
    }


    vector<string>::const_iterator I=scfnames.begin();
    Read r;
    while(I!=scfnames.end()){
      r.discard();
      try {
	r.loadDataFromSCF(*I);
      }
      catch(Notify n) {
	I++;
	n.setGravity(Notify::WARNING);
	n.handleError("main - scf loader");
	continue;
      }

      switch(opt_targetformat) {
      case 'f': {
	Read::setCoutType(Read::AS_FASTA);
	fout1 << r;
	Read::setCoutType(Read::AS_FASTAQUAL);
	fout2 << r;

	fout3 << dec;
	fout3 << ">" << r.getName() << endl;

	SCF thescf;
	string emsg;
	if(thescf.load(I->c_str(),emsg)>0){
	  MIRANOTIFY(Notify::FATAL,emsg);
	}

	uint32 cpl=0;
	uint32 numbases=thescf.getNumBases();
	for(uint32 i=0; i<numbases; i++){
	  fout3 << thescf.getPeakIndex(i) << " ";
	  if(cpl++==14){
	    cpl=0;
	    fout3 << "\n";
	  }
	}
	if(cpl!=0) fout3 << "\n";

	break;
      }
      case 'e': {
	Read::setCoutType(Read::AS_GAP4DA);
	string path=outfile+"/"+r.getName()+".exp";
	fout1.open(path.c_str(), ios::out);
	if(!fout1){
	  MIRANOTIFY(Notify::FATAL, "Could not open file for saving: " << path);
	}
	fout1 << r;
	fout1.close();
	break;
      }
      case 'c': {
	Read::setCoutType(Read::AS_CAF);
	fout1 << r;
	break;
      }
      case 'a': {
	Read::setCoutType(Read::AS_ACE);
	fout1 << r;
	break;
      }
      case 's': {
	Read::setCoutType(Read::AS_TEXTSHORT);
	fout1 << r;
	break;
      }
      case 't': {
	Read::setCoutType(Read::AS_TEXT);
	fout1 << r;
	break;
      }
      }
      I++;
    }

    if(opt_targetformat=='e') {
    } else {
      fout1.close();
      if(opt_targetformat=='f') {
	fout2.close();
	fout3.close();
      }
    }
  }
  catch(Notify n){
    n.handleError("mainConvert");
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

  if(toolchosen=="remix"){
    mainRemix(argc,argv);
  }else if(toolchosen=="cut"){
    mainCut(argc,argv);
  }else if(toolchosen=="convert"){
    mainConvert(argc,argv);
    //}else if(toolchosen=="bla"){
  }else{
    cout << "ERROR: unknown tool '" << toolchosen << "'\n\n";
    usageTool();
    exit(1);
  }


  return 0;
}
