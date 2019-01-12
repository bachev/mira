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

#include <iostream>
#include <iterator>
#include <string>


#include "modules/mod_memestim.H"
#include "modules/misc.H"


#include "mira/assembly.H"
#include "util/stlimprove.H"

#include "version.H"


using std::cout;
using std::cin;
using std::cerr;
using std::endl;


void mme_askChar(const std::string & question, const std::string & possibilities, char & answer, const char defchar)
{
  bool doloop=true;
  while(doloop){
    cout << question << ' ';
    if(!possibilities.empty()){
      cout <<"(";
      mstd::copy(possibilities,
		 std::ostream_iterator<char> (cout, "/"));
      cout << ") ";
    }
    if(defchar!=0){
      cout << "[" << defchar << "] ";
    }
    std::string input;
    getline(cin,input);

    // empty input, try to get defchar from possibilities if it exists
    if(input.empty() && defchar!=0) input=defchar;
    input.resize(1);
    for(uint32 i=0; i<possibilities.size(); i++){
      if(input[0] == possibilities[i]) {
	doloop=false;
	answer=input[0];
	break;
      }
    }
  }
  cout << answer << endl;
}

void mme_askDoubleNP(const std::string & question, double & answer, const std::string & defd)
{
  bool doloop=true;
  while(doloop){
    cout << question << ' ';
    if(!defd.empty()){
      cout << "[" << defd << "] ";
    }
    std::string input;
    getline(cin,input);

    // empty input, try to get def from possibilities if it exists
    if(input.empty() && !defd.empty()) {
      input=defd;
    }
    char * pend;
    answer=strtod(input.c_str(),&pend);
    doloop=false;

    // try to parse kilo, mega, giga
    while(*pend != 0 && isspace(*pend)) pend++;
    switch(toupper(*pend)){
    case 0 : break;
    case 'K' :{
      answer*=1000;
      break;
    }
    case 'M' :{
      answer*=1000000;
      break;
    }
    case 'G' :{
      answer*=1000000000;
      break;
    }
    default : {
      cout << "Please only use k, m, g as modifiers.\n";
      doloop=true;
    }
    }
  }
}

void mme_askDouble(const std::string & question, double & answer, const std::string & defd)
{
  mme_askDoubleNP(question, answer, defd);
  cout << answer << endl;
}

void mme_askInt(const std::string & question, int64 & answer, const std::string & defint)
{
  double tmp;
  mme_askDoubleNP(question, tmp, defint);
  answer=static_cast<int64>(tmp);
  cout << answer << endl;
}


void miraMemEstimate(int argc, char ** argv)
{
  int c;
  extern char *optarg;
  extern int optind;

  while (1){
    c = getopt(argc, argv, "v");
    if(c == -1) break;

    switch (c) {
    case 'v':
      cout << miraversion << endl;
      exit(0);
    default : {}
    }
  }


  cout << "This is MIRA " << miraversion << ".\n\n";

  cout << "Please cite: Chevreux, B., Wetter, T. and Suhai, S. (1999), Genome Sequence\nAssembly Using Trace Signals and Additional Sequence Information.\nComputer Science and Biology: Proceedings of the German Conference on\nBioinformatics (GCB) 99, pp. 45-56.\n\n";
  dumpStdMsg();

  cout << "\n\nmiraMEM helps you to estimate the memory needed to assemble a project.\n"
    "Please answer the questions below.\n\n"
    "Defaults are give in square brackets and chosen if you just press return.\n"
    "Hint: you can add k/m/g modifiers to your numbers to say kilo, mega or giga.\n\n";

  char yesno;
  char ptype=' ';
  char denovomapping;
  int64 seqsize=0;

  int64 numsanreads=0;
  int64 num454gs20reads=0;
  int64 num454flxreads=0;
  int64 num454titaniumreads=0;
  int64 numsxareads=0;
  int64 avgsxalen=0;
  int64 numpbsreads=0;
  int64 avgpbslen=0;
  int64 largestcontigexpected=0;

  // computed
  int64 totalexpectedbases=0;
  int64 totalreads=0;
  int64 readsinlargestcontig=0;
  int64 readbasesinlargestcontig=0;

  mme_askChar("Is it a genome or transcript (EST/tag/etc.) project?",
	      "ge",
	      ptype,
	      'g');
  if(ptype=='g'){
    mme_askInt("Size of genome?",
	       seqsize,
	       "4.5m");
    if(seqsize<100) {
      cout << "A sequence size of less than 100 bases is pretty improbable.\n"
	   << "Did you forget a modifier (k, m or g) to the number you gave?\n";
      exit(10);
    }
    largestcontigexpected=seqsize;
    if(largestcontigexpected>30*1000*1000){
      cout << "Looks like a larger eukaryote, guessing largest chromosome size: 30m\nChange if needed!\n";
      largestcontigexpected=30*1000*1000;
    }

    {
      std::string tmplc;
      std::ostringstream ostr;
      ostr << largestcontigexpected;
      tmplc=ostr.str();
      mme_askInt("Size of largest chromosome?",
		 largestcontigexpected,
		 tmplc);
    }

    mme_askChar("Is it a denovo or mapping assembly?",
		"dm",
		denovomapping,
		'd');
  }


  mme_askInt("Number of Sanger reads?",
	     numsanreads,
	     "0");
  mme_askChar("Are there 454 reads?",
	      "yn",
	      yesno,
	      'n');
  if(yesno=='y'){
    mme_askInt("Number of 454 GS20 reads?",
	       num454gs20reads,
	       "0");
    mme_askInt("Number of 454 FLX reads?",
	       num454flxreads,
	       "0");
    mme_askInt("Number of 454 Titanium reads?",
	       num454titaniumreads,
	       "0");
  }
  mme_askChar("Are there PacBio reads?",
	      "yn",
	      yesno,
	      'n');
  if(yesno=='y'){
    mme_askInt("Number of PacBio reads?",
	       numpbsreads,
	       "0");
    mme_askInt("Average PacBio length?",
	       avgpbslen,
	       "1100");
  }
  mme_askChar("Are there Solexa reads?",
	      "yn",
	      yesno,
	      'n');
  if(yesno=='y'){
    mme_askInt("Number of Solexa reads?",
	       numsxareads,
	       "0");
    mme_askInt("Average Solexe length?",
	       avgsxalen,
	       "75");
  }

  totalexpectedbases=numsanreads*1000;
  totalexpectedbases+=num454gs20reads*120;
  totalexpectedbases+=num454flxreads*260;
  totalexpectedbases+=num454titaniumreads*460;
  totalexpectedbases+=numsxareads*avgsxalen;
  totalexpectedbases+=numpbsreads*avgpbslen;

  totalreads=numsanreads;
  totalreads+=num454gs20reads;
  totalreads+=num454flxreads;
  totalreads+=num454titaniumreads;
  totalreads+=numsxareads;
  totalreads+=numpbsreads;

  if(ptype=='g'){
    if(denovomapping=='d'){
      readsinlargestcontig=totalreads/2;
      readbasesinlargestcontig=totalexpectedbases/2;
    }else{
      largestcontigexpected=seqsize;
      readsinlargestcontig=totalreads;
      readbasesinlargestcontig=totalexpectedbases;

      // if solexa is mapped, there are less reads due to
      //  coverage equivalent mapping and virtual long reads
      // be conservative, reduce only by 50%
      if(numsxareads>0){
	readsinlargestcontig-=numsxareads/2;
      }
    }
  }else{
    seqsize=50000;
    largestcontigexpected=seqsize;
    readsinlargestcontig=50000;
    readbasesinlargestcontig=readsinlargestcontig*1000; //10k reads times sanger length
  }

  // account for gaps with 454 reads
  if(num454flxreads>0 || num454gs20reads>0){
    largestcontigexpected+=largestcontigexpected/10;
    readbasesinlargestcontig+=readbasesinlargestcontig/10;
  }

  //cout << "totalreads: " << totalreads
  //     << "\nreadsinlargestcontig: " << readsinlargestcontig
  //     << "\ntotalexpectedbases: " << totalexpectedbases
  //     << "\nreadbasesinlargestcontig: " << readbasesinlargestcontig
  //     << endl;

  int64 livereads=totalreads+readsinlargestcontig;
  int64 livebases=totalexpectedbases+readbasesinlargestcontig;

  //cout << "livereads: " << livereads
  //     << "\nlivebases: " << livebases << endl;

  double avgcov=static_cast<double>(totalexpectedbases/seqsize);
  avgcov-=avgcov/8; // in general we have 12% loss of usable data

  int64 numskimhits=static_cast<int64>(avgcov*850000);  // estimate skim hits, very rough

  int64 memneeded=0;

  // what do the reads need?
  memneeded=
    livereads*sizeof(Read)       // class size
    +livereads*200               // additional strings etc.
    +livereads*4*sizeof(tag_t)   // on average 4 tags per read
    +livebases*8;                // sequences, qualities, adjustments, base flags

  // new: solexa reads don't have adjustments
  // yeah, but estimate is already small enough, keep it
  //memneeded-=(numsxareads*avgsxalen) * 2;


  // what does a contig need?
  //  (note: the needs for the reads themselves are already
  //   accounted for in the section above)
  memneeded+=
    //readsinlargestcontig*sizeof(Contig::contigread_t)
    //+readsinlargestcontig*sizeof(Contig::out_order_t)

    readsinlargestcontig*40       // 40 == rough guesstimate for PlacedContigReads
    +totalreads*9                              /* templates, mapping
						 allowedrefids */
    +largestcontigexpected*sizeof(Contig::consensus_counts_t)
    +largestcontigexpected*10;              // adjustments and some reserve

  int64 memforlargetables=0;
  // some more overhead by the assembly class
  memforlargetables+= totalreads*20;

  // get the skim edges accounted
  int64 skimhitsmem=numskimhits*2*sizeof(skimedges_t);
  // since 2.9.40 there's the possibility to cap that memory
  // use default value
  if(skimhitsmem>1024L*1024*1024){
    skimhitsmem=2LL*1024*1024*1024;
    if(numsxareads>0) skimhitsmem*=2;
  }

  // mem needed for temporary skim need
  int64 tmpskim=500*1000*1000;

  memforlargetables+=std::max(skimhitsmem,tmpskim);

  // possible vector leftover clip
  int64 memforpvlc=0;
  {
    // AS_readhitmiss & AS_readhmcovered
    memforpvlc=totalexpectedbases*8;
    // overhead of the structures above
    memforpvlc+=sizeof(std::vector<uint32>)*totalreads*2;

    // AS_count_rhm
    memforpvlc+=totalreads*4;
  }

  // ok, 1MB of additional small things
  int64 memneededfordata=memneeded+(1024*1024);

  // experience shows that not all has been accounted for
  //  and internal mem caching of memory allocators add another
  //  layer of RAM needs
  //
  // add 40% to estimates
  //  but not if whe have mapping with Solexas
  if(denovomapping!='m' && numsxareads==0){
    memneededfordata+=memneededfordata/100*40;
    memforlargetables+=memforlargetables/100*40;
  }

  cout.setf(std::ios::fixed, std::ios::floatfield);
  //cout.setf(std::ios::showpoint);
  cout.precision(1);

  cout << "\n\n************************* Estimates *************************\n\n";
  // last, if it's an EST assembly, there is no seqsize

  if(ptype=='e'){
    cout << "EST assembly, cannot give coverage estimate. Also, estimates"
      "\nmay be way off for pathological cases.\n";

  }else{
    cout << "The contigs will have an average coverage of ~ " << avgcov
	 << " (+/- 10%)"
      "\nEstimates may be way off for pathological cases.\n";
  }

  cout << "\nRAM estimates:"
    "\n" << std::setw(40) << "reads+contigs (unavoidable): ";
  byteToHumanReadableSize(memneededfordata,cout);
  cout << "\n" << std::setw(40) << "large tables (tunable): ";
  byteToHumanReadableSize(memforlargetables,cout);
  cout << "\n" << std::setw(40) << "" << "---------"
    "\n" << std::setw(40) << "total (peak): ";
  byteToHumanReadableSize(memforlargetables+memneededfordata,
			  cout);
  cout << "\n\n" << std::setw(40) << "add if using -CL:pvlc=yes : ";
  byteToHumanReadableSize(memforpvlc,cout);
  if(denovomapping=='m' && numsxareads>0){
    int64 notusingmerge=memneededfordata/100*40;
    cout << "\n" << std::setw(40) << "add if setting -CO:msr=no : ";
    byteToHumanReadableSize(notusingmerge,cout);
  }
  cout << "\n\n"
    "Note that some algorithms might try to grab more memory if"
    "\nthe need arises and the system has enough RAM. The options"
    "\nfor automatic memory management control this:"
    "\n  -AS:amm, -AS:kpmf, -AS:mps"
    "\nFurther switches that might reduce RAM (at cost of run time"
    "\nor accuracy):"
    "\n  -SK:mhim, -SK:mchr (both runtime); -SK:mhpr (accuracy)\n"
    "*************************************************************\n";

}


