/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2002 and later by Bastien Chevreux
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

#include "mira/assembly.H"
#include "version.H"


struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

typedef hash_multimap<const char*, int, hash<const char *>, eqstr> stringhash_t;
typedef stringhash_t::value_type stringhash_entry_t;




void uncover(Read & masked, Read & unmasked, int32 minrunlen=6, int32 minuncoverlen=7, bool uncoverall=false, bool alsorecover=false, bool ignorenbasemismatch=false);

void uncover(Read & masked, Read & unmasked, int32 minrunlen, int32 minuncoverlen, bool uncoverall, bool alsorecover, bool ignorenbasemismatch){

  if(masked.getLenSeq() != unmasked.getLenSeq()){
    cerr << "Mismatch in number of bases for sequence " << masked.getName() << " between masked and unmasked sequence, leaving sequence untouched.\n";
    return;
  }
  {
    for(uint32 i=0; i<masked.getLenSeq(); i++){
      char mc=tolower(masked.getBaseInSequence(i));
      char uc=tolower(unmasked.getBaseInSequence(i));
      if(mc != uc) {
	if(mc != 'x' && uc != 'x'){
	  if(ignorenbasemismatch && (mc != 'n' && uc != 'n')) {
	    cerr << "Mismatch in base number " << i << " for sequence " << masked.getName() << " between masked and unmasked sequence, leaving sequence untouched. (" << mc << " against " << uc << ")\n";
	    return;
	  }
	}
      }
    }
  }

  int32 ml=masked.getLMClipoff()-1;
  uint32 numuncover=0;
  uint32 countA=0;
  uint32 countC=0;
  uint32 countG=0;
  uint32 countT=0;
  //cout << ml << endl;
  while(ml>0 && numuncover <minuncoverlen) {
    switch(tolower(unmasked.getBaseInSequence(ml))) {
    case 'a': {
      countA++;
      break;
    }
    case 'c': {
      countC++;
      break;
    }
    case 'g': {
      countG++;
      break;
    }
    case 't': {
      countT++;
      break;
    }
    default: {}
    }
    ml--;
    numuncover++;
  }
  //cout << "countA: " << countA << endl;
  //cout << "countC: " << countC << endl;
  //cout << "countG: " << countG << endl;
  //cout << "countT: " << countT << endl;

  bool mustuncover=false;
  char runbase=' ';
  if(countA>=minrunlen) {
    mustuncover=true;
    runbase='a';
  }
  if(countT>=minrunlen) {
    mustuncover=true;
    runbase='t';
  }

  if(mustuncover) {
    ml=masked.getLMClipoff()-1;
    for(uint32 i=0; i<numuncover; i++, ml--) {
      // don't care about the quality
      masked.changeBaseInSequence(unmasked.getBaseInSequence(ml), 0, ml);
    }
    if(uncoverall) {
      while(ml>=0) {
	if(tolower(unmasked.getBaseInSequence(ml)) != runbase) break;
	masked.changeBaseInSequence(unmasked.getBaseInSequence(ml), 0, ml);
	ml--;
	numuncover++;
      }
    }
    cout << "Uncovered " << numuncover << " at start of " << masked.getName() << endl;

    if(alsorecover) {
      ml++;
      int32 ml2=-minuncoverlen;
      while(ml<unmasked.getLenSeq()){
	if(tolower(unmasked.getBaseInSequence(ml))!=runbase) break;
	if(ml2>=0) masked.changeBaseInSequence('x', 0, ml-minuncoverlen);
	ml++;ml2++;
      }
      if(ml2>0) {
	cout << "Recovered " << ml2 << " at start of " << masked.getName() << endl;
      }
    }
    masked.setLMClipoff(ml-minuncoverlen);
  }

  int32 mr=masked.getRMClipoff();
  numuncover=0;
  countA=0;
  countC=0;
  countG=0;
  countT=0;
  //cout << mr << endl;
  while(mr<unmasked.getLenSeq() && numuncover <minuncoverlen) {
    switch(tolower(unmasked.getBaseInSequence(mr))) {
    case 'a': {
      countA++;
      break;
    }
    case 'c': {
      countC++;
      break;
    }
    case 'g': {
      countG++;
      break;
    }
    case 't': {
      countT++;
      break;
    }
    default: {}
    }
    mr++;
    numuncover++;
  }
  //cout << "countA: " << countA << endl;
  //cout << "countC: " << countC << endl;
  //cout << "countG: " << countG << endl;
  //cout << "countT: " << countT << endl;

  mustuncover=false;
  runbase=' ';
  if(countA>=minrunlen) {
    mustuncover=true;
    runbase='a';
  }
  if(countT>=minrunlen) {
    mustuncover=true;
    runbase='t';
  }

  if(mustuncover) {
    mr=masked.getRMClipoff();
    for(uint32 i=0; i<numuncover; i++, mr++) {
      // don't care about the quality
      masked.changeBaseInSequence(unmasked.getBaseInSequence(mr), 0, mr);
    }

    if(uncoverall) {
      while(mr<unmasked.getLenSeq()) {
	if(tolower(unmasked.getBaseInSequence(mr)) != runbase) break;
	masked.changeBaseInSequence(unmasked.getBaseInSequence(mr), 0, mr);
	mr++;
	numuncover++;
      }
    }
    cout << "Uncovered " << numuncover << " at end of " << masked.getName() << endl;

    if(alsorecover) {
      mr--;
      int32 mr2=-minuncoverlen;
      while(mr>=0) {
	if(tolower(unmasked.getBaseInSequence(mr))!=runbase) break;
	if(mr2>=0) masked.changeBaseInSequence('x', 0, mr+minuncoverlen);
	mr--;mr2++;
      }
      if(mr2>0) {
	cout << "Recovered " << mr2 << " at end of " << masked.getName() << endl;
      }
    }
    masked.setRMClipoff(mr);
  }

}

void usage()
{
  cout << "uncover_at V1.0.1 (MIRALIB version "MIRALIBVERSION")\n\n";
  cout << "Written by Bastien Chevreux (bach@chevreux.org)\n\n";
  cout << "Uncovers Poly-AT stretches in masked parts of a fasta file,\n\
              writes a masked file with uncovered stretches to outfile." << endl;
  cout << "Usage:\n";
  cout << "\tuncover_at [-i number] [-l number] [-a] [-r] maskedfile unmaskedfile outfile\n\n";
  cout << "Options:\n\n\
\t-i\tinteger\tinitial length of stretch to be a poly-A/poly-T\n\t\t\t\t (defaults to 6)\n\
\t-l\tinteger\tlength of anchor (poly-A/poly-T stretch uncovered)\n\t\t\t\t (defaults to 7)\n\
\t-a\tswitch\tuncover all of the poly-A/poly-T stretch (even if\n\t\t\t\t longer than anchor)\n\
\t-r\tswitch\trecover poly-A/poly-T stretch, only leaving an anchor\n\t\t\t\t uncovered\n\
\t-n\tswitch\tignore mismatches between 'N' and a base\n\
\n\
\t-P\tstring\tParameters for the clipping, syntax like for the\n\
\t\t\t mira(1) assembler. If given, then quality clip and/or\n\
\t\t\t masked bases clips are computed for the reads and set\n\
\t\t\t depending on which options were used.\n\
\t\t\tUseful settings: '-CL:mbc*:qc*'.\n\
\t\t\t E.g. '-CL:mbc=on:mbcmeg=240:qc=on:qcwl=25'\n\
";
}


int main(int argc, char ** argv)
{
  FUNCSTART("int main(int argc, char ** argv)");

  cout << "uncover_at must be reworked and is not available at this time, sorry.\n";
  abort();

  int c;
  extern char *optarg;
  extern int optind;

  string opt_miraoptions="";
  bool opt_doclips=false;

  bool opt_uncoverall=false;
  bool opt_recover=false;
  bool opt_ignorenbasemismatch=false;

  int32 opt_minrunlen=6;
  int32 opt_minuncoverlen=7;

  while (1){
    c = getopt(argc, argv, "anri:l:P:h");
    if(c == -1) break;

    switch (c) {
    case 'P': {
      opt_miraoptions=optarg;
      opt_doclips=true;
      break;
    }
    case 'a': {
      opt_uncoverall=true;
      break;
    }
    case 'n': {
      opt_ignorenbasemismatch=true;
      break;
    }
    case 'r': {
      opt_recover=true;
      break;
    }
    case 'i': {
      opt_minrunlen=atoi(optarg);
      break;
    }
    case 'l': {
      opt_minuncoverlen=atoi(optarg);
      break;
    }
    case 'h':
    case '?': {
      usage();
      exit(0);
    }
    default:
      cerr << "Unknown option " << c << "\n";
      usage();
      exit(10);
    }
  }

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing all files (masked/unmasked/out) as arguments!\n";
    usage();
    exit(1);
  }

  if(argc-optind < 2) {
    cerr << argv[0] << ": " << "Missing two files from (masked/unmasked/out) as arguments!\n";
    usage();
    exit(1);
  }

  if(argc-optind < 3) {
    cerr << argv[0] << ": " << "Missing one file from (masked/unmasked/out) as argument!\n";
    usage();
    exit(1);
  }

  if(argc-optind > 3) {
    cerr << argv[0] << ": " << "Whoops, found more than (masked/unmasked/out) files as arguments left on the command line!\n";
    cerr << "Unparsed command line: ";
    for(;optind<argc;optind++) cerr <<argv[optind] << " ";
    cerr << endl;
    usage();
    exit(1);
  }

  string maskedfile=argv[optind++];
  string unmaskedfile=argv[optind++];
  string outfile=argv[optind];

  //cout << maskedfile << endl;
  //cout << unmaskedfile << endl;
  //cout << outfile << endl;

  //if(infile=="--help"){
  //  usage();
  //  exit(0);
  //}

  MIRAParameters P;
  // "-CL:mbc=1:mbcmeg=240"
  P.parse("-CL:mbc=0:qc=0");
  P.parse(opt_miraoptions.c_str());
  //cout << P;
  //p.getAssemblyParams().as_clip_maskedbase_gapsize=20;

  ReadPool maskedpool(&P), unmaskedpool(&P);

  try{
    maskedpool.loadDataFromFASTA(maskedfile);
    unmaskedpool.loadDataFromFASTA(unmaskedfile);

    assembly_parameters const & as_params= P.getAssemblyParams();

    if(opt_doclips) {
      cout << "Starting clips:";
      if(as_params.as_clip_quality) {
	cout << " quality";
	cout.flush();
	maskedpool.performQualClips(as_params.as_clip_quality_minqual,
				    as_params.as_clip_quality_winlen);
      }
      if(as_params.as_clip_quality && as_params.as_clip_maskedbases) cout << ",";
      if(as_params.as_clip_maskedbases) {
	cout << " masked characters ";
	cout.flush();
	maskedpool.performMaskCharClips(as_params.as_clip_maskedbase_gapsize,
					as_params.as_clip_maskedbase_maxfrontgap,
					as_params.as_clip_maskedbase_maxendgap);
      }
      cout << "done\n";
    }

    uint32 u_numreads=unmaskedpool.size();
    uint32 m_numreads=maskedpool.size();


    vector<string> u_names;
    for(uint32 i=0; i<u_numreads; i++) {
      if(!unmaskedpool.getRead(i).getName().empty()) {
	u_names.push_back(unmaskedpool.getRead(i).getName());
	//cout << unmaskedpool.getRead(i).getName() << endl;
      }
    }

    bool stopprocessing=false;
    stringhash_t M;
    pair<stringhash_t::const_iterator, stringhash_t::const_iterator> p;
    for(uint32 i=0; i< u_names.size(); i++){
      p=M.equal_range(u_names[i].c_str());
      if(p.first!=M.end()){
	cerr << "Sequence " << u_names[i] << " is present more than once in unmasked file " << unmaskedfile << endl;
	stopprocessing=true;
      }else{
	M.insert(stringhash_entry_t(u_names[i].c_str(), i));
      }
    }

    if(stopprocessing) {
      cerr << "Unrecoverable error (see above), stopping the processing.\n";
      exit(10);
    }

    ofstream fout;
    fout.open(outfile.c_str(), ios::out);
    if(!fout){
      throw Notify(Notify::FATAL, THISFUNC, outfile.c_str(), ": Could not open file for saving.");
    }

    //Read::setCoutType(Read::AS_FASTA);
    //Read::setCoutType(Read::AS_MASKEDFASTA);
    Read::setCoutType(Read::AS_FASTA);
    if(opt_doclips) Read::setCoutType(Read::AS_MASKEDMASKFASTA);
    for(uint32 mi=0; mi<m_numreads; mi++){
      p=M.equal_range(maskedpool.getRead(mi).getName().c_str());
      if(p.first!=M.end()){
	maskedpool.getRead(mi).setClipoffsToMaskedChars(25,40,80);
	uncover(maskedpool.getRead(mi),
		unmaskedpool.getRead((p.first)->second),
		opt_minrunlen,
		opt_minuncoverlen,
		opt_uncoverall,
		opt_recover,
		opt_ignorenbasemismatch);
	//Read::setCoutType(Read::AS_TEXTCLIPS);
	//fout << maskedpool.getRead(mi);
	//Read::setCoutType(Read::AS_MASKEDMASKFASTA);
	fout << maskedpool.getRead(mi);
      } else {
	fout << maskedpool.getRead(mi);
      }
    }

    fout.close();
    //thepool.dumpAsFASTA(fout,false,false);
  }
  catch(Notify n){
    n.handleError("main");
  }
  catch(Flow f){
    cerr << "Unexpected exception: Flow()\n";
  }

  FUNCEND();
  return 0;
}
