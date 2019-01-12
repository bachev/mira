#include <iostream>

#include "stdinc/stlincludes.H"

#include "mira/assembly.H"
#include "mira/readpool.H"
#include "mira/contig.H"
#include "EdIt/hypothesen.H"
#include "examine/scf_look.H"
#include "caf/caf.H"

#include "version.H"


MIRAParameters P;
ReadPool       readpool(&P);
list<Contig>   contigs;

void usage()
{
  cerr << "tagsnp\t(MIRALIB version " << MIRALIBVERSION << ")\n\n";
  cerr << "Usage:\n";

  cerr << "Options:\n";
  cerr << "\t-s <filename>\tload strain data from file\n";
  cerr << "\t-a\t\tassume SNPs instead of repeats\n";
  cerr << "\t-r <int>\tminimum reads per group (default: 1)\n";
  cerr << "\t-q <int>\tminimum qual for tagging (default: 30)\n";
  cerr << "\t-n <int>\tminimum neighbour qual for tagging (default: 20)\n";
  cerr << "\t-e <int>\tend read exclusion area (default: 25)\n";
  cerr << "\t-g\t\talso mark gap bases\n";
  cerr << "\t-m\t\talso mark multicolumn gap bases\n";

}


void save(string & cafout)
{

  if(!cafout.empty()){
    assout::saveAsCAF(contigs,cafout);
  } else {
    assout::dumpAsCAF(contigs,cout);
  }

  //
  //filename="out.tcs";
  //assout::saveAsTCS(contigs,filename);

  //filename="tagsnp_out.gap4da";
  //assout::saveAsGAP4DA(contigs,filename);

  //filename="featureanalysis.txt";
  //assout::saveFeatureAnalysis(400,100,contigs,readpool,
  //				filename,
  //				"featuresummary.txt",
  //				"featureprot.txt");

  //{
  //  string filename="out.html";
  //  cout << "Saving contigs to htmlfile: " << filename << endl;
  //  ofstream out(filename.c_str(), ios::out | ios::trunc);
  //  assout::dumpContigListAsHTML(contigs, "Super project", out);
  //  out.close();
  //}
}

void load (MIRAParameters * mp, string & cafin, string & strainin)
{
  cerr << "Loading project from CAF file: " << cafin << endl;

  CAF tcaf(readpool, contigs, mp);
  tcaf.load(cafin);

  if(!strainin.empty()){
    cerr << "Loading strain data";
    readpool.loadStrainData(strainin);
  }

  Assembly::refreshContigAndReadpoolValuesAfterLoading(readpool,contigs);
}

void doit()
{
  cout << "Tagging reads ..." << endl;
  list<Contig>::iterator I=contigs.begin();
  while(I!=contigs.end()){
    I->setParams(&P);

    uint32 numSRMB=0;
    uint32 numWRMB=0;
    uint32 numSNP=0;
    I->transposeReadSRMTagsToContig();
    //I->markPossibleRepeats(numSRMB, numWRMB, numSNP);
    vector<bool> readsmarkedsrm;
    I->newMarkPossibleRepeats(numSRMB,readsmarkedsrm);
    I++;
  }
}

void doit2()
{
  cout << "Tagging reads ..." << endl;
  list<Contig>::iterator I=contigs.begin();
  for(;I!=contigs.end(); I++){
    I->trashConsensusCache();
    I->markMissedSNPs();
  }
}


//            Assume SNP instead of repeats (asir)                 : No
//            Minimum reads per group needed for tagging (mrpg)    : 1
//            Minimum neighbour quality needed for tagging (mnq)   : 3
//            Minimum Group Quality needed for RMB Tagging (mgqrt) : 3
//            End-read Marking Exclusion Area in bases (emea)      : 25
//            Also mark gap bases (amgb)                           : Yes
//                Also mark gap bases - even multicolumn (amgbemc) : Yes


int main(int argc, char ** argv)
{
  FUNCSTART("int main(int argc, char ** argv)");

  int c;
  extern char *optarg;
  extern int optind;

  P.setContigAssumeSNPInsteadofRepeats(false);
  P.setContigMinReadsPerGroup(1);
  P.setContigMinNeighbourQuality(20);
  P.setContigMinGroupQuality(30);
  P.setContigEndReadMarkExclusionArea(25);
  P.setContigMarkGapBases(false);
  P.setContigMarkMulticolumnGapBases(false);

  P.setContigDisregardSpuriousRMBMismatches(false);

  string cafin="";
  string strainin="";

  while (1){
    c = getopt(argc, argv, "agms:r:n:q:e:");
    if(c == -1) break;

    switch (c) {
    case 's': {
      strainin=optarg;
      break;
    }
    case 'a': {
      P.setContigAssumeSNPInsteadofRepeats(true);
      break;
    }
    case 'r': {
      P.setContigMinReadsPerGroup(atoi(optarg));
      break;
    }
    case 'n': {
      P.setContigMinNeighbourQuality(atoi(optarg));
      break;
    }
    case 'q': {
      P.setContigMinGroupQuality(atoi(optarg));
      break;
    }
    case 'e': {
      P.setContigEndReadMarkExclusionArea(atoi(optarg));
      break;
    }
    case 'g': {
      P.setContigMarkGapBases(true);
      break;
    }
    case 'm': {
      P.setContigMarkMulticolumnGapBases(true);
      break;
    }
    case 'h':
    case '?': {
      usage();
      exit(0);
    }
    default : {}
    }
  }

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing at least infile as argument!\n";
    usage();
    exit(1);
  }

  string infile=argv[optind++];
  string outfile="";

  if(argc-optind > 1) {
    cerr << argv[0] << ": " << "More than infile or outfile as arguments?\n";
    usage();
    exit(1);
  }

  if(argc-optind == 1) {
    outfile=argv[optind];
  }

  try{
    load(&P, infile, strainin);

    doit();
    save(outfile);

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
