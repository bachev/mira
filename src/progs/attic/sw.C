#include <iostream>

//#include "readpool.H"
#include "mira/parameters.H"
#include "mira/assembly.H"
#include "caf/caf.H"

#include "version.H"



void usage()
{
}


int main(int argc, char ** argv)
{
  FUNCSTART("int main(int argc, char ** argv)");



  int c;
  extern char *optarg;
  extern int optind;

  while (1){
    c = getopt(argc, argv, "h");
    if(c == -1) break;

    switch (c) {
    case 'h':
    case '?':
      usage();
      exit(0);
    }
  }

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing dbfile and queryfile as arguments!\n";
    usage();
    exit(1);
  }

  if(argc-optind < 2) {
    cerr << argv[0] << ": " << "Missing either dbfile or queryfile as arguments!\n";
    usage();
    exit(1);
  }

  if(argc-optind > 2) {
    cerr << argv[0] << ": " << "Whoops, found more than dbfile and queryfile as arguments left on the command line!\n";
    cerr << "Unparsed command line: ";
    for(;optind<argc;optind++) cerr <<argv[optind] << " ";
    cerr << endl;
    usage();
    exit(1);
  }

  string dbfile=argv[optind++];
  string queryfile=argv[optind];

  if(dbfile=="--help"){
    usage();
    exit(0);
  }

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  ReadPool dbpool(&Pv), querypool(&Pv);

  align_parameters & ap= (align_parameters &) Pv[0].getAlignParams();
  assembly_parameters &as= (assembly_parameters &) Pv[0].getAssemblyParams();

  as.as_clip_possible_vectors=false;
  ap.al_kmin=100000;
  ap.al_kmax=100000;
  ap.al_min_score=100;

  try{
    uint32 dummy=0;
    dbpool.loadDataFromFASTA(dbfile,1,dummy, false);
    dummy=0;
    querypool.loadDataFromFASTA(queryfile,1,dummy, false);

    Align bla(&Pv[0]);
    list<AlignedDualSeq> madsl;

    for(int32 iq=0; iq<querypool.size(); iq++) {
      ProgressIndicator<int32> P (0, dbpool.size()-1);
      for(int32 id=0; id<dbpool.size(); id++) {
	bla.acquireSequences((char *) dbpool.getRead(id).getSeqAsChar(),
			     dbpool.getRead(id).getLenSeq(),
			     (char *) querypool.getRead(iq).getSeqAsChar(),
			     querypool.getRead(iq).getLenSeq(),
			     id, iq, 1, 1, true, 0);

	//scheisse, banded macht mist

	madsl.clear();
	bla.fullAlign(&madsl,false,false);
	list<AlignedDualSeq>::iterator I;
	cout << id << "\t" << madsl.size() << endl;
	{
	  int32 bestweight=0;
	  for(I=madsl.begin(); I!=madsl.end(); I++){
	    //CEBUG("MADSL entry:\n" << *I << endl);
	    if(I->isValid()==true){
	      cout << *I;
	    }
	  }
	}
	madsl.clear();
	//P.progress(id);
	cout.flush();
      }
    }

    cout << " done.\n";

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
