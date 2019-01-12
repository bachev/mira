#include <climits>
#include <iostream>

#include "io/fasta.H"
#include "mira/parameters.H"
#include "mira/align.H"
#include "mira/readgrouplib.H"

using namespace std;




//#include "s1s2_pacbio.inc"
//#include "s1s2_gaps.inc"
//include "s1s2_ooops.inc"
#include "s1s2_ooops2.inc"


#define TCSIZE 4*1024*1024
uint32 trashcache1[TCSIZE];
uint32 trashcache2[TCSIZE];

void trashCache()
{
  uint32 * sptr=trashcache1;
  uint32 * dptr=trashcache2;
  for(uint32 i=0; i<TCSIZE; ++i, ++sptr, ++dptr){
    *dptr=*sptr;
  }
}

int main(int argc, char ** argv)
{
  FUNCSTART("main");

  try{
    vector<MIRAParameters> Pv;
    MIRAParameters::setupStdMIRAParameters(Pv);

    //MIRAParameters::parse(argc, argv, Pv);
    string miraparams("--job=denovo,genome,accurate,pcbiolq PCBIOLQ_SETTINGS -AL:mrs=50");
    MIRAParameters::parse(miraparams.c_str(),Pv,false);

    uint8 st=ReadGroupLib::SEQTYPE_PACBIOLQ;

    auto & ap = Pv[st].getNonConstAlignParams();


    list<AlignedDualSeq> madsl;
    Align bla(&Pv[st]);
    //bla.setUseRLE(true);

    for(int ci=0; ci<1; ++ci){
      if(ci%10==0) cout << ci << endl;
      madsl.clear();
      trashCache();

      bool enforce_clean_ends=false;
      bool dontpenalisengaps=true;
      ap.al_kmin=7040;
      ap.al_kmax=7040;
      //ap.al_kmin=394;
      //ap.al_kmax=394;
      bla.acquireSequences((const char*) seq1, strlen(seq1),
      			   (const char*) seq2, strlen(seq2),
      			   1, 2, 1, 1,
      			   true, 0);

      bla.fullAlign(&madsl,enforce_clean_ends,dontpenalisengaps);
      //if(ci==0) bla.resetTimings();
    }

    bla.coutWhatWasGiven();

    list<AlignedDualSeq>::iterator I;
    cout << "madsl.size(): " << madsl.size() << endl;
    cout << "band hit: " << bla.wasBandHit() << endl;

    {
      for(I= madsl.begin(); I!=madsl.end(); I++){
	if(I->isValid()==true){
	  cout << *I;
	}
      }
    }

    cout << "\n\n\nmadsl has " << madsl.size() << " elements.\n" <<endl;

    bla.dumpTimings();

  }
  catch(Notify n){
    n.handleError("main");
  }
  cout << "\n\n";

  cout << "Tests ended.\n";

  FUNCEND();
  return 0;
}
