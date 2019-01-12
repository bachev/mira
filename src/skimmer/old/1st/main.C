#include <iostream.h>

#include "stdinc/stlincludes.H"
#include "mira/readpool.H"
#include "mira/parameters.H"

#include "skim.H"

uint32 CLASS_trace=0;
uint32 CLASS_debug=0;

int main(int argc, char ** argv)
{

  do{
    try{



      MIRAParameters P;

      P.setAlignMinScore(1);
      P.setAlignMinOverlap(1);

      P.cmdLine2Stream(argc, argv);
      cout << P;

#if 1
      ReadPool thepool;
      thepool.loadEXPs("fofn.txt");

      Skim theskim(thepool);

      theskim.go();
#endif

    }
    catch(Notify n){
      n.handleError("main");
    }
    catch(Flow f){
      cerr << "Unexpected exception: Flow()\n";
    }
  }while(0);

  return 0;
}




