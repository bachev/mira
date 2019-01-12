#include <iostream>

#include <limits.h>

#include "io/scf.H"
#include "io/exp.H"

#include "stdinc/stlincludes.H"

char * scfname = "5.scf";

int main(void)
{
  SCF * test1;
  test1=new SCF;

  try{
    test1->load(scfname);
    //    test1->save("good.scf");
  }
  catch(Notify n){
    n.handleError("main()");
  }
  test1->dump();
  //  test1->discard();


  return 0;
}

