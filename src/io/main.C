#include "scf.H"
#include "exp.H"

#include "stdinc/stlincludes.H"

char * scfname = "5.scf";
//char * scfname = "U13a010a06.t1SCF";
//char * expname = "U60f02a09.t1";
char * expname = "test.exp";
//char * filename = "xaa";
//char & filename=scffilename;


int main(void)
{
#if 0
  EXP * test2;
  test2=new EXP;
  try{
    //  SCF * test1;
    //  test1=new SCF;
    //    test1->load(scfname);
    test2->load(expname);

    test2->dump();
  }
  catch(Notify n){
    n.handleError("main()");
  }
  //  cout << test2->search('S','Q');

#else
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
#endif

  return 0;
}
