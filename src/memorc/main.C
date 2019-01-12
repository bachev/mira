#include "memorc.H"

#include <boost/algorithm/string.hpp>

#include "util/misc.H"

using std::cout;

// in MIRA: in preventinifiasco, not in library!
#ifdef MIRAMEMORC
MemORC::mbmap_t MemORC::MOC_memblocks;
MemORC::mbmap_t MemORC::MOC_hotblocks;
vector<uint64> MemORC::MOC_hotaidsrequested;

MemORC MemORC::MOC_semaphore; // keep last for memorc: when instantiated, sets readytouse
                   //  when destructed, clears readytouse
#endif


int main2(int argc, char ** argv)
{
#ifdef MIREMEMORC
  cout << "Start 1\n";
  std::string xx("sjkfghasfjkh");
  cout << "Start 2\n";
  std::string xxx("sjkfghasfjkh1");
  cout << "Start 3\n";
  std::string xxxx("sjkfghasfjkh2");
  cout << "Start 4\n";
  std::string xxxxx("sjkfghasfjkh3");
  cout << "Start 5\n";
  std::string y=xx;
  cout << "Start 6\n";
  y+=xxx;
  cout << "Start 7\n";
  y+=xxxx;
  cout << "Start 8\n";
  y+=xxxxx;
  cout << "Start 9\n";

  char * bla;
  printf("let's go\n");
  bla=new char[10];
  printf("got %p\n", bla);
  fflush(stdout);

  char * ptr=bla;

  for(int32 i=0; i<10; i++, ptr++){
    *ptr=1;
  }

  cout << "#### tag" << endl;
  //multitag_t t;

  cout << "#### std::string" << endl;
  std::string path;
  std::string miraprog;

  cout << "#### split" << endl;


  // bang almost immediately, at next new/delete call
  //*(ptr-40)=0;
  // bang, but later, on the delete[]
  bla+=2;


  splitFullPathAndFileName(argv[0],path,miraprog);
  cout << "#### done split" << endl;

  cout << "#### lower" << endl;
  boost::to_lower(miraprog);
  cout << "#### done lower" << endl;

  cout << "#### cout" << endl;
//  cout << t;


  char * ble;
  ble=new char[10];
  delete [] ble;
  printf("got %p\n", ble);

  delete [] bla;
#endif
  return 0;
}

int main(int argc, char ** argv)
{
#ifdef MIRAMEMORC
  MemORC::setChecking(true);
  cout << "start 000\n";
  auto x=main2(argc,argv);
#else
  printf("Not compiled with MemORC, cannot test.\n");
  int x=0;
#endif
  return x;
}
