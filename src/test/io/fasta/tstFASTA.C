#include <limits.h>
#include <iostream>

#include "io/fasta.H"
#include "stdinc/stlincludes.H"


#define NAMELEN 40


char singletest[][NAMELEN]= { 
  {"good.fasta"}, 
  {"emptyline.fasta"},
  {"multiple.fasta"},
  {"nogreater.fasta"},
  {"noname.fasta"},
  {"sillyname.fasta"},
  {"badspace.fasta"},
  {"wrongbase.fasta"},
  {"yyy.fasta"}
};

char qualtest[][NAMELEN]={
  {"good.fasta.qual"},
  {"basemissing.fasta.qual"},
  {"basemore.fasta.qual"},
  {"emptyline.fasta.qual"},
  {"multiple.fasta.qual"},
  {"nogreater.fasta.qual"},
  {"nonumber.fasta.qual"},
  {"tohigh.fasta.qual"},
  {"namemismatch.fasta.qual"},
  {"yyy.fasta.qual"}
};

char multipletest[][NAMELEN]={
  {"multiple.fasta.qual"},
  {"multiple-bad1.fasta.qual"},
  {"multiple-bad2.fasta.qual"},
  {"multiple-bad3.fasta.qual"},
  {"multiple-bad4.fasta.qual"},
};

int main(void)
{
  FUNCSTART("main");

  FASTA test;
  
  cout << "Testing empty object\n------------------------------------------------" << endl;
  try{
    test.dumpSequence(cout);
    test.dumpQuality(cout);
  }
  catch(Notify n){
    n.setGravity(Notify::WARNING);
    n.handleError("main()");
  }
  cout << "\n\n";

  for(uint32 i=0; i < sizeof(singletest)/NAMELEN; i++ ) {
    cout << "############ 1-" << i << endl;
    try{
      cout << "Testing " << singletest[i] << "\n------------------------------------------------" << endl;
      test.load(singletest[i]);
      test.dumpSequence(cout);
      test.dumpQuality(cout);
    }
    catch(Notify n){
      n.setGravity(Notify::WARNING);
      n.handleError("main()");
    }
    cout << "\n\n";
  }

  for(uint32 i=0; i < sizeof(qualtest)/NAMELEN; i++ ) {
    cout << "############ 2-" << i << endl;
    try{
      cout << "Testing " << qualtest[i] << "\n------------------------------------------------" << endl;
      test.load("good.fasta", qualtest[i]);
      test.dumpSequence(cout);
      test.dumpQuality(cout);
    }
    catch(Notify n){
      n.setGravity(Notify::WARNING);
      n.handleError("main()");
    }
    cout << "\n\n";
  }

  cout << "Testing load qual only\n------------------------------------------------" << endl;
  try{
    test.loadQual("good1.fasta.qual");
    test.dumpSequence(cout);
    test.dumpQuality(cout);
  }
  catch(Notify n){
    n.setGravity(Notify::WARNING);
    n.handleError("main()");
  }
  cout << "\n\n";

  cout << "Testing testIfSeqAndQualMatch\n------------------------------------------------" << endl;
  try{
    test.loadQual("good_bad.fasta.qual");
    test.dumpSequence(cout);
    test.dumpQuality(cout);
    test.testIfSeqAndQualMatch();
  }
  catch(Notify n){
    n.setGravity(Notify::WARNING);
    n.handleError("main()");
  }
  cout << "\n\n";

  cout << "Testing discard object\n------------------------------------------------" << endl;
  try{
    test.discard();
    test.dumpSequence(cout);
    test.dumpQuality(cout);
  }
  catch(Notify n){
    n.setGravity(Notify::WARNING);
    n.handleError("main()");
  }
  cout << "\n\n";


  for(uint32 i=0; i < sizeof(multipletest)/NAMELEN; i++ ) {
    cout << "############ 3-" << i << endl;
    try{
      cout << "Testing " << multipletest[i] << "\n------------------------------------------------" << endl;
      string fastain="multiple.fasta";
      ifstream fin;
      fin.open(fastain.c_str(), ios::in|ios::ate);
      if(!fin){
	MIRANOTIFY(Notify::WARNING, "File not found: " << fastain);
      }
      if(!fin.tellg()){
	MIRANOTIFY(Notify::FATAL, "Zero length file: " << fastain);
      }
      fin.seekg(0, ios::beg);
      
      ifstream qin;
      qin.open(multipletest[i], ios::in|ios::ate);
      if(!qin){
	MIRANOTIFY(Notify::WARNING, "File not found: " << multipletest[i]);
      }
      if(!qin.tellg()){
	MIRANOTIFY(Notify::FATAL, "Zero length file:"  << multipletest[i]);
      }
      qin.seekg(0, ios::beg);
      while(1){
	cout << "Loading ... ";
	cout.flush();
	test.loadNext(fin, qin);
	if(test.testIfEmpty()) {
	  cout << "no more sequences.\n";
	  break;
	}
	cout << "done." << endl;
	test.dumpSequence(cout);
	test.dumpQuality(cerr);
    }
    }
    catch(Notify n){
      n.setGravity(Notify::WARNING);
      n.handleError("main()");
    }
    cout << "\n\n";
  }



  cout << "Tests ended.\n";

  FUNCEND();
  return 0;
}


