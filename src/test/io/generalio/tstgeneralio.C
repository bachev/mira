#include <limits.h>
#include <iostream>

#include "io/generalio.H"

#if __GNUC__ >= 3
#define NO_NOCREATE
#endif 



#define NAMELEN 40


char normaltests[][NAMELEN]= { 
  {"kvin.txt"}
};

int main(int argc, char ** argv)
{
  FUNCSTART("main");
  try{
    cout << "Testing normal operation\n";
    
    for(uint32 i=0; i < sizeof(normaltests)/NAMELEN; i++ ) {
      cout << "############ 1-" << i << endl;
      try{
        cout << "Testing " << normaltests[i] << "\n------------------------------------------------" << endl;
	
	ifstream fin;
#ifdef NO_NOCREATE
	fin.open(normaltests[i], ios::in|ios::ate);
#else 
	fin.open(normaltests[i], ios::in|ios::ate|ios::nocreate);
#endif 
	if(!fin){
	  cerr << "File not found\n";
	}
	if(!fin.tellg()){
	  cerr << "Zero length file\n";
	}
	fin.seekg(0, ios::beg);

	string key, value;
	while(GeneralIO::readKeyValue(fin, key, value)){
	  cout << "key \"" << key << "\"\tvalue \"" << value << "\"" << endl;
	}
	cout << "\n\nKV list ended.\n";
      }
      catch(...){
	cerr << "Unknown exception trapped.\n";
      }
    }

  }
  catch(...){
    cerr << "Unknown exception trapped.\n";
  }


  cout << "\n\n";
  cout << "Tests ended.\n";

  FUNCEND();
  return 0;
}

