#include <limits.h>
#include <iostream>

#include "stdinc/stlincludes.H"
#include "mira/parameters.H"
#include "mira/readpool.H"


#define NAMELEN 40


char error[][NAMELEN]= { 
  {"fofnerror.txt"}
};
char normal[][NAMELEN]= { 
  {"fofn.txt"}
};
char straindata[][NAMELEN]= { 
  {"fofn.txt"}
};

int main(int argc, char ** argv)
{
  FUNCSTART("main");

  try{
    vector<MIRAParameters> Pv;
    MIRAParameters::setupStdMIRAParameters(Pv);
    MIRAParameters::parse(argc, argv, Pv);


   {
     ReadPool rp(&Pv);
     
     cout << "Testing error in fofn  operations\n";
     for(uint32 i=0; i < sizeof(error)/NAMELEN; i++ ) {
       cout << "############ 1-" << i << endl;
       try{
	 cout << "Testing " << error[i] << "\n------------------------------------------------" << endl;
	 uint32 dummy=0;
	 rp.loadEXPs(error[i],1,dummy);
	 cout << "Done.";
	 cout << "Readpool contains " << rp.size() << " reads." << endl;
       }
       catch(Notify n){
	 n.setGravity(Notify::WARNING);
	 n.handleError("main()");
       }
       cout << "\n\n";
     }
   }

    {
      ReadPool rp(&Pv);
      
      cout << "Testing normal EXP load operations\n";
      for(uint32 i=0; i < sizeof(normal)/NAMELEN; i++ ) {
	cout << "############ 2-" << i << endl;
	try{
	  cout << "Testing " << normal[i] << "\n------------------------------------------------" << endl;
	  uint32 dummy=0;
	  rp.loadEXPs(normal[i],1,dummy);
	  cout << "Done.";
	  cout << "Readpool contains " << rp.size() << " reads." << endl;
	}
	catch(Notify n){
	  n.setGravity(Notify::WARNING);
	  n.handleError("main()");
	}
	cout << "\n\n";
      }
    }

    {
      ReadPool rp(&Pv);
      
      cout << "Testing straindata load operations\n";
      for(uint32 i=0; i < sizeof(straindata)/NAMELEN; i++ ) {
	cout << "############ 3-" << i << endl;
	try{
	  cout << "Testing " << straindata[i] << "\n------------------------------------------------" << endl;
	  uint32 dummy=0;
	  rp.loadEXPs(straindata[i],1,dummy);
	  cout << "Done.";
	  cout << "Readpool contains " << rp.size() << " reads." << endl;
	}
	catch(Notify n){
	  n.setGravity(Notify::WARNING);
	  n.handleError("main()");
	}
	cout << "\n\n";
      }
    }
  }
  catch(Notify n){
    n.handleError("main");
  }
  cout << "\n\n";

  cout << "Tests ended.\n";

  FUNCEND();
  return 0;
}


