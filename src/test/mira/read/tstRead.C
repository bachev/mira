#include <limits.h>
#include <iostream>

#include "stdinc/stlincludes.H"
#include "mira/parameters.H"
#include "mira/read.H"


#define NAMELEN 40


char transferSVTagsToClip_tests[][NAMELEN]= { 
  {"rat1.exp"},
  {"rat2.exp"},
  {"rat3.exp"},
  {"rat4.exp"},
  {"rat5.exp"},
  {"rat7.exp"}
};

char brokenEXPs_tests[][NAMELEN]= {
  {"rat6.exp"}
};

char coutFASTA_tests[][NAMELEN]= { 
  {"rat1.exp"}
};

int main(int argc, char ** argv)
{
  FUNCSTART("main");

  try{
    vector<MIRAParameters> Pv;
    MIRAParameters::setupStdMIRAParameters(Pv);

    Read r;
    r.setCoutType(Read::AS_TEXTSHORT);

    cout << "Testing Read::transferSVTagsToClip\n";

    for(uint32 i=0; i < sizeof(transferSVTagsToClip_tests)/NAMELEN; i++ ) {
      cout << "############ 1-" << i << endl;
      try{
	cout << "Testing " << transferSVTagsToClip_tests[i] << "\n------------------------------------------------" << endl;
	cout << "Original:\n";
	r.discard();
	r.loadDataFromEXP(transferSVTagsToClip_tests[i]);
	cout << r;
	r.transferSVTagsToClip(10,50);
	cout << "\n\n\nAfter transfer:\n";
	cout << r;
      }
      catch(Notify n){
	n.setGravity(Notify::WARNING);
	n.handleError("main()");
      }
      cout << "\n\n";
    }


    cout << "Testing broken EXPs\n";
    
    for(uint32 i=0; i < sizeof(brokenEXPs_tests)/NAMELEN; i++ ) {
      cout << "############ 2-" << i << endl;
      try{
	cout << "Testing " << brokenEXPs_tests[i] << "\n------------------------------------------------" << endl;
	r.discard();
	r.loadDataFromEXP(brokenEXPs_tests[i]);
	cout << r;
      }
      catch(Notify n){
	n.setGravity(Notify::WARNING);
	n.handleError("main()");
      }
      cout << "\n\n";
    }


    cout << "Testing output as FASTA\n";
    
    for(uint32 i=0; i < sizeof(coutFASTA_tests)/NAMELEN; i++ ) {
      cout << "############ 3-" << i << endl;
      try{
	cout << "Testing " << coutFASTA_tests[i] << "\n------------------------------------------------" << endl;
	r.discard();
	r.loadDataFromEXP(coutFASTA_tests[i]);
	cout << "Complete read:\n";
	Read::setCoutType(Read::AS_FASTA);
	cout << r;

	cout << "Clipped read:\n";
	Read::setCoutType(Read::AS_CLIPPEDFASTA);
	cout << r;
      }
      catch(Notify n){
	n.setGravity(Notify::WARNING);
	n.handleError("main()");
      }
      cout << "\n\n";
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


