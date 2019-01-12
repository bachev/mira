#include <iostream>

#include "stdinc/stlincludes.H"
#include "mira/readpool.H"
#include "mira/parameters.H"

#include "skim.H"


int main(int argc, char ** argv)
{

  do{
    try{

      //{
      //	uint32 start=11248;
      //	uint32 end=19462;
      //	start=11000;
      //	end=20000;
      //	ProgressIndicator P(start, end);
      //	P.progress(start);
      //	for(uint32 i=start; i< end; i++) P.progress(i);
      //	P.progress(end);
      //	exit(0);
      //}

      //uint32 i=0;
      //uint32 num_colors=1025;
      //while((num_colors-1) >> ++i);
      //cout << "i " << i << endl;
      //cout << "2 " << (1<<i) << endl;
      //exit(0);

      MIRAParameters P;

      P.setAlignMinScore(1);
      P.setAlignMinOverlap(1);

      P.parse(argc, argv);
      cout << P;

#if 1
      skim_parameters const & skim_params= P.getSkimParams();

      ReadPool thepool(&P);
      //thepool.loadEXPs("fofn.txt");
      //thepool.loadDataFromFASTA("xaa");
      //thepool.loadDataFromFASTA("gbam.fasta");
      //thepool.loadDataFromFASTA("gbam_san454.fasta");
      //thepool.loadDataFromFASTA("gll454Reads.fna");
      //thepool.loadDataFromFASTA("rs454Reads.fna");
      //thepool.loadDataFromFASTA("rs_ref_mut_454Reads.fna");
      thepool.loadDataFromFASTA("good.fasta");

      for(uint32 i=0;i<thepool.size();i++){
	thepool.getRead(i).setUsedInAssembly(true);

	if(i==17086) cout << thepool.getRead(i).getName() << endl;

      }

      vector<uint32> overlapcounter;
      overlapcounter.resize(thepool.size(),0);
      vector<uint32> extends;
      bannedoverlappairs_t permanent_overlap_bans;
      permanent_overlap_bans.resize(thepool.size());
      
      possible_overlaps_t posfmatch;
      possible_overlaps_t poscmatch;

      Skim s2;
      
//  void skimGo (ReadPool & rp, 
//	       possible_overlaps_t  & posfmatch,
//	       possible_overlaps_t  & poscmatch,
//	       bannedoverlappairs_t & bannedoverlaps,
//	       vector<uint32>       & extends,
//	       vector<uint32>       & overlapcounter,
//
//	       uint32 maxmemusage = 500000000,
//	       
//	       bool verifyhashes=false, 
//	       bool takeextalso=false,
//	       
//	       uint32 initialmatch    = 4,
//	       uint32 matchconfirm    = 10,
//	       int32  percentrequired = 25,
//	       uint32 maxhitsperread  = 200);
  
      string pf="posfmatch.txt";
      string pc="poscmatch.txt";
      s2.skimGo(thepool,
		pf,
		pc,
		permanent_overlap_bans,
//		extends,
		overlapcounter,
		15000000,
//		1,
//		false,
//		false,
		false,
		16,
		4,
		50,
		skim_params.sk_maxhitsperread
	);

     //for(uint32 runningid=0;runningid<thepool.size();){
      //	cout << "Runningid: " << runningid << endl;
      //	runningid=s2.prepareSkim(thepool,runningid);
      //	s2.go(posfmatch, 
      //	      poscmatch,
      //	      permanent_overlap_bans,
      //	      extends,
      //	      overlapcounter,
      //	      skim_params.sk_initialmatch,
      //	      skim_params.sk_matchconfirm,
      //	      skim_params.sk_percentrequired,
      //	      skim_params.sk_maxhitsperread
      //	  );
      //}

      //Skim theskim(thepool);
      //
      //possible_overlaps_t posfmatch;
      //possible_overlaps_t poscmatch;
      //
      //theskim.go(posfmatch, poscmatch);

      if(0) {
	{
	  possible_overlaps_t::const_iterator I=posfmatch.begin();
	  while(I!=posfmatch.end()){
	    cout << I->first << " " << I->second.otherid << "\t" << I->second.eoffset << "\n";
	    I++;
	  }
	}
	{
	  possible_overlaps_t::const_iterator I=poscmatch.begin();
	  while(I!=poscmatch.end()){
	    cout << I->first << " -" << I->second.otherid << "\t" << I->second.eoffset << "\n";
	    I++;
	  }
	}
      }

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




