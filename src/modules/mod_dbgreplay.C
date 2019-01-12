#include <iostream>
#include <string>

#include "mira/assembly.H"
#include "mira/align.H"
#include "mira/parameters.H"
#include "mira/maf_parse.H"



using namespace std;


int dbgReplayMain(int argc, char ** argv)
{
  FUNCSTART("dbgReplayMain");

  string datafile="e.maf";
  string refrname="HUZ85:2416:1763";
  string newrname="HUZ85:1578:1024";
  string paramstring="--job=denovo,genome,accurate,iontor";

  vector<MIRAParameters> MPv;
  MIRAParameters::setupStdMIRAParameters(MPv);

  if(!paramstring.empty()){
    MIRAParameters::parse(paramstring.c_str(),MPv);
  }

  ReadPool       greadpool;
  list<Contig>   gcontigs;

  // not used, but needed for linking *sigh*
  Manifest manifest;
  Assembly as(manifest,MPv,false);

  MAFParse tmaf(&greadpool, &gcontigs, &MPv);
  vector<uint32> dummy;
  tmaf.load(datafile, ReadGroupLib::SEQTYPE_SANGER, 1, dummy, false);
  greadpool.makeTemplateIDs(NWWARN);

  greadpool.allowNameIndex(true);
  int32 refid=greadpool.getReadIndex(refrname);
  BUGIFTHROW(refid<0,"Ref " << refrname << " not in readpool?");
  int32 newid=greadpool.getReadIndex(newrname);
  BUGIFTHROW(newid<0,"New " << newrname << " not in readpool?");

  AlignedDualSeqFacts adsf;
  // refid, newid,
  //   delta,id1_rightdelta, id2_rightdelta, totallen, id1_direction, id2_direction, score_ratio
  //adsf.publicinit(refid,newid,2,0,9,91,1,1,100);
  int8 dirnewid=-1;
  bool newid_ismulticopy=true;
  adsf.publicinit(newid,refid,14,41,0,242,1,1,94,4,0,0,28,0);

  cout << adsf << endl;

  Contig::errorstatus_t errstat;
  vector<Align> aligncache;
  for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++) {
    Align a(&MPv[i]);
    aligncache.push_back(a);
  }

  cout << "before\n";
  Contig::setCoutType(Contig::AS_TEXT);
  cout << gcontigs.back();
  Contig::templateguessinfo_t tguess;
  gcontigs.back().addRead(aligncache,
			  &adsf,
			  refid,newid,
			  dirnewid,
			  newid_ismulticopy,
			  0,
			  tguess,
			  errstat);

  if(errstat.code == Contig::ENOERROR) {
    cout << "read added\n\n";
    Contig::setCoutType(Contig::AS_TEXT);
    cout << gcontigs.back();
  }else{
    cout << "read not added\n\n";
    Contig::setCoutType(Contig::AS_TEXT);
    cout << gcontigs.back();
  }


  return 0;
}
