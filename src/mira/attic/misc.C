
This file contains miscellaneous code that might serve again some day.


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////        Obsolete         ///////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


/*************************************************************************
 *
 * expects reads to have baseflags set  (by performHashAnalysis())
 *
 * doesn't seem to be a good idea
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

/*

void Assembly::performHashEditing()
{
  FUNCSTART("void Assembly::performHashEditing()");

  cout << "Hash analysis for editing:";

  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  uint32 basesperhash=as_fixparams.as_clip_pec_basesperhash;
  if(sizeof(uint64) < 8 && basesperhash > 15) basesperhash=15;
  {

    Skim s3;

    s3.setHashFrequencyRatios(skim_params.hs_freqest_minnormal,
			      skim_params.hs_freqest_maxnormal,
			      skim_params.hs_freqest_repeat,
			      skim_params.hs_freqest_heavyrepeat,
			      skim_params.hs_freqest_crazyrepeat,
			      skim_params.hs_nastyrepeatratio);

    s3.analyseHashes(AS_miraparams[0].getDirectoryParams().dir_tmp,
		     AS_readpool,
		     true,
		     false,
		     false,
		     true,
		     1,
		     basesperhash,
		     1,
		     false);
  }

  if(as_fixparams.as_dateoutput) dateStamp(cout);
  cout << '\n';

  cout << "Looking for proposed edits ... "; cout.flush();

  vector<uint8> maxhf;
  maxhf.reserve(10000);

  uint64 numbaseschanged=0;
  uint64 numreadschanged=0;

  for(uint32 actid=0; actid<AS_readpool.size(); actid++){
    Read & r=AS_readpool.getRead(actid);

    if(r.hasValidData()
       && r.hasBaseHashStats()
       && !(r.isBackbone()
	    || r.isRail())){

      maxhf.clear();
      maxhf.resize(r.getLenSeq(),0);

      bool wasedited=false;

      {
	int32 lpos=r.getLeftClipoff();
	vector<Read::bposhashstat_t>::const_iterator bhsI=r.getBPosHashStats().begin();
	vector<uint8>::iterator mhfI=maxhf.begin();
	advance(bhsI,lpos);
	advance(mhfI,lpos);

	uint32 counter=basesperhash;
	for(; lpos<static_cast<int32>(r.getLenSeq()); lpos++, bhsI++, mhfI++) {
	  *mhfI=(bhsI->fwd.getFrequency())>1;
	  if(*mhfI) counter=basesperhash;
	  if(counter) {
	    *mhfI=4;
	    --counter;
	  }
	}

	lpos=r.getLeftClipoff();
	mhfI=maxhf.begin();
	advance(mhfI,lpos);

	//for(; lpos<static_cast<int32>(r.getLenSeq()); lpos++) {
	//  cout << (uint16) maxhf[lpos] << ' ';
	//}
	//cout << endl;
	//lpos=r.getLeftClipoff();
	//for(; lpos<static_cast<int32>(r.getLenSeq()); lpos++) {
	//  cout << r.getBaseInSequence(lpos) << ' ';
	//}
	//cout << endl;
	//Read::setCoutType(Read::AS_TEXT);
	//cout << r;

	lpos=r.getLeftClipoff();
	for(; lpos<static_cast<int32>(r.getLenSeq()); lpos++, mhfI++) {
	  if(*mhfI) break;
	}

	int32 editstart=-1;
	for(; lpos<static_cast<int32>(r.getLenSeq()); lpos++, mhfI++) {
	  if(editstart<0){
	    if(*mhfI==0) {
	      editstart=lpos;
	    }
	  }else{
	    if(*mhfI) {
	      for(int32 ii=editstart; ii<lpos; ii++) {
		//editpositions.push_back(ii);
		r.changeBaseInSequence('n',0,ii);
		numbaseschanged++;
		wasedited=true;
	      }
	      editstart=-1;
	    }
	  }
	}

      }
      if(wasedited) numreadschanged++;

      //if(editpositions.size()){
      //	cout << r.getName() << ": wants to edit " << editpositions.size() << " positions\n";
      //}
    }
  }

  cout << "changed " << numbaseschanged << " bases to 'n' in " << numreadschanged << " reads.\n";

  FUNCEND();

  return;
}
//#define CEBUG(bla)
*/


/*************************************************************************
 *
 * BaCh 31.12.2012: Errrrm ... what's that function for? Was probably a
 *  quick hack for something or I planed some "extra" file for very
 *  special cases. Should probably be removed.
 *
 * REMOVEME!
 *
 * ugly and slow, but works and is fast enough
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
/*
void Assembly::mergeTemplateInfo(const string & tifile, const string & logname, const string & logprefix)
{
  FUNCSTART("void Assembly::mergeTemplateInfo(const string & tifile, const string & logname, const string & logprefix)");

  cout << "Merging template info from " << tifile << ":\n";

  CEBUG("Building hash table ... "); cout.flush();

  typedef boost::unordered_map<std::string, int32> strmap;
  strmap rnmap;
  strmap::iterator rnI;

  for(uint32 i=0; i<AS_readpool.size();i++){
    if(!AS_readpool[i].getName().empty()) {
      rnmap[AS_readpool[i].getName()]=i;
    }
  }
  CEBUG("done." << endl);

  ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname.c_str(), ios::out|ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  ifstream tifin;
  tifin.open(tifile.c_str(), ios::in|ios::ate);
  if(!tifin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << tifile);
  }
  streampos tifsize=tifin.tellg();
  tifin.seekg(0, ios::beg);

  ProgressIndicator<streamsize> P (0, tifsize,1000);

  string token;

  while(!tifin.eof()){
    tifin >> token;
    if(tifin.eof()) break;
    if(P.delaytrigger()) P.progress(tifin.tellg());

    //tifin >> sd_score >> sd_readname;

    if(tifin.eof()) break;

    if(token[0]=='+'){
      // new lib
    }else{
      // existing name
      bool foundname=false;
      rnI=rnmap.find(token);
      if(rnI==rnmap.end()) {
	CEBUG("Not found: " << token << endl);
	continue;
      }
      uint32 foundreadid=rnI->second;
      if(!AS_readpool[foundreadid].hasValidData()) continue;

      Read actread(AS_readpool[foundreadid]);
      assembly_parameters const & as_params= AS_miraparams[actread.getSequencingType()].getAssemblyParams();
    }
  }
  P.finishAtOnce();

  tifin.close();

  if(!logname.empty()){
    logfout.close();
  }

  cout << "\nDone." << endl;


  FUNCEND();
  return;
}
//#define CEBUG(bla)
*/




/*************************************************************************
 *
 * splits a sequence into overlapping subsequences
 *
 * AND
 *
 * saves pre-computed adsfacts file into log directory for later
 *  later reading
 * number of generated adsfacts is put in AS_numADSFacts_fromshreds
 *
 *
 * This saves enormous amount of time, but is not the "real" thing:
 *  matches between shreds that are non-overlapping from the start on are
 *  not made
 *
 *************************************************************************/

/*
void Assembly::shredReadsIntoReadPool(ReadPool & sourcepool, uint32 shredlen, uint32 shredoffsetinc, uint8 shredreadtype, const string & shredstrain)
{
  FUNCSTART("void Assembly::shredReadsIntoReadPool(ReadPool & sourcepool, uint32 shredlen, uint32 shredoffsetinc, uint8 shredreadtype, const string & shredstrain)");

  AS_numADSFacts_fromshreds=0;
  string adsfshredsfilename=AS_miraparams[0].getDirectoryParams().dir_tmp+"/shred.adsfacts";
  ofstream adsfout;
  adsfout.open((adsfshredsfilename+".adsfacts").c_str(), ios::out|ios::trunc);

  deque<uint32> overlapfifo;

  string shredseq;
  shredseq.reserve(shredlen);
  vector<base_quality_t> shredqual;
  shredqual.reserve(shredlen+10);
  string shredname;

  for(uint32 actsourceid=0; actsourceid < sourcepool.size(); actsourceid++){
    Read & sourceread = sourcepool.getRead(actsourceid);
    if(!sourceread.hasValidData()) continue;
    if(sourceread.getLenSeq() < shredlen) continue;

    uint32 actoffset=0;
    uint32 shredcounter=0;
    for(bool doloop=true; doloop; actoffset+=shredoffsetinc){
      uint32 fromi=actoffset;
      uint32 toi=actoffset+shredlen;
      if(toi>=sourceread.getLenSeq()) {
	toi=sourceread.getLenSeq();
	doloop=false;
      }
      shredseq.clear();
      shredqual.clear();
      for(; fromi<toi; fromi++){
	shredseq+=sourceread.getBaseInSequence(fromi);
	shredqual.push_back(sourceread.getQualityInSequence(fromi));
      }

      // if wished: lower quals to max as_cap454consensusqual
      if(AS_miraparams[0].getAssemblyParams().as_cap454consensusqual>0){
	vector<base_quality_t>::iterator qI=shredqual.begin();
	base_quality_t maxqual=AS_miraparams[0].getAssemblyParams().as_cap454consensusqual;
	for(;qI != shredqual.end(); qI++){
	  if(*qI>maxqual) *qI=maxqual;
	}
      }

      ostringstream ostr;
      ostr << "shred_" << shredcounter << "_" << sourceread.getName();
      shredname=ostr.str();

      AS_readpool.addNewEmptyRead();
      uint32 newreadid=AS_readpool.size()-1;
      Read & newread=AS_readpool.getRead(newreadid);
      newread.setName(shredname);
      newread.setSequenceFromString(shredseq);
      newread.setQualities(shredqual);
      newread.setStrain(shredstrain.c_str());
      newread.setSequencingType(shredreadtype);

      //cout << "\n----------------------------------------\nAdded " << shredname << '\n';
      // now insert the weights
      {
	overlapfifo.push_front(newreadid);
	deque<uint32>::iterator OFI=overlapfifo.begin();
	OFI++;
	int32 overlaplen=shredlen-shredoffsetinc;
	int32 totalshredoffset=shredoffsetinc;
	uint32 numelements=1;
	while(OFI != overlapfifo.end()) {
	  if(overlaplen<=0) break;

	  AlignedDualSeqFacts tmpadsf;
	  tmpadsf.publicinit(
	    *OFI,
	    newreadid,
	    static_cast<uint16>(totalshredoffset),
	    static_cast<uint16>(totalshredoffset
				-(AS_readpool.getRead(*OFI).getLenSeq()
				  -AS_readpool.getRead(newreadid).getLenSeq())),
	    0,
	    static_cast<uint16>((AS_readpool.getRead(*OFI).getLenSeq()+
				 AS_readpool.getRead(newreadid).getLenSeq()-overlaplen)),
	    1,
	    1,
	    100);

	  // output of the ADSfacts to file
	  // TODO: real ouput
	  // first weight and direction
	  // TODO: reduce weight to favorise real reads in assembly???
	  adsfout << overlaplen*10000 << "\t1\t";
	  tmpadsf.serialiseOut(adsfout);
	  adsfout << '\n';

	  AS_numADSFacts_fromshreds++;

	  OFI++;
	  overlaplen-=shredoffsetinc;
	  totalshredoffset+=shredoffsetinc;
	  numelements++;
	}
	if(overlapfifo.size()>numelements) overlapfifo.resize(numelements);
      }
      shredcounter++;
    }
    cout << "Shredded " << sourceread.getName() << " into " << shredcounter << " pieces.\n";
  }

  adsfout.close();

  FUNCEND();
}
*/


void HashStatistics::old_loadHashStatistics(string & hashstatfilename, uint8 basesperhash)
{
  FUNCSTART("void HashStatistics::old_loadHashStatisticsFile(string & hashstatfilename, uint8 basesperhash)");

  HS_hs_basesperhash=basesperhash;

  HS_hsv_hashstats.clear();
  HS_hsv_hsshortcuts.clear();

  HS_avg_freq_corrected=0;
  HS_avg_freq_raw=0;
  HS_avg_freq_taken=0;

  BUGIFTHROW(!fileExists(hashstatfilename),"No hash statistics file " << hashstatfilename << " to load data from?");
  auto fs=boost::filesystem::file_size(hashstatfilename);
  BUGIFTHROW(fs==0,"Empty file " << hashstatfilename << " ?");
  if(fs%sizeof(hashstat_t)){
    MIRANOTIFY(Notify::FATAL, "File probably not a hash stat: " << hashstatfilename);
  }

  HS_hsv_hashstats.resize(fs/sizeof(hashstat_t));

  CEBUG("Loading hash stats: " << HS_hsv_hashstats.size() << endl);

  FILE * fin;
  fin=fopen(hashstatfilename.c_str(), "r");
  size_t numread=myFRead(&HS_hsv_hashstats[0],sizeof(hashstat_t),HS_hsv_hashstats.size(),fin);
  if(numread != HS_hsv_hashstats.size()){
    MIRANOTIFY(Notify::FATAL, "Expected to read " << HS_hsv_hashstats.size() << " elements in hashfile " << hashstatfilename << " but read less (" << numread << "). Was the file deleted? Disk full?");
  }
  fclose(fin);

  if(HS_hsv_hashstats.begin() == HS_hsv_hashstats.end()) return;

  if(HS_logflag_hashcount){
    string logfile=hashstatfilename+".shouldneverbeseen.hashcount.usort";
    ofstream fout;
    fout.open(logfile.c_str(), ios::out);

    for(auto & hs : HS_hsv_hashstats){
      if(hs.hsc.hasfwdrevthresholdok) fout << hs.hsc.getCount() << "\n";
    }
  }

  sort(HS_hsv_hashstats.begin(),HS_hsv_hashstats.end(),sortHashStatComparatorByCountUp);

  if(HS_logflag_hashcount){
    string logfile=hashstatfilename+".shouldneverbeseen.hashcount.sort";
    ofstream fout;
    fout.open(logfile.c_str(), ios::out);

    for(auto & hs : HS_hsv_hashstats){
      if(hs.hsc.hasfwdrevthresholdok) fout << hs.hsc.getCount() << "\n";
    }
  }

  // do this before makeHashStatArrayShortcuts() as mHSAS() resorts
  //  HS_hsv_hashstats, but not by count!
  priv_calcAvgHashFreq();

  priv_makeHashStatArrayShortcuts();
}
