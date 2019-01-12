/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 * Copyright (C) 2000 and later by Bastien Chevreux
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */


#include "mira/assembly.H"

#include "mira/align.H"
#include "mira/ppathfinder.H"

#include "util/progressindic.H"
#include "util/stlimprove.H"

// BOOST
//#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/filesystem.hpp>

#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

#if 0
#include <valgrind/memcheck.h>
#define VALGRIND_LEAKCHECK
#endif


using std::cout;
using std::cerr;
using std::endl;

#define CEBUG(bla)

// cs1 for normal clocking ('user compatible' as is does not disturb)
//  cs2 for extensive clocking output, more for analysis of MIRA behaviour
//  during development

#ifndef PUBLICQUIET
#define CLOCK_STEPS2
#endif
#define CLOCK_STEPS2


//#define TRACKMEMUSAGE 1
#define TRACKMEMUSAGE 0




/*************************************************************************
 *
 * returns whether new strong repeat markers (SRMs) were found for
 *  any contig built in any stage
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }

bool Assembly::buildFirstContigs(const int32 passnr, const EDITParameters & eparams, const bool lastpass)
{
  FUNCSTART("void Assembly::buildFirstContigs()");

  CEBUG("BFC: " << passnr << "\t" << lastpass << endl);

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();
  edit_parameters const & ed_params= AS_miraparams[0].getEditParams();

  std::ofstream fout_povl(buildDefaultCheckpointFileName(AS_miraparams[0].getFileParams().chkpt_persistentoverlaps),std::ios::out|std::ios::app|std::ios::ate);

  AS_deleteoldresultfiles=true;

  AS_bfcstats.clear();
  AS_bfcstats.resize(2); // 0 for non-rep contigs, 1 for rep contigs

  AS_contigs.clear();
  Contig::setIDCounter(1);
  Contig::resetCERNumbering();

  AS_assemblyinfo.zeroInfo();

  AS_templateguesses.clear();
  if(AS_maxtemplateid>=0){
    AS_templateguesses.resize(AS_maxtemplateid+1);
    cout << "TGS: " << AS_templateguesses.size() << endl;
  }

  std::vector<Align> aligncache;
  setupAlignCache(aligncache);

  // AS_hasmcoverlaps will be initialised by Pathfinder
  AS_hasmcoverlaps.clear();

  // Keep track of coverage information for contigs so to be able to make
  //  average contig coverage predictions
  std::vector<std::vector<uint32> > covperstpercon(ReadGroupLib::SEQTYPE_END);
  std::vector<uint32> covtotalpercon;

  // initially, overlaps of every read can be reduced
  if(AS_needalloverlaps.empty()){
    AS_needalloverlaps.resize(AS_readpool.size(),false);
  }

  // filter out all reads that have no overlap at this stage
  // for this, first set all reads to "used" and "singlet", then
  //   re-allow those that have at least one SW ovelap

  AS_used_ids.clear();
  AS_used_ids.resize(AS_readpool.size(),1);
  AS_isdebris.clear();
  AS_isdebris.resize(AS_readpool.size(),DEBRIS_NOOVERLAP);

  // overlay AS_debrisreason
  {
    auto dstI=AS_isdebris.begin();
    auto srcI=AS_debrisreason.begin();
    for(; srcI!=AS_debrisreason.end(); ++dstI, ++srcI){
      if(*srcI) *dstI=*srcI;
    }
  }

  // take back all reads with overlaps
  for(const auto & cee : AS_confirmed_edges){
    AS_isdebris[cee.rid1]=0;
    AS_used_ids[cee.rid1]=0;
  }

  // Now, go through all the singlets and take them back into
  //  assembly if they have special MIRA tags attached
  // Also remove backbones & rails from the debris list
  //  and set them to unused
  {
    for(uint32 i=0; i< AS_readpool.size(); i++){
      if(AS_readpool[i].isUsedInAssembly()
	 && AS_readpool[i].getLenClippedSeq()>0
	 && (AS_readpool[i].hasTag(Read::REA_tagentry_idSRMr)
	     || AS_readpool[i].hasTag(Read::REA_tagentry_idCRMr)
	     || AS_readpool[i].hasTag(Read::REA_tagentry_idWRMr)
	     || AS_readpool[i].hasTag(Read::REA_tagentry_idSAOr)
	     || AS_readpool[i].hasTag(Read::REA_tagentry_idSROr)
	     || AS_readpool[i].hasTag(Read::REA_tagentry_idSIOr))) {
	AS_isdebris[i]=0;
	AS_used_ids[i]=0;
      }
      if(AS_readpool[i].isBackbone()
	 || AS_readpool[i].isRail()) {
	AS_used_ids[i]=0;
	AS_isdebris[i]=0;
      }
    }
  }

  // rollback digitally normalised reads to "used, debris"
  for(uint32 rpi=0; rpi< AS_readpool.size(); ++rpi){
    if(AS_debrisreason[rpi]==DEBRIS_DIGITAL_NORMALISATION){
      AS_isdebris[rpi]=DEBRIS_DIGITAL_NORMALISATION;
      AS_used_ids[rpi]=1;
    }
  }


  // get total length of backbones (if any)
  uint32 totalbblen=0;
  for(const auto & cle : AS_bbcontigs){
    totalbblen+=cle.getContigLength();
  }

  std::ofstream fout;
  if(as_fixparams.as_tmpf_unused_ids.size()!=0){
    fout.open((dir_params.dir_tmp+"/"+as_fixparams.as_tmpf_unused_ids), std::ios::out);
    fout.close();
  }

  // outside precomputed lowerbound of oedges, for PathFinder
  // is used lateron in constructStepByStep() of the Pathfinder
  //
  // however, every once in a while (every 10k to 100k reads used from
  //  the pool), the overlap edges vector will be compressed (unecessary
  //  edges thrown out) and therefore this vector will be emptied (to be reconstructed
  //  by pathfinder)
  //
  // On mapping 8m Solexas to 70 contigs (4MB), time goes down from 234 minutes to
  //  221 minutes (Core i940)
  //
  // effect of compression not as big as hoped for, but every little bit helps

  std::vector<necontainer_t::iterator> tmp_lowerbound_oedges;

  // TODO: currently filled by pathfinder, maybe good to have that already in skim?
  AS_hasreptoverlap.clear();
  AS_hasnoreptoverlap.clear();

  // this one really filled by pathfinder
  AS_maybespoilsport.clear();

  // PathFinder object can be created ouside loop and re-used

  PPathfinder qaf(&AS_miraparams,
		  &AS_readpool,
		  &AS_confirmed_edges,
		  &AS_adsfacts,
		  &aligncache,
		  &AS_used_ids,
		  &AS_multicopies,
		  &AS_hasmcoverlaps,
		  &AS_hasreptoverlap,
		  &AS_hasnoreptoverlap,
		  &AS_istroublemaker,
		  &AS_incorchim,
		  &AS_maybespoilsport,
		  &AS_wellconnected,
		  &tmp_lowerbound_oedges,
		  &AS_templateguesses);

  // this vector will hold the read IDs added by pathfinder to contig
  // + backbone + rail IDs
  // used after call to *_constructStepByStep() to set AS_used_ids[] elements
  std::vector<int32> tmp_ids_in_contig;

  uint32 overlapcompressstepping=AS_readpool.size()/5;
  if(overlapcompressstepping<10000) overlapcompressstepping=10000;
  //if(overlapcompressstepping<100) overlapcompressstepping=100;
  if(overlapcompressstepping>100000) overlapcompressstepping=100000;
  if(overlapcompressstepping>AS_readpool.size()) overlapcompressstepping=AS_readpool.size();
  uint32 nextoverlapcompress=overlapcompressstepping;
  cout << "overlapcompressstepping: " << overlapcompressstepping
       << "\nnextoverlapcompress: " << nextoverlapcompress << endl;

  bool foundSRMs=false;
  uint32 numsingletssincecleanup=0;

  bool maykillintermediatesinmglets=true;
  bool shouldmovesmallclusterstodebris=false;
  for(uint32 st=0;st < AS_miraparams.size(); st++){
    if(AS_miraparams[st].getAssemblyParams().as_savesimplesingletsinproject) maykillintermediatesinmglets=false;
    if(AS_miraparams[st].getAssemblyParams().as_minimum_readspercontig>1) shouldmovesmallclusterstodebris=true;
  }

  if(shouldmovesmallclusterstodebris) bfc_moveSmallClustersToDebris();


#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvloop;
  timeval tvtotal;
  gettimeofday(&tvtotal,nullptr);
#endif

  //uint32 countingunused=AS_used_ids.size();

  uint32 trackingunused=AS_used_ids.size();
  for(uint32 i=0; i<AS_used_ids.size(); i++){
    if(AS_used_ids[i]) --trackingunused;
  }

  uint32 numcontigs=1;
  // bug: if someone specifically sets as_maxcontigsperpass to 2^32-1, then
  //  this loop never runs.
  for(;trackingunused>0; ++numcontigs){
    CEBUG("bfc 1\n");
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    cout << '\n';

    bfc_sanityCheckASUSEDIDS(trackingunused,numcontigs);


//    // REMOVEME: paranoia check
//    for(uint32 ri=0; ri<AS_used_ids.size(); ++ri){
//      BUGIFTHROW(AS_isdebris[ri] && AS_used_ids[ri]==0, "AS_isdebris[ri] && AS_used_ids[ri]==0  for ri " << ri << "\t" << AS_readpool[ri].getName());
//    }



    // jump out of for loop if max number of contigs was reached
    if(as_fixparams.as_maxcontigsperpass>0 && numcontigs==as_fixparams.as_maxcontigsperpass+1) break;

    CEBUG("bfc 2\n");
    if(trackingunused>0){

#ifdef CLOCK_STEPS2
      gettimeofday(&tv,nullptr);
      tvloop=tv;
#endif

      //// compress the overlap edges when needed
      //if(AS_readpool.size()-unused > nextoverlapcompress){
      //	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      //	cout << "Compressing overlap edges ..."; cout.flush();
      //	necontainer_t::const_iterator srcI=AS_confirmed_edges.begin();
      //	necontainer_t::iterator dstI=AS_confirmed_edges.begin();
      //	for(; srcI != AS_confirmed_edges.end(); ++srcI){
      //	  if(AS_used_ids[srcI->rid1] == 0
      //	     || AS_used_ids[srcI->linked_with] == 0){
      //	    *dstI=*srcI;
      //	    ++dstI;
      //	  }
      //	}
      //	AS_confirmed_edges.resize(dstI-AS_confirmed_edges.begin());
      //	nextoverlapcompress=AS_readpool.size()-unused+overlapcompressstepping;
      //	cout << "done.\n";
      //	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      //}

      Contig::setIDCounter(numcontigs);
      // TODO: change wrt multiple MIRAparams
      Contig buildcon(&AS_miraparams, AS_readpool);

      //std::vector<int8> tmpused=AS_used_ids;

#ifdef CLOCK_STEPS2
      cout << "Timing BFC prelim1: " << diffsuseconds(tv) << endl;
#endif

      CEBUG("bfc 3\n");
      auto bbContigI=AS_bbcontigs.cbegin();
      if(numcontigs <= AS_bbcontigs.size()){
	advance(bbContigI,numcontigs-1);
      } else {
	bbContigI=AS_bbcontigs.end();
      }

      CEBUG("bfc 4\n");
      // set newreptmarked to true just to get into the for pass
      //  value is changed within pass
      bool newreptmarked=true;
      bool wasovercalledited=false;
      bool wasmajorovercalledited=false;

      bool mastermayeditovercalls=ed_params.ed_mira_automatic_contic_editing; // may be changed in one of the iter
      bool mayeditovercalls=mastermayeditovercalls; // recomputed after each contig build, just default init here

      // build a contig, repetitively until
      //  ... maximum number of iterations has passed
      //  ... or no now repeats were marked
      //  ... or no 454 edits happened in last iteration

      // Note: if we´re in last pass, iterate at least once if
      //  454 edits were made! (disregarding as_fixparams.as_numrmbbreakloops,
      //  maxiter gets adapted in loop then

      uint32 maxiter=as_fixparams.as_numrmbbreakloops;

      // Note: assemblies with a lot of passes (>=4) will
      //  get only one loop in the first pass. The reason: first pass already
      //  discovers a lot of repeats that are dealt with later on
      //  in makeAlignments(). It's faster to let makeAlignments() deal
      //  with "discovering" and marking reads than to loop here.

      if(as_fixparams.as_mark_repeats){
	if(as_fixparams.as_numpasses >= 4 && passnr==1){
	  maxiter=1;
	}
	// same thing: passes >= 6, pass 2 leads to maxiter = 2
	if(as_fixparams.as_numpasses >= 6 && passnr==2
	  && maxiter > 2){
	  maxiter=2;
	}
      }

      std::list<Contig::pbdse_t> pbdsev;

      bool markrepeatsduringstore=true;
      bool contignotok=false;
      bool continueiter=true;

#ifdef CLOCK_STEPS2
      gettimeofday(&tv,nullptr);
#endif
      buildcon.discard();
#ifdef CLOCK_STEPS2
      cout << "Timing BFC discard con: " << diffsuseconds(tv) << endl;
#endif

      std::string basename_forextsave;
      {
	std::ostringstream ostr;
	ostr << dir_params.dir_tmp << "/miratmp.pass_" << passnr << "_cb" << numcontigs << "_";

	basename_forextsave=ostr.str();
      }

      bool contigtrimmed=false;
      CEBUG("bfc 5\n");
      //for(uint32 iter=0; iter < maxiter && continueiter; ++iter){
      uint32 iter=0;
      do{
	std::string basename_forextsave_iter;
	{
	  std::ostringstream ostr;
	  ostr << basename_forextsave << "i" << iter << "_";

	  basename_forextsave_iter=ostr.str();
	}

	CEBUG("bfc 6/"<<iter << '\n');
	if(iter>0) bfc_sanityCheckASUSEDIDS(trackingunused,numcontigs);

	if(AS_hasbackbones
	   && passnr>=as_fixparams.as_startbackboneusage_inpass
	   && bbContigI != AS_bbcontigs.end()) {

#ifdef CLOCK_STEPS2
	  gettimeofday(&tv,nullptr);
#endif

	  // if we're in an iter>0, the contig is not empty
	  // as, we're re-initialising with fresh bb contig,
	  //  we need to do some housekeeping
	  if(buildcon.getContigLength()>0){
	    // release all reads
	    //Contig::setCoutType(Contig::AS_TEXT);
	    //cout << "Throwing away\n" << buildcon;
	    auto & cr=buildcon.getContigReads();
	    for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	      if(pcrI.getORPID()>=0){
		BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"bbrebuild ! AS_used_ids[...] ???");
		AS_used_ids[pcrI.getORPID()]=0;
		++trackingunused;
	      }
	    }
	  }

	  // new contig: initialise what's needed
	  buildcon=*bbContigI;
	  buildcon.setContigID(numcontigs);

	  {
	    //Contig::setCoutType(Contig::AS_TEXT);
	    //cout << "registering \n" << buildcon;

	    // track the reads
	    auto & cr=buildcon.getContigReads();
	    for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	      if(pcrI.getORPID()>=0){
		BUGIFTHROW(AS_used_ids[pcrI.getORPID()],"register AS_used_ids[" << pcrI.getORPID() << "]==" << static_cast<uint16>(AS_used_ids[pcrI.getORPID()]) << " ???");
		AS_used_ids[pcrI.getORPID()]=1;
		//--trackingunused;
	      }
	    }
	  }

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC copy bbcon: " << diffsuseconds(tv) << endl;
	  gettimeofday(&tv,nullptr);
#endif

	  // re-initialising the baselocks is necessary as reads might
	  //  have got new SRMc/WRMc tags in previous iterations
	  //  these are not known in the initial backbone contig, so
	  //  they must be made known
	  // TODO: 11.10.2012 not sure whether still good
	  buildcon.initialiseBaseLocks();

	  // tell contig to use backbone characters when possible for
	  //  tmp consensus
	  // TODO: make configurable?
	  buildcon.useBackbone4TmpConsensus(true);

	  // and by default try to merge short reads (Solexa, SOLiD)
	  buildcon.mergeNewSRReads(true);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC bbsetup remain: " << diffsuseconds(tv) << endl;
#endif

	} else {
	  // do nothing
	}

	CEBUG("bfc 7/"<<iter << '\n');

	buildcon.setContigNamePrefix(
	  AS_miraparams[0].getContigParams().con_nameprefix);

	if(iter==0){
	  cout << "Building new contig " << numcontigs;
	  if(buildcon.getNumBackbones()) cout << " from backbone " << buildcon.getContigName();
	  cout << endl;
	}else{
	  cout << "Rebuilding contig " << numcontigs << endl;
	}

	//Contig::setCoutType(Contig::AS_DEBUG);
	//cout << con;

	if(!AS_coverageperseqtype.empty()){
	  cout << "Setting contig coverage targets to: ";
	  for(uint8 ii=0; ii<AS_coverageperseqtype.size(); ii++){
	    cout << '\t' << AS_coverageperseqtype[ii];
	  }
	  cout << endl;
	  buildcon.setContigCoverageTarget(AS_coverageperseqtype);
	}

	if(as_fixparams.as_dateoutput) dateStamp(cout);

	cout << "Unused reads: " << trackingunused << endl;
	cout.flush();

	CEBUG("bfc 8/"<<iter << '\n');

	bfc_callPathfinder(passnr,iter,trackingunused,shouldmovesmallclusterstodebris,
			   buildcon,qaf);

	CEBUG("bfc 9/"<<iter << '\n');

	if(buildcon.getNumReadsInContig()==0){
	  cout << "WARNING WARNING WARNING! no reads in contig?!?!?!" << endl;
	  trackingunused=0;
	  break;
	}

	bool contigmeetsrequirements=bfc_checkIfContigMeetsRequirements(buildcon);
	if(!contigmeetsrequirements){

	  // we're going to throw away this contig
	  // however, save the small overlaps or else they might be completely lost!
#if 1
#ifdef CLOCK_STEPS2
	  gettimeofday(&tv,nullptr);
#endif
	  bfc_savePersistentSmallOverlaps(buildcon,passnr,fout_povl);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC persistent small overlaps1: " << diffsuseconds(tv) << endl;
#endif
#endif

	  newreptmarked=false;

	  // this contig won't be taken, but contig numbering would
	  //  go up in next loop.
	  // Countermeasure: decrease contig number by one to have the
	  //  loop increase it afterwards, effectivly re-using same contig number
	  --numcontigs;

	  uint32 numdebs=0;
	  for(auto crI=buildcon.getContigReads().begin();crI!=buildcon.getContigReads().end();++crI){
	    if(crI.getORPID() >= 0){
	      BUGIFTHROW(AS_isdebris[crI.getORPID()]>0,"Ooooops, read is already debris? " << crI.getORPID() << " " << AS_readpool[crI.getORPID()].getName() << endl);
	      BUGIFTHROW(!AS_used_ids[crI.getORPID()]>0,"Debris from an unused read in contig? " << crI.getORPID() << " " << AS_readpool[crI.getORPID()].getName() << endl);
	      ++numdebs;
	      AS_isdebris[crI.getORPID()]=DEBRIS_TINYCONTIG;
	    }
	  }
	  cout << "Discarding" << endl;
	  buildcon.discard();

	  cout << "\nContig does not meet requirement of minimum reads per contig."
	    "\nMoved " << numdebs << " reads to debris." << endl;

	  // break out of the iter loop
	  break;
	}

	cout << "\n\nFinished building." << endl;

	if(buildcon.getContigLength()>5 && buildcon.getNumReadsInContig()>1){
	  if(as_fixparams.as_put_asswithmira_tags){
	    buildcon.addTagToConsensus(0,
				       4,
				       '=',
				       "MIRA",
				       "Assembled with MIRA",
				       false);
	  }
	}

	try {
	  if(buildcon.getNumReadsInContig()>1){
	    if(as_fixparams.as_dateoutput) dateStamp(cout);
	    if(buildcon.getContigLength()>100000){
	      cout << "Calculating statistics (this may take a while)." << endl;
	    }

#ifdef CLOCK_STEPS2
	    gettimeofday(&tv,nullptr);
#endif
	    // see whether we can define new multicopies
	    if(as_fixparams.as_automatic_repeat_detection){
	      buildcon.analyseReadCoverage(AS_maxcoveragereached,
				      AS_multicopies,
				      AS_coverageperseqtype);
	    }
#ifdef CLOCK_STEPS2
	    cout << "Timing BFC analysereadcov: " << diffsuseconds(tv) << endl;
	    gettimeofday(&tv,nullptr);
#endif

	    buildcon.setCoutType(Contig::AS_TEXT);
	    buildcon.dumpStats(cout);

#ifdef CLOCK_STEPS2
	    cout << "Timing BFC cout constats: " << diffsuseconds(tv) << endl;
#endif
	    if(as_fixparams.as_dateoutput) dateStamp(cout);
	  }
	}
	catch (Notify n) {
	  n.setGravity(Notify::FATAL);
	  n.handleError(n.tif);
	}
	catch (...) {
	  cerr << "Darn, error with that contig. See darn.fasta.\n";
	  Read::setCoutType(Read::AS_CLIPPEDFASTA);
	  for(auto & cre : buildcon.getContigReads()) {
	    cout << cre;
	  }
	  abort();
	}

	CEBUG("bfc 10/"<<iter << '\n');

	// saving pre-edit
	saveExtTmpContig(buildcon,(basename_forextsave_iter+"pre"));

	bool wasedited=false;
	newreptmarked=false;

#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif


	bool goforit=true;
/*
  PacBio trials with strobed reads ... MIRA throws on real, non-strobed data
  furthermore, not sure whether this is really doable (accuracy ~85% *sigh*)

	pbdsev.clear();
	if(buildcon.hasSeqTypeData(ReadGroupLib::SEQTYPE_PACBIOHQ) && !(lastpass && iter==maxiter-1)){
	  uint32 maxcorrect=buildcon.createPacBioDarkStrobeEdits(pbdsev);

	  CEBUG("PacBio maxcorrect: " << maxcorrect << endl);

	  if(maxcorrect>5) {
	    goforit=false;
	    wasmajorovercalledited=true;
	  }

	  // another iteration is needed to get a valid contig
	  //if(maxiter<2) maxiter=2;
	}
*/

	mayeditovercalls= mastermayeditovercalls & buildcon.hasEditableOvercallData();

	if(goforit){
	  // handling of tricky overcalls (normally 454 & Ion)
	  std::vector<bool> readsmarkedsrm;

	  // mark areas with tricky overcalls
	  // in those areas, no repeat marker may be set
	  // do this always
	  if(mayeditovercalls && buildcon.getNumReadsInContig() >1){
	    cout << "P " << passnr << ", marked " << buildcon.editTrickyOvercalls(true,false,readsmarkedsrm) << " reads.\n";
	  }

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC edit tricky1: " << diffsuseconds(tv) << endl;
	  gettimeofday(&tv,nullptr);
#endif

	  markrepeatsduringstore=true;
	  if(as_fixparams.as_mark_repeats
	     && !as_fixparams.as_mark_repeats_onlyinresult
	     && buildcon.getNumReadsInContig() >1){
	    Contig::repeatmarker_stats_t contigrepstats;
	    newreptmarked=markRepeats(buildcon, readsmarkedsrm, contigrepstats);
	    AS_bfcstats[buildcon.getLongRepeatStatus()].numnewsrm+=contigrepstats.numSRMs;
	    foundSRMs|=newreptmarked;
	    markrepeatsduringstore=false;
	  }

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC mark reps: " << diffsuseconds(tv) << endl;
	  gettimeofday(&tv,nullptr);
#endif

	  CEBUG("bfc 11/"<<iter << '\n');

	  // edit only when no misassembled repeats found
	  //  this is to prevent the editor to try something foolish on
	  //  misassembled things
	  // same applies to overcall editing
	  //

	  wasovercalledited=false;

//	  if(buildcon.getNumReadsInContig() <= 8) mayeditovercalls=false;
//	  cout << "\nXXXXXXXXXXXXXX " << newreptmarked
//	       << " " << mayeditovercalls
//	       << " " << buildcon.getNumReadsInContig() << endl;
	  if(!newreptmarked && mayeditovercalls && buildcon.getNumReadsInContig() >1){
	    uint32 nummarks=buildcon.editTrickyOvercalls(false,false,readsmarkedsrm);
	    cout << "Edited " << nummarks << " reads.\n";
	    AS_bfcstats[buildcon.getLongRepeatStatus()].numeditovercall+=nummarks;
	    wasovercalledited=nummarks>0;

	    // major edit is edits in >=5% of the reads
	    wasmajorovercalledited=(nummarks>0 && nummarks >= buildcon.getNumReadsInContig()/20);

	    if(wasovercalledited){
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      cout << "Deleting superfluous gap columns (1) ... ";
	      cout.flush();
	      auto tmpnumdel=buildcon.deleteStarOnlyColumns(0, buildcon.getContigLength()-1);
	      cout << "done, deleted " << tmpnumdel << " columns.\n";
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	    }
#ifdef CLOCK_STEPS2
	    cout << "Timing BFC edit tricky2: " << diffsuseconds(tv) << endl;
	    gettimeofday(&tv,nullptr);
#endif
	  }

	  // BaCh 29.04.2011
	  // However, editSingleDiscrepancyNoHAFTag() is conservative enough to allow
	  //  that kind of edits all the time
	  //

	  if(ed_params.ed_mira_automatic_contic_editing
	     && ed_params.ed_kmer_singlets
	     && buildcon.getNumReadsInContig() >4){
	    // editmode: increasingly less conservative (but still conservative enough)
	    uint32 editmode=passnr-1;
	    if(editmode>3) editmode=3;
	    // careful, editmode 3 not suited for EST/RNASeq!
	    if(!AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
	       && editmode>2) {
	      editmode=2;
	    }
	    uint32 numedits=buildcon.editSingleDiscrepancyNoHAFTag(readsmarkedsrm,editmode);
	    cout << "\nEdited " << numedits << " positions.\n";
	    if(numedits>0){
	      AS_bfcstats[buildcon.getLongRepeatStatus()].numedithashfreq+=numedits;
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      cout << "Deleting superfluous gap columns (2) ... ";
	      cout.flush();
	      auto tmpnumdel=buildcon.deleteStarOnlyColumns(0, buildcon.getContigLength()-1);
	      cout << "done, deleted " << tmpnumdel << " columns.\n";
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	    }
#ifdef CLOCK_STEPS2
	    cout << "Timing BFC edit single discrepancy, no HAF: " << diffsuseconds(tv) << endl;
	    gettimeofday(&tv,nullptr);
#endif
	  }


/*
  Trials with original, raw PB reads (non-CCS, non-CLR), currently on hold

	  // PacBio LQ edits
	  if(buildcon.getNumReadsInContig() >4 && buildcon.hasSeqTypeData(ReadGroupLib::SEQTYPE_PACBIOLQ)){
	    uint32 numcoledits=0;
	    uint32 numreadedits=0;
	    buildcon.editPBSledgeHammer(readsmarkedsrm,numcoledits,numreadedits);
	    //buildcon.deleteStarOnlyColumns(0,I->getContigLength());
	    if(numreadedits){
	      cout << "PacBio low quality SledgeHammer edits: " << numcoledits << " c-edits, " << numreadedits << " r-edits.\n";
	    }
	  }
*/

	  // if we're in last pass and 454 / pacbio dark strobe edits were made, loop at least
	  //  once (switching off the 454 / pacbio edits for next pass)
	  if(lastpass
	     && wasovercalledited
	     && iter==maxiter-1){
	    if(maxiter>=1) {
	      mastermayeditovercalls=false;
	    }
	    maxiter++;
	  }

	  CEBUG("bfc 12/"<<iter << '\n');

#ifdef CLOCK_STEPS2
	  gettimeofday(&tv,nullptr);
#endif
	  // get rid of all PSHP tags
	  buildcon.deleteTagsInReads(Read::REA_defaulttag_PSHP.identifier);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC delPSHP: " << diffsuseconds(tv) << endl;
#endif

	  CEBUG("bfc 13/"<<iter << '\n');

#ifdef MIRA_HAS_EDIT
	  if(!newreptmarked && ed_params.ed_automatic_contic_editing!=0){
	    cout << "Editing temporary contig: ";
	    if (buildcon.getNumReadsInContig() >1) {
	      cout << endl;
#ifdef CLOCK_STEPS2
	      gettimeofday(&tv,nullptr);
#endif
	      wasedited=true;

	      editContigBack(buildcon, const_cast<EDITParameters &>(eparams));

	      //ScfBuffer::statistics();
	      //ScfBuffer::show();
	      //ScfBuffer::discard();
	      //ScfBuffer::statistics();
	      //ScfBuffer::show();

	      cout << "done.\n";
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      cout << "Deleting superfluous gap columns ... ";
	      cout.flush();
	      buildcon.deleteStarOnlyColumns(0, buildcon.getContigLength()-1);
	      cout << "done.\n";

#ifdef CLOCK_STEPS2
	      cout << "Timing BFC editconback and more: " << diffsuseconds(tv) << endl;
#endif
	      Contig::setCoutType(Contig::AS_TEXT);
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      buildcon.stats(cout);
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	    } else {
	      cout << "(only 1 read in contig, no editing)\n";
	    }
	  }
#endif


	  CEBUG("bfc 14/"<<iter << '\n');

	  // saving again if rept-marked or edited
	  if(newreptmarked || wasedited || wasovercalledited) {
	    if(wasovercalledited) {
	      saveExtTmpContig(buildcon,(basename_forextsave_iter+"post454"));
	    } else {
	      saveExtTmpContig(buildcon,(basename_forextsave_iter+"post"));
	    }
	  } else {
	    if ( as_fixparams.as_output_exttmp_fasta
		 || as_fixparams.as_output_exttmp_ace
		 || as_fixparams.as_output_exttmp_gap4da
		 || as_fixparams.as_output_exttmp_caf) {
	      cout << "No edit and no new repeat found, not saving extra temporary contig again.\n";
	    }
	  }
	}


	CEBUG("bfc 15/"<<iter << '\n');

	CEBUG("bfc 16/"<<iter << '\n');

	// Transfer all the reads fron the new contig into readpool
	//  if no misassembly was detected (or 454 editing happened)
#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	if(!newreptmarked || wasovercalledited || !pbdsev.empty()) {
	  transferContigReadsToReadpool(buildcon, pbdsev, passnr);

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC rp transfer: " << diffsuseconds(tv) << endl;
#endif
	} else {
	  cout << " Transfering read tags to readpool." << endl;
	  transferContigReadTagsToReadpool(buildcon, bbContigI);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC crtag2rp transfer: " << diffsuseconds(tv) << endl;
#endif
	}
	cout << "Done." << endl;

	CEBUG("bfc 17/"<<iter << '\n');

	contigtrimmed=false;
	if(passnr>1){
	  // reason: on highly repetitive data, reads trimmed away have SRMr tags,
	  // leading to long assembly times in later contigs. Not trimming away the
	  // reads in first pass makes the 1st pass finish quicker and arrive to the
	  // Smith-Waterman repeat disambiguation stage (SRMr->CRMr tags)
	  contignotok=bfc_trimDenovoIfNecessary(buildcon,foundSRMs,basename_forextsave_iter,trackingunused);
	  if(contignotok) ++AS_bfcstats[buildcon.getLongRepeatStatus()].numdisassemblies;
	  contigtrimmed=true;
	}

	++iter;

	continueiter=(iter < maxiter) && (newreptmarked | wasmajorovercalledited | contignotok);

	CEBUG("I have newreptmarked " << newreptmarked << " wasmajorovercalledited " << wasmajorovercalledited << " contignotok " << contignotok << "\tcontinueiter: " << continueiter << endl);

	bool nukeexistingcontig=false;
	if(newreptmarked){
	  cout << "Identified misassembled reads in contig.\n";
	}
	if(wasmajorovercalledited){
	  cout << "Had many overcall edits in reads.\n";
	  nukeexistingcontig=true;
	}
	if(!contigmeetsrequirements){
	  BUGIFTHROW(true,"We should never be here: contigmeetsrequirements");
	}

	if(continueiter && buildcon.getNumBackbones()>0){
	  cout << "Has backbones, no additional iteration allowed.\n";
	  continueiter=false;
	}

	if(lastpass) nukeexistingcontig=false;

	if(continueiter){
	  cout << "Need to loop contig building\n";
	  if(nukeexistingcontig){
	    cout << "Iteration will nuke contig built so far.\n";
	    //auto & cr=buildcon.getContigReads();
	    //for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	    //  if(pcrI.getORPID()>=0){
	    //	BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"nec ! AS_used_ids[...] ???");
	    //	AS_used_ids[pcrI.getORPID()]=0;
	    //	trackingunused+=1;
	    //  }
	    //}

	    for(auto & rid : qaf.getRIDsKnownInContig()){
	      if(rid>=0
		 && !AS_readpool[rid].isBackbone()
		 && !AS_readpool[rid].isRail()){
		if(AS_used_ids[rid]){
		  trackingunused+=1;
		  AS_used_ids[rid]=0;
		}
	      }
	    }

	    buildcon.discard();
	  }else{
	    cout << "Iteration will keep contig built so far.\n";
	  }
	}

      }while(continueiter);

      bfc_sanityCheckASUSEDIDS(trackingunused,numcontigs);

      // no contig? Then it was discarded, restart building one completely anew
      if(buildcon.getNumReadsInContig() == 0) continue;

      // this here to handle cases repeats should be marked only in results
      // TODO: check whether not to remove this parameter / option at all
      if(pbdsev.empty()
	 && lastpass
	 && as_fixparams.as_mark_repeats
	 && as_fixparams.as_mark_repeats_onlyinresult
	 && buildcon.getNumReadsInContig() >1){
	CEBUG("bfc 18\n");
#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	std::vector<bool> dummy;

	Contig::repeatmarker_stats_t contigrepstats;
	newreptmarked=markRepeats(buildcon, dummy,contigrepstats);
	AS_bfcstats[buildcon.getLongRepeatStatus()].numnewsrm+=contigrepstats.numSRMs;
	foundSRMs|=newreptmarked;
	markrepeatsduringstore=false;
#ifdef CLOCK_STEPS2
	cout << "Timing BFC markrep during store: " << diffsuseconds(tv) << endl;
#endif
      }

      if(buildcon.getNumReadsInContig() >0){
	if(newreptmarked){
	  cout << "\nAccepting probably misassembled contig, ";
	  if(contigtrimmed) {
	    cout << " but kept only best, non-problematic part.\n";
	  }else{
	    cout << " keeping as is.\n";
	  }
	}

	CEBUG("bfc 19\n");

#if 1
#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	bfc_savePersistentSmallOverlaps(buildcon,passnr,fout_povl);
#ifdef CLOCK_STEPS2
	cout << "Timing BFC persistent small overlaps2: " << diffsuseconds(tv) << endl;
#endif
#endif

	CEBUG("bfc 20\n");

	// If Illumina data present, see whether we find the contig to contain
	//  a sizeable amount of PhiX174 sequence. If yes, rename that contig.
#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	if(!AS_miraparams[0].getNagAndWarnParams().nw_warn_x174prefix.empty()
	   && buildcon.getNumReadsInContig()>1
	   && buildcon.getContigLength() > 50){
	  Read tmpc;
	  tmpc.disallowAdjustments();
	  tmpc.setName("dummy");
	  auto & tmpcons=buildcon.getTmpConsensus();
	  //cout << tmpcons << endl;
	  tmpc.setSequenceFromString(tmpcons);

	  auto numbaithits=AS_phix174hashstatistics.checkBaitHit(tmpc,false,0);
	  decltype(numbaithits) minhits=(buildcon.getContigLength()-30) / 3;   // 1/3rd of contig? -> phiX
	  if(numbaithits>=minhits){
	    buildcon.setContigNamePrefix(AS_miraparams[0].getNagAndWarnParams().nw_warn_x174prefix
					 +AS_miraparams[0].getContigParams().con_nameprefix);
	  }
	}
#ifdef CLOCK_STEPS2
	cout << "Timing BFC phix174 con: " << diffsuseconds(tv) << endl;
#endif

#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	bfc_storeContig(buildcon,numcontigs,markrepeatsduringstore,passnr,lastpass,true);
#ifdef CLOCK_STEPS2
	cout << "Timing BFC store con: " << diffsuseconds(tv) << endl;
#endif

	{
	  if(buildcon.getContigLength()>=5000){
	    cout << "Contig coverage analysis ";
	    if(buildcon.getContigLength()>=100000){
	      cout << "(this may take a while) ";
	    }
	    cout << "... "; cout.flush();
	    const Contig::constats_t & cs=buildcon.getStats();
	    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	      covperstpercon[st].push_back(static_cast<uint32>(cs.avg_covperst[st]+.5));
	      //cout << "cccccccccccc: " << covperstpercon[st].back() << endl;
	    }
	    covtotalpercon.push_back(static_cast<uint32>(cs.avg_coverage+.5));
	  }
	}

	// did the caller to assemble() ask for a callback for each contig built?
	// if we're on the last past, then call it
	if(lastpass && AS_contigbuilt_callbackfunc!=nullptr){
	  (*AS_contigbuilt_callbackfunc)(buildcon, AS_readpool);
	}

	if(as_fixparams.as_spoilerdetection) {
	  if(!as_fixparams.as_spdetect_lastpassonly
	     || passnr==as_fixparams.as_numpasses-1) {
	    huntSpoilSports(buildcon);
	  }
	}
      }

      // if we were in mapping mode, and on the last contig
      //  move all remaining reads to debris

      if(bbContigI != AS_bbcontigs.end()){
	++bbContigI;  // don't bother to decrease after using this ... will be re-init in next loop anyway
	if(bbContigI == AS_bbcontigs.end()){
	  uint32 adddebris=0;
	  for(size_t uid=0; uid<AS_used_ids.size(); ++uid){
	    if(AS_used_ids[uid]==0){
	      ++adddebris;
	      --trackingunused;
	      AS_isdebris[uid]=DEBRIS_NOTMAPPED;
	      AS_used_ids[uid]=1;
	    }
	  }
	  cout << "Last backbone mapped, not building new ones.\nMoved " << adddebris << " remaining reads to debris.\n";
	}
      }


    }  // if(trackingunused>0) ...


    //unused=0;
#ifdef CLOCK_STEPS2
    cout << "Timing BFC loop total: " << diffsuseconds(tvloop) << endl;
#endif
  }   //for(;trackingunused>0; ++numcontigs) ...


  // force a sanity check
  bfc_sanityCheckASUSEDIDS(trackingunused,0);

  // unused reads? we jumped out of contig creation
  // define unused reads as debris

  if(trackingunused>0){
    for(size_t uid=0; uid<AS_used_ids.size(); ++uid){
      if(AS_used_ids[uid]==0){
	AS_isdebris[uid]=DEBRIS_ABORTEDCONTIGCREATION;
	AS_used_ids[uid]=1;
      }
    }
  }

  // if the user wants "singlets", this will save the
  //  non-completely clipped debris as contigs
  for(uint32 st=1; st<ReadGroupLib::getNumSequencingTypes(); ++st){
    if(AS_miraparams[st].getAssemblyParams().as_savesimplesingletsinproject){
      saveDebrisAsContigs(numcontigs,passnr,lastpass);
      break; // yeah, but do this only once, will you?
    }
  }

  // Adapt debris
  // - DEBRIS_NOOVERLAP with MNRr tag changed to DEBRIS_MASKEDNASTYREPEAT
  // - DEBRIS_NOOVERLAP with megahub flag changed to DEBRIS_MEGAHUB
  //
  // Also create warning if needed
  {
    uint32 numnasty=0;
    uint32 nummega=0;
    for(size_t uid=0; uid<AS_isdebris.size(); ++uid){
      if(AS_isdebris[uid]==DEBRIS_NOOVERLAP){
	if(AS_readpool[uid].hasTag(Read::REA_tagentry_idMNRr)){
	  AS_isdebris[uid]=DEBRIS_MASKEDNASTYREPEAT;
	  ++numnasty;
	}else if(AS_readpool[uid].hasTag(Read::REA_tagentry_idHAF7)){
	  AS_isdebris[uid]=DEBRIS_MASKEDHAF7REPEAT;
	}else if(AS_readpool[uid].hasTag(Read::REA_tagentry_idHAF6)){
	  AS_isdebris[uid]=DEBRIS_MASKEDHAF6REPEAT;
	}else if(!AS_skimmegahubs.empty() && AS_skimmegahubs[uid]==1){
	  AS_isdebris[uid]=DEBRIS_MEGAHUB;
	  ++nummega;
	}
      }
    }
    if(nummega){
      std::string wstr("There are ");
      auto dummypercent=static_cast<uint32>((static_cast<double>(100.0/AS_readpool.size())*nummega)+.5);
      wstr+=boost::lexical_cast<std::string>(nummega)+" of your reads (~" + boost::lexical_cast<std::string>(dummypercent) + "%) declared as megahubs (too many possible connections to other reads) and where the SKIM fast overlap finder did not return enough results.\nIf you want those reads to be found (and subsequently assembled), you need to either\n- switch off filtering of megahubs (-SK:fmh=no)\n- increase the size of the megahub cap (-SK:mhc=...)\nBeware! Megahubs are the first defense mechanism of MIRA to keep assembly times reasonable.";
      if(nummega<1000) {
	wstr+=" Fortunately, the number of meguhbs is low enough to safely switch off filtering.";
      }else{
	wstr+=" With this number of megahubs, the assembly is expected to ";
	if(nummega<10000){
	  wstr+="take a couple of minutes longer per pass, switching off filtering should be OK.";
	}else if(nummega<100000){
	  wstr+="take some 30 minutes to 1 hour longer per pass when switching off filtering, if you want to wait for that, OK.";
	}else if(nummega<500000){
	  wstr+="take some 1 to 5 hours longer per pass when switching off filtering, that can get long.";
	}else{
	  wstr+="take >= 5 hours to several days longer per pass when switching off filtering, only for die-hards.";
	}
      }
      AS_warnings.setWarning("UNASSEMBLED_MEGAHUBS",2,"Megahub reads landed in debris file",wstr);
    }else{
      AS_warnings.clearWarning("UNASSEMBLED_MEGAHUBS");
    }
    if(numnasty){
      std::string wstr("There are ");
      auto dummypercent=static_cast<uint32>((static_cast<double>(100.0/AS_readpool.size())*numnasty)+.5);
      wstr+=boost::lexical_cast<std::string>(numnasty)+" of your reads (~" + boost::lexical_cast<std::string>(dummypercent) + "%) containing so many bases masked as nasty repeats that the SKIM fast overlap finder did not return enough results.\nIf you want those reads to be found (and subsequently assembled), you need to either\n- switch off masking of nasty repeats (-HS:mnr=no)\n- increase the the nasty repeat ratio (-HS:nrr=...) (consult the hash statistics in the output log of MIRA to find a good one)";
      if(!AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	wstr+="\n- as you do EST/RNASeq assembly, lossless digital normalisation (-HS:ldn=yes) may also be an option";
	if(AS_miraparams[0].getHashStatisticsParams().hs_apply_digitalnormalisation) {
	  wstr+=", but you already had it on";
	}
      }
      wstr+="\nBeware! Masking of nasty repeats is the second defense mechanism of MIRA (after megahubs) to keep assembly times reasonable. Changing the defaults will increase the time needed.";
      if(numnasty<1000) {
	wstr+=" Fortunately, the number of completely masked reads is low enough to safely switch off masking.";
      }else{
	wstr+=" With this number of masked reads, the assembly is expected to ";
	if(numnasty<10000){
	  wstr+="take a couple of minutes longer per pass, switching off masking should be OK.";
	}else if(numnasty<100000){
	  wstr+="take some 30 minutes to 1 hour longer per pass when switching off masking, if you want to wait for that, OK.";
	}else if(numnasty<500000){
	  wstr+="take some 1 to 5 hours longer per pass when switching off masking, that can get long.";
	}else{
	  wstr+="take >= 5 hours to several days longer per pass when switching off masking, only for die-hards.";
	}
      }
      AS_warnings.setWarning("UNASSEMBLED_NASTYREPEATS",2,"Reads with masked nasty repeats landed in debris file",wstr);
    }else{
      AS_warnings.clearWarning("UNASSEMBLED_NASTYREPEATS");
    }
  }

  cout << "\nBuildstats - RM positions        :\t"
       << AS_bfcstats[0].numnewsrm << '\t' << AS_bfcstats[1].numnewsrm
       << "\nBuildstats - overcall edits      :\t"
       << AS_bfcstats[0].numeditovercall << '\t' << AS_bfcstats[1].numeditovercall
       << "\nBuildstats - hash edits          :\t"
       << AS_bfcstats[0].numedithashfreq << '\t' << AS_bfcstats[1].numedithashfreq
       << "\nBuildstats - contig disassemblies:\t"
       << AS_bfcstats[0].numdisassemblies << '\t' << AS_bfcstats[1].numdisassemblies
       << endl << endl;

  analyseTemplateGuesses();
  if(lastpass && passnr==1){
    analyseTemplateGuesses();
  }

  if(lastpass) {
    saveDebrisList();
  }else{
    saveDebrisList(passnr, "", "_pass");
  }

  // reads that are debris or singlets apparently need every chance
  //  they can get to align, therefore subsequent passes should not
  //  reduce the overlaps
  // Change: for Solexa and Ion reads ... sorry, we'll just miss out
  {
    for(uint32 i=0; i<AS_needalloverlaps.size(); i++){
      if(AS_isdebris[i] && AS_readpool[i].getSequencingType() != ReadGroupLib::SEQTYPE_SOLEXA
	 && AS_isdebris[i] && AS_readpool[i].getSequencingType() != ReadGroupLib::SEQTYPE_IONTORRENT) AS_needalloverlaps[i]=true;
    }

    std::vector<int32> contigids;
    for(const auto & cle : AS_contigs){
      if(cle.getNumReadsInContig()==1){
	cle.getReadORPIDsAtContigPosition(contigids,0,0);
	AS_needalloverlaps[contigids[0]]=true;
      }
    }
  }

  // calculate median of average contig coverage
  // TODO: this is calculated only once because of problems
  //  with decreasing average/median in subsequent passes
  //  see if it can be improved
  if(as_fixparams.as_uniform_read_distribution
     && passnr+1>=as_fixparams.as_urd_startinpass
     && !covperstpercon[0].empty() && AS_coverageperseqtype.empty()){
    AS_coverageperseqtype.clear();

    cout << "Setting coverage analysis values for uniform read distribution:\n";
    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
      // use non-parallel sort, this is a small thing (right?)
      mstd::ssort(covperstpercon[st]);
      AS_coverageperseqtype.push_back(covperstpercon[st][covperstpercon[st].size()/2]);
      cout << '\t' << ReadGroupLib::getNameOfSequencingType(st) << " coverage:\t" << AS_coverageperseqtype.back() << '\n';
    }
    mstd::ssort(covtotalpercon);
    AS_coveragetotal=covtotalpercon[covtotalpercon.size()/2];
  }

  AS_steps[ASCONTIGSOK]=1;

  AS_used_ids.clear();

  //  saveAsCAF();

  FUNCEND();

  return foundSRMs;
}
#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::saveDebrisAsContigs(uint32 & numcontigs, const int32 passnr, const bool lastpass)
{
  std::vector<bool> savesinglet{false};

  for(uint32 st=1; st<ReadGroupLib::getNumSequencingTypes(); ++st){
    savesinglet.push_back(AS_miraparams[st].getAssemblyParams().as_savesimplesingletsinproject);
  }

  Contig buildcon(&AS_miraparams, AS_readpool);
  buildcon.setVerbose(false);

  cout << "Saving selected debris as singlets:\n";
  ProgressIndicator<int64> P(0, AS_isdebris.size());
  for(size_t uid=0; uid<AS_isdebris.size(); ++uid){
    P.increaseprogress();
    if(AS_isdebris[uid]==DEBRIS_NOOVERLAP
       && savesinglet[AS_readpool[uid].getSequencingType()]){
      AS_isdebris[uid]=DEBRIS_NOTDEBRIS;
      buildcon.discard();
      buildcon.setContigNamePrefix(AS_miraparams[0].getContigParams().con_nameprefix);
      buildcon.setContigID(numcontigs);
      buildcon.resetContigName();
      buildcon.addFirstRead(uid,1);
      bfc_storeContig(buildcon,numcontigs,false,passnr,lastpass,false);
      ++numcontigs;
    }
  }
  P.finishAtOnce();
  cout << "\n";
}


/*************************************************************************
 *
 * numexpected: number of reads expected to be unused
 * numcontigs is a simple "don't do that all the time" switch, but every 100th contig
 *
 *************************************************************************/

void Assembly::bfc_sanityCheckASUSEDIDS(uint32 numexpected, uint32 numcontigs)
{
  // REMOVEME from production code once stable

  FUNCSTART("void Assembly::bfc_sanityCheckASUSEDIDS(uint32 numexpected)");

  if(AS_miraparams[0].getSpecialParams().mi_extra_flag1 || numcontigs%100==0){
#ifdef CLOCK_STEPS2
    timeval tv;
    gettimeofday(&tv,nullptr);
#endif

    uint32 count=0;
    for(auto a : AS_used_ids) if(!a) ++count;

#ifdef CLOCK_STEPS2
    cout << "Timing BFC unused: " << diffsuseconds(tv) << endl;
    cout << "CUnused: " << count << endl;
    cout << "TUnused: " << numexpected << endl;
    cout << "AS_used_ids.size(): " << AS_used_ids.size() << endl;
#endif

    BUGIFTHROW(count!=numexpected,"count " << count << " != expected " << numexpected);
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::bfc_callPathfinder(const int32 passnr, const uint32 iter, uint32 & trackingunused, bool shouldmovesmallclusterstodebris,Contig & buildcon, PPathfinder & qaf)
{
  FUNCSTART("void Assembly::bfc_callPathfinder(const int32 passnr, const uint32 iter, uint32 trackingunused, bool shouldmovesmallclusterstodebris,Contig & buildcon, PPathfinder & qaf)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

#ifdef CLOCK_STEPS2
  timeval tv;
#endif

  // in which mode are we? de-novo or mapping?
  bool assemblymode_mapping=false;
  if(AS_hasbackbones
     && passnr >= AS_miraparams[0].getAssemblyParams().as_startbackboneusage_inpass){
    assemblymode_mapping=true;
  }

  cout << iter << "\tKnown 1: " << qaf.getRIDsKnownInContig().size() << endl;
  if(iter==0
     || buildcon.getContigReads().size()==0
     || buildcon.getContigLength()==0){
    CEBUG("iter 0, PF init new contig\n");
    qaf.prepareForNewContig(buildcon);
  }else if(assemblymode_mapping){
    CEBUG("iter n, mapping, PF resync contig\n");
    qaf.resyncContig();
  }
  cout << "Known 2: " << qaf.getRIDsKnownInContig().size() << endl;

  CEBUG("assemblymode_mapping: " << assemblymode_mapping << '\n');

  bool wantbootstrap=false;
  //if(assemblymode_mapping){
  //  for(uint32 st=0; st<AS_seqtypespresent.size(); ++st){
  //    if(AS_seqtypespresent[st] && AS_miraparams[st].getAssemblyParams().as_backbone_bootstrapnewbackbone) wantbootstrap=true;
  //  }
  //}
  if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA] && AS_miraparams[ReadGroupLib::SEQTYPE_SOLEXA].getAssemblyParams().as_backbone_bootstrapnewbackbone) wantbootstrap=true;

  if(assemblymode_mapping && AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]){
    if(wantbootstrap){
      CEBUG("mapping & solexa bootstrap\n");
      bfc_cp_mapWithSolexa(buildcon,qaf);
      buildcon.clipHomoPolyWithGapsAtEnds();


      //{
      //	std::ofstream fout("aftermap.maf");
      //	fout << "@Version\t2\t0\n";
      //	Contig::setCoutType(Contig::AS_MAF);
      //	Contig::setOutputRails(true);
      //	ReadGroupLib::dumpAllReadGroupsAsMAF(fout);
      //	fout << buildcon;
      //	Contig::setOutputRails(false);
      //}

      cout << "Looking at what to throw away ... "; cout.flush();
      priv_removePotentiallyWrongBaseInserts(buildcon);

      cout << "stripping ... ";cout.flush();
      buildcon.stripToBackbone();
      for(auto & rid : qaf.getRIDsKnownInContig()){
	if(rid>=0
	   && !AS_readpool[rid].isBackbone()
	   && !AS_readpool[rid].isRail()){
	  AS_used_ids[rid]=0;
	}
      }
      // chomp is needed here:
      //  removing reads from the contig may leave overhangs at the ends which are not
      //  covered by backbone. The alignment routines in Contig::addRead_wrapped() will not
      //  cope well with that as the calculation of the expected offset is then wrong
      //  (correct for indels *in* the reference, but not made for pseudo-indels at the
      //  ends of the contig)
      cout << "done, chomping ... ";cout.flush();
      buildcon.chompFront(-1);
      buildcon.chompBack(-1);

      cout << "done, synching ... ";cout.flush();
      qaf.resyncContig();
      cout << "done\n";
    }

    bfc_cp_mapWithSolexa(buildcon,qaf);
    buildcon.clipHomoPolyWithGapsAtEnds();

    // TODO: hack until there's a routine that only clears tags set by the
    //  makeIntelligentConsensus() functions.
    buildcon.clearConsensusTags();

    CEBUG("TU before " << trackingunused << endl);
    CEBUG("RIDs known " << qaf.getRIDsKnownInContig().size() << endl);
    trackingunused-=qaf.getRIDsKnownInContig().size();
    CEBUG("TU after " << trackingunused << endl);

  }else if(assemblymode_mapping){
    CEBUG("bfccp2" << endl);
    qaf.map();
    buildcon.coutAddReadTimings();
    cout << "Known 3: " << qaf.getRIDsKnownInContig().size() << endl;
    trackingunused-=qaf.getRIDsKnownInContig().size();
  }else{
#ifdef CLOCK_STEPS2
    gettimeofday(&tv,nullptr);
#endif

    CEBUG("use general pathfinder: " << AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms << '\n');

    bfc_sanityCheckASUSEDIDS(trackingunused, buildcon.getContigID());

    auto hacktracking=buildcon.getNumReadsInContig();
    if(buildcon.getNumReadsInContig()==0){
      CEBUG("bfccp3" << endl);
    }else{
      CEBUG("no bfccp3" << endl);
	  // so if the contig has some reads already? We're in an iter loop
	  // bump trackingunused by number of reads because after calling PPathfinder,
	  //  they'll be in the number of reads the pathfinder knows of.
    }

    qaf.denovo();

#ifdef CLOCK_STEPS2
    cout << "Timing BFC paf construct: " << diffsuseconds(tv) << endl;
#endif
    buildcon.coutAddReadTimings();

    cout << iter << "\tKnown 3: " << qaf.getRIDsKnownInContig().size() << endl;
    trackingunused+=hacktracking;
    trackingunused-=qaf.getRIDsKnownInContig().size();

    bfc_sanityCheckASUSEDIDS(trackingunused, buildcon.getContigID());

    // Trigger moving small clusters to debris if the pathfinder had to refill its
    //  start cache or when the start cache containes only singlets
    if(shouldmovesmallclusterstodebris
       && (qaf.startCacheRanDry() || qaf.startCacheHasSinglets())){
      cout << "Triggering additional cluster check:";
      if(shouldmovesmallclusterstodebris) cout << " shouldmovesmallclusterstodebris";
      if(qaf.startCacheRanDry()) cout << " startCacheRanDry";
      if(qaf.startCacheHasSinglets()) cout << " startCacheHasSinglets";
      cout << '\n';
      trackingunused-=bfc_moveSmallClustersToDebris();
    }
  }

  bfc_sanityCheckASUSEDIDS(trackingunused, buildcon.getContigID());

  // REMOVEME: paranoia check
  for(auto crI=buildcon.getContigReads().begin();crI!=buildcon.getContigReads().end();++crI){
    if(crI.getORPID() >= 0){
      BUGIFTHROW(!AS_used_ids[crI.getORPID()],"Ooooops, read is in contig but not used? " << crI.getORPID() << " " << AS_readpool[crI.getORPID()].getName() << endl);
    }
  }
}
#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::bfc_cp_mapWithSolexa(Contig & buildcon, PPathfinder & qaf)
{
  FUNCSTART("void Assembly::bfc_cp_mapWithSolexa(Contig & buildcon, PPathfinder & qaf)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  bool hasotherst=false;
  //AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]){
  for(uint8 st=0; st<AS_seqtypespresent.size(); ++st){
    if(st!=ReadGroupLib::SEQTYPE_SOLEXA && AS_seqtypespresent[st]) hasotherst=true;
  }

  CEBUG("HASOTHERST: " << hasotherst << endl);

  if(AS_miraparams[0].getFinalMappingParameters().fm_active){
    // -FM:act is set ... do it the (hopefully better) long and painful way

    auto & fmparams=AS_miraparams[ReadGroupLib::SEQTYPE_SOLEXA].getFinalMappingParameters();

    if(fmparams.fm_maxtotalerrors>0){
      qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_SOLEXA);
      auto mincoel=fmparams.fm_clean_end_dist;
      if (mincoel < 4) mincoel = 4;
      for(uint32 coel=28; coel>=mincoel; coel-=4){
	cout << "Gogo: coel " << coel << endl;
	buildcon.setSpecialSRAddConditions(fmparams.fm_maxtotalerrors,
					   fmparams.fm_maxgaps,
					   fmparams.fm_maxmismatches,
					   coel
	  );
	qaf.setWantsCleanOverlapEnds(coel);
	qaf.setMinTotalNonMatches(1);
	qaf.map();
	buildcon.updateBackboneConsensus();
      }
    }

    qaf.setWantsCleanOverlapEnds(0);
    qaf.setMinTotalNonMatches(0);

    CEBUG("bfccp1" << endl);

    if(fmparams.fm_mapperfect){
      cout << "Gogo: 100% mapping\n";
      buildcon.setSpecialSRAddConditions(0,0,0,0);
      qaf.map();
      buildcon.coutAddReadTimings();
    }

    if(hasotherst){
      cout << "Gogo: add others clean ends\n";
      qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_END);
      qaf.setWantsCleanOverlapEnds(16);
      qaf.setMinTotalNonMatches(0);
      // SRAddcondition is still 0,0,0 and cleanOverlap of 16 was also
      //  already done for Solexas, so no new Solexa should get added here
      qaf.map();
      buildcon.coutAddReadTimings();
    }

    qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_SOLEXA);
    qaf.setWantsCleanOverlapEnds(fmparams.fm_clean_end_dist);
    if(!fmparams.fm_mapperfect){
      qaf.setMinTotalNonMatches(1);
    }else{
      qaf.setMinTotalNonMatches(0);
    }
    if(fmparams.fm_maxtotalerrors>0 && qaf.getReadAddAttempts()>0) {
      if(fmparams.fm_maxmismatches != 0) {
	cout << "Gogo: mapping 1 mismatch\n";
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	buildcon.setSpecialSRAddConditions(1,0,1,
					   fmparams.fm_clean_end_dist
	  );
	qaf.map();
	buildcon.coutAddReadTimings();
	buildcon.updateBackboneConsensus();
      }

      if(fmparams.fm_maxgaps != 0) {
	cout << "Gogo: mapping 1 gap\n";
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	buildcon.setSpecialSRAddConditions(1,1,0,
					   fmparams.fm_clean_end_dist
	  );
	qaf.map();
	buildcon.coutAddReadTimings();
	buildcon.updateBackboneConsensus();
      }
    }

    if(fmparams.fm_maxtotalerrors>1 && qaf.getReadAddAttempts()>0) {
      if(fmparams.fm_maxmismatches >= 2){
	cout << "Gogo: mapping 2 mismatches\n";
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	buildcon.setSpecialSRAddConditions(2,0,2,
					   fmparams.fm_clean_end_dist
	  );
	qaf.map();
	buildcon.coutAddReadTimings();
	buildcon.updateBackboneConsensus();
      }
      if(fmparams.fm_maxgaps != 0 && fmparams.fm_maxmismatches != 0 ) {
	cout << "Gogo: mapping 1 gap, 1 mismatch\n";
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	buildcon.setSpecialSRAddConditions(2,1,1,
					   fmparams.fm_clean_end_dist
	  );
	qaf.map();
	buildcon.coutAddReadTimings();
	buildcon.updateBackboneConsensus();
      }

      if(fmparams.fm_maxgaps >= 2){
	cout << "Gogo: mapping 2 gaps\n";
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	buildcon.setSpecialSRAddConditions(2,2,0,
					   fmparams.fm_clean_end_dist
	  );
	qaf.map();
	buildcon.coutAddReadTimings();
	buildcon.updateBackboneConsensus();
      }
    }

    for(uint32 numerr=3; ; ){
      if(qaf.getReadAddAttempts()==0) break;
      cout << "Gogo: mapping all " << numerr << " errors\n";
      if(fmparams.fm_maxgaps>=0) {
	cout << "(max " << fmparams.fm_maxgaps << " gaps)\n";
      }
      if(fmparams.fm_maxmismatches>=0) {
	cout << "(max " << fmparams.fm_maxmismatches << " mismatches)\n";
      }
      if (fmparams.fm_maxgaps>=0 && fmparams.fm_maxgaps < numerr
	  && fmparams.fm_maxmismatches>=0 && fmparams.fm_maxmismatches < numerr){
	cout << "Max gaps and max mismatches are both <" << numerr << ", we're done.\n";
	break;
      }
      if(as_fixparams.as_dateoutput) dateStamp(cout);
      buildcon.setSpecialSRAddConditions(numerr,
					 fmparams.fm_maxgaps,
					 fmparams.fm_maxmismatches,
					 fmparams.fm_clean_end_dist
	);
      qaf.map();
      buildcon.coutAddReadTimings();
      buildcon.updateBackboneConsensus();
      if(numerr==fmparams.fm_maxtotalerrors) break;
      if(numerr<16){
	++numerr;
      }else if(numerr<32){
	numerr+=2;
      }else{
	numerr+=4;
      }
      if(numerr>fmparams.fm_maxtotalerrors) numerr=fmparams.fm_maxtotalerrors;
    }
  }else{
    // -AL:shme=0 means "quick mapping"
    // we did not go through the long and painful process above but now
    // have to force to map "whatever left" (which is all)
    hasotherst=true;
  }

  if(hasotherst){
    // TODO: no, not really the Solexa params. But probably never really used anyway, so ...
    auto & fmparams=AS_miraparams[ReadGroupLib::SEQTYPE_SOLEXA].getFinalMappingParameters();

    cout << "Gogo: mapping whatever left\n";
    qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_END);
    buildcon.setSpecialSRAddConditions(fmparams.fm_maxtotalerrors,
				       fmparams.fm_maxgaps,
				       fmparams.fm_maxmismatches,
				       0
      );
    qaf.setWantsCleanOverlapEnds(0);  // no clean end wanted
    qaf.setMinTotalNonMatches(0);  // no min total matches

    qaf.map();
    buildcon.updateBackboneConsensus();
  }

}
#define CEBUG(bla)



/*************************************************************************
 *
 * Move clusters smaller than wished minimum number of reads per contig
 *  to debris
 * However, not rails and not backbones!
 *
 * returns number of reads pushed to debris
 *
 *************************************************************************/

 //#define CEBUG(bla)   {cout << bla; cout.flush(); }
uint32 Assembly::bfc_moveSmallClustersToDebris()
{

  cout << "Moving small clusters to debris:\n";

  uint32 totaldebris=0;

  std::vector<uint32> numreadsperst(ReadGroupLib::getNumSequencingTypes(),0);
  std::vector<int32> clusteridperread;
  std::vector<std::list<int32> > readinclusterlist;

  clusterUnassembledReads(clusteridperread,readinclusterlist, AS_used_ids);

  uint32 clustercount=0;

  for(size_t ricli=0; ricli<readinclusterlist.size(); ++ricli){
    if(!readinclusterlist[ricli].empty()){
      uint32 totalreadsincluster=0;
      mstd::fill(numreadsperst,0);
      for(const auto & ricle : readinclusterlist[ricli]){
	++numreadsperst[AS_readpool[ricle].getSequencingType()];
	++totalreadsincluster;
      }
      bool takecluster=false;
      for(size_t st=0; st<numreadsperst.size(); ++st){
	if(numreadsperst[st]>0
	   && totalreadsincluster>=AS_miraparams[st].getAssemblyParams().as_minimum_readspercontig){
	  takecluster=true;
	}
      }
      // the above also kill mapping reads
      // therefore, first check whether an unsued rail is part of that cluster
      //  if yes, then don't kill cluster
      if(!takecluster){
	for(const auto & ricle : readinclusterlist[ricli]){
	  if(AS_readpool[ricle].isRail()
	     || AS_readpool[ricle].isBackbone()){

	    // should have this ... need to rework how as_usedids is filled (only after contig is made)
	    //&& !AS_used_ids[ricle]){

	    takecluster=false;
	    break;
	  }
	}
      }
      if(!takecluster){
	CEBUG("Killing cluster: " << ricli);
	for(const auto & ricle : readinclusterlist[ricli]){
	  CEBUG(" " << ricle << AS_readpool[ricle].getName());
	  AS_isdebris[ricle]=DEBRIS_TINYCLUSTER;
	  AS_used_ids[ricle]=1;
	  ++totaldebris;
	}
	CEBUG('\n');
      }
    }
  }

  // clusterUnassembledReads() did not return orphans as cluster list
  // therefore, look for unused read ids with no cluster number and also
  //  put them into debris

  for(size_t uid=0; uid<AS_used_ids.size(); ++uid){
    if(AS_used_ids[uid]==0 && clusteridperread[uid]==-1
       && !AS_readpool[uid].isRail()
       && !AS_readpool[uid].isBackbone()){
      CEBUG("Killing orphan: " << uid << endl);
      AS_isdebris[uid]=DEBRIS_TINYCLUSTERORPHAN;
      AS_used_ids[uid]=1;
      ++totaldebris;
    }
  }

  cout << "\nDone. " << totaldebris << " reads moved to debris.\n";

  return totaldebris;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * Checks if Contig meets specified requirements (atm: num of reads)
 *
 * If not, mark the reads in the contig as debris and empty the contig
 *
 *************************************************************************/

bool Assembly::bfc_checkIfContigMeetsRequirements(Contig & con)
{
  FUNCSTART("bool Assembly::bfc_checkIfContigMeetsRequirements(Contig & con)");

  bool contigok=false;

  std::vector<uint32> numreadsperst(ReadGroupLib::getNumSequencingTypes(),0);
  uint32 totalreadsincon=0;

  auto crI=con.getContigReads().begin();
  for(; crI!=con.getContigReads().end(); ++crI){
    if(crI.getORPID() >= 0){
      if(crI->isBackbone()){
	contigok=true;
	break;
      }
      ++numreadsperst[crI->getSequencingType()];
      ++totalreadsincon;
    }
  }

  if(!contigok){
    for(size_t st=0; st<numreadsperst.size(); ++st){
      if(numreadsperst[st]>0
	 && totalreadsincon>=AS_miraparams[st].getAssemblyParams().as_minimum_readspercontig){
	contigok=true;
      }
    }
  }

  return contigok;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::bfc_markRepReads(Contig & con)
{
  multitag_t tagtstR;
  tagtstR.setIdentifierStr("tstR");
  tagtstR.source=multitag_t::MT_tagsrcentry_idMIRA;

  auto & conreads=con.getContigReads();
  auto crI = conreads.begin();

  for(;crI != conreads.end(); crI++){
    if(crI.getORPID() >= 0
       && AS_multicopies[crI.getORPID()]) {
      cout << "xxxxxxxxxxxx mark " << crI.getORPID() << endl;
      Read & nonconstread = const_cast<Read &>(*crI);
      int32 rc=nonconstread.getRightClipoff()-1;
      if(rc<nonconstread.getLeftClipoff()) rc=nonconstread.getLeftClipoff();
      tagtstR.from=nonconstread.getLeftClipoff();
      tagtstR.to=rc;
      nonconstread.addTagO(tagtstR);
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::priv_tmpcheckroutine(Contig & buildcon)
{
  uint32 index=1;
  cout << "ptcr AS_used_ids[" << index << "]=" << static_cast<uint16>(AS_used_ids[index]) << endl;
  auto & cr=buildcon.getContigReads();
  bool foundidx=false;
  for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
    if(pcrI.getORPID()==index) {
      foundidx=true;
      break;
    }
  }
  cout << "Found idx: " << foundidx << endl;
}

/*************************************************************************
 *
 * return true if trimmed due to misassembly
 *
 *************************************************************************/

bool Assembly::bfc_trimDenovoIfNecessary(Contig & buildcon, bool foundSRMs, const std::string & basename_forextsave, uint32 & trackingunused)
{
  FUNCSTART("bool Assembly::bfc_trimDenovoIfNecessary(Contig & buildcon, bool foundSRMs, uint32 & trackingunused)");

  bool trimmedmisassembly=false;

  // if we have not been mapping:
  //  1) get misassembled parts out by looking at SRMc tags and trimming back to best range
  //  2) for genome assemblies with pairs: pair analysis and break contig at misassembled sites
  //  3) for genome assemblies: coverage analysis and remove reads in overcovered areas
  if(buildcon.getNumBackbones()==0){

    //priv_tmpcheckroutine(buildcon);
    if(foundSRMs){
      auto range=buildcon.findBestNonMisassembledRange();
      if(range.first>=0){
	cout<<"Found misassembly by repeat marker. Best range: " << range.first << ".." << range.second << '\t' << range.second-range.first << endl;
	trimmedmisassembly=true;
	auto & cr=buildcon.getContigReads();
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"srm ! AS_used_ids[...] ???");
	    AS_used_ids[pcrI.getORPID()]=0;
	  }
	}

	saveExtTmpContig(buildcon,(basename_forextsave+"_pretrimsrm"));

	cout << "Old trackingunused: " << trackingunused<< endl;
	trackingunused+=cr.size();
	cout << "Intermediate trackingunused: " << trackingunused<< endl;

	buildcon.trimContigToRange(range.first,range.second);

	// rewrite the AS_used_ids with the shortened
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    AS_used_ids[pcrI.getORPID()]=1;
	    --trackingunused;
	  }
	}
	cout << "New trackingunused: " << trackingunused<< endl;

	saveExtTmpContig(buildcon,(basename_forextsave+"_posttrimsrm"));
      }
    }

    //priv_tmpcheckroutine(buildcon);

    // pair consistency analysis and contig breaking for genome assemblies, denovo
    if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
       && AS_miraparams[0].getSpecialParams().mi_extra_flag3
       && !AS_hasbackbones){
      auto range=buildcon.findBestPairConsistencyRange();
      if(range.first>0 || range.second < buildcon.getContigLength()){
	cout<<"Found misassembly by pair consistency. Best range: " << range.first << ".." << range.second << '\t' << range.second-range.first << endl;
	trimmedmisassembly=true;
	auto & cr=buildcon.getContigReads();
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"pair consistency ! AS_used_ids[...] ???");
	    AS_used_ids[pcrI.getORPID()]=0;
	  }
	}

	saveExtTmpContig(buildcon,(basename_forextsave+"_pretrimpair"));

	cout << "Old trackingunused: " << trackingunused<< endl;
	trackingunused+=cr.size();
	cout << "Intermediate trackingunused: " << trackingunused<< endl;

	buildcon.trimContigToRange(range.first,range.second);

	// rewrite the AS_used_ids with the shortened
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    AS_used_ids[pcrI.getORPID()]=1;
	    --trackingunused;
	  }
	}
	cout << "New trackingunused: " << trackingunused<< endl;

	saveExtTmpContig(buildcon,(basename_forextsave+"_posttrimpair"));
      }
    }

    // priv_tmpcheckroutine(buildcon);

    // coverage analysis and coverage reduction for genome assemblies, denovo
    // note:
    // has a bug for data which was digitally normalised ...
    //  ... and besides, does absolutely not make sense to do this on
    //  digitally normalised data. Therefore, not done if diginorm active.
    if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
       && AS_miraparams[0].getSpecialParams().mi_extra_flag2
       && !AS_miraparams[0].getHashStatisticsParams().hs_apply_digitalnormalisation
       && !AS_hasbackbones){
      saveExtTmpContig(buildcon,(basename_forextsave+"_prered"));

      Contig::ccctype_t avgcovused=AS_coveragetotal;  // may be 0
      coverageinfo_t cinfo;
      std::vector<uint64> covvals;
      buildcon.collectCoverage(covvals);
      buildcon.calcStatsOnContainer(cinfo,covvals);
      cout << "1st covnum: " << cinfo << endl;

      // TODO: perhaps make this dependend of ratio mean vs stddev ?
      buildcon.calcSecondOrderStatsOnContainer(cinfo,covvals);
      cout << "2nd covnum: " << cinfo << endl;
      if(cinfo.median>2*avgcovused) avgcovused=cinfo.median;
      cout << "Using: " << avgcovused << endl;

      std::vector<uint8> peakindicator;
      buildcon.findPeaks(avgcovused,peakindicator);
      std::unordered_set<readid_t> readsremoved;
      buildcon.reduceReadsAtCoveragePeaks(avgcovused,peakindicator,readsremoved);
      cout << "Coverageremove: " << readsremoved.size() << endl;

      // if reads were removed, get the tracking corrected
      if(!readsremoved.empty()){
	for(auto & rid : readsremoved){
	  BUGIFTHROW(!AS_used_ids[rid],"rere ! AS_used_ids[rid] ???");
	  AS_used_ids[rid]=0;
	}
	cout << "Old trackingunused: " << trackingunused<< endl;
	trackingunused+=readsremoved.size();
	cout << "Intermediate trackingunused: " << trackingunused<< endl;

	saveExtTmpContig(buildcon,(basename_forextsave+"_postred"));
      }
    }
    //priv_tmpcheckroutine(buildcon);
  }

  return trimmedmisassembly;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Assembly::bfc_savePersistentSmallOverlaps(Contig & thiscon, const int32 passnr, std::ostream & fout)
{
  FUNCSTART("void Assembly::bfc_savePersistentSmallOverlaps(Contig & thiscon, const int32 passnr, std::ostream & fout)");

  CEBUG("bfc_savePersistentSmallOverlaps\n");

  // we'll save all overlaps between min and max
  auto minbph=priv_calcBasesPerHashOnPass(passnr);
  auto maxbph=priv_calcBasesPerHashOnPass(passnr+1);

  CEBUG("min " << minbph << " max " << maxbph << endl);

  if(minbph==maxbph) return;

  CEBUG("Needs persistent SO consideration. min " << minbph << " max " << maxbph << "\n");
  newedges_t tmpedge;

  auto pcrE=thiscon.getContigReads().end();
  for(auto pcrI=thiscon.getContigReads().begin(); pcrI!=pcrE; ++pcrI){
    auto orpid=pcrI.getOriginalReadPoolID();
    CEBUG("VVV orpid " << orpid << endl);
    tmpedge.rid1=orpid;
    auto runI=lower_bound(AS_confirmed_edges.begin(),
			  AS_confirmed_edges.end(),
			  tmpedge,
			  newedges_t::sortComparatorByRIDUp);

    // rare, but may happen: singlets with no overlaps
    // user may have asked to include them, but they did not have any overlaps
    // so, need to account for this.
    if(runI==AS_confirmed_edges.end()
       || runI->rid1 > orpid) continue;

    if(runI->rid1!=orpid){
      cout << "\n\nGaaaaaahhhhhhhh!\n"
	   << "runI->rid1 " << runI->rid1 << " != orpid " << orpid << endl;
      cout << "confedges.size() " << AS_confirmed_edges.size() << endl;
      cout << "runI dist: " << runI-AS_confirmed_edges.begin() << endl;
      cout << "runI -1 ->rid1: " << (runI-1)->rid1 << endl;
      BUGIFTHROW(runI->rid1!=orpid,"runI->rid1!=orpid ... really???");
    }

    for(;runI!=AS_confirmed_edges.end() && runI->rid1==orpid; ++runI){
      // save each pair of read-links only once, leave out the other half
      // (it's duplicated again in Assembly::loadAlignmentsFromFile())
      CEBUG("VVV " << AS_adsfacts[runI->adsfindex] << "VVV minmax " << minbph << " " << maxbph << "\n");
      if(runI->rid1 <= runI->linked_with){
	CEBUG("VVV consider\n");
	auto & actadsf = AS_adsfacts[runI->adsfindex];
	CEBUG(*runI << "\tovl: " << actadsf.getOverlapLen() << endl);
	if(actadsf.getOverlapLen() >= minbph
	   && actadsf.getOverlapLen() < maxbph
	  ){
	  CEBUG("VVV save\n");
	  fout << runI->best_weight
	       << '\t' << runI->direction
	       << '\t' << runI->ol_stronggood
	       << '\t' << runI->ol_weakgood
	       << '\t' << runI->ol_belowavgfreq
	       << '\t' << runI->ol_norept
	       << '\t' << runI->ol_rept
	       << '\t';
	  actadsf.serialiseOut(fout);
	  fout << '\n';
	}
      }
    }
  }

  fout.flush();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define VCOUT(bla)   {if(verbose) {cout << bla;} }

void Assembly::bfc_storeContig(Contig & con, uint32 & numcontigs, const bool mustmarkrepeats, const int32 passnr, const bool lastpass, const bool verbose)
{
  FUNCSTART("void Assembly::bfc_storeContig(Contig & con, uint32 & numcontigs, const bool mustmarkrepeats, const int32 passnr, const bool lastpass)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  auto & conreads=con.getContigReads();

  // look whether we store singlets in contig or not

  bool storecontig=true;

  if(con.getNumReadsInContig()==1){
    const assembly_parameters & as_rt_params = AS_miraparams[conreads.begin()->getSequencingType()].getAssemblyParams();

    storecontig=as_rt_params.as_savesimplesingletsinproject;
    if(conreads.begin()->hasTag(Read::REA_tagentry_idSRMr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idCRMr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idWRMr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idSROr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idSAOr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idSIOr)) {
      if(conreads.begin().getORPID() >= 0) AS_needalloverlaps[conreads.begin().getORPID()]=true;
      storecontig|=as_fixparams.as_savetaggedsingletsinproject;
     }
  }

  if(as_fixparams.as_backbone_trimoverhangingreads){
    con.trimMapOverhang();
  }

  // store contig names only for de-novo genome assemblies
  // used for listing large contigs at end of assembly
  std::string cnameforstore;
  if(!as_fixparams.as_assemblyjob_mapping
     && AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
    cnameforstore=con.getContigName();
  }

  // TODO: U13 hat 0 werte im ASCII?
  //Contig::setCoutType(Contig::AS_DEBUG);
  //cout << "Debug in bfc_storeContig()\n" << con;

  if(storecontig){
    if(as_fixparams.as_mark_repeats&&mustmarkrepeats){
      std::vector<bool> dummy;
      Contig::repeatmarker_stats_t contigrepstats;
      markRepeats(con, dummy,contigrepstats);
      AS_bfcstats[con.getLongRepeatStatus()].numnewsrm+=contigrepstats.numSRMs;
    }

    //bfc_markRepReads(con);

    VCOUT("Storing contig ... " << as_fixparams.as_mark_repeats << mustmarkrepeats);
    if(AS_hasbackbones){
      con.removeRails();

      // TODO: ask contig whether it has mappings
      //if(as_fixparams.as_loadSOLEXA || as_fixparams.as_loadSOLID){
      if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]
	 || AS_seqtypespresent[ReadGroupLib::SEQTYPE_ABISOLID]) {
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	cout << "Transforming CER mappings." << endl;
	con.transformCERMappingsToCoverageReads();
	cout << "done transforming CER mappings." << endl;
	if(as_fixparams.as_dateoutput) dateStamp(cout);
	//assout::saveAsMAF(con, getMAFFilename()+".bla", AS_deleteoldresultfiles);
      }
    }

    con.markFeaturesByConsensus(true, true, true);
    // transfer important tags to readpool
    VCOUT("Transfering tags to readpool." << endl);
    transferContigReadTagsToReadpool(con,AS_bbcontigs.end());

    con.updateStatsFromConsensusTags(true,true,true,true,true);

    // remove HAF tags if not wished
    if(!as_fixparams.as_buntify_reads){
      cout << "Removing HAF tags (set -GE:crkf if you want to change that)" << endl;
      for(auto & pcre : con.getContigReads()){
	for(uint32 i=0; i<Read::REA_allhaftags.size(); ++i){
	  const_cast<Read &>(pcre).deleteTag(Read::REA_allhaftags[i]);    // we know what we do ... *cough*
	}
      }
    }

    // store the contig information
    AS_assemblyinfo.storeContigStats(con.getStats(),cnameforstore);


    if(lastpass) {
      assout::saveStatistics(con,
			     getStatisticsFilename(),
			     AS_deleteoldresultfiles);
      assout::saveReadTagList(con,
			      getReadTagListFilename(),
			      AS_deleteoldresultfiles);
      assout::saveConsensusTagList(con,
				   getConsensusTagListFilename(),
				   AS_deleteoldresultfiles);
      assout::saveContigReadList(con,
				 getContigReadListFilename(),
				 AS_deleteoldresultfiles);
      if(as_fixparams.as_output_gff3){
	// store sequence for later
//	  AS_gff3defer_names.push_back(con.getContigName());
//TODO: weiterhier
//	  assout::saveTagsAsGFF3(con, getGFF3Filename(), AS_deleteoldresultfiles);
      }
      if(as_fixparams.as_output_caf){
	VCOUT("Saving CAF ... "; cout.flush());
	assout::saveAsCAF(con, getCAFFilename(), AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_maf){
	VCOUT("Saving MAF ... "; cout.flush());
	assout::saveAsMAF(con, getMAFFilename(), AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_gap4da){
	VCOUT("Saving gap4 direct assembly ... "; cout.flush());
	assout::saveAsGAP4DA(con,getGAP4DAFilename(),AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_fasta) {
	if(ReadGroupLib::getNumOfStrains()>1){
	  VCOUT("Saving strains as FASTA ... "; cout.flush());
	  assout::saveStrainsAsFASTAQ(con, AS_readpool,
				      buildDefaultResultsFileName(
					-1,"","", "",
					AS_miraparams[0].getAssemblyParams().as_outfile_FASTA,
					""),
				      false,
				      0,0,
				      AS_deleteoldresultfiles,false);
	}else{
	  VCOUT("Saving FASTA ... "; cout.flush());
	  assout::saveAsFASTA(con,
			      getFASTAFilename(),
			      getFASTAPaddedFilename(),
			      AS_deleteoldresultfiles);
	}
	VCOUT("done.\n");

      }
      if(as_fixparams.as_output_tcs) {
	VCOUT("Saving TCS ... "; cout.flush());
	assout::saveAsTCS(con, getTCSFilename(),AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_wiggle) {
	VCOUT("Saving Wiggle ... "; cout.flush());
	assout::saveAsWiggle(con, getWiggleFilename(),AS_deleteoldresultfiles, false);
	VCOUT("done.\n");
      }
      // TODO: enable these functions for incremental write
      //saveSNPAnalysis();
      //saveFeatureAnalysis();
      if(as_fixparams.as_output_txt){
	VCOUT("Saving text ... "; cout.flush());
	assout::saveAsTXT(con,getTXTFilename(),AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_ace){
	VCOUT("Saving ACE ... "; cout.flush());
	assout::saveAsACE(con,getACEFilename(),AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_html) {
	VCOUT("Saving HTML ... "; cout.flush());
	assout::dumpContigAsHTML(con,
				 getHTMLFilename(),
				 AS_deleteoldresultfiles,
				 AS_miraparams[0].getAssemblyParams().as_projectname_out);
	VCOUT("done.\n");
      }
    }else{
      assout::saveStatistics(con,
			     getStatisticsFilename(passnr, "", "_pass"),
			     AS_deleteoldresultfiles);
      assout::saveReadTagList(con,
			      getReadTagListFilename(passnr),
			      AS_deleteoldresultfiles);
      assout::saveConsensusTagList(con,getConsensusTagListFilename(passnr),
				   AS_deleteoldresultfiles);
      assout::saveContigReadList(con,
				 getContigReadListFilename(passnr, "", "_pass"),
				 AS_deleteoldresultfiles);

      if(as_fixparams.as_output_tmp_caf) {
	VCOUT("Saving temp CAF ... "; cout.flush());
	assout::saveAsCAF(con,
			  getCAFFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_tmp_maf) {
	VCOUT("Saving temp MAF ... "; cout.flush());
	assout::saveAsMAF(con,
			  getMAFFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_tmp_gap4da) {
	VCOUT("Saving temp gap4 direct assembly ... "; cout.flush());
	assout::saveAsGAP4DA(con,
			     getGAP4DAFilename(passnr, "", "_pass"),
			     AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_tmp_fasta){
	VCOUT("Saving temp FASTA ... "; cout.flush());
	assout::saveAsFASTA(con,
			    getFASTAFilename(passnr, "", "_pass"),
			    getFASTAPaddedFilename(passnr, "", "_pass"),
			    AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_tmp_txt){
	VCOUT("Saving temp text ... "; cout.flush());
	assout::saveAsTXT(con,
			  getTXTFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_tmp_ace) {
	VCOUT("Saving temp ACE ... "; cout.flush());
	assout::saveAsACE(con,
			  getACEFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      if(as_fixparams.as_output_tmp_tcs) {
	VCOUT("Saving temp TCS ... "; cout.flush());
	assout::saveAsTCS(con,
			  getTCSFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	VCOUT("done.\n");
      }
      //if(as_fixparams.as_output_tmp_html) saveAsHTML(passnr, "", "_pass");
      if(as_fixparams.as_output_tmp_html) {
	VCOUT("Saving temp HTML ... "; cout.flush());
	assout::dumpContigAsHTML(con,
				 getHTMLFilename(passnr, "", "_pass"),
				 AS_deleteoldresultfiles,
				 AS_miraparams[0].getAssemblyParams().as_projectname_out);
	VCOUT("done.\n");
      }
    }
    VCOUT("done." << endl);
  }else{
    // store the contig information
    AS_assemblyinfo.storeContigStats(con.getStats(),cnameforstore);

    if(conreads.begin().getORPID() >= 0) AS_isdebris[conreads.begin().getORPID()]=DEBRIS_UNSAVEDSINGLET;
    numcontigs--;
  }

  AS_deleteoldresultfiles=false;

  return;
}
#define VCOUT(bla)
