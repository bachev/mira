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


#include "assembly.H"

#include "mira/newpathfinder.H"


using namespace std;



#define CEBUG(bla)
#define CEBUGF(bla)



#if 1


/*************************************************************************
 *
 *
 *
 *************************************************************************/
void Assembly::prework()
{
  FUNCSTART("Assembly::prework()");

  ensureStandardDirectories(false);

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  if (as_fixparams.as_numrmbbreakloops <1) {
    const_cast<assembly_parameters &>(as_fixparams).as_numrmbbreakloops=1;
    cout << "Number of RMB break loops <1, setting to 1\n";
  }

  // find out whether we have 454 reads in the data set
  // also find out whether there are SRMr or CRMr tags that need attention
  AS_seqtypespresent.clear();
  AS_seqtypespresent.resize(ReadGroupLib::SEQTYPE_END,false);
  for(uint32 i=0; i< AS_readpool.size(); i++){
    if(AS_readpool.getRead(i).isUsedInAssembly()){
      AS_seqtypespresent[AS_readpool.getRead(i).getSequencingType()]=true;
    }
  }

  // allocate or reserve memory that is quite static from the size,
  //  reducing memory fragmentations at least a bit
  {
    AS_permanent_overlap_bans.nuke();
    AS_permanent_overlap_bans.resize(AS_readpool.size());

    AS_used_ids.reserve(AS_readpool.size());
    AS_clipleft.reserve(AS_readpool.size());
    AS_clipright.reserve(AS_readpool.size());
    AS_multicopies.reserve(AS_readpool.size());
    AS_hasmcoverlaps.reserve(AS_readpool.size());
    AS_istroublemaker.reserve(AS_readpool.size());

    // these 3 only temp filled, but block anyway to reduce fragmentation
    AS_allrmbsok.reserve(AS_readpool.size());
    AS_probablermbsnotok.reserve(AS_readpool.size());
    AS_weakrmbsnotok.reserve(AS_readpool.size());
  }


  findPossibleOverlaps(0, "", "_preassembly");

  if(as_fixparams.as_clip_possible_vectors
    || (as_fixparams.as_use_read_extension && as_fixparams.as_readextension_firstpassnum == 0)){
    cout << "Pre-assembly alignment search for read extension and / or vector clipping:\n";
    makeAlignments(Assembly::ma_takeall,
		   false,
		   false,
		   0,
		   "",
		   "_preassembly");
    loadAlignmentsFromFile(0, "", "_preassembly");
    if(as_fixparams.as_use_read_extension) {
      cout << "Pre-assembly read extension:\n";
      extendADS(0, "", "_preassembly1");
    }
    if(as_fixparams.as_clip_possible_vectors) {
      cout << "Pre-assembly vector clipping\n";
      performSeqVectorClippings();
    }

    //nukeSTLContainer(AS_adsfacts);
    //nukeSTLContainer(AS_confirmed_edges);
    AS_adsfacts.clear();
    AS_confirmed_edges.clear();

    postLoad();
  }


  bool exitearly=false;
  for(uint32 actpass=1; !exitearly && actpass<=as_fixparams.as_numpasses; actpass++){
    cout << "\n\nPass: " << actpass << endl;

    AS_contigs.clear();

    findPossibleOverlaps(actpass, "", "_pass");

    makeAlignments(Assembly::ma_takeMCandOverlapWMC,
		   false,
		   false,
		   actpass,
		   "",
		   "_pass");

    bool newreptfound=true;
    // maximum of three subloops, might be useful for genome sized projects
    for(uint32 subloop=0; subloop<3 && newreptfound; subloop++){
      loadAlignmentsFromFile(actpass, "", "_pass");

      newreptfound=buildRepeatContigs(actpass);

      cout << "Done building contigs.\n";

      if(newreptfound){
	cout << "New repeats found during contig building, adding additional alignment iteration\nfor quick repeat resolving.\n";
	AS_steps[ASADSLISTOK]=0;

	size_t oldnumbans=0;
	size_t oldnumsets=0;
	AS_permanent_overlap_bans.getNumBans(oldnumbans,oldnumsets);

	cout << "Old number of bans: " << oldnumbans << '\n';
	makeAlignments(Assembly::ma_needSRMrOrTwoCRMr,
		       true,
		       false,
		       actpass,
		       "",
		       "",
		       "repeat_resolve");

	size_t newnumbans=0;
	size_t newnumsets=0;
	AS_permanent_overlap_bans.getNumBans(newnumbans,newnumsets);

	cout << "New number of bans: " << newnumbans << '\n';
	if(oldnumbans==newnumbans){
	  cout << "No new bans appeared, next pass.\n";
	  newreptfound=false;
	}
      }
    }

    if(!newreptfound){
      if(actpass>=3 && actpass <as_fixparams.as_numpasses) exitearly=true;

      //// if we have 454 data that can be edited, we want at least
      ////  3 passes (for the simple 454 editing)
      //if(AS_454dosimpleedit){
      //	if(actpass>=3 && actpass <as_fixparams.as_numpasses) exitearly=true;
      //}else{
      //	if(actpass<as_fixparams.as_numpasses) exitearly=true;
      //}
    }

    //saveStatistics(actpass, "", "_pass");
    //saveDebrisList(actpass, "", "_pass");
    //saveReadTagList(actpass);
    //saveContigReadList(actpass, "", "_pass");
    if(as_fixparams.as_output_tmp_caf) saveAsCAF(actpass, "", "_pass");
  }

  if(exitearly){
    cout << "No new repeats found in previous pass, stopping early.\n";
  }

  FUNCEND();
  return;
}


/*************************************************************************
 *
 * either a multicopy, troublemaker or has SRMr or CRMr
 *
 *************************************************************************/

bool Assembly::ma_takeMCandOverlapWMC(Assembly & as, int32 rid1, int32 rid2)
{
  if(as.AS_multicopies[rid1]
     || as.AS_istroublemaker[rid1]
     || as.AS_multicopies[rid2]
     || as.AS_istroublemaker[rid2]) return true;
  if(Assembly::ma_needSRMrOrTwoCRMr(as, rid1, rid2)) return true;
  return false;
}



/*************************************************************************
 *
 * returns whether new strong repeat markers (SRMs) were found for
 *  any contig built in any stage
 *
 *************************************************************************/

bool Assembly::buildRepeatContigs(const int32 passnr)
{
  FUNCSTART("void Assembly::buildRepeatContigs()");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  contig_parameters const & con_params= AS_miraparams[0].getContigParams();
  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();

  AS_contigs.clear();
  Contig::setIDCounter(1);

  vector<Align> aligncache;
  setupAlignCache(aligncache);


  // filter out all reads that have no overlap at this stage
  // for this, first set all reads to "used" and "singlet", then
  //   re-allow those that have at least one SW ovelap of a read
  //   which is MULTICOPY!

  AS_used_ids.clear();
  AS_used_ids.resize(AS_readpool.size(),1);
  AS_isdebris.clear();
  AS_isdebris.resize(AS_readpool.size(),1);

  if(!AS_confirmed_edges.empty()){
    // search for singlets
    vector<newedges_t>::const_iterator ceI=AS_confirmed_edges.begin();
    for(;ceI!=AS_confirmed_edges.end();){
      int32 actrid=ceI->rid1;
      uint32 numoverlaps=0;
      //cout << actrid << endl;

      bool hasmcoverlap= AS_multicopies[actrid]>0;
      for(;ceI!=AS_confirmed_edges.end() && ceI->rid1 == actrid;ceI++){
	numoverlaps++;
	hasmcoverlap|=AS_multicopies[ceI->linked_with]>0;
      }
      if(numoverlaps && hasmcoverlap){
	AS_isdebris[actrid]=0;
	AS_used_ids[actrid]=0;
      }
    }
  }

  // Now, go through all the singlets and take them back into
  //  assembly if they have special MIRA tags attached
  // Also remove backbones & rails from the debris list
  {
    for(uint32 i=0; i< AS_readpool.size(); i++){
      if(AS_readpool.getRead(i).hasTag(Read::REA_tagentry_idSRMr)
	 || AS_readpool.getRead(i).hasTag(Read::REA_tagentry_idCRMr)
	 || AS_readpool.getRead(i).hasTag(Read::REA_tagentry_idWRMr)
	 || AS_readpool.getRead(i).hasTag(Read::REA_tagentry_idSAOr)
	 || AS_readpool.getRead(i).hasTag(Read::REA_tagentry_idSROr)
	 || AS_readpool.getRead(i).hasTag(Read::REA_tagentry_idSIOr)) {
	AS_isdebris[i]=0;
	AS_used_ids[i]=0;
      }
      if(AS_readpool.getRead(i).isBackbone()
	 || AS_readpool.getRead(i).isRail()) {
	AS_isdebris[i]=0;
      }
    }
  }

  // get out reads that are
  //  - invalid or not used in assembly
  //  - or backbone or rail
  {
    for(uint32 i=0; i< AS_readpool.size(); i++){
      if(AS_readpool.getRead(i).hasValidData()== false
	|| AS_readpool.getRead(i).isUsedInAssembly()==false){
	AS_used_ids[i]=-1;
      } else if(AS_readpool.getRead(i).isBackbone()
		|| AS_readpool.getRead(i).isRail()) {
	AS_used_ids[i]=1;
      }
    }
  }

  uint32 unused;
  uint32 numcontigs=0;
  ofstream fout;
  if(as_fixparams.as_tmpf_unused_ids.size()!=0){
    fout.open((dir_params.dir_tmp+"/"+as_fixparams.as_tmpf_unused_ids).c_str(), ios::out);
    fout.close();
  }

  // outside precomputed lowerbound of oedges, for PathFinder
  // is used lateron in constructStepByStep()
  vector<vector<newedges_t>::iterator> tmp_lowerbound_oedges;

  // PathFinder object can be created ouside loop and re-used
  //  (constructor needs to do a lot of lower_bound() searches.
  // Saves a lot of time in the endgame (singlets) of big projects
  Pathfinder paf(&AS_miraparams,
		 AS_readpool,
		 AS_confirmed_edges,
		 AS_adsfacts);
  vector<int32> pf_ids_in_contig;

  bool foundSRMs=false;

  do{
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    cout << "\n";

    numcontigs++;
    if(as_fixparams.as_tmpf_unused_ids.size()!=0){
      fout.open((dir_params.dir_tmp+"/"+as_fixparams.as_tmpf_unused_ids).c_str(), ios::out|ios::app);
      fout << "\nUnused for contig" << numcontigs << endl;
    }
    unused=0;
    for(uint32 i=0; i<AS_used_ids.size(); i++){
      if(AS_used_ids[i]==0){
	unused++;
	if(as_fixparams.as_tmpf_unused_ids.size()!=0){
	  fout << AS_readpool.getRead(i).getName();
	  fout << "\t" << AS_readpool.getRead(i).getSCFName() << endl;;
	}
      }
    }
    if(as_fixparams.as_tmpf_unused_ids.c_str()!=0){
      fout.close();
    }

    if(unused>0){
      Contig::setIDCounter(numcontigs);
      Contig con(&AS_miraparams, AS_readpool);

      // set newreptmarked to true just to get into the for pass
      //  value is changed within pass
      bool newreptmarked=true;

      // in contrast to buildFirstContigs(), build a contig once

      con.setContigNamePrefix(con_params.con_nameprefix);

      cout << "Building new contig " << numcontigs << endl;

      if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

      cout << "Unused reads: " << unused << endl;
      cout.flush();

      paf.constructStepByStep(aligncache,
			      &AS_used_ids,
			      &pf_ids_in_contig,
			      &AS_multicopies,
			      &AS_hasmcoverlaps,
			      &AS_istroublemaker,
			      &tmp_lowerbound_oedges,
			      con);

      cout << "\n\nFinished building the contig." << endl;

      try {
	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
	if(con.getContigLength()>100000){
	  cout << "Calculating statistics (this may take a while)." << endl;
	}
	con.stats(cout);
	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      }
      catch (...) {
	cerr << "Darn, error with that contig. See darn.fasta.\n";
	Read::setCoutType(Read::AS_CLIPPEDFASTA);
	for(auto & cre : con.getContigReads()) {
	  cout << cre;
	}
	exit(0);
      }

      string basename_forextsave;
      {
	ostringstream ostr;
	ostr << dir_params.dir_tmp << "/miralog.pass_" << passnr << "_cb" << numcontigs << "_";

	basename_forextsave=ostr.str();
      }
      // saving pre-edit
      saveExtTmpContig(con,(basename_forextsave+"pre"));

      newreptmarked=false;

      {
	cout << "Marking possibly misassembled repeats ...";
	cout.flush();
	Contig::repeatmarker_stats_t repstats;
	vector<bool> readsmarkedsrm;
	con.newMarkPossibleRepeats(repstats, readsmarkedsrm);
	if(repstats.numSRMs>0 || repstats.numWRMs>0 || repstats.numSNPs>0){
	  if(repstats.numSRMs>0) {
	    newreptmarked=true;
	    foundSRMs=true;
	  }
	  cout << "done. Found\n";
	  cout << " - " << repstats.numSRMs << " new Strong RMB (SRMc)\n";
	  cout << " - " << repstats.numWRMs << " new Weak RMB (WRMc)\n";
	  cout << " - " << repstats.numSNPs << " SNP\npositions tagged.";
	}else{
	  cout << "done. Found none." << endl;
	}
      }

      // saving again if rept-marked or edited
      if(newreptmarked){
	saveExtTmpContig(con,(basename_forextsave+"post"));
      } else {
	if ( as_fixparams.as_output_exttmp_fasta
	     || as_fixparams.as_output_exttmp_ace
	     || as_fixparams.as_output_exttmp_gap4da
	     || as_fixparams.as_output_exttmp_caf) {
	  cout << "No edit and no new repeat found, not saving extra temporary contig again.\n";
	}
      }

      // Transfer all the reads fron the new contig into readpool
      //  (we can do that because no editing is done while assembling
      //  repeats)

      cout << "Transfering reads to readpool.\n";

      auto & cr = con.getContigReads();
      auto crI = cr.begin();

      for(;crI!=cr.end(); ++crI){
	if(crI.getORPID()>=0){
	  AS_readpool[crI.getORPID()]=*crI;
	  AS_readpool[crI.getORPID()].removeGapsFromRead();
	}
      }

      // don't store 454 singlets in contigs
      // except if they have special MIRA tags
      bool storecontig=true;
      if(con.getNumReadsInContig()==1){
	auto & conreads=con.getContigReads();
	if(conreads.begin()->getSequencingType()==ReadGroupLib::SEQTYPE_454GS20
	   && !conreads.begin()->hasTag(Read::REA_tagentry_idSRMr)
	   && !conreads.begin()->hasTag(Read::REA_tagentry_idCRMr)
	   && !conreads.begin()->hasTag(Read::REA_tagentry_idWRMr)
	   && !conreads.begin()->hasTag(Read::REA_tagentry_idSAOr)
	   && !conreads.begin()->hasTag(Read::REA_tagentry_idSROr)
	   && !conreads.begin()->hasTag(Read::REA_tagentry_idSIOr)) {
	  if(conreads.begin().getORPID()>=0) AS_isdebris[conreads.begin().getORPID()]=1;
	  storecontig=false;
	}
      }
      if(storecontig){
	AS_contigs.push_back(con);
	AS_contigs.back().saveMem();
      }else{
	numcontigs--;
      }
    }
    //unused=0;

  }while(unused>0);

  AS_steps[ASCONTIGSOK]=1;

  AS_used_ids.clear();

  //  saveAsCAF();

  FUNCEND();

  return foundSRMs;
}







#endif
