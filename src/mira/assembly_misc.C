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

// functions to process reads
// currently in namespace and object assembly


#include <regex>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/format.hpp>


#include "util/progressindic.H"

#include "mira/assembly.H"
#include "mira/align.H"
#include "mira/ads.H"
#include "mira/dataprocessing.H"
#include "mira/hashstats.H"

#include "util/stlimprove.H"


using std::cout;
using std::endl;


//#define CEBUG(bla)   {if(CEBUGFLAG) {cout << bla; cout.flush();}}
#define CEBUG(bla)




//#define CEBUG(bla)   {if(id1==2282 && id2==342) {cout << bla; cout.flush();}}
//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla;}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

template<class TVHASH_T>
void Assembly::priv_phahelper(const std::string & filenameforkms, const std::string & signalfile, uint32 basesperhash, bool rarekmerfinalkill, int32 version, const std::string prefix, const std::string postfix, const std::string logname)
{
  FUNCSTART("void Assembly::priv_phahelper(uint32 basesperhash, bool rarekmerfinalkill, int32 version, const std::string prefix, const std::string postfix, const std::string logname)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();
  hashstatistics_parameters const & hs_params= AS_miraparams[0].getHashStatisticsParams();

  // BaCh 22.04.2013: MNR tags are now always deleted by assignReadBaseStatistics(), so must be also re-set
  uint32 nastyrepeatratio=hs_params.hs_nastyrepeatratio;;
  uint32 nastyrepeatcoverage=hs_params.hs_nastyrepeatcoverage;
  bool masknastyrepeats=hs_params.hs_masknastyrepeats;
  // but no masking if diginorm had already been applied
  if(nastyrepeatratio>0
     && hs_params.hs_masknastyrepeats
     && hs_params.hs_apply_digitalnormalisation
     && version>AS_applydiginorminpass){
    nastyrepeatratio=0;
    nastyrepeatcoverage=0;
    masknastyrepeats=false;
    cout << "Masking of nasty repeats switched of as diginorm already ran.\n";
  }

  HashStatistics<TVHASH_T> s3;
  std::string loadhsfn;
  {
    bool dosdbgchimera=AS_miraparams[0].getAssemblyParams().as_clip_sdbg_chimeradetection && basesperhash>25;
    bool dosdbgedit=AS_miraparams[0].getEditParams().ed_sdbg_readedit;

    s3.setHashFrequencyRatios(hs_params.hs_freqest_minnormal,
			      hs_params.hs_freqest_maxnormal,
			      hs_params.hs_freqest_repeat,
			      hs_params.hs_freqest_heavyrepeat,
			      hs_params.hs_freqest_crazyrepeat,
			      2,                                  // rare kmer count: occurrence <= that can be masked away later by rare kmer masking
			      nastyrepeatratio,
			      nastyrepeatcoverage);

    s3.setAvgHashFreqMinimum(hs_params.hs_freq_covestmin);

    bool havestats=false;
    if(!signalfile.empty() && fileExists(signalfile)){
      //cout << "XXX signal file exists " << signalfile << ", loading " << filenameforkms << endl;
      s3.loadHashStatistics(filenameforkms);
      loadhsfn=filenameforkms;
      havestats=true;
    }

    if(!havestats){
      loadhsfn=filenameforkms;
      //cout << "XXX no stats, setting " << loadhsfn << endl;

      uint32 rkfkval=0;
      if(rarekmerfinalkill) rkfkval=hs_params.hs_rare_kmer_final_kill;

      // if we do a chimera search OR some SDBG edits, do NOT take rails in to the
      //  calculation of hash statistics! These rails could have come from a "wrong"
      //  backbone and therefore led to wrong chimera recognition or read edits!
      // if neither is done, we can take the rails as "maybe, possibly true" to help in
      //  areas difficult for sequencing (e.g. GGCxG in Illuminas)
      bool alsorailsinhs=!(dosdbgchimera | dosdbgedit);

      s3.computeHashStatistics(AS_readpool,
			       hs_params.hs_memtouse,
			       true,
			       alsorailsinhs,
			       true,
			       as_fixparams.as_clip_pec_mkfr,   // TODO: ok to misuse pec_mkfr or needs own parameter?
			       rkfkval, // TODO: setting it here also influences SDBG routines. TODO: check if not wanted: filter afterwards
			       basesperhash,
			       filenameforkms,
			       AS_miraparams[0].getDirectoryParams().dir_tmp
	);
      if(!signalfile.empty()){
	std::ofstream sout(signalfile);
      }
    }
    s3.showHashStatisticsInfo();

    if(dosdbgchimera){
      uint32 trimfreq=3;
      auto avghashfreq=s3.getAvgHashFreqRaw();
      if(avghashfreq>60) trimfreq=4;
      if(avghashfreq>80) trimfreq=5;
      // version is pass ... have at least trimfreq = 4 on third pass and after
      if(version>=3 && trimfreq<4) trimfreq=4;

      if(basesperhash<25 && trimfreq<5) trimfreq=1;

      std::string tmpfname;
      tmpfname=buildFileName(0,"","",
			    as_fixparams.as_tmpf_clippings,
			     ".txt","",false);
      AS_dataprocessing.startLogging(tmpfname,false);
      cout << "Hunting down chimeras:\n";
      auto hits=AS_dataprocessing.performSDBGChimeraSearch_Pool(AS_readpool,s3,trimfreq,&AS_debrisreason,"HA-"+boost::lexical_cast<std::string>(basesperhash));
      cout << "\nHA-" << basesperhash << " chimera trophy count: " << hits << endl;

      cout << "Reloading statistics 1 " << loadhsfn << " ... "; cout.flush();
      s3.discard();
      s3.loadHashStatistics(loadhsfn);
      cout << "done.\n";
    }
    if(dosdbgedit){
      cout << "Performing HA SDBG edits:\n";
      auto trimfreq=1;  // trimfreq==1 -> very conservative
      if(version>=2) trimfreq=2;

      // more liberal edits not for EST/RNASEq (low coverage variants!)
      if(!AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	if(version>=4) trimfreq=3;
	if(version>=5) trimfreq=4;
      }
      auto hits=AS_dataprocessing.performSDBGEdits_Pool(AS_readpool,s3,trimfreq);
      cout << "\nHA SDBG edited reads: " << hits << endl;
      cout << "Reloading statistics 2 " << loadhsfn << " ... "; cout.flush();
      s3.discard();
      s3.loadHashStatistics(loadhsfn);
      cout << "done.\n";
    }

    // for EST (especially RNASeq): way to be less harsh for chimerakilling
    if(dosdbgchimera
       && basesperhash>=17 && basesperhash<=64
       && !AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){

      if(AS_estnochimerakill.empty()){
	AS_estnochimerakill.resize(AS_readpool.size(),0);
      }

      auto oldhssize=s3.getNumHashEntries();
      s3.trimHashStatsByFrequencyANDOR(3,3,3); // trim for "at least three values" (fwd/rev don't care)

      cout << "Assigning statistics values (2):\n";
      if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      s3.assignReadBaseStatistics_MultiThread(AS_readpool,skim_params.sk_numthreads, false,
					      true,  // calc kmer forks
					      0,     // kmerforks: minkmer=0, so take everything!
					      false  // kmerforks: we don't need fwd/rev, relaxed
	);
      if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

      // method: look for gaps in valid status (excluding endgaps)
      // if none, perfect, set flag so that read will not be killed via invalid kmer ends
      for(size_t rpi=0; rpi<AS_readpool.size(); ++rpi){
	auto & actread=AS_readpool[rpi];
	CEBUG("noestkill check " << actread.getName() << endl);
	if(actread.getLenClippedSeq() < basesperhash) continue;

	int32 xpos=actread.getLeftClipoff();
	bool hasvalid=false;
	bool hasgap=false;
	auto bhsI=actread.getBPosHashStats().cbegin();
	bhsI+=xpos;
	auto bhsE=bhsI+actread.getLenClippedSeq()-basesperhash+1;
	// look for first valid
	for(; bhsI!=bhsE; ++bhsI){
	  if(bhsI->fwd.isValid()) {
	    hasvalid=true;
	    CEBUG("noestkill valid1 at " << (bhsI-actread.getBPosHashStats().cbegin()) << endl);
	    break;
	  }
	}
	// look for first non-valid
	for(; bhsI!=bhsE; ++bhsI){
	  if(!bhsI->fwd.isValid()) {
	    CEBUG("noestkill invalid at " << (bhsI-actread.getBPosHashStats().cbegin()) << endl);
	    ++bhsI;
	    break;
	  }
	}
	// look for first valid
	for(; bhsI!=bhsE; ++bhsI){
	  if(bhsI->fwd.isValid()) {
	    hasgap=true;
	    CEBUG("noestkill gap end detect at " << (bhsI-actread.getBPosHashStats().cbegin()) << endl);
	    break;
	  }
	}
	CEBUG("noestkill check result " << hasvalid << " " << hasgap << endl);
	CEBUG(""; Read::setCoutType(Read::AS_TEXT); cout << actread);
	if(hasvalid && !hasgap){
	  CEBUG("noestkill OK " << AS_readpool[rpi].getName() << endl);
	  AS_estnochimerakill[rpi]=1;
	}
      }

      if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

      if(s3.getNumHashEntries() != oldhssize){
	cout << "Reloading statistics 3 " << loadhsfn << " ... "; cout.flush();
	s3.discard();
	s3.loadHashStatistics(loadhsfn);
	cout << "done.\n";
	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      }
    }

    // very rough estimator of coverage for contigs
    // works well for kmer size between 17 and ~60 on medium well to well covered projects
    // But force the first pass to be taken nevertheless
    if(version == 1 || (basesperhash>=15 && basesperhash<=65)){
      AS_hashstat_avghashfreq=s3.getAvgHashFreqRaw();
      cout << "Estimator of average coverage: " << (AS_hashstat_avghashfreq+basesperhash)/3 << endl;
      AS_assemblyinfo.setLargeContigCovForStats((AS_hashstat_avghashfreq+basesperhash)/3);
    }

    cout << "Assigning statistics values (3):\n";
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    {
      uint32 mincountkmerforks=3;
      if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	mincountkmerforks=AS_hashstat_avghashfreq/6; // /6 to conservative?? TODO: check as HashStat now does fwd and rev separately!
	if(mincountkmerforks<3) mincountkmerforks=3;
      }
      s3.assignReadBaseStatistics_MultiThread(AS_readpool, skim_params.sk_numthreads, masknastyrepeats,
					      true,  // calc kmer forks
					      mincountkmerforks,
					      true                   // kmerforks: need fwd/rev
	);
    }
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

    AS_dataprocessing.buntifyReadsByHashFreq_Pool(AS_readpool,basesperhash);
    if(basesperhash>=17){
      AS_dataprocessing.addKMerForkTags_Pool(AS_readpool,basesperhash);
    }
  }

  if(0){
    cout << "AFTER\n";
    for(uint32 actid=0; actid<AS_readpool.size(); ++actid){
      Read & r=AS_readpool[actid];
      Read::setCoutType(Read::AS_TEXT);
      if(1 || r.getName()=="mictbac_bg:1:2101:20987:13472"
         ){
        cout << r;
      }
    }
    exit(101);
  }

  //if(nastyrepeatratio){
  if(hs_params.hs_repeatlevel_in_infofile){
    std::string filename;

    if(logname.size()){
      filename=buildFileName(version, prefix, postfix, logname, "");
    }else{
      //filename=buildFileName(version, prefix, postfix,
      //			     as_fixparams.as_outfile_stats_readrepeats,
      //			     ".lst");

      //filename=buildDefaultInfoFileName(version, prefix, postfix,
      filename=buildDefaultInfoFileName(-1, "", "",
					"",
					as_fixparams.as_outfile_stats_readrepeats,
					".lst",
					true);
    }

    cout << "Writing read repeat info to: " << filename << " ... ";
    cout.flush();

    uint32 howmanys=0;
    uint32 howmanyt=0;
    uint32 repanalysislevel=hs_params.hs_repeatlevel_in_infofile;
    if(repanalysislevel<5) repanalysislevel=5;
    if(repanalysislevel>8) repanalysislevel=8;

    std::ofstream fout(filename, std::ios::out|std::ios::trunc);
    for(uint32 rpi=0; rpi<AS_readpool.size(); rpi++){
      Read & actread= AS_readpool.getRead(rpi);
      if(!actread.hasValidData()
	 || !actread.isUsedInAssembly()) continue;
      bool mustshow=false;
      if(actread.hasTag(Read::REA_tagentry_idHAF5,-1)) {
	if(repanalysislevel==5) mustshow=true;
      }else if(actread.hasTag(Read::REA_tagentry_idHAF6,-1)) {
	if(repanalysislevel<=6) mustshow=true;
      }else if(actread.hasTag(Read::REA_tagentry_idHAF7,-1)) {
	if(repanalysislevel<=7) mustshow=true;
      }else if(actread.hasTag(Read::REA_tagentry_idMNRr,-1)) {
	if(repanalysislevel<=8) mustshow=true;
      }
      if(mustshow){
	bool countedthisseq=false;
	for(uint32 tn=0; tn<actread.getNumOfTags(); tn++){
	  const multitag_t & acttag=actread.getTag(tn);
	  if(acttag.to-acttag.from +1 >= basesperhash){
	    mustshow=false;
	    if(acttag.identifier==Read::REA_tagentry_idHAF5) {
	      if(repanalysislevel==5) mustshow=true;
	    }else if(acttag.identifier==Read::REA_tagentry_idHAF6) {
	      if(repanalysislevel<=6) mustshow=true;
	    }else if(acttag.identifier==Read::REA_tagentry_idHAF7) {
	      if(repanalysislevel<=7) mustshow=true;
	    }else if(acttag.identifier==Read::REA_tagentry_idMNRr) {
	      if(repanalysislevel<=8) mustshow=true;
	    }
	    if(mustshow){
	      if(!countedthisseq){
		countedthisseq=true;
		++howmanys;
	      }
	      ++howmanyt;
	      fout << actread.getName() << '\t'
		   << acttag.getIdentifierStr() << '\t';
	      for(uint32 readpos=acttag.from; readpos<=acttag.to; readpos++){
		fout << static_cast<char>(toupper(actread.getBaseInSequence(readpos)));
	      }
	      fout << '\n';
	    }
	  }
	}
      }
    }

    cout << howmanys << " sequences with " << howmanyt << " masked stretches." << endl;
  }

  // quick check for estimated coverage in genome data
  if(!AS_donequickdenovocoveragecheck){
    // atm only for de-novo, think about doing it for mapping
    // (though there may be good reasons for high coverage in mappings)
    if(basesperhash>=17
       && warnAtHighCoverages(AS_hashstat_avghashfreq)
       && AS_miraparams[0].getNagAndWarnParams().nw_check_coverage==NWSTOP){
      MIRANOTIFY(Notify::FOOL,"High average coverage detected, see output log above respectively the 'WARNING' files in the info directory for more information. In case you wish to force MIRA to disregard this safety check, consider using '-NW:cac=warn' or '-NW:cac=no'");
    }
    AS_donequickdenovocoveragecheck=true;
  }

  if(0){
    AS_dataprocessing.performKMERRepeatTagging_Pool(AS_readpool,basesperhash);
  }

  if(0){
    Read::setCoutType(Read::AS_TEXT);
    cout << "AFTER2\n";
    for(uint32 actid=0; actid<AS_readpool.size(); ++actid){
      Read & r=AS_readpool.getRead(actid);
      cout << r;
    }
    exit(101);
  }

  if(hs_params.hs_masknastyrepeats && hs_params.hs_apply_digitalnormalisation){
    if(version==AS_applydiginorminpass){
      if(as_fixparams.as_dateoutput) dateStamp(cout);
      cout << "Performing digital normalisation: "; cout.flush();
      AS_dataprocessing.performDigitalNormalisation_Pool(AS_readpool,s3,&AS_debrisreason);
      cout << "done\n";
      if(as_fixparams.as_dateoutput) dateStamp(cout);
    }
    // also for later steps: if a read has diginorm data, throw out MNRr tags
    for(uint32 rpi=0; rpi<AS_readpool.size(); ++rpi){
      auto & actread=AS_readpool[rpi];
      if(actread.hasTag(Read::REA_tagentry_idDGNr)) actread.deleteTag(Read::REA_tagentry_idMNRr);
    }
  }

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::priv_performHashAnalysis(const std::string & kmerfilename, bool usesignal, bool rarekmerfinalkill, int32 version, const std::string prefix, const std::string postfix, const std::string logname)
{
  FUNCSTART("void Assembly::performHashAnalysis()");

  uint32 basesperhash=AS_miraparams[0].getSkimParams().sk_basesperhash;

  std::string filenameforkms;
  std::string signalfile;
  if(kmerfilename.empty()){
    filenameforkms=buildFileName(version, prefix, postfix,
				 AS_miraparams[0].getAssemblyParams().as_tmpf_kmerstatistics,
				 ".mhs.gz");
  }else{
    filenameforkms=kmerfilename;
  }

  if(usesignal){
    signalfile=buildFileName(version, prefix, postfix,
			     AS_miraparams[0].getAssemblyParams().as_tmpf_signal_kmerstats,
			     ".ok");
  }

  if(basesperhash<=32){
    priv_phahelper<vhash64_t>(filenameforkms,signalfile,basesperhash,rarekmerfinalkill,version,prefix,postfix,logname);
  }else if(basesperhash<=64){
    priv_phahelper<vhash128_t>(filenameforkms,signalfile,basesperhash,rarekmerfinalkill,version,prefix,postfix,logname);
  }else if(basesperhash<=128){
    priv_phahelper<vhash256_t>(filenameforkms,signalfile,basesperhash,rarekmerfinalkill,version,prefix,postfix,logname);
  }else if(basesperhash<=256){
    priv_phahelper<vhash512_t>(filenameforkms,signalfile,basesperhash,rarekmerfinalkill,version,prefix,postfix,logname);
  }else{
    MIRANOTIFY(Notify::FATAL,"Cannot perform a hash analysis with -SK:bph > 256 (you used " << basesperhash << "), though MIRA should've failed earlier, I admit.\n");
  }

  //CEBUG("BEFORE\n");
  //for(uint32 actid=0; actid<AS_readpool.size(); actid++){
  //  Read & r=AS_readpool.getRead(actid);
  //  r.integrityCheck();
  //  Read::setCoutType(Read::AS_TEXT);
  //  cout << r;
  //}


  if(AS_logflag_dumphashanalysis){
    std::string logfilename=buildFileName(version, "", "",
				     "elog.dp.hashanalysis_pass",
				     ".lst");
    //std::string logfilename=AS_miraparams[0].getDirectoryParams().dir_tmp+"/elog.dp.hashanalysis.lst";

    cout << "elog hashan: " << logfilename << endl;
    std::ofstream logfout(logfilename, std::ios::out|std::ios::trunc);

    for(uint32 rpi=0; rpi<AS_readpool.size(); rpi++){
      Read::setCoutType(Read::AS_TEXT);
      logfout << AS_readpool[rpi];
    }
  }

  FUNCEND();
  return;
}









/*************************************************************************
 *
 * expects reads to have baseflags set  (by performHashAnalysis())
 *
 *
 *************************************************************************/

// switch left/right?
// perform twice, then look for better?
//rs_pecimp2 rs:2108:4766:166384/1


//#define CEBUG(bla)   {cout << bla; cout.flush();}

uint64 Assembly::performProposedEndClips(uint32 basesperhash, const std::string & logname, const std::string & logprefix)
{
  FUNCSTART("void Assembly::performProposedEndClips(const std::string & logname, const std::string & logprefix)");

  const static std::string ggcproblem("ggc");

  cout << "Looking for proposed cutbacks ... "; cout.flush();

  std::ofstream logfout;
  std::vector<int32> lclip;
  std::vector<int32> rclip;

  lclip.resize(AS_readpool.size());
  rclip.resize(AS_readpool.size());
  for(uint32 rpi=0; rpi<AS_readpool.size(); ++rpi){
    lclip[rpi] = AS_readpool[rpi].getLeftClipoff();
    rclip[rpi] = AS_readpool[rpi].getRightClipoff();
  }

  if(!logname.empty()){
    logfout.open(logname, std::ios::out|std::ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  AS_dataprocessing.proposedEndClipping_Pool(AS_readpool, basesperhash);

  uint64 numbasesclipped=0;

  uint32 cbleft=0;
  uint32 cbright=0;
  uint32 killed=0;

  for(uint32 rpi=0; rpi<AS_readpool.size(); ++rpi){
    Read & actread=AS_readpool[rpi];
    if(lclip[rpi] != actread.getLeftClipoff()
       || rclip[rpi] != actread.getRightClipoff()){

      numbasesclipped+=rclip[rpi]-lclip[rpi]-actread.getLenClippedSeq();

      if(lclip[rpi] != actread.getLeftClipoff()){
	++cbleft;
	logfout << logprefix << " left "
		<< actread.getName() << '\t'
		<< lclip[rpi]
		<< " -> "  << actread.getLeftClipoff()
		<< "\n";
      }

      if(rclip[rpi] != actread.getRightClipoff()){
	++cbright;
	logfout << logprefix << " right "
		<< actread.getName() << '\t'
		<< rclip[rpi]
		<< " -> " << actread.getRightClipoff()
		<< "\n";
      }

      if(actread.getLenClippedSeq() < AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength){
	++killed;
	logfout << logprefix << " "
		<< actread.getName() << " killed, remaining length ("
		<< actread.getLenClippedSeq() << ")\n";
	if(rpi < AS_debrisreason.size()
	   && AS_debrisreason[rpi]==DEBRIS_NOTDEBRIS){
	  AS_debrisreason[rpi]=DEBRIS_CLIP_PROPOSEDENDCLIP;
	}
      }
    }
  }

  //Read::setCoutType(Read::AS_TEXT);

  logfout.close();

  cout << "done.\nPerformed clips:"
       << "\n\tNum reads cliped left      : " << cbleft
       << "\n\tNum reads cliped right     : " << cbright
       << "\n\tNum reads completely killed: " << killed
       << "\n\tTotal bases clipped        : " << numbasesclipped
       << "\n\n";

  return numbasesclipped;
}
//#define CEBUG(bla)


uint64 Assembly::performPECCHIMSRERKM(const std::string & logname, const std::string & logprefix)
{
  FUNCSTART("void Assembly::performPECCHIMSRERKM(const std::string & logname, const std::string & logprefix)");

  bool dopec=false;
  bool dorarekmermask=false;
  for(auto mpi=1; mpi < AS_miraparams.size(); ++mpi){
    dopec|=AS_seqtypespresent[mpi] && AS_miraparams[mpi].getAssemblyParams().as_clip_proposeendclips;
    dorarekmermask|=AS_seqtypespresent[mpi] && AS_miraparams[mpi].getAssemblyParams().as_clipmask_rarekmers;
  }

  bool dosdbgchimera=AS_miraparams[0].getAssemblyParams().as_clip_sdbg_chimeradetection;
  bool dosdbgedit=AS_miraparams[0].getEditParams().ed_sdbg_readedit;

  if(!dopec && !dosdbgchimera&& !dosdbgedit && !dorarekmermask) return 0;

  cout << "Hash analysis for ";
  {
    bool needand=false;;
    if(dopec) {
      cout << "proposed cutbacks";
      needand=true;
    }
    if(dosdbgchimera) {
      if(needand) cout << " and ";
      needand=true;
      cout << "chimera search";
    }
    if(dosdbgedit) {
      if(needand) cout << " and ";
      needand=true;
      cout << "read editing";
    }
    if(dorarekmermask) {
      if(needand) cout << " and ";
      cout << "rare kmer masking";
    }
    cout << ":\n";
  }

  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  hashstatistics_parameters const & hs_params= AS_miraparams[0].getHashStatisticsParams();

  uint32 basesperhash=as_fixparams.as_clip_pec_basesperhash;
  if(basesperhash > 32) {
    // BaCh 11.05.2016: ummmm ... now I can do with larger, can't I?
    MIRANOTIFY(Notify::FATAL,"-CL:pecbph can run only with up to 32, you set it to " << basesperhash << ". MIRA should've warned earlier though.\n");
  }

  {
    std::string merfile(AS_miraparams[0].getDirectoryParams().dir_tmp+"/peckmerstat.mhs.gz");
    HashStatistics<vhash64_t> s3;

    s3.setHashFrequencyRatios(hs_params.hs_freqest_minnormal,
			      hs_params.hs_freqest_maxnormal,
			      hs_params.hs_freqest_repeat,
			      hs_params.hs_freqest_heavyrepeat,
			      hs_params.hs_freqest_crazyrepeat,
			      2,                                  // rare kmer count: occurrence <= that can be masked away later by rare kmer masking
			      hs_params.hs_nastyrepeatratio,
			      hs_params.hs_nastyrepeatcoverage);

    std::string filenameforks(merfile);

    // if we do a chimera search OR some SDBG edits, do NOT take rails in to the
    //  calculation of hash statistics! These rails could have come from a "wrong"
    //  backbone and therefore led to wrong chimera recognition or read edits!
    // if neither is done, we can take the rails as "maybe, possibly true" to help in
    //  areas difficult for sequencing (e.g. GGCxG in Illuminas)
    bool alsorailsinhs=!(dosdbgchimera | dosdbgedit);

    s3.computeHashStatistics(AS_readpool,
			     hs_params.hs_memtouse,
			     true,
			     alsorailsinhs,
			     true,
			     as_fixparams.as_clip_pec_mkfr,
			     as_fixparams.as_clip_pec_mtk, // rare kmer early kill
			     basesperhash,
			     filenameforks,
			     AS_miraparams[0].getDirectoryParams().dir_tmp
      );
    s3.showHashStatisticsInfo();

    auto avghashfreq=s3.getAvgHashFreqRaw();

    // do this before assigning read statistics
    // chimera would be OK after, but performing SDBG edits currently trashes *ALL* read tags!
    if(dosdbgchimera || dosdbgedit){
      if(dosdbgchimera){
	uint32 trimfreq=3;
	if(avghashfreq>60) trimfreq=4;
	if(avghashfreq>80) trimfreq=5;

	std::string tmpfname;
	tmpfname=buildFileName(0,"","",
			       as_fixparams.as_tmpf_clippings,
			       ".txt","",false);
	AS_dataprocessing.startLogging(tmpfname,false);
	cout << "Hunting down chimeras:\n";
	auto hits=AS_dataprocessing.performSDBGChimeraSearch_Pool(AS_readpool,s3,trimfreq,&AS_debrisreason,"PEC");
	cout << "\nPEC chimera trophy count: " << hits << endl;
	AS_dataprocessing.stopLogging();

	BUGIFTHROW(AS_current_rls_byrg.size() != ReadGroupLib::getNumReadGroups(), "AS_current_rls_byrg.size() " << AS_current_rls_byrg.size() << " != ReadGroupLib::getNumReadGroups() " << ReadGroupLib::getNumReadGroups() << " ???");
	BUGIFTHROW(AS_current_rls_bytype.size() != ReadGroupLib::getNumSequencingTypes(), "AS_current_rls_bytype.size() " << AS_current_rls_bytype.size() << " != ReadGroupLib::getNumSequencingTypes() " << ReadGroupLib::getNumSequencingTypes() << " ???");

	for(auto & rls : AS_current_rls_bytype){
	  rls.RLS_count_chimera=0;
	}
	for(auto & rls : AS_current_rls_byrg){
	  rls.RLS_count_chimera=0;
	}
	for(uint32 rid=0; rid<AS_readpool.size(); ++rid){
	  if(AS_debrisreason[rid]==DEBRIS_CLIP_CHIMERA){
	    auto rgid=AS_readpool[rid].getReadGroupID().getLibId();
	    ++AS_current_rls_byrg[rgid].RLS_count_chimera;
	    ++AS_current_rls_bytype[AS_readpool[rid].getSequencingType()].RLS_count_chimera;
	  }
	}

	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
	cout << "Reloading statistics 3 " << merfile << " ... "; cout.flush();
	s3.discard();
	s3.loadHashStatistics(merfile);
	cout << "done.\n";
	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      }
      if(dosdbgedit){
	cout << "Performing PEC SDBG edits:\n";
	auto hits=AS_dataprocessing.performSDBGEdits_Pool(AS_readpool,s3,1); // trimfreq==1 -> very conservative
	cout << "\nPEC SDBG edited reads: " << hits << endl;
	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
	cout << "Reloading statistics 4 " << merfile << " ... "; cout.flush();
	s3.discard();
	s3.loadHashStatistics(merfile);
	cout << "done.\n";
	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      }
    }

    cout << "Assigning statistics values (1):\n";
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    s3.assignReadBaseStatistics_MultiThread(AS_readpool,skim_params.sk_numthreads, false,
					    true,  // calc kmer forks
					    0,     // kmerforks: minkmer=0, so take everything!
					    false  // kmerforks: we don't need fwd/rev, relaxed
      );
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

    if(basesperhash>=17
       && as_fixparams.as_clip_pec_mkfr <2){
      if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
	 && AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA] && avghashfreq >=50){
	// BaCh 11.04.2015: bad idea for projects with quite uneven coverage
        //cout << "Detected probable higher coverage in Illumina genome project, setting: -CL:pmkfr=2\n";
	//const_cast<assembly_parameters &>(as_fixparams).as_clip_pec_mkfr=2;
      }
    }

    // quick check for estimated coverage in genome data
    if(!AS_donequickdenovocoveragecheck){
      // atm only for de-novo, think about doing it for mapping
      // (though there may be good reasons for high coverage in mappings)
      if(basesperhash>=17
	 && warnAtHighCoverages(avghashfreq)
	 && AS_miraparams[0].getNagAndWarnParams().nw_check_coverage==NWSTOP){
	MIRANOTIFY(Notify::FATAL,"High average coverage detected, see output log above respectively the 'WARNING' files in the info directory for more information. In case you wish to force MIRA to disregard this safety check, consider using '-NW:cac=warn' or '-NW:cac=no'");
      }
      // Nope, do not set that here. The hash statistics of the main loop have the final word!
      // AS_donequickdenovocoveragecheck=true;
    }
  }

  if(as_fixparams.as_dateoutput) dateStamp(cout);
  cout << '\n';

  uint64 numbasesclipped=0;
  if(dopec){
    numbasesclipped=performProposedEndClips(as_fixparams.as_clip_pec_basesperhash,logname,logprefix);
  }

  if(1){
    AS_dataprocessing.startLogging(logname,false);

    AS_dataprocessing.performRareKMERMasking_Pool(AS_readpool,basesperhash,logprefix);
    AS_dataprocessing.stopLogging();
  }

  AS_dataprocessing.clipPolyBaseAtEnd_Pool(AS_readpool,logprefix);

  FUNCEND();

  return numbasesclipped;
}
//#define CEBUG(bla)






/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::cutBackPossibleChimeras(const std::string & logname, const std::string & logprefix, const std::vector<int32> & chuntleftcut, const std::vector<int32> & chuntrightcut, std::vector<bool> & chimeracutflag)
{
  FUNCSTART("void Assembly::cutBackPossibleChimeras(const std::string & logname, const std::string & logprefix, const std::vector<int32> & chuntleftcut, const std::vector<int32> & chuntrightcut)");

  BUGIFTHROW(chuntleftcut.size()!=chuntrightcut.size() && chuntleftcut.size() != AS_readpool.size(),"Arrays mismatch? chuntleftcut.size()!=chuntrightcut.size && chuntleftcut.size() != AS_readpool.size()");

  std::ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname, std::ios::out|std::ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  cout << "Cutting back possible chimeras ... "; cout.flush();

  if(!chimeracutflag.empty()){
    chimeracutflag.clear();
    chimeracutflag.resize(chuntleftcut.size(),false);
  }

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  for(uint32 actreadid=0;actreadid<chuntleftcut.size();actreadid++){
    Read & actread=AS_readpool.getRead(actreadid);
    if(actread.hasValidData()
       && !(actread.isBackbone()
	    || actread.isRail())){
      bool didcut=false;
      if(as_fixparams.as_clip_skimchimeradetection
	 && (chuntleftcut[actreadid]>0
	     || chuntrightcut[actreadid]>0)){
	logfout << logprefix << " possible chimera: " << actread.getName()
		<< "\t["
		<< actread.getLeftClipoff()
		<< ","
		<< actread.getRightClipoff()
		<< "[ using cfrag " << chuntleftcut[actreadid] << ":" << chuntrightcut[actreadid]
		<< " cut back to ";

	actread.setLSClipoff(actread.getLeftClipoff()+chuntleftcut[actreadid]);
	actread.setRSClipoff(actread.getLeftClipoff()+(chuntrightcut[actreadid]-chuntleftcut[actreadid])+1);
	didcut=true;
	if(!chimeracutflag.empty()){
	  chimeracutflag[actreadid]=true;
	}

	logfout << '['
		<< actread.getLeftClipoff()
		<< ","
		<< actread.getRightClipoff()
		<< "[\n";
      }

      if(!didcut
	 && (chuntleftcut[actreadid]<0
	     || chuntrightcut[actreadid]<0)){
	if(as_fixparams.as_clip_skimjunkdetection){
	  logfout << logprefix << " removed possible junk: " ;
	}else{
	  logfout << logprefix << " untouched possible junk: " ;
	}
	logfout << actread.getName()
		<< "\t["
		<< -chuntleftcut[actreadid]
		<< ","
		<< -chuntrightcut[actreadid]
		<< '\n';
	if(as_fixparams.as_clip_skimjunkdetection){
	  actread.setLSClipoff(actread.getLeftClipoff()-chuntleftcut[actreadid]);
	  actread.setRSClipoff(actread.getRightClipoff()+chuntrightcut[actreadid]);
	  if(!chimeracutflag.empty()){
	    chimeracutflag[actreadid]=true;
	  }
	}
      }
    }
  }

  cout << "done.\n";
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
/*
void Assembly::performPool_AdaptorRightClip(const std::string & logname, const std::string & logprefix, const uint8 seqtype)
{
  FUNCSTART("void Assembly::performPool_SolexaAdaptorRightClip(const std::string & logname, const std::string & logprefix, const uint8 seqtype);)");

// BOOST: regex not compatible with _GLIBCXX_DEBUG
#ifdef _GLIBCXX_DEBUG
  cout << "_GLIBCXX_DEBUG not compatible with std::regex :-(\n";
  return;
#endif

  BUGIFTHROW(seqtype>=ReadGroupLib::SEQTYPE_END,"Unknown seqtype " << static_cast<uint16>(seqtype) << "given.");

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  struct masterslavere_t {
    std::regex masterre;
    std::list<std::regex> slaveres;
    bool hasmaster;

    masterslavere_t(): hasmaster(false) {};
  };

  // prepare regular expressions
  std::list<masterslavere_t> adapres;
  {
    istd::stringstream tmpis;
    if(seqtype==ReadGroupLib::SEQTYPE_SOLEXA){
      static const char regexfile[] = {
#include "adaptorsregex.solexa.xxd.H"
	,0
      };
      tmpis.str(regexfile);
    }else if(seqtype==ReadGroupLib::SEQTYPE_IONTORRENT){
      static const char regexfile[] = {
#include "adaptorsregex.iontor.xxd.H"
	,0
      };
      tmpis.str(regexfile);
    }

    masterslavere_t tmpmsre;
    std::string line;

    while(true){
      getline(tmpis,line);
      if(tmpis.eof()) break;
      if(line[0]=='>'){
	adapres.push_back(tmpmsre);
	line.erase(0,1);         // get away the ">"
	boost::trim(line);
	if(!line.empty()){
	  boost::to_upper(line);
	  adapres.back().masterre=std::regex(line);
	  adapres.back().hasmaster=true;
	}
      }else{
	BUGIFTHROW(adapres.empty(),"Oooops, found no master expression?");
	boost::to_upper(line);
	adapres.back().slaveres.push_back(std::regex(line));
      }
    }
  }

  ReadPool adappool(&AS_miraparams);
  {
    istd::stringstream tmpis;

    if(seqtype==ReadGroupLib::SEQTYPE_SOLEXA){
      static const char adapfile[] = {
#include "adaptorsforclip.solexa.xxd.H"
	,0
      };
      tmpis.str(adapfile);
    }else if(seqtype==ReadGroupLib::SEQTYPE_IONTORRENT){
      static const char adapfile[] = {
#include "adaptorsforclip.iontor.xxd.H"
	,0
      };
      tmpis.str(adapfile);
    }else if(seqtype==ReadGroupLib::SEQTYPE_454GS20){
      static const char adapfile[] = {
#include "adaptorsforclip.454.xxd.H"
	,0
      };
      tmpis.str(adapfile);
    }

    std::string line;
    while(true){
      getline(tmpis,line);
      if(tmpis.eof()) break;
      line.erase(0,1);         // get away the ">"
      if(!line.empty()){
	size_t ereadidx=adappool.provideEmptyRead();
	Read & actread=adappool[ereadidx];
	actread.disallowAdjustments();
	actread.setName(line);
	getline(tmpis,line);
	if(tmpis.eof()) break;
	actread.setSequenceFromStd::String(line);
      }
    }
  }

  //adappool.dumpPoolInfo(cout);

  // Go back if nothing to be searched
  if(adappool.size()==0 && adapres.size()==0) return;

  cout << "Starting " << ReadGroupLib::getNameOfSequencingType(seqtype) << " known adaptor right clip ... "; cout.flush();

  Skim adapskim;
  adapskim.skimStreamPrepare(adappool,7,1);

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "Searching multithread now ... \n"; cout.flush();

  cout << static_cast<int16>(AS_miraparams[0].getSkimParams().sk_numthreads) << endl;

  std::vector<int32> clipres;
  adapskim.findAdaptorRightClip(AS_readpool,clipres,seqtype,9,AS_miraparams[0].getSkimParams().sk_numthreads);
  //adapskim.findAdaptorRightClip(AS_readpool,clipres,seqtype,9,1);
  //adapskim.findAdaptorRightClip(AS_readpool,clipres,seqtype,9,8);

  BUGIFTHROW(clipres.size()!=AS_readpool.size(),"clipres.size()!=AS_readpool.size()???");

  ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname, std::ios::out|std::ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  uint32 numclipped=0;

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "Searching for " <<  ReadGroupLib::getNameOfSequencingType(seqtype) << " partial end adaptors ... \n"; cout.flush();
  ProgressIndicator<int64> P(0, AS_readpool.size());
  for(uint32 actid=0; actid < AS_readpool.size(); actid++){
    P.progress(actid);
    Read & actread = AS_readpool.getRead(actid);
    if(actread.hasValidData()
       && actread.getSequencingType()==seqtype
       && !(actread.isBackbone() || actread.isRail())){

      auto oldrsclip=actread.getRSClipoff();
      if(clipres[actid]>=0){
	if(clipres[actid] < oldrsclip){
	  ++numclipped;
	  actread.setRSClipoff(clipres[actid]);
	  logfout << logprefix << " "
		  << ReadGroupLib::getNameOfSequencingType(seqtype)
		  << " adaptor: " << actread.getName()
		  << " changed right clip from " << oldrsclip << " to " << clipres[actid] << "\n";
	}
      }else if(!adapres.empty()){
	std::string seq(actread.getSeqAsChar());
	boost::to_upper(seq);

	auto flags = boost::match_default;
	boost::match_results<std::string::const_iterator> what;
	std::string::const_iterator start, end;

	for(auto & msre : adapres){
	  bool dosearch=true;
	  if(msre.hasmaster){
	    if(!regex_search(start, end, what, msre.masterre, flags)) {
	      dosearch=false;
	    }
	  }
	  bool breakit=false;
	  if(dosearch){
	    for(auto & thisre : msre.slaveres){
	      start = seq.begin();
	      end = seq.end();
	      if(regex_search(start, end, what, thisre, flags)) {
		if(what.position()< oldrsclip){
		  actread.setRSClipoff(what.position());
		  logfout << logprefix << " "
			  << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
			  << " partial end adaptor: " << actread.getName()
			  << " changed right clip from " << oldrsclip << " to " << what.position() << "\n";
		  breakit=true;
		  break;
		}
	      }
	    }
	  }
	  if(breakit) break;
	}

      }
    }
  }

  P.finishAtOnce();

  cout << "done. Clipped " << numclipped << " reads.\n";

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

  FUNCEND();
  return;
}

*/


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::correctContigs()
{
#ifdef MIRA_HAS_EDIT
  FUNCSTART("void Assembly::correctContigs()");

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "\nEditing contigs:" << endl;

  EDITParameters eparams;

  //  eparams.setDoEval();
  eparams.setStrictEvaluation(false);
  eparams.setConfirmationThreshold(0.5);
  eparams.setShowProgress(true);
  eparams.setVerbose(0);
  eparams.setShowProgress(true);


  int32 ccounter=0;
  ProgressIndicator<int64> P(0, AS_contigs.size());
  for(auto & cle : AS_contigs){
    P.progress(ccounter);
    try {
      //	CEBUG("Editing contig:" << ccounter << endl);
      //	CEBUG(cle);
      cout << "Editing contig:" << ccounter << endl;
      editContigBack(cle, eparams);
      ScfBuffer::discard();
      cout << "deleting star columns" << ccounter << endl;
      cle.deleteStarOnlyColumns(0, cle.getContigLength()-1);
      cout << "marking repeats" << ccounter << endl;

      Contig::repeatmarker_stats_t repstats;
      std::vector<bool> readsmarkedsrm;
      cle.newMarkPossibleRepeats(repstats, readsmarkedsrm);

      //	CEBUG("Corrected contig:" << endl);
      //	CEBUG(cle);
    }
    catch(Notify n){
      n.handleError("Error while examining fault-region");
    }

    I++;ccounter++;
  }

  P.finishAtOnce();

  cout << endl;

  FUNCEND();
#endif
  return;
}






/*************************************************************************
 *
 * Calculates possible sequence vector leftovers at the left side of a read
 * Reads that get a clip must be of Sanger type
 *
 * Does not clip backbone reads, rail reads, multicopyreads
 *  AND not areas protected by Staden GenBank Feature tags
 *
 * Clipping itself must be done afterwards in the performSeqVectorClippings()
 *  function. This was split in two parts to allow releasing of the
 *  big memory chunks AS_readhmcovered, AS_readhitmiss, etc.
 *
 *************************************************************************/


void Assembly::calcPossibleSeqVectorClipoffs(int32 version, const std::string prefix, const std::string postfix, const std::string logname)
{
  FUNCSTART("void Assembly::calcPossibleSeqVectorClipoffs(int32 version, const std::string prefix, const std::string postfix, const std::string logname)");

  if(AS_readhmcovered.size()==0 || AS_readhitmiss.size()==0) {
    cout << "\nNo vector clipping information available, aborting vector clip.\n";
    FUNCEND();
    return;
  }

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "\nCalculating possible vector leftovers ... ";
  cout.flush();
  //ProgressIndicator P (0, AS_readhmcovered.size()-1);

  AS_clipleft.clear();
  AS_clipright.clear();
  AS_clipleft.resize(AS_readhmcovered.size(),-1);
  AS_clipright.resize(AS_readhmcovered.size(),-1);

  std::string filename;
  if(logname.size()){
    filename=buildFileName(version, prefix, postfix, logname, ".txt");
  }else{
    filename=buildFileName(version, prefix, postfix,
			   AS_miraparams[0].getAssemblyParams().as_tmpf_vectorclip,
			   ".txt");
  }

  std::ofstream logout(filename, std::ios::out | std::ios::trunc);

  for(uint32 id=0; id<AS_readhmcovered.size(); id++) {
    if(AS_readpool.getRead(id).getSequencingType() != ReadGroupLib::SEQTYPE_SANGER
       || AS_readpool.getRead(id).isBackbone()
       || AS_readpool.getRead(id).isRail()
       || AS_multicopies[id]>0
      ) continue;


    //P.progress(id);

    uint32 clippos=0;
    bool mustclip=false;
    for(uint32 actpos=0; actpos<AS_readhmcovered[id].size(); actpos++) {
      if(actpos-clippos > 5) break;
      if(AS_readhmcovered[id][actpos]>=4) {
	if(AS_readhitmiss[id][actpos]) {
	  if(100.0/static_cast<double>(AS_readhmcovered[id][actpos])*static_cast<double>(AS_readhitmiss[id][actpos]) >= 30.0) {
	    clippos=actpos;
	    mustclip=true;
	  }
	}
      }
    }
    clippos++;

    // check that no GenBank Feature tags protect the area, else clip less
    {

      // FIXME: put all checks for that into read.C (*sigh*)

      for(uint32 i=0; i<AS_readpool.getRead(id).getNumOfTags(); i++){
	const multitag_t & acttag=AS_readpool.getRead(id).getTag(i);
	if(!acttag.isSourceMIRA()){
	  if(acttag.from<clippos) clippos=acttag.from;
	  if(acttag.to<=clippos) clippos=0;
	}
      }
    }

    // auf clip verzichten wenn nur 1 base betroffen (sieht zu doof aus)
    if(mustclip && clippos>1) {
      uint32 maxcliplenallowed=AS_miraparams[AS_readpool.getRead(id).getSequencingType()].getAssemblyParams().as_clip_vector_maxlenallowed;
      if(maxcliplenallowed == 0 || clippos <= maxcliplenallowed) {
	//AS_readpool.getRead(id).setClipoffs(AS_readpool.getRead(id).getLeftClipoff()+clippos,
	//				    AS_readpool.getRead(id).getRightClipoff(),
	//				    false);

	//AS_clipleft[id]=AS_readpool.getRead(id).getLeftClipoff()+clippos;

	AS_clipleft[id]=clippos;

	logout << "Clipped " << clippos << " bases on the left of " << AS_readpool.getRead(id).getName() << "\n";

      } else {
	if(clippos > maxcliplenallowed) {
	  logout << "Not clipped " << clippos << " bases on the left of " << AS_readpool.getRead(id).getName() << " , too long.\n";
	}
      }
    }
  }

  cout << "done.\n";

  AS_steps[ASVECTORSCLIPPED]=1;
  AS_steps[ASADSLISTOK]=0;

  FUNCEND();
}




/*************************************************************************
 *
 * Reads must be Sanger type
 *
 *
 *************************************************************************/

void Assembly::performSeqVectorClippings()
{
  FUNCSTART("void Assembly::performSeqVectorClippings()");

  cout << "\nPerforming vector clipping ... ";
  cout.flush();

  for(uint32 id=0; id<AS_clipleft.size(); id++) {
    if(AS_clipleft[id]>=0
       && AS_readpool.getRead(id).isSequencingType(ReadGroupLib::SEQTYPE_SANGER)) {
      AS_readpool.getRead(id).setClipoffs(AS_readpool.getRead(id).getLeftClipoff()+AS_clipleft[id],
					  AS_readpool.getRead(id).getRightClipoff(),
					  false);
    }
  }
  FUNCEND();

  AS_clipleft.clear();

  cout << "done." << endl;

  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

struct cliplen_t{
  int32 len;
  bool changed;
};


//#define CEBUGFLAG 1
void Assembly::extendADS(int32 version, const std::string prefix, const std::string postfix, const std::string logname)
{
  FUNCSTART("void Assembly::extendADS(int32 version, const std::string prefix, const std::string postfix, const std::string logname)");

//  if(AS_steps[ASADSLISTOK]==0){
//    makeAlignments();
//  }


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  // TODO: change to use different Aligns / MIRAparams depending
  //   on Sanger / 454 (/ PacBio ???) reads

  // TODO: what about PacBio? currently not used, but should it?

  MIRAParameters tmpparams = AS_miraparams[0];

  const_cast<align_parameters &>(tmpparams.getAlignParams()).al_min_relscore=5;

  assembly_parameters const & as_params= tmpparams.getAssemblyParams();

  std::string filename;
  if(logname.size()){
    filename=buildFileName(version, prefix, postfix, logname, ".txt");
  }else{
    filename=buildFileName(version, prefix, postfix,
			   as_params.as_tmpf_adsextend,
			   ".txt");
  }

  std::ofstream logout(filename, std::ios::out | std::ios::trunc);

  std::vector<cliplen_t> clips(AS_readpool.size());
  for(uint32 i=0; i<clips.size(); i++){
    clips[i].len=0;
    clips[i].changed=false;
  }

  std::list<AlignedDualSeq> madsl;

  try{
    // true for using memcache
    Align bla(&tmpparams);

    cout << "\n";
    if(as_params.as_dateoutput) dateStamp(cout);
    cout << "\nSearching possible read extensions (for Sanger and/or 454):\n";

    ProgressIndicator<int32> P(0, static_cast<int32>(AS_adsfacts.size())-1);
    uint32 pindic=0;

    for(auto I = AS_adsfacts.cbegin(); I!=AS_adsfacts.cend(); ++I){
      P.progress(pindic++);
      // first try: prolongate to end.
      int32 id1=I->getID1();
      int32 id2=I->getID2();

      // no sense to calc read extensions for reads where both seqtypes are said
      //  not to use extensions
      if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension == false
	 && AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension == false) continue;

      if(AS_permanent_overlap_bans.checkIfBanned(id1,id2)) {
	CEBUG("PermBan for: " << id1 << " " << id2 <<"\tskipping\n");
	continue;
      }

      CEBUG("\n\nid1: " << id1 << "\t" << AS_readpool.getRead(id1).getName() <<endl);
      CEBUG("id2: " << id2 << "\t" << AS_readpool.getRead(id2).getName() <<endl);

      // normally the sequences should have a length >0
      // but due to some clipping being done after SKIM (chimera etc.), it
      //  may happen they are 0 now. If that's the case, don't bother
      //  looking at.
      if(AS_readpool[id1].getLenClippedSeq() == 0
	 || AS_readpool[id2].getLenClippedSeq() == 0) continue;

      // check for sequencing types
      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)) continue;
      // let's allow PacBio HQ

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)) continue;

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)) continue;

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_TEXT)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_TEXT)) continue;

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_ABISOLID)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_ABISOLID)) continue;

      //if(clips[id1].changed && clips[id2].changed){
      //	CEBUG(id1 << " and " << id2 <<" already changed.\n");
      //	continue;
      //}

      madsl.clear();

#if CEBUGFLAG > 0
      //Read::setCoutType(Read::AS_TEXT);
      Read::setCoutType(Read::AS_TEXTCLIPS);
      CEBUG(AS_readpool.getRead(id1));
      CEBUG(AS_readpool.getRead(id2));
#endif

      if(I->getSequenceDirection(id1) * I->getSequenceDirection(id2) > 0){

	CEBUG("doalign\n");

	// evil hack warning
	// the &(* ...) construction is needed for gcc3 as it cannot convert
	//  a std::vector<char> iterator to char *   (*sigh*)

	int32 extendlen1=AS_readpool.getRead(id1).getRightExtend();
	int32 extendlen2=AS_readpool.getRead(id2).getRightExtend();

	if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension == false) {
	  extendlen1=0;
	}
	if(AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension == false){
	  extendlen2=0;
	}

	CEBUG("l1: " <<AS_readpool.getRead(id1).getLenClippedSeq() << endl);
	CEBUG("e1: " <<extendlen1 << endl);
	CEBUG("l2: " <<AS_readpool.getRead(id2).getLenClippedSeq() << endl);
	CEBUG("e2: " <<extendlen2 << endl);

	if(extendlen1 >= 10 || extendlen2 >= 10){
	  bla.acquireSequences(
	    &(*AS_readpool.getRead(id1).getActualSequence().begin())
	    +AS_readpool.getRead(id1).getLeftClipoff(),
	    AS_readpool.getRead(id1).getLenClippedSeq()+extendlen1,
	    &(*AS_readpool.getRead(id2).getActualSequence().begin())
	    +AS_readpool.getRead(id2).getLeftClipoff(),
	    AS_readpool.getRead(id2).getLenClippedSeq()+extendlen2,
	    id1, id2, 1, 1, true, I->getOffsetInAlignment(id2));
	  bla.setEnforceCleanEnds(false);
	  bla.setAffineGapScore(false);
	  bla.fullAlign(&madsl);

	  if(madsl.size()==0){
	    CEBUG("No results, less radical try.\n");

	    int32 tryseqlen1=0;
	    if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension) {
	      if(clips[id1].changed){
		extendlen1-=clips[id1].len;
	      }
	      extendlen1/=2;
	      tryseqlen1=AS_readpool.getRead(id1).getLenClippedSeq()+extendlen1;
	      if(clips[id1].changed){
		tryseqlen1+=clips[id1].len;
	      }
	      if(AS_readpool.getRead(id1).getLeftClipoff()+tryseqlen1 >= static_cast<int32>(AS_readpool.getRead(id1).getLenSeq())) {
		CEBUG("t1o: " <<tryseqlen1 << endl);
		tryseqlen1=AS_readpool.getRead(id1).getLenClippedSeq()+AS_readpool.getRead(id1).getRightExtend();
		CEBUG("t1n: " <<tryseqlen1 << endl);
	      }
	    }

	    int32 tryseqlen2=0;
	    if(AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension) {
	      if(clips[id2].changed){
		extendlen2-=clips[id2].len;
	      }
	      extendlen2/=2;
	      tryseqlen2=AS_readpool.getRead(id2).getLenClippedSeq()+extendlen2;
	      if(clips[id2].changed){
		tryseqlen2+=clips[id2].len;
	      }
	      if(AS_readpool.getRead(id2).getLeftClipoff()+tryseqlen2 >= static_cast<int32>(AS_readpool.getRead(id2).getLenSeq())) {
		CEBUG("t2o: " <<tryseqlen2 << endl);
		tryseqlen2=AS_readpool.getRead(id2).getLenClippedSeq()+AS_readpool.getRead(id2).getRightExtend();
		CEBUG("t2n: " <<tryseqlen2 << endl);
	      }
	    }

	    CEBUG("cc1: " <<clips[id1].changed << endl);
	    CEBUG("cl1: " <<clips[id1].len << endl);
	    CEBUG("l1: " <<AS_readpool.getRead(id1).getLenClippedSeq() << endl);
	    CEBUG("t1: " <<tryseqlen1 << endl);
	    CEBUG("cc2: " <<clips[id2].changed << endl);
	    CEBUG("cl2: " <<clips[id2].len << endl);
	    CEBUG("l2: " <<AS_readpool.getRead(id2).getLenClippedSeq() << endl);
	    CEBUG("t2: " <<tryseqlen2 << endl);
	    if(extendlen1 < 5 && extendlen2 < 5) {
	      CEBUG("skip" << endl);
	      continue;
	    }

	    if(tryseqlen1>0 && tryseqlen2>0){
	      bla.acquireSequences(
		&(*AS_readpool.getRead(id1).getActualSequence().begin())
		+AS_readpool.getRead(id1).getLeftClipoff(),
		tryseqlen1,
		&(*AS_readpool.getRead(id2).getActualSequence().begin())
		+AS_readpool.getRead(id2).getLeftClipoff(),
		tryseqlen2,
		id1, id2, 1, 1, true, I->getOffsetInAlignment(id2));
	    }
	  }
	}
      }else{
	if(I->getSequenceDirection(id2)>0){
	}else{
	}
      }

      if(madsl.size()==0){
	CEBUG("No results\n");
      }else{
	int32 bestweight=0;
	auto J=madsl.begin();
	while(J!=madsl.end()){
	  if(J->isValid()==false){
	    J=madsl.erase(J);
	  }else{
	    if(J->getWeight()>bestweight) bestweight=J->getWeight();
	    ++J;
	  }
	}
	// take only the best
	for(J= madsl.begin(); J!=madsl.end();){
	  if(J->getWeight() != bestweight){
	    J=madsl.erase(J);
	  } else {
	    ++J;
	  }
	}
//	  cout << "Ext. 1st success: " << id1 << "\t" << id2 << "\n";
//	  cout << *I;
//	  cout << *(madsl.begin());

	int32 lens1=0;
	int32 lens2=0;
	if(madsl.begin()->clipper(as_params.as_readextension_window_len,
				  as_params.as_readextension_window_maxerrors,
				  lens1, lens2)){
//	    cout << "Lalala\n";

	  lens1-=AS_readpool.getRead(id1).getLenClippedSeq();
	  lens2-=AS_readpool.getRead(id2).getLenClippedSeq();
	  CEBUG("o1: " << AS_readpool.getRead(id1).getLenClippedSeq() << "\tn: " << lens1);
	  CEBUG("\no2: " << AS_readpool.getRead(id2).getLenClippedSeq() << "\tn: " << lens2<<endl);


	  if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension){
	    if(lens1>5 && lens1>clips[id1].len){
	      clips[id1].len=lens1;
	      clips[id1].changed=true;
	    }
	  }

	  if(AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension){
	    if(lens2>5 && lens2>clips[id2].len){
	      clips[id2].len=lens2;
	      clips[id2].changed=true;
	    }
	  }
	}
      }
    }
    P.finishAtOnce();
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  int32 lenplus=0;
  int32 numchanged=0;
  for(uint32 rid=0; rid<clips.size(); rid++){
    if(AS_readpool.getRead(rid).isBackbone()
       || AS_readpool.getRead(rid).isRail()) continue;
    // contig join spoiler! do not extend back again!
    if(AS_readpool.getRead(rid).hasTag(Read::REA_defaulttag_CJSP.identifier)) continue;
    if(AS_miraparams[AS_readpool.getRead(rid).getSequencingType()].getAssemblyParams().as_use_read_extension) continue;

    if(clips[rid].changed){
      CEBUG("ID: " << rid << "\t" << AS_readpool.getRead(rid).getName() << "\toldlen: " << AS_readpool.getRead(rid).getLenClippedSeq());
      CEBUG("\tgained: " << clips[rid].len << endl);
      numchanged++;
      lenplus+=clips[rid].len;

      logout << AS_readpool.getRead(rid).getName() << "\t" << clips[rid].len << "\n";

      AS_readpool.getRead(rid).setClipoffs(AS_readpool.getRead(rid).getLeftClipoff(),
					 AS_readpool.getRead(rid).getLeftClipoff()+AS_readpool.getRead(rid).getLenClippedSeq()+clips[rid].len-1,
					 false);

      if(AS_readpool.getRead(rid).checkRead()){
	cout << AS_readpool.getRead(rid);
	MIRANOTIFY(Notify::INTERNAL, AS_readpool.getRead(rid).checkRead()) ;
      }
    }
  }

  cout << "\nChanged length of " << numchanged << " sequences."<< endl;
  if(numchanged!=0){
    cout << "Mean length gained in these sequences: " << static_cast<double>(lenplus)/ static_cast<double>(numchanged) << " bases." << endl;
  }

  AS_steps[ASADSLISTOK]=0;

  FUNCEND();
  return;
}
//#define CEBUGFLAG 0




#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::analyseOverlapHashProfile(std::vector<uint8> & profile, std::vector<skimedges_t>::const_iterator seI, ADSEstimator & adse)
{
  std::vector<uint32> longeststretch(7,0);
  std::vector<uint32> currentstretch(7,0);

  for(size_t pi=0; pi<profile.size(); pi++){
    //CEBUG(pi << '\t' << static_cast<uint16>(profile[pi]) << '\n');
    for(size_t si=0; si<7; si++){
      if(si==profile[pi]){
	currentstretch[si]++;
	if(currentstretch[si]>longeststretch[si]) longeststretch[si]=currentstretch[si];
      }else{
	currentstretch[si]=0;
      }
    }
  }

  if(longeststretch[3]<5){
    if(AS_skimstaken[seI->skimindex]==true){
      cout << "Remove seI: " << *seI;
      cout << "stretches:\n";
      for(size_t si=0; si<7; si++){
	cout << si << ' ' << longeststretch[si] << endl;
      }

      AS_skimstaken[seI->skimindex]=false;
      AS_numskimoverlaps[seI->rid1]--;
      AS_numskimoverlaps[seI->linked_with]--;
    }
  }
}

#define CEBUG(bla)



/*************************************************************************
 *
 * sorter to sort Contig::templateguessinfo_t from low to high
 *  first on readgroup id
 *  then segment_placement
 *  then on template size
 *
 *************************************************************************/

inline bool Assembly__sortTemplateGuessInfo_(const Contig::templateguessinfo_t & a,
					     const Contig::templateguessinfo_t & b);
inline bool Assembly__sortTemplateGuessInfo_(const Contig::templateguessinfo_t & a,
					     const Contig::templateguessinfo_t & b)
{
  if(a.rgid==b.rgid){
    if(a.splace_seen==b.splace_seen){
      return a.tsize_seen<b.tsize_seen;
    }
    return a.splace_seen<b.splace_seen;
  }
  return a.rgid<b.rgid;
}

/*************************************************************************
 *
 * Destroys AS_templateguesses by sorting it, therefore cleared at end of func
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::analyseTemplateGuesses()
{
//  cout << "tgs: " << AS_templateguesses.size() << endl;
  if(AS_templateguesses.empty()) return;
  if(AS_templateguesses.size()<65536){
    mstd::ssort(AS_templateguesses,Assembly__sortTemplateGuessInfo_);
  }else{
    mstd::psort(AS_templateguesses,Assembly__sortTemplateGuessInfo_);
  }

  {
    uint32 ind=0;
    for(auto & tge : AS_templateguesses){
      CEBUG(ind
	    << "\trgid: " << tge.rgid.getLibId()
	    << "\tspl: " << static_cast<int16>(tge.splace_seen)
	    << "\tts: " << tge.tsize_seen
	    << '\n'; ++ind);
    }
  }


  std::vector<rgtguess_t> atgpred;

  // collapse AS_templateguesses into predictions grouped by rgid and segment_placement
  auto tgI=AS_templateguesses.cbegin();
  for(; tgI!=AS_templateguesses.cend() && tgI->rgid.isDefaultNonValidReadGroupID(); ++tgI) {}
  if(tgI!=AS_templateguesses.cend()){
    while(tgI!=AS_templateguesses.cend()){
      auto tgS=tgI;
      // search end iterator for that grouping
      for(; tgI!=AS_templateguesses.cend() && tgI->rgid==tgS->rgid && tgI->splace_seen == tgS->splace_seen; ++tgI) {};
      auto tgE=tgI;
      auto numvals=tgI-tgS;
      if(numvals>=100){
	atgpred.resize(atgpred.size()+1);
	atgpred.back().count=numvals;
	atgpred.back().tg=*tgS;

	uint64 sum=0;
	for(tgI=tgS; tgI!=tgE; ++tgI) sum+=tgI->tsize_seen;
	double mean=static_cast<double>(sum)/static_cast<double>(numvals);
	double sp2sum=0;
	for(tgI=tgS; tgI!=tgE; ++tgI) {
	  auto m1=static_cast<double>(tgI->tsize_seen - mean);
	  sp2sum+=m1*m1;    // for stdev
	}
	double stdev=sqrt(sp2sum/(numvals-1));
	// skewness, see
	//   http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
	//   http://www.tc3.edu/instruct/sbrown/stat/shape.htm
	//
	// skewness >0 left skew
	// skewness <0 right skew
	// abs skewness >= 1 very skewed
	//     skewness >= .5 mildly skewed
	//     skewness >= .25 slightly skewed (own definition as we're working with higher number of measurements
	//
	// but for large data sets with tens or hundreds of thousands of measurement, even a few outliers (<1%) on one side
	//  (like happens in misassembled contigs) can induce a skew >1
	// Solution: recalc skew on a subset of values which takes only values at 3*stdev
	// also recalc new stdev on this subset

	sp2sum=0;
	double sp3sum=0;
	double sxstdev=stdev*3;
	uint64 numvals2=0;
	for(tgI=tgS; tgI!=tgE; ++tgI) {
	  auto m1=static_cast<double>(tgI->tsize_seen - mean);
	  if(abs(m1)<=sxstdev){
	    sp2sum+=m1*m1;    // for stdev
	    sp3sum+=m1*m1*m1; // for skewness
	    ++numvals2;
	  }
	}
	stdev=sqrt(sp2sum/(numvals2-1));
	double skewness=0;
	if(stdev>0) skewness=sp3sum/((numvals2-1)*stdev*stdev*stdev);

	// homebrew
	double lefttailfactor=static_cast<double>(2.0);
	double righttailfactor=static_cast<double>(2.0);
	if(skewness<=static_cast<double>(-1.0)){
	  lefttailfactor=static_cast<double>(2.4);
	  righttailfactor=static_cast<double>(1.6);
	}else if(skewness<=static_cast<double>(-0.5)){
	  lefttailfactor=static_cast<double>(2.2);
	  righttailfactor=static_cast<double>(1.8);
	}else if(skewness<=static_cast<double>(-0.25)){
	  lefttailfactor=static_cast<double>(2.1);
	  righttailfactor=static_cast<double>(1.9);
	}else if(skewness>=static_cast<double>(1.0)){
	  lefttailfactor=static_cast<double>(1.6);
	  righttailfactor=static_cast<double>(2.4);
	}else if(skewness>=static_cast<double>(0.5)){
	  lefttailfactor=static_cast<double>(1.8);
	  righttailfactor=static_cast<double>(2.2);
	}else if(skewness>=static_cast<double>(0.25)){
	  lefttailfactor=static_cast<double>(1.9);
	  righttailfactor=static_cast<double>(2.1);
	}

	atgpred.back().count=numvals2;
	atgpred.back().mean=mean;
	atgpred.back().stdev=stdev;
	atgpred.back().skewness=skewness;
	if(stdev>0){
	  atgpred.back().deduced_min=(mean-lefttailfactor*stdev > 0) ? (mean-lefttailfactor*stdev) : 0;
	  atgpred.back().deduced_max=mean+righttailfactor*stdev;
	}else{
	  atgpred.back().deduced_min=mean-mean/10;
	  atgpred.back().deduced_max=mean+mean/10;
	}
      }
    }
  }

  cout << "ATG PREDICTIONS\n";
  cout << std::fixed << std::setprecision(10);
  for(auto & ae : atgpred){
    cout << ae << endl;
  }

  if(AS_rgstinfo.size()<=ReadGroupLib::getNumReadGroups()){
    AS_rgstinfo.resize(ReadGroupLib::getNumReadGroups()+1);
  }

  // per readgroup, search the prediction with the highest count: that's going to be our final prediction
  auto atgpI=atgpred.cbegin();
  while(atgpI!=atgpred.cend()){
    auto best=*atgpI;
    for(; atgpI!=atgpred.cend() && atgpI->tg.rgid==best.tg.rgid; ++atgpI){
      if(atgpI->count >best.count) best=*atgpI;
    }
    cout << "Final prediction: " << best << endl;
    bool changedsomething=false;
    if(best.tg.rgid.wantSegmentPlacementEstimate()){
      best.tg.rgid.setSegmentPlacementCode(best.tg.splace_seen);
      best.tg.rgid.setWantSegmentPlacementEstimate(false);
      cout << "Set segment placement code.\n";
      AS_guessedtemplatevalues=true;
      changedsomething=true;
    } else if(best.tg.rgid.wantTemplateInfoEstimate()){
      best.tg.rgid.setInsizeFrom(best.deduced_min);
      best.tg.rgid.setInsizeTo(best.deduced_max);
      best.tg.rgid.setWantTemplateSizeEstimate(false);
      cout << "Set template size.\n";
      AS_guessedtemplatevalues=true;
      changedsomething=true;
    }
    if(changedsomething){
      AS_rgstinfo[best.tg.rgid.getLibId()].rgtguess=best;
    }
  }

  AS_templateguesses.clear();
}
//#define CEBUG(bla)


bool Assembly::warnChimeraContent()
{
  FUNCSTART("bool Assembly::warnChimeraContent()");
  BUGIFTHROW(AS_current_rls_byrg.size() != ReadGroupLib::getNumReadGroups(), "AS_current_rls_byrg.size() " << AS_current_rls_byrg.size() << " != ReadGroupLib::getNumReadGroups() " << ReadGroupLib::getNumReadGroups() << " ???");

  bool retval=false;

  std::vector<double> ratios;
  ratios.resize(ReadGroupLib::getNumReadGroups());

  uint32 warnlevel=2;

  std::string maxverdict;
  double maxratio=0.0;

  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    AS_current_rls_byrg[rgi].RLS_count_chimera=0;
    if(AS_current_rls_byrg[rgi].RLS_count_all){
      for(uint32 rid=0; rid<AS_readpool.size(); ++rid){
	if(AS_debrisreason[rid]==DEBRIS_CLIP_CHIMERA) ++AS_current_rls_byrg[rgi].RLS_count_chimera;
      }

      ratios[rgi]=100.0/AS_current_rls_byrg[rgi].RLS_count_all*AS_current_rls_byrg[rgi].RLS_count_chimera;

      if(ratios[rgi] > 0.0){
	retval=true;
	if(ratios[rgi] < 0.1){
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="still excellent, no need to worry";
	}else if(ratios[rgi] < 0.25){
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="ok, but could be better";
	}else if(ratios[rgi] < 0.5){
	  warnlevel=1;
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="not good";
	}else if(ratios[rgi] < 1.0){
	  warnlevel=1;
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="mediocre";
	}else if(ratios[rgi] < 2.0){
	  warnlevel=0;
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="inacceptable";
	}else if(ratios[rgi] < 4.0){
	  warnlevel=0;
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="absolutely appalling";
	}else if(ratios[rgi] < 8.0){
	  warnlevel=0;
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="a complete catastrophe";
	}else{
	  AS_current_rls_byrg[rgi].RLS_verdict_chimera="a reason to shoot the harddisk containing this data";
	  warnlevel=0;
	}
      }
      if(ratios[rgi]>maxratio){
	maxratio=ratios[rgi];
	maxverdict=AS_current_rls_byrg[rgi].RLS_verdict_chimera;
      }
    }
  }

  if(maxratio>0.0){
    std::string wstr("MIRA detected chimeric sequences in (at least) one of your readgroups. The maximum percentage found was ");
    if(maxratio<0.005){
      wstr+="<0.005";
    }else{
      wstr+=boost::str(boost::format("%.2f") % maxratio);
    }
    wstr+=+"%, which is "+maxverdict+".";
    if(maxratio>=0.25){
      if(maxratio>=4.0){
	wstr+="\n\nYour sequencing provider absolutely needs to get lower numbers, talk to them about it.";
      }else{
	wstr+="\n\nI suggest you ask your sequencing provider about this.";
      }
    }
    if(maxratio>=0.5){
      if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	wstr+="\n\nAs this is a genome assembly, you should be able to get away with it. But should you use the same library protocols for sequencing of RNASeq/EST data, this will create problems.";
      }else{
	wstr+="\n\nUsing a library with this amount of chimeric reads is very, very dangerous in RNASeq/EST projects.";
      }
    }
    wstr+="\n\nThe reads detected as chimeric are denoted in the 'debris' file in the info directory, the code they are marked with is CLIP_CHIMERA";
    AS_warnings.setWarning("CHIMERIC_READS",warnlevel,"Readgroup with chimeric reads",wstr);
  }
  return retval;
}

bool Assembly::warnAtSmileCoverage()
{
  std::string wstr;
  if(AS_assemblyinfo.huntForSmileCoverage(wstr)){
    AS_warnings.setWarning("CONCOV_SUSPICIOUS_DISTRIBUTION",0,"Suspicious distribution of contig coverages",wstr);
  }
  return !wstr.empty();
}

bool Assembly::warnAtHighCoverages(uint32 measuredcov)
{
  bool retval=false;
  if(AS_miraparams[0].getNagAndWarnParams().nw_check_coverage!=NWNONE
     && AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
     && measuredcov>AS_miraparams[0].getNagAndWarnParams().nw_check_covvalue){
    std::string wstr("You are running a genome ");
    if(AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping){
      wstr+="mapping";
    }else{
      wstr+="de-novo";
    }
    wstr+=" assembly and the current best estimation for average coverage is "
      +boost::lexical_cast<std::string>(measuredcov)
      +"x (note that this number can be +/- 20% off the real value). This is ";
    if(measuredcov>=80) wstr+="a pretty high coverage, ";
    wstr+="higher than the current warning threshold of "
      +boost::lexical_cast<std::string>(AS_miraparams[0].getNagAndWarnParams().nw_check_covvalue)
      +"x."
      "\n\nYou should try to get the average coverage not higher than, say, 60x to 100x for Illumina data or 40x to 60x for 454 and Ion Torrent data. Hybrid assemblies should target a total coverage of 80x to 100x as upper bound. For that, please downsample your input data."
      "\n\nThis warning has two major reasons:"
      "\n- for MIRA and other overlap based assemblers, the runtime and memory requirements for ultra-high coverage projects grow exponentially, so reducing the data helps you there"
      "\n- for all assemblers, the contiguity of an assembly can also suffer if the coverage is too high, i.e. you get more contigs than you would otherwise. Causes for this effect can be non-random sequencing errors or low frequency sub-populations with SNPs which become strong enough to be mistaken for possible repeats.";
    if(measuredcov>=150){
      if(measuredcov>=300){
	wstr+="\nA coverage of >300x ... no really, are you kidding me? *sigh*";
      }
      wstr+="\nWith the coverage you currently have, you *really* should downsample your data. You. Have. Been. Warned!";
    }
    wstr+="\nOf course, you can always choose to ignore these warning by changing -NW:cac and -NW:acv.";
    AS_warnings.setWarning("ASCOV_VERY_HIGH",1,"Very high average coverage",wstr);
    retval=true;
  }
  return retval;
}



/*************************************************************************
 *
 * mean function, with side-effect and definetely to be used only where
 *  it is currently called from
 *
 * Used to guess which inserts in a contigs are wrong after a first
 *  mapping round and makes sure they get removed
 *
 * E.g.:
 *        bb    cagtcatga***ctgcatgca
 *        r1    cag*catga***ctgcatgca
 *        r2    cag*catga***ctgcatgca
 *        r3    cag*catga***ctgcatgca
 *        ...
 *        rX    cag*catgaTTTctgcatgca
 *
 * with rX being some weird read (maybe low frequency variant, or sequencing
 *  error, or ...)
 *
 * Normally, the two stage mapping would calculate an intermediate new backbone
 *  to be
 *
 *        bbi   cag*catga***ctgcatgca
 *
 * but that is obviously not really the best guess and currently leads to
 *  misalignments with reads ending in this area. E.g.:
 *
 *        bbi   cag*catga***ctgcatgca
 *        rY               actgcatgca
 *
 * instead of
 *
 *        bbi   cag*catga***ctgcatgca
 *        rY            a***ctgcatgca
 *
 * This function will fake the contig CON_counts structure to enable
 *  Contig::deleteStarOnlyColumns() to remove the probably wrong inserts, so
 *  that the new intermediate contig looks like this:
 *
 *        bbi   cag*catgactgcatgca
 *
 * Side-effects:
 *  - CON_counts of contig is probably not reflecting truth anymore
 *  - reads of the contig also get edited, but this routine
 *    is currently used just before all the reads are discarded anyway.
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::priv_removePotentiallyWrongBaseInserts(Contig & con)
{
  FUNCSTART("void Assembly::priv_removePotentiallyWrongBaseInserts(Contig & con)");

  // try to look for strains which are only not in backbone
  std::vector<std::string> straincons(ReadGroupLib::getNumOfStrains());
  std::vector<base_quality_t> dummyqual;

  for(uint32 sid=0; sid<ReadGroupLib::getNumOfStrains(); ++sid){
    bool takesid=false;
    for(uint32 rgid=1; rgid<ReadGroupLib::getNumReadGroups();++rgid){
      auto rg=ReadGroupLib::getReadGroupID(rgid);
      if(rg.getStrainID()==sid
	 && !rg.isBackbone()
	 && !rg.isRail()){
	takesid=true;
	break;
      }
    }
    if(takesid){
      CEBUG("Taking sid " << sid << " " << ReadGroupLib::getStrainOfStrainID(sid) << endl);
      con.newConsensusGet(straincons[sid],dummyqual,sid);
    }
  }

  // look if we have *any* consensus made (meaning we'd have same strain backbone AND mapped reads).
  //  If not, recalc consensi for all strains
  {
    bool needrecalc=true;
    for(auto & s : straincons){
      if(!s.empty()) {
	needrecalc=false;
	break;
      }
    }
    if(needrecalc){
      CEBUG("Recalc all consensi\n");
      for(uint32 sid=0; sid<ReadGroupLib::getNumOfStrains(); ++sid){
	con.newConsensusGet(straincons[sid],dummyqual,sid);
      }
    }
  }

  {
    bool stillnocons=true;
    for(auto & s : straincons){
      if(!s.empty()) {
	stillnocons=false;
	break;
      }
    }
    BUGIFTHROW(stillnocons,"Ummmm ... no cons built???");
  }

  // fill up empty consensi
  for(auto & s : straincons){
    if(s.empty()) s.resize(con.getContigLength(),'*');
  }

  Contig::cccontainer_t & cc = const_cast<Contig::cccontainer_t &>(con.getConsensusCounts());

  auto ccI=cc.begin();
  for(uint32 actcontigpos=0; actcontigpos<straincons[0].size(); ++actcontigpos, ++ccI){
    // obviously, we should remove a gap backbone position only if it is not a valid IUPAC base
    //  or else we would actually edit bases of the backbone away (not good)
    if(!dptools::isValidIUPACBase(ccI->i_backbonecharorig)){
      bool maydelete=true;
      CEBUG("acp: " << actcontigpos << " : ");
      for(auto & s : straincons){
	CEBUG(s[actcontigpos]<<'#');
	if(s[actcontigpos]!='*'){
	  maydelete=false;
	  //break;
	}
      }
      CEBUG("\n");
      if(maydelete){
	CEBUG("rpwbi will delete " << actcontigpos << endl);
	ccI->total_cov=65535;
	ccI->star=65535;
      }
    }
  }
  con.deleteStarOnlyColumns(0,con.getContigLength(),false,65535);
}
//#define CEBUG(bla)
