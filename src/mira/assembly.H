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

#ifndef _assembly_h_
#define _assembly_h_

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include "mira/assembly_info.H"
#include "mira/assembly_output.H"
#include "mira/dataprocessing.H"
#include "mira/manifest.H"
#include "mira/overlapedges.H"
#include "mira/parameters.H"
#include "mira/readpool.H"
#include "mira/readlotstat.H"
#include "mira/skim.H"              // only for bannedoverlappairs_t ??? see if this can be improved
#include "mira/warnings.H"

class PPathfinder;


#ifdef MIRA_HAS_EDIT
#include "EdIt/hypothesen.H"
#include "examine/scf_look.H"
#else
typedef bool EDITParameters;
#endif



enum { ASREADPOOLOK=1, ASMATCHTABLEOK, ASADSLISTOK, ASCONTIGSOK,
       ASUSEDIDSOK, ASVECTORSCLIPPED,
       ASNUMOFSTEPS};


class Assembly
{
public:
  enum { DEBRIS_NOTDEBRIS=0,
	 DEBRIS_UNSPECIFIED,
	 DEBRIS_NOOVERLAP,
	 DEBRIS_MEGAHUB,
	 DEBRIS_MASKEDNASTYREPEAT,
	 DEBRIS_MASKEDHAF7REPEAT,
	 DEBRIS_MASKEDHAF6REPEAT,
	 DEBRIS_NOTMAPPED,
	 DEBRIS_ABORTEDCONTIGCREATION,
	 DEBRIS_TINYCONTIG,
	 DEBRIS_TINYCLUSTER,
	 DEBRIS_TINYCLUSTERORPHAN,
	 DEBRIS_UNSAVEDSINGLET,
	 DEBRIS_DIGITAL_NORMALISATION,
	 DEBRIS_SHORTONLOAD,
	 DEBRIS_CLIP_BADSOLEXAEND,
	 DEBRIS_CLIP_KNOWNADAPTORRIGHT,
	 DEBRIS_CLIP_QUALMINTHRESHOLD,
	 DEBRIS_CLIP_LOWERCASEFRONT,
	 DEBRIS_CLIP_LOWERCASEBACK,
	 DEBRIS_CLIP_QUALCLIPS,
	 DEBRIS_CLIP_MASKEDBASES,
	 DEBRIS_CLIP_BADSEQUENCESERACH,
	 DEBRIS_CLIP_POLYBASEATEND,
	 DEBRIS_CLIP_POLYAT,
	 DEBRIS_CLIP_MINLEFTCLIP,
	 DEBRIS_CLIP_MINRIGHTCLIP,
	 DEBRIS_CLIP_PHIX174,
	 DEBRIS_CLIP_RRNA,
	 DEBRIS_CLIP_RRNA_PAIR,
	 DEBRIS_CLIP_PROPOSEDENDCLIP,
	 DEBRIS_CLIP_CHIMERA,
	 DEBRIS_CLIP_TERMINALLYINCORRECTIBLEORCHIMERA,
	 DEBRIS_CLIP_INCORRECTIBLEENDORCHIMERA
  };

private:

  struct usedtmpfiles_t {
    std::string basename;
    std::string filename;

    size_t estimateMemoryUsage() const
      {
	size_t components=0;
	size_t cnum,cbytes,freecap,clba;

	components+=estimateMemoryUsageOfContainer(basename,false,cnum,cbytes,freecap,clba);
	components+=estimateMemoryUsageOfContainer(filename,false,cnum,cbytes,freecap,clba);
	return components;
      }
  };

  // no dedicated save, but reconstruct from command line save / manifest?
  std::vector<MIRAParameters> AS_miraparams;
  Manifest AS_manifest;

  Warnings AS_warnings;

  AssemblyInfo AS_assemblyinfo;

  DataProcessing AS_dataprocessing;

  std::list<Contig>   AS_contigs;

  // backbonecontigs
  std::list<Contig>   AS_bbcontigs;

  void (*AS_contigbuilt_callbackfunc)(Contig &, const ReadPool &);


  /***********************************************************************************
   *
   * important for checkpointing and resuming
   *
   ***********************************************************************************/

  // if true, then tries to resume an assembly from checkpoint files.
  // set through class constructor
  bool AS_resumeasembly;

  // even more resume:
  // if true, loads existing temp data at different stages of a pass if
  // if false, does not load temp data even if it exists, recomputes
  // set to "true" when doing an assembly which "resumes", set to false as soon as
  //  in a resume assembly a data file is missing
  bool AS_resumeisok;

  /** Needs saving *******************************************************************/

  // this holds the reads
  ReadPool       AS_readpool;    // via MAF

  // which read pairs are banned from overlap
  bannedoverlappairs_t AS_permanent_overlap_bans;  // no load/save done yet

  std::vector<uint32> AS_maxcoveragereached; /* the max coverage each read has
					    attained throughout the whole
					    assembly.  Helps to refine the
					    definition of multicopies.
					 */

  /** Needs saving for resume of skim results ****************************************/
  // if a read is well connected, i.e. has long overlaps (>50% of read)
  // currently initialised with false
  // set to true if AS_overlapcritlevell & -r both < 255
  // cleared and rebuilt in each call of findPossibleOverlaps()
  //
  // this could be a std::vector<bool> if only I found a way to save that beast to disk efficiently
  std::vector<uint8> AS_wellconnected;

  /** May need saving ****************************************************************/
  // this one is slightly different to readpool.size():
  //  backbone are counted, but railreads not
  uint32     AS_num_reads_in_ass;
  uint32     AS_num_reads_valid;
  uint32     AS_num_reads_too_small;

  // reads categorised as having incorrectible end or chimera
  // currently 0 or 1;
  // set by data processing search for incorrectible ends or chimera
  std::vector<uint8>  AS_incorchim;

  // reads categorised by SKIM as megahubs
  // currently 0 or 1; maybe later also 2 for "true" megahubs (and 1 for saved megahubs);
  // vector is empty when skim should not filter for megahubs
  std::vector<uint8>  AS_skimmegahubs;

  // Why was a read sorted away
  // Never cleared, to be used as template for AS_isdebris
  std::vector<uint8>  AS_debrisreason;

  //////////// not a big problem if not saved, but results are different //////////
  // TODO: maybe consolidate into one byte, using defined enums & class for accessing?

  // has info whether a read is a troublemaker
  std::vector<uint8>  AS_istroublemaker;

  // has info whether a read once was not in contig
  //  (throughout all assembly) and hence all overlaps should be
  //  taken to ensure it has all possibilities to get back into a
  //  contig (no skim edge reduction)
  std::vector<uint8>  AS_needalloverlaps;

  std::vector<uint8>  AS_multicopies;   /* reads with more overlaps than
				      expected on average have 1 here
				      pathfinder will start building
				      elsewhere, and include those
				      last
				   */

  std::vector<uint8>  AS_hasmcoverlaps;  /* reads that overlap with a read
				       that is categorised as multi-
				       copy get 1 here
				       initialised by pathfinder if
				       vector is empty (==once
				       every pass of MIRA)
				    */
  std::vector<uint8> AS_hasreptoverlap; /* reads that have ol_rept overlaps
					   init by pathfinder if empty (once every pass)
					   maybe unused atm
					*/
  std::vector<uint8> AS_hasnoreptoverlap; /* reads that have ol_norept overlaps
					     init by pathfinder if empty (once every pass)
					  */

  std::vector<uint8> AS_maybespoilsport; /* reads that are not well connected and have no
					    valid kmers at one end (or both)
					     init by pathfinder if empty (once every pass)
					  */

  std::vector<uint8> AS_estnochimerakill; /* atm: reads that have kmer <=64 without gaps in
					     valid hash status get a true here.
					     do not clear between passes.
					     TODO: should be saved
					  */


  /***********************************************************************************/
  /***********************************************************************************/

  // these hold the filenames of the "full" skim searches
  std::string AS_posfmatch_full_filename;
  std::string AS_poscmatch_full_filename;

  // these hold the filenames of the working set of skim searches
  // these might be a full or a reduced set
  std::string AS_posfmatch_filename;
  std::string AS_poscmatch_filename;


  std::vector<skimedges_t> AS_skim_edges;  // block with x elements
  std::vector<bool> AS_skimstaken;         // size of total skims in file

  std::vector<bool>   AS_readmaytakeskim;     // size of readpool
  std::vector<uint32> AS_numskimoverlaps;     // size of readpool
  std::vector<uint32> AS_numleftextendskims;  // size of readpool
  std::vector<uint32> AS_numrightextendskims; // size of readpool

  // minimum score ratio hits must have to be taken
  std::vector<uint8> AS_skimleftextendratio; // size of readpool
  std::vector<uint8> AS_skimrightextendratio; // size of readpool

  std::ofstream AS_CUMADSLofstream;

  // adsfcontainer_t / necontainer_t defined in overlapedges.H
  adsfcontainer_t AS_adsfacts;
  necontainer_t   AS_confirmed_edges;

  uint32 AS_numADSFacts_fromalignments;
  uint32 AS_numADSFacts_fromshreds;

// temporarily used
// emptied at each Skim pass, set when cutBackPossibleChimeras() called
//  queried in ...
  std::vector<bool> AS_chimeracutflag;
  bool AS_doneskimchimera;

  // Skim will say how many hits it has written to disk per read id
  //  used for reducing skim hits afterwards
  std::vector<uint32> AS_writtenskimhitsperid;

  // Skim will also give back info what the best overlap criterion levels
  //  are for left and right overlaps
  // also used to reduce skim hits by getting rid of unnecessary matches
  // [0] is for norept overlaps, [1] for rept
  std::vector<std::vector<uint8>> AS_overlapcritlevelvl;
  std::vector<std::vector<uint8>> AS_overlapcritlevelvr;

  // Small hash statistics to searrch for phiX174
  HashStatistics<vhash64_t> AS_phix174hashstatistics;

  /* temporary use for vector clipping
     holds space for clipped! sequence, more is not needed */
  std::vector<std::vector<uint32> > AS_readhitmiss;
  std::vector<std::vector<uint32> > AS_readhmcovered;
  std::vector<uint32>          AS_count_rhm;


  std::vector<int32>           AS_clipleft;       // holds result of clipping
  std::vector<int32>           AS_clipright;      //  until can be applied


  std::vector<int8>   AS_used_ids;         /* IDs which have been already
					  put into a contig are marked
					  with 1
					 IDs wich could not have been loaded
					  or are too small are marked -1
					 IDs not used yet are marked 0 */


  uint32 AS_hashstat_avghashfreq; // rough coverage estimator; if >0 then was set by performHashAnalysis()

  std::vector<uint32> AS_coverageperseqtype; /* This stores the median of all
					   average coverages encountered in
					   contigs at the end of the current round (per
					   sequencing type) Used for contig
					   bulding in next round. Also used to
					   analyse contig coverage and refine
					   AS_multicopies
					*/
  uint32 AS_coveragetotal; // as above, but avg. total coverage of all seqtypes


  // has info whether a read is not in a contig
  //  (per pass: initialised with AS_debrisreason)
  std::vector<uint8>  AS_isdebris;


  // this vector is for the repeat resolver after the building
  //  of contigs. if true, all overlaps with the read in question will
  //  be computed
  // it's only temporarily filled
  std::vector<bool>   AS_readsforrepeatresolve;

  // the following three vectors are for computing troublemakers
  // they are only temporarily filled
  std::vector<uint32> AS_allrmbsok;
  std::vector<uint32> AS_probablermbsnotok;
  std::vector<uint32> AS_weakrmbsnotok; /* three vectors together define a ratio
				      with which a given read has problems
				      with weak rmb mismatches
				      if notok > ok, then the weights
				      of a read in the overlap graph will be
				      reduced */

  std::vector<int8>   AS_steps;    // 0 not done, 1 done, -1 need reeval


  std::vector<bool>   AS_seqtypespresent;
  bool                AS_hasbackbones;

  std::list<usedtmpfiles_t> AS_usedtmpfiles;

  ///////////////

  // stats for buildfirstcontigs;
  struct bfcstats_t {
    uint32 numnewsrm;
    uint32 numeditovercall;
    uint32 numedithashfreq;
    uint32 numdisassemblies;

    bfcstats_t() : numnewsrm(0), numeditovercall(0), numedithashfreq(0), numdisassemblies(0) {};
  };
  std::vector<bfcstats_t> AS_bfcstats;  // vector of size 2: 0 for non-rep contigs, 1 for rep

  ///////////////

  bool AS_donequickdenovocoveragecheck;

  ///////////////

  int32 AS_maxtemplateid;
  std::vector<Contig::templateguessinfo_t> AS_templateguesses;

  bool AS_guessedtemplatevalues;

  // readgroup template guess
  struct rgtguess_t {
    uint32 count=0;
    Contig::templateguessinfo_t tg;
    double mean=0;
    double stdev=0;
    double skewness=0;
    int32 deduced_min=-1;
    int32 deduced_max=-1;

    friend std::ostream & operator<<(std::ostream &ostr, const rgtguess_t & atgt){
      ostr << "rgid: " << atgt.tg.rgid.getLibId()
	   << "\tc: " << atgt.count
	   << "\tsp: " << static_cast<int16>(atgt.tg.splace_seen)
	   << " (" << ReadGroupLib::getNameOfSegmentplacement(atgt.tg.splace_seen)
	   << ")\tm: " << atgt.mean
	   << "\td: " << atgt.stdev
	   << "\ts: " << atgt.skewness
	   << "\t-: " << atgt.deduced_min
	   << "\t+: " << atgt.deduced_max;
      return ostr;
    }
  };

  // AAAAAAAAAAAAAAAAAAARRRRRRRRRRRRRRRGGGGGGGGGGGGGGGG!
  // rgstinfo_t and readlotstats are somehow duplicated functionality
  //
  // TODO: fix that
  //
  // rg statistics + template info
  struct rgstinfo_t {
    uint64 total_bases_unclipped=0;
    uint64 total_bases_clipped=0;
    uint32 numreads=0;
    uint32 numreadsclipped=0;
    rgtguess_t rgtguess;
  };

  std::vector<rgstinfo_t> AS_rgstinfo;

  std::vector<ReadLotStat> AS_current_rls_bytype;  // will be set to size ReadGroupLib::SEQTYPE_END
  std::vector<ReadLotStat> AS_current_rls_byrg;    // will be set to size ReadGroupLib::getNumReadGroups();

  // copy of the above as snapshot just before the assembly enters pass 1
  std::vector<ReadLotStat> AS_fixed_rls_bytype;
  std::vector<ReadLotStat> AS_fixed_rls_byrg;

  ///////////////

  int64 AS_systemmemory; // initialised by constructor: how much mem the system has
                          // is 0 if info not available

  uint32 AS_applydiginorminpass=0; // set dynamically


  //
  bool AS_deleteoldresultfiles;


  bool AS_shouldrun_nfs_check;

  bool AS_logflag_dumprailreads=false;
  bool AS_logflag_dumphashanalysis=false;
  bool AS_logflag_oclevel=false;
  bool AS_logflag_swbbcheck=false;
  bool AS_logflag_adsdump=false;

  bool AS_logflag_loadedoverlaps=false;

  bool AS_logflag_dumpusedids=false;


  // temp objects, so as not to create them a couple of million times
  multitag_t AS_tmptag_CRMr;

  // to be set by outside caller
  bool AS_everythingwentfine;


  //Functions

/*************************************************************************
 *
 *  assembly_misc
 *
 *************************************************************************/
private:
//  void mergeTemplateInfo(const std::string & tifile,
//			 const std::string & logname,
//			 const std::string & logprefix);
  template<class TVHASH_T>
  void priv_phahelper(const std::string & kmerfilename,
		      const std::string & signalfile,
		      uint32 basesperhash,
		      bool rarekmerfinalkill,
		      int32 version,
		      const std::string prefix,
		      const std::string postfix,
		      const std::string logname);
  void priv_performHashAnalysis(const std::string & kmerfilename,
				bool usesignal,
				bool rarekmerfinalkill,
				int32 version,
				const std::string prefix="",
				const std::string postfix="",
				const std::string logname="");
  uint64 performProposedEndClips(uint32 basesperhash,
				 const std::string & tmpname,
				 const std::string & tmpprefix);
  uint64 performPECCHIMSRERKM(const std::string & tmpname,
			      const std::string & tmpprefix);
  void cutBackPossibleChimeras(const std::string & tmpname,
			       const std::string & tmpprefix,
			       const std::vector<int32> & chuntleftcut,
			       const std::vector<int32> & chuntrightcut,
			       std::vector<bool> & chimeracutflag);

  void correctContigs();               // TODO: have a look whether here?
  void calcPossibleSeqVectorClipoffs(int32 version=-1,
				     const std::string prefix="",
				     const std::string postfix="",
				     const std::string tmpname="");
  void performSeqVectorClippings();

  void extendADS(int32 version=-1,
		 const std::string prefix="",
		 const std::string postfix="",
		 const std::string tmpname="");

  void shredReadsIntoReadPool(ReadPool & sourcepool,
			      uint32 shredlen,
			      uint32 shredoffset,
			      uint8 shredreadtype,
			      const std::string & shredstrain);
  void analyseOverlapHashProfile(std::vector<uint8> & profile,
				 std::vector<skimedges_t>::const_iterator seI,
				 ADSEstimator & adse);
  void analyseTemplateGuesses();

  bool warnAtSmileCoverage();
  bool warnChimeraContent();
  bool warnAtHighCoverages(uint32 measuredcov);

  void priv_removePotentiallyWrongBaseInserts(Contig & con);

  //void clipTo200();
  //void performHashEditing();


/*************************************************************************
 *
 *  estassembly
 *
 *************************************************************************/
public:
  std::vector<std::string> assembleESTs();


/*************************************************************************
 *
 *  assembly_reduceskimhits
 *
 *************************************************************************/
private:
  void reduceSkimHits4(int32 version=-1,
		       const std::string prefix="",
		       const std::string postfix="",
		       const std::string tmpname="");
  void rsh4_countTotalSkimsTaken();
  void rsh4_denormaliseSkimHits(const std::string & dnsfile,
				std::list<int64> & idblocks,
				std::vector<uint64> & blockpos,
				std::vector<size_t> & blocklen);
  size_t rsh4_getNextNormalisedSkimBlock(std::list<int64> & idblocks,
					 int64 & blockstartid,
					 int64 & blockendid);
  size_t rsh4_loadNormalisedSkimHitBlock(const std::string & filename,
					 uint64 skimindex,
					 int64 blockstartid, int64 blockendid,
					 int8 rid1dir, int8 rid2dir);
  void rsh4_filterSkimHits(const std::string & oldfilename,
			   const std::string & newfilename,
			   uint64 & skimindex);
  void rsh4_takeThisSkim(const skimedges_t & seI,
			 ADSEstimator & adse,
			 bool calcadse);
  void rsh4_getNextSkimBlock(const std::string & dnsfile,
			     uint32 blocki,
			     const std::vector<uint64> & blockpos,
			     const std::vector<size_t> & blocklen);
  void rsh4_purgeSkimsOfReadsCutByChimera(std::string & filename);
  void rsh4_takeRailHits(const std::string & dnsfile,
			       const std::vector<uint64> & blockpos,
			       const std::vector<size_t> & blocklen);
  void rsh4_flagMulticopyReads(const std::string & dnsfile,
			       const std::vector<uint64> & blockpos,
			       const std::vector<size_t> & blocklen);
  void rsh4_weedOutBadHits(const std::string & dnsfile,
			   const std::vector<uint64> & blockpos,
			   const std::vector<size_t> & blocklen);
  void rsh4_take100PCMappingHits(const std::string & dnsfile,
				 const std::vector<uint64> & blockpos,
				 const std::vector<size_t> & blocklen);
  void rsh4_takeNonReptHitsThatExtend(uint32 nbest,
				      std::vector<uint32> & nbestl,
				      std::vector<uint32> & nbestr,
				      uint8 minscoreratio,
				      const std::string & dnsfile,
				      const std::vector<uint64> & blockpos,
				      const std::vector<size_t> & blocklen);
  void rsh4_takeReptNonReptBorders(uint8 minscoreratio,
				   const std::string & dnsfile,
				   const std::vector<uint64> & blockpos,
				   const std::vector<size_t> & blocklen);
  void rsh4_takeReptPEPEHitsThatExtend(uint32 nbest,
				       uint8 minscoreratio,
				       const std::string & dnsfile,
				       const std::vector<uint64> & blockpos,
				       const std::vector<size_t> & blocklen);
  void rsh4_takeReptNPENPEHitsThatExtend(uint32 nbest,
					 uint8 minscoreratio,
					 const std::string & dnsfile,
					 const std::vector<uint64> & blockpos,
					 const std::vector<size_t> & blocklen);
  void rsh4_takeReptPENPEHitsThatExtend(uint32 nbest,
					uint8 minscoreratio,
					const std::string & dnsfile,
					const std::vector<uint64> & blockpos,
					const std::vector<size_t> & blocklen);
  void rsh4_takeNonTakenReadsWithHitsThatExtend(uint32 nbest,
						std::vector<uint32> & nbestl,
						std::vector<uint32> & nbestr,
						uint8 minscoreratio,
						const std::string & dnsfile,
						const std::vector<uint64> & blockpos,
						const std::vector<size_t> & blocklen);
  void rsh4_take100PC100bpHitsThatExtend(uint32 nbest,
					 std::vector<uint32> & nbestl,
					 std::vector<uint32> & nbestr,
					 const std::string & dnsfile,
					 const std::vector<uint64> & blockpos,
					 const std::vector<size_t> & blocklen);
  void rsh4_takeWellconnected(uint32 nbest,
			      std::vector<uint32> & nbestl,
			      std::vector<uint32> & nbestr,
			      const std::string & dnsfile,
			      const std::vector<uint64> & blockpos,
			      const std::vector<size_t> & blocklen);
  void rsh4_takeNonTakenSideExtends(uint32 nbest,
				    uint8 minscoreratio,
				    const std::string & dnsfile,
				    const std::vector<uint64> & blockpos,
				    const std::vector<size_t> & blocklen);
  void rsh4_takeIdenticalReadHits(uint32 nbest,
				  const std::string & dnsfile,
				  const std::vector<uint64> & blockpos,
				  const std::vector<size_t> & blocklen);
  void rsh4_takeNonTakenReadsWithHits(uint32 nbest,
				      const std::string & dnsfile,
				      const std::vector<uint64> & blockpos,
				      const std::vector<size_t> & blocklen);
  size_t rsh4_tNTSEhelper(uint32 nbest,
			  ADSEstimator & adse,
			  std::vector<std::vector<skimedges_t>::const_iterator> & sev);
  void rsh4_takeAll(const std::string & dnsfile,
		    const std::vector<uint64> & blockpos,
		    const std::vector<size_t> & blocklen);
  void rsh4_takeNeedAllOverlaps_weakgood(const std::string & dnsfile,
					 const std::vector<uint64> & blockpos,
					 const std::vector<size_t> & blocklen);
  void rsh4_takeSolexaByCritLevel(uint32 ocvi,
				  uint32 nbest,
				  std::vector<uint32> & nbestl,
				  std::vector<uint32> & nbestr,
				  const std::string & dnsfile,
				  const std::vector<uint64> & blockpos,
				  const std::vector<size_t> & blocklen);
  void rsh4_takeTemplateOverlaps(const std::string & dnsfile,
				 const std::vector<uint64> & blockpos,
				 const std::vector<size_t> & blocklen);


/*************************************************************************
 *
 *  assembly_io
 *
 *************************************************************************/
private:
  void loadSequenceData_resume();
  void loadSequenceData_new();
  void basicDataChecks();
  void basicReadGroupChecks();
  void checkForReadNameLength(uint32 stoplength, bool stop);
  //void loadExternalSCFQualities();

  void loadSequenceDataFromManifest();
  static void streamSequenceDataFromManifest(std::vector<MIRAParameters> & miraparams,
					     Manifest & man,
					     ReadPool & readpool,
					     std::list<Contig> * contigsptr,
					     std::vector<readid_t> * readsasbackbonecontigs
    );

  void postLoadBackbone();
  void clipsAfterLoad();
  void addRailsToBackbones();

  void loadStrainData(const std::string & sdfile);

  void postLoad();
  void dumpSomeStatistics();

  void saveExtTmpContig(Contig & con, std::string basename);



public:
  void dumpContigs();
  void loadSequenceData();
  void dumpRailReads(std::ofstream & fout);

  std::string buildDefaultCheckpointFileName(const std::string & filename);
  std::string buildDefaultInfoFileName(int32 version,
				       const std::string & prefix,
				       const std::string & postfix,
				       const std::string & basename,
				       const std::string & defaultname,
				       const std::string & defaultextension,
				       bool removeold=false);
  std::string buildDefaultResultsFileName(int32 version,
					  const std::string & prefix,
					  const std::string & postfix,
					  const std::string & basename,
					  const std::string & defaultname,
					  const std::string & defaultextension,
					  bool removeold=false);
  std::string getContigReadListFilename(int32 version=-1,
					const std::string & prefix="",
					const std::string & postfix="",
					const std::string & crlname="");
  void saveContigReadList(int32 version=-1,
			  const std::string & prefix="",
			  const std::string & postfix="",
			  const std::string & crlname="",
			  bool deleteoldfile=false);
  std::string getStatisticsFilename(int32 version=-1,
				    const std::string & prefix="",
				    const std::string & postfix="",
				    const std::string & statname="");
  void saveStatistics(int32 version=-1,
		      const std::string & prefix="",
		      const std::string & postfix="",
		      const std::string & statname="",
		      bool deleteoldfile=false);
  std::string getAssemblyInfoFilename(int32 version=-1,
				      const std::string & prefix="",
				      const std::string & postfix="",
				      const std::string & statname="");
  void saveAssemblyInfo(int32 version=-1,
			const std::string & prefix="",
			const std::string & postfix="",
			const std::string & statname="",
			bool deleteoldfile=false);
  std::string getLargeContigsInfoFilename(int32 version=-1,
					  const std::string & prefix="",
					  const std::string & postfix="",
					  const std::string & statname="");
  void saveLargeContigsInfo(int32 version=-1,
			    const std::string & prefix="",
			    const std::string & postfix="",
			    const std::string & statname="",
			    bool deleteoldfile=false);
  std::string getRGSTInfoFilename(int32 version=-1,
				  const std::string & prefix="",
				  const std::string & postfix="",
				  const std::string & statname="");
  void saveRGSTInfo(int32 version=-1,
		    const std::string & prefix="",
		    const std::string & postfix="",
		    const std::string & statname="",
		    bool deleteoldfile=false);
  void saveDebrisList(int32 version=-1,
		      const std::string & prefix="",
		      const std::string & postfix="",
		      const std::string & debrisname="");
  std::string getReadTagListFilename(int32 version=-1,
				     const std::string & prefix="",
				     const std::string & postfix="",
				     const std::string & rtlname="");
  void saveReadTagList(int32 version=-1,
		       const std::string & prefix="",
		       const std::string & postfix="",
		       const std::string & rtlname="",
		       bool deleteoldfile=false);
  std::string getConsensusTagListFilename(int32 version=-1,
					  const std::string & prefix="",
					  const std::string & postfix="",
					  const std::string & ctlname="");
  void saveConsensusTagList(int32 version=-1,
			    const std::string & prefix="",
			    const std::string & postfix="",
			    const std::string & ctlname="",
			    bool deleteoldfile=false);
  void saveSNPList(int32 version=-1,
		   const std::string & prefix="",
		   const std::string & postfix="",
		   const std::string & saname="",
		   bool deleteoldfile=false);
  void saveFeatureAnalysis(int32 version=-1,
			   const std::string & prefix="",
			   const std::string & postfix="",
			   const std::string & faname="",
			   const std::string & fsname="",
			   const std::string & fcname="",
			   bool deleteoldfile=false);
  std::string getFASTAFilename(int32 version=-1,
			       const std::string & prefix="",
			       const std::string & postfix="",
			       const std::string & fastaname="");
  std::string getFASTAPaddedFilename(int32 version=-1,
				     const std::string & prefix="",
				     const std::string & postfix="",
				     const std::string & fastaname="");
  void saveAsFASTA(int32 version=-1,
		   const std::string & prefix="",
		   const std::string & postfix="",
		   const std::string & fastaname="",
		   bool deleteoldfile=false);
  void saveStrainsAsFASTAQUAL(int32 version=-1,
			      const std::string & prefix="",
			      const std::string & postfix="",
			      const std::string & fastaname="",
			      bool deleteoldfile=false);
  std::string getTCSFilename(int32 version=-1,
			     const std::string & prefix="",
			     const std::string & postfix="",
			     const std::string & txtname="");
  void saveAsTCS(int32 version=-1,
		 const std::string & prefix="",
		 const std::string & postfix="",
		 const std::string & tcsname="",
		 bool deleteoldfile=false);
  std::string getCAFFilename(int32 version=-1,
			     const std::string & prefix="",
			     const std::string & postfix="",
			     const std::string & cafname="");
  void saveAsCAF(int32 version=-1,
		 const std::string & prefix="",
		 const std::string & postfix="",
		 const std::string & cafname="",
		 bool deleteoldfile=false);
  std::string getMAFFilename(int32 version=-1,
			     const std::string & prefix="",
			     const std::string & postfix="",
			     const std::string & cafname="");
  void saveAsMAF(int32 version=-1,
		 const std::string & prefix="",
		 const std::string & postfix="",
		 const std::string & cafname="",
		 bool deleteoldfile=false);
  std::string getTXTFilename(int32 version=-1,
			     const std::string & prefix="",
			     const std::string & postfix="",
			     const std::string & txtname="");
  void saveAsTXT(int32 version=-1,
		 const std::string & prefix="",
		 const std::string & postfix="",
		 const std::string & txtname="",
		 bool deleteoldfile=false);
  std::string getACEFilename(int32 version=-1,
			     const std::string & prefix="",
			     const std::string & postfix="",
			     const std::string & txtname="");
  void saveAsACE(int32 version=-1,
		 const std::string & prefix="",
		 const std::string & postfix="",
		 const std::string & acename="",
		 bool deleteoldfile=false);
  std::string getWiggleFilename(int32 version=-1,
				const std::string & prefix="",
				const std::string & postfix="",
				const std::string & txtname="");
  void saveAsWiggle(int32 version=-1,
		    const std::string & prefix="",
		    const std::string & postfix="",
		    const std::string & acename="",
		    bool deleteoldfile=false);
  std::string getGAP4DAFilename(int32 version=-1,
				const std::string & prefix="",
				const std::string & postfix="",
				const std::string & dirname="");
  void saveAsGAP4DA(int32 version=-1,
		    const std::string & prefix="",
		    const std::string & postfix="",
		    const std::string & dirname="",
		    bool deleteoldfile=false);
  std::string getHTMLFilename(int32 version=-1,
			      const std::string & prefix="",
			      const std::string & postfix="",
			      const std::string & txtname="");
  void saveAsHTML(int32 version=-1,
		  const std::string & prefix="",
		  const std::string & postfix="",
		  const std::string & htmlname="",
		  bool deleteoldfile=false);


/*************************************************************************
 *
 *  assembly_repeatresolve
 *
 *************************************************************************/
private:
  static bool ma_takeMCandOverlapWMC(Assembly & as, int32 rid1, int32 rid2);
  bool buildRepeatContigs(const int32 passnr);
public:
  void prework();


/*************************************************************************
 *
 *  assembly_swalign
 *
 *************************************************************************/
private:
  void setupAlignCache(std::vector<Align> & aligncache);
  void makeAlignmentsFromPosMatchFile(const std::string & filename,
				      const int32 version,
				      const int8 direction,
				      const bool trans100percent,
				      bool (* checkfunction)(Assembly&,int32,int32),
				      std::ofstream & matchfout,
				      std::ofstream & rejectfout);

  void computeSWAlign(std::list<AlignedDualSeq> & madsl,
		      uint32 rid1,
		      uint32 rid2,
		      int32 eoffset,
		      int8 direction,
		      std::vector<Align> & chkalign,
		      int32 hintbandwidth);

  static bool ma_takeall(Assembly & as, int32 rid1, int32 rid2);
  static bool ma_needRRFlag(Assembly & as, int32 rid1, int32 rid2);
  static bool ma_needRRFlagAndBothCRMr(Assembly & as, int32 rid1, int32 rid2);
  static bool ma_needSRMrOrTwoCRMr(Assembly & as, int32 rid1, int32 rid2);

  void makeAlignments(bool (* checkfunction)(Assembly&,int32,int32),
		      bool takefullskimfilenames,
		      const bool trans100percent,
		      int32 version=-1,
		      const std::string prefix="",
		      const std::string postfix="",
		      const std::string tmpname="");
  void priv_loadAlignmentsFromFile(int32 version=-1,
				   const std::string prefix="",
				   const std::string postfix="",
				   const std::string tmpname="");
  void priv_laffhelper(const std::string & fn,
		       ProgressIndicator<int64> & P,
		       uint64 & runningADSFactnumber,
		       uint64 & numbannedADSFacts);

  uint32 getOverlapMalusDivider(int32 id1, int32 id2);

  void clusterUnassembledReads(std::vector<int32> & clusteridperread,
			       std::vector<std::list<int32> > & readinclusterlist,
			       const std::vector<int8> & usedids);

  void minimiseMADSL(std::list<AlignedDualSeq> & madsl);
  void cleanupMADSL(std::list<AlignedDualSeq> & madsl,
		    uint32 i,
		    uint32 j,
		    int8 direction,
		    bool flag_stronggood,
		    bool flag_weakgood,
		    bool flag_belowavgfreq,
		    bool flag_norept,
		    bool flag_rept);

  int32 checkADSForRepeatMismatches(AlignedDualSeq & ads);
  int32 checkADSForRepeatMismatches_wrapped(AlignedDualSeq & ads, bool & need2ndpass);
  void transcribeHits(AlignedDualSeq & anads);

  void recalcNonPerfectSkimMappingsBySW(int32 version=-1);
  void rnpskmbs_helper(const std::string & filename, const int32 version, const int8 direction);


/*************************************************************************
 *
 *  assembly_pbcorrect
 *
 *************************************************************************/

public:
  uint32 correctPBMain();
private:
  struct pbcounts_t {
    struct insertafter_t {
      std::string what;
      uint32 count;
      insertafter_t() : count(0) {};
      insertafter_t(const std::string & s) : what(s), count(1) {};
    };

    std::vector<insertafter_t> cxia; // after this base, one might want to insert one of the vector what-strings
                                     //  (which preferably contains only one element)

    uint32 cxa;       // count how many times that base should rather be exchanged against an 'A'
    uint32 cxc;       //  ... 'C'
    uint32 cxg;       //  ... 'G'
    uint32 cxt;       //  ... 'T'
    uint32 cxgap;     //  ... a gap

    uint32 cqsafe;    // count how often base was in a safe align (alignment region flanked by k-error-free sub-aligns

    pbcounts_t() : cxa(0),cxc(0),cxg(0),cxt(0),cxgap(0),cqsafe(0) {};

    friend std::ostream & operator<<(std::ostream &ostr, const pbcounts_t & p){
      ostr << p.cxa
	   << '\t' << p.cxc
	   << '\t' << p.cxg
	   << '\t' << p.cxt
	   << '\t' << p.cxgap
	   << '\t' << p.cqsafe
	   << '\t' << p.cxia.size();
      for(auto & iae : p.cxia){
	ostr << '\t' << iae.what << '(' << iae.count << ')';
      }
      return ostr;
    }
  };

  struct pbc_timing_t {
    suseconds_t hashanalysis;
    suseconds_t killnonhaf;
    suseconds_t possibleoverlaps;
    suseconds_t trypbcorrect;
    suseconds_t tpbc_setup;
    suseconds_t tpbc_sw;

    pbc_timing_t() : hashanalysis(0), killnonhaf(0), possibleoverlaps(0),trypbcorrect(0), tpbc_setup(0), tpbc_sw(0) {};

    friend std::ostream & operator<<(std::ostream &ostr, const pbc_timing_t & t){
      ostr << "Timing PBC hasha: " << t.hashanalysis
	   << "\nTiming PBC killnonhaf: " << t.killnonhaf
	   << "\nTiming PBC possibleoverlaps: " << t.possibleoverlaps
	   << "\nTiming PBC trypbcorrect: " << t.trypbcorrect
	   << "\nTiming PBC tpbc setup: " << t.tpbc_setup
	   << "\nTiming PBC tpbc SW: " << t.tpbc_sw
	   << std::endl;
      return ostr;
    }
  };

  // sequence to sequence aligninfo
  // this stores for every alignment of a PB to a PB sequence the values used / computed
  //  by Dynamic && Align for effectively used bandwidth and the minimum distance the align
  //  was from the bands
  // Subsequent alignment of the same sequences can use this information to optimise speed
  //
  // The key of the unordered_map is an uint64 which is made by concatenation of the two read
  // ids (which should stay int32 for the foreseeable future)

  struct s2saligninfo_t {
    int32 bandwidthused;
    int32 minbanddistance;

    bool alreadyaligned;

    s2saligninfo_t(int32 bwu, int32 mbd) : bandwidthused(bwu), minbanddistance(mbd), alreadyaligned(false) {};
  };
  std::unordered_map<uint64_t,s2saligninfo_t>  AS_s2saligninfo;

  // vectors containing proposed clips (to be applied after a PB correct pass)
  std::vector<int32> AS_naclipl;
  std::vector<int32> AS_naclipr;



  uint32 correctPBMainLight(uint32 startpass);
  uint32 correctPBMainNormal(uint32 startpass);

  uint32 cpbm_passHelper(uint32 actpass,
			 float ratioacceptvalue,
			 int32 additionalbelieveborder,
			 uint8 minbk,
			 bool tweakbkmar,
			 bool killnonhafstretches,
			 bool generatenaclips,
			 bool generaterleedits,
			 bool generatebaseedits,
			 pbc_timing_t & timing);

  uint32 tryPBCorrect(uint32 actpass,
		      float ratioacceptvalue,
		      int32 additionalbelieveborder,
		      uint8 minbk,
		      bool tweakbkmar,
		      bool generatenaclips,
		      bool generaterleedits,
		      bool generatebaseedits,
		      pbc_timing_t & timing);
  void tpbc_checkRead(Read & actread);

  void tpbc_fillCorrector(readid_t actrid,
			  readid_t otherrid,
			  AlignedDualSeq & ads,
			  std::vector<pbcounts_t> & correctorcounts,
			  std::vector<uint8> & bkmar,
			  std::vector<uint8> & bkmor);
  void tpbc_fillRLECorrector(readid_t actrid,
			     readid_t otherrid,
			     AlignedDualSeq & ads,
			     std::vector<std::vector<uint32> > & correctorcounts,
			     std::vector<uint8> & bkmar,
			     std::vector<uint8> & bkmor);


  void   tpbc_fc_makeBelieveKMERMap(readid_t arid,
				    readid_t orid,
				    AlignedDualSeq & ads,
				    std::vector<uint8> & bkmm,
				    int32 kmersizeused,
				    int32 additionalbelieveborder);
  void tpbc_fc_minimumBelieveKMerMap(std::vector<uint8> & bkmm, uint8 minbk);
  void tpbc_tweakBKMAR(readid_t actrid,
		       readid_t otherrid,
		       AlignedDualSeq & ads,
		       std::vector<uint8> & bkmar,
		       std::vector<uint8> & bkmor);

  uint32 tpbc_generateBaseEdits(uint32 actpass,
				readid_t actrid,
				std::vector<pbcounts_t> & correctorcounts);
  uint32 tpbc_generateRLEEdits(uint32 actpass,
			       readid_t actrid,
			       std::vector<std::vector<uint32>> & rlecc);
  void tpbc_generateNonAlignClips(readid_t actrid,
			     std::vector<pbcounts_t> & correctorcounts);
  void killNonHAFCoveredStretchesInReads();
  void colourReadsByRLE();

public:
  static void findDegeneratePolymerases(ReadPool & rp);

/*************************************************************************
 *
 *  assembly_buildcontigs
 *
 *************************************************************************/

private:

  bool buildFirstContigs(const int32 loopnr,
			 const EDITParameters & eparams,
			 const bool lastloop);
  void saveDebrisAsContigs(uint32 & numcontigs,
			   const int32 passnr,
			   const bool lastpass);
  void bfc_sanityCheckASUSEDIDS(uint32 numexpected, uint32 numcontigs);
  void bfc_callPathfinder(const int32 passnr,
			  const uint32 iter,
			  uint32 & trackingunused,
			  bool shouldmovesmallclusterstodebris,
			  Contig & buildcon,
			  PPathfinder & qaf);
  void bfc_cp_mapWithSolexa(Contig & buildcon, PPathfinder & qaf);
  uint32 bfc_moveSmallClustersToDebris();
  bool bfc_checkIfContigMeetsRequirements(Contig & con);
  void bfc_markRepReads(Contig & con);
  bool bfc_trimDenovoIfNecessary(Contig & buildcon,
				 bool foundSRMs,
				 const std::string & basename_forextsave,
				 uint32 & trackingunused);
  void bfc_storeContig(Contig & con,
		       uint32 & numcontigs,
		       const bool mustmarkrepeats,
		       const int32 passnr,
		       const bool lastpass,
		       const bool verbose);

  void bfc_savePersistentSmallOverlaps(Contig & con, const int32 passnr, std::ostream & fout);

  void priv_tmpcheckroutine(Contig & bc);


/*************************************************************************
 *
 *  assembly
 *
 *************************************************************************/

private:

  void foolCompiler();
  void init();
  void zeroVars();


  void ensureStandardDirectories(bool purge);

//  void markNastyReapeatsWithSkim(int32 version=-1,
//			    const std::string prefix="",
//			    const std::string postfix="",
//			    const std::string tmpname="");
  void findPossibleOverlaps(int32 version=-1,
			    const std::string prefix="",
			    const std::string postfix="",
			    const std::string tmpname="");
  void fpo_main(int32 version,
		const std::string prefix,
		const std::string postfix,
		const std::string tmpname);
  void fpo_buildFileNames(int32 version,
			  const std::string prefix,
			  const std::string postfix,
			  const std::string tmpname);
  void flagMulticopyReads(const std::vector<uint32> & overlapcounter,
			  const std::string & tmpfilename);



  void performSnapshot(uint32 actpass);
  void ssdReadPool(const std::string & filename);
  void ssdPassInfo(const std::string & filename, uint32 actpass);
  void ssdMaxCovReached(const std::string & filename);
  void ssdBannedOverlaps(const std::string & filename);

  void loadSnapshotData(uint32 & actpass);
  uint32 lsdPassInfo(const std::string & filename);
  void lsdMaxCovReached(const std::string & filename);
  void lsdBannedOverlaps(const std::string & filename);

  void loadResumeDataFPO(int32 version, const std::string & prefix, const std::string & postfix);
  void saveResumeDataFPO(int32 version, const std::string & prefix, const std::string & postfix);
  void loadResumeDataMA(int32 version, const std::string & prefix, const std::string & postfix);
  void saveResumeDataMA(int32 version, const std::string & prefix, const std::string & postfix);

  void transferContigReadsToReadpool(const Contig & buildcon,
				     std::list<Contig::pbdse_t> & pbdsev,
				     int32 passnr);
  void transferContigReadTagsToReadpool(const Contig & con,
					const std::list<Contig>::const_iterator bbContigI);

  void makeNewReadPoolFromContigs();

  bool checkPossibleMatch(const possible_overlaps_t & poverl,
			  const uint32 i,
			  const uint32 j) const;

  void search454Overcalls(AlignedDualSeq & anads);
  void dump454OvercallsArrays();
  uint32 perform454OvercallEdits(bool qualonly);

  void huntSpoilSports(Contig & chkcon);

  bool isReadGrinch(int32 thereadid,
		    int32 incontigid,
		    std::vector<int32> & readstocheck,
		    std::vector<int32> & contigidofendreads);
  void transferContigRMBTagsToPermanentOverlapBans(Contig & con);
  void banReadPairGroups(const std::vector<int32> & g1, const std::vector<int32> & g2);


  void editSimple454Overcalls(const uint32 pass);

  uint32 priv_calcBasesPerHashOnPass(uint32 passnr) const;
  void   priv_hackMergeTwoResultMAFs();

  //-----------------------------------------------------------------


  void dmi_dumpALine(std::ostream & ostr, const char * desc,
		     size_t numelem,
		     size_t bytes_size,
		     size_t freecapacity,
		     size_t lostbyalign);



  void preassembleTasks(bool usereadextension, bool clipvectorleftovers);
  void priv_setupAutoPasses();

  bool checkTerminationRequest();

public:
  Assembly(Manifest & manifest, std::vector<MIRAParameters> & params, bool resumeassembly);
  ~Assembly();

  Assembly(Assembly const &other) = delete;
  Assembly const & operator=(Assembly const & other) = delete;

  void discard();
  void dumpMemInfo();

  inline void setContigBuiltCallback(void (*cbfunc)(Contig &, const ReadPool &) = nullptr){
    AS_contigbuilt_callbackfunc=cbfunc;
  }
  void assemble();
  void setEverythingWentFine(bool b) { AS_everythingwentfine=b; };

  static bool markRepeats(Contig & con,
			  std::vector<bool> & readsmarkedsrm,
			  Contig::repeatmarker_stats_t & repstats);

  static void refreshContigAndReadpoolValuesAfterLoading(
    ReadPool & rp,
    std::list<Contig> & contigs,
    uint8 nagwarn_templateproblems
    );

  const Read & getRead(uint32 index);
  std::list<Contig> & getContigs() {return AS_contigs;};
  ReadPool & getReadPool() {return AS_readpool;};

  void checkForNFSMountOnTmpDir();

  uint32 cleanupOldFile(const std::string & basename, const std::string & filename);
  std::string buildFileName(int32 version=-1,
			    const std::string & prefix="",
			    const std::string & postfix="",
			    const std::string & basename="",
			    const std::string & suffix="",
			    const std::string & dirname="",
			    bool removeold=true);


  void test();

  void setExtendedLog(bool f) {
    AS_logflag_dumprailreads=f;
    AS_logflag_dumphashanalysis=f;
    AS_logflag_oclevel=f;
    AS_logflag_swbbcheck=f;
    AS_logflag_adsdump=f;
    AS_logflag_dumpusedids=f;
    AS_logflag_loadedoverlaps=f;
  }

  void sortReadPool(std::vector<uint32> & dummy);
};


#endif
