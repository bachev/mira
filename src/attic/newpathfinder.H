/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2002 and later by Bastien Chevreux
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

#ifndef _bas_pathfinder_h_
#define _bas_pathfinder_h_


#include "stdinc/defines.H"

#include <iostream>

#include "mira/assembly.H"
#include "mira/readpool.H"
#include "mira/contig.H"



struct nextreadtoadd_t{
  int32 refid;
  int32 newid;                   // number of read in seqpool
  uint32 weight;
  int16  direction_newid;
  const AlignedDualSeqFacts * ads_node;

  uint8 foundmethod;

  friend std::ostream & operator<<(std::ostream &ostr, nextreadtoadd_t const &nrta);
};


class Pathfinder
{
public:
  enum {FOUND_NORMAL_AND_SRMB,
	FOUND_NORMAL,
	FOUND_SPEEDYMAP,
	FOUND_NONMULTICOPY,

	FOUND_TEMPLATE_AND_TPARTNER,
	//FOUND_TEMPLATE_AND_OVERLAPPINGTPARTNER,
	FOUND_TPARTNER_IN_CONTIG,
	FOUND_OVERLAPPING_PARTNER,

	FOUND_TEMPLATE_NONMULTICOPY_AND_TPARTNER,
	FOUND_RAIL_AND_NONMULTICOPY,
	FOUND_RAIL_AND_MULTICOPY,
	FOUND_ENUMEND
  };

  enum {FOUND_STAGE1=100,
	FOUND_STAGE2,
	FOUND_STAGE100
  };

private:
  //Structures

  struct beststartinfo_t {
    uint32 bsi_clustersize;
    uint32 bsi_rid;
    // don't need to add that one, just for debugging
    uint32 bsi_numconnects;

    friend std::ostream & operator<<(std::ostream &ostr, beststartinfo_t const &bsi){
      ostr << "BSI: rid " << bsi.bsi_rid
	   << "\tcs:" << bsi.bsi_clustersize
	   << "\tnc: " << bsi.bsi_numconnects << std::endl;
      return ostr;
    }
  };

  //Variables
private:
  uint32 PAF_valid;

  std::vector<MIRAParameters> * PAF_miraparams;
  pathfinder_parameters PAF_pafparams;


  ReadPool * PAF_readpool;
  std::vector<newedges_t>   * PAF_forward_edges;

  std::vector<int8>  * PAF_used_ids_ptr;
  std::vector<int32> * PAF_ids_in_contig_ptr;
  std::vector<uint8> * PAF_istroublemaker_ptr;

  // size of number of reads in assembly
  // is cleared on every constructStepByStep()
  std::vector<int32>   PAF_fresh_ids_in_contig;
  std::vector<int32>   PAF_cambium_ids;

  // temporary object, used to shuffle and compress cambiumids
  // as class variable to reduce memory allocs
  std::vector<int32>   PAF_tmpcambiumids;

  // size of number of reads in assembly
  // must contain all 0 at begin of constructStepByStep()
  std::vector<int8>    PAF_used_ids_in_this_run;
  bool            PAF_uiitr_is_clean; /* tells whether
				       * PAF_used_ids_in_this_run is clean
				       * (i.e. all 0). If not, the setup
				       * routine should initialise it
				       *
				       * used to save cleaning time for big
				       * assemblies
				       */

  std::vector<AlignedDualSeqFacts> * PAF_adsfacts;

  /* in this run read has a counter. Each alignment try decreases
     by 1, when 0 is reached ->  blacklisted
     blacklisting is erased in-between redo-loops

     Also: if a read in the contig has no more potential un-assembled
     partner-reads -> blacklisted

   */
  std::vector<uint8>   PAF_blacklisted_ids;

  // not used yet
  // std::vector<uint8>   PAF_hasSRMBtags;     /* cache for quick lookup,
  // 					  initialised the first time
  // 					  the pathfinder object is
  // 					  started for building a new contig
  // 				       */
  std::vector<uint8> * PAF_multicopies_ptr; /* reads with more overlaps than
					  expected on average have 1
					  provided by assembly class
					  pathfinder will start building
					  elsewhere, and include those
					  last
				       */

  std::vector<uint8> * PAF_wellconnected_ptr; /* vector of bool values
					    reads having a level 0 overlap
					    criterion with left extend and
					    level 0 for right extend are
					    well connected
					    (see skim and assembly for level
					    settings)
					 */

  std::vector<uint8> * PAF_hasmcoverlap_ptr; /* reads that overlap with a read
					   that is categorised as multi-
					   copy get 1 here
					   initialised by pathfinder if
					   vector is empty (==once
					   every pass of MIRA)
					*/

  /* this one is a hack to fasten up different functions
     lower_bound() is called very often, so caching all
     possible results (= number of reads in readpool)
     speeds the things up
     New from 2.7.4 on: moved outside class as needs "a lot" of time in
     constructor.
     Moved to call of constructStepByStep() so that it can be re-used
  */
  std::vector<std::vector<newedges_t>::iterator> * PAF_lowerbound_oedges_ptr;


  /*
     When in backbone assembly, this vector is initialised with
     all read ids that have an overlap with a rail
     The initial mapping is the to simply work through this vector
     and try mapping each read in there once
     Once the vector is empty, switches to normal backbone assembly.
   */
  std::vector<uint32> PAF_railoverlapcache;

  /*
    When in backbone assembly, this records the rails in the contig when starting
    the pathfinder
    Used by prepareRailOverlapCache() to save time in projects with many reference /
    backbone sequences
  */
  std::vector<uint32> PAF_railsincontig;


  /*
    For each contig, newpathfinder searches the whole remaining
    read space for the then best start site, i.e., having the best chance
    of having a really long, good and easy assembly stretch.

    For assemblies with millions of reads, this gets slow if the assembly
    is fragmented (i.e. 20 million reads in 20k contigs) as each search
    may take up to a second.

    If MIRA is told to use caching of all best chosen start sites (done in
    n4_searchBestStrongGoodStartEnvironment*(), this is the structure it uses.
   */
  std::vector<beststartinfo_t> PAF_beststartcache;;

  /* set by n4_searchBestStrongGoodStartEnvironment_sub()
     says what the startcache contains
     0 == best quality, strong good overlaps found, well connected reads
     1 == medium quality, overlaps but not strong good, not multicopy, well connected reads
     2 == medium quality, overlaps but not strong good. multicopy, well connected reads
     3 == not well connected
     4 == singlets
  */
  uint8 PAF_bsccontent;
  bool PAF_usesbeststartcache; // set to true if caching routines were used

  // counter

  //  when building in public quiet mode, counts how many
  //  success/failure characters have been printed out
  //  if a given limit is reached, a newline has to be printed
  // avoids overlong lines in the logfile
  uint16 PAF_buildcontig_newlinecounter;

  int16 PAF_skipwholecontigscan_counter;

  // count how many of the newedges_t elements in the vector are banned
  //  in the constructStepByStep() run
  // this allows to save time when starting the next contig
  //  if no overlaps were banned, the banned flag in the edges does
  //  not need to be reset
  // on de-novo 6m Solexa paired-end 36mers, this saves ~60 to 70% time
  //  in the contig building phase (and 45-50% overall)
  uint32 PAF_overlapsbanned;

  // small store is a store with a certain capacity (50k?) that takes
  //  up iterators to all banned overlaps until full.
  // So when it comes to reset pf_banned flag, if store is not full only
  //  the overlap edges whose iterators are in this store have to be cleared
  //  instead of iterating through the complete overlap edges. The latter
  //  is needed if store is full.
  // on de-novo 6m Solexa paired-end 36mers, this reduces another
  //  15% in the contig building phase (38 minutes to 33 minutes)
  // tradeoff of 5 minutes for 400k RAM is pretty good (160k RAM would
  //  suffice in that example, but let's plan for a bit more)
  // influence of mechanism grows with number of overlaps, so
  //  assemblies with more reads will profit exponentially
  std::vector<std::vector<newedges_t>::iterator> PAF_overlapsbanned_smallstore;


  // counts how many read add attempts there were (contig.addRead())
  //  in this constructStepByStep() (needed for Solexa mapping assembly,
  //  where loops can stop if no more adding attempt exists)
  uint32 PAF_readaddattempts;

  // counts how many non-multicopy reads have been added since the
  //  last reject from the contig
  uint32 PAF_nonmulticopiesadded;

  // flag: whether this contig may add multicopy/multicopy overlaps
  bool PAF_musthonourmcmcflag;
  bool PAF_mayaddmcmcoverlaps;


public:

  //Functions
private:
  void foolCompiler();

  uint32 searchBestEnvironment();
  void buildContig(std::vector<Align> & aligncache,int32 startid, Contig & con);

  void checkindex(int32);

  void findOverlappingTemplatePartnerID(const std::vector<int32> & readsincontig,
					nextreadtoadd_t & resultread);
  void findNextOverlap(const std::vector<int32> & readsincontig,
		       nextreadtoadd_t & resultread);
  void prepareRailOverlapCache(const std::vector<bool> & allowedrefids);
  void findNextBackboneOverlapNormal(const std::vector<int32> & readsincontig,
				     nextreadtoadd_t & resultread,
				     const std::vector<bool> & allowedrefids,
				     bool allowmulticopies,
				     bool allowtroublemakers);
  void findNextBackboneOverlapQuick(nextreadtoadd_t & resultread,
				    const std::vector<bool> & allowedrefids,
				    bool allowmulticopies,
				    bool allowtroublemakers,
				    bool allowsmallhits);


  void   n3_basicCSBSSetup();
  void   n3_handleReadNotAligned(Contig::errorstatus_t & contigerrstat,
				 nextreadtoadd_t const &nrta);

  bool n4_checkQuickRules(const AlignedDualSeqFacts * ads_node, uint32 newreadid);

  void   n4_basicCSBSSetup();
  void   n4_handleReadNotAligned(Contig::errorstatus_t & contigerrstat,
				 nextreadtoadd_t const &nrta);
  int32  n4_searchBestStrongGoodStartEnvironmentEST();
  int32  n4_searchBestStrongGoodStartEnvironmentGenome();
  int32  n4_searchBestStrongGoodStartEnvironment_sub();
  void   n4_searchBestStrongGoodStartEnvironment_subFillCache(bool wanttroublemakercheck,
							      bool wantstronggoodcheck,
							      bool wantmulticopycheck,
							      bool wantwellconnectedcheck);
  static bool sortbeststartinfo_t_(const beststartinfo_t & a,
				   const beststartinfo_t & b);

  void   n4_buildContigStartWithStrongGood(
    std::vector<Align> & aligncache,
    int32 startid,
    Contig & con);
  bool   n4_findNextGoodOverlap(
    bool wantfreqcheck,
    bool allowbelowavgfreq,
    bool allowspeedypair,
    std::vector<int32> & cambiumidsincontig,
    nextreadtoadd_t & resultread,
    clock_t & maxallowedclocks);
  void   n4_findNextPillarCounterpart(
    std::vector<int32> & cambiumidsincontig,
    nextreadtoadd_t & resultread,
    clock_t & maxallowedclocks);
  void   n4_setupAllowedReferences(
    std::vector<int32> & allowed_references,
    std::vector<int32> & ids_in_contig);
  void   n4_findNextAllowedRefOverlap(
    std::vector<int32> & cambiumidsincontig,
    nextreadtoadd_t & resultread,
    clock_t & maxallowedclocks);

  void   n4_buildESTContigStartWithStrongGood(
    std::vector<Align> & aligncache,
    int32 startid,
    Contig & con);


public:
//  Pathfinder(MIRAParameters* params,
//	     ReadPool & readpool,
//	     const overlap_edges_t & forward_edges,
////	     const slist<AlignedDualSeqFacts> & adslist,
//	     bannedoverlappairs_t & tmp_srmb_overlap_bans);

  Pathfinder(std::vector<MIRAParameters> * params,
	     ReadPool & readpool,
	     std::vector<newedges_t> & forward_edges,
	     std::vector<AlignedDualSeqFacts> & adsfacts);

  Pathfinder(Pathfinder const &other);
  ~Pathfinder();

  Pathfinder const & operator=(Pathfinder const & other);
  friend std::ostream & operator<<(std::ostream &ostr, Pathfinder const &soc);

  void discard();

  void constructStepByStep(std::vector<Align> & aligncache,
			   std::vector<int8> * used_ids,
			   std::vector<int32> * ids_in_contig, // this is a return value
			   std::vector<uint8> * multicopies,
			   std::vector<uint8> * hasmcoverlaps,
			   std::vector<uint8> * istroublemaker,
			   std::vector<std::vector<newedges_t>::iterator > * lowerbound_oedges_ptr,
			   Contig & con);

  uint32 getReadAddAttempts() const {return PAF_readaddattempts;}

//  void n3_constructStepByStep(
//    std::vector<Align> & aligncache,
//    std::vector<int8> * used_ids,
//    std::vector<uint8> * multicopies,
//    std::vector<uint8> * hasmcoverlaps,
//    std::vector<uint8> * istroublemaker,
//    std::vector<std::vector<newedges_t>::iterator > * lowerbound_oedges_ptr,
//    std::vector<int32> * proposedSRMrightclips,
//    Contig & con);


  void n4_constructStepByStep(std::vector<Align> & aligncache,
			      std::vector<int8> * used_ids,
			      std::vector<int32> * ids_in_contig, // this is a return value
			      std::vector<uint8> * multicopies,
			      std::vector<uint8> * hasmcoverlaps,
			      std::vector<uint8> * istroublemaker,
			      std::vector<std::vector<newedges_t>::iterator > * lowerbound_oedges_ptr,
			      std::vector<uint8> * wellconnected,
			      Contig & con);

  size_t n4_getNumElementsInStartCache();

};


#endif
