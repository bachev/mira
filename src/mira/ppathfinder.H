/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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

#ifndef _mira_ppathfinder_h_
#define _mira_ppathfinder_h_

#include <array>
#include <queue>
#include <utility>

#include "stdinc/defines.H"

#include "mira/overlapedges.H"
#include "mira/contig.H"


class PPathfinder
{
  // Types
private:
  struct beststartinfo_t {
    uint32 bsi_clustersize;
    readid_t bsi_rid;
    // don't need to add that one, just for debugging
    uint32 bsi_numconnects;

    friend std::ostream & operator<<(std::ostream &ostr, beststartinfo_t const &bsi){
      ostr << "BSI: rid " << bsi.bsi_rid
	   << "\tcs:" << bsi.bsi_clustersize
	   << "\tnc: " << bsi.bsi_numconnects << std::endl;
      return ostr;
    }

    inline static bool ltclustersize(const beststartinfo_t & a, const beststartinfo_t & b){
      return a.bsi_clustersize < b.bsi_clustersize;
    }

  };


  struct nextreadtoadd_t{
    readid_t refid;
    readid_t newid;                   // number of read in seqpool
    uint32 weight;
    int16  direction_newid;
    const AlignedDualSeqFacts * ads_node;

    uint8 foundqueuenum;
    friend std::ostream & operator<<(std::ostream &ostr, nextreadtoadd_t const &nrta);
  };


  typedef std::pair<uint32,necontainer_t::iterator> ppfweightelem_t;
  typedef std::priority_queue<ppfweightelem_t> ppfweightqueue_t;

  //Variables
private:
  enum queue_types {
    QTG_TPARTNERNOREPT_OLNOREPTSTRONG_WELLCONNECTED=0,     // 0    // TPARTNER: no rept; THIS OVERLAP: norept; NEW READ: well connected
    QTG_TPARTNERNOREPT_OLNOREPTWEAK_WELLCONNECTED,	    // 1
    QTG_TPARTNERNOREPT_OLNOREPTBELOWAVG_WELLCONNECTED,	    // 2

    QTG_TPARTNERNOREPT_OLNOREPTSTRONG,      // 3    // TPARTNER: no rept; THIS OVERLAP: norept
    QTG_TPARTNERNOREPT_OLNOREPTWEAK,	    // 4
    QTG_TPARTNERNOREPT_OLNOREPTBELOWAVG,    // 5

    QTG_TPARTNERNOTREPT_OLNOREPTSTRONG_WELLCONNECTED,      // 6    // TPARTNER: HAF4; THIS OVERLAP: norept; NEW READ: well connected
    QTG_TPARTNERNOTREPT_OLNOREPTWEAK_WELLCONNECTED,	    // 7
    QTG_TPARTNERNOTREPT_OLNOREPTBELOWAVG_WELLCONNECTED,    // 8

    QTG_TPARTNERNOTREPT_OLNOREPTSTRONG,      // 9    // TPARTNER: HAF4; THIS OVERLAP: norept
    QTG_TPARTNERNOTREPT_OLNOREPTWEAK,	    // 10
    QTG_TPARTNERNOTREPT_OLNOREPTBELOWAVG,    // 11

    QTG_OLNOREPTSTRONG_WELLCONNECTED,			    // 12
    QTG_OLNOREPTWEAK_WELLCONNECTED,			    // 13
    QTG_OLNOREPTBELOWAVG_WELLCONNECTED,		    // 14

    QTG_TPARTNERNOREPT_STRONG_WELLCONNECTED,		    // 14a, maybe also a 14b and c (weak/belowavg)?

    QTG_OLNOREPTSTRONG,			    // 15
    QTG_OLNOREPTWEAK,			    // 16
    QTG_OLNOREPTBELOWAVG,		    // 17

    QTG_TPARTNERNOREPT_OLNOREPTOTHER,	    // 18
    QTG_TPARTNERNOREPT_OLNOTREPTOTHER,	    // 19
    QTG_TPARTNERNOTREPT_OLNOREPTOTHER,	    // 20
    QTG_TPARTNERNOTREPT_OLNOTREPTOTHER,	    // 21
    QTG_OLNOREPTOTHER,			    // 22

    QTG_TPARTNERNOREPT_OLREPT,		    // 23
    QTG_TPARTNERNOTREPT_OLREPT,		    // 24

    QTG_OLNOTREPT,                           // 25   // norept==false, but also rept==false (HAF4)

    QTG_OLREPTSTRONG_WELLCONNECTED,	    // 26
    QTG_OLREPTWEAK_WELLCONNECTED,	    // 27
    QTG_OLREPTSTRONG,			    // 28
    QTG_OLREPTWEAK,			    // 29
    QTG_OLREPT,				    // 30

    QTG_MAYBESPOILSPORT,                    // 31

    QTG_INCORRECTIBLEORCHIMERA,             // 32 keep this as very, very last resort option

    QTG_OLOTHER,                             // 33   should stay empty??

    QTG_RELEGATEDBYPP,                      // 34

    QTG_END = 36,


    QTE_TPARTNERWCREPT6P_WCREPT6P_OLREPT_REPT6PWC=0,  //  0
    QTE_TPARTNERWCREPT5_WCREPT5_OLREPT_REPT6PWC,      //  1
    QTE_TPARTNERWCREPT5_WCREPT5_OLREPT_REPT5WC,	      //  2

    QTE_TPARTNERWC_WCREPT6P_OLREPT_REPT6PWC,	      //  3
    QTE_TPARTNERWC_WCREPT5_OLREPT_REPT6PWC,	      //  4
    QTE_TPARTNERWC_WCREPT5_OLREPT_REPT5WC,	      //  5

    QTE_TPARTNERWC_WC_OLREPT_REPT6PWC,		      //  6
    QTE_TPARTNERWC_WC_OLREPT_REPT5WC,		      //  7

    QTE_TPARTNERWC_WC_OLREPT_WC,		      //  8
 ///--- OL_REPT && TPARTNERWC ^^^
    QTE_WCREPT6P_OLREPT_REPT6PWC,		      //  9
    QTE_WCREPT5_OLREPT_REPT6PWC,		      // 10
    QTE_WCREPT5_OLREPT_REPT5WC,			      // 11

    QTE_REPT6P_OLREPT_REPT6PWC,			      // 12
    QTE_REPT5_OLREPT_REPT6PWC,			      // 13
    QTE_REPT5_OLREPT_REPT5WC,			      // 14

    QTE_TPARTNER_WC_OLREPT_WC,			      // 15
    QTE_TPARTNER_OLREPT_WC,			      // 16

    QTE_WC_OLREPT_WC,				      // 17
    QTE_OLREPT_WC,				      // 18
///--- OLREPT ^^^
    QTE_TPARTNERWC_WC_OL_WC,			      // 19
    QTE_TPARTNER_WC_OL_WC,			      // 20
    QTE_TPARTNER_OL_WC,				      // 21
    QTE_TPARTNER,				      // 22  // nope for EST ... intron still may spoil

    QTE_WC_OL_WC,				      // 23
    QTE_OL_WC,					      // 24

    QTE_REPT6P,					      // 25
    QTE_REPT5,					      // 26
    QTE_H4,					      // 27
    QTE_H3,					      // 28

    QTE_INCORRECTIBLEORCHIMERA,                       // 29 keep this as very, very last resort option

    QTE_OTHER,                                        // 30    should stay empty??

    QTE_RELEGATEDBYPP,                                // 31

    QTE_END
  };

  std::array<ppfweightqueue_t,QTG_END> PPF_queues; // careful in case QTE_END is bigger!

  std::vector<MIRAParameters> * PPF_miraparams_ptr;
  ReadPool * PPF_readpool_ptr;
  necontainer_t  * PPF_overlap_edges_ptr;
  adsfcontainer_t * PPF_adsfacts_ptr;
  std::vector<Align> * PPF_aligncache_ptr;

  std::vector<int8>  * PPF_used_ids_ptr;
  std::vector<uint8> * PPF_multicopies_ptr; /* reads with more overlaps than
					  expected on average have 1
					  provided by assembly class
					  pathfinder will start building
					  elsewhere, and include those
					  last
				       */
  std::vector<uint8> * PPF_hasmcoverlap_ptr; /* reads that overlap with a read
					   that is categorised as multi-
					   copy get 1 here
					   initialised by pathfinder if
					   vector is empty (==once
					   every pass of MIRA)
					   // TODO: unused now???
					*/
  std::vector<uint8> * PPF_hasreptoverlap_ptr; /* reads that have ol_rept overlaps
						  init by pathfinder if empty (once every pass)
						  maybe unused atm
					      */
  std::vector<uint8> * PPF_hasnoreptoverlap_ptr; /* reads that have ol_norept overlaps
						    init by pathfinder if empty (once every pass)
						 */
  std::vector<uint8> * PPF_istroublemaker_ptr;
  std::vector<uint8> * PPF_incorrectibleorchimera_ptr;
  std::vector<uint8> * PPF_maybespoilsport_ptr;
  std::vector<uint8> * PPF_wellconnected_ptr; /* vector of bool values
					    reads having a level 0 overlap
					    criterion with left extend and
					    level 0 for right extend are
					    well connected
					    (see skim and assembly for level
					    settings)
					 */

  /* lower_bound() is called very often, so caching/precomputing all
     possible results (= number of reads in readpool)
     speeds the things up
  */
  std::vector<necontainer_t::iterator> * PPF_lowerbound_oedges_ptr;
  std::vector<Contig::templateguessinfo_t> * PPF_astemplateguesses_ptr; // beware, may be rightfully empty


  // atm only for mapping
  uint32 PPF_wantscleanoverlapends;   // but this one could also be used for de-novo I think
  // special for mapping
  uint32 PPF_mintotalnonmatches;
  // special for mapping
  uint8 PPF_allowedseqtype;

  /*********************************************************************************
   *********************************************************************************/

  Contig * PPF_actcontig_ptr;

  /*********************************************************************************
   *********************************************************************************/

  const pathfinder_parameters * PPF_pafparams_ptr;
  std::vector<beststartinfo_t> PPF_beststartcache;;

  // size of number of reads in assembly
  // tells what overlap type led to a read being added to a contig
  // double function: lookup for which read is in contig
  // 4: added via ol_norept
  // 3: added via ol_norept==false && ol_rept==false
  // 2: added via other
  // 1: added by backbone
  // 0: not added
  enum {
    ADDED_NOTADDED=0,
    ADDED_BY_BACKBONE,
    ADDED_BY_OTHER,
    ADDED_BY_NOTREPT,
    ADDED_BY_NOREPT
  };
  std::vector<int8>    PPF_ids_added_oltype;

  // simple list of ids that were added to a contig, in the order of addition
  //  used as a quick lookup for cleanup routines which need to restore
  //  pristine, all-zero PPF_ids_added_oltype, PPF_blacklisted_ids etc vectors.
  std::vector<readid_t>    PPF_ids_in_contig_list;

  // simple list of rail reads in a contig when mapping
  std::vector<readid_t>    PPF_rails_in_contig_list;

  // A simple queue for rail overlaps, using deque as STL queue has no clear() :-(
  //  IDs contained here are from non-rails, but have an overlap to a rail
  std::deque<readid_t>     PPF_railoverlapcache;
  // bool yes/no whether a rail already used
  std::vector<uint8>       PPF_tmpproc_readalreadyrailed;

  Contig::errorstatus_t PPF_contigerrstat;

  /* set by n4_searchBestStrongGoodStartEnvironment_sub()
     says what the startcache contains
     0 == best quality, strong good overlaps found, well connected reads
     1 == medium quality, overlaps but not strong good, not multicopy, well connected reads
     2 == medium quality, overlaps but not strong good. multicopy, well connected reads
     3 == not well connected
     4 == singlets
  */
  enum {
    BSCC_GENOME_BESTQUAL=0,
    BSCC_GENOME_MEDQUAL_TROUBLEMAKER_NOTSTRONG_NOTMULTICOPY_WELLCONNECTED_NOKMERFORK,
    BSCC_GENOME_MEDQUAL_TROUBLEMAKER_NOTSTRONG_MULTICOPY_WELLCONNECTED_NOKMERFORK,
    BSCC_GENOME_BADQUAL_NOTWELLCONNECTED,
    BSCC_EST_WELLCONNECTED_REPT6P,
    BSCC_EST_WELLCONNECTED_REPT,
    BSCC_EST_WELLCONNECTED_H4,
    BSCC_EST_WELLCONNECTED_H3,
    BSCC_EST_WELLCONNECTED_H2,
    BSCC_EST_NOTWELLCONNECTED,
    BSCC_SINGLETS,
  };
  uint8 PPF_bsccontent;
  bool  PPF_bsrandry;    // whether startcache had to be refilled in last call to Pathfinder

  // small store is a store with a certain capacity (100k?) that takes
  //  up iterators to all banned overlaps until full.
  // So when it comes to reset pf_banned flag, if store is not full only
  //  the overlap edges whose iterators are in this store have to be cleared
  //  instead of iterating through the complete overlap edges. The latter
  //  is needed if store is full.
  // on de-novo 6m Solexa paired-end 36mers, this reduces another
  //  15% in the contig building phase (38 minutes to 33 minutes)
  // tradeoff of 5 minutes for 800k RAM is pretty good (160k RAM would
  //  suffice in that example, but let's plan for a bit more)
  // influence of mechanism grows with number of overlaps, so
  //  assemblies with more reads will profit exponentially
  std::vector<necontainer_t::iterator> PPF_overlapsbanned_smallstore;

  // blacklisting: queue to handle blacklist decay ...
  std::queue<std::vector<readid_t>> PPF_blacklist_queues;
  // ... and bool vector (size of readpool) for quick access
  // 0: normal; 1: blacklisted
  std::vector<int8> PPF_blacklisted_ids;

  // temporary vector to be used by routines which need an array of size of readpool
  // initiliased to all 0 by ppathfinder, but routines must take care themselves
  //  to clean up after they used it!
  std::vector<uint8> PPF_tmparray;
  std::vector<uint32> PPF_tmparray_idxused;

  // currently for EST data
  // max and min HAF level in reads
  std::vector<uint8> PPF_haflevel_max;
  std::vector<uint8> PPF_haflevel_min;
  uint8 PPF_haflevel_maxseen;
  uint8 PPF_haflevel_minseen;

  // counts how many read add attempts there were to the contig (contig.addRead())
  //  in the last call to map() or denovo() (needed for Solexa mapping assembly,
  //  where loops can stop if no more adding attempt exists)
  uint32 PPF_readaddattempts;

  /*********************************************************************************
   *********************************************************************************/

  uint32 PPF_buildcontig_newlinecounter;
  suseconds_t PPF_timing_pathsearch;
  suseconds_t PPF_timing_connadd;

  static bool PPF_staticinit;

public:


  //Functions
private:
  static bool staticInit();

  void priv_ppFillNoRept();
  void priv_ppFillRept();
  void priv_ppFillSpoilSport();

  void priv_initialiseLowerBoundOEdges();
  void priv_showProgress();
  void priv_basicSetup();

  void priv_fillDenovoStartCache() {
    if(PPF_beststartcache.empty()){
      PPF_bsrandry=true;
      if(PPF_pafparams_ptr->paf_use_genomic_algorithms) {
	priv_fdns_genome();
      }else{
	priv_fdns_est();
      }
    }
  }
  void priv_fdns_genome();
  void priv_fdns_g_subFillCache(bool wanttroublemakercheck,
				bool wantstronggoodcheck,
				bool wantmulticopycheck,
				bool wantwellconnectedcheck,
				bool wantnokmerfork);
  void priv_fdns_est();
  void priv_fdns_e_subFillCache(bool wantwellconnectedcheck,
				uint8 minallowedfreq);
  void priv_fdns_fillSinglets();

  readid_t priv_getNextStartID();
  readid_t priv_gnsi_helper();
  void priv_fillDenovoQueues(size_t maxdist);
  inline uint32 priv_insertRIDIntoDenovoQueues(readid_t rid){
    if(PPF_pafparams_ptr->paf_use_genomic_algorithms) {
      return priv_iridnq_genome(rid);
    }
    return priv_iridnq_est(rid);
  }
  uint32 priv_iridnq_genome(readid_t rid);
  uint32 priv_iridnq_est(readid_t rid);

  uint32 priv_getNextOverlapFromDenovoQueue(necontainer_t::iterator & oeI);

  void priv_loopDenovo(){
    priv_ld_genome_and_est();
  }
  void priv_ld_genome_and_est();

  void priv_handleReadNotAligned(necontainer_t::iterator oeI, nextreadtoadd_t const &nrta);
  void priv_munchBlacklist(bool force);

  void priv_fillHAFLevelInfo();

  void priv_prepareRailOverlapCache();
  void priv_prochelper1(uint8 seqtype);

  necontainer_t::iterator priv_findNextBackboneOverlapQuick(nextreadtoadd_t & resultread,
								      bool allowmulticopies,
								      bool allowtroublemakers,
								      bool allowsmallhits);

  void priv_storeTemplateGuess(readid_t newid, Contig::templateguessinfo_t & tguess);

public:
  PPathfinder(std::vector<MIRAParameters> * params,
	      ReadPool * readpool,
	      necontainer_t * overlap_edges,
	      adsfcontainer_t * adsfacts,
	      std::vector<Align> * aligncache,
	      std::vector<int8> * used_ids,
	      std::vector<uint8> * multicopies,
	      std::vector<uint8> * hasmcoverlaps,
	      std::vector<uint8> * hasreptoverlaps,
	      std::vector<uint8> * hasnoreptoverlaps,
	      std::vector<uint8> * istroublemaker,
	      std::vector<uint8> * incorrectibleorchimera,
	      std::vector<uint8> * maybespoilsport,
	      std::vector<uint8> * wellconnected,
	      std::vector<necontainer_t::iterator > * lowerbound_oedges_ptr,
	      std::vector<Contig::templateguessinfo_t> * astemplateguess
    );

  ~PPathfinder();

  PPathfinder(PPathfinder const &other) = delete;
  PPathfinder const & operator=(PPathfinder const & other) = delete;

  //friend std::ostream & operator<<(std::ostream &ostr, PPathfinder const &ppf);

  void discard();

  void prepareForNewContig(Contig & con);
  void resyncContig();

  void denovo();
  void map();
  void mapAndDenovo();

  const std::vector<readid_t> & getRIDsKnownInContig() const { return PPF_ids_in_contig_list;}
  bool startCacheRanDry() const { return PPF_bsrandry;}
  bool startCacheHasSinglets() const { return PPF_bsccontent==BSCC_SINGLETS;}
  uint32 getReadAddAttempts() const {return PPF_readaddattempts;}

  void setWantsCleanOverlapEnds(uint32 len) {PPF_wantscleanoverlapends=len;}
  void setMinTotalNonMatches(uint32 n) {PPF_mintotalnonmatches=n;}
  void setAllowedSeqTypeForMapping(uint8 st) {PPF_allowedseqtype=st;}
};


#endif
