/*
 * Written by Bastien Chevreux (BaCh)
 *
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


// 	$Id$


#ifndef _mira_structs_h_
#define _mira_structs_h_

#include <string>
#include <vector>

//#include "stdinc/stlincludes.H"
#include "stdinc/types.H"
#include "enums.H"


// ----------------- For parametrisation

struct special_parameters
{
  uint32 sp_est_startstep;

  bool mi_iknowwhatido;

  bool mi_extended_log;
  bool mi_extra_flag1;
  bool mi_extra_flag2;
  bool mi_extra_flag3;

  uint32 mi_as_largecontigsize;
  uint32 mi_as_largecontigsize4stats;

  bool sp_parse_checktechnologypresence;
};


enum { NWNONE=0, NWSTOP, NWWARN };
struct nagandwarn_parameters
{
  uint8  nw_check_nfs;
  uint8  nw_check_templateproblems;
  uint8  nw_check_multipassmapping;

  uint8  nw_check_maxreadnamelength;
  uint32 nw_check_mrnlvalue;
  uint8  nw_check_readnamesra;
  uint8  nw_check_duplicatereadnames;

  uint8  nw_check_illuminajunkinest;
  uint8  nw_check_proposedendclipping;

  uint8  nw_check_coverage;
  uint32 nw_check_covvalue;

  std::string nw_warn_x174prefix;
};

struct assembly_parameters
{
//  uint32 as_pattern_length;     // Length of the pattern to use for matching
//  uint32 as_num_end_errors;      // # of errors allowed when matching ends
//  uint32 as_num_middle_errors;   // # of errors allowed when matching middle parts
//  uint32 as_num_middle_runs;     // # of middle part runs to make. at this time only exactly 1 supported

  uint32 as_numthreads;

  uint64 as_amm_maxprocesssize;
  uint64 as_amm_keeppercentfree;

  uint32 as_cleanup_tmpfiles;  // clean up the temporary files after assembly?

  uint32 as_minimum_readlength;      // min. length of reads to be used in assembly
  uint32 as_minimum_readspercontig;      // min. reads per contig

  uint32 as_backbone_raillength;   // raillength for backbone
  uint32 as_backbone_railoverlap;   // railoverlap for backbone
  uint32 as_backbone_outlen;       // "rail"length for backbone in output
  int32 as_startbackboneusage_inpass; // in which pass to use backbones first

  bool   as_backbone_bootstrapnewbackbone;
  bool   as_backbone_trimoverhangingreads;
  bool   as_automemmanagement;

  bool   as_assemblyjob_accurate;
  bool   as_assemblyjob_mapping;
  bool   as_assemblyjob_preprocessonly;
  uint32  as_numpasses;
  uint32  as_maxcontigsperpass;
  uint32  as_numrmbbreakloops;
  bool   as_filecheck_only;

  std::vector<uint32> as_bphseries;

  bool   as_enforce_qualsinreads;
  bool   as_wants_qualityfile;

  bool   as_savesimplesingletsinproject;
  bool   as_savetaggedsingletsinproject;

  bool   as_output_html;
  bool   as_output_txt;
  bool   as_output_caf;
  bool   as_output_maf;
  bool   as_output_tcs;
  bool   as_output_ace;
  bool   as_output_gff3;
  bool   as_output_wiggle;
  bool   as_output_fasta;
  bool   as_output_gap4da;
  bool   as_output_tmp_html;
  bool   as_output_tmp_tcs;
  bool   as_output_tmp_txt;
  bool   as_output_tmp_caf;
  bool   as_output_tmp_maf;
  bool   as_output_tmp_ace;
  bool   as_output_tmp_fasta;
  bool   as_output_tmp_gap4da;
  bool   as_output_exttmp_alsosinglets;
  bool   as_output_exttmp_html;
  bool   as_output_exttmp_caf;
  bool   as_output_exttmp_ace;
  bool   as_output_exttmp_fasta;
  bool   as_output_exttmp_gap4da;

  bool   as_output_removerollovertmps;
  bool   as_output_removetmpdir;

  bool   as_clip_lowercase_front;
  bool   as_clip_lowercase_back;

  bool   as_clip_sdbg_chimeradetection;
  bool   as_clip_kmer_junkdetection;
  bool   as_clip_kmer_junkkill;

  bool   as_clip_skimchimeradetection;
  bool   as_clip_skimjunkdetection;

  bool   as_clip_knownadaptorsright;
  bool   as_search_phix174;
  bool   as_filter_phix174;
  bool   as_filter_rrna;
  bool   as_filter_rrna_pairs;
  uint32 as_filter_rrna_numkmers;

  bool   as_clip_proposeendclips;
  bool   as_clip_pec_continuous;
  bool   as_clip_pec_sxaggcxg;
  bool   as_clip_pec_ffr; // front forward-rev
  bool   as_clip_pec_bfr; // back  forward-rev
  bool   as_clip_pec_fcmst; // front confirmed multiple seq-types
  bool   as_clip_pec_bcmst; // back  confirmed multiple seq-types
  bool   as_clip_pec_fsalp; // front seen at low pos, 27.04.2012 not a good criterion, should not be used
  bool   as_clip_pec_bsalp; // back  seen at low pos, 27.04.2012 not a good criterion, should not be used
  uint32 as_clip_pec_basesperhash;
  uint16 as_clip_pec_ffreq;     // front freq
  uint16 as_clip_pec_bfreq;     // back  freq
  uint32 as_clip_pec_mkfr;  // min kmer for forward-rev criterion
  uint32 as_clip_pec_mtk;  // min total kmer (both directions)

  bool   as_clip_badsolexaends;

  bool   as_clip_possible_vectors;
  bool   as_clip_quality;
  bool   as_clip_badstretchquality;
  bool   as_clip_maskedbases;

  bool   as_put_asswithmira_tags;

  base_quality_t as_clip_quality_minqual;
  base_quality_t as_clip_badstretchquality_minqual;

  base_quality_t as_clip_quality_minthreshold;
  uint32         as_clip_quality_numminthreshold;

  int32  as_clip_quality_winlen;
  int32  as_clip_badstretchquality_winlen;

  int32  as_clip_vector_maxlenallowed;
  int32  as_clip_maskedbase_gapsize;
  int32  as_clip_maskedbase_maxfrontgap;
  int32  as_clip_maskedbase_maxendgap;
  int32  as_clip_ssahamerge_gapsize;
  int32  as_clip_ssahamerge_maxfrontgap;
  int32  as_clip_ssahamerge_maxendgap;
  int32  as_clip_ssahamerge_strictfrontclip;
  int32  as_clip_ssahamerge_strictendclip;
  int32  as_clip_minslrequired;
  int32  as_clip_minqlsetto;
  int32  as_clip_minsrrequired;
  int32  as_clip_minqrsetto;
  bool   as_clip_ensureminimumleftclipoff;
  bool   as_clip_ensureminimumrightclipoff;

  bool    as_clip_polyat;
  bool    as_clip_polyat_keeppolystretch;
  uint32  as_clip_polyat_len;
  uint32  as_clip_polyat_maxerrors;
  uint32  as_clip_polyat_maxgap;

  bool    as_clip_3ppolybase;
  uint32  as_clip_3ppolybase_len;
  uint32  as_clip_3ppolybase_maxerrors;
  uint32  as_clip_3ppolybase_maxgap;

  //uint32  as_clipmask_rarekmers;
  bool    as_clipmask_rarekmers;

  bool    as_use_read_extension;
  uint32  as_readextension_window_len;
  uint32  as_readextension_window_maxerrors;
  uint32  as_readextension_firstpassnum;
  uint32  as_readextension_lastpassnum;

  bool   as_mark_repeats;
  bool   as_mark_repeats_onlyinresult;
//  bool   as_skimeachpass;

  bool   as_spoilerdetection;
  int32  as_spdetect_lastpassonly;

  bool   as_automatic_repeat_detection;
  double as_ard_multicopythreshold;
  int32  as_ard_multicopyminlen;
  int32  as_ard_multicopygrace;
  bool   as_uniform_read_distribution;
  int32  as_urd_startinpass;
  double as_urd_cutoffmultiplier;

  bool   as_dateoutput;

  bool as_backbone_strainname_forceforall;

  bool as_buntify_reads;
  bool as_rle_reads;

  std::string as_user_kmerseries;

  std::string as_projectname_in;
  std::string as_projectname_out;
  std::string as_tmpf_unused_ids;
  std::string as_tmpf_ads;
  std::string as_tmpf_adsextend;
  std::string as_tmpf_scfreadfail;
  std::string as_tmpf_scfreadfatallywrong;
  std::string as_tmpf_vectorclip;
  std::string as_tmpf_spoiler;
  std::string as_tmpf_posmatch;
  std::string as_tmpf_poolinfo;
  std::string as_tmpf_normalisedskim;
  std::string as_tmpf_clippings;
  std::string as_tmpf_kmerstatistics;

  ////////// needed for extended resuming
  std::string as_tmpf_signal_kmerstats;
  std::string as_tmpf_signal_findpossibleoverlaps;
  std::string as_tmpf_signal_mainalignments;

  std::string as_tmpf_wellconnected;
  std::string as_tmpf_banned_overlaps;
  std::string as_tmpf_istroublemaker;
  std::string as_tmpf_needalloverlaps;
  std::string as_tmpf_multicopies;
  std::string as_tmpf_hasmcoverlap;
  std::string as_tmpf_debrisreason;
  std::string as_tmpf_skimmegahubs;

  //////////

  std::string as_backbone_rail_fromstrain;

  std::string as_outfile_FASTA;
  std::string as_outfile_FASTAUNPADDED;
  std::string as_outfile_FASTAUNPADDEDQUAL;
  std::string as_outfile_FASTAPADDED;
  std::string as_outfile_FASTAPADDEDQUAL;
  std::string as_outfile_CAF;
  std::string as_outfile_MAF;
  std::string as_outfile_TCS;
  std::string as_outfile_ACE;
  std::string as_outfile_WIGGLE;
  std::string as_outfile_HTML;
  std::string as_outfile_TXT;
  std::string as_outdir_GAP4DA;

  std::string as_outfile_stats_contigstats;
  std::string as_outfile_stats_contigtags;
  std::string as_outfile_stats_crlist;
  std::string as_outfile_stats_debrislist;
  std::string as_outfile_stats_featureanalysis;
  std::string as_outfile_stats_featuresummary;
  std::string as_outfile_stats_featuresequences;
  std::string as_outfile_stats_featurecoverage;
  std::string as_outfile_stats_info;
  std::string as_outfile_stats_warnings;
  std::string as_outfile_stats_largecontigs;
  std::string as_outfile_stats_snpanalysis;
  std::string as_outfile_stats_snpenvironment;
  std::string as_outfile_stats_reads_invalid;
  std::string as_outfile_stats_reads_tooshort;
  std::string as_outfile_stats_readrepeats;
  std::string as_outfile_stats_readtags;
  std::string as_outfile_stats_rgstinfo;
};

struct finalmapping_parameters
{
  bool  fm_active;
  bool  fm_mapperfect;
  int32 fm_maxtotalerrors; /* maximum errors an alignment of Solexa sequences
			      may have. If -FM:me is < 0, automatically
			      determining optimal value.
			    */
  int32 fm_maxmismatches; /* maximum mismatch errors an alignment of Solexa
			     sequences may have
			     if <0, maxtotalerrors is implicit
			  */
  int32 fm_maxgaps; /* maximum gap errors an alignment of Solexa sequences may
		       have
		       if <0, maxtotalerrors is implicit
		     */
  uint32 fm_clean_end_dist; /* like ads_clean_end_distance, just for final mapping */
};

struct directory_parameters
{
  std::string dir_top;
  std::string dir_tmp;
  std::string dir_tmp_redirectedto;
  std::string dir_results;
  std::string dir_info;
  std::string dir_checkpoint;
  std::string dir_checkpoint_tmp;

  std::string dir_tmp_symlink;
};


struct file_parameters
{
  std::string chkpt_persistentoverlaps;
  std::string chkpt_bannedoverlaps;
  std::string chkpt_maxcovreached;
  std::string chkpt_passinfo;
  std::string chkpt_readpool;
};


struct skim_parameters
{
  uint32 sk_numthreads;
  uint32 sk_basesperhash;
  uint32 sk_bph_max;
  uint32 sk_bph_increasestep;
  uint32 sk_hashsavestepping;
  int32  sk_percentrequired;
  uint32 sk_maxhitsperread;

  bool   sk_filtermegahubs;
  double sk_maxmegahubratio;
  uint32 sk_megahubcap;

  bool   sk_alsoskimrevcomp;

  bool   sk_swcheckonbackbones;

  //uint32 sk_maxhashesinmem;
  uint32 sk_maxhashesinmem;
  uint32 sk_memcaphitreduction;
};


struct hashstatistics_parameters
{
  double hs_freqest_minnormal;
  double hs_freqest_maxnormal;
  double hs_freqest_repeat;
  double hs_freqest_heavyrepeat;
  double hs_freqest_crazyrepeat;
  uint32 hs_freq_covestmin;

  bool   hs_masknastyrepeats;
  uint32 hs_nastyrepeatratio;
  uint32 hs_nastyrepeatcoverage;

  uint32 hs_repeatlevel_in_infofile;

  bool   hs_apply_digitalnormalisation;

  uint32 hs_memtouse;
  uint32 hs_rare_kmer_final_kill;
};



struct align_parameters
{
  uint32 al_max_cutoff;	        /* failsafe. if too many solutions exist,
				   stop calculations here */

  // The following values will all need to be reached if a solution
  //  has to be regarded as valid.
  int32 al_min_score;	        /* Minimum score for two aligned sequences
				   to be accepted as valid solution */
  int32 al_min_relscore;	/* Percentage minimum of the score of two
				   alignments to have to be recognised as
				   valid solution */
  int32 al_min_overlap;     // minimum overlap of two alignments in bases

  int32  al_kmin;           // minimum number of errors in banded Smith-Waterman
  int32  al_kmax;           // maximum ...
  int32  al_kpercent;       /* k errors deduced by the length of expected overlap in %
			       10% and 150 ovelap -> k= 15 */


  bool   ads_enforce_clean_ends; /* when true, ads will tolerate no ACGT/ACGT
				    mismatch within distance of ends of
				    overlap */
  int32  ads_clean_end_distance; // distance for above
  int32  ads_clean_end_mismatchallowed; // how many mismatches/gaps are allowed in ends


  bool   ads_extra_gap_penalty;   // in ADS, use extra penalty function for long gaps?
  std::string ads_gp_functionstring; // the below as given by user
  std::vector<uint32> ads_gp_function;  // gap penalty level;
  uint32 ads_max_gppercent;  // Max. gap penalty in percent

//  bool   ads_extra_mismatch_penalty;
//  uint32 ads_emp_windowlen;
//  uint32 ads_emp_maxmismatches;
};

struct dynamic_parameters
{
  int32 dyn_score_multiplier;   /* by which factor the score differs from a normal
				   +1/-1/-2 matrix */


  int32 dyn_score_match;		// Scores for a match
  int32 dyn_score_mismatch;      //   "  a mismatch
  int32 dyn_score_halfmatch;	 // Scores for a half match (e.g. M-A or M-V)
  int32 dyn_score_halfmismatch;  //   "  a half mismatch (e.g. M-D)
  int32 dyn_score_nmatch;        //   "  a match with N
  int32 dyn_score_npenaltymatch; //   "  a match with N (light penalty)
  int32 dyn_score_gap;           //   "  gap penalty
  int32 dyn_score_oldgap;        //   "  gap penalty
  int32 dyn_score_oldgapmatch;   //   "  gap penalty, but base "below" gap match
  int32 dyn_score_ltermgap;      //   "  left terminal gap penalties
  int32 dyn_score_rtermgap;      //   "  right terminal gap penalties
  int32 dyn_score_termgap;      //   "  (min from lterm and rterm)
};


struct pathfinder_parameters
{
  bool paf_use_genomic_algorithms;

  bool paf_use_emergency_blacklist;
  bool paf_use_emergency_search_stop;
  bool paf_use_max_contig_buildtime;

  bool   paf_use_quick_rule;
  int8   paf_quickrule_minsim1;
  int8   paf_quickrule_minsim2;
  int32  paf_quickrule_minlen1;
  int32  paf_quickrule_minlen2;

  uint32 paf_bbquickoverlap_minlen;

  uint32 paf_minreadstocachebuildstart;

  int32 paf_max_contig_buildtime;
  int32 paf_ess_depth;

  long int paf_max_startcache_filltime;

  int16 paf_skipwholecontigscan;

  clock_t paf_maxcontigclockticks;
  clock_t paf_nextread_maxcttoess;
};

struct contig_parameters
{
  std::string con_nameprefix;
  int32 con_reject_on_drop_in_relscore;
  int32 con_min_relscore;	/* Percentage minimum of the score of two
				   alignments to have to be recognised as
				   valid solution */
  int32 con_danger_analyse_mode;        // analysejob
  int32 con_danger_max_error_rate;      /* reject when more than x% bases
					   from danger zones in read to add
					   don't match
					   the actual consensus */
  bool  con_danger_analyse_Ns;          // true when N also have to be analysed
//  int32 con_rmb_numzone_trigger;       /* how many rmb zones are needed to
//					   be sure that it's really a rmb */
//  int32 con_min_groupqual_srmbwrmb_change;
  base_quality_t con_minrmbneighbourqual;
  base_quality_t con_mingroupqualforrmbtagging;
  uint32 con_mincoveragepercentage;
  int32 con_minreadspergroup;
  int32 con_endreadmarkexclusionarea;
  bool  con_emea_setzero_on_clipping_pec;

  int32 con_output_html_cpl;            // characters per line for HTML
  int32 con_output_text_cpl;            // characters per line for TEXT
  char  con_output_html_gapfill;            // endgap fill-character for TEXT
  char  con_output_text_gapfill;            // endgap fill-character for TEXT

  uint32 con_gap_override_ratio;

  bool  con_assume_snp_insteadof_rmb;
  bool  con_disregard_spurious_rmb_mismatches;
  bool  con_also_mark_gap_bases;
  bool  con_also_mark_gap_bases_needbothstrands;
  bool  con_also_mark_gap_bases_evenmc;   // even multicolumn
  bool  con_shorttagcomments;

  bool  con_force_nonIUPACconsensus; // for all seqtypes
  bool  con_force_nonIUPACconsensus_perseqtype;  // to be retired

  bool  con_mergeshortreads;
  int32 con_msr_keependsunmapped;
  int32 con_msr_maxerrors;
};

struct edit_parameters
{
  bool ed_mira_automatic_contic_editing;
  bool ed_homopolymer_overcalls;
  bool ed_kmer_singlets;

  bool ed_sdbg_readedit;

  bool ed_edit_automatic_contic_editing;
  bool ed_strict_editing_mode;
  uint32 ed_confirmation_threshold;
};


#endif
