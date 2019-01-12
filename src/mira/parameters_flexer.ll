option stack
%option noyywrap
%option c++
%option outfile="lex.yy.c"
%option prefix="MP"

%{
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

#include <fstream>
#include "parameters_tokens.h"

#include "enums.H"

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1110
#endif

  int filenameid=0;
  int singlechar_id=0;

%}

ANID [A-Za-z][A-Za-z0-9_]*
INT [+\-]?[0-9]+
FLOAT [0-9]+"."[0-9]+
FILENAME [/A-Za-z0-9_#.][\-/A-Za-z0-9_#.]+

%x JOB_MODE
%x GET_UNTIL_NEWLINE
%x GE_MODE
%x SB_MODE
%x LR_MODE
%x LR_EQMODE
%x LR_RNSMODE
%x CL_MODE
%x SR_MODE
%x FILE_MODE
%x OUTPUT_MODE
%x FN_MODE
%x ASK_YN_MODE
%x GET_SINGLECHAR_STRING_MODE
%x PAF_MODE
%x EDIT_MODE
%x SKIM_MODE
%x KMERSTAT_MODE
%x ALIGN_MODE
%x FINALMAP_MODE
%x AS_MODE
%x DP_MODE
%x CO_MODE
%x DIR_MODE
%x MI_MODE
%x NW_MODE
%x NW_CHOICEMODE
%x CO_VALMODE
%x CSV_NUMBERS_MODE
%x COMMENT_MODE
%x ASKFORTECH_NOCLIPPING
%%

<GE_MODE>"number_of_threads" |
<GE_MODE>"not"                   {return MP_as_numthreads;}
<GE_MODE>"keep_percent_memory_free" |
<GE_MODE>"kpmf"                   {return MP_as_amm_keeppercentfree;}
<GE_MODE>"max_process_size" |
<GE_MODE>"mps"                   {return MP_as_amm_maxprocesssize;}
<GE_MODE>"automatic_memory_management" |
<GE_MODE>"amm"                 { yy_push_state(ASK_YN_MODE); return MP_as_automemmanagement;}

<GE_MODE>"clean_tmp_files" |
<GE_MODE>"ctf"                 { yy_push_state(ASK_YN_MODE); return MP_as_cleanup_tmp_files;}

<GE_MODE>"print_date" |
<GE_MODE>"pd"              { yy_push_state(ASK_YN_MODE); return MP_as_nodateoutput;}
<GE_MODE>"est_snp_pipeline_step" |
<GE_MODE>"esps"               { return MP_sp_est_startstep;}
<GE_MODE>"colour_reads_by_kmer_frequency" |
<GE_MODE>"crkf"                 { yy_push_state(ASK_YN_MODE); return MP_as_buntify_reads;}
<GE_MODE>"preprocess_only" |
<GE_MODE>"ppo"                 { yy_push_state(ASK_YN_MODE); return MP_as_preprocess_only;}
<GE_MODE>"bang_on_throw" |
<GE_MODE>"bot"              { yy_push_state(ASK_YN_MODE); return MP_as_bangonthrow;}


<LR_MODE>"wants_quality_file" |
<LR_MODE>"wqf"              {yy_push_state(ASK_YN_MODE); return MP_as_wants_qualityfile;}

<LR_EQMODE>[nN][oO][nN][eE ]   { BEGIN(LR_MODE); return E_QUAL_NONE;}
<LR_EQMODE>[sS][cC][fF]        { BEGIN(LR_MODE); return E_QUAL_SCF;}
<LR_EQMODE>{ANID}              { BEGIN(0); return MP_UNRECOGNISED_STRING;}
<LR_EQMODE>{FLOAT}             { BEGIN(0); return MP_UNRECOGNISED_STRING;}
<LR_EQMODE>{INT}               { BEGIN(0); return MP_UNRECOGNISED_STRING;}
<LR_EQMODE>\n                  { BEGIN(0); return MP_UNRECOGNISED_STRING;}
<LR_MODE>"filecheck_only" |
<LR_MODE>"fo"                 { yy_push_state(ASK_YN_MODE); return MP_as_filecheck_only;}


<SB_MODE>"backbone_raillength" |
<SB_MODE>"brl"               { return MP_as_backbone_raillength;}
<SB_MODE>"backbone_railoverlap" |
<SB_MODE>"bro"               { return MP_as_backbone_railoverlap;}
<SB_MODE>"backbone_rail_fromstrain" |
<SB_MODE>"brfs"               { filenameid=MP_as_backbone_rail_fromstrain; yy_push_state(FN_MODE);}
<SB_MODE>"bsnffa" |
<SB_MODE>"backbone_strainname_forceforall"  { yy_push_state(ASK_YN_MODE); return MP_as_backbone_strainname_forceforall;}
<SB_MODE>"backbone_basequals" |
<SB_MODE>"bbq"               { return MP_as_backbone_basequals;}
<SB_MODE>"startbackboneusage_inpass" |
<SB_MODE>"sbuip"               { return MP_as_startbackboneusage_inpass;}
<SB_MODE>"trim_overhanging_reads" |
<SB_MODE>"tor"                 { yy_push_state(ASK_YN_MODE); return MP_as_backbone_trimoverhangingreads;}
<SB_MODE>"bootstrap_new_backbone" |
<SB_MODE>"bnb"                 { yy_push_state(ASK_YN_MODE); return MP_as_backbone_bootstrapnewbackbone;}


<CL_MODE>"msvs_strict_front_clip" |
<CL_MODE>"msvssfc"              { return MP_as_clip_ssahamerge_strictfrontclip;}
<CL_MODE>"msvs_strict_end_clip" |
<CL_MODE>"msvssec"              { return MP_as_clip_ssahamerge_strictendclip;}
<CL_MODE>"msvs_gap_size" |
<CL_MODE>"msvsgs"              { return MP_as_clip_ssahamerge_gapsize;}
<CL_MODE>"msvs_max_front_gap" |
<CL_MODE>"msvsmfg"             { return MP_as_clip_ssahamerge_maxfrontgap;}
<CL_MODE>"msvs_max_end_gap" |
<CL_MODE>"msvsmeg"             { return MP_as_clip_ssahamerge_maxendgap;}
<CL_MODE>"possible_vector_leftover_clip" |
<CL_MODE>"pvlc"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_possible_vectors;}
<CL_MODE>"pvc_maxlenallowed" |
<CL_MODE>"pvcmla"               { return MP_as_clip_vector_maxlenallowed;}
<CL_MODE>"bsqc_minimum_quality" |
<CL_MODE>"bsqcmq"               { return MP_as_clip_badstretchquality_minqual;}
<CL_MODE>"bsqc_window_length" |
<CL_MODE>"bsqcwl"               { return MP_as_clip_badstretchquality_winlen;}
<CL_MODE>"bad_stretch_quality_clip" |
<CL_MODE>"bsqc"                 { yy_push_state(ASK_YN_MODE); return MP_as_clip_badstretchquality;}
<CL_MODE>"qc_minimum_quality" |
<CL_MODE>"qcmq"               { return MP_as_clip_quality_minqual;}
<CL_MODE>"qc_window_length" |
<CL_MODE>"qcwl"               { return MP_as_clip_quality_winlen;}
<CL_MODE>"quality_clip" |
<CL_MODE>"qc"                 { yy_push_state(ASK_YN_MODE); return MP_as_clip_quality;}
<CL_MODE>"mbc_gap_size" |
<CL_MODE>"mbcgs"              { return MP_as_clip_maskedbases_gapsize;}
<CL_MODE>"mbc_max_front_gap" |
<CL_MODE>"mbcmfg"             { return MP_as_clip_maskedbases_maxfrontgap;}
<CL_MODE>"mbc_max_end_gap" |
<CL_MODE>"mbcmeg"             { return MP_as_clip_maskedbases_maxendgap;}
<CL_MODE>"maskedbase_clip" |
<CL_MODE>"mbc"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_maskedbases;}
<CL_MODE>"lowercase_clip_front" |
<CL_MODE>"lccf"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_lowercase_front;}
<CL_MODE>"lowercase_clip_back" |
<CL_MODE>"lccb"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_lowercase_back;}
<CL_MODE>"ensure_minimum_left_clip" |
<CL_MODE>"emlc"               { yy_push_state(ASK_YN_MODE); return MP_as_clip_ensureminimumleftclip;}
<CL_MODE>"minimum_left_clip_required" |
<CL_MODE>"mlcr"               { return MP_as_clip_minimumleftcliprequired;}
<CL_MODE>"set_minimum_left_clip_to" |
<CL_MODE>"smlc"               { return MP_as_clip_setminimumleftclip;}
<CL_MODE>"emrc"               { yy_push_state(ASK_YN_MODE); return MP_as_clip_ensureminimumrightclip;}
<CL_MODE>"minimum_right_clip_required" |
<CL_MODE>"mrcr"               { return MP_as_clip_minimumrightcliprequired;}
<CL_MODE>"set_minimum_right_clip_to" |
<CL_MODE>"smrc"               { return MP_as_clip_setminimumrightclip;}
<CL_MODE>"min_quality_threshold_for_entire_read_number_of_bases" |
<CL_MODE>"mqtfernob"               { return MP_as_clip_quality_numminthreshold;}
<CL_MODE>"min_quality_threshold_for_entire_read" |
<CL_MODE>"mqtfer"               { return MP_as_clip_quality_minthreshold;}

<CL_MODE>"clip_3ppolybase" |
<CL_MODE>"c3pp"              {yy_push_state(ASK_YN_MODE); return MP_as_clip_3ppolybase;}
<CL_MODE>"c3pp_min_signal_len" |
<CL_MODE>"c3ppmsl"              {return MP_as_clip_3ppolybase_len;}
<CL_MODE>"c3pp_max_errors_allowed" |
<CL_MODE>"c3ppmea"              {return MP_as_clip_3ppolybase_maxerrors;}
<CL_MODE>"c3pp_max_gap_from_ends" |
<CL_MODE>"c3ppmgfe"              {return MP_as_clip_3ppolybase_maxgap;}

<CL_MODE>"clip_polyat" |
<CL_MODE>"cpat"              {yy_push_state(ASK_YN_MODE); return MP_as_clip_polyat;}
<CL_MODE>"cp_keep_poly_stretch" |
<CL_MODE>"cpkps"              {yy_push_state(ASK_YN_MODE); return MP_as_clip_polyat_keeppolystretch;}
<CL_MODE>"cp_min_sequence_len" |
<CL_MODE>"cpmsl"              {return MP_as_clip_polyat_len;}
<CL_MODE>"cp_max_errors_allowed" |
<CL_MODE>"cpmea"              {return MP_as_clip_polyat_maxerrors;}
<CL_MODE>"cp_max_gap_from_ends" |
<CL_MODE>"cpmgfe"              {return MP_as_clip_polyat_maxgap;}

<CL_MODE>"pec_kmer_size" |
<CL_MODE>"peckms"             { return MP_as_clip_pec_basesperhash;}
<CL_MODE>"pec_bases_per_hash" |
<CL_MODE>"pecbph"                { return MP_ERROR_RENAMED_BPH_KMER;}
<CL_MODE>"propose_end_clip" |
<CL_MODE>"pec"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_proposeendclips;}
<CL_MODE>"pec_continuous" |
<CL_MODE>"pecc"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_continuous;}
<CL_MODE>"handle_solexa_ggcxg_problem" |
<CL_MODE>"pechsgp"            { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_sxaggcxg;}
<CL_MODE>"pffreq"            { return MP_as_clip_pec_ffreq;}
<CL_MODE>"pbfreq"            { return MP_as_clip_pec_bfreq;}
<CL_MODE>"pmkfr" |
<CL_MODE>"pec_minimum_kmer_forward_reverse"     { return MP_as_clip_pec_mkfr;}
<CL_MODE>"pmtk" |
<CL_MODE>"pec_minimum_total_kmer"     { return MP_as_clip_pec_mtk;}
<CL_MODE>"pffore"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_ffr;}
<CL_MODE>"pbfore"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_bfr;}
<CL_MODE>"pfcmst"            { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_fcmst;}
<CL_MODE>"pbcmst"            { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_bcmst;}
<CL_MODE>"pfsalp"            { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_fsalp;}
<CL_MODE>"pbsalp"            { yy_push_state(ASK_YN_MODE); return MP_as_clip_pec_bsalp;}
<CL_MODE>"gb_chimeradetectionclip" |
<CL_MODE>"gbcdc"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_sdbg_chimeradetection;}
<CL_MODE>"kmerjunk_detection" |
<CL_MODE>"kjd"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_kmer_junkdetection;}
<CL_MODE>"kmerjunk_completekill" |
<CL_MODE>"kjck"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_kmer_junkkill;}
<CL_MODE>"apply_skim_chimeradetectionclip" |
<CL_MODE>"ascdc"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_skimchimeradetection;}
<CL_MODE>"apply_skim_junkdetectionclip" |
<CL_MODE>"asjdc"              { yy_push_state(ASK_YN_MODE); return MP_as_clip_skimjunkdetection;}
<CL_MODE>"clip_bad_solexaends" |
<CL_MODE>"cbse"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_badsolexaends;}
<CL_MODE>"search_phix174" |
<CL_MODE>"spx174"               { yy_push_state(ASK_YN_MODE); return MP_as_search_phix174;}
<CL_MODE>"filter_phix174" |
<CL_MODE>"fpx174"               { yy_push_state(ASK_YN_MODE); return MP_as_filter_phix174;}
<CL_MODE>"filter_rrna_numkmers" |
<CL_MODE>"frrnank"              { return MP_as_filter_rrna_numkmers;}
<CL_MODE>"filter_rrna_pairs" |
<CL_MODE>"frrnap"               { yy_push_state(ASK_YN_MODE); return MP_as_filter_rrna_pairs;}
<CL_MODE>"filter_rrna" |
<CL_MODE>"frrna"                { yy_push_state(ASK_YN_MODE); return MP_as_filter_rrna;}
<CL_MODE>"clip_known_adaptorsright" |
<CL_MODE>"ckar"                { yy_push_state(ASK_YN_MODE); return MP_as_clip_knownadaptorsright;}
<CL_MODE>"rare_kmer_mask" |
<CL_MODE>"rkm"                { yy_push_state(ASK_YN_MODE); return MP_as_clipmask_rarekmers;}




<OUTPUT_MODE>"savesimplesingletsinproject"  |
<OUTPUT_MODE>"sssip"               { yy_push_state(ASK_YN_MODE); return MP_as_savesimplesingletsinproject;}
<OUTPUT_MODE>"savetaggedsingletsinproject"  |
<OUTPUT_MODE>"stsip"               { yy_push_state(ASK_YN_MODE); return MP_as_savesimplesingletsinproject;}
<OUTPUT_MODE>"remove_rollover_tmps"  |
<OUTPUT_MODE>"rrot"               { yy_push_state(ASK_YN_MODE); return MP_as_output_removerollovertmps;}
<OUTPUT_MODE>"remove_tmp_dir"  |
<OUTPUT_MODE>"rtd"               { yy_push_state(ASK_YN_MODE); return MP_as_output_removetmpdir;}

<OUTPUT_MODE>"output_tmpresult_html" |
<OUTPUT_MODE>"oth"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_html;}
<OUTPUT_MODE>"output_tmpresult_text" |
<OUTPUT_MODE>"ott"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_text;}
<OUTPUT_MODE>"output_tmpresult_tcs" |
<OUTPUT_MODE>"ots"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_tcs;}
<OUTPUT_MODE>"output_tmpresult_caf" |
<OUTPUT_MODE>"otc"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_caf;}
<OUTPUT_MODE>"output_tmpresult_maf" |
<OUTPUT_MODE>"otm"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_maf;}
<OUTPUT_MODE>"output_tmpresult_ace" |
<OUTPUT_MODE>"ota"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_ace;}
<OUTPUT_MODE>"output_tmpresult_gap4da" |
<OUTPUT_MODE>"otg"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_gap4da;}
<OUTPUT_MODE>"output_tmpresult_fasta" |
<OUTPUT_MODE>"otf"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_tmp_fasta;}

<OUTPUT_MODE>"output_exttmpresult_alsosinglets" |
<OUTPUT_MODE>"oetas"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_exttmp_alsosinglets;}
<OUTPUT_MODE>"output_exttmpresult_html" |
<OUTPUT_MODE>"oeth"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_exttmp_html;}
<OUTPUT_MODE>"output_exttmpresult_caf" |
<OUTPUT_MODE>"oetc"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_exttmp_caf;}
<OUTPUT_MODE>"output_exttmpresult_ace" |
<OUTPUT_MODE>"oeta"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_exttmp_ace;}
<OUTPUT_MODE>"output_exttmpresult_gap4da" |
<OUTPUT_MODE>"oetg"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_exttmp_gap4da;}
<OUTPUT_MODE>"output_exttmpresult_fasta" |
<OUTPUT_MODE>"oetf"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_exttmp_fasta;}

<OUTPUT_MODE>"output_result_caf" |
<OUTPUT_MODE>"orc"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_caf;}
<OUTPUT_MODE>"output_result_maf" |
<OUTPUT_MODE>"orm"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_maf;}
<OUTPUT_MODE>"output_result_html" |
<OUTPUT_MODE>"orh"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_html;}
<OUTPUT_MODE>"output_result_text" |
<OUTPUT_MODE>"ort"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_text;}
<OUTPUT_MODE>"output_result_ace" |
<OUTPUT_MODE>"ora"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_ace;}
<OUTPUT_MODE>"output_result_wiggle" |
<OUTPUT_MODE>"orw"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_wiggle;}
<OUTPUT_MODE>"output_result_tcs" |
<OUTPUT_MODE>"ors"                { yy_push_state(ASK_YN_MODE); return MP_as_output_tcs;}
<OUTPUT_MODE>"output_result_gff3" |
<OUTPUT_MODE>"org3"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_gff3;}
<OUTPUT_MODE>"output_result_gap4da" |
<OUTPUT_MODE>"org"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_gap4da;}
<OUTPUT_MODE>"output_result_fasta" |
<OUTPUT_MODE>"orf"                 { yy_push_state(ASK_YN_MODE); return MP_as_output_fasta;}

<OUTPUT_MODE>"hcpl" |
<OUTPUT_MODE>"html_chars_per_line" { return MP_con_output_html_cpl;}
<OUTPUT_MODE>"tcpl" |
<OUTPUT_MODE>"text_chars_per_line" { return MP_con_output_text_cpl;}
<OUTPUT_MODE>"hegfc" |
<OUTPUT_MODE>"html_endgap_fillchar" { singlechar_id=MP_con_output_html_gapfill; yy_push_state(GET_SINGLECHAR_STRING_MODE);}
<OUTPUT_MODE>"tegfc" |
<OUTPUT_MODE>"text_endgap_fillchar" { singlechar_id=MP_con_output_text_gapfill;yy_push_state(GET_SINGLECHAR_STRING_MODE);}

<GET_SINGLECHAR_STRING_MODE>[[:blank:]]*[=:]{1}[[:blank:]]*
<GET_SINGLECHAR_STRING_MODE>\'.\'          {yy_pop_state(); return singlechar_id;}
<GET_SINGLECHAR_STRING_MODE>\".\"          {yy_pop_state(); return singlechar_id;}
<GET_SINGLECHAR_STRING_MODE>[[:blank:]]
<GET_SINGLECHAR_STRING_MODE>.          {yy_pop_state(); return singlechar_id;}


<ASK_YN_MODE>[oO][nN]          {yy_pop_state(); return 1;}
<ASK_YN_MODE>[oO][fF][fF]      {yy_pop_state(); return 0;}
<ASK_YN_MODE>[tT][rR][uU][eE]     {yy_pop_state(); return 1;}
<ASK_YN_MODE>[fF][aA][lL][sS][eE] {yy_pop_state(); return 0;}
<ASK_YN_MODE>[yY][eE][sS]      {yy_pop_state(); return 1;}
<ASK_YN_MODE>[nN][oO]          {yy_pop_state(); return 0;}
<ASK_YN_MODE>[yYtT]      {yy_pop_state(); return 1;}
<ASK_YN_MODE>[nNfF]      {yy_pop_state(); return 0;}
<ASK_YN_MODE>[= \t]               {}
<ASK_YN_MODE>.               {yy_pop_state(); return MP_UNRECOGNISED_STRING;}

<DIR_MODE>"tmp_redirected_to" |
<DIR_MODE>"trt"               {filenameid=MP_dir_tmp_redirectedto; yy_push_state(FN_MODE); }
<DIR_MODE>"tmp"               {filenameid=MP_dir_tmp; yy_push_state(FN_MODE); }

<FN_MODE>{FILENAME}          { yy_pop_state(); return filenameid;}


<AS_MODE>"minimum_read_length" |
<AS_MODE>"mrl"                   { return MP_as_minimum_readlength;}
<AS_MODE>"minimum_reads_per_contig" |
<AS_MODE>"mrpc"                   { return MP_as_minimum_readspercontig;}
<AS_MODE>"num_of_passes" |
<AS_MODE>"nop"               { return(MP_as_numpasses);}
<AS_MODE>"kms" |
<AS_MODE>"kmer_series"       { yy_push_state(CSV_NUMBERS_MODE); return MP_as_kmerseries;}
<AS_MODE>"skim_each_pass" |
<AS_MODE>"sep"                 { return MP_ERROR_REMOVED;}
<AS_MODE>"rmb_break_loops" |
<AS_MODE>"rbl"                   { return MP_as_numrmbbreakloops;}
<AS_MODE>"max_contigs_per_pass" |
<AS_MODE>"mcpp"                 { return MP_as_maxcontigsperpass;}
<AS_MODE>"spoiler_detection" |
<AS_MODE>"sd"               { yy_push_state(ASK_YN_MODE); return MP_as_spoiler_detection;}
<AS_MODE>"sd_last_pass_only" |
<AS_MODE>"sdlpo"               { yy_push_state(ASK_YN_MODE); return MP_as_sd_lastpassonly;}
<AS_MODE>"use_run_length_encoding" |
<AS_MODE>"urle"                { yy_push_state(ASK_YN_MODE); return(MP_as_rle_reads);}
<AS_MODE>"automatic_repeat_detection" |
<AS_MODE>"ard"             { yy_push_state(ASK_YN_MODE); return MP_as_automatic_repeat_detection;}
<AS_MODE>"coverage_threshold" |
<AS_MODE>"ardct"             { return MP_as_ard_multicopythreshold;}
<AS_MODE>"min_length" |
<AS_MODE>"ardml"             { return MP_as_ard_multicopyminlen;}
<AS_MODE>"grace_length" |
<AS_MODE>"ardgl"             { return MP_as_ard_multicopygrace;}
<AS_MODE>"cutoff_multiplier" |
<AS_MODE>"urdcm"             { return MP_as_urd_cutoffmultiplier;}
<AS_MODE>"startinpass" |
<AS_MODE>"urdsip"             { return MP_as_urd_startinpass;}
<AS_MODE>"uniform_read_distribution" |
<AS_MODE>"urd"               { yy_push_state(ASK_YN_MODE); return MP_as_uniform_read_distribution;}
<AS_MODE>"use_genomic_pathfinder" |
<AS_MODE>"ugpf"               { yy_push_state(ASK_YN_MODE); return MP_paf_use_genomic_pathfinder;}
<AS_MODE>"use_emergency_search_stop" |
<AS_MODE>"uess"               { yy_push_state(ASK_YN_MODE); return MP_paf_use_emergency_search_stop;}
<AS_MODE>"use_emergency_blacklist" |
<AS_MODE>"uebl"               { yy_push_state(ASK_YN_MODE); return MP_paf_use_emergency_blacklist;}
<AS_MODE>"ess_partnerdepth" |
<AS_MODE>"esspd"              {return MP_paf_ess_partnerdepth;}
<AS_MODE>"use_max_contig_buildtime" |
<AS_MODE>"umcbt"              {yy_push_state(ASK_YN_MODE); return MP_paf_use_max_contig_buildtime;}
<AS_MODE>"buildtime_in_seconds" |
<AS_MODE>"bts"              {return MP_paf_buildtime_inseconds;}
<AS_MODE>"enforce_presence_of_qualities" |
<AS_MODE>"epoq"              {yy_push_state(ASK_YN_MODE); return MP_as_enforce_qualsinreads;}

<CSV_NUMBERS_MODE>[0-9, \t]+       { yy_pop_state(); return filenameid;}

<DP_MODE>"use_read_extension" |
<DP_MODE>"ure"                 { yy_push_state(ASK_YN_MODE); return MP_as_extend_reads;}
<DP_MODE>"read_extension_window_len" |
<DP_MODE>"rewl"              {return MP_as_readextension_window_len;}
<DP_MODE>"read_extension_window_maxerrors" |
<DP_MODE>"rewme"              {return MP_as_readextension_window_maxerrors;}
<DP_MODE>"first_extension_in_pass" |
<DP_MODE>"feip"              {return MP_as_readextension_firstpassnum;}
<DP_MODE>"last_extension_in_pass" |
<DP_MODE>"leip"              {return MP_as_readextension_lastpassnum;}


<EDIT_MODE>"gb_read_editing" |
<EDIT_MODE>"gbre"                 { yy_push_state(ASK_YN_MODE); return MP_ed_sdbg_readedit;}
<EDIT_MODE>"mira_automatic_contig_editing" |
<EDIT_MODE>"mace"                 { yy_push_state(ASK_YN_MODE); return MP_ed_mira_automatic_contig_editing;}
<EDIT_MODE>"edit_kmer_singlets"
<EDIT_MODE>"eks"                 { yy_push_state(ASK_YN_MODE); return MP_ed_kmer_singlets;}
<EDIT_MODE>"edit_homopolymer_overcalls"
<EDIT_MODE>"ehpo"                 { yy_push_state(ASK_YN_MODE); return MP_ed_homopolymer_overcalls;}

<EDIT_MODE>"edit_automatic_contig_editing" |
<EDIT_MODE>"eace"                 { yy_push_state(ASK_YN_MODE); return MP_ed_edit_automatic_contig_editing;}
<EDIT_MODE>"strict_editing_mode" |
<EDIT_MODE>"sem"                 { yy_push_state(ASK_YN_MODE); return MP_ed_strict;}
<EDIT_MODE>"confirmation_threshold" |
<EDIT_MODE>"ct"                   { return MP_ed_confirmation_threshold;}

<SKIM_MODE>"number_of_threads" |
<SKIM_MODE>"not"                   {return MP_sk_numthreads;}
<SKIM_MODE>"kmsmax"                   {return MP_sk_bph_max;}
<SKIM_MODE>"kmsaipp"                   {return MP_sk_bph_increasestep;}
<SKIM_MODE>"bases_per_hash" |
<SKIM_MODE>"bph"                    { return MP_ERROR_RENAMED_BPH_KMER;}
<SKIM_MODE>"hash_save_stepping" |
<SKIM_MODE>"hss"                   {return MP_ERROR_RENAMED_BPH_KMER;}
<SKIM_MODE>"kmer_size" |
<SKIM_MODE>"kms"                   {return MP_sk_basesperhash;}
<SKIM_MODE>"kmer_save_stepping" |
<SKIM_MODE>"kss"                   {return MP_sk_hashsavestepping;}
<SKIM_MODE>"percent_required" |
<SKIM_MODE>"pr"                   {return MP_sk_percentrequired;}
<SKIM_MODE>"maxhits_perread" |
<SKIM_MODE>"mhpr"                 {return MP_sk_maxhitsperread;}
<SKIM_MODE>"maxhashesinmemory" |
<SKIM_MODE>"mhim"                 {return MP_ERROR_RENAMED_BPH_KMER;}
<SKIM_MODE>"maxkmersinmemory" |
<SKIM_MODE>"mkim"                 {return MP_sk_maxhashesinmemory;}
<SKIM_MODE>"memcap_hitreduction" |
<SKIM_MODE>"mchr"                 {return MP_sk_memcaphitreduction;}
<SKIM_MODE>"also_compute_reverse_complements" |
<SKIM_MODE>"acrc"                {yy_push_state(ASK_YN_MODE); return MP_sk_alsoskimrevcomp;}
<SKIM_MODE>"sw_check_on_backbones" |
<SKIM_MODE>"swcob"                {yy_push_state(ASK_YN_MODE); return MP_sk_swcheckonbackbones;}
<SKIM_MODE>"filter_megahubs" |
<SKIM_MODE>"fmh"                {yy_push_state(ASK_YN_MODE); return MP_sk_filtermegahubs;}
<SKIM_MODE>"max_megahub_ratio" |
<SKIM_MODE>"mmhr"                 {return MP_sk_maxmegahubratio;}
<SKIM_MODE>"megahub_cap" |
<SKIM_MODE>"mhc"                 {return MP_sk_megahubcap;}

<SKIM_MODE>"nasty_repeat_ratio" |
<SKIM_MODE>"nrr"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"nasty_repeat_coverage" |
<SKIM_MODE>"nrc"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"repeatlevel_in_infofile" |
<SKIM_MODE>"rliif"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"mask_nasty_repeats" |
<SKIM_MODE>"mnr"                 { return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"freq_covestmin" |
<SKIM_MODE>"fcem"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"freq_est_minnormal" |
<SKIM_MODE>"fenn"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"freq_est_maxnormal" |
<SKIM_MODE>"fexn"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"freq_est_repeat" |
<SKIM_MODE>"fer"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"freq_est_heavyrepeat" |
<SKIM_MODE>"fehr"                 {return MP_ERROR_MOVED_SECTION_KS;}
<SKIM_MODE>"freq_est_crazyrepeat" |
<SKIM_MODE>"fecr"                 {return MP_ERROR_MOVED_SECTION_KS;}


<KMERSTAT_MODE>"nasty_repeat_ratio" |
<KMERSTAT_MODE>"nrr"                 {return MP_hs_nastyrepeatratio;}
<KMERSTAT_MODE>"nasty_repeat_coverage" |
<KMERSTAT_MODE>"nrc"                 {return MP_hs_nastyrepeatcoverage;}
<KMERSTAT_MODE>"repeatlevel_in_infofile" |
<KMERSTAT_MODE>"rliif"                 {return MP_hs_repeatlevel_in_infofile;}
<KMERSTAT_MODE>"mask_nasty_repeats" |
<KMERSTAT_MODE>"mnr"                 {yy_push_state(ASK_YN_MODE); return MP_hs_masknastyrepeats;}
<KMERSTAT_MODE>"freq_covestmin" |
<KMERSTAT_MODE>"fcem"                 {return MP_hs_freq_covestmin;}
<KMERSTAT_MODE>"freq_est_minnormal" |
<KMERSTAT_MODE>"fenn"                 {return MP_hs_freqest_minnormal;}
<KMERSTAT_MODE>"freq_est_maxnormal" |
<KMERSTAT_MODE>"fexn"                 {return MP_hs_freqest_maxnormal;}
<KMERSTAT_MODE>"freq_est_repeat" |
<KMERSTAT_MODE>"fer"                 {return MP_hs_freqest_repeat;}
<KMERSTAT_MODE>"freq_est_heavyrepeat" |
<KMERSTAT_MODE>"fehr"                 {return MP_hs_freqest_heavyrepeat;}
<KMERSTAT_MODE>"freq_est_crazyrepeat" |
<KMERSTAT_MODE>"fecr"                 {return MP_hs_freqest_crazyrepeat;}
<KMERSTAT_MODE>"lossless_digital_normalisation" |
<KMERSTAT_MODE>"ldn"                 {yy_push_state(ASK_YN_MODE); return MP_hs_applydigitalnormalisation;}
<KMERSTAT_MODE>"rare_kmer_final_kill" |
<KMERSTAT_MODE>"rkfk"                 {return MP_hs_rare_kmer_final_kill;}
<KMERSTAT_MODE>"memtouse" |
<KMERSTAT_MODE>"mtu"                  {return MP_hs_memtouse;}
<KMERSTAT_MODE>"million_hashes_per_buffer" |
<KMERSTAT_MODE>"mhpb"                 {return MP_ERROR_MOVED_SECTION_KS;}
<KMERSTAT_MODE>"million_kmers_per_buffer" |
<KMERSTAT_MODE>"mkpb"                 {return MP_ERROR_MOVED_SECTION_KS;}


<ALIGN_MODE>"bandwidth_in_percent" |
<ALIGN_MODE>"bip"                   {return MP_al_bip;}
<ALIGN_MODE>"bandwidth_max" |
<ALIGN_MODE>"bmax"                   {return MP_al_bmax;}
<ALIGN_MODE>"bandwidth_min" |
<ALIGN_MODE>"bmin"                   {return MP_al_bmin;}
<ALIGN_MODE>"min_score" |
<ALIGN_MODE>"ms"                   {return MP_al_min_score;}
<ALIGN_MODE>"min_overlap" |
<ALIGN_MODE>"mo"                   {return MP_al_min_overlap;}
<ALIGN_MODE>"min_relative_score" |
<ALIGN_MODE>"mrs"                  {return MP_al_min_relscore;}
<ALIGN_MODE>"egp_level" |
<ALIGN_MODE>"egpl"                 { yy_push_state(CSV_NUMBERS_MODE); return MP_ads_gp_functionstring;}
<ALIGN_MODE>"extra_gap_penalty" |
<ALIGN_MODE>"egp"                  { yy_push_state(ASK_YN_MODE); return MP_ads_extra_gap_penalty;}
<ALIGN_MODE>"max_egp_percent" |
<ALIGN_MODE>"megpp"                  {return MP_ads_max_gppercent;}
<ALIGN_MODE>"enforce_clean_ends" |
<ALIGN_MODE>"ece"                  { yy_push_state(ASK_YN_MODE); return MP_ads_enforce_clean_ends;}
<ALIGN_MODE>"clean_end_distance" |
<ALIGN_MODE>"ced"                  {return MP_ads_clean_end_distance;}
<ALIGN_MODE>"clean_end_mismatch_allowed" |
<ALIGN_MODE>"cema"                  {return MP_ads_clean_end_mismatchallowed;}

<FINALMAP_MODE>"active" |
<FINALMAP_MODE>"act"                   {yy_push_state(ASK_YN_MODE); return MP_fm_active;}
<FINALMAP_MODE>"map_perfect_matches" |
<FINALMAP_MODE>"mpm"                   {yy_push_state(ASK_YN_MODE); return MP_fm_mapperfect;}
<FINALMAP_MODE>"maxtotalerrors" |
<FINALMAP_MODE>"mte"                   {return MP_fm_maxtotalerrors;}
<FINALMAP_MODE>"maxmismatches" |
<FINALMAP_MODE>"mmm"                   {return MP_fm_maxmismatches;}
<FINALMAP_MODE>"maxgaps" |
<FINALMAP_MODE>"mg"                   {return MP_fm_maxgaps;}
<FINALMAP_MODE>"clean_end_distance" |
<FINALMAP_MODE>"ced"                   {return MP_fm_clean_end_dist;}

<CO_MODE>"analysis" |
<CO_MODE>"an"                      {BEGIN(CO_VALMODE); return MP_con_analyse_mode;}
<CO_VALMODE>"none"                 {BEGIN(CO_MODE); return 0;}
<CO_VALMODE>"text"                 {BEGIN(CO_MODE); return 1;}
<CO_VALMODE>"signal"               {BEGIN(CO_MODE); return 2;}
<CO_VALMODE>{ANID}              {BEGIN(CO_MODE); return MP_UNRECOGNISED_STRING;}
<CO_VALMODE>{FLOAT}             {BEGIN(CO_MODE); return MP_UNRECOGNISED_STRING;}
<CO_VALMODE>{INT}               {BEGIN(CO_MODE); return MP_UNRECOGNISED_STRING;}
<CO_VALMODE>\n               {BEGIN(CO_MODE); return MP_UNRECOGNISED_STRING;}
<CO_MODE>"rej_on_dropinrelscore" |
<CO_MODE>"rodirs"                  {return MP_con_rodirs;}
<CO_MODE>"cmin_relative_score" |
<CO_MODE>"cmrs"                  {return MP_con_min_relscore;}
<CO_MODE>"danger_max_error_rate" |
<CO_MODE>"dmer"                    {return MP_con_danger_max_error_rate;}
<CO_MODE>"mark_repeats" |
<CO_MODE>"mr"                 { yy_push_state(ASK_YN_MODE); return MP_as_mark_repeats;}
<CO_MODE>"only_in_result" |
<CO_MODE>"oir"            |
<CO_MODE>"mroir"                 { yy_push_state(ASK_YN_MODE); return MP_as_mark_repeats_only_in_result;}
<CO_MODE>"assume_snp_instead_repeat" |
<CO_MODE>"asir"                 { yy_push_state(ASK_YN_MODE); return MP_con_assume_snp_insteadof_rmb;}
<CO_MODE>"min_rmb_neighbour_qual" |
<CO_MODE>"min_rmb_neighbor_qual" |
<CO_MODE>"mnq"                    {return MP_con_min_rmb_neighbourqual;}
<CO_MODE>"min_reads_per_group" |
<CO_MODE>"mrpg"                    {return MP_con_min_readspergroup;}
<CO_MODE>"min_groupqual_for_rmb_tagging" |
<CO_MODE>"mgqrt"                    {return MP_con_min_groupqual_for_rmb_tagging;}
<CO_MODE>"min_coverage_percentage" |
<CO_MODE>"mcp"                    {return MP_con_min_coverage_percentage;}
<CO_MODE>"emea_set1_on_clipping_pec" |
<CO_MODE>"emeas1clpec"              {return MP_con_emea_setzero_on_clipping_pec;}
<CO_MODE>"endread_mark_exclusion_area" |
<CO_MODE>"emea"                    {return MP_con_endread_mark_exclusion_area;}
<CO_MODE>"also_mark_gap_bases_even_multicolumn" |
<CO_MODE>"amgbemc"                 { yy_push_state(ASK_YN_MODE); return MP_con_also_mark_gap_bases_evenmc;}
<CO_MODE>"also_mark_gap_bases_needbothstrands" |
<CO_MODE>"amgbnbs"                 { yy_push_state(ASK_YN_MODE); return MP_con_also_mark_gap_bases_needbothstrands;}
<CO_MODE>"also_mark_gap_bases" |
<CO_MODE>"amgb"                 { yy_push_state(ASK_YN_MODE); return MP_con_also_mark_gap_bases;}
<CO_MODE>"name_prefix" |
<CO_MODE>"name" |
<CO_MODE>"np"                    {filenameid=MP_con_name_prefix; yy_push_state(FN_MODE);}
<CO_MODE>"force_nonIUPACconsensus_perseqtype" |
<CO_MODE>"fnic"               { yy_push_state(ASK_YN_MODE); return MP_con_force_nonIUPACconsensus;}
<CO_MODE>"msr_maxerrors" |
<CO_MODE>"msrme"               { return MP_con_msr_maxerrors;}
<CO_MODE>"msr_keepcontigendsunmerged" |
<CO_MODE>"msrkceu"               { return MP_con_msr_keependsunmapped;}
<CO_MODE>"merge_short_reads" |
<CO_MODE>"msr"               { yy_push_state(ASK_YN_MODE); return MP_con_mergeshortreads;}
<CO_MODE>"gap_override_ratio" |
<CO_MODE>"gor"                    {return MP_con_gap_override_ratio;}

<MI_MODE>"extended_log" |
<MI_MODE>"el"                  { yy_push_state(ASK_YN_MODE); return MP_mi_extendedlog;}
<MI_MODE>"iknowwhatido" |
<MI_MODE>"ikwid"                  { yy_push_state(ASK_YN_MODE); return MP_mi_iknowwhatido;}
<MI_MODE>"large_contig_size_for_stats" |
<MI_MODE>"lcs4s"                  { return MP_mi_as_largecontigsize4stats;}
<MI_MODE>"large_contig_size" |
<MI_MODE>"lcs"                  { return MP_mi_as_largecontigsize;}
<MI_MODE>"ef1"                  { yy_push_state(ASK_YN_MODE); return MP_mi_extra_flag1;}
<MI_MODE>"ef2"                  { yy_push_state(ASK_YN_MODE); return MP_mi_extra_flag2;}
<MI_MODE>"ef3"                  { yy_push_state(ASK_YN_MODE); return MP_mi_extra_flag3;}
<MI_MODE>"sonfs"                { return MP_ERROR_MOVED_SECTION_NW;}
<MI_MODE>"somrnl"               { return MP_ERROR_MOVED_SECTION_NW;}

<NW_MODE>"check_nfs" |
<NW_MODE>"cnfs"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_nfs;}
<NW_MODE>"check_template_problems" |
<NW_MODE>"ctp"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_templateproblems;}
<NW_MODE>"check_duplicate_readnames" |
<NW_MODE>"cdrn"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_duplicatereadnames;}
<NW_MODE>"check_sra_readnames" |
<NW_MODE>"csrn"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_readnamesra;}
<NW_MODE>"check_maxreadnamelength" |
<NW_MODE>"cmrnl"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_maxreadnamelength;}
<NW_MODE>"maxreadnamelength" |
<NW_MODE>"mrnl"                  { return MP_nw_check_mrnlvalue;}
<NW_MODE>"check_multipassmapping" |
<NW_MODE>"cmpm"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_multipassmapping;}
<NW_MODE>"check_illuminajunkinest" |
<NW_MODE>"cijie"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_illuminajunkinest;}
<NW_MODE>"check_proposedendclip" |
<NW_MODE>"cpec"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_proposedendclip;}
<NW_MODE>"check_average_coverage" |
<NW_MODE>"cac"                  { yy_push_state(NW_CHOICEMODE); return MP_nw_check_coverage;}
<NW_MODE>"average_coverage_value" |
<NW_MODE>"acv"                  { return MP_nw_check_covvalue;}
<NW_MODE>"x174_contig_prefix" |
<NW_MODE>"x174cp"                  { filenameid=MP_nw_warn_x174prefix; yy_push_state(FN_MODE); }

<NW_CHOICEMODE>"no"          {yy_pop_state(); return 0;}
<NW_CHOICEMODE>"stop"          {yy_pop_state(); return 1;}
<NW_CHOICEMODE>"warn"          {yy_pop_state(); return 2;}
<NW_CHOICEMODE>[= \t]        {}
<NW_CHOICEMODE>{ANID}              { BEGIN(0); return MP_UNRECOGNISED_STRING;}
<NW_CHOICEMODE>{FLOAT}             { BEGIN(0); return MP_UNRECOGNISED_STRING;}
<NW_CHOICEMODE>{INT}               { BEGIN(0); return MP_UNRECOGNISED_STRING;}


<PAF_MODE>"skip_whole_contig_scan" |
<PAF_MODE>"swcs"                   { return MP_paf_skip_whole_contig_scan;}
<PAF_MODE>"use_quick_rule" |
<PAF_MODE>"uqr"                   { return MP_paf_use_quick_rule;}
<PAF_MODE>"quickrule_minlen1" |
<PAF_MODE>"qrml1"                   { return MP_paf_quickrule_minlen1;}
<PAF_MODE>"quickrule_minsim1" |
<PAF_MODE>"qrms1"                   { return MP_paf_quickrule_minsim1;}
<PAF_MODE>"quickrule_minlen2" |
<PAF_MODE>"qrml2"                   { return MP_paf_quickrule_minlen2;}
<PAF_MODE>"quickrule_minsim2" |
<PAF_MODE>"qrms2"                   { return MP_paf_quickrule_minsim2;}
<PAF_MODE>"backbonequickoverlap_minlen" |
<PAF_MODE>"bqoml"                   { return MP_paf_bbquickoverlap_minlen;}
<PAF_MODE>"max_startcache_filltime" |
<PAF_MODE>"mscft"              {return MP_paf_max_startcache_filltime;}


<JOB_MODE>draft      { return MP_jobdef_draft;}
<JOB_MODE>normal     { return MP_jobdef_normal   ;}
<JOB_MODE>accurate   { return MP_jobdef_accurate ;}
<JOB_MODE>genome     { return MP_jobdef_genome   ;}
<JOB_MODE>est        { return MP_jobdef_est      ;}
<JOB_MODE>fragments  { return MP_jobdef_fragments    ;}
<JOB_MODE>clustering { return MP_jobdef_clustering   ;}
<JOB_MODE>esps1      { return MP_jobdef_estsnppipeline1 ;}
<JOB_MODE>esps2      { return MP_jobdef_estsnppipeline2 ;}
<JOB_MODE>esps3      { return MP_jobdef_estsnppipeline3 ;}
<JOB_MODE>denovo     { return MP_jobdef_denovo   ;}
<JOB_MODE>mapping    { return MP_jobdef_mapping  ;}
<JOB_MODE>[Ss]anger     { return MP_jobdef_sanger   ;}
<JOB_MODE>454	     { return MP_jobdef_454	;}
<JOB_MODE>iontor |
<JOB_MODE>IonTor	     { return MP_jobdef_iontor	;}
<JOB_MODE>PacBioLQ |
<JOB_MODE>pacbiolq |
<JOB_MODE>PcBioLQ |
<JOB_MODE>pcbiolq     { return MP_jobdef_pacbiolq	;}
<JOB_MODE>PacBioHQ |
<JOB_MODE>pacbiohq |
<JOB_MODE>PcBioHQ |
<JOB_MODE>pcbiohq     { return MP_jobdef_pacbiohq	;}
<JOB_MODE>[Tt]ext     { return MP_jobdef_text   ;}
<JOB_MODE>[Ss]olexa     { return MP_jobdef_solexa   ;}
<JOB_MODE>[Ss]olid	     { return MP_jobdef_solid	;}
<JOB_MODE>[,:]
<JOB_MODE>{FILENAME} { return MP_UNRECOGNISED_STRING;}
<JOB_MODE>[\n\t ]    { BEGIN(0); return MP_jobdefend;}
<JOB_MODE>.  |
<JOB_MODE><<EOF>>    { BEGIN(0); return MP_jobdefend;}
<*>-{1,2}job[\t ]*=  { BEGIN(JOB_MODE); return MP_jobdefstart;}


<ASKFORTECH_NOCLIPPING>all      { return MP_quickmode_noclipping_all;}
<ASKFORTECH_NOCLIPPING>sanger   { return MP_quickmode_noclipping_sanger;}
<ASKFORTECH_NOCLIPPING>454      { return MP_quickmode_noclipping_454;}
<ASKFORTECH_NOCLIPPING>iontor   { return MP_quickmode_noclipping_iontor;}
<ASKFORTECH_NOCLIPPING>pcbiohq   { return MP_quickmode_noclipping_pacbiohq;}
<ASKFORTECH_NOCLIPPING>pcbiolq   { return MP_quickmode_noclipping_pacbiolq;}
<ASKFORTECH_NOCLIPPING>text     { return MP_quickmode_noclipping_text;}
<ASKFORTECH_NOCLIPPING>solexa   { return MP_quickmode_noclipping_solexa;}
<ASKFORTECH_NOCLIPPING>solid    { return MP_quickmode_noclipping_solid;}
<ASKFORTECH_NOCLIPPING>[ \t]    { BEGIN(0); }


<*>"-GENERAL" |
<*>"-GE"                       {BEGIN(GE_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-LOADREADS" |
<*>"-LR"                       {BEGIN(LR_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-STRAIN/BACKBONE" |
<*>"-SB"                       {BEGIN(SB_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-CLIPPING" |
<*>"-CL"                       {BEGIN(CL_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-OUTPUT" |
<*>"-OUT"    |
<*>"-OU"                       {BEGIN(OUTPUT_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-FILENAME" |
<*>"-FN"                          { BEGIN(FILE_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-ASSEMBLY" |
<*>"-AS"                           {BEGIN(AS_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-DATAPROCESSING" |
<*>"-DP"                           {BEGIN(DP_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-EDIT" |
<*>"-ED"                           {BEGIN(EDIT_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-SKIM" |
<*>"-SK"                           {BEGIN(SKIM_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-HASHSTATISTICS" |
<*>"-HASHSTAT" |
<*>"-HS"                           {return MP_ERROR_MOVED_SECTION_KS;}
<*>"-KMERSTATISTICS" |
<*>"-KMERSTAT" |
<*>"-KS"                           {BEGIN(KMERSTAT_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-ALIGN" |
<*>"-AL"                           {BEGIN(ALIGN_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-FINALMAP" |
<*>"-FM"                           {BEGIN(FINALMAP_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-CONTIG" |
<*>"-CO"                           {BEGIN(CO_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-DIRECTORY" |
<*>"-DIR" |
<*>"-DI"                           {BEGIN(DIR_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-MISC" |
<*>"-MI"                           {BEGIN(MI_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-NAGANDWARN" |
<*>"-NAG_AND_WARN" |
<*>"-NAW" |
<*>"-NW"                           {BEGIN(NW_MODE); return MP_PARSING_SECTIONCHANGE;}
<*>"-PATHFINDER" |
<*>"-PF"                         { BEGIN(PAF_MODE); return MP_PARSING_SECTIONCHANGE;}

<*>COMMON_SETTINGS    { return MP_params_for_common;}
<*>SANGER_SETTINGS    { return MP_params_for_sanger;}
<*>454_SETTINGS       { return MP_params_for_454;}
<*>IONTOR_SETTINGS    { return MP_params_for_iontor;}
<*>PCBIOLQ_SETTINGS    { return MP_params_for_pacbiolq;}
<*>PCBIOHQ_SETTINGS    { return MP_params_for_pacbiohq;}
<*>TEXT_SETTINGS    { return MP_params_for_text;}
<*>SOLEXA_SETTINGS    { return MP_params_for_solexa;}
<*>SOLID_SETTINGS     { return MP_params_for_solid;}

<*>_COMMON_SETTINGS    { return MP_silentparams_for_common;}
<*>_SANGER_SETTINGS    { return MP_silentparams_for_sanger;}
<*>_454_SETTINGS       { return MP_silentparams_for_454;}
<*>_IONTOR_SETTINGS    { return MP_silentparams_for_iontor;}
<*>_PCBIOLQ_SETTINGS    { return MP_silentparams_for_pacbiolq;}
<*>_PCBIOHQ_SETTINGS    { return MP_silentparams_for_pacbiohq;}
<*>_TEXT_SETTINGS    { return MP_silentparams_for_text;}
<*>_SOLEXA_SETTINGS    { return MP_silentparams_for_solexa;}
<*>_SOLID_SETTINGS     { return MP_silentparams_for_solid;}


<*>-{1,2}borg   { return MP_quickmode_borg;}

<*>-{1,2}parameterfile[\t ]*= |
<*>-{1,2}params[\t ]*= |
<*>-{1,2}parameters[\t ]*=   { filenameid=MP_quickmode_loadparam; yy_push_state(FN_MODE);}

<*>-{1,2}hirep_something  { return MP_quickmode_hirep_something;}
<*>-{1,2}hirep_good  { return MP_quickmode_hirep_good;}
<*>-{1,2}hirep_best  { return MP_quickmode_hirep_best;}
<*>-{1,2}highlyrepetitive  { return MP_quickmode_highlyrepetitive;}
<*>-{1,2}highqualitydata  { return MP_quickmode_highqualitydata;}
<*>-{1,2}lowqualitydata  { return MP_quickmode_lowqualitydata;}


<*>-{1,2}noclipping[s]{0,1}=     { BEGIN(ASKFORTECH_NOCLIPPING);}

<*>-{1,2}noclipping[s]{0,1}    { return MP_quickmode_noclipping_all;}

<*>-{1,2}noquality |
<*>-{1,2}noqualities    { return MP_quickmode_noquality_all;}

<*>"*=BEGIN0=*" {BEGIN(0); return MP_PARSING_SECTIONRESET;} /*Hack for
							     parsing command
							     lines
                                   where each blank resets the state */
<*>{ANID}  {return MP_ANID;}
<*>{FLOAT} {return MP_FLOAT;}
<*>{INT}   {return MP_INT;}

<*>":"
<*>"="
<*>";"

<*>#   {yy_push_state(COMMENT_MODE);} /* munch commentline, # is start of comment */
<COMMENT_MODE>[^\n]*\n  { yy_pop_state();}

<*>[\t\n ]   {/* munch these, if not recognised earlier */ }

<*>[-]{1,}[\t\n ] {return MP_ERROR_DASHES;}
<*>[-+] {return MP_ERROR;}

%%

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1110
#endif

/* deprecated

<*>-{1,2}horrid         { return MP_quickmode_horrid;}
<*>-{1,2}estmode        { return MP_quickmode_estmode;}
<*>-{1,2}genomedraft    { return MP_quickmode_genomedraft;}
<*>-{1,2}genomenormal   { return MP_quickmode_genomenormal;}
<*>-{1,2}genomeaccurate { return MP_quickmode_genomeaccurate;}
<*>-{1,2}mappingdraft   { return MP_quickmode_mappingdraft;}
<*>-{1,2}mappingnormal  { return MP_quickmode_mappingnormal;}
<*>-{1,2}mappingaccurate { return MP_quickmode_mappingaccurate;}

<GE_MODE>"projectinname" |
<GE_MODE>"projectin" |
<GE_MODE>"proin"                 { filenameid=MP_as_projectname_in; yy_push_state(FN_MODE);}
<GE_MODE>"projectoutname" |
<GE_MODE>"projectout" |
<GE_MODE>"proout"                 { filenameid=MP_as_projectname_out; yy_push_state(FN_MODE);}
<GE_MODE>"projectname" |
<GE_MODE>"project" |
<GE_MODE>"pro"                 { filenameid=MP_as_projectname; yy_push_state(FN_MODE);}

<*>-{1,2}proin[\t ]*= |
<*>-{1,2}projectin[\t ]*= |
<*>-{1,2}projectinname[\t ]*=  { filenameid=MP_as_projectname_in; yy_push_state(FN_MODE);}
<*>-{1,2}proout[\t ]*= |
<*>-{1,2}projectout[\t ]*= |
<*>-{1,2}projectoutname[\t ]*=  { filenameid=MP_as_projectname_out; yy_push_state(FN_MODE);}
<*>-{1,2}project[\t ]*= |
<*>-{1,2}projectname[\t ]*=  { filenameid=MP_as_projectname; yy_push_state(FN_MODE);}
 */
