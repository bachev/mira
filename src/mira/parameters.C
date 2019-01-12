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


#include <mira/parameters.H>

#include <boost/algorithm/string.hpp>   // trim, split etc.
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "mira/parameters_tokens.h"
#include "mira/read.H"
#include "util/machineinfo.H"
#include "util/fileanddisk.H"
#include "util/fmttext.H"


using std::cout;
using std::cerr;
using std::endl;

//#define CEBUG(bla) {cout << bla; cout.flush();}
#define CEBUG(bla)



bool   MIRAParameters::MP_errorinparams=false;
std::string MIRAParameters::MP_currentparametersection;
std::vector<std::string> MIRAParameters::MP_loadfilename;

std::string MIRAParameters::MP_binpath;
std::string MIRAParameters::MP_bindir;
std::string MIRAParameters::MP_homedir;
std::string MIRAParameters::MP_sharedir;
std::string MIRAParameters::MP_mhslibdir;



MIRAParameters::MIRAParameters()
{
  MP_errorinparams=false;
  MP_parsedsomesettings=false;
  MP_jobtechsettings=false;
  MP_loadfilename.clear();
  MP_loadfilename.reserve(12);


  // force the MIRAParameters object to check for presence of technologies
  //  in the manifest readgroups if a technology-dependend parameter was used
  mp_special_params.sp_parse_checktechnologypresence=true;

  // start step 0 is not a valid start step for miraSearchESTSNPs
  // but this at least catches errors if they occur
  mp_special_params.sp_est_startstep=0;
  mp_special_params.mi_extended_log=false;
  mp_special_params.mi_iknowwhatido=false;
  mp_special_params.mi_as_largecontigsize=500;
  mp_special_params.mi_as_largecontigsize4stats=5000;
  mp_special_params.mi_extra_flag1=false;
  mp_special_params.mi_extra_flag2=true;
  mp_special_params.mi_extra_flag3=false;

  mp_nagandwarn_params.nw_check_maxreadnamelength=NWSTOP;
  mp_nagandwarn_params.nw_check_mrnlvalue=40;
  mp_nagandwarn_params.nw_check_nfs=NWSTOP;
  mp_nagandwarn_params.nw_check_templateproblems=NWSTOP;
  mp_nagandwarn_params.nw_check_duplicatereadnames=NWSTOP;
  mp_nagandwarn_params.nw_check_readnamesra=NWSTOP;
  mp_nagandwarn_params.nw_check_illuminajunkinest=NWSTOP;
  mp_nagandwarn_params.nw_check_proposedendclipping=NWSTOP;
  mp_nagandwarn_params.nw_check_multipassmapping=NWSTOP;
  mp_nagandwarn_params.nw_check_coverage=NWSTOP;
  mp_nagandwarn_params.nw_check_covvalue=80;

  mp_assembly_params.as_numthreads=2;
  mp_assembly_params.as_automemmanagement=true;
  mp_assembly_params.as_amm_keeppercentfree=15;   // use all system mem minus 15%
  mp_assembly_params.as_amm_maxprocesssize=0;  // 0 = unlimited, use keep percent free

  mp_skim_params.sk_numthreads=mp_assembly_params.as_numthreads;
  mp_skim_params.sk_basesperhash=17;
  mp_skim_params.sk_bph_max=0;
  mp_skim_params.sk_bph_increasestep=0;
  mp_skim_params.sk_hashsavestepping=4;
  mp_skim_params.sk_percentrequired=50;
  mp_skim_params.sk_maxhitsperread=2000;
  mp_skim_params.sk_maxhashesinmem=15000000;
  mp_skim_params.sk_memcaphitreduction=1024;
  mp_skim_params.sk_alsoskimrevcomp=true;
  mp_skim_params.sk_swcheckonbackbones=true;
  mp_skim_params.sk_filtermegahubs=true;
  mp_skim_params.sk_maxmegahubratio=0.0;
  mp_skim_params.sk_megahubcap=150000;

  mp_hashstatistics_params.hs_freqest_minnormal=0.4;
  mp_hashstatistics_params.hs_freqest_maxnormal=1.6;
  mp_hashstatistics_params.hs_freqest_repeat=2.1;
  mp_hashstatistics_params.hs_freqest_heavyrepeat=8;
  mp_hashstatistics_params.hs_freqest_crazyrepeat=20;
  mp_hashstatistics_params.hs_freq_covestmin=0;
  mp_hashstatistics_params.hs_masknastyrepeats=true;
  mp_hashstatistics_params.hs_nastyrepeatratio=100;
  mp_hashstatistics_params.hs_nastyrepeatcoverage=0;
  mp_hashstatistics_params.hs_repeatlevel_in_infofile=6;
  mp_hashstatistics_params.hs_apply_digitalnormalisation=false;
  mp_hashstatistics_params.hs_rare_kmer_final_kill=0;
  mp_hashstatistics_params.hs_memtouse=75;

  mp_pathfinder_params.paf_use_genomic_algorithms=true;
  mp_pathfinder_params.paf_use_emergency_blacklist=true;
  mp_pathfinder_params.paf_use_emergency_search_stop=true;
  mp_pathfinder_params.paf_use_max_contig_buildtime=false;
  mp_pathfinder_params.paf_ess_depth=500;
  mp_pathfinder_params.paf_max_contig_buildtime=10000;
  mp_pathfinder_params.paf_max_startcache_filltime=6;
  mp_pathfinder_params.paf_use_quick_rule=true;
  mp_pathfinder_params.paf_quickrule_minlen1=200;
  mp_pathfinder_params.paf_quickrule_minsim1=90;
  mp_pathfinder_params.paf_quickrule_minlen2=100;
  mp_pathfinder_params.paf_quickrule_minsim2=95;
  mp_pathfinder_params.paf_bbquickoverlap_minlen=150;

  mp_pathfinder_params.paf_skipwholecontigscan=5;

  setPathfinderMaxContigTime(mp_pathfinder_params.paf_max_contig_buildtime);
  setPathfinderNextReadMaxTimeToESS(1);

  mp_edit_params.ed_sdbg_readedit=true;

  mp_edit_params.ed_mira_automatic_contic_editing=true;
  mp_edit_params.ed_kmer_singlets=true;
  mp_edit_params.ed_homopolymer_overcalls=false;

  // Thomas' editor, currently unused
  mp_edit_params.ed_edit_automatic_contic_editing=false;
  mp_edit_params.ed_strict_editing_mode=false;
  mp_edit_params.ed_confirmation_threshold=50;

  mp_assembly_params.as_dateoutput=true;

  mp_assembly_params.as_minimum_readlength=80;
  mp_assembly_params.as_minimum_readspercontig=1;

  mp_assembly_params.as_buntify_reads=true;
  mp_assembly_params.as_rle_reads=false;

  mp_assembly_params.as_backbone_rail_fromstrain="";
  mp_assembly_params.as_backbone_strainname_forceforall=false;
  mp_assembly_params.as_backbone_raillength=0;
  mp_assembly_params.as_backbone_railoverlap=0;
  mp_assembly_params.as_backbone_outlen=30000;
  mp_assembly_params.as_backbone_trimoverhangingreads=true;
  mp_assembly_params.as_backbone_bootstrapnewbackbone=true;

  mp_assembly_params.as_assemblyjob_accurate=true;
  mp_assembly_params.as_assemblyjob_mapping=false;
  mp_assembly_params.as_assemblyjob_preprocessonly=false;

  mp_assembly_params.as_enforce_qualsinreads=true;
  mp_assembly_params.as_wants_qualityfile=true;

  mp_assembly_params.as_cleanup_tmpfiles=0;
  mp_assembly_params.as_spoilerdetection=true;
  mp_assembly_params.as_spdetect_lastpassonly=true;
  mp_assembly_params.as_automatic_repeat_detection=true;
  mp_assembly_params.as_ard_multicopythreshold=2.0;
  mp_assembly_params.as_ard_multicopyminlen=400;
  mp_assembly_params.as_ard_multicopygrace=40;
  mp_assembly_params.as_uniform_read_distribution=true;
  mp_assembly_params.as_urd_startinpass=3;
  mp_assembly_params.as_urd_cutoffmultiplier=1.5;
  mp_assembly_params.as_numpasses=0;
  mp_assembly_params.as_maxcontigsperpass=0;
  mp_assembly_params.as_numrmbbreakloops=2;
  mp_assembly_params.as_startbackboneusage_inpass=3;
  mp_assembly_params.as_use_read_extension=true;
  mp_assembly_params.as_readextension_window_len=30;
  mp_assembly_params.as_readextension_window_maxerrors=2;
  mp_assembly_params.as_readextension_firstpassnum=0;
  mp_assembly_params.as_readextension_lastpassnum=0;
  //mp_assembly_params.as_skimeachpass=true;
  mp_assembly_params.as_mark_repeats=true;
  mp_assembly_params.as_mark_repeats_onlyinresult=false;

  mp_assembly_params.as_projectname_in="mira";
  mp_assembly_params.as_projectname_out=mp_assembly_params.as_projectname_in;

  mp_assembly_params.as_put_asswithmira_tags=true;

  mp_assembly_params.as_filecheck_only=false;

  mp_assembly_params.as_savesimplesingletsinproject=false;
  mp_assembly_params.as_savetaggedsingletsinproject=true;

  mp_assembly_params.as_output_caf=true;
  mp_assembly_params.as_output_maf=true;
  mp_assembly_params.as_output_fasta=true;
  mp_assembly_params.as_output_gap4da=false;
  mp_assembly_params.as_output_ace=false;
  mp_assembly_params.as_output_gff3=false;
  mp_assembly_params.as_output_wiggle=true;
  mp_assembly_params.as_output_html=false;
  mp_assembly_params.as_output_tcs=true;
  mp_assembly_params.as_output_txt=true;

  mp_assembly_params.as_output_tmp_html=false;
  mp_assembly_params.as_output_tmp_tcs=false;
  mp_assembly_params.as_output_tmp_txt=false;
  mp_assembly_params.as_output_tmp_ace=false;
  mp_assembly_params.as_output_tmp_fasta=false;
  mp_assembly_params.as_output_tmp_gap4da=false;
  mp_assembly_params.as_output_tmp_caf=true;
  mp_assembly_params.as_output_tmp_maf=false;

  mp_assembly_params.as_output_exttmp_html=false;
  mp_assembly_params.as_output_exttmp_ace=false;
  mp_assembly_params.as_output_exttmp_fasta=false;
  mp_assembly_params.as_output_exttmp_gap4da=false;
  mp_assembly_params.as_output_exttmp_caf=false;
  mp_assembly_params.as_output_exttmp_alsosinglets=false;

  mp_assembly_params.as_output_removerollovertmps=true;
  mp_assembly_params.as_output_removetmpdir=false;

  mp_assembly_params.as_clip_sdbg_chimeradetection=true;
  mp_assembly_params.as_clip_kmer_junkdetection=true;
  mp_assembly_params.as_clip_kmer_junkkill=true;

  mp_assembly_params.as_clip_skimchimeradetection=false;
  mp_assembly_params.as_clip_skimjunkdetection=false;
  mp_assembly_params.as_clip_proposeendclips=true;
  mp_assembly_params.as_clip_pec_sxaggcxg=true;
  mp_assembly_params.as_clip_pec_continuous=true;
  mp_assembly_params.as_clip_pec_ffr=false;
  mp_assembly_params.as_clip_pec_bfr=false;
  mp_assembly_params.as_clip_pec_fcmst=false;
  mp_assembly_params.as_clip_pec_bcmst=false;
  mp_assembly_params.as_clip_pec_fsalp=false;
  mp_assembly_params.as_clip_pec_bsalp=false;
  mp_assembly_params.as_clip_pec_basesperhash=17;
  mp_assembly_params.as_clip_pec_ffreq=0;
  mp_assembly_params.as_clip_pec_bfreq=0;
  mp_assembly_params.as_clip_pec_mkfr=1;
  mp_assembly_params.as_clip_pec_mtk=3;

  mp_assembly_params.as_clip_badsolexaends=false;
  mp_assembly_params.as_clip_knownadaptorsright=false;
  mp_assembly_params.as_search_phix174=false;
  mp_assembly_params.as_filter_phix174=false;
  mp_assembly_params.as_filter_rrna=false;
  mp_assembly_params.as_filter_rrna_pairs=false;
  mp_assembly_params.as_filter_rrna_numkmers=17;

  mp_assembly_params.as_clip_lowercase_front=false;
  mp_assembly_params.as_clip_lowercase_back=false;

  mp_assembly_params.as_clip_polyat=false;
  mp_assembly_params.as_clip_polyat_keeppolystretch=false;
  mp_assembly_params.as_clip_polyat_len=12;
  mp_assembly_params.as_clip_polyat_maxerrors=1;
  mp_assembly_params.as_clip_polyat_maxgap=9;

  mp_assembly_params.as_clip_3ppolybase=false;
  mp_assembly_params.as_clip_3ppolybase_len=15;
  mp_assembly_params.as_clip_3ppolybase_maxerrors=3;
  mp_assembly_params.as_clip_3ppolybase_maxgap=9;

  mp_assembly_params.as_clip_possible_vectors=true;
  mp_assembly_params.as_clip_quality=true;
  mp_assembly_params.as_clip_badstretchquality=true;
  mp_assembly_params.as_clip_maskedbases=true;

  mp_assembly_params.as_clip_quality_minthreshold=0;
  mp_assembly_params.as_clip_quality_numminthreshold=0;

  mp_assembly_params.as_clip_quality_minqual=20;
  mp_assembly_params.as_clip_quality_winlen=30;
  mp_assembly_params.as_clip_badstretchquality_minqual=5;
  mp_assembly_params.as_clip_badstretchquality_winlen=20;
  mp_assembly_params.as_clip_vector_maxlenallowed=18;
  mp_assembly_params.as_clip_maskedbase_gapsize=20;
  mp_assembly_params.as_clip_maskedbase_maxfrontgap=40;
  mp_assembly_params.as_clip_maskedbase_maxendgap=60;
  mp_assembly_params.as_clip_ssahamerge_gapsize=10;
  mp_assembly_params.as_clip_ssahamerge_maxfrontgap=60;
  mp_assembly_params.as_clip_ssahamerge_maxendgap=120;
  mp_assembly_params.as_clip_ssahamerge_strictfrontclip=0;
  mp_assembly_params.as_clip_ssahamerge_strictendclip=0;
  mp_assembly_params.as_clip_ensureminimumleftclipoff=true;
  mp_assembly_params.as_clip_ensureminimumrightclipoff=false;
  mp_assembly_params.as_clip_minslrequired=25;
  mp_assembly_params.as_clip_minqlsetto=30;
  mp_assembly_params.as_clip_minsrrequired=10;
  mp_assembly_params.as_clip_minqrsetto=20;

  mp_assembly_params.as_clipmask_rarekmers=false;

// Initialise the scoring scheme for the similarity matrix
  mp_dynamic_params.dyn_score_multiplier=10;
  mp_dynamic_params.dyn_score_match=1*mp_dynamic_params.dyn_score_multiplier;
  mp_dynamic_params.dyn_score_mismatch=-1*mp_dynamic_params.dyn_score_multiplier;
  mp_dynamic_params.dyn_score_halfmatch=static_cast<int32>(.5*mp_dynamic_params.dyn_score_multiplier);
  mp_dynamic_params.dyn_score_halfmismatch=static_cast<int32>(-.5*mp_dynamic_params.dyn_score_multiplier);
  mp_dynamic_params.dyn_score_nmatch=0;
  mp_dynamic_params.dyn_score_npenaltymatch=-1;
  mp_dynamic_params.dyn_score_gap=-2*mp_dynamic_params.dyn_score_multiplier;
  mp_dynamic_params.dyn_score_oldgap=static_cast<int32>((-1.5)*mp_dynamic_params.dyn_score_multiplier);
  mp_dynamic_params.dyn_score_oldgapmatch=static_cast<int32>((-1.1)*mp_dynamic_params.dyn_score_multiplier);
  //mp_dynamic_params.dyn_score_gap=-5;
  //mp_dynamic_params.dyn_score_oldgap=-4;
  mp_dynamic_params.dyn_score_ltermgap=0;      // DO NOT TOUCH
  mp_dynamic_params.dyn_score_rtermgap=0;      // DO NOT TOUCH
  mp_dynamic_params.dyn_score_termgap=std::min(mp_dynamic_params.dyn_score_ltermgap,
					  mp_dynamic_params.dyn_score_rtermgap);       // 0  DO NOT CHANGE THIS ONE

  mp_align_params.al_max_cutoff=1;
  mp_align_params.al_min_score=15*mp_dynamic_params.dyn_score_multiplier;
  mp_align_params.al_min_relscore=65;
  mp_align_params.al_min_overlap=17;
  mp_align_params.al_kmin=25;
  mp_align_params.al_kmax=100;
  mp_align_params.al_kpercent=15;

  mp_finalmap_params.fm_active=true;
  mp_finalmap_params.fm_mapperfect=true;
  mp_finalmap_params.fm_maxtotalerrors=-1;
  mp_finalmap_params.fm_maxmismatches=-1;
  mp_finalmap_params.fm_maxgaps=-1;
  mp_finalmap_params.fm_clean_end_dist=0;

  mp_align_params.ads_extra_gap_penalty=false;
  setAlignGapPenaltyLevel(std::string("0,5,10,20,40,80,100"));
  mp_align_params.ads_max_gppercent=100;
  mp_align_params.ads_enforce_clean_ends=false;
  mp_align_params.ads_clean_end_distance=0;
  mp_align_params.ads_clean_end_mismatchallowed=0;
  //mp_align_params.ads_extra_mismatch_penalty=true;
  //mp_align_params.ads_emp_windowlen=30;
  //mp_align_params.ads_emp_maxmismatches=15;

  mp_contig_params.con_reject_on_drop_in_relscore=15;
  mp_contig_params.con_min_relscore=-1;
  mp_contig_params.con_danger_analyse_mode=ANALYSE_SIGNAL;
  mp_contig_params.con_danger_max_error_rate=1;         // 1% error rate

  mp_contig_params.con_assume_snp_insteadof_rmb=false;
  mp_contig_params.con_disregard_spurious_rmb_mismatches=true;
  mp_contig_params.con_also_mark_gap_bases=true;
  mp_contig_params.con_also_mark_gap_bases_needbothstrands=true;
  mp_contig_params.con_also_mark_gap_bases_evenmc=true;
  mp_contig_params.con_shorttagcomments=true;
  mp_contig_params.con_minrmbneighbourqual=20;
  mp_contig_params.con_mingroupqualforrmbtagging=30;
  mp_contig_params.con_mincoveragepercentage=0;
  //mp_contig_params.con_min_groupqual_srmbwrmb_change=45;
  //mp_contig_params.con_rmb_numzone_trigger=1;
  mp_contig_params.con_endreadmarkexclusionarea=25;
  mp_contig_params.con_emea_setzero_on_clipping_pec=true;
  mp_contig_params.con_minreadspergroup=2;
  mp_contig_params.con_gap_override_ratio=66;

  mp_contig_params.con_force_nonIUPACconsensus=false;
  mp_contig_params.con_force_nonIUPACconsensus_perseqtype=false;
  mp_contig_params.con_mergeshortreads=false;
  mp_contig_params.con_msr_keependsunmapped=-1;
  mp_contig_params.con_msr_maxerrors=0;

  mp_contig_params.con_output_text_cpl=60;
  mp_contig_params.con_output_html_cpl=60;
  mp_contig_params.con_output_text_gapfill=' ';
  mp_contig_params.con_output_html_gapfill=' ';

  mp_directory_params.dir_tmp_redirectedto="";
}


MIRAParameters::~MIRAParameters()
{
}



void MIRAParameters::setupStdMIRAParameters(std::vector<MIRAParameters> & Pv, bool verbose)
{
  Pv.clear();
  Pv.resize(ReadGroupLib::SEQTYPE_END);

  generateProjectNames(Pv);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  {
    auto testPv=Pv;
    const char * ms_dn_f =
      "--job=denovo,fragments,accurate,sanger,454,iontor,pcbiolq,pcbiohq,text,solexa,solid"
      ;
    const char * ms_dn_c =
      "--job=denovo,clustering,accurate,sanger,454,iontor,pcbiolq,pcbiohq,text,solexa,solid"
      ;
    const char * ms_dn_g =
      "--job=denovo,genome,accurate,sanger,454,iontor,pcbiolq,pcbiohq,text,solexa,solid"
      ;
    const char ms_dn_e[] =
      "--job=denovo,est,accurate,sanger,454,iontor,pcbiolq,pcbiohq,text,solexa,solid"
      ;
    const char ms_m_g[] =
      "--job=mapping,genome,accurate,sanger,454,iontor,pcbiolq,pcbiohq,text,solexa,solid"
      ;
    const char ms_m_e[] =
      "--job=mapping,est,accurate,sanger,454,iontor,pcbiolq,pcbiohq,text,solexa,solid"
      ;

    const char * apptr=ms_dn_f;
    try {
      parseQuickmodeNoTechSettingsChange(apptr,
					 "Default settings",
					 testPv, verbose);
      apptr=ms_dn_c;
      parseQuickmodeNoTechSettingsChange(apptr,
					 "Default settings",
					 testPv, verbose);
      apptr=ms_dn_g;
      parseQuickmodeNoTechSettingsChange(apptr,
					 "Default settings",
					 testPv, verbose);
      apptr=ms_dn_e;
      parseQuickmodeNoTechSettingsChange(apptr,
					 "Default settings",
					 testPv, verbose);
      apptr=ms_m_g;
      parseQuickmodeNoTechSettingsChange(apptr,
					 "Default settings",
					 testPv, verbose);
      apptr=ms_m_e;
      parseQuickmodeNoTechSettingsChange(apptr,
					 "Default settings",
					 testPv, verbose);
    }
    catch (Notify n) {
      std::ostringstream my__emsg;
      my__emsg << "Error while parsing MIRA-internal standard parameters:\n"
	       << apptr
	       << "\n\nThere is nothing you can do about it, this is a blunder by the author."
	       << "\nPlease file a bug report immediately.";
      n.setMsg(my__emsg.str().c_str());
      n.setGravity(Notify::INTERNAL);
      n.handleError(n.tif);
    }
  }

  //for(auto & pve : Pv){
  //  pve.MP_parsedsomesettings=false;
  //}

}


// adapts a few values automatically after parsing
void MIRAParameters::postParsingChanges(std::vector<MIRAParameters> & Pv)
{
  FUNCSTART("void MIRAParameters::postParsingChanges(std::vector<MIRAParameters> & Pv)");
  BUGIFTHROW(Pv.empty(), "Empty MIRAParameters vector???");

  assembly_parameters const & as_fixparams= Pv[0].mp_assembly_params;

  if(as_fixparams.as_clip_proposeendclips
     && Pv[0].mp_contig_params.con_emea_setzero_on_clipping_pec){
    cout << "-CL:pec and -CO:emeas1clpec are set, setting -CO:emea values to 1.\n";
    for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++){
      Pv[i].mp_contig_params.con_endreadmarkexclusionarea=1;
    }
  }

  if(Pv[0].mp_nagandwarn_params.nw_check_multipassmapping!=NWNONE
     && as_fixparams.as_assemblyjob_mapping
     && (as_fixparams.as_numpasses>1
	 || as_fixparams.as_numrmbbreakloops>1)){
    std::string wmsg="You are running a mapping job with:";
    if(as_fixparams.as_numpasses>1){
      wmsg+="\n-AS:nop="+boost::lexical_cast<std::string>(as_fixparams.as_numpasses);
    }
    if(as_fixparams.as_numrmbbreakloops>1){
      wmsg+="\n-AS:rbl="+boost::lexical_cast<std::string>(as_fixparams.as_numrmbbreakloops);
    }
    wmsg+="\n\nA mapping assembly with the above parameter(s) >1 will probably not"
      " produce what you expect. For haploid organisms or for simple SNP search in"
      " multiploid organisms, maybe you should rethink your choice."
      "\nHowever, certain use cases do benefit from the parameters you have chosen."
      " E.g., doing clean mappings of multiploid organisms. If you think"
      " you want this, please use the following parameter to let MIRA continue:"
      "\n -NW:cmpm=warn"
      "\nor"
      "\n -NW:cmpm=no"
      ;
    if(Pv[0].mp_nagandwarn_params.nw_check_multipassmapping==NWSTOP){
      MIRANOTIFY(Notify::FATAL,wmsg);
    }else{
      cout << "\nWARNING!\n" << FmtText::wordWrap(wmsg,80) << endl;
    }
  }

  correctTmpDirectory(Pv);

  FUNCEND();
}

// the "intelligent" version of dumpAllParams()
// looks which data are to be loaded or where settings have been parsed for
//  and sets the vector with the needed data accordingly
void MIRAParameters::dumpAllParams(std::vector<MIRAParameters> & Pv, std::ostream & ostr)
{
  FUNCSTART("void MIRAParameters::dumpAllParams(std::vector<MIRAParameters> & Pv, std::ostream & ostr)");
  BUGIFTHROW(Pv.empty(), "Empty MIRAParameters vector???");

  std::vector<int> indexesInPv;
  //indexesInPv.push_back(0);

  for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++){
    if(Pv[i].MP_jobtechsettings){
      indexesInPv.push_back(i);
    }
  }

  // Uh ... can't be empty or we'll get a SEGFAULT. Simply print out all in that case.
  if(indexesInPv.empty()){
    for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++){
      indexesInPv.push_back(i);
    }
  }

  // What a cludge:
  // multiParamPrint() uses size of vector to decide whether it is a technology dependent parameter or not
  // Which, for projects with just one technology, would not print out the "[tech]" qualifier. Not good.
  // Therefore, add the SEQTYPE_END qualifier (and make sure it's caught later on)
  // TODO: rethink this whole thing.
  indexesInPv.push_back(ReadGroupLib::SEQTYPE_END);

  dumpAllParams(Pv, indexesInPv, ostr);

  FUNCEND();
}

void MIRAParameters::dumpAllParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  FUNCSTART("void MIRAParameters::dumpAllParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)");

  BUGIFTHROW(indexesInPv.empty(), "Trying to dump nothing?\n");

  ostr << "------------------------------------------------------------------------------\nParameter settings seen for:\n";

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  for(uint32 i=0; i<indexesInPv.size(); i++){
    if(indexesInPv[i] != SEQTYPE_END){
      if(i>0) ostr << ", ";
      ostr << ReadGroupLib::getNameOfSequencingType(i) << " data";
    }
  }

  ostr << "\n\nUsed parameter settings:\n";
  dumpGeneralParams(Pv, indexesInPv, ostr);
  dumpLoadParams(Pv, indexesInPv, ostr);
  dumpAssemblyParams(Pv, indexesInPv, ostr);
  dumpStrainBackboneParams(Pv, indexesInPv, ostr);
  dumpDataProcessingParams(Pv, indexesInPv, ostr);
  dumpClippingParams(Pv, indexesInPv, ostr);
  dumpSkimParams(Pv, indexesInPv, ostr);
  dumpHashStatisticsParams(Pv, indexesInPv, ostr);
  dumpPathfinderParams(Pv, indexesInPv, ostr);
  dumpAlignParams(Pv, indexesInPv, ostr);
  dumpFinalMappingParams(Pv, indexesInPv, ostr);
  dumpContigParams(Pv, indexesInPv, ostr);
  dumpEditParams(Pv, indexesInPv, ostr);
  dumpMiscParams(Pv, indexesInPv, ostr);
  dumpNagAndWarnParams(Pv, indexesInPv, ostr);
  dumpDirectoryParams(Pv, indexesInPv, ostr);
  dumpFileInParams(Pv, indexesInPv, ostr);
  dumpFileOutParams(Pv, indexesInPv, ostr);
  dumpFileTempParams(Pv, indexesInPv, ostr);
  dumpOutputCustomisationParams(Pv, indexesInPv, ostr);
  dumpFileDirectoryOutNamesParams(Pv, indexesInPv, ostr);

  ostr << "------------------------------------------------------------------------------\n";

  FUNCEND();
}


#define OUTSTRING(o,xx) {if(xx.size()!=0){o<<(xx)<<"\n";}else{o<<"(none)\n";}}

// static func
void MIRAParameters::dumpFileDirectoryOutNamesParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  ostr << "\n    File / directory output names:\n";
  ostr << "\tCAF             : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_CAF+".caf"));
  ostr << "\tMAF             : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_MAF+".maf"));
  ostr << "\tFASTA           : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_FASTAUNPADDED+".fasta"));
  ostr << "\tFASTA quality   : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_FASTAUNPADDED+".fasta.qual"));
  ostr << "\tFASTA (padded)  : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_FASTAPADDED+".fasta"));
  ostr << "\tFASTA qual.(pad): ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_FASTAPADDED+".fasta.qual"));
  ostr << "\tGAP4 (directory): ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outdir_GAP4DA+".gap4da"));
  ostr << "\tACE             : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_ACE+".ace"));
  ostr << "\tHTML            : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_HTML+".html"));
  ostr << "\tSimple text     : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_TXT+".txt"));
  ostr << "\tTCS overview    : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_TCS+".tcs"));
  ostr << "\tWiggle          : ";
  OUTSTRING(ostr, (Pv[0].mp_assembly_params.as_outfile_WIGGLE+".wig"));
}


// static func
void MIRAParameters::dumpOutputCustomisationParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=46;

  std::vector<int> singlePvIndex;

  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }


  ostr << "\n    Alignment output customisation:\n";

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_contig_params.con_output_text_cpl,
		  "\t", "TEXT characters per line (tcpl)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_contig_params.con_output_html_cpl,
		  "\t", "HTML characters per line (hcpl)",
		  fieldlength);

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_contig_params.con_output_text_gapfill,
		  "\t", "TEXT end gap fill character (tegfc)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_contig_params.con_output_html_gapfill,
		  "\t", "HTML end gap fill character (hegfc)",
		  fieldlength);

}


// static func
void MIRAParameters::dumpFileTempParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=46;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n    Temporary result files:\n";

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_caf,
		      "\t", "Saved as CAF                       (otc)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_maf,
		      "\t", "Saved as MAF                       (otm)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_fasta,
		      "\t", "Saved as FASTA                     (otf)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_gap4da,
		      "\t", "Saved as GAP4 (directed assembly)  (otg)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_ace,
		      "\t", "Saved as phrap ACE                 (ota)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_html,
		      "\t", "Saved as HTML                      (oth)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_tcs,
		      "\t", "Saved as Transposed Contig Summary (ots)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tmp_txt,
		      "\t", "Saved as simple text format        (ott)",
		      fieldlength);

  ostr << "\n    Extended temporary result files:\n";

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_exttmp_caf,
		      "\t", "Saved as CAF                      (oetc)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_exttmp_fasta,
		      "\t", "Saved as FASTA                    (oetf)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_exttmp_gap4da,
		      "\t", "Saved as GAP4 (directed assembly) (oetg)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_exttmp_ace,
		      "\t", "Saved as phrap ACE                (oeta)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_exttmp_html,
		      "\t", "Saved as HTML                     (oeth)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_exttmp_alsosinglets,
		      "\t", "Save also singlets               (oetas)",
		      fieldlength);
}


// static func
void MIRAParameters::dumpFileOutParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=46;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Output files (-OUTPUT/-OUT):\n";

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_savesimplesingletsinproject,
		      "\t", "Save simple singlets in project (sssip)",
		      fieldlength);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_savetaggedsingletsinproject,
		      "\t", "Save tagged singlets in project (stsip)",
		      fieldlength);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_removerollovertmps,
		      "\t", "Remove rollover tmps (rrot)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_removetmpdir,
		      "\t", "Remove tmp directory (rtd)",
		      fieldlength);


  ostr << "\n    Result files:\n";

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_caf,
		      "\t", "Saved as CAF                       (orc)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_maf,
		      "\t", "Saved as MAF                       (orm)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_fasta,
		      "\t", "Saved as FASTA                     (orf)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_gap4da,
		      "\t", "Saved as GAP4 (directed assembly)  (org)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_ace,
		      "\t", "Saved as phrap ACE                 (ora)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_gff3,
		      "\t", "Saved as GFF3                     (org3)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_html,
		      "\t", "Saved as HTML                      (orh)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_tcs,
		      "\t", "Saved as Transposed Contig Summary (ors)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_txt,
		      "\t", "Saved as simple text format        (ort)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_output_wiggle,
		      "\t", "Saved as wiggle                    (orw)",
		      fieldlength);
}


// static func
void MIRAParameters::dumpFileInParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  return;
  const int32 fieldlength=46;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  File names (-FN):\n";
  ostr << '\n';

}


// static func
void MIRAParameters::dumpDirectoryParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=35;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Directories (-DI):\n";

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_directory_params.dir_top,
		  "\t", "Top directory for writing files",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_directory_params.dir_results,
		  "\t", "For writing result files",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_directory_params.dir_info,
		  "\t", "For writing result info files",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_directory_params.dir_tmp,
		  "\t", "For writing tmp files",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_directory_params.dir_tmp_redirectedto,
		  "\t", "Tmp redirected to (trt)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_directory_params.dir_checkpoint,
		  "\t", "For writing checkpoint files",
		  fieldlength);

  // TODO: hide?
  //ostr << "\tFor writing gap4 DA res.: " << mp_assembly_params.as_outdir_GAP4DA << endl;

}


// static func
void MIRAParameters::dumpEditParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Edit options (-ED):\n";

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_edit_params.ed_sdbg_readedit,
		      "\t", "GB read editing (gbre)",
		      fieldlength);

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_edit_params.ed_mira_automatic_contic_editing,
		      "\t", "Mira automatic contig editing (mace)",
		      fieldlength);

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_edit_params.ed_kmer_singlets,
		      "\t    ", "Edit kmer singlets (eks)",
		      fieldlength-4);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_edit_params.ed_homopolymer_overcalls,
		      "\t    ", "Edit homopolymer overcalls (ehpo)",
		      fieldlength-4);


//  ostr << "     Sanger only:\n";
//
//  multiParamPrintBool(Pv, indexesInPv, ostr,
//		      Pv[0].mp_edit_params.ed_automatic_contic_editing,
//		      "\t", "EdIt automatic contig editing (eace)",
//		      fieldlength);
//  multiParamPrintBool(Pv, singlePvIndex, ostr,
//		      Pv[0].mp_edit_params.ed_strict_editing_mode,
//		      "\t", "Strict editing mode (sem)",
//		      fieldlength);
//  multiParamPrint(Pv, singlePvIndex, ostr,
//		  Pv[0].mp_edit_params.ed_confirmation_threshold,
//		  "\t", "Confirmation threshold in percent (ct)",
//		  fieldlength);
}


// static func
void MIRAParameters::dumpContigParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=58;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Contig parameters (-CO):\n";

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_contig_params.con_nameprefix,
		  "\t", "Name prefix (np)",
		  fieldlength);

// TODO: deprecate?
//  ostr << "\tError analysis (an)                                      : ";
//  switch(mp_contig_params.con_danger_analyse_mode){
//  case ANALYSE_SIGNAL:{
//    ostr << "SCF signal (signal)\n";
//    break;
//  }
//  case ANALYSE_TEXT:{
//    ostr << "text only (text)\n";
//    break;
//  }
//  case ANALYSE_NONE:{
//    ostr << "none (none)\n";
//    break;
//  }
//  default:{
//    ostr << "Unknown??? (please contact the authors)\n";
//  }
//  }


  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_reject_on_drop_in_relscore,
		  "\t", "Reject on drop in relative alignment score in % (rodirs)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_min_relscore,
		  "\t", "CMinimum relative score in % (cmrs)",
		  fieldlength);

// TODO: deprecate?
//  ostr << "\tMax. error rate in dangerous zones in % (dmer)           : " << mp_contig_params.con_danger_max_error_rate << "\n";


  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_mark_repeats,
		      "\t", "Mark repeats (mr)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_mark_repeats_onlyinresult,
		      "\t    ", "Only in result (mroir)",
		      fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_contig_params.con_assume_snp_insteadof_rmb,
		      "\t    ", "Assume SNP instead of repeats (asir)",
		      fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_minreadspergroup,
		  "\t    ", "Minimum reads per group needed for tagging (mrpg)",
		  fieldlength-4);
  multiParamPrintNumericChar(
    Pv, indexesInPv, ostr,
    Pv[0].mp_contig_params.con_minrmbneighbourqual,
    "\t    ", "Minimum neighbour quality needed for tagging (mnq)",
    fieldlength-4);
  multiParamPrintNumericChar(
    Pv, indexesInPv, ostr,
    Pv[0].mp_contig_params.con_mingroupqualforrmbtagging,
    "\t    ", "Minimum Group Quality needed for RMB Tagging (mgqrt)",
    fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_contig_params.con_mincoveragepercentage,
		  "\t    ", "Minimum coverage percentage (mcp)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_endreadmarkexclusionarea,
		  "\t    ", "End-read Marking Exclusion Area in bases (emea)",
		  fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_contig_params.con_emea_setzero_on_clipping_pec,
		      "\t        ", "Set to 1 on clipping PEC (emeas1clpec)",
		      fieldlength-8);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_contig_params.con_also_mark_gap_bases,
		      "\t    ", "Also mark gap bases (amgb)",
		      fieldlength-4);
  multiParamPrintBool(
    Pv, indexesInPv, ostr,
    Pv[0].mp_contig_params.con_also_mark_gap_bases_evenmc,
    "\t        ", "Also mark gap bases - even multicolumn (amgbemc)",
    fieldlength-8);
  multiParamPrintBool(
    Pv, indexesInPv, ostr,
    Pv[0].mp_contig_params.con_also_mark_gap_bases_needbothstrands,
    "\t        ", "Also mark gap bases - need both strands (amgbnbs)",
    fieldlength-8);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_contig_params.con_force_nonIUPACconsensus,
		      "\t", "Force non-IUPAC consensus (fnic)",
		      fieldlength);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_contig_params.con_mergeshortreads,
		      "\t", "Merge short reads (msr)",
		      fieldlength);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_msr_maxerrors,
		  "\t    ", "Max errors (msrme)",
		  fieldlength-4);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_msr_keependsunmapped,
		  "\t    ", "Keep ends unmerged (msrkeu)",
		  fieldlength-4);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_contig_params.con_gap_override_ratio,
		  "\t", "Gap override ratio (gor)",
		  fieldlength);
}



// static func
void MIRAParameters::dumpPathfinderParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Pathfinder options (-PF):\n";

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_pathfinder_params.paf_use_quick_rule,
		      "\t", "Use quick rule (uqr)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_pathfinder_params.paf_quickrule_minlen1,
		  "\t    ", "Quick rule min len 1 (qrml1)",
		  fieldlength-4);
  multiParamPrintNumericChar(Pv, indexesInPv, ostr,
			     Pv[0].mp_pathfinder_params.paf_quickrule_minsim1,
			     "\t    ", "Quick rule min sim 1 (qrms1)",
			     fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_pathfinder_params.paf_quickrule_minlen2,
		  "\t    ", "Quick rule min len 2 (qrml2)",
		  fieldlength-4);
  multiParamPrintNumericChar(Pv, indexesInPv, ostr,
			     Pv[0].mp_pathfinder_params.paf_quickrule_minsim2,
			     "\t    ", "Quick rule min sim 2 (qrms2)",
			     fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_pathfinder_params.paf_bbquickoverlap_minlen,
		  "\t", "Backbone quick overlap min len (bqoml)",
		  fieldlength);

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_pathfinder_params.paf_max_startcache_filltime,
		  "\t", "Max. start cache fill time (mscft)",
		  fieldlength);
}

// static func
void MIRAParameters::dumpSkimParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Parameters for SKIM algorithm (-SK):\n";

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_numthreads,
		  "\t", "Number of threads (not)",
		  fieldlength);
  ostr << '\n';
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_skim_params.sk_alsoskimrevcomp,
		      "\t", "Also compute reverse complements (acrc)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_basesperhash,
		  "\t", "Kmer size (kms)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_bph_increasestep,
		  "\t    ", "Automatic increase per pass (kmsaipp)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_bph_max,
		  "\t    ", "Kmer size max(kmsmax)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_hashsavestepping,
		  "\t", "Kmer save stepping (kss)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_skim_params.sk_percentrequired,
		  "\t", "Percent required (pr)",
		  fieldlength);
  ostr << '\n';
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_maxhitsperread,
		  "\t", "Max hits per read (mhpr)",
		  fieldlength);
  ostr << '\n';
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_skim_params.sk_filtermegahubs,
		      "\t", "Filter megahubs (fmh)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_megahubcap,
		  "\t    ", "Megahub cap (mhc)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_maxmegahubratio,
		  "\t    ", "Max megahub ratio (mmhr)",
		  fieldlength-4);

  ostr << '\n';
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_skim_params.sk_swcheckonbackbones,
		      "\t", "SW check on backbones (swcob)",
		      fieldlength);

  ostr << '\n';
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_maxhashesinmem,
		  "\t", "Max kmers in memory (mkim)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_skim_params.sk_memcaphitreduction,
		  "\t", "MemCap: hit reduction (mchr)",
		  fieldlength);

}



// static func
void MIRAParameters::dumpHashStatisticsParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Parameters for Kmer Statistics (-KS):\n";

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_freq_covestmin,
		  "\t", "Freq. cov. estim. min (fcem)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_freqest_minnormal,
		  "\t", "Freq. estim. min normal (fenn)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_freqest_maxnormal,
		  "\t", "Freq. estim. max normal (fexn)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_freqest_repeat,
		  "\t", "Freq. estim. repeat (fer)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_freqest_heavyrepeat,
		  "\t", "Freq. estim. heavy repeat (fehr)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_freqest_crazyrepeat,
		  "\t", "Freq. estim. crazy (fecr)",
		  fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_hashstatistics_params.hs_masknastyrepeats,
		      "\t", "Mask nasty repeats (mnr)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_nastyrepeatratio,
		  "\t    ", "Nasty repeat ratio (nrr)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_nastyrepeatcoverage,
		  "\t    ", "Nasty repeat coverage (nrc)",
		  fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_hashstatistics_params.hs_apply_digitalnormalisation,
		      "\t    ", "Lossless digital normalisation (ldn)",
		      fieldlength-4);

  ostr << '\n';
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_repeatlevel_in_infofile,
		  "\t", "Repeat level in info file (rliif)",
		  fieldlength);

  ostr << '\n';
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_memtouse,
		  "\t", "Memory to use (mtu)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_hashstatistics_params.hs_rare_kmer_final_kill,
		  "\t", "Rare kmer final kill (rkfk)",
		  fieldlength);
}



// static func
void MIRAParameters::dumpClippingParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Clipping options (-CL):\n";
  ostr << "\tSSAHA2 or SMALT clipping:\n";
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_ssahamerge_gapsize,
		  "\t    ", "Gap size (msvsgs)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_ssahamerge_maxfrontgap,
		  "\t    ", "Max front gap (msvsmfg)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_ssahamerge_maxendgap,
		  "\t    ", "Max end gap (msvsmeg)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_ssahamerge_strictfrontclip,
		  "\t    ", "Strict front clip (msvssfc)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_ssahamerge_strictendclip,
		  "\t    ", "Strict end clip (msvssec)",
		  fieldlength-4);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_possible_vectors,
		      "\t", "Possible vector leftover clip (pvlc)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_vector_maxlenallowed,
		  "\t    ", "maximum len allowed (pvcmla)",
		  fieldlength-4);

  multiParamPrintNumericChar(Pv, indexesInPv, ostr,
			     Pv[0].mp_assembly_params.as_clip_quality_minthreshold,
			     "\t", "Min qual. threshold for entire read (mqtfer)",
			     fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_quality_numminthreshold,
		  "\t    ", "Number of bases (mqtfernob)",
		  fieldlength-4);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_quality,
		      "\t", "Quality clip (qc)",
		      fieldlength);
  multiParamPrintNumericChar(Pv, indexesInPv, ostr,
			     Pv[0].mp_assembly_params.as_clip_quality_minqual,
			     "\t    ", "Minimum quality (qcmq)",
			     fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_quality_winlen,
		  "\t    ", "Window length (qcwl)",
		  fieldlength-4);


  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_badstretchquality,
		      "\t", "Bad stretch quality clip (bsqc)",
		      fieldlength);
  multiParamPrintNumericChar(Pv, indexesInPv, ostr,
			     Pv[0].mp_assembly_params.as_clip_badstretchquality_minqual,
			     "\t    ", "Minimum quality (bsqcmq)",
			     fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_badstretchquality_winlen,
		  "\t    ", "Window length (bsqcwl)",
		  fieldlength-4);



  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_maskedbases,
		      "\t", "Masked bases clip (mbc)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_maskedbase_gapsize,
		  "\t    ", "Gap size (mbcgs)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_maskedbase_maxfrontgap,
		  "\t    ", "Max front gap (mbcmfg)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_maskedbase_maxendgap,
		  "\t    ", "Max end gap (mbcmeg)",
		  fieldlength-4);


  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_lowercase_front,
		      "\t", "Lower case clip front (lccf)",
		      fieldlength);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_lowercase_back,
		      "\t", "Lower case clip back (lccb)",
		      fieldlength);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_polyat,
		      "\t", "Clip poly A/T at ends (cpat)",
		      fieldlength);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_polyat_keeppolystretch,
		      "\t    ", "Keep poly-a signal (cpkps)",
		      fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_polyat_len,
		  "\t    ", "Minimum signal length (cpmsl)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_polyat_maxerrors,
		  "\t    ", "Max errors allowed (cpmea)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_polyat_maxgap,
		  "\t    ", "Max gap from ends (cpmgfe)",
		  fieldlength-4);


  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_3ppolybase,
		      "\t", "Clip 3 prime polybase (c3pp)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_3ppolybase_len,
		  "\t    ", "Minimum signal length (c3ppmsl)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_3ppolybase_maxerrors,
		  "\t    ", "Max errors allowed (c3ppmea)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_3ppolybase_maxgap,
		  "\t    ", "Max gap from ends (c3ppmgfe)",
		  fieldlength-4);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_knownadaptorsright,
		      "\t", "Clip known adaptors right (ckar)",
		      fieldlength);

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_ensureminimumleftclipoff,
		      "\t", "Ensure minimum left clip (emlc)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_minslrequired,
		  "\t    ", "Minimum left clip req. (mlcr)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_minqlsetto,
		  "\t    ", "Set minimum left clip to (smlc)",
		  fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_ensureminimumrightclipoff,
		      "\t", "Ensure minimum right clip (emrc)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_minsrrequired,
		  "\t    ", "Minimum right clip req. (mrcr)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_minqrsetto,
		  "\t    ", "Set minimum right clip to (smrc)",
		  fieldlength-4);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_sdbg_chimeradetection,
		      "\t", "GB chimera detection clip (gbcdc)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_kmer_junkdetection,
		      "\t", "KMER junk detection (kjd)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_kmer_junkkill,
		      "\t    ", "KMER junk complete kill (kjck)",
		      fieldlength-4);

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_skimchimeradetection,
		      "\t", "DEPRECATED! Apply SKIM chimera detection clip (ascdc)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_skimjunkdetection,
		      "\t", "DEPRECATED! Apply SKIM junk detection clip (asjdc)",
		      fieldlength);
  ostr << '\n';

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_proposeendclips,
		      "\t", "Propose end clips (pec)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_clip_pec_basesperhash,
		  "\t    ", "Kmer size (peckms)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_clip_pec_mkfr,
		  "\t    ", "Minimum kmer for forward-rev (pmkfr)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_clip_pec_mtk,
		  "\t    ", "Minimum total kmer (pmtk)",
		  fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_sxaggcxg,
		      "\t    ", "Handle Solexa GGCxG problem (pechsgp)",
		      fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_continuous,
		      "\t    ", "Continuous (pecc)",
		      fieldlength-4);

  ostr << '\n';

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clipmask_rarekmers,
		      "\t", "Rare kmer mask (rkm)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_pec_ffreq,
		  "\t    ", "Front freq (pffreq)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_clip_pec_bfreq,
		  "\t    ", "Back freq (pbfreq)",
		  fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_ffr,
		      "\t    ", "Front forward-rev (pffore)",
		      fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_bfr,
		      "\t    ", "Back forward-rev (pbfore)",
		      fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_fcmst,
		      "\t    ", "Front conf. multi-seq type (pfcmst)",
		      fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_bcmst,
		      "\t    ", "Back conf. multi-seq type (pbcmst)",
		      fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_fsalp,
		      "\t    ", "Front seen at low pos (pfsalp)",
		      fieldlength-4);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_pec_bsalp,
		      "\t    ", "Back seen at low pos (pbsalp)",
		      fieldlength-4);

  ostr << '\n';

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_clip_badsolexaends,
		      "\t", "Clip bad solexa ends (cbse)",
		      fieldlength);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_search_phix174,
		      "\t", "Search PhiX174 (spx174)",
		      fieldlength);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_filter_phix174,
		      "\t    ", "Filter PhiX174 (fpx174)",
		      fieldlength-4);


  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_filter_rrna,
		      "\t", "Filter rRNA (frrna)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_filter_rrna_pairs,
		      "\t    ", "Pairs (frrnap)",
		      fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_filter_rrna_numkmers,
		  "\t    ", "Number of kmers (frrnank)",
		  fieldlength-4);

  return;
}


// static func
void MIRAParameters::dumpDataProcessingParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
    if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Dataprocessing options (-DP):\n";

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_use_read_extension,
		      "\t", "Use read extensions (ure)",
		      fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_readextension_window_len,
		  "\t    ", "Read extension window length (rewl)",
		  fieldlength-4);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_readextension_window_maxerrors,
		  "\t    ", "Read extension w. maxerrors (rewme)",
		  fieldlength-4);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_readextension_firstpassnum,
		  "\t    ", "First extension in pass (feip)",
		  fieldlength-4);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_readextension_lastpassnum,
		  "\t    ", "Last extension in pass (leip)",
		  fieldlength-4);
}

// static func
void MIRAParameters::dumpStrainBackboneParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
    if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Strain and backbone options (-SB):\n";

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_backbone_bootstrapnewbackbone,
		      "\t", "Bootstrap new backbone (bnb)",
		      fieldlength);

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_startbackboneusage_inpass,
		  "\t","Start backbone usage in pass (sbuip)",
		  fieldlength);
//  multiParamPrintBool(Pv, singlePvIndex, ostr,
//		      Pv[0].mp_assembly_params.as_backbone_strainname_forceforall,
//		      "\t    ", "Force for all (bsnffa)",
//		      fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_backbone_rail_fromstrain,
		  "\t", "Backbone rail from strain (brfs)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_backbone_raillength,
		  "\t", "Backbone rail length (brl)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_backbone_railoverlap,
		  "\t", "Backbone rail overlap (bro)",
		  fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_backbone_trimoverhangingreads,
		      "\t", "Trim overhanging reads (tor)",
		      fieldlength);
}

// static func
void MIRAParameters::dumpAssemblyParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
    if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Assembly options (-AS):\n";

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_pathfinder_params.paf_use_genomic_algorithms,
		      "\t", "Use genomic pathfinder (ugpf)",
		      fieldlength);

  ostr << '\n';

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_numpasses,
		  "\t", "Number of passes (nop)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_user_kmerseries,
		  "\t", "Kmer series (kms)",
		  fieldlength);
//  multiParamPrintBool(Pv, singlePvIndex, ostr,
//		      Pv[0].mp_assembly_params.as_skimeachpass,
//		      "\t    ", "Skim each pass (sep)",
//		      fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_numrmbbreakloops,
		  "\t", "Maximum number of RMB break loops (rbl)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_maxcontigsperpass,
		  "\t", "Maximum contigs per pass (mcpp)",
		  fieldlength);

  ostr << '\n';

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_minimum_readlength,
		  "\t", "Minimum read length (mrl)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_minimum_readspercontig,
		  "\t", "Minimum reads per contig (mrpc)",
		  fieldlength);
  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_enforce_qualsinreads,
		      "\t", "Enforce presence of qualities (epoq)",
		      fieldlength);


#ifndef PUBLICQUIET
  ostr << '\n';
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_rle_reads,
		  "\t", "Use run length encoding (urle)",
		  fieldlength);
#endif

  ostr << '\n';
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_automatic_repeat_detection,
		  "\t", "Automatic repeat detection (ard)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_ard_multicopythreshold,
		  "\t    ", "Coverage threshold (ardct)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_ard_multicopyminlen,
		  "\t    ", "Minimum length (ardml)",
		  fieldlength-4);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_ard_multicopygrace,
		  "\t    ", "Grace length (ardgl)",
		  fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_uniform_read_distribution,
		      "\t    ", "Use uniform read distribution (urd)",
		      fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_urd_startinpass,
		  "\t      ", "Start in pass (urdsip)",
		  fieldlength-6);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_assembly_params.as_urd_cutoffmultiplier,
		  "\t      ", "Cutoff multiplier (urdcm)",
		  fieldlength-6);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_spoilerdetection,
		      "\t", "Spoiler detection (sd)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_spdetect_lastpassonly,
		      "\t    ", "Last pass only (sdlpo)",
		      fieldlength-4);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_pathfinder_params.paf_use_emergency_search_stop,
		      "\t", "Use emergency search stop (uess)",
		      fieldlength);

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_pathfinder_params.paf_ess_depth,
		  "\t    ", "ESS partner depth (esspd)",
		  fieldlength-4);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_pathfinder_params.paf_use_emergency_blacklist,
		      "\t", "Use emergency blacklist (uebl)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_pathfinder_params.paf_use_max_contig_buildtime,
		      "\t", "Use max. contig build time (umcbt)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_pathfinder_params.paf_max_contig_buildtime,
		  "\t    ", "Build time in seconds (bts)",
		  fieldlength-4);

}



// static func
void MIRAParameters::dumpLoadParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
    if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Load reads options (-LR):\n";

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_assembly_params.as_wants_qualityfile,
		      "\t", "Wants quality file (wqf)",
		      fieldlength);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_filecheck_only,
		      "\t",
		      "Filecheck only (fo)",
		      fieldlength);

}

void MIRAParameters::dumpGeneralParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    cout << "Pushing back " << indexesInPv.front() << endl;
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "  General (-GE):\n";

  multiParamPrint(Pv, singlePvIndex, ostr, Pv[0].mp_assembly_params.as_projectname_out,
		  "\t", "Project name",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_numthreads,
		  "\t", "Number of threads (not)",
		  fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_automemmanagement,
		      "\t", "Automatic memory management (amm)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_amm_keeppercentfree,
		  "\t    ", "Keep percent memory free (kpmf)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_assembly_params.as_amm_maxprocesssize,
		  "\t    ", "Max. process size (mps)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_special_params.sp_est_startstep,
		  "\t",
		  "EST SNP pipeline step (esps)",
		  fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_buntify_reads,
		      "\t",
		      "Colour reads by kmer frequency (crkf)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_assembly_params.as_assemblyjob_preprocessonly,
		      "\t",
		      "Preprocess only (ppo)",
		      fieldlength);
}

void MIRAParameters::dumpMiscParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Misc (-MI):\n";

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_special_params.mi_as_largecontigsize,
		  "\t", "Large contig size (lcs)",
		  fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_special_params.mi_as_largecontigsize4stats,
		  "\t", "Large contig size for stats (lcs4s)",
		  fieldlength);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_special_params.mi_iknowwhatido,
		      "\t", "I know what I do (ikwid)",
		      fieldlength);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_special_params.mi_extra_flag1,
		      "\t", "Extra flag 1 / sanity track check (ef1)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_special_params.mi_extra_flag2,
		      "\t", "Extra flag 2 / dnredreadsatpeaks (ef2)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_special_params.mi_extra_flag3,
		      "\t", "Extra flag 3 / pelibdisassemble (ef3)",
		      fieldlength);
  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_special_params.mi_extended_log,
		      "\t", "Extended log (el)",
		      fieldlength);
}

void MIRAParameters::dumpNagAndWarnParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=45;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Nag and Warn (-NW):\n";

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_nfs,
			 "\t", "Check NFS (cnfs)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_multipassmapping,
			 "\t", "Check multi pass mapping (cmpm)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_templateproblems,
			 "\t", "Check template problems (ctp)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_readnamesra,
			 "\t", "Check SRA read names (csrn)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_duplicatereadnames,
			 "\t", "Check duplicate read names (cdrn)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_illuminajunkinest,
			 "\t", "Check Illumina junk in EST (cijie)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_proposedendclipping,
			 "\t", "Check proposed end clipping (cpec)",
			 fieldlength);

  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_maxreadnamelength,
			 "\t", "Check max read name length (cmrnl)",
			 fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_nagandwarn_params.nw_check_mrnlvalue,
		  "\t    ", "Max read name length (mrnl)",
		  fieldlength-4);
  multiParamPrintNagWarn(Pv, singlePvIndex, ostr,
			 Pv[0].mp_nagandwarn_params.nw_check_coverage,
			 "\t", "Check average coverage (cac)",
			 fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_nagandwarn_params.nw_check_covvalue,
		  "\t    ", "Average coverage value (acv)",
		  fieldlength-4);
  multiParamPrint(Pv, singlePvIndex, ostr, Pv[0].mp_nagandwarn_params.nw_warn_x174prefix,
		  "\t", "(Phi)X174 contig prefix (x174cp)",
		  fieldlength);
}

// static func
void MIRAParameters::dumpAlignParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=40;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Align parameters for Smith-Waterman align (-AL):\n";

  multiParamPrint(Pv, indexesInPv, ostr, Pv[0].mp_align_params.al_kpercent,
		  "\t", "Bandwidth in percent (bip)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_align_params.al_kmax,
		  "\t",
		  "Bandwidth max (bmax)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_align_params.al_kmin,
		  "\t",
		  "Bandwidth min (bmin)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_align_params.al_min_score,
		  "\t",
		  "Minimum score (ms)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_align_params.al_min_overlap,
		  "\t",
		  "Minimum overlap (mo)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_align_params.al_min_relscore,
		  "\t",
		  "Minimum relative score in % (mrs)",
		  fieldlength);

  ostr << '\n';

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_align_params.ads_enforce_clean_ends,
		      "\t",
		      "Enforce clean ends (ece)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_align_params.ads_clean_end_distance,
		  "\t    ",
		  "Clean end distance (ced)",
		  fieldlength-4);

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_align_params.ads_clean_end_mismatchallowed,
		  "\t    ",
		  "Clean end mismatch allowed (cema)",
		  fieldlength-4);

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_align_params.ads_extra_gap_penalty,
		      "\t",
		      "Extra gap penalty (egp)",
		      fieldlength);
  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_align_params.ads_gp_functionstring,
		  "\t    ",
		  "extra gap penalty levels (egpl)",
		  fieldlength-4);

  multiParamPrint(Pv, singlePvIndex, ostr,
		  Pv[0].mp_align_params.ads_max_gppercent,
		  "\t    ",
		  "Max. egp in percent (megpp)",
		  fieldlength-4);
}


void MIRAParameters::dumpFinalMappingParams(std::vector<MIRAParameters> & Pv, const std::vector<int> & indexesInPv, std::ostream & ostr)
{
  const int32 fieldlength=40;

  std::vector<int> singlePvIndex;
  if(indexesInPv.size()==1){
    singlePvIndex.push_back(indexesInPv.front());
  }else{
    singlePvIndex.push_back(0);
  }

  ostr << "\n  Final mapping parameters (-FM):\n";

  multiParamPrintBool(Pv, singlePvIndex, ostr,
		      Pv[0].mp_finalmap_params.fm_active,
		      "\t",
		      "Active (act)",
		      fieldlength);

  ostr << '\n';

  multiParamPrintBool(Pv, indexesInPv, ostr,
		      Pv[0].mp_finalmap_params.fm_mapperfect,
		      "\t",
		      "Map perfect matches (mpm)",
		      fieldlength);

  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_finalmap_params.fm_maxtotalerrors,
		  "\t",
		  "Max total errors (mte)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_finalmap_params.fm_maxmismatches,
		  "\t",
		  "Max mismatches (mmm)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_finalmap_params.fm_maxgaps,
		  "\t",
		  "Max total errors (mg)",
		  fieldlength);
  multiParamPrint(Pv, indexesInPv, ostr,
		  Pv[0].mp_finalmap_params.fm_clean_end_dist,
		  "\t",
		  "Clean end distance (ced)",
		  fieldlength);

}


void MIRAParameters::consistencyCheck(bool verbose)
{
  checkMin("-GE:not",mp_assembly_params.as_numthreads,
	   static_cast<uint32>(1),
	   static_cast<uint32>(1), verbose);
  checkMax("-GE:not",mp_assembly_params.as_numthreads,
	   static_cast<uint32>(256),
	   static_cast<uint32>(256), verbose);
  checkMin("-AS:mrl",mp_assembly_params.as_minimum_readlength,
	   static_cast<uint32>(20),
	   static_cast<uint32>(20), verbose);
  checkMin("-AS:bts",mp_pathfinder_params.paf_max_contig_buildtime,60,120, verbose);
  checkMin("-AS:rbl",mp_assembly_params.as_numrmbbreakloops,
	   static_cast<uint32>(1),
	   static_cast<uint32>(1), verbose);

  checkMin("-CL:peckms",mp_assembly_params.as_clip_pec_basesperhash,
	   static_cast<uint32>(10),
	   static_cast<uint32>(10), verbose);

  checkMax("-SK:kms",mp_skim_params.sk_basesperhash,
	   static_cast<uint32>(256),
	   static_cast<uint32>(256), verbose);
//  if(sizeof(unsigned long) == 4) {
//    if(mp_skim_params.sk_basesperhash>=16) {
//      cout << "We're on a 32 bit machine, must adapt -SK:bph.\n";
//    }
//    checkMax("-SK:bph",mp_skim_params.sk_basesperhash,
//	     static_cast<uint32>(15),
//	     static_cast<uint32>(15), verbose);
//  }
  checkMax("-SB:sbuip",mp_assembly_params.as_startbackboneusage_inpass,
	   static_cast<int32>(mp_assembly_params.as_numpasses),
	   static_cast<int32>(mp_assembly_params.as_numpasses), verbose);

  // TODO: more consistency checks down the parameter list

  checkMin("-OUT:tcpl",mp_contig_params.con_output_text_cpl,10,60, verbose);
  checkMin("-OUT:hcpl",mp_contig_params.con_output_html_cpl,10,60, verbose);
}

void MIRAParameters::generateProjectNames(std::vector<MIRAParameters> & Pv, std::string name)
{
  generateProjectInNames(Pv, name);
  generateProjectOutNames(Pv, name);

  return;
}


// TODO: to be retired as soon as ssaha2 / smalt is integrated into manifest loading system
void MIRAParameters::generateProjectInNames(std::vector<MIRAParameters> & Pv, std::string name) {
  if(name.empty()) name=Pv[0].mp_assembly_params.as_projectname_in;

  Pv[0].mp_assembly_params.as_projectname_in=name;

  return;
}

void MIRAParameters::generateProjectOutNames(std::vector<MIRAParameters> & Pv, std::string name) {
  if(name.empty()) name=Pv[0].mp_assembly_params.as_projectname_out;
  Pv[0].mp_assembly_params.as_projectname_out=name;
  Pv[0].mp_assembly_params.as_outfile_FASTAUNPADDED=name+"_out.unpadded";
  Pv[0].mp_assembly_params.as_outfile_FASTAUNPADDEDQUAL=name+"_out.unpadded";
  Pv[0].mp_assembly_params.as_outfile_FASTAPADDED=name+"_out.padded";
  Pv[0].mp_assembly_params.as_outfile_FASTAPADDEDQUAL=name+"_out.padded";
  Pv[0].mp_assembly_params.as_outfile_FASTA=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_CAF=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_MAF=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_HTML=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_TXT=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_ACE=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_WIGGLE=name+"_out";
  Pv[0].mp_assembly_params.as_outfile_TCS=name+"_out";
  Pv[0].mp_assembly_params.as_outdir_GAP4DA=name+"_out";

  Pv[0].mp_file_params.chkpt_persistentoverlaps="povl.txt";
  Pv[0].mp_file_params.chkpt_bannedoverlaps="bannedOverlaps.txt";
  Pv[0].mp_file_params.chkpt_maxcovreached="maxCovReached.txt";
  Pv[0].mp_file_params.chkpt_passinfo="passInfo.txt";
  Pv[0].mp_file_params.chkpt_readpool="readpool.maf";

  Pv[0].mp_assembly_params.as_outfile_stats_reads_invalid=name+"_info_reads_invalid";
  Pv[0].mp_assembly_params.as_outfile_stats_reads_tooshort=name+"_info_reads_tooshort";
  Pv[0].mp_assembly_params.as_outfile_stats_contigstats=name+"_info_contigstats";
  Pv[0].mp_assembly_params.as_outfile_stats_info=name+"_info_assembly";
  Pv[0].mp_assembly_params.as_outfile_stats_warnings=name+"_info_WARNINGS";
  Pv[0].mp_assembly_params.as_outfile_stats_debrislist=name+"_info_debrislist";
  Pv[0].mp_assembly_params.as_outfile_stats_crlist=name+"_info_contigreadlist";
  Pv[0].mp_assembly_params.as_outfile_stats_readtags=name+"_info_readtaglist";
  Pv[0].mp_assembly_params.as_outfile_stats_contigtags=name+"_info_consensustaglist";
  Pv[0].mp_assembly_params.as_outfile_stats_snpanalysis=name+"_info_snplist";
  Pv[0].mp_assembly_params.as_outfile_stats_snpenvironment=name+"_info_snpenvironment";
  Pv[0].mp_assembly_params.as_outfile_stats_featureanalysis=name+"_info_featureanalysis";
  Pv[0].mp_assembly_params.as_outfile_stats_featuresummary=name+"_info_featuresummary";
  Pv[0].mp_assembly_params.as_outfile_stats_featuresequences=name+"_info_featuresequences";
  Pv[0].mp_assembly_params.as_outfile_stats_featurecoverage=name+"_info_featurecoverage";
  Pv[0].mp_assembly_params.as_outfile_stats_readrepeats=name+"_info_readrepeats";
  Pv[0].mp_assembly_params.as_outfile_stats_largecontigs=name+"_info_largecontigs";
  Pv[0].mp_assembly_params.as_outfile_stats_rgstinfo=name+"_info_readgroups";

  Pv[0].mp_assembly_params.as_tmpf_spoiler=name+"_int_contigjoinspoiler";
  Pv[0].mp_assembly_params.as_tmpf_adsextend=name+"_int_alignextends";
  Pv[0].mp_assembly_params.as_tmpf_vectorclip=name+"_int_vectorclip";
  Pv[0].mp_assembly_params.as_tmpf_posmatch=name+"_int_posmatch";
  //Pv[0].mp_assembly_params.as_tmpf_overlap_criterion_levels=name+"_int_oclevel";
  Pv[0].mp_assembly_params.as_tmpf_normalisedskim=name+"_int_normalisedskims";
  Pv[0].mp_assembly_params.as_tmpf_clippings=name+"_int_clippings";
  Pv[0].mp_assembly_params.as_tmpf_ads=name+"_int_ads";
  Pv[0].mp_assembly_params.as_tmpf_poolinfo=name+"_log_readpoolinfo";
  Pv[0].mp_assembly_params.as_tmpf_kmerstatistics=name+"_int_kmerstats";

  Pv[0].mp_assembly_params.as_tmpf_signal_findpossibleoverlaps=name+"_signal_findpossibleoverlaps";
  Pv[0].mp_assembly_params.as_tmpf_signal_mainalignments=name+"_signal_mainalignments";
  Pv[0].mp_assembly_params.as_tmpf_signal_kmerstats=name+"_signal_kmerstats";

  Pv[0].mp_assembly_params.as_tmpf_wellconnected=name+"_int_wellconnected";
  Pv[0].mp_assembly_params.as_tmpf_banned_overlaps=name+"_int_banned_overlaps";
  Pv[0].mp_assembly_params.as_tmpf_istroublemaker=name+"_int_istroublemaker";
  Pv[0].mp_assembly_params.as_tmpf_needalloverlaps=name+"_int_needalloverlaps";
  Pv[0].mp_assembly_params.as_tmpf_multicopies=name+"_int_multicopies";
  Pv[0].mp_assembly_params.as_tmpf_hasmcoverlap=name+"_int_hasmcoverlap";
  Pv[0].mp_assembly_params.as_tmpf_debrisreason=name+"_int_debrisreason";
  Pv[0].mp_assembly_params.as_tmpf_skimmegahubs=name+"_int_skimmegahubs";

  std::string topdir=name+"_assembly";
  Pv[0].mp_directory_params.dir_top=topdir;
  Pv[0].mp_directory_params.dir_tmp_symlink.clear();
  Pv[0].mp_directory_params.dir_tmp=topdir+"/"+name+"_d_tmp";
  Pv[0].mp_directory_params.dir_results=topdir+"/"+name+"_d_results";
  Pv[0].mp_directory_params.dir_info=topdir+"/"+name+"_d_info";
  Pv[0].mp_directory_params.dir_checkpoint=topdir+"/"+name+"_d_chkpt";
  Pv[0].mp_directory_params.dir_checkpoint_tmp=topdir+"/"+name+"_d_chkpt_tmp";

  //Pv[0].mp_assembly_params.as_tmpf_unused_ids="miratmp.unused_ids";
  Pv[0].mp_assembly_params.as_tmpf_unused_ids="";

  Pv[0].mp_assembly_params.as_tmpf_scfreadfail=name+"_info_scfreadfail";
  Pv[0].mp_assembly_params.as_tmpf_scfreadfatallywrong=name+"_error_scfreadfatallywrong";

  Pv[0].mp_contig_params.con_nameprefix=name;

  return;
}

void MIRAParameters::correctTmpDirectory(std::vector<MIRAParameters> & Pv)
{
  if(Pv[0].mp_directory_params.dir_tmp_redirectedto.size()){
    Pv[0].mp_directory_params.dir_tmp_symlink=Pv[0].mp_directory_params.dir_tmp;
    Pv[0].mp_directory_params.dir_tmp=
      Pv[0].mp_directory_params.dir_tmp_redirectedto
      +"/"
      +Pv[0].mp_assembly_params.as_projectname_out
      +"_d_tmp";
  }
}

//#define OUTYESNO(o, yn) { if(yn){o << "Yes\n";}else{ostr << "No\n";}}
//#define OUTCHARSTRING(o,xx) {o<<"'"<<(xx)<<"'\n";}


std::ostream & operator<<(std::ostream &ostr, MIRAParameters const  &mp)
{

  std::vector<MIRAParameters> Pv;
  Pv.push_back(mp);
  std::vector<int> bla;
  bla.push_back(0);

  MIRAParameters::dumpAllParams(
    Pv,
    bla,
    ostr);

  return ostr;
}



void MIRAParameters::setPathfinderMaxContigTime(uint32 t)
{
  mp_pathfinder_params.paf_maxcontigclockticks=t*sysconf(_SC_CLK_TCK);
}

void MIRAParameters::setPathfinderNextReadMaxTimeToESS(uint32 t)
{
  mp_pathfinder_params.paf_nextread_maxcttoess=t*sysconf(_SC_CLK_TCK);
}

void MIRAParameters::setContigForceNonIUPAC(bool perseq, bool amongseq)
{
  mp_contig_params.con_force_nonIUPACconsensus_perseqtype=perseq;
  mp_contig_params.con_force_nonIUPACconsensus=amongseq;
}

//void MIRAParameters::setAlignGapPenaltyLevel(uint32 level)
bool MIRAParameters::setAlignGapPenaltyLevel(const std::string & gplevels, std::stringstream * errstreamptr)
{
  FUNCSTART("bool MIRAParameters::setAlignGapPenaltyLevel(const std::string & gplevels, std::stringstream * errstream)");

  static boost::char_separator<char> separator(",");

  bool retval_error=false;

  mp_align_params.ads_gp_functionstring=gplevels;
  mp_align_params.ads_gp_function.clear();
  // 0 stars: no additional penalty
  mp_align_params.ads_gp_function.push_back(0);

  boost::tokenizer<boost::char_separator<char> > tok(gplevels,separator);
  for(auto & te : tok){
    std::string tmps(te);
    boost::trim(tmps);
    auto tmp=atoi(tmps.c_str());
    if(tmp<0 || tmp>100){
      std::stringstream msg;
      msg << "ERROR -AS:kmerseries: "
	  << gplevels << "\n" << tmps << " is <= 0. kmer values must be >0\n";
      if(errstreamptr!=nullptr) {
	*errstreamptr << msg.str();
	retval_error=true;
      }else{
	MIRANOTIFY(Notify::FATAL,msg.str());
      }
    }
    mp_align_params.ads_gp_function.push_back(tmp);
  }

  return retval_error;
}

/*
  switch(level){
  case 0:{
    mp_align_params.ads_gp_function.push_back(0);     // 1 gap
    mp_align_params.ads_gp_function.push_back(5);     // 2 gaps
    mp_align_params.ads_gp_function.push_back(10);    // 3 gaps etc.
    mp_align_params.ads_gp_function.push_back(20);
    mp_align_params.ads_gp_function.push_back(40);
    mp_align_params.ads_gp_function.push_back(80);
    mp_align_params.ads_gp_function.push_back(100);
    break;
  }
  case 1:{
    mp_align_params.ads_gp_function.push_back(0);
    mp_align_params.ads_gp_function.push_back(10);
    mp_align_params.ads_gp_function.push_back(25);
    mp_align_params.ads_gp_function.push_back(50);
    mp_align_params.ads_gp_function.push_back(100);
    break;
  }
  case 2:{
    mp_align_params.ads_gp_function.push_back(0);
    mp_align_params.ads_gp_function.push_back(10);
    mp_align_params.ads_gp_function.push_back(50);
    mp_align_params.ads_gp_function.push_back(100);
    break;
  }
  case 10:{
    // jump function for codon sized gaps
    mp_align_params.ads_gp_function.push_back(0);     // 1 gap
    mp_align_params.ads_gp_function.push_back(5);     // 2 gaps
    mp_align_params.ads_gp_function.push_back(100);   // 3 or more gaps
    break;
  }
  default:{
    mp_align_params.ads_gp_function.push_back(0);
    mp_align_params.ads_gp_function.push_back(10);
    mp_align_params.ads_gp_function.push_back(50);
    mp_align_params.ads_gp_function.push_back(100);
  }
  }
}
*/


void MIRAParameters::loadParams(const std::string & pfile, std::vector<MIRAParameters> & Pv)
{
  FUNCSTART("void MIRAParameters::loadParams(const std::string & pfile)");

  cout << "Loading parameters from file: " << pfile << endl;

  if(MP_loadfilename.size()>=10) {
    MIRANOTIFY(Notify::FATAL, "Already loading 10 other files ... there's something unusual about that ... really: " << pfile);
  }

  for(uint32 i=0; i<MP_loadfilename.size(); i++) {
    if(MP_loadfilename[i] == pfile) {
      MIRANOTIFY(Notify::FATAL, "Already loading that file, recursion, sorry: " << pfile);
    }
  }

  MP_loadfilename.push_back(pfile);

  {
    std::ifstream fin(pfile, std::ios::in);
    if(!fin){
      MIRANOTIFY(Notify::FATAL, "File not found: " << pfile);
    }
    parse(fin, Pv, false);
  }

  MP_loadfilename.pop_back();

  Pv[0].consistencyCheck(true);

  FUNCEND();
  return;
}

void MIRAParameters::parse(int argc, char **argv, std::vector<MIRAParameters> & Pv, bool verbose)
{
  std::stringstream tss;

  cout << "Parsing parameters:";
  for(int32 i=1; i<argc; i++){
    cout << " " << argv[i];
    //    cout << "Doing " << i << "\t" << argv[i] << endl;
    tss << argv[i] << "  *=BEGIN0=*";
  }
  cout << "\n\n";

  parse(tss, Pv, verbose);
}

void MIRAParameters::parse(const char * params, std::vector<MIRAParameters> & Pv, bool verbose)
{
  std::stringstream tss;
  tss << params;
  parse(tss, Pv, verbose);
}

void MIRAParameters::parseQuickmode(const char * params, const char * qm, std::vector<MIRAParameters> & Pv, bool verbose)
{
  std::stringstream tss;

  //verbose=true;

  if(verbose){
    if(strlen(qm)>0) {
      cout << "Using quickmode switch " << qm << " : ";
    }
    cout << params << endl;
  }
  tss << params;
  parse(tss, Pv, verbose);
}

void MIRAParameters::parseQuickmodeNoTechSettingsChange(const char * params, const char * qm, std::vector<MIRAParameters> & Pv, bool verbose)
{
  std::vector<bool> save;
  saveParsedSettingsValues(Pv, save);
  parseQuickmode(params, qm, Pv, verbose);
  restoreParsedSettingsValues(Pv, save);
}

int32 MIRAParameters::gimmeAnInt(FlexLexer * lexer, std::stringstream & errstream)
{
  std::string currenttoken=lexer->YYText();
  if(lexer->yylex() != MP_INT){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\t\tToken '" << currenttoken << "'\n*\tExpected a number after this, not '" << lexer->YYText() <<"'\n\n";
    MP_errorinparams=true;
  }

  CEBUG("\t\tInt: " << lexer->YYText() << endl);
  return atoi(lexer->YYText());
}

double MIRAParameters::gimmeADouble(FlexLexer * lexer, std::stringstream & errstream)
{
  std::string currenttoken=lexer->YYText();
  int yyretcode=lexer->yylex();
  if(yyretcode != MP_INT && yyretcode != MP_FLOAT ){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\t\tToken '" << currenttoken << "'\n*\tExpected a number (int or float) after this, not '" << lexer->YYText() <<"'\n\n";
    MP_errorinparams=true;
  }
  CEBUG("\t\tDouble: " << lexer->YYText() << endl);
  return atof(lexer->YYText());
}

int32 MIRAParameters::getFixedStringMode(FlexLexer * lexer, std::stringstream & errstream)
{
  std::string currenttoken=lexer->YYText();
  int32 tmp=lexer->yylex();
  if(tmp==MP_UNRECOGNISED_STRING){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\t\tToken '" << currenttoken << "'\n*\tNon recognised string '" << lexer->YYText() << "', probably expected something\n*\tlike yes|no|on|off|true|false|y|n|t|f or other fixed strings.\n\n";
    MP_errorinparams=true;
  }
  CEBUG("\t\tFixed string: " << lexer->YYText() << " -> " << tmp << endl);
  return tmp;
}

//void MIRAParameters::checkCOMMON(MIRAParameters * actpar, std::vector<MIRAParameters> & Pv, FlexLexer * lexer, stringstream & errstream, const string & laststset){
//  if(!Pv.empty()){
//    if(actpar!=&Pv[0]){
//      errstream << "* Parameter section: '" << MP_currentparametersection << "'\n*\tParameter '" << lexer->YYText() << "' can only be set as COMMON_SETTINGS, not individually\n*\tfor a specific sequencing type (" << laststset << ").\n\n";
//      MP_errorinparams=true;
//    }
//  }
//}

void MIRAParameters::checkCOMMON(const std::string & currentst, FlexLexer * lexer, std::stringstream & errstream){
  if(currentst != "COMMON_SETTINGS"
     && currentst != "_COMMON_SETTINGS"){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\n*\tParameter '" << lexer->YYText() << "' can only be set as COMMON_SETTINGS, not individually\n*\tfor a specific sequencing type (" << currentst << ").\n\n";

    CEBUG("Argh!!! " << MP_currentparametersection << " " << lexer->YYText() << " not COMMON but " << currentst << endl);

    MP_errorinparams=true;
  }
}
void MIRAParameters::checkNONCOMMON(const std::string & currentst, FlexLexer * lexer, std::stringstream & errstream){
  if(currentst == "COMMON_SETTINGS"
     || currentst == "_COMMON_SETTINGS"){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\n*\tParameter '" << lexer->YYText() << "' can only be set as sequencing type specific"
      "\n*\tparameter (SANGER_SETTINGS, 454_SETTINGS, etc.pp)\n*\tand not for " << currentst << ".\n\n";

    CEBUG("Argh!!! " << MP_currentparametersection << " " << lexer->YYText() << " not NONCOMMON but " << currentst << endl);

    MP_errorinparams=true;
  }
}


std::string MIRAParameters::createAllTechString(const std::string & settings)
{
  std::string ret=
    "\n_SANGER_SETTINGS\n\t    " + settings
    + "\n_454_SETTINGS\n\t    " + settings
    + "\n_IONTOR_SETTINGS\n\t    " + settings
    + "\n_PCBIOHQ_SETTINGS\n\t    " + settings
    + "\n_PCBIOLQ_SETTINGS\n\t    " + settings
    + "\n_TEXT_SETTINGS\n\t    " + settings
    + "\n_SOLEXA_SETTINGS\n\t    " + settings
    + "\n_SOLID_SETTINGS\n\t    " + settings
    + "\n";
    ;
  return ret;
}

static const std::string noclipping_string = "-CL:pec=no:pvlc=no:qc=no:bsqc=no:mbc=no:lccf=no:lccb=no:emlc=no:emrc=no:c3pp=no:cpat=no:mqtfer=0:ckar=no:rkm=no:cbse=no";
static const std::string noquality_string = "-LR:wqf=no -AS:epoq=no -CL:qc=no:bsqc=no";


// parses either into MIRAParameters Pv vector
//  or into single MIRAParameters object
void MIRAParameters::parse(std::istream & is, std::vector<MIRAParameters> & Pv, bool verbose)
{
  FUNCSTART("void MIRAParameters::parse(std::istream & is, std::vector<MIRAParameters> & Pv)");

  std::stringstream errstream;
  FlexLexer* lexer = new MPFlexLexer(&is);

  std::string currentseqtypesettings="COMMON_SETTINGS";

  MP_errorinparams=false;
  MP_currentparametersection="(none)";

  MIRAParameters * actpar;
  if(!Pv.empty()){
    actpar=&Pv[ReadGroupLib::SEQTYPE_SANGER];
  }else{
    MIRANOTIFY(Notify::INTERNAL, "No parameter object to parse into? Not good.");
  }

  // vector to hold quickswitch --job definitions
  std::vector<std::vector<uint32> > jobdefs(4);

  int yyretcode=-1;
  while(yyretcode!=0){
    yyretcode=lexer->yylex();
    CEBUG("******: " << currentseqtypesettings << " " << yyretcode << "\t" << lexer->YYText() << endl);
    switch(yyretcode){
    case 0: {break;}                              // do nothing, eof
    case MP_PARSING_SECTIONCHANGE: {
      MP_currentparametersection=lexer->YYText();
      break;
    }
    case MP_PARSING_SECTIONRESET: {
      MP_currentparametersection="(none)";
      break;
    }
    case MP_as_numthreads:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_numthreads=gimmeAnInt(lexer,errstream);
      if(actpar->mp_assembly_params.as_numthreads==0){
	actpar->mp_assembly_params.as_numthreads=MachineInfo::getCoresTotal();
	if(actpar->mp_assembly_params.as_numthreads==0) actpar->mp_assembly_params.as_numthreads=2;
      }
      actpar->mp_skim_params.sk_numthreads=actpar->mp_assembly_params.as_numthreads;
      break;
    }
    case MP_as_automemmanagement:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_automemmanagement=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_amm_keeppercentfree:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_amm_keeppercentfree=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_amm_maxprocesssize:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_amm_maxprocesssize=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_ed_sdbg_readedit:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_edit_params.ed_sdbg_readedit=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_ed_mira_automatic_contig_editing:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_edit_params.ed_mira_automatic_contic_editing=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_ed_kmer_singlets:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_edit_params.ed_kmer_singlets=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_ed_homopolymer_overcalls:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_edit_params.ed_homopolymer_overcalls=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_ed_strict:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_edit_params.ed_strict_editing_mode=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_ed_confirmation_threshold:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_edit_params.ed_confirmation_threshold=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sp_est_startstep:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.sp_est_startstep=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_nw_check_nfs:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_nfs=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_multipassmapping:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_multipassmapping=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_templateproblems:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_templateproblems=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_duplicatereadnames:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_duplicatereadnames=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_illuminajunkinest:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_illuminajunkinest=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_proposedendclip:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_proposedendclipping=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_readnamesra:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_readnamesra=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_maxreadnamelength:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_maxreadnamelength=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_mrnlvalue:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_mrnlvalue=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_nw_check_coverage:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_coverage=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_nw_check_covvalue:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_nagandwarn_params.nw_check_covvalue=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_nw_warn_x174prefix:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      if(lexer->YYText()==nullptr){
	errstream << "ERROR Nag and Warn PhiX contig prefix not found?\n";
	MP_errorinparams=true;
	break;
      }
      actpar->mp_nagandwarn_params.nw_warn_x174prefix=lexer->YYText();
      break;
    }
    case MP_mi_extendedlog:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_extended_log=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_mi_iknowwhatido:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_iknowwhatido=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_mi_as_largecontigsize:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_as_largecontigsize=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_mi_as_largecontigsize4stats:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_as_largecontigsize4stats=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_mi_extra_flag1:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_extra_flag1=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_mi_extra_flag2:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_extra_flag2=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_mi_extra_flag3:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_special_params.mi_extra_flag3=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_kmerseries:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      int32 tmp=lexer->yylex();
      if(lexer->YYText()==nullptr){
	errstream << "ERROR -AS:kmerseries: ahum?\n";
	MP_errorinparams=true;
	break;
      }
      actpar->mp_assembly_params.as_user_kmerseries=lexer->YYText();
      actpar->mp_assembly_params.as_bphseries.clear();

      boost::char_separator<char> separator(",");
      boost::tokenizer<boost::char_separator<char> > tok(actpar->mp_assembly_params.as_user_kmerseries,separator);
      for(auto & te : tok){
	std::string tmps(te);
	boost::trim(tmps);
	int64 tmp=atoi(tmps.c_str());
	if(tmp<=0){
	  errstream << "ERROR -AS:kmerseries: " << actpar->mp_assembly_params.as_user_kmerseries << "\n" << tmps << " is <= 0. kmer values must be >0\n";
	  MP_errorinparams=true;
	  break;
	}
	actpar->mp_assembly_params.as_bphseries.push_back(static_cast<uint32>(tmp));
      }
      if(!actpar->mp_assembly_params.as_bphseries.empty()){
	actpar->mp_assembly_params.as_numpasses=0;
      }

      break;
    }
    case MP_as_extend_reads:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_use_read_extension=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_cleanup_tmp_files:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_cleanup_tmpfiles=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_ssahamerge_gapsize:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ssahamerge_gapsize=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_ssahamerge_maxfrontgap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ssahamerge_maxfrontgap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_ssahamerge_maxendgap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ssahamerge_maxendgap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_ssahamerge_strictfrontclip:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ssahamerge_strictfrontclip=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_ssahamerge_strictendclip:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ssahamerge_strictendclip=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_possible_vectors: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_possible_vectors=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_vector_maxlenallowed: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_vector_maxlenallowed=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_ensureminimumleftclip: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ensureminimumleftclipoff=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_minimumleftcliprequired: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_minslrequired=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_setminimumleftclip: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_minqlsetto=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_ensureminimumrightclip: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_ensureminimumrightclipoff=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_minimumrightcliprequired: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_minsrrequired=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_setminimumrightclip: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_minqrsetto=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_quality_minthreshold:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_quality_minthreshold=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_quality_numminthreshold:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_quality_numminthreshold=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_quality: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_quality=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_quality_minqual:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_quality_minqual=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_quality_winlen:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_quality_winlen=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_badstretchquality: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_badstretchquality=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_badstretchquality_minqual:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_badstretchquality_minqual=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_badstretchquality_winlen:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_badstretchquality_winlen=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_maskedbases: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_maskedbases=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_maskedbases_gapsize:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_maskedbase_gapsize=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_maskedbases_maxfrontgap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_maskedbase_maxfrontgap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_maskedbases_maxendgap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_maskedbase_maxendgap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_lowercase_front: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_lowercase_front=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_lowercase_back: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_lowercase_back=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_sdbg_chimeradetection: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_sdbg_chimeradetection=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_kmer_junkdetection: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_kmer_junkdetection=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_kmer_junkkill: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_kmer_junkkill=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_skimchimeradetection: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_skimchimeradetection=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_skimjunkdetection: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_skimjunkdetection=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_proposeendclips: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_proposeendclips=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_continuous: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_continuous=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_sxaggcxg: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_sxaggcxg=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_basesperhash:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_basesperhash=gimmeAnInt(lexer,errstream);
      break;
    }

    case MP_as_clip_pec_ffreq:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_ffreq=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_bfreq:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_bfreq=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_mkfr: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_mkfr=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_mtk: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_mtk=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_ffr: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_ffr=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_bfr: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_bfr=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_fcmst: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_fcmst=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_bcmst: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_bcmst=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_fsalp: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_fsalp=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_pec_bsalp: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_pec_bsalp=getFixedStringMode(lexer,errstream);
      break;
    }

    case MP_as_clip_badsolexaends: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_badsolexaends=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_search_phix174: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_search_phix174=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_filter_phix174: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_filter_phix174=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_filter_rrna: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_filter_rrna=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_filter_rrna_pairs: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_filter_rrna_pairs=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_filter_rrna_numkmers: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_filter_rrna_numkmers=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_knownadaptorsright: {
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_knownadaptorsright=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_name_prefix: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      if(lexer->YYText()==nullptr){
	errstream << "ERROR -CO:name_prefix: string not found?\n";
	MP_errorinparams=true;
	break;
      }
      actpar->mp_contig_params.con_nameprefix=lexer->YYText();
      break;
    }
    case MP_as_mark_repeats:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_mark_repeats=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_mark_repeats_only_in_result:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_mark_repeats_onlyinresult=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_assume_snp_insteadof_rmb:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_assume_snp_insteadof_rmb=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_spoiler_detection:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_spoilerdetection=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_sd_lastpassonly:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_spdetect_lastpassonly=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_uniform_read_distribution:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_uniform_read_distribution=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_urd_startinpass:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_urd_startinpass=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_urd_cutoffmultiplier:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_urd_cutoffmultiplier=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_as_rle_reads:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_rle_reads=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_automatic_repeat_detection:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_automatic_repeat_detection=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_ard_multicopythreshold:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_ard_multicopythreshold=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_as_ard_multicopyminlen:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_ard_multicopyminlen=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_ard_multicopygrace:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_ard_multicopygrace=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_polyat:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_polyat=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_polyat_keeppolystretch:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_polyat_keeppolystretch=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_polyat_len:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_polyat_len=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_polyat_maxerrors:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_polyat_maxerrors=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_polyat_maxgap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_polyat_maxgap=gimmeAnInt(lexer,errstream);
      break;
    }

    case MP_as_clip_3ppolybase:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_3ppolybase=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_clip_3ppolybase_len:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_3ppolybase_len=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_3ppolybase_maxerrors:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_3ppolybase_maxerrors=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clip_3ppolybase_maxgap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clip_3ppolybase_maxgap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_clipmask_rarekmers:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_clipmask_rarekmers=getFixedStringMode(lexer,errstream);
      break;
    }

    case MP_as_readextension_window_len:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_readextension_window_len=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_readextension_window_maxerrors:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_readextension_window_maxerrors=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_readextension_firstpassnum:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_readextension_firstpassnum=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_readextension_lastpassnum:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_readextension_lastpassnum=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_filecheck_only:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_filecheck_only=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_backbone_trimoverhangingreads:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_backbone_trimoverhangingreads=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_backbone_bootstrapnewbackbone:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_backbone_bootstrapnewbackbone=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_backbone_strainname_forceforall: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_backbone_strainname_forceforall=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_backbone_rail_fromstrain: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      if(lexer->YYText()==nullptr){
	errstream << "ERROR backbone rail from strain: name not found?\n";
	MP_errorinparams=true;
	break;
      }
      actpar->mp_assembly_params.as_backbone_rail_fromstrain=lexer->YYText();
      break;
    }
    case MP_as_startbackboneusage_inpass:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_startbackboneusage_inpass=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_backbone_raillength:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_backbone_raillength=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_backbone_railoverlap:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_backbone_railoverlap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_backbone_outlen:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_backbone_outlen=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_buntify_reads:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_buntify_reads=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_preprocess_only:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_assemblyjob_preprocessonly=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_savesimplesingletsinproject:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_savesimplesingletsinproject=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_savetaggedsingletsinproject:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_savetaggedsingletsinproject=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_enforce_qualsinreads:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_enforce_qualsinreads=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_wants_qualityfile:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_wants_qualityfile=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_output_html_cpl:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_output_html_cpl=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_output_text_cpl:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_output_text_cpl=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_output_html_gapfill: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      if(lexer->YYText()==nullptr){
	errstream << "-OUT:hegfc: character not found?\n";
	MP_errorinparams=true;
	break;
      }
      std::string tmp=lexer->YYText();
      if(tmp.size()>1) tmp=tmp[tmp.size()-2];
      actpar->mp_contig_params.con_output_html_gapfill=tmp[0];
      break;
    }
    case MP_con_output_text_gapfill: {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      if(lexer->YYText()==nullptr){
	errstream << "-OUT:tegfc: character not found?\n";
	MP_errorinparams=true;
	break;
      }
      std::string tmp=lexer->YYText();
      if(tmp.size()>1) tmp=tmp[tmp.size()-2];
      actpar->mp_contig_params.con_output_text_gapfill=tmp[0];
      break;
    }
    case MP_as_output_removerollovertmps:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_removerollovertmps=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_removetmpdir:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_removetmpdir=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_html:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_html=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_tcs:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_tcs=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_text:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_txt=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_caf:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_caf=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_maf:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_maf=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_ace:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_ace=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_gap4da:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_gap4da=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tmp_fasta:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tmp_fasta=getFixedStringMode(lexer,errstream);
      break;
    }

    case MP_as_output_exttmp_html:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_exttmp_html=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_exttmp_caf:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_exttmp_caf=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_exttmp_ace:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_exttmp_ace=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_exttmp_gap4da:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_exttmp_gap4da=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_exttmp_fasta:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_exttmp_fasta=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_exttmp_alsosinglets:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_exttmp_alsosinglets=getFixedStringMode(lexer,errstream);
      break;
    }

    case MP_as_output_html:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_html=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_text:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_txt=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_ace:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_ace=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_gff3:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_gff3=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_wiggle:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_wiggle=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_tcs:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_tcs=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_gap4da:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_gap4da=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_caf:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_caf=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_maf:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_maf=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_output_fasta:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_output_fasta=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_numpasses:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_numpasses=gimmeAnInt(lexer,errstream);
      break;
    }
//    case MP_as_skimeachpass:{
//      checkCOMMON(currentseqtypesettings, lexer, errstream);
//      actpar->mp_assembly_params.as_skimeachpass=getFixedStringMode(lexer,errstream);
//      break;
//    }
    case MP_as_numrmbbreakloops:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_numrmbbreakloops=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_maxcontigsperpass:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_maxcontigsperpass=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_minimum_readlength:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_minimum_readlength=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_as_minimum_readspercontig:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_minimum_readspercontig=gimmeAnInt(lexer,errstream);
      break;
    }
//    case MP_ads_extra_mismatch_penalty:{
//      int32 tmp=getFixedStringMode(lexer,errstream);
//      if(tmp==MP_UNRECOGNISED_STRING){
//	errstream << "-AL:emp" << endl;
//	break;
//      }
//      actpar->mp_align_params.ads_extra_mismatch_penalty=(tmp>0);
//      break;
//    }
//    case MP_ads_emp_windowlen:{
//      actpar->mp_align_params.ads_emp_windowlen=gimmeAnInt(lexer,errstream);
//      break;
//    }
//    case MP_ads_emp_maxmismatches:{
//      actpar->mp_align_params.ads_emp_maxmismatches=gimmeAnInt(lexer,errstream);
//      break;
//    }
    case MP_ads_enforce_clean_ends:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.ads_enforce_clean_ends=(getFixedStringMode(lexer,errstream)>0);
      break;
    }
    case MP_ads_clean_end_distance:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.ads_clean_end_distance=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_ads_clean_end_mismatchallowed:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.ads_clean_end_mismatchallowed=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_ads_extra_gap_penalty:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.ads_extra_gap_penalty=(getFixedStringMode(lexer,errstream)>0);
      break;
    }
    case MP_ads_max_gppercent:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.ads_max_gppercent=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_ads_gp_functionstring:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      int32 tmp=lexer->yylex();
      if(lexer->YYText()==nullptr){
	errstream << "ERROR -AL:egpl: ahum?\n";
	MP_errorinparams=true;
	break;
      }
      MP_errorinparams|=actpar->setAlignGapPenaltyLevel(lexer->YYText(),&errstream);
      break;
    }
    case MP_sk_numthreads:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_numthreads=gimmeAnInt(lexer,errstream);
      if(actpar->mp_skim_params.sk_numthreads==0){
	actpar->mp_skim_params.sk_numthreads=MachineInfo::getCoresTotal();
	if(actpar->mp_skim_params.sk_numthreads==0) actpar->mp_skim_params.sk_numthreads=2;
      }
      break;
    }
    case MP_sk_alsoskimrevcomp:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_alsoskimrevcomp=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_sk_swcheckonbackbones:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_swcheckonbackbones=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_sk_basesperhash:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_basesperhash=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_bph_increasestep:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_bph_increasestep=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_bph_max:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_bph_max=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_hashsavestepping:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_hashsavestepping=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_percentrequired:{
      actpar->mp_skim_params.sk_percentrequired=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_maxhitsperread:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_maxhitsperread=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_maxhashesinmemory:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_maxhashesinmem=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_memcaphitreduction:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_memcaphitreduction=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_filtermegahubs:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_filtermegahubs=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_sk_megahubcap:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_maxmegahubratio=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_sk_maxmegahubratio:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_skim_params.sk_maxmegahubratio=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_hs_masknastyrepeats:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_masknastyrepeats=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_hs_nastyrepeatratio:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_nastyrepeatratio=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_hs_nastyrepeatcoverage:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_nastyrepeatcoverage=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_hs_repeatlevel_in_infofile:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_repeatlevel_in_infofile=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_hs_freq_covestmin:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_freq_covestmin=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_hs_freqest_minnormal:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_freqest_minnormal=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_hs_freqest_maxnormal:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_freqest_maxnormal=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_hs_freqest_repeat:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_freqest_repeat=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_hs_freqest_heavyrepeat:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_freqest_heavyrepeat=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_hs_freqest_crazyrepeat:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_freqest_crazyrepeat=gimmeADouble(lexer,errstream);
      break;
    }
    case MP_hs_applydigitalnormalisation:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_apply_digitalnormalisation=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_hs_memtouse:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_memtouse=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_hs_rare_kmer_final_kill:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_hashstatistics_params.hs_rare_kmer_final_kill=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_al_bip:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.al_kpercent=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_al_bmin:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.al_kmin=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_al_bmax:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.al_kmax=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_al_min_overlap:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.al_min_overlap=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_al_min_score:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.al_min_score=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_al_min_relscore:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_align_params.al_min_relscore=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_fm_active:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_finalmap_params.fm_active=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_fm_mapperfect:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_finalmap_params.fm_mapperfect=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_fm_maxtotalerrors:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      int32 tmp=gimmeAnInt(lexer,errstream);
      for(uint32 p=0; p< Pv.size(); p++) {
	Pv[p].mp_finalmap_params.fm_maxtotalerrors=tmp;
      }
      break;
    }
    case MP_fm_maxgaps:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      int32 tmp=gimmeAnInt(lexer,errstream);
      for(uint32 p=0; p< Pv.size(); p++) {
	Pv[p].mp_finalmap_params.fm_maxgaps=tmp;
      }
      break;
    }
    case MP_fm_maxmismatches:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      int32 tmp=gimmeAnInt(lexer,errstream);
      for(uint32 p=0; p< Pv.size(); p++) {
	Pv[p].mp_finalmap_params.fm_maxmismatches=tmp;
      }
      break;
    }
    case MP_fm_clean_end_dist:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      int32 tmp=gimmeAnInt(lexer,errstream);
      for(uint32 p=0; p< Pv.size(); p++) {
	Pv[p].mp_finalmap_params.fm_clean_end_dist=tmp;
      }
      break;
    }
    case MP_con_analyse_mode:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_danger_analyse_mode=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_rodirs:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_reject_on_drop_in_relscore=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_min_relscore:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_min_relscore=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_also_mark_gap_bases:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_also_mark_gap_bases=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_also_mark_gap_bases_evenmc:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_also_mark_gap_bases_evenmc=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_also_mark_gap_bases_needbothstrands:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_also_mark_gap_bases_needbothstrands=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_min_groupqual_for_rmb_tagging:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_mingroupqualforrmbtagging=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_min_coverage_percentage:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_mincoveragepercentage=gimmeAnInt(lexer,errstream);
      break;
    }
    //case MP_con_min_groupqual_for_wrmbsrmb_change:{
    //  actpar->mp_contig_params.con_min_groupqual_srmbwrmb_change=gimmeAnInt(lexer,errstream);
    //  break;
    //}
    //case MP_con_num_rmb_zone_trigger:{
    //  actpar->mp_contig_params.con_rmb_numzone_trigger=gimmeAnInt(lexer,errstream);
    //  break;
    //}
    case MP_con_min_readspergroup:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_minreadspergroup=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_min_rmb_neighbourqual:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_minrmbneighbourqual=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_endread_mark_exclusion_area:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_endreadmarkexclusionarea=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_emea_setzero_on_clipping_pec:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_emea_setzero_on_clipping_pec=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_gap_override_ratio:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_gap_override_ratio=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_danger_max_error_rate:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_danger_max_error_rate=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_force_nonIUPACconsensus:{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_force_nonIUPACconsensus=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_mergeshortreads:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_mergeshortreads=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_con_msr_maxerrors:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_msr_maxerrors=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_con_msr_keependsunmapped:{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_contig_params.con_msr_keependsunmapped=gimmeAnInt(lexer,errstream);
      break;
    }

    case MP_paf_use_genomic_pathfinder :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_use_genomic_algorithms=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_paf_use_emergency_blacklist :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_use_emergency_blacklist=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_paf_use_emergency_search_stop :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_use_emergency_search_stop=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_paf_ess_partnerdepth :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_ess_depth=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_paf_use_max_contig_buildtime :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_use_max_contig_buildtime=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_paf_buildtime_inseconds :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_max_contig_buildtime=gimmeAnInt(lexer,errstream);
      actpar->setPathfinderMaxContigTime(actpar->mp_pathfinder_params.paf_max_contig_buildtime);
      break;
    }
    case MP_paf_max_startcache_filltime :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_max_startcache_filltime=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_paf_skip_whole_contig_scan :{
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_skipwholecontigscan=static_cast<int16>(getFixedStringMode(lexer,errstream));
      break;
    }
    case MP_paf_use_quick_rule :{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_use_quick_rule=static_cast<int16>(getFixedStringMode(lexer,errstream));
      break;
    }
    case MP_paf_quickrule_minlen1 :{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_quickrule_minlen1=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_paf_quickrule_minsim1 :{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_quickrule_minsim1=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_paf_quickrule_minlen2 :{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_quickrule_minlen2=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_paf_quickrule_minsim2 :{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_quickrule_minsim2=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_paf_bbquickoverlap_minlen :{
      checkNONCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_pathfinder_params.paf_bbquickoverlap_minlen=gimmeAnInt(lexer,errstream);
      break;
    }
    case MP_fn_cafout: {
      if(lexer->YYText()==nullptr){
	errstream << "ERROR Filename cafout: filename not found?\n";
	MP_errorinparams=true;
	break;
      }
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_outfile_CAF=lexer->YYText();
      break;
    }
    case MP_dir_tmp: {
      if(lexer->YYText()==nullptr){
	errstream << "ERROR Directory name tmp: name not found?\n";
	MP_errorinparams=true;
	break;
      }
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_directory_params.dir_tmp=lexer->YYText();
      break;
    }
    case MP_dir_tmp_redirectedto: {
      if(lexer->YYText()==nullptr){
	errstream << "ERROR Directory name tmp_redirected_to: name not found?\n";
	MP_errorinparams=true;
	break;
      }
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_directory_params.dir_tmp_redirectedto=lexer->YYText();
      break;
    }
    case MP_quickmode_loadparam: {
      if(lexer->YYText()==nullptr){
	errstream << "ERROR --params=: filename not found?\n";
	MP_errorinparams=true;
	break;
      }
      loadParams(lexer->YYText(), Pv);
      break;
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    case MP_quickmode_noquality_all: {
      std::string modestring = createAllTechString(noquality_string);
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noquality", Pv, verbose);
      break;
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    case MP_quickmode_noclipping_all: {
      std::string modestring = createAllTechString(noclipping_string);
      modestring+="\nCOMMON_SETTINGS -CL:ascdc=no:asjdc=no:gbcdc=no:kjd=no:kjck=no:frrna=no";
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    case MP_quickmode_noclipping_sanger: {
      std::string modestring =
	"\n_SANGER_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_454: {
      std::string modestring =
	"\n_454_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_iontor: {
      std::string modestring =
	"\n_IONTOR_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_pacbiohq: {
      std::string modestring =
	"\n_PCBIOHQ_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_pacbiolq: {
      std::string modestring =
	"\n_PCBIOLQ_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_text: {
      std::string modestring =
	"\n_TEXT_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_solexa: {
      std::string modestring =
	"\n_SOLEXA_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
    case MP_quickmode_noclipping_solid: {
      std::string modestring =
	"\n_SOLID_SETTINGS\n\t    " + noclipping_string;
      parseQuickmodeNoTechSettingsChange(modestring.c_str(), "-noclipping", Pv, verbose);
      break;
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    case MP_quickmode_hirep_something:
    case MP_quickmode_hirep_good:
    case MP_quickmode_hirep_best: {
      if(verbose) cout << "Changing some parameters for highly repetitive data, mode ";
      std::string modestring("\n_COMMON_SETTINGS"
			"\n\t-CO:mr=yes:mroir=false");
      if(yyretcode==MP_quickmode_hirep_something){
	if(verbose) cout << "'something'";
	modestring+=
	  "\n\t-KS:mnr=yes:nrr=10:nrc=0"
	  "\n";
      }else if(yyretcode==MP_quickmode_hirep_good){
	if(verbose) cout << "'good'";
	modestring+=
	  "\n\t-KS:mnr=yes:nrr=200:nrc=0"
	  "\n";
      }else if(yyretcode==MP_quickmode_hirep_best){
	if(verbose) cout << "'best'";
	modestring+=
	  "\n\t-KS:mnr=no"
	  "\n";
      }else{
	MIRANOTIFY(Notify::INTERNAL,"Huh? Unknown highly repetitive flag?");
      }
      parseQuickmode(modestring.c_str(), "", Pv, verbose);
      if(verbose) cout << "  - fixed settings\n";

      MIRAParameters * tmpactpar;
      if(!Pv.empty()) {
	tmpactpar=&Pv[0];
      }else{
	tmpactpar=actpar;
      }
      if(tmpactpar->mp_assembly_params.as_numpasses<6){
	if(verbose) cout << "  - increassing number of passes (-AS:nop) ";
	tmpactpar->mp_assembly_params.as_numpasses++;
	if(tmpactpar->mp_assembly_params.as_numpasses<=5){
	  tmpactpar->mp_assembly_params.as_numpasses++;
	  if(verbose) cout << "by two.\n";
	}else{
	  if(verbose) cout << "by one.\n";
	}
      }
      if(tmpactpar->mp_assembly_params.as_numrmbbreakloops<3){
	if(verbose) cout << "  - increasing maximum of RMB break loop (-AS:rbl).\n";
	tmpactpar->mp_assembly_params.as_numrmbbreakloops++;
      }
      if(verbose) cout << "Done\n";

      break;
    }
    case MP_quickmode_highlyrepetitive: {
      errstream << "I'm sorry, '--highlyrepetitive' has been replaced by either"
	"\n  1) '--hirep_something' for \"give me something\" (fast)"
	"\n  2) '--hirep_good' for \"give me good resolve\" (slow)"
	"\n  3) '--hirep_best' for \"give me best resolve\" (very slow)"
	"\n";
      break;
    }
    case MP_quickmode_lowqualitydata: {
      if(verbose) cout << "Adjusting parameters for low quality data (-lowqualitydata):\n"
		    "  Increassing (-CO:mrpg) by 1.\n"
		    "  Switching on -CO:amgbnbs=yes\n";
      //"  Switching on read extension (-DP:ure=yes)"

      if(!Pv.empty()) {
	for(uint32 i=0; i<Pv.size(); i++){
	  Pv[i].mp_contig_params.con_minreadspergroup+=1;
	  Pv[i].mp_contig_params.con_also_mark_gap_bases_needbothstrands=true;
	  if(Pv[i].mp_contig_params.con_reject_on_drop_in_relscore>20){
	    Pv[i].mp_contig_params.con_reject_on_drop_in_relscore-=10;
	  }
	  if(Pv[i].mp_contig_params.con_min_relscore>20){
	    Pv[i].mp_contig_params.con_reject_on_drop_in_relscore-=10;
	  }
	  Pv[i].mp_contig_params.con_mingroupqualforrmbtagging-=5;

	  Pv[i].mp_align_params.ads_extra_gap_penalty=true;
	  //Pv[i].setAlignGapPenaltyLevel(0);
	  Pv[i].mp_align_params.ads_max_gppercent=100;
	}
      }else{
	actpar->mp_contig_params.con_minreadspergroup+=1;
	actpar->mp_contig_params.con_also_mark_gap_bases_needbothstrands=true;
	if(actpar->mp_contig_params.con_reject_on_drop_in_relscore>20){
	  actpar->mp_contig_params.con_reject_on_drop_in_relscore-=10;
	}
	actpar->mp_contig_params.con_mingroupqualforrmbtagging-=5;

	actpar->mp_align_params.ads_extra_gap_penalty=true;
	//actpar->setAlignGapPenaltyLevel(0);
	actpar->mp_align_params.ads_max_gppercent=100;
      }
      break;
    }
    case MP_quickmode_highqualitydata: {
      if(verbose) cout << "Adjusting parameters for high quality data (-highqualitydata):\n"
		    "  Increassing (-CL:qcmq) by 4.\n"
		    "  Increassing (-CO:mnq) by 4.\n"
		    "  Increassing (-CO:mgqrt) by 2.\n";
      // why does gcc 4.3.2 warn
      //  warning: conversion to ‘base_quality_t’ from ‘int’ may alter its value [-Wconversion]
      // in the next 6 lines with += ???
      if(!Pv.empty()) {
	for(uint32 i=0; i<Pv.size(); i++){
	  Pv[i].mp_assembly_params.as_clip_quality_minqual+=4;
	  Pv[i].mp_contig_params.con_minrmbneighbourqual+=4;
	  Pv[i].mp_contig_params.con_mingroupqualforrmbtagging+=2;
	}
      }else{
	actpar->mp_assembly_params.as_clip_quality_minqual+=4;
	actpar->mp_contig_params.con_minrmbneighbourqual+=4;
	actpar->mp_contig_params.con_mingroupqualforrmbtagging+=2;
      }
      break;
    }
    case MP_as_nodateoutput : {
      checkCOMMON(currentseqtypesettings, lexer, errstream);
      actpar->mp_assembly_params.as_dateoutput=getFixedStringMode(lexer,errstream);
      break;
    }
    case MP_as_bangonthrow : {
      getFixedStringMode(lexer,errstream);
      Notify::setBangOnThrow(true);
      cout << "\n############# Bang on throw: MIRA will raise a SigTrap when encountering INTERNAL or FATAL errors.\n\n";
      break;
    }

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    case MP_silentparams_for_common :
    case MP_silentparams_for_sanger :
    case MP_silentparams_for_454 :
    case MP_silentparams_for_iontor :
    case MP_silentparams_for_pacbiohq :
    case MP_silentparams_for_pacbiolq :
    case MP_silentparams_for_text :
    case MP_silentparams_for_solexa :
    case MP_silentparams_for_solid :
    case MP_params_for_common :
    case MP_params_for_sanger :
    case MP_params_for_454 :
    case MP_params_for_iontor :
    case MP_params_for_pacbiohq :
    case MP_params_for_pacbiolq :
    case MP_params_for_text :
    case MP_params_for_solexa :
    case MP_params_for_solid : {
      BUGIFTHROW(Pv.empty(),"Current parsing mode does not allow for changing the sequencing type.\n");
      currentseqtypesettings=lexer->YYText();

      uint32 newst=SEQTYPE_END;
      bool realst=true;
      switch(yyretcode){
      case MP_silentparams_for_common :
      case MP_params_for_common :{
	newst=SEQTYPE_SANGER;
	realst=false;
	break;
      }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      case MP_silentparams_for_sanger :
      case MP_params_for_sanger : {
      	newst=SEQTYPE_SANGER;
      	break;
      }
      case MP_silentparams_for_454 :
      case MP_params_for_454 : {
      	newst=SEQTYPE_454GS20;
      	break;
      }
      case MP_silentparams_for_iontor :
      case MP_params_for_iontor : {
      	newst=SEQTYPE_IONTORRENT;
      	break;
      }
      case MP_silentparams_for_pacbiolq :
      case MP_params_for_pacbiolq : {
      	newst=SEQTYPE_PACBIOLQ;
      	break;
      }
      case MP_silentparams_for_pacbiohq :
      case MP_params_for_pacbiohq : {
      	newst=SEQTYPE_PACBIOHQ;
      	break;
      }
      case MP_silentparams_for_solexa :
      case MP_params_for_solexa : {
      	newst=SEQTYPE_SOLEXA;
      	break;
      }
      case MP_silentparams_for_text :
      case MP_params_for_text : {
      	newst=SEQTYPE_TEXT;
      	break;
      }
      case MP_silentparams_for_solid :
      case MP_params_for_solid : {
      	newst=SEQTYPE_ABISOLID;
      	break;
      }
      default: {
	// do nothing
      }
      }
      BUGIFTHROW(newst>=Pv.size(),"Trying to switch to " << ReadGroupLib::getNameOfSequencingType(newst) << " settings, but not enough parameter objects exist?");

      actpar=&Pv[newst];

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      switch(yyretcode){
      case MP_params_for_sanger :
      case MP_params_for_454 :
      case MP_params_for_iontor :
      case MP_params_for_pacbiohq :
      case MP_params_for_pacbiolq :
      case MP_params_for_text :
      case MP_params_for_solexa :
      case MP_params_for_solid : {
	CEBUG("\nParsed some settings!\n");
	actpar->MP_parsedsomesettings=true;
	break;
      }
      default : {
	// do nothing
      }
      }
      break;
    }
    case MP_jobdefstart : {
      MP_currentparametersection=lexer->YYText();
      // do nothing else?
      break;
    }
    case MP_jobdef_draft:
    case MP_jobdef_normal:
    case MP_jobdef_accurate: {
      jobdefs[JA_QUALITY].push_back(yyretcode);
      if(yyretcode==MP_jobdef_normal){
	if(verbose) cout << "--job=normal is deprecated and will be removed in later versions of MIRA"
		      "\nPlease use only 'draft' or 'accurate'\n";
      }
      break;
    }
    case MP_jobdef_genome:
    case MP_jobdef_est:
    case MP_jobdef_fragments:
    case MP_jobdef_clustering:
    case MP_jobdef_estsnppipeline1:
    case MP_jobdef_estsnppipeline2:
    case MP_jobdef_estsnppipeline3: {
      jobdefs[JA_TYPE].push_back(yyretcode);
      break;
    }
    case MP_jobdef_denovo:
    case MP_jobdef_mapping: {
      jobdefs[JA_METHOD].push_back(yyretcode);
      break;
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    case MP_jobdef_sanger:
    case MP_jobdef_454:
    case MP_jobdef_iontor:
    case MP_jobdef_pacbiolq:
    case MP_jobdef_pacbiohq:
    case MP_jobdef_text:
    case MP_jobdef_solexa:
    case MP_jobdef_solid: {
      jobdefs[JA_TECH].push_back(yyretcode);
      break;
    }
    case MP_jobdefend : {
      interpretJobDefs(Pv, jobdefs, errstream);

      // clean the jobdefs array
      {
	size_t tmp=jobdefs.size();
	jobdefs.clear();
	jobdefs.resize(tmp);
      }
      MP_currentparametersection="(none)";
      break;
    }
    case MP_ERROR : {
      errstream << "* Parameter section: '" << MP_currentparametersection
		<< "'\n*\tunrecognised string or unexpected character: " << lexer->YYText();
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_ERROR_RENAMED_BPH_KMER : {
      errstream << "* Parameter section: '" << MP_currentparametersection
		<< "'\n*\tParameter: " << lexer->YYText()
		<< "\n*\tParameters with a name 'bph' (bases_per_hash) in MIRA 4.0.x and before were renamed to use 'kms' (kmer_size) as this term is what the world now uses. Sorry about that.";
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_ERROR_REMOVED : {
      errstream << "* Parameter section: '" << MP_currentparametersection
		<< "'\n*\tParameter: " << lexer->YYText()
		<< "\n*\tThis parameter was removed completely, sorry about that.";
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_ERROR_MOVED_SECTION_NW : {
      errstream << "* Parameter section: '" << MP_currentparametersection
		<< "'\n*\tParameter: " << lexer->YYText()
		<< "\n*\tThis was moved to section '-NW' (eventually also renamed), sorry about that.";
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_ERROR_MOVED_SECTION_KS : {
      errstream << "* Parameter section: '" << MP_currentparametersection
		<< "'\n*\tParameter: " << lexer->YYText()
		<< "\n*\tSection -HASHSTATISTICS was renamed to section '-KMERSTATISTICS'. Furthermore, parameters with the term 'hash' were renamed with the term 'kmer'. Also, a parameters were renamed or changed meaning. Please look up the new name in the documentation, sorry about that.";
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_ERROR_DASHES : {
      errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
      errstream <<  "*\tone or several dashes with a blank too much behind them: " << lexer->YYText();
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_UNRECOGNISED_STRING : {
      errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
      errstream << "*\texpected something else than: " << lexer->YYText();
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    case MP_FLOAT :
    case MP_INT :
    case MP_ANID : {
      errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
      errstream << "*\tunrecognised string or unexpected character: " << lexer->YYText();
      if(MP_errorinparams){
	errstream << "\n*\t(may be due to previous errors)";
      }
      errstream <<"\n\n";
      MP_errorinparams=true;
      break;
    }
    default:{
      errstream << "* Oooooops, internal error!\n";
      errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
      errstream <<  "*\tparsed valid parameter '" << lexer->YYText() << "' which is not handled by Parameter object\n";
      errstream << "*\tyyretcode: " << yyretcode << "\n\n";
      MP_errorinparams=true;
    }
    }
  }

  delete lexer;

  if(Pv[0].mp_special_params.sp_parse_checktechnologypresence){
    for(uint32 st=0; st<Pv.size(); ++st){
      if(Pv[st].MP_parsedsomesettings && !Pv[st].MP_jobtechsettings){
	errstream << "* Seen parameters for " << ReadGroupLib::getNameOfSequencingType(st) << ", but no readgroup has that sequencing technology?\n\n";
      }
    }
  }

  if(errstream.str().size()){
    std::vector<int> bla;
    for(uint32 i=0; i<Pv.size(); i++) bla.push_back(i);

    cout << "\n\n";

    //MIRAParameters::dumpAllParams(
    //  Pv,
    //  bla,
    //  cerr);

    cout << "========================= Parameter parsing error(s) ==========================\n\n";
    cout << errstream.str();
    cout << "===============================================================================\n";
    MIRANOTIFY(Notify::FATAL, "Error while parsing parameters, sorry.");
  }

  //cout << endl;
  Pv[0].consistencyCheck(verbose);

  for(uint32 st=1; st<Pv.size(); ++st){
    if(Pv[st].MP_parsedsomesettings
       && Pv[0].mp_skim_params.sk_basesperhash>Pv[st].mp_align_params.al_min_overlap){
      if(verbose) cout << "WARNING: -SK:kms=" << Pv[0].mp_skim_params.sk_basesperhash
		       << " is larger than -AL:mo=" << Pv[st].mp_align_params.al_min_overlap
		       << " for " << ReadGroupLib::getShortNameOfSequencingType(st) << ". Some potential overlaps will not be found.\n";
    }
  }

  // hack to get
  //   -AL:egp:egpl:megpp:ece:ced:cema
  //  copied to every seqtype from the common settings as assembly_swalign.C and contig.C
  //  use seqtype dependent params for -AL:bip:bmax:bmin:ms:mo:mrs)
  //
  // Above can make sense for read/read SW calc, but would break down in hybrid de-novo
  //  assemblies during contig build phase.
  // E.g.: 454 & Solexa. Long 454 with homopolymers error gets added (lenient on extra gaps),
  //  but that would stop Solexas being added (if strict on extra gaps).
  for(uint32 st=1; st<Pv.size(); ++st){
    Pv[st].mp_align_params.ads_extra_gap_penalty=Pv[0].mp_align_params.ads_extra_gap_penalty;
    Pv[st].mp_align_params.ads_gp_function=Pv[0].mp_align_params.ads_gp_function;
    Pv[st].mp_align_params.ads_max_gppercent=Pv[0].mp_align_params.ads_max_gppercent;
    Pv[st].mp_align_params.ads_enforce_clean_ends=Pv[0].mp_align_params.ads_enforce_clean_ends;
    Pv[st].mp_align_params.ads_clean_end_distance=Pv[0].mp_align_params.ads_clean_end_distance;
    Pv[st].mp_align_params.ads_clean_end_mismatchallowed=Pv[0].mp_align_params.ads_clean_end_mismatchallowed;
  }

  FUNCEND();
}





void MIRAParameters::interpretJobDefs(std::vector<MIRAParameters> & Pv, std::vector<std::vector<uint32> > & jobdefs, std::stringstream & errstream)
{
  FUNCSTART("void MIRAParameters::interpretJobDefs(std::vector<std::vector<uint16> > & jobdefs, std::stringstream & errstream)");

  BUGIFTHROW(Pv.size()!=SEQTYPE_END,"Pv.size()!=SEQTYPE_END ???");

  for(auto pvi=0; pvi<Pv.size(); ++pvi){
    Pv[pvi].MP_jobtechsettings=false;
  }

  if(jobdefs[JA_QUALITY].size()==0){
    cout << "Seen no assembly quality in job definition, assuming 'accurate'.\n";
    jobdefs[JA_QUALITY].push_back(MP_jobdef_accurate);
  }else if(jobdefs[JA_QUALITY].size()>1){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
    errstream <<  "*\tSeen multiple assembly qualities in job definition, pick only one\n";
    MP_errorinparams=true;
  }
  if(jobdefs[JA_TYPE].size()==0){
    cout << "Seen no assembly type in job definition, assuming 'genome'.\n";
    jobdefs[JA_TYPE].push_back(MP_jobdef_genome);
  }else if(jobdefs[JA_TYPE].size()>1){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
    errstream <<  "*\tSeen multiple assembly types in job definition, pick only one\n";
    MP_errorinparams=true;
  }
  if(jobdefs[JA_METHOD].size()==0){
    cout << "Seen no assembly method in job definition, assuming 'denovo'.\n";
    jobdefs[JA_METHOD].push_back(MP_jobdef_denovo);
  }else if(jobdefs[JA_METHOD].size()>1){
    errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
    errstream <<  "*\tSeen multiple assembly methods in job definition, pick only one\n";
    MP_errorinparams=true;
  }

  Pv[0].getNonConstAssemblyParams().as_assemblyjob_accurate=(jobdefs[JA_QUALITY].front() == MP_jobdef_accurate);
  if(jobdefs[JA_METHOD].front()==MP_jobdef_mapping){
    Pv[0].getNonConstAssemblyParams().as_assemblyjob_mapping=true;
  }

  bool isfragments=false;
  if(jobdefs[JA_TYPE].back()==MP_jobdef_fragments){
    // fragments is a specialised mode of EST
    // i.e., has all data reduction security features (diginorm etc.) etc. switched off.
    isfragments=true;
    jobdefs[JA_TYPE].back()=MP_jobdef_est;
  }

  bool isclustering=false;
  if(jobdefs[JA_TYPE].back()==MP_jobdef_clustering){
    // clustering is a specialised mode of fragments is a specialised mode of EST
    // i.e., fragments + has some clustering settings
    isfragments=true;
    isclustering=true;
    jobdefs[JA_TYPE].back()=MP_jobdef_est;
  }


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  {
    for(uint32 i=0; i<jobdefs[JA_TECH].size(); i++){
      uint32 actst=SEQTYPE_END;
      switch(jobdefs[JA_TECH][i]){
      case MP_jobdef_sanger : {
	actst=SEQTYPE_SANGER;
	break;
      }
      case MP_jobdef_454 : {
	actst=SEQTYPE_454GS20;
	break;
      }
      case MP_jobdef_iontor : {
	actst=SEQTYPE_IONTORRENT;
	break;
      }
      case MP_jobdef_pacbiolq : {
	actst=SEQTYPE_PACBIOLQ;
	break;
      }
      case MP_jobdef_pacbiohq : {
	actst=SEQTYPE_PACBIOHQ;
	break;
      }
      case MP_jobdef_text : {
	actst=SEQTYPE_TEXT;
	break;
      }
      case MP_jobdef_solexa : {
	actst=SEQTYPE_SOLEXA;
	break;
      }
      case MP_jobdef_solid : {
	actst=SEQTYPE_ABISOLID;
	break;
      }
      default: {
	// do nothing
      }
      }

      if(Pv[actst].MP_jobtechsettings){
	errstream << "* Parameter section: '" << MP_currentparametersection << "'\n";
	errstream << "*\tSeen '" << ReadGroupLib::getNameOfSequencingType(actst) << "' more than once, did you mean a different sequencing\n\ttechnology?\n";
	MP_errorinparams=true;
      }
      Pv[actst].MP_jobtechsettings=true;
    }
  }

  bool hasSHORTREADS=Pv[SEQTYPE_SOLEXA].MP_jobtechsettings | Pv[SEQTYPE_ABISOLID].MP_jobtechsettings;

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  CEBUG("Pv[SEQTYPE_SANGER].MP_jobtechsettings: " << Pv[SEQTYPE_SANGER].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_454GS20].MP_jobtechsettings: " << Pv[SEQTYPE_454GS20].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings: " << Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings: " << Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings: " << Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_SOLEXA].MP_jobtechsettings: " << Pv[SEQTYPE_TEXT].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_SOLEXA].MP_jobtechsettings: " << Pv[SEQTYPE_SOLEXA].MP_jobtechsettings << endl);
  CEBUG("Pv[SEQTYPE_ABISOLID].MP_jobtechsettings: " << Pv[SEQTYPE_ABISOLID].MP_jobtechsettings << endl);
  CEBUG("hasSHORTREADS: " << hasSHORTREADS << endl);

  // complete: 1.6 2.1

  ///////////////////////////////////////////////////////////////////////////////////
  ////  The following qual settings are for:  genome denovo

  std::string modestring=
    "\nCOMMON_SETTINGS"
    "\n\t-GE:not=0:crkf=yes:amm=yes:kpmf=15:mps=0:ppo=no"
    "\n\t-MI:ikwid=no:el=false:lcs=500:lcs4s=5000"
    "\n\t-NW:cnfs=stop:ctp=stop:cdrn=stop:cijie=stop:cpec=stop:csrn=stop:cmrnl=stop:mrnl=40:cac=stop:acv=80:x174cp=WarnILMNPhiX174_"
    "\n\t-LR:fo=no"
    "\n\t-AS:nop=0;urle=no:ugpf=yes:urd=no:urdsip=3:mcpp=0"
    "\n\t-SB:bsnffa=no"
    "\n\t-CL:peckms=17:pechsgp=yes:ascdc=yes:asjdc=no:gbcdc=yes:kjd=yes:kjck=no"
    "\n\t-CO:mr=yes:mroir=false:asir=no:mcp=0"
    "\n\t    emeas1clpec=yes:fnic=no"
    "\n\t-SK:kmsaipp=0:acrc=yes:pr=80:swcob=no:fmh=yes:mhc=150000:mmhr=0:mchr=2048"
    "\n\t-KS:ldn=no:mnr=yes:nrr=100:nrc=0:fenn=0.4:fexn=1.6:fer=1.9:fehr=8:fecr=20:fcem=0"
    "\n\t-KS:mtu=75:rkfk=1"
    "\n\t-PF:mscft=5"
    "\n\t-ED:mace=yes:eks=yes"
    "\n\t-AL:egp=yes:egpl=0,5,10,20,40,80,100"  // "1,2,4,8,16,32,64,100"
    "\n\t"
    "\n\t-OUT:rrot=yes:rtd=no"
    "\n\t-OUT:ors=no:orm=yes:ora=no:ort=no"
    "\n\t-OUT:otm=yes:otf=yes:otc=no"
    ;
  if(jobdefs[JA_METHOD].back()==MP_jobdef_denovo) {
    modestring+=
      "\n\t-CO:mcp=10"
      ;
  }


  if(Pv[SEQTYPE_SANGER].MP_jobtechsettings){
    // nothing atm
    modestring+=
      "\nCOMMON_SETTINGS"
      "\n\t-CL:pmkfr=1"
      "\n\t-OUT:org=no"
      "\nSANGER_SETTINGS"
      "\n\t-AS:mrl=80:mrpc=2"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=400:ardgl=40"
      "\n\t-DP:ure=yes:rewl=30:rewme=2:feip=0;leip=0"
      "\n\t-CL:bsqc=yes:bsqcmq=20:bsqcwl=30:mbc=yes:mbcgs=20:mbcmfg=40:mbcmeg=60:mqtfer=0:mqtfernob=0"
      "\n\t    emlc=yes:mlcr=25:smlc=30:qc=no:cpat=no:c3pp=no:lccf=no:lccb=no:ckar=no"
      "\n\t    pec=yes:pffreq=1:pbfreq=1:pffore=false:pbfore=false:pfcmst=false:pbcmst=false:pfsalp=false:pbsalp=false"
      "\n\t-CO:mrpg=2:emea=25"
      "\n\t    amgb=yes:amgbemc=yes:amgbnbs=yes"
      "\n\t-ED:ehpo=yes"
      "\n\t-AL:bip=15:bmin=25:bmax=70:mo=17:ms=30:mrs=65"
      "\n\t-PF:uqr=yes:qrml1=200:qrms1=90:qrml2:100:qrms2=95:bqoml=150"
      "\n"
      ;
  }

  if(Pv[SEQTYPE_TEXT].MP_jobtechsettings){
    modestring+=
      "\nCOMMON_SETTINGS"
      "\n\t-CL:pmkfr=1"
      "\nTEXT_SETTINGS"
      "\n\t-AS:mrl=80:mrpc=2"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=400:ardgl=40"
      "\n\t-DP:ure=no"
      "\n\t-CO:mrpg=2:emea=1"
      "\n\t    amgb=yes:amgbemc=yes:amgbnbs=yes"
      "\n\t-CL:pec=no:pffreq=0:pbfreq=0:pffore=false:pbfore=false:pfcmst=false:pbcmst=false:pfsalp=false:pbsalp=false"
      "\n\t-ED:ehpo=yes"
      "\n\t-AL:bip=15:bmin=25:bmax=70:mo=17:ms=30:mrs=65"
      "\n\t-PF:uqr=yes:qrml1=200:qrms1=90:qrml2:100:qrms2=95:bqoml=150"
      "\nTEXT_SETTINGS\n\t"
      +
      noclipping_string
      +
      "\n"
      ;
  }

  if(Pv[SEQTYPE_454GS20].MP_jobtechsettings){
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-CL:pmkfr=1:peckms=27"
      "\n\t-OUT:org=no"
      "\n\t-KS:fenn=0.4:fexn=1.6:fer=1.9:fehr=8:fecr=20"
      "\n\t-KS:rkfk=0"
      "\n\t-AL:egp=yes:egpl=0,0,100"
      "\n\t"
      "\n454_SETTINGS"
      "\n\t-AS:mrl=40:mrpc=5"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=200:ardgl=20"
      "\n\t-DP:ure=no:rewl=15:rewme=2:feip=0;leip=0"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:emlc=no:mlcr=4:smlc=4:emrc=no:mrcr=10:smrc=15:mbc=yes:mbcgs=5:mbcmfg=12:mbcmeg=12:mqtfer=0:mqtfernob=0"
      "\n\t    msvsgs=8:msvsmfg=8:msvsmeg=12:msvssfc=0:msvssec=0:c3pp=no:cpat=no:lccf=yes:lccb=yes:ckar=yes"
      // TODO: pbsalp=true is not good for Solexa due to GGCxG errors also very early in some reads
      //       therefore transplanted this to 454, but would need some checks.
      "\n\t    pec=yes:pffreq=1:pbfreq=1:pffore=false:pbfore=yes:pfcmst=false:pbcmst=yes:pfsalp=false:pbsalp=false"
      "\n\t-AL:ms=15:mo=17:mrs=70:bip=20:bmin=20:bmax=80"
      "\n\t-CO:rodirs=30:mrpg=4:mnq=20:mgqrt=25:emea=10:amgb=no"
      "\n\t    gor=66"
      "\n\t-ED:ehpo=yes"
      "\n\t-PF:uqr=yes:qrml1=80:qrms1=90:qrml2=60:qrms2=95:bqoml=80"
      "\n\t-SK:pr=80"
      ;
  }

  if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings){
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-CL:peckms=27:pmkfr=1"
      "\n\t-OUT:org=no"
      "\n\t-KS:fenn=0.4:fexn=1.6:fer=1.9:fehr=8:fecr=20"
      "\n\t-KS:rkfk=0"
      "\n\t-AL:egp=yes:egpl=0,0,100"
      "\n\t"
      "\nIONTOR_SETTINGS"
      "\n\t-AS:mrl=40:mrpc=5"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=200:ardgl=20"
      "\n\t-DP:ure=no:rewl=15:rewme=2:feip=0;leip=0"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:emlc=no:mlcr=4:smlc=4:emrc=no:mrcr=10:smrc=15:mbc=yes:mbcgs=5:mbcmfg=12:mbcmeg=12:mqtfer=0:mqtfernob=0"
      "\n\t    msvsgs=8:msvsmfg=8:msvsmeg=12:msvssfc=0:msvssec=0:c3pp=no:cpat=no:lccf=yes:lccb=yes:ckar=yes"
      //TODO: check whether like Sanger (now) or like 454 would be better
      "\n\t    pec=yes:pffreq=1:pbfreq=1:pffore=false:pbfore=false:pfcmst=false:pbcmst=false:pfsalp=false:pbsalp=false"
      "\n\t-AL:ms=15:mo=17:mrs=70:bip=20:bmin=20:bmax=80"
      "\n\t-CO:rodirs=25:mrpg=4:mnq=20:mgqrt=25:emea=10:amgb=no"
      "\n\t    gor=66"
      "\n\t-ED:ehpo=yes"
      "\n\t-PF:uqr=yes:qrml1=80:qrms1=90:qrml2=60:qrms2=95:bqoml=80"
      "\n\t-SK:pr=80"
      ;
  }

  // TODO: adapt once PacBio data analysed
  // template from 454
  // mrpg=3 for less coverage expected.
  if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings){
    // TODO: PacBioLQ and combo
    // This is for HQ only
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-CL:peckms=27:pmkfr=1"
      "\n\t-OUT:org=no"
      "\n\t-KS:fenn=0.4:fexn=1.6:fer=1.9:fehr=8:fecr=20"
      "\n\t-AL:egp=yes:egpl=0,0,100"
      "\n\t"
      "\nPCBIOHQ_SETTINGS"
      "\n\t-AS:mrl=40:mrpc=5"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=200:ardgl=20"
      "\n\t-DP:ure=no:rewl=15:rewme=2:feip=0;leip=0"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:emlc=no:mlcr=4:smlc=4:emrc=no:mrcr=10:smrc=15:mbc=yes:mbcgs=5:mbcmfg=12:mbcmeg=12"
      "\n\t    msvsgs=8:msvsmfg=8:msvsmeg=12:msvssfc=0:msvssec=0:c3pp=no:cpat=no:lccf=yes:lccb=yes:ckar=no"
      "\n\t    pec=no:pffreq=1:pbfreq=1:pffore=false:pbfore=false:pfcmst=false:pbcmst=false:pfsalp=false:pbsalp=false"
      "\n\t-AL:ms=15:mo=17:mrs=70:bip=20:bmin=20:bmax=200"
      "\n\t-CO:rodirs=30:mrpg=3:mnq=20:mgqrt=25:emea=10:amgb=no"
      "\n\t    gor=66"
      "\n\t-ED:ehpo=yes"
      "\n\t-PF:uqr=yes:qrml1=80:qrms1=90:qrml2=60:qrms2=95:bqoml=80"
      "\n\t-SK:pr=80"
      ;
  }

  if(Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings){
    // note: -CL:pec and -SK:kms set at end again to override all for PacBio LQ data
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-CL:pmkfr=1"
      "\n\t-OUT:org=no"
      "\n\t-SK:kms=6:mmhr=90"
      "\n\t-KS:fenn=0.4:fexn=1.6:fer=1.9:fehr=8:fecr=20"
      "\n\t-ED:mace=no"
      "\n\t-AL:egp=no"
      "\n\t"
      "\nPCBIOLQ_SETTINGS"
      "\n\t-AL:bip=50:bmax=1000:mrs=10"
      "\n\t-AS:mrl=100:mrpc=2"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=200:ardgl=20"
      "\n\t-DP:ure=no:rewl=15:rewme=2:feip=0;leip=0"
      "\n\t-CO:rodirs=30:mrpg=6:mnq=20:mgqrt=25:emea=10:amgb=no"
      "\n\t    gor=66"
      "\n\t-CL"
      "\n\t    pec=no:pffreq=1:pbfreq=1:pffore=false:pbfore=false:pfcmst=false:pbcmst=false:pfsalp=false:pbsalp=false"
      "\n\t-ED:ehpo=no"
      "\n\t-PF:uqr=yes:qrml1=80:qrms1=90:qrml2=60:qrms2=95:bqoml=80"
      "\n\t-SK:pr=80"
      "\n\t"
      +noclipping_string;
      ;
  }


  if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings){
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-AS:sd=no:rbl=2:urd=no:ard=yes"
      "\n\t-OUT:org=no"
      "\n\t-SK:kss=1:mchr=4096"
      "\n\t-CO:mr=yes"
      "\n\t-PF:swcs=12"
      "\n\t-CL:peckms=31:pmkfr=1"
      "\n\t-KS:fenn=0.4:fexn=1.6:fer=1.9:fehr=8:fecr=20"
      "\n\t-KS:rkfk=0"
      "\n\t-AL:egp=yes:egpl=0,0,100"
      "\n\t"
      "\nSOLEXA_SETTINGS"
      "\n\t-OUT:sssip=no:stsip=no"
      "\n\t-AS:mrl=20:mrpc=10"
      "\n\t    urdcm=1.5:ardct=2.5:ardml=300:ardgl=20"
      "\n\t-DP:ure=no"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:emlc=no:mlcr=0:smlc=0:mbc=no:mbcgs=5:mbcmfg=12:mbcmeg=12:mqtfer=5:mqtfernob=15"
      "\n\t    spx174=true"
      "\n\t    fpx174=no"
      "\n\t    msvsgs=1:msvsmfg=2:msvsmeg=2:msvssfc=0:msvssec=0:cpat=no:lccf=no:lccb=no:ckar=yes"
      "\n\t    c3pp=yes:c3ppmsl=15:c3ppmea=3:c3ppmgfe=9:cbse=yes"
      // 27.04.2012 pbsalp=true is not good for Solexa due to GGCxG errors also very early in some reads
      // 27.04.2012 pbfreq=1 is not good for Solexa due to identical GGCxG error patterns in multiple reads
      //  i.e., for front and back clipping, only forward/reverse is good
      //  (maybe pbfreq=1 for low coverage data?)
      "\n\t    pec=yes:rkm=no:pffreq=0:pbfreq=0:pffore=yes:pbfore=yes:pfcmst=yes:pbcmst=yes:pfsalp=false:pbsalp=false"
      // 29.08.2014: why do I have mo=25 here? I might want mo=17 ???
      "\n\t-AL:ms=15:mo=17:mrs=90:bip=20:bmin=20:bmax=80"
      "\n\t-CO:rodirs=30:mrpg=4:mnq=20:mgqrt=30:emea=4:amgb=yes"
      "\n\t    msr=yes:msrme=0:msrkceu=-1"
      "\n\t-ED:ehpo=no"
      "\n\t-PF:uqr=yes:qrml1=-95:qrms1=100:qrml2=-85:qrms2=100:bqoml=20"
      "\n\t-SK:pr=90"
      ;
    if(Pv[SEQTYPE_454GS20].MP_jobtechsettings){
      modestring+=
	"\n\t"
	"\nCOMMON_SETTINGS"
	"\n\t-CL:peckms=27"
	;
    }
    if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings){
      modestring+=
	"\n\t"
	"\nCOMMON_SETTINGS"
	"\n\t-CL:peckms=27"
	;
    }
    // TODO: PacBio LQ ?
    if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings){
      modestring+=
	"\n\t"
	"\nCOMMON_SETTINGS"
	"\n\t-CL:peckms=27"
	;
    }
  }

  if(Pv[SEQTYPE_ABISOLID].MP_jobtechsettings){
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t"
      "\nSOLID_SETTINGS"
      "\n\t-DP:ure=no"
      "\n\t-AS:mrl=30"
      "\n\t    urdcm=1.5:ardct=2.0:ardml=200:ardgl=20"
      ;
  }


  if(jobdefs[JA_QUALITY].front()==MP_jobdef_draft){
    CEBUG("MP_jobdef_draft\n");
    static const char g_draftonly[] =
      "\nCOMMON_SETTINGS"
      "\n\t-AS:rbl=1:sd=no:ard=yes:urd=no"
      "\n\t-SK:kss=8:pr=70:mhpr=200"
      "\n\t"
      "\n_SANGER_SETTINGS"
      "\n\t-DP:ure=no"
      "\n\t-CL:pvlc=no"
      "\n\t-AL:mo=17:ms=30:mrs=65"
      "\n\t-CO:rodirs=15"
      ;
    modestring+=g_draftonly;

    if(Pv[SEQTYPE_454GS20].MP_jobtechsettings) {
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-AS:rbl=2"
	"\n"
	"\n_454_SETTINGS"
	"\n\t-AL:mrs=70"
	"\n\t-DP:ure=no"
	"\n\t-CL:emrc=no"
	"\n\t-SK:pr=80"
	;
    }

    if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-AS:rbl=2"
	"\n"
	"\n_IONTOR_SETTINGS"
	"\n\t-AL:mrs=70:mo=21"
	"\n\t-DP:ure=no"
	"\n\t-CL:emrc=no"
	"\n\t-SK:pr=70"
	;
    }

  // TODO: adapt once PacBio data analysed
  // template from 454
    if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings || Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings) {
      // TODO: PacBio HQ only atm
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-AS:rbl=2"
	"\n"
	;
    }

    if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings) {
      // TODO: PacBio ??? LQ?
      if(Pv[SEQTYPE_SANGER].MP_jobtechsettings || Pv[SEQTYPE_454GS20].MP_jobtechsettings || Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings || Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-SK:pr=85:kss=4"
	  "\n_454_SETTINGS"
	  "\n\t-SK:pr=85"
	  "\n_IONTOR_SETTINGS"
	  "\n\t-SK:pr=85"
	  "\n_PCBIOHQ_SETTINGS"
	  "\n\t-SK:pr=85"
	  ;
      }else{
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=2"
	  "\n\t-SK:kss=4"
	  "\n_SOLEXA_SETTINGS"
	  "\n\t-SK:pr=90"
	  ;
      }
    }

  }else{
    CEBUG("not MP_jobdef_draft\n");
    static const char g_level23[]=
      "\nCOMMON_SETTINGS"
      "\n\t-AS:sd=yes:sdlpo=yes"
      "\n\t"
      "\n_SANGER_SETTINGS"
      "\n\t-DP:ure=yes:feip=0;leip=0"
      "\n\t-CL:pvlc=yes:pvcmla=18"
      ;
    modestring+=g_level23;

    if(jobdefs[JA_QUALITY].front()==MP_jobdef_normal){
      CEBUG("MP_jobdef_normal\n");
      static const char g_normalonly[] =
	"\nCOMMON_SETTINGS"
	"\n\t-AS:rbl=2:urdsip=3"
	"\n\t-SK:kss=4:pr=70:mhpr=2000"
	"\n\t"
	"\n_SANGER_SETTINGS"
	"\n\t-AL:bip=15:bmin=25:bmax=100"
	"\n\t-CO:rodirs=20"
	;
      modestring+=g_normalonly;

      if(Pv[SEQTYPE_454GS20].MP_jobtechsettings) {
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=2:urdsip=3"
	  "\n"
	  "\n_SANGER_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n"
	  "\n_454_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=80"
	  ;
      }

      if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=2:urdsip=3"
	  "\n"
	  "\n_SANGER_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n"
	  "\n_IONTOR_SETTINGS"
	  "\n\t-AL:mrs=70:mo=19"
	  "\n\t-SK:pr=50"
	  ;
      }


      // TODO: adapt once PacBio data analysed
      // template from 454
      if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings) {
	// TODO: PacBio HQ atm
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=2:urdsip=3"
	  "\n"
	  "\n_SANGER_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n"
	  "\n_PCBIOHQ_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=80"
	  ;
      }

      if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings) {
	if(Pv[SEQTYPE_SANGER].MP_jobtechsettings || Pv[SEQTYPE_454GS20].MP_jobtechsettings || Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings) {
	// TODO: PacBio HQ atm
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-SK:pr=85:kss=1"
	    "\n_454_SETTINGS"
	    "\n\t-SK:pr=85"
	    "\n_PCBIOHQ_SETTINGS"
	    "\n\t-SK:pr=85"
	    ;
	}else{
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-AS:rbl=2"
	    "\n\t-SK:kss=1:mhpr=100"
	    "\n"
	    "\n_SOLEXA_SETTINGS"
	    "\n\t-SK:pr=90"
	    ;
	}
      }

    }else{
      CEBUG("not MP_jobdef_draft ie MP_jobdef_accurate\n");
      static const char g_accurateonly[] =
	"\nCOMMON_SETTINGS"
	"\n\t-AS:rbl=2:urdsip=3"
	"\n\t-SK:kss=1:pr=65:mhpr=2000"
	"\n_SANGER_SETTINGS"
	"\n\t-AL:bip=20:bmin=25:bmax=130"
	"\n\t-CO:rodirs=25"
	;

      modestring+=g_accurateonly;

      if(Pv[SEQTYPE_454GS20].MP_jobtechsettings) {
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=3:urdsip=4"
	  "\n"
	  "\n_SANGER_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=70"
	  "\n"
	  "\n_454_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=80"
	  ;
      }


      if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=3:urdsip=4"
	  "\n"
	  "\n_SANGER_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=70"
	  "\n"
	  "\n_IONTOR_SETTINGS"
	  "\n\t-AL:mrs=70:mo=19"
	  "\n\t-SK:pr=50"
	  ;
      }

      // TODO: adapt once PacBio data analysed
      // template from 454
      if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings) {
	// TODO: PacBio HQ atm
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:rbl=3:urdsip=4"
	  "\n"
	  "\n_SANGER_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=70"
	  "\n"
	  "\n_PCBIOHQ_SETTINGS"
	  "\n\t-AL:mrs=70"
	  "\n\t-SK:pr=80"
	  ;
      }

      if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings) {
	if(Pv[SEQTYPE_SANGER].MP_jobtechsettings || Pv[SEQTYPE_454GS20].MP_jobtechsettings || Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings || Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
	  // TODO: PacBio HQ atm
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-SK:pr=85:kss=1"
	    "\n_454_SETTINGS"
	    "\n\t-SK:pr=85"
	    "\n_IONTOR_SETTINGS"
	    "\n\t-SK:pr=85"
	    "\n_PCBIOHQ_SETTINGS"
	    "\n\t-SK:pr=85"
	    ;
	}else{
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-AS:rbl=2"
	    "\n\t-SK:kss=1:mhpr=2000"
	    "\n_SOLEXA_SETTINGS"
	    "\n\t-SK:pr=95"
	    ;
	}
      }
    }
  }

  if(!jobdefs[JA_METHOD].empty() && jobdefs[JA_METHOD].front()==MP_jobdef_mapping){
    // modifier for mapping assembly

    CEBUG("MP_jobdef_mapping\n");

    addModifiersForMapping(modestring, jobdefs, Pv, hasSHORTREADS);
  }

  //
  // settings for EST / RNASeq
  //

  std::string estbasesettings=
    "\nCOMMON_SETTINGS"
    "\n\t-GE:crkf=no"
    "\n\t-MI:lcs=500:lcs4s=1000"
    "\n\t-AS:sd=no:ard=no:urd=no:ugpf=no:uess=yes:esspd=500:umcbt=yes:bts=360"
    "\n\t-CL:peckms=25:ascdc=no:asjdc=no:gbcdc=yes:kjd=yes:kjck=yes"
    "\n\t    frrna=yes:frrnap=yes:frrnank=20"
    "\n\t-SK:mhpr=30:mhc=250000"
    "\n\t-KS:ldn=yes:mnr=yes:fcem=30:nrc=200"
    "\n\t-CO:fnic=yes"               // force non-IUPAC results
    "\n\t-OUT:orw=no"
    "\n\t-AL:egp=yes"
    ;
  if(Pv[SEQTYPE_SANGER].MP_jobtechsettings){
    estbasesettings+=
      "\nSANGER_SETTINGS"
      "\n\t-AS:mrpc=2"
      "\n\t-AL:mrs=85"
      "\n\t-DP:ure=no"
      "\n\t-CL:pvlc=no:qc=yes:bsqc=no:mbc=yes:c3pp=no:emlc=no:emrc=no:mqtfer=0:mqtfernob=0"
      "\n\t    cpat=yes:cpkps=yes:cpmsl=12:cpmea=1:cpmgfe=20000"
      "\n\t    pec=no"
      "\n\t-CO:rodirs=10"
      ;
  }
  if(Pv[SEQTYPE_454GS20].MP_jobtechsettings) {
    estbasesettings+=
      "\n454_SETTINGS"
      "\n\t-AS:mrpc=2"
      "\n\t-AL:mrs=80"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:mbc=yes:c3pp=no:emlc=no:emrc=no:mqtfer=0:mqtfernob=0"
      "\n\t    cpat=yes:cpkps=yes:cpmsl=12:cpmea=1:cpmgfe=20000"
      "\n\t    pec=yes:pecc=yes"
      "\n\t-CO:rodirs=15"
      ;
  }
  if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
    estbasesettings+=
      "\nIONTOR_SETTINGS"
      "\n\t-AS:mrpc=2"
      "\n\t-AL:mrs=80"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:mbc=yes:c3pp=no:emlc=no:emrc=no:mqtfer=0:mqtfernob=0"
      "\n\t    cpat=yes:cpkps=yes:cpmsl=12:cpmea=1:cpmgfe=20000"
      "\n\t    pec=yes:pecc=yes"
      "\n\t-CO:rodirs=15"
      ;
  }
  if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings) {
    estbasesettings+=
      "\nCOMMON_SETTINGS"
      "\n\t-CL:peckms=31:pmkfr=1:pmtk=3:ascdc=yes"
      "\nSOLEXA_SETTINGS"
      "\n\t-AS:mrpc=4"
      "\n\t-AL:mrs=90"
      "\n\t-CL:pvlc=no:qc=no:bsqc=no:mbc=yes:emlc=no:emrc=no:mqtfer=5:mqtfernob=15"
      "\n\t    spx174=yes:fpx174=yes"
      "\n\t    cpat=yes:cpkps=yes:cpmsl=15:cpmea=1:cpmgfe=20000"
      "\n\t    c3pp=yes:c3ppmsl=15:c3ppmea=3:c3ppmgfe=9"
      "\n\t    pec=yes:pecc=yes:rkm=no"
      "\n\t-CO:rodirs=15"
      ;
  }
  if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings) {
    estbasesettings+=
      "\nPCBIOHQ_SETTINGS"
      "\n\t-AS:mrpc=2"
      "\n\t-AL:mrs=85"
      "\n\t-DP:ure=no"
      "\n\t-CL:pvlc=no:qc=yes:bsqc=no:mbc=yes:c3pp=no:cpat=yes:cpkps=yes:cpmsl=12:cpmea=1:cpmgfe=20000:emlc=no:emrc=no:mqtfer=0:mqtfernob=0:pec=yes:pecc=yes"
      "\n\t-CO:rodirs=10"
      ;
  }
  // TODO: PacBio LQ?
  if(Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings) {
    estbasesettings+=
      "\nPCBIOLQ_SETTINGS"
      "\n\t-AS:mrpc=2"
      "\n\t-DP:ure=no"
      "\n\t-CL:pvlc=no:qc=yes:bsqc=no:mbc=yes:c3pp=no:cpat=yes:cpkps=yes:cpmsl=12:cpmea=1:cpmgfe=20000:emlc=no:emrc=no:mqtfer=0:mqtfernob=0"
      "\n\t-CO:rodirs=30"
      ;
  }
  if(!jobdefs[JA_METHOD].empty() && jobdefs[JA_METHOD].front()==MP_jobdef_mapping){
    addModifiersForMapping(estbasesettings, jobdefs, Pv, hasSHORTREADS);
  }

  switch(jobdefs[JA_TYPE].front()){
    case MP_jobdef_genome : {
      // genome is handled in the whole part above
      // no need to set anything here
      break;
    }
    case MP_jobdef_est : {
      CEBUG("MP_jobdef_est\n");
      modestring+=estbasesettings;
      // give Sanger and 454 Sequences a bit more maximum time (longer sequences)
      if(Pv[SEQTYPE_SANGER].MP_jobtechsettings || Pv[SEQTYPE_454GS20].MP_jobtechsettings){
	modestring+=
	  "\nCOMMON_SETTINGS"
	  "\n\t-AS:bts=720"
	  ;
      }
      break;
    }
    case MP_jobdef_estsnppipeline1 : {
      CEBUG("MP_jobdef_estsnppipeline1\n");
      modestring+=estbasesettings;
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-GE:proout=step1:esps=1"
	"\n\t-AS:rbl=4"
	"\n\t-CO:mr=yes:mroir=no:asir=no"
	"\n\t-OUT:orc=yes:orf=yes:org=no:ora=no:ort=no:orh=no"
	"\n\t-ED:mace=no"
	;

      if(Pv[SEQTYPE_SANGER].MP_jobtechsettings){
	modestring+=
	  "\nSANGER_SETTINGS"
	  "\n\t-AL:ms=15:mo=17:mrs=65"
	  "\n\t-CO:emea=15:rodirs=20:amgb=yes:amgbemc=yes:amgbnbs=yes"
	  ;
      }
      // TODO: PacBio LQ?
      if(Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings) {
	modestring+=
	  "\nPCBIOLQ_SETTINGS"
	  "\n\t-AL:ms=15:mo=17:mrs=10"
	  "\n\t-CO:emea=15:rodirs=20:amgb=yes:amgbemc=yes:amgbnbs=yes"
	  ;
      }
      if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings) {
	modestring+=
	  "\nPCBIOHQ_SETTINGS"
	  "\n\t-AL:ms=15:mo=17:mrs=65"
	  "\n\t-CO:emea=15:rodirs=20:amgb=yes:amgbemc=yes:amgbnbs=yes"
	  ;
      }
      if(Pv[SEQTYPE_454GS20].MP_jobtechsettings) {
	modestring+=
	  "\n454_SETTINGS"
	  "\n\t-AL:ms=15:mo=17:mrs=65"
	  "\n\t-CO:emea=5:rodirs=25:amgb=no"
	  ;
      }
      if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) {
	modestring+=
	  "\nIONTOR_SETTINGS"
	  "\n\t-AL:ms=15:mo=17:mrs=65"
	  "\n\t-CO:emea=5:rodirs=25:amgb=no"
	  ;
      }
      break;
    }
    case MP_jobdef_estsnppipeline2 : {
      CEBUG("MP_jobdef_estsnppipeline2\n");
      modestring+=estbasesettings;
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-GE:proout=step2:esps=2"
	"\n\t-AS:rbl=7"
	"\n\t-CO:mr=yes:mroir=no:asir=no:fnic=yes"
	"\n\t-OUT:orc=yes:orf=yes:org=no:ora=no:ort=no:orh=no"
	"\n\t-ED:mace=no"
	"\n\t-AL:egp=yes:egpl=0,0,100"
	;
      if(Pv[SEQTYPE_SANGER].MP_jobtechsettings){
	modestring+=
	  "\nSANGER_SETTINGS"
	  "\n\t" + noclipping_string +
	  "\n\t-AL:ms=30:mo=30:mrs=75"
	  "\n\t-CO:emea=15:rodirs=10:amgb=yes:amgbemc=yes:amgbnbs=yes"
	  ;
      }
      if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings){
	modestring+=
	  "\nPCBIOHQ_SETTINGS"
	  "\n\t" + noclipping_string +
	  "\n\t-AL:ms=30:mo=30:mrs=75"
	  "\n\t-CO:emea=15:rodirs=10:amgb=yes:amgbemc=yes:amgbnbs=yes"
	  ;
      }
      // TODO: PacBio LQ ?
      if(Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings){
	modestring+=
	  "\nPCBIOLQ_SETTINGS"
	  "\n\t" + noclipping_string +
	  "\n\t-AL:ms=30:mo=30:mrs=10"
	  "\n\t-CO:emea=15:rodirs=10:amgb=yes:amgbemc=yes:amgbnbs=yes"
	  ;
      }
      if(Pv[SEQTYPE_454GS20].MP_jobtechsettings){
	modestring+=
	  "\n454_SETTINGS"
	  "\n\t" + noclipping_string +
	  "\n\t-AL:ms=30:mo=30:mrs=75"
	  "\n\t-CO:emea=5:rodirs=15:amgb=no"
	  ;
      }
      if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings){
	modestring+=
	  "\nIONTOR_SETTINGS"
	  "\n\t" + noclipping_string +
	  "\n\t-AL:ms=30:mo=30:mrs=75"
	  "\n\t-CO:emea=5:rodirs=15:amgb=no"
	  ;
      }
      break;
    }
    case MP_jobdef_estsnppipeline3 : {
      CEBUG("MP_jobdef_estsnppipeline3\n");
      modestring+=estbasesettings;
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-GE:proout=step3:esps=3"
	"\n\t-FN:cafin=step2_reads.caf:straindatain=step2_straindata_in.txt"
	"\n\t-AS:rbl=1"
	"\n\t-CO:mr=yes:mroir=no:asir=yes"
	"\n\t-OUT:orc=yes:orf=yes:org=no:ora=no:ort=no:orh=no"
	"\n\t-ED:mace=no"
	"\n\t-AL:egp=yes:egpl=0,0,0,0,0,0,0,0,0,0,0,0,0,0,100"  // 15 gaps, reject if 5 codons or more as gap
	;
	if(Pv[SEQTYPE_SANGER].MP_jobtechsettings){
	  modestring+=
	    "\nSANGER_SETTINGS"
	    "\n\t" + noclipping_string +
	    "\n\t-AL:ms=30:mo=30:mrs=70"
	    "\n\t-CO:rodirs=12:mrpg=1:emea=3:amgb=yes:amgbemc=yes:amgbnbs=yes"
	    ;
	}
	// TODO: PacBio LQ ?
	if(Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings){
	  modestring+=
	    "\nPCBIOLQ_SETTINGS"
	    "\n\t" + noclipping_string +
	    "\n\t-AL:ms=30:mo=30:mrs=70"
	    "\n\t-CO:rodirs=12:mrpg=1:emea=3:amgb=yes:amgbemc=yes:amgbnbs=yes"
	    ;
	}
	if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings){
	  modestring+=
	    "\nPCBIOHQ_SETTINGS"
	    "\n\t" + noclipping_string +
	    "\n\t-AL:ms=30:mo=30:mrs=70"
	    "\n\t-CO:rodirs=12:mrpg=1:emea=3:amgb=yes:amgbemc=yes:amgbnbs=yes"
	    ;
	}
	if(Pv[SEQTYPE_454GS20].MP_jobtechsettings){
	  modestring+=
	    "\n454_SETTINGS"
	    "\n\t" + noclipping_string +
	    "\n\t-AL:ms=30:mo=30:mrs=70"
	    "\n\t-CO:rodirs=12:mrpg=1:emea=3:amgb=yes:amgbemc=yes:amgbnbs=yes"
	    ;
	}
	if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings){
	  modestring+=
	    "\nIONTOR_SETTINGS"
	    "\n\t" + noclipping_string +
	    "\n\t-AL:ms=30:mo=30:mrs=70"
	    "\n\t-CO:rodirs=12:mrpg=1:emea=3:amgb=yes:amgbemc=yes:amgbnbs=yes"
	    ;
	}
      break;
    }
  default : {
  }
  }

  // fragments at end for override of diginorm etc.
  if(isfragments){
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-KS:mnr=no:ldn=no"
      "\n\t-SK:mhpr=20000:fmh=no"
      "\n\t-CO:fnic=yes"               // force non-IUPAC results
      "\n\t-CL:frrna=no"               // do not filter away rRNA
      "\n\t";
  }

  // clustering at end
  if(isclustering){
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-AS:kms=17"                          // one pass only, 17mers
      "\n\t-ED:gbre=no:mace=no"                 // no de-bruijn editing, no MIRA editing
      "\n\t-CO:asir=yes:fnic=yes"               // assume SNPs, force non-IUPAC results
      "\n\t--noclipping"                        // absolutely no clipping
      "\n\t-NW:cijie=no:cpec=no"                // turn off nag warnings for missing clipping options
      ;


    std::string alltech(
      "-AS:mrpc=1 -OUT:sssip=yes:stsip=yes"      // make sure singlets are made
      " -AS:mrl=50"                               // Assembly: min 50 read length
      " -CO:mrpg=1:mnq=1:mgqrt=1"                 // tag everything by setting min qual to 1
      " -CO:amgb=yes:amgbemc=yes:amgbnbs=no"      // tag really all gaps
      );
    if(jobdefs[JA_QUALITY].back() == MP_jobdef_accurate){
      modestring+="\n\t-AL:egp=yes:egpl=100:ece=no";   // extra gap: reject all gaps larger than 2 bp. No clean ends.
      alltech+=" -AL:mo=50:mrs=96 -CO:cmrs=96";             // alignment: min 50 overlap, (c)mrs at ~2% diff
    }else{
      modestring+="\n\t-AL:egp=yes:egpl=0,0,0,0,0,0,0,0,0,0,0,0,100:ece=no";   // extra gap: reject all gaps larger than 12 bp. No clean ends.
      alltech+=" -AL:mo=50:mrs=85 -CO:cmrs=85";             // alignment: min 50 overlap, (c)mrs at ~7.5% diff
    }

    modestring+=createAllTechString(alltech);
  }

  if(Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings){
    // note: -CL:pec and -SK:kms set at end again to override all for PacBio LQ data
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-SK:kms=6:mmhr=90"
      "\n\t";
  }


  CEBUG("Standard settings after quality, method and type lookup:" << modestring << endl);

  parseQuickmode(modestring.c_str(), "JobDefs", Pv, false);

  FUNCEND();
}


void MIRAParameters::addModifiersForMapping(std::string & modestring, std::vector<std::vector<uint32> > & jobdefs, std::vector<MIRAParameters> & Pv, bool hasSHORTREADS)
{
    modestring+=
      "\n\t"
      "\nCOMMON_SETTINGS"
      "\n\t-GE:crkf=no"
      "\n\t-AS:ard=yes:urd=no:nop=1:rbl=1"
      "\n\t-CL:ascdc=no"
      "\n\t-SB:sbuip=1:brl=0:bro=0:tor=yes"
      "\n\t-SK:kms=16:kss=4:pr=60:mhpr=1000"
      "\n\t-KS:ldn=no:mnr=no"
      "\n\t-CO:mr=yes"
      "\n\t-PF:swcs=12"
      "\n\t-AL:egp=no"
      "\n\t-FM:active=yes"
      ;

    if(jobdefs[JA_QUALITY].front()==MP_jobdef_accurate){
      modestring+=
	"\n\t"
	"\nCOMMON_SETTINGS"
	"\n\t-SK:swcob=yes";
    }

    if(Pv[SEQTYPE_SANGER].MP_jobtechsettings){
      // empty atm
      modestring+=
	"\nSANGER_SETTINGS"
	"\n\t-AS:mrl=20"
	"\n\t-SB:bnb=no"
	;
    }
    if(Pv[SEQTYPE_454GS20].MP_jobtechsettings){
      // empty atm
      modestring+=
	"\n454_SETTINGS"
	"\n\t-AS:mrl=20"
	"\n\t-SB:bnb=no"
	;
    }
    if(Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings){
      modestring+=
	"\n\t"
	"\nCOMMON_SETTINGS"
	"\n\t-NW:acv=120"
	"\nIONTOR_SETTINGS"
	"\n\t-AS:mrl=20"
	"\n\t-SB:bnb=no"
	;
    }
    if(Pv[SEQTYPE_TEXT].MP_jobtechsettings){
      modestring+=
	"\nTEXT_SETTINGS"
	"\n\t-AS:mrl=20"
	;
    }
    if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings){
      modestring+=
	"\n\t"
	"\nCOMMON_SETTINGS"
	"\n\t-SK:kms=17:kss=1:pr=60"
	"\n\t-CO:mr=yes"
	"\n\t-PF:swcs=12"
	"\n\t-NW:acv=160"
	"\n\t"
	"\nSOLEXA_SETTINGS"
	"\n\t-AS:mrl=20"
	"\n\t    urdcm=1.5:ardct=2.0:ardml=200:ardgl=20"
	"\n\t-DP:ure=no"
	"\n\t-CL:pvlc=no:qc=no:bsqc=no:emlc=no:mlcr=0:smlc=0:mbc=no:mbcgs=5:mbcmfg=12:mbcmeg=12"
	"\n\t    msvsgs=1:msvsmfg=2:msvsmeg=2:msvssfc=0:msvssec=0:cpat=no:mqtfer=0"
	"\n\t-AL:ms=15:mo=17:mrs=80:bip=20:bmin=20:bmax=80"
	"\n\t-FM:mpm=true:mte=-1:mmm=-1:mg=-1:ced=0"
	"\n\t-CO:rodirs=30:cmrs=-1:mrpg=3:mnq=20:mgqrt=30:emea=4:amgb=yes"
	"\n\t-PF:uqr=yes:qrml1=-90:qrms1=100:qrml2=-80:qrms2=100:bqoml=20"
	"\n\t-SK:pr=60"
	"\n\t-SB:bnb=yes"
	;
      if(Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings){
	// empty atm
      }
    }

    modestring+=
      "\nCOMMON_SETTINGS";

    // also take care of brl by checking which sequencing types are used
    if(hasSHORTREADS){
      // TODO: check good brl length
      modestring+=
	"\nCOMMON_SETTINGS"
	"\n\t-SB:sbuip=0:brl=0";
      // BaCh: 27.04.2013
      // mroir=yes heavily interferes with decision making on whether
      //  overcalls may be edited or not (assembly.C), as this is decided
      //  upon existence of SRM tags. And these do not get set if
      //  mroir is true, which then sometimes leads to wrong overcall editing
      //"\n\t-CO:mroir=yes";
      if(Pv[SEQTYPE_SOLEXA].MP_jobtechsettings){
	if(jobdefs[JA_QUALITY].front()==MP_jobdef_draft){
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-SB:bro=0"
	    "\n\t-SK:mhpr=1000"
	    "\n_SANGER_SETTINGS"
	    "\n\t-SK:pr=70"
	    "\n_454_SETTINGS"
	    "\n\t-SK:pr=70"
	    "\n_IONTOR_SETTINGS"
	    "\n\t-SK:pr=70"
	    "\n_PCBIOHQ_SETTINGS"
	    "\n\t-SK:pr=70"
	    // TODO:PacBio LQ ?
	    "\n_SOLEXA_SETTINGS"
	    "\n\t-AL:mrs=70"
	    "\n\t-CO:msr=yes"
	    "\n\t-SK:pr=90"
	    ;
	}else if(jobdefs[JA_QUALITY].front()==MP_jobdef_normal){
	  if(!Pv[SEQTYPE_SANGER].MP_jobtechsettings && !Pv[SEQTYPE_454GS20].MP_jobtechsettings && !Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) modestring+="\nCOMMON_SETTINGS\n-SK:kms=12:kss=1:mhpr=1500";
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-SB:brl=0:bro=0"
	    "\n_SANGER_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_454_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_IONTOR_SETTINGS"
	    "\n\t-SK:pr=60"
	    // TODO:PacBio LQ ?
	    "\n_PCBIOHQ_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_SOLEXA_SETTINGS"
	    "\n\t-AL:mrs=75"
	    "\n\t-CO:msr=yes"
	    "\n\t-SK:pr=75"
	    ;
	}else{
	  if(!Pv[SEQTYPE_SANGER].MP_jobtechsettings && !Pv[SEQTYPE_454GS20].MP_jobtechsettings && !Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings && !Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings) modestring+="\nCOMMON_SETTINGS\n-SK:kms=10:kss=1:mhpr=2000";
	  modestring+=
	    "\nCOMMON_SETTINGS"
	    "\n\t-SB:brl=0:bro=0"
	    "\n_SANGER_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_454_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_IONTOR_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_PCBIOHQ_SETTINGS"
	    "\n\t-SK:pr=60"
	    "\n_SOLEXA_SETTINGS"
	    "\n\t-AL:mrs=60"
	    "\n\t-CO:msr=yes"
	    "\n\t-SK:pr=60"
	    ;
	}
      }else{
	if(jobdefs[JA_QUALITY].front()==MP_jobdef_draft){
	  modestring+="\n\t-SB:bro=20";
	}else if(jobdefs[JA_QUALITY].front()==MP_jobdef_normal){
	  modestring+="\n\t-SB:bro=35";
	  if(!Pv[SEQTYPE_SANGER].MP_jobtechsettings && !Pv[SEQTYPE_454GS20].MP_jobtechsettings && !Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings && !Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings && !Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) modestring+="\n\n-SK:kms=12";
	}else{
	  modestring+="\n\t-SB:bro=40";
	  if(!Pv[SEQTYPE_SANGER].MP_jobtechsettings && !Pv[SEQTYPE_454GS20].MP_jobtechsettings && !Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings && !Pv[SEQTYPE_PACBIOLQ].MP_jobtechsettings && !Pv[SEQTYPE_IONTORRENT].MP_jobtechsettings) modestring+="\n\n-SK:kms=10";
	}
      }
    }else if(Pv[SEQTYPE_454GS20].MP_jobtechsettings || Pv[SEQTYPE_PACBIOHQ].MP_jobtechsettings){
      // TODO: PacBio LQ?
      modestring+="\n\t-SB:brl=0:bro=0"
	"\n\t-SB:sbuip=0";
    }else{
      modestring+="\n\t-SB:brl=0:bro=0"
	"\n\t-SB:sbuip=0";
    }
}


void MIRAParameters::saveParsedSettingsValues(std::vector<MIRAParameters> & Pv, std::vector<bool> & saved)
{
  saved.clear();
  if(!Pv.empty()) {
    saved.resize(Pv.size());
    for(uint32 i=0; i<Pv.size(); i++){
      saved[i]=Pv[i].MP_parsedsomesettings;
    }
  }
}

void MIRAParameters::restoreParsedSettingsValues(std::vector<MIRAParameters> & Pv, std::vector<bool> & saved)
{
  if(saved.size()==Pv.size()) {
    for(uint32 i=0; i<saved.size(); i++){
      Pv[i].MP_parsedsomesettings=saved[i];
    }
  }
}


/*************************************************************************
 *
 * This is so depressing: this cannot be called from a static initialiser
 *  as findLocationOfSelfBinary() calls boost::read_symlink and that
 *  crashes during __static_initialization_and_destruction_0() of the C/C++
 *  environment (at least with GCC)
 *
 *************************************************************************/

void MIRAParameters::initStdLocations()
{
  // find out where the called binary is installed
  std::string calledbin;
  findLocationOfSelfBinary(calledbin);

  MP_binpath=calledbin;

  std::string path;
  std::string miraprog;
  splitFullPathAndFileName(calledbin,path,miraprog);
  MP_bindir=path;
  MP_homedir=path+"/..";
  MP_sharedir=MP_homedir+"/share/mira";
  MP_mhslibdir=MP_sharedir+"/mhs";
  return;
}


const std::string & MIRAParameters::getBinPath() {
  if(MP_binpath.empty()) initStdLocations();
  return MP_binpath;
}
const std::string & MIRAParameters::getBinDir(){
  if(MP_binpath.empty()) initStdLocations();
  return MP_bindir;
}
const std::string & MIRAParameters::getHomeDir(){
  if(MP_binpath.empty()) initStdLocations();
  return MP_homedir;
}
const std::string & MIRAParameters::getShareDir()
{
  if(MP_binpath.empty()) initStdLocations();
  return MP_sharedir;
}
const std::string & MIRAParameters::getMHSLibDir()
{
  if(MP_binpath.empty()) initStdLocations();
  return MP_mhslibdir;
}
