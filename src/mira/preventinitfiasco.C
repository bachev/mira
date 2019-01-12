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


#include <memorc/memorc.H>

#include <io/annotationmappings.H>
#include <mira/readgrouplib.H>
#include <mira/multitag.H>
#include <mira/read.H>
#include <mira/contig.H>
#include <mira/gff_parse.H>
#include <util/misc.H>




// simplest way to avoid the static initialisation fiasco:
//  define static class variables for libraries here

/*************************************************************************
 *
 * MemORC
 *
 *************************************************************************/
#ifdef MIRAMEMORC
MemORC::mbmap_t MemORC::MOC_memblocks;
MemORC::mbmap_t MemORC::MOC_hotblocks;
std::vector<uint64> MemORC::MOC_hotaidsrequested;

MemORC MemORC::MOC_semaphore; // keep last for memorc: when instantiated, sets readytouse
                   //  when destructed, clears readytouse
#endif

/*************************************************************************
 *
 * multitag_t
 *
 *************************************************************************/

bool multitag_t::MT_fastnew=false;
const std::vector<uint32> multitag_t::MT_zerosizeuint32v;

StringContainer<uint8> multitag_t::MT_sc_mttagsource("multitags tagsource");
StringContainer<uint16> multitag_t::MT_sc_mtidentifier("multitags identifier");
StringContainer<uint32> multitag_t::MT_sc_mtcomment("multitags comment");


std::vector<int8> multitag_t::MT_cache_identifier_isgff3so_entry;


#ifdef MIRA_HAS_EDIT
#include <EdIt/hypothesen.H>
multitag_t fault_region::fr_tagED_C(Read::REA_defaulttag_ED_C);
multitag_t fault_region::fr_tagED_D(Read::REA_defaulttag_ED_D);
multitag_t fault_region::fr_tagED_I(Read::REA_defaulttag_ED_I);
#endif



// multitag_t::newIdentifier() and newComment() below will need initialised
//  static entries in the multitag class. Therefore, make sure we initialised
//   that before using newIdentifier() or newComment()

const multitag_t::mte_src_t multitag_t::MT_tagsrcentry_idEmpty=multitag_t::newSource("");
const multitag_t::mte_src_t multitag_t::MT_tagsrcentry_idMIRA=multitag_t::newSource("MIRA");
const multitag_t::mte_src_t multitag_t::MT_tagsrcentry_idGenBank=multitag_t::newSource("GenBank");
const multitag_t::mte_src_t multitag_t::MT_tagsrcentry_idGFF3=multitag_t::newSource("GFF3");


/*************************************************************************
 *
 * Contig
 *
 *************************************************************************/

std::unordered_set<std::string> Contig::CON_cebugnames;
std::vector<multitag_t::mte_id_t> Contig::CON_danger_zones_ids;
std::vector<multitag_t::mte_id_t> Contig::CON_baselock_ids;
std::vector<multitag_t::mte_id_t> Contig::CON_snplock_ids;

const multitag_t::mte_id_t Contig::CON_tagentry_idEmpty=multitag_t::newIdentifier("");
const multitag_t::mte_id_t Contig::CON_tagentry_idALUS=multitag_t::newIdentifier("ALUS");
const multitag_t::mte_id_t Contig::CON_tagentry_idREPT=multitag_t::newIdentifier("REPT");
const multitag_t::mte_id_t Contig::CON_tagentry_idSRMc=multitag_t::newIdentifier("SRMc");
const multitag_t::mte_id_t Contig::CON_tagentry_idWRMc=multitag_t::newIdentifier("WRMc");
const multitag_t::mte_id_t Contig::CON_tagentry_idSAOc=multitag_t::newIdentifier("SAOc");
const multitag_t::mte_id_t Contig::CON_tagentry_idSROc=multitag_t::newIdentifier("SROc");
const multitag_t::mte_id_t Contig::CON_tagentry_idSIOc=multitag_t::newIdentifier("SIOc");
const multitag_t::mte_id_t Contig::CON_tagentry_idSGPc=multitag_t::newIdentifier("SGPc");
const multitag_t::mte_id_t Contig::CON_tagentry_idPSHP=multitag_t::newIdentifier("PSHP");
const multitag_t::mte_id_t Contig::CON_tagentry_idED_D=multitag_t::newIdentifier("ED_D");
const multitag_t::mte_id_t Contig::CON_tagentry_idED_C=multitag_t::newIdentifier("ED_C");
const multitag_t::mte_id_t Contig::CON_tagentry_idED_I=multitag_t::newIdentifier("ED_I");

const multitag_t::mte_id_t Contig::CON_tagentry_idESDN=multitag_t::newIdentifier("ESDN");

const multitag_t::mte_id_t Contig::CON_tagentry_idSTMS=multitag_t::newIdentifier("STMS");
const multitag_t::mte_id_t Contig::CON_tagentry_idSTMU=multitag_t::newIdentifier("STMU");
const multitag_t::mte_id_t Contig::CON_tagentry_idUNSc=multitag_t::newIdentifier("UNSc");   // UNSure, contig

const multitag_t::mte_id_t Contig::CON_tagentry_idIUPc=multitag_t::newIdentifier("IUPc");   // IUPAC in consensus

const multitag_t::mte_id_t Contig::CON_tagentry_idMCVc=multitag_t::newIdentifier("MCVc");   // missing coverage in consensus
const multitag_t::mte_id_t Contig::CON_tagentry_idDGPc=multitag_t::newIdentifier("DGPc");   // Dubious Gap Position

const multitag_t::mte_id_t Contig::CON_tagentry_idSOFApolyA_signal_sequence=multitag_t::newIdentifier("polyA_signal_sequence");

const multitag_t::mte_co_t Contig::CON_tagentry_coEmpty=multitag_t::newComment("");


/*************************************************************************
 *
 * Read
 *
 *************************************************************************/

StringContainer<uint32> Read::REA_sc_readname("Read:: read name");
StringContainer<uint8> Read::REA_sc_processstatus("Read:: process status");
StringContainer<uint32> Read::REA_sc_asped("asped");

const multitag_t::mte_id_t Read::REA_tagentry_idEmpty=multitag_t::newIdentifier("");

const multitag_t::mte_id_t Read::REA_tagentry_idMINF=multitag_t::newIdentifier("MINF");
const multitag_t::mte_id_t Read::REA_tagentry_idMIT2=multitag_t::newIdentifier("MIT2");
const multitag_t::mte_id_t Read::REA_tagentry_idCOMM=multitag_t::newIdentifier("COMM");

const multitag_t::mte_id_t Read::REA_tagentry_idSRMr=multitag_t::newIdentifier("SRMr");
const multitag_t::mte_id_t Read::REA_tagentry_idCRMr=multitag_t::newIdentifier("CRMr");
const multitag_t::mte_id_t Read::REA_tagentry_idWRMr=multitag_t::newIdentifier("WRMr");
const multitag_t::mte_id_t Read::REA_tagentry_idSAOr=multitag_t::newIdentifier("SAOr"); //  SNP intrA Organism in Read
const multitag_t::mte_id_t Read::REA_tagentry_idSROr=multitag_t::newIdentifier("SROr"); //  SNP inteR Organism in Read
const multitag_t::mte_id_t Read::REA_tagentry_idSIOr=multitag_t::newIdentifier("SIOr"); //  SNP Intra- and inter Organism in Read
const multitag_t::mte_id_t Read::REA_tagentry_idSAOm=multitag_t::newIdentifier("SAOm"); // 'm': Marker versions for GFF3
const multitag_t::mte_id_t Read::REA_tagentry_idSROm=multitag_t::newIdentifier("SROm");
const multitag_t::mte_id_t Read::REA_tagentry_idSIOm=multitag_t::newIdentifier("SIOm");
const multitag_t::mte_id_t Read::REA_tagentry_idMCVm=multitag_t::newIdentifier("MCVm");
const multitag_t::mte_id_t Read::REA_tagentry_idSRMm=multitag_t::newIdentifier("SRMm");
const multitag_t::mte_id_t Read::REA_tagentry_idWRMm=multitag_t::newIdentifier("WRMm");

const multitag_t::mte_id_t Read::REA_tagentry_idESDN=multitag_t::newIdentifier("ESDN");

const multitag_t::mte_id_t Read::REA_tagentry_idUNSr=multitag_t::newIdentifier("UNSr");
const multitag_t::mte_id_t Read::REA_tagentry_idMNRr=multitag_t::newIdentifier("MNRr");      // Masked Nasty Repeat

const multitag_t::mte_id_t Read::REA_tagentry_idHAF0=multitag_t::newIdentifier("HAF0");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF1=multitag_t::newIdentifier("HAF1");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF2=multitag_t::newIdentifier("HAF2");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF3=multitag_t::newIdentifier("HAF3");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF4=multitag_t::newIdentifier("HAF4");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF5=multitag_t::newIdentifier("HAF5");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF6=multitag_t::newIdentifier("HAF6");
const multitag_t::mte_id_t Read::REA_tagentry_idHAF7=multitag_t::newIdentifier("HAF7");

const multitag_t::mte_id_t Read::REA_tagentry_idRLE1=multitag_t::newIdentifier("RLE1");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE2=multitag_t::newIdentifier("RLE2");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE3=multitag_t::newIdentifier("RLE3");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE4=multitag_t::newIdentifier("RLE4");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE5=multitag_t::newIdentifier("RLE5");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE6=multitag_t::newIdentifier("RLE6");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE7=multitag_t::newIdentifier("RLE7");
const multitag_t::mte_id_t Read::REA_tagentry_idRLE8=multitag_t::newIdentifier("RLE8");

const multitag_t::mte_id_t Read::REA_tagentry_idKMRF=multitag_t::newIdentifier("KMRF");


const multitag_t::mte_id_t Read::REA_tagentry_idDGNr=multitag_t::newIdentifier("DGNr");

const multitag_t::mte_id_t Read::REA_tagentry_idMFSM=multitag_t::newIdentifier("MFSM");  // MIRA Force Short-Read Merge


const multitag_t::mte_id_t Read::REA_tagentry_idALUS=multitag_t::newIdentifier("ALUS");
const multitag_t::mte_id_t Read::REA_tagentry_idREPT=multitag_t::newIdentifier("REPT");
const multitag_t::mte_id_t Read::REA_tagentry_idSVEC=multitag_t::newIdentifier("SVEC");

const multitag_t::mte_id_t Read::REA_tagentry_idSOFAdatabank_entry	      =multitag_t::newIdentifier("databank_entry");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAcontig		      =multitag_t::newIdentifier("contig");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAgene		      =multitag_t::newIdentifier("gene");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFACDS		              =multitag_t::newIdentifier("CDS");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAexon		      =multitag_t::newIdentifier("exon");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAintron		      =multitag_t::newIdentifier("intron");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFApolyA_sequence            =multitag_t::newIdentifier("polyA_sequence");

const multitag_t::mte_id_t Read::REA_tagentry_idSOFAmRNA		      =multitag_t::newIdentifier("mRNA");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAtranscript	              =multitag_t::newIdentifier("transcript");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAprimary_transcript        =multitag_t::newIdentifier("primary_transcript");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFArRNA		      =multitag_t::newIdentifier("rRNA");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAscRNA		      =multitag_t::newIdentifier("scRNA");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAsnRNA		      =multitag_t::newIdentifier("snRNA");
const multitag_t::mte_id_t Read::REA_tagentry_idSOFAtRNA                      =multitag_t::newIdentifier("tRNA");

const multitag_t::mte_id_t Read::REA_tagentry_idSOFAmiraitagstore=multitag_t::newIdentifier("experimental_feature");

const multitag_t Read::REA_defaulttag_SRMr("SRMr","","MIRA");
const multitag_t Read::REA_defaulttag_CRMr("CRMr","","MIRA");
const multitag_t Read::REA_defaulttag_WRMr("WRMr","","MIRA");
const multitag_t Read::REA_defaulttag_SAOr("SAOr","","MIRA");
const multitag_t Read::REA_defaulttag_SROr("SROr","","MIRA");
const multitag_t Read::REA_defaulttag_SIOr("SIOr","","MIRA");

const multitag_t Read::REA_defaulttag_ED_C("ED_C","","EdIt");
const multitag_t Read::REA_defaulttag_ED_D("ED_D","","EdIt");
const multitag_t Read::REA_defaulttag_ED_I("ED_I","","EdIt");

const multitag_t Read::REA_defaulttag_UNSr("UNSr","","MIRA");
const multitag_t Read::REA_defaulttag_MNRr("MNRr","","MIRA");

const multitag_t Read::REA_defaulttag_PSHP("PSHP","","MIRA");
const multitag_t Read::REA_defaulttag_CJSP("CJSP","","MIRA");

const multitag_t Read::REA_defaulttag_SOFApolyA_sequence("polyA_sequence","","MIRA");


const multitag_t::mte_co_t Read::REA_tagentry_coEmpty=multitag_t::newComment("");
const multitag_t::mte_co_t Read::REA_tagentry_coUnknown=multitag_t::newComment("UNKNOWN??? Please contact author.");
const multitag_t::mte_co_t Read::REA_tagentry_coSRMr=multitag_t::newComment("Strong Repeat Marker base");
const multitag_t::mte_co_t Read::REA_tagentry_coCRMr=multitag_t::newComment("Carbon-copy Repeat Marker base");
const multitag_t::mte_co_t Read::REA_tagentry_coWRMr=multitag_t::newComment("Weak Repeat Marker base");
const multitag_t::mte_co_t Read::REA_tagentry_coSAOr=multitag_t::newComment("SNP intrA Organism");
const multitag_t::mte_co_t Read::REA_tagentry_coSROr=multitag_t::newComment("SNP inteR Organism");
const multitag_t::mte_co_t Read::REA_tagentry_coSIOr=multitag_t::newComment("SNP Intra- and inter Organism");
const multitag_t::mte_co_t Read::REA_tagentry_coPSHP=multitag_t::newComment("Pyrosequencing Suspicious HomoPolymer");
const multitag_t::mte_co_t Read::REA_tagentry_coUNSr=multitag_t::newComment("Unsure, read");

const std::vector<multitag_t::mte_id_t> Read::REA_allhaftags={Read::REA_tagentry_idHAF0,
							      Read::REA_tagentry_idHAF1,
							      Read::REA_tagentry_idHAF2,
							      Read::REA_tagentry_idHAF3,
							      Read::REA_tagentry_idHAF4,
							      Read::REA_tagentry_idHAF5,
							      Read::REA_tagentry_idHAF6,
							      Read::REA_tagentry_idHAF7};

/*************************************************************************
 *
 * ReadGroupLib
 *
 *************************************************************************/

std::vector<ReadGroupLib::rginfo_t> ReadGroupLib::RG_static_infolib;

bool ReadGroupLib::RG_strainids_clean;
int8 ReadGroupLib::RG_numstrains;

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

const std::vector<std::string> ReadGroupLib::RG_namesofseqtypes = {
  "Sanger",
  "454",
  "IonTor",
  "PcBioHQ",
  "PcBioLQ",
  "Text",
  "Solexa",
  "Solid",
  "UNDEFINED!!! SHOULD NEVER BE SEEN"
};
const std::vector<std::string> ReadGroupLib::RG_lcnamesofseqtypes;
const std::vector<std::string> ReadGroupLib::RG_shortnamesofseqtypes = {
  "san",
  "454",
  "ion",
  "pbh",
  "pbl",
  "txt",
  "sxa",
  "sid",
  "NDF!!! SHOULD NEVER BE SEEN!!!"
};
const std::vector<std::string> ReadGroupLib::RG_samnamesofseqtypes = {
  "CAPILLARY",
  "LS454",
  "IONTORRENT",
  "PACBIO",
  "PACBIO",
  "TEXT",
  "ILLUMINA",
  "SOLID",
  "UNDEFINED!!! SHOULD NEVER BE SEEN"
};

const std::vector<std::string> ReadGroupLib::RG_namesofnamingschemes = {
  "unknown",
  "Sanger",
  "TIGR",
  "fr",
  "Solexa",
  "StLouis",
  "SRA",
  "none",
  "UNDEFINED!!! SHOULD NEVER BE SEEN"
};

const std::vector<std::string> ReadGroupLib::RG_namesofsegmentplacements = {
  "outies",
  "innies",
  "unknown",
  "lefties",
  "righties",
  "samedir",
  "UNDEFINED!!! SHOULD NEVER BE SEEN"
};


const std::string RG_emptystring;

// keep this last!
const bool ReadGroupLib::RG_initialisedstatics=ReadGroupLib::staticInitialiser();



/*************************************************************************
 *
 * gff_parse
 *
 *************************************************************************/

const std::string GFFParse::GFFP_emptystring;

const std::vector<std::string> GFFParse::GFFP_gff3scankeys= {
  "ID=",       // having this first second speeds up for standard GFF3 tags
  "gff3str=",  // and this second speeds up scanning for MIRA tags (MAF V1 files)
  "Name=",
  "Alias=",
  "Parent=",
  "Note=",
  "Target=",
  "Dbxref=",
  "miragff3=",
  "miraitag=",
  "gff3sco=",
  "gff3pha=",
  "gff3src="
};

std::unordered_map<std::string,std::regex> GFFParse::GFFP_regex_extractkeys;
