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


#include "util/fileanddisk.H"

#include "mira/assembly.H"
#include "mira/maf_parse.H"
#include "mira/readpool_io.H"


#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <regex>

using std::cout;
using std::cerr;
using std::endl;



#define CEBUG(bla)
#define CEBUGF(bla)





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpContigs()
{
  cout << "The assembled project has " << AS_contigs.size() << " objects.\n";

  Contig::setCoutType(Contig::AS_TEXT);
  for(const auto & ce : AS_contigs){
    cout << ce << "\n";
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Assembly::loadSequenceData()
{
  FUNCSTART("void Assembly::loadSequenceData()");

#ifdef MIRAMEMORC
  MemORC::statistics();
#endif

  discard();

#ifdef MIRAMEMORC
  MemORC::statistics();
#endif

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "\n";

  if(AS_resumeasembly){
    if(AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping){
      MIRANOTIFY(Notify::FATAL,"Resume is currently possible only for de-novo assemblies, sorry.");
    }
    loadSequenceData_resume();
  }else{
    //loadSequenceData_new();
    loadSequenceDataFromManifest();
    if(AS_miraparams[0].getAssemblyParams().as_rle_reads){
      cout << "RLE reads ... "; cout.flush();
      for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
	AS_readpool[ri].rleRead();
	//cout << AS_readpool[ri] << endl;
      }
      cout << "done." << endl;
      cout << "Tagging reads by RLE ... "; cout.flush();
      //colourReadsByRLE();
      cout << "done." << endl;
    }
  }

  AS_warnings.dumpWarnings();

#ifdef MIRAMEMORC
  MemORC::statistics();
#endif
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadSequenceData_resume()
{
  FUNCSTART("void Assembly::loadSequenceData_resume()");

//  cout << "Assembly::loadSequenceData_resume() must be adapted to readgroups\n"; exit(0);

  ReadGroupLib::discard();

  MAFParse mafp(&AS_readpool, nullptr, &AS_miraparams);
  std::vector<uint32> dummy;
  mafp.load(buildDefaultCheckpointFileName(AS_miraparams[0].getFileParams().chkpt_readpool),
	    ReadGroupLib::SEQTYPE_SANGER,
	    1,
	    dummy,
	    false,
	    nullptr
	  );

  bool templatesusable=AS_readpool.makeTemplateIDs(AS_miraparams[0].getNagAndWarnParams().nw_check_templateproblems);
  if(!templatesusable) {
    cout << "No useful template information found.\n";
  }

//  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
//
//  std::vector<uint32> lrperseqtype(ReadGroupLib::SEQTYPE_END,0);
//
//  size_t seqsloaded=loadMAF(buildDefaultCheckpointFileName(as_fixparams.as_infile_chkptMAF),
//			    ReadGroupLib::SEQTYPE_SANGER,
//			    1,                    // load directly
//			    lrperseqtype);
//

  dumpSomeStatistics();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::loadSequenceDataFromManifest()
{
  FUNCSTART("void Assembly::loadSequenceDataFromManifest()");

  std::vector<Manifest::manifestloadentry_t> & manifestdata2load=AS_manifest.MAN_manifestdata2load;

  // quick consistency check for mapping
  {
    bool seenbb=false;
    bool seennotbb=false;
    for(const auto & mle : manifestdata2load){
      if(mle.loadasbackbone) {
	seenbb=true;
      }else{
	seennotbb=true;
      }
    }
    if(AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping){
      if(!seenbb){
	MIRANOTIFY(Notify::FATAL,"The \"job=...\" definition of the manifest says you want a mapping assembly, but no backbone sequence was given in the readgroups.\n\nDid you forget an 'is_reference' in one of the readgroups?");
      }
    }else{
      if(seenbb){
	MIRANOTIFY(Notify::FATAL,"The \"job=...\" definition of the manifest says you want a de-novo assembly, but there is a backbone/reference sequence given in the readgroups.\n\nThis is a slight contradiction, make up your mind please.");
      }
    }
  }

  // denovo non-genome (EST/RNASeq, fragments etc) using Illumina but without -CL:kjd:kjck is not a good idea
  if(!AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping
     && !AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
     && (!AS_miraparams[0].getAssemblyParams().as_clip_kmer_junkdetection
	 || !AS_miraparams[0].getAssemblyParams().as_clip_kmer_junkkill)) {
    bool needwarn=false;
    for(const auto & mle : manifestdata2load){
      if(mle.rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA) {
	needwarn=true;
      }
    }

    if(needwarn){
      std::string wstr("You started a de-novo assembly in non-genome mode (est, fragments or clustering) using Illumina technology reads with the following parameters: ");
      if(!AS_miraparams[0].getAssemblyParams().as_clip_kmer_junkdetection) wstr+="-CL:kjd=no";
      if(!AS_miraparams[0].getAssemblyParams().as_clip_kmer_junkkill){
	if(!wstr.empty()) wstr+=' ';
	wstr+="-CL:kjck=no";
      }
      wstr+="\nThis is a pretty bad idea. I strongly suggest you turn the parameter(s) above to 'on'. Really, I mean it.";
      if(AS_miraparams[0].getNagAndWarnParams().nw_check_illuminajunkinest==NWSTOP){
	wstr+=" Alternatively: switch off this warning via -NW:cijie=no or warn only with -NW:cijie=warn'";
	MIRANOTIFY(Notify::FOOL,wstr);
      }else{
	wstr+=" Well, you chose to forego this warning with -NW:cijie=WARN";
      }
      if(AS_miraparams[0].getNagAndWarnParams().nw_check_illuminajunkinest!=NWNONE){
	std::string tmpstr("PARAMETER_WARN_ILLUMINA_JUNK");
	AS_warnings.setWarning(tmpstr.c_str(),1,"Parameters switched off junk detection in Illumina non-genome assembly",wstr);
      }
    }
  }

  if(AS_miraparams[0].getNagAndWarnParams().nw_check_proposedendclipping != NWNONE) {
    // for any sequencing technolgy except TEXT, check for proposed end clipping.
    bool needwarn=false;
    std::vector<uint8> seqseen(ReadGroupLib::SEQTYPE_END,0);
    for(const auto & mle : manifestdata2load){
      if(!mle.loadasbackbone && mle.rgid.getSequencingType()!=ReadGroupLib::SEQTYPE_TEXT) {
	seqseen[mle.rgid.getSequencingType()]=1;
	needwarn=true;
      }
    }
    if(needwarn){
      std::string tstr;
      std::string pstr;
      for(uint32 st=1; st<seqseen.size(); ++st){
	if(seqseen[st] && AS_miraparams[st].getAssemblyParams().as_clip_proposeendclips==false){
	  if(!tstr.empty()) {
	    tstr+=',';
	  }
	  tstr+=ReadGroupLib::getNameOfSequencingType(st);
	  if(!AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
	     &&  AS_miraparams[st].getAssemblyParams().as_clip_pec_continuous==false
	     && (st==ReadGroupLib::SEQTYPE_SOLEXA
		 || st==ReadGroupLib::SEQTYPE_IONTORRENT)){
	  }
	}
      }
      if(!tstr.empty()){
	std::string wstr("The proposed end clipping (-CL:pec) has been switched off.");
	wstr+="\nThis is a pretty bad idea. I strongly suggest you turn it back on with the following line in the manifest file:\n  parameters= COMMON_SETTINGS -CL:pec=yes\nReally, I mean it.";
	if(AS_miraparams[0].getNagAndWarnParams().nw_check_proposedendclipping==NWSTOP){
	  wstr+=" Alternatively: switch off this warning via -NW:cpec=no or warn only with -NW:cpec=warn'";
	  MIRANOTIFY(Notify::FOOL,wstr);
	}else{
	  wstr+=" Well, you chose to forego this warning with -NW:cpec=warn";
	  std::string tmpstr("PARAMETER_WARN_PROPOSED_END_CLIPPING");
	  AS_warnings.setWarning(tmpstr.c_str(),1,"Parameters switched off proposed end clipping",wstr);
	}
      }
    }
  }


  // for Illumina & Ion, check for continuous proposed end clipping.
  // in de-novo EST/fragment/clustering mode
  if(AS_miraparams[0].getNagAndWarnParams().nw_check_proposedendclipping != NWNONE
     && !AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
     && !AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping){
    bool needwarn=false;
    std::vector<uint8> seqseen(ReadGroupLib::SEQTYPE_END,0);
    for(const auto & mle : manifestdata2load){
      if(!mle.loadasbackbone
	 && (mle.rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA
	     || mle.rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT)) {
	seqseen[mle.rgid.getSequencingType()]=1;
	needwarn=true;
      }
    }
    if(needwarn){
      std::string tstr;
      std::string pstr;
      std::string tech;
      for(uint32 st=1; st<seqseen.size(); ++st){
	if(seqseen[st] && AS_miraparams[st].getAssemblyParams().as_clip_pec_continuous==false){
	  if(!tstr.empty()) {
	    tstr+=',';
	    pstr+=' ';
	  }else{
	    pstr+="parameters= ";
	  }
	  tech=ReadGroupLib::getNameOfSequencingType(st);
	  tstr+=tech;
	  boost::to_upper(tech);
	  pstr+=tech+"_SETTINGS -CL:pecc=yes";
	}
      }
      if(!tstr.empty()){
	std::string wstr("The continuous proposed end clipping has been switched off for the following technology: "+tstr);
	wstr+="\nIn de-novo, non-genome assembly modes, this is a pretty bad idea. I strongly suggest you turn it back on with the following line in the manifest file:\n  "+pstr+"\nReally, I mean it.";
	if(AS_miraparams[0].getNagAndWarnParams().nw_check_proposedendclipping==NWSTOP){
	  wstr+=" Alternatively: switch off this warning via -NW:cpec=no or warn only with -NW:cpec=warn'";
	  MIRANOTIFY(Notify::FOOL,wstr);
	}else{
	  wstr+=" Well, you chose to forego this warning with -NW:cpec=warn";
	  std::string tmpstr("PARAMETER_WARN_PROPOSED_END_CLIPPING");
	  AS_warnings.setWarning(tmpstr.c_str(),1,"Parameters switched off proposed end clipping",wstr);
	}
      }
    }
  }


  // for delayed transformation reads to contigs
  // needed to simplify life when loading .fna and .gff3 with annotations
  // this gives the Readpool loader the opportunity to load the sequence,
  //  then the annotation, map the annotation to the sequence
  // then only make the contigs
  //
  // if this were not done this way, one would have to annotate both the readpool read
  //   AND search contigs for this read and annotate that one too
  // TODO: might be an idea for "re-annotation"
  std::vector<readid_t> readsasbackbonecontigs;

  // this loads the data completely (no contig or read callbacks given)
  streamSequenceDataFromManifest(AS_miraparams,
				 AS_manifest,
				 AS_readpool,
				 &AS_bbcontigs,
				 &readsasbackbonecontigs);

  checkForReadNameLength(AS_miraparams[0].getNagAndWarnParams().nw_check_mrnlvalue,
			 AS_miraparams[0].getNagAndWarnParams().nw_check_maxreadnamelength==NWSTOP);

  bool hassomeerror=false;
  {
    std::string logname(buildFileName(0,"","",
				 AS_miraparams[0].getAssemblyParams().as_tmpf_clippings,
				 ".txt","",false));
    std::string logprefix("load ancillary: ");

  for(const auto & mle : manifestdata2load){
    for(auto & fnfte : mle.ancillaryfilesfoundfordata){
      cout << "Merging ancillary data " << fnfte.fn << " type " << fnfte.ft << endl;
	if(fnfte.ft=="xml"){
	  AS_readpool.mergeXMLTraceInfo(fnfte.fn,AS_miraparams[0].getNagAndWarnParams().nw_check_templateproblems);
	}else if(fnfte.ft=="ssaha2"){
	  AS_readpool.mergeSSAHA2SMALTVecScreenData(fnfte.fn,
						    false,
						    AS_miraparams,
						    logname,
						    logprefix);
	}else if(fnfte.ft=="smalt"){
	  AS_readpool.mergeSSAHA2SMALTVecScreenData(fnfte.fn,
						    true,
						    AS_miraparams,
						    logname,
						    logprefix);
	}else{
	  BUGIFTHROW(true, "Unknown ancillary file type '" << fnfte.fn << "' for data '" << fnfte.ft << "'\n");
	}
      }
    }
  }

  if(hassomeerror){
    MIRANOTIFY(Notify::FATAL,"While looking at or loading ancillary files named in the manifest, some errors occurred. Please check the log.");
  }

  // delayed creation of backbone contigs from single reads
  if(!readsasbackbonecontigs.empty()){
    Contig con(&AS_miraparams, AS_readpool);
    for(auto & rpi : readsasbackbonecontigs){
      AS_bbcontigs.push_back(con);
      AS_bbcontigs.back().addFirstRead(rpi,1);
    }
  }

  if(!AS_bbcontigs.empty()){
    postLoadBackbone();
  }


  // small sanity check (because things may go wrong somewhere without me knowing ...)
  for(uint32 rgi=1; rgi< ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    BUGIFTHROW(rgid.getSequencingType() >= ReadGroupLib::SEQTYPE_END, "readgroup with impossible sequencing type found: st="<<static_cast<uint16>(rgid.getSequencingType()) << " for readgroup\n" << rgid);
  }

  // special handling for 454 data (20.09.2008)
  //  as there's no separate seqvec clip in the XML files (until now)
  //  and MIRA now uses a "clip back, extend later" strategy for 454
  //  reads, the 454 adaptor must be protected from extension as it
  //  happens often enough that two reads start the sequencing adaptor
  //  at the same time ... WHICH WOULD THEN BE UNCOVERED!
  //
  // E.g.
  //  ACCGTCAGTCAGTCAGTGTTGACGTGTCAccctgagacacgcaacaggggatagacaaggca
  //  ACCGTACGTCAG*CAGTGTTGACGTGTCAccctgagacacgcaacaggggatagacaaggca
  //
  // Two possible worarounds
  // 1) instruct extendADS() not to extend into lower case (bad, relies
  //    on case information)
  // 2) transform the right qual clip into a vec clip if there is no
  //    vec clip
  //
  // we'll do number 2 here

  // find out what we have in the pool
  AS_seqtypespresent.clear();
  AS_seqtypespresent.resize(ReadGroupLib::SEQTYPE_END,false);
  for(uint32 rgi=1; rgi< ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(!rgid.isRail() && !rgid.isBackbone()){
      AS_seqtypespresent[rgid.getSequencingType()]=true;
    }
  }

#if CPP_READ_SEQTYPE_END != 8
#error "Check if new seqtype needs same workaround."
#endif
  if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_454GS20]){
    uint32 changecount=0;
    for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
      if(!AS_readpool.getRead(rnr).isBackbone() && !AS_readpool.getRead(rnr).isRail()){
	if(AS_readpool.getRead(rnr).hasValidData()){
	  // if no right seq vec but a right clip
	  if(AS_readpool.getRead(rnr).getRSClipoff() == static_cast<int32>(AS_readpool.getRead(rnr).getLenSeq())
	     && AS_readpool.getRead(rnr).getRQClipoff() != static_cast<int32>(AS_readpool.getRead(rnr).getLenSeq())){
	    // make right seq vec = right clip
	    AS_readpool.getRead(rnr).setRSClipoff(AS_readpool.getRead(rnr).getRQClipoff());
	    changecount++;
	  }
	}
      }
    }
    if(changecount){
      cout << "Note: " << changecount << " reads with 454 data had quality clips given, but no sequencing vector clip.\n"
	   << "For MIRA to run properly with read extension, those quality clips have been\n"
	   << "changed to sequencing vector clips.\n\n";
    }
  }

  postLoad();
  cout << "Have read pool with " << AS_readpool.size() << " reads.\n";

  dumpSomeStatistics();

  basicDataChecks();
  basicReadGroupChecks();

  clipsAfterLoad();

  dumpSomeStatistics();

  std::vector<uint32> dummy;
  if(AS_bbcontigs.empty()) {
    sortReadPool(dummy);
  }else{
    addRailsToBackbones();
    sortReadPool(dummy);

    std::vector<uint32> reversemap;
    for(auto & bbce : AS_bbcontigs){
      bbce.exchangeReadIDs(dummy,reversemap);
      // a few things to be redone
      //bbcI->setupAsBackBoneContig();
    }
    dumpSomeStatistics();
  }

  std::ofstream fout((AS_miraparams[0].getDirectoryParams().dir_tmp+'/'+AS_miraparams[0].getAssemblyParams().as_tmpf_poolinfo+".lst"), std::ios::out| std::ios::trunc);
  AS_readpool.dumpPoolInfo(fout);

  FUNCEND();
  return;
}


void Assembly::sortReadPool(std::vector<uint32> & dummy)
{
  FUNCSTART("void Assembly::sortReadPool()");

  AS_readpool.checkTemplateIDs("Before sorting: template partners have template mismatch!");

  dummy.clear();
  cout << "Sorting reads ... "; cout.flush();
  AS_readpool.sortPoolToMIRAStandard(dummy);
  cout << "done."<<endl;

  if(AS_debrisreason.size()){
    BUGIFTHROW(AS_debrisreason.size() != dummy.size(),"AS_debrisreason.size() != dummy.size()");
    std::vector<uint8> tmp;
    tmp.resize(AS_debrisreason.size());
    auto dstI=tmp.begin();
    auto srcI=dummy.begin();
    for(; dstI!=tmp.end(); ++srcI, ++dstI){
      *dstI=AS_debrisreason[*srcI];
    }
    AS_debrisreason.swap(tmp);
  }

  // The reads have Template Partner IDs internally which are now wrong
  // Get them corrected
  {
    std::vector<uint32> backmap(dummy.size(),-1);
    for(uint32 ri=0;ri<backmap.size(); ++ri){
      backmap[dummy[ri]]=ri;
    }

    for(uint32 ri=0;ri<AS_readpool.size(); ++ri){
      CEBUG(cout << AS_readpool[ri].getName() << "\t" << AS_readpool[ri].getTemplatePartnerID());
      if(AS_readpool[ri].getTemplatePartnerID()>=0){
	CEBUG("\t" << AS_readpool[AS_readpool[ri].getTemplatePartnerID()].getName());
	AS_readpool[ri].setTemplatePartnerID(backmap[AS_readpool[ri].getTemplatePartnerID()]);
	CEBUG("\t" << AS_readpool[AS_readpool[ri].getTemplatePartnerID()].getName());
      }else{
	CEBUG("\t-");
      }
      CEBUG("\n");
    }

    AS_readpool.checkTemplateIDs("After sorting: template partners have template mismatch!");
  }
}

/*************************************************************************
 *
 * All pointer parameters may be nullptr
 * - if contigsptr is nullptr, then the CAF and MAF loading routines will
 *   know they do not need to build contig structures, saves time and
 *   memory
 * - readsasbackbonecontigs will contain a list of readpool ids of those
 *   *single reads!* which were not embedded in a contig but should be
 *   seen as backbone sequence (like e.g., loaded from FASTA, GFF3, FASTQ,
 *   etc.pp
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::streamSequenceDataFromManifest(std::vector<MIRAParameters> & miraparams, Manifest & man, ReadPool & readpool, std::list<Contig> * contigsptr, std::vector<readid_t> * readsasbackbonecontigs)
{
  FUNCSTART("void Assembly::streamSequenceDataFromManifest(std::vector<MIRAParameters> & miraparams, Manifest & man, ReadPool & readpool, std::list<Contig> * contigsptr, std::vector<readid_t> * readsasbackbonecontigs)");

  ReadPoolIO rpio(readpool);
  rpio.setAttributeProgressIndicator(true);
  rpio.setAttributeFASTQQualOffset(0); // in case we load FASTQs, we want to adapt ourselves
  rpio.setAttributeFASTQTransformName(true); // get those names from Illumina transformed
  rpio.setAttributeFASTQAPreserveComment(false); // but throw away all comments

  auto & manifestdata2load=man.MAN_manifestdata2load;

  bool hassomeerror=false;
  for(const auto & mle : manifestdata2load){
    for(auto & fnfte : mle.mainfilesfoundfordata){
      // should have already been set by the manifest parser, just in case
      BUGIFTHROW(mle.loadasbackbone && !mle.rgid.isBackbone(),"manifest entry says to load as backbone, but readgroup is not backbone?");

      std::string fn2;

      if(fnfte.ft=="fasta"){
	fn2=fnfte.fn+".qual";
	rpio.setAttributeFASTAQualFileWanted(true);
      }else if(fnfte.ft=="fna"){
	rpio.setAttributeFASTAQualFileWanted(false);
      }

      size_t oldrpsize=readpool.size();
      size_t oldbbsize=0;
      if(contigsptr!=nullptr) oldbbsize=contigsptr->size();

      if(((mle.loadasbackbone && contigsptr!=nullptr) || contigsptr!=nullptr)
	 && (fnfte.ft=="maf" || fnfte.ft=="caf")){
	// Note: this branch really only if we want contigs or load as backbone (may be
	//  two different usages from the caller function)
	rpio.setAttributesForContigs(contigsptr,&miraparams);
	cout << "Loading reads or assembled contigs ";
	if(mle.loadasbackbone) cout << "as reference backbone ";
	cout << "from " << fnfte.fn << " type " << fnfte.ft << endl;
      }else{
	rpio.setAttributesForContigs(nullptr,&miraparams);
	if(mle.loadasbackbone){
	  cout << "Loading reference backbone from " << fnfte.fn << " type " << fnfte.ft << endl;
	  rpio.setAttributeFASTAQualFileWanted(false);
	}else{
	  cout << "Loading reads from " << fnfte.fn << " type " << fnfte.ft << endl;
	}
      }
      rpio.registerFile(fnfte.ft,fnfte.fn,fn2,mle.rgid,false);
      rpio.loadNextSeqs(-1,-1); // load all
      if(fnfte.ft=="fastq"){
	readpool.adaptFASTQQualValues(oldrpsize,readpool.size(),0,true);
      }

      if(contigsptr!=nullptr && contigsptr->size() != oldbbsize){
	cout << "contained " << contigsptr->size()-oldbbsize << " contigs. Only the contigs will be added as backbone.\n";
      }else if(mle.loadasbackbone){
	// there were no contigs in this file: set qualities (if needed) and mark the reads to be added as contigs later
	for(size_t rpi=oldrpsize; rpi<readpool.size(); ++rpi) {
	  if(readpool[rpi].hasValidData()){
	    if(!readpool[rpi].hasQuality()){
	      readpool[rpi].setQualities(readpool[rpi].getReadGroupID().getDefaultQual());
	      readpool[rpi].setQualityFlag(false);
	    }
	    if(readsasbackbonecontigs!=nullptr) readsasbackbonecontigs->push_back(static_cast<readid_t>(rpi));
	  }
	}
      }
    }
  }

  if(hassomeerror){
    MIRANOTIFY(Notify::FATAL,"While looking at or loading data files named in the manifest, some errors occurred. Please check the log.");
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::addRailsToBackbones()
{
  FUNCSTART("void Assembly::addRailsToBackbones()");
  std::vector<uint32> lrperseqtype(ReadGroupLib::SEQTYPE_END,0);
  uint32 longestread=0;

  for(uint32 i=0; i<AS_readpool.size(); ++i){
    if(AS_readpool[i].isRail()
       || AS_readpool[i].isBackbone()
       || AS_readpool[i].isCoverageEquivalentRead()) continue;
    if(AS_readpool[i].getLenClippedSeq()>longestread) longestread=AS_readpool[i].getLenClippedSeq();
    if(AS_readpool[i].getLenClippedSeq()>lrperseqtype[AS_readpool[i].getSequencingType()]) lrperseqtype[AS_readpool[i].getSequencingType()]=AS_readpool[i].getLenClippedSeq();
  }

  if(longestread==0){
    MIRANOTIFY(Notify::FATAL, "No read with sequence length >0 present? Did you provide data to load?");
  }

  // if wanted, determine -FM:me parameter automatically
  {
    auto & fmparams=AS_miraparams[ReadGroupLib::SEQTYPE_SOLEXA].getFinalMappingParameters();
    if(AS_miraparams[0].getFinalMappingParameters().fm_active
       && fmparams.fm_maxtotalerrors < 0
       && lrperseqtype[ReadGroupLib::SEQTYPE_SOLEXA]){
      cout << "-FM:me for Solexa is < 0, automatically determining optimal value.\n";
      AS_miraparams[ReadGroupLib::SEQTYPE_SOLEXA].getNonConstFinalMappingParameters().fm_maxtotalerrors=
	lrperseqtype[ReadGroupLib::SEQTYPE_SOLEXA]*15/100;
      cout << "set -FM:me "
	   << fmparams.fm_maxtotalerrors
	   << '\n';
    }
  }

  // if wanted, determine rail length and overlap automatically
  if(AS_miraparams[0].getAssemblyParams().as_backbone_raillength == 0){
    cout << "-SB:brl is 0, automatically determining optimal value.\n";

    // add 15% to longest read (so accomodate insertion), then times 2
    uint32 newraillength=(longestread*115/100) * 2;
    if(newraillength > 32760){
      cout << "Optimal rail would be longer than 32k, adjusting down to 32k.\n";
      newraillength=32760;
    }
    AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength=newraillength;
    cout << "brl: "
	 << AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength
	 << '\n';
  }
  if(AS_miraparams[0].getAssemblyParams().as_backbone_railoverlap == 0){
    cout << "-SB:bro is 0, automatically determining optimal value.\n";
    AS_miraparams[0].getNonConstAssemblyParams().as_backbone_railoverlap=
      AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength/2;
    cout << "bro: "
	 << AS_miraparams[0].getNonConstAssemblyParams().as_backbone_railoverlap
	 << '\n';
  }
  if(AS_miraparams[0].getAssemblyParams().as_backbone_railoverlap >=
     AS_miraparams[0].getAssemblyParams().as_backbone_raillength){
    cout << "-SB:bro is >= -SB:brl ... adjusting -SB:bro to (-SB:brl)-1\n";
    AS_miraparams[0].getNonConstAssemblyParams().as_backbone_railoverlap=
      AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength-1;
  }

  // TODO: see whether rails need strains!!!
  // Answer: yes, it helps. E.g.: markFeaturesByConsensus() and other routines do not
  //  need to work around a strain which says "I have 3 strains" (backbone, rails (==empty) and
  //  reads).
  // Set the strain name to be equal to the strain name of the first backbone read encountered

  size_t numrailscreated=0;

  ProgressIndicator<int32> P (0, AS_bbcontigs.size());
  for(auto & bbc : AS_bbcontigs){
    P.increaseprogress();
    bbc.recalcTemplateIDsAndStrainPresent();

    std::string strainname;
    bool foundbbread=false;
    for(auto & cr : bbc.getContigReads()){
      if(cr.isBackbone()){
	strainname=cr.getStrainName();
	foundbbread=true;
	break;
      }
    }
    BUGIFTHROW(foundbbread==false,"no backbone read found in backbone???");

    bool bbvalue=true;

    // add the rails
    if(bbvalue) {
      auto verbose=bbc.getVerbose();
      bbc.setVerbose(false);
      numrailscreated+=bbc.addRails(
	AS_miraparams[0].getAssemblyParams().as_backbone_raillength,
	AS_miraparams[0].getAssemblyParams().as_backbone_railoverlap,
	strainname,
	AS_miraparams[0].getAssemblyParams().as_backbone_strainname_forceforall,
	AS_miraparams[0].getAssemblyParams().as_backbone_rail_fromstrain,
	false);
      bbc.setVerbose(verbose);
    }
  }
  P.finishAtOnce();
  cout << endl;

  AS_debrisreason.resize(AS_readpool.size(),DEBRIS_NOTDEBRIS);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::clipsAfterLoad()
{
  AS_debrisreason.clear();
  AS_debrisreason.resize(AS_readpool.size(),DEBRIS_NOTDEBRIS);

  std::vector<std::unique_ptr<DataProcessing>> dpv(AS_miraparams[0].getAssemblyParams().as_numthreads);

  for(uint32 ti=0; ti<dpv.size(); ++ti){
    dpv[ti]=std::unique_ptr<DataProcessing>(new DataProcessing(&AS_miraparams));
    std::string logname(buildFileName(0,"","",
				 AS_miraparams[0].getAssemblyParams().as_tmpf_clippings
				 + "_t" + boost::lexical_cast<std::string>(ti),
				 ".txt","",false));
    cout << logname << endl;
    dpv[ti]->startLogging(logname,false);
  }

  std::string logprefix("loadclip: ");
  cout << "Post-load clips:\n";

  DataProcessing::stdTreatmentPool_MultiThread(AS_miraparams,AS_dataprocessing,dpv,AS_readpool,&AS_debrisreason,logprefix,true);
  cout << endl;

  if(AS_dataprocessing.DP_stats.cphix174){
    cout << "SEARCH MSG: PhiX 174 found: " << AS_dataprocessing.DP_stats.cphix174 << endl;
  }
  if(AS_dataprocessing.DP_stats.crrna){
    cout << "CLIP MSG: rRNA found: " << AS_dataprocessing.DP_stats.crrna << endl;
  }
  if(AS_dataprocessing.DP_stats.cadapright){
    cout << "CLIP MSG: Adaptor right found: " << AS_dataprocessing.DP_stats.cadapright << endl;
  }
  if(AS_dataprocessing.DP_stats.cadaprightpartial){
    cout << "CLIP MSG: Partial adaptor right found: " << AS_dataprocessing.DP_stats.cadaprightpartial << endl;
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpRailReads(std::ofstream & fout)
{
  Read::setCoutType(Read::AS_FASTA);
  for(uint32 i=0; i<AS_readpool.size(); ++i){
    if(AS_readpool[i].isRail()){
      fout << AS_readpool[i].isRail();
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::basicDataChecks()
{
  FUNCSTART("void Assembly::basicDataChecks()");

  uint32 nummsglr=0;
  uint32 nummsgqual=0;
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    if(AS_readpool[ri].getLenClippedSeq() > MIRASRLEN
       && !AS_readpool[ri].isBackbone()){
      if(++nummsglr<=10){
	cout << "Read " << AS_readpool[ri].getName() << " has " << AS_readpool[ri].getLenClippedSeq() << " bases (clipped). Too long (>" << MIRASRLEN << "\n";
	if(nummsglr==10) cout << "More long reads may exist, but stopping output here.\n";
      }
    }
    bool bahqual=false;
    for(auto q : AS_readpool[ri].getQualities()){
      if(q>=100){
	bahqual=true;
	break;
      }
    }
    if(bahqual){
      if(++nummsgqual<=10){
	cout << "Read " << AS_readpool[ri].getName() << " has quality values >100, this is illegal.\n";
	if(nummsgqual==10) cout << "More reads with bad quals may exist, but stopping output here.\n";
      }
    }
  }

  if(nummsglr){
    MIRANOTIFY(Notify::FATAL,"MIRA found " << nummsglr << " very long reads (too long for normal reads), see log above.");
  }
  if(nummsgqual){
    MIRANOTIFY(Notify::FATAL,"MIRA found " << nummsgqual << " reads with illegal qualities, see log above.");
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::basicReadGroupChecks()
{
  FUNCSTART("void Assembly::basicReadGroupChecks()");

  bool haserror_ep=false; // expected pair
  bool haserror_sra=false; // unexpected sra names

  for(uint32 rglid=1; rglid<ReadGroupLib::getNumReadGroups(); ++rglid){
    auto rgid = ReadGroupLib::getReadGroupID(rglid);

    uint32 numsra=0;
    uint32 numnotsra=0;
    for(uint32 rid=0; rid<AS_readpool.size(); ++rid){
      if(AS_readpool[rid].getReadGroupID()!=rgid) continue;
      auto & rname=AS_readpool[rid].getName();
      if(Read::checkStringForSRANamingScheme(rname)){
	++numsra;
	if(numsra<=10
	   && rgid.getReadNamingScheme() != ReadGroupLib::SCHEME_SRARAW
	   && AS_miraparams[0].getNagAndWarnParams().nw_check_readnamesra != NWNONE){
	  cout <<"Read " << rname << " is in SRA naming scheme, that is unexpected.\n";
	  if(numsra==10) cout << "More may be present but will not be shown anymore.\n";
	}
      }else{
	++numnotsra;
      }
    }
    if(numsra && rgid.getReadNamingScheme() != ReadGroupLib::SCHEME_SRARAW
      && AS_miraparams[0].getNagAndWarnParams().nw_check_readnamesra != NWNONE){
      std::string wstr("In readgroup ");
      wstr+=boost::lexical_cast<std::string>(rgid.getLibId())+" (named: '"+rgid.getGroupName()+"'), reads names were found which seem to follow the SRA naming scheme but the readgroup does have the naming scheme 'sra'. In case your reads were downloaded from the SRA, it is strongly suggested to tell MIRA that your read names use the corresponding naming scheme (use 'segment_naming=sra' in the manifest).";
      std::string tmpstr("READGROUP_UNEXPECTED_SRANAMES_"+boost::lexical_cast<std::string>(rgid.getLibId()));
      AS_warnings.setWarning(tmpstr.c_str(),2,"Read names in SRA naming scheme unexpected in readgroup",wstr);
      haserror_sra=true;
    }

    cout << "Checking pairs of readgroup " << rglid << " (named: '" << rgid.getGroupName() << "'): ";
    cout.flush();
    uint32 numreadswithtpids=0;
    uint32 numreads=0;
    for(uint32 rid=0; rid<AS_readpool.size(); ++rid){
      if(AS_readpool[rid].getReadGroupID()==rgid){
	++numreads;
	if(AS_readpool[rid].getTemplatePartnerID()>=0) ++numreadswithtpids;
      }
    }
    cout << " found " << numreadswithtpids<<endl;
    if(rgid.expectsReadPairs()){
      if(numreads && numreadswithtpids==0) {
	cout << "WARNING: in the above readgroup, no read is paired although the manifest says there should be pairs. This is fishy!\n";
	haserror_ep=true;
      }
    }else{
      if(numreadswithtpids>0) {
	std::string wstr("In readgroup ");
	wstr+=boost::lexical_cast<std::string>(rgid.getLibId())+" (named: '"+rgid.getGroupName()+"') paired reads were found but no pairing information given in the manifest. MIRA will estimate 'template_size' and 'segment_placement'.\nYou can suppress this warning by using the keyword 'autopairing' in the readgroup definition of the manifest file.";
	if(numsra && rgid.getReadNamingScheme() != ReadGroupLib::SCHEME_SRARAW){
	  wstr+="\nNote that this may be caused by reads following the SRA naming scheme but you having not set the naming scheme of the readgroup to be 'sra'.";
	}
	AS_warnings.setWarning("READGROUP_UNEXPECTED_PAIRS",2,"Unexpected pairs in readgroup",wstr);
      }
    }
  }

  if(haserror_ep){
    MIRANOTIFY(Notify::FATAL,"MIRA found readgroups where pairs are expected but no read has a partner. See log above and then check your input please (either manifest file or data files loaded or segment_naming scheme).");
  }

  if(haserror_sra
     && AS_miraparams[0].getNagAndWarnParams().nw_check_readnamesra == NWSTOP){
    MIRANOTIFY(Notify::FATAL,"MIRA found readgroups where read names follow the SRA naming scheme, but that was not expected, please see log above.\nRemedies: use 'segment_naming=sra' in manifest file for this readgroup. Alternatively: switch off this warning via -NW:csrn=no or warn only with -NW:csrn=warn'");
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::checkForReadNameLength(uint32 stoplength, bool stop)
{
  FUNCSTART("void Assembly::checkForReadNameLength(uint32 stoplength)");

  if(stoplength==0) return;

  // if names are too long, print out only the first 20 of each read group
  std::vector<uint32> countperrg(ReadGroupLib::getNumReadGroups(),20);

  uint32 count=0;
  bool maybemore=false;
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    if(AS_readpool[ri].getName().size()>stoplength){
      if(countperrg[AS_readpool[ri].getReadGroupID().getLibId()]){
	if(--countperrg[AS_readpool[ri].getReadGroupID().getLibId()] == 0) maybemore=true;
	if(count==0) {
	  cout << "List of read names which have problems with name length:\n";
	}
	cout << "Name too long: " << AS_readpool[ri].getName() << '\n';
      }
      ++count;
    }
  }
  if(maybemore){
    cout << count << " reads had a long name length, for brevity's sake not all were listed.\n";
  }
  if(count>0){
    std::string emsg=boost::lexical_cast<std::string>(count)+
      " reads were detected with names longer than "
      +boost::lexical_cast<std::string>(stoplength)+
      " characters (see output log for more details).\n\n"
      "While MIRA and many other programs have no problem with that, some older "
      "programs have restrictions concerning the length of the read name.\n"
      "\nExample given: the pipeline\n"
      "     CAF -> caf2gap -> gap2caf\n"
      "will stop working at the gap2caf stage if there are read names having > 40 characters "
      "where the names differ only at >40 characters.\n"
      "\nThis is a warning only, but as a couple of people were bitten by this, the default "
      "behaviour of MIRA is to stop when it sees that potential problem.\n"
      "\nYou might want to rename your reads to have <= "
      +boost::lexical_cast<std::string>(stoplength) +
      " characters. Instead of renaming reads in the input files, maybe the 'rename_prefix' functionality of manifest files is useful for you there.\n"
      "\nOn the other hand, you also can ignore this potential problem and force MIRA to "
      "continue by using the parameter: '-NW:cmrnl=warn'  or  '-NW:cmrnl=no'\n";
    if(stop){
      MIRANOTIFY(Notify::FATAL,emsg);
    }else{
      cout << "WARNING!\n" << endl;
      AS_warnings.setWarning("READ_NAME_TOO_LONG",2,"Long read names",emsg);
    }
  }
  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::postLoadBackbone()
{
  FUNCSTART("void Assembly::postLoadBackbone()");

  AS_hasbackbones=true;

  cout << "Deleting gap columns in backbones ... "; cout.flush();
  for(auto I=AS_bbcontigs.begin(); I!=AS_bbcontigs.end(); I++){
    I->deleteStarOnlyColumns(0,I->getContigLength());
  }

  cout << "Postprocessing backbone(s) ... this may take a while."<< endl;

  // mark all reads loaded in backbone as backbone
  // check that they are not named "ContigX"
  // set the strain to "backbone"
  // Backbones will not be included is Skim, makeAlignments etc.


  static const std::regex badseqnameexp("^Contig[0-9]+$");
  //return regex_match(s, e);

  // set MFSM tags
  if(0){
    std::vector<multitag_t::mte_id_t> idstoreplace;
    {
      std::string tmp="FLTR";
      idstoreplace.push_back(Read::getTagID(tmp));
      tmp="FrRN";
      idstoreplace.push_back(Read::getTagID(tmp));
    }
    std::string mfsm="MFSM";
    multitag_t::mte_id_t mtid=Read::getTagID(mfsm);
    for(auto & ce : AS_bbcontigs){
      auto & conreads=ce.getContigReads();
      for(auto & cre : conreads){
	for(uint32 it=0; it<cre.getNumOfTags(); ++it){
	  multitag_t tmp=cre.getTag(it);
	  if(find(idstoreplace.begin(),idstoreplace.end(),tmp.identifier)!=idstoreplace.end()){
	    tmp.identifier=mtid;
	    tmp.source=multitag_t::MT_tagsrcentry_idMIRA;
	    Read & ncr=const_cast<Read &>(cre);
	    ncr.addTagO(tmp);
	  }
	}
      }
    }
  }


  bool contignamesok=true;
  {
    uint32 bbnum=0;
    cout << AS_bbcontigs.size() << " to process\n";
    for(auto & ce : AS_bbcontigs){
      ++bbnum;
      // first, find a name for that contig

      auto & conreads=ce.getContigReads();

      // if it is a single read contig,
      //  set the name for that contig to be the name of the read
      if(conreads.size()==1) {
	if(conreads.begin()->getName().size()){
	  if(conreads.begin()->getName()[0]=='C'){
// BOOST: regex not compatible with _GLIBCXX_DEBUG
#ifndef _GLIBCXX_DEBUG
	    if(std::regex_match(conreads.begin()->getName(), badseqnameexp)){
	      cout << "Bad name for backbone sequence " << bbnum << ": " << conreads.begin()->getName() << '\n';
	      cout << "Backbone sequences may NOT be name 'ContigX' with 'X' being any number.\n";
	      contignamesok=false;
	    }
#endif
	  }
	}else{
	  cout << "There's a backbone sequence (number " << bbnum << ") without a name? Too bad, not allowed.\n";
	  contignamesok=false;
	}

	/* BaCh 14.09.2010: why did I add _bb? Let's remove and see whether it's better.
	   BaCh 26.10.2010: I remembered. A contig may not be named like a sequence, CAF will
	                    dump an error.
			    Solution: ...?
	*/
	ce.setContigName(conreads.begin()->getName()+"_bb");

      }

      cout << ce.getContigName() << "\t" << ce.getContigLength() << endl;

      //bool bbvalue=true;

      ////  except singlets?!
      //if(conreads.size()==1 && I->getContigLength()<4000) {
      //	bbvalue=false;
      //}

      // now let the contig do the rest of the setup
      ce.setupAsBackBoneContig();
    }
  }

  if(!contignamesok){
    MIRANOTIFY(Notify::FATAL,"Some backbones had either no names or a bad name (see log above). Stopping here, fix your sequence names.\n")
  }

  ReadGroupLib::dumpStrainIDSummary();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::postLoad()
{
  FUNCSTART("void Assembly::postLoad()");

  if(!AS_readpool.checkForDuplicateReadNames(AS_miraparams[0].getNagAndWarnParams().nw_check_duplicatereadnames)){
    MIRANOTIFY(Notify::FATAL,"MIRA found duplicate read names in your data (see log above for more info).\n\nThis should never, never be!\n\nYou may have loaded a file more than once in the manifest or\nreads may be present multiple times across your input file(s).\nEither way: fix that!\n");
  }

  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  // how many sequences in assembly


  // count how many have valid data, clean up stars in reads
  // count how many have SCF data: if none, switch off editing.

  {
    cout << "Checking reads for trace data (loading qualities if needed):\n";

    ProgressIndicator<int32> P (0, AS_readpool.size());

    std::ofstream fout((dir_params.dir_tmp+'/'+as_fixparams.as_outfile_stats_reads_invalid), std::ios::out| std::ios::trunc);

    bool can_edit=false;
    AS_num_reads_valid=0;
    for(uint32 i=0;i<AS_readpool.size();i++){
      P.progress(i);
      if(AS_readpool.getRead(i).hasValidData()){
	AS_num_reads_valid++;
	if(AS_readpool[i].isBackbone()){
	  AS_hasbackbones=true;
	}
	AS_readpool.getRead(i).removeGapsFromRead();
	if(AS_readpool.getRead(i).hasSCFData(true)){
	  can_edit=true;
	}
      } else {
	if(!AS_readpool.getRead(i).getName().empty()) {
	  fout << AS_readpool.getRead(i).getName() << "\n";
	//} else if (!AS_readpool.getRead(i).getEXPName().empty()) {
	//  fout << AS_readpool.getRead(i).getEXPName() << "\n";
	} else {
	  fout << "Unknown read (loaded as number: " << i << ")\n";
	}
      }
    }
    P.finishAtOnce();
    cout << endl;

    if(!can_edit) {
      cout << "No SCF data present in any read, EdIt automatic contig editing for Sanger data is now switched off.\n";
      const_cast<edit_parameters &>(AS_miraparams[0].getEditParams()).ed_edit_automatic_contic_editing=false;
    }
  }

  cout << AS_num_reads_valid << " reads with valid data for assembly.\n";

  if(!AS_num_reads_valid){
    MIRANOTIFY(Notify::FATAL, "No valid read in assembly?");
  }


  bool templatesusable=AS_readpool.makeTemplateIDs(AS_miraparams[0].getNagAndWarnParams().nw_check_templateproblems);
  if(!templatesusable) {
    cout << "No useful template information found.\n";
  }


  //re-adjust bbcontigs template and strain ids in contig reads
  if(AS_hasbackbones){
    for(auto & ce : AS_bbcontigs){
      ce.recalcTemplateIDsAndStrainPresent();
    }
  }

  ReadGroupLib::dumpStrainIDSummary();

  // look for quality values in reads
  {
    bool stopall=false;
    for(uint32 i=0;i<AS_readpool.size();i++){
      Read & actread=AS_readpool.getRead(i);
      if(actread.isBackbone()
	 || actread.isRail()) continue;
      if(AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_enforce_qualsinreads
	 && actread.hasValidData()
	 && actread.hasQuality()==false
	 && !actread.hasUserDefaultQuality()){
	cout << "No quality data found: (" << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType()) << ") " << actread.getName() << '\n';
	stopall=true;
      }
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    if(stopall) {
      MIRANOTIFY(Notify::FATAL,"Some reads had no quality values given (see log above),\nplease check your input data.\nIf sure that this is ok for your data, switch off this check with -AS:epoq=no\nfor any sequencing type you wish (Sanger, 454, IonTorrent, PacbioLQ, PacbioHQ, Text, Solexa, ...).\nAlso consider the '--noqualities' parameter setting.\nAlternatively, you can force to switch off this check for specific readgroups by using the 'default_qual' setting in the manifest.")
    }

  }


  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpSomeStatistics()
{
  FUNCSTART("void Assembly::dumpSomeStatistics()");
  // initialise the assembly_structure and do some statistics

  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  uint32 backbonereads=0;
  uint32 railreads=0;
  {
    std::ofstream fout((dir_params.dir_tmp+'/'+as_fixparams.as_outfile_stats_reads_tooshort), std::ios::out| std::ios::trunc);

    AS_num_reads_too_small=0;
    for(uint32 i=0;i<AS_readpool.size();i++){
      Read & actread=AS_readpool.getRead(i);

      //AS_ok_for_assembly[i]=1;
      actread.setUsedInAssembly(true);

      if(actread.hasValidData()==false){
	cout << actread.getName() << ": unable to load or other reason for invalid data.\n";
	//AS_ok_for_assembly[i]=0;
	actread.setUsedInAssembly(false);
      }else{
	if(actread.isBackbone()){
	  actread.setUsedInAssembly(false);
	  backbonereads++;
	} else if(actread.isRail()) {
	  railreads++;
	}else{
	  // throw out on minumum length if no template partner is present
	  if(actread.getLenClippedSeq() < AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength){
	    //cout << "Short length: "
	    //	 << actread.getName() << " ("
	    //	 << actread.getShortNameOfSequencingType(actread.getSequencingType())
	    //	 << "): only " << actread.getLenClippedSeq()
	    //	 << " good bases, ";
	    fout << actread.getName();
	    if(actread.getTemplatePartnerID() == -1){
	      //cout << "need: "
	      //	   << AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength
	      //	   << ". No paired end partner, rejected.\n";
	      fout << " too small and no paired end\n";
	      AS_num_reads_too_small++;
	      //AS_ok_for_assembly[i]=0;
	      actread.setUsedInAssembly(false);
	    }else{
	      if(actread.getLenClippedSeq() < 20){
		//cout << "really too small, rejected.\n";
		fout << " too small even with paired end\n";
		AS_num_reads_too_small++;
		//AS_ok_for_assembly[i]=0;
		actread.setUsedInAssembly(false);
	      }else{
		fout << " saved by paired-end\n";
		//cout << "accepted as paired-end partner is present.\n";
	      }
	    }
	  }
	}
      }
    }
  }

  if(AS_logflag_dumpusedids){
    std::ofstream fout((dir_params.dir_tmp+"/elog.usedids.lst"), std::ios::out | std::ios::trunc);
    for(uint32 i=0; i<AS_used_ids.size(); i++){
      if(AS_readpool[i].isUsedInAssembly()){
	fout << AS_readpool[i].getName()<< '\n';
      }
    }
  }

  // TODO: also take reads too short into statistics
  //


  AS_current_rls_byrg.clear();
  AS_current_rls_byrg.resize(ReadGroupLib::getNumReadGroups());
  AS_current_rls_bytype.clear();
  AS_current_rls_bytype.resize(ReadGroupLib::SEQTYPE_END);

  //remove("log.noqualities");
  std::ofstream fout((dir_params.dir_tmp+"/miralog.noqualities"), std::ios::out| std::ios::trunc);
  for(uint32 i=0;i<AS_readpool.size();i++){
    Read & actread=AS_readpool.getRead(i);
    if(actread.isBackbone() == false
       && actread.isRail() == false){

      auto rgi=actread.getReadGroupID().getLibId();

      BUGIFTHROW(rgi==0,"Read " << actread.getName() << " has rgid of 0??? Ouch.");
      BUGIFTHROW(rgi>=ReadGroupLib::getNumReadGroups(),"Read " << actread.getName() << " has rgid of " << rgi << " and num of readgroups is " << ReadGroupLib::getNumReadGroups() << "??? Ouch.");

      auto st=actread.getSequencingType();

      if(actread.hasQuality()==false){
	++AS_current_rls_bytype[st].RLS_count_noqual;
	++AS_current_rls_byrg[rgi].RLS_count_noqual;
	fout << actread.getName() << endl;
      }
      ++AS_current_rls_bytype[st].RLS_count_all;
      ++AS_current_rls_byrg[rgi].RLS_count_all;
      AS_current_rls_bytype[st].RLS_len_all+=actread.getLenSeq();
      AS_current_rls_byrg[rgi].RLS_len_all+=actread.getLenSeq();
      if(actread.isUsedInAssembly()){
	++AS_current_rls_bytype[st].RLS_count_used;
	++AS_current_rls_byrg[rgi].RLS_count_used;
	AS_current_rls_bytype[st].RLS_len_used+=actread.getLenClippedSeq();
	AS_current_rls_byrg[rgi].RLS_len_used+=actread.getLenClippedSeq();
      }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      switch(actread.getSequencingType()){
      case ReadGroupLib::SEQTYPE_TEXT :
      case ReadGroupLib::SEQTYPE_SANGER :
      case ReadGroupLib::SEQTYPE_PACBIOHQ :
      case ReadGroupLib::SEQTYPE_PACBIOLQ :
      case ReadGroupLib::SEQTYPE_IONTORRENT :
      case ReadGroupLib::SEQTYPE_454GS20 :
      case ReadGroupLib::SEQTYPE_SOLEXA : {
	if(actread.getLenClippedSeq() == actread.getLenSeq()){
	  ++AS_current_rls_bytype[st].RLS_count_withoutclips;
	  ++AS_current_rls_byrg[rgi].RLS_count_withoutclips;
	}
	break;
      }
      case ReadGroupLib::SEQTYPE_ABISOLID : {
	MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 20a.");
	break;
      }
      default : {
	cerr << "Sequencing type " << actread.getSequencingType() << " unknown?\n";
	MIRANOTIFY(Notify::FATAL, "Found unknown sequencing type in read.");
      }
      }
    }
  }
  fout.close();

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  cout << "\n===========================================================================\n";
  cout << "Backbones: " << backbonereads << "\tBackbone rails: " << railreads << "\n";

  static std::vector<uint8> displayst= {
    ReadGroupLib::SEQTYPE_SANGER,
    ReadGroupLib::SEQTYPE_454GS20,
    ReadGroupLib::SEQTYPE_IONTORRENT,
    ReadGroupLib::SEQTYPE_PACBIOHQ,
    ReadGroupLib::SEQTYPE_PACBIOLQ,
    ReadGroupLib::SEQTYPE_TEXT,
    ReadGroupLib::SEQTYPE_SOLEXA,
    ReadGroupLib::SEQTYPE_ABISOLID
  };

  cout << "Sequencing technology statistics:\n";
  cout << "\n\t";
  //cout << "\tSanger\t454\tIonTor\tPcBioHQ\tPcBioLQ\tText\tSolexa\tSOLiD\n";
  for(auto st : displayst){
    cout << "\t" << ReadGroupLib::getNameOfSequencingType(st);
  }
  cout << '\n';

  cout << "\t\t------------------------------------------------------------\n";
  cout << "Total reads";
  for(auto st : displayst){
    cout << "\t" << AS_current_rls_bytype[st].RLS_count_all;
  }
  cout << '\n';

  cout << "Reads wo qual";
  for(auto st : displayst){
    cout << "\t" << AS_current_rls_bytype[st].RLS_count_noqual;
  }
  cout << '\n';

  cout << "Used reads";
  for(auto st : displayst){
    cout << "\t" << AS_current_rls_bytype[st].RLS_count_used;
  }
  cout << '\n';

  cout << "Avg. tot rlen";
  for(auto st : displayst){
    cout << "\t" << AS_current_rls_bytype[st].getAvgLenAll();
  }
  cout << '\n';

  cout << "Avg. used rlen";
  for(auto st : displayst){
    cout << "\t" << AS_current_rls_bytype[st].getAvgLenUsed();
  }
  cout << '\n';

  cout << "W/o clips";
  for(auto st : displayst){
    cout << "\t" << AS_current_rls_bytype[st].RLS_count_withoutclips;
  }
  cout << '\n';

  cout << "\n\nReadgroup statistics:\n";
  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    if(AS_current_rls_byrg[rgi].RLS_count_all){
      cout << "RG " << rgi << "\t" << ReadGroupLib::getReadGroupID(rgi).getNameOfSequencingType()
	   << "\tavg total len: " << AS_current_rls_byrg[rgi].getAvgLenAll()
	   << "\tavg clip len: " << AS_current_rls_byrg[rgi].getAvgLenUsed()
	   << "\ttotal bases: " << AS_current_rls_byrg[rgi].RLS_len_all
	   << "\tused bases: " << AS_current_rls_byrg[rgi].RLS_len_used
	   << '\n';
    }
  }

  cout << "===========================================================================\n\n";


  if(AS_readpool.size()-AS_num_reads_too_small-backbonereads-railreads<=0) {
    //MIRANOTIFY(Notify::FATAL, "No read can be used for assembly.");
  }

  cout << endl;

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveExtTmpContig(Contig & con, std::string basename)
{
  FUNCSTART("void Assembly::saveExtTmpContig(std::string prepost)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  //directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();

  if(con.getNumReadsInContig() > 1
     ||   as_fixparams.as_output_exttmp_alsosinglets){

    if (as_fixparams.as_output_exttmp_caf) {
      std::string filename=basename+".caf";

      cout << "Logging this contig to file: " << filename << endl;

      std::ofstream cafout(filename, std::ios::out | std::ios::trunc);
      Contig::setCoutType(Contig::AS_CAF);
      cafout << con;
    }

    if (as_fixparams.as_output_exttmp_ace) {
      std::string filename=basename+".ace";
      cout << "Logging this contig to file: " << filename << endl;

      std::ofstream aceout(filename, std::ios::out | std::ios::trunc);
      Contig::setCoutType(Contig::AS_ACE);
      aceout << con;
    }


    if (as_fixparams.as_output_exttmp_fasta) {
      std::string filename=basename+".fasta";
      std::string qualname=filename+".qual";

      cout << "Logging this contig to files: " << filename << "  and  " << qualname << endl;

      {
	std::ofstream fastaout(filename, std::ios::out | std::ios::trunc);
	Contig::setCoutType(Contig::AS_FASTA);
	fastaout << con;
      }
      {
	std::ofstream qualout(qualname, std::ios::out | std::ios::trunc);
	Contig::setCoutType(Contig::AS_FASTAQUAL);
	qualout << con;
      }
    }


    if (as_fixparams.as_output_exttmp_gap4da) {
      std::string dirname=basename+".gap4da";

      cout << "Logging this contig to directory: " << dirname << endl;
      if(ensureDirectory(dirname,true)){
	MIRANOTIFY(Notify::FATAL, "Cannot make sure the directory exist? Aborting.");
      }

      Contig::setCoutType(Contig::AS_GAP4DA);
      std::ofstream fofnout((dirname+"/fofn"), std::ios::out | std::ios::trunc);
      con.saveAsGAP4DA(dirname, fofnout);
    }
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::buildDefaultCheckpointFileName(const std::string & filename)
{
  return AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/"+filename;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::buildDefaultInfoFileName(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, const std::string & defaultname, const std::string & defaultextension, bool removeold)
{
  std::string dirname;
  if(version>=0){
    dirname=AS_miraparams[0].getDirectoryParams().dir_tmp;
  }else{
    dirname=AS_miraparams[0].getDirectoryParams().dir_info;
  }

  std::string filename;
  if(basename.size()){
    filename=buildFileName(version, prefix, postfix,
			   basename, defaultextension,
			   "",
			   removeold);
  }else{
    filename=buildFileName(version, prefix, postfix,
			   defaultname, defaultextension,
			   dirname,
			   removeold);
  }

  return filename;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::buildDefaultResultsFileName(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, const std::string & defaultname, const std::string & defaultextension, bool removeold)
{
  std::string dirname;
  if(version>=0){
    dirname=AS_miraparams[0].getDirectoryParams().dir_tmp;
  }else{
    dirname=AS_miraparams[0].getDirectoryParams().dir_results;
  }

  std::string filename;
  if(basename.size()){
    filename=buildFileName(version, prefix, postfix,
			   basename, defaultextension,
			   "",
			   removeold);
  }else{
    filename=buildFileName(version, prefix, postfix,
			   defaultname, defaultextension,
			   dirname,
			   removeold);
  }

  return filename;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getContigReadListFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_crlist,
    ".txt");
}

void Assembly::saveContigReadList(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getContigReadListFilename(version, prefix, postfix, basename));
  assout::saveContigReadList(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getStatisticsFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_contigstats,
    ".txt");
}

void Assembly::saveStatistics(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getStatisticsFilename(version, prefix, postfix, basename));
  assout::saveStatistics(AS_contigs,filename, deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getAssemblyInfoFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_info,
    ".txt");
}

void Assembly::saveAssemblyInfo(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getAssemblyInfoFilename(version, prefix, postfix, basename));
  assout::saveAssemblyInfo(AS_assemblyinfo,filename, deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getLargeContigsInfoFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_largecontigs,
    ".txt");
}

void Assembly::saveLargeContigsInfo(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getLargeContigsInfoFilename(version, prefix, postfix, basename));
  assout::saveLargeContigsInfo(AS_assemblyinfo,filename, deleteoldfile);
}

/*************************************************************************
 *
 * readgroup stat & template info
 *
 *
 *************************************************************************/

std::string Assembly::getRGSTInfoFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_rgstinfo,
    ".txt");
}

void Assembly::saveRGSTInfo(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getRGSTInfoFilename(version, prefix, postfix, basename));

  cout << "Saving readgroup info to file: " << filename << endl;
  std::ofstream fout(filename, std::ios::out | std::ios::trunc);

  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    fout << "Readgroup " << rgi;
    if(rgid.getGroupName().empty()){
      fout << " (unnamed)\n";
    }else{
      fout << " (" << rgid.getGroupName() << ")\n";
    }
    fout << "------------------------------------------\n";
    fout << "Technology: " << ReadGroupLib::getNameOfSequencingType(rgid.getSequencingType()) << '\n';
    fout << "Num. reads: " << AS_current_rls_byrg[rgi].RLS_count_all << '\n';
    fout << "  thereof chimera: "
	 << AS_current_rls_byrg[rgi].RLS_count_chimera
	 << " (";
    if(AS_current_rls_byrg[rgi].RLS_count_chimera==0){
      fout << "maybe not searched for?";
    }else{
      auto tmpratio=100.0/AS_current_rls_byrg[rgi].RLS_count_all*AS_current_rls_byrg[rgi].RLS_count_chimera;
      fout << boost::str(boost::format("%.2f") % tmpratio);
      fout << "%, " << AS_current_rls_byrg[rgi].RLS_verdict_chimera;
    }
    fout << ")\n";
    fout << "\nLibrary type: ";
    if(rgid.hasTemplateInfo()) {
      fout << "paired "
	   << ReadGroupLib::getNameOfSegmentplacement(AS_rgstinfo[rgi].rgtguess.tg.splace_seen)
	   << "\n    mean: " << AS_rgstinfo[rgi].rgtguess.mean
	   << " stdev: " << AS_rgstinfo[rgi].rgtguess.stdev
	   << " skewness: " << AS_rgstinfo[rgi].rgtguess.skewness
	   << "\n    min: " << rgid.getInsizeFrom()
	   << "\tmax: " << rgid.getInsizeTo()
	   << endl;
    }else{
      fout << "shotgun\n";
    }
    fout << "\n\n";
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveDebrisList(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  FUNCSTART("void Assembly::saveDebrisInfo(int32 version, const std::string & prefix, const std::string & postfix, const std::string & debrisfilename)");

  std::string filename(buildDefaultInfoFileName(
		    version, prefix, postfix, basename,
		    AS_miraparams[0].getAssemblyParams().as_outfile_stats_debrislist,
		    ".txt"));

  cout << "Saving debris list to file: " << filename << endl;
  std::ofstream fout(filename, std::ios::out | std::ios::trunc);

  for(uint32 i=0; i< AS_readpool.size(); i++){
    if(AS_isdebris[i]) {
      fout << AS_readpool.getRead(i).getName();
      switch(AS_isdebris[i]) {
      case DEBRIS_NOTDEBRIS : {
	// NOP, can never happen but quieten compiler warning
	break;
      }
      case DEBRIS_SHORTONLOAD : {
	fout << "\tSHORTONLOAD\n";
	break;
      }
      case DEBRIS_UNSPECIFIED : {
	fout << "\tUNSPECIFIED\n";
	break;
      }
      case DEBRIS_NOOVERLAP : {
	fout << "\tNO_OVERLAP\n";
	break;
      }
      case DEBRIS_MEGAHUB : {
	fout << "\tFILTER_MEGAHUB\n";
	break;
      }
      case DEBRIS_MASKEDNASTYREPEAT : {
	fout << "\tMASKED_NASTY_REPEAT\n";
	break;
      }
      case DEBRIS_MASKEDHAF7REPEAT : {
	fout << "\tMASKED_HAF7_REPEAT\n";
	break;
      }
      case DEBRIS_MASKEDHAF6REPEAT : {
	fout << "\tMASKED_HAF6_REPEAT\n";
	break;
      }
      case DEBRIS_NOTMAPPED : {
	fout << "\tNOT_MAPPED\n";
	break;
      }
      case DEBRIS_ABORTEDCONTIGCREATION : {
	fout << "\tABORTED_CONTIG_CREATION\n";
	break;
      }
      case DEBRIS_TINYCONTIG : {
	fout << "\tTINY_CONTIG\n";
	break;
      }
      case DEBRIS_TINYCLUSTER : {
	fout << "\tTINY_CLUSTER\n";
	break;
      }
      case DEBRIS_TINYCLUSTERORPHAN : {
	fout << "\tTINY_CLUSTER_ORPHAN\n";
	break;
      }
      case DEBRIS_UNSAVEDSINGLET : {
	fout << "\tUNSAVED_SINGLET\n";
	break;
      }
      case DEBRIS_DIGITAL_NORMALISATION : {
	fout << "\tDIGITAL_NORMALISATION\n";
	break;
      }
      case DEBRIS_CLIP_BADSOLEXAEND : {
	fout << "\tCLIP_BAD_SOLEXA_END\n";
	break;
      }
      case DEBRIS_CLIP_KNOWNADAPTORRIGHT : {
	fout << "\tCLIP_KNOWNADAPTORRIGHT\n";
	break;
      }
      case DEBRIS_CLIP_QUALMINTHRESHOLD : {
	fout << "\tCLIP_QUALMINTHRESHOLD\n";
	break;
      }
      case DEBRIS_CLIP_LOWERCASEFRONT : {
	fout << "\tCLIP_LOWERCASEFRONT\n";
	break;
      }
      case DEBRIS_CLIP_LOWERCASEBACK : {
	fout << "\tCLIP_LOWERCASEBACK\n";
	break;
      }
      case DEBRIS_CLIP_QUALCLIPS : {
	fout << "\tCLIP_QUALCLIPS\n";
	break;
      }
      case DEBRIS_CLIP_MASKEDBASES : {
	fout << "\tCLIP_MASKEDBASES\n";
	break;
      }
      case DEBRIS_CLIP_BADSEQUENCESERACH : {
	fout << "\tCLIP_BADSEQUENCESERACH\n";
	break;
      }
      case DEBRIS_CLIP_POLYBASEATEND : {
	fout << "\tCLIP_POLYBASEATEND\n";
	break;
      }
      case DEBRIS_CLIP_POLYAT : {
	fout << "\tCLIP_POLYAT\n";
	break;
      }
      case DEBRIS_CLIP_MINLEFTCLIP : {
	fout << "\tCLIP_MINLEFTCLIP\n";
	break;
      }
      case DEBRIS_CLIP_MINRIGHTCLIP : {
	fout << "\tCLIP_MINRIGHTCLIP\n";
	break;
      }
      case DEBRIS_CLIP_PHIX174 : {
	fout << "\tCLIP_PHIX174\n";
	break;
      }
      case DEBRIS_CLIP_RRNA : {
	fout << "\tCLIP_RRNA\n";
	break;
      }
      case DEBRIS_CLIP_RRNA_PAIR : {
	fout << "\tCLIP_RRNA_PAIR\n";
	break;
      }
      case DEBRIS_CLIP_PROPOSEDENDCLIP : {
	fout << "\tCLIP_PROPOSEDENDCLIP\n";
	break;
      }
      case DEBRIS_CLIP_CHIMERA : {
	fout << "\tCLIP_CHIMERA\n";
	break;
      }
      case DEBRIS_CLIP_TERMINALLYINCORRECTIBLEORCHIMERA : {
	fout << "\tCLIP_TERMINALLYINCORRECTIBLEORCHIMERA\n";
	break;
      }
      case DEBRIS_CLIP_INCORRECTIBLEENDORCHIMERA : {
	fout << "\tCLIP_INCORRECTIBLEENDORCHIMERA\n";
	break;
      }
      default : {
	fout << "\tNO_CODE_YET?_" << static_cast<uint16>(AS_isdebris[i]) << "\n";
      }
      }
    }
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getReadTagListFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_readtags,
    ".txt");
}

void Assembly::saveReadTagList(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getReadTagListFilename(version, prefix, postfix, basename));
  assout::saveReadTagList(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getConsensusTagListFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
   return buildDefaultInfoFileName(
     version, prefix, postfix, basename,
     AS_miraparams[0].getAssemblyParams().as_outfile_stats_contigtags,
     ".txt");
}

void Assembly::saveConsensusTagList(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
  {
  std::string filename(getConsensusTagListFilename(version, prefix, postfix, basename));
  assout::saveConsensusTagList(AS_contigs,filename,deleteoldfile);
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveSNPList(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(buildDefaultInfoFileName(
		    version, prefix, postfix, basename,
		    AS_miraparams[0].getAssemblyParams().as_outfile_stats_snpanalysis,
		    ".txt"));
  assout::saveSNPList(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveFeatureAnalysis(int32 version, const std::string & prefix, const std::string & postfix, const std::string & faname, const std::string & fsname, const std::string & fcname, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveFeatureAnalysis(int32 version, const std::string & prefix, const std::string & postfix, const std::string & faname, const std::string & faname, bool deleteoldfile)");

  std::string dirname;
  if(version>=0){
    dirname=AS_miraparams[0].getDirectoryParams().dir_tmp;
  }else{
    dirname=AS_miraparams[0].getDirectoryParams().dir_info;
  }

  std::string filenamea;
  if(faname.size()){
    filenamea=buildFileName(version, prefix, postfix, faname, ".txt");
  }else{
    filenamea=buildFileName(version, prefix, postfix,
			    AS_miraparams[0].getAssemblyParams().as_outfile_stats_featureanalysis,
			    ".txt",
			   dirname);
  }

  std::string filenames;
  if(fsname.size()){
    filenames=buildFileName(version, prefix, postfix, fsname, ".txt");
  }else{
    filenames=buildFileName(version, prefix, postfix,
			    AS_miraparams[0].getAssemblyParams().as_outfile_stats_featuresummary,
			    ".txt",
			   dirname);
  }

  std::string filenamec;
  if(fcname.size()){
    filenamec=buildFileName(version, prefix, postfix, fcname, ".txt");
  }else{
    filenamec=buildFileName(version, prefix, postfix,
			    AS_miraparams[0].getAssemblyParams().as_outfile_stats_featuresequences,
			    ".txt",
			   dirname);
  }

  assout::saveFeatureAnalysis(AS_contigs,AS_readpool,
			      filenamea,filenames,filenamec,
			      deleteoldfile);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getFASTAFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_FASTAUNPADDED,
    ".fasta");
}
std::string Assembly::getFASTAPaddedFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_FASTAPADDED,
    ".fasta");
}

void Assembly::saveAsFASTA(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getFASTAFilename(version, prefix, postfix, basename));
  std::string paddedfilename(getFASTAPaddedFilename(version, prefix, postfix, basename));
  assout::saveAsFASTA(AS_contigs,filename,paddedfilename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveStrainsAsFASTAQUAL(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveStrainsAsFASTAQUAL(int32 version, const std::string & prefix, const std::string & postfix, const std::string & fastaname)");

  std::string filename(buildDefaultResultsFileName(
		    version, prefix, postfix, basename,
		    AS_miraparams[0].getAssemblyParams().as_outfile_FASTAPADDED,
		    ""));
  assout::saveStrainsAsFASTAQ(AS_contigs,AS_readpool,
			      filename,
			      false,0,0,
			      deleteoldfile);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getTCSFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_TCS,
    ".tcs");
}
void Assembly::saveAsTCS(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsTCS(int32 version, const std::string & prefix, const std::string & postfix, const std::string & tcsname)");

  std::string filename(getTCSFilename(version, prefix, postfix, basename));
  assout::saveAsTCS(AS_contigs,filename,deleteoldfile);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getCAFFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_CAF,
    ".caf");
}

void Assembly::saveAsCAF(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getCAFFilename(version, prefix, postfix, basename));
  assout::saveAsCAF(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getMAFFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_MAF,
    ".maf");
}

void Assembly::saveAsMAF(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getMAFFilename(version, prefix, postfix, basename));
  assout::saveAsMAF(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getTXTFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_TXT,
    ".txt");
}

void Assembly::saveAsTXT(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string filename(getTXTFilename(version, prefix, postfix, basename));
  assout::saveAsTXT(AS_contigs,filename,deleteoldfile);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getACEFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_ACE,
    ".ace");
}
void Assembly::saveAsACE(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsACE(int32 version, const std::string & prefix, const std::string & postfix, const std::string & acename)");

  std::string filename(getACEFilename(version, prefix, postfix, basename));
  assout::saveAsACE(AS_contigs,filename,deleteoldfile);

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getWiggleFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_WIGGLE,
    ".wig");
}
void Assembly::saveAsWiggle(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsWiggle(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)");

  std::string filename(getWiggleFilename(version, prefix, postfix, basename));
  assout::saveAsWiggle(AS_contigs,filename,deleteoldfile,false);

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getGAP4DAFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outdir_GAP4DA,
    ".gap4da");
}

void Assembly::saveAsGAP4DA(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  std::string subdirname(getGAP4DAFilename(version, prefix, postfix, basename));
  assout::saveAsGAP4DA(AS_contigs,subdirname,deleteoldfile);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Assembly::getHTMLFilename(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_HTML,
    ".html");
}

void Assembly::saveAsHTML(int32 version, const std::string & prefix, const std::string & postfix, const std::string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsHTML(int32 version, const std::string & prefix, const std::string & postfix, const std::string & htmlname, bool deleteoldfile)");

  std::string filename(getHTMLFilename(version, prefix, postfix, basename));

  std::string projectname(AS_miraparams[0].getAssemblyParams().as_projectname_out);

  cout << "Saving contigs to file: " << filename << endl;

  assout::dumpContigListAsHTML(AS_contigs, filename, deleteoldfile, projectname);

  FUNCEND();
}
