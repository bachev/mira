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


#include "contig.H"

#include "io/annotationmappings.H"
#include "mira/gff_parse.H"

using std::cout;
using std::cerr;
using std::endl;


/*************************************************************************
 *
 * copy first all Genbank features of all reads into one big vector
 *  if allowedfeatures not empty: only the features named there
 *  if forbiddenfeatures not empty: do not take features named there
 *
 * reads that have  Fsrc, Fgen, FCDS, Fexn, or Fint can get additional
 *  intergenic regions simulated if wished
 *
 * list is sorted ascending before given back by contigfrom, to and lastly
 *  by identifier (Fsrc, Fgen, FCDS, rest)
 *
 *************************************************************************/
#define CEBUGF2(bla)
//#define CEBUGF2(bla)  {cout << bla; cout.flush();}

void ___bad_extStr(std::string & what, std::string & ext){
  if(!what.empty()) what+="; ";
  what+=ext;
}

void Contig::getGBFSummary(std::list<gbfsummary_t> & allGBfeatures, const std::vector<multitag_t::mte_id_t> & allowedfeatures, const std::vector<multitag_t::mte_id_t> & forbiddenfeatures, bool simulateintergenics, bool alsomiraconstags) const
{
  FUNCSTART("void Contig::getGBFSummary(std::list<gbfsummary_t> & allGBfeatures, const std::vector<multitag_t::mte_id_t> & allowedfeatures, const std::vector<multitag_t::mte_id_t> & forbiddenfeatures, bool simulateintergenics) const");

  allGBfeatures.clear();

  std::vector<int8> featuremask;
  uint32 igrcounter=0;

  gff3attributes_t gff3a;

  for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
    featuremask.clear();

    //CEBUGF2(pcrI->read);
    bool mustsimulateigr=false;
    if(simulateintergenics){
      for(auto & thetag : pcrI->getTags()){
	// with or without databank entry???
	if(thetag.identifier==Read::REA_tagentry_idSOFAdatabank_entry
	   || thetag.identifier==Read::REA_tagentry_idSOFAgene
	   || thetag.identifier==Read::REA_tagentry_idSOFACDS
	   || thetag.identifier==Read::REA_tagentry_idSOFAexon
	   || thetag.identifier==Read::REA_tagentry_idSOFAintron
	   || thetag.identifier==Read::REA_tagentry_idSOFAmRNA
	   || thetag.identifier==Read::REA_tagentry_idSOFAtranscript
	   || thetag.identifier==Read::REA_tagentry_idSOFAprimary_transcript
	   || thetag.identifier==Read::REA_tagentry_idSOFArRNA
	   || thetag.identifier==Read::REA_tagentry_idSOFAscRNA
	   || thetag.identifier==Read::REA_tagentry_idSOFAsnRNA
	   || thetag.identifier==Read::REA_tagentry_idSOFAtRNA
	  ){
	  mustsimulateigr=true;
	  break;
	}
      }
    }

    for(auto & thetag : pcrI->getTags()){
      CEBUGF2("la this: " << thetag << endl);
      if(thetag.isSourceMIRA()) continue;

      // next tag if tag is completely outside of the good read part
      if(!((thetag.from >= pcrI->getLeftClipoff()
	    && thetag.from < pcrI->getRightClipoff())
	   || (thetag.to >= pcrI->getLeftClipoff()
	       && thetag.to < pcrI->getRightClipoff()))) continue;

      bool musttakeit=true;
      if(!allowedfeatures.empty()) {
	musttakeit=false;
	for(uint32 i=0;i<allowedfeatures.size();i++){
	  if(thetag.identifier == allowedfeatures[i]){
	    musttakeit=true;
	    break;
	  }
	}
      }


      if(!forbiddenfeatures.empty()) {
	for(uint32 i=0;i<forbiddenfeatures.size();i++){
	  if(thetag.identifier == forbiddenfeatures[i]){
	    musttakeit=false;
	    break;
	  }
	}
      }

      bool isgene=false;
      if(thetag.identifier==Read::REA_tagentry_idSOFAgene
	 || thetag.identifier==Read::REA_tagentry_idSOFACDS
	 || thetag.identifier==Read::REA_tagentry_idSOFAexon
	 || thetag.identifier==Read::REA_tagentry_idSOFAintron
	 || thetag.identifier==Read::REA_tagentry_idSOFAmRNA
	 || thetag.identifier==Read::REA_tagentry_idSOFAtranscript
	 || thetag.identifier==Read::REA_tagentry_idSOFAprimary_transcript
	 || thetag.identifier==Read::REA_tagentry_idSOFArRNA
	 || thetag.identifier==Read::REA_tagentry_idSOFAscRNA
	 || thetag.identifier==Read::REA_tagentry_idSOFAsnRNA
	 || thetag.identifier==Read::REA_tagentry_idSOFAtRNA
	){
	isgene=true;
      }

      if(musttakeit) {
	gbfsummary_t onefeature(pcrI);
	std::string extracted;

	CEBUGF2("Have this:\n"<<thetag<<endl);

	onefeature.identifier=thetag.getIdentifierStr();

	onefeature.mustbetranslated=false;
	if(thetag.identifier==Read::REA_tagentry_idSOFACDS){
	  onefeature.mustbetranslated=true;
	}
	onefeature.cfrom=pcrI.unclippedReadPos2ContigPos(thetag.from);
	onefeature.cto=pcrI.unclippedReadPos2ContigPos(thetag.to);

	onefeature.direction=thetag.getStrandDirection();
	onefeature.strainid=pcrI->getStrainID();

	if(thetag.commentisgff3){
	  CEBUGF2("Tag is GFF3, parsing attributes\n");
	  GFFParse::parseGFF3Attributes(thetag.getCommentStr(),gff3a);
	  onefeature.translationtable=11;

	  int8_t codonstart=-1; // means: not set

	  for(auto & gfe : gff3a){
	    if(!gfe.values.empty()){
	      if(gfe.tag=="locus_tag"){
		onefeature.locustag=gfe.values[0];
	      }else if(gfe.tag=="old_locus_tag"){
		onefeature.oldlocustag=gfe.values[0];
	      }else if(gfe.tag=="function"){
		___bad_extStr(onefeature.function,gfe.values[0]);
	      }else if(gfe.tag=="eC_number" || gfe.tag=="EC_number"){
		___bad_extStr(onefeature.ecnumber,gfe.values[0]);
	      }else if(gfe.tag=="product"){
		___bad_extStr(onefeature.product,gfe.values[0]);
	      }else if(gfe.tag=="note" or gfe.tag=="Note"){
		___bad_extStr(onefeature.note,gfe.values[0]);
	      }else if(gfe.tag=="inference"){
		___bad_extStr(onefeature.inference,gfe.values[0]);
	      }else if(gfe.tag=="gene_synonym"){
		___bad_extStr(onefeature.genesynonym,gfe.values[0]);
	      }else if(gfe.tag=="regulatory_class"){
		___bad_extStr(onefeature.regulatoryclass,gfe.values[0]);
	      }else if(gfe.tag=="gO_process" || gfe.tag=="GO_process"){
		___bad_extStr(onefeature.goprocess,gfe.values[0]);
	      }else if(gfe.tag=="gO_function" || gfe.tag=="GO_function"){
		___bad_extStr(onefeature.gofunction,gfe.values[0]);
	      }else if(gfe.tag=="gO_component" || gfe.tag=="GO_component"){
		___bad_extStr(onefeature.gocomponent,gfe.values[0]);
	      }else if(gfe.tag=="transl_table"){
		onefeature.translationtable=static_cast<int8>(atoi(gfe.values[0].c_str()));
	      }else if(gfe.tag=="codon_start"){
		codonstart=static_cast<int8>(atoi(gfe.values[0].c_str()));
	      }
	    }
	  }
	  if(thetag.phase!=3) codonstart=thetag.phase+1;  // MIRA tag V2 phase takes precedence (should be same anyway)
	  if(codonstart>0) onefeature.codonstart=codonstart;
	  onefeature.gene=GFFParse::extractCommonName(thetag.getCommentStr(),true);
	}else{
	  CEBUGF2("Tag is not GFF3, simply adding\n");
	  onefeature.note="Note="+thetag.getCommentStr();
	}

	onefeature.isgene=isgene;

	CEBUGF2("Extracted feature:\n" << onefeature);

	allGBfeatures.push_back(onefeature);
      }

      // preparing for intergenic regions: mask away this feature from read
      if(isgene){
	if(featuremask.empty()) featuremask.resize(pcrI->getLenSeq(),1);
	for(uint32 pos=thetag.from; pos<=thetag.to; pos++){
	  if(pos<pcrI->getLenSeq()) featuremask[pos]=0;
	}
      }
    }

    // if there were GenBank features in this read, simulate intergenic
    //  features
    if(!featuremask.empty() && mustsimulateigr){
      CEBUGF2("mustsimulateigr: " << mustsimulateigr << '\n');
      CEBUGF2("Simulate IGR for " << pcrI->getName() << '\n');
      bool inintergenic=false;
      uint32 runstart=0;
      for(uint32 pos=0; pos <= pcrI->getLenSeq(); pos++){
	if(pos < pcrI->getLenSeq()){
	  CEBUGF2("featuremask[" << pos << "]:\t" << static_cast<uint16>(featuremask[pos]) << endl);
	}
	if(pos == pcrI->getLenSeq() || featuremask[pos]==0){
	  if(inintergenic){
	    gbfsummary_t onefeature(pcrI);

	    onefeature.identifier="Figr";

	    onefeature.mustbetranslated=false;
	    int32 tmp=pcrI.unclippedReadPos2ContigPos(runstart);
	    if(tmp<0) tmp=0;
	    onefeature.cfrom=tmp;
	    if(pos == pcrI->getLenSeq()){
	      // carefull: if at end of sequence, we must take -1
	      //  or else the tag will be one too long
	      onefeature.cto=pcrI.unclippedReadPos2ContigPos(pos-1);
	    }else{
	      onefeature.cto=pcrI.unclippedReadPos2ContigPos(pos);
	    }
	    onefeature.direction=1;
	    onefeature.strainid=pcrI->getStrainID();
	    onefeature.isgene=false;
	    allGBfeatures.push_back(onefeature);

	    CEBUGF2("Made IGR:\n" << onefeature);

	    inintergenic=false;
	    runstart=0;
	  }
	}else{
	  if(!inintergenic){
	    inintergenic=true;
	    runstart=pos;
	  }
	}
      }
    }
  }

  CEBUGF2("going to sort 1\n");
  allGBfeatures.sort(Contig::gbfsummary_t_comparator);
  CEBUGF2("sorted\n");

  // assign names to the intergenic regions
  if(simulateintergenics){
    std::string gene;
    std::string function;
    std::string ecnumber;
    std::string product;
    std::string note;

    std::string concatstring="; ";

    for(auto gbfsI=allGBfeatures.begin(); gbfsI != allGBfeatures.end(); ++gbfsI){
      if(gbfsI->identifier == "Figr"){
      	std::ostringstream newlocusostr;
      	std::ostringstream newnoteostr;

      	newlocusostr << "IGR_";
      	newnoteostr << "Intergenic between ";

CEBUGF2("#1\n");

      	bool foundprevious=false;
      	if(gbfsI != allGBfeatures.begin()){
	  auto gbfsImm=gbfsI;

CEBUGF2("#2\n");
	  while (gbfsImm != allGBfeatures.begin()
		 && gbfsI->pcrI == gbfsImm->pcrI
		 && (gbfsImm == gbfsI || gbfsImm->isgene==false)) gbfsImm--;

CEBUGF2("#3\n");
	  if(gbfsI->pcrI == gbfsImm->pcrI) {
	    newlocusostr << gbfsImm->locustag << '_';
CEBUGF2("#4\n");
	    concatAllGBFInfoForLocus(allGBfeatures, gbfsImm, concatstring, gene, function, ecnumber, product, note);
CEBUGF2("#5\n");
	    if(gene.empty()){
	      newnoteostr << gbfsImm->locustag << " and ";
	    }else{
	      newnoteostr << gene << " and ";
	    }
CEBUGF2("#6\n");
	    foundprevious=true;
	  }
	}
CEBUGF2("#7\n");
      	if(!foundprevious){
      	  newlocusostr << "seqstart_";
      	  newnoteostr << "sequence start and ";
      	}

      	if(gbfsI != allGBfeatures.end()) {
CEBUGF2("#8\n");
	  auto gbfsImm=gbfsI;
	  while (gbfsImm != allGBfeatures.end()
		 && gbfsI->pcrI == gbfsImm->pcrI
		 && (gbfsImm == gbfsI || gbfsImm->isgene==false)) gbfsImm++;

CEBUGF2("#9\n");
	  foundprevious=false;
	  if(gbfsImm != allGBfeatures.end()
	     && gbfsI->pcrI == gbfsImm->pcrI) {
CEBUGF2("#a\n");
	    newlocusostr << gbfsImm->locustag;
CEBUGF2("#b\n");
	    concatAllGBFInfoForLocus(allGBfeatures, gbfsImm, concatstring, gene, function, ecnumber, product, note);
CEBUGF2("#c\n");
	    if(gene.empty()){
	      newnoteostr << gbfsImm->locustag;
	    }else{
	      newnoteostr << gene;
	    }
CEBUGF2("#d\n");
	    foundprevious=true;
	  }
	}
	if(!foundprevious){
      	  newlocusostr << "seqend";
      	  newnoteostr << "sequence end";
      	}

CEBUGF2("#e\n");
      	gbfsI->locustag=newlocusostr.str();
      	gbfsI->note = newnoteostr.str();
CEBUGF2("#f\n");
      }
    }
    CEBUGF2("going to sort 2\n");
    allGBfeatures.sort(Contig::gbfsummary_t_comparator);
    CEBUGF2("sorted\n");
  }

  {
    // silly gcc 4.6.x does not do "gbfsummary_t onefeature(CON_reads.end());" *sigh*
    auto tmp=CON_reads.end();
    bool hasnewtag=false;
    for(auto & cte : CON_consensus_tags){
      if((!cte.getSourceStr().empty() && !cte.isSourceMIRA())
	 || (alsomiraconstags)){
//	     (AnnotationMappings::isMIRAEntry(cte.getIdentifierStr())
//	      || AnnotationMappings::isValidGAP4Entry(cte.getIdentifierStr())))){
	gbfsummary_t onefeature(tmp);
	if(cte.identifier == CON_tagentry_idSAOc
	   || cte.identifier == CON_tagentry_idSIOc
	   || cte.identifier == CON_tagentry_idSROc
	   || cte.identifier == CON_tagentry_idSRMc
	   || cte.identifier == CON_tagentry_idIUPc
	   || cte.identifier == CON_tagentry_idUNSc) {
	  onefeature.identifier="sequence_variant";
	}else if(cte.identifier == CON_tagentry_idMCVc) {
	  onefeature.identifier="gap";
	}else{
	  //onefeature.identifier="region";
	}
	if(onefeature.identifier.empty()) continue;
	onefeature.cfrom=cte.from;
	onefeature.cto=cte.to;
	onefeature.direction=1;
	if(cte.commentisgff3){
	  onefeature.note=GFFParse::extractKeytag("Note",cte.getCommentStr());
	}else{
	  onefeature.note=cte.getCommentStr();
	}
	allGBfeatures.push_back(onefeature);
	hasnewtag=true;

	//cout << "Seen: " << cte << endl;
	//cout << "Pushed back: " << onefeature << endl;
      }
    }

    if(hasnewtag) {
      CEBUGF2("going to sort 3\n");
      allGBfeatures.sort(Contig::gbfsummary_t_comparator);
      CEBUGF2("sorted\n");
    }
  }
}
//#define CEBUGF2(bla)


/*************************************************************************
 *
 * append b to a if b != a;
 * if a was not empty, also add a concat string inbetween
 *
 *************************************************************************/

void Contig::myappend(std::string & a, const std::string & b, const std::string & concatstring) const
{
  if(!b.empty() && b != a) {
    if(!a.empty()) a+=concatstring;
    a+=b;
  }
}


/*************************************************************************
 *
 * GenBank annotations are sometimes scattered across different identifiers
 *  of a locus. Starting from a known locus, this function collects
 *  all information for gene, function, ecnumber, product and note
 *
 * Expects the list to be sorted
 *
 *************************************************************************/

void Contig::concatAllGBFInfoForLocus(const std::list<gbfsummary_t> & allGBfeatures, std::list<gbfsummary_t>::const_iterator gbelementI, const std::string & concatstring, std::string & gene, std::string & function, std::string & ecnumber, std::string & product, std::string & note) const
{
  gene.clear();
  function.clear();
  ecnumber.clear();
  product.clear();
  note.clear();

  if(gbelementI == allGBfeatures.end()) return;

  auto gbfsI=gbelementI;
  // search first element with that locus
  while (gbfsI != allGBfeatures.begin() && gbfsI->cfrom >= gbelementI->cfrom) --gbfsI;
  ++gbfsI;
  auto gbfsE=gbfsI;
  while (gbfsE != allGBfeatures.end() && gbfsE->cfrom <= gbelementI->cto) ++gbfsE;

  for(; gbfsI != gbfsE; ++gbfsI){
    if(gbfsI->locustag == gbelementI->locustag
      && gbfsI->cfrom == gbelementI->cfrom
      && gbfsI->cto == gbelementI->cto) {
      myappend(gene, gbfsI->gene, concatstring);
      myappend(function, gbfsI->function, concatstring);
      myappend(ecnumber, gbfsI->ecnumber, concatstring);
      myappend(product, gbfsI->product, concatstring);
      myappend(note, gbfsI->note, concatstring);
    }
  }

  return;
}



/*************************************************************************
 *
 * go through reads of a contig and and fetch all SequenceOntology
 *  tags to the return list
 *
 * position of tags in list will be positions in the contig
 *
 * // HMMM list is sorted ascending before given back by contigfrom, to and lastly
 *  by identifier (Fsrc, Fgen, FCDS, rest)
 *
 *************************************************************************/
//#define CEBUGF2(bla)  {cout << bla; cout.flush();}

#define CEBUGF2(bla)
void Contig::getSeqOntTags(std::list<contigSOtag_t> & allSOfeatures)
{
  FUNCSTART("void Contig::getSeqOntTags(std::list<contigSOtag_t> & allSOfeatures)");

  allSOfeatures.clear();

  contigSOtag_t emptycsot;

  for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
    if(pcrI.getORPID() < 0) continue;
    // for testing:
    // pcrI->sortTagsForGFF3();
    for(uint32 ti=0; ti < pcrI->getNumOfTags(); ++ti){
      const multitag_t & acttag=pcrI->getTag(ti);

      if(!acttag.identifierIsValidGFF3SOEntry()) continue;

      // next tag if tag is completely outside of the good read part
      if(!((acttag.from >= pcrI->getLeftClipoff()
	    && acttag.from < pcrI->getRightClipoff())
	   || (acttag.to >= pcrI->getLeftClipoff()
	       && acttag.to < pcrI->getRightClipoff()))) continue;
      allSOfeatures.push_back(emptycsot);
      allSOfeatures.back().multitag=acttag;
      allSOfeatures.back().multitag.from=pcrI.unclippedReadPos2ContigPos(acttag.from);
      allSOfeatures.back().multitag.to=pcrI.unclippedReadPos2ContigPos(acttag.to);
      // our tags might be on reads which are reversed (e.g. when on assembled reads instead the backbone),
      //  which would lead to the multitag "to" < "from". Not good.
      // Therefore take care that the multitag gets those positions right
      // TODO: this is too fragile. Multitag should have accessor set/get functions which check that
      if(allSOfeatures.back().multitag.to < allSOfeatures.back().multitag.from){
	std::swap(allSOfeatures.back().multitag.to,allSOfeatures.back().multitag.from);
      }
      BUGIFTHROW(allSOfeatures.back().multitag.from>=CON_counts.size(),"allSOfeatures.back().multitag.from>=CON_counts.size() ???");
      BUGIFTHROW(allSOfeatures.back().multitag.to>=CON_counts.size(),"allSOfeatures.back().multitag.from>toCON_counts.size() ???");
    }
  }
}
#define CEBUGF2



/*************************************************************************
 *
 * return in "result":
 *   first in list: values for whole contig
 *   then:          values for each feature
 *
 *
 *************************************************************************/

void Contig::calcSOTagCoverage(const std::list<contigSOtag_t> & features, std::list<tagcoverageinfo_t> & result)
{
  static tagcoverageinfo_t emptytci;

  tagcoverageinfo_t wholecontig_tci;
  std::vector<uint64> covvals;

  result.clear();
  covvals.reserve(100000);

  {
    // calc values for whole contig;
    contigSOtag_t tmp;
    tmp.multitag.from=0;
    tmp.multitag.to=static_cast<uint32>(CON_counts.size()-1);
    tmp.multitag.identifier=Read::REA_tagentry_idSOFAcontig;
    tmp.multitag.source=multitag_t::MT_tagsrcentry_idMIRA;
    tmp.multitag.setCommentStr("Note=MIRA: second order coverage values for whole contig");
    calcSOTagCoverage_helper(tmp, wholecontig_tci, covvals, wholecontig_tci);
  }
  calcSecondOrderStatsOnContainer(wholecontig_tci.ccinfo, covvals);
  result.push_back(wholecontig_tci);

  for(auto fI=features.begin(); fI!=features.end(); ++fI){
    result.push_back(emptytci);
    calcSOTagCoverage_helper(*fI, result.back(), covvals, wholecontig_tci);
  }

  return;
}


/*************************************************************************
 *
 * ftci : feature tci
 * ctci : comparator tci (usually whole contig)
 *
 *
 *************************************************************************/

void Contig::calcSOTagCoverage_helper(const contigSOtag_t & csot, tagcoverageinfo_t & ftci, std::vector<uint64> & covvals, const tagcoverageinfo_t & ctci)
{
  ftci.csot=csot;

  uint64 numvals=csot.multitag.to - csot.multitag.from+1;
  if(numvals>0){
    collectCoverage(csot.multitag.from, csot.multitag.to, covvals,false);
    calcStatsOnContainer(ftci.ccinfo, covvals);

    if(ctci.ccinfo.median > ftci.ccinfo.median){
      ftci.comparator_factor=-ctci.ccinfo.median/ftci.ccinfo.median;
    }else{
      ftci.comparator_factor=ftci.ccinfo.median/ctci.ccinfo.median;
    }

    if(static_cast<int64>(ftci.ccinfo.median + 1*ftci.ccinfo.stddev)
       < static_cast<int64>(ctci.ccinfo.median)-1*static_cast<int64>(ctci.ccinfo.stddev)){
      if(static_cast<int64>(ftci.ccinfo.median + 2*ftci.ccinfo.stddev)
	 < static_cast<int64>(ctci.ccinfo.median)-2*static_cast<int64>(ctci.ccinfo.stddev)){
	ftci.comparator_text="probable";
      }else{
	ftci.comparator_text="possible";
      }
      if(ftci.ccinfo.median < ctci.ccinfo.stddev) {
	ftci.comparator_text+="_deletion";
      }else{
	ftci.comparator_text+="_CNV_down";
      }
    }else if(static_cast<int64>(ftci.ccinfo.median - 1*ftci.ccinfo.stddev)
	     > static_cast<int64>(ctci.ccinfo.median)+1*static_cast<int64>(ctci.ccinfo.stddev)){
      if(static_cast<int64>(ftci.ccinfo.median - 2*ftci.ccinfo.stddev)
	 > static_cast<int64>(ctci.ccinfo.median)+1*static_cast<int64>(ctci.ccinfo.stddev)){
	ftci.comparator_text="probable";
      }else{
	ftci.comparator_text="possible";
      }
      ftci.comparator_text+="_CNV_up";
    }else{
      ftci.comparator_text="normal";
    }
  }

  return;
}
