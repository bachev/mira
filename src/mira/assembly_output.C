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

/*
 * functions to output results of assembly
 * sourced out from assembly.C to get smaller file sizes and better modularity
 */


#include "mira/assembly_output.H"


#include <boost/lexical_cast.hpp>
#include <functional>                  // for str_hash

#include "io/annotationmappings.H"
#include "util/fileanddisk.H"
#include "util/dptools.H"
#include "util/prettyprint_container.H"
#include "mira/gff_parse.H"

#include "util/stlimprove.H"


using std::cout;
using std::cerr;
using std::endl;


#define CEBUG(bla)
#define CEBUGF(bla)

#ifndef PUBLICQUIET
//#define CEBUGNPQ(bla)  {cout << bla; cout.flush();}
#define CEBUGNPQ(bla)
#else
#define CEBUGNPQ(bla)
#endif

//#define CEBUGNPQ(bla)  {cout << bla; cout.flush();}
//#define CEBUG(bla)  {cout << bla; cout.flush();}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveContigReadList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveContigReadList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  std::ofstream fout;
  if(!openFileForAppend(filename,fout,deleteoldfile)){
    Contig::dumpContigReadList_Head(fout);
  }

  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	try{
	  cle.dumpContigReadList_Body(fout);
	}
	catch (Notify n) {
	  cerr << "Error while dumping " << cle.getContigName() << ".\n";
	  n.handleError(THISFUNC);
	}
      }
    }
  }

  FUNCEND();
}

void assout::saveContigReadList(Contig & con, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveContigReadList(Contig & con, const std::string & filename)");

  std::ofstream fout;
  if(!openFileForAppend(filename,fout,deleteoldfile)){
    Contig::dumpContigReadList_Head(fout);
  }
  con.dumpContigReadList_Body(fout);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveStatistics(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveStatistics(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  //cout << "Saving project statistics to file: " << filename << endl;

  std::ofstream statout;
  if(!openFileForAppend(filename,statout, deleteoldfile)){
    Contig::dumpContigStatistics_Head(statout);
  }

  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	try{
	  cle.dumpContigStatistics_Body(statout);
	}
	catch (Notify n) {
	  cerr << "Error while dumping " << cle.getContigName() << ".\n";
	  n.handleError(THISFUNC);
	}
      }
    }
  }

  FUNCEND();
}

void assout::saveStatistics(Contig & con, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveStatistics(Contig & con, const std::string & filename)");

  std::ofstream statout;
  try{
    if(!openFileForAppend(filename,statout, deleteoldfile)){
      Contig::dumpContigStatistics_Head(statout);
    }
    con.dumpContigStatistics_Body(statout);
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveAssemblyInfo(AssemblyInfo & asi, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveAssemblyInfo(AssemblyInfo & asi, const std::string & filename, bool deleteoldfile)");

  (void) deleteoldfile;

  std::ofstream infoout;
  openFileForAppend(filename,infoout, true);
  dateStamp(infoout);
  infoout << '\n' << asi;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveLargeContigsInfo(AssemblyInfo & asi, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveLargeContigsInfo(AssemblyInfo & asi, const std::string & filename, bool deleteoldfile)");

  (void) deleteoldfile;

  std::ofstream infoout;
  openFileForAppend(filename,infoout, true);
  asi.dumpLargeContigNames(infoout);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveReadTagList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveReadTagList(std::list<Contig> & clist, const std::string & filename)");

  //cout << "Saving read tag list to file: " << filename << endl;

  std::ofstream rtout;
  if(!openFileForAppend(filename,rtout, deleteoldfile)){
    Contig::dumpReadTagList_Head(rtout);
  }

  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	try{
	  cle.dumpReadTagList_Body(rtout);
	}
	catch (Notify n) {
	  cerr << "Error while dumping " << cle.getContigName() << ".\n";
	  n.handleError(THISFUNC);
	}
      }
    }
  }

  FUNCEND();
}

void assout::saveReadTagList(Contig & con, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveReadTagList(Contig & con, const std::string & filename)");

  std::ofstream rtout;
  try{
    if(!openFileForAppend(filename,rtout,deleteoldfile)){
      Contig::dumpReadTagList_Head(rtout);
    }
    con.dumpReadTagList_Body(rtout);
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveConsensusTagList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveConsensusTagList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  //cout << "Saving contig tag list to file: " << filename << endl;

  std::ofstream ctout;
  if(!openFileForAppend(filename,ctout,deleteoldfile)){
    Contig::dumpConsensusTagList_Head(ctout);
  }

  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	try{
	  cle.dumpConsensusTagList_Body(ctout);
	}
	catch (Notify n) {
	  cerr << "Error while dumping " << cle.getContigName() << ".\n";
	  n.handleError(THISFUNC);
	}
      }
    }
  }

  FUNCEND();
}

void assout::saveConsensusTagList(Contig & con, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveConsensusTagList(Con & con, const std::string & filename)");

  std::ofstream ctout;
  try{
    if(!openFileForAppend(filename,ctout,deleteoldfile)){
      Contig::dumpConsensusTagList_Head(ctout);
    }
    con.dumpConsensusTagList_Body(ctout);
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::getPreviousLocusTag(const std::list<gbfsummary_t> & gbfs, std::list<gbfsummary_t>::const_iterator gbfsI, std::string & locustag, std::string & gene)
{
  locustag=gbfsI->locustag;
  gene.clear();
  bool found=false;
  for(;!found;){
    if(gbfsI == gbfs.begin()) break;
    gbfsI--;
    if(gbfsI->locustag != locustag){
      found=true;
      locustag=gbfsI->locustag;
      gene=gbfsI->gene;
      while(gene.empty()){
	if(gbfsI == gbfs.begin()) break;
	gbfsI--;
	if(gbfsI->locustag==locustag) {
	  gene=gbfsI->gene;
	}else{
	  break;
	}
      }
    }
  }
  if(!found) locustag.clear();
}

void assout::getNextLocusTag(const std::list<gbfsummary_t> & gbfs, std::list<gbfsummary_t>::const_iterator gbfsI, std::string & locustag, std::string & gene)
{
  locustag=gbfsI->locustag;
  gene.clear();
  bool found=false;
  for(;!found;){
    gbfsI++;
    if(gbfsI == gbfs.end()) break;
    if(gbfsI->locustag != locustag){
      found=true;
      locustag=gbfsI->locustag;
      gene=gbfsI->gene;
      while(gene.empty()){
	gbfsI++;
	if(gbfsI == gbfs.end() || gbfsI->locustag != locustag) break;
	gene=gbfsI->gene;
      }
    }
  }
  if(!found) locustag.clear();
}

/*************************************************************************
 *
 * this function relies on the fact that Contig::getGBFSummary()
 *  is used with the option to fake intergenic features
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla; cout.flush();}
void assout::saveSNPList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveSNPList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  cout << "Saving SNP list to file: " << filename << endl;

  std::ofstream saout;
  if(!openFileForAppend(filename,saout,deleteoldfile)){
    saout << "#Note: this file is only useful when having assembled with data"
      "\n#containing GenBank features."
      "\n#A list of all SNP tags (SROc, SIOc, SAOc and eventually MCVc is"
      "\n#in the 'consensustaglist' file."
      "\n#Contig\t"
      "Contig position\t"
      "Reference\t"
      "Reference position\t"
      "Reference map position\t"
      "SNP type\t"
      "countAf\t"
      "countAr\t"
      "countCf\t"
      "countCr\t"
      "countGf\t"
      "countGr\t"
      "countTf\t"
      "countTr\t"
      "count*f\t"
      "count*r\t"
      "Feature type\t"
      "Hit locus\t"
      "Gene name\t"
      "EC num\t"
      "Product\t"
      "Function\t"
      "Note\n";
  }

  saout << std::setprecision(1) << std::fixed;

  // to get reference positions right for consensus bases, we need to get
  // the reference sequence
  std::string refseq;
  std::vector<base_quality_t> refqual;

  // search for the rg id with the refernce in it (may be not present)
  uint32 refid=1;
  for(; refid<ReadGroupLib::getNumReadGroups(); ++refid){
    auto rgid=ReadGroupLib::getReadGroupID(refid);
    if(rgid.isBackbone() && !rgid.isRail()) break;
  }
  // and refpos as re-usable vector for the loop further down;
  std::vector<uint32> refpos;

  for(auto & cle : clist){
    cle.sortConsensusTags();
    const PlacedContigReads & pcr = cle.getContigReads();
    for(auto & cr : pcr) {
      const_cast<Read &>(cr).sortTags();
    }
  }

  // we will look at all GBF features (empty allowed)
  std::vector<multitag_t::mte_id_t> allowedfeatures;
  // but not at Fsrc
  std::vector<multitag_t::mte_id_t> forbiddenfeatures;
  forbiddenfeatures.push_back(Read::REA_tagentry_idSOFAdatabank_entry);
  //forbiddenfeatures.push_back("FCDS");

  std::string gene;
  std::string function;
  std::string ecnumber;
  std::string product;
  std::string note;

  for(auto & cle : clist){
    try{
      // generate contig pos to reference translations
      // remember, the output is with 1 as first base!
      refpos.clear();
      refpos.resize(cle.getContigLength());
      if(refid==ReadGroupLib::getNumReadGroups()){
	// default translations
	uint32 count=1;
	for(auto & x : refpos) x=count++;
      }else{
	// translations if we HAVE a reference
	cle.newConsensusGet(refseq, refqual,
			    ReadGroupLib::getReadGroupID(refid).getStrainID());
	uint32 count=1;
	{
	  auto rsI=refseq.begin();
	  auto rpI=refpos.begin();
	  // jump over first Xs which mark non-covered reference sequence
	  for(;rsI!=refseq.end() && toupper(*rsI)=='X'; ++rsI, ++rpI){
	    *rpI=count;
	  }
	  // calc positions, gaps do not increase pos-count
	  for(;rsI!=refseq.end(); ++rsI, ++rpI){
	    *rpI=count;
	    if(*rsI!='*') ++count;
	  }
	}
	{
	  // and simulate a jump at the end of the sequence which might consist of Xs.
	  --count;
	  auto rsI=refseq.rbegin();
	  auto rpI=refpos.rbegin();
	  for(;rsI!=refseq.rend() && toupper(*rsI)=='X'; ++rsI, ++rpI){
	    *rpI=count;
	  }
	}
      }

      std::list<gbfsummary_t> allGBfeatures;
      cle.getGBFSummary(allGBfeatures,allowedfeatures,forbiddenfeatures, true, true);

      for(auto & cont : cle.getConsensusTags()) {
	if(!(cont.identifier==Contig::CON_tagentry_idSAOc
	     || cont.identifier==Contig::CON_tagentry_idSIOc
	     || cont.identifier==Contig::CON_tagentry_idSROc
	     )) continue;
	bool infeature=false;

	CEBUG("cont ID: " << cont.getIdentifierStr() << '\t' << cont.from << '\t' << cont.to << endl);

	auto lastfrom = cle.getContigLength()+1;
	auto lastto = lastfrom;
	for(auto gbfsI=allGBfeatures.cbegin(); gbfsI != allGBfeatures.cend(); ++gbfsI){
	  if(gbfsI->identifier == multitag_t::getIdentifierStr(Read::REA_tagentry_idSOFACDS)) continue;
	  CEBUG(' ' << gbfsI->identifier << '\t' << gbfsI->locustag << '\t' << gbfsI->cfrom << '\t' << gbfsI->cto << '\n');

	  // Added this hacky thing to suppress multiple output of the same sequence variant.
	  // MIRA currently tags SIOc and SAOc at the same places in the genome as
	  //  SAOc contains info about read counts while SIOc and SROc do not
	  // But this leads to 2 sequence_variant entries in the GBF list, meaning that a
	  //  SIOc and SAOc entry at one position each get double outputs, i.e., four lines
	  //  for one position. Not good.
	  // This gets it back to two lines.
	  if (gbfsI->cfrom == lastfrom && gbfsI->cto == lastto) {
	    CEBUG("Already worked on that, no need to redo.\n");
	    continue;
	  }
	  lastfrom=gbfsI->cfrom;
	  lastto=gbfsI->cto;


	  if(gbfsI->cfrom <= cont.from && gbfsI->cto >= cont.to){
	    CEBUG("Hit! gbfsI: " << *gbfsI
		  << "\ncont ID: " << cont.getIdentifierStr() << '\t' << cont.from << '\t' << cont.to << endl);

	    int32 adjposfroml=cont.from;
	    int32 adjposfromr=adjposfroml;

	    std::string wherename;
	    int32 wherelen=cle.getContigLength();
	    if(gbfsI->pcrI != cle.getContigReads().end()){
	      const Read & featureread=*(gbfsI->pcrI);
	      CEBUG("is in read: " << featureread.getName() << '\n');
	      //CEBUG(featureread);
	      wherename=featureread.getName();
	      int32 rreadposfrom=cle.getRealReadPos(cont.from,gbfsI->pcrI);
	      CEBUG("rreadposfrom: " << rreadposfrom << endl);
	      if(rreadposfrom>=0 && rreadposfrom<featureread.getLenClippedSeq()){
		adjposfroml=featureread.getLowerNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
		adjposfromr=featureread.getUpperNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
	      }
	    }else{
	      CEBUG("Is in consensus\n");
	      wherename="on consensus";
	      // get the positions from the reference reconstruction done earlier
	      for(; adjposfroml>0 && refseq[adjposfroml]=='*'; --adjposfroml);
	      adjposfroml=refpos[adjposfroml];
	      for(; adjposfromr<refseq.size()-1 && refseq[adjposfromr]=='*'; ++adjposfromr);
	      adjposfromr=refpos[adjposfromr];
	    }

	    std::string serialc;
	    infeature=true;

	    CEBUG("#1" << endl);

	    cle.concatAllGBFInfoForLocus(allGBfeatures, gbfsI, "; ", gene, function, ecnumber, product, note);

	    CEBUG("#2" << endl);

	    //<< '\t' << I->paddedPos2UnpaddedPos(cont.from)
	    saout << cle.getContigName()
		  << '\t' << cont.from
		  << '\t' << wherename
		  << '\t' << adjposfroml;
	    if(adjposfroml!=adjposfromr) saout << ':' << adjposfromr;
	    saout << '\t' << (1.0/wherelen*adjposfroml*360)
		  << '\t' << cont.getIdentifierStr()
		  << '\t' << cont.countACGTGapf[0]
		  << '\t' << cont.countACGTGapr[0]
		  << '\t' << cont.countACGTGapf[1]
		  << '\t' << cont.countACGTGapr[1]
		  << '\t' << cont.countACGTGapf[2]
		  << '\t' << cont.countACGTGapr[2]
		  << '\t' << cont.countACGTGapf[3]
		  << '\t' << cont.countACGTGapr[3]
		  << '\t' << cont.countACGTGapf[4]
		  << '\t' << cont.countACGTGapr[4]
		  << '\t' << gbfsI->identifier
		  << '\t' << gbfsI->locustag
		  << '\t' << gene
		  << '\t' << ecnumber
		  << '\t' << product
		  << '\t' << function
		  << '\t' << note
		  << endl;
	    //<< '\n';
	  }
	}

	if(!infeature) {
	  saout << cle.getContigName()
		<< '\t' << cont.from
		<< "\tn/a"
		<< "\tn/a"
		<< "\tn/a"
		<< '\t' << cont.getIdentifierStr()
		<< '\t' << cont.countACGTGapf[0]
		<< '\t' << cont.countACGTGapr[0]
		<< '\t' << cont.countACGTGapf[1]
		<< '\t' << cont.countACGTGapr[1]
		<< '\t' << cont.countACGTGapf[2]
		<< '\t' << cont.countACGTGapr[2]
		<< '\t' << cont.countACGTGapf[3]
		<< '\t' << cont.countACGTGapr[3]
		<< '\t' << cont.countACGTGapf[4]
		<< '\t' << cont.countACGTGapr[4]
		<< '\t'
		<< '\t'
		<< '\t'
		<< '\t'
		<< '\t'
		<< '\t'
		<< "\tNot covered by any GenBank feature"
		<< endl;
	}
      }
    }
    catch (Notify n) {
      cerr << "Error while dumping " << cle.getContigName() << ".\n";
      n.handleError(THISFUNC);
    }
    FUNCEND();
  }
}
//#define CEBUG(bla) {cout << bla; cout.flush();}





/*************************************************************************
 *
 * this function relies on tags having SO feature names
 *
 *************************************************************************/

void assout::saveCoverageInfo(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveCoverageInfo(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  cout << "Saving coverage info to file: " << filename << endl;

  std::ofstream saout;
  if(!openFileForAppend(filename,saout,deleteoldfile)){
    saout << "#Note: this file is only useful when having assembled with data"
      "\n#containing SequenceOntology features."
      "\n#Contig\t"
      "Feature ID\t"
      "Feature name\t"
      "Feature type\t"
      "Feature from\t"
      "Feature to\t"
      "Cov Min\t"
      "Cov Max\t"
      "Cov Mean\t"
      "Cov Median\t"
      "Cov StdDev\t"
      "PseudoReadCount\t"
      "Cov Status\t"
      "Cov Factor\t"
      "Product\t"
      "Function\t"
      "GO component\t"
      "GO function\t"
      "GO process\t"
      "Note\t"
      "GFF3-attributes\t"
      "\n"
      ;
  }

  std::list<contigSOtag_t> allSOtags;
  std::list<tagcoverageinfo_t> tagcovinfo;

  static const std::string idstring("ID");
  static const std::string locusstring("locus_tag");
  static const std::string prodstring("product");
  static const std::string funcstring("function");
  static const std::string gocostring("gO_component");
  static const std::string gofustring("gO_function");
  static const std::string goprstring("gO_process");
  static const std::string notestring("Note");
  static const std::string notestring2("note");

  for(auto & cle : clist){
    cle.getSeqOntTags(allSOtags);
    cle.calcSOTagCoverage(allSOtags,tagcovinfo);
    const tagcoverageinfo_t & contig_tci=tagcovinfo.front();

    for(auto & tce : tagcovinfo){
      saout << std::setprecision(1) << std::fixed;

      std::string id(GFFParse::extractKeytag(idstring,tce.csot.multitag.getCommentStr()));
      if(id.empty()){
	// no ID? Maybe because we loaded a GenBank file. Take locus_tag instead.
	id=GFFParse::extractKeytag(locusstring,tce.csot.multitag.getCommentStr());
      }
      std::string cname(GFFParse::extractCommonName(tce.csot.multitag.getCommentStr()));
      if(cname==id) cname.clear();
      std::string thisnote(GFFParse::extractKeytag(notestring,tce.csot.multitag.getCommentStr()));
      if(thisnote.empty()){
	// *sigh* Some files have "note=" instead of "Note="
	thisnote=GFFParse::extractKeytag(notestring2,tce.csot.multitag.getCommentStr());
      }

      // pseudo read count
      // how many perfect 100bp reads hit this area.
      auto pseudorc=static_cast<uint64>(tce.ccinfo.mean+0.5);
      {
	auto flen=tce.csot.multitag.to-tce.csot.multitag.from+1;
	if(flen >= 100){
	  pseudorc=(tce.ccinfo.sum+50)/100;
	}else{
	  // is there something to do?
	  //pseudorc*=flen;
	}
      }

      saout << cle.getContigName()
	    << "\t" << replaceEmptyString(id)
	    << "\t" << replaceEmptyString(cname)
	    << "\t" << replaceEmptyString(tce.csot.multitag.getIdentifierStr())
	    << "\t" << tce.csot.multitag.from+1
	    << "\t" << tce.csot.multitag.to+1
	    << "\t" << tce.ccinfo.min
	    << "\t" << tce.ccinfo.max
	    << "\t" << tce.ccinfo.mean
	    << "\t" << tce.ccinfo.median
	    << "\t" << tce.ccinfo.stddev
	    << "\t" << pseudorc
	    << "\t" << replaceEmptyString(tce.comparator_text)
	    << "\t" << std::setprecision(2) << std::fixed << tce.comparator_factor
	    << "\t" << htmlCode(replaceEmptyString(GFFParse::extractKeytag(prodstring,tce.csot.multitag.getCommentStr())))
	    << "\t" << htmlCode(replaceEmptyString(GFFParse::extractKeytag(funcstring,tce.csot.multitag.getCommentStr())))
	    << "\t" << htmlCode(replaceEmptyString(GFFParse::extractKeytag(gocostring,tce.csot.multitag.getCommentStr())))
	    << "\t" << htmlCode(replaceEmptyString(GFFParse::extractKeytag(gofustring,tce.csot.multitag.getCommentStr())))
	    << "\t" << htmlCode(replaceEmptyString(GFFParse::extractKeytag(goprstring,tce.csot.multitag.getCommentStr())))
	    << "\t" << htmlCode(replaceEmptyString(thisnote))
	    << "\t" << tce.csot.multitag.getCommentStr()
	    << "\n";
    }
  }
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla; cout.flush();}

void assout::saveSNPSurroundingAsHTML(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveSNPList(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  cout << "Saving SNP surroundings as HTML to file: " << filename << endl;

  std::ofstream htmlout;
  if(!openFileForAppend(filename,htmlout,deleteoldfile)){
    dumpHTMLHeader("SOME PROJECT NAME",htmlout);
  }
  htmlout << "<p>";

  //htmlout << std::setprecision(1) << std::fixed;

  for(auto & cle : clist){
    cle.sortConsensusTags();
    for(auto & cr : cle.getContigReads()) {
      const_cast<Read &>(cr).sortTags();
    }
  }

  CEBUG("Num contigs: " << clist.size() << '\n');

  uint32 surrounding=39;
  uint32 maxdistformerge=5;


  // we will look at all GBF features (empty allowed)
  std::vector<multitag_t::mte_id_t> allowedfeatures;
  // but not at Fsrc
  std::vector<multitag_t::mte_id_t> forbiddenfeatures;
  forbiddenfeatures.push_back(Read::REA_tagentry_idSOFAdatabank_entry);

  std::list<gbfsummary_t> allGBfeatures;

  // two loops: first to collect all names, then to
  //  print out navigation and snps

  std::list<snpenv_t> snpenvironments;
  uint32 snpjumpnamenr=0;

  snpenv_t dummysnpenv;
  dummysnpenv.from=0;
  dummysnpenv.to=0;
  dummysnpenv.numsnps=1;
  for(auto cI=clist.begin(); cI!=clist.end(); ++cI){
    cI->getGBFSummary(allGBfeatures,allowedfeatures,forbiddenfeatures, true, true);
    for(auto & cont : cI->getConsensusTags()) {
      if(!(cont.identifier==Contig::CON_tagentry_idSAOc
	   || cont.identifier==Contig::CON_tagentry_idSIOc
	   || cont.identifier==Contig::CON_tagentry_idSROc
	   || cont.identifier==Contig::CON_tagentry_idMCVc
	   )) continue;
      //
      CEBUG("Found tag: " << cont.getIdentifierStr() << "\t" << cont.from << '\n');
      if(snpenvironments.empty()){
	dummysnpenv.cI=cI;
	dummysnpenv.from=cont.from;
	dummysnpenv.to=cont.to;
	snpenvironments.push_back(dummysnpenv);
      }else{
	if(snpenvironments.back().cI == cI
	   && cont.to <= snpenvironments.back().to+maxdistformerge){
	  snpenvironments.back().to=cont.to;
	  snpenvironments.back().numsnps+=1;
	}else{
	  dummysnpenv.cI=cI;
	  dummysnpenv.from=cont.from;
	  dummysnpenv.to=cont.to;
	  snpenvironments.push_back(dummysnpenv);
	}
      }
    }
  }

  CEBUG("SNP environments: " << snpenvironments.size() << '\n');

  for(auto & see : snpenvironments){
    CEBUG("Checking seI: " << see);

    std::ostringstream snpcoordostr;
    std::string firstlocus="";
    std::string lastlocus="";
    int32 lowestpos=0x7fffffff;
    int32 highestpos=0;

    bool featurefound=false;
    for(auto gbfsI=allGBfeatures.begin(); gbfsI != allGBfeatures.end() && !featurefound; ++gbfsI){
      // careful, snpenv from/to is [...[ while gbfsI is [...] !
      if(see.to<=gbfsI->cfrom
	 || see.from > gbfsI->cto) continue;
      {
	CEBUG(' ' << gbfsI->identifier << '\t' << gbfsI->locustag << '\t' << gbfsI->cfrom << '\t' << gbfsI->cto << ' ' << gbfsI->note << endl);
	CEBUG("gbfsI: " << *gbfsI);
	CEBUG("Hit!\n");

	//// we can stop if we find a Fgen
	//if(gbfsI->identifier == "Fgen"){
	//  featurefound=true;
	//}

	see.gbfslist.push_back(*gbfsI);

	shortgbfinfo_t sgitmp;
	sgitmp.identifier=gbfsI->identifier;
	sgitmp.locustag=gbfsI->locustag;
	auto cI=see.cI;
	cI->concatAllGBFInfoForLocus(allGBfeatures, gbfsI, "; ", sgitmp.gene, sgitmp.function, sgitmp.ecnumber, sgitmp.product, sgitmp.note);
	see.sgilist.push_back(sgitmp);

	int32 adjposfroml=see.from;
	int32 adjposfromr=adjposfroml;

	// careful, gbfsummary_t may have been made by a consensus tag,
	//  therefore having a pcrI of end()
	if(gbfsI->pcrI != cI->getContigReads().end()){
	  const Read & featureread=*(gbfsI->pcrI);
	  int32 rreadposfrom=cI->getRealReadPos(see.from,gbfsI->pcrI);
	  CEBUG("rreadposfrom: " << rreadposfrom << endl);
	  if(rreadposfrom>=0 && rreadposfrom<featureread.getLenClippedSeq()){
	    adjposfroml=featureread.getLowerNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
	    adjposfromr=featureread.getUpperNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
	  }
	  CEBUG("yo" << endl);
	}
	if(adjposfroml<lowestpos){
	  lowestpos=adjposfroml;
	}
	if(adjposfromr>highestpos){
	  highestpos=adjposfromr;
	}
	if(!gbfsI->locustag.empty()){
	  if(firstlocus.empty()) firstlocus=gbfsI->locustag;
	  lastlocus=gbfsI->locustag;
	}
      }
    }
    if(firstlocus.empty()){
      snpcoordostr << "nodesc_" << lowestpos;
    }else{
      if(firstlocus==lastlocus){
	snpcoordostr << firstlocus << '_' << lowestpos;
      }else{
	snpcoordostr << firstlocus << '-' << lastlocus << '_' << lowestpos;
      }
    }
    if(lowestpos!=highestpos) snpcoordostr << ':' << highestpos;
    see.snpjumpname=snpcoordostr.str();
    CEBUG("jumpname: " << see.snpjumpname << endl);

    if(see.gbfslist.empty()) {
      CEBUG("WARNING: No feature found???\n" << see << endl);
    }
  }

  // for each snpenvironment:
  //  if there were gene or CDS: keep just those entries with the same locus tags
  //  else if there were entries with locus tags, just keep those entries

  for(auto & see : snpenvironments){
    bool hasGENE=false;
    bool hasCDS=false;
    bool hasLOCUS=false;
    auto gbfslI=see.gbfslist.cbegin();
    for(; gbfslI != see.gbfslist.end(); ++gbfslI){
      if(!gbfslI->locustag.empty()){
	hasLOCUS=true;
	if(gbfslI->identifier=="Fgen"){
	  hasGENE=true;
	}else if(gbfslI->identifier=="FCDS") {
	  hasCDS=true;
	}
      }
    }

    if(hasGENE || hasCDS || hasLOCUS){
      CEBUG("SOMETHING TO DELETE? " << hasGENE << hasCDS << hasLOCUS << see);
      auto sgiI=see.sgilist.cbegin();
      for(gbfslI=see.gbfslist.cbegin(); gbfslI != see.gbfslist.cend(); ++gbfslI, ++sgiI){
	bool deleteit=false;
	if(hasGENE){
	  if(gbfslI->identifier != "Fgen") deleteit=true;
	}else if(hasCDS){
	  if(gbfslI->identifier != "FCDS") deleteit=true;
	}else if(hasLOCUS){
	  if(gbfslI->locustag.empty()) deleteit=true;
	}
	if(deleteit){
	  CEBUG("DELETE gbfs: " << *gbfslI);
	  gbfslI=see.gbfslist.erase(gbfslI);
	  CEBUG("DELETE sgi: " << *sgiI);
	  sgiI=see.sgilist.erase(sgiI);
	}
      }
      CEBUG("AFTER DELETE: " << see);
    }else{
      CEBUG("NOTHING TO DELETE: " << see);
    }
  }


  // perfect, create jumpcoord strings etc.


  {
    std::string currentcontig;
    for(auto & see : snpenvironments){
      auto cI=see.cI;
      CEBUG("New seI. Contig: " << cI->getContigName() << endl);
      CEBUG(see);
      if(currentcontig != cI->getContigName()) {
	if(!currentcontig.empty()){
	  htmlout << "</table><p>";
	}
	currentcontig = cI->getContigName();
	htmlout << "<H1>Difference regions for " << currentcontig << "</H1>\n";

	{
	  uint32 numbbstrains=0;
	  std::vector<bool> seenstrains(cI->getContigReads().size(),false);
	  auto crI=cI->getContigReads().begin();
	  auto crE=cI->getContigReads().end();
	  for(; crI != crE; crI++){
	    if(crI.getORPID()>=0 && crI->isBackbone()){
	      if(!seenstrains[crI->getStrainID()]){
		numbbstrains++;
		seenstrains[crI->getStrainID()]=true;
	      }
	    }
	  }
	  if(numbbstrains>0){
	    if(numbbstrains==1) {
	      htmlout << "Backbone strain: ";
	    }else{
	      htmlout << "Backbone strains: ";
	    }
	    crI=cI->getContigReads().begin();
	    for(; crI != crE; crI++){
	      if(crI.getORPID()>=0){
		if(seenstrains[crI->getStrainID()]){
		  seenstrains[crI->getStrainID()]=false;
		  htmlout << crI->getStrainName() << ' ';
		}
	      }
	    }
	    htmlout << "<p>\n";
	  }

	  uint32 numstrains=0;
	  crI=cI->getContigReads().begin();
	  for(; crI != crE; crI++){
	    if(crI.getORPID()>=0 && crI->isBackbone()==false){
	      if(!seenstrains[crI->getStrainID()]){
		numstrains++;
		seenstrains[crI->getStrainID()]=true;
	      }
	    }
	  }
	  if(numstrains>0){
	    if(numbbstrains==1) {
	      htmlout << "Mapped strain: ";
	    }else{
	      htmlout << "Mapped strains: ";
	    }
	    crI=cI->getContigReads().begin();
	    for(; crI != crE; crI++){
	      if(crI.getORPID()>=0){
		if(seenstrains[crI->getStrainID()]){
		  seenstrains[crI->getStrainID()]=false;
		  htmlout << crI->getStrainName() << ' ';
		}
	      }
	    }
	    htmlout << "<p>\n";
	  }
	}


	htmlout << "<table BORDER=\"0\" CELLSPACING=\"10\" CELLPADDING=\"0\">\n";
      }

      htmlout << "<tr align=\"left\">"
	      << "<td CLASS=\"jtable1\"><a href=\"#"
	      << see.snpjumpname << "\">"
	      << see.snpjumpname << "</a></td>"
	      << "<td  CLASS=\"jtable2\">";
      uint32 num=0;
      for(auto sgiI=see.sgilist.begin(); sgiI!=see.sgilist.end(); sgiI++, num++){
	if(!(sgiI->gene).empty()){
	  if(num>0) htmlout << "; ";
	  htmlout << sgiI->gene;
	}
      }
      htmlout << "</td>"
	      << "<td  CLASS=\"jtable3\">";
      if(see.sgilist.size()==0) {
	// valid if left or right of annotated backbone
	//htmlout << "0 locus hit??? Strange.";
      } else if(see.sgilist.size()>1) {
	htmlout << "Multiple locus hit";
      }else{
	bool mustbr=false;
	if(!see.sgilist.begin()->function.empty()){
	  htmlout << see.sgilist.begin()->function;
	  mustbr=true;
	}
	if(!see.sgilist.begin()->product.empty()){
	  if(mustbr) htmlout << "<br>";
	  htmlout << see.sgilist.begin()->product;
	  mustbr=true;
	}
	if(!see.sgilist.begin()->note.empty()){
	  if(mustbr) htmlout << "<br>";
	  htmlout << see.sgilist.begin()->note;
	}

	//if(!see.sgilist.begin()->function.empty()){
	//  htmlout << see.sgilist.begin()->function;
	//}else{
	//  htmlout << see.sgilist.begin()->note;
	//}
      }
      htmlout << "</td></tr>\n";
    }
    if(!currentcontig.empty()){
      htmlout << "</table><p>";
    }
  }

  for(auto seI=snpenvironments.begin(); seI!=snpenvironments.end(); ++seI){
    auto cI=seI->cI;

    htmlout << "<H1>"
	    << "<a NAME=\"" << seI->snpjumpname << "\"></a>"
	    << seI->snpjumpname << "</H1>\n";
    if(seI != snpenvironments.begin()){
      seI--;
      htmlout << "<a href=\"#" << seI->snpjumpname << "\">";
      seI++;
    }
    htmlout << "previous";
    if(seI != snpenvironments.begin()){
      htmlout << "</a>";
    }
    htmlout << "  ";
    seI++;
    if(seI != snpenvironments.end()){
      htmlout << "<a href=\"#" << seI->snpjumpname << "\">";
    }
    htmlout << "next";
    if(seI != snpenvironments.end()){
      htmlout << "</a>";
    }
    seI--;

    htmlout << "<br/><br/>";

    uint32 num=0;
    for(auto sgiI=seI->sgilist.begin(); sgiI!=seI->sgilist.end(); sgiI++, num++){
      if(sgiI->identifier != "Figr"
	 && (!sgiI->locustag.empty() || !sgiI->gene.empty())){
	htmlout << "Hitting: " << sgiI->locustag;
	if(!sgiI->gene.empty()){
	  htmlout << " (" << sgiI->gene << ')';
	}else{
	  htmlout << " (" << *sgiI << ')';
	}
	htmlout << "<br>\n";
      }
      if(!sgiI->product.empty()){
	htmlout << "Product: " << sgiI->product << "<br>\n";
      }
      if(!sgiI->function.empty()){
	htmlout << "Function: " << sgiI->function << "<br>\n";
      }
      if(!sgiI->note.empty()){
	htmlout << "Note: " << sgiI->note << "<br>\n";
      }
    }


    cI->dumpAsHTML(htmlout, seI->from - surrounding, seI->to + surrounding);
  }

  CEBUG("\nDone with this contig" << endl);
}

#define CEBUG(bla)



/*************************************************************************
 *
 * result may contain '@' for positions not covered by a strain.
 *  if you do not want that, set fillholesinstraingenomes to true
 *
 *************************************************************************/

void assout::makeAllStrainGenomes(Contig & contigI, base_quality_t minqual, std::string & consseq, std::vector<base_quality_t> & consqual, std::vector<std::string> & strain_consseq, std::vector< std::vector<base_quality_t> > & strain_consqual, strainid2name_t & strainnames_in_contig, bool fillholesinstraingenomes)
{
  FUNCSTART("void makeAllStrainGenomes(Contig & contigI, std::string & consseq, std::vector<base_quality_t> & consqual, std::vector<std::string> & strain_consseq, std::vector< std::vector<base_quality_t> > & strain_consqual, strainid2name_t & strainnames_in_contig, bool fillholesinstraingenomes)");

  CEBUGNPQ("Making main genome consensus\n");

  // make a common consensus for all strains: mincoverage=0, minqual given
  contigI.calcConsensi(0,minqual,'X');

  // fetch the combined consensus
  contigI.newConsensusGet(consseq, consqual);

  CEBUGNPQ("Making strain genome consensi\n");

  // and make a specific consensus for every strain present in
  //  this contig
  // We'll do by fetching the strain specific consensus and apply, if wanted,
  //  the function to fill in gaps and Ns from the combined consensus
  strainnames_in_contig.clear();

  //// search the strains present in the contig
  //// TODO: contig now track strains present ... should use that!
  //for(auto & cre : contigI.getContigReads()){
  //  if(strainnames_in_contig.find(cre.getStrainID())==strainnames_in_contig.end()){
  //	strainnames_in_contig.insert(s_idname_pair_t(cre.getStrainID(),cre.getStrainName()));
  //	CEBUGNPQ(static_cast<int32>(cre.getStrainID()) << "\t-- " << cre.getStrain() << endl);
  //  }
  //}

  // No, take strains from ReadGroupLib
  for(uint32 id=1; id<ReadGroupLib::getNumReadGroups(); ++id){
    auto rgid=ReadGroupLib::getReadGroupID(id);
    if(strainnames_in_contig.find(rgid.getStrainID())==strainnames_in_contig.end()){
      strainnames_in_contig.insert(s_idname_pair_t(rgid.getStrainID(),rgid.getStrainName()));
      CEBUGNPQ(static_cast<int32>(rgid.getStrainID()) << "\t-- " << rgid.getStrainName() << endl);
    }
  }

  CEBUGNPQ("Number of strains present: " << strainnames_in_contig.size() << endl);
  strain_consseq.clear();
  strain_consqual.clear();

  if(strainnames_in_contig.size()>0) {
    auto SI=strainnames_in_contig.end();
    --SI;
    int32 largest=(SI->first)+1;
    strain_consseq.resize(largest);
    strain_consqual.resize(largest);

    // fill the consensus of the strains present
    SI=strainnames_in_contig.begin();
    for(;SI != strainnames_in_contig.end(); SI++) {
      contigI.newConsensusGet(strain_consseq[SI->first],
			      strain_consqual[SI->first],
			      SI->first);

      if(fillholesinstraingenomes) {
	// now merge the strain consensus:
	//  fill up N or holes in the strain consensus with bases from
	//  the main consensus
	BUGIFTHROW(consseq.size() != strain_consseq[SI->first].size(),"consseq.size() != strain_consseq[SI->first].size())?");
	//CEBUGNPQ("strain_consseq[" << (uint16) *SI << "]: " <<  strain_consseq[*SI] << endl);
	for(uint32 i=0; i<consseq.size(); i++){
	  auto up = static_cast<char>(toupper(strain_consseq[SI->first][i]));
	  if(up == 'X' || up == 'N' || up == '@'){
	    strain_consseq[SI->first][i]=static_cast<char>(tolower(consseq[i]));
	    strain_consqual[SI->first][i]=consqual[i];
	  }
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
struct codoninfo_t {
  std::string codon;
  uint32 count=0;
  double frequency=0.0;
  char   aa='?';
  uint8  index=0;

  friend std::ostream & operator<<(std::ostream &ostr, const codoninfo_t & ci){
    ostr << static_cast<uint16>(ci.index)
	 << '\t' << ci.codon
	 << '\t' << ci.aa
	 << '\t' << ci.count
	 << '\t' << ci.frequency;
    return ostr;
  }
  static inline bool sortByAAUpCountDown(const codoninfo_t & i1, const codoninfo_t & i2){
    if(i1.aa==i2.aa){
      return i1.count>i2.count;
    }
    return i1.aa<i2.aa;
  }
  static inline bool sortByIndexUp(const codoninfo_t & i1, const codoninfo_t & i2){
    return i1.index<i2.index;
  }
};


//#define CEBUGNPQ(bla) {cout << bla; cout.flush();}
void assout::saveFeatureAnalysis(std::list<Contig> & clist, ReadPool & rpool, const std::string & faname, const std::string & fsname, const std::string & fcname, bool deleteoldfile)
{
  FUNCSTART("saveFeatureAnalysis(std::list<Contig> & clist, ReadPool & rpool, const std::string & faname, const std::string & faname)");

  std::hash<std::string> str_hash;

  cout << "Saving feature analysis to file: " << faname << endl;
  std::ofstream faout;
  if(!openFileForAppend(faname,faout,deleteoldfile)){
    faout << "#Contig\t"
      "Reference strain\t"
      "Mutant strain\t"
      "Internal contig pos\t"
      "Position in reference\t"
      "Map pos in reference\t"
      "FType\t"
      "Locus\t"
      "Direction\t"
      "Gene name\t"
      "Distance of SNP\t"
      "DNA effect\t"
      "Nucleotide change\t"
      "Codon in reference\t"
      "Codon in mutant\t"
      "Codonfreq in reference\t"
      "Codonfreq in mutant\t"
      "Quality reference\t"
      "Quality mutant\t"
      "AA change\t"
      "AA position reference\t"
      "Codon position reference\t"
      "AA position mutant\t"
      "Codon position mutant\t"
      "Effect of mutation on protein\t"
      "GO process\t"
      "GO function\t"
      "Product\t"
      "Function\t"
      "Note\n";
  }

  faout << std::setprecision(1) << std::fixed;

  cout << "Saving impact summary to file: " << fsname << endl;
  std::ofstream fsout;
  if(!openFileForAppend(fsname,fsout,deleteoldfile)){
    fsout << "#Locus\t"
      "Gene name\t"
      "Feature type\t"
      "Feature checksum\t"
      "Genome map pos\t"

      "Interesting?\t"
      "Your own filter\t"

      "Coverage status\t"

      "First codon is start\t"

      "Changed start codon\t"
      "Destroyed start codon\t"
      "Changed stop codon\t"
      "Destroyed stop codon\t"
      "Premature stop codon\t"

      "Intergenic mutation\t"

      "Insertion in locus\t"
      "Deletion in locus\t"
      "Silent in locus\t"
      "AA change in locus\t"

      "Insertion untranslated\t"
      "Deletion untranslated\t"
      "Silent untranslated\t"
      "AA change untranslated\t"

      "GO process\t"
      "GO function\t"

      "Product\tFunction\tNote"
      "\n";
  }

  fsout << std::setprecision(1) << std::fixed;

  cout << "Saving feature sequences to file: " << fcname << endl;
  std::ofstream fcout;
  if(!openFileForAppend(fcname,fcout,deleteoldfile)){
    fcout << "#Contig\t"
      "FType\t"
      "Locus\t"
      "Gene name\t"
      "Feature strain name\t"
      "Mutant strain name 1 ...\t"
      "Protein in feature strain\t"
      "Protein in mutant strain 1 ...\t"
      "DNA in feature strain\t"
      "DNA in mutant strain 1 ...\n";
  }

  cout << "Collecting strain sequences and creating codon usage table" << endl;

  uint8 codontableused=0;
  std::vector<codoninfo_t> codoninfo(64);

  uint32 refstrainid=100;
  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rg=ReadGroupLib::getReadGroupID(rgi);
    if(rg.isBackbone() && !rg.isRail()){
      refstrainid=rg.getStrainID();
      break;
    }
  }
  // we will look at all GBF features (empty == all allowed)
  std::vector<multitag_t::mte_id_t> allowedfeatures;
  allowedfeatures.push_back(Read::REA_tagentry_idSOFACDS);
  std::vector<multitag_t::mte_id_t> forbiddenfeatures;
  for(auto & cle : clist){
    cle.sortConsensusTags();

    // copy first all Genbank features of all reads into one big vector
    //  (makes the search loop afterwards easier)

    std::list<gbfsummary_t> allGBfeatures;
    cle.getGBFSummary(allGBfeatures,allowedfeatures,forbiddenfeatures,true, false);

    // continue with next contig if there are no features in here
    if(allGBfeatures.empty()) continue;

    std::string consseq;
    std::vector<base_quality_t> consqual;
    std::vector<std::string> strain_consseq;
    std::vector< std::vector<base_quality_t> > strain_consqual;
    strainid2name_t strainnames_in_contig;

    makeAllStrainGenomes(cle,
			 0,
			 consseq,
			 consqual,
			 strain_consseq,
			 strain_consqual,
			 strainnames_in_contig,
			 false);

    if(refstrainid!=100){
      std::string featuredna;
      std::string featureprot;

      bool hascodontablewarning=false;
      for(auto & feature : allGBfeatures){
	// unsigned, always true: if(feature.cfrom<0 || feature.cto<0
	if(feature.cfrom >= consseq.size() || feature.cto >= consseq.size()){
	  CEBUGNPQ("Out of bounds, next." << endl);
	  continue;
	}

	if(codontableused==0){
	  codontableused=feature.translationtable;

	  std::vector<std::string> dummy1;
	  std::vector<bool> dummy2;
	  std::vector<char> aaresult;

	  char acgt[]="ACGT";
	  double unused=-1;
	  std::string codon("xxx");
	  for(int b1=0; b1<4; ++b1){
	    codon[0]=acgt[b1];
	    for(int b2=0; b2<4; ++b2){
	      codon[1]=acgt[b2];
	      for(int b3=0; b3<4; ++b3){
		codon[2]=acgt[b3];
		auto codonindex=dptools::getIndexOfCodon(codon);
		codoninfo[codonindex].codon=codon;
		codoninfo[codonindex].index=codonindex;

		dptools::codon2AminoAcids(codontableused,
					  codon[0],
					  codon[1],
					  codon[2],
					  dummy1,
					  aaresult,
					  dummy2);
		codoninfo[codonindex].aa=aaresult[0];
	      }
	    }
	  }
	}else{
	  if(codontableused!=feature.translationtable){
	    if(!hascodontablewarning){
	      hascodontablewarning=true;
	      cout << "WARNING! Contig " << cle.getContigName() << ": first features in annotation used a codon table of " << static_cast<uint16>(codontableused) << ", but now the feature " << " defines a codon table of " << static_cast<uint16>(feature.translationtable) << ". The codon frequency analysis will be completely wrong for this contig.\n";
	    }
	  }
	}

	// look how feature is translated
	if(feature.direction>0) {
	  CEBUGNPQ("#x1\n");
	  dptools::dnaToProtein(strain_consseq[refstrainid],
				featureprot,
				featuredna,
				feature.cfrom,
				feature.cto,
				feature.direction,
				feature.translationtable,
				feature.codonstart,
				true);
	}else{
	  CEBUGNPQ("#x2\n");
	  dptools::dnaToProtein(strain_consseq[refstrainid],
				featureprot,
				featuredna,
				feature.cto,
				feature.cfrom,
				feature.direction,
				feature.translationtable,
				feature.codonstart,
				true);
	}
	CEBUGNPQ("Protein for ref " << refstrainid << ": " << featureprot << endl);
	CEBUGNPQ("DNA for ref " << refstrainid << ": " << featuredna << endl);

	for(uint32 protpos=0; protpos<featureprot.size(); ++protpos){
	  auto codonindex=dptools::getIndexOfCodon(featuredna[3*protpos],featuredna[3*protpos+1],featuredna[3*protpos+2]);
	  if(codonindex<64){
	    codoninfo[codonindex].count+=1;
	  }
	}
      }
    }
  }

  // use non-parallel sort, this is a tiny thing
  mstd::ssort(codoninfo, codoninfo_t::sortByAAUpCountDown);
  {
    auto ciS=codoninfo.begin();
    while(ciS!=codoninfo.end()){
      auto ciE=ciS;
      uint64 tcount=0;
      for(; ciE!=codoninfo.end() && ciS->aa == ciE->aa; ++ciE) {tcount+=ciE->count;}
      if(tcount>0){
	for(auto ciI=ciS; ciI!=ciE; ++ciI){
	  ciI->frequency=static_cast<double>(100)/tcount*ciI->count;
	}
      }
      ciS=ciE;
    }
  }
  cout << "Codon table:\n";
  for(auto & cie : codoninfo){
    cout << cie << endl;
  }
  // use non-parallel sort, this is a tiny thing
  mstd::ssort(codoninfo, codoninfo_t::sortByIndexUp);

  // we will look at all GBF features (empty == all allowed)
  allowedfeatures.clear();
  //allowedfeatures.push_back("Fgen");
  allowedfeatures.push_back(Read::REA_tagentry_idSOFACDS);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAexon);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAintron);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAmRNA);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAtranscript);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAprimary_transcript);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFArRNA);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAscRNA);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAsnRNA);
  allowedfeatures.push_back(Read::REA_tagentry_idSOFAtRNA);

  // but not at Fsrc and Fgen
  forbiddenfeatures.clear();
  //forbiddenfeatures.push_back("Fsrc");
  //forbiddenfeatures.push_back("Fgen");
  //forbiddenfeatures.push_back("Fcon");
  //forbiddenfeatures.push_back("Fmod");
  //forbiddenfeatures.push_back("Fold");
  //forbiddenfeatures.push_back("Frpr");
  //forbiddenfeatures.push_back("Frpu");
  //forbiddenfeatures.push_back("F???");
  //forbiddenfeatures.push_back("Fvar");

  // now, go through each contig and analyse the effect of SNPs on feature
  //  sequences

  for(auto & cle : clist){
    cle.sortConsensusTags();

    // copy first all Genbank features of all reads into one big vector
    //  (makes the search loop afterwards easier)

    std::list<gbfsummary_t> allGBfeatures;
    // BaCh 11.09.2013: don't take consensus features into analysis
    // MIRA tags showing up (don't want to filter them atm) and gene /CDS features
    //  in consensus are not set by MIRA anyway. Therefore: save work, implement
    //  only if need arises (i.e. probably never)
    cle.getGBFSummary(allGBfeatures,allowedfeatures,forbiddenfeatures,true, false);

    // continue with next contig if there are no features in here
    if(allGBfeatures.empty()) continue;

    CEBUGNPQ("Preparing statistics\n");

    std::vector<proteinchangesummary_t> allproteinchanges;
    {
      static const proteinchangesummary_t pcs = {
	false,
	false,
	false,
	false,
	false,
	false,
	false,
	0,
	0,0,0,0,
	0,0,0,0,
	0,"",
	""
      };
      allproteinchanges.resize(allGBfeatures.size(),pcs);
    }

    // idee: nicht nur einen Consensus generieren, sondern für jeden
    //       strain einen eigenen. Dann kann man einfach strainweise
    //       die auswirkungen bei Proteinen abschätzen
    //  contig.C umändern.
    //  dann nicht abgedeckte stellen in einem strain (N's) jeweils durch
    //    IUPAC consensus anderer strains auffüllen

    std::string consseq;
    std::vector<base_quality_t> consqual;
    std::vector<std::string> strain_consseq;
    std::vector< std::vector<base_quality_t> > strain_consqual;
    strainid2name_t strainnames_in_contig;

    makeAllStrainGenomes(cle,
			 0,
			 consseq,
			 consqual,
			 strain_consseq,
			 strain_consqual,
			 strainnames_in_contig,
			 false);

    CEBUGNPQ("Going through SNPs\n");

    // make a quick lookup vector for which position is intergenic
    std::vector<uint8> intergeniclookup(consseq.size(),0);
    {
      //std::list<gbfsummary_t>::const_iterator F=allGBfeatures.begin();
      for(auto & f : allGBfeatures){
	if(f.identifier == "Figr"){
	  for(uint32 iglpos=f.cfrom; iglpos <= f.cto; iglpos++){
	    if(iglpos<intergeniclookup.size()) intergeniclookup[iglpos]=1;
	  }
	}
      }
    }

    // Now, go through each GB features and analyse the impact
    //  of SNPs to it
    // Also calculate the "real" dna and protein sequence for every
    //  strain on the fly

    std::vector< std::vector<std::string> > strain_featuredna;
    std::vector< std::vector<std::string> > strain_featureprot;
    strain_featuredna.resize(allGBfeatures.size());
    strain_featureprot.resize(allGBfeatures.size());
    for(uint32 i=0; i<allGBfeatures.size(); ++i){
      strain_featuredna[i].resize(strain_consseq.size());
      strain_featureprot[i].resize(strain_consseq.size());
    }

    auto fI=allGBfeatures.cbegin();
    auto pI=allproteinchanges.begin();
    for(uint32 featurecount=0; fI!=allGBfeatures.end(); ++fI, ++pI, ++featurecount){
      if(cle.getConsensusTags().empty()) continue;

      //CEBUGNPQ("We're looking at " << fI->identifier << "\t" << fI->locustag << '\t' << fI->cfrom << '\t' << fI->cto << endl);
      CEBUGNPQ("We're looking at " << *fI);
      // unsigned, always true: if(fI->cfrom<0 || fI->cto<0
      if(fI->cfrom >= consseq.size() || fI->cto >= consseq.size()){
	CEBUGNPQ("Out of bounds, next." << endl);
	pI->coveragestatus="out of contig";
	continue;
      }

      // featuredir is also the increment
      int8 featuredir=fI->direction;;
      CEBUGNPQ("Increment/featuredir is: " << static_cast<int16>(featuredir) << endl);

      pI->coveragestatus="ok";
      pI->firstcodonisstart=true;

      // collect information on the DNA and protein sequences
      for(uint32 strainid=0; strainid<strain_consseq.size(); ++strainid){
	if(strain_consseq[strainid].size()){

	  // check whether we have non covered bases in this strain
	  bool missingcompletely=false;
	  {
	    uint32 basesnocov=0;

	    const char * dnap=strain_consseq[strainid].c_str();
	    dnap+=fI->cfrom;
	    for(uint32 ii=fI->cfrom; ii<= fI->cto; dnap++, ii++){
	      if(toupper(*dnap)=='X') basesnocov++;
	    }

	    CEBUGNPQ("basesnocov: " << basesnocov << '\n');
	    if(basesnocov >= (fI->cto - fI->cfrom)){
	      pI->coveragestatus="completely missing?";
	      missingcompletely=true;
	    }else if(basesnocov > 0
		     && pI->coveragestatus =="ok"){
	      pI->coveragestatus="partly missing?";
	    }
	  }

	  // look how feature is translated
	  if(featuredir>0) {
	    CEBUGNPQ("#x1\n");
	    uint32 translateto=fI->cto;
	    if(fI->mustbetranslated) translateto=fI->cfrom;
	    dptools::dnaToProtein(strain_consseq[strainid],
				  strain_featureprot[featurecount][strainid],
				  strain_featuredna[featurecount][strainid],
				  fI->cfrom,
				  translateto,
				  featuredir,
				  fI->translationtable,
				  fI->codonstart,
				  fI->mustbetranslated);
	  }else{
	    CEBUGNPQ("#x2\n");
	    uint32 translateto=fI->cfrom;
	    if(fI->mustbetranslated) translateto=fI->cto;
	    dptools::dnaToProtein(strain_consseq[strainid],
				  strain_featureprot[featurecount][strainid],
				  strain_featuredna[featurecount][strainid],
				  fI->cto,
				  translateto,
				  featuredir,
				  fI->translationtable,
				  fI->codonstart,
				  fI->mustbetranslated);
	  }
	  CEBUGNPQ("Protein for strain " << strainid << ": " << strain_featureprot[featurecount][strainid] << endl);
	  CEBUGNPQ("DNA for strain " << strainid << ": " << strain_featuredna[featurecount][strainid] << endl);

	  // clear "protein" of features that are not to be translated
	  //  or where the coverage is completely missing
	  if(!fI->mustbetranslated
	     || missingcompletely){
	    strain_featureprot[featurecount][strainid].clear();
	    pI->firstcodonisstart=false;
	  }else if(strain_featuredna[featurecount][strainid].size() >= 3
		   && pI->firstcodonisstart){
	    pI->firstcodonisstart=dptools::isCodonStart(fI->translationtable,
						       strain_featuredna[featurecount][strainid][0],
						       strain_featuredna[featurecount][strainid][1],
						       strain_featuredna[featurecount][strainid][2]);
	    CEBUGNPQ("firstcodonisstart: " << pI->firstcodonisstart << '\n');
	  }
	}
      }



      // go through each SNP consensustag
      auto CTact=cle.getConsensusTags().cbegin();
      auto CTlast=cle.getConsensusTags().cend()-1;

      if(featuredir<0) {
	std::swap(CTact,CTlast);
      }

      bool lastconstagreached=false;
      for(; !lastconstagreached; CTact+=featuredir){
	if(CTact==CTlast) lastconstagreached=true;

	if(CTact->identifier == Contig::CON_tagentry_idSAOc
	   || CTact->identifier == Contig::CON_tagentry_idSIOc
	   || CTact->identifier == Contig::CON_tagentry_idSROc) {

	  bool mustconsider=false;
	  int32 reldistance=0;
	  //int8 updownlocation=0;
	  // look where the SNP is located: in the feature (location =0),
	  //  upstream (-1) or downstream (1)

	  //CEBUGNPQ("Testing F " << fI->cfrom << " " << fI->cto << "\t");
	  //CEBUGNPQ(CTact->getIdentifierStr() << " " << CTact->from << " " << CTact->to << endl);

	  // In dubio pro reo. To counter
	  //   - wrong annotation in initial genome or
	  //   - SNPs in mutant strain
	  // we look at the coordinates as dictated by the protein lengths
	  //  and decide upon these whether a SNP is inlocus,
	  //  up- or downstream


	  int32 realfcfrom=fI->cfrom;
	  int32 realfcto=fI->cto;
	  {
	    // first, get max protein length

	    int32 newgenelen=0;
	    {
	      size_t maxprotlen=0;
	      for(uint32 strainid=0; strainid<strain_consseq.size(); strainid++){
		maxprotlen=std::max(maxprotlen,strain_featureprot[featurecount][strainid].size());
	      }
	      newgenelen=static_cast<int32>(maxprotlen*3);
	    }

	    if(newgenelen > static_cast<int32>(fI->cto - fI->cfrom)){
	      CEBUG("Must adapt gene length to: " << newgenelen << '\n');
	      if(featuredir>0){
		realfcto=fI->cfrom+newgenelen;
	      }else{
		realfcfrom=fI->cto-newgenelen;
	      }
	    }
	  }

	  bool taginintergenic=false;
	  // always >=0: CTact->from >=0
	  if(CTact->from < consseq.size()
	     && intergeniclookup[CTact->from]){
	    taginintergenic=true;
	  }

	  if(static_cast<int32>(CTact->from) >= realfcfrom
	     && static_cast<int32>(CTact->from) <= realfcto){
	    mustconsider=true;
	    if(featuredir>0){
	      reldistance=CTact->from - (realfcfrom);
	    }else{
	      reldistance=(realfcto)-(CTact->from);
	    }
	  }

	  if(mustconsider){
	    CEBUGNPQ("Possible candidate for: " << *fI << endl);
	    CEBUGNPQ("CTact: " << *CTact << endl);

	    auto * constouse = &consseq;
	    auto * qualtouse = &consqual;
	    if(fI->strainid >= 0
	       && fI->strainid < strain_consseq.size()){
	      constouse=&strain_consseq[fI->strainid];
	      qualtouse=&strain_consqual[fI->strainid];
	    }

	    char baseinfeaturestrain='?';
	    if(fI->pcrI != cle.getContigReads().end()){
	      baseinfeaturestrain=static_cast<char>(toupper(cle.getBaseInRead(CTact->from,fI->pcrI)));
	      if(featuredir<0) baseinfeaturestrain=dptools::getComplementIUPACBase(baseinfeaturestrain);
	    }else{
	      baseinfeaturestrain=(*constouse)[CTact->from];
	    }

	    CEBUGNPQ("  baseinfeaturestrain: " << baseinfeaturestrain << "\n");

	    // compile some information on the position in the strain
	    //  containing the feature
	    std::string codoninfeature;
	    std::vector<std::string> codonvariantsinfeature;
	    std::vector<char> aainfeature;
	    std::vector<bool> isstartinfeature;
	    int32 aanumberinfeature=0;
	    int8  posinaainfeature=0;

	    if(fI->mustbetranslated){
	      if(featuredir>0) {
		dptools::infoOnAAatDNAPos(*constouse,
					  CTact->from,
					  realfcfrom,
					  featuredir,
					  fI->translationtable,
					  fI->codonstart,
					  codoninfeature,
					  codonvariantsinfeature,
					  aainfeature,
					  isstartinfeature,
					  aanumberinfeature,
					  posinaainfeature);
	      }else{
		dptools::infoOnAAatDNAPos(*constouse,
					  CTact->from,
					  realfcto,
					  featuredir,
					  fI->translationtable,
					  fI->codonstart,
					  codoninfeature,
					  codonvariantsinfeature,
					  aainfeature,
					  isstartinfeature,
					  aanumberinfeature,
					  posinaainfeature);
	      }

	      //BUGIFTHROW(aainfeature.size()==0,"aainfeature.size()==0?");
	      //BUGIFTHROW(isstartinfeature.size()==0,"isstartinfeature.size()==0?");
	      if(aainfeature.empty()) aainfeature.push_back('?');
	      if(isstartinfeature.empty()) isstartinfeature.push_back(false);
	    } else {
	      aainfeature.push_back('/');
	      isstartinfeature.push_back(false);
	    }

	    CEBUGNPQ("codoninfeature        : " << codoninfeature << endl);
	    CEBUGNPQ("codonvariantsinfeature: " << codonvariantsinfeature << endl);
	    CEBUGNPQ("aainfeature     : " << aainfeature << endl);
	    CEBUGNPQ("isstartinfeature: " << isstartinfeature << endl);

	    double reffreq=-1.0;
	    if(codonvariantsinfeature.size()==1){
	      reffreq=codoninfo[dptools::getIndexOfCodon(codonvariantsinfeature[0])].frequency;
	    }

	    CEBUGNPQ("\nLooking for strains ...");

	    for(auto & snic : strainnames_in_contig){
	      int32 actstrainid=snic.first;
	      CEBUGNPQ("Strainid: " << actstrainid << " ... ");
	      // do not look in own strain of feature
	      if(actstrainid==fI->strainid) continue;
	      // do not look in empty strains when inlocus
	      if(strain_featuredna[featurecount][actstrainid].size()==0) continue;

	      CEBUGNPQ("CHECK!\n");

	      char baseinstrain='?';
	      if(featuredir>0){
		baseinstrain=static_cast<char>(toupper(strain_consseq[actstrainid][CTact->from]));
	      }else{
		baseinstrain=static_cast<char>(toupper(dptools::getComplementIUPACBase(strain_consseq[actstrainid][CTact->from])));
	      }
	      CEBUGNPQ("baseinfeaturestrain : " << baseinfeaturestrain<<endl);
	      CEBUGNPQ("baseinstrain " << actstrainid << ": " << baseinstrain<<endl);
	      if(baseinstrain==baseinfeaturestrain) continue;

	      CEBUGNPQ("Ok, base in strain different from feature.\n");

	      pI->isaffected=true;

	      bool makesframeshift=false;
	      std::string dnaeffect="basechange";
	      if(baseinfeaturestrain=='*' || baseinstrain=='*'){
		makesframeshift=true;
		if(baseinstrain=='*'){
		  dnaeffect="deletion";
		}else{
		  dnaeffect="insertion";
		}
	      }

	      CEBUGNPQ("Makes frame shift: " << makesframeshift << '\n');
	      CEBUGNPQ("Effect on DNA: " << dnaeffect << '\n');

	      // compile some information on the position in the strain
	      //  we are looking at
	      std::string codoninstrain;
	      std::vector<std::string> codonvariantsinstrain;
	      std::vector<std::string> codonvariantsinstrainnewonly;

	      // aainstrain is original. e.g.: "TI"
	      std::vector<char> aainstrain;

	      // aainstrainnewonly is the above MINUS amino acids which were also in feature strain
	      // E.g.: feature has "T" and aainstrain has "TI", then aainstrainnewonly has "I"
	      std::vector<char> aainstrainnewonly;  // only "new" aminoacids
	      std::vector<bool> isstartinstrain;
	      int32 aanumberinstrain=0;
	      int8  posinaainstrain=0;
	      if(fI->mustbetranslated){
		CEBUGNPQ("must be translated.\n");
		if(featuredir>0) {
		  dptools::infoOnAAatDNAPos(strain_consseq[actstrainid],
					    CTact->from,
					    realfcfrom,
					    featuredir,
					    fI->translationtable,
					    fI->codonstart,
					    codoninstrain,
					    codonvariantsinstrain,
					    aainstrain,
					    isstartinstrain,
					    aanumberinstrain,
					    posinaainstrain);
		}else{
		  dptools::infoOnAAatDNAPos(strain_consseq[actstrainid],
					    CTact->from,
					    realfcto,
					    featuredir,
					    fI->translationtable,
					    fI->codonstart,
					    codoninstrain,
					    codonvariantsinstrain,
					    aainstrain,
					    isstartinstrain,
					    aanumberinstrain,
					    posinaainstrain);
		}
		//BUGIFTHROW(aainstrain.size()==0,"aainstrain.size()==0?");
		//BUGIFTHROW(isstartinstrain.size()==0,"isstartinstrain.size()==0?");
		if(aainstrain.empty()) aainstrain.push_back('?');
		if(isstartinstrain.empty()) isstartinstrain.push_back(false);
	      }else{
		CEBUGNPQ("not translated.\n");
		aainstrain.push_back('/');
		isstartinstrain.push_back(false);
	      }

	      // create aainstrainnewonly
	      {
		for(auto & aais : aainstrain){
		  bool found=false;
		  for(auto & aaif : aainfeature){
		    if(aaif == aais) {
		      found=true;
		      break;
		    }
		  }
		  if(!found) {
		    aainstrainnewonly.push_back(aais);
		  }
		}
	      }

	      CEBUGNPQ("\naainstrain     : " << aainstrain << '\n');
	      CEBUGNPQ("aainstrainnewonly: " << aainstrainnewonly << '\n');

	      // create aainstrainnewonly
	      {
		for(auto & cvis : codonvariantsinstrain){
		  bool found=false;
		  for(auto & cvif : codonvariantsinfeature){
		    if(cvif == cvis) {
		      found=true;
		      break;
		    }
		  }
		  if(!found) {
		    codonvariantsinstrainnewonly.push_back(cvis);
		  }
		}
	      }
	      CEBUGNPQ("cv strain       : " << codonvariantsinstrain << '\n');
	      CEBUGNPQ("cv strainnewonly: " << codonvariantsinstrainnewonly << '\n');

	      double strainfreq=-1.0;
	      if(codonvariantsinstrainnewonly.empty() && codonvariantsinstrain.size()==1){
		strainfreq=codoninfo[dptools::getIndexOfCodon(codonvariantsinstrain[0])].frequency;
	      }else if(codonvariantsinstrainnewonly.size()==1){
		strainfreq=codoninfo[dptools::getIndexOfCodon(codonvariantsinstrainnewonly[0])].frequency;
	      }

	      // untranslated SNPs are SNPs that occur within the
	      //  feature, but where the protein still ended earlier
	      std::string effectonprot="none";
	      if(fI->mustbetranslated){
		CEBUGNPQ("#0\n");
		bool untranslated=false;
		if(aanumberinstrain>=static_cast<int32>(strain_featureprot[featurecount][actstrainid].size())) {
		  untranslated=true;

		  CEBUGNPQ("Untranslated!\tAAnumberinstrain: " << aanumberinstrain << "\tAAnumberinfeatureprot: " << strain_featureprot[featurecount][actstrainid].size() << endl);

		}

		// frameshift, aa changes and silents
		if(makesframeshift) {
		  CEBUGNPQ("#1\n");
		  effectonprot="frameshift";
		  // corrective measures for "untranslated"
		  if(untranslated){
		    CEBUGNPQ("#2a\n");
		    if(baseinstrain=='*') {
		      CEBUGNPQ("#3a\n");
		      pI->deletionuntranslated+=1;
		    }else{
		      CEBUGNPQ("#4a\n");
		      pI->insertionuntranslated+=1;
		    }
		  }else{
		    CEBUGNPQ("#2b\n");
		    if(baseinstrain=='*') {
		      CEBUGNPQ("#3b\n");
		      pI->deletioninlocus+=1;
		    }else{
		      CEBUGNPQ("#4b\n");
		      pI->insertioninlocus+=1;
		    }
		  }
		}else if(aainstrainnewonly.empty()){
		  // now new aa created? that's a silent thing
		  effectonprot="none";
		  if(untranslated){
		    CEBUGNPQ("#5a\n");
		    pI->silentuntranslated+=1;
		  }else{
		    CEBUGNPQ("#5b\n");
		    pI->silentinlocus+=1;
		  }
		}else{
		  CEBUGNPQ("#6\n");
		  effectonprot.clear();
		  for(auto & aais : aainstrainnewonly){
		    if(!effectonprot.empty()) effectonprot+=" / ";
		    if(aainfeature[0]!=aais){
		      CEBUGNPQ("#7\n");
		      effectonprot+="aa change";
		      if(untranslated){
			CEBUGNPQ("#7a\n");
			pI->aachangeuntranslated+=1;
		      }else{
			CEBUGNPQ("#7b\n");
			pI->aachangeinlocus+=1;
		      }
		    }else{
		      CEBUGNPQ("#8\n");
		      effectonprot+="none";
		      if(untranslated){
			CEBUGNPQ("#9\n");
			pI->silentuntranslated+=1;
		      }else{
			CEBUGNPQ("#10\n");
			pI->silentinlocus+=1;
		      }
		    }
		  }
		}
		if(aainfeature[0]!='*'){
		  CEBUGNPQ("#11\n");
		  for(auto & aais : aainstrainnewonly){
		    if(aais=='*') {
		      effectonprot+=" / creates new stop codon";
		      pI->prematurestop=true;
		      break;
		    }
		  }
		}

		// destroy/change Start/Stop
		if(aanumberinfeature==0 && isstartinfeature[0]) {
		  CEBUGNPQ("#12\n");
		  if(!isstartinstrain[0]){
		    CEBUGNPQ("#13\n");
		    effectonprot+=" / destroys start codon";
		    pI->destroyedstart=true;
		  }else{
		    CEBUGNPQ("#14\n");
		    effectonprot+=" / switches to another start codon";
		    pI->changedstart=true;
		  }
		}
		if(aainfeature[0]=='*'){
		  CEBUGNPQ("#15\n");
		  for(auto & aais : aainstrainnewonly){
		    if(aais!='*'){
		      CEBUGNPQ("#16\n");
		      effectonprot+=" / destroys stop codon";
		      pI->destroyedstop=true;
		    }else{
		      CEBUGNPQ("#17\n");
		      effectonprot+=" / switches to another stop codon";
		      pI->changedstop=true;
		    }
		  }
		}
		if(untranslated){
		  CEBUGNPQ("#18\n");
		  effectonprot+=" :: untranslated, protein ended earlier";
		}
	      }else{
		CEBUGNPQ("#19\n");
		effectonprot="unknown";

		if(taginintergenic){
		  CEBUGNPQ("#20\n");
		  pI->mutinintergenic+=1;
		}
	      }


	      // for the feature analysis, only print out inlocus
	      int32 adjposfroml=CTact->from;
	      int32 adjposfromr=adjposfroml;
	      if(fI->pcrI != cle.getContigReads().end()){
		const Read & featureread=*(fI->pcrI);
		int32 rreadposfrom=cle.getRealReadPos(CTact->from,fI->pcrI);
		if(rreadposfrom>=0){
		  CEBUGNPQ("#c1" << endl);
		  adjposfroml=featureread.getLowerNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
		  CEBUGNPQ("#c2" << endl);
		  adjposfromr=featureread.getUpperNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
		}
	      }
	      CEBUGNPQ("#21\n");
	      faout << cle.getContigName() << "\t";
	      if(fI->pcrI != cle.getContigReads().end()){
		faout << ReadGroupLib::getStrainOfStrainID(fI->pcrI->getStrainID()) << "\t";
	      }else{
		if(fI->strainid>=0 && fI->strainid<strain_consseq.size()){
		  faout << ReadGroupLib::getStrainOfStrainID(fI->strainid) << "\t";
		}else{
		  faout << "NA\t";
		}
	      }
	      CEBUGNPQ("#22\n");
	      faout << ReadGroupLib::getStrainOfStrainID(actstrainid) << "\t";
	      faout << CTact->from << "\t" << adjposfroml;
	      if(adjposfroml!=adjposfromr) faout << ':' << adjposfromr;
	      faout << '\t' << (1.0/cle.getContigLength()*adjposfroml*360);
	      faout << '\t' << fI->identifier;
	      faout << '\t' << fI->locustag << '\t';
	      if(featuredir<0) {
		faout << "complement\t";
	      }else{
		faout << "forward\t";
	      }
	      faout << fI->gene << "\t";
	      faout.flush();

	      CEBUGNPQ("#23\n");

	      //hitother not used as this is pure biology
	      //faout << "hitsother?" << "\t";
	      //faout.flush();
	      faout << reldistance << "\t";
	      faout.flush();
	      faout << dnaeffect << "\t";
	      faout.flush();
	      faout << baseinfeaturestrain << "->"; // << baseinstrain
	      for(auto & tmpc : dptools::getNucleicAcidFromIUPAC(baseinstrain)) faout << tmpc;
	      faout << '\t' << codoninfeature << "\t" << codoninstrain << "\t";
	      faout.flush();

	      if(reffreq>=0) faout << reffreq;
	      faout << "\t";
	      if(strainfreq>=0) faout << strainfreq;
	      faout << "\t";

	      CEBUGNPQ("#24\n");
	      if(fI->pcrI != cle.getContigReads().end()){
		faout << static_cast<uint16>(cle.getQualityInRead(CTact->from,fI->pcrI)) << "\t";
	      }else{
		faout << static_cast<uint16>((*qualtouse)[CTact->from]) << "\t";
	      }
	      CEBUGNPQ("#25\n");
	      faout << static_cast<uint16>(consqual[CTact->from]) << "\t";

	      faout.flush();

	      for(auto & tmpf : aainfeature) faout << tmpf;
	      faout << "->";
	      for(auto & tmps : aainstrain) {
		bool found=false;
		for(auto & tmpsn : aainstrainnewonly) {
		  if(tmps==tmpsn){
		    found=true;
		    break;
		  }
		}
		if(found) {
		  faout << tmps;
		}else{
		  faout << static_cast<char>(tolower(tmps));
		}
	      }
	      faout << "\t";

	      CEBUGNPQ("#26\n");

	      faout << aanumberinfeature << "\t";
	      faout << static_cast<int16>(posinaainfeature) << "\t";
	      faout << aanumberinstrain << "\t";
	      faout << static_cast<int16>(posinaainstrain) << "\t";

	      faout.flush();

	      faout << effectonprot << "\t";

	      faout << fI->goprocess << "\t";
	      faout << fI->gofunction << "\t";

	      faout << fI->product << "\t";
	      faout << fI->function << "\t";
	      faout << fI->note << "\t";

	      faout << endl;
	      CEBUGNPQ("#27\n");
	    }
	  }
	}
      }
    }

    // went through all features of contig
    // now print summary for each interesting feature

    cout << "Summary of changes:\n";

    fI=allGBfeatures.begin();
    pI=allproteinchanges.begin();
    for(uint32 featurecount=0; fI!=allGBfeatures.end(); ++fI, ++pI, ++featurecount){
      if(AnnotationMappings::isMIRAEntry(fI->identifier)) continue;

      CEBUGNPQ("Summary for " << fI->locustag << '\t' << fI->cfrom << '\t' << fI->cto << endl);

      fcout << cle.getContigName() << '\t' <<  fI->identifier << '\t' << fI->locustag << '\t' << fI->gene;
      if(fI->strainid>=0 && fI->strainid<strain_consseq.size()){
	fcout << "\t" << ReadGroupLib::getStrainOfStrainID(fI->strainid);
      }else{
	fcout << "\tNA";
      }

      for(uint32 strainid=0; strainid<strain_consseq.size(); strainid++){
	if(strainid!=static_cast<uint32>(fI->strainid)){
	  fcout << "\t" << ReadGroupLib::getStrainOfStrainID(strainid);
	}
      }

      uint64 dnahash=0;
      if(fI->strainid>=0 && fI->strainid<strain_featureprot.size()){
	std::string reference=strain_featureprot[featurecount][fI->strainid];
	adjustCaseOfSequences(reference, strain_featureprot[featurecount]);

	fcout << "\t" << reference;

	for(uint32 strainid=0; strainid<strain_consseq.size(); strainid++){
	  if(strainid!=static_cast<uint32>(fI->strainid)){
	    fcout << "\t" << strain_featureprot[featurecount][strainid];
	  }
	}

	reference=strain_featuredna[featurecount][fI->strainid];
	dnahash=str_hash(reference);
	for(uint32 strainid=0; strainid<strain_consseq.size(); strainid++){
	  if(strainid!=static_cast<uint32>(fI->strainid)){
	    uint64 tmphash=str_hash(strain_featuredna[featurecount][strainid]);
	    if(tmphash!=dnahash) dnahash=tmphash;
	  }
	}
	dnahash&=0x3fffffff; // only 30 bits please

	adjustCaseOfSequences(reference, strain_featuredna[featurecount]);

	fcout << "\t" << reference;

	for(uint32 strainid=0; strainid<strain_consseq.size(); strainid++){
	  if(strainid!=static_cast<uint32>(fI->strainid)){
	    fcout << "\t" << strain_featuredna[featurecount][strainid];
	  }
	}
      }
      fcout << endl;


      int32 adjposfroml=fI->cfrom;
      if(fI->pcrI != cle.getContigReads().end()){
	int32 rreadposfrom=cle.getRealReadPos(fI->cfrom,fI->pcrI);
	CEBUG("#d1" << endl);
	if(rreadposfrom>=0){
	  adjposfroml=fI->pcrI->getLowerNonGapAdjustmentPosOfReadPos(rreadposfrom)+1;
	  CEBUG("#d2" << endl);
	}
      }

      fsout << fI->locustag << '\t' << fI->gene
	    << '\t' << fI->identifier
	    << "\tz" << toOwnBase32(dnahash,30)
	    << '\t' << (1.0/cle.getContigLength()*adjposfroml*360) << '\t';

      // make a quick summary whether one should look at this protein
      {
	// interestlevel 0 = no, 1 = perhaps, 2 = yes
	uint32 interestlevel=0;
	if(
	  pI->changedstart
	  || pI->changedstop

	  || pI->silentinlocus > 0

	  || pI->silentuntranslated    > 0

	  ){
	  CEBUGNPQ("Interest 1\n");
	  interestlevel=1;
	}
	if(
	  pI->destroyedstart
	  || pI->destroyedstop
	  || pI->prematurestop

	  || pI->insertioninlocus > 0
	  || pI->deletioninlocus  > 0
	  || pI->aachangeinlocus  > 0

	  || pI->insertionuntranslated > 0           // in dubio pro reo
	  || pI->deletionuntranslated  > 0           // in dubio pro reo
	  || pI->aachangeuntranslated  > 0           // in dubio pro reo

	  || (pI->coveragestatus!="ok"
	      && pI->coveragestatus!="out of contig")
	  ){
	  CEBUGNPQ("Interest 2/std\n");
	  interestlevel=2;
	}

	// account for wrong annotation
	if(!pI->firstcodonisstart
	   &&	(pI->silentinlocus >0
		 || pI->insertionuntranslated > 0           // in dubio pro reo
		 || pI->deletionuntranslated  > 0           // in dubio pro reo
		 || pI->aachangeuntranslated  > 0           // in dubio pro reo
	     )){
	  CEBUGNPQ("Interest 2/account wrong anno\n");
	  interestlevel=2;
	}

	// intergenics are always interesting
	if(fI->identifier == "Figr"
	   && pI->mutinintergenic){
	  interestlevel=2;
	  CEBUGNPQ("Interest 2/intergenic\n");
	}


	if(interestlevel==0){
	  fsout << "no\t";
	}else if(interestlevel==1){
	  fsout << "perhaps\t";
	}else{
	  fsout << "yes\t";
	}
      }

      fsout << '\t';       // own filter

      fsout << pI->coveragestatus << '\t';

      if(fI->identifier == "FCDS"){
	if(pI->firstcodonisstart) {fsout << "yes\t";}else{fsout << "no\t";}
	if(pI->changedstart) {fsout << "yes\t";}else{fsout << "no\t";}
	if(pI->destroyedstart) {fsout << "yes\t";}else{fsout << "no\t";}
	if(pI->changedstop) {fsout << "yes\t";}else{fsout << "no\t";}
	if(pI->destroyedstop) {fsout << "yes\t";}else{fsout << "no\t";}
	if(pI->prematurestop) {fsout << "yes\t";}else{fsout << "no\t";}
      }else{
	fsout << "n/a\t"
	      << "n/a\t"
	      << "n/a\t"
	      << "n/a\t"
	      << "n/a\t"
	      << "n/a\t";
      }

      fsout << pI->mutinintergenic << "\t";

      fsout << pI->insertioninlocus << "\t";
      fsout << pI->deletioninlocus << "\t";
      fsout << pI->silentinlocus << "\t";
      fsout << pI->aachangeinlocus << "\t";

      fsout << pI->insertionuntranslated << "\t";
      fsout << pI->deletionuntranslated << "\t";
      fsout << pI->silentuntranslated << "\t";
      fsout << pI->aachangeuntranslated << "\t";

      fsout << fI->goprocess << "\t";
      fsout << fI->gofunction << "\t";

      fsout << fI->product << "\t";
      fsout << fI->function << "\t";
      fsout << fI->note;

      fsout << endl;


    }
  }

  cout << "Done with feature analysis." << endl;
}
//#define CEBUGNPQ(bla)



/*************************************************************************
 *
 * Lowercase all strings, then compare each mutant string to reference string
 *
 * If there's a difference, uppercase reference string at that pos
 *  as well as the given mutant string at that pos
 *
 *************************************************************************/

void assout::adjustCaseOfSequences(std::string & reference, std::vector<std::string> & mutants)
{
  size_t maxstring=reference.size();
  for(uint32 m=0; m<mutants.size(); m++) maxstring=std::max(maxstring,mutants[m].size());

  CEBUG("maxstring: " << maxstring << endl);

  for(size_t i=0; i<maxstring; i++){
    if(i<reference.size()) {
      reference[i]=static_cast<char>(tolower(reference[i]));
      bool ismutant=false;
      for(size_t m=0; m<mutants.size(); m++) {
	if(i<mutants[m].size()){
	  mutants[m][i]=static_cast<char>(tolower(mutants[m][i]));
	  CEBUG("r[i] / m[m][i] 1: " << reference[i] << ' ' << mutants[m][i] << '\n');

	  if(reference[i] != mutants[m][i]){
	    CEBUG("r[i] / m[m][i] 1.6: " << reference[i] << ' ' << mutants[m][i] << '\n');
	    mutants[m][i]=static_cast<char>(toupper(mutants[m][i]));
	  }
	  CEBUG("r[i] / m[m][i] 2: " << reference[i] << ' ' << mutants[m][i] << '\n');
	}else if(m!=0 || mutants[m].size()!=0){
	  // the above is for catching the "default" strain which might be
	  //  completely empty in assemblies where each read has a strain attached

	  // this mutant has ended, uppercase reference
	  ismutant=true;
	}
      }
      if(ismutant) reference[i]=static_cast<char>(toupper(reference[i]));
    }else{
      // reference string has already ended, uppercase mutants
      for(size_t m=0; m<mutants.size(); m++) {
	if(i<mutants[m].size()) mutants[m][i]=static_cast<char>(toupper(mutants[m][i]));
      }
    }
  }

  CEBUG("refstr: " << reference << '\n');
  for(uint32 m=0; m<mutants.size(); m++) {
    CEBUG("mutstr: " << mutants[m] << '\n');
  }

}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveAsFASTA(std::list<Contig> & clist, const std::string & filename, const std::string & paddedfilename, bool deleteoldfile)
{
  FUNCSTART("void saveAsFASTA(std::list<Contig> & clist, const std::string & filename, const std::string & paddedfilename)");

  std::ofstream fastaout;
  std::ofstream fastapaddedout;
  std::ofstream qualout;
  std::ofstream qualpaddedout;

  std::string qualname;
  std::string paddedqualname;
  if(filename.size()) {
    qualname=filename+".qual";
    paddedqualname=paddedfilename+".qual";
  }

  openFileForAppend(filename, fastaout, deleteoldfile);
  openFileForAppend(paddedfilename, fastapaddedout, deleteoldfile);
  openFileForAppend(qualname, qualout, deleteoldfile);
  openFileForAppend(paddedqualname, qualpaddedout, deleteoldfile);

  // Save contigs first, then singlets
  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	try{
	  Contig::setCoutType(Contig::AS_FASTAPADDED);
	  fastapaddedout << cle;
	  Contig::setCoutType(Contig::AS_FASTAPADDEDQUAL);
	  qualpaddedout << cle;
	  Contig::setCoutType(Contig::AS_FASTA);
	  fastaout << cle;
	  Contig::setCoutType(Contig::AS_FASTAQUAL);
	  qualout << cle;
	}
	catch (Notify n) {
	  cerr << "Error while dumping " << cle.getContigName() << ".\n";
	  n.handleError(THISFUNC);
	}
      }
    }
  }

  FUNCEND();
}

void assout::saveAsFASTA(Contig & con, const std::string & filename, const std::string & paddedfilename, bool deleteoldfile)
{
  FUNCSTART("void saveAsFASTA(Contig & con, const std::string & filename, const std::string & paddedfilename)");

  std::ofstream fastaout;
  std::ofstream fastapaddedout;
  std::ofstream qualout;
  std::ofstream qualpaddedout;

  std::string qualname;
  std::string paddedqualname;
  if(filename.size()) {
    qualname=filename+".qual";
    paddedqualname=paddedfilename+".qual";
  }

  openFileForAppend(filename, fastaout, deleteoldfile);
  openFileForAppend(paddedfilename, fastapaddedout, deleteoldfile);
  openFileForAppend(qualname, qualout, deleteoldfile);
  openFileForAppend(paddedqualname, qualpaddedout, deleteoldfile);

  try{
    Contig::setCoutType(Contig::AS_FASTAPADDED);
    fastapaddedout << con;
    Contig::setCoutType(Contig::AS_FASTAPADDEDQUAL);
    qualpaddedout << con;
    Contig::setCoutType(Contig::AS_FASTA);
    fastaout << con;
    Contig::setCoutType(Contig::AS_FASTAQUAL);
    qualout << con;
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveStrainsAsFASTAQ(std::list<Contig> & clist, const ReadPool & rp, const std::string & paddedfilename, bool asfastq, uint32 mincoverage, base_quality_t minqual, bool deleteoldfile, bool fillholesinstrain)
{
  FUNCSTART("void saveAsFASTAQ(std::list<Contig> & clist, const std::string & paddedfilename)");

  //cout << "Saving padded strain contigs to FASTA file: " << paddedfilename << "_<strainname>.fasta" << endl;
  //cout << "(qualities have .qual appended)" << endl;

  // Save contigs first, then singlets
  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	saveStrainsAsFASTAQ(cle,rp,paddedfilename,asfastq,mincoverage,minqual,deleteoldfile,fillholesinstrain);
	deleteoldfile=false;
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

void assout::saveStrainsAsFASTAQ(Contig & outcon, const ReadPool & rp, const std::string & paddedfilename, bool asfastq, uint32 mincoverage, base_quality_t minqual, bool deleteoldfile, bool fillholesinstrain)
{
  FUNCSTART("void assout::saveStrainsAsFASTAQ(Contig & outcon, const ReadPool & rp, const std::string & paddedfilename, bool asfastq, uint32 mincoverage, base_quality_t minqual, bool deleteoldfile, bool fillholesinstrain)");

  // make sure to dump "AllStrains" only if we have more than one strain
  int32 startid=-1;
  if(ReadGroupLib::getNumOfStrains()==1) startid=0;

  for(int32 strainid=startid; strainid<static_cast<int32>(ReadGroupLib::getNumOfStrains()); ++strainid){
    //if(outcon.getNumReadsPerStrain(strainid)==0) continue;
    std::string strainfilename;
    if(strainid>=0){
      strainfilename=ReadGroupLib::getStrainOfStrainID(strainid);
    }else{
      strainfilename="AllStrains";
    }
    if(outcon.getVerbose()) cout << "### calc " << strainfilename << " with strainid " << strainid << endl;
    try{
      if(asfastq){
	std::ofstream fastqpaddedout;
	openFileForAppend((paddedfilename+"_"+strainfilename+".padded.fastq"), fastqpaddedout, deleteoldfile);
	outcon.dumpStrainAsFASTQ(fastqpaddedout,
				 mincoverage,
				 minqual,
				 strainid,
				 true,
				 fillholesinstrain);

	fastqpaddedout.close();
	openFileForAppend((paddedfilename+"_"+strainfilename+".unpadded.fastq"), fastqpaddedout, deleteoldfile);
	outcon.dumpStrainAsFASTQ(fastqpaddedout,
				 mincoverage,
				 minqual,
				 strainid,
				 false,
				 fillholesinstrain);

	fastqpaddedout.close();
      }else{
	std::ofstream fastapaddedout;
	std::ofstream qualpaddedout;
	openFileForAppend((paddedfilename+"_"+strainfilename+".padded.fasta"), fastapaddedout, deleteoldfile);
	openFileForAppend((paddedfilename+"_"+strainfilename+".padded.fasta.qual"), qualpaddedout, deleteoldfile);

	outcon.dumpStrainAsFASTAQUAL(fastapaddedout,
				     qualpaddedout,
				     mincoverage,
				     minqual,
				     strainid,
				     true,
				     fillholesinstrain);

	fastapaddedout.close();
	qualpaddedout.close();

	openFileForAppend((paddedfilename+"_"+strainfilename+".unpadded.fasta"), fastapaddedout,deleteoldfile);
	openFileForAppend((paddedfilename+"_"+strainfilename+".unpadded.fasta.qual"), qualpaddedout, deleteoldfile);

	outcon.dumpStrainAsFASTAQUAL(fastapaddedout,
				     qualpaddedout,
				     mincoverage,
				     minqual,
				     strainid,
				     false,
				     fillholesinstrain);
	fastapaddedout.close();
	qualpaddedout.close();
      }
    }
    catch (Notify n) {
      cerr << "Error while dumping " << outcon.getContigName() << " for strain " << strainfilename << ".\n";
      n.handleError(THISFUNC);
    }
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveStrainsAsGBF(std::list<Contig> & clist, const ReadPool & rp, const std::string & gbfbasename, base_quality_t minqual, bool fillholesinstraingenomes, bool deleteoldfile)
{
  // Save contigs first, then singlets
  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	std::string consseq;
	std::vector<base_quality_t> consqual;
	std::vector<std::string> strain_consseq;
	std::vector< std::vector<base_quality_t> > strain_consqual;
	strainid2name_t strainnames_in_contig;

	makeAllStrainGenomes(cle,
			     minqual,
			     consseq,
			     consqual,
			     strain_consseq,
			     strain_consqual,
			     strainnames_in_contig,
			     fillholesinstraingenomes);

	if(!fillholesinstraingenomes){
	  consseq.clear();
	  consqual.clear();
	}

	// we will look at all GBF features (empty allowed)
	std::vector<multitag_t::mte_id_t> allowedfeatures;
	// but not at Fsrc and Fgen
	std::vector<multitag_t::mte_id_t> forbiddenfeatures;
	forbiddenfeatures.push_back(Read::REA_tagentry_idSOFAdatabank_entry);

	// copy first all Genbank features of all reads into one big vector
	//  (makes the search loop afterwards easier)

	std::list<gbfsummary_t> allGBfeatures;
	cle.getGBFSummary(allGBfeatures,
			  allowedfeatures,
			  forbiddenfeatures,
			  false,
			  true);

	//// dump contigs only for strains that are in this contig

	for(auto & snic : strainnames_in_contig){
	  int32 actstrainid=snic.first;
	  CEBUGNPQ("Strainid: " << actstrainid << "\tStrainname:" << snic.second << endl);

	  std::ofstream gbfout;
	  openFileForAppend((gbfbasename+"_"+ReadGroupLib::getStrainOfStrainID(actstrainid)+".gbf"),
			    gbfout,
			    deleteoldfile);

	  // beware! strain_consseq[actstrainid] will be different after call
	  dumpSequenceAsGBF_priv(cle.getContigName(),
				 snic.second,
				 strain_consseq[actstrainid],
				 consseq,
				 allGBfeatures,
				 gbfout);

	}
      }
    }
  }

}


/*************************************************************************
 *
 *  beware! UGLY! (but faster than first copying the whole std::string):
 *   "seq" will be different after call (lower case and eventually
 *   filled up with altseq if a character in seq is == 'n'
 *
 *************************************************************************/

//#define CEBUGF(bla)  {cout << bla; cout.flush();}

void assout::dumpSequenceAsGBF_priv(const std::string & seqname, const std::string & strainname, std::string & seq, const std::string & consseq, const std::list<gbfsummary_t> & allGBfeatures, std::ostream & fout)
{
  // position mapping of consensus: padded (in memory) vs unpadded
  // positions with a * get position of last non-* assigned
  uint32 actdepadpos=0;
  std::vector<uint32> depadposmap(seq.size(),0);
  for(uint32 seqi=0; seqi < seq.size();++seqi){
    depadposmap[seqi]=actdepadpos;
    if(seq[seqi]!='*') ++actdepadpos;
  }

  CEBUGF("depadposmap.size(): " << depadposmap.size() << endl);

  fout << "LOCUS       ";
  // dummy
  fout.width(16);
  fout << std::left << seqname << "  ";
  fout << actdepadpos << " bp    DNA     ENV    ";

  {
    char timestr[80];
    struct tm *timepoint;
    time_t t;

    time(&t);
    timepoint=localtime(&t);

    strftime(&timestr[0],60,"%d-%b-%Y",timepoint);

    fout << timestr << "\n";
  }

  fout << "SOURCE      " <<
    "DNA, strain: "<< strainname <<
    "\nACCESSION   " << seqname <<
    "\nCOMMENT     Generated by the MIRA assembler\n"
    "FEATURES             Location/Qualifiers\n"
    "     source          1.." << actdepadpos << endl <<
    "                     /strain=\"" << strainname << "\"\n";

  for(auto & gbkf : allGBfeatures){
    CEBUGF("###gbkf:\n" << gbkf);
    //CEBUGF("\n" << gbkf.locustag << "\t" << gbkf.cfrom << "\t" << gbkf.cto << "\t");
    if(gbkf.cfrom >= depadposmap.size() || gbkf.cto >= depadposmap.size()) continue;

    const std::string & gbfstring=AnnotationMappings::translateSOfeat2GenBankfeat(gbkf.identifier);
    if(gbfstring.empty()) continue;

    fout << "     ";
    fout.width(16);

    fout <<std::left << gbfstring;
    CEBUGF("0");
    if(gbkf.direction >=0) {
      fout << depadposmap[gbkf.cfrom]+1 << ".." << depadposmap[gbkf.cto]+1 << "\n";
    }else{
      fout << "complement(" << depadposmap[gbkf.cfrom]+1 << ".." << depadposmap[gbkf.cto]+1 << ")\n";
    }

    CEBUGF("1");

    dumpTextAsGBFValueLine_priv("/gene",gbkf.gene,fout);
    dumpTextAsGBFValueLine_priv("/locus_tag",gbkf.locustag,fout);
    dumpTextAsGBFValueLine_priv("/old_locus_tag",gbkf.oldlocustag,fout);
    dumpTextAsGBFValueLine_priv("/gene_synonym",gbkf.genesynonym,fout);
    dumpTextAsGBFValueLine_priv("/regulatory_class",gbkf.regulatoryclass,fout);
    dumpTextAsGBFValueLine_priv("/function",gbkf.function,fout);
    dumpTextAsGBFValueLine_priv("/EC_number",gbkf.ecnumber,fout);
    dumpTextAsGBFValueLine_priv("/product",gbkf.product,fout);
    dumpTextAsGBFValueLine_priv("/GO_process",gbkf.goprocess,fout);
    dumpTextAsGBFValueLine_priv("/GO_function",gbkf.gofunction,fout);
    dumpTextAsGBFValueLine_priv("/GO_component",gbkf.gocomponent,fout);
    dumpTextAsGBFValueLine_priv("/inference",gbkf.inference,fout);
    dumpTextAsGBFValueLine_priv("/note",gbkf.note,fout);

    if(gbfstring == "CDS" || gbfstring == "gene"){
      if(gbkf.translationtable > 0){
	fout << "                     /transl_table=" << static_cast<uint16>(gbkf.translationtable) << '\n';
      }
      if(gbkf.codonstart > 0 && gbkf.codonstart < 4){
	fout << "                     /codon_start=" << static_cast<uint16>(gbkf.codonstart) << '\n';
      }
    }

    CEBUGF("5");
    if(gbkf.mustbetranslated){
      // dump translation
      std::string featureprot;
      std::string featuredna;
      if(gbkf.direction >=0) {
	dptools::dnaToProtein(seq,
			      featureprot,
			      featuredna,
			      gbkf.cfrom,
			      gbkf.cfrom,
			      1,
			      gbkf.translationtable,
			      gbkf.codonstart,
			      gbkf.mustbetranslated);
      }else{
	dptools::dnaToProtein(seq,
			      featureprot,
			      featuredna,
			      gbkf.cto,
			      gbkf.cto,
			      -1,
			      gbkf.translationtable,
			      gbkf.codonstart,
			      gbkf.mustbetranslated);
      }
      CEBUGF("6");
      if(!featureprot.empty()){
	// dnaToProtein includes the * stop codon, remove that before output
	featureprot.resize(featureprot.size()-1);
	dumpTextAsGBFValueLine_priv("/translation",featureprot,fout);
      }
      CEBUGF("7");
    }
  }

  fout << "ORIGIN\n";

  uint8 posinline=0;
  uint32 depadpos=1;
  bool isnewline=true;
  bool addspace=false;
  for(uint32 seqi=0; seqi<seq.size(); ++seqi){
    if(isnewline) {
      isnewline=false;
      fout.width(9);
      fout << std::right << depadpos;
    }
    if(posinline%10 == 0){
      addspace=true;
    }
    if(seq[seqi]!='*'){
      if(addspace) {
	addspace=false;
	fout << " ";
      }
      fout << seq[seqi];
      ++depadpos;
      ++posinline;
      if(posinline==60) {
	fout << "\n";
	posinline=0;
	isnewline=true;
      }
    }
  }
  if(!isnewline){
    fout << "\n";
  }
  fout << "//\n";

}

//#define CEBUGF(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void assout::dumpTextAsGBFValueLine_priv(const std::string & key, const std::string & value, std::ostream & fout)
{
  if(value.empty()) return;

  fout << "                     " << key << "=\"";
  fout.flush();

  std::string::size_type posinline=21+key.size()+2;

  std::string::size_type fpos=0;
  std::string::size_type tpos=0;

  bool firstwordinline=true;

  int32 runawaystop=100000;

  // print out "words" until the end of the string
  while(fpos<value.size()){

    if(--runawaystop==0) {
      cout << "Abort: runawaystop\n";
      exit(1);
    }

    bool splitit=false;
    bool onnextline=false;

    // find next space (or end of string)
    tpos=value.find(" ",fpos);
    if(tpos==std::string::npos){
      tpos=value.size();
    }

    auto spaceneeded=tpos-fpos;
    if(!firstwordinline) spaceneeded++;

    // see whether there is enough place remaining in the line
    //  (careful for space when not firstword??)
    if(posinline+spaceneeded > 79) onnextline=true;

    // max 58 chars in one line, if greater, we´ll have to split anyways
    if(spaceneeded > 58) splitit=true;

    if(splitit) {
      for(; fpos<value.size() && posinline<79; fpos++,posinline++) {
	fout << static_cast<char>(value[fpos]);
      }
      if(posinline==79){
	fout << "\n                     ";
	posinline=21;
	firstwordinline=true;
      }
    } else {
      if (onnextline) {
	fout << "\n                     ";
	posinline=21;
	firstwordinline=true;
      }
      // there is enough place for this word, simply print it
      if(!firstwordinline) {
	fout << " ";
	posinline++;
      }
      fout << value.substr(fpos, tpos-fpos);
      posinline+=tpos-fpos;
      fpos+=tpos-fpos+1;  // NOT spaceneeded, +1 to jump over space
      firstwordinline=false;
    }

  }
  if(posinline>=79) fout << "\n                     ";

  fout << "\"\n";
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::dumpContigs(std::list<Contig> & clist, std::ostream & fout)
{
  FUNCSTART("void dumpContigs(std::list<Contig> & clist, std::ostream & fout)");
  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      try{
	if(saveme) {
	  fout << cle;
	}
      }
      catch (Notify n) {
	cerr << "Error while dumping " << cle.getContigName() << ".\n";
	n.handleError(THISFUNC);
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

void assout::saveAs_TYPE(std::list<Contig> & clist, const std::string & filename, const uint8 type, bool deleteoldfile)
{
  FUNCSTART("void saveAs_TYPE(std::list<Contig> & clist, const std::string & filename, const uint8 type, bool deleteoldfile)");

  std::ofstream fout;

  // nope, this doesn't speed up simple output
  //
  // char mybuf[1024*1024];
  // fout.rdbuf()->pubsetbuf(mybuf,1024*1024);

  if(!openFileForAppend(filename,fout, deleteoldfile)){
    if(type==Contig::AS_TCS) Contig::dumpTCS_Head(fout);
    if(type==Contig::AS_MAF) {
      ReadGroupLib::resetSaveStatus();
      Contig::dumpMAF_Head(fout);
    }
  }
  Contig::setCoutType(type);
  dumpContigs(clist,fout);

  FUNCEND();
}
void assout::saveAs_TYPE(Contig & con, const std::string & filename, const uint8 type, bool deleteoldfile)
{
  FUNCSTART("void saveAs_TYPE(Contig & con, const std::string & filename, const uint8 type)");

  std::ofstream fout;

  // nope, this doesn't speed up simple output
  //
  // char mybuf[1024*1024];
  // fout.rdbuf()->pubsetbuf(mybuf,1024*1024);

  if(!openFileForAppend(filename,fout, deleteoldfile)){
    if(type==Contig::AS_TCS) Contig::dumpTCS_Head(fout);
    if(type==Contig::AS_MAF) {
      ReadGroupLib::resetSaveStatus();
      Contig::dumpMAF_Head(fout);
    }
  }
  Contig::setCoutType(type);
  try{
    fout << con;
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::dumpAsACE(std::list<Contig> & clist, std::ostream & aceout)
{
  FUNCSTART("void dumpAsACE(std::list<Contig> & clist, std::ostream & aceout)");

  aceout << "AS " << clist.size();
  {
    uint32 sum=0;
    for(auto & cle : clist){
      sum+=cle.getNumReadsInContig();
    }
    aceout << " " << sum << endl << endl;
  }

  Contig::setCoutType(Contig::AS_ACE);
  dumpContigs(clist, aceout);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveAsACE(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)
{
  //declare a new file object for appending and overwriting
  std::fstream fio;

  uint32 oldnumcontigs=0;
  uint32 oldnumreads=0;
  saveAsACE_openACE(fio,
		    filename,
		    deleteoldfile,
		    oldnumcontigs,
		    oldnumreads);

  Contig::setCoutType(Contig::AS_ACE);
  dumpContigs(clist, fio);

  uint32 newnumreads=oldnumreads;
  for(auto & cle : clist){
    newnumreads+=cle.getNumReadsInContig();
  }

  saveAsACE_rewriteHeader(fio,
			  oldnumcontigs+static_cast<uint32>(clist.size()),
			  newnumreads);
}

void assout::saveAsACE(Contig & con, const std::string & filename, bool deleteoldfile)
{
  FUNCSTART("void saveAsACE(Contig & con, const std::string & filename, bool deleteoldfile)");

  //declare a new file object for appending and overwriting
  std::fstream fio;

  uint32 oldnumcontigs=0;
  uint32 oldnumreads=0;

  saveAsACE_openACE(fio,
		    filename,
		    deleteoldfile,
		    oldnumcontigs,
		    oldnumreads);
  Contig::setCoutType(Contig::AS_ACE);
  try{
    fio << con;
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
  }
  saveAsACE_rewriteHeader(fio,
			  oldnumcontigs+1,
			  oldnumreads+con.getNumReadsInContig());
  FUNCEND();
}

void assout::saveAsACE_openACE(std::fstream & fio, const std::string & filename, bool deleteoldfile, uint32 & numcontigs, uint32 & numreads)
{
  FUNCSTART("void saveAsACE_openACE(std::fstream & fio, const std::string & filename, uint32 & numcontigs, uint32 & numreads)");

  numcontigs=0;
  numreads=0;

  struct stat st;
  if(deleteoldfile || stat(filename.c_str(),&st)) {
    // oh, new file
    fio.open(filename, std::ios::out | std::ios::in | std::ios::trunc);
    // write a long empty line to reserve space for header
    fio << "                                                                                                                                                                           \n\n";
    return;
  }
  fio.open(filename, std::ios::ate | std::ios::out | std::ios::in);

  long fiosize = fio.tellp();
  fio.seekp(0);

  std::string dummy;
  if(fio.eof()){
    MIRANOTIFY(Notify::INTERNAL, "the ace file is present but seems to be empty: " << filename);
  }
  getline(fio,dummy);
  if(dummy.size()<50){
    MIRANOTIFY(Notify::INTERNAL, "first line is too short for rewriting: " << filename);
  }
  fio.seekp(0);
  fio >> dummy >> numcontigs >> numreads;
  fio.seekp(fiosize);
}

void assout::saveAsACE_rewriteHeader(std::fstream & fio, const uint32 numcontigs, const uint32 numreads)
{
  fio.seekp(0);
  //fio << "AS " << numcontigs << ' ' << numreads <<
  //  "                                                  ";

  std::string tmp="AS ";
  tmp+=boost::lexical_cast<std::string>(numcontigs);
  tmp+=" ";
  tmp+=boost::lexical_cast<std::string>(numreads);
  while(tmp.size()<50) tmp+=" ";
  fio << tmp;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::saveAsGAP4DA(std::list<Contig> & clist, const std::string & dirname, bool deleteolddir)
{
  FUNCSTART("void saveAsGAP4DA(std::list<Contig> & clist, const std::string & dirname)");

  if(ensureDirectory(dirname, deleteolddir)){
    MIRANOTIFY(Notify::FATAL, "Could not make sure that directory '" << dirname << "' exists, aborting MIRA.");
  }

  Contig::setCoutType(Contig::AS_GAP4DA);
  std::ofstream fofnout((dirname+"/fofn"), std::ios::out | std::ios::app);

  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      if(saveme){
	try{
	  cle.saveAsGAP4DA(dirname, fofnout);
	}
	catch (Notify n) {
	  cerr << "Error while dumping " << cle.getContigName() << ".\n";
	  n.handleError(THISFUNC);
	}
      }
    }
  }

  FUNCEND();
  return;
}

void assout::saveAsGAP4DA(Contig & con, const std::string & dirname, bool deleteolddir)
{
  FUNCSTART("void saveAsGAP4DA(Contig & con, const std::string & dirname)");

  if(ensureDirectory(dirname, deleteolddir)){
    MIRANOTIFY(Notify::FATAL, "Could not make sure that directory '" << dirname << "' exists, aborting MIRA.");
  }

  Contig::setCoutType(Contig::AS_GAP4DA);
  std::ofstream fofnout((dirname+"/fofn"), std::ios::out | std::ios::app);
  try{
    con.saveAsGAP4DA(dirname, fofnout);
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << ".\n";
    n.handleError(THISFUNC);
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

void assout::saveAsWiggle(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile, bool gcinsteadcov)
{
  FUNCSTART("void saveAsWiggle(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile)");

  {
    std::ofstream fout;
    openFileForAppend(filename,fout, deleteoldfile);
  }

  for(uint32 savewhat=0; savewhat<2; savewhat++){
    for(auto & cle : clist){
      bool saveme=false;
      if(savewhat==0 && cle.getNumReadsInContig()>1) saveme=true;
      if(savewhat==1 && cle.getNumReadsInContig()==1) saveme=true;
      try{
	if(saveme){
	  saveAsWiggle(cle, filename, false, gcinsteadcov);
	}


      }
      catch (Notify n) {
	cerr << "Error while dumping " << cle.getContigName() << ".\n";
	n.handleError(THISFUNC);
      }
    }
  }

  FUNCEND();
}

void assout::saveAsWiggle(Contig & con, const std::string & filename, bool deleteoldfile, bool gcinsteadcov)
{
  FUNCSTART("void saveAsWiggle(Contig & con, const std::string & filename, bool deleteoldfile)");

  std::ofstream fout;
  openFileForAppend(filename,fout, deleteoldfile);

  try{
    std::vector<int32> strainidsofbackbone;
    con.getStrainsOfBackbone(strainidsofbackbone);
    int32 bbstrainid=-1;
    if(strainidsofbackbone.size()==1) bbstrainid=strainidsofbackbone.front();

    std::string consseq;
    std::vector<base_quality_t> dummy;

    con.calcConsensi(0,0,'X');
    con.newConsensusGet(consseq,
			dummy,
			bbstrainid);

    if(gcinsteadcov){
      con.dumpGCWiggle_Body(fout, consseq);
    }else{
      con.dumpWiggle_Body(fout, consseq, true);
    }
  }
  catch (Notify n) {
    cerr << "Error while saving " << con.getContigName() << " as wiggle.\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::dumpHTMLHeader(const std::string & projectname, std::ostream & htmlout)
{
  htmlout << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\">\n\
<html>\n\
<head>\n\
   <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n\
   <meta name=\"GENERATOR\" content=\"MIRA (c) Bastien Chevreux & EdIt (c) Thomas Pfisterer;\">\n\
   <meta name=\"Author\" content=\"";

  //char * logname= getenv("LOGNAME");
  //htmlout << getlogin() << "\">\n

  {
    // REMARK:
    // GCC 3.4 tells this
    // : warning: Using 'getpwuid' in statically linked applications
    //   requires at runtime the shared libraries from the glibc version
    //   used for linking
    //
    //  This might also be the reason for a reported crash of
    //   miraconvert on a 2.6 kernel (with another glibc than my
    //   home machine
    //
    // Resolve: get back to getlogin() even if manual says it can be
    //  easily "fooled". This application is not security critical.

    //struct passwd * pws = getpwuid(getuid());
    //bool noname=true;
    //if(pws != nullptr) {
    //	if(pws->pw_name != nullptr && pws->pw_gecos != nullptr) {
    //	  htmlout << pws->pw_name << " (" << pws->pw_gecos << ")";
    //	  noname=false;
    //	}
    //}

    char * namestr=getlogin();
    if(namestr==nullptr) {
      namestr=getenv("LOGNAME");
      if(namestr==nullptr) {
	namestr=getenv("USER");
      }
    }

    if(namestr!=nullptr){
      htmlout << namestr;
    }else{
      htmlout << "unknown";
    }
  }

  htmlout << "\">\n<meta name=\"Description\" content=\"Assembled shotgun project\">\n\
   <title>";
  htmlout << "Project " << projectname << " </title>\n\
  <STYLE TYPE=\"text/css\">\n\
  <!--\n\
  \n\
   .FCDS {color:black;  background-color:#4AA090;}\n\
   .FrRN {color:black;  background-color:#f41e8e;}\n\
   .FtRN {color:black;  background-color:#736cdc;}\n\
   .FmxR {color:black;  background-color:#653BD9;}\n\
   .MISM {color:red;  background-color:#dddddd;}\n\
   .SRMr {color:black;  background-color:#ff5050;}\n\
   .SRMc {color:black;  background-color:#ff5050;}\n\
   .WRMr {color:black;  background-color:orange;}\n\
   .WRMc {color:black;  background-color:orange;}\n\
   .SROr {color:black;  background-color:#00ced1;}\n\
   .SROc {color:black;  background-color:#00ced1;}\n\
   .SAOr {color:black;  background-color:#2e8b57;}\n\
   .SAOc {color:black;  background-color:#2e8b57;}\n\
   .SIOr {color:black;  background-color:#98fb98;}\n\
   .SIOc {color:black;  background-color:#98fb98;}\n\
   .MCVc {color:black;  background-color:#cc3333;}\n\
   .POLY {color:black;  background-color:#ffff99;}\n\
   .EDxD {color:black;  background-color:#db7093;}\n\
   .EDxI {color:black;  background-color:#db7093;}\n\
   .EDxC {color:black;  background-color:#db7093;}\n\
   .IUPC {color:black;  background-color:#cccccc;}\n\
\n\
BODY  { font-family: sans-serif;\n\
  color: #000000 ;\n\
}\n\
\n\
   .jtable1 {\n\
     color : black; \n\
     background-color : #cccccc ;\n\
     font-size: normal ;\n\
     font-style: normal ;\n\
     font-family: sans-serif ; \n\
     font-weight: normal ;\n\
     text-align: left ; \n\
     vertical-align: top ;\n\
     padding: 10px;\n\
   }\n\
   .jtable2 {\n\
     color : black; \n\
     background-color : #eeeeee ;\n\
     font-size: normal ;\n\
     font-style: normal ;\n\
     font-family: sans-serif ; \n\
     font-weight: normal ;\n\
     text-align: left ; \n\
     vertical-align: top ;\n\
     padding: 10px;\n\
   }\n\
   .jtable3 {\n\
     color : black; \n\
     background-color : white ;\n\
     font-size: normal ;\n\
     font-style: normal ;\n\
     font-family: sans-serif ; \n\
     font-weight: normal ;\n\
     text-align: left ; \n\
     vertical-align: top ;\n\
     padding: 10px;\n\
   }\n\
\n\
  -->\n\
</STYLE>\n\
</head>\n\
<body TEXT=\"#000000\" BGCOLOR=\"#FFFFFF\" LINK=\"#FF0000\" VLINK=\"#551A8B\" ALINK=\"#000088\">\n";

  //   .ALUS {color:black;  background-color:#90ee90;}\n
  // #66ffff = azure (pi*daumen)

  htmlout << "<h1><center>Tag legend</center></h1>\n";

  htmlout << "<center>\n";
  htmlout << "<table CELLSPACING=0 CELLPADDING=0 NOSAVE >\n";

  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"FCDS\">&nbsp;</SPAN> = FCDS;</tt></td><td>Feature CDS (coding sequence)</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"FtRN\">&nbsp;</SPAN> = FtRN;</tt></td><td>tRNA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"FrRN\">&nbsp;</SPAN> = FrRN;</tt></td><td>rRNA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"FmxR\">&nbsp;</SPAN> = Fm-R;</tt></td><td>misc. RNA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"MISM\">&nbsp;</SPAN> = MISM;</tt></td><td>Mismatch (discrepancy) between reads and consensus</td></tr>\n";
  //htmlout << "<tr><td><tt><SPAN CLASS=\"ALUS\">&nbsp;</SPAN> = ALUS;</tt> Repetitive ALU sequence</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"SRMr\">&nbsp;</SPAN> = SRMx;</tt></td><td>Strong Repeat Marker Base set by MIRA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"WRMr\">&nbsp;</SPAN> = WRMx;</tt></td><td>Weak Repeat Marker Base set by MIRA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"SROr\">&nbsp;</SPAN> = SROx;</tt></td><td>SNP inteR Organism (Read/Consensus) set by MIRA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"SAOr\">&nbsp;</SPAN> = SAOx;</tt></td><td>SNP intrA Organism (Read/Consensus) set by MIRA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"SIOr\">&nbsp;</SPAN> = SIOx;</tt></td><td>SNP Inter- and intra-Organism (Read/Consensus) set by MIRA</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"MCVc\">&nbsp;</SPAN> = MCVc;</tt></td><td>Missing CoVerage in Consensus (set by MIRA)</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"POLY\">&nbsp;</SPAN> = POLY;</tt></td><td>Poly-A signal</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"EDxD\">&nbsp;</SPAN> = EDxD;</tt></td><td>Delete operation set by EdIt</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"EDxI\">&nbsp;</SPAN> = EDxI;</tt></td><td>Insert operation set by EdIt</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"EDxC\">&nbsp;</SPAN> = EDxC;</tt></td><td>Change operation set by EdIt</td></tr>\n";
  htmlout << "<tr align=\"left\"><td><tt><SPAN CLASS=\"IUPC\">&nbsp;</SPAN> = IUPAC;</tt></td><td> IUPAC base (shows only in HTML output)</td></tr>\n";

  htmlout<< "</table></center>\n";
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void assout::dumpContigListAsHTML(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile, const std::string & projectname)
{
  FUNCSTART("void dumpContigListAsHTML(std::list<Contig> & clist, const std::string & filename, bool deleteoldfile, const std::string & projectname)");

  std::ofstream fout;
  if(!openFileForAppend(filename,fout, deleteoldfile)){
    dumpHTMLHeader(projectname, fout);
  }

  // A contig list at the top of the file is currently not possible
  //  anymore *sigh*
  // maybe splitting into two files and then concatenate at the of the process?
  //
  //std::list<Contig>::iterator I=clist.begin();
  //htmlout << "<h1><center>Contig Std::List</center></h1>\n";
  //while(I!=clist.end()){
  //
  //  if(I->getContigReads().size() > 1) {
  //	htmlout << "<a href=\"#" << I->getContigName() << "\">Contig " <<  I->getContigID();
  //  }else{
  //	htmlout << "<a href=\"#" << I->getContigName() << "\">Singlet " <<  I->getContigID();
  //  }
  //  htmlout << " (" << I->getContigLength() << ")</a>";
  //  I++;
  //  if(I!=clist.end()){
  //	htmlout << ", ";
  //  }else{
  //	htmlout << "\n<p>\n";
  //  }
  //}

  Contig::setCoutType(Contig::AS_HTML);
  dumpContigs(clist,fout);

  // This is also bad ... when should the HTML be closed?
  //fout << "\n</body></html>";

  return;

  FUNCEND();
}


void assout::dumpContigAsHTML(Contig & con, const std::string & filename, bool deleteoldfile, const std::string & projectname)
{
  FUNCSTART("void dumpContigAsHTML(Contig & con, const std::string & filename, bool deleteoldfile, const std::string & projectname)");

  std::ofstream fout;
  if(!openFileForAppend(filename,fout, deleteoldfile)){
    dumpHTMLHeader(projectname, fout);
  }

  // A contig list at the top of the file is currently not possible
  //  anymore *sigh*
  // maybe splitting into two files and then concatenate at the of the process?

  Contig::setCoutType(Contig::AS_HTML);
  try{
    fout << con;
  }
  catch (Notify n) {
    cerr << "Error while dumping " << con.getContigName() << " as HTML.\n";
    n.handleError(THISFUNC);
  }

  // This is also bad ... when should the HTML be closed?
  //fout << "\n</body></html>";

  return;

  FUNCEND();
}
