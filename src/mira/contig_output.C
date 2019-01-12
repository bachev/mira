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


#include "mira/contig.H"

#include "util/dptools.H"
#include "util/misc.H"
#include "mira/ads.H"
#include "mira/gff_parse.H"
#include "mira/sam_collect.H"

#include "util/stlimprove.H"

#include <boost/format.hpp>

using std::cout;
using std::cerr;
using std::endl;

//#define CEBUGFLAG
#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpStats(std::ostream &ostr)
{
  FUNCSTART("std::ostream & Contig::dumpStats(std::ostream &ostr, Contig const  &con)");

  calcStats();

  // shouldn't be needed as setprecision() works only only once (Josuttis)
  //  but apparently gcc does not honour this
  auto oldprecision=cout.precision();

  if(CON_outtype==AS_HTML){
    ostr << "<h2>Statistics</h2><p>\n";
    ostr << "To be reworked!\n";
    //ostr << "Contig id: " << CON_id << "<br>" << '\n';
    //ostr << "Contig length: " << CON_counts.size() << "<br>" << '\n';
    //ostr << "Avg. contig coverage: " << CON_stats.avg_coverage << "<br>" << '\n';
    //
    //ostr << "Consensus contains: A: " << CON_stats.AinC << "\tC: " << CON_stats.CinC << "\tG: " << CON_stats.GinC << "\tT: " << CON_stats.TinC << "\tN: " << CON_stats.NinC << "\tIUPAC: " << CON_stats.IUPACinC << "\tFunny: " << CON_stats.FunnyInC << "\t*: " << CON_stats.starInC << "<p>" << "\n";
    //
    //
    //ostr << "\nNum reads: " << CON_stats.total_reads << "<br>" << '\n';
    //
    //ostr << "Avg. read length: " << CON_stats.avg_readlen << "<br>" << '\n';
    //
    //ostr << "Reads contain " << CON_stats.total_readlen-CON_stats.starInR-CON_stats.NinR << " bases, " << CON_stats.NinR << " Ns and " << CON_stats.starInR << " gaps." << "<p>" << "\n";
    //
    //ostr << "<p>\n";
  }else{
    ostr << "\n-------------- Contig statistics ----------------\n";
    ostr << "Contig id: " << CON_id << '\n';;

    ostr << "Contig length: " << CON_counts.size() << "\n\n";

    ostr << "\t\t";
    for(uint32 i=0;i<ReadGroupLib::SEQTYPE_END; i++){
      ostr << std::setw(12) << ReadGroupLib::getNameOfSequencingType(i);
    }
    ostr << "\nNum. reads\t";
    for(uint32 i=0;i<ReadGroupLib::SEQTYPE_END; i++){
      ostr << std::setw(12) << CON_stats.readsperst[i];
    }

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

    ostr << "\n100% merged reads\t"
	 << std::setw(4) << "-" << std::setw(12) << "-" << std::setw(12) << "-" << std::setw(12) << "-" << std::setw(12) << "-" << std::setw(12) << "-"
	 << std::setw(12) << CON_nummergedreads_perseqtype[0]
      ;

    ostr << "\nMax. coverage\t";
    for(uint32 i=0;i<ReadGroupLib::SEQTYPE_END; i++){
      ostr << std::setw(12) << CON_stats.max_covperst[i];
    }

    ostr << "\nAvg. coverage\t";
    for(uint32 i=0;i<ReadGroupLib::SEQTYPE_END; i++){
      ostr << std::setw(12) << std::setprecision(3) << CON_stats.avg_covperst[i];
    }

    ostr << "\n\nMax. contig coverage: " << CON_stats.max_coverage << '\n';
    ostr << "Avg. contig coverage: " << std::setprecision(3) << CON_stats.avg_coverage << '\n';

    ostr << "\nConsensus contains:\tA: " << CON_stats.AinC << "\tC: " << CON_stats.CinC << "\tG: " << CON_stats.GinC << "\tT: " << CON_stats.TinC << "\tN: " << CON_stats.NinC << "\n\t\t\tIUPAC: " << CON_stats.IUPACinC << "\tFunny: " << CON_stats.FunnyInC << "\t*: " << CON_stats.starInC << "\n";

    ostr << "GC content: ";
    if(CON_stats.CinC+CON_stats.GinC>0){
      ostr << 100.0/(CON_stats.AinC+CON_stats.CinC+CON_stats.GinC+CON_stats.TinC)*(CON_stats.CinC+CON_stats.GinC) << "%\n";
    }else{
      ostr << "0%\n";
    }

    if(CON_stats.numnocoverage>=0){
      ostr << "\nContig positions with coverage: " << CON_stats.numnocoverage << '\n';
    }

    //ostr << "\nReads contain " << CON_stats.total_readlen-CON_stats.starInR-CON_stats.NinR << " bases, " << CON_stats.NinR << " Ns and " << CON_stats.starInR << " gaps.\n";

    ostr << "-------------------------------------------------\n";
  }

  cout.precision(oldprecision);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::setCoutType(uint8 type)
{
  FUNCSTART("void Contig::setCoutType(uint8 type)");
  switch(type){
  case AS_TEXT:
  case AS_DEBUG:
  case AS_FASTA:
  case AS_FASTAQUAL:
  case AS_FASTAPADDED:
  case AS_FASTAPADDEDQUAL:
  case AS_HTML:
  case AS_CAF:
  case AS_MAF:
  case AS_ACE:
  case AS_TCS:
  case AS_GAP4DA: {
    CON_outtype=type;
    break;
  }
  default:{
    MIRANOTIFY(Notify::INTERNAL, "Wrong type is not one of TEXT, HTML, CAF. ACE or GAP4DA.");
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

std::ostream & operator<<(std::ostream &ostr, Contig const  &con)
{
  FUNCSTART("friend std::ostream & Contig::operator<<(std::ostream &ostr, Contig const  &con)");

  Contig & nonconstcon = const_cast<Contig &>(con);

  //nonconstcon.definalise();
  nonconstcon.finalise();

  if(con.CON_outtype==Contig::AS_DEBUG){
    nonconstcon.priv_dumpAsDebug(ostr);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_CAF){
    nonconstcon.priv_dumpAsCAF(ostr);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_MAF){
    nonconstcon.priv_dumpAsMAF(ostr);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_ACE){
    nonconstcon.priv_dumpAsACE(ostr);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_TCS){
    nonconstcon.priv_dumpAsTCS(ostr);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_FASTA){
    nonconstcon.priv_dumpAsFASTA(ostr,false);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_FASTAQUAL){
    nonconstcon.priv_dumpAsFASTAQUAL(ostr,false);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_FASTAPADDED){
    nonconstcon.priv_dumpAsFASTA(ostr,true);
    FUNCEND();
    return ostr;
  }
  if(con.CON_outtype==Contig::AS_FASTAPADDEDQUAL){
    nonconstcon.priv_dumpAsFASTAQUAL(ostr,true);
    FUNCEND();
    return ostr;
  }


  if(con.CON_outtype==Contig::AS_HTML){
    ostr << "<a NAME=\"" << con.getContigName() << "\"></a>\n";
    ostr << "<h1><center>" << con.getContigName() << "</center></h1>\n";
  }

  nonconstcon.dumpStats(ostr);

  std::string consseq;
  std::vector<base_quality_t> consqual;
  nonconstcon.newConsensusGet(consseq, consqual);

  // TODO: ??? why ???
  //for(uint32 i=0; i<consseq.size(); i++){
  //  nonconstcon.CON_tmpcons[i]=consseq[i];
  //}


  if(con.CON_finalised==true){
    nonconstcon.priv_dumpAsTextOrHTML(ostr, con.CON_outtype, consseq, consqual, 0, nonconstcon.getContigLength());
  }else{
    ostr << "Consensus not finalised, no more information to output.\n";
  }

  FUNCEND();

  return ostr;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpAsText(std::ostream &ostr, int32 frompos, int32 topos)
{
  if(topos < frompos || frompos > static_cast<int32>(CON_counts.size())) return;

  std::string consseq;
  std::vector<base_quality_t> consqual;
  newConsensusGet(consseq, consqual);
  priv_dumpAsTextOrHTML(ostr, Contig::AS_TEXT, consseq, consqual, frompos, topos);
}

void Contig::dumpAsHTML(std::ostream &ostr, int32 frompos, int32 topos)
{
  if(topos < frompos || frompos > static_cast<int32>(CON_counts.size())) return;

  std::string consseq;
  std::vector<base_quality_t> consqual;
  newConsensusGet(consseq, consqual);
  priv_dumpAsTextOrHTML(ostr, Contig::AS_HTML, consseq, consqual, frompos, topos);
}

void Contig::priv_dumpAsTextOrHTML(std::ostream &ostr, const uint8 outtype, const std::string & consseq, const std::vector<base_quality_t> & consqual, int32 frompos, int32 topos)
{
  FUNCSTART("void Contig::priv_dumpAsTextOrHTML(std::ostream &ostr, const uint8 outtype, const std::string & consseq, const std::vector<base_quality_t> & consqual, int32 frompos, int32 topos)");

  (void) consqual;

  if(topos < frompos || frompos > static_cast<int32>(CON_counts.size())) return;
  if(frompos<0) frompos=0;
  if(topos<0 || topos > static_cast<int32>(CON_counts.size())) topos=static_cast<int32>(CON_counts.size());

  if(outtype==Contig::AS_TEXT){
    ostr << "Sequence:\n";
  }else{
    ostr << "<h2>Sequence:</h2><p>\n";
  }

  const contig_parameters & con_params = (*(CON_miraparams))[0].getContigParams();

  int32 cpl=con_params.con_output_text_cpl;
  char egfillchar=con_params.con_output_text_gapfill;
  if(outtype==Contig::AS_HTML){
    cpl=con_params.con_output_html_cpl;
    egfillchar=con_params.con_output_html_gapfill;
  }
  if(cpl>topos-frompos) cpl=topos-frompos;
  if(egfillchar=='\t') egfillchar=' ';

  size_t longestnamesize=0;
  for(auto & cre : CON_reads){
    longestnamesize=std::max(longestnamesize,cre.getName().size());
  }
  ++longestnamesize;

  std::string spanid; // for HTML output

  // initialise the iterator for getting through the contig
  rcci_t rcci(this,
	      nullptr,  // all strainids
	      nullptr,  // all seqtypes
	      CON_outputrails,            // take rails only on demand
	      true,     // take backbones
	      true);   // take reads without readpool-reads
  rcci.advance(frompos);

  for(int32 conpos=frompos; conpos < topos; conpos+=cpl){
    // collect all pcrIs which will be in the view of cpl columns
    // take the pcrIs in the current column ...
    auto pcrIsinthisview=rcci.getPCRIsInCol();
    // ... then advance cpl-1 columns and add to this the newly arrived pcrIs
    for(int32 aci=0; aci<cpl-1; ++aci){
      rcci.advance();
      for(auto pI=rcci.getItToNewPCRIs(); pI!=rcci.getPCRIsInCol().end(); ++pI){
	pcrIsinthisview.push_back(*pI);
      }
    }
    // we advanced rcci by cpl-1, must got one further for the view
    //  (next iteration of for loop)
    rcci.advance();

    if(outtype==Contig::AS_HTML) ostr << "<table CELLSPACING=0 CELLPADDING=0 NOSAVE >\n";
    if(outtype==Contig::AS_HTML){
      ostr << "<tr NOSAVE>\n<td NOSAVE><div align=right>";
    }
    if(outtype==Contig::AS_TEXT){
      ostr << std::right << std::setw(static_cast<int>(longestnamesize));
    }
    ostr << conpos << ' ';
    if(outtype==Contig::AS_HTML) ostr << "</div></td>\n";
    CEBUG("\t\t\t\t\t\t");
    if(outtype==Contig::AS_TEXT){
      //ostr << "\t\t\t";
      ostr << ' ';
      for(int32 i=0; i<((cpl-1)/10)+1; i++) {
	ostr << "|    .    ";
      }
      ostr << '\n';
    }else if(outtype==Contig::AS_HTML){
      ostr << "<td><tt>";
      for(int32 i=0; i<((cpl-1)/10)+1; i++) {
	ostr << "|&nbsp;&nbsp;&nbsp;&nbsp;.&nbsp;&nbsp;&nbsp;&nbsp;";
      }
      ostr << "</tt></td>\n</tr>\n";
    }

    //for(uint32 i=0; i<CON_reads.size(); i++){
    for(auto & pcrI : pcrIsinthisview){
      if(!CON_outputrails && pcrI->isRail()) continue;

      //Read & actread=CON_reads[origindex].read;

      //if((conpos >= CON_outputorder[i].offset_start
      //	  || conpos+cpl >= CON_outputorder[i].offset_start)
      //	 && conpos<static_cast<int32>(CON_outputorder[i].offset_start+actread.getLenClippedSeq())){
      if(1){


#ifdef CEBUGFLAG
	ostr << "offset: " << pcrI.getReadStartOffset()
	ostr << "\tsize: " << pcrI->getLenClippedSeq();
	ostr << "\torpid: " << pcrI.getORPID()
	ostr << "\t";
#endif

	if(outtype==Contig::AS_HTML) {
	  ostr << "<tr><td>";
	}else if(outtype==Contig::AS_TEXT){
	  ostr << std::left << std::setw(static_cast<int>(longestnamesize));
	}
	ostr << pcrI->getName();

	std::vector<char>::const_iterator oI;
	if(pcrI.getReadDirection() > 0){
	  oI=pcrI->getClippedSeqIterator();
	  ostr << '+';
	}else{
	  oI=pcrI->getClippedComplementSeqIterator();
	  ostr << '-';
	}
	if(outtype==Contig::AS_TEXT){
	  ostr << ' ';
	}else if(outtype==Contig::AS_HTML){
	  ostr << "</td>\n<td><tt>";
	}

	int32 index = conpos-pcrI.getReadStartOffset();
	advance(oI, index);

	for(int32 linej=0; linej<cpl; ++oI, ++linej){
	  if(index+linej >= 0
	     && index+linej < static_cast<int32>(pcrI->getLenClippedSeq())){
	    char tmp=static_cast<char>(tolower(*oI));
	    if(tmp!=consseq[conpos+linej]) tmp=static_cast<char>(toupper(*oI));
	    bool setspan=false;
	    if(outtype==Contig::AS_HTML){
	      int32 rp=pcrI.contigPos2UnclippedReadPos(conpos+linej);

	      if(pcrI.getReadDirection() < 0){
		rp=pcrI->getLenSeq()-rp-1;
	      }

	      spanid.clear();
	      int32 spanprecedence=9999;
	      for(uint32 tagi=0; tagi < pcrI->getNumOfTags() && spanprecedence != 0; ++tagi){
		const multitag_t & thetag=pcrI->getTag(tagi);
		if(rp >= static_cast<int32>(thetag.from)
		   && rp <= static_cast<int32>(thetag.to)){
		  if(thetag.identifier==Read::REA_defaulttag_ED_D.identifier){
		    spanid="EDxD";
		    setspan=true;
		    spanprecedence=5;
		  }
		  if(thetag.identifier==Read::REA_defaulttag_ED_C.identifier){
		    spanid="EDxC";
		    setspan=true;
		    spanprecedence=5;
		  }
		  if(thetag.identifier==Read::REA_defaulttag_ED_I.identifier){
		    spanid="EDxI";
		    setspan=true;
		    spanprecedence=5;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSRMr){
		    if(dptools::isValidStarBase(consseq[conpos+linej])
		       && toupper(consseq[conpos+linej])!=toupper(tmp)){
		      spanid="MISM";
		    }else{
		      spanid="SRMr";
		    }
		    setspan=true;
		    spanprecedence=1;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idWRMr){
		    spanid="WRMr";
		    setspan=true;
		    spanprecedence=1;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSROr){
		    spanid="SROr";
		    setspan=true;
		    spanprecedence=0;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSAOr){
		    spanid="SAOr";
		    setspan=true;
		    spanprecedence=0;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSIOr){
		    spanid="SIOr";
		    setspan=true;
		    spanprecedence=0;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSOFApolyA_sequence){
		    spanid="POLY";
		    setspan=true;
		    spanprecedence=3;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSOFACDS){
		    spanid="FCDS";
		    setspan=true;
		    spanprecedence=3;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSOFAtRNA){
		    spanid="FtRNA";
		    setspan=true;
		    spanprecedence=3;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSOFArRNA){
		    spanid="FrRN";
		    setspan=true;
		    spanprecedence=3;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idSOFAtranscript){
		    spanid="FmxR";
		    setspan=true;
		    spanprecedence=3;
		  }
		  if(thetag.identifier==Read::REA_tagentry_idALUS){
		    spanid="ALUS";
		    setspan=true;
		    spanprecedence=4;
		  }
		}
	      }
	      if(setspan==false){
		if(toupper(consseq[conpos+linej])!=toupper(tmp)){
		  spanid="MISM";
		  setspan=true;
		}
	      }
	      if(setspan){
		ostr << "<SPAN CLASS=\"" << spanid << "\">";
	      }
	    }

	    ostr << tmp;

	    if(setspan){
	      ostr << "</SPAN>";
	    }
	  }else{
	    if(outtype==Contig::AS_TEXT){
	      ostr << egfillchar;
	    }else if(outtype==Contig::AS_HTML){
	      switch(egfillchar) {
	      case ' ' : {
		ostr << "&nbsp;";
		break;
	      }
	      default : {
		ostr << egfillchar;
	      }
	      }
	    }

	    // we're already over the limit, oI should point to .end()
	    // this has absolutely no relevance in normal code, but when compiling
	    //  in _GLIBCXX_DEBUG mode, we get an error complaining about the
	    //  fact that we try to increase the iterator past .end()
	    // therefore, reduce it before the next loop
	    if(index+linej >= static_cast<int32>(pcrI->getLenClippedSeq())){
	      --oI;
	    }
	  }
	}

	if(outtype==Contig::AS_TEXT){
	  ostr << '\n';
	}else if(outtype==Contig::AS_HTML){
	  ostr << "</tt></td></tr>\n";
	}
      }
    }
    CEBUG("\t\t\t\t\t\t");
    if(outtype==Contig::AS_TEXT){
      ostr << std::setw(static_cast<int>(longestnamesize+2)) << ' ';
      for(int32 j=0; j<cpl; ++j){
	ostr << "-";
      }
      ostr << '\n' << std::left << std::setw(static_cast<int>(longestnamesize+2)) << "Consensus:";
    }else if(outtype==Contig::AS_HTML){
      ostr << "<tr><td>Consensus:</td>\n<td><tt>";
    }

    CEBUG("\t\t\t\t\t\t");

    for(int32 conj=0; conj<cpl; ++conj){
      if(conpos+conj<topos){
	bool setspan=false;
	if(outtype==Contig::AS_HTML){
	  for(size_t ctagi=0; ctagi < CON_consensus_tags.size(); ++ctagi){
	    if(conpos+conj>=static_cast<int32>(CON_consensus_tags[ctagi].from) && conpos+conj<=static_cast<int32>(CON_consensus_tags[ctagi].to)){
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idMCVc){
		ostr << "<SPAN CLASS=\"MCVc\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idED_D){
		ostr << "<SPAN CLASS=\"EDxD\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idED_C){
		ostr << "<SPAN CLASS=\"EDxC\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idED_I){
		ostr << "<SPAN CLASS=\"EDxI\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idSRMc){
		ostr << "<SPAN CLASS=\"SRMc\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idWRMc){
		ostr << "<SPAN CLASS=\"WRMc\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idSROc){
		ostr << "<SPAN CLASS=\"SROc\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idSAOc){
		ostr << "<SPAN CLASS=\"SAOc\">";
		setspan=true;
		break;
	      }
	      if(CON_consensus_tags[ctagi].identifier==CON_tagentry_idSIOc){
		ostr << "<SPAN CLASS=\"SIOc\">";
		setspan=true;
		break;
	      }
	    }
	  }
	}
	if(outtype==Contig::AS_HTML
	   && !setspan
	   && !dptools::isValidBase(consseq[conpos+conj])
	   && dptools::isValidIUPACBase(consseq[conpos+conj])) {
	  setspan=true;
	  ostr << "<SPAN CLASS=\"IUPC\">";
	}
	ostr << consseq[conpos+conj];
	if(setspan){
	  ostr << "</SPAN>";
	}
      }else{
	if(outtype==Contig::AS_TEXT){
	  ostr << egfillchar;
	}else if(outtype==Contig::AS_HTML){
	  switch(egfillchar) {
	  case ' ' : {
	    ostr << "&nbsp;";
	    break;
	  }
	  default : {
	    ostr << egfillchar;
	  }
	  }
	}
      }
    }

    if(outtype==Contig::AS_TEXT){
      ostr << "\n\n";
    }else if(outtype==Contig::AS_HTML){
      ostr << "</tt>\n</td></tr>\n";
    }
    if(outtype==Contig::AS_HTML) ostr << "\n</table>\n<p>\n";
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpConReads()
{
  cout << "%%% dumping contig reads" << '\n';
  auto crI=CON_reads.begin();
  for(uint32 i=0; crI != CON_reads.end(); ++crI, ++i){
    cout << i << "\t" << crI.getORPID() << "\t" << crI.getReadStartOffset();
    cout << "\t" << crI.getReadDirection() << "\t" << crI->getName();
    if(crI->isBackbone()) cout << "\tbb";
    if(crI->isRail()) cout << "\trail";
    cout << '\n';
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsDebug(std::ostream & ostr)
{
  FUNCSTART("void Contig::priv_dumpAsDebug(std::ostream & ostr)");

  // TODO: adapt to new struct

  uint32 count=0;
  for(auto ccI=CON_counts.cbegin(); ccI != CON_counts.cend(); ++count, ++ccI){
    ostr << count << ":\t" << ccI->A
	 << '\t' << ccI->C
	 << '\t' << ccI->G
	 << '\t' << ccI->T
	 << '\t' << ccI->N
	 << '\t' << ccI->X
	 << '\t' << ccI->star
	 << '\t' << ccI->total_cov
	 << '\t' << ccI->baselock
	 << '\t' << ccI->snplock
	 << '\t' << static_cast<char>(ccI->i_backbonecharorig)
	 << '\t' << static_cast<char>(ccI->i_backbonecharupdated)
	 << '\t' << static_cast<uint16>(ccI->i_backbonequalorig)
	 << "\tf " << ccI->bbcountsf[0]
	 << "\tr " << ccI->bbcountsr[0]
	 << '\n';
  }

  ostr.flush();
  FUNCEND();

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsTCS(std::ostream & ostr)
{
  FUNCSTART("void Contig::priv_dumpAsTCS(std::ostream & ostr)");

  std::string consseq;
  std::vector<base_quality_t> consqual;

  //OLDgetConsensus1(consseq, consqual, false, 0, 0, -1, '@', &ostr);

  ostr << "TCS output currently not available, please contact the author.\n";

  ostr.flush();
  FUNCEND();

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

std::string Contig::priv_buildFASTAQHeaderComment(uint32 len, bool alsoavgcov)
{
  std::string ret("    ");
  if(!CON_stats.statsvalid) calcStats();
  if(alsoavgcov) {
    ret+="cov="+boost::str(boost::format("%.2f") % CON_stats.avg_coverage) + " ";
  }
  ret+="len="+boost::str(boost::format("%1i") % len) + " ";
  ret+="gc="+boost::str(boost::format("%.2f") % CON_stats.gccontent) + " ";
  ret+="nseq="+boost::str(boost::format("%1i") % CON_reads.size());

  return ret;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsFASTA(std::ostream & ostr, bool padded)
{
  FUNCSTART("void Contig::priv_dumpAsFASTA(std::ostream & ostr, bool padded)");

  std::string consseq;
  std::vector<base_quality_t> consqual;

  newConsensusGet(consseq, consqual);

  // remove pads if wanted
  if(!padded){
    priv_depadSeqQual(consseq,nullptr);
  }

  ostr << ">" << getContigName() << priv_buildFASTAQHeaderComment(consseq.size(),true) << '\n';

  uint32 cpl=0;
  for(auto & se : consseq){
    ostr << se;
    if(++cpl==80){
      cpl=0;
      ostr << "\n";
    }
  }
  if(cpl!=0) ostr << "\n";

  ostr.flush();
  FUNCEND();

  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Contig::priv_depadSeqQual(std::string & seq, std::vector<base_quality_t> * qual)
{
  auto dstsI=seq.begin();
  auto srcsI=dstsI;
  if(qual==nullptr){
    for(; srcsI!=seq.end(); ++srcsI){
      *dstsI=*srcsI;
      if(*srcsI!='*'){
	++dstsI;
      }
    }
  }else{
    auto dstqI=qual->begin();
    auto srcqI=dstqI;
    for(; srcsI!=seq.end(); ++srcsI, ++srcqI){
      *dstqI=*srcqI;
      *dstsI=*srcsI;
      if(*srcsI!='*'){
	++dstsI;
	++dstqI;
      }
    }
    qual->resize(dstsI-seq.begin());
  }
  seq.resize(dstsI-seq.begin());
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsFASTAQUAL(std::ostream & ostr, bool padded)
{
  FUNCSTART("void Contig::priv_dumpAsFASTAQUAL(std::ostream & ostr, bool padded)");

  std::string consseq;
  std::vector<base_quality_t> consqual;

  // don't need to pass padded here, we just need info on stars *
  // BUT: it is probable that we can use the cached temporary cheat consensus!
  newConsensusGet(consseq, consqual);

  BUGIFTHROW(consseq.size() != consqual.size(), "sequence size != qual size ???");

  // remove pads if wanted
  if(!padded){
    priv_depadSeqQual(consseq, &consqual);
  }

  ostr << ">" << getContigName() << priv_buildFASTAQHeaderComment(consseq.size(),true) << '\n';

  uint32 cpl=0;
  for(auto & qe : consqual){
    ostr << static_cast<uint16>(qe) << " ";
    if(++cpl==20){
      cpl=0;
      ostr << "\n";
    }
  }
  if(cpl!=0) ostr << "\n";

  ostr.flush();
  FUNCEND();

  return;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_getStrainSeqQualForFASTAQ(std::string & strainseq, std::vector<base_quality_t> & strainqual, uint32 mincoverage, base_quality_t minqual, int32 strainidtotake, bool padded, bool fillholesinstrain)
{
  FUNCSTART("void Contig::priv_getStrainSeqQualForFASTAQ(std::string & strainseq, std::vector<base_quality_t> & strainqual, uint32 mincoverage, base_quality_t minqual, int32 strainidtotake, bool padded, bool fillholesinstrain)");

  calcConsensi(mincoverage,minqual,'X');
  newConsensusGet(strainseq, strainqual, strainidtotake);

  // fill holes of strainseq with data from consseq if wanted
  if(fillholesinstrain){
    calcConsensi(mincoverage,0,'X');
    std::string consseq;
    std::vector<base_quality_t> consqual;
    newConsensusGet(consseq, consqual);
    BUGIFTHROW(consseq.size() != strainseq.size(), "cons size != strain size ???");

    auto csI=consseq.cbegin();
    auto cqI=consqual.cbegin();
    auto ssI=strainseq.begin();
    auto sqI=strainqual.begin();
    for(; csI!=consseq.cend(); ++csI, ++cqI, ++ssI, ++sqI){
      auto up=static_cast<char>(toupper(*ssI));
      if(up=='N' || up=='X' || *ssI=='@') {
	*ssI=*csI;
	*sqI=*cqI;
      }
    }
  }

  if(!padded){
    priv_depadSeqQual(strainseq, &strainqual);
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpStrainAsFASTAQUAL(std::ostream & fout,std::ostream & qout, uint32 mincoverage, base_quality_t minqual, int32 strainidtotake, bool padded, bool fillholesinstrain)
{
  FUNCSTART("void Contig::dumpAsFASTA(std::ostream & fout, std::ostream & qout, base_quality_t minqual, int32 strainidtotake, bool padded, bool fillholesinstrain)");

  std::string strainseq;
  std::vector<base_quality_t> strainqual;

  priv_getStrainSeqQualForFASTAQ(strainseq,strainqual,mincoverage, minqual, strainidtotake, padded, fillholesinstrain);

  {
    auto tmp=std::string(">"+getContigName())+priv_buildFASTAQHeaderComment(strainseq.size(),true);
    fout << tmp << '\n';
    qout << tmp << '\n';
  }

  uint32 cpl=0;
  uint32 qpl=0;
  auto qI=strainqual.cbegin();
  for(auto * sptr=strainseq.c_str(); *sptr; ++sptr,++qI){
    fout << static_cast<char>(*sptr);
    qout << static_cast<uint16>(*qI) << ' ';
    if(++cpl==80){
      cpl=0;
      fout << "\n";
    }
    if(++qpl==20){
      qpl=0;
      qout << "\n";
    }
  }
  if(cpl!=0) fout << "\n";
  if(qpl!=0) qout << "\n";

  fout.flush();
  qout.flush();

  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpStrainAsFASTQ(std::ostream & fout, uint32 mincoverage, base_quality_t minqual, int32 strainidtotake, bool padded, bool fillholesinstrain)
{
  FUNCSTART("void Contig::dumpAsFASTq(std::ostream & fout, base_quality_t minqual, int32 strainidtotake, bool padded, bool fillholesinstrain)");

  std::string strainseq;
  std::vector<base_quality_t> strainqual;

  priv_getStrainSeqQualForFASTAQ(strainseq,strainqual,mincoverage, minqual, strainidtotake, padded, fillholesinstrain);

  fout << '@' << getContigName() << priv_buildFASTAQHeaderComment(strainseq.size(),true) << '\n';

  fout << strainseq << "\n+\n";

  for(auto & qe : strainqual){
    fout << static_cast<char>((qe)+33);
  }
  fout <<'\n';
  fout.flush();

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsCAF(std::ostream & ostr)
{
  FUNCSTART("void Contig::priv_dumpAsCAF(std::ostream & ostr)");

  finalise();

  // make the consensus here as the routines also may set tags (STRM)
  //  to the contig
  // but first delete the old tags

  deleteTagsInReads(CON_tagentry_idSTMS);
  deleteTagsInReads(CON_tagentry_idSTMU);
  std::string consstring;
  std::vector<base_quality_t> consqual;

  newConsensusGet(consstring, consqual);

  Read::setCoutType(Read::AS_CAF);

  // output rails first
  if(CON_outputrails){
    for(auto & cre : CON_reads){
      if(cre.isRail()) ostr << cre;
    }
  }

  //ostr << "// dubidu" << '\n';

  // Reads
  for(auto & cre : CON_reads){
    if(!cre.isRail()) ostr << cre;
  }


  ostr << "Sequence : " << getContigName() << '\n';
  ostr << "Is_contig\nPadded\n";

  // Assembled from
  for(auto crI=CON_reads.begin(); crI!=CON_reads.end(); ++crI){
    if(!CON_outputrails && crI->isRail()) continue;
    ostr << "Assembled_from " << crI->getName() << " ";
    if(crI.getReadDirection()>0){
      ostr << crI.getReadStartOffset()+1 << " " << crI.getReadStartOffset()+crI->getLenClippedSeq() << " ";
    }else{
      ostr << crI.getReadStartOffset()+crI->getLenClippedSeq() << " " << crI.getReadStartOffset()+1 << " ";
    }
    ostr << crI->getLeftClipoff()+1 << " " << crI->getRightClipoff() << '\n';
  }

  // Tags
  for(auto & cte : CON_consensus_tags) cte.dumpAsCAF(ostr);

  ostr << "\n";

  {
    ostr << "DNA : " << getContigName() << '\n';

    const char * cptr = consstring.c_str();

    //makeTmpConsensus(0, CON_counts.size()+1);
    //char * cptr = CON_tmpcons;

    uint32 cpl=0;
    while(*cptr){
      if(*cptr=='*'){
	ostr << "-";
      }else{
	ostr << *cptr;
      }
      cptr++;
      if(cpl++==59){
	cpl=0;
	ostr << "\n";
      }
    }
    ostr << "\n\n";
  }

  {
    ostr << "BaseQuality : " << getContigName() << '\n';

    //makeTmpConsensus(0, CON_counts.size()+1);
    //char * cptr = CON_tmpcons;

    uint32 cpl=0;
    for(auto qI=consqual.begin(); qI!=consqual.end(); ++qI){
      if(cpl++==24){
	ostr << static_cast<uint16>(*qI) << '\n';
	cpl=0;
      } else {
	ostr << static_cast<uint16>(*qI) << ' ';
      }
    }
    ostr << "\n\n";
  }

  ostr.flush();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsMAF(std::ostream & ostr)
{
  FUNCSTART("void Contig::priv_dumpAsMAF(std::ostream & ostr)");

  finalise();

  // make the consensus here as the routines also may set tags (STRM)
  //  to the contig
  // but first delete the old tags

  deleteTagsInReads(CON_tagentry_idSTMS);
  deleteTagsInReads(CON_tagentry_idSTMU);
  std::string consstring;
  std::vector<base_quality_t> consqual;

  newConsensusGet(consstring, consqual);

  // for each readgroup present, dump the readgrouplib info
  // using the "save" so as not to dump readgroups multiple times
  //
  // because of SAM conversion woes (the MIRA sam converter does not like
  //  definition of new readgroups when reads were already defined), do
  //  this for all know readgroups, not only the ones we have in this contig

  ReadGroupLib::saveAllReadGroupsAsMAF(ostr);

  Read::setCoutType(Read::AS_MAF);

  // CO = contig name
  // NR = num reads (optional)
  // LC = contig length
  // CT = Consensus Tag (identifier from to comment).

  // scraped AF = assembled from (like in CAF: readname from to)
  // AT = assemble to (like in CAF AF, but no readname: from to)
  //      furthermore, they follow RD/ER block immediately so that
  //      the parser does not need to perform time consuming string lookups

  // CS = sequence (in one line, gaps as '*')
  // CQ = base quality (in one line, FASTQ-33 format)

  // \\ = start of read data
  // // = end of read data

  // EC = End Contig (marker for MAF parsing, mandatory)

  ostr << "CO\t" << getContigName()
       << "\nNR\t" << CON_reads.size()
       << "\nLC\t" << consstring.size() << '\n';

  // Tags
  for(auto cte : CON_consensus_tags) cte.dumpAsMAF(ostr,"CT");

  ostr << "CS\t";
  ostr << consstring << '\n';

  ostr << "CQ\t";
  for(auto qv : consqual) ostr << static_cast<char>(qv+33);

  ostr << "\n\\\\\n";

  // output first rails (if wanted), then other reads
  // it==0 -> rails
  // it==1 -> all other reads
  int8 it=1;
  if(CON_outputrails) it=0;
  for(; it<2; ++it){
    for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
      bool outit=false;
      if(it==0 && pcrI->isRail()) outit=true;
      if(it==1 && !pcrI->isRail()) outit=true;
      if(outit){
	ostr << *pcrI;

	// and the AT line
	ostr << "AT\t";
	if(pcrI.getReadDirection() > 0){
	  ostr << pcrI.getReadStartOffset()+1 << '\t' << pcrI.getReadStartOffset()+pcrI->getLenClippedSeq() << '\t';
	}else{
	  ostr << pcrI.getReadStartOffset()+pcrI->getLenClippedSeq() << '\t' << pcrI.getReadStartOffset()+1 << '\t';
	}
	ostr << pcrI->getLeftClipoff()+1 << '\t' << pcrI->getRightClipoff() << '\n';
      }
    }
  }


  //ostr << "// dubidu" << '\n';

  ostr << "//\n";

  ostr << "EC\n";

  ostr.flush();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpAsSAM(std::ostream & ostr, const SAMCollect & samc, bool alsobackbone)
{
  FUNCSTART("void Contig::dumpAsSAM(std::ostream & ostr)");

  finalise();

/*
  currently tg_index does not understand this

  // find out which strain should be used for the consensus
  // i.e., strains which are not backbone nor rails
  try{
    std::vector<int32> sidsnotbackbone;
    for(uint32 rgint=0; rgint<CON_readsperreadgroup.size(); ++rgint){
      if(CON_readsperreadgroup[rgint]>0){
	auto rgid=ReadGroupLib::getReadGroupID(rgint);
	if(!rgid.isBackbone() && !rgid.isRail()){
	  //cout << "### Push back sid " << static_cast<int32>(rgid.getStrainID()) << endl;
	  sidsnotbackbone.push_back(rgid.getStrainID());
	}
      }
    }
    sort(sidsnotbackbone.begin(),sidsnotbackbone.end());
    auto lI=unique(sidsnotbackbone.begin(),sidsnotbackbone.end());
    sidsnotbackbone.resize(lI-sidsnotbackbone.begin());
    BUGIFTHROW(sidsnotbackbone.empty(),"Ooooops, found no strain? Maybe the input data is at fault, but MIRA should have caught that earlier.");
    if(sidsnotbackbone.size()>1){
      cout << "For contig " << getContigName() << ": found " << sidsnotbackbone.size() << " strains with data, but SAM can embedd only one consensus. Will not embedd consensus for this contig.\n";
    }else{
      // make the consensus here as the routines also may set tags (STRM)
      //  to the contig
      // but first delete the old tags

      deleteTagsInReads(CON_tagentry_idSTMS);
      deleteTagsInReads(CON_tagentry_idSTMU);
      string consstring;
      std::vector<base_quality_t> consqual;
      newConsensusGet(consstring, consqual,sidsnotbackbone[0]);

      ostr << getContigName() << "\t516\t" << getContigName() << "\t1\t50\t";

      // CIGAR
      uint32 mcount=0;
      uint32 pcount=0;
      for(auto seqptr=consstring.c_str(); *seqptr; ++seqptr){
	if(*seqptr != '*'){
	  if(pcount){
	    ostr << pcount << "P";
	    pcount=0;
	  }
	  ++mcount;
	}else{
	  if(mcount){
	    ostr << mcount << "M";
	    mcount=0;
	  }
	  ++pcount;
	}
      }
      if(pcount){
	ostr << pcount << "P";
      }
      if(mcount){
	ostr << mcount << "M";
      }

      ostr << "\t*\t0\t0\t";

      // SEQ
      for(auto seqptr=consstring.c_str(); *seqptr; ++seqptr){
	if(*seqptr!='*') ostr << *seqptr;
      }
      // QUAL
      auto qI=consqual.cbegin();
      for(auto seqptr=consstring.c_str(); *seqptr; ++seqptr, ++qI){
	if(*seqptr!='*') ostr << (*qI + 33);
      }
    }
  }
  catch(Notify n){
    ReadGroupLib::debugDumpReadGroupInfo(cout);
    cout << "\nOooops, error while looking for SAM consensus?\n\n";
    n.handleError(THISFUNC);
  }
*/

  auto ctI=CON_consensus_tags.begin();
  try{
    std::string seqstr;
    std::string qualstr;
    SAMCollect::samrinfo_t samri(false);
    for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
      if(pcrI->isBackbone() && !alsobackbone) continue;

      // first, contig tags up till position of current read
      for(; ctI!=CON_consensus_tags.end() && ctI->from <= pcrI.getReadStartOffset(); ++ctI){
	ctI->dumpAsSAM(ostr,getContigName());
      }

      if(pcrI->getTemplate().empty()){
	MIRANOTIFY(Notify::FATAL,"The read " << pcrI->getName() << " is without a template? This should not be!");
	ostr << pcrI->getName() << '\t';
      }else{
	ostr << pcrI->getTemplate() << '\t';
      }

      if(!samc.getSAMRInfo(pcrI->getName(),samri)){
	MIRANOTIFY(Notify::INTERNAL,"Could not collect samrinfo_t for read " << pcrI->getName() << " ???");
      }

      //cout << "actSI: " << samri;

      // OUCH, the FLAG!
      uint32 flag=samri.samflags;
      // reverse flag, easier to set here than in SAMCollect for, e.g., single reads
      if(pcrI.getReadDirection()<0){
	flag|=0x10;
      }
      ostr << flag;

      ostr << '\t' << samc.getContigName(samri) << '\t';

      ostr << pcrI.getReadStartOffset()+1
	   << '\t' << "255"                  // mapping quality: none, so 255
	   << '\t';

      // CIGAR, seq & qual string
      {
	// Calculate preceeding 'S' statement (if any)
	// Cannot use the left and right clips as these may contain gaps!
	//  E.g.: through post alignment clipping
	// Need to do this the hard way, i.e., walk and count.

	const char * seqptr=nullptr;
	int32 needwalk=0;
	if(pcrI.getReadDirection()>0){
	  seqptr=pcrI->getSeqAsChar();
	  if(pcrI->getLeftClipoff()){
	    needwalk=pcrI->getLeftClipoff();
	  }
	}else{
	  seqptr=pcrI->getComplementSeqAsChar();
	  if(pcrI->getLenSeq()-pcrI->getRightClipoff()>0){
	    needwalk=pcrI->getLenSeq()-pcrI->getRightClipoff();
	  }
	}
	const char * endptr=seqptr+pcrI->getLenSeq();

	uint32 scount=0;
	for(; needwalk>0; --needwalk, ++seqptr){
	  if(*seqptr != '*') ++scount;
	}
	if(scount>0) ostr << scount << "S";

	uint32 mcount=0;
	uint32 dcount=0;
	for(uint32 counter=0; counter<pcrI->getLenClippedSeq(); ++counter, ++seqptr){
	  if(*seqptr != '*'){
	    if(dcount){
	      ostr << dcount << "D";
	      dcount=0;
	    }
	    ++mcount;
	  }else{
	    if(mcount){
	      ostr << mcount << "M";
	      mcount=0;
	    }
	    ++dcount;
	  }
	}
	if(dcount){
	  ostr << dcount << "D";
	}
	if(mcount){
	  ostr << mcount << "M";
	}

	// trailing S (if needed)
	// again, cannot use clips, must walk.
	scount=0;
	for(; seqptr!=endptr; ++seqptr){
	  if(*seqptr != '*') ++scount;
	}
	if(scount>0) ostr << scount << "S";
      }

      ostr << '\t' << samc.getRNextEntry(samri)
	   << '\t' << samri.pnext
	   << '\t' << samri.tlen
	   << '\t';

      {
	seqstr.clear();
	qualstr.clear();

	const char * seqptr=pcrI->getSeqAsChar();
	int32 qualindex=0;
	int32 qualincr=1;
	if(pcrI.getReadDirection()<0){
	  seqptr=pcrI->getComplementSeqAsChar();
	  qualindex=pcrI->getLenSeq()-1;
	  qualincr=-1;
	}


	for(uint32 counter=0; counter<pcrI->getLenSeq(); ++counter, ++seqptr, qualindex+=qualincr){
	  if(*seqptr!='*'){
	    seqstr+=*seqptr;
	    qualstr+=static_cast<char>(pcrI->getQualityInSequence(qualindex)+33);
	  }
	}
      }

      ostr << seqstr
	   << '\t' << qualstr
	   << "\tRG:Z:" << pcrI->getReadGroupID().getLibId();

      if(pcrI->getNumOfTags()){
	ostr << "\tPT:Z:";
	bool wantpipe=false;
	for(auto & readtag : pcrI->getTags()){
	  if(wantpipe) ostr << '|';
	  wantpipe=true;
	  if(pcrI.getReadDirection()>0){
	    readtag.dumpAsSAM(ostr);
	  }else{
	    readtag.dumpAsSAM(ostr,pcrI->getLenSeq());
	  }
	}
      }
      ostr << '\n';
    }

    // we're done, except maybe some trailing consensus tags
    for(; ctI!=CON_consensus_tags.end(); ++ctI) ctI->dumpAsSAM(ostr,getContigName());
  }
  catch(Notify n){
    cout << "Oooops, error while writing SAM?\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}



/*************************************************************************
 *
 * I need to count BS lines before output, so this function
 *  is called twice: once for counting only, once for output
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
int32 Contig::priv_helper_dumpAsACE_BSLines(std::string & consseq, bool alsodump, std::ostream & ostr)
{
  FUNCSTART("int32 Contig::priv_helper_dumpAsACE_BSLines(std::string & consseq, bool alsodump, std::ostream & ostr)");

  int32 numBSlines=0;

  CEBUG("consseq: " << consseq << endl);

  const char * cptr = consseq.c_str();
  int32 conpos = 0;
  int32 stretch_begin=1;

  //std::vector<out_order_t>::const_iterator outorderI = CON_outputorder.begin();

  // initialise the iterator for getting through the contig
  rcci_t rcci(this,
	      nullptr,  // all strainids
	      nullptr,  // all seqtypes
	      false,    // no rails
	      true,     // take backbones
	      true);   // take reads without readpool-reads

  auto currentpcrI=rcci.getPCRIsInCol().front();
  for(;*cptr; ++conpos, rcci.advance(), ++cptr){
    char consbase = static_cast<char>(toupper(*cptr));
    base_quality_t bestqual=0;
    auto bqpcrI=rcci.getPCRIsInCol().front();

    CEBUG("caa: conpos: " << conpos << '\t' << consbase << '\n');

    // look for best qual equal to cons base
    bool foundcurrindexincol=false;
    bool foundaread=false;
    //for(uint32 i= 0; i<rcci.read_ids_in_col.size(); i++) {
    for(auto & pcrI : rcci.getPCRIsInCol()){
      CEBUG("caa: read: " << pcrI->getName());
      if(pcrI==currentpcrI && toupper(getBaseInRead(conpos,pcrI)) == consbase){
	foundcurrindexincol=true;
      }
      if(toupper(getBaseInRead(conpos,pcrI)) == consbase
	 || (toupper(getBaseInRead(conpos,pcrI)) == 'X' && consbase=='N')) {
	auto q=getQualityInRead(conpos,pcrI);
	if(q >= bestqual) {
	  bestqual=q;
	  bqpcrI=pcrI;
	  foundaread=true;
	  CEBUG("\tbq");
	}
      }
      CEBUG('\n');
    }

    if(!foundaread){
      bqpcrI=rcci.getPCRIsInCol().front();
      bestqual=getQualityInRead(conpos,bqpcrI);
    }

    CEBUG("bq: " << static_cast<uint16>(bestqual) << "\nbqpcrI: " << bqpcrI->getName() << '\n');

    // look if current read index (if set) is not within 5 of best qual
    //  then output bs line
    if((!foundcurrindexincol || bestqual - getQualityInRead(conpos,currentpcrI) > 80)
      && bqpcrI != currentpcrI) {
      CEBUG("Jump new read: " << bqpcrI->getName() << '\n');
      ++numBSlines;
      if(alsodump){
	ostr << "BS " << stretch_begin << " " << conpos << " ";
	BUGIFTHROW(currentpcrI==CON_reads.end(),"currentpcrI==CON_reads.end() ???");
	ostr << currentpcrI->getName() << '\n';
      }
      stretch_begin=conpos+1;
      currentpcrI=bqpcrI;
    }
  }

  ++numBSlines;
  if(alsodump) {
    ostr << "BS " << stretch_begin << " " << conpos << " ";
    BUGIFTHROW(currentpcrI==CON_reads.end(),"currentpcrI==CON_reads.end() ???");
    ostr << currentpcrI->getName() << "\n\n";
  }

  FUNCEND();
  return numBSlines;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpAsACE(std::ostream & ostr)
{
  FUNCSTART("void Contig::priv_dumpAsACE(std::ostream & ostr)");

  // 0 basesegments, don't know if i use the golden path
  //if(CON_reads.size() > 1){
  //  ostr << "CO C" << CON_id << " " << CON_counts.size() << " ";
  //}else{
  //  ostr << "CO S" << CON_id << " " << CON_counts.size() << " ";
  //}
  //ostr << CON_reads.size() << " 0 U" << '\n';

  std::string consseq;
  std::vector<base_quality_t> consqual;

  finalise();
  newConsensusGet(consseq, consqual);

  ostr << "CO " << getContigName() << " " << CON_counts.size() << " ";

  // TODO: check, isn't there some tracking now by the contig?
  if(!CON_outputrails){
    uint32 numreads=0;
    for(auto read : CON_reads){
      if(!read.isRail()) ++numreads;
    }
    ostr << numreads << " ";
  }else{
    ostr << CON_reads.size() << " ";
  }
  ostr << priv_helper_dumpAsACE_BSLines(consseq,false,ostr) << " U\n";

  {
    const char * cptr = consseq.c_str();
    uint32 cpl=0;
    while(*cptr){
      ostr << *cptr;
      ++cptr;
      if(cpl++==59){
	cpl=0;
	ostr << "\n";
      }
    }
    ostr << "\n\n";
  }

  {
    // qual for pads MUST NOT be output
    ostr << "BQ\n";
    const char * cptr = consseq.c_str();
    uint32 cpl=0;
    //while(I!=consqual.end()){
    for(auto cqI : consqual){
      if(*cptr !='*') {
	ostr << static_cast<uint16>(cqI) << " ";
	if(cpl++==19){
	  cpl=0;
	  ostr << "\n";
	}
      }
      ++cptr;
    }
    ostr << "\n\n";
  }

  // AF lines
  for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
    if(!CON_outputrails && pcrI->isRail()) continue;
    ostr << "AF " << pcrI->getName();
    if(pcrI.getReadDirection() >0) {
      ostr << " U ";
      ostr << (static_cast<int32>(pcrI.getReadStartOffset())+1 - static_cast<int32>(pcrI->getLeftClipoff()) )<< '\n';
    } else {
      ostr << " C ";
      ostr << (static_cast<int32>(pcrI.getReadStartOffset()) + static_cast<int32>(pcrI->getRightClipoff()) - static_cast<int32>(pcrI->getLenSeq()) + 1 )<< '\n';
    }
  }
  ostr << "\n";

  // BS lines
  priv_helper_dumpAsACE_BSLines(consseq,true,ostr);

  // Tags
  for(auto contag : CON_consensus_tags){
    ostr << "CT{\n";
    //ostr << "C" << CON_id;
    ostr << getContigName() << ' ' << contag.getIdentifierStr();
    if(contag.getStrandDirection()==-1) {
      ostr << " MIRA " << (contag.to)+1 << ' ' << (contag.from)+1;
    }else{
      ostr << " MIRA " << (contag.from)+1 << ' ' << (contag.to)+1;
    }
    ostr << " 020202:121212 NoTrans\n";
    if(contag.comment != CON_tagentry_coEmpty){
      ostr << "COMMENT{\n" << contag.getCommentStr() << "\nC}\n";
    }
    ostr << "}\n\n";
  }

  // Reads
  for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
    if(!CON_outputrails && pcrI->isRail()) continue;
    if (pcrI.getReadDirection() >0 ) {
      Read::setCoutType(Read::AS_ACE);
      ostr << *pcrI;
    } else {
      Read::setCoutType(Read::AS_ACE_COMPLEMENT);
      ostr << *pcrI;
    }
  }

  ostr << endl;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::saveAsGAP4DA(const std::string & dirname, std::ostream & fofnstr)
{
  FUNCSTART("void Contig::saveAsGAP4DA()");

  finalise();

  // Reads
  {
    std::string anchor="";
    std::string APline="";
    //for(uint32 cooi=0; cooi<CON_outputorder.size(); cooi++){
    for(auto pcrI=CON_reads.begin(); pcrI != CON_reads.end(); ++pcrI){
      if(!CON_outputrails && pcrI->isRail()) continue;
      if(APline.empty()){
	if(pcrI.getReadDirection() > 0){
	  APline="*new* +";
	} else {
	  APline="*new* -";
	}
	anchor=pcrI->getName();
      } else {
	std::ostringstream tmp;
	if(pcrI.getReadDirection() > 0){
	  tmp << anchor << " + " << pcrI.getReadStartOffset() << " -1";
	}else{
	  tmp << anchor << " - " << pcrI.getReadStartOffset() << " -1";
	}
	APline=tmp.str();
      }

      // we really don't want EXP filenames with a | as part of the filename
      std::string expfilename=pcrI->getName()+".exp";
      for(uint32 ii=0; ii<expfilename.size(); ii++){
	if(expfilename[ii]=='|') expfilename[ii]='_';
      }

      std::ofstream expout((dirname+"/"+expfilename), std::ios::out | std::ios::trunc);
      const_cast<Read &>(*pcrI).dumpAsGAP4DA(expout, APline);


      // We're putting consensus tags all in the first exp file
      // james bonfield suggested the last file, but this lead to unexplainable
      //  off by 1,2,3,... shifts.
      if(pcrI==CON_reads.begin()) {
	for(const auto & contag : CON_consensus_tags){
	  expout << "TC   " << contag.getIdentifierStr();

	  if(pcrI.getReadDirection() > 0) {
	    // if the first read of that contig is saved in + direction,
	    // we must add the left clippoff of the first read that
	    //  was written to the file
	    int32 lq=pcrI->getLeftClipoff();
	    expout << " = " << (contag.from)+1+lq << ".." << (contag.to)+1+lq;
	  } else {
	    // if the first read of that contig is saved complemented,
	    //  save the tag positions also complemented
	    //  (and take right extend)
	    int32 lq=pcrI->getRightClipoff();
	    expout << " = " << lq-static_cast<int32>(contag.to) << ".." << lq-static_cast<int32>(contag.from);
	  }
	  const char * cs=contag.getCommentStr().c_str();
	  if (*cs!=0) expout << "\nTC        ";
	  while(*cs) {
	    if(*cs == '\n') expout << "\nTC        ";
	    expout << *cs;
	    cs++;
	  }
	  expout << '\n';
	}
      }

      expout.close();

      fofnstr << expfilename << '\n';
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

void Contig::dumpContigReadList_Head(std::ostream &ostr)
{
  ostr << "#\n#conName\treadName\n#\n";
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpContigReadList_Body(std::ostream &ostr)
{
  for(auto & cre : CON_reads){
    ostr << getContigName() << "\t";
    ostr << cre.getName() << "\n";
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpContigStatistics_Head(std::ostream &ostr)
{
  ostr << std::setw(20) << std::left << "#name";
  ostr << "\tlength\tav.qual\t#-reads\tmx.cov.\tav.cov\tGC%\tCnIUPAC\tCnFunny\tCnN\tCnX\tCnGap\tCnNoCov\n";
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpContigStatistics_Body(std::ostream &ostr)
{
  Contig::constats_t s=getStats();

  ostr << std::setw(20) << std::left << getContigName();
  ostr << '\t';
  ostr << s.conlength_nogap << '\t';
  ostr << s.avg_conqual << '\t';
  ostr << getNumReadsInContig() << '\t';
  ostr << s.max_coverage << '\t';
  auto f = ostr.flags();
  ostr.setf(std::ios::fixed, std::ios::floatfield);
  ostr.setf(std::ios::showpoint);
  ostr.precision(2);
  ostr << s.avg_coverage << '\t';
  ostr << s.gccontent << '\t';
  ostr.flags(f);

  //ostr << s.avg_readlen << '\t';
  //ostr << "--\t";

  //ostr << s.AinC << '\t';
  //ostr << s.CinC << '\t';
  //ostr << s.GinC << '\t';
  //ostr << s.TinC << '\t';
  ostr << s.IUPACinC << '\t'
       << s.FunnyInC << '\t'
       << s.NinC << '\t'
       << s.XinC << '\t'
       << s.starInC << '\t'
       << s.numnocoverage
       << '\n';
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpReadTagList_Head(std::ostream &ostr)
{
  ostr << "#\n"
       << "#conName\tcFromPadded\tcToPadded\tcFromUnpadded\tcToUnpadded\ttype\trName\trFromPadded\trToPadded\trFromUnpadded\trToUnpadded\tcomment\n"
       << "#\n";
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpReadTagList_Body(std::ostream &ostr)
{
  FUNCSTART("void Contig::dumpReadTagList_Body(std::ostream &ostr)");
  std::string serialc;

  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI) {
    if(!pcrI->isBackbone()) const_cast<Read &>(*pcrI).sortTags();

    for(auto & thetag : pcrI->getTags()){
      //const multitag_t & thetag= cre.getTag(t);
      int32 cpf=pcrI.unclippedReadPos2ContigPos(thetag.from);
      int32 cpt=pcrI.unclippedReadPos2ContigPos(thetag.to);
      ostr << getContigName() << "\t"
	   << cpf << "\t"
	   << cpt << "\t";
      if(cpf>=0) ostr << paddedPos2UnpaddedPos(cpf);
      ostr << "\t";
      if(cpt>=0) ostr << paddedPos2UnpaddedPos(cpt);
      ostr << "\t"
	   << thetag.getIdentifierStr() << "\t"
	   << pcrI->getName() << "\t"
	   << thetag.from << "\t"
	   << thetag.to << "\t"
	   << pcrI->getLowerNonGapAdjustmentPosOfReadPos(thetag.from) << "\t"
	   << pcrI->getUpperNonGapAdjustmentPosOfReadPos(thetag.to) << "\t";
      thetag.serialiseComment(serialc);
      ostr << serialc << "\n";
    }
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpConsensusTagList_Head(std::ostream &ostr)
{
  ostr << "#\n"
    "#conName\tsizePadded\tfromPadded\ttoPadded\tfromUnPadded\ttoUnPadded\tlen\ttype\tSNPqual\tqualA\tqualC\tqualG\tqualT\tqual*\tcomment\n"
    "#" << endl;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpConsensusTagList_Body(std::ostream &ostr)
{
  FUNCSTART("void Contig::dumpConsensusTagList_Body(std::ostream &ostr)");

  std::string serialc;

  // FIXME:
  //  very tricky. paddedPos2UnpaddedPos() needs an existing consensus. If it does not exist,
  //   it will recalculate one
  // But this recalculation may add tags to the consensus!!!
  // Therefore, force this possible recalc to be done before the anything else, or else hilarity ensues
  //  in the for loop below as it gets its elements changed away under its ass.
  paddedPos2UnpaddedPos(0);

  sortConsensusTags();
  for(auto & cte : getConsensusTags()){
    ostr << getContigName() << "\t";
    ostr << getContigLength() << "\t";
    ostr << cte.from << "\t";
    ostr << cte.to << "\t";
    ostr << paddedPos2UnpaddedPos(cte.from) << "\t";
    ostr << paddedPos2UnpaddedPos(cte.to) << "\t";
    ostr << cte.to-cte.from+1 << '\t';
    ostr << cte.getIdentifierStr() << "\t";
    if(cte.additionalinfo_initialised) {
      base_quality_t snpQual=255;
      //if((cte.qualA)>0 && cte.qualA<snpQual) snpQual=cte.qualA;
      //if((cte.qualC)>0 && cte.qualC<snpQual) snpQual=cte.qualC;
      //if((cte.qualG)>0 && cte.qualG<snpQual) snpQual=cte.qualG;
      //if((cte.qualT)>0 && cte.qualT<snpQual) snpQual=cte.qualT;
      //if((cte.qualStar)>0 && cte.qualStar<snpQual) snpQual=cte.qualStar;

      for(uint8 i=0; i<5; i++) {
	if((cte.qualACGTGap[i])>0 && cte.qualACGTGap[i]<snpQual) snpQual=cte.qualACGTGap[i];
      }

      if(snpQual==255) snpQual=0;

      ostr << static_cast<uint16>(snpQual) << "\t";
      for(uint8 i=0; i<5; i++) {
	ostr << static_cast<uint16>(cte.qualACGTGap[i]) << "\t";
      }
      //ostr << static_cast<uint16>(cte.qualACGTGap[0]) << "\t";
      //ostr << static_cast<uint16>(cte.qualACGTGap[1]) << "\t";
      //ostr << static_cast<uint16>(cte.qualACGTGap[2]) << "\t";
      //ostr << static_cast<uint16>(cte.qualACGTGap[3]) << "\t";
      //ostr << static_cast<uint16>(cte.qualACGTGap[4]) << "\t";
    }else{
      ostr << "\t\t\t\t\t\t";
    }
    cte.serialiseComment(serialc);
    ostr << serialc << "\n";
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpMAF_Head(std::ostream &ostr)
{
  ostr << "@Version\t2\t0\n";
  ostr << "@Program\tMIRALIB\n";
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpTCS_Head(std::ostream &ostr)
{
  ostr << "#TCS V1.0\n";
  ostr << "#\n";
  ostr << std::setw(20) << std::left << "#conName"
       << std::setw(9) << std::right << "padPos" << std::setw(9) << "upadPos"
       << " | B" << "  Q |"
       << std::setw(5) << "tcov"
       << std::setw(5) << "covA"
       << std::setw(5) << "covC"
       << std::setw(5) << "covG"
       << std::setw(5) << "covT"
       << std::setw(5) << "cov*"
       << " |"
       << std::setw(3) << "qA"
       << std::setw(3) << "qC"
       << std::setw(3) << "qG"
       << std::setw(3) << "qT"
       << std::setw(3) << "q*"
       << " |"
       << std::setw(3) << "S"
       << " | Tags\n"
       << "#" << endl;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpTCS_Body(std::ostream &ostr)
{
  priv_dumpAsTCS(ostr);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::dumpWiggle_Head(std::ostream &ostr)
{
  ostr << "#TCS V0.9\n"
    "#\n"
    "#conName paddedPos unpaddedPos Base bastyp bqual totcov qualA qualC qualG qualT qual*   covA covC covG covT cov*   Tags\n"
    "#" << endl;
}

/*************************************************************************
 *
 * dumps the coverage of a strain as wiggle file
 * must give padded consensus to routine (save recalculation)
 *
 *************************************************************************/

void Contig::dumpWiggle_Body(std::ostream &ostr, std::string & cons, bool aspadded)
{
  FUNCSTART("void Contig::dumpWiggle_Body(std::ostream &ostr, std::string & cons)");
  BUGIFTHROW(cons.size() != CON_counts.size(), getContigName() << ": consensus.size() != CON_counts.size() ???");

  ccctype_t maxcov=0;
  for(auto & cce : CON_counts){
    maxcov=std::max(maxcov, cce.total_cov);
  }

  uint32 maxstrains=ReadGroupLib::getNumOfStrains();
  std::vector<bool> strainlookup(maxstrains,false);


  ostr << "track type=wiggle_0 name=\"" << getContigName()
       << "cov\" description=\""
       << getContigName() << " ";

  auto crI=CON_reads.begin();
  uint32 numbackbonestrains=0;
  for(;crI != CON_reads.end(); ++crI){
    if(crI.getORPID()>=0 && crI->isBackbone()){
      if(!strainlookup[crI->getStrainID()]){
	strainlookup[crI->getStrainID()]=true;
	if(numbackbonestrains==0){
	  ostr << "BB strain:";
	}
	ostr << ' ' << crI->getStrainName();
	++numbackbonestrains;
      }
    }
  }

  mstd::fill(strainlookup,false);
  uint32 numstrains=0;
  crI=CON_reads.begin();
  for(;crI != CON_reads.end(); ++crI){
    if(crI.getORPID()>=0 && !crI->isBackbone() && !crI->isRail()){
      if(!strainlookup[crI->getStrainID()]){
	strainlookup[crI->getStrainID()]=true;
	if(numstrains==0){
	  if(numbackbonestrains>0){
	    ostr << " Mapped strain(s):";
	  }else{
	    ostr << " Strain(s):";
	  }
	}
	if(numstrains<10){
	  ostr << ' ' << crI->getStrainName();
	} else if(numstrains==10){
	  ostr << "...";
	}
	++numstrains;
      }
    }
  }

  uint32 span=4;

  BUGIFTHROW(span==0, "span for wiggle file cannot be 0");

  ostr << "\" visibility=full autoScale=off "
       << "viewLimits=0:" << maxcov << " color=0,200,100 "
       << "maxHeightPixels=100:50:20 graphType=bar priority=30\n"
       << "fixedStep chrom=" << getContigName()
       << " start=1 step=" << span << " span=" << span << endl;

  auto ccI=CON_counts.begin();
  for(uint32 cpos=0; cpos<cons.size();){
    uint32 tcov=0;
    uint32 counter=0;
    bool haszerovalue=false;
    for(; counter<span && ccI!=CON_counts.end(); cpos++, ccI++){
      if(aspadded
	 || (cons[cpos] != '*'
	     && (cons[cpos] != 'X'))){ // || (cons[cpos] == 'X' && ccI->total_cov!=0))){
	if(ccI->total_cov){
	  tcov+=ccI->total_cov;
	}else{
	  haszerovalue=true;
	}
	counter++;
      }
    }
    if(haszerovalue || counter==0){
      ostr << "0\n";
    }else if(counter!=0) {
      ostr << tcov/counter << '\n';
    }
  }

  FUNCEND();
}



/*************************************************************************
 *
 * dumps the GC content of a strain as wiggle file
 * must give padded consensus to routine (save recalculation)
 *
 *************************************************************************/

void Contig::dumpGCWiggle_Body(std::ostream &ostr, std::string & cons)
{
  FUNCSTART("void Contig::dumpGCWiggle_Body(std::ostream &ostr, std::string & cons)");
  BUGIFTHROW(cons.size() != CON_counts.size(), getContigName() << ": consensus.size() != CON_counts.size() ???");

  uint32 span=4;

  BUGIFTHROW(span==0, "span for wiggle file cannot be 0");

  ostr << "track type=wiggle_0 name=\"" << getContigName()
       << "gccont\" description=\""
       << getContigName() << " GC content"
       << "\" visibility=full autoScale=off "
       << "viewLimits=0:100 color=0,100,200 "
       << "maxHeightPixels=100:50:20 graphType=line priority=30\n"
       << "fixedStep chrom=" << getContigName()
       << " start=1 step=" << span << " span=" << span << endl;

  int32 window=50;
  int32 topos=window/2;
  int32 frompos=-topos;

  std::vector<int32> counts(255,0);

  for(int32 i=frompos; i<topos; ++i){
    if(i>=0 && i<cons.size()){
      ++counts[toupper(cons[i])];
    }
  }

  for(int32 cpos=0; cpos<cons.size();){
    uint32 counter=0;
    for(; counter<span && cpos < cons.size(); cpos++){
      if(cons[cpos] != '*'
	 && (cons[cpos] != 'X')){ // || (cons[cpos] == 'X' && ccI->total_cov!=0))){
	if(frompos>=0) --counts[toupper(cons[frompos])];
	++frompos;
	++topos;
	if(topos < cons.size()) ++counts[toupper(cons[topos])];
	++counter;
      }
    }
    uint32 allcounts=counts['A']+counts['C']+counts['G']+counts['T'];
    if(allcounts>0){
      double gc=100*static_cast<double>(counts['C']+counts['G'])/static_cast<double>(allcounts);
      ostr << static_cast<uint16>(gc) << "\n";
    }else{
      ostr << "0\n";
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

void Contig::dumpTagsAsGFF3(std::ostream & ostr, std::string & padcons)
{
  FUNCSTART("void Contig::dumpTagsAsGFF3(std::ostream & ostr, std::string & padcons)");

  static const std::string namestr="Name"; // for comparison later
  static const std::string idstr="ID"; // for comparison later

  //  build a mapping for positions from padded sequence to unpadded sequence
  std::vector<uint32> depadposmap;
  if(padcons.empty()){
    cout << "no padcons\n";
  }else{
    BUGIFTHROW(padcons.size() != CON_counts.size(),"padcons.size() != CON_counts.size() ???");
    cout << "havepadcons " << padcons.size() << endl;
    uint32 actdepadpos=0;
    depadposmap.resize(padcons.size(),0);
    for(uint32 seqi=0; seqi < padcons.size();seqi++){
      depadposmap[seqi]=actdepadpos;
      char actbase=static_cast<char>(tolower(padcons[seqi]));
      if(actbase!='*') actdepadpos++;
    }
  }

  std::deque<multitag_t> tagswithconpos;

  // dump tags from backbone reads ( == initial annotation)
  multitag_t tmptag;
  std::string tmps;
  for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
    if(CON_reads.size()==1 || pcrI->isBackbone()){
      for(uint32 ti=0; ti < pcrI->getNumOfTags(); ++ti){
	if(pcrI->getTag(ti).identifier == Read::REA_tagentry_idSRMr
	   || pcrI->getTag(ti).identifier == Read::REA_tagentry_idWRMr
	   || pcrI->getTag(ti).identifier == Read::REA_tagentry_idSROr
	   || pcrI->getTag(ti).identifier == Read::REA_tagentry_idSIOr
	   || pcrI->getTag(ti).identifier == Read::REA_tagentry_idSAOr
	   || pcrI->getTag(ti).identifier == Read::REA_tagentry_idUNSr
	   || pcrI->getTag(ti).identifier == Read::REA_tagentry_idESDN
	  ){
	  // output them or not? at the moment rather not
	}else{
	  if(pcrI->getTag(ti).identifier == Read::REA_tagentry_idSOFAcontig){
	    // need to rewrite "Name=..." to match the one of the contig
	    gff3attributes_t tmpa;
	    GFFParse::parseGFF3Attributes(pcrI->getTag(ti).getCommentStr(),
					  tmpa);
	    for(auto & ae : tmpa){
	      if(ae.tag==namestr
		 || ae.tag==idstr){
		ae.values.resize(1);
		ae.values[0]=getContigName();
	      }
	    }
	    GFFParse::createGFF3AttributeString(tmpa, tmps);
	    tmptag=pcrI->getTag(ti);
	    tmptag.setCommentStr(tmps);

	    tagswithconpos.push_back(tmptag);
	  }else{
	    tagswithconpos.push_back(pcrI->getTag(ti));
	  }
	  // the tags contain read positions. We need to convert these into contig positions,
	  //  then into depadded sequence positions.
	  // However, there may be tags (starting/ending/completely) outside the good range of the contig
	  // resolution:
	  //  completely outside -> drop
	  //  partly inside -> adjust start/stop
	  tagswithconpos.back().from=pcrI.unclippedReadPos2ContigPos(tagswithconpos.back().from);
	  tagswithconpos.back().to=pcrI.unclippedReadPos2ContigPos(tagswithconpos.back().to);

	  if(!depadposmap.empty()){
	    // completely outside?
	    if(tagswithconpos.back().from >= depadposmap.size() &&
	       tagswithconpos.back().to >= depadposmap.size()) continue;

	    // partly?
	    if(tagswithconpos.back().from >= depadposmap.size()) tagswithconpos.back().from = 0;
	    tagswithconpos.back().from=depadposmap[tagswithconpos.back().from];
	    if(tagswithconpos.back().to >= depadposmap.size()) tagswithconpos.back().to = depadposmap.size()-1;
	    tagswithconpos.back().to=depadposmap[tagswithconpos.back().to];
	  }
	}
      }
    }
  }

  // dump consensus tags ( == SNP and other markers )

  {
    for(auto & ct : CON_consensus_tags){
      tmptag=ct;
      bool saveit=false;
      if(tmptag.identifier==Contig::CON_tagentry_idSROc){
	tmptag.identifier=Read::REA_tagentry_idSROm;
	saveit=true;
      }else if(tmptag.identifier==Contig::CON_tagentry_idSIOc){
	tmptag.identifier=Read::REA_tagentry_idSIOm;
	saveit=true;
      }else if(tmptag.identifier==Contig::CON_tagentry_idSAOc){
	tmptag.identifier=Read::REA_tagentry_idSAOm;
	saveit=true;
      }else if(tmptag.identifier==Contig::CON_tagentry_idMCVc){
	tmptag.identifier=Read::REA_tagentry_idMCVm;
	saveit=true;
      }else if(tmptag.identifier==Contig::CON_tagentry_idSRMc){
	tmptag.identifier=Read::REA_tagentry_idSRMm;
	saveit=true;
      }else if(tmptag.identifier==Contig::CON_tagentry_idWRMc){
	tmptag.identifier=Read::REA_tagentry_idWRMm;
	saveit=true;
      }
      // hmmm, the above is a bit restrictive
      // let's save everything and see whether problems occur
      saveit=true;
      if(saveit) {
	if(!depadposmap.empty()) {
	  tmptag.from=depadposmap[tmptag.from];
	  tmptag.to=depadposmap[tmptag.to];
	}
	tagswithconpos.push_back(tmptag);
      }
    }
  }

  for(auto & tag : tagswithconpos){
    tag.dumpAsGFF3(ostr,getContigName());
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::priv_dumpReplay(std::ofstream & eout, const AlignedDualSeqFacts * initialadsf, int32 refid, int32 newid, int32 direction_frnid, bool newid_ismulticopy, int32 forcegrow)
{
  Read::setCoutType(Read::AS_MAF);
  setCoutType(AS_MAF);
  setOutputRails(true);

  eout << "This file allows one to replay an error that MIRA just encountered\n\n"
    "Please do not delete but contact the author (bach@chevreux.org) immediately\n\n"
    "addRead_Error!\n\n\n\n";
  if(initialadsf!=nullptr) eout << "\ninitialadsf:\n" << *initialadsf;
  eout << "\n\nrefid: " << refid << " (" << CON_readpool->getRead(refid).getName()
       << ")\nnewid: " << newid << " (" << CON_readpool->getRead(newid).getName()
       << ")\ndirection_frnid: " << direction_frnid
       << "\nnewid_ismulticopy: " << newid_ismulticopy
       << "\nforcegrow: " << forcegrow;

  if(initialadsf!=nullptr){
    eout << "\n\nint8 dirnewid=" << direction_frnid
	 << ";  // " << static_cast<int16>(initialadsf->getSequenceDirection(newid))
	 << " " << static_cast<int16>(initialadsf->getSequenceDirection(refid))
	 << "\nadsf.publicinit(" << initialadsf->getID1()
	 << "," << initialadsf->getID2()
	 << "," << initialadsf->getDelta()
	 << "," << initialadsf->getRightDelta(initialadsf->getID1())
	 << "," << initialadsf->getRightDelta(initialadsf->getID2())
	 << "," << initialadsf->getTotalLen()
	 << "," << static_cast<int16>(initialadsf->getSequenceDirection(initialadsf->getID1()))
	 << "," << static_cast<int16>(initialadsf->getSequenceDirection(initialadsf->getID2()))
	 << "," << static_cast<uint16>(initialadsf->getScoreRatio())
	 << "," << initialadsf->getTotalNonMatches()
	 << "," << initialadsf->get5pLenContiguousMatch(initialadsf->getID1())
	 << "," << initialadsf->get3pLenContiguousMatch(initialadsf->getID1())
	 << "," << initialadsf->get5pLenContiguousMatch(initialadsf->getID2())
	 << "," << initialadsf->get3pLenContiguousMatch(initialadsf->getID2())
	 << ");\n";
  }

  eout << "readgroups known:\n";
  ReadGroupLib::dumpAllReadGroupsAsMAF(eout);

  eout << "\n\nOffending read:\n" << CON_readpool->getRead(newid)
       << endl
       << "\nOffending contig:\n"
       << "@Version\t2\t0\n" << *this
       << endl;
}


/*************************************************************************
 *
 * Dump some internal things on the contig
 *
 * note: incomplete
 *
 *************************************************************************/

void Contig::dumpStatus(std::ostream & ostr)
{
  ostr << "ContigDump: " << CON_id << " " << CON_nameprefix << " " << CON_name << endl
       << "CON_finalised: " << CON_finalised << endl;
  if(CON_readpool==nullptr){
    ostr << "CON_readpool: nullptr !\n";
  }else{
    ostr << "CON_readpool: set\n";
  }
  ostr << "CON_reads: " << CON_reads.size() << endl
       << "CON_counts: " << CON_counts.size() << endl
       << "CON_templates_present: " << CON_templates_present.size() << endl
       << "CON_consensus_tags: " << CON_consensus_tags.size() << endl
       << "CON_targetcoverageperst: " << CON_targetcoverageperst.size();
  for(uint32 i=0; i<CON_targetcoverageperst.size();++i){
    ostr << "\t" << CON_targetcoverageperst[i] << endl;
  }

  ostr << "CON_allowedrefids: " << CON_allowedrefids.size() << endl
       << "CON_2tmpcons: " << CON_2tmpcons.size() << endl
       << "CON_tmpcons_from_backbone: " <<  CON_tmpcons_from_backbone << endl
       << "CON_specialsraddconditions: " << CON_specialsraddconditions << endl
       << "CON_ssrc_maxtotalerrors: " << CON_ssrc_maxtotalerrors << endl
       << "CON_ssrc_maxgaps: " << CON_ssrc_maxgaps << endl
       << "CON_ssrc_maxmismatches: " << CON_ssrc_maxmismatches << endl
       << "num backbones: " << getNumBackbones() << endl;

  for(uint32 i=0; i<NUMMERGESEQTYPES; ++i){
    ostr << "CON_nummergedreads_perseqtype[" << i << "]: " << CON_nummergedreads_perseqtype[i] << endl;
  }

  ostr << "CON_fixedconsseq: " << CON_fixedconsseq.size() << endl
       << "CON_fixedconsqual: " << CON_fixedconsqual.size() << endl
       << "CON_conscalc_mincov: " << CON_conscalc_mincov << endl
       << "CON_strainconsqual: " <<  CON_strainconsqual.size() << endl;

  for(uint32 i=0; i<CON_strainconsqual.size(); ++i){
    ostr << "\t" << i << " " << CON_strainconsqual[i].size() << endl;
  }

  ostr << "CON_strainadjustments: " << CON_strainadjustments.size() << endl;
  for(uint32 i=0; i<CON_strainadjustments.size(); ++i){
    ostr << "\t" << i << " " << CON_strainadjustments[i].size() << endl;
  }

  ostr << "CON_readsperstrain: " << CON_readsperstrain.size() << endl;
  for(uint32 i=0; i<CON_readsperstrain.size(); ++i){
    ostr << "\t" << i << " " << ReadGroupLib::getStrainOfStrainID(i) << "\t" << CON_readsperstrain[i] << endl;
  }

  ostr << "CON_readsperreadgroup: " << CON_readsperreadgroup.size() << endl;
  for(uint32 i=0; i<CON_readsperreadgroup.size(); ++i){
    ostr << "\t" << i << " " << ReadGroupLib::getStrainOfStrainID(i) << "\t" << CON_readsperreadgroup[i] << endl;
  }

  ostr << "CON_last_dangerous_overlaps: " << CON_last_dangerous_overlaps.size() << endl
       << "CON_contains_long_repeats_only: " << CON_contains_long_repeats_only << endl;

}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::debugDump()
{
  FUNCSTART("void Contig::debugDump()");

  cout << "Contig::debugDump()" << endl;
  uint32 cci=0;
  for(auto & cce : CON_counts){
    cout << "cci: " << cci << cce << '\n';
    ++cci;
  }

  FUNCEND();
}
