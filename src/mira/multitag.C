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


#include <io/annotationmappings.H>
#include <mira/multitag.H>
#include <mira/read.H>
#include <mira/contig.H>
#include <mira/gff_parse.H>
#include <util/misc.H>


using std::cout;
using std::cerr;
using std::endl;


/*************************************************************************
 *
 * trashes/releases all string containers of multitag
 *
 *************************************************************************/

void multitag_t::trashContainers()
{
  MT_sc_mttagsource.trash();
  MT_sc_mtidentifier.trash();
  MT_sc_mtcomment.trash();
}

/*************************************************************************
 *
 * GFF3
 *
 *
 *************************************************************************/

void multitag_t::dumpAsGFF3(std::ostream & ostr, const char * seqid) const
{
  FUNCSTART("void multitag_t::dumpAsGFF3(std::ostream & ostr, const char * seqid) const");

  std::string gff3mnemonic(getIdentifierStr());
  bool addmiraitag=false;
  if(!AnnotationMappings::isValidGFF3SOEntry(gff3mnemonic)) {
    gff3mnemonic="experimental_feature";
    addmiraitag=true;
  }
  ostr << seqid;
  if(getSourceStr().empty()){
    ostr << "\tunknown";
  }else{
    ostr << '\t' << getSourceStr();
  }
//	     << '\t' << getIdentifierStr()

  uint32 ofrom=from;
  uint32 oto=to;
//  if(!poslookup.empty()){
//    BUGIFTHROW(ofrom >= poslookup.size(),"ofrom >= poslookup.size() ???");
//    BUGIFTHROW(oto >= poslookup.size(),"oto >= poslookup.size() ???");
//    ofrom=poslookup[ofrom];
//    oto=poslookup[oto];
//  }
  char ostrand=getStrand();
  switch(ostrand){
  case '=' : {
    ostrand='.';
    break;
  }
  case '+' :
  case '-' :
  case '?' : {
    break;
  }
  default : {
    MIRANOTIFY(Notify::FATAL,"unknown strand character " << ostrand << " (" << static_cast<uint16>(ostrand) << ")");
  }
  }
  ostr << '\t' << gff3mnemonic
       << '\t' << ofrom+1
       << '\t' << oto+1
       << "\t."                           // TODO: get score out of comment string and store here
       <<  "\t" << ostrand;
  if(phase==3){
    ostr << "\t.";
  }else{
    ostr << "\t" << static_cast<uint16>(phase);
  }
  ostr << "\t";
  bool needsemicolon=false;
  if(!getCommentStr().empty()){
    if(!commentisgff3){
      std::string tmpstr;
      gff3Code(getCommentStr(),tmpstr);
      ostr << "Note="<<tmpstr;
      needsemicolon=true;
    }else{
      ostr << getCommentStr();  // TODO: without score!
      needsemicolon=true;
    }
  }
  if(addmiraitag){
    if(needsemicolon) ostr << ";";
    ostr << "miraitag=" << getIdentifierStr();
    needsemicolon=true;
  }
  ostr << '\n';
  FUNCEND();
}


/*************************************************************************
 *
 * CAF
 *
 *
 *************************************************************************/

void multitag_t::dumpAsCAF(std::ostream & ostr) const
{
  std::string xgap4(AnnotationMappings::translateSOfeat2GAP4feat(getIdentifierStr()));
  if(xgap4.empty()){
    xgap4=AnnotationMappings::translateSOfeat2XGAP4feat(getIdentifierStr());
    if(xgap4.empty()){
      xgap4=getIdentifierStr();
    }
  }

  if(getStrandDirection()==-1) {
    ostr << "Tag " << xgap4 << ' ' << (to)+1 << ' ' << (from)+1;
  }else{
    ostr << "Tag " << xgap4 << ' ' << (from)+1 << ' ' << (to)+1;
  }

  ostr << " \"";
  dumpCommentAsAttributes(ostr);
  ostr << "\"\n";
}


void multitag_t::dumpDebug(std::ostream & ostr) const
{
  ostr << "From: " << from <<  endl;
  ostr << "To: " << to << endl;
  ostr << "Strand: " << getStrand() << endl;
  ostr << "Phase: " << static_cast<uint16>(phase) << endl;
  ostr << "Identifier: (" << static_cast<uint64>(identifier.getSCID());
  ostr.flush();
  ostr << ")\t" << getIdentifierStr() << endl;
  ostr << "Comment: (" << static_cast<uint64>(comment.getSCID());
  ostr.flush();
  ostr << ")\t" << getCommentStr() << endl;
  ostr << "Source: (" << static_cast<uint64>(source.getSCID());
  ostr.flush();
  ostr << ")\t" << getSourceStr() << endl;
  ostr << "Comm is GFF3: " << commentisgff3 << endl;
}



void multitag_t::dumpAsMAF(std::ostream & ostr, const char * type) const
{
  ostr << type << "\t" << getIdentifierStr() << '\t' << (from)+1 << '\t' << (to)+1
       << '\t' << getStrand()
       << '\t' << getSourceStr();
  if(phase==3){
    ostr << "\t.";
  }else{
    ostr << '\t' << static_cast<uint16>(phase);
  }

  if(getCommentStr().empty()){
    ostr << '\n';
  }else{
    if(!commentisgff3){
      std::string tmpstr;
      gff3Code(getCommentStr(),tmpstr);
      ostr << "\tNote="<<tmpstr << '\n';
    }else{
      ostr << '\t' << getCommentStr() << '\n';
    }
  }
}


/*************************************************************************
 *
 * dumps as SAM contig tag line
 *
 *
 *************************************************************************/

void multitag_t::dumpAsSAM(std::ostream & ostr, const std::string & contigname) const
{
  std::string xgap4(AnnotationMappings::translateSOfeat2GAP4feat(getIdentifierStr()));
  if(xgap4.empty()){
    xgap4=AnnotationMappings::translateSOfeat2XGAP4feat(getIdentifierStr());
    if(xgap4.empty()){
      xgap4=getIdentifierStr();
    }
  }
  ostr << "*\t768\t" << contigname
       << '\t' << from+1
       << "\t255"
       << '\t' << to+1-from << "M\t*\t0\t0\t*\t*\tCT:Z:";
  if(getStrand()=='='){
    ostr << '.';
  }else{
    ostr << getStrand();
  }
  ostr << ";" << xgap4
       << ";";

  if(!getCommentStr().empty()){
    if(!commentisgff3){
      std::string tmpstr;
      gff3Code(getCommentStr(),tmpstr);
      ostr << "Note="<<tmpstr;
    }else{
      ostr << getCommentStr();
    }
  }
  ostr << '\n';
}


/*************************************************************************
 *
 * dumps as SAM read tag
 *
 * rlen>0 means to dump tags in reverse direction, the positions
 *  then being rlen-pos
 *
 *************************************************************************/

void multitag_t::dumpAsSAM(std::ostream & ostr, int32 rlen) const
{
  if(rlen==0){
    ostr << from+1 << ';' << to+1 << ';';
  }else{
    ostr << rlen-to << ';' << rlen-from << ';';
  }

  if(getStrand()=='='){
    ostr << '.';
  }else{
    ostr << getStrand();
  }
  std::string xgap4(AnnotationMappings::translateSOfeat2GAP4feat(getIdentifierStr()));
  if(xgap4.empty()){
    xgap4=AnnotationMappings::translateSOfeat2XGAP4feat(getIdentifierStr());
    if(xgap4.empty()){
      xgap4=getIdentifierStr();
    }
  }
  ostr << ";" << xgap4
       << ";";

  if(!getCommentStr().empty()){
    if(!commentisgff3){
      std::string tmpstr;
      gff3Code(getCommentStr(),tmpstr);
      ostr << "Note="<<tmpstr;
    }else{
      ostr << getCommentStr();
    }
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void multitag_t::dumpCommentAsAttributes(std::ostream & ostr) const
{
  bool needsemicolon=false;
  if(!getCommentStr().empty()){
    if(!commentisgff3){
      std::string tmpstr;
      gff3Code(getCommentStr(),tmpstr);
      ostr << "Note="<<tmpstr;
    }else{
      ostr << getCommentStr();
    }
    needsemicolon=true;
  }

  if(getStrand()!='='){
    if(needsemicolon) ostr << ";";
    ostr << "gff3str=" << getStrand();
    needsemicolon=true;
  }else if(from != to){
    if(needsemicolon) ostr << ";";
    ostr << "gff3str=.";
    needsemicolon=true;
  }
  if(phase!=3){
    if(needsemicolon) ostr << ";";
    ostr << "gff3pha=" << static_cast<uint16>(phase);
    needsemicolon=true;
  }
  if(!getSourceStr().empty()){
    if(needsemicolon) ostr << ";";
    ostr << "gff3src=" << getSourceStr();
    needsemicolon=true;
  }
}
