/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ) and Bastien Chevreux
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
 * Routines fro simple handling of EXP files
 *
 */


#include "tokens.h"
#include "exp.H"

#include <cstdlib>     // atoi() under Cygwin
#include <memory>

#define CEBUG(bla)

using std::cout;
using std::cerr;
using std::endl;


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void tag_t::serialiseComment(const std::string & comment, std::string & result)
{
  result.clear();
  result.reserve(comment.size()+40);
  for(uint32 i=0; i< comment.size(); i++){
    if(comment[i]!='\n'){
      if(comment[i]=='\t'){
	result+=' ';
      }else{
	result+=comment[i];
      }
    }else{
      result+=" :: ";
    }
  }
  return;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void tag_t::serialiseComment(std::string & result) const
{
  serialiseComment(comment,result);
  return;
}



/*************************************************************************
 *
 * works for serialised comments (where the keys are separated by " :: "
 *  as well as for keys that are in line by line setup
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool tag_t::extractGenBankKeyValueFromComment(const std::string & whatkey, std::string & result) const
{
  return extractGenBankKeyValueFromComment(comment, whatkey, result);
}

bool tag_t::extractGenBankKeyValueFromComment(const std::string & comment, const std::string & whatkey, std::string & result)
{
  CEBUG(endl << whatkey << '\t' << comment << endl);
  result.clear();
  if(comment.empty()) return false;
  auto fpos=comment.find(whatkey,0);
  CEBUG("1: " << fpos << endl);
  if(fpos==std::string::npos) return false;
  CEBUG("2: " << fpos << endl);
  fpos=comment.find_first_not_of(" =\t\"",fpos+whatkey.size());
  CEBUG("3: " << fpos << endl);
  if(fpos==std::string::npos) return false;

  CEBUG("4: " << fpos << endl);
  auto tpos=comment.find(" :: ",fpos);
  if(tpos==std::string::npos){
    CEBUG("5: " << tpos << endl);
    tpos=comment.find_first_of("\n\"",fpos);
    CEBUG("6: " << tpos << endl);
    if(tpos==std::string::npos){
      CEBUG("7: " << tpos << endl);
      tpos=comment.size()-1;
    }
  }
  CEBUG("8: " << tpos << endl);
  tpos=comment.find_last_not_of(" =\t\"",tpos);
  CEBUG("9: " << tpos << endl);
  if(tpos==std::string::npos
    || tpos<fpos) return false;

  CEBUG("10: " << tpos << endl);
  result=comment.substr(fpos, tpos-fpos+1);
  return true;
}
//#define CEBUG(bla)



/*************************************************************************
 *************************************************************************
 *************************************************************************
 *************************************************************************
 *************************************************************************
 *************************************************************************
 *************************************************************************
 *************************************************************************/




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
EXP::EXP()
{
  FUNCSTART("EXP::EXP()");

  // 2000 is a good upper bound for sequence lengths
  //  should sequences be longer, the STL-behind-the-scene doubling
  //  will take care of it.
  EXP_sq.reserve(2000);

  zeroVars();

  FUNCEND();
}

EXP::~EXP()
{
  FUNCSTART("EXP::~EXP()");

  discard();

  FUNCEND();
}


void EXP::discard()
{
  FUNCSTART("void EXP::discard()");

  EXP_qualities.clear();
  EXP_tags.clear();
  EXP_on_ranges.clear();

  zeroVars();

  FUNCEND();
}



void EXP::zeroVars()
{
  EXP_sq.clear();
  EXP_id.clear();
  EXP_le.clear();
  EXP_cn.clear();
  EXP_cf.clear();
  EXP_dt.clear();
  EXP_mt.clear();
  EXP_li.clear();
  EXP_ma.clear();
  EXP_mc.clear();
  EXP_mn.clear();
  EXP_en.clear();
  EXP_ln.clear();
  EXP_lt.clear();
  EXP_ps.clear();
  EXP_sf.clear();
  EXP_sq.clear();
  EXP_cv.clear();
  EXP_sv.clear();
  EXP_ss.clear();
  EXP_tn.clear();
  EXP_op.clear();
  EXP_bc.clear();
  EXP_pn.clear();

  EXP_aq=-1;
  EXP_pr=-1;
  EXP_sc=-1;
  EXP_sp=-1;
  EXP_sl=-1;
  EXP_sr=-1;
  EXP_ql=-1;
  EXP_qr=-1;
  EXP_cl=-1;
  EXP_cr=-1;
  EXP_ch=-1;
  EXP_si_from=-1;
  EXP_si_to=-1;

  EXP_len_seq=0;
}


//char * EXP::conditionalStrCpy(const char * src)
//{
//  if(src == nullptr) return nullptr;
//  uint32 len=strlen(src);
//  char * dest= new char[len+1];
//  strcpy(dest, src);
//  return dest;
//}


int32 EXP::gimmeAnInt(FlexLexer * lexer, const char * filename)
{
  FUNCSTART("int32 EXP::gimmeAnInt(FlexLexer * lexer)");

  int yyretcode=lexer->yylex();
  int pm=1;

  CEBUG(lexer->YYText());
  if(yyretcode==EXPT_PLUS){
    yyretcode=lexer->yylex();
  }
  if(yyretcode==EXPT_MINUS){
    yyretcode=lexer->yylex();
    pm=-1;
  }
  if(yyretcode==EXPT_INT || yyretcode==EXPT_FLOAT){
  }else{
    cerr << "yyret:" << yyretcode << "\t";
    cerr << lexer->YYText();
    MIRANOTIFY(Notify::FATAL, "expected a number: " << filename);
  }

  FUNCEND();
  return pm*atoi(lexer->YYText());
}


// once the SQ has been recognized, get the sequence
//  copies the sequence without \n and blanks into a 0 terminated
//  string
void EXP::getSequence(FlexLexer * lexer, const char * filename)
{
  FUNCSTART("int32 EXP::gimmeSequence(FlexLexer * lexer)");

  int yyretcode=lexer->yylex();

  while(yyretcode == EXPT_SQseq) {

    const char * pos=lexer->YYText();
    for(;;pos++){
      if(!*pos) break;
      if(*pos=='/') break;
      if(*pos==' ') continue;
      if(*pos=='\n') continue;
      EXP_sq.push_back(*pos);
    }
    yyretcode=lexer->yylex();
  }

  if(yyretcode!=EXPT_EOSQ){
    if(yyretcode==EXPT_ILLBASEINSQ){
      cerr << lexer->YYText() << endl;
      MIRANOTIFY(Notify::FATAL, "Illegal base SQ: " << filename);
    }
    MIRANOTIFY(Notify::FATAL, "error while reading SQ: " << filename);
  }

  EXP_len_seq=EXP_sq.size();

  FUNCEND();
}

void EXP::load(const char * filename, bool warning)
{
  FUNCSTART("void EXP::load()");

  std::ifstream fin(filename, std::ios::in|std::ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::WARNING, "EXP file not found for loading: " << filename);
  }
  if(!fin.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length EXP file: " << filename);
  }
  fin.seekg(0, std::ios::beg);

  //FlexLexer* lexer = new EXPFlexLexer(&fin);
  std::unique_ptr<FlexLexer> lexer(new EXPFlexLexer(&fin));


  int yyretcode=-1;
  bool dontread=false;
  while(yyretcode!=0){
    if(dontread==false) yyretcode=lexer->yylex();
    dontread=false;
    switch(yyretcode){
    case 0: break;                              // do nothing, eof
    case EXPT_EOL: break;                       // do nothing, eol that isn't handled otherwise
    case EXPT_UNKNOWN:{
      if(warning){
	cout << "WARNING (in " << filename << "): unknown identifier " << lexer->YYText() << ", assuming oneliner with content:  ";
      }
      lexer->yylex();
      if(warning){
	cout << lexer->YYText() << endl;
      }
      break;
    }

    case EXPT_CF: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_cf=lexer->YYText(); break;}
    case EXPT_CN: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_cn=lexer->YYText(); break;}
    case EXPT_CV: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_cv=lexer->YYText(); break;}
    case EXPT_DT: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_dt=lexer->YYText(); break;}
    case EXPT_EN: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_en=lexer->YYText(); break;}
    case EXPT_ID: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_id=lexer->YYText(); break;}
    case EXPT_LE: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_le=lexer->YYText(); break;}
    case EXPT_LI: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_li=lexer->YYText(); break;}
    case EXPT_LN: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_ln=lexer->YYText(); break;}
    case EXPT_LT: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_lt=lexer->YYText(); break;}
    case EXPT_MA: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_ma=lexer->YYText(); break;}
    case EXPT_MC: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_mc=lexer->YYText(); break;}
    case EXPT_MN: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_mn=lexer->YYText(); break;}
    case EXPT_MT: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_mt=lexer->YYText(); break;}
    case EXPT_PS: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_ps=lexer->YYText(); break;}
    case EXPT_SF:{
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_sf=lexer->YYText(); break;}
    case EXPT_SS: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_ss=lexer->YYText(); break;}
    case EXPT_SV: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_sv=lexer->YYText(); break;}
    case EXPT_TN: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_tn=lexer->YYText(); break;}
    case EXPT_OP: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_op=lexer->YYText(); break;}
    case EXPT_BC: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_bc=lexer->YYText(); break;}
    case EXPT_PN: {
      lexer->yylex(); CEBUG(lexer->YYText()); EXP_pn=lexer->YYText(); break;}


    case EXPT_CL: { EXP_cl=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_CR: { EXP_cr=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_QL: { EXP_ql=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_QR: { EXP_qr=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_SL: { EXP_sl=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_SR: { EXP_sr=gimmeAnInt(lexer.get(), filename); break;}

    case EXPT_AQ: { EXP_aq=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_CH: { EXP_ch=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_PR: { EXP_pr=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_SC: { EXP_sc=gimmeAnInt(lexer.get(), filename); break;}
    case EXPT_SP: { EXP_sp=gimmeAnInt(lexer.get(), filename); break;}

    case EXPT_SQ: { getSequence(lexer.get(), filename); break;}

    case EXPT_AV:{
      CEBUG("Found AV:");
      while((yyretcode=lexer->yylex())==EXPT_INT){
	EXP_qualities.push_back(atoi(lexer->YYText()));
      }
      break;
    }
    case EXPT_CS:{
      CEBUG("Found CS:");
      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_INT){
	//	cout << lexer->YYText() << yyretcode;
	MIRANOTIFY(Notify::FATAL, "expected int after CS: " << filename);
      }
      EXP_cl=atoi(lexer->YYText());
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_PP){
	MIRANOTIFY(Notify::FATAL, "expected .. after CS int: " << filename);
      }
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_INT){
	MIRANOTIFY(Notify::FATAL, "expected int after CS int..: " << filename);
      }
      EXP_cr=atoi(lexer->YYText());
      CEBUG(lexer->YYText() << "found");
      break;
    }
    case EXPT_SI:{
      CEBUG("Found SI:");
      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_INT){
	//	cout << lexer->YYText() << yyretcode;
	MIRANOTIFY(Notify::FATAL, "expected int after SI: " << filename);
      }
      EXP_si_from=atoi(lexer->YYText());
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_PP){
	MIRANOTIFY(Notify::FATAL, "expected .. after SI int: " << filename);
      }
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_INT){
	MIRANOTIFY(Notify::FATAL, "expected int after SI int..: " << filename);
      }
      EXP_si_to=atoi(lexer->YYText());
      CEBUG(lexer->YYText() << "found");
      break;
    }
    case EXPT_TG:{

      const char * cptr;
      tag_t  tmptag;
      EXP_tags.push_back(tmptag);
      tag_t & newtag=EXP_tags.back();

      CEBUG(EXP_tags.size());

      CEBUG("Found TG: ");
      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_ANID){
	cout << lexer->YYText() << yyretcode;
	MIRANOTIFY(Notify::FATAL, "expected alphanum id after TG: " << filename);
      }
      cptr=lexer->YYText();
      while(*cptr) {
	//	cout << *cptr;cout.flush();
       	if(*cptr!='\n') newtag.identifier+=*cptr;
	cptr++;
      }
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_PME){
	MIRANOTIFY(Notify::FATAL, "expected strand direction (+-=) after TG id: " << filename);
      }
      newtag.strand=(lexer->YYText())[0];

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_INT){
	MIRANOTIFY(Notify::FATAL, "expected int after TG: " << filename);
      }
      newtag.from=atoi(lexer->YYText());
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_PP){
	MIRANOTIFY(Notify::FATAL, "expected .. after TG id int: " << filename);
      }
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();
      if(yyretcode!=EXPT_INT){
	MIRANOTIFY(Notify::FATAL, "expected int after TG id int..:" << filename);
      }
      newtag.to=atoi(lexer->YYText());
      CEBUG(lexer->YYText() << "found");

      yyretcode=lexer->yylex();

      if(yyretcode==EXPT_TAGONELINER){
	CEBUG("search TAG-ONELINER:");
	bool needcopy=false;
	cptr=lexer->YYText();
	while(*cptr){
	  if(*cptr++!=' '){
	    needcopy=true;
	    break;
	  }
	}
	if(needcopy==true){
	  CEBUG("found:");
	  cptr=lexer->YYText();
	  while(*cptr){
	    if(*cptr!='\n') newtag.comment+=*cptr;
	    cptr++;
	  }
	}
	yyretcode=lexer->yylex();
      }

      uint32 numofcomments=0;
      while(yyretcode==EXPT_TGC){
	CEBUG("Found TGC:");
	lexer->yylex();
	if(numofcomments>0) newtag.comment+='\n';
	numofcomments++;
	CEBUG(lexer->YYText() << "found");
	cptr=lexer->YYText();
	while(*cptr){
	  if(*cptr!='\n') newtag.comment+=*cptr;
	  cptr++;
	}
	yyretcode=lexer->yylex();
      }
      dontread=true;
      break;
    }

    case EXPT_ON:{
      CEBUG("Found ON:");
      EXP_on_ranges.clear();
      yyretcode=lexer->yylex();
      int32 val;
      while(yyretcode!=EXPT_EOL){
	if(yyretcode==EXPT_INT){
	  // a single number (0) showing an insert
	  //val=atoi(lexer->YYText());
	  //EXP_on_ranges.push_back(val);
	  //EXP_on_ranges.push_back(val);
	  EXP_on_ranges.push_back(0);
	  EXP_on_ranges.push_back(0);
	}else if(yyretcode==EXPT_RANGELOW){
	  val=atoi(lexer->YYText());
	  EXP_on_ranges.push_back(val);
	  yyretcode=lexer->yylex();
	  if(yyretcode!=EXPT_RANGEHIGH){
	    cerr << "Error while reading: " << lexer->YYText() << "(code: " << yyretcode << ")" << endl;
	    MIRANOTIFY(Notify::FATAL, "expected int after .. in ON tag: " << filename);
	  }
	  val=atoi(lexer->YYText());
	  EXP_on_ranges.push_back(val);
	}else{
	  cerr << "Error while reading: " << lexer->YYText() << "(code: " << yyretcode << ")" << endl;
	  MIRANOTIFY(Notify::FATAL, "expected int after ON: " << filename);
	}
	yyretcode=lexer->yylex();
      }
      break;
    }

    case EXPT_JENAQUIRK: break;

    default:{
      cerr << lexer->YYText() << endl;
      cerr << yyretcode << endl;
      MIRANOTIFY(Notify::FATAL, "illegal token, this looks like an error in the EXP file: " << filename);
    }
    }
  }

  if(EXP_sl<0) EXP_sl=0;
  if(EXP_sl>static_cast<int32> (EXP_len_seq)) {
    //and warning
    EXP_sl=EXP_len_seq-1;
  }
  if(EXP_ql<0) EXP_ql=0;
  if(EXP_ql>static_cast<int32> (EXP_len_seq)) {
    //and warning
    EXP_ql=EXP_len_seq-1;
  }
  if(EXP_sr<0) EXP_sr=EXP_len_seq;
  if(EXP_sr>static_cast<int32> (EXP_len_seq)){
    //and warning
    EXP_sr=EXP_len_seq;
  }
  if(EXP_qr<0) EXP_qr=EXP_len_seq;
  if(EXP_qr>static_cast<int32> (EXP_len_seq)){
    //and warning
    EXP_qr=EXP_len_seq;
  }

  if(EXP_cl> static_cast<int32> (EXP_len_seq)) {
    //and warning
    EXP_cl=EXP_len_seq-1;
  }
  if(EXP_cl<0) EXP_cl=0;

  if(EXP_cr>static_cast<int32> (EXP_len_seq)){
    //and warning
    EXP_cr=EXP_len_seq;
  }
  if(EXP_cr<0) EXP_cr=EXP_len_seq;

  FUNCEND();
  return;
}


void EXP::dump()
{


  cout << "EXP_bc      " << EXP_bc     << endl;
  cout << "EXP_cf      " << EXP_cf     << endl;
  cout << "EXP_cn      " << EXP_cn     << endl;
  cout << "EXP_cv      " << EXP_cv     << endl;
  cout << "EXP_dt      " << EXP_dt     << endl;
  cout << "EXP_en      " << EXP_en     << endl;
  cout << "EXP_id      " << EXP_id     << endl;
  cout << "EXP_le      " << EXP_le     << endl;
  cout << "EXP_li      " << EXP_li     << endl;
  cout << "EXP_ln      " << EXP_ln     << endl;
  cout << "EXP_lt      " << EXP_lt     << endl;
  cout << "EXP_ma      " << EXP_ma     << endl;
  cout << "EXP_mc      " << EXP_mc     << endl;
  cout << "EXP_mn      " << EXP_mn     << endl;
  cout << "EXP_mt      " << EXP_mt     << endl;
  cout << "EXP_ps      " << EXP_ps     << endl;
  cout << "EXP_sf      " << EXP_sf     << endl;
  cout << "EXP_sv      " << EXP_sv     << endl;
  cout << "EXP_ss      " << EXP_ss     << endl;
  cout << "EXP_tn      " << EXP_tn     << endl;
  cout << "EXP_op      " << EXP_op     << endl;
  cout << "EXP_pn      " << EXP_pn     << endl;

  cout << "EXP_sq      " << EXP_sq     << endl;

  cout << "EXP_aq      " << EXP_aq     << endl;
  cout << "EXP_ch      " << EXP_ch     << endl;
  cout << "EXP_pr      " << EXP_pr     << endl;
  cout << "EXP_sc      " << EXP_sc     << endl;
  cout << "EXP_si_from " << EXP_si_from<< endl;
  cout << "EXP_si_to   " << EXP_si_to  << endl;
  cout << "EXP_sp      " << EXP_sp     << endl;
  cout << "EXP_cl      " << EXP_cl     << endl;
  cout << "EXP_cr      " << EXP_cr     << endl;
  cout << "EXP_ql      " << EXP_ql     << endl;
  cout << "EXP_qr      " << EXP_qr     << endl;
  cout << "EXP_sl      " << EXP_sl     << endl;
  cout << "EXP_sr      " << EXP_sr     << endl;

  cout << "EXP_len_seq " << EXP_len_seq<< endl;

}
