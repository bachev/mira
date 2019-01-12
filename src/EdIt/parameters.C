/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
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

#include "EdIt/parameters.H"


EDITParameters::EDITParameters()
{

  tag_parameters.tag_insert    = true;
  tag_parameters.tag_delete    = true;
  tag_parameters.tag_alter     = true;
  tag_parameters.tag_asterisk  = true;
  tag_parameters.tag_consensus = true;


  file_parameters.caf_input   = conditionalStrCpy("EdIt.in.caf");
  file_parameters.caf_output  = conditionalStrCpy("EdIt.out.caf");
  file_parameters.log_output  = conditionalStrCpy("EdIt.log");
  file_parameters.text_output = nullptr;
  file_parameters.html_output = nullptr;
  file_parameters.caf_void = false;

  strand_parameters.make_double_stranded     = true;
  strand_parameters.doublestrand_machinetype = false;
  strand_parameters.doublestrand_minscore = 85;
  strand_parameters.doublestrand_umgebung = 2;
  strand_parameters.mincoverage_single = 4;

  eval_parameters.min_solutions_checked  = 3;
  eval_parameters.max_score_over_best    = 20;
  eval_parameters.confirmation_threshold = 0.6;

  // eval_parameters.strict_evaluation      = true;
  eval_parameters.strict_evaluation      = false;
  eval_parameters.edit_quality_below     = -1;

  setEditByColumn(false);
  setEditByRegion(true);
  setEditBySmallRegions(false);

  command_edit   = true;   // Perform edit operations
  command_eval   = true;   // Evaluate hypotheses
  command_hypo   = true;   // Generate hypotheses

  eval_range.contig_id_von = -1;     // for all Contigs
  eval_range.contig_id_bis = -1;     // for all Contigs
  eval_range.contig_pos_von = -1;    // for all bases
  eval_range.contig_pos_bis = -1;    // for all bases

  verbose       = 0;
  show_progress = true;
  debug         = 0;

  help = false;   // print usage() instead of doing something

  log_cleared=false;
}


EDITParameters::~EDITParameters()
{
  delete [] file_parameters.caf_input;
  delete [] file_parameters.caf_output;
  delete [] file_parameters.log_output;
  delete [] file_parameters.html_output;
  delete [] file_parameters.text_output;
}


// return 0 if ok
int32 EDITParameters::setTagParams(ep_tag_parameters & tp)
{
  tag_parameters = tp;
  return 0;
}


ep_tag_parameters const & EDITParameters::getTagParams() const
{
  return tag_parameters;
}



void EDITParameters::logOperation(char *s) const
{
  static int32 count = 0;
  ofstream x;

  if (file_parameters.caf_void == false) {
    if(!log_cleared){
      x.open(file_parameters.log_output);
      x.close();
      bool & lc=const_cast<bool &>(log_cleared);
      lc=true;
    }
    x.open(file_parameters.log_output, ios::app);
    x << ++count << "\t" << s << endl;
    x.close();
  }
}




void EDITParameters::cmdLine2Stream(int argc, char **argv)
{
  stringstream tss;

  for(int32 i=1; i<argc; i++){
    //cout << "Doing " << i << "\t" << argv[i] << endl;
    tss << argv[i] << "-=BEGIN0=-";
  }

  parse(tss);
}


int32 EDITParameters::gimmeAnInt(FlexLexer * lexer)
{
  if(lexer->yylex() != EP_INT){
    cout << "Expected an int";
  }
  return atoi(lexer->YYText());
}


int32 EDITParameters::getConAnalyseMode(FlexLexer * lexer)
{
  int32 tmp=lexer->yylex();
  if(tmp==666){
    cout << "Non recognised am option";
    tmp=0;
  }
  return tmp;
}


bool EDITParameters::doContig(int32 id)
{
  // evaluate all contigs
  if (eval_range.contig_id_von < 0) return true;

  // evaluate contigs starting with...
  if (id >= eval_range.contig_id_von) {
    if (eval_range.contig_id_bis < 0) {
      return true;
    } else {
      return (id <= eval_range.contig_id_bis);
    }
  } else {
    return false;
  }
}


void EDITParameters::setFilenames(char *input, char *output, char *logfile)
{
  if (input != nullptr) {
    delete [] file_parameters.caf_input;
    file_parameters.caf_input = new char[strlen(input) + 1];
    strcpy(file_parameters.caf_input, input);
  }

  if (output != nullptr) {
    delete [] file_parameters.caf_output;
    file_parameters.caf_output = new char[strlen(output) + 1];
    strcpy(file_parameters.caf_output, output);
  }

  if (logfile != nullptr) {
    delete [] file_parameters.log_output;
    file_parameters.log_output = new char[strlen(logfile) + 1];
    strcpy(file_parameters.log_output, logfile);
  }

  file_parameters.caf_void = (logfile == nullptr);
}



void EDITParameters::parse(istream & is)
{
  FlexLexer* lexer = new EPFlexLexer(&is);

  int yyretcode=-1;

  while(yyretcode!=0){
    yyretcode=lexer->yylex();

    switch(yyretcode){
    case 0: break;                              // do nothing, eof

    case EP_ALTER_TAG :{
      tag_parameters.tag_alter = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: alter tag set to "
	   << tag_parameters.tag_alter << endl;
      break;
    }

    case EP_DELETE_TAG: {
      tag_parameters.tag_delete = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: delete tag set to "
	   << tag_parameters.tag_alter << endl;
      break;
    }

    case EP_INSERT_TAG: {
      tag_parameters.tag_insert = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: insert tag set to "
	   << tag_parameters.tag_alter << endl;
      break;
    }

    case EP_CONSENSUS_TAG: {
      tag_parameters.tag_consensus = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: consensus tag set to "
	   << tag_parameters.tag_consensus << endl;
      break;
    }

    case EP_VERBOSE: {
      verbose = gimmeAnInt(lexer);
      cout << "PARAMETER: verbose mode set to " << verbose << endl;
#ifdef RUNONLY
      cout << "           compiled as run-only; some output ignored" << endl;
#endif
      break;
    }

    case EP_DEBUG: {
      debug = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: debug mode set to " << debug << endl;
      break;
    }

    case EP_ALL_TAGS: {
      tag_parameters.tag_insert    = true;
      tag_parameters.tag_delete    = true;
      tag_parameters.tag_alter     = true;
      tag_parameters.tag_consensus = true;
      cout << "PARAMETER: insert, delete, alter, consensus tag set to 1."
	   << endl;
      break;
    }

   case EP_NO_TAGS: {
      tag_parameters.tag_insert    = false;
      tag_parameters.tag_delete    = false;
      tag_parameters.tag_alter     = false;
      tag_parameters.tag_consensus = false;
      cout << "PARAMETER: insert, delete, alter, consensus tag set to 0."
	   << endl;
      break;
    }

    case EP_DELETE_ASTERISK: {
      tag_parameters.tag_asterisk = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: delete asterisk columns set to "
	   << tag_parameters.tag_asterisk << endl;
      break;
    }

    case EP_caf_in_file: {
      delete [] file_parameters.caf_input;
      file_parameters.caf_input = new char[lexer->YYLeng() + 1];
      strcpy(file_parameters.caf_input, lexer->YYText());
      cout << "PARAMETER: input file : "
	   << file_parameters.caf_input << endl;
      break;
    }

    case EP_caf_out_file: {
      int32 l;

      delete [] file_parameters.caf_output;
      file_parameters.caf_output = new char[lexer->YYLeng() + 1];
      strcpy(file_parameters.caf_output, lexer->YYText());
      cout << "PARAMETER: output file : "
	   << file_parameters.caf_output << endl;
      break;
    }

   case EP_log_out_file: {
      int32 l;

      delete [] file_parameters.log_output;
      file_parameters.log_output = new char[lexer->YYLeng() + 1];
      strcpy(file_parameters.log_output, lexer->YYText());
      cout << "PARAMETER: log file : "
	   << file_parameters.log_output << endl;
      break;
    }

    case EP_html_out_file: {
      int32 l;

      delete [] file_parameters.html_output;
      file_parameters.html_output = new char[lexer->YYLeng() + 1];
      strcpy(file_parameters.html_output, lexer->YYText());
      cout << "PARAMETER: html file : "
	   << file_parameters.html_output << endl;
      break;
    }


    case EP_text_out_file: {
      int32 l;

      delete [] file_parameters.text_output;
      file_parameters.text_output = new char[lexer->YYLeng() + 1];
      strcpy(file_parameters.text_output, lexer->YYText());
      cout << "PARAMETER: text file : "
	   << file_parameters.text_output << endl;
      break;
    }


    case EP_outfile_void: {
      file_parameters.caf_void = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: create caf output: "
	   << file_parameters.caf_void << endl;
      break;
    }


    case EP_make_double: {
      strand_parameters.make_double_stranded = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: make double-stranded "
	   << strand_parameters.make_double_stranded << endl;
      break;
    }

    case EP_use_machinetype: {
      strand_parameters.doublestrand_machinetype = ((gimmeAnInt(lexer) == 1));
      cout << "PARAMETER: make double-stranded if different sequencers "
	   << strand_parameters.doublestrand_machinetype << endl;
      break;
    }

    case EP_doublestrand_minscore: {
      int32 score =  gimmeAnInt(lexer);
      if (score <= 100 && score >= 0) {
	strand_parameters.doublestrand_minscore = score;
      } else {
	cerr << "ERROR: minscore doublestranding outside [0,100]!\n";
      }
      break;
    }

    case EP_mincoverage_single: {
      int32 coverage = gimmeAnInt(lexer);
      if (coverage <= 100 && coverage >= -1) {
	strand_parameters.mincoverage_single = coverage;
	cout << "PARAMETER: minimum coverage for single stranded regions: "
	     << coverage << endl;
      } else {
	cerr << "ERROR: parameter 'cover' outside [-1,100]!\n";
      }
      break;
    }

    case EP_ANID: {
      cout <<  "PARAMETER: ERROR - unknown option: "
	   << lexer->YYText() <<endl;
      cout.flush();
      break;
    }

    case EP_DO_NOP: {
      cout << "PARAMETER: do nothing (read the caf-file)" << endl;
      command_hypo = command_eval = command_edit = false;
      break;
    }

    case EP_DO_ALL: {
      cout << "PARAMETER: do all (default)" << endl;
      command_hypo = command_eval = command_edit = true;
      break;
    }

    case EP_DO_HYPO: {
      cout << "PARAMETER: generate hypotheses only (no eval, no edit)"
	   << endl;
      command_hypo = true;
      command_eval = command_edit = false;
      break;
    }

    case EP_DO_EVAL: {
      cout << "PARAMETER: evaluate hypotheses but do not edit" << endl;
      command_hypo = command_eval = true;
      command_edit = false;
      break;
    }

    case EP_DO_CONTIG: {
      eval_range.contig_id_von = gimmeAnInt(lexer);
      eval_range.contig_id_bis = eval_range.contig_id_von ;

      cout << "PARAMETER: edit single contig - "
	   << eval_range.contig_id_von << endl;

      setDoAll();

      break;
    }

    case EP_DO_RANGE: {
      eval_range.contig_pos_von = gimmeAnInt(lexer);
      eval_range.contig_pos_bis = gimmeAnInt(lexer) ;

      cout << "PARAMETER: eval range von "
	   << eval_range.contig_pos_von << " - "
	   << eval_range.contig_pos_bis << endl;

      break;
    }

    case EP_DO_STRICT: {
      cout << "PARAMETER: strict evaluation ON" << endl;
      eval_parameters.strict_evaluation = true;
      break;
    }

    case EP_DO_NOSTRICT: {
      cout << "PARAMETER: strict evaluation OFF" << endl;
      eval_parameters.strict_evaluation = false;
      break;
    }

    case EP_DO_THRESHOLD: {
      setConfirmationThreshold((float)(gimmeAnInt(lexer))/100.0) ;
      cout << "PARAMETER: confirmation threshold "
	   << getConfirmationThreshold() << endl;
      break;
    }

    case EP_DO_LOWQUAL: {
      setEditLowQualityBases(gimmeAnInt(lexer));
      cout << "PARAMETER: edit if base quality below "
	   << getLowQualityThreshold() << endl;
      break;
    }

    case EP_DO_REGIONS: {
      setEditByRegion(gimmeAnInt(lexer) != 0);
      cout << "PARAMETER: edit by fault region is "
	   << getEditByRegion() << endl;
      break;
    }

    case EP_DO_COLUMNS: {
      setEditByColumn(gimmeAnInt(lexer) != 0);
      cout << "PARAMETER: edit by column is "
	   << getEditByColumn() << endl;
      break;
    }

    case EP_DO_SMALL_REGIONS: {
      setEditBySmallRegions(gimmeAnInt(lexer) != 0);
      cout << "PARAMETER: edit by small regions is "
	   << getEditBySmallRegions() << endl;
      break;
    }

    case EP_HELP: {
      setHelp(true);
      break;
    }

    default:{
	cout <<  "PARAMETER: ERROR - recognised but not caught: "
	     << lexer->YYText() <<endl;
	cout.flush();
    }
    }

  }

  delete lexer;

}



ostream &operator<<(ostream &ostr, EDITParameters const &i)
{
  ostr << "\nEdIt parameters : " << endl;
  ostr << "  Edit by:";

  if (i.getEditByRegion()) { ostr << " 'regions'"; }
  if (i.getEditBySmallRegions()) { ostr << " 'small regions'"; }
  if (i.getEditByColumn()) { ostr << "'columns'"; }
  ostr << "." << endl;

  ostr << "  Edit single contig: " << i.eval_range.contig_id_von << endl;
  ostr << "  Check at least " << i.eval_parameters.min_solutions_checked
       << " solutions until " << i.eval_parameters.max_score_over_best
       << " % over best score." << endl;

  ostr << "  Threshold for network evaluation: "
       << i.eval_parameters.confirmation_threshold
       << "." << endl;
  ostr << "  Edit quality below "
       << i.eval_parameters.edit_quality_below << endl;
  ostr << "  Strict evaluation :" ;
  if (i.eval_parameters.strict_evaluation) {
    ostr << " YES."  << endl;
  } else {
    ostr << " NO." << endl;
  }

  ostr << "  Tags : ";
  if (i.tag_parameters.tag_insert)    { ostr << " 'insert'"; }
  if (i.tag_parameters.tag_delete)    { ostr << " 'delete'"; }
  if (i.tag_parameters.tag_alter)     { ostr << " 'alter'"; }
  if (i.tag_parameters.tag_consensus) { ostr << " 'consensus'"; }
  if (i.tag_parameters.tag_asterisk)  { ostr << " 'asterisk'"; }
  ostr << "." << endl;

  ostr << "  Files: ";
  ostr << " " << i.file_parameters.caf_input << " ==> "
       << i.file_parameters.caf_output
       << "     log: " << i.file_parameters.log_output << endl;

  ostr << endl;
  return ostr;
}
