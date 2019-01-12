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


#ifndef _edit_parameters_H
#define _edit_parameters_H

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include <stdlib.h>

#include "errorhandling/errorhandling.H"
#include "util/misc.H"

#include <sstream>

#include <fstream>
//#include <algo.h>                     // STL min max template
#include "EdIt/parameters_tokens.h"

#undef yyFlexLexer
#define yyFlexLexer EPFlexLexer
#include <FlexLexer.h>



// 	$Id$



struct ep_tag_parameters
{
  bool tag_alter;
  bool tag_insert;
  bool tag_delete;
  bool tag_asterisk;
  bool tag_consensus;
};


struct ep_file_parameters
{
  char *caf_input;
  char *caf_output;
  char *log_output;
  char *text_output;
  char *html_output;
  bool caf_void;
};



struct ep_strand_parameters
{
  bool  make_double_stranded;
  bool  doublestrand_machinetype;   // verschiedene sequenziermaschinen = doublestranded
  int32 doublestrand_minscore;
  int32 doublestrand_umgebung;
  int32 mincoverage_single;
};


struct ep_evaluation_parameters
{
  int32  min_solutions_checked;    // minimum number of solutions tested
  // threshold: evaluate solutions if their score is less than x(%) over
  // higher than the score of the best solution
  int32  max_score_over_best;
  float  confirmation_threshold;   // KNN result threshold to confirm hypo.
  bool   strict_evaluation;        // Edit only when completely confirmed.
  // Edit all bases below this quality value even without trace confirmation
  // if set to -1 (default) no automatic editing of low quality bases is done
  int32  edit_quality_below;
  bool   edit_by_columns;
  bool   edit_by_regions;
  bool   edit_by_small_regions;
};



struct ep_evaluation_range
{
  int32 contig_id_von;          // edit contigs starting from
  int32 contig_id_bis;          // until...
  int32 contig_len_von;         // edit contigs longer than...
  int32 contig_len_bis;         //              shorter than...
  int32 contig_pos_von;         // edit contig between
  int32 contig_pos_bis;         //              and...
};



class EDITParameters
{
private:

  ep_tag_parameters          tag_parameters;
  ep_file_parameters         file_parameters;
  ep_strand_parameters       strand_parameters;
  ep_evaluation_parameters   eval_parameters;
  ep_evaluation_range        eval_range;

  bool  command_hypo;
  bool  command_eval;
  bool  command_edit;

  int32 verbose;
  bool  show_progress;
  bool  debug;

  bool  log_cleared;

  bool  help;    // parameter "help" was set

  //Functions
private:
  void  foolCompiler();
  int32 gimmeAnInt(FlexLexer * lexer);
  int32 getConAnalyseMode(FlexLexer * lexer);

public:
  EDITParameters();
  ~EDITParameters();

  ep_tag_parameters  const & getTagParams() const ;

  void cmdLine2Stream(int argc, char **argv);
  void parse(istream & is);

  int32 setTagParams(ep_tag_parameters & ap);
  void  logOperation(char *s) const;

  void  setFilenames(char *input, char *output, char *logfile);
  float getConfirmationThreshold() const
    { return eval_parameters.confirmation_threshold; }
  void  setConfirmationThreshold(const float f)
    { eval_parameters.confirmation_threshold = f; }
  bool getStrictEvaluation()
    { return eval_parameters.strict_evaluation; }
  void  setStrictEvaluation(const bool mode)
    { eval_parameters.strict_evaluation = mode; }
  bool  getEditLowQualityBases()
    { return (eval_parameters.edit_quality_below > 0); }
  int32 getLowQualityThreshold()
    { return eval_parameters.edit_quality_below; }
  void  setEditLowQualityBases(int32 quality_threshold)
    { eval_parameters.edit_quality_below = quality_threshold; }
  bool  getEditByRegion() const
    { return eval_parameters.edit_by_regions; }
  void  setEditByRegion(bool b)
    { eval_parameters.edit_by_regions = b; }
  bool  getEditByColumn() const
    { return eval_parameters.edit_by_columns; }
  void  setEditByColumn(bool b)
    { eval_parameters.edit_by_columns = b; }
  bool  getEditBySmallRegions() const
    { return eval_parameters.edit_by_small_regions; }
  void  setEditBySmallRegions(bool b)
    { eval_parameters.edit_by_small_regions = b; }

  bool   tagInsertOperations() const
    { return tag_parameters.tag_insert; };
  bool   tagDeleteOperations() const
    { return tag_parameters.tag_delete; };
  bool   tagAlterOperations() const
    { return tag_parameters.tag_alter; };
  bool   tagDeleteAsteriskColumns() const
    { return tag_parameters.tag_asterisk; };
  bool   tagConsensus() const
    { return tag_parameters.tag_consensus; };

  char * getCafInputFileName()  const
    { return file_parameters.caf_input; };
  char * getCafOutputFileName() const
    { return file_parameters.caf_output; };
  bool   createCafOutput() const
    { return file_parameters.caf_void == false; }
  bool   createHtmlOutput() const
    { return file_parameters.html_output != nullptr; }
  char*  getHtmlOutputFileName() const
    { return file_parameters.html_output; }
  bool   createTextOutput() const
    { return file_parameters.text_output != nullptr; }
  char*  getTextOutputFileName() const
    { return file_parameters.text_output; }

  bool   makeDoubleStranded() const
    { return strand_parameters.make_double_stranded; };
  void   makeDoubleStranded(const bool ds)
    { strand_parameters.make_double_stranded = ds; };
  int32  getDoubleStrandMinScore() const
    { return strand_parameters.doublestrand_minscore; };
  int32  getDoublestrandUmgebung() const
    { return strand_parameters.doublestrand_umgebung; };
  int32  getMinCoverage() const
    { return strand_parameters.mincoverage_single; };

  int32  getMinSolutionsChecked() const
    { return eval_parameters.min_solutions_checked; };
  int32  getMaxScoreOverBest() const
    { return eval_parameters.max_score_over_best; };

  bool doSomething() const
    { return command_hypo == true; }
  bool doHypo() const
    { return command_hypo == true; }
  bool doEval() const
    { return command_eval == true; }
  bool doEdit() const
    { return command_edit == true; }


  bool  doContig(int32 id);
  int32 evalAt(int32 pos)
    {
      if (eval_range.contig_pos_von == -1) return true;
      return (pos >= eval_range.contig_pos_von &&
	      pos <= eval_range.contig_pos_bis);
    }
  void  setEvalRange(int32 pos_von, int32 pos_bis)
    {
      eval_range.contig_pos_von = pos_von;
      eval_range.contig_pos_bis = pos_bis;
    }

  bool isVerbose() const
    { return verbose > 0; }

  bool isVerbose(int32 i) const
    { return verbose > i; }

  bool isDebug() const
    { return debug == true; }

  void setVerbose(bool isitVerbose)
    { if (isitVerbose) { verbose = 1; } else { verbose = 0; }}

  void setVerbose(int32 verbmode)
    { verbose = verbmode; }

  void setShowProgress(bool b)
    { show_progress = b; }

  bool getShowProgress()
    { return show_progress; }



  void setDoNothing() {
    command_hypo = false;
    command_eval = false;
    command_edit = false;
  }

  void setDoHypo() {
    command_hypo = true;
    command_eval = true;
    command_edit = false;
  }

  void setDoEval() {
    command_hypo = true;
    command_eval = true;
    command_edit = false;
  }

  void setDoAll() {
    command_hypo = true;
    command_eval = true;
    command_edit = true;
  }

  bool getHelp() { return help; }
  void setHelp(bool h) { help = h; }

  friend ostream &operator<<(ostream &ostr, EDITParameters const &i);
};

#endif
