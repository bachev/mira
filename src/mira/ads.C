/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1998-2000 by the German Cancer Research Center (Deutsches
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
 * 2 Sequences aligned routines
 * AlignedDualSeq manages the data produced by the alignement routines. Two
 *  sequences have been aligned and are stored in an ADS object.
 * Routines for computing scores and some other classification number are
 *  provided, too.
 *
 */


#include "ads.H"

#include "errorhandling/errorhandling.H"
#include "util/dptools.H"
#include "mira/ads.H"

#include <cmath>

using std::cout;
using std::endl;


#define CEBUG(bla)

std::vector<double> AlignedDualSeq::ADS_powofstarcounter;

uint32 AlignedDualSeq::ADS_c_matrix_valid=0;
uint32 AlignedDualSeq::ADS_s_matrix_valid=0;

int8 AlignedDualSeq::ADS_realscore_matrix[ADS_MATSIZE][ADS_MATSIZE];
int8 AlignedDualSeq::ADS_expectedscore_matrix[ADS_MATSIZE][ADS_MATSIZE];
char AlignedDualSeq::ADS_consensus_matrix[ADS_MATSIZE][ADS_MATSIZE];
char AlignedDualSeq::ADS_consensus_gap_matrix[ADS_MATSIZE][ADS_MATSIZE];




const bool AlignedDualSeq::ADS_initialisedstatics=AlignedDualSeq::staticInitialiser();



/*************************************************************************
 *
 * a static class initialiser
 *
 * initialises a few static variables / arrays _before_ main is called
 *
 *
 *************************************************************************/

bool AlignedDualSeq::staticInitialiser()
{
  ADS_powofstarcounter.resize(50,1);
  for(uint32 i=1; i<ADS_powofstarcounter.size(); i++){
    ADS_powofstarcounter[i]=pow(static_cast<double>(.9),static_cast<double>(i));
  }
  return true;
}


// Constructor with sequences not known
AlignedDualSeq::AlignedDualSeq(MIRAParameters * params) : AlignedDualSeqFacts ()
{
  FUNCSTART("AlignedDualSeq::AlignedDualSeq(MIRAParameters * params)");

  ADS_valid=0;
  ADS_miraparams=params;

  init();

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void AlignedDualSeq::init()
{
  ADS_aligned_seq1=nullptr;
  ADS_aligned_seq2=nullptr;
  ADS_consensus_seq=nullptr;
  ADS_consensus_gap_seq=nullptr;
  ADS_cursize=0;
  ADS_seq1=nullptr;
  ADS_seq2=nullptr;

  ADS_contained=0;
  ADS_score=0;
  ADS_expected_score=0;
  ADS_weight=0;
  ADS_nummismatches=0;
  ADS_numgaps=0;
  ADS_maxcontiguousgaps=0;
  ADS_len1=0;
  ADS_len2=0;

  ADS_affinegapscore=false;

  ADS_minbanddistance=0;
  ADS_bandwidthused=0;

  // initialise the following match matrices
  //  ADS_realscore_matrix
  //  ADS_expectedscore_matrix,
  //  ADS_consensus_matrix
  //
  // can't do that in static initialiser because we're using dyp.*
  // TODO: perhaps rethink

  if(!ADS_s_matrix_valid){
    dynamic_parameters dyp= ADS_miraparams->getDynamicParams();
    // these should not happen
    matinit((char *) &ADS_realscore_matrix,' ',' ',static_cast<char>(dyp.dyn_score_nmatch));
    matinit((char *) &ADS_realscore_matrix,' ','*',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,' ','1',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,' ','2',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,' ','3',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,' ','4',(char) dyp.dyn_score_nmatch);
    // these ok
    matinit((char *) &ADS_realscore_matrix,' ','N',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','X',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','A',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','C',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','G',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','T',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','M',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','R',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','W',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','S',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','Y',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','K',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','V',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','H',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','D',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_realscore_matrix,' ','B',(char) dyp.dyn_score_termgap);

    // these should not happen
    matinit((char *) &ADS_realscore_matrix,'*','*',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','1',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','2',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','3',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','4',(char) dyp.dyn_score_gap);
    // these ok
    matinit((char *) &ADS_realscore_matrix,'*','N',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','X',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','A',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','C',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','G',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','T',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','M',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','R',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','W',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','S',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','Y',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','K',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','V',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','H',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','D',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'*','B',(char) dyp.dyn_score_gap);

    matinit((char *) &ADS_realscore_matrix,'#','*',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'#','N',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','X',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','A',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','C',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','G',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','T',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','M',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','R',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','W',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','S',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','Y',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','K',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','V',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','H',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','D',(char) dyp.dyn_score_gap);
    matinit((char *) &ADS_realscore_matrix,'#','B',(char) dyp.dyn_score_gap);

    matinit((char *) &ADS_realscore_matrix,'1','N',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','X',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','A',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'1','C',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','G',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','T',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','M',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','R',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','W',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','S',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','Y',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','K',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','V',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','H',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','D',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'1','B',(char) dyp.dyn_score_oldgap);

    matinit((char *) &ADS_realscore_matrix,'2','N',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','X',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','A',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','C',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'2','G',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','T',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','M',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','R',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','W',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','S',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','Y',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','K',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','V',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','H',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','D',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'2','B',(char) dyp.dyn_score_oldgap);

    matinit((char *) &ADS_realscore_matrix,'3','N',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','X',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','A',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','C',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','G',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'3','T',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','M',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','R',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','W',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','S',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','Y',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','K',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','V',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','H',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','D',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'3','B',(char) dyp.dyn_score_oldgap);

    matinit((char *) &ADS_realscore_matrix,'4','N',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','X',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','A',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','C',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','G',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','T',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'4','M',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','R',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','W',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','S',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','Y',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','K',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','V',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','H',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','D',(char) dyp.dyn_score_oldgap);
    matinit((char *) &ADS_realscore_matrix,'4','B',(char) dyp.dyn_score_oldgap);

    matinit((char *) &ADS_realscore_matrix,'1','1',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'1','2',(char) dyp.dyn_score_oldgapmatch);
    matinit((char *) &ADS_realscore_matrix,'1','3',(char) dyp.dyn_score_oldgapmatch);
    matinit((char *) &ADS_realscore_matrix,'1','4',(char) dyp.dyn_score_oldgapmatch);
    matinit((char *) &ADS_realscore_matrix,'1','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'2','2',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'2','3',(char) dyp.dyn_score_oldgapmatch);
    matinit((char *) &ADS_realscore_matrix,'2','4',(char) dyp.dyn_score_oldgapmatch);
    matinit((char *) &ADS_realscore_matrix,'2','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'3','3',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'3','4',(char) dyp.dyn_score_oldgapmatch);
    matinit((char *) &ADS_realscore_matrix,'3','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'4','4',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'4','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'N','N',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','X',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','A',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','C',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','G',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','T',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','M',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','R',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','V',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','H',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','D',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'N','B',(char) dyp.dyn_score_nmatch);


    matinit((char *) &ADS_realscore_matrix,'X','X',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','N',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','A',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','C',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','G',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','T',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','M',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','R',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','V',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','H',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','D',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'X','B',(char) dyp.dyn_score_nmatch);

    matinit((char *) &ADS_realscore_matrix,'A','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'A','C',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'A','G',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'A','T',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'A','M',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'A','R',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'A','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'A','S',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'A','Y',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'A','K',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'A','V',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'A','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'A','D',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'A','B',(char) dyp.dyn_score_mismatch);

    matinit((char *) &ADS_realscore_matrix,'C','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'C','G',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'C','T',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'C','M',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'C','R',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'C','W',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'C','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'C','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'C','K',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'C','V',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'C','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'C','D',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'C','B',(char) dyp.dyn_score_halfmismatch);

    matinit((char *) &ADS_realscore_matrix,'G','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'G','T',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'G','M',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'G','R',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'G','W',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'G','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'G','Y',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'G','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'G','V',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'G','H',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'G','D',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'G','B',(char) dyp.dyn_score_halfmismatch);

    matinit((char *) &ADS_realscore_matrix,'T','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'T','M',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'T','R',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'T','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'T','S',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'T','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'T','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'T','V',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'T','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'T','D',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'T','B',(char) dyp.dyn_score_halfmismatch);

    matinit((char *) &ADS_realscore_matrix,'M','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'M','R',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'M','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'M','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'M','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'M','K',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'M','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'M','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'M','D',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'M','B',(char) dyp.dyn_score_halfmismatch);

    matinit((char *) &ADS_realscore_matrix,'R','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'R','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'R','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'R','Y',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'R','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'R','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'R','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'R','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'R','B',(char) dyp.dyn_score_halfmismatch);

    matinit((char *) &ADS_realscore_matrix,'W','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'W','S',(char) dyp.dyn_score_mismatch);
    matinit((char *) &ADS_realscore_matrix,'W','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'W','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'W','V',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'W','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'W','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'W','B',(char) dyp.dyn_score_halfmismatch);

    matinit((char *) &ADS_realscore_matrix,'S','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'S','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'S','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'S','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'S','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'S','D',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'S','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'Y','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'Y','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_realscore_matrix,'Y','V',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'Y','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'Y','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'Y','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'K','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'K','V',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'K','H',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_realscore_matrix,'K','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'K','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'V','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'V','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'V','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'V','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'H','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'H','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_realscore_matrix,'H','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'D','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_realscore_matrix,'D','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_realscore_matrix,'B','B',(char) dyp.dyn_score_mismatch);




   // shouldn't happen
    matinit((char *) &ADS_expectedscore_matrix,' ',' ',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,' ','*',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,' ','1',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,' ','2',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,' ','3',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,' ','4',(char) dyp.dyn_score_nmatch);
    // these ok
    matinit((char *) &ADS_expectedscore_matrix,' ','N',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','X',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','A',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','C',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','G',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','T',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','M',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','R',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','W',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','S',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','Y',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','K',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','V',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','H',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','D',(char) dyp.dyn_score_termgap);
    matinit((char *) &ADS_expectedscore_matrix,' ','B',(char) dyp.dyn_score_termgap);

      // shouldn't happen
    matinit((char *) &ADS_expectedscore_matrix,'*','*',(char) dyp.dyn_score_match);
    // these ok
    matinit((char *) &ADS_expectedscore_matrix,'*','1',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'*','2',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'*','3',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'*','4',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'*','N',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'*','X',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'*','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'*','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'#','*',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','N',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'#','X',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'#','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'#','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'1','N',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','X',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','A',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'1','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'2','N',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','X',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','C',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'2','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'3','N',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','X',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','G',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'3','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'4','N',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','X',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','T',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'4','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'1','1',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'1','2',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_expectedscore_matrix,'1','3',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_expectedscore_matrix,'1','4',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_expectedscore_matrix,'1','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'2','2',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'2','3',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_expectedscore_matrix,'2','4',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_expectedscore_matrix,'2','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'3','3',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'3','4',(char) dyp.dyn_score_halfmismatch);
    matinit((char *) &ADS_expectedscore_matrix,'3','#',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'4','4',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'4','#',(char) dyp.dyn_score_halfmatch);


    matinit((char *) &ADS_expectedscore_matrix,'N','N',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','X',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','A',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','C',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','G',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','T',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','M',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','R',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','W',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','S',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','Y',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','K',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','V',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','H',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','D',(char) dyp.dyn_score_npenaltymatch);
    matinit((char *) &ADS_expectedscore_matrix,'N','B',(char) dyp.dyn_score_npenaltymatch);

    matinit((char *) &ADS_expectedscore_matrix,'X','X',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','N',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','A',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','C',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','G',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','T',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','M',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','R',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','W',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','S',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','Y',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','K',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','V',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','H',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','D',(char) dyp.dyn_score_nmatch);
    matinit((char *) &ADS_expectedscore_matrix,'X','B',(char) dyp.dyn_score_nmatch);

    matinit((char *) &ADS_expectedscore_matrix,'A','A',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'A','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'C','C',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'C','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'G','G',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'G','B',(char) dyp.dyn_score_match);

    matinit((char *) &ADS_expectedscore_matrix,'T','T',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'T','B',(char) dyp.dyn_score_match);

    // best thing we could hope for with different bases is a halfmatch
    matinit((char *) &ADS_expectedscore_matrix,'M','M',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'M','R',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','W',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','S',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','Y',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','K',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'M','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'R','R',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'R','W',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','S',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','Y',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','K',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'R','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'W','W',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'W','S',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'W','Y',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'W','K',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'W','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'W','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'W','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'W','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'S','S',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'S','Y',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'S','K',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'S','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'S','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'S','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'S','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'Y','Y',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'Y','K',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'Y','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'Y','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'Y','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'Y','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'K','K',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'K','V',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'K','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'K','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'K','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'V','V',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'V','H',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'V','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'V','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'H','H',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'H','D',(char) dyp.dyn_score_halfmatch);
    matinit((char *) &ADS_expectedscore_matrix,'H','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'D','D',(char) dyp.dyn_score_match);
    matinit((char *) &ADS_expectedscore_matrix,'D','B',(char) dyp.dyn_score_halfmatch);

    matinit((char *) &ADS_expectedscore_matrix,'B','B',(char) dyp.dyn_score_match);



    ADS_s_matrix_valid=1;
  }

  // build up a matrix once for all objects
  if(!ADS_c_matrix_valid){

    matinit((char *) &ADS_consensus_matrix,' ',' ',' ');
    matinit((char *) &ADS_consensus_matrix,' ','*','*');
    matinit((char *) &ADS_consensus_matrix,' ','#','*');
    matinit((char *) &ADS_consensus_matrix,' ','N','N');
    matinit((char *) &ADS_consensus_matrix,' ','X','X');
    matinit((char *) &ADS_consensus_matrix,' ','A','A');
    matinit((char *) &ADS_consensus_matrix,' ','C','C');
    matinit((char *) &ADS_consensus_matrix,' ','G','G');
    matinit((char *) &ADS_consensus_matrix,' ','T','T');
    matinit((char *) &ADS_consensus_matrix,' ','M','M');
    matinit((char *) &ADS_consensus_matrix,' ','R','R');
    matinit((char *) &ADS_consensus_matrix,' ','W','W');
    matinit((char *) &ADS_consensus_matrix,' ','S','S');
    matinit((char *) &ADS_consensus_matrix,' ','Y','Y');
    matinit((char *) &ADS_consensus_matrix,' ','K','K');
    matinit((char *) &ADS_consensus_matrix,' ','V','V');
    matinit((char *) &ADS_consensus_matrix,' ','H','H');
    matinit((char *) &ADS_consensus_matrix,' ','D','D');
    matinit((char *) &ADS_consensus_matrix,' ','B','B');

    matinit((char *) &ADS_consensus_matrix,'*','*','*');
    matinit((char *) &ADS_consensus_matrix,'*','N','N');
    matinit((char *) &ADS_consensus_matrix,'*','X','X');
    matinit((char *) &ADS_consensus_matrix,'*','A','A');
    matinit((char *) &ADS_consensus_matrix,'*','C','C');
    matinit((char *) &ADS_consensus_matrix,'*','G','G');
    matinit((char *) &ADS_consensus_matrix,'*','T','T');
    matinit((char *) &ADS_consensus_matrix,'*','M','M');
    matinit((char *) &ADS_consensus_matrix,'*','R','R');
    matinit((char *) &ADS_consensus_matrix,'*','W','W');
    matinit((char *) &ADS_consensus_matrix,'*','S','S');
    matinit((char *) &ADS_consensus_matrix,'*','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'*','K','K');
    matinit((char *) &ADS_consensus_matrix,'*','V','V');
    matinit((char *) &ADS_consensus_matrix,'*','H','H');
    matinit((char *) &ADS_consensus_matrix,'*','D','D');
    matinit((char *) &ADS_consensus_matrix,'*','B','B');

    matinit((char *) &ADS_consensus_matrix,'#','#','*');
    matinit((char *) &ADS_consensus_matrix,'#','*','*');
    matinit((char *) &ADS_consensus_matrix,'#','N','N');
    matinit((char *) &ADS_consensus_matrix,'#','X','X');
    matinit((char *) &ADS_consensus_matrix,'#','A','A');
    matinit((char *) &ADS_consensus_matrix,'#','C','C');
    matinit((char *) &ADS_consensus_matrix,'#','G','G');
    matinit((char *) &ADS_consensus_matrix,'#','T','T');
    matinit((char *) &ADS_consensus_matrix,'#','M','M');
    matinit((char *) &ADS_consensus_matrix,'#','R','R');
    matinit((char *) &ADS_consensus_matrix,'#','W','W');
    matinit((char *) &ADS_consensus_matrix,'#','S','S');
    matinit((char *) &ADS_consensus_matrix,'#','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'#','K','K');
    matinit((char *) &ADS_consensus_matrix,'#','V','V');
    matinit((char *) &ADS_consensus_matrix,'#','H','H');
    matinit((char *) &ADS_consensus_matrix,'#','D','D');
    matinit((char *) &ADS_consensus_matrix,'#','B','B');

    matinit((char *) &ADS_consensus_matrix,'1','#','*');
    matinit((char *) &ADS_consensus_matrix,'1','*','*');
    matinit((char *) &ADS_consensus_matrix,'1','N','N');
    matinit((char *) &ADS_consensus_matrix,'1','X','X');
    matinit((char *) &ADS_consensus_matrix,'1','A','A');
    matinit((char *) &ADS_consensus_matrix,'1','C','C');
    matinit((char *) &ADS_consensus_matrix,'1','G','G');
    matinit((char *) &ADS_consensus_matrix,'1','T','T');
    matinit((char *) &ADS_consensus_matrix,'1','M','M');
    matinit((char *) &ADS_consensus_matrix,'1','R','R');
    matinit((char *) &ADS_consensus_matrix,'1','W','W');
    matinit((char *) &ADS_consensus_matrix,'1','S','S');
    matinit((char *) &ADS_consensus_matrix,'1','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'1','K','K');
    matinit((char *) &ADS_consensus_matrix,'1','V','V');
    matinit((char *) &ADS_consensus_matrix,'1','H','H');
    matinit((char *) &ADS_consensus_matrix,'1','D','D');
    matinit((char *) &ADS_consensus_matrix,'1','B','B');

    matinit((char *) &ADS_consensus_matrix,'2','#','*');
    matinit((char *) &ADS_consensus_matrix,'2','*','*');
    matinit((char *) &ADS_consensus_matrix,'2','N','N');
    matinit((char *) &ADS_consensus_matrix,'2','X','X');
    matinit((char *) &ADS_consensus_matrix,'2','A','A');
    matinit((char *) &ADS_consensus_matrix,'2','C','C');
    matinit((char *) &ADS_consensus_matrix,'2','G','G');
    matinit((char *) &ADS_consensus_matrix,'2','T','T');
    matinit((char *) &ADS_consensus_matrix,'2','M','M');
    matinit((char *) &ADS_consensus_matrix,'2','R','R');
    matinit((char *) &ADS_consensus_matrix,'2','W','W');
    matinit((char *) &ADS_consensus_matrix,'2','S','S');
    matinit((char *) &ADS_consensus_matrix,'2','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'2','K','K');
    matinit((char *) &ADS_consensus_matrix,'2','V','V');
    matinit((char *) &ADS_consensus_matrix,'2','H','H');
    matinit((char *) &ADS_consensus_matrix,'2','D','D');
    matinit((char *) &ADS_consensus_matrix,'2','B','B');

    matinit((char *) &ADS_consensus_matrix,'3','#','*');
    matinit((char *) &ADS_consensus_matrix,'3','*','*');
    matinit((char *) &ADS_consensus_matrix,'3','N','N');
    matinit((char *) &ADS_consensus_matrix,'3','X','X');
    matinit((char *) &ADS_consensus_matrix,'3','A','A');
    matinit((char *) &ADS_consensus_matrix,'3','C','C');
    matinit((char *) &ADS_consensus_matrix,'3','G','G');
    matinit((char *) &ADS_consensus_matrix,'3','T','T');
    matinit((char *) &ADS_consensus_matrix,'3','M','M');
    matinit((char *) &ADS_consensus_matrix,'3','R','R');
    matinit((char *) &ADS_consensus_matrix,'3','W','W');
    matinit((char *) &ADS_consensus_matrix,'3','S','S');
    matinit((char *) &ADS_consensus_matrix,'3','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'3','K','K');
    matinit((char *) &ADS_consensus_matrix,'3','V','V');
    matinit((char *) &ADS_consensus_matrix,'3','H','H');
    matinit((char *) &ADS_consensus_matrix,'3','D','D');
    matinit((char *) &ADS_consensus_matrix,'3','B','B');

    matinit((char *) &ADS_consensus_matrix,'4','#','*');
    matinit((char *) &ADS_consensus_matrix,'4','*','*');
    matinit((char *) &ADS_consensus_matrix,'4','N','N');
    matinit((char *) &ADS_consensus_matrix,'4','X','X');
    matinit((char *) &ADS_consensus_matrix,'4','A','A');
    matinit((char *) &ADS_consensus_matrix,'4','C','C');
    matinit((char *) &ADS_consensus_matrix,'4','G','G');
    matinit((char *) &ADS_consensus_matrix,'4','T','T');
    matinit((char *) &ADS_consensus_matrix,'4','M','M');
    matinit((char *) &ADS_consensus_matrix,'4','R','R');
    matinit((char *) &ADS_consensus_matrix,'4','W','W');
    matinit((char *) &ADS_consensus_matrix,'4','S','S');
    matinit((char *) &ADS_consensus_matrix,'4','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'4','K','K');
    matinit((char *) &ADS_consensus_matrix,'4','V','V');
    matinit((char *) &ADS_consensus_matrix,'4','H','H');
    matinit((char *) &ADS_consensus_matrix,'4','D','D');
    matinit((char *) &ADS_consensus_matrix,'4','B','B');



    matinit((char *) &ADS_consensus_matrix,'N','N','N');
    matinit((char *) &ADS_consensus_matrix,'N','X','N');
    matinit((char *) &ADS_consensus_matrix,'N','A','A');
    matinit((char *) &ADS_consensus_matrix,'N','C','C');
    matinit((char *) &ADS_consensus_matrix,'N','G','G');
    matinit((char *) &ADS_consensus_matrix,'N','T','T');
    matinit((char *) &ADS_consensus_matrix,'N','M','M');
    matinit((char *) &ADS_consensus_matrix,'N','R','R');
    matinit((char *) &ADS_consensus_matrix,'N','W','W');
    matinit((char *) &ADS_consensus_matrix,'N','S','S');
    matinit((char *) &ADS_consensus_matrix,'N','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'N','K','K');
    matinit((char *) &ADS_consensus_matrix,'N','V','V');
    matinit((char *) &ADS_consensus_matrix,'N','H','H');
    matinit((char *) &ADS_consensus_matrix,'N','D','D');
    matinit((char *) &ADS_consensus_matrix,'N','B','B');

    matinit((char *) &ADS_consensus_matrix,'X','X','X');
    matinit((char *) &ADS_consensus_matrix,'X','N','N');
    matinit((char *) &ADS_consensus_matrix,'X','A','A');
    matinit((char *) &ADS_consensus_matrix,'X','C','C');
    matinit((char *) &ADS_consensus_matrix,'X','G','G');
    matinit((char *) &ADS_consensus_matrix,'X','T','T');
    matinit((char *) &ADS_consensus_matrix,'X','M','M');
    matinit((char *) &ADS_consensus_matrix,'X','R','R');
    matinit((char *) &ADS_consensus_matrix,'X','W','W');
    matinit((char *) &ADS_consensus_matrix,'X','S','S');
    matinit((char *) &ADS_consensus_matrix,'X','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'X','K','K');
    matinit((char *) &ADS_consensus_matrix,'X','V','V');
    matinit((char *) &ADS_consensus_matrix,'X','H','H');
    matinit((char *) &ADS_consensus_matrix,'X','D','D');
    matinit((char *) &ADS_consensus_matrix,'X','B','B');

    matinit((char *) &ADS_consensus_matrix,'A','A','A');
    matinit((char *) &ADS_consensus_matrix,'A','C','M');
    matinit((char *) &ADS_consensus_matrix,'A','G','R');
    matinit((char *) &ADS_consensus_matrix,'A','T','W');
    matinit((char *) &ADS_consensus_matrix,'A','M','M');
    matinit((char *) &ADS_consensus_matrix,'A','R','R');
    matinit((char *) &ADS_consensus_matrix,'A','W','W');
    matinit((char *) &ADS_consensus_matrix,'A','S','V');
    matinit((char *) &ADS_consensus_matrix,'A','Y','H');
    matinit((char *) &ADS_consensus_matrix,'A','K','D');
    matinit((char *) &ADS_consensus_matrix,'A','V','V');
    matinit((char *) &ADS_consensus_matrix,'A','H','H');
    matinit((char *) &ADS_consensus_matrix,'A','D','D');
    matinit((char *) &ADS_consensus_matrix,'A','B','N');

    matinit((char *) &ADS_consensus_matrix,'C','C','C');
    matinit((char *) &ADS_consensus_matrix,'C','G','S');
    matinit((char *) &ADS_consensus_matrix,'C','T','Y');
    matinit((char *) &ADS_consensus_matrix,'C','M','M');
    matinit((char *) &ADS_consensus_matrix,'C','R','V');
    matinit((char *) &ADS_consensus_matrix,'C','W','H');
    matinit((char *) &ADS_consensus_matrix,'C','S','S');
    matinit((char *) &ADS_consensus_matrix,'C','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'C','K','B');
    matinit((char *) &ADS_consensus_matrix,'C','V','V');
    matinit((char *) &ADS_consensus_matrix,'C','H','H');
    matinit((char *) &ADS_consensus_matrix,'C','D','N');
    matinit((char *) &ADS_consensus_matrix,'C','B','B');

    matinit((char *) &ADS_consensus_matrix,'G','G','G');
    matinit((char *) &ADS_consensus_matrix,'G','T','K');
    matinit((char *) &ADS_consensus_matrix,'G','M','V');
    matinit((char *) &ADS_consensus_matrix,'G','R','R');
    matinit((char *) &ADS_consensus_matrix,'G','W','D');
    matinit((char *) &ADS_consensus_matrix,'G','S','S');
    matinit((char *) &ADS_consensus_matrix,'G','Y','B');
    matinit((char *) &ADS_consensus_matrix,'G','K','K');
    matinit((char *) &ADS_consensus_matrix,'G','V','V');
    matinit((char *) &ADS_consensus_matrix,'G','H','N');
    matinit((char *) &ADS_consensus_matrix,'G','D','D');
    matinit((char *) &ADS_consensus_matrix,'G','B','B');

    matinit((char *) &ADS_consensus_matrix,'T','T','T');
    matinit((char *) &ADS_consensus_matrix,'T','M','H');
    matinit((char *) &ADS_consensus_matrix,'T','R','D');
    matinit((char *) &ADS_consensus_matrix,'T','W','W');
    matinit((char *) &ADS_consensus_matrix,'T','S','B');
    matinit((char *) &ADS_consensus_matrix,'T','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'T','K','K');
    matinit((char *) &ADS_consensus_matrix,'T','V','N');
    matinit((char *) &ADS_consensus_matrix,'T','H','H');
    matinit((char *) &ADS_consensus_matrix,'T','D','D');
    matinit((char *) &ADS_consensus_matrix,'T','B','B');

    matinit((char *) &ADS_consensus_matrix,'M','M','M');
    matinit((char *) &ADS_consensus_matrix,'M','R','V');
    matinit((char *) &ADS_consensus_matrix,'M','W','H');
    matinit((char *) &ADS_consensus_matrix,'M','S','V');
    matinit((char *) &ADS_consensus_matrix,'M','Y','H');
    matinit((char *) &ADS_consensus_matrix,'M','K','N');
    matinit((char *) &ADS_consensus_matrix,'M','V','V');
    matinit((char *) &ADS_consensus_matrix,'M','H','H');
    matinit((char *) &ADS_consensus_matrix,'M','D','N');
    matinit((char *) &ADS_consensus_matrix,'M','B','N');

    matinit((char *) &ADS_consensus_matrix,'R','R','R');
    matinit((char *) &ADS_consensus_matrix,'R','W','D');
    matinit((char *) &ADS_consensus_matrix,'R','S','V');
    matinit((char *) &ADS_consensus_matrix,'R','Y','N');
    matinit((char *) &ADS_consensus_matrix,'R','K','D');
    matinit((char *) &ADS_consensus_matrix,'R','V','V');
    matinit((char *) &ADS_consensus_matrix,'R','H','N');
    matinit((char *) &ADS_consensus_matrix,'R','D','D');
    matinit((char *) &ADS_consensus_matrix,'R','B','N');

    matinit((char *) &ADS_consensus_matrix,'W','W','W');
    matinit((char *) &ADS_consensus_matrix,'W','S','N');
    matinit((char *) &ADS_consensus_matrix,'W','Y','H');
    matinit((char *) &ADS_consensus_matrix,'W','K','D');
    matinit((char *) &ADS_consensus_matrix,'W','V','N');
    matinit((char *) &ADS_consensus_matrix,'W','H','H');
    matinit((char *) &ADS_consensus_matrix,'W','D','D');
    matinit((char *) &ADS_consensus_matrix,'W','B','A');

    matinit((char *) &ADS_consensus_matrix,'S','S','S');
    matinit((char *) &ADS_consensus_matrix,'S','Y','B');
    matinit((char *) &ADS_consensus_matrix,'S','K','B');
    matinit((char *) &ADS_consensus_matrix,'S','V','S');
    matinit((char *) &ADS_consensus_matrix,'S','H','N');
    matinit((char *) &ADS_consensus_matrix,'S','D','N');
    matinit((char *) &ADS_consensus_matrix,'S','B','B');

    matinit((char *) &ADS_consensus_matrix,'Y','Y','Y');
    matinit((char *) &ADS_consensus_matrix,'Y','K','B');
    matinit((char *) &ADS_consensus_matrix,'Y','V','N');
    matinit((char *) &ADS_consensus_matrix,'Y','H','H');
    matinit((char *) &ADS_consensus_matrix,'Y','D','N');
    matinit((char *) &ADS_consensus_matrix,'Y','B','B');

    matinit((char *) &ADS_consensus_matrix,'K','K','K');
    matinit((char *) &ADS_consensus_matrix,'K','V','N');
    matinit((char *) &ADS_consensus_matrix,'K','H','N');
    matinit((char *) &ADS_consensus_matrix,'K','D','D');
    matinit((char *) &ADS_consensus_matrix,'K','B','B');

    matinit((char *) &ADS_consensus_matrix,'V','V','V');
    matinit((char *) &ADS_consensus_matrix,'V','H','N');
    matinit((char *) &ADS_consensus_matrix,'V','D','N');
    matinit((char *) &ADS_consensus_matrix,'V','B','N');

    matinit((char *) &ADS_consensus_matrix,'H','H','H');
    matinit((char *) &ADS_consensus_matrix,'H','D','N');
    matinit((char *) &ADS_consensus_matrix,'H','B','N');

    matinit((char *) &ADS_consensus_matrix,'D','D','D');
    matinit((char *) &ADS_consensus_matrix,'D','B','B');

    matinit((char *) &ADS_consensus_matrix,'B','B','B');

    memcpy(&ADS_consensus_gap_matrix, &ADS_consensus_matrix, sizeof(ADS_consensus_matrix));

    matinit((char *) &ADS_consensus_gap_matrix,'*','#','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','*','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','A','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','C','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','G','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','T','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','N','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','X','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','M','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','R','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','W','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','S','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','Y','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','K','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','V','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','H','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','D','*');
    matinit((char *) &ADS_consensus_gap_matrix,'*','B','*');

    matinit((char *) &ADS_consensus_gap_matrix,'#','#','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','*','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','A','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','C','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','G','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','T','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','N','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','X','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','M','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','R','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','W','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','S','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','Y','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','K','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','V','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','H','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','D','*');
    matinit((char *) &ADS_consensus_gap_matrix,'#','B','*');

    matinit((char *) &ADS_consensus_gap_matrix,'1','#','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','*','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','A','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','C','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','G','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','T','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','N','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','X','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','M','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','R','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','W','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','S','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','Y','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','K','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','V','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','H','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','D','*');
    matinit((char *) &ADS_consensus_gap_matrix,'1','B','*');

    matinit((char *) &ADS_consensus_gap_matrix,'2','#','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','*','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','A','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','C','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','G','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','T','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','N','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','X','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','M','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','R','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','W','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','S','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','Y','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','K','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','V','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','H','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','D','*');
    matinit((char *) &ADS_consensus_gap_matrix,'2','B','*');

    matinit((char *) &ADS_consensus_gap_matrix,'3','#','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','*','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','A','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','C','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','G','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','T','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','N','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','X','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','M','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','R','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','W','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','S','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','Y','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','K','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','V','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','H','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','D','*');
    matinit((char *) &ADS_consensus_gap_matrix,'3','B','*');

    matinit((char *) &ADS_consensus_gap_matrix,'4','#','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','*','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','A','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','C','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','G','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','T','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','N','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','X','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','M','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','R','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','W','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','S','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','Y','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','K','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','V','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','H','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','D','*');
    matinit((char *) &ADS_consensus_gap_matrix,'4','B','*');

    ADS_c_matrix_valid=1;
  }

  return;
}


// Copy operator, needed by copy-constructor
AlignedDualSeq const & AlignedDualSeq::operator=(AlignedDualSeq const &other)
{
  FUNCSTART("AlignedDualSeq const & AlignedDualSeq::operator=(AlignedDualSeq const &other)");

  if(this != &other){
    discard();
    ADS_valid=other.ADS_valid;
    if(ADS_valid){
      ADS_miraparams=other.ADS_miraparams;
      ADS_len1 = other.ADS_len1 ;
      ADS_len2 = other.ADS_len2 ;
      ADS_contained= other.ADS_contained;
      ADS_score = other.ADS_score;
      ADS_expected_score = other.ADS_expected_score ;
      ADS_weight=other.ADS_weight;

      ADS_nummismatches=other.ADS_nummismatches;
      ADS_numgaps=other.ADS_numgaps;
      ADS_maxcontiguousgaps=other.ADS_maxcontiguousgaps;

      ADS_affinegapscore=other.ADS_affinegapscore;

      ADS_minbanddistance=other.ADS_minbanddistance;
      ADS_bandwidthused=other.ADS_bandwidthused;

      if(other.ADS_aligned_seq1!=nullptr){
	ADS_cursize=other.ADSF_total_len+1;
	ADS_aligned_seq1= new char[ADS_cursize];
	ADS_aligned_seq2= new char[ADS_cursize];
	ADS_consensus_seq= new char[ADS_cursize];
	ADS_consensus_gap_seq= new char[ADS_cursize];
	strcpy(ADS_aligned_seq1, other.ADS_aligned_seq1);
	strcpy(ADS_aligned_seq2, other.ADS_aligned_seq2);
	strcpy(ADS_consensus_seq, other.ADS_consensus_seq);
	strcpy(ADS_consensus_gap_seq, other.ADS_consensus_gap_seq);

	ADS_seq1=ADS_aligned_seq1;
	ADS_seq2=ADS_aligned_seq2+ADSF_delta;;
      }else{
	ADS_cursize=0;
	ADS_aligned_seq1=nullptr;
	ADS_aligned_seq2=nullptr;
	ADS_consensus_seq=nullptr;
	ADS_consensus_gap_seq=nullptr;
	ADS_seq1=nullptr;
	ADS_seq2=nullptr;
      }

    }
  }

  FUNCEND();
  return *this;
}


// Copy constructor
AlignedDualSeq::AlignedDualSeq(AlignedDualSeq const &other) : AlignedDualSeqFacts (other)
{
  FUNCSTART("AlignedDualSeq::AlignedDualSeq(AlignedDualSeq const &other)");

  ADS_valid=0;

  *this=other;

  FUNCEND();
}


// Destructor
AlignedDualSeq::~AlignedDualSeq()
{
  FUNCSTART("AlignedDualSeq::~AlignedDualSeq()");

  discard();              // free used memory

  FUNCEND();
}



// Discard
// Free used spaced and mark object as invalid (without destructing object)
void AlignedDualSeq::discard()
{
  FUNCSTART("void AlignedDualSeq::discard()");

  if(ADS_valid){
    if(ADS_aligned_seq1!=nullptr) delete [] ADS_aligned_seq1;
    ADS_aligned_seq1=nullptr;
    if(ADS_aligned_seq2!=nullptr) delete [] ADS_aligned_seq2;
    ADS_aligned_seq2=nullptr;
    if(ADS_consensus_seq!=nullptr) delete [] ADS_consensus_seq;
    ADS_consensus_seq=nullptr;
    if(ADS_consensus_gap_seq!=nullptr) delete [] ADS_consensus_gap_seq;
    ADS_consensus_gap_seq=nullptr;
  }
  ADS_cursize=0;

  FUNCEND();
}



// get the aligned sequences and store them while upcasing all
//  all characters. Only IUPACs are allowed for sequence characters,
//  blanks (' ') are endgaps. Nothing else allowed.
// The first sequence (seq1) always begins with a non endgap character.
// No errorchecking is done for endgaps within the sequences itself.
// Compute some classification numbers like score etc.
void AlignedDualSeq::acquireSequences(const char * seq1, const char * seq2, readid_t id1, readid_t id2, int8 id1dir, int8 id2dir,
				      bool enforce_clean_ends,
				      bool affinegapscore,
				      int32 bandwidthused, int32 minbanddistance)
{
  FUNCSTART("void AlignedDualSeq::acquireSequences(const char * seq1, const char * seq2, readid_t id1, readid_t id2, int8 id1dir, int8 id2dir, bool enforce_clean_ends, bool affinegapscore, int32 bandwidthused, int32 minbanddistance)");

  {
    auto tmplen=strlen(seq1);
    BUGIFTHROW(tmplen!=strlen(seq2)," Sequence lengths unequal: " << tmplen << " vs. " << strlen(seq2));
    BUGIFTHROW(tmplen==0, "Sequence lengths 0?");
    BUGIFTHROW(tmplen>65530, "Alignment of two sequences: " << tmplen << " > 65530 bases, cannot handle, shouldn't have happened!");
    ADSF_total_len= static_cast<uint16>(tmplen);
  }
  CEBUG("ADSF_total_len: " << ADSF_total_len << endl);

  ADS_minbanddistance=minbanddistance;
  ADS_bandwidthused=bandwidthused;

  if(ADS_aligned_seq1==nullptr
     || ADS_cursize<ADSF_total_len+1){
    discard();
    ADS_cursize=ADSF_total_len+1;                          // 0 character!
    ADS_aligned_seq1=new char[ADS_cursize];
    ADS_aligned_seq2=new char[ADS_cursize];
    ADS_consensus_seq=new char[ADS_cursize];
    ADS_consensus_gap_seq=new char[ADS_cursize];
  }

  ADSF_id1=id1;
  ADSF_id2=id2;
  ADSF_id1and2_directions=0;
  if(id1dir>0) ADSF_id1and2_directions=1;
  if(id2dir>0) ADSF_id1and2_directions|=2;
  // which sequence first?
  if(*seq1==' '){
    std::swap(seq1, seq2);
    ADSF_id1=id2;
    ADSF_id2=id1;
    //ADSF_id1_direction=id2dir;
    //ADSF_id2_direction=id1dir;

    ADSF_id1and2_directions=0;
    if(id2dir>0) ADSF_id1and2_directions=1;
    if(id1dir>0) ADSF_id1and2_directions|=2;
  }

  ADS_affinegapscore=affinegapscore;

  // copy sequences (upcase them) and compute some values
  uint32 delta_trigger=1;
  {
    const char * src1ptr=seq1;
    const char * src2ptr=seq2;
    char * dst1ptr=ADS_aligned_seq1;
    char * dst2ptr=ADS_aligned_seq2;
    for(uint32 runi=0; runi<ADSF_total_len; ++runi, ++src1ptr, ++src2ptr, ++dst1ptr, ++dst2ptr){
      uint32 overlap=1;

      CEBUG("" << *src1ptr << *src2ptr << " ");

      *dst1ptr=static_cast<char>(toupper(*src1ptr));
      switch(*dst1ptr){
      case ' ': {overlap=0; break;}       // no overlap if a ' '
      case '-':	{*dst1ptr='N'; break;}
      default: break;
      }
      CEBUG("!"<<  overlap);

      overlap++;
      *dst2ptr=static_cast<char>(toupper(*src2ptr));
      switch(*dst2ptr){
      case ' ': {overlap=0; break;}       // no overlap if a ' '
      case '-': {*dst2ptr='N'; }     // fall-through!
      default:{
	if(delta_trigger){
	  delta_trigger=0;
	  ADSF_delta=static_cast<uint16>(runi);
	}
	break;
      }
      }

      CEBUG(overlap);
    }
    *dst1ptr=0;                     // terminate sequence with 0
    *dst2ptr=0;                     // terminate sequence with 0
  }

  ADS_seq1=ADS_aligned_seq1;
  ADS_seq2=ADS_aligned_seq2+ADSF_delta;

  CEBUG("delta" << ADSF_delta);
  CEBUG("seq2" << *ADS_seq2<< "\n");

  // compute the sequence lengths and the right deltas
  {
    char * ptr= ADS_seq1;
    uint32 i=0;

    while(*ptr != 0 && *ptr!=' '){
      ++i; ++ptr;
    }
    ADS_len1=static_cast<uint16>(i);
    i=0;
    while(*ptr++ != 0) ++i;
    ADSF_id1_rightdelta=static_cast<uint16>(i);
  }

  {
    char * ptr= ADS_seq2;
    uint32 i=0;

    while(*ptr != 0 && *ptr!=' '){
      ++i; ++ptr;
    }
    ADS_len2=static_cast<uint16>(i);
    i=0;
    while(*ptr++ != 0) ++i;
    ADSF_id2_rightdelta=static_cast<uint16>(i);
  }

  CEBUG("len1" << ADS_len1);
  CEBUG("len2" << ADS_len2);

  // look if one of the alignements is contained in the other
  ADS_contained=0;
  if(ADSF_delta==0){
    if(ADS_len1==ADS_len2){
      ADS_contained=2;           // both alignments contain each other
    }else{
      // find out which sequence must be sequence 1
      if(ADS_len1<ADS_len2){
	std::swap(ADS_aligned_seq1,ADS_aligned_seq2);
	std::swap(ADS_seq1, ADS_seq2);
	std::swap(ADS_len1, ADS_len2);
	std::swap(ADSF_id1, ADSF_id2);
	std::swap(ADSF_id1_rightdelta, ADSF_id2_rightdelta);
	if(ADSF_id1and2_directions==1){
	  ADSF_id1and2_directions=2;
	}else if(ADSF_id1and2_directions==1){
	  ADSF_id1and2_directions=1;
	}
      }
      ADS_contained=1;
    }
  }else{
    if(ADS_len1==ADSF_total_len){
      ADS_contained=1;
    }
  }

  ADS_valid=1;

  // now compute the consensus first, then the score
  // (score needs consenus for simple gap evaluations)
  consensus();
  score(enforce_clean_ends, affinegapscore);

  if(ADSF_score_ratio>100
     || (ADSF_score_ratio==100
	 && (ADS_nummismatches>0 || ADS_numgaps>0))) ADSF_score_ratio=99;

  ADS_weight=ADS_score*ADSF_score_ratio*ADSF_score_ratio;

  // remaining ADSF_ entries ...
  if(ADS_numgaps+ADS_nummismatches>1023){
    ADSF_totalnonmatches=1023;
  }else{
    ADSF_totalnonmatches=static_cast<uint16>(ADS_numgaps+ADS_nummismatches);
  }
  // TODO: can this be deferred somehow until "really needed[tm]"?
  calcEndsLenContiguousMatch();


  //cout << *this;

  FUNCEND();
}

void AlignedDualSeq::matinit(char * arr, char a, char b, char value)
{
  a=static_cast<char>(toupper(a));
  b=static_cast<char>(toupper(b));
  arr[a+b*ADS_MATSIZE]=value;
  arr[b+a*ADS_MATSIZE]=value;
  a=static_cast<char>(tolower(a));
  arr[a+b*ADS_MATSIZE]=value;
  arr[b+a*ADS_MATSIZE]=value;
  b=static_cast<char>(tolower(b));
  arr[a+b*ADS_MATSIZE]=value;
  arr[b+a*ADS_MATSIZE]=value;
  a=static_cast<char>(toupper(a));
  arr[a+b*ADS_MATSIZE]=value;
  arr[b+a*ADS_MATSIZE]=value;
}




/*************************************************************************
 *
 * compute the score of the aligned sequences
 * uses scoring scheme defined in DynamicParams section
 * put in penalty for more than 1 star in a row
 *
 *************************************************************************/

//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif

int32 AlignedDualSeq::score(bool enforce_clean_ends, bool affinegapscore)
{
  FUNCSTART("int32 AlignedDualSeq::score(bool enforce_clean_ends, bool affinegapscore)");

  dynamic_parameters dyp= ADS_miraparams->getDynamicParams();
  align_parameters alpar= ADS_miraparams->getAlignParams();

  if(!ADS_valid || !ADS_s_matrix_valid){
    MIRANOTIFY(Notify::FATAL, " Object not initialised.");
  }

  double realscore=0.0;
  double expectedscore=0.0;
  int32 addpenpercent=0;

  uint32 badendcounter1=0;
  uint32 badendcounter2=0;

  ADS_maxcontiguousgaps=0;

  if(getOverlapLen()){
    CEBUG("OVERLAP SCORE CALC\n");
    char * ptr1=ADS_aligned_seq1+ADSF_delta;
    char * ptr2=ADS_aligned_seq2+ADSF_delta;
    char * ptrconsgap=ADS_consensus_gap_seq+ADSF_delta;
    uint32 starcounter=0;

    uint8 a,b;

    // this loop goes until oli==1 in a normal way, checking each character
    //  from the ads
    // I append an additional loop for oli==0, to get the case where
    //  one of the sequences ends in stars (can happen)
    //  For this last loop, the characters are set manually to ' '
    // Reason for the extra loop included in loop:
    //  I don't want the code for addpenpercent doubled after the loop
    for(int32 oli=getOverlapLen(); oli>=0;--oli, ++ptr1, ++ptr2, ++ptrconsgap){
      if(likely(oli>0)) {
	a=static_cast<uint8>(*ptr1);
	b=static_cast<uint8>(*ptr2);
      } else {
	a=' ';
	b=' ';
      }

      CEBUG(a << b << " ");

      if(enforce_clean_ends){
	if(unlikely(toupper(a) != toupper(b))) {
	  if(likely(dptools::isValidACGTStarBase(a) && dptools::isValidACGTStarBase(b))){
	    if(oli <alpar.ads_clean_end_distance) {
	      ++badendcounter1;
	      CEBUG("BADEND1!!! ");
	    }
	    if(getOverlapLen()-oli < alpar.ads_clean_end_distance){
	      ++badendcounter2;
	      CEBUG("BADEND2!!! ");
	    }
	  }
	}
      }

      if(*ptrconsgap=='*'){
	if (affinegapscore) {
	  if (starcounter == 0) {
	    realscore += dyp.dyn_score_gap;
	  } else {
	    //realscore += dyp.dyn_score_mismatch;
	    // no bonus or malus
	  }
	} else {
	  // if starcounter < size of ADS_powofstarcounter, add to real score
	  // if not, the value would be 0 anyway
	  // This makes a increasing penalty for inreasing gap lengths
	  if(starcounter<ADS_powofstarcounter.size()){
	    realscore+=ADS_realscore_matrix[a][b]*ADS_powofstarcounter[starcounter];
	  }
	}
	expectedscore+=ADS_expectedscore_matrix[a][b];
	starcounter++;
      }else{
	if(starcounter>0){
	  if(starcounter>ADS_maxcontiguousgaps) ADS_maxcontiguousgaps=static_cast<uint16>(starcounter);
	  if(starcounter>=alpar.ads_gp_function.size()){
	    addpenpercent+=alpar.ads_gp_function[alpar.ads_gp_function.size()-1];
	  }else{
	    addpenpercent+=alpar.ads_gp_function[starcounter];
	  }
	}
	starcounter=0;
	realscore+=ADS_realscore_matrix[a][b];
	expectedscore+=ADS_expectedscore_matrix[a][b];
      }
      CEBUG("\tes: " << expectedscore << "\trs: " << realscore<< "\n");
    }
  }

  if(addpenpercent>static_cast<int32>(alpar.ads_max_gppercent)){
    addpenpercent=alpar.ads_max_gppercent;
  }

  CEBUG("\nExpected score " << expectedscore << endl);
  CEBUG("Real score " << realscore << endl);
  CEBUG("Max contiguous gaps: " << ADS_maxcontiguousgaps << endl);
  CEBUG("Clean end killer: " << enforce_clean_ends << " " << alpar.ads_clean_end_distance << " " << alpar.ads_clean_end_mismatchallowed << endl);
  CEBUG("Additional penalty (%) " << addpenpercent << endl);

  if(alpar.ads_extra_gap_penalty){
    realscore-=realscore*addpenpercent/100;
    CEBUG("Adjusted real score " << realscore << endl);
  }

  if(realscore<0.0){
    realscore=0.0;
  }

  if(badendcounter1>alpar.ads_clean_end_mismatchallowed
     || badendcounter2>alpar.ads_clean_end_mismatchallowed){
    CEBUG("Killer: " << badendcounter1 << ' ' << badendcounter2 << '\n');
    realscore=0.0;
  }

  ADS_score=static_cast<int32>(realscore/dyp.dyn_score_multiplier);
  ADS_expected_score=static_cast<int32>(expectedscore/dyp.dyn_score_multiplier);
  if(expectedscore>0 && realscore>0){
    ADSF_score_ratio=
      static_cast<int8>((100.0 / expectedscore * realscore)+0.5);
  }else{
    ADSF_score_ratio=0;
  }

  CEBUG(*this);

  FUNCEND();

  return static_cast<int32>(realscore);
}


/*************************************************************************
 *
 * Computes two kinds of consensus: one where gap against a base gives a
 *  base, once where it gives a gap
 *
 * Additionally, counts how many mismatches or real gaps are present and
 *  stores those numbers
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void AlignedDualSeq::consensus()
{
  FUNCSTART("void AlignedDualSeq::consensus()");

  char *src1ptr=ADS_aligned_seq1;
  char *src2ptr=ADS_aligned_seq2;
  char *dstptr1=ADS_consensus_seq;
  char *dstptr2=ADS_consensus_gap_seq;
  ADS_nummismatches=0;
  ADS_numgaps=0;
  for(uint32 i=0; i<ADSF_total_len; i++, src1ptr++, src2ptr++, dstptr1++, dstptr2++){
    // generating the two consensus sequences

    CEBUG(i << "  " << *src1ptr << " with " << (uint16) *src2ptr << " is ");
    CEBUG((uint16) ADS_consensus_matrix[static_cast<uint8>(*src1ptr)][static_cast<uint8>(*src2ptr)]);
    CEBUG(endl);
    *dstptr1=ADS_consensus_matrix[static_cast<uint8>(*src1ptr)][static_cast<uint8>(*src2ptr)];
    *dstptr2=ADS_consensus_gap_matrix[static_cast<uint8>(*src1ptr)][static_cast<uint8>(*src2ptr)];

    // counting mismatches and gaps
    // count errors only if none is a blank (endgap)
    if(*src1ptr!=' ' && *src2ptr!=' '){
      // if the consensus gap result is not [ACGTNX], we need to seriously
      //  check for mismatches
      if(!dptools::isValidBase(*dstptr2)){
	CEBUG("Check for mismatch: ");
	// first, check whether we have a gap situation
	bool foundreason=false;
	if(*src1ptr=='*' || *src1ptr=='#'
	   || *src2ptr=='*' || *src2ptr=='#') {

	  CEBUG("gap? ");

	  // if one of the two is a IUPAC, we have a gap/base situation
	  //  and must count a gap
	  if(dptools::isValidIUPACBase(*src1ptr)
	     || dptools::isValidIUPACBase(*src2ptr)){
	    CEBUG("yes!\n");
	    foundreason=true;
	    ADS_numgaps++;
	  }else if((*src1ptr=='*' || *src1ptr=='#')
		   && (*src2ptr=='*' || *src2ptr=='#')) {
	    // a gap/gap (or oldgap) which we don't count as error
	    foundreason=true;
	  }
	}
	if(!foundreason){
	  CEBUG("no! mismatch. ");
	  // so, we have a potential mismatch situation without gaps
	  //
	  // if both are valid bases, it's a normal mismatch involving
	  //  [ACGT] vs another [ACGT] (in a [ACGT vs [NX] situation,
	  //  we wouldn't be in this branch anyway
	  if(dptools::isValidBase(*src1ptr) && dptools::isValidBase(*src2ptr)){
	    CEBUG("normal base/base mismatch.\n");
	    ADS_nummismatches++;
	  } else {
	    // we are in a IUPAC vs IUPAC situation
	    if(dptools::isValidBase(*src1ptr)){
	      CEBUG("src1 base. ");
	      // valid base vs. IUPAC
	      // count mismatch if the base is not contained in the IUPAC
	      if(!dptools::hasNucleicAcidInIUPAC(*src1ptr,*src2ptr)){
		CEBUG("src2 IUPAC mismatch!");
		ADS_nummismatches++;
	      }
	      CEBUG('\n');
	    }else if(dptools::isValidBase(*src2ptr)){
	      CEBUG("src2 base. ");
	      // valid base vs. IUPAC
	      // count mismatch if the base is not contained in the IUPAC
	      if(!dptools::hasNucleicAcidInIUPAC(*src2ptr,*src1ptr)){
		CEBUG("src1 IUPAC mismatch!");
		ADS_nummismatches++;
	      }
	      CEBUG('\n');
	    }else{
	      CEBUG("IUPAC vs IUPAC. ");
	      // IUPAC vs. IUPAC
	      // if they're unequal, it's a mismatch
	      // (we don't check further on partially equality or the like)
	      if(*src1ptr != *src2ptr){
		CEBUG("mismatch. ");
		ADS_nummismatches++;
	      }
	      CEBUG('\n');
	    }
	  }
	}
      }
    }
  }
  *dstptr1=0;
  *dstptr2=0;

  CEBUG("Gaps: " << ADS_numgaps);
  CEBUG("\nMismatches: " << ADS_nummismatches << '\n');

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

// dump some values to stream
std::ostream & operator<<(std::ostream & ostr, const AlignedDualSeq &ads)
{
  FUNCSTART("void AlignedDualSeq::dump()");

  if(!ads.ADS_valid){
    ostr << "Invalid / Not initialised.\n";
  }else{
    ostr << "Aligned (1) Seq1 and (2) Seq2 with (3) difference monitor and (4) consensus:\n";
    if(ads.ADS_expected_score!=ads.ADS_score){
      ostr << "Not equal.\n";
    }
    if( ads.ADS_aligned_seq1 != nullptr
        && ads.ADS_aligned_seq2 != nullptr
	&& ads.ADS_consensus_seq != nullptr
        && ads.ADS_consensus_gap_seq != nullptr){
      const uint32 cpl=60;
      uint32 cut =cpl;

      bool addtab=false;
      if(ads.ADSF_id1 >999 || ads.ADSF_id2 >999) addtab=true;

      char *s1ptr= ads.ADS_aligned_seq1;
      char *s2ptr= ads.ADS_aligned_seq2;
      char *cptr= ads.ADS_consensus_seq;
      char *cgptr= ads.ADS_consensus_gap_seq;
      for(uint32 i=0; i<ads.ADSF_total_len; i+=cut, s1ptr+=cut, s2ptr+=cut, cptr+=cut, cgptr+=cut){
	if(i+cut>ads.ADSF_total_len) cut=ads.ADSF_total_len%cpl;
	char savetmp;

	if(addtab) ostr << "\t";
	ostr << "\t|    .    |    .    |    .    |    .    |    .    |    .    ";
	ostr << endl;
	ostr << "ID1:" << ads.ADSF_id1 << "\t";
	if(addtab && ads.ADSF_id1 <1000) ostr << "\t";
	savetmp=*(s1ptr+cut);
	*(s1ptr+cut)=0;
	ostr << s1ptr << endl;
	*(s1ptr+cut)=savetmp;

	ostr << "ID2:" << ads.ADSF_id2 << "\t";
	if(addtab && ads.ADSF_id2 <1000) ostr << "\t";
	savetmp=*(s2ptr+cut);
	*(s2ptr+cut)=0;
	ostr << s2ptr << endl;
	*(s2ptr+cut)=savetmp;

	ostr << i << "-\t";
	if(addtab) ostr << "\t";

	for(uint32 j=0; j<cut; j++){
	  char a=*(s1ptr+j);
	  char b=*(s2ptr+j);
	  if(a==b || a==' ' || b==' '){
	    ostr << " ";
	  }else{
	    if(a=='*' || b=='*'){
	      if(a=='#' || b=='#'){
		ostr << " ";
	      }else{
		ostr << "*";
	      }
	    }else{
	      if(a=='N' || b=='N' || a=='X' || b=='X'){
		ostr << "!";
	      }else{
		ostr << "X";
	      }
	    }
	  }
	}
	ostr << endl;

	ostr << i << "C\t";
	if(addtab) ostr << "\t";

	savetmp=*(cptr+cut);
	*(cptr+cut)=0;
	ostr << cptr << "\n";
	*(cptr+cut)=savetmp;

	ostr << i << "CG\t";
	if(addtab) ostr << "\t";

	savetmp=*(cgptr+cut);
	*(cgptr+cut)=0;
	ostr << cgptr << "\n\n";
	*(cgptr+cut)=savetmp;

      }
    }else{
      cout << "ADS object apparently in saveMem mode, cannot show alignment.\n";
    }

    ostr << "Len1: " << ads.ADS_len1 << '\n';
    ostr << "Len2: " << ads.ADS_len2 << '\n';

    ostr << static_cast<AlignedDualSeqFacts>(ads);
    ostr << "Expected Score: " << ads.ADS_expected_score << '\n';
    ostr << "Score: " << ads.ADS_score << '\n';
    ostr << "Weight: " << ads.ADS_weight << '\n';
    ostr << "Contained: " << static_cast<uint16>(ads.ADS_contained) << endl;

    ostr << "Mismatches: " << ads.ADS_nummismatches << '\n';
    ostr << "Gaps: " << ads.ADS_numgaps << '\n';
    ostr << "Max contiguous gaps: " << ads.ADS_maxcontiguousgaps << '\n';
    ostr << "Affine gap score: " << ads.ADS_affinegapscore << '\n';
  }

  FUNCEND();
  return ostr;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

const char * AlignedDualSeq::getAlignedSequence(readid_t id) const
{
  FUNCSTART("char * getAlignedSequence(readid_t id)");
  if(!ADS_valid){
    MIRANOTIFY(Notify::FATAL, " Object not initialised.");
  }
  if(ADS_aligned_seq1==nullptr || ADS_aligned_seq2==nullptr){
    MIRANOTIFY(Notify::FATAL, " Tried to get aligned sequences on object where saveMem() had been called previously!");
  }
  if(id==ADSF_id1){
    FUNCEND();
    return ADS_seq1;
  }else if(id==ADSF_id2){
    FUNCEND();
    return ADS_seq2;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

const char * AlignedDualSeq::getSequenceAtOverlapStart(readid_t id) const
{
  FUNCSTART("char * getSequenceAtOverlapStart(readid_t id)");
  if(!ADS_valid){
    MIRANOTIFY(Notify::FATAL, " Object not initialised.");
  }
  if(ADS_aligned_seq1==nullptr || ADS_aligned_seq2==nullptr){
    MIRANOTIFY(Notify::FATAL, " Tried to get aligned sequences on object where saveMem() had been called previously!");
  }
  if(id==ADSF_id1){
    FUNCEND();
    return ADS_seq1+ADSF_delta;
  }else if(id==ADSF_id2){
    FUNCEND();
    return ADS_seq2;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

uint32 AlignedDualSeq::getLenOfAlignedSequence(readid_t id) const
{
  FUNCSTART("uint32 getLenOfAlignedSequence(readid_t id)");
  if(!ADS_valid){
    MIRANOTIFY(Notify::FATAL, " Object not initialised.");
  }
  if(id==ADSF_id1){
    FUNCEND();
    return ADS_len1;
  }else if(id==ADSF_id2){
    FUNCEND();
    return ADS_len2;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool AlignedDualSeq::isContained(readid_t id) const
{
  if(!ADS_contained) return false;
  if(ADS_contained==2) return true;
  if(id==ADSF_id2) return true;
  return false;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUG(bla)
void AlignedDualSeq::calcEndsLenContiguousMatch()
{
  FUNCSTART("uint32 AlignedDualSeq::calcEndsLenContiguousMatch()");
  if(!ADS_valid){
    MIRANOTIFY(Notify::FATAL, " Object not initialised.");
  }
  BUGIFTHROW(ADS_aligned_seq1==nullptr,"ADS_aligned_seq1==nullptr ???");
  BUGIFTHROW(ADS_aligned_seq2==nullptr,"ADS_aligned_seq2==nullptr ???");

  CEBUG("CCCCCCCCCCCCCCCCCC\n");

  int32 local_5pconmatch1=0;
  int32 local_3pconmatch1=0;
  int32 local_5pconmatch2=0;
  int32 local_3pconmatch2=0;

  // forward
  auto as1=ADS_aligned_seq1;
  auto as2=ADS_aligned_seq2;
  bool b1=true; // seq 1 is clean from blank
  bool b2=true; // seq 2 is clean from blank
  int32 runl=0;
  CEBUG("FORWARD S\n");
  for(; (runl<ADS_cursize && (*as1==' ' || *as2==' ')); ++as1, ++as2, ++runl){
    if(*as1!=' ') b1=false;
    if(*as2!=' ') b2=false;
    CEBUG(runl << '\t' << *as1 << '\t' << *as2 << endl);
  }
  CEBUG("FORWARD B " << b1 << ' ' << b2 << endl);
  {
    auto left=runl;
    CEBUG("FORWARD M\n");
    for(; runl<ADS_cursize && toupper(*as1)==toupper(*as2); ++as1, ++as2, ++runl){
      CEBUG(runl << '\t' << *as1 << '\t' << *as2 << endl);
    }
    if(b1) local_5pconmatch1=runl-left;
    if(b2) local_5pconmatch2=runl-left;
    CEBUG("local_5pconmatch1 " << local_5pconmatch1 << endl);
    CEBUG("local_5pconmatch2 " << local_5pconmatch2 << endl);
  }

  // backward
  b1=true; // seq 1 is clean from blank
  b2=true; // seq 2 is clean from blank
  as1=ADS_aligned_seq1+ADSF_total_len-1; // there's the final 0 to be accounted for
  as2=ADS_aligned_seq2+ADSF_total_len-1; // there's the final 0 to be accounted for
  runl=ADS_cursize;
  CEBUG("BACKWARD S\n");
  for(; (runl>0 && (*as1==' ' || *as2==' ')); --as1, --as2, --runl){
    if(*as1!=' ') b1=false;
    if(*as2!=' ') b2=false;
    CEBUG(runl << '\t' << *as1 << '\t' << *as2 << endl);
  }
  CEBUG("BACKWARD B " << b1 << ' ' << b2 << endl);
  {
    auto right=runl;
    CEBUG("BACKWARD M\n");
    for(; runl>0 && toupper(*as1)==toupper(*as2); --as1, --as2, --runl){
      CEBUG(runl << '\t' << *as1 << '\t' << *as2 << endl);
    }
    if(b1) local_3pconmatch1=right-runl;
    if(b2) local_3pconmatch2=right-runl;
    CEBUG("local_3pconmatch1 " << local_3pconmatch1 << endl);
    CEBUG("local_3pconmatch2 " << local_3pconmatch2 << endl);
  }

  // swap if needed
  if(getSequenceDirection(ADSF_id1)<0) std::swap(local_5pconmatch1,local_3pconmatch1);
  if(getSequenceDirection(ADSF_id2)<0) std::swap(local_5pconmatch2,local_3pconmatch2);

  // fill in ADSF_... pendants

  if(local_5pconmatch1>=28){
    ADSF_5pconmatch1=7;
  }else{
    ADSF_5pconmatch1=static_cast<uint16>(local_5pconmatch1/4);
  }
  if(local_3pconmatch1>=28){
    ADSF_3pconmatch1=7;
  }else{
    ADSF_3pconmatch1=static_cast<uint16>(local_3pconmatch1/4);
  }
  if(local_5pconmatch2>=28){
    ADSF_5pconmatch2=7;
  }else{
    ADSF_5pconmatch2=static_cast<uint16>(local_5pconmatch2/4);
  }
  if(local_3pconmatch2>=28){
    ADSF_3pconmatch2=7;
  }else{
    ADSF_3pconmatch2=static_cast<uint16>(local_3pconmatch2/4);
  }

  CEBUG("ADSF_5pconmatch1 " << ADSF_5pconmatch1 << endl);
  CEBUG("ADSF_3pconmatch1 " << ADSF_3pconmatch1 << endl);
  CEBUG("ADSF_5pconmatch2 " << ADSF_5pconmatch2 << endl);
  CEBUG("ADSF_3pconmatch2 " << ADSF_3pconmatch2 << endl);


  return;
}
#define CEBUG(bla)



/*************************************************************************
 *
 * Save memory by deleting non-necessary info once an alignment is computed
 *
 *
 *************************************************************************/

void AlignedDualSeq::saveMem(bool delete_seq, bool delete_consseq)
{
  FUNCSTART("void AlignedDualSeq::saveMem(bool delete_seq, bool delete_consseq)");

  if(!ADS_valid){
    MIRANOTIFY(Notify::FATAL, " Object not initialised.");
  }

  if(delete_seq){
    if(ADS_aligned_seq1!=nullptr) delete [] ADS_aligned_seq1;
    ADS_aligned_seq1=nullptr;
    if(ADS_aligned_seq2!=nullptr) delete [] ADS_aligned_seq2;
    ADS_aligned_seq2=nullptr;

    // avoid dangling pointers
    ADS_seq1=nullptr;
    ADS_seq2=nullptr;
  }
  if(delete_consseq){
    if(ADS_consensus_seq!=nullptr) delete [] ADS_consensus_seq;
    ADS_consensus_seq=nullptr;
    if(ADS_consensus_gap_seq!=nullptr) delete [] ADS_consensus_gap_seq;
    ADS_consensus_gap_seq=nullptr;
  }

  FUNCEND();
}



/*************************************************************************
 *
 * Starting from two aligned sequences, computes the maximum good
 *  length of these sequences in the alignment
 *
 * In a given window, a maxiumum number of allowed errors may be
 *  encountered before the search algorithm stops.
 *
 * return length for each sequence
 *
 * currently used in extendADS()
 *
 *************************************************************************/

bool AlignedDualSeq::clipper(uint32 winlen, int32 allerr, int32 & retgoodlen1, int32 & retgoodlen2) const
{
  FUNCSTART("void AlignedDualSeq::clipper(int32 winlen, int32 numerr, int32 retries, int32 & retgoodlen1, int32 & retgoodlen2) const");

  BUGIFTHROW(!ADS_valid, " Object not initialised.");

  BUGIFTHROW((ADS_aligned_seq1==nullptr || ADS_aligned_seq2==nullptr)," Tried to use function on object where saveMem() had been called previously!");

  BUGIFTHROW(winlen==0,"winlen == 0 ?");

  if(getOverlapLen()<=winlen){
    // nicht verlängerbar;
    //cout << "too short";
    return false;
    FUNCEND();
  }

  //  cout << "hipping";

  dynamic_parameters dyp= ADS_miraparams->getDynamicParams();


  // Precompute initial window
  int32 errors=0;
  char * s1r=ADS_aligned_seq1+ADSF_delta;
  char * s2r=ADS_aligned_seq2+ADSF_delta;
  char * s1l=s1r;
  char * s2l=s2r;
  uint32 nORx=0;
  for(uint32 i=0; i<winlen-1; i++, s1r++, s2r++){
    int32 comp=ADS_realscore_matrix[static_cast<uint8>(*s1r)][static_cast<uint8>(*s2r)];
    //    cout << *s1r << *s2r << endl;
    if(comp == dyp.dyn_score_mismatch
       || comp == dyp.dyn_score_gap){
      errors++;
    }
    if(toupper(*s1r)=='N' || toupper(*s1r)=='X'
       || toupper(*s2r)=='N' || toupper(*s2r)=='X') {
	nORx++;
    }else{
	nORx=0;
    }
  }


  // now slide window over alignment

  int32 acceptedoverlap=0;
  int32 lastnoerror=0;
  for(; acceptedoverlap < static_cast<int32>(getOverlapLen())-static_cast<int32>(winlen); acceptedoverlap++, s1l++, s2l++, s1r++, s2r++){
    int32 comp=ADS_realscore_matrix[static_cast<uint8>(*s1r)][static_cast<uint8>(*s2r)];
    if(comp == dyp.dyn_score_mismatch
       || comp == dyp.dyn_score_gap){
      errors++;
    }

    if(toupper(*s1r)=='N' || toupper(*s1r)=='X'
       || toupper(*s2r)=='N' || toupper(*s2r)=='X') {
	nORx++;
    }else{
	nORx=0;

	if(errors==0) lastnoerror=acceptedoverlap;
    }

    //cout << *s1r << *s2r << " " <<  errors << " " << nORx << endl;

    if(nORx==winlen*2) break;

    if(errors>allerr) break;

    comp=ADS_realscore_matrix[static_cast<uint8>(*s1l)][static_cast<uint8>(*s2l)];
    if(comp == dyp.dyn_score_mismatch
       || comp == dyp.dyn_score_gap){
      errors--;
    }
    //    cout << "\t"<<*s1l << *s2l << errors<< endl;
  }

  //acceptedoverlap-=nORx;
  //acceptedoverlap-=5;
  //
  //acceptedoverlap--;
  //if(acceptedoverlap<5) {
  //  FUNCEND();
  //  return false;
  //}
  //
  //retgoodlen2=acceptedoverlap+winlen;
  //retgoodlen1=ADSF_delta+acceptedoverlap;

  if(lastnoerror<5) {
    FUNCEND();
    return false;
  }

  retgoodlen2=lastnoerror+winlen;
  retgoodlen1=ADSF_delta+lastnoerror+winlen;

  FUNCEND();

  return true;
}





ADSEstimator::ADSEstimator()
{
  ADSE_id1=0xffffffff;
  ADSE_id2=0xffffffff;
  ADSE_len1=0;
  ADSE_len2=0;
  ADSE_lexpandof1=0;
  ADSE_rexpandof1=0;
  ADSE_lexpandof2=0;
  ADSE_rexpandof2=0;
  ADSE_dir1=0;
  ADSE_dir2=0;
}



/*************************************************************************
 *
 * This class is not for storage, but meant as simple helper to estimate
 *  different values of an alignment of two sequences.
 *
 * e.g.  len1 = 13; len2 = 12
 *
 *  s1  ------------>               s1      ------------>
 *  s2      ----------->	    s2  ----------->
 *
 * init with offsets1tos2 = 4	   init with offsets1tos2 = -4
 *
 * lexpandof1=0			   lexpandof1=4
 * rexpandof1=3			   rexpandof1=0
 * lexpandof2=4			   lexpandof2=0
 * rexpandof2=0                    rexpandof2=5
 *
 *
 *  s1  ------------>               s1      ------------>
 *  s2      <-----------	    s2  <-----------
 *
 * init with offsets1tos2 = 4      init with offsets1tos2 = -4
 *
 * lexpandof1=0			   lexpandof1=4
 * rexpandof1=3			   rexpandof1=0
 * lexpandof2=0			   lexpandof2=5
 * rexpandof2=4			   rexpandof2=0
 *
 *************************************************************************/

// handling of s1 in - and s2 in + is more a makeshift ... improve when time

void ADSEstimator::calcNewEstimate(int32 offsets1tos2, uint32 len1, uint32 len2, readid_t id1, readid_t id2, int8 id1dir, int8 id2dir)
{
  FUNCSTART("void ADSEstimator::calcNewEstimate(int32 offsets1tos2, uint32 len1, uint32 len2, readid_t id1, readid_t id2, int8 id1dir, int8 id2dir)");

  BUGIFTHROW(len1==0,"len1==0 ? id1= " << id1 << " id2= " << id2);
  BUGIFTHROW(len2==0,"len2==0 ? id2= " << id2 << " id1= " << id1);

  if(id1dir<0 && id2dir<0){
    id1dir=1;
    id2dir=1;
  }

  bool mustswap=false;

  ADSE_len1=len1;
  ADSE_len2=len2;
  ADSE_id1=id1;
  ADSE_id2=id2;
  ADSE_dir1=id1dir;
  ADSE_dir2=id2dir;

  if(id1dir<0){
    mustswap=true;
  }

  uint32 absoffsets1tos2=abs(offsets1tos2);

  if(id2dir>0) {
    if(offsets1tos2 >= 0) {
//cout << "adse 1\n";
      ADSE_lexpandof1=0;
      ADSE_lexpandof2=absoffsets1tos2;
      ADSE_totallen=absoffsets1tos2+len2;
      if(ADSE_totallen<len1) ADSE_totallen=len1;
      ADSE_rexpandof1=ADSE_totallen-len1;
      ADSE_rexpandof2=ADSE_totallen-absoffsets1tos2-len2;
      ADSE_overlaplen=ADSE_totallen-ADSE_lexpandof2-ADSE_rexpandof1-ADSE_rexpandof2;
    }else{
//cout << "adse 2\n";
      ADSE_lexpandof1=absoffsets1tos2;
      ADSE_lexpandof2=0;
      ADSE_totallen=absoffsets1tos2+len1;
      if(ADSE_totallen<len2) ADSE_totallen=len2;
      ADSE_rexpandof1=ADSE_totallen-absoffsets1tos2-len1;
      ADSE_rexpandof2=ADSE_totallen-len2;
      ADSE_overlaplen=ADSE_totallen-ADSE_lexpandof1-ADSE_rexpandof1-ADSE_rexpandof2;
    }
  }else{
    if(offsets1tos2 >= 0) {
//cout << "adse 3\n";
      ADSE_lexpandof1=0;
      ADSE_rexpandof2=absoffsets1tos2;
      ADSE_totallen=absoffsets1tos2+len2;
      if(ADSE_totallen<len1) ADSE_totallen=len1;
      ADSE_rexpandof1=ADSE_totallen-len1;
      ADSE_lexpandof2=ADSE_totallen-absoffsets1tos2-len2;
      ADSE_overlaplen=ADSE_totallen-ADSE_rexpandof2-ADSE_lexpandof2-ADSE_rexpandof1;
    }else{
//cout << "adse 4\n";
      ADSE_lexpandof1=absoffsets1tos2;
      ADSE_rexpandof2=0;
      ADSE_totallen=absoffsets1tos2+len1;
      if(ADSE_totallen<len2) ADSE_totallen=len2;
      ADSE_rexpandof1=ADSE_totallen-absoffsets1tos2-len1;
      ADSE_lexpandof2=ADSE_totallen-len2;
      ADSE_overlaplen=ADSE_totallen-ADSE_lexpandof1-ADSE_lexpandof2-ADSE_rexpandof1;
    }
  }

  // is this necessary?
  if(mustswap){
    std::swap(ADSE_id1,ADSE_id2);
    std::swap(ADSE_len1,ADSE_len2);
    std::swap(ADSE_dir1,ADSE_dir2);
    std::swap(ADSE_lexpandof1,ADSE_lexpandof2);
    std::swap(ADSE_rexpandof1,ADSE_rexpandof2);
    std::swap(ADSE_rexpandof2,ADSE_lexpandof2);
  }

  if(ADSE_overlaplen>1073741824
     || ADSE_overlaplen>ADSE_totallen
     || (offsets1tos2>=0 && offsets1tos2>=len1)
     || (offsets1tos2<0 && abs(offsets1tos2)>=len2)){
    cout << "something is really wrong with that ADSE!"
	 << "\n" << offsets1tos2
	 << "\n" << len1
	 << "\n" << len2
	 << "\n" << id1
	 << "\n" << id2
	 << "\n" << static_cast<int16>(id1dir)
	 << "\n" << static_cast<int16>(id2dir)
	 << "\n" << *this << endl;
    BUGIFTHROW(offsets1tos2>=0 && offsets1tos2>=len1,"offsets1tos2>=0 && offsets1tos2>=len1 ?");
    BUGIFTHROW(offsets1tos2<0 && abs(offsets1tos2)>=len2,"offsets1tos2<0 && abs(offsets1tos2)>=len2 ?");
    BUGIFTHROW(ADSE_overlaplen>ADSE_totallen,"ADSE_overlaplen>ADSE_totallen");
    BUGIFTHROW(true,"Some other reason ...");
  }

  return;
}


// like above, but here expect values as deliver by skim
// i.e., if one of the two direction is rever, then the
//  eoffset is "left" only when id1<id2 && dir2<0
//                           or id1>id2 && dir1<0
// else it's the *right* offset ... in which case one must
//  swap a few values
void ADSEstimator::calcNewEstimateFromSkim(int32 offsets1tos2, uint32 len1, uint32 len2, readid_t id1, readid_t id2, int8 id1dir, int8 id2dir)
{
  bool itsright=false;
  if(id1dir*id2dir<0){
    if((id1<id2 && id1dir<0)
       || (id1>id2 && id2dir<0)) itsright=true;
  }
  calcNewEstimate(offsets1tos2,len1,len2,id1,id2,id1dir,id2dir);
  if(itsright){
    std::swap(ADSE_lexpandof1,ADSE_rexpandof1);
    std::swap(ADSE_lexpandof2,ADSE_rexpandof2);
  }
}


uint32 ADSEstimator::getEstimatedLeftExpand(readid_t id) const
{
  FUNCSTART("uint32 ADSEstimator::getEstimatedLeftExpand(readid_t id) const");

  if(ADSE_id1==id){
    return ADSE_lexpandof1;
  }else if(ADSE_id2==id){
    return ADSE_lexpandof2;
  }
  MIRANOTIFY(Notify::INTERNAL,"id " << id << " not in ADSEstimator.\n");
}


uint32 ADSEstimator::getEstimatedRightExpand(readid_t id) const
{
  FUNCSTART("uint32 ADSEstimator::getEstimatedRightExpand(readid_t id) const");

  if(ADSE_id1==id){
    return ADSE_rexpandof1;
  }else if(ADSE_id2==id){
    return ADSE_rexpandof2;
  }
  MIRANOTIFY(Notify::INTERNAL,"id " << id << " not in ADSEstimator.\n");
}

uint32 ADSEstimator::getLen(readid_t id) const
{
  FUNCSTART("uint32 ADSEstimator::getLen(readid_t id) const");

  if(ADSE_id1==id){
    return ADSE_len1;
  }else if(ADSE_id2==id){
    return ADSE_len2;
  }
  MIRANOTIFY(Notify::INTERNAL,"id " << id << " not in ADSEstimator.\n");
}


uint32 ADSEstimator::getContainmentLevel() const
{
  FUNCSTART("uint32 ADSEstimator::getContainmentLevel() const");

  bool s1isexpanded=true;
  bool s2isexpanded=true;
  if(ADSE_lexpandof1 == 0
     && ADSE_rexpandof1 ==0){
    s1isexpanded=false;
  }
  if(ADSE_lexpandof2 == 0
     && ADSE_rexpandof2 ==0){
    s2isexpanded=false;
  }
  if(s1isexpanded) {
    if(s2isexpanded) return 0;
    return 1;
  }else if(s2isexpanded) {
    return 1;
  }
  return 2;
}

// returns id of sequence that is completely contained in other
// if none is contained, return -1
// if both sequences are contained in each other, always returns id2
int32 ADSEstimator::getIDOfContained() const
{
  FUNCSTART("int32 ADSEstimator::getIDOfContained() const");

  if(ADSE_lexpandof1 == 0
     && ADSE_rexpandof1 ==0){
    return ADSE_id2;
  } else if(ADSE_lexpandof2 == 0
     && ADSE_rexpandof2 ==0){
    return ADSE_id1;
  }

  return -1;
}


void ADSEstimator::getPositionsInForwardSequenceofAlignmentStart(readid_t id1, readid_t id2, int32 & offset1, int32 & offset2) const
{
  FUNCSTART("void ADSEstimator::getStartOffsetPositions(readid_t id1, readid_t id2, int32 & offset1, int32 & offset2) const");

  bool needswap=false;
  if(id1==ADSE_id2){
    std::swap(id1, id2);
    needswap=true;
  }

  if(ADSE_lexpandof1 == 0){
    if(ADSE_dir2>0){
      offset1=ADSE_lexpandof2;
      offset2=0;
    }else{
      offset1=ADSE_rexpandof2;
      offset2=ADSE_len2;
    }
  }else{
    if(ADSE_dir2>0){
      offset1=0;
      offset2=ADSE_lexpandof1;
    }else{
      offset1=0;
      offset2=ADSE_len2-ADSE_lexpandof1;
    }
  }

  if(needswap) std::swap(offset1,offset2);

  return;
}
