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


#include "dynamic.H"

#include <mira/parameters.H>

#include <climits>

using std::cout;
using std::cerr;
using std::endl;


//#define CEBUGFLAG
#define CLOCK_STEPS1

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif



#if __GNUC__ >= 3
#define prefetchrs(p)     __builtin_prefetch((p), 0, 0)
#define prefetchws(p)     __builtin_prefetch((p), 1, 0)
#define prefetchrl(p)     __builtin_prefetch((p), 0, 3)
#define prefetchwl(p)     __builtin_prefetch((p), 1, 3)
#else
#define prefetchrs(p)
#define prefetchws(p)
#define prefetchrl(p)
#define prefetchwl(p)
#endif


uint64 Dynamic::DYN_alloccounts=0;
uint64 Dynamic::DYN_alloccountm=0;
int16 Dynamic::DYN_matvalid=0;
int32 Dynamic::DYN_match_matrix[DYN_MATSIZE][DYN_MATSIZE];


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Constructor with init
Dynamic::Dynamic(MIRAParameters * params, const char * seq1, uint32 len1, const char * seq2, uint32 len2, bool calcwithoffset, int32 expectedoffset)
{
  FUNCSTART("Dynamic::Dynamic(MIRAParameters * params, const char * seq1, uint32 len1, const char * seq2, uint32 len2, bool calcwithoffset, int32 expectedoffset)");

  DYN_miraparams=params;

  init();

  setSequences(seq1, len1, seq2, len2, calcwithoffset, expectedoffset);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Constructor without direct init.
// Before using this instance later, one must call setSequences or else
//  it will throw an error
Dynamic::Dynamic(MIRAParameters * params)
{
  FUNCSTART("Dynamic::Dynamic(MIRAParameters * params)");

  DYN_miraparams=params;

  init();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Dynamic::zeroVars()
{
  FUNCSTART("void Dynamic::zeroVars()");

  DYN_maxscore=INT_MIN;

  DYN_sequence1=nullptr;
  DYN_sequence2=nullptr;
  DYN_s1size=0;
  DYN_s2size=0;
  DYN_knowngaps=0;
  DYN_simmatrix=nullptr;
  DYN_smsize=0;

  DYN_validseq=false;
  DYN_matrixcalculated=false;
  DYN_calcwithoffset=false;
  DYN_eoffset=0;

  DYN_bandwidth=0;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Dynamic::matinit(char a, char b, int32 score)
{
  a=static_cast<char>(toupper(a));
  b=static_cast<char>(toupper(b));
  DYN_match_matrix[static_cast<uint8>(a)][static_cast<uint8>(b)]=score;
  DYN_match_matrix[static_cast<uint8>(b)][static_cast<uint8>(a)]=score;
  a=static_cast<char>(tolower(a));
  DYN_match_matrix[static_cast<uint8>(a)][static_cast<uint8>(b)]=score;
  DYN_match_matrix[static_cast<uint8>(b)][static_cast<uint8>(a)]=score;
  b=static_cast<char>(tolower(b));
  DYN_match_matrix[static_cast<uint8>(a)][static_cast<uint8>(b)]=score;
  DYN_match_matrix[static_cast<uint8>(b)][static_cast<uint8>(a)]=score;
  a=static_cast<char>(toupper(a));
  DYN_match_matrix[static_cast<uint8>(a)][static_cast<uint8>(b)]=score;
  DYN_match_matrix[static_cast<uint8>(b)][static_cast<uint8>(a)]=score;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Dynamic::init()
{
  FUNCSTART("Dynamic::init()");

  zeroVars();

  if(DYN_matvalid==0){
    dynamic_parameters const & dp = DYN_miraparams->getDynamicParams();

    matinit('X','X',dp.dyn_score_nmatch);
    matinit('X','N',dp.dyn_score_nmatch);
    matinit('X','A',dp.dyn_score_nmatch);
    matinit('X','C',dp.dyn_score_nmatch);
    matinit('X','G',dp.dyn_score_nmatch);
    matinit('X','T',dp.dyn_score_nmatch);
    matinit('X','M',dp.dyn_score_nmatch);
    matinit('X','R',dp.dyn_score_nmatch);
    matinit('X','W',dp.dyn_score_nmatch);
    matinit('X','S',dp.dyn_score_nmatch);
    matinit('X','Y',dp.dyn_score_nmatch);
    matinit('X','K',dp.dyn_score_nmatch);
    matinit('X','V',dp.dyn_score_nmatch);
    matinit('X','H',dp.dyn_score_nmatch);
    matinit('X','D',dp.dyn_score_nmatch);
    matinit('X','B',dp.dyn_score_nmatch);

#if 0
    // prior to PACBIO, this was standard
    // however, with strobe sequencing the other way round
    //  makes much more sense
    // perhaps also not only with strobes ... hmmm.
    matinit('N','N',dp.dyn_score_nmatch);
    matinit('N','A',dp.dyn_score_npenaltymatch);
    matinit('N','C',dp.dyn_score_npenaltymatch);
    matinit('N','G',dp.dyn_score_npenaltymatch);
    matinit('N','T',dp.dyn_score_npenaltymatch);
    matinit('N','M',dp.dyn_score_npenaltymatch);
    matinit('N','R',dp.dyn_score_npenaltymatch);
    matinit('N','W',dp.dyn_score_npenaltymatch);
    matinit('N','S',dp.dyn_score_npenaltymatch);
    matinit('N','Y',dp.dyn_score_npenaltymatch);
    matinit('N','K',dp.dyn_score_npenaltymatch);
    matinit('N','V',dp.dyn_score_npenaltymatch);
    matinit('N','H',dp.dyn_score_npenaltymatch);
    matinit('N','D',dp.dyn_score_npenaltymatch);
    matinit('N','B',dp.dyn_score_npenaltymatch);
#else
    matinit('N','N',dp.dyn_score_npenaltymatch);
    matinit('N','A',dp.dyn_score_nmatch);
    matinit('N','C',dp.dyn_score_nmatch);
    matinit('N','G',dp.dyn_score_nmatch);
    matinit('N','T',dp.dyn_score_nmatch);
    matinit('N','M',dp.dyn_score_nmatch);
    matinit('N','R',dp.dyn_score_nmatch);
    matinit('N','W',dp.dyn_score_nmatch);
    matinit('N','S',dp.dyn_score_nmatch);
    matinit('N','Y',dp.dyn_score_nmatch);
    matinit('N','K',dp.dyn_score_nmatch);
    matinit('N','V',dp.dyn_score_nmatch);
    matinit('N','H',dp.dyn_score_nmatch);
    matinit('N','D',dp.dyn_score_nmatch);
    matinit('N','B',dp.dyn_score_nmatch);
#endif

    matinit('X','#',dp.dyn_score_oldgap);
    matinit('N','#',dp.dyn_score_oldgap);
    matinit('A','#',dp.dyn_score_oldgap);
    matinit('C','#',dp.dyn_score_oldgap);
    matinit('G','#',dp.dyn_score_oldgap);
    matinit('T','#',dp.dyn_score_oldgap);
    matinit('M','#',dp.dyn_score_oldgap);
    matinit('R','#',dp.dyn_score_oldgap);
    matinit('W','#',dp.dyn_score_oldgap);
    matinit('S','#',dp.dyn_score_oldgap);
    matinit('Y','#',dp.dyn_score_oldgap);
    matinit('K','#',dp.dyn_score_oldgap);
    matinit('V','#',dp.dyn_score_oldgap);
    matinit('H','#',dp.dyn_score_oldgap);
    matinit('D','#',dp.dyn_score_oldgap);
    matinit('B','#',dp.dyn_score_oldgap);

    // 1 == gap, but A below
    matinit('X','1',dp.dyn_score_oldgap);
    matinit('N','1',dp.dyn_score_oldgap);
    matinit('A','1',dp.dyn_score_halfmatch);
    matinit('C','1',dp.dyn_score_oldgap);
    matinit('G','1',dp.dyn_score_oldgap);
    matinit('T','1',dp.dyn_score_oldgap);
    matinit('M','1',dp.dyn_score_oldgap);
    matinit('R','1',dp.dyn_score_oldgap);
    matinit('W','1',dp.dyn_score_oldgap);
    matinit('S','1',dp.dyn_score_oldgap);
    matinit('Y','1',dp.dyn_score_oldgap);
    matinit('K','1',dp.dyn_score_oldgap);
    matinit('V','1',dp.dyn_score_oldgap);
    matinit('H','1',dp.dyn_score_oldgap);
    matinit('D','1',dp.dyn_score_oldgap);
    matinit('B','1',dp.dyn_score_oldgap);

    // 2 == gap, but C below
    matinit('X','2',dp.dyn_score_oldgap);
    matinit('N','2',dp.dyn_score_oldgap);
    matinit('A','2',dp.dyn_score_oldgap);
    matinit('C','2',dp.dyn_score_halfmatch);
    matinit('G','2',dp.dyn_score_oldgap);
    matinit('T','2',dp.dyn_score_oldgap);
    matinit('M','2',dp.dyn_score_oldgap);
    matinit('R','2',dp.dyn_score_oldgap);
    matinit('W','2',dp.dyn_score_oldgap);
    matinit('S','2',dp.dyn_score_oldgap);
    matinit('Y','2',dp.dyn_score_oldgap);
    matinit('K','2',dp.dyn_score_oldgap);
    matinit('V','2',dp.dyn_score_oldgap);
    matinit('H','2',dp.dyn_score_oldgap);
    matinit('D','2',dp.dyn_score_oldgap);
    matinit('B','2',dp.dyn_score_oldgap);

    // 3 == gap, but G below
    matinit('X','3',dp.dyn_score_oldgap);
    matinit('N','3',dp.dyn_score_oldgap);
    matinit('A','3',dp.dyn_score_oldgap);
    matinit('C','3',dp.dyn_score_oldgap);
    matinit('G','3',dp.dyn_score_halfmatch);
    matinit('T','3',dp.dyn_score_oldgap);
    matinit('M','3',dp.dyn_score_oldgap);
    matinit('R','3',dp.dyn_score_oldgap);
    matinit('W','3',dp.dyn_score_oldgap);
    matinit('S','3',dp.dyn_score_oldgap);
    matinit('Y','3',dp.dyn_score_oldgap);
    matinit('K','3',dp.dyn_score_oldgap);
    matinit('V','3',dp.dyn_score_oldgap);
    matinit('H','3',dp.dyn_score_oldgap);
    matinit('D','3',dp.dyn_score_oldgap);
    matinit('B','3',dp.dyn_score_oldgap);

    // 4 == gap, but T below
    matinit('X','4',dp.dyn_score_oldgap);
    matinit('N','4',dp.dyn_score_oldgap);
    matinit('A','4',dp.dyn_score_oldgap);
    matinit('C','4',dp.dyn_score_oldgap);
    matinit('G','4',dp.dyn_score_oldgap);
    matinit('T','4',dp.dyn_score_halfmatch);
    matinit('M','4',dp.dyn_score_oldgap);
    matinit('R','4',dp.dyn_score_oldgap);
    matinit('W','4',dp.dyn_score_oldgap);
    matinit('S','4',dp.dyn_score_oldgap);
    matinit('Y','4',dp.dyn_score_oldgap);
    matinit('K','4',dp.dyn_score_oldgap);
    matinit('V','4',dp.dyn_score_oldgap);
    matinit('H','4',dp.dyn_score_oldgap);
    matinit('D','4',dp.dyn_score_oldgap);
    matinit('B','4',dp.dyn_score_oldgap);

    matinit('1','1',dp.dyn_score_match);
    matinit('1','2',dp.dyn_score_oldgapmatch);
    matinit('1','3',dp.dyn_score_oldgapmatch);
    matinit('1','4',dp.dyn_score_oldgapmatch);

    matinit('2','2',dp.dyn_score_match);
    matinit('2','3',dp.dyn_score_oldgapmatch);
    matinit('2','4',dp.dyn_score_oldgapmatch);

    matinit('3','3',dp.dyn_score_match);
    matinit('3','4',dp.dyn_score_oldgapmatch);

    matinit('4','4',dp.dyn_score_match);


    matinit('A','A',dp.dyn_score_match);
    matinit('A','C',dp.dyn_score_mismatch);
    matinit('A','G',dp.dyn_score_mismatch);
    matinit('A','T',dp.dyn_score_mismatch);
    matinit('A','M',dp.dyn_score_nmatch);
    matinit('A','R',dp.dyn_score_nmatch);
    matinit('A','W',dp.dyn_score_nmatch);
    matinit('A','S',dp.dyn_score_mismatch);
    matinit('A','Y',dp.dyn_score_mismatch);
    matinit('A','K',dp.dyn_score_mismatch);
    matinit('A','V',dp.dyn_score_halfmismatch);
    matinit('A','H',dp.dyn_score_halfmismatch);
    matinit('A','D',dp.dyn_score_halfmismatch);
    matinit('A','B',dp.dyn_score_mismatch);

    matinit('C','C',dp.dyn_score_match);
    matinit('C','G',dp.dyn_score_mismatch);
    matinit('C','T',dp.dyn_score_mismatch);
    matinit('C','M',dp.dyn_score_nmatch);
    matinit('C','R',dp.dyn_score_mismatch);
    matinit('C','W',dp.dyn_score_mismatch);
    matinit('C','S',dp.dyn_score_nmatch);
    matinit('C','Y',dp.dyn_score_nmatch);
    matinit('C','K',dp.dyn_score_mismatch);
    matinit('C','V',dp.dyn_score_halfmismatch);
    matinit('C','H',dp.dyn_score_halfmismatch);
    matinit('C','D',dp.dyn_score_mismatch);
    matinit('C','B',dp.dyn_score_halfmismatch);

    matinit('G','G',dp.dyn_score_match);
    matinit('G','T',dp.dyn_score_mismatch);
    matinit('G','M',dp.dyn_score_mismatch);
    matinit('G','R',dp.dyn_score_nmatch);
    matinit('G','W',dp.dyn_score_mismatch);
    matinit('G','S',dp.dyn_score_nmatch);
    matinit('G','Y',dp.dyn_score_mismatch);
    matinit('G','K',dp.dyn_score_nmatch);
    matinit('G','V',dp.dyn_score_halfmismatch);
    matinit('G','H',dp.dyn_score_mismatch);
    matinit('G','D',dp.dyn_score_halfmismatch);
    matinit('G','B',dp.dyn_score_halfmismatch);

    matinit('T','T',dp.dyn_score_match);
    matinit('T','M',dp.dyn_score_mismatch);
    matinit('T','R',dp.dyn_score_mismatch);
    matinit('T','W',dp.dyn_score_nmatch);
    matinit('T','S',dp.dyn_score_mismatch);
    matinit('T','Y',dp.dyn_score_nmatch);
    matinit('T','K',dp.dyn_score_nmatch);
    matinit('T','V',dp.dyn_score_mismatch);
    matinit('T','H',dp.dyn_score_halfmismatch);
    matinit('T','D',dp.dyn_score_halfmismatch);
    matinit('T','B',dp.dyn_score_halfmismatch);

    matinit('M','M',dp.dyn_score_match);
    matinit('M','R',dp.dyn_score_nmatch);
    matinit('M','W',dp.dyn_score_nmatch);
    matinit('M','S',dp.dyn_score_nmatch);
    matinit('M','Y',dp.dyn_score_nmatch);
    matinit('M','K',dp.dyn_score_mismatch);
    matinit('M','V',dp.dyn_score_halfmatch);
    matinit('M','H',dp.dyn_score_halfmatch);
    matinit('M','D',dp.dyn_score_halfmismatch);
    matinit('M','B',dp.dyn_score_halfmismatch);

    matinit('R','R',dp.dyn_score_match);
    matinit('R','W',dp.dyn_score_nmatch);
    matinit('R','S',dp.dyn_score_nmatch);
    matinit('R','Y',dp.dyn_score_mismatch);
    matinit('R','K',dp.dyn_score_nmatch);
    matinit('R','V',dp.dyn_score_halfmatch);
    matinit('R','H',dp.dyn_score_halfmismatch);
    matinit('R','D',dp.dyn_score_halfmatch);
    matinit('R','B',dp.dyn_score_halfmismatch);

    matinit('W','W',dp.dyn_score_match);
    matinit('W','S',dp.dyn_score_mismatch);
    matinit('W','Y',dp.dyn_score_nmatch);
    matinit('W','K',dp.dyn_score_nmatch);
    matinit('W','V',dp.dyn_score_halfmismatch);
    matinit('W','H',dp.dyn_score_halfmatch);
    matinit('W','D',dp.dyn_score_halfmatch);
    matinit('W','B',dp.dyn_score_halfmismatch);

    matinit('S','S',dp.dyn_score_match);
    matinit('S','Y',dp.dyn_score_nmatch);
    matinit('S','K',dp.dyn_score_nmatch);
    matinit('S','V',dp.dyn_score_halfmatch);
    matinit('S','H',dp.dyn_score_halfmismatch);
    matinit('S','D',dp.dyn_score_halfmismatch);
    matinit('S','B',dp.dyn_score_halfmatch);

    matinit('Y','Y',dp.dyn_score_match);
    matinit('Y','K',dp.dyn_score_nmatch);
    matinit('Y','V',dp.dyn_score_halfmismatch);
    matinit('Y','H',dp.dyn_score_halfmatch);
    matinit('Y','D',dp.dyn_score_halfmismatch);
    matinit('Y','B',dp.dyn_score_halfmatch);

    matinit('K','K',dp.dyn_score_match);
    matinit('K','V',dp.dyn_score_halfmismatch);
    matinit('K','H',dp.dyn_score_halfmismatch);
    matinit('K','D',dp.dyn_score_halfmatch);
    matinit('K','B',dp.dyn_score_halfmatch);

    matinit('V','V',dp.dyn_score_match);
    matinit('V','H',dp.dyn_score_halfmatch);
    matinit('V','D',dp.dyn_score_halfmatch);
    matinit('V','B',dp.dyn_score_halfmatch);

    matinit('H','H',dp.dyn_score_match);
    matinit('H','D',dp.dyn_score_halfmatch);
    matinit('H','B',dp.dyn_score_halfmatch);

    matinit('D','D',dp.dyn_score_match);
    matinit('D','B',dp.dyn_score_halfmatch);

    matinit('B','B',dp.dyn_score_match);


    DYN_matvalid=1;
  }

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
// Destructor to clean up the memory
// It just calls discard(), as everything that might need to be cleaned
//  up is done there
Dynamic::~Dynamic()
{
  FUNCSTART("Dynamic::~Dynamic()");

  discard();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// setSequences lets this instance now which sequences it should
//  make an dynamic programming matrix of.
// The sequences need not to be terminated with 0. Instead, the length
//  must be given as parameter (this allows live 'cutouts' from longer
//  longer sequences.
// Once acquired, the dynamic matrix is immediately computed.
// TODO:
//  It's bad to compute it immediately. If one wants to change parameters,
//  it must be recomputed.
void Dynamic::setSequences(const char * seq1, uint32 len1, const char * seq2, uint32 len2, bool calcwithoffset, int32 expectedoffset)
{
  FUNCSTART("Dynamic::setSequences(const char * seq1, uint32 len1, const char * seq2, uint32 len2, bool calcwithoffset, int32 expectedoffset)");

  BUGIFTHROW(len1==0 || len2==0,"len1 (" << len1 << ") or len2 (" << len2 << ") == 0 ?");
  BUGIFTHROW(len1>200000 || len2>200000,"len1 (" << len1 << ") or len2 (" << len2 << ") >200000 ? Seems unbelievable.");

  DYN_calcwithoffset=calcwithoffset;
  DYN_eoffset=expectedoffset;

  DYN_maxscore=INT_MIN;
  DYN_lastrc_maxscore=INT_MIN;
  DYN_knowngaps=0;

  if(DYN_sequence1 == nullptr || DYN_s1size<len1+1){
    if(DYN_sequence1 != nullptr) delete [] DYN_sequence1;
    DYN_s1size=len1+1;
    if(DYN_s1size<2000) DYN_s1size=2000;
    DYN_sequence1= new char[DYN_s1size];
    ++DYN_alloccounts;
  }
  if(DYN_sequence2 == nullptr || DYN_s2size<=len2+1){
    if(DYN_sequence2 != nullptr) delete [] DYN_sequence2;
    DYN_s2size=len2+1;
    if(DYN_s2size<2000) DYN_s2size=2000;
    DYN_sequence2= new char[DYN_s2size];
    ++DYN_alloccounts;
  }

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  // Copy first sequence
  try{
    //DYN_sequence1= new char[len1+1];
    DYN_knowngaps=sequenceCopy(DYN_sequence1, seq1, len1);
    DYN_len_seq1=len1;
  }
  catch(Notify n){
    cout << "Error while copying sequence 1" << endl;
    //    n.handleError(THISFUNC);
    throw(n);
  }

  // Copy second sequence
  try{
    //DYN_sequence2= new char[len2+1];
    DYN_knowngaps-=sequenceCopy(DYN_sequence2, seq2, len2);
    DYN_len_seq2=len2;
  }
  catch(Notify n){
    cout << "Error while copying sequence 2" << endl;
    //    n.handleError(THISFUNC);
    throw(n);
  }

#ifdef CLOCK_STEPS1
  DYN_timing_seqcopy+=diffsuseconds(tv);
#endif


  if(DYN_knowngaps<0) DYN_knowngaps=-DYN_knowngaps;

  // len+1 to get the 0 row and 0 column into the matrix
  //  +1 at the end ta allow the matrix computation algorithm to be easier
#if 1

  {
    uint32 sizeneeded=(len1+1)*(len2+1) +1;

    //cout << "sizeneeded1: " << sizeneeded << endl;

    if (DYN_simmatrix==nullptr) {
      // ok, on first use, we're taking at least 1024^2+1 bytes
      sizeneeded=std::max(sizeneeded,static_cast<uint32>(1024*1024+1));

      //cout << "DYN_simmatrix: " << DYN_simmatrix << "\tsizeneeded2: " << sizeneeded << endl;

      DYN_simmatrix= new int32[sizeneeded];
      DYN_smsize=sizeneeded;
    } else {
      if(DYN_smsize < sizeneeded) {
	delete [] DYN_simmatrix;
	DYN_simmatrix= new int32[sizeneeded];
	DYN_smsize=sizeneeded;
	++DYN_alloccountm;
      }
    }
  }

#else
#define SPECIALCEBUG
  if (DYN_simmatrix!=nullptr) delete [] DYN_smsize;

  DYN_simmatrix= new int32[(len1+1)*(len2+1) +1 +2000];
  // Nur wenn +2000 bei Speicheralloziierung
  {
    int32 * ptr= DYN_simmatrix+(len1+1)*(len2+1) +1;
    for(uint32 i=0; i<2000; i++, ptr++){
      *ptr=0xDEADBEEF;
    }
  }
#endif

  DYN_validseq=true;
  DYN_matrixcalculated=false;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Discard
// Frees the used memory
void Dynamic::discard()
{
  FUNCSTART("Dynamic::discard()");

  if(DYN_sequence1!=nullptr) delete [] DYN_sequence1;
  if(DYN_sequence2!=nullptr) delete [] DYN_sequence2;
  if(DYN_simmatrix!=nullptr) delete [] DYN_simmatrix;

  FUNCEND();
}




/*************************************************************************
 *
 * Copy the given sequence (modifying it) and give back number of gaps
 *
 * This makes this instance independant from the source sequence as we
 *  modify it!
 *
 *  translate all into uppercase bases
 *  translate undefined bases to N,
 *  translate gaps to "#"
 *  add 0 terminator to the sequence
 *
 *************************************************************************/

int32 Dynamic::sequenceCopy(char * to, const char * from, uint32 len)
{
  FUNCSTART("void Dynamic::sequenceCopy(char * to, const char * from, uint32 len)");

  int32 numgaps=0;
  auto ofrom=from;
  for(;from<ofrom+len; ++from, ++to){
    char base=static_cast<char>(toupper(*from));
    if(likely(dptools::isValidIUPACBase(base))){
      *to=base;
    }else{
      switch(base){
      case '*':{
	*to='#';
	++numgaps;
	break;
      }
      case '1':
      case '2':
      case '3':
      case '4':{
	*to=base;
	++numgaps;
	break;
      }
      case '-':
	*to='N';
	break;
      default:{
	cout << "Position: " << from-ofrom << "\t" << std::hex << static_cast<uint16>(*from) << "\t" << *from << endl;
	MIRANOTIFY(Notify::FATAL, "Unknown base in read: " << from);
      }
      }
    }
  }
  *to=0;

  FUNCEND();

  return numgaps;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  // The sequences are known, the memory for the similarity matrix
  //  allocated, now compute this matrix

void Dynamic::computeMatrix()
{
  FUNCSTART("void Dynamic::computeMatrix()");

  //if an expected offset of sequence 1 and 2 has been given, use this and
  // compute a banded matrix, else whole matrix
  if(DYN_calcwithoffset){
    //    computeSimMatrix();
    computeBSimMatrix();
  }else{
    BUGIFTHROW(true,"computeSimMatrix() not available atm, sorry. use calcwithoffset and large band.");
    //computeSimMatrix();
  }

  FUNCEND();
}
/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

/* expected offset of seq2 to seq1
       s1 -------                                   s1 -------
         s2 --------     eoffset =2	          s2 --------     eoffset =-2

            s2				                s2
          0------			             0-^----
         0				            0  xpos=2
	 |				       	    |
         ypos=2				            |
	 |				       	    |
      s1 |				         s1 |
	 |				       	    |
	 |				       	    |
	 |				       	    |
          -------                                    -------

*/

void Dynamic::computeBSimMatrix()
{
  FUNCSTART("Dynamic::computeBSimMatrix()");

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  // xpos and ypos are the offsets computed from the expected overlap
  // xrun and yrun point to the cell in the matrix to be computed next
  //  and to the left and above of the cell in the simmatrix
  // eoverlap is the expected overlap
  int32 xpos=0;
  int32 ypos=0;
  int32 xrun=0;
  int32 yrun=0;

  int32 eoverlap;
  if(DYN_eoffset>=0){
    eoverlap=static_cast<int32>(DYN_len_seq1)-DYN_eoffset;
    if(static_cast<int32>(DYN_len_seq2)<eoverlap) eoverlap=DYN_len_seq2;
    ypos=DYN_eoffset;
  }else{
    eoverlap=static_cast<int32>(DYN_len_seq2)+DYN_eoffset;
    if(static_cast<int32>(DYN_len_seq1)<eoverlap) eoverlap=DYN_len_seq1;
    xpos=-DYN_eoffset;
  }


  // I allow % inserts/deletes, so make the
  //  band % bases of the expected overlap
  // If band > length of one of the sequences reduce it
  //
  // kerrors is the number of cells left and right from the
  //  idel diagonal that are computed

  int32 kerrors=DYN_knowngaps+static_cast<int32>(eoverlap*static_cast<float>(DYN_miraparams->getAlignParams().al_kpercent)/100.0f);

  if(kerrors<DYN_miraparams->getAlignParams().al_kmin) kerrors=DYN_miraparams->getAlignParams().al_kmin;
  if(kerrors>DYN_miraparams->getAlignParams().al_kmax) {
    kerrors=DYN_miraparams->getAlignParams().al_kmax;
  }

  if(kerrors>static_cast<int32>(std::max(DYN_len_seq1,DYN_len_seq2))) kerrors=static_cast<int32>(std::max(DYN_len_seq1,DYN_len_seq2));
  if(kerrors<5) kerrors=5;

  if(kerrors-DYN_knowngaps < 0) kerrors=DYN_miraparams->getAlignParams().al_kmin+DYN_knowngaps;

  const int32 bandwidth=2*kerrors+1;

  DYN_bandwidth=bandwidth;

  //cout << "KGAP: " << DYN_knowngaps << endl;
  //cout << "KERR: " << kerrors << endl;

  if(DYN_eoffset>=0){
    DYN_leftbandx=-ypos-kerrors;
    DYN_rightbandx=-ypos+kerrors;
  }else{
    DYN_leftbandx=xpos-kerrors;
    DYN_rightbandx=xpos+kerrors;
  }
  DYN_bandsafety=eoverlap/100;
  if(DYN_bandsafety<1){
    DYN_bandsafety=1;
  }else if(DYN_bandsafety>15){
    DYN_bandsafety=15;
  }

  // xpos and ypos are the offsets computed from the expected overlap
  if(ypos>0){
    yrun=ypos-kerrors;
    if(yrun<0) yrun=0;
  }
  if(xpos>0){
    xrun=xpos-kerrors;
    if(xrun<0) xrun=0;
  }

  CEBUG("lens1: " << DYN_len_seq1 << endl);
  CEBUG("lens2: " << DYN_len_seq2 << endl);
  CEBUG("eoffset: " << DYN_eoffset << endl);
  CEBUG("eoverlap: " << eoverlap << endl);
  CEBUG("kerrors: " << kerrors << endl);

  dynamic_parameters const & DYN_params = DYN_miraparams->getDynamicParams();

#ifdef MATRIXDEBUG
  {
    // Fill the matrix with 666 for testing and debugging purposes
    int32 * ptr=DYN_simmatrix;
    for(uint32 i=0; i<(DYN_len_seq1+1)*(DYN_len_seq2+1); ++i, ++ptr) *ptr=666;
  }
#endif

#ifdef CLOCK_STEPS1
  timeval tl;
  gettimeofday(&tl,nullptr);
#endif

  // fill last row with -1
  {
    int32 * ptrr=DYN_simmatrix+(DYN_len_seq1)*(DYN_len_seq2+1);
    for(uint32 i=0; i<=DYN_len_seq2; ++i){
      *ptrr++= -1;
    }
    // no measurable effect on PacBio 6k overlap
    //memset(ptrr,-1,DYN_len_seq2*sizeof(int32));
  }
  //cout << "###1 " << diffsuseconds(tl) << endl;

  // fill first row terminal gap penalties.
  {
    int32 * ptrr=DYN_simmatrix;
    for(uint32 i=0; i<=DYN_len_seq2; ++i){
      *ptrr++= DYN_params.dyn_score_ltermgap*i;
    }
  }
  //cout << "###2 " << diffsuseconds(tl) << endl;

  // filling first and last column almost certainly will lead to cache misses
  // on large sequences (pacbio, with 7 or 9 Kb), this is responsible for 1/8
  //  to 1/4 of the total time spent in banded Smith-Waterman, depending on band-width
  //
  // therefore, perform a preloading beforehand so that the values will be
  //  at least in usual L3 caches, probably even in L2
  // remember, last column + 1 == first column (next row)

  // Apple MacAir corei7: this is *SO* absolutely, totally, utterly frustrating
  //  the prefetch below has exactly 0 effect.
  // putting a prefetch in the main loop also 0 effect, regardless on
  //  how far I let it prefetch
  // TODO: check on home Linux box

//  {
//    int32 * ptrlc=DYN_simmatrix+DYN_len_seq2;
//    for(uint32 i=0; i<=DYN_len_seq1; ++i, ptrlc+=DYN_len_seq2+1){
//      __builtin_prefetch (ptrlc);
//      //prefetchwl(ptrlc);
//    }
//  }
  //cout << "### " << diffsuseconds(tl) << endl;

  // fill first column with terminal gap penalties, last column with -1
  {
    int32 * ptrfc=DYN_simmatrix;
    int32 * ptrlc=DYN_simmatrix+DYN_len_seq2;
    for(uint32 i=0; i<=DYN_len_seq1; ++i, ptrfc+=DYN_len_seq2+1, ptrlc+=DYN_len_seq2+1){
      *ptrfc= DYN_params.dyn_score_rtermgap*i;
      *ptrlc= -1;
      //prefetchwl(ptrlc+100*(DYN_len_seq2+1));
    }
  }
  //cout << "###4 " << diffsuseconds(tl) << endl;

#ifdef CLOCK_STEPS1
  DYN_timing_bswm_setup+=diffsuseconds(tl);
#endif

  if((xpos==0 && ypos-kerrors>static_cast<int32>(DYN_len_seq1))
     || (ypos==0 && xpos-kerrors>static_cast<int32>(DYN_len_seq2))){
    DYN_lastrc_maxscore=0;
    DYN_maxscore=0;
    FUNCEND();
    return;
  }else{


#ifdef MATRIXDEBUG
    int32 bandlimit=DYN_BANDLIMIT;
#endif

    // Local variables are faster than class variables
    const int32 s_sgap=DYN_params.dyn_score_gap;


#ifdef CLOCK_STEPS1
    gettimeofday(&tl,nullptr);
#endif
    {
      // Neu Teil 1

      int32 dolen=0;
      int32 doheight=0;

      if(ypos>0){
	CEBUG("Runter");
	if(ypos-kerrors>=0){
	  // spitzes Dreieck
	  //
	  //  .....
	  //  X....
	  //  XX...
	  // etc
	  CEBUG("Spitz\n");
	  xrun=xpos;
	  yrun=ypos-kerrors;
	  dolen=0;
	  doheight=2*kerrors;
	  if(yrun+doheight>static_cast<int32>(DYN_len_seq1)) doheight=DYN_len_seq1-yrun;
	  // xpos=0!
	  if(xrun+doheight>static_cast<int32>(DYN_len_seq2)) doheight=DYN_len_seq2-xrun;
	}else{
	  // oben flaches Dreieck oder recht
	  xrun=0;
	  yrun=0;
	  // ypos-kerrors < 0 !
	  dolen=-(ypos-kerrors);
	  if(dolen>=static_cast<int32>(DYN_len_seq2)){
	    CEBUG("Rechteck\n");
	    // kein Dreieck, Rechteck. volle Zeile muss berechnet werden, später
	    doheight=0;
	  }else{
	    CEBUG("Flach\n");
	    doheight=2*kerrors+(ypos-kerrors);

	    if(yrun+doheight>static_cast<int32>(DYN_len_seq1)) doheight=DYN_len_seq1-yrun;
	    if(xrun+dolen+doheight>static_cast<int32>(DYN_len_seq2)) doheight=DYN_len_seq2-xrun-dolen;
	  }
	}
      }else{
	CEBUG("Rechts");
	if(xpos-kerrors<0){
	  dolen=xpos+kerrors;
	  CEBUG("\ndolentmp: " << dolen<< endl);
	  if(dolen>=static_cast<int32>(DYN_len_seq2)){
	    CEBUG("Rechteck\n");
	    // kein Dreieck, Rechteck. volle Zeile muss berechnet werden, später
	    doheight=0;
	  } else {
	    doheight=-(xpos-kerrors);
	    CEBUG("doheighttmp: " << doheight<< endl);
	    if(yrun+doheight>static_cast<int32>(DYN_len_seq1)) {
	      CEBUG("p1\n");
	      doheight=DYN_len_seq1-yrun;
	    }
	    if(xrun+dolen+doheight>static_cast<int32>(DYN_len_seq2)) {
	      CEBUG("p2\n");
	      doheight=DYN_len_seq2-xrun-dolen;
	    }
	  }
	}
      }

      CEBUG("\nTeil 1\n");
      CEBUG("xpos: " << xpos<< endl);
      CEBUG("ypos: " << ypos<< endl);
      CEBUG("xrun: " << xrun<< endl);
      CEBUG("yrun: " << yrun<< endl);
      CEBUG("dolen: " << dolen<< endl);
      CEBUG("doheight: " << doheight<< endl);

#ifdef MATRIXDEBUG
      if(yrun>0) *(DYN_simmatrix+yrun*(DYN_len_seq2+1)+1)=bandlimit;
#endif

      int32 * ptrt=DYN_simmatrix+(yrun+1)*(DYN_len_seq2+1)+1;
      for(int32 zeile=0; zeile<doheight; ++zeile, ++yrun, ++dolen){
	const int32 * ptrla=DYN_simmatrix+yrun*(DYN_len_seq2+1);
	const int32 * ptra=ptrla+1;
	const int32 * ptrl=ptrla+DYN_len_seq2+1;
	ptrt=const_cast<int32 *>(ptrl+1);

	const char * s_ds2s=DYN_sequence2;
	CEBUG("1 " << DYN_sequence1[yrun] << ":");

	const int32 * mmp= (int32 *) &DYN_match_matrix+DYN_sequence1[yrun]*DYN_MATSIZE;
	for(int32 spalte=0; spalte<dolen; ++spalte,
	      ++s_ds2s,++ptrla,++ptra,++ptrl,++ptrt){
	  //CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")");
	  CEBUG("1");
	  *ptrt= std::max((*ptra)+s_sgap,
		     std::max((*ptrla)+mmp[*s_ds2s], (*ptrl)+s_sgap ));
	  //prefetchrl(ptra+1);
	  //prefetchrl(ptrla+1);
	  //prefetchrl(ptrl+1);
	  //prefetchwl(ptrt+1);
	}
	//CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")" << endl);
	CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")" << endl);
	*ptrt= std::max((*ptrla)+mmp[*s_ds2s], (*ptrl)+s_sgap );

#ifdef MATRIXDEBUG
	*(ptrt+1)=bandlimit;
#endif
      }

      //if(doheight >0 && xrun+dolen>=static_cast<int32>(DYN_len_seq2)) {
      //	BUGIFTHROW(true,"Is this dead code? Hopefully, but it was triggered anyway :-( It was");
      //	*(ptrt+1)=-212;
      //}
    }

#ifdef CLOCK_STEPS1
    DYN_timing_bswm_p1+=diffsuseconds(tl);
#endif

#ifdef MATRIXDEBUG
    cout << "Mathead\n";
    dump();
#endif

#ifdef CLOCK_STEPS1
    gettimeofday(&tl,nullptr);
#endif

    {
/*
  Teil 2a: mittelteil voll?
  takes care of eventually full square covered by the band.
  e.g. (marked with X)

    ....
    \...
     \..
      \.
    XXXX
    XXXX
    XXXX
    \...
     \..

 */
      CEBUG("\nTeil 2a\n");
      int32 runtoline=-DYN_leftbandx;
      if(runtoline>static_cast<int32>(DYN_len_seq1)) runtoline=DYN_len_seq1;
      if(bandwidth<static_cast<int32>(DYN_len_seq2)) runtoline=-1;

      CEBUG("xpos: " << xpos<< endl);
      CEBUG("ypos: " << ypos<< endl);
      CEBUG("xrun: " << xrun<< endl);
      CEBUG("yrun: " << yrun<< endl);
      CEBUG("runtoline: " << runtoline<< endl);

      const int32 * ptrla=DYN_simmatrix+yrun*(DYN_len_seq2+1)+xrun;
      const int32 * ptra=ptrla+1;
      const int32 * ptrl=ptrla+DYN_len_seq2+1;
      int32 * ptrt=const_cast<int32 *>(ptrl+1);
      for(; yrun<runtoline; ++yrun,
	    ++ptrla, ++ptra, ++ptrl, ++ptrt){
	const char * s_ds2s=DYN_sequence2;
	const int32 * mmp= (int32 *) &DYN_match_matrix+DYN_sequence1[yrun]*DYN_MATSIZE;

	CEBUG("2a " << DYN_sequence1[yrun] << ":");
	for(uint32 spalte=0; spalte<DYN_len_seq2; ++spalte,
	      ++s_ds2s,++ptrla,++ptra,++ptrl,++ptrt){
	  //CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")");
	  CEBUG("a");
	  *ptrt= std::max((*ptra)+s_sgap,
		     std::max((*ptrla)+mmp[*s_ds2s], (*ptrl)+s_sgap ));
	  //prefetchrl(ptra+1);
	  //prefetchrl(ptrla+1);
	  //prefetchrl(ptrl+1);
	  //prefetchwl(ptrt+1);
	}
	CEBUG(endl);
      }
    }

#ifdef CLOCK_STEPS1
  DYN_timing_bswm_p2a+=diffsuseconds(tl);
#endif

#ifdef MATRIXDEBUG
    cout << "Matmid 2a\n";
    dump();
#endif

#ifdef CLOCK_STEPS1
    gettimeofday(&tl,nullptr);
#endif
    {
      // Teil 2b
/*
  the running middle

     \.......
  XXXXX......
  .XXXXX.....
  ..XXXXX....
  ...XXXXX...
  ....XXXXX..
  .....XXXXX.
  ......XXXXX
  .......\
 */
      CEBUG("\nTeil 2b\n");

      int32 doheight=DYN_len_seq1-yrun;
      CEBUG("doheight-pre: " << doheight<< endl);
      if(xrun+doheight-1+bandwidth>=static_cast<int32>(DYN_len_seq2)){
	doheight= DYN_len_seq2-xrun-bandwidth+1;
      }
      if(yrun>static_cast<int32>(DYN_len_seq1)) doheight=0;
      CEBUG("xpos: " << xpos<< endl);
      CEBUG("ypos: " << ypos<< endl);
      CEBUG("xrun: " << xrun<< endl);
      CEBUG("yrun: " << yrun<< endl);
      CEBUG("doheight: " << doheight<< endl);

      const int32 * ptrla=DYN_simmatrix+yrun*(DYN_len_seq2+1)+xrun;
      const int32 * ptra=ptrla+1;
      const int32 * ptrl=ptrla+DYN_len_seq2+1;
      int32 * ptrt=const_cast<int32 *>(ptrl+1);
      for(int32 zeile=0; zeile<doheight; ++zeile, ++xrun, ++yrun){

	// no effect, maybe even slightly detrimental
	//{
	//  auto ptrp=ptrl+2*(DYN_len_seq2+1);
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//  prefetchrl(ptrp);
	//  ptrp+=16;
	//}
	const char * s_ds2s=DYN_sequence2+xrun;
	const int32 * mmp= (int32 *) &DYN_match_matrix+DYN_sequence1[yrun]*DYN_MATSIZE;

	CEBUG("2b " << DYN_sequence1[yrun] << ":");
	CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")");

#ifdef MATRIXDEBUG
	*(ptrt-1)=bandlimit;
#endif
	*ptrt= std::max((*ptrla)+mmp[*s_ds2s], (*ptra)+s_sgap);
	++s_ds2s,++ptrla,++ptra,++ptrl,++ptrt;

//*
	// BaCh 08.08.2013
	// manual loop unrolling has no measurable effect
	for(int32 spalte=0; spalte<bandwidth-2; ++spalte,
	      ++s_ds2s,++ptrla,++ptra,++ptrl,++ptrt){
	  //CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")");
	  CEBUG("b");
	  *ptrt= std::max((*ptra)+s_sgap,
		     std::max((*ptrla)+mmp[*s_ds2s], (*ptrl)+s_sgap ));
	  // BaCh 08.08.2013
	  // bah! Having these prefetches here makes BSW 3.5 % slower
	  //  on 6000bp overlap of PacBio raw data
	  //prefetchrl(ptrla+1);
	  //prefetchrl(ptra+1);
	  //prefetchrl(ptrl+1);
	  //prefetchwl(ptrt+1);
	  //
	  // also, prefetching one line in advance has zero effect, two lines is even detrimental
	  //prefetchrl(ptrla+2*(DYN_len_seq2+1));
	}
//*/

	CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")" << endl);
	*ptrt= std::max((*ptrla)+mmp[*s_ds2s], (*ptrl)+s_sgap );

#ifdef MATRIXDEBUG
	*(ptrt+1)=bandlimit;
#endif

	ptrla+=(DYN_len_seq2+1)-2*kerrors+1;
	// no effect? maybe a very, very slight effect
	//prefetchrl(ptrla);
	ptra=ptrla+1;
	ptrl=ptrla+DYN_len_seq2+1;
	ptrt=const_cast<int32 *>(ptrl+1);
	// no effect?
	//prefetchwl(ptrt-1);
      }
    }

#ifdef CLOCK_STEPS1
  DYN_timing_bswm_p2b+=diffsuseconds(tl);
#endif

#ifdef MATRIXDEBUG
    cout << "Matmid 2b\n";
    dump();
#endif

#ifdef CLOCK_STEPS1
    gettimeofday(&tl,nullptr);
#endif
    {
      // Teil 3
      int32 dolen=DYN_len_seq2-xrun-1;
      int32 doheight=dolen+1;
      // TODO: check     12.08.1999 >=DYN_len_seq1+1 changed to > DYN_len_seq1
      if(yrun+doheight>static_cast<int32>(DYN_len_seq1)) doheight=DYN_len_seq1-yrun;
      if(yrun>=static_cast<int32>(DYN_len_seq1)) doheight=0;

      CEBUG("\nTeil 3\n");
      CEBUG("dolen: " << dolen<< endl);
      CEBUG("doheight: " << doheight<< endl);
      CEBUG("xrun: " << xrun<< endl);
      CEBUG("yrun: " << yrun<< endl);

      for(int32 zeile=0; zeile<doheight; ++zeile, ++xrun, ++yrun, --dolen){
	const int32 * ptrla=DYN_simmatrix+yrun*(DYN_len_seq2+1)+xrun;
	const int32 * ptra=ptrla+1;
	const int32 * ptrl=ptrla+DYN_len_seq2+1;
	int32 * ptrt=const_cast<int32 *>(ptrl+1);

	const char * s_ds2s=DYN_sequence2+xrun;
	const int32 * mmp= (int32 *) &DYN_match_matrix+DYN_sequence1[yrun]*DYN_MATSIZE;

	CEBUG("3 " << DYN_sequence1[kerrors+zeile] << ":");
	CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")");

#ifdef MATRIXDEBUG
	*(ptrt-1)=bandlimit;
#endif
	*ptrt= std::max((*ptrla)+mmp[*s_ds2s], (*ptra)+s_sgap);
	++s_ds2s,++ptrla,++ptra,++ptrl,++ptrt;

	for(int32 spalte=0; spalte<dolen; ++spalte,
	      ++s_ds2s,++ptrla,++ptra,++ptrl,++ptrt){
	  //CEBUG(*s_ds2s << "(" << mmp[*s_ds2s] << ")");
	  CEBUG("3");
	  *ptrt= std::max((*ptra)+s_sgap,
		     std::max((*ptrla)+mmp[*s_ds2s], (*ptrl)+s_sgap ));
	  //prefetchrl(ptra+1);
	  //prefetchrl(ptrla+1);
	  //prefetchrl(ptrl+1);
	  //prefetchwl(ptrt+1);
	}
	CEBUG(endl);
      }
    }

#ifdef CLOCK_STEPS1
  DYN_timing_bswm_p3+=diffsuseconds(tl);
#endif


#ifdef MATRIXDEBUG
    cout << "Mattail\n";
    dump();
#endif

#ifdef CLOCK_STEPS1
    gettimeofday(&tl,nullptr);
#endif

    //  search for maximum in last column and row

    int32 tmpmax=0;
    {
      int32 * ptrr=DYN_simmatrix+(DYN_len_seq1)*(DYN_len_seq2+1);
      for(uint32 i=0; i<=DYN_len_seq2; ++i, ++ptrr){
	tmpmax= std::max(tmpmax, *ptrr);
      }
    }
    {
      int32 * ptrc=DYN_simmatrix+DYN_len_seq2;
      for(uint32 i=0; i<=DYN_len_seq1; ++i, ptrc+=DYN_len_seq2+1){
	tmpmax= std::max(tmpmax, *ptrc);
      }
    }

#ifdef CLOCK_STEPS1
    DYN_timing_bswm_cleanband+=diffsuseconds(tl);
#endif

    DYN_lastrc_maxscore=tmpmax;
    DYN_maxscore=tmpmax;

    //    cout << DYN_maxscore << endl << DYN_lastrc_maxscore << endl;

  }


#ifdef SPECIALCEBUG
  // Special debug: activate only when +2000 in memory alloc has been done
  {
    int32 * ptr= DYN_simmatrix+(DYN_len_seq1+1)*(DYN_len_seq2+1) +1;
    bool bangit=false;
    for(uint32 i=0; i<2000; i++, ptr++){
	  if(*ptr!=0xDEADBEEF){
	    cout << i << ": " << hex << *ptr<< dec << endl;
	    bangit=true;
	  }
    }
    if(bangit){
	  cout << "Seq1:\n";
	  for(int32 i=0; i< DYN_len_seq1; i++){
	    cout << DYN_sequence1[i];
	  }
	  cout << endl;
	  cout << "Seq2:\n";
	  for(int32 i=0; i< DYN_len_seq2; i++){
	    cout << DYN_sequence2[i];
	  }
	  cout << endl;
	  cout << "Expected offset: " << DYN_eoffset << endl;
	  MIRANOTIFY(Notify::INTERNAL, "Bang, crossed border.") ;
    }
  }
#undef SPECIALCEBUG
#endif

  // Special debugon();
//  if(DYN_len_seq1==560
//     || DYN_len_seq2==560){
//    {
//	cout << "Matlastrow:\n";
//	int32 * ptrr=DYN_simmatrix+(DYN_len_seq1)*(DYN_len_seq2+1);
//	for(int32 i=0; i<=DYN_len_seq2;i++, ptrr++){
//	  cout << *ptrr << "\t";
//	}
//	cout << "\n";
//    }
//    {
//	cout << "Matlastcol:\n";
//	int32 * ptrc=DYN_simmatrix+DYN_len_seq2;
//	for(int32 i=0; i<=DYN_len_seq1;i++, ptrc+=DYN_len_seq2+1){
//	  cout << *ptrc << "\t";
//	}
//	cout << endl;
//    }
//  }

#ifdef CLOCK_STEPS1
  DYN_timing_bswmatrix+=diffsuseconds(tv);
#endif

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Dump the dynamic programming matrix to stdout
void Dynamic::dump()
{
  if(!DYN_matrixcalculated) return;
  const int32 * ptr=DYN_simmatrix;
  for(uint32 i=0; i<DYN_len_seq1+1; ++i){
    for(uint32 j=0; j<DYN_len_seq2+1; ++j){
      cout << *ptr++ << "\t";
    }
    cout << endl;
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// Compute the dynamic programming matrix.
//
// Warning: don't use that one anymore!

/*
void Dynamic::computeSimMatrix()
{
  FUNCSTART("Dynamic::computeSimMatrix()");

  dynamic_parameters const & DYN_params = DYN_miraparams->getDynamicParams();

  // fill first row and first column with terminal gap penalties.
  {
    int32 * ptrr=DYN_simmatrix;
    for(int32 i=0; i<DYN_len_seq2+1;i++){
	*ptrr++= DYN_params.dyn_score_ltermgap*i;
    }
  }

  {
    int32 * ptrc=DYN_simmatrix;
    for(int32 i=0; i<DYN_len_seq1+1;i++){
	*ptrc= DYN_params.dyn_score_rtermgap*i;
	ptrc+=DYN_len_seq2+1;
    }
  }


  {
    // Local variables are faster than class variables
    int32   s_sgap=DYN_params.dyn_score_gap;
    int32   s_smatch=DYN_params.dyn_score_match;
    int32   s_smismatch=DYN_params.dyn_score_mismatch;
    int32   s_snmatch=DYN_params.dyn_score_nmatch;
    int32   s_maxscore=INT_MIN;

    if(DYN_params.dyn_score_termgap==0) s_maxscore=0;

    uint32  s_dls2=DYN_len_seq2;

    // Somewhere between 20%-30% faster that the old version
    // Using pointers is much faster than arrays
    int32 * ptrla=DYN_simmatrix;
    int32 * ptra=ptrla+1;
    int32 * ptrl=DYN_simmatrix+DYN_len_seq2+1;
    int32 * ptrt=ptrl+1;
    for(int32 i=0; i<DYN_len_seq1; i++){
	char   s_ds1i=DYN_sequence1[i];
	char * s_ds2j=DYN_sequence2;
	for(uint32 j=0; j<s_dls2; j++
	    , s_ds2j++, ptrla++,ptra++,ptrl++,ptrt++){
	  int32 vgl;
	  if((s_ds1i== 'N') || (*s_ds2j=='N') || (s_ds1i== 'X') || (*s_ds2j=='X')){
	    vgl=s_snmatch;
	  }else{
	    if(s_ds1i== *s_ds2j){
	      vgl=s_smatch;
	    }else{
	      vgl=s_smismatch;
	    }
	  }
	  *ptrt= std::max((*ptra+s_sgap), std::max((*ptrla+vgl), (*ptrl+s_sgap) ));
	  if(*ptrt>s_maxscore) s_maxscore=*ptrt;

	}
	// Move the pointers one element further to jump over the
	// additional first column which has already been pre-'calculated'
	ptrla++;ptra++;ptrl++;ptrt++;
    }

    DYN_maxscore=s_maxscore;

    // Now search for the maximum in the last row and last column
    int32 tmpmax=0;
    for(int32 i=0;i<=DYN_len_seq1;i++){
	tmpmax= std::max(tmpmax, DYN_simmatrix[i*(DYN_len_seq2+1)+DYN_len_seq2]);
    }
    for(int32 i=0;i<DYN_len_seq2;i++){
	tmpmax= std::max(tmpmax, DYN_simmatrix[((DYN_len_seq1)*(DYN_len_seq2+1)+i)]);
    }
    DYN_lastrc_maxscore=tmpmax;

    //    cout << DYN_maxscore << endl << DYN_lastrc_maxscore << endl;

  }

  // Special debugon();
//  if(DYN_len_seq1==560
//     || DYN_len_seq2==560){
//    {
//	cout << "Matlastrow:\n";
//	int32 * ptrr=DYN_simmatrix+(DYN_len_seq1)*(DYN_len_seq2+1);
//	for(int32 i=0; i<=DYN_len_seq2;i++, ptrr++){
//	  cout << *ptrr << "\t";
//	}
//	cout << "\n";
//    }
//    {
//	cout << "Matlastcol:\n";
//	int32 * ptrc=DYN_simmatrix+DYN_len_seq2;
//	for(int32 i=0; i<=DYN_len_seq1;i++, ptrc+=DYN_len_seq2+1){
//	  cout << *ptrc << "\t";
//	}
//	cout << endl;
//    }
//  }

  FUNCEND();
}

//*/


void Dynamic::coutWhatWasGiven()
{
  cout << "Dynamic\n--------\nUh oh ... hunting a bug, aren't you?\n";
  if(DYN_sequence1 != nullptr) {
    cout << "Seq1: " << DYN_sequence1;
  }else{
    cout << "Seq1: nullptr";
  }
  if(DYN_sequence1 != nullptr) {
    cout << "\nSeq2: " << DYN_sequence2;
  }else{
    cout << "\nSeq2: nullptr";
  }
  cout << "\nvalidseq: " << DYN_validseq
       << "\nmat calc: " << DYN_matrixcalculated
       << "\nUse offset? " << DYN_calcwithoffset
       << "\nExp offset: " << DYN_eoffset
       << "\nlbx: " << DYN_leftbandx
       << "\nrbx: " << DYN_rightbandx
       << "\nbs: " << DYN_bandsafety
       << "\nbw: " << DYN_bandwidth
       << "\n";
}
