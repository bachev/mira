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


#include "mira/align.H"

#include "errorhandling/errorhandling.H"

#include "mira/ads.H"


using std::cout;
using std::cerr;
using std::endl;


#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif

#define CLOCK_STEPS1


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


uint64 Align::AL_alloccount=0;



/*************************************************************************
 *
 *
 *
 *************************************************************************/

Align::Align(MIRAParameters * params): Dynamic(params)
{
  FUNCSTART("Align::Align(MIRAParameters * params): Dynamic(params)");

  AL_miraparams=params;

  AL_tmpads=nullptr;
  AL_alseq1=nullptr;
  AL_alseq2=nullptr;
  AL_as12size=0;

  AL_userle=false;
  AL_rle_create_packed_align=false;
  AL_rle_create_stdleft_align=true;
  AL_rle_create_nonstdright_align=false;

  init();
  resetTimings();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::init()
{
  FUNCSTART("void Align::init()");

  AL_no_solutions=0;
  AL_no_diff_solutions=0;
  AL_max_relscore=0;

  AL_mpcache_dyn_score_multiplier=0;
  AL_mpcache_dyn_score_gap=0;
  AL_mpcache_al_max_cutoff=0;
  AL_mpcache_al_min_score=0;
  AL_mpcache_al_min_overlap=0;
  AL_mpcache_al_min_relscore=0;

  AL_mpset_al_min_overlap=0;
  AL_mpset_al_min_relscore=0;

  AL_valid=131;

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

Align::~Align()
{
  FUNCSTART("Align::~Align()");

  if(AL_tmpads != nullptr) delete AL_tmpads;
  if(AL_alseq1 != nullptr) delete [] AL_alseq1;
  if(AL_alseq2 != nullptr) delete [] AL_alseq2;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::acquireSequences(const char * seq1, uint32 len1, const char * seq2, uint32 len2, int32 id1, int32 id2, int8 id1dir, int8 id2dir, bool calcwithoffset, int32 expectedoffset)
{
  FUNCSTART("Align::acquireSequences(const char * seq1, uint32 len1, const char * seq2, uint32 len2, int32 id1, int32 id2, int8 id1dir, int8 id2dir, bool calcwithoffset, int32 expectedoffset)");

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  if(AL_userle){
    auto newoffset=expectedoffset;
    pa_packSeqToRLE(seq1,len1,AL_rles1, AL_rlev1, newoffset);
    if(newoffset>=0) expectedoffset=newoffset;
    newoffset=-expectedoffset;
    pa_packSeqToRLE(seq2,len2,AL_rles2, AL_rlev2, newoffset);
    if(newoffset>=0) expectedoffset=-newoffset;
    seq1=AL_rles1.c_str();
    len1=static_cast<uint32>(AL_rles1.size());
    seq2=AL_rles2.c_str();
    len2=static_cast<uint32>(AL_rles2.size());
  }else{
    AL_rles1.clear();
    AL_rles2.clear();
    AL_rlev1.clear();
    AL_rlev2.clear();
  }

  Dynamic::setSequences(seq1, len1, seq2, len2, calcwithoffset, expectedoffset);

  AL_id1=id1;
  AL_id2=id2;
  AL_id1dir=id1dir;
  AL_id2dir=id2dir;

  init();

#ifdef CLOCK_STEPS1
  AL_timing_acquires+=diffsuseconds(tv);
#endif

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::discard()
{
  FUNCSTART("Align::discard()");

  if(AL_valid != 131) {
    cerr << "AL_valid not valid?!\n";
    exit(0);
  }

  Dynamic::discard();

  if(AL_alseq1!=nullptr) delete [] AL_alseq1;
  if(AL_alseq2!=nullptr) delete [] AL_alseq2;
  AL_alseq1=nullptr;
  AL_alseq2=nullptr;
  AL_as12size=0;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::prepareAlign(std::list<AlignedDualSeq> * adslist)
{
  FUNCSTART("void Align::prepareAlign(std::list<AlignedDualSeq> * adslist)");

  BUGIFTHROW(!DYN_validseq,": Programming error. Tried to align without proper initialisation of Align-object.");

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  AL_adslist=adslist;
  if(AL_tmpads==nullptr) AL_tmpads= new AlignedDualSeq(AL_miraparams);

  AL_align_maxlen=DYN_len_seq1+DYN_len_seq2+1;

  if(AL_as12size < AL_align_maxlen){
    if(AL_alseq1 != nullptr) delete [] AL_alseq1;
    if(AL_alseq2 != nullptr) delete [] AL_alseq2;

    AL_as12size=AL_align_maxlen;
    if(AL_as12size<2000) AL_as12size=2000;

    AL_alseq1= new char[AL_as12size];
    AL_alseq2= new char[AL_as12size];
    AL_alloccount+=2;
  }

  AL_no_solutions=0;
  AL_no_diff_solutions=0;
  AL_max_relscore=0;

  AL_error_hit_band=false;

  // clear the aligned sequences
  // TODO: check: do I need this???
  memset(AL_alseq1,0,AL_align_maxlen);
  memset(AL_alseq2,0,AL_align_maxlen);

  AL_allen=AL_align_maxlen-1;

  AL_seq1ptr=DYN_sequence1+DYN_len_seq1-1;
  AL_seq2ptr=DYN_sequence2+DYN_len_seq2-1;


#ifdef CLOCK_STEPS1
  AL_timing_prepalign+=diffsuseconds(tv);
#endif

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::setRAlignParams()
{
  align_parameters const & AL_params = AL_miraparams->getAlignParams();
  dynamic_parameters const & DYN_params = DYN_miraparams->getDynamicParams();

  if(AL_mpset_al_min_relscore!=0){
    AL_mpcache_al_min_relscore=AL_mpset_al_min_relscore;
  }else{
    AL_mpcache_al_min_relscore=AL_params.al_min_relscore;
  }
  if(AL_mpset_al_min_overlap!=0){
    AL_mpcache_al_min_overlap=AL_mpset_al_min_overlap;
  }else{
    AL_mpcache_al_min_overlap=AL_params.al_min_overlap;
  }

  AL_mpcache_al_min_score=AL_params.al_min_score;
  AL_mpcache_al_max_cutoff=AL_params.al_max_cutoff;

  AL_mpcache_dyn_score_gap=DYN_params.dyn_score_gap;
  AL_mpcache_dyn_score_multiplier=DYN_params.dyn_score_multiplier;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::simpleAlign(std::list<AlignedDualSeq> * adslist)
{
  FUNCSTART("void Align::simpleAlign(std::list<AlignedDualSeq> * adslist)");

  if(!DYN_matrixcalculated) Dynamic::computeMatrix();

  prepareAlign(adslist);

  AL_new_solution=1;
  AL_cutoff_counter=0;
  AL_minbanddistance=1000000;

  setRAlignParams();

  rAlign(DYN_len_seq1, DYN_len_seq2,'d',false);

  //delete AL_tmpads;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::fullAlign(std::list<AlignedDualSeq> * adslist)
{
  FUNCSTART("void Align::fullAlign(std::list<AlignedDualSeq> * adslist)");

  if(!DYN_matrixcalculated) Dynamic::computeMatrix();

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  try {
    prepareAlign(adslist);

    setRAlignParams();

    termAlignFancy();
  }
  catch(Notify n){
    cout << "Full align failed!"
	 << endl;
    coutWhatWasGiven();
    n.handleError(THISFUNC);
  }


#ifdef CLOCK_STEPS1
  AL_timing_fullalign+=diffsuseconds(tv);
#endif

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define AL_TATRACE
//#define CEBUG(bla)   {cout << bla; cout.flush();}

//#define CEBUG(bla)   {if(AL_id1==-1 && AL_id2==14200) cout << bla; cout.flush();}


// looks like multiplier of 10000 rather than 100 is more sure

void Align::termAlign()
{
  FUNCSTART("void Align::termAlign()");

  align_parameters const & AL_params = AL_miraparams->getAlignParams();

  // FIXME: bad fix. if both are 1, termAlign will read across memborders
  // Albeit this should never happen, I'll fix termAlign later
  if(AL_params.al_min_score==1 && AL_params.al_min_overlap == 1){
    align_parameters & AL_ncparams = const_cast<align_parameters &>(AL_miraparams->getAlignParams());
    AL_ncparams.al_min_score=2;
    AL_ncparams.al_min_overlap=2;
  }

  // simple tactics now: just take the best value from last column
  //  or row

  int32 maxvalrow=-1;
  int32 maxvalcol=maxvalrow;

  uint32 maxvalrowpos=0;
  uint32 maxvalcolpos=0;
  {
    CEBUG("\n\nLast row (left to right):\n");
    for(uint32 i=1; i<=DYN_len_seq2; ++i) {
      int32 vAL_actpos=DYN_simmatrix[((DYN_len_seq1)*(DYN_len_seq2+1)+i)];
      CEBUG("R-Looking: " << i << "\tval:" << vAL_actpos << '\n');
      int32 relscore=(vAL_actpos*10000)/std::min(i,DYN_len_seq1);
      if(relscore >= AL_params.al_min_relscore
	 && vAL_actpos>=AL_params.al_min_score
	 && vAL_actpos>=maxvalrow){
	maxvalrow=vAL_actpos;
	maxvalrowpos=i;
      }
    }

    CEBUG("Last column:\n");
    for(uint32 i=1; i<=DYN_len_seq1; ++i) {
      int32 vAL_actpos=(DYN_simmatrix[i*(DYN_len_seq2+1)+DYN_len_seq2]);
      CEBUG("C-Looking: " << i << "\tval: " << vAL_actpos << '\n');
      int32 relscore=(vAL_actpos*10000)/std::min(i,DYN_len_seq2);
      if(relscore >= AL_params.al_min_relscore
	 && vAL_actpos>=AL_params.al_min_score
	 && vAL_actpos >= maxvalcol){
	maxvalcol=vAL_actpos;
	maxvalcolpos=i;
      }
    }
  }

  CEBUG("Maxvalrow: " << maxvalrow << "\tpos" << maxvalrowpos << "\n");
  CEBUG("Maxvalcol: " << maxvalcol << "\tpos" << maxvalcolpos << "\n");

  AL_new_solution=1;
  AL_cutoff_counter=0;
  AL_minbanddistance=1000000;

  //if(maxvalrow==maxvalcol
  //   && maxvalrowpos+1==DYN_len_seq2
  //   && maxvalcolpos+1==DYN_len_seq1){
  //  // we start in the corner
  //  rAlign(DYN_len_seq1, DYN_len_seq2,'d');
  //}else
  if(maxvalrow>maxvalcol){
    if(maxvalrowpos>0) {
      CEBUG("maxvalrowpos: " << maxvalrowpos << '\t' << "DYN_len_seq2: " << DYN_len_seq2 << endl);
      for(uint32 i=maxvalrowpos; i<DYN_len_seq2; ++i){
	AL_alseq1[--AL_allen]=' ';
	AL_alseq2[AL_allen]=*AL_seq2ptr--;
      }
      CEBUG("DoAlign");
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      rAlign(DYN_len_seq1, maxvalrowpos,'d',false);
#ifdef CLOCK_STEPS1
      AL_timing_raligntot+=diffsuseconds(tv);
#endif
    }
  } else {
    if(maxvalcolpos>0) {
      CEBUG("maxvalcolpos: " << maxvalcolpos << '\t' << "DYN_len_seq1: " << DYN_len_seq1 << endl);
      for(uint32 i=maxvalcolpos; i<DYN_len_seq1; ++i){
	AL_alseq1[--AL_allen]=*AL_seq1ptr--;
	AL_alseq2[AL_allen]=' ';
      }
      CEBUG("DoAlign");
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      rAlign(maxvalcolpos, DYN_len_seq2,'d',false);
#ifdef CLOCK_STEPS1
      AL_timing_raligntot+=diffsuseconds(tv);
#endif
    }
  }
  CEBUG("I'm out! AL_adslist->size(): " << AL_adslist->size() << endl;);

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * Fancy version of termAlign()
 * Like above, starts alignment from best value in last column or row
 *  But when it's through with initial alignment, looks left and right
 *  for "next best" secondary starting points defined as in:
 *  - within x% of the maximum
 *  - breaking monotonic descent of values.
 *
 *  1420                 #
 *                      # #
 *                    ##   #   #
 *  1390           # #      ### #
 *  1375          # #            #
 *               #                #
 *            ...                  ...
 *                       ^
 *                      initial
 *                 ^           ^
 *               secondary   secondary
 *
 * Reason:
 * I did (and will not) not implement affine gap cost. But in mapping
 *  applications, reads with longer deletions near end get a mathematically
 *  'better' value for a collapsed alignment than for the potential
 *  'true' alignment.
 *
 *  E.g.:
 *   ID1:-1  CTCCTACTTCTACTTCT424423TCTTCTTCTTCGTCCCCGACA
 *   ID2:10  CTCCTACTTCTACTTCTTCTTCTTCTTCGTCC
 *   Diff                     XXXXXX     X  X
 *
 *  instead of
 *   ID1:-1  CTCCTACTTCTACTTCT424423TCTTCTTCTTCGTCCCCGACA
 *   ID2:10  CTCCTACTTCTACTTCT******TCTTCTTCTTCGTCC
 *   Diff                     XXXXXX
 *
 * The fancy version will find alternative starts (max one left and one right)
 *  and let the ads::score() sort through the possibilities, as ads::score()
 *  has a modifier for affine gap costs.
 *
 * IMPORTANT NOTE: given I want this only in mapping, I implement only part of it
 *   Only looks in last column because
 *     1. reference cut-out within a contigm will be longer and in ID1
 *     2. alternative start sites will be to the "increasing-i" side
 *
 * This function is such an ugly hack ...
 *
 *************************************************************************/

//#define AL_TATRACE
//#define CEBUG(bla)   {cout << bla; cout.flush();}

//#define CEBUG(bla)   {if(AL_id1==-1 && AL_id2==14200) cout << bla; cout.flush();}


// looks like multiplier of 10000 rather than 100 is more sure

void Align::termAlignFancy()
{
  FUNCSTART("void Align::termAlignFancy()");

  align_parameters const & AL_params = AL_miraparams->getAlignParams();

  // FIXME: bad fix. if both are 1, termAlign will read across memborders
  // Albeit this should never happen, I'll fix termAlign later
  if(AL_params.al_min_score==1 && AL_params.al_min_overlap == 1){
    align_parameters & AL_ncparams = const_cast<align_parameters &>(AL_miraparams->getAlignParams());
    AL_ncparams.al_min_score=2;
    AL_ncparams.al_min_overlap=2;
  }

  // simple tactics now: just take the best value from last column
  //  or row

  int32 maxvalrow=-1;
  int32 maxvalcol=maxvalrow;

  uint32 maxvalrowpos=0;
  uint32 maxvalcolpos=0;
  {
    CEBUG("\n\nLast row (left to right):\n");
    for(uint32 i=1; i<=DYN_len_seq2; ++i) {
      int32 vAL_actpos=DYN_simmatrix[((DYN_len_seq1)*(DYN_len_seq2+1)+i)];
      CEBUG("R-Looking: " << i << "\tval:" << vAL_actpos << '\n');
      int32 relscore=(vAL_actpos*10000)/std::min(i,DYN_len_seq1);
      if(relscore >= AL_params.al_min_relscore
	 && vAL_actpos>=AL_params.al_min_score
	 && vAL_actpos>=maxvalrow){
	maxvalrow=vAL_actpos;
	maxvalrowpos=i;
      }
    }

    CEBUG("Last column:\n");
    for(uint32 i=1; i<=DYN_len_seq1; ++i) {
      int32 vAL_actpos=(DYN_simmatrix[i*(DYN_len_seq2+1)+DYN_len_seq2]);
      CEBUG("C-Looking: " << i << "\tval: " << vAL_actpos << '\n');
      int32 relscore=(vAL_actpos*10000)/std::min(i,DYN_len_seq2);
      if(relscore >= AL_params.al_min_relscore
	 && vAL_actpos>=AL_params.al_min_score
	 && vAL_actpos >= maxvalcol){
	maxvalcol=vAL_actpos;
	maxvalcolpos=i;
      }
    }
  }

  CEBUG("Maxvalrow: " << maxvalrow << "\tpos" << maxvalrowpos << "\n");
  CEBUG("Maxvalcol: " << maxvalcol << "\tpos" << maxvalcolpos << "\n");

  AL_new_solution=1;
  AL_cutoff_counter=0;
  AL_minbanddistance=1000000;

  //if(maxvalrow==maxvalcol
  //   && maxvalrowpos+1==DYN_len_seq2
  //   && maxvalcolpos+1==DYN_len_seq1){
  //  // we start in the corner
  //  rAlign(DYN_len_seq1, DYN_len_seq2,'d');
  //}else
  if(maxvalrow>maxvalcol){
    if(maxvalrowpos>0) {
      CEBUG("maxvalrowpos: " << maxvalrowpos << '\t' << "DYN_len_seq2: " << DYN_len_seq2 << endl);
      for(uint32 i=maxvalrowpos; i<DYN_len_seq2; ++i){
	AL_alseq1[--AL_allen]=' ';
	AL_alseq2[AL_allen]=*AL_seq2ptr--;
      }
      CEBUG("DoAlign");
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      rAlign(DYN_len_seq1, maxvalrowpos,'d',false);
#ifdef CLOCK_STEPS1
      AL_timing_raligntot+=diffsuseconds(tv);
#endif
    }
  } else {
    if(maxvalcolpos>0) {
      CEBUG("maxvalcolpos: " << maxvalcolpos << '\t' << "DYN_len_seq1: " << DYN_len_seq1 << endl);

      // atm 10%
      auto threshold = maxvalcol - maxvalcol / 10;

      // search alternative start col
      uint32 altcolpos1=0;
      auto oldval=maxvalcol;
      for(uint32 i=maxvalcolpos; i<=DYN_len_seq1; ++i) {
	int32 vAL_actpos=(DYN_simmatrix[i*(DYN_len_seq2+1)+DYN_len_seq2]);
	CEBUG("C-Looking a1: " << i << "\tval: " << vAL_actpos << '\n');
	if (vAL_actpos < threshold) break;
	if (vAL_actpos > oldval) {
	  oldval = vAL_actpos;
	  while (true) {
	    if (i > DYN_len_seq1) break;
	    ++i;
	    vAL_actpos=(DYN_simmatrix[i*(DYN_len_seq2+1)+DYN_len_seq2]);
	    CEBUG("C-Looking a1+: " << i << "\tval: " << vAL_actpos << '\n');
	    if (vAL_actpos < oldval) break; // found second max
	    oldval = vAL_actpos;
	  }
	  altcolpos1 = i - 1;
	  CEBUG("Found altcolpos1: " << altcolpos1 << '\n');
	  break;
	}
	oldval=vAL_actpos;
      }

      for(uint32 i=maxvalcolpos; i<DYN_len_seq1; ++i){
	AL_alseq1[--AL_allen]=*AL_seq1ptr--;
	AL_alseq2[AL_allen]=' ';
      }
      CEBUG("DoAlign");
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      rAlign(maxvalcolpos, DYN_len_seq2,'d',false);

      if (altcolpos1) {
	AL_allen=AL_align_maxlen-1;
	AL_seq1ptr=DYN_sequence1+DYN_len_seq1-1;
	AL_seq2ptr=DYN_sequence2+DYN_len_seq2-1;
	AL_new_solution=1;
	AL_cutoff_counter=0;
	AL_minbanddistance=1000000;
	for(uint32 i=altcolpos1; i<DYN_len_seq1; ++i){
	  AL_alseq1[--AL_allen]=*AL_seq1ptr--;
	  AL_alseq2[AL_allen]=' ';
	}
	CEBUG("DoAlignAlt");
	rAlign(altcolpos1, DYN_len_seq2,'d',false);
      }
#ifdef CLOCK_STEPS1
      AL_timing_raligntot+=diffsuseconds(tv);
#endif
    }
  }
  CEBUG("I'm out! AL_adslist->size(): " << AL_adslist->size() << endl;);

  FUNCEND();
}
//#define CEBUG(bla)





/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define FUNCSTART(bla)  static const char * THISFUNC = bla"  ";
//#define FUNCTRACE(bla) { cout << THISFUNC << bla; cout.flush();}
//#define FUNCEND()

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Align::rAlign(uint32 irow, uint32 jcol, char lastdir, bool hadn)
{
  FUNCSTART("void Align::rAlign(uint32 irow, uint32 jcol, char lastdir, bool hadn)");

//  align_parameters const & AL_params = AL_miraparams->getAlignParams();
//  dynamic_parameters const & DYN_params = DYN_miraparams->getDynamicParams();

  uint32 mll=DYN_len_seq2+1;

  //  prefetchrl(&DYN_match_matrix[static_cast<uint8>(*AL_seq1ptr)][static_cast<uint8>(*AL_seq2ptr)]);
  //  prefetchrl(&DYN_simmatrix[i*mll+j-1]);
  //  prefetchrl(&DYN_simmatrix[i*mll+j]);
  //  prefetchrl(&DYN_simmatrix[(i-1)*mll+j-1]);
  //  prefetchrl(&DYN_simmatrix[(i-1)*mll+j]);

  CEBUG("Dong!\n");
  CEBUG("irow: " << irow << "\tjcol: " << jcol << "\tAL_allen: " << AL_allen<< endl);
  CEBUG("s[i,j]:" << DYN_simmatrix[irow*(DYN_len_seq2+1)+jcol]<< endl);
  if(irow>0 && jcol>0){
    CEBUG("s[i-1,j-1]:" << DYN_simmatrix[((irow-1)*(DYN_len_seq2+1))+jcol-1] << endl);
  }

  bool hasn=false;
  if(AL_seq1ptr>=DYN_sequence1 && *AL_seq1ptr){
    CEBUG("seq1ptr points on: " << *AL_seq1ptr << endl);
    hasn=(*AL_seq1ptr=='N');
  }else{
    CEBUG("seq1ptr points on: nullptr\n");
  }
  if(AL_seq2ptr>=DYN_sequence2 && *AL_seq2ptr){
    CEBUG("seq2ptr points on: " << *AL_seq2ptr << endl);
    hasn=hasn | (*AL_seq2ptr=='N');
  }else{
    CEBUG("seq2ptr points on: nullptr\n");
  }
  if(unlikely(AL_allen>AL_align_maxlen)) {
    cerr << "allen:: "<< AL_allen;
    MIRANOTIFY(Notify::INTERNAL, ": FOOOOOO!.") ;
  }

  if(unlikely(AL_cutoff_counter==AL_mpcache_al_max_cutoff)) {
    CEBUG("back... (because of cutoff)\n");
    return;
  }
  if(unlikely(AL_error_hit_band)) {
    CEBUG("back... (because of band hit)\n");
    return;
  }

  if(unlikely(irow==0 && jcol==0)) {
    if(AL_seq1ptr==DYN_sequence1-1 && AL_seq2ptr==DYN_sequence2-1){
      const char * adss1ptr=AL_alseq1+AL_allen;
      const char * adss2ptr=AL_alseq2+AL_allen;
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      AL_tmpads->acquireSequences(adss1ptr, adss2ptr,
				  AL_id1, AL_id2,
				  AL_id1dir, AL_id2dir,
				  AL_enforce_clean_ends,
				  AL_affine_gap_scorees,
				  DYN_bandwidth, AL_minbanddistance);
#ifdef CLOCK_STEPS1
      AL_timing_ra_adsacquire+=diffsuseconds(tv);
#endif
      CEBUG("Solution\n");
      CEBUG(*AL_tmpads);

      if(AL_tmpads->getScore() >= static_cast<int32>(AL_mpcache_al_min_score / AL_mpcache_dyn_score_multiplier)
	 && ((AL_tmpads->getOverlapLen() >= static_cast<uint32>(AL_mpcache_al_min_overlap)
	      && AL_tmpads->getScoreRatio() >= static_cast<int32>(AL_mpcache_al_min_relscore))
	     || (true
		 && AL_tmpads->getOverlapLen() >= 17
		 && AL_tmpads->getScoreRatio() == 100))) {

#ifdef CLOCK_STEPS1
	gettimeofday(&tv,nullptr);
#endif
	AL_adslist->push_back(*AL_tmpads);
#ifdef CLOCK_STEPS1
	AL_timing_ra_adslist+=diffsuseconds(tv);
#endif

	if(!AL_rlev1.empty()){
	  if(AL_rle_create_stdleft_align){
	    ra_expandRLEAlignments(adss1ptr,adss2ptr,true);
#ifdef CLOCK_STEPS1
	    gettimeofday(&tv,nullptr);
#endif
	    AL_tmpads->acquireSequences(AL_unrles1.c_str(), AL_unrles2.c_str(),
					AL_id1, AL_id2,
					AL_id1dir, AL_id2dir,
					AL_enforce_clean_ends,
					AL_affine_gap_scorees,
					DYN_bandwidth, AL_minbanddistance);
#ifdef CLOCK_STEPS1
	    AL_timing_ra_adsacquire+=diffsuseconds(tv);
#endif
#ifdef CLOCK_STEPS1
	    gettimeofday(&tv,nullptr);
#endif
	    AL_adslist->push_back(*AL_tmpads);
#ifdef CLOCK_STEPS1
	    AL_timing_ra_adslist+=diffsuseconds(tv);
#endif
	  }
	  if(AL_rle_create_nonstdright_align){
	    ra_expandRLEAlignments(adss1ptr,adss2ptr,false);
#ifdef CLOCK_STEPS1
	    gettimeofday(&tv,nullptr);
#endif
	    AL_tmpads->acquireSequences(AL_unrles1.c_str(), AL_unrles2.c_str(),
					AL_id1, AL_id2,
					AL_id1dir, AL_id2dir,
					AL_enforce_clean_ends,
					AL_affine_gap_scorees,
					DYN_bandwidth, AL_minbanddistance);
#ifdef CLOCK_STEPS1
	    AL_timing_ra_adsacquire+=diffsuseconds(tv);
#endif
#ifdef CLOCK_STEPS1
	    gettimeofday(&tv,nullptr);
#endif
	    AL_adslist->push_back(*AL_tmpads);
#ifdef CLOCK_STEPS1
	    AL_timing_ra_adslist+=diffsuseconds(tv);
#endif
	  }
	}


      }

      ++AL_cutoff_counter;
      ++AL_no_solutions;
      if(AL_new_solution){
	AL_new_solution=0;
	++AL_no_diff_solutions;
      }
    }
    CEBUG("back...\n");
    return;
  }

  if(irow==0){
    CEBUG("irow=0 left ...\n");

    AL_alseq1[--AL_allen]=' ';
    AL_alseq2[AL_allen]=*AL_seq2ptr--;
    rAlign(irow,jcol-1,'l',hasn);
//    AL_alseq1[AL_allen]='^';
//    AL_alseq2[AL_allen++]='^';
    AL_allen++;
    AL_seq2ptr++;
  }else if(jcol==0){

    CEBUG("jcol=0 up ...\n");

    AL_alseq1[--AL_allen]=*AL_seq1ptr--;
    AL_alseq2[AL_allen]=' ';
    rAlign(irow-1,jcol,'u',hasn);
//    AL_alseq1[AL_allen]='%';
//    AL_alseq2[AL_allen++]='%';
    AL_allen++;
    AL_seq1ptr++;
  } else {
    // minimum distance from left band ...
    AL_minbanddistance= std::min(AL_minbanddistance,static_cast<int32>(jcol)-(DYN_leftbandx+static_cast<int32>(irow)));
    BUGIFTHROW(AL_minbanddistance<0,"left AL_minbanddistance " << AL_minbanddistance << " < 0 ???\nirow: " << irow << " jcol: " << jcol << " lbx: " << DYN_leftbandx);
    // ... and right band
    AL_minbanddistance= std::min(AL_minbanddistance,DYN_rightbandx+static_cast<int32>(irow)-static_cast<int32>(jcol));
    BUGIFTHROW(AL_minbanddistance<0,"right AL_minbanddistance " << AL_minbanddistance << " < 0 ???\nirow: " << irow << " jcol: " << jcol << " rbx: " << DYN_rightbandx);

    if(AL_minbanddistance <= DYN_bandsafety){
      AL_error_hit_band=true;
      CEBUG("HIT BAND LIMIT " << DYN_bandsafety << " : " << static_cast<int32>(jcol)-(DYN_leftbandx+static_cast<int32>(irow)) << " " << DYN_rightbandx+static_cast<int32>(irow)-static_cast<int32>(jcol) << endl);
    }else{
      int32 vgl=DYN_match_matrix[static_cast<uint8>(*AL_seq1ptr)][static_cast<uint8>(*AL_seq2ptr)];

      CEBUG("vgl: " << vgl << endl);

      // precompute possibilities, optimised for memory access
      // at the same time, check whether we hit a band limit
      if(DYN_simmatrix[irow*mll+jcol-1]==DYN_BANDLIMIT) AL_error_hit_band=true;
      bool leftok=(DYN_simmatrix[irow*mll+jcol-1]+AL_mpcache_dyn_score_gap==DYN_simmatrix[irow*mll+jcol]);
      bool diagok=(DYN_simmatrix[(irow-1)*mll+jcol-1]+vgl==DYN_simmatrix[irow*mll+jcol]);
      if(DYN_simmatrix[(irow-1)*mll+jcol]==DYN_BANDLIMIT) AL_error_hit_band=true;
      bool upok=(DYN_simmatrix[(irow-1)*mll+jcol]+AL_mpcache_dyn_score_gap==DYN_simmatrix[irow*mll+jcol]);

      // This is a special rule for strobe sequencing
      // If a stretch of 'N' occurs and there's the choice between
      //  a match and a gap, choose gap
      // this prevents things like:
      //    ..nnnnnnnaaaaaaaaat.....
      //    ..nnnnnnn*****nnnnt.....
      // and make it
      //    ..nnnnnnnaaaaaaaaat.....
      //    ..nn*****nnnnnnnnnt.....
      //
      // i.e., align gaps to N

      // BaCh 21.04.2010: not that good this idea.
      // - new PacBio dark strobe editing strategy doesn't need this
      // - no directly visible improvement
      // -> remove.
      //if(*AL_seq1ptr=='N'
      //   && *AL_seq2ptr=='N'
      //   && diagok
      //   && (leftok || upok)){
      //  diagok=false;
      //}

      // prevent
      //    ..cccctccaccgaatgcctaa
      //    ..cccctcc****a-----------------
      // and make it
      //    ..cccctccaccgaatgcctaa
      //    ..cccctcca****-----------------

      if(hadn
	 && diagok
	 && (leftok || upok)){
	diagok=false;
      }

      // now, look which possibilities we have
      //  if a diagonal and (up or left) are equal, then take
      //  the same direction as last time
      // this prevents things like:
      //    ..tggaaaaaaaaat.....
      //    ..t*g*********t.....
      // and make it
      //    ..tggaaaaaaaaat.....
      //    ..tg**********t.....


      //if(diagok){
      if(upok){
	if(lastdir=='u') {
	  diagok=false;
	}
      }
      if(leftok){
	if(lastdir=='l') {
	  diagok=false;
	  upok=false;
	}
      }
      //}

      if(diagok){

	CEBUG("diagonal ...\n");

	AL_alseq1[--AL_allen]=*AL_seq1ptr--;
	AL_alseq2[AL_allen]=*AL_seq2ptr--;
	rAlign(irow-1,jcol-1,'d',hasn);
//	AL_alseq1[AL_allen]='/';
//	AL_alseq2[AL_allen++]='/';
	AL_allen++;
	AL_seq1ptr++;
	AL_seq2ptr++;
      }

      if(upok && !AL_error_hit_band){

	CEBUG("up ...\n");

	AL_alseq1[--AL_allen]=*AL_seq1ptr--;
	AL_alseq2[AL_allen]='*';
	rAlign(irow-1,jcol,'u',hasn);
//	AL_alseq1[AL_allen]='$';
//	AL_alseq2[AL_allen++]='$';
	AL_allen++;
	AL_seq1ptr++;
      }

      if(leftok && !AL_error_hit_band){

	CEBUG("left ...\n");

	AL_alseq1[--AL_allen]='*';
	AL_alseq2[AL_allen]=*AL_seq2ptr--;
	rAlign(irow,jcol-1,'l',hasn);
//	AL_alseq1[AL_allen]='#';
//	AL_alseq2[AL_allen++]='#';
	AL_allen++;
	AL_seq2ptr++;
      }
    }
  }

  CEBUG("back...\n");


  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::pa_packSeqToRLE(const char * seq, uint32 len, std::string & rles, std::vector<uint32> & rlev, int32 & eoffset)
{
  FUNCSTART("void Align::pa_packSeqToRLE(const char * seq, uint32 len, std::string & rles, std::vector<uint32> rlev)");

  rles.clear();
  rlev.clear();

  uint32 run=0;
  char runchar=*seq;
  for(uint32 ri=0; ri<len; ++seq, ++ri){
    if(*seq==runchar){
      ++run;
    }else{
      rles.push_back(runchar);
      runchar=*seq;
      rlev.push_back(run);
      run=1;
    }
    if(static_cast<int32>(ri)==eoffset) eoffset=static_cast<uint32>(rles.size());
  }
  rles.push_back(runchar);
  rlev.push_back(run);
}


/*************************************************************************
 *
 * Transform packed RLE alignment
 *    s1        AAACGTC
 *    s2        AAA**TC
 * to (shakedirstdleft==true)
 *    s1        **AAACGTTTC
 *    s2        AAAAA**TTTC
 * or (shakedirstdleft==false)
 *    s1        AAA**CGTTTC
 *    s2        AAAAA**TTTC
 *
 * according to RLE values stored for reads
 *
 *************************************************************************/

void Align::ra_expandRLEAlignments(const char * s1, const char * s2, bool shakedirstdleft)
{
  FUNCSTART("void Align::ra_expandRLEAlignments(const char * s1, const char * s2, bool shakedirstdleft)");

  AL_unrles1.clear();
  AL_unrles2.clear();
  auto vs1I=AL_rlev1.cbegin();
  auto vs2I=AL_rlev2.cbegin();
  for(; *s1 ; ++s1, ++s2){
    uint32 lens1=1;
    uint32 lens2=1;
    if(*s1==' '){
      if(*s2!='*'){
	lens2=*vs2I;
	++vs2I;
      }
      lens1=lens2;
    }else if(*s2==' '){
      if(*s1!='*'){
	lens1=*vs1I;
	++vs1I;
      }
      lens2=lens1;
    }else{
      if(*s1!='*'){
	lens1=*vs1I;
	++vs1I;
      }
      if(*s2!='*'){
	lens2=*vs2I;
	++vs2I;
      }
    }
    uint32 maxlen=std::max(lens1,lens2);
    uint32 run=0;
    if(shakedirstdleft){
      for(; run<maxlen-lens1; ++run) AL_unrles1.push_back('*');
      for(; run<maxlen; ++run) AL_unrles1.push_back(*s1);
      run=0;
      for(; run<maxlen-lens2; ++run) AL_unrles2.push_back('*');
      for(; run<maxlen; ++run) AL_unrles2.push_back(*s2);
    }else{
      for(; run<lens1; ++run) AL_unrles1.push_back(*s1);
      for(; run<maxlen; ++run) AL_unrles1.push_back('*');
      run=0;
      for(; run<lens2; ++run) AL_unrles2.push_back(*s2);
      for(; run<maxlen; ++run) AL_unrles2.push_back('*');
    }
  }
  BUGIFTHROW(vs1I!=AL_rlev1.cend(),"vs1I!=AL_rlev1.cend()");
  BUGIFTHROW(vs2I!=AL_rlev2.cend(),"vs2I!=AL_rlev2.cend()");
}

/*************************************************************************
 *
 * Transform
 *    s1        AAA**CGT
 *    s2        AAAAA**T
 * to
 *    s1        AAACGT
 *    s2        AAAAAT
 *
 * Really? Not sure I want that for PacBio.
 *
 *************************************************************************/

/*
void Align::ra_fitRLEIndents()
{
}
//*/

void Align::coutWhatWasGiven()
{
  Dynamic::coutWhatWasGiven();

  cout << "Align\n------"
       << "\nAL_id1: " << AL_id1
       << "\nAL_id2: " << AL_id2
       << "\nAL_id1dir: " << static_cast<int16>(AL_id1dir)
       << "\nAL_id2dir: " << static_cast<int16>(AL_id2dir)
       << endl;
}


void Align::resetTimings()
{
  DYN_timing_seqcopy=0;
  DYN_timing_bswmatrix=0;
  DYN_timing_bswm_setup=0;
  DYN_timing_bswm_p1=0;
  DYN_timing_bswm_p2a=0;
  DYN_timing_bswm_p2b=0;
  DYN_timing_bswm_p3=0;
  DYN_timing_bswm_cleanband=0;
  AL_timing_acquires=0;
  AL_timing_fullalign=0;
  AL_timing_prepalign=0;
  AL_timing_raligntot=0;
  AL_timing_ra_adsacquire=0;
  AL_timing_ra_adslist=0;
}

void Align::dumpTimings()
{
  cout << "Align timing DYN seqcpy : " << DYN_timing_seqcopy << endl;
  cout << "Align timing DYN bsw su : " << DYN_timing_bswm_setup << endl;
  cout << "Align timing DYN bsw p1 : " << DYN_timing_bswm_p1 << endl;
  cout << "Align timing DYN bsw p2a: " << DYN_timing_bswm_p2a << endl;
  cout << "Align timing DYN bsw p2b: " << DYN_timing_bswm_p2b << endl;
  cout << "Align timing DYN bsw p3 : " << DYN_timing_bswm_p3 << endl;
  cout << "Align timing DYN bsw cb : " << DYN_timing_bswm_cleanband << endl;

  cout << "Align timing DYN bsw    : " << DYN_timing_bswmatrix << endl;
  cout << "Align timing AL acqu s  : " << AL_timing_acquires << endl;
  cout << "Align timing AL full    : " << AL_timing_fullalign << endl;
  cout << "Align timing AL prep    : " << AL_timing_prepalign << endl;
  cout << "Align timing AL ralignt : " << AL_timing_raligntot << endl;
  cout << "Align timing AL ralignc : " << AL_timing_raligntot-AL_timing_ra_adsacquire-AL_timing_ra_adslist << endl;
  cout << "Align timing AL ads a   : " << AL_timing_ra_adsacquire << endl;
  cout << "Align timing AL ads s   : " << AL_timing_ra_adslist << endl;
}
