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
 * 2 Sequences aligned, fact storage
 * AlignedDualSeq manages the data produced by the alignement routines. Two
 *  sequences have been aligned and are stored in an ADS object.
 *
 * The ADSFacts object is a smaller storage unit for the plain facts,
 *  smaller and therefore less memory consuming (454 data *sigh*)
 *
 */


#include "adsfacts.H"

#include "errorhandling/errorhandling.H"


std::ostream & operator<<(std::ostream & ostr, const AlignedDualSeqFacts &adsf)
{
  ostr << "ID1:" << adsf.ADSF_id1 << '\n';
  ostr << "ID2:" << adsf.ADSF_id2 << '\n';
  ostr << "Direction1: " << static_cast<int16>(adsf.getSequenceDirection(adsf.ADSF_id1)) << '\n';
  ostr << "Direction2: " << static_cast<int16>(adsf.getSequenceDirection(adsf.ADSF_id2)) << '\n';
  ostr << "Delta Seq2 to Seq1: " << adsf.ADSF_delta << '\n';
  ostr << "ID1 right delta: " << adsf.ADSF_id1_rightdelta << '\n';
  ostr << "ID2 right delta: " << adsf.ADSF_id2_rightdelta << '\n';
  ostr << "ID1 5p clean: " << adsf.ADSF_5pconmatch1*4
       << "\nID1 3p clean: " << adsf.ADSF_3pconmatch1*4
       << "\nID2 5p clean: " << adsf.ADSF_5pconmatch2*4
       << "\nID2 3p clean: " << adsf.ADSF_3pconmatch2*4 << '\n';
  ostr << "Overlap length: " << adsf.getOverlapLen() << '\n';
  ostr << "Total length: " << adsf.ADSF_total_len << '\n';
  ostr << "Score ratio: " << static_cast<uint16>(adsf.ADSF_score_ratio) << '\n';
  return ostr;
}

void AlignedDualSeqFacts::serialiseOut(std::ostream & ostr)
{
  ostr << ADSF_id1
       << '\t' << ADSF_id2
       << '\t' << static_cast<int16>(getSequenceDirection(ADSF_id1))
       << '\t' << static_cast<int16>(getSequenceDirection(ADSF_id2))
       << '\t' << ADSF_delta
       << '\t' << ADSF_id1_rightdelta
       << '\t' << ADSF_id2_rightdelta
       << '\t' << getOverlapLen()
       << '\t' << ADSF_total_len
       << '\t' << static_cast<uint16>(ADSF_score_ratio)
       << '\t' << static_cast<uint16>(ADSF_totalnonmatches)
       << '\t' << ADSF_5pconmatch1
       << '\t' << ADSF_3pconmatch1
       << '\t' << ADSF_5pconmatch2
       << '\t' << ADSF_3pconmatch2;

  return;
}

void AlignedDualSeqFacts::serialiseIn(std::istream & istr)
{
  istr >> ADSF_id1;
  istr >> ADSF_id2;

  ADSF_id1and2_directions=0;
  // no idea why I can't directly stream into an int8 *sigh*
  int32 tmp;
  istr >> tmp;
  if(tmp>0) ADSF_id1and2_directions=1;
  istr >> tmp;
  if(tmp>0) ADSF_id1and2_directions|=2;
  istr >> ADSF_delta;
  istr >> ADSF_id1_rightdelta;
  istr >> ADSF_id2_rightdelta;
  istr >> tmp;
  istr >> ADSF_total_len;
  istr >> tmp;
  ADSF_score_ratio=static_cast<int8>(tmp);

  // bitfields need temp variable
  istr >> tmp;
  ADSF_totalnonmatches=tmp;
  istr >> tmp;
  ADSF_5pconmatch1=tmp;
  istr >> tmp;
  ADSF_3pconmatch1=tmp;
  istr >> tmp;
  ADSF_5pconmatch2=tmp;
  istr >> tmp;
  ADSF_3pconmatch2=tmp;


  //cout << *this;
  return;
}


int8 AlignedDualSeqFacts::getSequenceDirection(readid_t id) const
{
  FUNCSTART("int8 AlignedDualSeqFacts::getSequenceDirection(readid_t id) const");
  if(id==ADSF_id1){
    FUNCEND();
    if(ADSF_id1and2_directions & 0x1) return 1;
    return -1;
  }else if(id==ADSF_id2){
    FUNCEND();
    if(ADSF_id1and2_directions & 0x2) return 1;
    return -1;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}

int32 AlignedDualSeqFacts::getOtherID(readid_t id) const
{
  FUNCSTART("uint32 ADS::getOtherID(ureadid_t id)");
  if(id==ADSF_id1){
    FUNCEND();
    return ADSF_id2;
  }else if(id==ADSF_id2){
    FUNCEND();
    return ADSF_id1;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}


uint32 AlignedDualSeqFacts::getOffsetInAlignment(readid_t id) const
{
  FUNCSTART("uint32 AlignedDualSeqFacts::getOffsetInAlignment(readid_t id) const");
  if(id==ADSF_id1){
    FUNCEND();
    return 0;
  }else if(id==ADSF_id2){
    FUNCEND();
    return ADSF_delta;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}

uint32 AlignedDualSeqFacts::getRightOffsetInAlignment(readid_t id) const
{
  FUNCSTART("uint32 ADS::getOffsetInAlignment(readid_t id)");
  if(id==ADSF_id1){
    FUNCEND();
    return ADSF_id1_rightdelta;
  }else if(id==ADSF_id2){
    FUNCEND();
    return ADSF_id2_rightdelta;
  }else{
    MIRANOTIFY(Notify::FATAL, "ID not in alignment.");
  }
}

void AlignedDualSeqFacts::publicinit(readid_t id1, readid_t id2, uint16 delta, uint16 id1_rightdelta, uint16 id2_rightdelta, uint16 total_len, int8 id1_direction, int8 id2_direction, int8 score_ratio, uint16 totalnonmatches, uint16 s5pcm1, uint16 s3pcm1, uint16 s5pcm2, uint16 s3pcm2)
{
  FUNCSTART("void AlignedDualSeqFacts::publicinit(readid_t id1, readid_t id2, uint16 delta, uint16 id1_rightdelta, uint16 id2_rightdelta, uint16 total_len, int8 id1_direction, int8 id2_direction, int8 score_ratio, uint16 totalnonmatches, uint16 s5pcm1, uint16 s3pcm1, uint16 s5pcm2, uint16 s3pcm2)");
  ADSF_id1=id1;
  ADSF_id2=id2;
  ADSF_delta=delta;
  ADSF_id1_rightdelta=id1_rightdelta;
  ADSF_id2_rightdelta=id2_rightdelta;
  ADSF_total_len=total_len;
  ADSF_5pconmatch1=(s5pcm1>=28) ? 7 : s5pcm1/4;
  ADSF_3pconmatch1=(s3pcm1>=28) ? 7 : s3pcm1/4;
  ADSF_5pconmatch2=(s5pcm2>=28) ? 7 : s5pcm2/4;
  ADSF_3pconmatch2=(s3pcm2>=28) ? 7 : s3pcm2/4;

  ADSF_id1and2_directions=0;
  if(id1_direction>0) ADSF_id1and2_directions=1;
  if(id2_direction>0) ADSF_id1and2_directions|=2;
  ADSF_score_ratio=score_ratio;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

uint32 AlignedDualSeqFacts::get5pLenContiguousMatch(readid_t id) const
{
  if(ADSF_id1==id) return ADSF_5pconmatch1*4;
  return ADSF_5pconmatch2*4;
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

uint32 AlignedDualSeqFacts::get3pLenContiguousMatch(readid_t id) const
{
  if(ADSF_id1==id) return ADSF_3pconmatch1*4;
  return ADSF_3pconmatch2*4;
}


uint16 AlignedDualSeqFacts::getRightDelta(readid_t id) const
{
  FUNCSTART("uint16 AlignedDualSeqFacts::getRightDelta(readid_t id) const");
  if(id==ADSF_id1){
    return ADSF_id1_rightdelta;
  }else if(unlikely(id!=ADSF_id2)){
    MIRANOTIFY(Notify::FATAL, "ID " << id << " not in alignment.");
  }
  return ADSF_id2_rightdelta;
}
