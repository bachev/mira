/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2012 and later by Bastien Chevreux
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
#include "mira/simple_2Dsignalprocessing.H"

#include "util/stlimprove.H"

#include <boost/lexical_cast.hpp>

using std::cout;
using std::cerr;
using std::endl;


template <class PI>
void mergePeaks(PI hbegin, PI hend, PI lbegin, PI lend)
{
  while(hbegin!=hend){
    if(*hbegin){
      while(hbegin != hend && *hbegin){
	++hbegin;
	++lbegin;
      }
      while(lbegin != lend && *lbegin){
	*hbegin=1;
	++hbegin;
	++lbegin;
      }
    }else{
      ++hbegin;
      ++lbegin;
    }
  }
}

// simple, dumb downslope
template <class PI, class CI>
void extendPeaks(PI pbegin, PI pend, CI cbegin, CI cend, Contig::ccctype_t minthresh, Contig::ccctype_t maxthresh)
{
  while(pbegin!=pend){
    if(*pbegin){
      while(pbegin!=pend && *pbegin) {
	++pbegin;
	++cbegin;
      }
      if(pbegin!=pend) {
	auto oldval=(cbegin-1)->total_cov;
	while(pbegin!=pend && ! *pbegin) {
	  if(cbegin->total_cov < minthresh
	     || cbegin->total_cov <= oldval) break;
	  oldval=cbegin->total_cov;
	  *pbegin=1;
	  ++pbegin;
	  ++cbegin;
	}
      }
    }else{
      ++pbegin;
      ++cbegin;
    }
  }
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   //{cout << bla; cout.flush();}
#define CEBUG(bla)
void Contig::findPeaks(ccctype_t avgcov, std::vector<uint8> & peakindicator200)
{
  peakindicator200.clear();
  if(CON_counts.empty()) return;

  std::vector<uint8> peakindicator150;
  findPeaks_helper(avgcov,avgcov*2,peakindicator200);
  findPeaks_helper(avgcov,avgcov+avgcov/2,peakindicator150);
  mergePeaks(peakindicator200.begin(),peakindicator200.end(),
	     peakindicator150.begin(),peakindicator150.end());
  CEBUG(""; dbgContainerToWiggle(peakindicator200,getContigName(),"07a_merged"));
  mergePeaks(peakindicator200.rbegin(),peakindicator200.rend(),
	     peakindicator150.rbegin(),peakindicator150.rend());
  CEBUG(""; dbgContainerToWiggle(peakindicator200,getContigName(),"07b_merged"));
}
#undef CEBUG

/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   //{cout << bla; cout.flush();}
#define CEBUG(bla)
void Contig::findPeaks_helper(ccctype_t avgcov, ccctype_t threshold, std::vector<uint8> & peakindicator)
{
  peakindicator.clear();
  peakindicator.resize(CON_counts.size(),0);

  if(CON_counts.empty()) return;

  std::string tstr=boost::lexical_cast<std::string>(threshold);
  //CEBUG(threshold << "\t" << tstr << endl);

  auto piI=peakindicator.begin();
  for(auto & cce: CON_counts){
    if(cce.total_cov>=threshold) *piI=1;
    ++piI;
  }
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"01_piraw_"+tstr));

  erode(peakindicator.begin(),peakindicator.end(),20);
  erode(peakindicator.rbegin(),peakindicator.rend(),20);
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"02_erode1_"+tstr));

  dilate(peakindicator.begin(),peakindicator.end(),40);
  dilate(peakindicator.rbegin(),peakindicator.rend(),40);
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"03_dilate_"+tstr));

  erode(peakindicator.begin(),peakindicator.end(),20);
  erode(peakindicator.rbegin(),peakindicator.rend(),20);
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"04_erode2_"+tstr));

  removeSmallPeaks(peakindicator,50);
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"05_remsmallp_"+tstr));

  // peak >= threshold found, but extend the peak area downslope to ~1x (avgcov)
  extendPeaks(peakindicator.begin(),peakindicator.end(),
	      CON_counts.begin(), CON_counts.end(),
	      avgcov, threshold);
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"06a_extendp_"+tstr));
  // hdeque (CON_counts) has no reverse iterators yet ... :-(
  // workaround (if I don't want to duplicate extendPeaks()):
  //   reverse the container
  // what a waste
  mstd::reverse(CON_counts);
  extendPeaks(peakindicator.rbegin(),peakindicator.rend(),
	      CON_counts.begin(), CON_counts.end(),
	      avgcov, threshold);
  mstd::reverse(CON_counts);
  CEBUG(""; dbgContainerToWiggle(peakindicator,getContigName(),"06b_extendp_"+tstr));

  return;
}
#undef CEBUG



/*************************************************************************
 *
 * Note: currently cannot constify as CON_counts, when using hdeque, has no
 *   "real" const_iterator
 *
 * From are GenBank compliant boundaries, i.e., they are [...] (including) and
 *   not C-like boundaries [...[ (excluding)
 *
 * countgaps: if true, gaps are counted toward coverage. If not, not.
 *   (countgaps==false atm used only by ...)
 *
 *************************************************************************/

void Contig::collectCoverage(uint32 from, uint32 to, std::vector<uint64> & covvals, bool countgaps)
{
  FUNCSTART("void Contig::collectCoverage(uint32 from, uint32 to, std::vector<uint64> & covvals)");

  covvals.clear();
  if(CON_counts.empty()) return;

  BUGIFTHROW(from>to,"from " << from << " > to " << to << " ???");
  BUGIFTHROW(from>=CON_counts.size(),"from " << from << " >= CON_counts.size() " << CON_counts.size() << " ???");
  BUGIFTHROW(to>=CON_counts.size(),"to " << to << " >= CON_counts.size() " << CON_counts.size() << " ???");
  BUGIFTHROW(CON_counts.empty() && (from!=0 || to!=0),"empty contig, but from " << from << " and to " << to);

  uint64 numvals=to - from+1;
  if(numvals>0){
    covvals.resize(numvals,0);
    auto cvI=covvals.begin();
    auto ccI=CON_counts.cbegin()+from;
    for(auto ccE=ccI+numvals; ccI!=ccE; ++ccI, ++cvI){
      if(countgaps) {
	*cvI=ccI->total_cov;
      }else if(ccI->total_cov > ccI->star) {   // should always be true, but ... oh well
	*cvI=ccI->total_cov-ccI->star;
      } // else 0, but pre-initialised
    }
  }

}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::calcStatsOnContainer(coverageinfo_t & cinfo, std::vector<uint64> & covvals) const
{
  FUNCSTART("void Contig::calcStatsOnContainer(coverageinfo_t & cinfo, std::vector<uint64> & covvals)");

  cinfo.clear();

  if(covvals.empty()) return;

  uint64 sum=0;
  for(auto cve : covvals){
    sum+=cve;
    cinfo.min=std::min(cinfo.min,static_cast<uint64>(cve));
    cinfo.max=std::max(cinfo.max,static_cast<uint64>(cve));
  }
  cinfo.sum=sum;
  cinfo.mean=static_cast<double>(sum)/static_cast<double>(covvals.size());
  // use non-parallel sort, this is more or less a small to medium thing
  /// then again ...
  if(covvals.size()<65536){
    mstd::ssort(covvals);
  }else{
    mstd::psort(covvals);
  }

  if(covvals.size()%2==0){
    cinfo.median=static_cast<double>(covvals[covvals.size()/2]+covvals[covvals.size()/2-1])/2;
  }else{
    cinfo.median=static_cast<double>(covvals[covvals.size()/2]);
  }

  double sqsum=0.0;
  for(auto cvI=covvals.begin(); cvI!=covvals.end(); ++cvI){
    sqsum+=(static_cast<double>(*cvI)-cinfo.mean)*(static_cast<double>(*cvI)-cinfo.mean);
  }

  cinfo.stddev=sqrt(static_cast<double>(sqsum)/static_cast<double>(covvals.size()));
}

/*************************************************************************
 *
 * precondition: covvals must be sorted ascending
 *
 *
 *
 *************************************************************************/

void Contig::calcSecondOrderStatsOnContainer(coverageinfo_t & tci, const std::vector<uint64> & covvals) const
{

  int64 lowerb=static_cast<int64>(tci.mean/3);
  int64 upperb=static_cast<int64>(tci.mean*3);

  auto cvS=covvals.begin();

  // advance to beginning of area
  for(; cvS != covvals.end() && static_cast<int64>(*cvS) < lowerb; ++cvS) {};

  auto cvE=cvS;

  uint64 sum=0;
  for(; cvE != covvals.end() && static_cast<int64>(*cvE) <= upperb; ++cvE){
    sum+=*cvE;
  }

  int64 numvals=cvE-cvS;

  if(numvals>0){
    tci.mean=static_cast<double>(sum)/static_cast<double>(numvals);
    if(numvals%2==0){
      tci.median=static_cast<double>(*(cvS+numvals/2)+(*(cvS+(numvals/2-1))))/2.0;
    }else{
      tci.median=static_cast<double>(*(cvS+numvals/2));
    }

    double sqsum=0.0;
    for(auto cvI=cvS; cvI!=cvE; ++cvI){
      sqsum+=(static_cast<double>(*cvI)-tci.mean)*(static_cast<double>(*cvI)-tci.mean);
    }

    tci.stddev=sqrt(static_cast<double>(sqsum)/static_cast<double>(covvals.size()));
  }

  return;
}
