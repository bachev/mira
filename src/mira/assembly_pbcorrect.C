/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2013 and later by Bastien Chevreux
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


#include "mira/assembly.H"

// BOOST
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "errorhandling/errorhandling.H"
#include "util/progressindic.H"
#include "util/dptools.H"
#include "util/fileanddisk.H"
#include "caf/caf.H"
#include "mira/align.H"

#include "mira/ads.H"
#include "mira/simple_2Dsignalprocessing.H"

using std::cout;
using std::endl;

#define CEBUG(bla)

// cs1 for normal clocking ('user compatible' as is does not disturb)
//  cs2 for extensive clocking output, more for analysis of MIRA behaviour
//  during development

#ifndef PUBLICQUIET
#define CLOCK_STEPS2
#endif
#define CLOCK_STEPS2

//ecopbownH5_31561


#define CEBUG(bla) { if(docebug) {cout << bla; cout.flush();}}


std::string readofinterest1="ecopbownH5_38629a";
std::string readofinterest2="ecopbownH5_25312a";

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint32 Assembly::correctPBMain()
{
  FUNCSTART("Assembly::correctPBMain()");

  uint32 startpass=1;

  uint32 totaleditsmade=0;

  basicDataChecks();
  ensureStandardDirectories(false);

  AS_miraparams[0].getNonConstAssemblyParams().as_dateoutput=false;
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_skimchimeradetection=false;  // needs to be off?
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_skimjunkdetection=false;  // needs to be off?

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  // fill quick lookup which sequencing types are present as well as whether backbones are there
  AS_seqtypespresent.clear();
  AS_seqtypespresent.resize(ReadGroupLib::SEQTYPE_END,false);
  for(uint32 rgi=0; rgi< ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(rgid.isBackbone()){
      AS_hasbackbones=true;
    }
    if(!rgid.isRail() && !rgid.isBackbone()){
      AS_seqtypespresent[rgid.getSequencingType()]=true;
    }
  }

  AS_permanent_overlap_bans.resize(AS_readpool.size());
  AS_used_ids.resize(AS_readpool.size(),0);
  AS_needalloverlaps.resize(AS_readpool.size(),false);   // not really used, but computeSWAlign() needs it

  if(startpass==1){
    findDegeneratePolymerases(AS_readpool);
  }

  cout << "RLEing reads\n";
  Read::setCoutType(Read::AS_TEXT);
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    AS_readpool[ri].rleRead();
    //cout << AS_readpool[ri] << endl;
  }

/*
  {
    std::string tmpfname;
    tmpfname=buildFileName(0,"","",
			  as_fixparams.as_tmpf_clippings,
			  ".txt");

    // doing it twice catches a few outliers missed the first time
    std::string logprefix="proposed cutback 1a: ";
    uint64 numclipped=performNewProposedCutbackClips(tmpfname,logprefix);
    if(numclipped>0){
      logprefix="proposed cutback 1b: ";
      performNewProposedCutbackClips(tmpfname,logprefix);
    }else{
      cout << "No bases clipped in first pec round, skipping second round.\n";
    }
    dumpSomeStatistics();
  }
*/

  totaleditsmade+=correctPBMainNormal(startpass);
  //totaleditsmade+=correctPBMainLight(startpass);

  cout << "Hard trim" << endl;
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    AS_readpool[ri].performHardTrim();
  }
  cout << "Hard trim done" << endl;

  {
    std::ofstream ostr("corrected_intermediate_f.fastq",std::ios::out);
    AS_readpool.dumpAs(ostr,Read::AS_FASTQ,true);
  }

  //colourReadsByRLE();

  cout << "derle" << endl;
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    // carefull: hard trim will have completely discarded empty reads
    if(AS_readpool[ri].getRLEValues() != nullptr){
      AS_readpool[ri].deRLERead();
    }
  }
  cout << "derle done" << endl;

  {
    std::ofstream ostr("corrected.maf",std::ios::out);
    AS_readpool.dumpAs(ostr,Read::AS_MAF,true);
  }
  {
    std::ofstream ostr("corrected.fastq",std::ios::out);
    AS_readpool.dumpAs(ostr,Read::AS_FASTQ,true);
  }

  AS_naclipl.clear();
  AS_naclipr.clear();

  return totaleditsmade;
}

//ecopbownH5_31561



uint32 Assembly::correctPBMainLight(uint32 startpass)
{
  FUNCSTART("Assembly::correctPBMainLight()");
  uint32 totaleditsmade=0;

  pbc_timing_t timing;

  // pass 1
  uint32 actpass=startpass;
  int32 additionalbelieveborder=0;
  // later passes: already larger kmer and borders, but still problems at microrepeat sites
  //  so, in later passes ... very rare things: don't believe
  uint8 minbk=0;
  bool tweakbkmar=false;
  bool killnonhafstretches=false;
  bool generatenaclips=false;
  bool generaterleedits=false;
  bool generatebaseedits=true;

  float ratioacceptvalue=0.0072; // ratio < .72 % not accepted

  // using mkfr=4 (instead of 1) here considerably increases the quality of the resulting sequences
  //  (assemblies have less homopolymer problems) at the expense of worse resolving of low coverage
  //  regions, hence more data clipped
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_pec_mkfr=4;

  if(actpass==1){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

  // pass 2
  ratioacceptvalue=0.013;
  AS_miraparams[0].getNonConstSkimParams().sk_basesperhash=29;
  if(actpass==2){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

  timeval tv;
  gettimeofday(&tv,nullptr);

  priv_performHashAnalysis("",false,actpass==1, actpass, "", "_pass");

  timing.hashanalysis+=diffsuseconds(tv);

  gettimeofday(&tv,nullptr);
  killNonHAFCoveredStretchesInReads();
  timing.killnonhaf+=diffsuseconds(tv);

  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    auto & actread=AS_readpool[ri];
    if(actread.getLenClippedSeq()>0){
      auto cI=actread.getClippedSeqAsChar();
      uint32 numn=0;
      for(uint32 csl=0; csl<actread.getLenClippedSeq(); ++csl, ++cI){
	if(toupper(*cI)=='N') ++numn;
      }
      float maskratio=100.0/actread.getLenClippedSeq()*numn;
      cout << "LPerc: " << maskratio << " " << actread.getName() << endl;
    }else{
      cout << "LPerc: out " << actread.getName() << endl;
    }
  }

  cout << "Saving final results ..." << endl;
  // TODO: gcc version 4.6.1 (GCC)
  // oooooops: src/tcmalloc.cc:390] Attempt to free invalid pointer: 0x7a54a0
  //std::ofstream ostr("corrected_intermediate_"+boost::lexical_cast<std::string>(actpass)+".fastq",std::ios::out);
  std::string fname("corrected_intermediate_f.fastq");
  std::ofstream ostr(fname,std::ios::out);
  Read::setCoutType(Read::AS_FASTQ);
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    auto tmpread=AS_readpool[ri];
    tmpread.performHardTrim();
    if(tmpread.getRLEValues() != nullptr){
      tmpread.deRLERead();
    }
    ostr << tmpread;
  }



  return totaleditsmade;
}


uint32 Assembly::correctPBMainNormal(uint32 startpass)
{
  FUNCSTART("Assembly::correctPBMainNormal()");

  uint32 totaleditsmade=0;

  pbc_timing_t timing;

  // pass 1
  uint32 actpass=startpass;
  int32 additionalbelieveborder=0;
  // later passes: already larger kmer and borders, but still problems at microrepeat sites
  //  so, in later passes ... very rare things: don't believe
  uint8 minbk=0;
  bool tweakbkmar=false;
  bool killnonhafstretches=false;
  bool generatenaclips=false;
  bool generaterleedits=false;
  bool generatebaseedits=true;

  float ratioacceptvalue=0.0072; // ratio < .72 % not accepted

  // using mkfr=4 (instead of 1) here considerably increases the quality of the resulting sequences
  //  (assemblies have less homopolymer problems) at the expense of worse resolving of low coverage
  //  regions, hence more data clipped
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_pec_mkfr=4;

  if(actpass==1){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

  // pass 2
  ratioacceptvalue=0.013;
  if(actpass==2){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

  // mkfr=1 at this place counterbalances mkfr=4 above (if set to that high) by resolving again a little bit
  //  better areas with less coverage. Less data clipped. But seems to have more indels even in good coverage
  //AS_miraparams[0].getNonConstAssemblyParams().as_clip_pec_mkfr=1;

  // pass 3
  ratioacceptvalue=0.033;
  if(actpass==3){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }


  // optional pass 4
  // see above, to get better resolution in low cov areas
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_pec_mkfr=1;
  if(actpass==4){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

  // ?
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_pec_mkfr=4;


  // pass 5:
  //  - replace with N (which will hopefully be edited away lateron)
  //  - non-align clips
  ratioacceptvalue=0.036;
  AS_miraparams[0].getNonConstSkimParams().sk_basesperhash=29;
  additionalbelieveborder=6;
  killnonhafstretches=true;
  generatenaclips=true;

  if(actpass==5){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

  // ?
  AS_miraparams[0].getNonConstAssemblyParams().as_clip_pec_mkfr=1;

  // pass 6 polishing
  ratioacceptvalue=0.0625;
  minbk=4;
  tweakbkmar=true;
  killnonhafstretches=false;
  //generatenaclips=false; // ?
  generatenaclips=false;
  if(actpass==6){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

//*
  // pass 7: RLE edits
  ratioacceptvalue=0.073;
  generaterleedits=true;
  generatebaseedits=false;
  additionalbelieveborder=0;
  minbk=0;
  if(actpass==7){
    totaleditsmade+=cpbm_passHelper(actpass,ratioacceptvalue,additionalbelieveborder,minbk,tweakbkmar,killnonhafstretches,
				    generatenaclips,generaterleedits,generatebaseedits,timing);
    ++actpass;
  }

//*/

  cout << "Total edits made: " << totaleditsmade << '\n';

  cout << timing;

  return totaleditsmade;
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

uint32 Assembly::cpbm_passHelper(uint32 actpass, float ratioacceptvalue, int32 additionalbelieveborder, uint8 minbk, bool tweakbkmar, bool killnonhafstretches, bool generatenaclips, bool generaterleedits, bool generatebaseedits, pbc_timing_t & timing)
{
  cout << "\n\nPB correct pass: " << actpass << " / " << AS_miraparams[0].getAssemblyParams().as_numpasses << '\n';

  if(generatenaclips){
    AS_naclipl.clear();
    AS_naclipl.resize(AS_readpool.size(),-1);
    AS_naclipr.clear();
    AS_naclipr.resize(AS_readpool.size(),-1);
  }

  timeval tv;
  gettimeofday(&tv,nullptr);

  priv_performHashAnalysis("",false,actpass==1, actpass, "", "_pass");

  timing.hashanalysis+=diffsuseconds(tv);

  if(killnonhafstretches) {
    gettimeofday(&tv,nullptr);
    killNonHAFCoveredStretchesInReads();
    timing.killnonhaf+=diffsuseconds(tv);

    for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
      auto & actread=AS_readpool[ri];
      if(actread.getLenClippedSeq()>0){
	auto cI=actread.getClippedSeqAsChar();
	uint32 numn=0;
	for(uint32 csl=0; csl<actread.getLenClippedSeq(); ++csl, ++cI){
	  if(toupper(*cI)=='N') ++numn;
	}
	float maskratio=100.0/actread.getLenClippedSeq()*numn;
	cout << "KPerc: " << maskratio << " " << actread.getName() << endl;
      }else{
	cout << "KPerc: out " << actread.getName() << endl;
      }
    }

  }

  gettimeofday(&tv,nullptr);
  findPossibleOverlaps(actpass, "", "_pass");
  timing.possibleoverlaps+=diffsuseconds(tv);

  gettimeofday(&tv,nullptr);
  auto editsmade=tryPBCorrect(actpass,
			      ratioacceptvalue,
			      additionalbelieveborder,minbk,tweakbkmar,
			      generatenaclips,generaterleedits,generatebaseedits,
			      timing);
  timing.trypbcorrect+=diffsuseconds(tv);

  cout << "Edits made in pass " << actpass << ": " << editsmade << '\n';

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);


  // apply non-align clips
  // IMPORTANT: do this *before* somehow changing the length of the sequences,
  //            e.g. before removing gap bases (below)
  if(generatenaclips){
    for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
      if(AS_naclipl[ri]>0){
	AS_readpool[ri].setLQClipoff(AS_naclipl[ri]);
      }
      if(AS_naclipr[ri]>0){
	AS_readpool[ri].setRQClipoff(AS_naclipr[ri]);
      }
    }
    AS_naclipl.clear();
    AS_naclipr.clear();
  }

  cout << "Removing gap bases ..."; cout.flush();
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    AS_readpool[ri].removeGapsFromRead();
  }
  cout << "done.\n";

  cout << "Saving intermediate results ..." << endl;
  // TODO: gcc version 4.6.1 (GCC)
  // oooooops: src/tcmalloc.cc:390] Attempt to free invalid pointer: 0x7a54a0
  //std::ofstream ostr("corrected_intermediate_"+boost::lexical_cast<std::string>(actpass)+".fastq",std::ios::out);
  std::string fname("corrected_intermediate_");
  fname+=boost::lexical_cast<std::string>(actpass);
  fname+=".fastq";
  std::ofstream ostr(fname,std::ios::out);
  Read::setCoutType(Read::AS_FASTQ);
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    auto tmpread=AS_readpool[ri];
    tmpread.performHardTrim();
    if(tmpread.getRLEValues() != nullptr){
      tmpread.deRLERead();
    }
    ostr << tmpread;
  }

  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    auto & actread=AS_readpool[ri];
    if(actread.getLenClippedSeq()>0){
      auto cI=actread.getClippedSeqAsChar();
      uint32 numn=0;
      for(uint32 csl=0; csl<actread.getLenClippedSeq(); ++csl, ++cI){
	if(toupper(*cI)=='N') ++numn;
      }
      float maskratio=100.0/actread.getLenClippedSeq()*numn;
      cout << "NPerc: " << maskratio << " " << actread.getName();
      if(maskratio>25.0){
	cout << " complete kill";
	//for(uint32 hcpos=0; hcpos<hafcov.size(); ++hcpos){
	//	actread.changeBaseInClippedSequence('n',0,hcpos);
	//}
	actread.setRQClipoff(actread.getLeftClipoff());
      }
      cout << endl;
     }else{
      cout << "NPerc: out " << actread.getName() << endl;
    }
  }

  return editsmade;
}

/*************************************************************************
 *
 * Degenerate polymerases lead to bases being sequenced multiple times
 *  instead of the correct number of times, leading to sequences like this
 *  xxxxxxxxxxxAAAAAAAACCCGGGGGGGGTTTTTTTTCCGTTTTTTTTAAATTGGGGGGGGGGGG...
 *
 * While in RLE space the bases may match quite well to other reads, they
 *  obviously have totally wrong base counts which could make RLE correction
 *  less accurate. Therefore, simply cut away the degenerate part of a read
 *
 * Method: Count number of base changes of a non-RLE(!) read in a given
 *  window (I'm using 96 bases). Then a simple cutoff at 48 to create a
 *  basic map of what is 'good', then a bit of 2D signal signal processing,
 *  then look for large 'non-good' areas, done.
 *
 *************************************************************************/

void Assembly::findDegeneratePolymerases(ReadPool & rp)
{
  FUNCSTART("void Assembly::findDegeneratePolymerases()");

  //bool docebug=true;
  bool docebug=false;

  std::vector<uint32> bchanges;
  bchanges.reserve(50000);

  std::vector<uint8> goodmap;
  goodmap.reserve(50000);

  uint32 bandsize=48;

  cout << "Looking for degenerate polymerases:\n";

  ProgressIndicator<int64> P(0,rp.size());
  uint32 numdegen=0;
  for(uint32 ri=0; ri<rp.size(); ++ri){
    P.increaseprogress();
    auto & actread=rp[ri];

    CEBUG("FDP " << actread.getName() << endl);
    if(actread.getLenClippedSeq()<bandsize+10) continue;

    bchanges.clear();
    bchanges.resize(actread.getLenClippedSeq(),0);
    goodmap.clear();
    goodmap.resize(actread.getLenClippedSeq(),255);


    auto bI=bchanges.begin()+bandsize;
    auto * seq = actread.getClippedSeqAsChar()+bandsize;
    auto * send = seq + actread.getLenClippedSeq()-bandsize;
    for(;seq < send; ++seq, ++bI){
      auto * wseq = seq-bandsize;
      auto * eseq = seq+bandsize;
      char currbase=0;
      char numchanges=0;
      for(; wseq<eseq; ++wseq){
	if(toupper(*wseq)!=currbase) ++numchanges;
	currbase=toupper(*wseq);
      }
      *bI=numchanges;
    }

    bI=bchanges.begin();
    auto gI=goodmap.begin();
    for(; bI!=bchanges.end(); ++bI, ++gI){
      if(*bI<48) *gI=0;
    }
    // close small holes len 40
    dilate(goodmap.begin(),goodmap.end(),20);
    dilate(goodmap.rbegin(),goodmap.rend(),20);
    // undo the above and additionally eliminate small peaks len 40
    erode(goodmap.begin(),goodmap.end(),20+20);
    erode(goodmap.rbegin(),goodmap.rend(),20+20);
    // back to start
    dilate(goodmap.begin(),goodmap.end(),20);
    dilate(goodmap.rbegin(),goodmap.rend(),20);

    bI=bchanges.begin();
    gI=goodmap.begin();
    for(; bI!=bchanges.end(); ++bI, ++gI){
      CEBUG("fdpi: " << bI-bchanges.begin() << "\t" << *bI << "\t" << static_cast<uint16>(*gI) << endl);
    }

    gI=goodmap.begin()+bandsize;
    while(gI!=goodmap.end()){
      if(*gI==0){
	auto hI=gI;
	for(; hI!=goodmap.end() && *hI==0; ++hI){}
	if(hI-gI>100){
	  CEBUG(actread.getName() << " found degen polymerase: " << gI-goodmap.begin() << endl);
	  actread.setRQClipoff(gI-goodmap.begin());
	  ++numdegen;
	  break;
	}else{
	  gI=hI;
	}
      }else{
	++gI;
      }
    }
  }

  P.finishAtOnce();

  cout << "\nClipped " << numdegen << " reads with degenerated polymerases.\n";

  CEBUG("FDPOUT" << endl);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::tpbc_checkRead(Read & actread)
{
  FUNCSTART("void Assembly::tpbc_checkRead(Read & actread)");

  auto * ptr= actread.getClippedComplementSeqAsChar();
  for(uint32 ri=0; ri<actread.getLenClippedSeq(); ++ri, ++ptr){
    if(!dptools::isValidIUPACStarBase(*ptr)){
      Read::setCoutType(Read::AS_MAF);
      cout << actread;
      BUGIFTHROW(true,"Ouchi X for " << actread.getName() << " with " << *ptr << " (" << static_cast<uint16>(*ptr) << ")");
    }
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
uint32 Assembly::tryPBCorrect(uint32 actpass, float ratioacceptvalue, int32 additionalbelieveborder, uint8 minbk, bool tweakbkmar, bool generatenaclips, bool generaterleedits, bool generatebaseedits, pbc_timing_t & timing)
{
  FUNCSTART("uint32 Assembly::tryPBCorrect(uint32 actpass, int32 additionalbelieveborder)");

  bool docebug=false;

  uint32 editsmade=0;

  timeval tv;
  gettimeofday(&tv,nullptr);

  uint64 tsizeneeded=getFileSize(AS_posfmatch_filename)+getFileSize(AS_poscmatch_filename);
  tsizeneeded/=sizeof(skimhitforsave_t);
  tsizeneeded*=2;

  std::vector<skimhitforsave_t> vsh;
  vsh.reserve(tsizeneeded);


//*
  loadVector(vsh,AS_posfmatch_filename,0);

  // denormalise this part (forward / forward)
  // also set fwd / rev directions (not done atm by skim)
  {
    auto vI=vsh.begin();
    for(auto nume=vsh.size(); nume>0; --nume, ++vI){
      vI->setRID1Dir(1);
      vI->setRID2Dir(1);
      vsh.push_back(*vI);
      std::swap(vsh.back().rid1,vsh.back().rid2);
      vsh.back().eoffset=-vsh.back().eoffset;
    }
  }
//*/

//*
  auto oldsize=vsh.size();
  loadVector(vsh,AS_poscmatch_filename,0);

  // denormalise this part (forward / reverse)
  // some bits of magic here:
  //   r1    ------>
  //   r2        <-----
  // first becomes
  //   r2        <-----
  //   r1    ------>
  // and this is conceptually identical to
  //   r2    ------>
  //   r1        <-----
  // The added bonus is that the align object will always get the "interesting" read
  //  in the same, forward direction

  ADSEstimator adse;
  auto vI=vsh.begin()+oldsize;
  for(auto nume=vsh.size()-oldsize; nume>0; --nume, ++vI){
    adse.calcNewEstimateFromSkim(
      vI->eoffset,
      AS_readpool[vI->rid1].getLenClippedSeq(),
      AS_readpool[vI->rid2].getLenClippedSeq(),
      vI->rid1,
      vI->rid2,
      1,
      -1);

    vI->setRID1Dir(1);
    vI->setRID2Dir(-1);
    vsh.push_back(*vI);

    vsh.back().eoffset=adse.getEstimatedRightExpand(vsh.back().rid1);
    std::swap(vsh.back().rid1,vsh.back().rid2);
  }

  //BUGIFTHROW(vsh.size() != vsh.capacity(),"vsh.size() != vsh.capacity() ???");

//*/

  // TODO: maybe sort differently: have longest reads sorted to front? Hmmmm.
  sort(vsh.begin(),vsh.end(), skimhitforsave_t::stdSortCmp);

  cout << "VSH:\n";
  for(auto & vshe : vsh){
    CEBUG(vshe.rid1 << '\t' << vshe.rid2 << '\n');
  }

  for(auto & s2se : AS_s2saligninfo){
    s2se.second.alreadyaligned=false;
  }

  auto lmax=AS_readpool[0].getLenClippedSeq();
  for(auto ri=1; ri<AS_readpool.size(); ++ri){
    lmax=std::max(lmax,AS_readpool[ri].getLenClippedSeq());
  }

  std::vector<pbcounts_t> correctorcounts;
  std::vector<std::vector<uint32> > rlecorrectorcounts;
  correctorcounts.reserve(lmax);
  rlecorrectorcounts.reserve(lmax);

  std::vector<Align> chkalign;
  setupAlignCache(chkalign);

  std::list<AlignedDualSeq> madsl;

  // these two taken out of tpbc_fillCorrector() as also used in new RLE corrector routines
  std::vector<uint8> bkmar;  // believe kmer actrid
  std::vector<uint8> bkmor;  // believe kmer otherrid
  int32 kmersizeused=AS_miraparams[0].getNonConstSkimParams().sk_basesperhash;

  timing.tpbc_setup+=diffsuseconds(tv);


  cout << "Correcting " << AS_readpool.size() << " reads via " << vsh.size() << " potential matches:\n";
  ProgressIndicator<int64> P(0,vsh.size());
  readid_t actrid=-1;
  for(auto vI=vsh.cbegin(); vI!=vsh.cend(); ++vI){
    P.increaseprogress();
    if(vI->rid1 != actrid){
      if(actrid != -1){
	if(generatenaclips) tpbc_generateNonAlignClips(actrid,correctorcounts);
	if(generaterleedits) tpbc_generateRLEEdits(actpass,actrid,rlecorrectorcounts);
	if(generatebaseedits) editsmade+=tpbc_generateBaseEdits(actpass,actrid,correctorcounts);
	tpbc_checkRead(AS_readpool[actrid]);
      }
      actrid=vI->rid1;
      correctorcounts.clear();
      correctorcounts.resize(AS_readpool[actrid].getLenClippedSeq());
      rlecorrectorcounts.clear();
      rlecorrectorcounts.resize(AS_readpool[actrid].getLenClippedSeq());
    }
    BUGIFTHROW(vI->getRID1Dir()<0,"vI->getRID1Dir()<0 ???");

    int32 hintbandwidth=-1;
    uint64 uomapkey=vI->rid1;
    uomapkey<<=32;
    uomapkey+=vI->rid2;

    int32 safetydist=35;
    auto s2sI=AS_s2saligninfo.find(uomapkey);
    if(s2sI!=AS_s2saligninfo.end()){
      if(s2sI->second.minbanddistance>=safetydist){
	hintbandwidth=s2sI->second.bandwidthused/2 - (s2sI->second.minbanddistance-safetydist);
      }
    }

    adse.calcNewEstimateFromSkim(
      vI->eoffset,
      AS_readpool[vI->rid1].getLenClippedSeq(),
      AS_readpool[vI->rid2].getLenClippedSeq(),
      vI->rid1,
      vI->rid2,
      vI->getRID1Dir(),
      vI->getRID2Dir());

    auto estimovl=adse.getEstimatedOverlap();
    CEBUG("RL pos: " << vI-vsh.begin() << endl);
    CEBUG("astats bsw " << AS_readpool[vI->rid1].getName() << "\t" << AS_readpool[vI->rid2].getName() << '\t');

    if(actpass>1 && hintbandwidth<0) hintbandwidth=200;
    if(hintbandwidth<0) hintbandwidth=400;

    if(hintbandwidth>=0) CEBUG("hint ");
    CEBUG(hintbandwidth << '\t' << estimovl << '\t' << *vI);
    if(s2sI != AS_s2saligninfo.end() && s2sI->second.alreadyaligned){
      CEBUG("astats already aligned\n");
    }else{
      gettimeofday(&tv,nullptr);
      computeSWAlign(madsl,vI->rid1,vI->rid2,vI->eoffset,vI->getRID1Dir()*vI->getRID2Dir(),chkalign,hintbandwidth);
      timing.tpbc_sw+=diffsuseconds(tv);

      if(!madsl.empty()){
	if(s2sI!=AS_s2saligninfo.end()){
	  s2sI->second.bandwidthused=madsl.front().getBandwidthUsed();
	  s2sI->second.minbanddistance=madsl.front().getMinBandDistance();
	}else{
	  // gaaaaaaaah, insert() does not give back the correct iterator, need to ...
	  AS_s2saligninfo.insert(
	    std::pair<uint64,s2saligninfo_t>(uomapkey,
					     s2saligninfo_t(madsl.front().getBandwidthUsed(),
							    madsl.front().getMinBandDistance())));
	  // ... find it again
	  s2sI=AS_s2saligninfo.find(uomapkey);
	}
	s2sI->second.alreadyaligned=true;
	CEBUG("astats aar " << vI->rid1 << " " << vI->rid2 << "\tbwu: " << madsl.front().getBandwidthUsed()
	      << "\tmbd: " << madsl.front().getMinBandDistance()
	      << endl);
      }else{
	// ooooooops? what to do? Let's double the hint via bandwidth used
	// and set the minbanddistance to a low value so that in the next pass
	//  we maybe get again an alignment
	CEBUG("astats missed!\n");
	s2sI->second.bandwidthused*=2;
	s2sI->second.minbanddistance=5;
      }

      if(!madsl.empty()){
	bkmar.clear();
	bkmar.resize(madsl.front().getOverlapLen(),0);
	tpbc_fc_makeBelieveKMERMap(actrid,vI->rid2, madsl.front(), bkmar, kmersizeused, additionalbelieveborder);
	if(minbk>1) tpbc_fc_minimumBelieveKMerMap(bkmar,minbk);

	bkmor.clear();
	bkmor.resize(madsl.front().getOverlapLen(),0);
	tpbc_fc_makeBelieveKMERMap(vI->rid2,actrid, madsl.front(), bkmor, kmersizeused, additionalbelieveborder);
	if(minbk>1) tpbc_fc_minimumBelieveKMerMap(bkmor,minbk);

	if(tweakbkmar) tpbc_tweakBKMAR(actrid,vI->rid2,madsl.front(),bkmar,bkmor);

	tpbc_fillCorrector(actrid,vI->rid2,madsl.front(),correctorcounts,bkmar,bkmor);
	tpbc_fillRLECorrector(actrid,vI->rid2,madsl.front(),rlecorrectorcounts,bkmar,bkmor);
      }
    }
  }
  // catch edits of last read
  if(actrid != -1){
    if(generatenaclips) tpbc_generateNonAlignClips(actrid,correctorcounts);
    if(generaterleedits) tpbc_generateRLEEdits(actpass,actrid,rlecorrectorcounts);
    if(generatebaseedits) editsmade+=tpbc_generateBaseEdits(actpass,actrid,correctorcounts);
    tpbc_checkRead(AS_readpool[actrid]);
  }

  P.finishAtOnce();
  cout << endl;

  return editsmade;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
void Assembly::tpbc_fillRLECorrector(readid_t actrid, readid_t otherrid, AlignedDualSeq & ads, std::vector<std::vector<uint32> > & rlecc,std::vector<uint8> & bkmar,std::vector<uint8> & bkmor)
{
  FUNCSTART("void Assembly::tpbc_fillRLECorrector(readid_t actrid, readid_t otherrid, AlignedDualSeq & ads, std::vector<std::vector<uint32> > & rlecc,std::vector<uint8> & bkmar,std::vector<uint8> & bkmor)");

  BUGIFTHROW(AS_readpool[otherrid].getRLEValues()==nullptr,"AS_readpool[otherrid].getRLEValues()==nullptr");

  bool docebug=false;
  if(false
     || AS_readpool[actrid].getName()==readofinterest1
     || AS_readpool[actrid].getName()==readofinterest2) docebug=true;
  if(false
     || AS_readpool[otherrid].getName()==readofinterest1
     || AS_readpool[otherrid].getName()==readofinterest2) docebug=true;

  CEBUG("fillrle " << AS_readpool[actrid].getName() << " " << AS_readpool[otherrid].getName() << endl);

  auto & oread=AS_readpool[otherrid];

  auto & orlev= *(oread.getRLEValues());

  //if(AS_readpool[otherrid].getName()!="ecopbownH5_15338"){
  //  CEBUG("SHOWME RLE:\n");
  //  for(uint32 pi=0; pi<orlev.size(); ++pi){
  //    cout << "pi: " << pi << "\tv: " << static_cast<uint16>(orlev[pi]) << endl;
  //  }
  //}

  bool waaaaahseen=false;

  auto posor=oread.getLeftClipoff();
  int32 deltaposor=1;
  if(ads.getSequenceDirection(otherrid)<0){
    posor=oread.getRightClipoff()-1;
    deltaposor=-1;
  }
  posor+=deltaposor*ads.getOffsetInAlignment(actrid);
  auto aseq = ads.getSequenceAtOverlapStart(actrid);
  auto oseq = ads.getSequenceAtOverlapStart(otherrid);
  auto rccI = rlecc.begin() + ads.getOffsetInAlignment(otherrid);
  for(uint32 ovli=0; ovli<ads.getOverlapLen(); ++ovli, ++aseq, ++oseq){
    BUGIFTHROW(rccI==rlecc.end(),"rccI==rlecc.end() ???");
    BUGIFTHROW(posor<0,"posor<0");
    BUGIFTHROW(posor>=oread.getLenSeq(),"posor " << posor << " >= oread.getLenSeq() " << oread.getLenSeq());
    CEBUG("rp: " << rccI-rlecc.begin() << "\tarp: " << AS_readpool[actrid].getAdjustmentPosOfReadPos(rccI-rlecc.begin()) << "\tovli: " << ovli << "\tposor: " << posor << '\t' << *aseq << '\t' << *oseq << '\t');

    char check='/';
    if(deltaposor>0){
      check=toupper(oread.getBaseInSequence(posor));
    }else{
      check=toupper(dptools::getComplementIUPACBase(oread.getBaseInSequence(posor)));
    }
    CEBUG(check);
    if(*oseq!='*' && *oseq!='#'){
      if(*oseq==check){
	CEBUG(" OK");
      }else{
	waaaaahseen=true;
	CEBUG(" waaaaaah!");
      }
    }


    if(*aseq!=*oseq){
      CEBUG("\td");
    }else{
      CEBUG("\t ");
    }


    CEBUG('\t' << static_cast<int16>(bkmar[ovli]));
    CEBUG('\t' << static_cast<int16>(bkmor[ovli]));
    if(*aseq!=*oseq
       && bkmar[ovli]>0
       && bkmor[ovli]>0){
      CEBUG("\touch");
    }

    CEBUG("\t" << rlecc[rccI-rlecc.begin()].size());

    if(*aseq==*oseq
       && dptools::isValidACGTBase(*oseq)
       && bkmar[ovli]>1
       && bkmor[ovli]>1){
      // check case and qual! no lower case, no 0 qual as these would be from already corrected/inserted bases

      if(isupper(oread.getBaseInSequence(posor))
	 && (oread.getQualities())[posor]>0){
	 // adjustment vector must match!!!   {pos-1, pos, pos+1}
	 //       if deletion occurred next to it, don't take. E.g. {pos-1, pos, pos+2}
	 //
	 //       Rationale: e.g.     xATAAAAx   ->  xATAx  (RLE: 114) ->   xAx    (RLE: 1)     -> ouch!
	 //                           xAAAAAAx   ->  xA**x  (RLE: 5--) ->   xAx    (RLE: 5)     -> OK
	 //
	 //       TODO: check what to do on inserted bases: E.g. {pos-1, pos, -1}

	// so, first check whether we're not at the border of the other read ...
	if(posor>0 && posor<oread.getLenSeq()-1){
	  // and now for the adjacent adjustment values
	  if(oread.getAdjustmentPosOfReadPos(posor)-1 == oread.getAdjustmentPosOfReadPos(posor-1)
	     && oread.getAdjustmentPosOfReadPos(posor)+1 == oread.getAdjustmentPosOfReadPos(posor+1)){
	    rlecc[rccI-rlecc.begin()].push_back(orlev[posor]);
	    CEBUG(" plus1 " << static_cast<uint16>(orlev[posor]));
	    if(orlev[posor]>1) {
	      CEBUG(" yeah");
	    }
	  }
	}
      }
    }

    CEBUG(endl);

    if(*aseq != '*') ++rccI;               // careful: not '#' as this is an gap edited false base of this round
    if(*oseq != '*') posor+=deltaposor;    // careful: not '#' as this is an gap edited false base of this round
  }

  BUGIFTHROW(waaaaahseen,"waaaaahseen");
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
void Assembly::tpbc_fillCorrector(readid_t actrid, readid_t otherrid, AlignedDualSeq & ads, std::vector<pbcounts_t> & correctorcounts,std::vector<uint8> & bkmar,std::vector<uint8> & bkmor)
{
  FUNCSTART("void Assembly::tpbc_fillCorrector(readid_t actrid, AlignedDualSeq & ads, std::vector<pbcounts_t> & correctorcounts)");
  BUGIFTHROW(ads.getSequenceDirection(actrid)<0,"ads.getSequenceDirection(actrid)<0");

  bool docebug=false;
  if(false
     || AS_readpool[actrid].getName()==readofinterest1
     || AS_readpool[actrid].getName()==readofinterest2) docebug=true;
  if(false
     || AS_readpool[otherrid].getName()==readofinterest1
     || AS_readpool[otherrid].getName()==readofinterest2) docebug=true;

  CEBUG("belief " << AS_readpool[actrid].getName() << " " << AS_readpool[otherrid].getName() << '\n');
  auto aseq = ads.getSequenceAtOverlapStart(actrid);
  auto oseq = ads.getSequenceAtOverlapStart(otherrid);
  auto ccI = correctorcounts.begin() + ads.getOffsetInAlignment(otherrid);
  auto ccssI = ccI;  // cc safe start
  for(uint32 ovli=0; ovli<ads.getOverlapLen(); ++ovli, ++aseq, ++oseq){
    // beware: no 'continue' in this loop! because of  "if(*aseq != '*') ++ccI;" at end

    BUGIFTHROW(ccI==correctorcounts.end(),"ccI==correctorcounts.end() ???");
    CEBUG("rp: " << ccI-correctorcounts.begin() << "\tarp: " << AS_readpool[actrid].getAdjustmentPosOfReadPos(ccI-correctorcounts.begin()) << "\tovli: " << ovli << '\t' << *aseq << '\t' << *oseq << '\t');
    if(*aseq!=*oseq){
      CEBUG('d');
    }else{
      CEBUG(' ');
    }

    CEBUG('\t' << static_cast<int16>(bkmar[ovli]));
    CEBUG('\t' << static_cast<int16>(bkmor[ovli]));
    if(*aseq!=*oseq
       && bkmar[ovli]>0
       && bkmor[ovli]>0){
      CEBUG("\touch");
    }

    if(*aseq!=*oseq || ovli+1==ads.getOverlapLen()){
      // yeah, well, the ovli+1 above will be, strictly speaking, wrong if the last characters
      //  of an alignment did not match. Then again: does not really matter for the purpose
      //  whether the length of the match was n or n-1
      // but at leas we save ourselves a post-loop check
      if(ccI-ccssI>=40){
	for(; ccssI!=ccI; ++ccssI) ++ccssI->cqsafe;
      }
      ccssI=ccI+1;
    }


    bool cancorrect=true;
    bool isbasechange=false;
    bool isinsert=false;
    if(*aseq!=*oseq
       && bkmar[ovli]==0
       && bkmor[ovli]>0){
      // make sure that if we are in an indel situation, that the bkmor are set for the whole
      //  run +/-1 position
      if(dptools::isValidACGTNBase(*aseq) && dptools::isValidACGTBase(*oseq)){
	isbasechange=true;
      }else if(*oseq=='*' && dptools::isValidACGTNBase(*aseq)){
	cancorrect=false;
	auto * checkptr=oseq;
	uint32 dx=0;
	// move to begin of run
	while(ovli-dx>=0 && (*checkptr=='*' || *checkptr==*aseq)) {--dx; --checkptr;}
	// check front set?
	if(ovli-dx>=0 && bkmor[ovli-dx]){
	  checkptr=oseq;
	  dx=0;
	  // move to end of run
	  while(ovli+dx<ads.getOverlapLen() && (*checkptr=='*' || *checkptr==*aseq)) {++dx; ++checkptr;}
	  if(ovli+dx<ads.getOverlapLen() && bkmor[ovli+dx]){
	    cancorrect=true;
	  }
	}
	if(cancorrect){
	  isbasechange=true;
	}else{
	  CEBUG("\trunnocorrect");
	}
      }else if(ovli>0 && *aseq=='*' && *(aseq-1)!='*' && dptools::isValidACGTNBase(*oseq)){
	// above ... ovli >0 means: even if we had a gap at the first place (which will hopefully never
	//  happen anyway, do not try to correct
	// above ... *(aseq-1)!='*' means: if we already were in a run, do not try to correct (first gap
	//  of run will have done it's job already. Can use aseq-1 directly because of ovli>0, so we're
	//  not at the beginning of the sequence anyway

	// TODO: check bkmor in acteq gap-run for full belief (bit like for dx check above)

	isinsert=true;
	std::string potential_insert;
	auto * iaptr=aseq;
	auto * ioptr=oseq;
	for(; *iaptr=='*'; ++iaptr, ++ioptr){
	  // careful here: earlier reads may have been corrected toward a gap, in the Align sequences
	  //  those show up as '#' and we do not want these!
	  if(*ioptr!='#') potential_insert+=*ioptr;
	}
	// should always be true, but you never know ...
	if(!potential_insert.empty()){
	  bool foundinsert=false;
	  for(auto & ie : ccI->cxia){
	    if(ie.what==potential_insert){
	      foundinsert=true;
	      ++ie.count;
	      break;
	    }
	  }
	  if(!foundinsert){
	    ccI->cxia.push_back(pbcounts_t::insertafter_t(potential_insert));
	  }
	}
      }

      if(cancorrect){
	if(isbasechange){
	  CEBUG("\tcancorrect_bc");
	  if(*aseq != '*'){
	    switch(*oseq){
	    case 'A' :
	      ccI->cxa+=1;
	      break;
	    case 'C' :
	      ccI->cxc+=1;
	      break;
	    case 'G' :
	      ccI->cxg+=1;
	      break;
	    case 'T' :
	      ccI->cxt+=1;
	      break;
	    case '*' :
	      ccI->cxgap+=1;
	      break;
	    case '#' :
	      // 'oldgap' marker from AlignedDualSeq ... well, rather newly inserted gap in already edited read
	      // intentionally do nothing atm
	      break;
	    default :
	      // bug if not IUPAC (for IUPAC we don't do anything)
	      BUGIFTHROW(!dptools::isValidIUPACBase(*oseq),"Uh oh ... " << *oseq << "(" << static_cast<uint16>(*oseq) << ") is not a IUPAC base. At this stage???\n" << ads);
	    }
	  }
	}else if(isinsert){
	  CEBUG("\tcancorrect_is");
	}
      }
    }

    CEBUG(endl);

    if(*aseq != '*') ++ccI;  // careful: not '#' as this is an gap edited false base of this round
  }
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * ecopbownH5_38193, around RLE pos 280 (nonRLE 390)
 *
 * CTGCTCACGCATACGCACG   reality
 * CTGCTC  GCATA GCACG   read
 *
 * CTGCT*  *CA*AcGCACG   after pass 1
 *
 * Gaaah ... the deleted bases are from "spurious" misalignments with other reads, reaching
 * the threshold level of 2 minedits to support the edit. This continues to:
 *
 * CTGCTCACGCATACGCACG   reality
 * CTGCTC******AcGCACG   pass 5
 *       rrrrr rrrrr     <--- microrepeat
 *
 * The microrepeat ACGAtACGAc then creates the belief map (against a 'true' read)
 *
 * rp: 233 arp: 226        ovli: 233       C       C               5       14
 * rp: 234 arp: 227        ovli: 234       T       T               5       14
 * rp: 235 arp: 228        ovli: 235       G       G               5       14
 * rp: 236 arp: 229        ovli: 236       C       C               5       14
 * rp: 237 arp: 230        ovli: 237       T       T               5       14
 * rp: 238 arp: 231        ovli: 238       C       C               5       14
 * rp: 239 arp: 232        ovli: 239       *       A       d       4       14      ouch
 * rp: 239 arp: 232        ovli: 240       *       C       d       3       14      ouch
 * rp: 239 arp: 232        ovli: 241       *       G       d       2       14      ouch
 * rp: 239 arp: 232        ovli: 242       *       C       d       1       14      ouch
 * rp: 239 arp: 232        ovli: 243       *       A       d       0       14
 * rp: 239 arp: 232        ovli: 244       *       T       d       0       14
 * rp: 239 arp: 232        ovli: 245       A       A               0       14
 * rp: 240 arp: 233        ovli: 246       C       C               0       14
 *
 * respectively
 *
 * rp: 233 arp: 226        ovli: 233       C       C               5       14
 * rp: 234 arp: 227        ovli: 234       T       T               5       14
 * rp: 235 arp: 228        ovli: 235       G       G               5       14
 * rp: 236 arp: 229        ovli: 236       C       C               5       14
 * rp: 237 arp: 230        ovli: 237       T       T               5       14
 * rp: 238 arp: 231        ovli: 238       C       C               5       14
 * rp: 239 arp: 232        ovli: 239       *       A       d       4       14      ouch
 * rp: 239 arp: 232        ovli: 240       *       C       d       0       14
 * rp: 239 arp: 232        ovli: 241       *       G       d       0       14
 * rp: 239 arp: 232        ovli: 242       *       C       d       0       14
 * rp: 239 arp: 232        ovli: 243       *       A       d       0       14
 * rp: 239 arp: 232        ovli: 244       *       T       d       0       14
 * rp: 239 arp: 232        ovli: 245       A       A               0       14
 * rp: 240 arp: 233        ovli: 246       C       C               0       14
 *
 * when using a minimum belief value of 4 in pass 5.
 *
 * minedits 3 solves the problem, but at the price of
 * - slightly more coverage needed
 * - and worse performance in lower coverage areas (14x and lower)
 *
 * minimum belief value of 5 would also solve the problem (in this case),
 *  but not for slightly longer microrepeats. Also would worsen performance in
 *  lower coverage areas.
 *
 *
 * So, try the following: in polishing passes (>=5), search for gap / base
 *  discrepancies where at least one position in the bkmar is 0 (maybe also
 *  one pos before or behind the stretch?) and all the  bkmor are >= 6 (or
 *  other, higher value?)
 * If that is the case, set all bkmar positions in the affected gap / base stretch to 0
 * Like so:
 *
 * rp: 237 arp: 230        ovli: 237       T       T               5       14
 * rp: 238 arp: 231        ovli: 238       C       C               5       14
 * rp: 239 arp: 232        ovli: 239       *       A       d       0       14      // newly set to 0, no ouch!
 * rp: 239 arp: 232        ovli: 240       *       C       d       0       14
 * rp: 239 arp: 232        ovli: 241       *       G       d       0       14
 * rp: 239 arp: 232        ovli: 242       *       C       d       0       14
 * rp: 239 arp: 232        ovli: 243       *       A       d       0       14
 * rp: 239 arp: 232        ovli: 244       *       T       d       0       14
 * rp: 239 arp: 232        ovli: 245       A       A               0       14
 * rp: 240 arp: 233        ovli: 246       C       C               0       14
 *
 * Then the normal fillCorrector routines should be able to handle this correctly
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
void Assembly::tpbc_tweakBKMAR(readid_t actrid, readid_t otherrid, AlignedDualSeq & ads, std::vector<uint8> & bkmar,std::vector<uint8> & bkmor)
{
  FUNCSTART("void Assembly::tpbc_tweakBKMAR(readid_t actrid, readid_t otherrid, AlignedDualSeq & ads, std::vector<uint8> & bkmar,std::vector<uint8> & bkmor)");
  BUGIFTHROW(ads.getSequenceDirection(actrid)<0,"ads.getSequenceDirection(actrid)<0");

  bool docebug=false;
  if(false
     || AS_readpool[actrid].getName()==readofinterest1
     || AS_readpool[actrid].getName()==readofinterest2) docebug=true;
  if(false
     || AS_readpool[otherrid].getName()==readofinterest1
     || AS_readpool[otherrid].getName()==readofinterest2) docebug=true;

  CEBUG("tweak " << AS_readpool[actrid].getName() << " " << AS_readpool[otherrid].getName() << '\n');
  auto aseq = ads.getSequenceAtOverlapStart(actrid);
  auto oseq = ads.getSequenceAtOverlapStart(otherrid);
  for(uint32 ovli=0; ovli<ads.getOverlapLen(); ++ovli, ++aseq, ++oseq){
    // beware: no 'continue' in this loop! because of  "if(*aseq != '*') ++ccI;" at end

    CEBUG("ovli: " << ovli << '\t' << *aseq << '\t' << *oseq << '\t');
    if(*aseq!=*oseq){
      CEBUG('d');
    }else{
      CEBUG(' ');
    }

    CEBUG('\t' << static_cast<int16>(bkmar[ovli]));
    CEBUG('\t' << static_cast<int16>(bkmor[ovli]));

    const uint8 minbkmor=6;

    // this will run for every pos of a gap / base stretch ... I don't care atm
    if(*aseq!=*oseq
       && *aseq=='*'
       && bkmar[ovli]==0
       && bkmor[ovli]>=minbkmor){
      CEBUG("\tptweak");
      auto * checkptr=aseq;
      uint32 dx=0;
      // move to begin of run
      while(ovli-dx>=0 && *checkptr=='*' && bkmor[ovli-dx]>=minbkmor) {--dx; --checkptr;}
      checkptr=aseq;
      uint32 dy=0;
      // move to end of run
      while(ovli+dy<ads.getOverlapLen() && *checkptr=='*' && bkmor[ovli+dy]>=minbkmor) {++dy; ++checkptr;}

      CEBUG(" " << ovli-dx << " - " << ovli+dy);
      // eliminate bkmar values
      for(uint32 ei=ovli-dx; ei<ovli+dy; ++ei) bkmar[ei]=0;
    }

    CEBUG(endl);
  }
}
//#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
void Assembly::tpbc_fc_makeBelieveKMERMap(readid_t arid, readid_t orid, AlignedDualSeq & ads, std::vector<uint8> & bkmm, int32 kmersizeused, int32 additionalbelieveborder)
{
  FUNCSTART("void Assembly::tpbc_fc_makeBelieveKMERMap(readid_t arid, readid_t orid, AlignedDualSeq & ads, std::vector<uint8> & bkmm, int32 kmersizeused, int32 additionalbelieveborder)");

  // may have been reason for an error dumped:
  // because of for loop below: for(int32 ovli=0; ovli<ads.getOverlapLen()-kmersizeused; ...)
  if(kmersizeused>=ads.getOverlapLen()) return;

  try{
    bool docebug=false;
    if(false
       || AS_readpool[arid].getName()==readofinterest1
       || AS_readpool[arid].getName()==readofinterest2) docebug=true;
    if(false
       || AS_readpool[orid].getName()==readofinterest1
       || AS_readpool[orid].getName()==readofinterest2) docebug=true;

    CEBUG("makeBelieve " << AS_readpool[arid].getName() << "\t" << arid << " " << orid << '\n' << ads << '\n' << "adsend" << '\n');

    auto offsettooverlap=ads.getOffsetInAlignment(orid);
    if(offsettooverlap==0) offsettooverlap=ads.getOffsetInAlignment(arid);

    CEBUG("offsettooverlap: " << offsettooverlap << '\n');

    int32 offsetintoread=0;
    int32 deltax=1;
    bool isfwd=ads.getSequenceDirection(arid)>0;
    if(isfwd){
      offsetintoread=AS_readpool[arid].getLeftClipoff()+ads.getOffsetInAlignment(orid);
    }else{
      deltax=-1;
      offsetintoread=AS_readpool[arid].getRightClipoff()-1-ads.getOffsetInAlignment(orid);
    }

    auto & actread=AS_readpool[arid];

    //Read::setCoutType(Read::AS_TEXT);
    //CEBUG(actread << '\n');

    auto * seq=ads.getSequenceAtOverlapStart(arid);
    for(int32 ovli=0; ovli<ads.getOverlapLen()-kmersizeused; ++ovli, ++seq){
      CEBUG("ovli " << ovli << "\toir " << offsetintoread << "\tarp: " << AS_readpool[arid].getAdjustmentPosOfReadPos(offsetintoread) << "\tseq: " << *seq);
      if(*seq != '*'){
	bool isok=true;
	if(isfwd){
	  if(!actread.getBPosHashStats(offsetintoread).fwd.isValid()
	     || actread.getBPosHashStats(offsetintoread).fwd.getFrequency()<2
	     || !actread.getBPosHashStats(offsetintoread).fwd.hasConfirmedFwdRev()){
	    isok=false;
	  }
	}else{
	  if(!actread.getBPosHashStats(offsetintoread).rev.isValid()
	     || actread.getBPosHashStats(offsetintoread).rev.getFrequency()<2
	     || !actread.getBPosHashStats(offsetintoread).rev.hasConfirmedFwdRev()){
	    isok=false;
	  }
	}
	if(isok){
	  // first need to find out how many bases we will really need to cover as
	  //  the overlap string can be AA*********ACGATACATATATCCTA which is significantly
	  //  different from 'kmersizeused'
	  uint32 checklen=0;
	  {
	    auto * cseq=seq;
	    for(uint32 seenbases=0; seenbases<kmersizeused && ovli+checklen < ads.getOverlapLen(); ++checklen, ++cseq){
	      if(*cseq != '*' && *cseq != '#') {
		++seenbases;
	      }
	    }
	  }

	  if(checklen>1){
	    uint32 leftborder=1;
	    for(;leftborder<checklen && (*(seq+leftborder)=='*' || *(seq+leftborder)=='#'); ++leftborder) {}
	    char bseq=*(seq+leftborder);
	    for(; leftborder<checklen && *(seq+leftborder) == bseq; ++leftborder) {}
	    uint32 rightborder=checklen-1;
	    for(;rightborder>0 && (*(seq+rightborder)=='*' || *(seq+rightborder)=='#'); --rightborder) {}
	    bseq=*(seq+rightborder);
	    for(; rightborder>0 && *(seq+rightborder) == bseq; --rightborder) {}
	    CEBUG("\tlen " << rightborder-leftborder);
	    CEBUG("\tfrom " << leftborder+ovli << "\tto " << rightborder+ovli);
	    leftborder+=additionalbelieveborder;
	    rightborder-=additionalbelieveborder;
	    for(; leftborder<=rightborder; ++leftborder){
	      BUGIFTHROW(ovli+leftborder>=bkmm.size(),"ovli " << ovli << " + leftborder " << leftborder << " >= bkmm.size() " << bkmm.size() << " ???");
	      bkmm[ovli+leftborder]+=1;
	    }
	  }
	}
	offsetintoread+=deltax;
      }
      CEBUG(endl);
    }
  }
  catch(Notify n){
    cout << "Houston, we have a problem.\n";
    cout << "arid: " << arid << "\torid: " << orid << endl;
    cout << "bkmm.size(): " << bkmm.size()
	 << "\tkmersizeused:" << kmersizeused
	 << "\tadditionalbelieveborder:" << additionalbelieveborder << endl;
    Read::setCoutType(Read::AS_FASTQ);
    cout << AS_readpool[arid] << endl;
    cout << AS_readpool[orid] << endl;
    cout << "ADS:\n" << ads << endl;
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::tpbc_fc_minimumBelieveKMerMap(std::vector<uint8> & bkmm, uint8 minbk)
{
  for(auto & bke : bkmm){
    if(bke<minbk) bke=0;
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
uint32 Assembly::tpbc_generateRLEEdits(uint32 actpass, readid_t actrid, std::vector<std::vector<uint32>> & rlecc)
{
  FUNCSTART("uint32 Assembly::tpbc_generateRLEEdits(uint32 actpass, readid_t actrid, std::vector<std::vector<uint32>> & rlecc)");

  BUGIFTHROW(AS_readpool[actrid].getRLEValues()==nullptr,"AS_readpool[actrid].getRLEValues()==nullptr");

  bool docebug=false;
  //bool docebug=true;
  if(false
     || AS_readpool[actrid].getName()==readofinterest1
     || AS_readpool[actrid].getName()==readofinterest2) docebug=true;

  auto & actread=AS_readpool[actrid];

  CEBUG("ccr " << actread.getName() << "\tlc: " << actread.getLeftClipoff() << "\trc: " << actread.getRightClipoff() << '\n');

  // first, add our own RLE values to the lot
  {
    auto * aseq=actread.getClippedSeqAsChar();
    auto * eptr=aseq+actread.getLenClippedSeq();
    auto rvI=(actread.getRLEValues())->begin()+actread.getLeftClipoff();
    auto rqI=actread.getQualities().begin()+actread.getLeftClipoff();
    auto rccI=rlecc.begin();
    uint32 actpos=actread.getLeftClipoff();
    for(; aseq!=eptr; ++aseq, ++rvI, ++rqI, ++rccI, ++actpos){
      BUGIFTHROW(rccI==rlecc.end(),"rccI==rlecc.end()");
      if(dptools::isValidACGTBase(*aseq) && isupper(*aseq) && *rqI!=0){
	// same rule as for filling RLE corrector: need to have adjustments left and right concur
	if(actpos>0 && actpos<actread.getLenClippedSeq()-1){
	  if(actread.getAdjustmentPosOfReadPos(actpos)-1 == actread.getAdjustmentPosOfReadPos(actpos-1)
	     && actread.getAdjustmentPosOfReadPos(actpos)+1 == actread.getAdjustmentPosOfReadPos(actpos+1)){
	    rccI->push_back(*rvI);
	  }
	}
      }
    }
    BUGIFTHROW(rccI!=rlecc.end(),"rccI!=rlecc.end()");
  }

  // get to the fun stuff: if all RLE values the same, perfect. If not, then let's make a choice.
  {
    auto * aseq=actread.getClippedSeqAsChar();
    uint32 actpos=actread.getLeftClipoff();
    auto rvI=(const_cast<std::vector<uint8> *>(actread.getRLEValues()))->begin()+actread.getLeftClipoff();
    for(auto rccI=rlecc.begin(); rccI!=rlecc.end(); ++rccI, ++actpos, ++aseq, ++rvI){
      CEBUG("rpos: " << rccI-rlecc.begin() << '\t' << rccI->size() << '\t' << *aseq << "\tRLE" << static_cast<uint16>(*rvI));
      if(rccI->empty()){
	CEBUG("\tempty\n");
      }else{
	bool allsame=true;
	auto firstval=rccI->front();
	for(auto & rval : *rccI){
	  if(rval!=firstval){
	    allsame=false;
	    break;
	  }
	}
	uint32 choice=0;
	if(allsame){
	  CEBUG("\tall same: " << firstval);
	  choice=firstval;
	}else{
	  CEBUG("\tnot same {");
	  sort(rccI->begin(), rccI->end());
	  uint32 currval=rccI->front();
	  uint32 maxval=currval;
	  uint32 currcount=0;
	  uint32 maxcount=0;
	  uint64 totalval=0;
	  bool istie=false;
	  for(auto & rval : *rccI){
	    CEBUG(" " << rval);
	    if(rval==currval){
	      ++currcount;
	    }else{
	      if(maxcount==currcount){
		istie=true;
	      }else if(currcount>maxcount){
		maxcount=currcount;
		maxval=currval;
		istie=false;
	      }
	      currval=rval;
	      currcount=1;
	    }
	    totalval+=rval;
	  }
	  if(maxcount==currcount){
	    istie=true;
	  }else if(currcount>maxcount){
	    maxcount=currcount;
	    maxval=currval;
	    istie=false;
	  }
	  float avg=static_cast<float>(totalval)/rccI->size();
	  uint32 avgi=avg+0.5f;
	  uint32 median=0;
	  if(rccI->size() % 2){
	    median=(*rccI)[rccI->size()/2];
	  }else{
	    median=static_cast<uint32>(static_cast<float>((*rccI)[rccI->size()/2]+(*rccI)[rccI->size()/2-1])/2+0.5f);
	  }
	  CEBUG(" }");
	  CEBUG("\tmxc: " << maxcount << "/" << maxval << " med: " << median << " avg: " << avg);

	  if(avgi==median && median == maxval){
	    CEBUG(" uni");
	    choice=avgi;
	  }else{
	    CEBUG(" shoot");
	    if(avgi==median ||
	       avgi==maxval){
	      choice=avgi;
	    }else if(median==maxval){
	      choice=median;
	    }else{
	      // ooooops, all three are different???
	      // atm take median
	      CEBUG("wtf");
	      choice=median;
	    }
	  }
	}
	CEBUG("\tch: " << choice);
	if(choice!=*rvI){
	  CEBUG(" rlech!");
	  *rvI=choice;
	}
	CEBUG(endl);
      }
    }
  }

  return 0;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) { cout << bla; cout.flush();}
uint32 Assembly::tpbc_generateBaseEdits(uint32 actpass, readid_t actrid, std::vector<pbcounts_t> & correctorcounts)
{
  FUNCSTART("uint32 Assembly::generateBaseEdits(readid_t actrid, std::vector<pbcounts_t> & correctorcounts)");

  bool docebug=false;
  //bool docebug=true;
  if(false
     || AS_readpool[actrid].getName()==readofinterest1
     || AS_readpool[actrid].getName()==readofinterest2) docebug=true;

  uint32 editsmade=0;

  auto & actread=AS_readpool[actrid];

//*
  // get rid of conflict positions (+/-1) as this could be a clear sign of 'problems'
  // well, except in n-stretches, where alignment uncertainties may lead to disagreements

  //if(actpass!=5)

  uint32 numminedits=3;
  bool prekillconflict=true;
  bool preremovelower=false;

  if(actpass>=5) {
    preremovelower=true;
    numminedits=4;
  }

  CEBUG("ccu vector for " << actrid << '\t' << actread.getName() << '\n');
  {
    pbcounts_t pbczero;
    auto * aseq=actread.getClippedSeqAsChar();
    uint32 actcpos=0;
    for(auto ccI=correctorcounts.begin(); ccI != correctorcounts.end(); ++ccI, ++aseq, ++actcpos){
      CEBUG("acp: " << actcpos << '\t' << actread.getAdjustmentPosOfReadPos(actread.calcClippedPos2RawPos(actcpos)) << '\t' << actread.getBaseInClippedSequence(actcpos) << '\t' << *ccI << '\t');
      if(toupper(*aseq)!='N'){

	if(preremovelower){
	  if(ccI->cxa && ccI->cxa<numminedits){
	    ccI->cxa=0;
	    CEBUG(" prlA");
	  }
	  if(ccI->cxc && ccI->cxc<numminedits){
	    ccI->cxc=0;
	    CEBUG(" prlC");
	  }
	  if(ccI->cxg && ccI->cxg<numminedits){
	    ccI->cxg=0;
	    CEBUG(" prlG");
	  }
	  if(ccI->cxt && ccI->cxt<numminedits){
	    ccI->cxt=0;
	    CEBUG(" prlT");
	  }
	  if(ccI->cxgap && ccI->cxgap<numminedits){
	    ccI->cxgap=0;
	    CEBUG(" prl*");
	  }
	  if(ccI->cxia.size()>1){
	    auto iI=ccI->cxia.begin();
	    while(iI!=ccI->cxia.end()){
	      if(iI->count < 4){
		CEBUG(" prlI " << iI->what);
		iI=ccI->cxia.erase(iI);
	      }else{
		++iI;
	      }
	    }
	  }
	}

	uint8 epos=0;
	epos+=ccI->cxa>=numminedits;
	epos+=ccI->cxc>=numminedits;
	epos+=ccI->cxg>=numminedits;
	epos+=ccI->cxt>=numminedits;
	epos+=ccI->cxgap>=numminedits;
	// Nope, can have good inserts with good edits!
	// epos+=(ccI->cxia.size()==1 && ccI->cxia.front().count>4);
	if(epos>1){
	  if(ccI!=correctorcounts.begin()) *(ccI-1)=pbczero;
	  *ccI=pbczero;
	  if(ccI+1!=correctorcounts.end()) *(ccI+1)=pbczero;
	}
      }
      CEBUG(endl);
    }
  }
//*/

  CEBUG("cc vector for " << actrid << '\t' << actread.getName() << '\n');

  uint32 actcpos=correctorcounts.size()-1;
  for(auto ccI=correctorcounts.rbegin(); ccI != correctorcounts.rend(); ++ccI, --actcpos){
    CEBUG("acp: " << actcpos << '\t' << actread.getAdjustmentPosOfReadPos(actread.calcClippedPos2RawPos(actcpos)) << '\t' << actread.getBaseInClippedSequence(actcpos) << '\t' << *ccI);

    uint8 epos=0;
    char cchar=toupper(actread.getBaseInClippedSequence(actcpos));

    // due to the 'N'-stretch uncertainty, we need to act with a detour
    auto maxvotes=std::max(ccI->cxa,std::max(ccI->cxc,std::max(ccI->cxg,std::max(ccI->cxt,ccI->cxgap))));
    if(maxvotes>=numminedits){
      if(ccI->cxa == maxvotes){
	++epos;
	actread.changeBaseInClippedSequence('a',0,actcpos);

	// TODO: bang here, kill the RLE values, too!

	CEBUG("\tccedit A");
      }
      if(ccI->cxc == maxvotes){
	++epos;
	actread.changeBaseInClippedSequence('c',0,actcpos);
	CEBUG("\tccedit C");
      }
      if(ccI->cxg == maxvotes){
	++epos;
	actread.changeBaseInClippedSequence('g',0,actcpos);
	CEBUG("\tccedit G");
      }
      if(ccI->cxt == maxvotes){
	++epos;
	actread.changeBaseInClippedSequence('t',0,actcpos);
	CEBUG("\tccedit T");
      }
      if(ccI->cxgap >= maxvotes){
	++epos;
	actread.changeBaseInClippedSequence('*',0,actcpos);
	CEBUG("\tccedit newgap");
      }
      if(epos>1){
	if(cchar=='N'){
	  // whoooops, shouldn't happen all to often. But if it does, undo the edit
	  actread.changeBaseInClippedSequence('n',0,actcpos);
	  CEBUG(" ccedit UNDO");
	  // make sure the potential insert does not happen
	  cchar='Z';
	}else{
	  BUGIFTHROW(true,"epos>1, no N ???");
	}
      }
    }

    // do we need to insert something?
    // it is perfectly valid to have a change of a base plus insert,
    //  but for everything else, no inserts please
    if(ccI->cxia.size() == 1 && ccI->cxia.front().count >= 4){
      CEBUG("\tccedit insert " << ccI->cxia.front().what);
      for(auto sI=ccI->cxia.front().what.rbegin(); sI!=ccI->cxia.front().what.rend(); ++sI){
	actread.insertBaseInClippedSequence(tolower(*sI),0,actcpos,false);
      }
      ++epos;
    }
    editsmade+=epos;

    CEBUG(endl);
  }

  return editsmade;
}

// ecopbownH5_26908  at ~ 550 ... probably due to ecopbownH5_15286 low complexity errors at end

// small_assembly_orig: gap4 @ 13340 ... some TGGT reads miscorrected to T*GT after 5 rounds :-(((
//    +15 ecopbownH5_25906
//    +16 ecopbownH5_29223
//    -17 ecopbownH5_38640




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::killNonHAFCoveredStretchesInReads()
{
  FUNCSTART("void Assembly::killNonHAFCoveredStretchesInReads()");

  std::vector<uint8> hafcov;
  for(uint32 ri=0; ri<AS_readpool.size();++ri){
    auto & actread=AS_readpool[ri];
    hafcov.clear();
    hafcov.resize(actread.getLenClippedSeq(),0);
    for(auto & te : actread.getTags()){
      if(te.identifier==Read::REA_tagentry_idHAF2
	 || te.identifier==Read::REA_tagentry_idHAF3
	 || te.identifier==Read::REA_tagentry_idHAF4
	 || te.identifier==Read::REA_tagentry_idHAF5
	 || te.identifier==Read::REA_tagentry_idHAF6
	 || te.identifier==Read::REA_tagentry_idHAF7){
	int32 from=te.from;
	from-=actread.getLeftClipoff();
	int32 to=te.to;
	to-=actread.getLeftClipoff();
	for(; from<=to; ++from){
	  if(from>=0 && from<hafcov.size()) hafcov[from]=1;
	}
      }
    }

    uint32 nummasked=0;
    for(uint32 hcpos=0; hcpos<hafcov.size(); ++hcpos){
      if(hafcov[hcpos]==0) {
	++nummasked;
	actread.changeBaseInClippedSequence('n',0,hcpos);

	// TODO: bang here, kill the RLE values, too!

      }
    }
    float maskratio=100.0/hafcov.size()*nummasked;
    cout << "NMask: " << ri << '\t' << actread.getName() << '\t' << nummasked << " / " << hafcov.size()
	 << " (" << maskratio << "%)";
    if(maskratio>50.0){
      cout << " complete kill";
      //for(uint32 hcpos=0; hcpos<hafcov.size(); ++hcpos){
      //	actread.changeBaseInClippedSequence('n',0,hcpos);
      //}
      actread.setRQClipoff(actread.getLeftClipoff());
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

//#define CEBUG(bla) { cout << bla; cout.flush();}
void Assembly::tpbc_generateNonAlignClips(readid_t actrid, std::vector<pbcounts_t> & correctorcounts)
{
  FUNCSTART("uint32 Assembly::generateNonAlignClips(readid_t actrid, std::vector<pbcounts_t> & correctorcounts)");


  BUGIFTHROW(AS_naclipl.size()!=AS_readpool.size(),"AS_naclipl.size() " << AS_naclipl.size() << " != AS_readpool.size() " << AS_readpool.size());
  BUGIFTHROW(AS_naclipr.size()!=AS_readpool.size(),"AS_naclipr.size() " << AS_naclipr.size() << " != AS_readpool.size() " << AS_readpool.size());

  // *sigh* through indels occurring in the generatePBEdits() call prior to the call to this function,
  //  things may be slightly out of sync
  // just make sure this is called late in the whole correction proces, ideally (next to) last pass

  //bool docebug=true;
  bool docebug=false;

  CEBUG("generateNonAlignClips GNA\n");

  auto & actread=AS_readpool[actrid];

  // signal processing tmp vector
  // we need this tmp vectory as the erode/dilate operations further down below may
  //  temporarily expand into "non-existing" space left and right of a vector
  //  just of size correctorcounts.size(), and then this gives wrong results at
  //  the borders in the subsequent operation
  //
  // we could get away with just using one vector without copying, but atm this is
  //  prototype code anyway. Should this ever become bottleneck in the PacBio correction
  //  process (I don't think so), re-engineer.

  uint32 borderpad=100;
  std::vector<uint32> sptmpv(correctorcounts.size()+2*borderpad,0);

//#define GIMMEWIG
#ifdef GIMMEWIG
  std::string wigname=actread.getName();
  std::vector<uint32> cqtmpv(correctorcounts.size()+2*borderpad,0);
#endif

  {

    {
      auto scI=sptmpv.begin()+borderpad;
      auto ccI=correctorcounts.begin();
#ifdef GIMMEWIG
      auto qcI=cqtmpv.begin();
      for(; ccI!=correctorcounts.end(); ++scI, ++ccI, ++qcI){
	if(ccI->cqsafe>=2) *scI=1;
	*qcI=ccI->cqsafe;
      }
#else
      for(; ccI!=correctorcounts.end(); ++scI, ++ccI){
	if(ccI->cqsafe>=2) *scI=1;
      }
#endif
    }

#ifdef GIMMEWIG
    dbgContainerToWiggle(cqtmpv,wigname+"_cq",wigname+"_cq");
    dbgContainerToWiggle(sptmpv,wigname+"_raw",wigname+"_raw");
#endif

    erode(sptmpv.begin(),sptmpv.end(),20);
    erode(sptmpv.rbegin(),sptmpv.rend(),20);
    dilate(sptmpv.begin(),sptmpv.end(),70);
    dilate(sptmpv.rbegin(),sptmpv.rend(),70);
    erode(sptmpv.begin(),sptmpv.end(),200);
    erode(sptmpv.rbegin(),sptmpv.rend(),200);
    dilate(sptmpv.begin(),sptmpv.end(),150);
    dilate(sptmpv.rbegin(),sptmpv.rend(),150);

#ifdef GIMMEWIG
    dbgContainerToWiggle(sptmpv,wigname+"_sp",wigname+"_sp");
#endif
  }

  // create safecounts from really used subpart of sptmpv
  // atm as copy, maybe need sptmpv later?
  std::vector<uint32> safecounts(sptmpv.begin()+borderpad,sptmpv.end()-borderpad);
  sptmpv.clear();

  {
    CEBUG("dsp for " << actread.getName() << endl);
    auto scI=safecounts.begin();
    auto ccI=correctorcounts.begin();
    for(; ccI!=correctorcounts.end(); ++scI, ++ccI){
      CEBUG(scI-safecounts.begin()
	    << '\t' << ccI->cqsafe
	    << '\t' << *scI
	    << endl);
    }
  }

  int32 numlclips=0;
  {
    auto scI=safecounts.begin();
    for(; scI!=safecounts.end() && *scI<1; ++scI) {}

    numlclips=scI-safecounts.begin();
    AS_naclipl[actrid]=actread.getLeftClipoff()+numlclips;
    CEBUG("GNA clip left " << actread.getName() << ": " << numlclips << "\tnl: " << AS_naclipl[actrid] << endl);
  }

  int32 numrclips=0;
  {
    auto scI=safecounts.rbegin();
    for(; scI!=safecounts.rend() && *scI<1; ++scI) {}
    numrclips=scI-safecounts.rbegin();
    AS_naclipr[actrid]=actread.getRightClipoff()-numrclips;
    CEBUG("GNA clip right " << actread.getName() << ": " << numrclips  << "\tnr: " << AS_naclipr[actrid] << endl);
  }

  // can happen, if whole read is potentially killed above
  if(AS_naclipr[actrid]<AS_naclipl[actrid]){
    AS_naclipl[actrid]=0;
    AS_naclipr[actrid]=0;
    CEBUG("GNA killed " << actread.getName() << endl);
  }


  // check for holes remaining and eventually choose longest stretch
  // TODO: later, not choose longest, but create subreads for stretches >=500bp
  if(AS_naclipl[actrid]<AS_naclipr[actrid]){
    auto scI=safecounts.begin()+numlclips;
    auto scE=safecounts.end()-numrclips;

    CEBUG("Checking " << scI-safecounts.begin() << ' ' << scE-safecounts.begin() << endl);
    uint32 longesthole=0;
    uint32 holelen=0;
    uint32 longeststretch=0;
    uint32 stretchlen=0;
    auto lsscE=scE;

    bool instretch=*scI>0;
    for(; scI<scE; ++scI){
      if(*scI>0){
	if(instretch){
	  // stretchcont
	  ++stretchlen;
	}else{
	  // hole2stretch
	  longesthole=std::max(longesthole,holelen);
	  stretchlen=1;
	  instretch=true;
	}
      }else{
	if(!instretch){
	  // hole cont
	  ++holelen;
	}else{
	  // stretch2hole
	  longeststretch=std::max(longeststretch,stretchlen);
	  if(stretchlen==longeststretch){
	    lsscE=scI;
	  }
	  holelen=1;
	  instretch=false;
	}
      }
    }
    if(instretch) {
      longeststretch=std::max(longeststretch,stretchlen);
      lsscE=scI;
    }else{
      longesthole=std::max(longesthole,holelen);
    }

    CEBUG("LH: " << longesthole << "\tLS: " << longeststretch << endl);

    if(longesthole>250){
      AS_naclipr[actrid]=actread.getLeftClipoff()+(lsscE-safecounts.begin());
      AS_naclipl[actrid]=AS_naclipr[actrid]-longeststretch;
      CEBUG("GNA: Has hole, using subread " << AS_naclipl[actrid] << ' ' << AS_naclipr[actrid] << endl);
      CEBUG("raw: " << lsscE-safecounts.begin()-longeststretch << " - " << lsscE-safecounts.begin() << endl);
    }
  }
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::colourReadsByRLE()
{
  FUNCSTART("void Assembly::colourReadsByRLE()");

  multitag_t mt;
  for(uint32 rid=0; rid<AS_readpool.size(); ++rid){
    auto & actread=AS_readpool[rid];
    if(actread.getRLEValues()==nullptr) continue;

    auto rvI=actread.getRLEValues()->begin();
    for(uint32 rp=0; rp<actread.getLenSeq(); ++rp, ++rvI){
      mt.identifier=Read::REA_tagentry_idMINF;
      switch(*rvI){
      case 0:
	break;
      case 1:
	mt.identifier=Read::REA_tagentry_idRLE1;
	break;
      case 2:
	mt.identifier=Read::REA_tagentry_idRLE2;
	break;
      case 3:
	mt.identifier=Read::REA_tagentry_idRLE3;
	break;
      case 4:
	mt.identifier=Read::REA_tagentry_idRLE4;
	break;
      case 5:
	mt.identifier=Read::REA_tagentry_idRLE5;
	break;
      case 6:
	mt.identifier=Read::REA_tagentry_idRLE6;
	break;
      case 7:
	mt.identifier=Read::REA_tagentry_idRLE7;
	break;
      default:
	// 8 or higher
	mt.identifier=Read::REA_tagentry_idRLE8;
	break;
      }
      if(mt.identifier!=Read::REA_tagentry_idMINF){
	mt.from=rp;
	mt.to=rp;
	actread.addTagO(mt);
      }
    }
  }
}

/*

PacBio Bell adaptor:
ATCTCTCTCttttcctcctcctccgttgttgttgttGAGAGAGAT

---------

!!! ecopbownH5_19195: typical case of sequenced around the adapter ... without adapter! *sigh*

Blasted at NCBI:

Xuzhou 21: 498120 - 497938
... TACGGTAATTCAAGGGGTGAGACCGAATGAGTCCCTTGCCCGCAGTCAGTCGGAACGTTAAATGGCAAATCAAGAACAATCAGGACAACCTTCACGGAAAGTCAATGTGCAGATCCAATGATAATGCCCCTGCCAACGTCCCTCCGCGCCCAGCGGAGTTGTGATATTCGTTACAGTTAAAGGAGTTTCGGTG

Xuzhou 21: 497913 - 497973
GGATACGTACCAAGCCTAACTGATCCACAGAAATCCCCTTTAACGTAACGGAAATCAACATCGCGGGCGAGACCGGCAGGGGCCATT
AATCAATTGATTCTGCAACATGATTTCCTGAATTGTCCTATTGCCCTGAATTTGCCCATATTCAACAGTTCCGCGTGCCGGTGAATGCGGGCAAGATCATACGACTCACCCTGAAACGTAGT...

I.e., a simple jump on the genome (25 bases ... length of inner core/loop of adapter?) and started in the other direction:
                         497938 <------------------------------------------------------------------------------
497913 ---------------------------------------------------------------------------------->

---------

//*/
