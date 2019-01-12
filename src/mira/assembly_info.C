/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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


#include "version.H"

#include "mira/assembly_info.H"
#include "errorhandling/errorhandling.H"

#include "util/stlimprove.H"


// Plain vanilla constructor
AssemblyInfo::AssemblyInfo()
{
  FUNCSTART("AssemblyInfo::AssemblyInfo()");

  zeroInfo();

  ASI_conf_large_mincontigsize=0;
  ASI_conf_large_consize4stats=5000;
  ASI_conf_large_totalcov4stats=0;
  ASI_conf_large_mintotalcov=0;
  ASI_conf_large_minavgcov_perst.resize(ReadGroupLib::SEQTYPE_END,0);

  FUNCEND();
}

AssemblyInfo::~AssemblyInfo()
{
  FUNCSTART("AssemblyInfo::~AssemblyInfo()");

  discard();

  FUNCEND();
}


// nothing to do at the moment
void AssemblyInfo::discard()
{
  FUNCSTART("AssemblyInfo::discard()");

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//AssemblyInfo::AssemblyInfo(AssemblyInfo const &other)
//{
//  FUNCSTART("AssemblyInfo::AssemblyInfo(AssemblyInfo const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//AssemblyInfo const & AssemblyInfo::operator=(AssemblyInfo const & other)
//{
//  FUNCSTART("AssemblyInfo const & AssemblyInfo::operator=(AssemblyInfo const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

std::ostream & operator<<(std::ostream &ostr, AssemblyInfo & asi)
{
  FUNCSTART("friend std::ostream & AssemblyInfo::operator<<(std::ostream &ostr, const  &???)");

  asi.dumpCurrentInfo(ostr);

  FUNCEND();
  return ostr;
}


/*************************************************************************
 *
 * reset all statistics about the assembly
 *
 *
 *************************************************************************/

void AssemblyInfo::zeroInfo()
{
  ASI_contigstats.clear();
  ASI_nameslargecontigs.clear();
  ASI_contignames.clear();

  zeroStats();
}

void AssemblyInfo::zeroStats()
{
  ASI_numreads_total=0;
  ASI_numbases_total=0;

  ASI_numreads_assembled=0;
  ASI_numbases_assembled=0;

  ASI_numreads_singlets=0;
  ASI_numbases_singlets=0;

  for(uint32 i=0; i<2; i++){
    ASI_numcontigs[i]=0;
    ASI_sizeconsensus[i]=0;
    ASI_numIUPACs[i]=0;
    ASI_numSRMc[i]=0;
    ASI_numWRMc[i]=0;
    ASI_numSNPs[i]=0;
    ASI_numSTMS[i]=0;
    ASI_numSTMU[i]=0;

    ASI_maxcoverage[i]=0;
    ASI_maxcoverage_perst[i].clear();
    ASI_maxcoverage_perst[i].resize(ReadGroupLib::SEQTYPE_END,0);
    ASI_avgcoverage[i]=0;
    ASI_avgcoverage_perst[i].clear();
    ASI_avgcoverage_perst[i].resize(ReadGroupLib::SEQTYPE_END,0);

    ASI_largestconsize[i]=0;
    ASI_n50consize[i]=0;
    ASI_n90consize[i]=0;
    ASI_n95consize[i]=0;
    ASI_avgconsize[i]=0;
    ASI_avgconqual[i]=0;

    ASI_numcon_noqualread[i]=0;
    ASI_numcon_somequalreadmissing[i]=0;
  }

}


void AssemblyInfo::storeContigStats(const Contig::constats_t & cs, const std::string & cname)
{
  ASI_contigstats.push_back(cs);
  if(!cname.empty()){
    ASI_contignames.push_back(cname);
  }
}

void AssemblyInfo::setLargeContigCovPerST(uint32 cov, uint8 seqtype)
{
  FUNCSTART("void AssemblyInfo::setLargeContigCovPerST(uint32 cov, uint8 seqtype)");

  BUGIFTHROW(seqtype>=ReadGroupLib::SEQTYPE_END, "Illegal sequence type?");

  ASI_conf_large_minavgcov_perst.resize(ReadGroupLib::SEQTYPE_END);
  ASI_conf_large_minavgcov_perst[seqtype]=cov;

  FUNCEND();
}


void AssemblyInfo::calcCurrentInfo()
{
  FUNCSTART("void AssemblyInfo::calcCurrentInfo()");

  zeroStats();

  std::vector<double> avgtotalcov[2];
  std::vector<std::vector<std::vector<double> > > avgcovperst(2);
  avgcovperst[0].resize(ReadGroupLib::SEQTYPE_END);
  avgcovperst[1].resize(ReadGroupLib::SEQTYPE_END);

  std::vector<uint32> consizes[2];
  std::vector<uint32> conquals[2];

  ASI_nameslargecontigs.clear();
  auto cnI=ASI_contignames.cend();
  if(ASI_contignames.size()==ASI_contigstats.size()){
    cnI=ASI_contignames.cbegin();
  }

  auto csI=ASI_contigstats.begin();
  for(; csI != ASI_contigstats.end(); ++csI){
    // 0 == all contigs, 1 == large contigs: when minsize and cov have been reached
    uint32 coninfooffset=0;
    if(csI->conlength_nogap >= ASI_conf_large_mincontigsize){
      bool avgcovfound=false;
      if(ASI_conf_large_mintotalcov>0
	 && csI->avg_coverage >= ASI_conf_large_mintotalcov){
	avgcovfound=true;
      }else{
	for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	  if(ASI_conf_large_minavgcov_perst[st]>0
	     && csI->avg_covperst[st] >= ASI_conf_large_minavgcov_perst[st]){
	    avgcovfound=true;
	    break;
	  }
	}
      }
      if(avgcovfound) coninfooffset=1;
    }

    // for bacteria in growth phase / projects with very uneven coverage across genome,
    //  the above may fail (contigs with several kb below 50% of avg. coverage)
    // try a rescue strategy:
    // - non-rep contigs with twice the needed min size and with 2/3 of needed coverage are also taken
    // TODO: this could probably be merged with loop above, but kept apart atm in case logic needs to change
    if(coninfooffset==0                     // not taken yet
       && !csI->contains_long_repeats_only  // non reps
       && csI->conlength_nogap >= 2*ASI_conf_large_mincontigsize){   // 2x min length large contigs
      bool avgcovfound=false;
      if(ASI_conf_large_mintotalcov>0
	 && csI->avg_coverage>= ASI_conf_large_mintotalcov*(static_cast<double>(2.0)/3)){
	avgcovfound=true;
      }else{
	for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	  if(ASI_conf_large_minavgcov_perst[st]>0
	     && csI->avg_covperst[st] >= ASI_conf_large_minavgcov_perst[st]*(static_cast<double>(2.0)/3)){
	    avgcovfound=true;
	    break;
	  }
	}
      }
      if(avgcovfound) coninfooffset=1;
    }

    if(coninfooffset && cnI!=ASI_contignames.end()){
      ASI_nameslargecontigs.push_back(*cnI);
      csI->islargecontig=1;
    }else{
      csI->islargecontig=-1;
    }

    if(cnI!=ASI_contignames.end()){
      ++cnI;
    }

    ASI_numreads_assembled+=csI->total_reads;
    //ASI_numbases_assembled+=...

    // store stats for all contigs and also for large contigs
    //  in second loop if we have a large contig
    for(uint32 i=0; i<=coninfooffset; ++i){
      consizes[i].push_back(csI->conlength_nogap);
      conquals[i].push_back(csI->avg_conqual);

      ASI_numIUPACs[i]+=csI->IUPACinC;
      ASI_numSRMc[i]+=csI->numSRMc;
      ASI_numWRMc[i]+=csI->numWRMc;

      if(csI->total_reads>1) {
	ASI_numcontigs[i]++;
      }else if(i==0){
	// increase singlets only once (in all contigs loop)
	ASI_numreads_singlets++;
      }

      if(csI->total_reads==csI->numreads_noqual) {
	ASI_numcon_noqualread[i]++;
      }else if(csI->total_reads!=csI->numreads_withqual) {
	ASI_numcon_somequalreadmissing[i]++;
      }

      ASI_sizeconsensus[i]+=csI->conlength_nogap;
      ASI_largestconsize[i]=std::max(ASI_largestconsize[i],csI->conlength_nogap);
      ASI_maxcoverage[i]=std::max(ASI_maxcoverage[i],csI->max_coverage);
      if(csI->conlength_nogap>=ASI_conf_large_consize4stats && csI->avg_coverage>= ASI_conf_large_totalcov4stats) {
	avgtotalcov[i].push_back(csI->avg_coverage);
      }
      for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; ++st){
	ASI_maxcoverage_perst[i][st]=std::max(
	  ASI_maxcoverage_perst[i][st],csI->max_covperst[st]
	  );
	if(csI->conlength_nogap>=ASI_conf_large_consize4stats && csI->avg_coverage>= ASI_conf_large_totalcov4stats) {
	  avgcovperst[i][st].push_back(csI->avg_covperst[st]);
	}
      }
    }
  }

  for(uint32 i=0; i<2; i++){
    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
      if(avgcovperst[i].size()){
	// use non-parallel sort, this is a tiny thing
	mstd::ssort(avgcovperst[i][st]);
	auto acI=avgcovperst[i][st].begin();
	for(; acI!=avgcovperst[i][st].end() && *acI==0.0; ++acI);
	if(acI != avgcovperst[i][st].end()){
	  ASI_avgcoverage_perst[i][st]=
	    avgcovperst[i][st][(avgcovperst[i][st].end()-acI)/2];
	}
      }
    }

    if(!consizes[i].empty()){
    // use non-parallel sort, this is a tiny thing
      mstd::ssort(consizes[i], std::greater<uint32>());

      ASI_n50consize[i]=0;
      ASI_n90consize[i]=0;
      ASI_n95consize[i]=0;

      uint32 totalconsize=0;
      for(auto & cse : consizes[i]){
	totalconsize+=cse;
	if(ASI_n50consize[i]==0
	   && totalconsize >= (0.5 * ASI_sizeconsensus[i])){
	  ASI_n50consize[i]=cse;
	}
	if(ASI_n90consize[i]==0
	   && totalconsize >= (0.9 * ASI_sizeconsensus[i])){
	  ASI_n90consize[i]=cse;
	}
	if(ASI_n95consize[i]==0
	   && totalconsize >= (0.95 * ASI_sizeconsensus[i])){
	  ASI_n95consize[i]=cse;
	}
      }
    }

    if(!conquals[i].empty()){
      // use non-parallel sort, this is a tiny thing
      mstd::ssort(conquals[i], std::greater<uint32>());
      ASI_avgconqual[i]=conquals[i][(conquals[i].size())/2];
    }
    if(!avgtotalcov[i].empty()){
      // use non-parallel sort, this is a tiny thing
      mstd::ssort(avgtotalcov[i], std::greater<uint32>());
      ASI_avgcoverage[i]=avgtotalcov[i][(avgtotalcov[i].size())/2];
    }
  }

  FUNCEND();
}


void AssemblyInfo::dumpLargeContigNames(std::ostream & ostr)
{
  for(auto & lce : ASI_nameslargecontigs){
    ostr << lce << '\n';
  }
}

void AssemblyInfo::dumpCurrentInfo(std::ostream & ostr)
{
  calcCurrentInfo();

  ostr << "Assembly information:\n"
       << "=====================\n\n";

  dateStamp(ostr);
  ostr << "MIRA version: " << miraversion << "\n\n";

  ostr << "Num. reads assembled: " << ASI_numreads_assembled << '\n';
  ostr << "Num. singlets: " << ASI_numreads_singlets << '\n';

  ostr.setf(std::ios::fixed, std::ios::floatfield);
  ostr.setf(std::ios::showpoint);
  ostr.precision(2);

  for(int32 i=1; i>=0; i--){
    if(i==1){
      ostr << "\n\nCoverage assessment (calculated from contigs >= "
	   << ASI_conf_large_consize4stats
	   << " with coverage >= " << ASI_conf_large_totalcov4stats
	   << "):\n=========================================================\n";
      ostr << "  Avg. total coverage: " << ASI_avgcoverage[i];
      ostr << "\n  Avg. coverage per sequencing technology";
      for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	ostr << "\n\t" << ReadGroupLib::getNameOfSequencingType(st)
	     << ":\t" << ASI_avgcoverage_perst[i][st];
      }

      ostr << "\n\n\nLarge contigs (makes less sense for EST assemblies):\n"
	   << "====================================================\nWith\tContig size\t\t>= " << ASI_conf_large_mincontigsize
	   << "\n\tAND (Total avg. Cov\t>= " << ASI_conf_large_mintotalcov;
      for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	ostr << "\n\t     OR Cov(" << ReadGroupLib::getShortNameOfSequencingType(st)
	     << ")\t>= " << ASI_conf_large_minavgcov_perst[st];
      }
      ostr << "\n\t    )\n\n";
    }else{
      ostr << "\nAll contigs:\n"
	   << "============\n";
    }
    ostr << "  Length assessment:\n  ------------------\n"
	 << "  Number of contigs:\t" << ASI_numcontigs[i]
	 << "\n  Total consensus:\t" << ASI_sizeconsensus[i]
	 << "\n  Largest contig:\t" << ASI_largestconsize[i]
	 << "\n  N50 contig size:\t" << ASI_n50consize[i]
	 << "\n  N90 contig size:\t" << ASI_n90consize[i]
	 << "\n  N95 contig size:\t" << ASI_n95consize[i];

    ostr << "\n\n  Coverage assessment:\n  --------------------\n"
	 << "  Max coverage (total):\t"
	 << ASI_maxcoverage[i];
    ostr << "\n  Max coverage per sequencing technology";
    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
      ostr << "\n\t" << ReadGroupLib::getNameOfSequencingType(st)
	   << ":\t" << ASI_maxcoverage_perst[i][st];
    }
//    ostr << "\n  Avg. total coverage (size >= "
//	 << ASI_conf_large_consize4stats << "): " << ASI_avgcoverage[i];
//    ostr << "\n  Avg. coverage (contig size >= "<< ASI_conf_large_consize4stats << ")";
//    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
//      ostr << "\n\t" << Read::getNameOfSequencingType(st)
//	   << ":\t" << ASI_avgcoverage_perst[i][st];
//    }

    ostr << "\n\n  Quality assessment:\n  -------------------";
    ostr << "\n  Average consensus quality:\t\t\t" << ASI_avgconqual[i];
    checkThesePrinter(ostr, ASI_numIUPACs[i],
		      "\n  Consensus bases with IUPAC:\t\t\t");
    checkThesePrinter(ostr, ASI_numSRMc[i],
		      "\n  Strong unresolved repeat positions (SRMc):\t");
    checkThesePrinter(ostr, ASI_numWRMc[i],
		      "\n  Weak unresolved repeat positions (WRMc):\t");
    checkThesePrinter(ostr, ASI_numSTMU[i],
		      "\n  Sequencing Type Mismatch Unsolved (STMU):\t");
    checkThesePrinter(ostr, ASI_numcon_noqualread[i],
		      "\n  Contigs having only reads wo qual:\t\t");
    checkThesePrinter(ostr, ASI_numcon_somequalreadmissing[i],
		      "\n  Contigs with reads wo qual values:\t\t");
    ostr << "\n\n";
  }
}

void AssemblyInfo::checkThesePrinter(std::ostream & ostr, uint32 val, const char * s)
{
  ostr << s << val;
  if(val==0){
    ostr << "\t(excellent)";
  }else{
    ostr << "\t(you might want to check these)";
  }
}



/*************************************************************************
 *
 * return warning level of smile coverage, 0== no warning
 *
 *************************************************************************/
uint32 AssemblyInfo::huntForSmileCoverage(std::string & warntext)
{
  FUNCSTART("uint32 AssemblyInfo::huntForSmileCoverage(std::string & warntext)");
  uint32 retvalue=0;

  auto covtotake=ASI_avgcoverage[1];
  if(covtotake==0.0) covtotake=ASI_avgcoverage[0];

  auto thmin=covtotake*static_cast<double>(.75);
  auto thmax=covtotake*static_cast<double>(1.25);

  uint32 num_lnr=0;   // lnr= large, non-rep contigs
  uint32 num_below=0;
  uint32 num_above=0;

  uint64 sumsize_lnr=0;
  uint64 sumsize_below=0;
  uint64 sumsize_above=0;

  for(const auto & cse : ASI_contigstats){
    BUGIFTHROW(cse.islargecontig==0,"cse.islargecontig==0 ???");
    if(!cse.contains_long_repeats_only && cse.islargecontig>0){
      ++num_lnr;
      sumsize_lnr+=cse.conlength_nogap;
      if(cse.avg_coverage<thmin){
	++num_below;
	sumsize_below+=cse.conlength_nogap;
      }else if(cse.avg_coverage>thmax){
	++num_above;
	sumsize_above+=cse.conlength_nogap;
      }
    }
  }

  std::ostringstream wmsg;

  // test for proportion of contigs below ...
  if(sumsize_below >= sumsize_lnr*static_cast<double>(0.15)){
    ++retvalue;
    wmsg << "- " << num_below << " contig(s) with a total of " << sumsize_below << " bases (= " << 100.0/sumsize_lnr*sumsize_below << "% of bases in all non-repetitive large contigs) have an average coverage less than 75% of the average coverage of all non-repetitive large contigs.\n";
  }
  // test for proportion of contigs above ...
  if(sumsize_above >= sumsize_lnr*static_cast<double>(0.15)){
    ++retvalue;
    wmsg << "- " << num_above << " contig(s) with a total of " << sumsize_above << " bases (= " << 100.0/sumsize_lnr*sumsize_above << "% of bases in all non-repetitive contigs) have an average coverage more than 125% of the average coverage of all non-repetitive large contigs.\n";
  }

  // test for proportion of contigs above ...
  if(sumsize_above+sumsize_below >= sumsize_lnr*static_cast<double>(0.3)){
    ++retvalue;
    wmsg << "- " << num_above+num_below << " contig(s) with a total of " << sumsize_above+sumsize_below << " bases (= " << 100.0/sumsize_lnr*(sumsize_above+sumsize_below) << "% of bases in all non-repetitive contigs) have an average coverage 25% above or below the average coverage of all non-repetitive large contigs.\n";
  }

  warntext.clear();
  if(retvalue){
    wmsg << "Summary: found " << retvalue << " indicator(s) for coverage problem(s).\n\nIf the DNA you are assembling is bacterial, this could indicate that you sampled and sequenced DNA from exponential or late exponential phase of a bacterial population. This leads to a coverage bias toward the origin of replication, hence false positive detection of repeats, hence an assembly which is more fragmented than it could be or may have misassemblies in regions located toward the opposite of the origin of replication."
      "\nOnly available countermeasure: for your next sequencing project, do not sample in exponential phase but sample in stationary phase (if possible).\n";
    warntext=wmsg.str();
  }

  return retvalue;
}
