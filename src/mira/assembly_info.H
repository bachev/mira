/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2008 and later by Bastien Chevreux
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

#ifndef _mira_assemblyinfo_h_
#define _mira_assemblyinfo_h_

#include <iostream>

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include "mira/contig.H"

class AssemblyInfo
{
//Variables
public:
  uint32         ASI_conf_large_consize4stats; /* min size to calc statistics
						  for large contigs */
  uint32         ASI_conf_large_totalcov4stats; /* min total cov to calc statistics
						   for large contigs */

  // this is what user considers as large contig
  uint32         ASI_conf_large_mincontigsize; // minimum size of large contigs
  uint32         ASI_conf_large_mintotalcov; // minimum totcov of large contigs
  std::vector<uint32> ASI_conf_large_minavgcov_perst; /* minimum coverage per
						       sequencing type for
						       large contigs */
  std::deque<Contig::constats_t> ASI_contigstats;
  std::deque<std::string> ASI_contignames;

  std::deque<std::string> ASI_nameslargecontigs;

  uint32 ASI_numreads_total;
  uint32 ASI_numbases_total;

  uint32 ASI_numreads_assembled;     // all reads in contigs
  uint32 ASI_numbases_assembled;     // all bases of reads in contrigs

  uint32 ASI_numreads_singlets;      // numsinglets (+debris)
  uint32 ASI_numbases_singlets;

  // 0 = all contigs, 1 = large contigs
  // some values not used for all contigs
  uint32 ASI_numcontigs[2];
  uint32 ASI_sizeconsensus[2];
  uint32 ASI_numIUPACs[2];
  uint32 ASI_numSRMc[2];
  uint32 ASI_numWRMc[2];
  uint32 ASI_numSNPs[2];            // sum SROc, SAOc, SIOc
  uint32 ASI_numSTMU[2];
  uint32 ASI_numSTMS[2];

  uint32 ASI_maxcoverage[2];
  std::vector<uint32> ASI_maxcoverage_perst[2];
  double ASI_avgcoverage[2];
  std::vector<double> ASI_avgcoverage_perst[2];

  uint32 ASI_largestconsize[2];
  uint32 ASI_n50consize[2];
  uint32 ASI_n90consize[2];
  uint32 ASI_n95consize[2];
  uint32 ASI_avgconsize[2];
  uint32 ASI_avgconqual[2];

  uint32 ASI_numcon_noqualread[2];
  uint32 ASI_numcon_somequalreadmissing[2];

  //Functions
private:
  void foolCompiler();
  void discard();

  void zeroStats();
  void checkThesePrinter(std::ostream & ostr, uint32 val, const char * s);
  void dumpCurrentInfo(std::ostream & ostr);

public:
  AssemblyInfo();
  AssemblyInfo(AssemblyInfo const &other);
  ~AssemblyInfo();

  AssemblyInfo const & operator=(AssemblyInfo const & other);
  friend std::ostream & operator<<(std::ostream &ostr, AssemblyInfo &asi);

  void zeroInfo();
  void calcCurrentInfo();

  inline void setLargeContigSize(uint32 s) { ASI_conf_large_mincontigsize=s;};
  inline void setLargeTotalCov(uint32 s) { ASI_conf_large_mintotalcov=s;};
  void setLargeContigCovPerST(uint32 cov, uint8 seqtype);

  void setLargeContigSizeForStats(uint32 s) { ASI_conf_large_consize4stats=s;};
  void setLargeContigCovForStats(uint32 s) { ASI_conf_large_totalcov4stats=s;};

  void storeContigStats(const Contig::constats_t & cs, const std::string & cname);
  void dumpLargeContigNames(std::ostream & ostr);

  uint32 huntForSmileCoverage(std::string & warntext);
};


#endif
