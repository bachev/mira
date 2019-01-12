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


#include <fstream>

#include <string>
#include <vector>
#include <list>

#include "mira/parameters.H"
#include "mira/readpool_io.H"
#include "mira/hashstats.H"
#include "mira/vhash.H"


class MiraMer
{
private:
  std::string MER_job{"create"};
  std::string MER_outmhs{"miramer.mhs.gz"};

  uint32 MER_basesperhash=31;
  uint32 MER_rarekmerearlykill=0;
  uint32 MER_optthreads=0;

  HashStatistics<vhash64_t> MER_hs64;
  HashStatistics<vhash128_t> MER_hs128;
  HashStatistics<vhash256_t> MER_hs256;
  HashStatistics<vhash512_t> MER_hs512;

private:
  template<typename TVHASH_T>
  void mer_fhs_helper1(int argc, char ** argv,HashStatistics<TVHASH_T> & hs);
  template<typename TVHASH_T>
  void mer_bdbg_helper1(int argc, char ** argv,HashStatistics<TVHASH_T> & hs);
  template<typename TVHASH_T>
  void mer_diff_helper1(int argc, char ** argv);

private:
  void usage();

  void merCreateHashStats(int argc, char ** argv);
  void merFilterHashStats(int argc, char ** argv);
  void merInfoHashStats(int argc, char ** argv);
  void merSortHashStats(int argc, char ** argv);
  void merDumpHashStats(int argc, char ** argv);
  void merDumpDebug(int argc, char ** argv);
  void merDumpHashDistrib(int argc, char ** argv);
  void merDiffHashStats(int argc, char ** argv);
  void merBuildDBGHashStats(int argc, char ** argv);

  void merDeltaTest(int argc, char ** argv);

public:
  MiraMer() : MER_job("create") {};

  int mainMiraMer(int argc, char ** argv);
};
