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

#include "util/fileanddisk.H"

#include "mira/parameters.H"
#include "mira/readpool_io.H"
#include "mira/hashstats.H"
#include "mira/vhash.H"


class MiraBait
{
private:
  enum { PS_NOPAIR=0, PS_INTERLEAVE, PS_2FILES};
  struct wqueueunit_t {
    uint8 wqu_status; // loading, baiting, saving etc.

    uint8 pairstatus; // 0==no pair; 1==interleave; 2==paired2files

    ReadPool rp1;
    ReadPool rp2;
    std::vector<uint8> take1;
    std::vector<uint8> take2;
  };

  struct files_t {
    std::string infilename1;
    std::string infilename2;

    uint8 intype1; // Read::AS_FASTQ etc.
    uint8 intype2; // Read::AS_FASTQ etc.
    uint8 writetype; // Read::AS_FASTQ etc.

    std::ofstream hitfout1;
    std::ofstream missfout1;
    std::ofstream hitfout2;
    std::ofstream missfout2;
  };

  // TODO: get the statics out into class variables

  static std::unordered_map<std::string,uint8> MB_fromtypemap;
  static std::unordered_map<std::string,uint8> MB_totypemap;

  static std::list<wqueueunit_t> MB_workqueue;
  static files_t MB_files;

  static std::vector<MIRAParameters> MB_Pv;

  //static std::string MB_fromtype;
  //static std::string MB_baitfromtype;
  //static std::string MB_totype;

  static uint8 MB_tortype;

  static std::list<std::string> MB_baitfiles;
  static std::list<std::string> MB_infiles;

  static std::string MB_tmpdirname; // for hashstat tmp files

  static std::string MB_hashstatfname;
  static std::string MB_savehashstatfname;

  static bool MB_testflag;

  static bool MB_mergeoutput;

  static bool MB_wantbaithits;
  static bool MB_wantbaitmiss;
  static bool MB_fwdandrev;
  static bool MB_dustfilter;

  static bool MB_changeseqcase;
  static char MB_maskchar;

  static int16 MB_dustperc;

  static std::string MB_hitpath;
  static std::string MB_misspath;
  static std::string MB_nameprefix;

  static int32 MB_numbaithits;

  static uint32 MB_optthreads;
  static int32  MB_optmbtouse;

  static bool MB_mustdeletetargetfiles;

  static HashStatistics<vhash512_t> MB_hashstatistics;

  static std::list<char> MB_filepairinfo;    // p = 2 files, P = 1 file interleave

  static uint32 MB_basesperhash;
  static uint64 MB_baitpoolsize;
  static uint64 MB_numreadsread;
  static uint64 MB_numpairsbaited;
  static uint64 MB_numpairsmissed;
  static uint64 MB_numunpairedbaited;
  static uint64 MB_numunpairedmissed;

  static uint64 MB_numclippedreadsinload; // number of reads which have clips already when loaded (CAF/MAF)

  static bool   MB_signal_ctrlc;

private:
  static void ctrlCHandler(int s);
  static void usage();
  static uint8 checkFromType(std::string & fromtype);
  static uint8 checkToType(std::string & totype);

  static void saveReadPool(ReadPool & rp);

  static uint8 setupRPIO(std::string filename, ReadGroupLib::ReadGroupID rgid, ReadPoolIO & rpio, uint8 & ziptype);

  template<typename TVHASH_T>
  static void baitReads(HashStatistics<TVHASH_T> & hs, const ReadPool & rp, std::vector<uint8> & take);
  template<typename TVHASH_T>
  static void parallelBaitReads(HashStatistics<TVHASH_T> & hs, const ReadPool & rp, std::vector<uint8> & take);
  template<typename TVHASH_T>
  static void doBaitWithHS(HashStatistics<TVHASH_T> & hs);

  static void saveWQueueElement(wqueueunit_t & wqu);

  static void setupOutfiles(const std::string & fname, uint8 rtype, uint8 ziptype, std::ofstream & hitfout, std::ofstream & missfout);

public:
  ~MiraBait();

  int mainMiraBait(int argc, char ** argv);

};
