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
#include "mira/dataprocessing.H"


class MiraSQT
{
private:
  enum { PS_NOPAIR=0, PS_INTERLEAVE, PS_2FILES};
  struct wqueueunit_t {
    uint8 wqu_status; // loading, baiting, saving etc.

    uint8 pairstatus; // 0==no pair; 1==interleave; 2==paired2files

    ReadPool rp1;
    ReadPool rp2;
    //std::vector<uint8> take1;
    //std::vector<uint8> take2;
  };

  struct files_t {
    std::string infilename1;
    std::string infilename2;

    uint8 intype1; // Read::AS_FASTQ etc.
    uint8 intype2; // Read::AS_FASTQ etc.
    uint8 writetype; // Read::AS_FASTQ etc.

    std::ofstream tfout1;
    std::ofstream tfout2;
  };

  // TODO: get the statics out into class variables

  static std::unordered_map<std::string,uint8> MS_fromtypemap;
  static std::unordered_map<std::string,uint8> MS_totypemap;

  static std::list<wqueueunit_t> MS_workqueue;
  static files_t MS_files;

  static std::vector<MIRAParameters> MS_Pv;
  static DataProcessing MS_dataprocessing;
  static std::vector<std::unique_ptr<DataProcessing>> MS_dpv;

  //static std::string MS_fromtype;
  //static std::string MS_baitfromtype;
  //static std::string MS_totype;

  static uint8 MS_tortype;

  static std::list<std::string> MS_infiles;

  static std::string MS_tmpdirname; // for hashstat tmp files

  static std::string MS_hashstatfname;

  static bool MS_testflag;

  static bool MS_fwdandrev;
  static bool MS_changeseqcase;
  static bool MS_dustfilter;

  static int16 MS_dustperc;

  static std::string MS_hitpath;
  static std::string MS_nameprefix;

  static uint32 MS_optthreads;
  static int32  MS_optmbtouse;

  static bool MS_mustdeletetargetfiles;

  static HashStatistics<vhash512_t> MS_hashstatistics;

  static std::list<char> MS_filepairinfo;    // p = 2 files, P = 1 file interleave

  static uint32 MS_basesperhash;
  static uint64 MS_numreadsread;

  static bool   MS_signal_ctrlc;

private:
  static void ctrlCHandler(int s);
  static void usage();
  static uint8 checkFromType(std::string & fromtype);
  static uint8 checkToType(std::string & totype);

  static void saveReadPool(ReadPool & rp);

  static uint8 setupRPIO(std::string filename, ReadGroupLib::ReadGroupID rgid, ReadPoolIO & rpio, uint8 & ziptype);

  static void setFakeFlags(ReadPool & rp);

  static void parallelSimpleClip(ReadPool & rp);
  template<typename TVHASH_T>
  static void doTreatWithHS(HashStatistics<TVHASH_T> & hs);

  template<typename TVHASH_T>
  static void parallelProposedEndClip(HashStatistics<TVHASH_T> & mbhs, ReadPool & rp);
  template<typename TVHASH_T>
  static void pecSingleRead(HashStatistics<TVHASH_T> & mbhs, Read & actread);

  static void saveWQueueElement(wqueueunit_t & wqu);

  static void setupOutfiles(const std::string & fname, uint8 rtype, uint8 ziptype, std::ofstream & hitfout);

public:
  ~MiraSQT();

  int mainMiraSQT(int argc, char ** argv);

};
