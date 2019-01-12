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

#ifndef _bas_dataprocessing_h_
#define _bas_dataprocessing_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "mira/hashstats.H"
#include "mira/multitag.H"
#include "mira/vhash.H"

class MIRAParameters;
class ReadPool;
template<typename TVHASH_T> class Skim;

#include <regex>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>



// DataProcessing
// One instatiated object should be used in one thread only, so do not share a DataProcessing
//  object between several threads unless access is serialised by the threads
// But each thread having it's own object is fine ... and meant to be that way.

class DataProcessing
{
public:
  struct dpstats_t {
    uint32 cphix174;
    uint32 crrna;
    uint32 cadapright;
    uint32 cadaprightpartial;

    dpstats_t() : cphix174(0), crrna(0), cadapright(0), cadaprightpartial(0) {};
  };

  dpstats_t DP_stats;

  //Variables
private:
  struct poolskim_t {
    ReadPool * poolptr;
    Skim<vhash64_t> * skimptr;

    poolskim_t() : poolptr(nullptr), skimptr(nullptr) {};
    ~poolskim_t();
  };

  static std::vector<poolskim_t> DP_adapskims;
  static boost::mutex DP_ps_changemutex;         // exclusive mutex for write access to DP_*skims, we need that as it's a static variable

  static HashStatistics<vhash64_t> DP_phix174hashstatistics;
  static bool DP_px174hs_init;
  static boost::mutex DP_px174hs_changemutex;         // exclusive mutex for write access to hashstatistics, we need that as it's a static variable

  static HashStatistics<vhash64_t> DP_rrnahashstatistics;
  static bool DP_rrnahs_init;
  static boost::mutex DP_rrnahs_changemutex;         // exclusive mutex for write access to hashstatistics, we need that as it's a static variable

  // while Skim is thread safe in the adaptor search, the boost RegEx are not when using the same objects
  // therefore, keep those per DataProcessing object and not static
  // no mutex here as one DataProcessing object is meant to be used by one thread only
  struct masterslavere_t {
    std::regex masterre;
    std::list<std::regex> slaveres;
    bool hasmaster;

    masterslavere_t(): hasmaster(false) {};
  };
  struct adapres_t {
    std::list<masterslavere_t> adapres;
    bool areinit;

    adapres_t() : areinit(false) {}
  };
  std::vector<adapres_t> DP_adapres;



  std::vector<MIRAParameters> * DP_miraparams_ptr;
  int32 DP_threadid; // when using multi-threaded versions, this is the thread-id this object uses (-1==non-multi-threaded)

  std::ofstream DP_logfout;
  std::string DP_logname;

  static std::string DP_ggcstring;

  struct threadsharecontrol_t {
    boost::mutex accessmutex;
    int32 from;
    int32 to;
    int32 todo;
    int32 done;
    int32 stepping;
  };

  // For multithreaded baiting within a single HashStatistics object (e.g. Phi X 174 search)

  //typename std::vector<typename HashStatistics<vhash64_t>::vhrap_t> DP_baiting_singlereadvhraparray;
  //std::vector<uint8>   DP_baiting_tagmaskvector;


  //
  // TEMPS
  //

  // uint8 vector for temporary use
  std::vector<uint8> DP_tmpvu8;

  // uint32 array with 128 elements for temporary use
  uint32 DP_tmpau32_128[128];

  // multitag with SOFA poly AT signature for temporary use
  multitag_t DP_tmpmtpolyAT;


public:


  //functions
private:
  void init();
  void zeroVars();

  static bool priv_staticInitialiser();

  void priv_EnsureAdapRegexes(ReadGroupLib::ReadGroupID rgid);
  void addAdapRegexes(ReadGroupLib::ReadGroupID rgid, const char * regexfile);

  void priv_EnsureAdapSkims(ReadGroupLib::ReadGroupID rgid);
  void priv_EnsurePhiX174Statistics();
  void priv_EnsureRRNAStatistics();
  void priv_constructorSkimPool(ReadGroupLib::ReadGroupID rgid, std::vector<poolskim_t> & skimpool, const uint32 basesperhash, const char * adapfile);

  void priv_buntifyHelper(uint8 allowedfreq,
			  uint32 basesperhash,
			  std::vector<Read::bposhashstat_t>::const_iterator bhsI,
			  std::vector<Read::bposhashstat_t>::const_iterator bhsE,
			  std::vector<uint8>::iterator tfI,
			  std::vector<uint8>::iterator tfE);

  bool searchPolyBaseFrom5Prime(Read & actread,
				const char polybase,
				int32 & lpolystart,
				int32 & rpolyend,
				const uint32 mincount,
				const uint32 maxbad,
				int32 grace);
  bool searchPolyBaseFrom3Prime(Read & actread,
				const char polybase,
				int32 & lpolystart,
				int32 & rpolyend,
				const uint32 mincount,
				const uint32 maxbad,
				int32 grace);
  static void priv_stdTreatmentThread(uint32 threadnum,
				      threadsharecontrol_t * tscptr,
				      std::vector<MIRAParameters> * mpptr,
				      DataProcessing * dpptr,
				      ReadPool * rpoolptr,
				      std::vector<uint8> * debrisreasonptr,
				      std::string * logprefixptr);

  static void priv_stp_helperDebris(std::vector<MIRAParameters> & mp,
				    ReadPool & rpool,
				    Read & actread,
				    int32 rid,
				    std::vector<uint8> * debrisreasonptr,
				    uint8 reason);

  void pechelper_fwd(Read & actread, uint32 basesperhash);
  void pechelper_rev(Read & actread, uint32 basesperhash);

public:
  DataProcessing(std::vector<MIRAParameters> * params);
  ~DataProcessing();

  DataProcessing(DataProcessing const &other)=delete;
  DataProcessing const & operator=(DataProcessing const & other)=delete;

  friend std::ostream & operator<<(std::ostream &ostr, DataProcessing const &dp);

  void discard();

  void setThreadID(int32 id) {DP_threadid=id;};

  void startLogging(const std::string & filename, bool newfile);
  void stopLogging();

  // Pool versions
  // Take care of everything from MIRA parameters
  //
  // Read versions
  // Perform the task unconditionally on the reads

  void performRareKMERMasking_Pool(ReadPool & rpool, uint32 basesperhash, const std::string & logprefix);
  void performRareKMERMasking_Read(Read & actread, uint32 basesperhash, const std::string & logprefix);

  void performKMERRepeatTagging_Pool(ReadPool & rpool, uint32 basesperhash);
  void performKMERRepeatTagging_Read(Read & actread, uint32 basesperhash);

  void buntifyReadsByHashFreq_Pool(ReadPool & rpool, uint32 basesperhash);
  void buntifyReadsByHashFreq_Read(Read & actread, uint32 basesperhash);

  void addKMerForkTags_Pool(ReadPool & rp, uint32 basesperhash);

  void proposedEndClipping_Pool(ReadPool & rpool, uint32 basesperhash);
  void proposedEndClipping_Read(Read & actread, uint32 basesperhash);

  void clipBadSolexaEnds_Pool(ReadPool & rpool, const std::string & logprefix);
  void clipBadSolexaEnds_Read(Read & actread, const std::string & logprefix);

  void lowerCaseClipping_Pool(ReadPool & rpool, const std::string & logprefix);
  void lowerCaseClippingFront_Read(Read & actread, const std::string & logprefix);
  void lowerCaseClippingBack_Read(Read & actread, const std::string & logprefix);

  void qualClips_Pool(ReadPool & rpool, const std::string & logprefix);
  void qualClips_Read(Read & actread, const std::string & logprefix);

  void maskClips_Pool(ReadPool & rpool, const std::string & logprefix);
  void maskClips_Read(Read & actread, const std::string & logprefix, int32 gapsize=-1, int32 maxfrontgap=-1, int32 maxendgap=-1);

  void minimumQualityThreshold_Pool(ReadPool & rpool, const std::string & logprefix);
  bool minimumQualityThreshold_Read(Read & actread, const std::string & logprefix);

  void minimumLeftClip_Pool(ReadPool & rp, const std::string & logprefix,
			    bool qual, bool seqvec, bool mask);
  void minimumLeftClip_Read(Read & actread, const std::string & logprefix,
			    bool qual, bool seqvec, bool mask);

  void minimumRightClip_Pool(ReadPool & rp, const std::string & logprefix,
			     bool qual, bool seqvec, bool mask);
  void minimumRightClip_Read(Read & actread, const std::string & logprefix,
			     bool qual, bool seqvec, bool mask);

  void badSequenceSearch_Pool(ReadPool & rpool, const std::string & logprefix);
  void badSequenceSearch_Read(Read & actread, const std::string & logprefix);

  void clipPolyATAtEnds_Pool(ReadPool & rpool, const std::string & logprefix);
  void clipPolyATAtEnds_Read(Read & actread, const std::string & logprefix);

  void clipPolyBaseAtEnd_Pool(ReadPool & rpool, const std::string & logprefix);
  void clipPolyBaseAtEnd_Read(Read & actread, const std::string & logprefix);

  void adaptorRightClip_Pool(ReadPool & rpool, const std::string & logprefix);
  void adaptorRightClip_Read(Read & actread, const std::string & logprefix);

  static bool addPhiX174ToReadpool(ReadPool & baitrp);
  void seqMatchPhiX174_Read(Read & actread, const std::string & logprefix, bool filter);
  void seqMatchRRNA_Read(Read & actread, const std::string & logprefix, uint32 numkmers);


  uint32 markReadsWithInvalidKMerEndsAsChimeras_Pool(ReadPool & rp,
						     uint32 basesperhash,
						     bool killreads,
						     std::vector<uint8> * debrisreasonptr,
						     std::vector<uint8> & ischimera,
						     std::vector<uint8> & nochimerakill,
						     const std::string & logprefix);
  bool checkReadForInvalidKMerEndsAsChimera(Read & actread,
					    uint32 basesperhash);
  uint32 markReadsWithRareKMersAsChimeras_Pool(ReadPool & rp,
					       uint32 basesperhash,
					       bool killreads,
					       std::vector<uint8> * debrisreasonptr,
					       std::vector<uint8> & ischimera,
					       std::vector<uint8> & nochimerakill,
					       const std::string & logprefix);
  bool checkReadForRareKMersAsChimera(Read & actread,
				      uint32 basesperhash);

//  static void stdTreatmentPool_MultiThread(ReadPool & rpool);
  static void stdTreatmentPool_SingleThread(std::vector<MIRAParameters> & mp,
					    DataProcessing & dp,
					    ReadPool & rpool,
					    std::vector<uint8> * debrisreason,
					    std::string & logprefix,
					    bool oprogress,
					    int32 fromid=-1,
					    int32 toid=-1);
  static void stdTreatmentPool_MultiThread(std::vector<MIRAParameters> & mp,
					   DataProcessing & dpcollector,
					   std::vector<std::unique_ptr<DataProcessing>> & dpv,
					   ReadPool & rpool,
					   std::vector<uint8> * debrisreasonptr,
					   std::string & logprefix,
					   bool progress,
					   int32 fromid=-1,
					   int32 toid=-1);

  template<typename TVHASH_T>
  static void performDigitalNormalisation_Pool(ReadPool & rp,
					       HashStatistics<TVHASH_T> & hsd,
					       std::vector<uint8> * debrisreasonptr);

  template<typename TVHASH_T>
  uint32 performSDBGChimeraSearch_Pool(ReadPool & rp,
				       HashStatistics<TVHASH_T> & hsd,
				       uint32 trimfreq,
				       std::vector<uint8> * debrisreasonptr,
				       const std::string & logprefix);
  template<typename TVHASH_T>
  bool performSDBGChimeraSearch_Read(Read & actread, HashStatistics<TVHASH_T> & hsd, const std::string & logprefix);

  template<typename TVHASH_T>
  uint32 performSDBGEdits_Pool(ReadPool & rp,
			       HashStatistics<TVHASH_T> & hsd,
			       uint32 trimfreq);
};


#endif
