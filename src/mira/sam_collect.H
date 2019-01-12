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
 *
 */


#ifndef _samcollect_h
#define _samcollect_h


#include <string>
#include <vector>

#include <unordered_map>
#include <unordered_set>

#include "stdinc/types.H"
#include "stdinc/stlincludes.H"

#include "mira/readgrouplib.H"


class SAMCollect
{
public:
  struct samrinfo_t {
    size_t initial_samri_id;     // 0.. size of SAMC_samris -1
    size_t contigid;     // in which contig is the read
    size_t templateid;   // part of which template

    int32 clippedstartpos;
    int32 clippedlen;
    int8  dir;
    uint8  template_segment;
    ReadGroupLib::ReadGroupID rgid;

    uint32 samflags;

    int32 rnext_conid;
    int32 pnext;
    int32 tlen;

    samrinfo_t(bool dummy) {(void) dummy;};  // lazy initialiser, we know what we do
    samrinfo_t(size_t sid, size_t cid, size_t tid, uint8 ts, ReadGroupLib::ReadGroupID rg, int32 csp, int32 cl, int8 di, uint32 sf) :
      initial_samri_id(sid),
      contigid(cid), templateid(tid),
      clippedstartpos(csp), clippedlen(cl),
      dir(di),
      template_segment(ts),
      rgid(rg),
      samflags(sf),
      rnext_conid(-1),
      pnext(0),
      tlen(0)
      {};
    // sort ascending by template, then by template-segment
    static inline bool lt_templateid(const samrinfo_t & a, const samrinfo_t & b) {
      if(a.templateid == b.templateid){
	return a.template_segment < b.template_segment;
      }
      return a.templateid < b.templateid;};

    friend std::ostream & operator<<(std::ostream &ostr, const samrinfo_t & si){
      ostr << "SI: "
	   << "isid: " << si.initial_samri_id
	   << "\tcid: " << si.contigid
	   << "\ttid: " << si.templateid
	   << "\tcsp: " << si.clippedstartpos
	   << "\tcln: " << si.clippedlen
	   << "\tdir: " << static_cast<int16>(si.dir)
	   << "\tts: " << static_cast<uint16>(si.template_segment)
	   << "\trgid: " << si.rgid.getLibId()
	   << "\tsfl: " << si.samflags
	   << "\trnxt_id: " << si.rnext_conid
	   << "\tpnxt: " << si.pnext
	   << "\ttlen: " << si.tlen
	   << "\n";
      return ostr;
    }
  };

  std::vector<samrinfo_t> SAMC_samris;

  std::unordered_set<std::string> SAMC_namecheck_contig;

  std::vector<std::string> SAMC_contignames;
  std::vector<size_t> SAMC_contiglengths;
  std::vector<std::string> SAMC_templatenames;
  std::unordered_map<std::string, size_t> SAMC_tname2tid;
  std::unordered_map<std::string, size_t> SAMC_rname2samriid;

  std::string SAMC_headerstring;

private:
  void errorMsgMAFFormat(const std::string & filename,
			 size_t linenumber,
			 const std::string & line,
			 const char * msg);

  void collectInfoFromMAF(const std::string & mafname);
  void processSAMRIs();

public:
  void processMAF(const std::string & mafname);
  void createSAMHeader();

  bool getSAMRInfo(const std::string & readname, samrinfo_t & samri) const;
  const std::string & getContigName(samrinfo_t & samri) const;
  const std::string & getRNextEntry(samrinfo_t & samri) const;
};





#endif
