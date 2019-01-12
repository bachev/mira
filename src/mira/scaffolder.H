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

#ifndef _scaffolder_h_
#define _scaffolder_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include <unordered_map>

// fwd declaration of Contig so that contig.H does not need to be included
//  by this .H
class Contig;

class Scaffolder
{
  //struct
private:
  struct contiginfo_t {
    std::string name;    // contigname;
  };

  struct placedcontig_t {
    size_t contigid=0;
    int8 dir=0;
    int32 disttonextcontig=0;
  };

  struct scaffold_t {
    std::vector<placedcontig_t> placedcontigs;
  };

  struct readlinkinfo_t {
    size_t contigid=0;
    size_t templateid=0;       // template id as created via MAF/CAF reading function
    int32  distance=0;     // distance from border of contig
    int8   dir=0;
    int8   wantdir=0;

    friend std::ostream & operator<<(std::ostream &ostr, readlinkinfo_t const &rli){
      ostr << "cid: " << rli.contigid << "\ttid: " << rli.templateid
	   << "\td: " << rli.distance
	   << "\tdr: " << static_cast<int16>(rli.dir)
	   << "\twdr: " << static_cast<int16>(rli.wantdir)
	;
      return ostr;
    }
  };

  //Variables
private:
  std::vector<contiginfo_t> SCA_contiginfo;     // index == contigid
  std::vector<scaffold_t>   SCA_scaffolds;      // index == scaffoldid;
  std::vector<readlinkinfo_t> SCA_readlinks;

  std::unordered_map<std::string, uint32> SCA_tname2tid;

public:


  //Functions
private:
  void foolCompiler();
  void init();
  void zeroVars();

public:
  Scaffolder();
  ~Scaffolder();

  Scaffolder(Scaffolder const &other) = delete;
  Scaffolder const & operator=(Scaffolder const & other) = delete;

  friend std::ostream & operator<<(std::ostream &ostr, Scaffolder const &sca);

  void discard();

  void storeInfoFreshContig(Contig & con);

  void dumpDebug();
};


#endif
