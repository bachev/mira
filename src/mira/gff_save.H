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

#ifndef _bas_gffsave_h_
#define _bas_gffsave_h_

#include <iostream>

#include "stdinc/defines.H"

class Read;
class ReadPool;
class Contig;

class GFFSave
{
  //Variables
private:
  std::string   GFFS_foutname;
  std::ofstream GFFS_foutp;      // file out for padded

  std::vector<std::string> GFFS_gff3_strainnames;
  std::vector<std::ofstream *> GFFS_fouts;      // file out for unpadded strains, size of strainnames, one for each strain

  // saved data for deferred save in GFF3 format
  std::list<std::string> GFFS_gff3_names;           // contig names
  std::list<std::vector<std::string>> GFFS_gff3_paddedcons;      // padded contig cons, one list-element per contig, and per contig one sequence for each strain


public:


  //Functions
private:
  void foolCompiler();
  void zeroVars();

public:
  GFFSave();
  GFFSave(GFFSave const &other);
  ~GFFSave();

  GFFSave const & operator=(GFFSave const & other);
  friend std::ostream & operator<<(std::ostream &ostr, GFFSave const &gffs);

  void discard();

  void open(const char * filename);
  void open(const std::string & filename) {open(filename.c_str());}
  inline bool is_open() const {return GFFS_foutp.is_open();};
  void close();

  void acquireContig(Contig & con, const ReadPool & rp);
  void acquireRead(Read & read);
};


#endif
