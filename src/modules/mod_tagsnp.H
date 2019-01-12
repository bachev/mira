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

#include <iostream>
#include <fstream>

#include <list>

#include "mira/readpool.H"
#include "mira/contig.H"


class tagsnp
{
private:
  static std::ofstream TS_fout;

  static std::list<Contig> TS_clist;   // needed for CAF conversion (and GBF)

private:
  void usage();


  static void doit2(std::list<Contig> & contigs);
  static void saveCList(std::list<Contig> & clist, ReadPool & rp);
  static void cafload_callback(std::list<Contig> & clist, ReadPool & rp);

public:
  int mainTagSNP(int argc, char ** argv);
};

