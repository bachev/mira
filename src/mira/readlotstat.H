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

#ifndef _bas_readlotstat_h_
#define _bas_readlotstat_h_

#include <iostream>
#include <string>

#include "stdinc/defines.H"


class ReadLotStat
{
  //Variables
private:

public:
  uint64 RLS_count_all=0;
  uint64 RLS_count_withoutclips=0;
  uint64 RLS_count_used=0;
  uint64 RLS_count_noqual=0;
  uint64 RLS_count_clipleft=0;
  uint64 RLS_count_clipright=0;
  uint64 RLS_count_chimera=0;

  uint64 RLS_len_all=0;
  uint64 RLS_len_used=0;

  uint64 getAvgLenAll() const { return RLS_count_all==0 ? 0 : RLS_len_all/RLS_count_all; }
  uint64 getAvgLenUsed() const { return RLS_count_used==0 ? 0 : RLS_len_used/RLS_count_used; }

  std::string RLS_verdict_chimera;

  //ReadLotStat();
  //ReadLotStat(ReadLotStat const &other);
  //~ReadLotStat();

  friend std::ostream & operator<<(std::ostream &ostr, ReadLotStat const &daq);
};


#endif
