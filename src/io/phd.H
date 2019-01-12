/*
 * Written by Bastien Chevreux (BaCh)
 *
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


#ifndef _phd_h_
#define _phd_h_

#include <fstream>
#include <string>

#include "stdinc/stlincludes.H"

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"


class PHD
{
private:
  std::string         PHD_name;
  std::string         PHD_sequence;

  std::vector<base_quality_t>  PHD_qualities;
  std::vector<uint32>          PHD_peakindex;

private:
  void   zeroVars();

public:
  PHD();
  virtual ~PHD();
  friend std::ostream & operator<<(std::ostream &ostr, PHD const &thephd);
  PHD const & operator=(PHD const & other);

  void load(const char * filename);
  void loadNextSeq(std::ifstream & phdfin);
  bool testIfEmpty();

  void discard();

  inline const std::string   & getSequence() const { return PHD_sequence;}
  inline const std::string   & getName() const { return PHD_name;}
  inline const std::vector<base_quality_t>  & getQualities() const { return PHD_qualities;}
  inline const std::vector<uint32> & getPeakIndex() const { return PHD_peakindex;}

};

#endif
