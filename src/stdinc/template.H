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

#ifndef _bas_???_h_
#define _bas_???_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"



class CLASS
{
  //Variables
private:

public:


  //Functions
private:
  void init();
  void zeroVars();

public:
  CLASS();
  CLASS(CLASS const &other);
  ~CLASS();

  CLASS const & operator=(CLASS const & other);
  friend std::ostream & operator<<(std::ostream &ostr, CLASS const &???);

  void discard();
};


#endif
