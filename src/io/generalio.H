/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2002 and later by Bastien Chevreux
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

#ifndef _bas_gio_h_
#define _bas_gio_h_

#include <iostream>
#include <fstream>

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"


#if __GNUC__ >= 3
#include <sstream>
#include <istream>

#else

#include <strstream>
#endif

#include <string>

#include "errorhandling/errorhandling.H"


namespace GeneralIO
{
  bool readKeyValue(std::istream & sin, std::string & key, std::string & value);
}


#endif
