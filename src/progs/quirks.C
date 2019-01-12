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

#include <boost/filesystem.hpp>

#include <iostream>

#include <cstdlib>

// make the "tcmalloc: large alloc" messages from TCMallom disappear
// by setting the reporting environment variable to a very large value
// see: http://groups.google.com/group/google-perftools/browse_thread/thread/24a003fc35f3d470?pli=1
void quietenTCMALLOC()
{
  setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD",
	 "1099511627776", // 1 TiB, should be enough to quieten almost everything
	 0); // "0" means "do not overwrite if already existing"
}

/*
  On some systems, Boost::filesystem (at least up to 1.50) throws a standard
  exception when using some functions: locale::facet::_S_create_c_locale name not valid
  Only countermeasure possible: setenv
  Using setlocale or std::locale::global does *NOT* work as workaround (tried and tested)
 */
void fixLocaleQuirk()
{
  try{
    // this must work
    boost::filesystem::path fp = boost::filesystem::current_path();
  }
  catch(...){
    // if not, we're on a system with quirks
    setenv("LC_ALL",
	   "C",
	   1);

    std::cout << "Your system seems to be older or have some quirks with locale settings."
      "\nUsing the LC_ALL=C workaround."
      "\nIf you don't want that, fix your system ;-)\n";
  }
}


void fixQuirks()
{
  quietenTCMALLOC();
  fixLocaleQuirk();
}

