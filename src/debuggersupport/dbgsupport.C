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
 */

#include "dbgsupport.H"

#include <cstdio>

bool mira___seendebugger=miraDetectDebugger();

bool seenDebugger()
{
  return mira___seendebugger;
}

void miraDetectGDB()
{
  mira___seendebugger=false;

  // fileno() is not std C, does not work, e.g., on Windows
#if _POSIX_C_SOURCE >= 1 || _XOPEN_SOURCE || _POSIX_SOURCE
  // take /bin instead of /tmp, some machines apparently have access to /tmp denied
  FILE *fd = fopen("/bin", "r");
  // careful, the above may fail! (bug report from Chris Hoefler)
  if(fd!=nullptr){
    if(fileno(fd) > 5){
      mira___seendebugger=true;
    }
    fclose(fd);
  }
#endif

}

bool miraDetectDebugger()
{
  miraDetectGDB();
  return mira___seendebugger;
}
