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

#ifndef _bas_machineinfo_h_
#define _bas_machineinfo_h_

#include <string>
#include <iostream>

#include "stdinc/defines.H"
#include "stdinc/types.H"

class MachineInfo
{
private:
  //Variables
  static bool MI_staticinitialised;

  static uint64 MI_memtotal;
  static uint64 MI_corestotal;


public:


  //Functions
private:
  static bool staticInit();

  static uint64 grepMemSizeFromProcFS(const char * fname, const std::string & identifier);
  static uint64 computeMemTotal();
  static uint64 computeMemAvail();
  static uint64 computeNumCores();

public:
  MachineInfo();

  MachineInfo const & operator=(MachineInfo const & other);
  friend std::ostream & operator<<(std::ostream &ostr, MachineInfo const &mi);

  inline static uint64 getCoresTotal() {return MI_corestotal;}
  inline static uint64 getMemTotal() {return MI_memtotal;}
  static uint64 getVMSize() {return grepMemSizeFromProcFS("/proc/self/status","VmSize:"); }
  inline static uint64 getMemAvail() {return computeMemAvail();}

};


#endif
