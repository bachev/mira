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

#include "machineinfo.H"

#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach_host.h>
#endif

#include <fstream>

#include <set>
#include <regex>


uint64 MachineInfo::MI_memtotal=0;
uint64 MachineInfo::MI_corestotal=0;

bool   MachineInfo::MI_staticinitialised=staticInit();


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool MachineInfo::staticInit()
{
  MI_memtotal=computeMemTotal();
  MI_corestotal=computeNumCores();
  return true;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


uint64 MachineInfo::computeNumCores()
{
#ifdef __APPLE__
  size_t uints=sizeof(uint64);
  uint64 numc=0;
  sysctlbyname("hw.physicalcpu", &numc, &uints, nullptr, 0);
  if(numc<1) numc=1;

  return numc;
#else

  static const std::regex re_processor("^processor[\\s]*:[\\s]*([0-9]+)");
  static const std::regex re_physicalid("^physical id[\\s]*:[\\s]*([0-9]+)");
  static const std::regex re_coreid("^core id[\\s]*:[\\s]*([0-9]+)");

  std::set<uint32> pidcidset; // should be enough for up to 2^16 CPUs with 2^16 cores each
  std::string tmpline;

  std::ifstream fin("/proc/cpuinfo", std::ios::in);

  if(fin){
    int16 physicalid=-1;
    int16 coreid=-1;
    std::smatch what;
    while(getline(fin,tmpline)){
      if(std::regex_search(tmpline.cbegin(), tmpline.cend(),
			   what, re_processor, std::regex_constants::match_default)){
	if(physicalid>=0 && coreid>=0){
	  uint32 pidcid=physicalid;
	  pidcid<<=16;
	  pidcid+=coreid;
	  pidcidset.insert(pidcid);
	  physicalid=-1;
	  coreid=-1;
	}
      }
      if(std::regex_search(tmpline.cbegin(), tmpline.cend(),
			   what, re_physicalid, std::regex_constants::match_default)){
	if(what.size()>=2){
	  std::string tmps(what[1].first, what[1].second);
	  physicalid=static_cast<int16>(atoi(tmps.c_str()));
	}
      }
      if(std::regex_search(tmpline.cbegin(), tmpline.cend(),
			   what, re_coreid, std::regex_constants::match_default)){
	if(what.size()>=2){
	  std::string tmps(what[1].first, what[1].second);
	  coreid=static_cast<int16>(atoi(tmps.c_str()));
	}
      }
    }
  }

  uint64 retval=1;
  if(!pidcidset.empty()) retval=pidcidset.size();
  return retval;
#endif
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint64 MachineInfo::computeMemTotal()
{
#ifdef __APPLE__
  size_t uints=sizeof(uint64);
  uint64 mem=0;
  //sysctlbyname("hw.memsize", &mem, &uints, nullptr, 0);
  sysctlbyname("hw.usermem", &mem, &uints, nullptr, 0);

  return mem;
#else
  return grepMemSizeFromProcFS("/proc/meminfo","MemTotal:");
#endif
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint64 MachineInfo::computeMemAvail()
{
#ifdef __APPLE__
  uint64 retval=0;

  // from http://stackoverflow.com/questions/6094444/how-can-i-programmatically-check-free-system-memory-on-mac-like-the-activity-mon
  int mib[6];
  mib[0] = CTL_HW;
  mib[1] = HW_PAGESIZE;

  int pagesize;
  size_t length;
  length = sizeof (pagesize);
  if (sysctl (mib, 2, &pagesize, &length, nullptr, 0) < 0){
    fprintf (stderr, "getting page size");
  }

  mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
  vm_statistics_data_t vmstat;
  if(KERN_SUCCESS == host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vmstat, &count)){
    retval=vmstat.free_count + vmstat.inactive_count;
  }

  return retval*pagesize;
#else
  return grepMemSizeFromProcFS("/proc/meminfo","MemFree:")+grepMemSizeFromProcFS("/proc/meminfo","Cached:");
#endif
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint64 MachineInfo::grepMemSizeFromProcFS(const char * fname, const std::string & identifier)
{
  uint64 retval=0;

  std::ifstream fin(fname, std::ios::in);
  std::string dummy;
  if(fin){
    while(!fin.eof()){
      fin >> dummy;
      if(dummy==identifier){
	fin >> retval;
	fin >> dummy;
	if(dummy=="kB" || dummy=="kiB"){
	  retval*=1024;
	}else if(dummy=="mB" || dummy=="miB"){
	  retval*=1024*1024;
	}else if(dummy=="gB" || dummy=="giB"){
	  retval*=1024*1024*1024;
	}
	break;
      }
      if(fin.eof()) break;
    }
    fin.close();
  }
  return retval;
}
