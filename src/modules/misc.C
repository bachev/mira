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

#include "modules/misc.H"

#include "stdinc/defines.H"
#include "stdinc/types.H"
#include "io/generalio.H"
#include "errorhandling/errorhandling.H"

#include <iostream>
#include <fstream>

#include <cstdlib>      // exit() on Cygwin

using std::cout;
using std::cin;
using std::cerr;
using std::endl;


extern const char * compileinfo;


void dumpStdMsg()
{
  cout <<
    "To (un-)subscribe the MIRA mailing lists, see:\n"
    "\thttp://www.chevreux.org/mira_mailinglists.html\n\n"
    "After subscribing, mail general questions to the MIRA talk mailing list:\n"
    "\tmira_talk@freelists.org\n\n"
    "\nTo report issues or ask for features, please use the GitHub issue\nsystem at:\n"
    "\thttps://github.com/bachev/mira/issues\n"
    "This ensures that requests do not get lost.\n\n\n";

  bool addnl=false;

  cout << compileinfo;
#ifdef CEBUGFLAG
  cout << "Compiled in debug output mode.\n";
  addnl=true;
#endif
#ifdef TRACEFLAG
  cout << "Compiled with trace mode.\n";
  addnl=true;
#endif
#ifdef BOUNDTRACKFLAG
  cout << "Compiled in boundtracking mode.\n";
  addnl=true;
#endif
#ifdef BUGTRACKFLAG
  cout << "Compiled in bugtracking mode.\n";
  addnl=true;
#endif
#ifdef PARANOIABUGTRACKFLAG
  cout << "Compiled in paranoia bugtracking mode.\n";
  addnl=true;
#endif
#ifdef ENABLE64
  cout << "Compiled with ENABLE64 activated.\n";
  addnl=true;
#else
  cout << "Compiled with ENABLE64 de-activated.\n";
  addnl=true;
#endif
#ifdef MIRAMEMORC
  cout << "Compiled with memory overrun checks, MIRA *will* be slower.\n";
  addnl=true;
#endif

  cout << "Runtime settings (sorry, for debug):"
       << "\n\tSize of size_t  : " << sizeof(size_t)
       << "\n\tSize of uint32  : " << sizeof(uint32)
       << "\n\tSize of uint32_t: " << sizeof(uint32_t)
       << "\n\tSize of uint64  : " << sizeof(uint64)
       << "\n\tSize of uint64_t: " << sizeof(uint64_t)
       << "\nCurrent system: ";
  {
    cout.flush();
    int tmp=system("uname -a");
    // don't complain about unused variable
    (void) tmp;
  }

  if(addnl) cout << endl;
}


General::strintmap General::GE_nameselectionmap;
bool General::GE_namesread=false;

void General::makeSelectionStringSet(std::string & filename)
{
  FUNCSTART("void makeSelectionStringSet(std::string & filename)");

  std::ifstream fin(filename, std::ios::in);
  if(!fin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }
  fin.seekg(0, std::ios::beg);

  std::string elemname, dummy;
  strintmap::iterator nI;
  uint32 numread=0;
  while(GeneralIO::readKeyValue(fin, elemname,dummy)){
    nI=GE_nameselectionmap.find(elemname);
    if(nI==GE_nameselectionmap.end()) {
      GE_nameselectionmap[elemname]=numread;
      numread++;
    }
  }
  GE_namesread=true;

  if(numread>0 && GE_nameselectionmap.empty()) {
    cerr << "ehhh?";
    exit(10);
  }

  FUNCEND();
}



bool General::checkNamePresence(std::string & name)
{
  if(!GE_namesread) return true;
  return (GE_nameselectionmap.find(name) != GE_nameselectionmap.end());
}

bool General::hasNames()
{
  return GE_namesread;
}

size_t General::getNameOrder(const std::string & name)
{
  if(GE_nameselectionmap.empty()) return (0-1);
  strintmap::iterator nI=GE_nameselectionmap.find(name);
  if (nI == GE_nameselectionmap.end()) return (0-1);
  return nI->second;
}

