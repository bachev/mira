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


#ifndef _mira_std_defines_h_
#define _mira_std_defines_h_

#include <csignal>

#include "stdinc/types.H"
#include "debuggersupport/dbgsupport.H"


#if __GNUC__ >= 4
#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)
#else
#define likely(x) (x)
#define unlikely(x) (x)
#endif

//#define MIRA_HAS_EDIT
#define BITSET(n,f) ((f)|=(1<<(n)))
#define BITTEST(n,f) (((f)&(1<<(n)))!=0)

#if 0
#define TIMERESTRICTED
#endif

//#define BANGONTHROW




#ifdef __APPLE__
// somehow, the code produced by gcc on OSX does not catch all exceptions?!
// make sure the user gets to read the most important things
#ifndef COUTONTHROW
#define COUTONTHROW
#endif
#endif

#ifdef BANGONTHROW
#define BANGNOW {raise(SIGTRAP);}
#ifndef COUTONTHROW
#define COUTONTHROW
#endif
#else
#define BANGNOW
#endif

#ifdef COUTONTHROW
#define COUTERRMSG {std::cout << std::endl << "In func: " << THISFUNC << std::endl << "Throw message: " << my__emsg.str() << std::endl;}
#else
#define COUTERRMSG
#endif


#define ERROR(bla)  { std::cout.flush(); std::cerr << "\nFatal error in: " << THISFUNC << ":  " << bla; exit(1);}
#define WARNING(bla)  { std::cout.flush(); std::cerr << "Warning in: " << THISFUNC << ":  " << bla << endl;}

#define FUNCSTART(bla) static const char * THISFUNC = bla; (void) THISFUNC;
#define FUNCEND()

#define MIRANOTIFY(nlevel, throwmsg) {std::ostringstream my__emsg; my__emsg << throwmsg; COUTERRMSG; if(seenDebugger()) {raise(SIGTRAP);}; if(nlevel>2){BANGNOW;}; throw Notify(nlevel, THISFUNC, my__emsg.str().c_str());}



#ifdef BUGTRACKFLAG
#define BUGIFTHROW(ifcond, throwmsg) { if(unlikely(ifcond)) {std::ostringstream my__emsg; my__emsg << throwmsg; COUTERRMSG; if(seenDebugger()) {raise(SIGTRAP);};BANGNOW; throw Notify(Notify::INTERNAL, THISFUNC, my__emsg.str().c_str());}}
#define BUGIF(ifcond, statement) { if(ifcond) {statement;}}
#define BUGSTAT(statement) {statement;}
#else
#define BUGIFTHROW(ifcond, throwmsg)
#define BUGIF(ifcond, statement)
#define BUGSTAT(statement)
#endif

#ifdef BOUNDTRACKFLAG
//#define BOUNDCHECK(x, l, h) assert((x) >= (l) && (x) < (h));
#define BOUNDCHECK(x, l, h) if(unlikely(!((x) >= (l) && (x) < (h)))) MIRANOTIFY(Notify::INTERNAL, "Oooops, boundcheck failed: " << x << ' ' << l << ' ' << h );
#else
#define BOUNDCHECK(x, l, h)
#endif


#define MIRASRLEN 32760

#endif
