/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ) and Bastien Chevreux
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

// 	$Id$

#ifndef _errorclasses_h_
#define _errorclasses_h_

#include <string>
#include <sstream>


#include "stdinc/defines.H"
#include "util/misc.H"





class Notify
{
public:
  enum {SILENT=0, REMARK, WARNING, FATAL, FOOL, EXTERNAL, INSTALL, INTERNAL};

public:
  std::string tif;               // thrown_in_function;
  std::string msg1;              // normally: Errordescription
  uint8  gravity;

  // variables
private:
  static bool NO_bangonthrow;

  //Functions
private:

  void foolCompiler();

  void zeroPointers();

  void discard();

public:
  Notify(uint8 gravity, const char * _tif);
  Notify(uint8 gravity, const char * _tif, const char * _msg1);

  Notify(Notify const &other);
  virtual ~Notify();

  friend std::ostream & operator<<(std::ostream &ostr, Notify const &n);
  Notify const & operator=(Notify const & other);

  void setGravity(uint8 _gravity);
  void setMsg(const char * msg) {msg1=msg;}

  void handleError(const char * _cif);     // caught in function
  void handleError(const std::string & _cif) {handleError(_cif.c_str());};     // caught in function

  static void setBangOnThrow(bool b) {NO_bangonthrow=b;}
};


// Flow is for non-fatal, non-usernotification flow control
// Its sole purpose is to be thrown in a function
//  and to be caught by the caller to show that there has been
//  a minor quirk.
class Flow
{
public:
  enum {UNDEFINED, POSSIBLE_DEAD_END, DEAD_END, TEMPLATE, PRMB_MISMATCH};
  int32 fl_flowtype;
  int32 fl_data;

private:

public:
  Flow(int32 flowtype=UNDEFINED);
  Flow(int32 flowtype, int32 data);

  friend std::ostream & operator<<(std::ostream &ostr, Flow const &f);

  void exitWhenUndefined();
};

#endif
