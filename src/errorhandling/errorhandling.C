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


#include "errorhandling.H"
#include "util/fmttext.H"

#include <cstdlib>     // exit(), free() under Cygwin
#include <csignal>
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;

bool Notify::NO_bangonthrow=false;

Notify::Notify(uint8 _gravity, const char * _tif)
{
  FUNCSTART("Notify::Notify(uint8 _gravity, const char * _tif)");
  zeroPointers();

  cerr.flush();
  gravity=_gravity;
  tif=_tif;

  if(NO_bangonthrow && _gravity>WARNING){
    cout << *this;
    raise(SIGTRAP);
  }

  FUNCEND();
}

Notify::Notify(uint8 _gravity, const char * _tif, const char * _msg1)
{
  FUNCSTART("Notify::Notify(uint8 _gravity, const char * _tif, const char * _msg1)");

  zeroPointers();

  cerr.flush();
  gravity=_gravity;
  tif=_tif;
  msg1= _msg1;

  if(NO_bangonthrow && _gravity>WARNING){
    cout << *this;
    raise(SIGTRAP);
  }

  FUNCEND();
}


Notify::Notify(Notify const &other)
{
  FUNCSTART("Notify::Notify(Notify const &other)");

  zeroPointers();

  *this=other;                               // call the copy operator

  FUNCEND();
}


// Copy operator, needed by copy-constructor
Notify const & Notify::operator=(Notify const & other)
{
  FUNCSTART("Notify const & Contig::operator=(Notify const & other)");

  if(this != &other){
    discard();

    tif=other.tif;
    msg1=other.msg1;
    gravity=other.gravity;
  }

  FUNCEND();
  return *this;
}


void Notify::zeroPointers()
{
  FUNCSTART("void Notify::zeroPointers()");

  gravity=FATAL;
  tif.clear();
  msg1.clear();

  FUNCEND();
}

Notify::~Notify()
{
  FUNCSTART("Notify::~Notify()");

  discard();

  FUNCEND();
}

void Notify::discard()
{
  FUNCSTART("Notify::discard()");

  zeroPointers();
  FUNCEND();
}

void Notify::setGravity(uint8 _gravity)
{
  if(gravity!=INTERNAL)  gravity=_gravity;
}

std::ostream & operator<<(std::ostream &ostr, Notify const &n)
{
  // flush stdout and stderr, there might be still some thing in it
  ostr.flush();
  cerr.flush();

  ostr << "\n";
  switch(n.gravity){
  case Notify::INTERNAL: {
    ostr << FmtText::wordWrap("Internal logic/programming/debugging error (*sigh* this should not have happened)\n\n",80);
    break;
  }
  case Notify::FATAL: {
    ostr << "Fatal error (may be due to problems of the input data or parameters):\n\n";
    break;
  }
  case Notify::WARNING: {
    ostr << "Warning:\n\n";
    break;
  }
  case Notify::REMARK: {
    ostr << "Remark:\n\n";
    break;
  }
  }
  if(!n.msg1.empty()){
    ostr << FmtText::makeTextSign(n.msg1,80);
  }
  if(!n.tif.empty()){
    ostr << "\n->Thrown: " << n.tif;
  }
  ostr << endl;

  return ostr;
}

void Notify::handleError(const char * _cif)
{
  FUNCSTART("void Notify::handleError(const char * _cif)");

  if(gravity!=SILENT) {

    cout << *this;
    if(_cif!=nullptr){
      cout << "->Caught: " << _cif << endl << endl;
    }

    std::string dueto;
    if(gravity==INTERNAL) {
      cout << FmtText::wordWrap("Aborting process, probably due to an internal error.\n\n"
				"If you want to report the error, please do so on\n"
				"\thttp://sourceforge.net/p/mira-assembler/tickets/\n"
				"and also give a short notice on the mira talk mailing list.\n\n"
				"If reporting, please do not delete the log and checkpoint directories, there may be files in them which could be needed to find the problem.",80) << endl;
    }else if(gravity==FATAL) {
      dueto="probably due to an error in the input data or parameters in the manifest file.";
    }else if(gravity==FOOL) {
      dueto="to prevent you from doing something considered foolish.";
    }else if(gravity==INSTALL) {
      dueto="due to an error in the installation of MIRA.";
    }else if(gravity==EXTERNAL) {
      dueto="due to an external cause.";
    }
    if(!dueto.empty()){
      std::string wws("Aborting process "+dueto+" Please check the lines above in this ouput for more information.\n");
      cout << FmtText::wordWrap(wws,80);
    }
    if(gravity==FATAL || gravity==INTERNAL || gravity==FOOL) {
      cout << FmtText::wordWrap("\nFor information on subscribing or unsubscribing to mira talk, see:\n  http://www.freelists.org/list/mira_talk\n",80) << endl;
	char *cwd;
	if ((cwd = getcwd(nullptr, 4096)) == nullptr) {
	    perror("pwd error");
	    exit(2);
	}
	cout << "CWD: " << cwd << endl;
	free(cwd);

	cout << FmtText::wordWrap("\nThank you for noticing that this is *NOT* a crash, but a controlled program stop.\n",80);

	exit(100);
    }
  }

  FUNCEND();
}



//Flow::Flow()
//{
//}
Flow::Flow(int32 flowtype)
{
  fl_flowtype=flowtype;
  fl_data=0;
}

Flow::Flow(int32 flowtype, int32 data)
{
  fl_flowtype=flowtype;
  fl_data=data;
}

std::ostream & operator<<(std::ostream &ostr, Flow const  &f)
{
  ostr << "Flowtype: ";
  switch(f.fl_flowtype){
  case Flow::UNDEFINED:{
    ostr << "undefined";
    break;
  }
  case Flow::POSSIBLE_DEAD_END:{
    ostr << "possible dead end ahead";
    break;
  }
  case Flow::DEAD_END:{
    ostr << "encountered dead end";
    break;
  }
  default:{
    ostr << "Unknown to me.";
  }
  }

  ostr << "\nData: " << f.fl_data << endl;

  return ostr;
}

void Flow::exitWhenUndefined()
{
  FUNCSTART("void Flow::exitWhenUndefined()");

  if(fl_flowtype==Flow::UNDEFINED){
    throw Notify(Notify::INTERNAL, THISFUNC, "Object not valid.");
  }

  FUNCEND();
  return;
}
