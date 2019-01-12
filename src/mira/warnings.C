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

#include "mira/warnings.H"

#include "errorhandling/errorhandling.H"
#include "util/fmttext.H"

using std::cout;
using std::cerr;
using std::endl;

// Plain vanilla constructor
Warnings::Warnings()
{
  FUNCSTART("Warnings::Warnings()");

  zeroVars();
  init();

  FUNCEND();
}

void Warnings::zeroVars()
{
  FUNCSTART("void Warnings::zeroVars()");
  FUNCEND();
}

void Warnings::init()
{
  FUNCSTART("void Warnings::init()");
  FUNCEND();
}



Warnings::~Warnings()
{
  FUNCSTART("Warnings::~Warnings()");

  discard();

  FUNCEND();
}


void Warnings::discard()
{
  FUNCSTART("Warnings::discard()");

  zeroVars();

  FUNCEND();
}

void Warnings::priv_setWarning(std::string & shortcode,uint32 level,std::string & title,std::string & message,bool donotdump,bool add)
{
  FUNCSTART("void Warnings::priv_setWarning(std::string & shortcode,uint32 level,std::string & title,std::string & message,bool add)");

  BUGIFTHROW(level>2,"level>2 ???");
  if(add){
    auto mI=WA_messages.find(shortcode);
    if(mI==WA_messages.end()){
      WA_messages[shortcode]=warnmsg_t(title,message,level);
    }else{
      mI->second.warnlevel=level;
      mI->second.title=title;
      mI->second.message+=message;
    }
  }else{
    WA_messages[shortcode]=warnmsg_t(title,message,level);
  }
  priv_dumpSingleWarning(*(WA_messages.find(shortcode)),true,cout);
}


bool Warnings::priv_clearWarning(std::string & shortcode)
{
  FUNCSTART("void Warnings::priv_clearWarning(std::string & shortcode)");

  auto num=WA_messages.erase(shortcode);
  if(num){
    cout << "Removed warning: " << shortcode << '\n';
  }
  return num!=0;
}


std::ostream & operator<<(std::ostream &ostr, Warnings const &war)
{
  FUNCSTART("friend std::ostream & Warnings::operator<<(std::ostream &ostr, const  &war)");

  for(uint32 wl=0; wl<3; ++wl){
    for(auto & wme : war.WA_messages){
      if(wme.second.warnlevel==wl){
	war.priv_dumpSingleWarning(wme,true,ostr);
      }
    }
  }

  FUNCEND();
  return ostr;
}

void Warnings::priv_dumpSingleWarning(const std::pair<const std::string, warnmsg_t> & wm, bool withheader, std::ostream & ostr) const
{
  if(withheader){
    ostr << "-------- ";
    if(wm.second.warnlevel==0){
      ostr << "CRITICAL";
    }else if(wm.second.warnlevel==1){
      ostr << "MEDIUM";
    }else if(wm.second.warnlevel==2){
      ostr << "MINOR";
    }
    ostr << " warning --------\n\n";
  }
  ostr << "MIRA warncode: " << wm.first << "\nTitle: " << wm.second.title << "\n\n";
  ostr << FmtText::wordWrap(wm.second.message,80) << '\n';
}


void Warnings::dumpWarnings() const
{
  FUNCSTART("void Warnings::dumpWarnings() const");

  if(WA_filebasename.empty()){
    cout << *this;
  }else{
    for(uint32 wl=0; wl<3; ++wl){
      std::string fname(WA_filebasename);
      if(wl==0){
	fname+="_critical";
      }else if(wl==1){
	fname+="_medium";
      }else if(wl==2){
	fname+="_minor";
      }
      fname+=".txt";
      std::ofstream fout(fname, std::ios::out | std::ios::trunc);
      bool hadoutput=false;
      for(auto & wme : WA_messages){
	if(wme.second.warnlevel==wl){
	  if(hadoutput) fout << "\n";
	  hadoutput=true;
	  priv_dumpSingleWarning(wme,false,fout);
	  fout << "\n--------------------------------------------------------------------------------\n";
	}
      }
      fout.close();
      if(fout.fail()){
	MIRANOTIFY(Notify::FATAL,"Could not write to " << fname << " ?");
      }
    }
  }
}

void Warnings::dumpWarning(std::string & shortcode, std::ostream & ostr) const
{
  auto mI=WA_messages.find(shortcode);
  if(mI!=WA_messages.end()){
    priv_dumpSingleWarning(*mI,true,ostr);
  }
}
