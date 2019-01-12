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

#ifndef _mira_warnings_h_
#define _mira_warnings_h_

#include <ostream>

#include "stdinc/defines.H"

#include <unordered_map>


class Warnings
{
  //Variables
private:

  struct warnmsg_t {
    std::string title;
    std::string message;
    uint32 warnlevel;

    warnmsg_t(std::string & t,std::string & m,uint32 wl) : title(t), message(m), warnlevel(wl) {};
    warnmsg_t() {}; // need this for emplace, should however never be used "normally"
  };

  std::unordered_map<std::string,warnmsg_t> WA_messages;

  std::string WA_filebasename;

public:


  //Functions
private:
  void init();
  void zeroVars();

  void priv_setWarning(std::string & shortcode,
		       uint32 level,
		       std::string & title,
		       std::string & message,
		       bool donotdump,
		       bool add);
  bool priv_clearWarning(std::string & shortcode);

  void priv_dumpSingleWarning(const std::pair<const std::string, warnmsg_t> & wm,
			      bool withheader, std::ostream & ostr) const;

public:
  Warnings();
  Warnings(Warnings const &other);
  ~Warnings();

  Warnings const & operator=(Warnings const & other);
  friend std::ostream & operator<<(std::ostream &ostr, Warnings const &war);

  void discard();

  inline void setWarning(std::string & shortcode,
			 uint32 level,
			 std::string & title,
			 std::string & message) {priv_setWarning(shortcode,level,title,message,false,false);}
  inline void setWarning(const char * shortcode,
			 uint32 level,
			 const char * title,
			 std::string & message) {
    std::string sc(shortcode);
    std::string ti(title);
    priv_setWarning(sc,level,ti,message,false,false);
  }
  inline void setWarning(const char * shortcode,
			 uint32 level,
			 const char * title,
			 const char * message) {
    std::string sc(shortcode);
    std::string ti(title);
    std::string mg(message);
    priv_setWarning(sc,level,ti,mg,false,false);
  }
  inline void addToWarning(std::string & shortcode,
			   uint32 level,
			   std::string & title,
			   std::string & message) {priv_setWarning(shortcode,level,title,message,false,true);}
  inline bool clearWarning(std::string & shortcode) {return priv_clearWarning(shortcode);}
  inline bool clearWarning(const char * shortcode) {
    std::string sc(shortcode);
    return priv_clearWarning(sc);
  }

  void dumpWarnings() const;
  void dumpWarning(std::string & shortcode, std::ostream & ostr) const;
  void dumpWarning(const char * shortcode, std::ostream & ostr) const {
    std::string sc(shortcode); dumpWarning(sc,ostr);
  }
  void setOutputPath(const std::string & p) {WA_filebasename=p;}
  void setOutputPath(const char * p) {WA_filebasename=p;}
};


#endif
