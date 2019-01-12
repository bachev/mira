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

#ifndef _util_fmttext_h
#define _util_fmttext_h

#include <sstream>
#include <string>

namespace FmtText
{
  std::string wordWrap(const char * intxt, uint32_t llen);
  inline std::string wordWrap(const std::string & intxt, uint32_t llen) { return wordWrap(intxt.c_str(),llen);}
  inline std::string wordWrap(const std::ostringstream & intxt, uint32_t llen) { return wordWrap(intxt.str().c_str(),llen);}
  inline std::string wordWrap(const char * intxt) { return wordWrap(intxt,78);}
  inline std::string wordWrap(const std::string & intxt) { return wordWrap(intxt,78);}
  inline std::string wordWrap(const std::ostringstream & intxt) { return wordWrap(intxt,78);}

  std::string makeTextSign(const char * intxt, uint32_t llen);
  inline std::string makeTextSign(const std::string & intxt, uint32_t llen) { return makeTextSign(intxt.c_str(),llen);}
  inline std::string makeTextSign(const std::ostringstream & intxt, uint32_t llen) { return makeTextSign(intxt.str().c_str(),llen);}
  inline std::string makeTextSign(const char * intxt) { return makeTextSign(intxt,80);}
  inline std::string makeTextSign(const std::string & intxt) { return makeTextSign(intxt.c_str(),80);}
  inline std::string makeTextSign(const std::ostringstream & intxt) { return makeTextSign(intxt.str().c_str(),80);}

}

#endif
