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

#ifndef _bas_boostiostrutil_h_
#define _bas_boostiostrutil_h_

#include <string>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <boost/algorithm/string.hpp>


template <class FS, class IS>
void setupStreamFilterByName(FS & fs, IS & is, const std::string & filename)
{
  boost::filesystem::path runfn(filename);

  std::list<std::string> extlist;
  while(true){
    if(runfn.has_extension()){
      // native does not compile on Cygwin, gcc 4.8.3 and BOOST 1.56???
      // extlist.push_back(runfn.extension().native();
      extlist.push_back(runfn.extension().string());
    }
    if(boost::filesystem::is_symlink(runfn)){
      runfn=boost::filesystem::read_symlink(runfn);
    }else{
      break;
    }
  }


  for(auto & ext : extlist){
    if(boost::iequals(ext,".gz")){
      fs.push(boost::iostreams::gzip_decompressor());
      break;
    }else if(boost::iequals(ext,".bz2")){
      fs.push(boost::iostreams::bzip2_decompressor());
      break;
    }else if(boost::iequals(ext,".z")){
      fs.push(boost::iostreams::zlib_decompressor());
      break;
    }
  }

  fs.push(is);

}



#endif
