/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2002 and later by Bastien Chevreux
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


#include "generalio.H"



namespace GeneralIO {

/*************************************************************************
 *
 * reads next key value pair from an istream
 *
 * return true if successful, false if end of istream
 * lines that have # as first non-whitespace are comments and jumped over
 * key may NOT have blanks, value may
 * key and value are trimmed off leading and trailing whitespaces.
 *
 *************************************************************************/

  bool readKeyValue(std::istream & sin, std::string & key, std::string & value)
  {
    std::string tmpline;
    value.resize(0);
    key.resize(0);
    while(key.size()==0){
      // get out of here if were at the end of a stream
      if(sin.eof()) break;

      // read a new line
      while(!sin.rdstate()) {
	char tmpc= static_cast<char>(sin.get());
        if(sin.rdstate()) break;
	if(tmpc=='\n' || tmpc== -1) break;
	tmpline+=tmpc;
      }

      if(sin.rdstate()) return 0;
      //cout << "Line: " << tmpline << endl;;

      // use the istrstream facility to parse that :)
      std::istringstream in(tmpline);
      in >> key;

      // read the value string only if line is not an empty or commented line
      if(key.size()>0 && key[0]!='#') {
	// throw away leading blanks
	while(!in.eof() && isspace(in.peek())) in.get();
	// ket the value string
	while(!in.eof() && in.peek()!=-1) value+=static_cast<char>(in.get());
	// throw away trailing blanks
	while(value.size() > 0 && isspace(value[value.size()-1])) value.resize(value.size()-1);
	//cout << "--- key \"" << key << "\"\tvalue \"" << value << "\"" << endl;
      } else {
	//cout << "Comment or empty: " << key << "\t" << tmpline << endl;

	// throw away the key so that we can loop
	key.resize(0);

	// clear tmpline for next loop
	tmpline.resize(0);
      }
    }

    //cout << "KV: " << key << "\t" << value << endl;

    return (key.size()>0);
  }
}
