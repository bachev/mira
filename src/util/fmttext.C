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


#include "util/fmttext.H"

namespace FmtText
{
  std::string wordWrap(const char * intxt, uint32_t llen)
  {
    std::string result;
    if(intxt!=nullptr){
      uint32_t colcount=0;
      bool addinsert=false;
      std::string word;
      while(true){
	word.clear();
	// blanks at beginning of line are taken verbatim
	// tabs are converted to 8 blanks
	if(colcount==0){
	  for(; *intxt!=0; ++intxt){
	    if(*intxt==' '){
	      word+=*intxt;
	    }else if(*intxt=='\t'){
	      word+="        ";
	    }else{
	      break;
	    }
	  }
	}else{
	  // munch blanks
	  for(; *intxt!=0 && (*intxt==' ' || *intxt=='\t'); ++intxt);
	}
	// build word
	for(; *intxt!=0 && *intxt!=' ' && *intxt!='\t' && *intxt!='\n'; ++intxt) word+=*intxt;
	if(colcount==0 && word=="-"){
	  addinsert=true;
	  llen-=2;
	}
	if(colcount>0 && !word.empty() && colcount+word.size()+1 > llen){
	  result+='\n';
	  colcount=0;
	  if(addinsert) result+="  ";
	}
	if(colcount>0) {
	  result+=' ';
	  ++colcount;
	}
	colcount+=static_cast<uint32_t>(word.size());
	result+=word;
	if(*intxt=='\n'){
	  result+='\n';
	  colcount=0;
	  if(addinsert){
	    addinsert=false;
	    llen+=2;
	  }
	}
	if(*intxt==0) break;
	++intxt;
      }
    }
    return result;
  }

  std::string makeTextSign(const char * intxt, uint32_t llen)
  {
    std::string result;
    auto wrapped=wordWrap(intxt,llen-4);
    if(!wrapped.empty() && wrapped.back()!='\n') wrapped+='\n';
    for(uint32_t xi=0; xi<llen; ++xi) result+='*';
    result+='\n';
    uint32_t colcount=0;
    for(auto & se : wrapped){
      if(colcount==0){
	result+="* ";
	colcount=2;
      }
      if(se=='\n'){
	for(; colcount< llen-1; ++colcount) result+=' ';
	result+="*\n";
	colcount=0;
      }else{
	result+=se;
	++colcount;
      }
    }
    for(uint32_t xi=0; xi<llen; ++xi) result+='*';
    return result;
  }

}
