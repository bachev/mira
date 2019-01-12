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

#include "...H"

// Plain vanilla constructor
CLASS::CLASS()
{
  FUNCSTART("CLASS::CLASS()");

  zeroVars();
  init();

  ERROR("Not implemented yet.");
  FUNCEND();
}

void CLASS::zeroVars()
{
  FUNCSTART("void CLASS::zeroVars()");
  FUNCEND();
}

void CLASS::init()
{
  FUNCSTART("void CLASS::init()");
  FUNCEND();
}



CLASS::~CLASS()
{
  FUNCSTART("CLASS::~CLASS()");

  discard();

  ERROR("Not implemented yet.");
  FUNCEND();
}


void CLASS::discard()
{
  FUNCSTART("CLASS::discard()");

  if(??_valid==0){
    ERROR("Not implemented yet.");
  }

  zeroVars();

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//CLASS::CLASS(CLASS const &other)
//{
//  FUNCSTART("CLASS::CLASS(CLASS const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//CLASS const & CLASS::operator=(CLASS const & other)
//{
//  FUNCSTART("CLASS const & CLASS::operator=(CLASS const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, CLASS const &???)
//{
//  FUNCSTART("friend ostream & CLASS::operator<<(ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}
