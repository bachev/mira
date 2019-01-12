/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
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
 *
 * Buffer for SCF_look objects
 *
 * Written by Thomas Pfisterer
 *
 * 16.12.99 bufferReadCopy: returns a copy of the object in the buffer.
 * 22.09.99 Using stl to store the data; change of storage concept.
 *
 *
 * bufferRead:    ==> pointer to SCF_look object
 *                    to each call of bufferRead you have to call:
 * bufferDelete   ==> release the object, buffer can remove it.
 *
 * bufferReadCopy ==> copy of the object. You have to delete the object
 *                    for memory deallocation.
 */


#ifndef _buffer_buffer_H_
#define _buffer_buffer_H_


#include "examine/scf_look.H"

// 	$Id$

class SCF_look;


class ScfBuffer {
private:
  static list<SCF_look*> scfBufferRead;
  static list<int32>    readCount;
  static int32    scfBufferSize;
  static int32    bufReads, fileReads;
  static bool     isValid;

  static SCF_look *search_insert(const char * fName);
  static SCF_look *search_delete(const char * fName);
  static bool deleteIfUnused(const int32 anzahl);
public:
  ScfBuffer();
  ~ScfBuffer();

  static void initialize();
  static void resize(int32 buffersize);
  static void discard();
  static SCF_look *bufferRead(const Read & aRead, int32 richtung);
  static SCF_look &bufferReadCopy(const Read & aRead, int32 richtung);
  static void  bufferDelete(SCF_look *aRead);
  static void  statistics();

  static void show();
};





#endif
