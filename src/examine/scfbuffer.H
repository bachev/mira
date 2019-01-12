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
 * SCF examination methods
 *
 */


#ifndef _scfbuffer_scfbuffer_H_
#define _scfbuffer_scfbuffer_H_


class SCF_look;


class SCF_buffer {
private:
  //  static int32    scfBufferSize;
  enum { scfBufferSize = 40 };

  static SCF_look *scfBufferRead[scfBufferSize];
  static int32    scfBufferIndex[scfBufferSize];
  //static SCF_look **scfBufferRead;
  //static int32    *scfBufferIndex;
  static int32    scfBufferLast;
  static int32    bufReads, fileReads;
  static bool     isValid;
public:
  SCF_buffer();
  ~SCF_buffer();

  static void initialize();
  static void discard();
  static SCF_look *search(const char * fName);
  static SCF_look *bufferRead(const char * fname);
  static SCF_look *bufferRead(const Read & aRead, int richtung);
  static void  bufferDelete(SCF_look *aRead);
  static void  statistics();

  static void show();
};





#endif
