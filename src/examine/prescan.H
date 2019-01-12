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


#ifndef _prescan_prescan_H_
#define _prescan_prescan_H_


#include <iostream>
#include "mira/read.H"
#include "examine/scf_look.H"



struct prescan_parameter{
  int32    offset;
  bool     relativ;
  bool     liste;
  int32    lock_back_bases;
};



class prescan {
  enum { arraysize = 65};

  static int32 triplett_ranking[arraysize];
  static bool  ready;
  static bool  verbose;

public:
  prescan();
  ~prescan();

  int32 initialize();
  void  out(const char *aText);
  int32 load(char *filename);
  int32 getRanking(int32 index);
};



class ringbuffer {
  enum { buffersize = 10 };

  int32 buffer_index;
  int32 buffer_usage;

  int32 *buffer;
  int32 buffer_size;

public:
  ringbuffer();
  ringbuffer(int32 size);
  ~ringbuffer();
  void   add(int32 element);
  bool   getMean(int32 &mean);
};



int32 baseToIndex1(char base);
int32 baseToIndex2(char prev_base, char base);
int32 baseToIndex3(char prev_base, char base, char next_base);









#endif
