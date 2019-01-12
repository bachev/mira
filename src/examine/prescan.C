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
 * Hilfsfunktionen zur Analyse der Chemieeigenschaften
 *
 */


#include "examine/prescan.H"
#include <math.h>
#include "stdlib.h"



bool  prescan::ready   = false;
bool  prescan::verbose = true;
int32 prescan::triplett_ranking[65];




prescan::prescan()
{
  initialize();
}



prescan::~prescan()
{
}




int32 prescan::initialize()
{
  if (!ready) {
    return load("prescan.par");
  }
  return 1;
}



void prescan::out(const char *aText)
{
  if (true == verbose) {
    cout << aText;
  }
}


int32 prescan::getRanking(int32 index)
{
  if (index >= 0 && index < arraysize) {
    return triplett_ranking[index];
  } else {
    return (arraysize / 2);
  }
}


int32 prescan::load(char *filename) {
  ifstream file;
  int32    index, freq, val;
  char     code[10];
  int32    loop = 0;

  ready=true;

  try {
    file.open(filename, ios::in);

    if (!file) {
      // FIXME   geht so nicht; ~ wird offensichtlich nicht expandiert
      out("Try to open ~/prescan.par.\n");
      file.open("~/prescan.par", ios::in);
    }

    // Initialize with identical values (e.g. no effect of sequence)
    for (int l=0; l < arraysize; l++) {
      triplett_ranking[l] = (arraysize / 2);
    }

    if (!file) {
      out("Couldn't open prescan information!\n");
      return -1;
    }

    loop = 0;
    do {
      file >> index >> code >> freq >> val;
      triplett_ranking[index-1] = loop;

      loop++;
    } while (!file.eof() && loop < arraysize);

    file.close();
  }
  catch(Notify n) {
    file.close();
    return -1;
  }

  out("Prescan information loaded.\n");
  return 1;
}





ringbuffer::ringbuffer()
{
  buffer_usage = 0;
  buffer_index = 0;
  buffer_size  = buffersize;
  buffer = new int32 [buffersize];
}


ringbuffer::ringbuffer(int32 size)
{
  buffer_usage = 0;
  buffer_index = 0;
  buffer_size  = size;
  buffer = new int32 [buffer_size];
}


ringbuffer::~ringbuffer()
{
  delete buffer;
}



void ringbuffer::add(int32 element)
{
  buffer[buffer_index] = element;
  buffer_index = (buffer_index+1) % buffer_size;
  buffer_usage++;
}



bool ringbuffer::getMean(int32 &mean)
{
  int32 sum = 0;

  if (buffer_usage < buffer_size) {
    mean = 0;
    return false;
  }

  for (int32 i=0; i < buffer_size; i++) {
    sum += buffer[i];
  }

  mean = sum / buffer_size;
  return true;
}



int32 baseToIndex1(char base)
{
  switch(toupper(base)) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  default: return 4;
  }
}


int32 baseToIndex2(char prev_base, char base)
{
  int i = 0;

  switch(toupper(prev_base)) {
  case 'A': i = 0;  break;
  case 'C': i = 4;  break;
  case 'G': i = 8;  break;
  case 'T': i = 12; break;
  default: return 16;
  }

  switch(toupper(base)) {
  case 'A': return i;
  case 'C': return i+1;
  case 'G': return i+2;
  case 'T': return i+3;
  default: return 16;
  }
}


int32 baseToIndex3(char prev_base, char base, char next_base)
{
  int i = 0;

  switch(toupper(prev_base)) {
  case 'A': i = 0;   break;
  case 'C': i = 16;  break;
  case 'G': i = 32;  break;
  case 'T': i = 48;  break;
  default: return 64;
  }

  switch(toupper(base)) {
  case 'A': break;
  case 'C': i += 4;  break;
  case 'G': i += 8;  break;
  case 'T': i += 12; break;
  default: return 64;
  }

  switch(toupper(next_base)) {
  case 'A': return i;
  case 'C': return i+1;
  case 'G': return i+2;
  case 'T': return i+3;
  default: return 64;
  }
}
