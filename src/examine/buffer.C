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
 * SCF examination methods: Buffer for SCF-Files
 *
 * Written by Thomas Pfisterer
 *
 * 22.09.99 Using stl for keeping the date; change of storage concept.
 *          resize is now possible.
 * 16.12.99 bufferReadCopy implemented.
 */


#include "examine/buffer.H"



bool  ScfBuffer::isValid   = false;
int32 ScfBuffer::bufReads  = 0;
int32 ScfBuffer::fileReads = 0;
int32 ScfBuffer::scfBufferSize = 25;
list<SCF_look*> ScfBuffer::scfBufferRead;
list<int32> ScfBuffer::readCount;


void ScfBuffer::initialize()
{
  scfBufferRead.clear();
  readCount.clear();

  bufReads  = 0;
  fileReads = 0;
  isValid = true;
}


void ScfBuffer::resize(int32 buffersize)
{
  scfBufferSize = buffersize;
  //  scfBufferRead.resize(scfBufferSize);
}



void ScfBuffer::discard()
{
  if (!isValid) return;    // do not discard an uninitialized buffer

  list<SCF_look*>::iterator I = scfBufferRead.begin();
  while (I != scfBufferRead.end()) {
    delete *I;
    I++;
  }

  scfBufferRead.clear();
  readCount.clear();
}



void ScfBuffer::statistics() {
  cerr << "\nSCF-Buffer Statistics" << endl;
  cerr << "---------------------" << endl;
  cerr << "Buffersize : " << scfBufferSize << endl;
  cerr << "ListBuffersize : " << scfBufferRead.size() << endl;
  cerr << "Read Operations: " << (bufReads + fileReads)<<endl;
  if (bufReads > 0) {
    cerr << "Buffered : " << 100 * bufReads /
      (bufReads + fileReads)
	 << "%" << endl;
  }
  cerr << endl;
}



SCF_look* ScfBuffer::search_insert(const char *fname)
{
  list<SCF_look*>::iterator I = scfBufferRead.begin();
  list<int32>::iterator C = readCount.begin();

  //cout << "Searching (insert) " << fname << endl;

  while (I != scfBufferRead.end()) {
    if (0 == strcmp(fname, (*I)->getFileName())) {
      (*C)++;

      bufReads++;
      return *I;
    }
    I++;
    C++;
  }
  return nullptr;
}



SCF_look* ScfBuffer::search_delete(const char *fname)
{
  list<SCF_look*>::iterator I = scfBufferRead.begin();
  list<int32>::iterator C = readCount.begin();

  //cout << "Searching (delete) " << fname << endl;
  //
  //if(fname == nullptr) {
  //  cout << *(fname-1) ;
  //}

  while (I != scfBufferRead.end()) {
    if (0 == strcmp(fname, (*I)->getFileName())) {
      (*C)--;  // decrease counter

      if (scfBufferRead.size() > scfBufferSize) {
	deleteIfUnused(scfBufferSize - scfBufferRead.size());
      }
      return *I;
    }
    I++;
    C++;
  }

  return nullptr;  // this should not happen
}



bool ScfBuffer::deleteIfUnused(const int32 anzahl)
{
  list<SCF_look*>::iterator I = scfBufferRead.begin();
  list<int32>::iterator C = readCount.begin();
  int32 count = anzahl;

  while(count>0 && I != scfBufferRead.end()) {
    if (*C < 1) {
      //cout << "Remove from buffer " << (*I)->getFileName() << endl;

      delete *I;
      I = scfBufferRead.erase(I);
      C = readCount.erase(C);
      count--;
    } else {
      I++;
      C++;
    }
  }
  //cout << "Objects in use: " << scfBufferRead.size() << endl;

  return (count == 0);
}



SCF_look* ScfBuffer::bufferRead(const Read &aRead, int32 richtung)
{
  FUNCSTART("SCF_bufferRead(const char * fname)");

  if (!isValid) {
    ScfBuffer::initialize();
  }

  // BaCh: added hasSCFData(), we know what we're doing, so the const cast is ok
  if(!const_cast<Read &>(aRead).hasSCFData()) {
    MIRANOTIFY(Notify::WARNING, "no SCF Data here " << aRead.getName() << "... is this a SCF file?.");
  }

  // cout << "BufferRead " << aRead.getName() << endl;

  SCF_look *aScfRead = search_insert(aRead.getName().c_str());

  if (aScfRead == nullptr) {
    /* Nichts gefunden, aber Read meint, da gaebe es was. Dann müssen wir wohl laden */
    try {
      aScfRead = new SCF_look;
      aScfRead->load(aRead, richtung);
      scfBufferRead.push_back(aScfRead);
      readCount.push_back(1);

      fileReads++; // count filereads

      //      cout << "using " << scfBufferRead.size() << " of "
      //	   << scfBufferSize << " objects " << endl;

      if (scfBufferRead.size() > scfBufferSize) {

	if (!deleteIfUnused(1)) {
	  cout << "Warning: Buffer overflow (no real problem)" << endl;
	}
      }
    }

    catch(Notify n) {
      cout << "ScfBuffer: Error loading SCF-File " << aRead.getName() << endl;
      delete aScfRead;
      aScfRead = nullptr;

      n.handleError(THISFUNC);

      throw Notify(Notify::WARNING, THISFUNC, "");
    }
  }

  return aScfRead;
}



SCF_look& ScfBuffer::bufferReadCopy(const Read &aRead, int32 richtung)
{
  // We read the pointer from the buffer and copy the object to result
  SCF_look *result = new SCF_look;
  SCF_look *bufferRead = nullptr;

  //Read::setCoutType(Read::AS_TEXTSHORT);
  //cout << aRead;
  bufferRead = ScfBuffer::bufferRead(aRead, richtung);

  if (bufferRead != nullptr) {
     *result = *bufferRead;
  }

  // We remove the object from the buffer: we have a copy, we do not need
  // it any more.
  ScfBuffer::bufferDelete(result);

  return *result;
}



void ScfBuffer::bufferDelete(SCF_look *aRead) {
  list<SCF_look*>::iterator I = scfBufferRead.begin();

  if (!isValid) {
    cerr << "Deleting read from uninititialize buffer.... ignored" << endl;
    return;      // do not delete from an uninitialized buffer
  }

  if (aRead == nullptr) return;

  //cout << "BufferDelete " << aRead->getFileName() << endl;

  if (search_delete(aRead->getFileName()) == nullptr) {
    //cerr << "Removing SCF-File from buffer that is not in the buffer!"<<endl;
  }
  aRead = nullptr;
}



void ScfBuffer::show()
{
  list<SCF_look*>::iterator I = scfBufferRead.begin();
  list<int32>::iterator C = readCount.begin();

  SCF_look *aRead = nullptr;

  cout << "\nBuffer contains: " << endl;
  while (I != scfBufferRead.end()) {
    cout << "\t" << (*I)->getFileName()
	 << "\t" << *C  //<< "\t" << hex << &(*C) << dec
	 << endl;

    I++;
    C++;
  }
}
