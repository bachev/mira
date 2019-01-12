/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ) and Bastien Chevreux
 * Copyright (C) 2000 and later by Bastien Chevreux
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


#include "scf.H"


#include <cstdlib>     // mkstemp() under Cygwin, also needs -U__STRICT_ANSI__


#define CEBUG(bla)

#include <cstring>


using std::cout;
using std::cerr;
using std::endl;

std::vector<std::string> SCF::SCF_suffixalternatives;



SCF::SCF()
{
  FUNCSTART("SCF::SCF()");

  zeroVars();

  if(SCF_suffixalternatives.size()==0) {
    SCF_suffixalternatives.push_back(".gz");
    SCF_suffixalternatives.push_back(".scf");
    SCF_suffixalternatives.push_back(".scf.gz");
    SCF_suffixalternatives.push_back(".Z");
    SCF_suffixalternatives.push_back(".scf.Z");
  }

  FUNCEND();
}

SCF::~SCF()
{
  FUNCSTART("SCF::~SCF()");

  discard();

  FUNCEND();
}


//void SCF::discard()
//Discard this SCF
//Free the memory that could have been allocated for the SCF.
void SCF::discard()
{
  FUNCSTART("SCF::discard()");
  if(SCF_tmpfnameforload) {
    remove(SCF_tmpfnameforload);
    delete [] SCF_tmpfnameforload;
    SCF_tmpfnameforload=nullptr;
  }

  if(SCF_samples_A!=nullptr) delete [] SCF_samples_A;
  if(SCF_samples_C!=nullptr) delete [] SCF_samples_C;
  if(SCF_samples_G!=nullptr) delete [] SCF_samples_G;
  if(SCF_samples_T!=nullptr) delete [] SCF_samples_T;

  if(SCF_peak_index) delete [] SCF_peak_index;
  if(SCF_prob_A)     delete [] SCF_prob_A;
  if(SCF_prob_C)     delete [] SCF_prob_C;
  if(SCF_prob_G)     delete [] SCF_prob_G;
  if(SCF_prob_T)     delete [] SCF_prob_T;
  if(SCF_prob_sub)     delete [] SCF_prob_sub;
  if(SCF_prob_ins)     delete [] SCF_prob_ins;
  if(SCF_prob_del)     delete [] SCF_prob_del;
  if(SCF_bases)      delete [] SCF_bases;

  if(SCF_comments) delete [] SCF_comments;

  if(SCF_private_data) delete [] SCF_private_data;

  zeroVars();

  FUNCEND();
}

void SCF::zeroVars()
{
  FUNCSTART("void SCF::zeroVars()");

  SCF_tmpfnameforload=nullptr;

  SCF_samples_A=nullptr;
  SCF_samples_C=nullptr;
  SCF_samples_G=nullptr;
  SCF_samples_T=nullptr;

  SCF_peak_index=nullptr;
  SCF_prob_A=nullptr;
  SCF_prob_C=nullptr;
  SCF_prob_G=nullptr;
  SCF_prob_T=nullptr;
  SCF_prob_sub=nullptr;
  SCF_prob_ins=nullptr;
  SCF_prob_del=nullptr;
  SCF_bases=nullptr;

  SCF_comments=nullptr;
  SCF_private_data=nullptr;

  SCF_valid_data=0;
  SCF_header.magic_number=0;

  FUNCEND();
}


// Copy constructor
//  no discard needed as this object will be freshly created when
//  called through this constructor
SCF::SCF(SCF const &other)
{
  FUNCSTART("SCF::SCF(SCF const &other)");

  zeroVars();

  *this=other;                               // call the copy operator

  FUNCEND();
}


// Copy operator, needed by copy-constructor
SCF const & SCF::operator=(SCF const & other)
{
  FUNCSTART("SCF const & SCF::operator=(SCF const & other)");

  if(this != &other){
    discard();

    SCF_header=other.SCF_header;
    if(SCF_header.private_size!=0){
      SCF_private_data= new char[SCF_header.private_size];
      memcpy(SCF_private_data, other.SCF_private_data, SCF_header.private_size);
    }

    if(other.SCF_comments!=nullptr){
      size_t len=strlen(other.SCF_comments);
      SCF_comments= new char [len+1];
      strcpy(SCF_comments, other.SCF_comments);
    }

    if(SCF_header.samples!=0){
      SCF_samples_A=new uint16 [SCF_header.samples];
      SCF_samples_C=new uint16 [SCF_header.samples];
      SCF_samples_G=new uint16 [SCF_header.samples];
      SCF_samples_T=new uint16 [SCF_header.samples];

      memcpy(SCF_samples_A, other.SCF_samples_A, SCF_header.samples*sizeof(uint16));
      memcpy(SCF_samples_C, other.SCF_samples_C, SCF_header.samples*sizeof(uint16));
      memcpy(SCF_samples_G, other.SCF_samples_G, SCF_header.samples*sizeof(uint16));
      memcpy(SCF_samples_T, other.SCF_samples_T, SCF_header.samples*sizeof(uint16));
    }

    if(SCF_header.bases!=0){
      SCF_prob_A=new uint8 [SCF_header.bases];
      SCF_prob_C=new uint8 [SCF_header.bases];
      SCF_prob_G=new uint8 [SCF_header.bases];
      SCF_prob_T=new uint8 [SCF_header.bases];
      SCF_prob_sub=new uint8 [SCF_header.bases];
      SCF_prob_ins=new uint8 [SCF_header.bases];
      SCF_prob_del=new uint8 [SCF_header.bases];

      SCF_bases=new char [SCF_header.bases];
      SCF_peak_index=new uint32 [SCF_header.bases];

      memcpy(SCF_prob_A, other.SCF_prob_A, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_prob_C, other.SCF_prob_C, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_prob_G, other.SCF_prob_G, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_prob_T, other.SCF_prob_T, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_bases, other.SCF_bases, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_prob_sub, other.SCF_prob_sub, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_prob_ins, other.SCF_prob_ins, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_prob_del, other.SCF_prob_del, SCF_header.bases*sizeof(uint8));
      memcpy(SCF_peak_index, other.SCF_peak_index, SCF_header.bases*sizeof(uint32));

    }

    SCF_valid_data=other.SCF_valid_data;
  }

  FUNCEND();
  return *this;
}

// return
//  if ok: 0 and emsg is empty
//  if non-fatal error: 1 and emsg is not empty
int8 SCF::load(const char * givenname, std::string & emsg)
{
  FUNCSTART("SCF::load(const char * filename)");

  emsg.clear();

  // free memory if previously used
  discard();

  if(givenname==nullptr){
    MIRANOTIFY(Notify::INTERNAL, "got nullptr as SCF givenname???");
  }

  if(strlen(givenname)==0){
    MIRANOTIFY(Notify::INTERNAL, "givenname length is 0???");
  }

  std::string filename=givenname;
  std::ifstream fin(filename, std::ios::in|std::ios::binary);
  if(!fin){
    bool filefound=false;
    for(uint32 ai=0; ai<SCF_suffixalternatives.size(); ai++) {
      auto compfn=filename+SCF_suffixalternatives[ai];
      fin.clear();
      fin.open(compfn, std::ios::in|std::ios::binary);
      if(fin){
	filefound=true;
	filename=compfn;
	break;
      }
    }

    if(!filefound) {
      discard();
      emsg="SCF file not found for loading, even tried .Z, .gz, .scf, .scf.gz and .scf.Z extension: ";
      emsg+=givenname;
      return 1;
    }
  }
  cout << "Loading " << filename << endl;
  SCF_lastloadedscfname=filename;

  fin.read((char *) &SCF_header, sizeof(SCF_header));
  if(fin.fail()){
    discard();
    MIRANOTIFY(Notify::FATAL, "Unexpected EOF while reading SCF header: " << filename);
  }

  convertHeaderByteOrderToHost();

  if(SCF_header.magic_number!=SCF_MAGIC){
    //uh oh ... now we have a problem.
    // this file could still be a compressed one,
    // check for gzip, compress and pack. bzip is not recognised.

    CEBUG("Not magic!\n");

    if(static_cast<uint16>(SCF_header.magic_number>>16)==0x1f8b \
       ||static_cast<uint16>(SCF_header.magic_number>>16)==0x1f9d \
       ||static_cast<uint16>(SCF_header.magic_number>>16)==0x1f1e){

      CEBUG("Compressed magic!\n");
      //Ok, seems to be a packed file. As gzip can handle all these
      // methods, use only gzip to unpack.

      // first, close the original file
      fin.close();

      SCF_tmpfnameforload=new char[32];
      strcpy(SCF_tmpfnameforload,"/tmp/mirascftmpXXXXXX");

      if(!mkstemp(SCF_tmpfnameforload)){
	discard();
	MIRANOTIFY(Notify::FATAL, "Couldn't create tmpfname. Unable to uncompress: " << filename);
      }
      char * command_buffer= new char[4096];
      sprintf(command_buffer,"gzip -d -c <%s 1>%s 2>/dev/null", filename.c_str(), SCF_tmpfnameforload);
      CEBUG(command_buffer);
      if(system(command_buffer)){
	if(system(command_buffer)){
	  if(system(command_buffer)){
	    cerr << "Load retries failed, giving up.\n";
	    delete [] command_buffer;
	    discard();
	    MIRANOTIFY(Notify::FATAL, "Unpacking of file failed. Not enough space in tmp-filesystem, corrupt file or process table full? : " << filename);
	  }
	}
      }

      delete [] command_buffer;

      //File has been uncompressed and seems ok. Open this
      // temporary file and read the header.
      fin.clear();
      fin.open(SCF_tmpfnameforload, std::ios::in|std::ios::binary);
      if(!fin){
	discard();
	MIRANOTIFY(Notify::FATAL, "Could not open temporary unpacked file for read. This should not have happened! : " << filename);
      }
      fin.read((char *) &SCF_header, sizeof(SCF_header));
      if(fin.fail()){
	discard();
	MIRANOTIFY(Notify::FATAL, "Unexpected EOF while reading header from temporarily uncompressed file. This should not have happened: " << filename);
      }

      convertHeaderByteOrderToHost();

      if(SCF_header.magic_number!=SCF_MAGIC){

	//cout << hex << SCF_header.magic_number << endl;
	//cout << SCF_tmpfnameforload << endl;
	//exit(0);

	discard();
	MIRANOTIFY(Notify::FATAL, "File is not an SCF: " << filename);
      }
    }else{
	discard();
	MIRANOTIFY(Notify::FATAL, "File is not an SCF and does not appear to be compressed: " << filename);
    }
  }

  //Before reading from file, allocate the memory we need
  if(SCF_header.samples) {
    SCF_samples_A=new uint16[SCF_header.samples];
    SCF_samples_C=new uint16[SCF_header.samples];
    SCF_samples_G=new uint16[SCF_header.samples];
    SCF_samples_T=new uint16[SCF_header.samples];
  }

  if(SCF_header.bases) {
    SCF_peak_index=new uint32[SCF_header.bases];
    SCF_prob_A= new uint8[SCF_header.bases];
    SCF_prob_C= new uint8[SCF_header.bases];
    SCF_prob_G= new uint8[SCF_header.bases];
    SCF_prob_T= new uint8[SCF_header.bases];
    SCF_bases=new char[SCF_header.bases];
    SCF_prob_sub= new uint8[SCF_header.bases];
    SCF_prob_ins= new uint8[SCF_header.bases];
    SCF_prob_del= new uint8[SCF_header.bases];
  }

  if(SCF_header.comments_size){
    SCF_comments=new SCF_Comments[SCF_header.comments_size];
  }

  if(SCF_header.private_size){
    SCF_private_data=new char[SCF_header.private_size];
  }


  char version='?';

  {
    char *ptr;
    uint32 *cptr;
    ptr=(char *)&SCF_header.version;
    cptr=(uint32 *) &SCF_header.version;
    *cptr=htonl(*cptr);
    version= ptr[0];
    *cptr=ntohl(*cptr);
  }

  if(version < '1' && version >'3'){
    discard();
    MIRANOTIFY(Notify::FATAL, "unknown SCF Version, file not supported: " << filename);
  }

  if(version=='1'){
    discard();
    MIRANOTIFY(Notify::FATAL, "SCF Version 1 files not supported: " << filename);
  }
  if(version>'3'){
    discard();
    MIRANOTIFY(Notify::FATAL, "SCF Version >3 unknown and therefore not supported: " << filename);
  }

  if(version=='2'){
    //read an SCF version 2 file
    if(SCF_header.sample_size==1){
      discard();
      MIRANOTIFY(Notify::FATAL, "Sorry, 8 bit samples in SCF V2 not supported yet, as it was not tested: " << filename);
    }


    //the temporary SCF version 2 structure gets filled with arrays.
    if(SCF_header.samples) {
      SCF_tmp_scf2.samples.samples2=new SCF_Samples2[SCF_header.samples];
      fin.read((char *) SCF_tmp_scf2.samples.samples2, SCF_header.samples*sizeof(SCF_Samples2));
    }
    if(fin.fail()){
      delete [] SCF_tmp_scf2.samples.samples2;
      discard();
      MIRANOTIFY(Notify::FATAL, "Unexpected EOF while reading samples in SCF: " << filename);
    }
    if(SCF_header.bases) {
      SCF_tmp_scf2.bases=new SCF_Bases[SCF_header.bases];
      fin.read((char *) SCF_tmp_scf2.bases, SCF_header.bases*sizeof(SCF_Bases));
    }
    if(fin.fail()){
      delete [] SCF_tmp_scf2.bases;
      delete [] SCF_tmp_scf2.samples.samples2;
      discard();
      MIRANOTIFY(Notify::FATAL, "Unexpected EOF while reading bases in SCF: " << filename);
    }

    //We now have the traces and bases in memory, scf version 1 or 2
    // Now translate these data into the class format

    for(uint32 i=0;i<SCF_header.samples;i++){
      SCF_samples_A[i]=SCF_tmp_scf2.samples.samples2[i].sample_A;
      SCF_samples_C[i]=SCF_tmp_scf2.samples.samples2[i].sample_C;
      SCF_samples_G[i]=SCF_tmp_scf2.samples.samples2[i].sample_G;
      SCF_samples_T[i]=SCF_tmp_scf2.samples.samples2[i].sample_T;
    }
    for(uint32 i=0;i<SCF_header.bases;i++){
      SCF_peak_index[i]=SCF_tmp_scf2.bases[i].peak_index;
      SCF_prob_A[i]=SCF_tmp_scf2.bases[i].prob_A;
      SCF_prob_C[i]=SCF_tmp_scf2.bases[i].prob_C;
      SCF_prob_G[i]=SCF_tmp_scf2.bases[i].prob_G;
      SCF_prob_T[i]=SCF_tmp_scf2.bases[i].prob_T;
      SCF_bases[i]=SCF_tmp_scf2.bases[i].base;
      SCF_prob_sub[i]=SCF_tmp_scf2.bases[i].prob_sub;
      SCF_prob_ins[i]=SCF_tmp_scf2.bases[i].prob_ins;
      SCF_prob_del[i]=SCF_tmp_scf2.bases[i].prob_del;
    }

    if(SCF_header.bases) delete [] SCF_tmp_scf2.bases;
    if(SCF_header.samples) delete [] SCF_tmp_scf2.samples.samples2;

    // convert the samples if needed
    if(SCF_header.sample_size==2){
      convert2ByteOrderToHost(SCF_samples_A, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_C, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_G, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_T, SCF_header.samples);
    }

  }else{
    //read an SCF version 3 file

    if(SCF_header.samples) {
      fin.read((char *) SCF_samples_A, SCF_header.samples*SCF_header.sample_size);
      fin.read((char *) SCF_samples_C, SCF_header.samples*SCF_header.sample_size);
      fin.read((char *) SCF_samples_G, SCF_header.samples*SCF_header.sample_size);
      fin.read((char *) SCF_samples_T, SCF_header.samples*SCF_header.sample_size);
    }

    if(SCF_header.bases) {
      fin.read((char *) SCF_peak_index, SCF_header.bases*4);
      fin.read((char *) SCF_prob_A, SCF_header.bases);
      fin.read((char *) SCF_prob_C, SCF_header.bases);
      fin.read((char *) SCF_prob_G, SCF_header.bases);
      fin.read((char *) SCF_prob_T, SCF_header.bases);
      fin.read((char *) SCF_bases, SCF_header.bases);
      fin.read((char *) SCF_prob_sub, SCF_header.bases);
      fin.read((char *) SCF_prob_ins, SCF_header.bases);
      fin.read((char *) SCF_prob_del, SCF_header.bases);
    }

    // convert the samples if needed
    if(SCF_header.sample_size==2){
      convert2ByteOrderToHost(SCF_samples_A, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_C, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_G, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_T, SCF_header.samples);
    }

    if(fin.fail()){
      discard();
      MIRANOTIFY(Notify::FATAL, "Unexpected EOF, data missing in SCF: " << filename);
    }

    if(SCF_header.sample_size==1){
      // *sigh* SCF V3 with 1 byte samples: SCF_samples_? are filled
      //  with 8 bit instead of 16 bit values.

      // samples in SCF V3 are stored as delta values, undelta them
      // I don't like reinterpret_cast in this case. Paranoia?
      undelta8((uint8 *)SCF_samples_A);
      undelta8((uint8 *)SCF_samples_C);
      undelta8((uint8 *)SCF_samples_G);
      undelta8((uint8 *)SCF_samples_T);

      //  copy and convert 8 bit to 16 bit
      // I don't like reinterpret_cast in this case. Paranoia?
      uint8  * src=((uint8 *)(SCF_samples_A))+SCF_header.samples;
      //      src+=SCF_header.samples;
      uint16 * dst=&SCF_samples_A[SCF_header.samples];
      for(; dst!=SCF_samples_A;){
	*(--dst)=static_cast<uint16>(*(--src));
      }

      src=((uint8 *)(SCF_samples_C))+SCF_header.samples;
      dst=&SCF_samples_C[SCF_header.samples];
      for(; dst!=SCF_samples_C;){
	*(--dst)=static_cast<uint16>(*(--src));
      }

      src=((uint8 *)(SCF_samples_G))+SCF_header.samples;
      dst=&SCF_samples_G[SCF_header.samples];
      for(; dst!=SCF_samples_G;){
	*(--dst)=static_cast<uint16>(*(--src));
      }

      src=((uint8 *)(SCF_samples_T))+SCF_header.samples;
      dst=&SCF_samples_T[SCF_header.samples];
      for(; dst!=SCF_samples_T;){
	*(--dst)=static_cast<uint16>(*(--src));
      }
    }else{
      // samples in SCF V3 are stored as delta values, undelta them
      undelta16(SCF_samples_A);
      undelta16(SCF_samples_C);
      undelta16(SCF_samples_G);
      undelta16(SCF_samples_T);
    }
  }

  //good, now read the comments and private data (if any)

  if(SCF_header.comments_size){
    fin.read((char *) SCF_comments, SCF_header.comments_size);
  }
  if(fin.fail()){
    discard();
    MIRANOTIFY(Notify::FATAL,"Unexpected EOF while reading comments in SCF: " << filename);
  }

  if(SCF_header.private_size){
    fin.read((char *) SCF_private_data, SCF_header.private_size);
  }
  if(fin.fail()){
    discard();
    MIRANOTIFY(Notify::FATAL, "Unexpected EOF while reading private data SCF: " << filename);
  }

  // Convert the byte order (no effect on machines where the architecture
  //  already works in network byte order)
  // convert the peak indices
  convert4ByteOrderToHost(SCF_peak_index, SCF_header.bases);

  // check if the base peak positions and the number of samples disagree
  if(SCF_header.bases && (SCF_header.samples <= SCF_peak_index[SCF_header.bases-1])){
    // Damn, they disagree: allocate new samples area and
    //  fill the tail with 0 (if not too much)

    WARNING("SCF file \""<< filename << "\" has inconsistent peak positions and number of samples.\nError corrected internally, but you DO want to have a look at that file and find out how it was corrupted!");

    uint16 * src;
    uint16 * dst;
    uint32 new_num_samples=SCF_peak_index[SCF_header.bases-1]+1;

    uint16 * tmp_scf_samples= new uint16[new_num_samples];
    src=SCF_samples_A;
    dst=tmp_scf_samples;
    //    memset(tmp_scf_samples, 0, new_num_samples*2);
    for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
      *dst=*src;
    }
    delete [] SCF_samples_A;
    SCF_samples_A=tmp_scf_samples;

    tmp_scf_samples= new uint16[new_num_samples];
    src=SCF_samples_C;
    dst=tmp_scf_samples;
    //    memset(tmp_scf_samples, 0, new_num_samples*2);
    for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
      *dst=*src;
    }
    delete [] SCF_samples_C;
    SCF_samples_C=tmp_scf_samples;

    tmp_scf_samples= new uint16[new_num_samples];
    src=SCF_samples_G;
    dst=tmp_scf_samples;
    //    memset(tmp_scf_samples, 0, new_num_samples*2);
    for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
      *dst=*src;
    }
    delete [] SCF_samples_G;
    SCF_samples_G=tmp_scf_samples;

    tmp_scf_samples= new uint16[new_num_samples];
    src=SCF_samples_T;
    dst=tmp_scf_samples;
    //    memset(tmp_scf_samples, 0, new_num_samples*2);
    for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
      *dst=*src;
    }
    delete [] SCF_samples_T;
    SCF_samples_T=tmp_scf_samples;

    SCF_header.samples=new_num_samples;
  }

  // if a temporary file was used, remove it
  if(SCF_tmpfnameforload) {
    remove(SCF_tmpfnameforload);
    delete [] SCF_tmpfnameforload;
    SCF_tmpfnameforload=nullptr;
  }

  //looks like everything's ok
  SCF_valid_data=1;

  //// but we don't trust all this, so make a thorough check
  //if(checkSCFDataOK(true)==false){
  //  if(checkSCFDataOK(false)==false){
  //    discard();
  //    throw Notify(Notify::WARNING, THISFUNC, filename.c_str(), ": Broken SCF file will not be used and you DO want to have a look at that file and find out how it was corrupted!");
  //  }
  //}

  FUNCEND();

  return 0;
}



// trytocorrect: phred has got a problem with SCF files, sometimes Peak indices
//   are wrongly placed.
bool SCF::checkSCFDataOK()
{
  FUNCSTART("bool SCF::checkSCFData()");

  if(SCF_valid_data==0) return false;

  // darn, it might be, that the peak positions are _not_ in order due to
  //  compression
  // we can't check it.

  return true;

//  {
//    // Check peak indices
//    uint32 * ptr=SCF_peak_index;
//    for(uint32 i=0; i<SCF_header.bases-1; i++, ptr++){
//	if(*ptr > *(ptr+1)){
//	  if(trytocorrect){
//	    *ptr=*(ptr+1);
//	  }else{
//	    WARNING("Inconsistent peak positions indices.\n");
//	    return false;
//	  }
//	}
//    }
//  }

  FUNCEND();

  return true;
}

void SCF::correctPeakIndices()
{
  if(SCF_header.bases==0) return;

  // Check peak indices
  uint32 * ptr=SCF_peak_index;
  for(uint32 i=0; i<SCF_header.bases-1; i++, ptr++){
    if(*ptr > *(ptr+1)){
      *ptr=*(ptr+1);
    }
  }
}


inline void SCF::convert4ByteOrderToHost(uint32 *ptr, uint32 num)
{
  for(;num;num--,ptr++){
    *ptr=ntohl(*ptr);
  }
}

inline void SCF::convert2ByteOrderToHost(uint16 *ptr, uint32 num)
{
  for(;num;num--,ptr++){
    *ptr=ntohs(*ptr);
  }
}

inline void SCF::convertHeaderByteOrderToHost()
{
  convert4ByteOrderToHost((uint32 *) &SCF_header,
			  (sizeof(SCF_header)/sizeof(uint32)));
}


inline void SCF::convert4ByteOrderToNetwork(uint32 *ptr, uint32 num)
{
  for(;num;num--,ptr++){
    *ptr=htonl(*ptr);
  }
}

inline void SCF::convert2ByteOrderToNetwork(uint16 *ptr, uint32 num)
{
  for(;num;num--,ptr++){
    *ptr=htons(*ptr);
  }
}

inline void SCF::convertHeaderByteOrderToNetwork()
{
  convert4ByteOrderToNetwork((uint32 *) &SCF_header,
			  (sizeof(SCF_header)/sizeof(uint32)));
}



//undelta SCF-V3 Samples
//No Errorchecking, expects valid data
//Algorithm as described in Staden manual and adapted
// TODO as template
void SCF::undelta16(uint16 * ptr){
  FUNCSTART("void SCF::undelta16(uint16 * ptr)");

  uint16 p_sample=0;
  uint16 * optr=ptr;
  for(uint32 i=0;i<SCF_header.samples;i++){
    *ptr=*ptr+p_sample;
    p_sample=*ptr++;
  }
  p_sample=0;
  for(uint32 i=0;i<SCF_header.samples;i++){
    *optr=*optr+p_sample;
    p_sample=*optr++;
  }
  FUNCEND();
}

void SCF::delta16(uint16 * ptr)
{
  FUNCSTART("void SCF::delta16(uint16 * ptr)");

  uint16 p_delta=0;
  uint16 p_sample=0;
  uint16 * optr=ptr;
  for(uint32 i=0; i<SCF_header.samples; i++){
    p_sample=*ptr;
    *ptr=*ptr-p_delta;
    ptr++;
    p_delta=p_sample;
  }
  p_delta=0;
  for(uint32 i=0; i<SCF_header.samples; i++){
    p_sample=*optr;
    *optr=*optr-p_delta;
    optr++;
    p_delta=p_sample;
  }

  FUNCEND();
}

void SCF::undelta8(uint8 * ptr){
  FUNCSTART("void SCF::undelta8(uint8 * ptr)");

  uint8 p_sample=0;
  uint8 * optr=ptr;
  for(uint32 i=0;i<SCF_header.samples;i++){
    *ptr=*ptr+p_sample;
    p_sample=*ptr++;
  }
  p_sample=0;
  for(uint32 i=0;i<SCF_header.samples;i++){
    *optr=*optr+p_sample;
    p_sample=*optr++;
  }
  FUNCEND();
}

void SCF::delta8(uint8 * ptr)
{
  FUNCSTART("void SCF::delta8(uint16 * ptr)");

  uint8 p_delta=0;
  uint8 p_sample=0;
  uint8 * optr=ptr;
  for(uint32 i=0; i<SCF_header.samples; i++){
    p_sample=*ptr;
    *ptr=*ptr-p_delta;
    ptr++;
    p_delta=p_sample;
  }
  p_delta=0;
  for(uint32 i=0; i<SCF_header.samples; i++){
    p_sample=*optr;
    *optr=*optr-p_delta;
    optr++;
    p_delta=p_sample;
  }

  FUNCEND();
}


int8 SCF::dump()
{
  FUNCSTART("SCF::dump()");

  if(!SCF_valid_data){
    cout << "SCF Object not initialised.\n";
    return 1;
  }
  cout << "[Header]\n";
  char *ptr;
  uint32 *cptr;

  ptr=(char *)&SCF_header.magic_number;
  cptr=&SCF_header.magic_number;

  *cptr=htonl(*cptr);
  cout << ptr[0] << ptr[1]<<ptr[2]<<ptr[3] << "\t\t# magic_number" << endl;
  *cptr=ntohl(*cptr);

  cout << SCF_header.samples << "\t\t# samples" << endl;
  cout << SCF_header.samples_offset << "\t\t# samples_offset" << endl;
  cout << SCF_header.bases << "\t\t# bases" << endl;
  cout << "\t\t# (OBSOLETE)" << endl;
  cout << "\t\t# (OBSOLETE)" << endl;
  cout << SCF_header.bases_offset << "\t\t# bases_offset" << endl;
  cout << SCF_header.comments_size << "\t\t# comments_size" << endl;
  cout << SCF_header.comments_offset << "\t\t# comments_offset" << endl;
  ptr=(char *)&SCF_header.version;
  cptr=(uint32 *) &SCF_header.version;
  *cptr=htonl(*cptr);
  cout << ptr[0] << ptr[1]<<ptr[2]<<ptr[3] << "\t\t# magic_number" << endl;
  *cptr=ntohl(*cptr);
  cout << SCF_header.sample_size << "\t\t# sample_size" << endl;
  cout << SCF_header.code_set << "\t\t# code_set" << endl;
  cout << SCF_header.private_size << "\t\t# private_size" << endl;
  cout << SCF_header.private_offset << "\t\t# private_offset" << endl;
  for(int i=0;i<18;i++){
    cout << SCF_header.spare[i] << "\t\t# spare["<<i<<"]" << endl;
  }

  cout << "\n\n[Bases]\n";
  //flush the cout buffer (ther still could be something in it)
  cout.flush();

  // use printf for this as it is a _lot_ faster
  for(uint32 i=0;i<SCF_header.bases;i++){
    printf("%c %6i %3i %3i %3i %3i # %i\n",SCF_bases[i], SCF_peak_index[i], \
	   SCF_prob_A[i] ,SCF_prob_C[i] ,SCF_prob_G[i] ,SCF_prob_T[i], i);
  }
  // flush the stdout
  fflush(stdout);


  cout << "\n\n[A-Trace]\n";
  dumpSample(&SCF_samples_A[0]);
  cout << "\n\n[C-Trace]\n";
  dumpSample(&SCF_samples_C[0]);
  cout << "\n\n[G-Trace]\n";
  dumpSample(&SCF_samples_G[0]);
  cout << "\n\n[T-Trace]\n";
  dumpSample(&SCF_samples_T[0]);

  if(SCF_header.comments_size){
    cout << "\n\n[Comments]\n" << SCF_comments;
  }

  FUNCEND();
  return 0;
}


// print the trace-samples to stdout
// uses c-printf because it's _a lot_ faster
void SCF::dumpSample(const uint16 * ptr)
{
  FUNCSTART("void SCF::dumpSample(const uint16 * ptr)");

  // flush the cout buffer (ther still could be something in it)
  cout.flush();

  for(uint32 i=0;i<SCF_header.samples;){
    //cout <<hex<< *ptr++ << "\t# " << i++ << endl;
    printf("%d\t# %i\n",*ptr++,i++);
  }

  // flush the stdout
  fflush(stdout);
  FUNCEND();
}




void SCF::save(const char * filename)
{
  FUNCSTART("void SCF::save(const char * filename)");

  if(!SCF_valid_data){
    MIRANOTIFY(Notify::WARNING, "object not valid/initialised: " << filename);
  }

  uint32 *cptr;
  cptr=(uint32*)&SCF_header.version;
  *cptr=htonl(*cptr);
  if(SCF_header.version[0]=='2'){
    //WARNING("File read was SCF version 2, will now write as SCF version 3.");
    SCF_header.version[0]='3';
  }
  SCF_header.version[2]='1';
  SCF_header.version[3]='0';
  *cptr=ntohl(*cptr);

  recalcHeader();

  std::ofstream fout(filename, std::ios::out|std::ios::binary);
  if(!fout){
    MIRANOTIFY(Notify::WARNING, "Could not open file for saving: " << filename);
  }

  convertHeaderByteOrderToNetwork();
  fout.write((char *) &SCF_header, sizeof(SCF_header));
  convertHeaderByteOrderToHost();

  if(SCF_header.samples) {
    // is sample precision 1 byte or 2 byte?
    if(SCF_header.sample_size==2) {
      // 2 byte

      // delta the SCF samples temporarily for saving
      delta16(SCF_samples_A);
      delta16(SCF_samples_C);
      delta16(SCF_samples_G);
      delta16(SCF_samples_T);

      convert2ByteOrderToNetwork(SCF_samples_A, SCF_header.samples);
      convert2ByteOrderToNetwork(SCF_samples_C, SCF_header.samples);
      convert2ByteOrderToNetwork(SCF_samples_G, SCF_header.samples);
      convert2ByteOrderToNetwork(SCF_samples_T, SCF_header.samples);

      fout.write((char *) SCF_samples_A, SCF_header.samples*SCF_header.sample_size);
      fout.write((char *) SCF_samples_C, SCF_header.samples*SCF_header.sample_size);
      fout.write((char *) SCF_samples_G, SCF_header.samples*SCF_header.sample_size);
      fout.write((char *) SCF_samples_T, SCF_header.samples*SCF_header.sample_size);

      convert2ByteOrderToHost(SCF_samples_A, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_C, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_G, SCF_header.samples);
      convert2ByteOrderToHost(SCF_samples_T, SCF_header.samples);

      // undelta again
      undelta16(SCF_samples_A);
      undelta16(SCF_samples_C);
      undelta16(SCF_samples_G);
      undelta16(SCF_samples_T);
    } else {
      // *sigh* 1 byte

      uint8 * tmpsamplesA= new uint8[SCF_header.samples];
      uint8 * tmpsamplesC= new uint8[SCF_header.samples];
      uint8 * tmpsamplesG= new uint8[SCF_header.samples];
      uint8 * tmpsamplesT= new uint8[SCF_header.samples];

      uint8 * dst=tmpsamplesA;
      uint16 * src=SCF_samples_A;
      for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
	*dst=static_cast<uint8>(*src);
      }
      dst=tmpsamplesC;
      src=SCF_samples_C;
      for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
	*dst=static_cast<uint8>(*src);
      }
      dst=tmpsamplesG;
      src=SCF_samples_G;
      for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
	*dst=static_cast<uint8>(*src);
      }
      dst=tmpsamplesT;
      src=SCF_samples_T;
      for(uint32 i=0; i<SCF_header.samples; i++, src++, dst++){
	*dst=static_cast<uint8>(*src);
      }

      delta8(tmpsamplesA);
      delta8(tmpsamplesC);
      delta8(tmpsamplesG);
      delta8(tmpsamplesT);

      fout.write((char *) tmpsamplesA, SCF_header.samples*SCF_header.sample_size);
      fout.write((char *) tmpsamplesC, SCF_header.samples*SCF_header.sample_size);
      fout.write((char *) tmpsamplesG, SCF_header.samples*SCF_header.sample_size);
      fout.write((char *) tmpsamplesT, SCF_header.samples*SCF_header.sample_size);

      delete [] tmpsamplesA;
      delete [] tmpsamplesC;
      delete [] tmpsamplesG;
      delete [] tmpsamplesT;
    }
  }

  if(SCF_header.bases) {
    convert4ByteOrderToNetwork(SCF_peak_index, SCF_header.bases);
    fout.write((char *) SCF_peak_index, SCF_header.bases*4);
    convert4ByteOrderToHost(SCF_peak_index, SCF_header.bases);

    fout.write((char *) SCF_prob_A, SCF_header.bases);
    fout.write((char *) SCF_prob_C, SCF_header.bases);
    fout.write((char *) SCF_prob_G, SCF_header.bases);
    fout.write((char *) SCF_prob_T, SCF_header.bases);
    fout.write((char *) SCF_bases, SCF_header.bases);

    // create a temp array filed with zeros to write 3*SCF_header.bases
    //  as space which is not used yet
    {
      char * bla=new char[SCF_header.bases*3];
      memset(bla,0,SCF_header.bases*3);             // clear that block
      fout.write((char *) bla, SCF_header.bases*3);
      delete [] bla;
    }
  }

  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL, "Error while writing data into SCF: " << filename);
  }

  //good, now write the comments and private data (if any)

  if(SCF_header.comments_size){
    fout.write((char *) SCF_comments, SCF_header.comments_size);
  }
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL, "Error while writing comments into SCF: " << filename);
  }

  if(SCF_header.private_size){
    fout.write((char *) SCF_private_data, SCF_header.private_size);
  }
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL, "Error while writing private data SCF: " << filename);
  }

  FUNCEND();
}


// Cuts out the part of the SCF between frombase (including) and tobase
//  (excluding)
// After the operation, the SCF consists only of the bases between those
//  two bounds, the rest is discarded.
void SCF::cutBases(uint32 frombase, uint32 tobase)
{
  FUNCSTART("void SCF::cutBases(uint32 frombase, uint32 tobase)");

  if(!SCF_valid_data){
    MIRANOTIFY(Notify::WARNING, "object not valid/initialised.");
  }

  if(tobase>SCF_header.bases){
    MIRANOTIFY(Notify::FATAL, "Right bound > number of bases in SCF.");
  }
  if(frombase>tobase){
    MIRANOTIFY(Notify::FATAL, "Left bound > right bound?");
  }

  // search for how many samples on the left side to cut
  // if frombase = 0 then none
  uint32 lsamples_to_cut=0;

  if(frombase){
    lsamples_to_cut=(SCF_peak_index[frombase]-SCF_peak_index[frombase-1])/2+SCF_peak_index[frombase-1];
  }

  // search for how many samples on the right side to cut
  // if frombase = SCF_header.bases then none
  uint32 rsamples_to_cut=SCF_header.samples;

  if(tobase!=SCF_header.bases){
    rsamples_to_cut=(SCF_peak_index[tobase]-SCF_peak_index[tobase-1])/2+SCF_peak_index[tobase-1];
  }

  uint32 samples_to_copy=rsamples_to_cut-lsamples_to_cut;

  {
    // cut the samples

    uint16 * tmpsamples=new uint16[samples_to_copy];
    memcpy(tmpsamples, SCF_samples_A+lsamples_to_cut, samples_to_copy*2);
    delete [] SCF_samples_A;
    SCF_samples_A=tmpsamples;

    tmpsamples=new uint16[samples_to_copy];
    memcpy(tmpsamples, SCF_samples_C+lsamples_to_cut, samples_to_copy*2);
    delete [] SCF_samples_C;
    SCF_samples_C=tmpsamples;

    tmpsamples=new uint16[samples_to_copy];
    memcpy(tmpsamples, SCF_samples_G+lsamples_to_cut, samples_to_copy*2);
    delete [] SCF_samples_G;
    SCF_samples_G=tmpsamples;

    tmpsamples=new uint16[samples_to_copy];
    memcpy(tmpsamples, SCF_samples_T+lsamples_to_cut, samples_to_copy*2);
    delete [] SCF_samples_T;
    SCF_samples_T=tmpsamples;
  }

  uint32 bases_to_copy=tobase-frombase;

  {
    // cut the base-probabilities and the bases

    uint8 * tmpprob= new uint8[bases_to_copy];
    memcpy(tmpprob, SCF_prob_A+frombase, bases_to_copy);
    delete [] SCF_prob_A;
    SCF_prob_A=tmpprob;

    tmpprob= new uint8[bases_to_copy];
    memcpy(tmpprob, SCF_prob_C+frombase, bases_to_copy);
    delete [] SCF_prob_C;
    SCF_prob_C=tmpprob;

    tmpprob= new uint8[bases_to_copy];
    memcpy(tmpprob, SCF_prob_G+frombase, bases_to_copy);
    delete [] SCF_prob_G;
    SCF_prob_G=tmpprob;

    tmpprob= new uint8[bases_to_copy];
    memcpy(tmpprob, SCF_prob_T+frombase, bases_to_copy);
    delete [] SCF_prob_T;
    SCF_prob_T=tmpprob;

    char * tmpbases= new char[bases_to_copy];
    memcpy(tmpbases, SCF_bases+frombase, bases_to_copy);
    delete [] SCF_bases;
    SCF_bases=tmpbases;
  }

  {
    // cut the peak index and adjust it
    uint32 * tmpindex= new uint32[bases_to_copy];
    memcpy(tmpindex, SCF_peak_index+frombase, bases_to_copy*4);
    delete [] SCF_peak_index;
    SCF_peak_index=tmpindex;

    {
      for(uint32 i=0;i<bases_to_copy;i++){
	SCF_peak_index[i] -= lsamples_to_cut;
      }
    }
  }

  // Now adjust the header values
  SCF_header.samples=samples_to_copy;
  SCF_header.bases=bases_to_copy;

  // and recalculate the Header
  recalcHeader();


  FUNCEND();
}



void SCF::recalcHeader()
{
  uint32 total_offset=128;

  SCF_header.samples_offset= total_offset;
  total_offset+=SCF_header.samples*SCF_header.sample_size*4;
  SCF_header.bases_offset=total_offset;
  total_offset+=SCF_header.bases*12;
  SCF_header.comments_offset=total_offset;
  total_offset+=SCF_header.comments_size;
  SCF_header.private_offset=total_offset;
}


//void SCF::setBases(uint32 numbases, const char * newbases, const uint8 * newprobs, const uint32 * newpeakindex)
void SCF::setBases(const std::string & newbases, const std::vector<uint8> & newprobs, const std::vector<uint32> & newpeakindex)
{
  FUNCSTART("void SCF::setBases(uint32 numbases, const char * bases, const char * prob, const uint32 * peakpos)");

  if(newbases.size() != newprobs.size()
     || newprobs.size() != newpeakindex.size()) {
    cerr << "Number of bases: " << newbases.size() << endl;
    cerr << "Number of quals: " << newprobs.size() << endl;
    cerr << "Number of peaks: " << newpeakindex.size() << endl;
    MIRANOTIFY(Notify::FATAL, "Number of elements in bases, qualities and peakindex arrays are not equal!");
  }

  if(SCF_peak_index!=nullptr) delete [] SCF_peak_index;
  if(SCF_prob_A!=nullptr) delete [] SCF_prob_A;
  if(SCF_prob_C!=nullptr) delete [] SCF_prob_C;
  if(SCF_prob_G!=nullptr) delete [] SCF_prob_G;
  if(SCF_prob_T!=nullptr) delete [] SCF_prob_T;
  if(SCF_bases!=nullptr) delete [] SCF_bases;

  SCF_peak_index=nullptr;
  SCF_prob_A=nullptr;
  SCF_prob_C=nullptr;
  SCF_prob_G=nullptr;
  SCF_prob_T=nullptr;
  SCF_bases=nullptr;

  uint32 numbases=static_cast<uint32>(newbases.size());
  SCF_header.bases=numbases;

  if(numbases!=0) {

    SCF_peak_index= new uint32[numbases];
    SCF_prob_A= new uint8[numbases];
    SCF_prob_C= new uint8[numbases];
    SCF_prob_G= new uint8[numbases];
    SCF_prob_T= new uint8[numbases];
    SCF_prob_sub= new uint8[numbases];
    SCF_prob_ins= new uint8[numbases];
    SCF_prob_del= new uint8[numbases];
    SCF_bases= new char[numbases];

    memset(SCF_prob_A, 0, numbases*sizeof(uint8));
    memset(SCF_prob_C, 0, numbases*sizeof(uint8));
    memset(SCF_prob_G, 0, numbases*sizeof(uint8));
    memset(SCF_prob_T, 0, numbases*sizeof(uint8));
    memset(SCF_prob_sub, 0, numbases*sizeof(uint8));
    memset(SCF_prob_ins, 0, numbases*sizeof(uint8));
    memset(SCF_prob_del, 0, numbases*sizeof(uint8));

    //memcpy(SCF_peak_index, newpeakindex, numbases*sizeof(uint32));
    //memcpy(SCF_bases, newbases, numbases*sizeof(uint8));
    {
      uint32 * target=SCF_peak_index;
      auto I=newpeakindex.cbegin();
      for(uint32 i=0; I != newpeakindex.end(); I++, target++, i++){
	if(*I >= SCF_header.samples) {
	  cerr << "Pos: " << i << "\t illegal peak index: " << *I << "\tas SCF has only " << SCF_header.samples << " samples." << endl;
	  MIRANOTIFY(Notify::FATAL, "Illegal base1.");
	}

	*target=*I;
      }
    }
    {
      char * target=SCF_bases;
      for(uint32 i=0; i!= newbases.size(); i++, target++){
	*target=newbases[i];
      }
    }

    uint8 * probA=SCF_prob_A;
    uint8 * probC=SCF_prob_C;
    uint8 * probG=SCF_prob_G;
    uint8 * probT=SCF_prob_T;
    auto theprob=newprobs.cbegin();
    for(uint32 i=0; i<numbases; i++, theprob++, probA++, probC++, probG++, probT++) {
      switch(tolower(newbases[i])){
      case 'a': {
	*probA=*theprob;
	break;
      }
      case 'c': {
	*probC=*theprob;
	break;
      }
      case 'g': {
	*probG=*theprob;
	break;
      }
      case 't': {
	*probT=*theprob;
	break;
      }
      case 'n': {
	*probA=*theprob;
	*probC=*theprob;
	*probG=*theprob;
	*probT=*theprob;
	break;
      }
      case 'r': {
	*probA=*theprob;
	*probG=*theprob;
	break;
      }
      case 'y': {
	*probC=*theprob;
	*probT=*theprob;
	break;
      }
      case 'm': {
	*probA=*theprob;
	*probC=*theprob;
	break;
      }
      case 's': {
	*probC=*theprob;
	*probG=*theprob;
	break;
      }
      case 'k': {
	*probG=*theprob;
	*probT=*theprob;
	break;
      }
      case 'w': {
	*probA=*theprob;
	*probT=*theprob;
	break;
      }
      case 'h': {
	*probA=*theprob;
	*probC=*theprob;
	*probT=*theprob;
	break;
      }
      case 'b': {
	*probC=*theprob;
	*probG=*theprob;
	*probT=*theprob;
	break;
      }
      case 'v': {
	*probA=*theprob;
	*probC=*theprob;
	*probG=*theprob;
	break;
      }
      case 'd': {
	*probA=*theprob;
	*probG=*theprob;
	*probT=*theprob;
	break;
      }
      default:{
	cerr << "Pos1: " << i << "\t illegal base: '" << newbases[i] << "'" << endl;
	MIRANOTIFY(Notify::FATAL, "Illegal base1.");
      }
      }
    }
  }

  recalcHeader();

  FUNCEND();
}



//{\em Ambiguity Codes} R (G oder A), Y (T oder C), M (A oder C), S (C oder
//G), K (G oder T), W (T oder A), H (not G), B (not A), V (not T), D (not
//C), N (A, C, G oder T)  (z.B. \cite{ALLEX_97})

void SCF::transposeAmbiguityCodes()
{
  FUNCSTART("void SCF::transposeAmbiguityCodes()");

  if(SCF_valid_data){
    char * ptr = SCF_bases;
    for(uint32 i=0; i< SCF_header.bases; i++, ptr++){
      switch(*ptr){
      case 'A':
      case 'a':
      case 'C':
      case 'c':
      case 'G':
      case 'g':
      case 'T':
      case 't':
      case 'N':
      case 'n': break;

      case '-': {
	*ptr='N';
	break;
      }

      case 'X':

      case 'R':
      case 'Y':
      case 'M':
      case 'S':
      case 'K':
      case 'W':
      case 'H':
      case 'B':
      case 'V':
      case 'D': {
	*ptr='N';
	break;
      }

      case 'x':

      case 'r':
      case 'y':
      case 'm':
      case 's':
      case 'k':
      case 'w':
      case 'h':
      case 'b':
      case 'v':
      case 'd':{
	*ptr='n';
	break;
      }
      default: {
	cerr << "Pos2: " << i << "\t illegal base: '" << *ptr << "'" << endl;
	MIRANOTIFY(Notify::FATAL, "Illegal base2.");
      }
      }
    }
  }
  FUNCEND();
}

char  SCF::getBase(uint32 index) const
{
  FUNCSTART("");
  BUGIFTHROW(index>=SCF_header.bases, "index > number of bases");
  FUNCEND();
  return SCF_bases[index];
}

uint8 SCF::getCalledBaseProb(uint32 index) const
{
  FUNCSTART("");
  BUGIFTHROW(index>=SCF_header.bases, "index > number of bases");
  FUNCEND();

  switch(toupper(SCF_bases[index])){
  case 'A':{
    return SCF_prob_A[index];
  }
  case 'C':{
    return SCF_prob_C[index];
  }
  case 'G':{
    return SCF_prob_G[index];
  }
  case 'T':{
    return SCF_prob_T[index];
  }
  default:{
    return std::max(SCF_prob_T[index], std::max(SCF_prob_G[index], std::max(SCF_prob_A[index],SCF_prob_C[index])));
  }
  }
}

uint8 SCF::getAProb(uint32 index) const
{
  FUNCSTART("uint8 SCF::getAProb(uint32 index) const");
  BUGIFTHROW(index>=SCF_header.bases, "index > number of bases");
  FUNCEND();
  return SCF_prob_A[index];
}

uint8 SCF::getCProb(uint32 index) const
{
  FUNCSTART("uint8 SCF::getCProb(uint32 index) const");
  BUGIFTHROW(index>=SCF_header.bases, "index > number of bases");
  FUNCEND();
  return SCF_prob_C[index];
}

uint8 SCF::getGProb(uint32 index) const
{
  FUNCSTART("uint8 SCF::getGProb(uint32 index) const");
  BUGIFTHROW(index>=SCF_header.bases, "index > number of bases");
  FUNCEND();
  return SCF_prob_G[index];
}

uint8 SCF::getTProb(uint32 index) const
{
  FUNCSTART("uint8 SCF::getTProb(uint32 index) const");
  BUGIFTHROW(index>=SCF_header.bases, "index > number of bases");
  FUNCEND();
  return SCF_prob_T[index];
}


uint16 SCF::getASample(uint32 samplepos) const
{
  FUNCSTART("uint16 getASample(uint32 samplepos) const");
  BUGIFTHROW(samplepos>=SCF_header.samples, "samplepos > number of samples");
  FUNCEND();
  return SCF_samples_A[samplepos];
}

uint16 SCF::getCSample(uint32 samplepos) const
{
  FUNCSTART("uint16 getCSample(uint32 samplepos) const");
  BUGIFTHROW(samplepos>=SCF_header.samples, "samplepos > number of samples");
  FUNCEND();
  return SCF_samples_C[samplepos];
}

uint16 SCF::getGSample(uint32 samplepos) const
{
  FUNCSTART("uint16 getGSample(uint32 samplepos) const");
  BUGIFTHROW(samplepos>=SCF_header.samples, "samplepos > number of samples");
  FUNCEND();
  return SCF_samples_G[samplepos];
}

uint16 SCF::getTSample(uint32 samplepos) const
{
  FUNCSTART("uint16 getTSample(uint32 samplepos) const");
  BUGIFTHROW(samplepos>=SCF_header.samples, "samplepos > number of samples");
  FUNCEND();
  return SCF_samples_T[samplepos];
}


uint32 SCF::getPeakIndex(uint32 basepos) const
{
  FUNCSTART("uint8 SCF::getPeakIndex(uint32 basepos) const");
  BUGIFTHROW(basepos>=SCF_header.bases, "basepos > number of bases");
  FUNCEND();
  return SCF_peak_index[basepos];
}

void SCF::addComment(const char * comment)
{
  FUNCSTART("void SCF::addComment(char * newcomment)");

  if(comment==nullptr) return;
  uint32 newcomlen=SCF_header.comments_size+static_cast<uint32>(strlen(comment)+1);

  char * newcom= new char[newcomlen];

  strcpy(newcom, SCF_comments);
  strcat(newcom, comment);
  delete [] SCF_comments;
  SCF_comments=newcom;
  SCF_header.comments_size=newcomlen-1;

  recalcHeader();

  FUNCEND();
}
