/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
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

//#include <valgrind/callgrind.h>

#include <iostream>
#include <string>

#include <boost/bind.hpp>
#include <boost/crc.hpp>      // for boost::crc_basic, boost::augmented_crc

#include "version.H"

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"
#include "util/dptools.H"
#include "util/misc.H"
#include "io/fastq-mira.H"



using namespace std;

class MiraZip
{
private:
  typedef uint32 readid_t;
  typedef uint16 readlength_t;
  typedef uint64 baithash_t;

  /**********************************************************************
   *
   * for baits
   *
   **********************************************************************/

  struct baitinfo_t {
    baithash_t   bait;
    readid_t     rid;
    readlength_t pos;

    friend ostream & operator<<(ostream &ostr, const baitinfo_t & b){
      ostr << "hash: " << hex << b.bait
	   << "\trid: " << dec << b.rid
	   << "\thpos: " << b.pos;
      return ostr;
    }

  };

  uint32 MZ_baitlength;
  vector<baitinfo_t>   MZ_allbaitinfo;

  /**********************************************************************
   *
   * All reads:
   * reading file: offsets in file for each read
   * basic read info: length
   *
   **********************************************************************/

  vector<streampos>    MZ_readfileoffsets;
  vector<readlength_t> MZ_readlengths;
  vector<uint32> MZ_readcrc32;
  vector<bool>   MZ_readhascrc32;
  vector<uint16> MZ_readcrc16;
  vector<bool>   MZ_readhascrc16;

  /**********************************************************************
   *
   *
   *
   **********************************************************************/

  struct basicread_t {
    string name;
    string comment;
    string seq;
    string quals;
  };

  struct readcache_t {
    vector<basicread_t> cachereads;
    list<size_t>        cachefreelist;     // currently unused elements in cache
    list<size_t>        cacheagelist;      // oldest accessed reads at end, newest front
    list<size_t>        mynullptr;

    // two functions:
    //  a) whether read is cached (non nullptr) and
    //  b) pointer to iterator in Age(!) list so that it can be
    //     quickly moved (phuh)
    vector<list<size_t>::iterator>   cacheageiterators;

    readcache_t(readid_t totalreadnum) {
      cacheageiterators.resize(totalreadnum,mynullptr.end());
      cachereads.resize(100000);
      for(size_t ci=0; ci<cachereads.size(); ++ci) cachefreelist.push_back(ci);
    }
    basicread_t & loadRead(FastQ & fqloader, readid_t rid, streampos offset);
  };



  /**********************************************************************
   *
   * other
   *
   **********************************************************************/

  FastQ     MZ_fqloader;
  ofstream  MZ_fout;


private:
  // functions
  void getFileStats(const string & filename,
		    readid_t & numreads,
		    size_t & totalseqlen,
		    size_t & totalnamlen,
		    size_t & totalcomlen);
  void loadBaits(const string & filename);
  void createAllBaitsForSequence(vector<baitinfo_t> & biv,
				 const char * seq,
				 size_t slen,
				 readid_t rid,
				 const uint32 basesperhash);
  bool createCRC16(const string & seq, uint16 & crc16val) const;
  bool createCRC32(const string & seq, uint32 & crc32val) const;
  void chooseBaits(vector<baitinfo_t> & biv);
  void makeOrder(const string & filename, readid_t numreads);
  bool sortABIByBaitAndPosComp(const baitinfo_t & a, const baitinfo_t & b);
  static bool sortABIByBaitComp(const baitinfo_t & a, const baitinfo_t & b);
  size_t compareBaitInfoStringLengths(readcache_t & rc, const baitinfo_t & a, const baitinfo_t & b);
  void dumpGivenFASTQ(baitinfo_t & bi, ostream & fout);

public:
  MiraZip();
  int main(int argc, char ** argv);
  void usage();
};



//#define CEBUG(bla)   {cout << bla; cout.flush();}
MiraZip::basicread_t & MiraZip::readcache_t::loadRead(FastQ & fqloader, readid_t rid, streampos offset)
{
  FUNCSTART("void loadRead(FastQ & fqloader, readid_t rid, streampos offset)");
  size_t crindex=0;
  bool mustload=true;
  CEBUG("loadRead(...,"<<rid<<","<<offset<<")" << endl);
  if(cacheageiterators[rid]!=mynullptr.end()){
    list<size_t>::iterator tmpit=cacheageiterators[rid];
    CEBUG("ok1"<<endl);
    crindex=*tmpit;
    CEBUG("crindex1: " << crindex << endl);
    CEBUG("ok2"<<endl);
    cacheagelist.splice(cacheagelist.begin(),cacheagelist,tmpit);
    mustload=false;
  }else if(!cachefreelist.empty()){
    list<size_t>::iterator tmpit=cachefreelist.end();
    --tmpit;
    crindex=*tmpit;
    CEBUG("crindex2: " << crindex << endl);
    cacheagelist.splice(cacheagelist.begin(),cachefreelist,tmpit);
    cacheageiterators[rid]=cacheagelist.begin();
  }else if(!cacheagelist.empty()){
    list<size_t>::iterator tmpit=cacheagelist.end();
    --tmpit;
    cacheageiterators[*tmpit]=mynullptr.end();
    crindex=*tmpit;
    CEBUG("crindex3: " << crindex << endl);
    cacheagelist.splice(cacheagelist.begin(),cacheagelist,tmpit);
    cacheageiterators[rid]=cacheagelist.begin();
  }else{
    MIRANOTIFY(Notify::FATAL,"I should never have arrived here!");
  }
  CEBUG("cf: " << cachefreelist.size() << endl);
  CEBUG("ca: " << cacheagelist.size() << endl);
  if(mustload){
    if(fqloader.loadNext(offset) < 0) {
      MIRANOTIFY(Notify::FATAL,"Oooops? Seeking read " << rid << " at streampos " <<  offset << " yielded no result? Probably internal error, but maybe the file has changed in the mean time?");
    }
    cachereads[crindex].name=fqloader.getName();
    cachereads[crindex].comment=fqloader.getComment();
    cachereads[crindex].seq=fqloader.getSequence();
    cachereads[crindex].quals=fqloader.getQuals();
  }
  FUNCEND();

  return cachereads[crindex];
}
#define CEBUG(bla)


MiraZip::MiraZip()
{
  MZ_baitlength=32;
}

void MiraZip::usage()
{
  cout << "mirazip ...\n";
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/
void MiraZip::chooseBaits(vector<baitinfo_t> & biv)
{
  if(!biv.empty()){
    MZ_allbaitinfo.push_back(biv.front());
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/
#define MAXREADSIZEALLOWED 65534

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void MiraZip::createAllBaitsForSequence(vector<baitinfo_t> & biv, const char * seq, size_t slen, readid_t rid, const uint32 basesperhash)
{
  FUNCSTART("void MiraZip::createAllBaitsForSequence(vector<baitinfo_t> & biv, const char * seq, size_t slen, readid_t rid, const uint8 basesperhash)");

  BUGIFTHROW(basesperhash==0, "basesperhash == 0 ?");
  BUGIFTHROW(basesperhash>32, "basesperhash > 32 ?");
  if(slen>MAXREADSIZEALLOWED){
    MIRANOTIFY(Notify::FATAL,"Read " << rid << " is longer than MAXREADSIZEALLOWED (" << MAXREADSIZEALLOWED << ") bases. mirazip cannot handle than, sorry.");
  }

  biv.clear();
  biv.resize(slen);

//  Read::setCoutType(Read::AS_TEXT);
  CEBUG("rid: " << rid << endl);
  CEBUG("seq: " << seq << endl);
  CEBUG("strlen(seq): " << strlen(seq) << endl);
  CEBUG("slen: " << slen << endl);

  baithash_t acthash =0;
  baithash_t hashmask=1;

  // *grml* undefined behaviour of left shift for 64 shifts makes this cludge necessary
  if(basesperhash==32){
    hashmask=0;
  }else{
    hashmask<<=basesperhash*2;
  }
  CEBUG("hash mask: " << hex << hashmask << dec << '\n');
  hashmask--;

  CEBUG("sizeof baithash_t: " << sizeof(baithash_t) << '\n');
  CEBUG("bases per hash: " << static_cast<uint16>(basesperhash) << '\n');
  CEBUG("hash mask: " << hex << hashmask << dec << '\n');

  uint32 goods=0;
  uint32 bads=0;
  uint32 baseok=0;

  char actbase;

  vector<baitinfo_t>::iterator biI=biv.begin();

  for(uint16 seqi=0; seqi<slen; seqi++, seq++){
    acthash<<=2;
    acthash&=hashmask;
    baseok++;

    actbase=static_cast<char>(toupper(*seq));

    CEBUG(seqi << '\t' << actbase << endl);

    switch (actbase) {
    case 'A' : break;
    case 'C' : {
      acthash+=1;
      break;
    }
    case 'G' : {
      acthash+=2;
      break;
    }
    case 'T' : {
      acthash+=3;
      break;
    }
    default : {
      // the IUPAC bases are treated like N and X
      if(dptools::isValidIUPACStarBase(actbase)) {
	// break hash making (which is actually better than behaving
	//  like another character in case of multiple bases with
	//  IUPAC or '*')
	acthash=0;
	baseok=0;
      } else {
	cout << "Unknown base '" << *seq << "' (ASCII " << static_cast<uint16>(*seq) << ") at position " << seqi << " in read " << rid << endl;
	exit(100);
      }
    }
    }

    CEBUG(seqi << ' ' << *seq << ' ');
    if(baseok >= basesperhash) {
      goods++;
      //biv.resize(1+biv.size());
      biI->bait=acthash;
      biI->rid=rid;
      biI->pos=seqi-basesperhash+1;
      CEBUG("saved hash: " << *biI << '\n');
      if(biI->pos > 65000){
	MIRANOTIFY(Notify::INTERNAL,"oooopsi: " << *biI);
      }
      ++biI;
    } else {
      //cout << "Missed hash" << endl;
      if(seqi>=basesperhash) {
	bads++;
      }
    }
    CEBUG('\n');
  }

#ifdef CEBUGFLAG
  CEBUG("goods: " << goods << endl);
  CEBUG("bads: " << bads << endl);
#endif

  biv.resize(goods);
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool MiraZip::createCRC32(const string & seq, uint32 & crc32val) const
{
  crc32val=0;
  if(seq.size() <= sizeof(baithash_t)*4+16) return false;
  const char * sI=seq.data();
  sI+=sizeof(baithash_t)*4;

  boost::crc_32_type result;
  result.process_bytes(sI, seq.size()-(sizeof(baithash_t)*4+16));
  crc32val=result.checksum();

  return true;
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

bool MiraZip::createCRC16(const string & seq, uint16 & crc16val) const
{
  crc16val=0;
  if(seq.size() <= 16) return false;

  const char * sI=seq.data();
  sI+=seq.size()-16;

  boost::crc_16_type result;
  result.process_bytes(sI, 16);
  crc16val=result.checksum();

  return true;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void MiraZip::loadBaits(const string & filename)
{
  FUNCSTART("void MiraZip::loadBaits(const string & filename)");

  readid_t numreads=0;
  vector<baitinfo_t> tmpbi;
  tmpbi.reserve(1000);

  // tell reader to track tellg()
  MZ_fqloader.openFile(filename,true);

  while (MZ_fqloader.loadNext() >= 0) {
    createAllBaitsForSequence(tmpbi,
			      MZ_fqloader.getSequence().c_str(),
			      MZ_fqloader.getSequence().size(),
			      numreads,
			      MZ_baitlength);
    //cout << "tmpbi size: " << tmpbi.size()<<endl;
    chooseBaits(tmpbi);
    ++numreads;
    MZ_readfileoffsets.push_back(MZ_fqloader.getStreamposOfRead());
    MZ_readlengths.push_back(MZ_fqloader.getSequence().size());
    uint32 crc32val=0;
    uint16 crc16val=0;
    MZ_readhascrc32.push_back(createCRC32(MZ_fqloader.getSequence(),crc32val));
    MZ_readcrc32.push_back(crc32val);
    MZ_readhascrc16.push_back(createCRC16(MZ_fqloader.getSequence(),crc16val));
    MZ_readcrc16.push_back(crc16val);
  }

  FUNCEND();
}


/*************************************************************************
 *
 * Careful: stats are not zeroed when entering the function, they rather
 *  add up
 *
 *************************************************************************/

void MiraZip::getFileStats(const string & filename, readid_t & numreads, size_t & totalseqlen, size_t & totalnamlen, size_t & totalcomlen)
{
  FUNCSTART("void MiraZip::getFileStats(string filename, readid_t & numreads, size_t & totalseqlen, size_t & totalcomlen)");

#if 0
  gzFile fp;
  kseq_t *seq;

  fp = gzopen(filename.c_str(), "r");
  if(fp==Z_NULL) {
    MIRANOTIFY(Notify::FATAL,"Could not open FASTQ file '" << filename << "'. Is it present? Is it readable? Did you want to load your data in another format?");
  }
  seq = kseq_init(fp);

  cout << "Counting sequences in file: "; cout.flush();

  // well, count sequences and and total lengths
  int dummy;
  while ((dummy = kseq_read(seq)) >= 0) {
    ++numreads;
    totalseqlen+=seq->seq.l;
    totalnamlen+=seq->name.l;
    totalcomlen+=seq->comment.l;
    cout << "@" << seq->name.s;
    if(seq->comment.l){
      cout << ' ' << seq->comment.s;
    }
    cout << '\t' << seq->seq.s << '\t';
    if(seq->qual.l){
      cout << "+\t" << seq->qual.s << '\t';
    }
    cout << '\n';
  }
#endif

  MZ_fqloader.openFile(filename);

  while (MZ_fqloader.loadNext() >= 0) {
    ++numreads;
    totalseqlen+=MZ_fqloader.getSequence().size();
    totalnamlen+=MZ_fqloader.getName().size();
    totalcomlen+=MZ_fqloader.getComment().size();
    //cout << "\"@" << MZ_fqloader.getName() << "\"";
    //if(MZ_fqloader.getComment().size()){
    //  cout << " \"" << MZ_fqloader.getComment() << "\"";
    //}
    //cout << '\t' << MZ_fqloader.getSequence() << '\t';
    //if(MZ_fqloader.getQuals().size()){
    //  cout << "+\t" << MZ_fqloader.getQuals() << '\t';
    //}
    //cout << '\n';
  }

  cout << "#r: " << numreads
       << "\ts: " << totalseqlen
       << "\tn: " << totalnamlen
       << "\tc: " << totalcomlen
       << endl;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/
void MiraZip::dumpGivenFASTQ(baitinfo_t & bi, ostream & fout)
{
  FUNCSTART("void MiraZip::dumpGivenFASTQ(baitinfo_t & bi)");
  if(MZ_fqloader.loadNext(MZ_readfileoffsets[bi.rid]) < 0) {
    MIRANOTIFY(Notify::FATAL,"Oooops? Seeking read at streampos " <<  MZ_readfileoffsets[bi.rid] << " yielded no result? Probably internal error, but maybe the file has changed in the mean time?");
  }

  fout << "@" << MZ_fqloader.getName();
  if(MZ_fqloader.getComment().size()){
    fout << " " << MZ_fqloader.getComment();
  }
  fout << '\n' << MZ_fqloader.getSequence() << '\n';
  if(MZ_fqloader.getQuals().size()){
    fout << "+\n" << MZ_fqloader.getQuals() << '\n';
  }
  FUNCEND();
}

/*************************************************************************
 *
 * sort by bait first, on equality by position, both ascending
 *
 *************************************************************************/

bool MiraZip::sortABIByBaitAndPosComp(const baitinfo_t & a, const baitinfo_t & b)
{
  if(a.bait==b.bait){
    if(a.pos==b.pos){
      if(MZ_readcrc32[a.rid]==MZ_readcrc32[b.rid]){
	return MZ_readcrc16[a.rid]<MZ_readcrc16[b.rid];
      }
      return MZ_readcrc32[a.rid]<MZ_readcrc32[b.rid];
    }
    return a.pos<b.pos;
  }
  return a.bait<b.bait;
}

bool MiraZip::sortABIByBaitComp(const baitinfo_t & a, const baitinfo_t & b)
{
  return a.bait<b.bait;
}

//bool MiraZip::compareReadEquality(readcache_t & rc,const baitinfo_t & a, const baitinfo_t & b)
//{
//  FUNCSTART("size_t MiraZip::compareReadEquality(readcache_t & rc,const baitinfo_t & a, const baitinfo_t & b)");
//
//  if(a.bait!=b.bait) return false;
//
//  basicread_t & lread=rc.loadRead(MZ_fqloader,a.rid,MZ_readfileoffsets[a.rid]);
//  string aso=lread.seq;
//
//  basicread_t & mread=rc.loadRead(MZ_fqloader,b.rid,MZ_readfileoffsets[b.rid]);
//  string bso=mread.seq;
//}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
size_t MiraZip::compareBaitInfoStringLengths(readcache_t & rc,const baitinfo_t & a, const baitinfo_t & b)
{
  FUNCSTART("size_t MiraZip::compareBaitInfoStringLengths(const baitinfo_t & a, const baitinfo_t & b)");

  basicread_t & lread=rc.loadRead(MZ_fqloader,a.rid,MZ_readfileoffsets[a.rid]);
  string aso=lread.seq;

  basicread_t & mread=rc.loadRead(MZ_fqloader,b.rid,MZ_readfileoffsets[b.rid]);
  string bso=mread.seq;

  string::const_iterator aI=aso.begin();
  advance(aI,a.pos+MZ_baitlength);
  string::const_iterator bI=bso.begin();
  advance(bI,b.pos+MZ_baitlength);

  //cout << "Comp " << a.rid << "," << b.rid << ":\n" << aso << '\n' << bso << endl;

  size_t matchlen=MZ_baitlength;
  for(; aI!=aso.end() && bI!=bso.end() && *aI== *bI; ++aI, ++bI) ++matchlen;

  FUNCEND();
  return matchlen;
}
//#define CEBUG(bla)

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void MiraZip::makeOrder(const string & filename, readid_t numreads)
{
  FUNCSTART("void MiraZip::makeOrder(const string & filename, readid_t numreads)");

  if(MZ_allbaitinfo.empty()) return;

  cout << "now sorting" << endl;
  //sort(MZ_allbaitinfo.begin(),MZ_allbaitinfo.end(),MiraZip::sortABIByBaitAndPosComp);
  sort(MZ_allbaitinfo.begin(),
       MZ_allbaitinfo.end(),
       boost::bind( &MiraZip::sortABIByBaitAndPosComp, this, _1, _2 ));

  MZ_fqloader.openFile(filename);

  vector<uint8> taken(numreads,0);
  readcache_t readcache(numreads);


  int64_t referenceid=-1;
  readid_t numtodump=numreads;
  baitinfo_t refbi=MZ_allbaitinfo.back();
  refbi.rid=-1;

  vector<baitinfo_t> tmpbiv;
  tmpbiv.reserve(1000);

  cout << "dumping" << endl;

  while(numtodump){
    if(refbi.rid==-1){
      cout << "Need new\n";
      while(MZ_allbaitinfo.size() && taken[MZ_allbaitinfo.back().rid]) MZ_allbaitinfo.pop_back();
      refbi=MZ_allbaitinfo.back();
      MZ_allbaitinfo.pop_back();
      dumpGivenFASTQ(refbi,MZ_fout);
      taken[refbi.rid]=1;
      --numtodump;
    }else{
      CEBUG("Have ref " << refbi << endl);

      basicread_t & refread=readcache.loadRead(MZ_fqloader,refbi.rid,MZ_readfileoffsets[refbi.rid]);
      //if(MZ_fqloader.loadNext(MZ_readfileoffsets[refbi.rid]) < 0) {
      //	MIRANOTIFY(Notify::FATAL,"Oooops? Seeking read at streampos " <<  MZ_readfileoffsets[refbi.rid] << " yielded no result? Probably internal error, but maybe the file has changed in the mean time?");
      //}
      createAllBaitsForSequence(tmpbiv,
				refread.seq.c_str(),
				refread.seq.size(),
				refbi.rid,
				MZ_baitlength);

      uint64_t allchecks=0;
      uint64_t longest=0;
      baitinfo_t longestbi=refbi;
      longestbi.rid=-1;
      vector<baitinfo_t>::iterator actbiI=tmpbiv.begin();
      for(; actbiI != tmpbiv.end(); ++actbiI){
	CEBUG("Check1: " << *actbiI << endl);
	// if longest we found is longer than any other possibility, we do not need to check further
	if(actbiI->pos + longest >= MZ_readlengths[actbiI->rid]) break;
	vector<baitinfo_t>::iterator lbiI=lower_bound(MZ_allbaitinfo.begin(),MZ_allbaitinfo.end(),*actbiI,MiraZip::sortABIByBaitComp);
	if(lbiI!=MZ_allbaitinfo.end()){
	  CEBUG("Check2: " << *actbiI << endl);
	  uint64_t innerchecks=0;
	  for(; lbiI != MZ_allbaitinfo.end() && lbiI->bait == actbiI->bait; ++lbiI){
	    CEBUG("Check3: " << *actbiI << "\t" << *lbiI << endl);
	    if(lbiI->rid != actbiI->rid
	       && taken[lbiI->rid]==0
	       && MZ_readlengths[lbiI->rid]-lbiI->pos > longest){
	      int64 thiscmp=compareBaitInfoStringLengths(readcache,*actbiI,*lbiI);
	      ++allchecks;
	      ++innerchecks;
	      CEBUG("Check4 (" << allchecks << ","<<innerchecks<<")\t"<< thiscmp << '\t' << longest);
	      if(thiscmp>longest){
		CEBUG(" is longest");
		longest=thiscmp;
		longestbi=*lbiI;
		innerchecks=0;
	      }
	      CEBUG("\n");
	      // stop runaway comaprisons which did not yield better result for some time
	      if(innerchecks>=200) break;
	    }
	  }
	}
      }
      CEBUG("I'm at end" << endl);
      if(longestbi.rid!=-1){
	refbi=longestbi;
	CEBUG("Take this (" << longest << ","<<allchecks<<"): " << longestbi<< endl);
	dumpGivenFASTQ(refbi,MZ_fout);
	taken[refbi.rid]=1;
	--numtodump;
      }else{
	refbi.rid=-1;
      }
    }
  }

  vector<baitinfo_t>::const_iterator biI=MZ_allbaitinfo.begin();
  for(; biI != MZ_allbaitinfo.end(); ++biI){
    //cout << biI->rid << " seeking " << MZ_readfileoffsets[biI->rid] << endl;
    if(MZ_fqloader.loadNext(MZ_readfileoffsets[biI->rid]) < 0) {
      MIRANOTIFY(Notify::FATAL,"Oooops? Seeking read at streampos " <<  MZ_readfileoffsets[biI->rid] << " yielded no result? Probably internal error, but maybe the file has changed in the mean time?");
    }

  }

  FUNCEND();
}
#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

int MiraZip::main(int argc, char ** argv)
{
  readid_t numreads=0;
  size_t totalseqlen=0;
  size_t totalnamlen=0;
  size_t totalcomlen=0;

  string filename("mziptest.fastq");

  //string filename("bs4511.fastq");

  getFileStats(filename,numreads,totalseqlen,totalnamlen,totalcomlen);

  MZ_allbaitinfo.clear();
  MZ_allbaitinfo.reserve(numreads);
  MZ_readfileoffsets.clear();
  MZ_readfileoffsets.reserve(numreads+1);
  MZ_readlengths.clear();
  MZ_readlengths.reserve(numreads+1);
  MZ_readcrc16.clear();
  MZ_readcrc16.reserve(numreads+1);
  MZ_readcrc32.clear();
  MZ_readcrc32.reserve(numreads+1);
  MZ_readhascrc16.clear();
  MZ_readhascrc16.reserve(numreads+1);
  MZ_readhascrc32.clear();
  MZ_readhascrc32.reserve(numreads+1);

  loadBaits(filename);

  MZ_fout.open("sorted", ios::out);
  makeOrder(filename,numreads);
  //dumpReads();

  return 0;
}


int main(int argc, char ** argv)
{
  //CALLGRIND_STOP_INSTRUMENTATION;

  FUNCSTART("int main(int argc, char ** argv)");

  string path;
  string convertprog;
  splitFullPathAndFileName(argv[0],path,convertprog);

  std::transform(convertprog.begin(),
		 convertprog.end(),
		 convertprog.begin(),
		 (int(*)(int))std::tolower); // now, that's what I call ugly

  try {
    if(convertprog=="mirazip"){
      MiraZip mz;
      mz.main(argc, argv);
    }
  }
  catch(Notify n){
    n.handleError("main");
  }
  catch(Flow f){
    cout << "INTERNAL ERROR: Unexpected exception: Flow()\n";
    exit(100);
  }
  catch(const std::bad_alloc & e){
    cout << "Out of memory detected, exception message is: ";
    cout << e.what() << endl;

    if(sizeof(size_t) == sizeof(int32)){
      cout << "\nYou are running a 32 bit executable. Please note that the maximum"
	"\ntheoretical memory a 32 bit programm can use (be it in Linux, Windows or"
	"\nother) is 4 GiB, in practice less: between 2.7 and 3.3 GiB. This is valid"
	"\neven if your machine has hundreds of GiB."
	"\nShould your machine have more that 4 GiB, use a 64 bit OS and a 64 bit"
	"\nversion of MIRA.";
    }

    cout << "\n\nIf you have questions on why this happened, please send the last 1000"
      "\nlines of the output log (or better: the complete file) to the author"
      "\ntogether with a short summary of your assembly project.\n\n";

    exit(100);
  }
  catch(const ios_base::failure & e){
    cout << "Failure in IO stream detected, exception message is: "
	 << e.what() << endl
	 << "\nWe perhaps ran out of disk space or hit a disk quota?\n";
    exit(100);
  }
  catch (exception& e)
  {
    cout << "A 'standard' exception occurred (that's NOT normal):\n" << e.what() << "\n\nIf the cause is not immediatly obvious, please contact: bach@chevreux.org\n\n";
    exit(100);
  }
  catch(...){
    cout << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    exit(100);
  }

  return 0;
}
