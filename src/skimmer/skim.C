/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 and later by Bastien Chevreux
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 *
 */

// 	$Id$

#ifndef lint
static char vcid[] = "$Id$";
#endif /* lint */

#include "skim.H"





//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif




#define MAXVHASHMASK 0xFFFFFFL





//hier dann weiter wieso kommt bei partitionen was anderes raus als ohne?



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Skim::foolCompiler()
{
#include "stdinc/foolcompiler.C"
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

// Plain vanilla constructor
Skim::Skim()
{
  FUNCSTART("Skim::Skim(ReadPool & rp)");

  init();

  FUNCEND();
}


void Skim::init()
{
  FUNCSTART("Skim::init()");

  SKIM3_possiblehits=0;
  SKIM3_acceptedhits=0;


  SKIM3_basesperhash=16;
  SKIM3_hashsavestepping=4;
  SKIM3_percentrequired=50;

  SKIM_totalpermbans=0;
  SKIM3_totalhitschosen=0;

  FUNCEND()
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Skim::~Skim()
{
  FUNCSTART("Skim::~Skim()");
  //  ERROR("Not implemented yet.");

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//Skim::Skim(Skim const &other)
//{
//  FUNCSTART("Skim::Skim(Skim const &other)");
//
//  SKIM_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//
//// Copy operator, needed by copy-constructor
//Skim const & Skim::operator=(Skim const & other)
//{
//  FUNCSTART("Skim const & Skim::operator=(Skim const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}




///*************************************************************************
// *
// *
// *
// *
// *************************************************************************/
//
//ostream & operator<<(ostream &ostr, Skim const &theskim)
//{
//  FUNCSTART("friend ostream & Skim::operator<<(ostream &ostr, const  &theskim)");
//  //  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Skim::discard()
{
  FUNCSTART("Skim::discard()");

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}

void Skim::skimGo(ReadPool & rp,
		  string               & posfmatchname,
		  string               & poscmatchname,
		  bannedoverlappairs_t & bannedoverlaps,
		  vector<uint32>       & overlapcounter,
		  vector<uint32>       & rawhashhitcounter,
		  uint32 maxmemusage,
		  bool onlyagainstrails,
		  uint8  bph,
		  uint8  hss,
		  int32  percentrequired,
		  uint32 maxhitsperread)
{
  dateStamp(cout);
  Read::setCoutType(Read::AS_CLIPPEDFASTA);

  init();

  SKIM3_readpool=&rp;

  // if extalso is true, the take not only the clipped sequence, but
  //  also the extends to it

  SKIM3_onlyagainstrails=onlyagainstrails;

  SKIM3_basesperhash=bph;
  SKIM3_hashsavestepping=hss;
  SKIM3_percentrequired=percentrequired;
  SKIM3_maxhitsperread=maxhitsperread;

  SKIM3_overlapcounter=&overlapcounter;
  SKIM3_rawhashitcounter=&rawhashhitcounter;
  SKIM3_bannedoverlaps=&bannedoverlaps;

  SKIM_partfirstreadid=0;         // partition first read id
  SKIM_partlastreadid=0;         // partition last read id

  SKIM3_posfmatchfout.open(posfmatchname.c_str(), ios::out| ios::trunc);
  SKIM3_poscmatchfout.open(poscmatchname.c_str(), ios::out| ios::trunc);

  uint32 numpartitions=computePartition(maxmemusage*SKIM3_hashsavestepping,true);

  CEBUG("We will get " << numpartitions << " partitions.\n");

  //cout << "\nSKIMMER: Using maximum of " << maxmemusage << " hashes stored in memory, " << numpartitions << " partitions will be computed." << endl << endl;

  CEBUG("Progressend: " << SKIM_progressend << endl);

  SKIM3_megahubs.resize(SKIM3_readpool->size(),0);

  cout << "Now running partitioned skimmer with " << numpartitions << " partitions:" << endl;

  SKIM_progressindicator= new ProgressIndicator(0,SKIM_progressend);

  for(uint32 actpartition=1; actpartition<=numpartitions; actpartition++){
    CEBUG("\nWorking on partition " << actpartition << "/" << numpartitions << endl);

    computePartition(maxmemusage*SKIM3_hashsavestepping,false);

    CEBUG("Will contain read IDs " << SKIM_partfirstreadid << " to " << SKIM_partlastreadid-1 << endl);

    prepareSkim();

    SKIM_partfirstreadid=SKIM_partlastreadid;
  }

  SKIM_progressindicator->finishAtOnce();

  cout << "\n";

  SKIM3_posfmatchfout.close();
  SKIM3_poscmatchfout.close();

  {
    uint32 megahubs=0;
    for(uint32 i=0; i<SKIM3_megahubs.size(); i++){
      if(SKIM3_megahubs[i]>0) {
	megahubs++;
	cout << "Megahub: " << i << '\n';
      }
    }
    cout << "Total megahubs: " << megahubs << endl;
  }

  //cout << "\nSkim summary:\n\texcellent: " << SKIM_totalexcellenthits << "\n\te-rejected: " << SKIM_totalexcellentrejectedhits << "\n\tgood: " << SKIM_totalgoodhits << "\n\tpossible: " << SKIM_totalpossiblehits << "\n\tpermbans: " << SKIM_totalpermbans << "\n\tredohash: " << SKIM_totalredohash;
  cout << "\nSkim summary:\n\taccepted: " << SKIM3_acceptedhits << "\n\tpossible: " << SKIM3_possiblehits  << "\n\tpermbans: " << SKIM_totalpermbans;
  cout << "\n\nHits chosen: " << SKIM3_totalhitschosen << "\n\n";

  dateStamp(cout);

  cout << endl;

  delete SKIM_progressindicator;

}
#define CEBUG(bla)
#define CEBUGF(bla)




/*************************************************************************
 *
 * computes either
 *  - starting from the first read id, the last read id to use for the
 *    next partition (SKIM_partfirstreadid and SKIM_partlastreadid)
 * or
 *  - starting from the first read id, the total number of partitions
 * to use in the next skim run.(returned)
 *
 * The maxmemusage is the dominating factor: each partition may not use
 *  (much) more than that.
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

uint32 Skim::computePartition(uint32 maxmemusage, bool computenumpartitions)
{
  FUNCSTART("uint32 Skim::computePartition(uint32 maxmemusage, bool computenumpartitions)");

  uint32 numpartitions=0;
  uint32 totalseqlen=0;
  uint32 maxseqlen=0;

  SKIM_partlastreadid=SKIM_partfirstreadid;

  if(computenumpartitions) {
    SKIM_progressend=SKIM3_readpool->size()-SKIM_partlastreadid;
  }

  for(; SKIM_partlastreadid<SKIM3_readpool->size(); SKIM_partlastreadid++) {
    if(!SKIM3_readpool->getRead(SKIM_partlastreadid).hasValidData()
       || !SKIM3_readpool->getRead(SKIM_partlastreadid).isUsedInAssembly()) continue;

    if(SKIM3_readpool->getRead(SKIM_partlastreadid).getLenClippedSeq() >29900) {
      NOTIFYMSG("Read " << SKIM3_readpool->getRead(SKIM_partlastreadid).getName() << " is longer than 29900 bases. SKIM cannot handle this, aborting.\n");
      throw Notify(Notify::FATAL, THISFUNC, notifymsg) ;
    }

    maxseqlen=max(maxseqlen,SKIM3_readpool->getRead(SKIM_partlastreadid).getLenClippedSeq());

    totalseqlen+=SKIM3_readpool->getRead(SKIM_partlastreadid).getLenClippedSeq();
    //if(SKIM_takeextalso){
    //  totalseqlen+=SKIM3_readpool->getRead(SKIM_partlastreadid).getRightExtend();
    //}

    if(totalseqlen>maxmemusage) {
      if(computenumpartitions){
	totalseqlen=0;
	numpartitions++;

	SKIM_progressend+=SKIM3_readpool->size()-SKIM_partlastreadid;

      }else{
	SKIM_partlastreadid++;
	break;
      }
    }
  }

  if(computenumpartitions) {
    if(totalseqlen>0) {
      numpartitions++;
      SKIM_progressend+=SKIM3_readpool->size()-SKIM_partlastreadid;
    }
    SKIM_progressend*=2;
  }

  //// Compute SKIM_maxoffsets
  //if(maxseqlen>0){
  //  // this beauty is slow, but a cute hack to find out the number of
  //  //  the highest bit which power of two fits the given value
  //  // e.g. 1024 -> 10 -> 2^10 = 1024 fits 1024
  //  //      1025 -> 11 -> 2^11 = 2048 fits 1025
  //  SKIM_mo_shiftmultiplier=0;
  //  while((maxseqlen-1) >> ++SKIM_mo_shiftmultiplier);
  //}else{
  //  SKIM_mo_shiftmultiplier=4;
  //}
  //
  //// restrict shift multiplier to 11 (and therefore maxoffsets to 2048)
  //if(SKIM_mo_shiftmultiplier>11) SKIM_mo_shiftmultiplier=11;
  //SKIM_maxoffsets=(1<<SKIM_mo_shiftmultiplier);

  FUNCEND();

  return numpartitions;
}





/*************************************************************************
 *
 * sorter to sort from low to high, but lower 24bit grouped
 *
 *
 *************************************************************************/

inline bool Skim__sortVHRAPArray_(const vhrap_t & a,
			    const vhrap_t & b);
inline bool Skim__sortVHRAPArray_(const vhrap_t & a, const vhrap_t & b)
{
  ////if(a.vhash == b.vhash) {
  ////  return a.readid < b.readid;
  ////}
  //return a.vhash < b.vhash;

  if((a.vhash & 0xFFFFFFL) != (b.vhash & 0xFFFFFFL)) {
    return (a.vhash & 0xFFFFFFL) < (b.vhash & 0xFFFFFFL);
  }
  return a.vhash < b.vhash;
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

void Skim::prepareSkim()
{
  FUNCSTART("void Skim::prepareSkim()");

  uint32 totalseqlen=0;
  uint32 totalseqs=0;

  for(uint32 i=SKIM_partfirstreadid; i<SKIM_partlastreadid; i++) {
    if(!SKIM3_readpool->getRead(i).hasValidData()
      || !SKIM3_readpool->getRead(i).isUsedInAssembly()) continue;
    totalseqlen+=SKIM3_readpool->getRead(i).getLenClippedSeq();
    totalseqs++;
    //if(SKIM_takeextalso) totalseqlen+=SKIM3_readpool->getRead(i).getRightExtend();
  }

  //dateStamp(cout);
  CEBUG("\nPreparing skim data: "  << SKIM_partfirstreadid << " to " << SKIM_partlastreadid << endl);
  CEBUG(totalseqs << " sequences to skim, totalling " << totalseqlen << " bases." << endl);


  // next loop:
  //  transform each read into a series of forward
  //   hashes, store them into the array along with info whether
  //   each hash position is valid or not,
  //  also store the info how many hashes each sequence produced

  uint32 totalhashes=0;

  vector<vhrap_t> vhraparray(totalseqlen/SKIM3_hashsavestepping);
  vector<vhrap_t>::iterator vhraparrayI=vhraparray.begin();

  //ProgressIndicator P(partfirstreadid, partlastreadid);
  for(uint32 seqnr=SKIM_partfirstreadid; seqnr<SKIM_partlastreadid; seqnr++){
    //P.progress(seqnr);
    Read & actread= SKIM3_readpool->getRead(seqnr);
    if(!actread.hasValidData()
      || !actread.isUsedInAssembly()) continue;

    uint32 slen=actread.getLenClippedSeq();
    //if(SKIM_takeextalso) slen+=actread.getRightExtend();

    if(slen>=8) {
      uint32 hashesmade= transformSeqToVariableHash(
	seqnr,
	actread.getClippedSeqAsChar(),
	slen,
	actread.getName(),
	SKIM3_basesperhash,
	vhraparrayI,
	false,
	SKIM3_hashsavestepping
	);

      //CEBUG(seqnr << "\t" << totalhashes << "\t" << slen << endl);
      totalhashes+=hashesmade;
    }
  }

  //P.progress(partlastreadid);
  CEBUG("Totalseqlen " << totalseqlen << endl);
  CEBUG("Computed " << totalhashes << " linkpoints." << endl);

  //{
  //  CEBUG("Partition unsorted:\n");
  //  cout << "Partition unsorted:\n";
  //  vector<vhrap_t>::const_iterator vaI=vhraparray.begin();
  //  for(; vaI!=vhraparray.end(); vaI++){
  //    cout << "###\t" << hex << vaI->vhash << dec << '\t' << vaI->readid << '\t' << vaI->hashpos << '\n';
  //  }
  //  cout << "###########################\n";
  //}

  CEBUG("Resizing array" << endl);
  vhraparray.resize(totalhashes);

  CEBUG("Sorting array" << endl);
  sort(vhraparray.begin(), vhraparray.end(), Skim__sortVHRAPArray_);

  //{
  //  CEBUG("Partition sorted:\n");
  //  vector<vhrap_t>::const_iterator vaI=vhraparray.begin();
  //  for(; vaI!=vhraparray.end(); vaI++){
  //    cout << "###\t" << hex << vaI->vhash << dec << '\t' << vaI->readid << '\t' << vaI->hashpos << '\n';
  //  }
  //  cout << "###########################" << endl;
  //}

  makeVHRAPArrayShortcuts(vhraparray, SKIM3_basesperhash);
  checkForHashes(1);
  checkForHashes(-1);

  FUNCEND();
}
//#define CEBUG(bla)
//#define CEBUGF(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

void Skim::makeVHRAPArrayShortcuts(vector<vhrap_t> & vhraparray, const uint8 basesperhash)
{
  //cout << "Making VHRAPArrayShortcuts" << endl;

  SKIM3_vashortcuts_begin.clear();
  SKIM3_vashortcuts_end.clear();
  vector<vhrap_t>::const_iterator vaI=vhraparray.begin();
  if(vaI==vhraparray.end()) return;

  SKIM3_vashortcuts_begin.resize(
    1<<(min(static_cast<uint8>(12),basesperhash)*2),
    static_cast<vector<vhrap_t>::iterator>(nullptr)
    );

  SKIM3_vashortcuts_end.resize(
    1<<(min(static_cast<uint8>(12),basesperhash)*2),
    static_cast<vector<vhrap_t>::iterator>(nullptr)
    );

  vhash_t acthash= (vaI->vhash & MAXVHASHMASK);
  while(vaI != vhraparray.end()){
    SKIM3_vashortcuts_begin[acthash]=vaI;
    for(;vaI != vhraparray.end() && (vaI->vhash & MAXVHASHMASK) == acthash; vaI++);
    SKIM3_vashortcuts_end[acthash]=vaI;
    //cout << "vhash: " << hex << acthash << "\t" << dec << SKIM3_vashortcuts_end[acthash]-SKIM3_vashortcuts_begin[acthash] << '\n';
    if(vaI != vhraparray.end()) acthash= vaI->vhash & MAXVHASHMASK;
  }
}

//#define CEBUG(bla)
//#define CEBUGF(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}


uint32 Skim::transformSeqToVariableHash (const uint32 readid, const char * seq, uint32 slen, const string & readname, const uint8 basesperhash, vector<vhrap_t>::iterator & vhraparrayI, const bool countonly, const uint8 hashsavestepping)
{
  FUNCSTART("void Skim::transformSeqToVariableHash (const char * seq, uint32 slen, const string & readname, const uint8 basesperhash)");


  BUGIFTHROW(basesperhash>32, "basesperhash > 32 ?");
  BUGIFTHROW(hashsavestepping<1, "hashsavestepping < 1 ?");

  vhash_t acthash =0;
  vhash_t hashmask=1;
  hashmask<<=(basesperhash*2);
  hashmask--;

  CEBUG("sizeof vhash_t: " << sizeof(vhash_t) << '\n');
  CEBUG("bases per hash: " << (uint16) basesperhash << '\n');
  CEBUG("hash mask: " << hex << hashmask << dec << '\n');

  // first hash made must also be saved
  uint8 hashsavecounter=1;

  uint32 goods=0;
  uint32 bads=0;
  vector<vhrap_t>::iterator initial_vaI=vhraparrayI;

  CEBUG("Hashing " << readname << "\n");

  uint32  baseok=0;

  for(uint16 i=0; i<slen; i++, seq++){
    acthash<<=2;
    acthash&=hashmask;
    baseok++;

    char actbase=toupper(*seq);

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
      if(dptools::isValidIUPACStarBase(actbase)) {
      // the IUPAC bases are treated like N and X

      // do nothing: behave like an 'a'   (wrong in 75%, but right in 25%
      //  and doesn't break the hash making!)
      // feel free to make addhash+1, +2 or +3 for c, g ot t instead

      // or

      // uncomment this to break hash making (which is actually better,
      //  in case of multiple bases with 'n' or '*')
	acthash=0;
	baseok=0;

	// set hashsavecounter to 1 so that the next good hash generated is saved!
	hashsavecounter=1;
      } else {
	cout << "Unknown base '" << *seq << "' (ASCII " << static_cast<uint16>(*seq) << ") at position " << i << " in _CLIPPED_ sequence " << readname << endl;
	exit(100);
      }
    }
    }

    CEBUG(*seq << " ");
    if(baseok >= basesperhash) {
      CEBUG("Hash: " << hex << acthash << dec << endl);
      goods++;
      if(!countonly) {
	if(--hashsavecounter == 0){
	  vhraparrayI->vhash=acthash;
	  vhraparrayI->readid=readid;
	  vhraparrayI->hashpos=i;
	  vhraparrayI++;

	  hashsavecounter=hashsavestepping;
	}
      }

    } else {
      //cout << "Missed hash" << endl;
      if(i>=basesperhash) {
	bads++;
      }

    }
  }

  //for(uint32 i=0; i<basesperhash; i++, hashp++, hashokp++) {
  //  *hashp=0;
  //  *hashokp=0;
  //}

#ifdef CEBUGFLAG
  //CEBUG(readname << " generated " << hashokp-orghashokp << " hashes.\n");
  CEBUG("goods: " << goods << endl);
  CEBUG("bads: " << bads << endl);
#endif

  return (vhraparrayI-initial_vaI);
}



/*************************************************************************
 *
 * sorter to sort from low to high
 *
 *
 *************************************************************************/

inline bool Skim__compareVHRAPArrayElem_(const vhrap_t & one, const vhrap_t & other)
{
  return one.vhash < other.vhash;
};


bool Skim__sortreadhashmatch_t_(const readhashmatch_t & a,
			    const readhashmatch_t & b);
bool Skim__sortreadhashmatch_t_(const readhashmatch_t & a, const readhashmatch_t & b)
{
  if(a.rid2 == b.rid2){
    if(a.eoffset == b.eoffset) return a.hashpos1 < b.hashpos1;
    return a.eoffset < b.eoffset;
  }
  return a.rid2 < b.rid2;
}

/*************************************************************************
 *
 * sorter to sort from high to low
 *
 *
 *************************************************************************/

bool Skim__sortMWByPercent_(const matchwithsorter_t & a,
			    const matchwithsorter_t & b);
bool Skim__sortMWByPercent_(const matchwithsorter_t & a, const matchwithsorter_t & b)
{
  if(a.percent_in_overlap == b.percent_in_overlap) {
    return a.numhashes > b.numhashes;
  }
  return a.percent_in_overlap > b.percent_in_overlap;
}

bool Skim__sortMWByNumHashes_(const matchwithsorter_t & a,
			    const matchwithsorter_t & b);
bool Skim__sortMWByNumHashes_(const matchwithsorter_t & a, const matchwithsorter_t & b)
{
  if(a.numhashes == b.numhashes){
     return a.percent_in_overlap > b.percent_in_overlap;
  }
  return a.numhashes > b.numhashes;
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

void Skim::checkForHashes(const int8 direction)
{
  uint32 allseqnum=SKIM3_readpool->size();
  uint32 reducedseqnum=SKIM_partlastreadid-SKIM_partfirstreadid;

  vector<vhrap_t> singlereadvhraparray;
  singlereadvhraparray.reserve(5000);
  SKIM3_readhashmatches.clear();

  vector<matchwithsorter_t> tmpmatchwith;
  tmpmatchwith.reserve(2000);

  CEBUG("Checking " << reducedseqnum << " sequences\n");

  //possible_overlaps_t * posmatch=SKIM3_posfmatch;
  //if(direction<0) posmatch=SKIM3_poscmatch;
  ofstream * posmatchfout=&SKIM3_posfmatchfout;
  if(direction<0) posmatchfout=&SKIM3_poscmatchfout;

  for(uint32 actreadid=SKIM_partfirstreadid; actreadid<allseqnum; actreadid++){
    SKIM_progressindicator->increaseprogress();
    //P.progress(actreadid);

    //if(actreadid>100) return;

    // don't need to go through identified megahubs again
    if(SKIM3_megahubs[actreadid]>0) continue;

    Read & actread= SKIM3_readpool->getRead(actreadid);
    if(!actread.hasValidData()
      || !actread.isUsedInAssembly()) continue;

    uint32 slen=actread.getLenClippedSeq();

    if(slen<SKIM3_basesperhash) continue;

    singlereadvhraparray.resize(slen);

    vector<vhrap_t>::iterator srvaI=singlereadvhraparray.begin();
    uint32 hashesmade;
    if(direction>0) {
      hashesmade=transformSeqToVariableHash(
	actreadid,
	actread.getClippedSeqAsChar(),
	slen,
	actread.getName(),
	SKIM3_basesperhash,
	srvaI,
	false,
	1
	);
    }else{
      hashesmade=transformSeqToVariableHash(
	actreadid,
	actread.getClippedComplementSeqAsChar(),
	slen,
	actread.getName(),
	SKIM3_basesperhash,
	srvaI,
	false,
	1
	);
    }

    singlereadvhraparray.resize(hashesmade);

    srvaI=singlereadvhraparray.begin();
    vector<vhrap_t>::const_iterator lowerbound;
    vector<vhrap_t>::const_iterator upperbound;
    uint32 truetestsm2hits=0;
    for(; srvaI != singlereadvhraparray.end(); srvaI++){
      lowerbound=SKIM3_vashortcuts_begin[srvaI->vhash & MAXVHASHMASK];
      upperbound=SKIM3_vashortcuts_end[srvaI->vhash & MAXVHASHMASK];

      if(static_cast<vector<vhrap_t>::const_iterator>(nullptr)!=lowerbound){
	if(SKIM3_basesperhash>12){
	  // with more than 12 bases in a hash, the vhrap array is
	  //  subdivided
	  pair<vector<vhrap_t>::const_iterator, vector<vhrap_t>::const_iterator>
	    p=equal_range(lowerbound,
			  upperbound,
			  *srvaI,
			  Skim__compareVHRAPArrayElem_);
	  lowerbound=p.first;
	  upperbound=p.second;
	}

	for(;lowerbound!=upperbound; lowerbound++){
	  truetestsm2hits++;

//seufz irgendwie noch nicht ausgereift

	  if(actreadid > lowerbound->readid){
	    SKIM3_readhashmatches.resize(SKIM3_readhashmatches.size()+1);
	    SKIM3_readhashmatches.back().rid2=lowerbound->readid;
	    SKIM3_readhashmatches.back().eoffset=srvaI->hashpos - lowerbound->hashpos;
	    SKIM3_readhashmatches.back().hashpos1=srvaI->hashpos;
	    //SKIM3_readhashmatches.back().hashpos2=lowerbound->hashpos;
	  }
	}
      }
    }

    if(actreadid % 1 == 0) CEBUG("actreadid: " << actreadid << "\tSKIM3_readhashmatches.size(): " << SKIM3_readhashmatches.size() << "\ttruetestsm2hits: " << truetestsm2hits << endl);

    // Hmmm ... this does not represent the full truth, but without "if" there
    // is a partition effect in the data. Bad, but cannot be helped (except
    // going through all reads in all partitions, which is unnecessary for the
    // search itself and effectively.doubles the SKIM time)
    if(SKIM_partfirstreadid==0)(*SKIM3_rawhashitcounter)[actreadid]+=truetestsm2hits;

    if(SKIM3_readhashmatches.size()>0){
      if(truetestsm2hits<100000) {
	checkForPotentialHits(actreadid, tmpmatchwith);


	if(tmpmatchwith.size()<SKIM3_maxhitsperread) {
	  SKIM3_totalhitschosen+=tmpmatchwith.size();
	  for(uint32 i=0; i<tmpmatchwith.size(); i++) {
	    *posmatchfout << tmpmatchwith[i].otherid
			  << '\t' << actreadid
			  << '\t' << -tmpmatchwith[i].eoffset
			  << '\t' << tmpmatchwith[i].percent_in_overlap
			  << '\t' << tmpmatchwith[i].numhashes
			  << '\n';
	  }
	} else {
	  uint32 incounter=0;

	  // take best matches by percent
	  sort(tmpmatchwith.begin(), tmpmatchwith.end(), Skim__sortMWByPercent_);
	  vector<matchwithsorter_t>::iterator I=tmpmatchwith.begin();
	  for(;incounter<SKIM3_maxhitsperread/2 && I!=tmpmatchwith.end(); I++) {
	    if(!(I->taken)){
	      I->taken=true;
	      *posmatchfout << I->otherid
			    << '\t' << actreadid
			    << '\t' << -I->eoffset
			    << '\t' << I->percent_in_overlap
			    << '\t' << I->numhashes
			    << '\n';

	      SKIM3_totalhitschosen++;
	      incounter++;
	    }
	  }

	  // take best matches by num hashes
	  sort(tmpmatchwith.begin(), tmpmatchwith.end(), Skim__sortMWByNumHashes_);
	  I=tmpmatchwith.begin();
	  for(;incounter<SKIM3_maxhitsperread && I!=tmpmatchwith.end(); I++) {
	    if(!(I->taken)){
	      I->taken=true;
	      *posmatchfout << I->otherid
			    << '\t' << actreadid
			    << '\t' << -I->eoffset
			    << '\t' << I->percent_in_overlap
			    << '\t' << I->numhashes
			    << '\n';

	      SKIM3_totalhitschosen++;
	      incounter++;
	    }
	  }
	}


      }else{
	CEBUG("Identified megahub: " << actreadid << "\tSKIM3_readhashmatches.size(): " << SKIM3_readhashmatches.size() << endl);
	SKIM3_megahubs[actreadid]=1;
      }
      SKIM3_readhashmatches.clear();
    }
  }
}

//#define CEBUG(bla)
//#define CEBUGF(bla)







//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUGF(bla)  {cout << bla; cout.flush();}

void Skim::checkForPotentialHits(const uint32 actreadid, vector<matchwithsorter_t> & tmpmatchwith)
{
  CEBUG("Potential hits of " << actreadid << " (" << SKIM3_readpool->getRead(actreadid).getLenClippedSeq() << ")\n----------------\n");
  CEBUG(SKIM3_readpool->getRead(actreadid) << endl);

  tmpmatchwith.clear();

  sort(SKIM3_readhashmatches.begin(), SKIM3_readhashmatches.end(), Skim__sortreadhashmatch_t_);

  bool actreadisrail=SKIM3_readpool->getRead(actreadid).isRail();

  vector<readhashmatch_t>::const_iterator sI=SKIM3_readhashmatches.begin();
  uint32 possiblehits=0;
  uint32 acceptedhits=0;

  uint32 countid=sI->rid2;
  while(sI != SKIM3_readhashmatches.end()){
    possiblehits++;

    uint32 rid2=sI->rid2;
    uint16 hp1min=0xffff;
    uint16 hp1max=0;
    //uint16 hp2min=0xffff;
    //uint16 hp2max=0;
    int32  eoffsetmin=0x7fffffff;
    int32  eoffsetmax=0x80000000;
    int32 oldeoffset=sI->eoffset;

    uint32 numhashes=0;

    for(;sI != SKIM3_readhashmatches.end() && sI->rid2 == countid; sI++){
      // this ensures that the eoffset between two following
      //  entries may not differ by too much (10 bases here)
      // IF they do, then this is treated like a different hit
      //  by breaking the loop
      if(abs(sI->eoffset - oldeoffset) > 10){
	CEBUG("BREAKER!\n");
	break;
      }
      numhashes++;
      hp1min=min(hp1min,sI->hashpos1);
      hp1max=max(hp1max,sI->hashpos1);
      //hp2min=min(hp2min,sI->hashpos2);
      //hp2max=max(hp2max,sI->hashpos2);
      eoffsetmin=min(eoffsetmin,sI->eoffset);
      eoffsetmax=max(eoffsetmax,sI->eoffset);
      oldeoffset=sI->eoffset;

      CEBUG(sI->rid2);
      CEBUG("\t" << SKIM3_readpool->getRead(sI->rid2).getLenClippedSeq());
      CEBUG("\t" << sI->eoffset);
      CEBUG("\t" << sI->hashpos1);
      //CEBUG("\t" << sI->hashpos2);
      CEBUG('\n');
    }

    // disregard this potential match if
    //  1) we scan only against rails and both reads are non-rails
    //  2) this read pair has been banned previously
    if((SKIM3_onlyagainstrails
	&& (!actreadisrail && !SKIM3_readpool->getRead(rid2).isRail()))
       ||
       (*SKIM3_bannedoverlaps)[actreadid].count(rid2) > 0){
      if(sI!=SKIM3_readhashmatches.end()) countid=sI->rid2;
      continue;
    }

    int32 maxoverlap;

    // adjust min positions for the hash length
    hp1min-=(SKIM3_basesperhash-1);
    //hp2min-=(SKIM3_basesperhash-1);

    int32 eoffsetmean=eoffsetmin+(eoffsetmax-eoffsetmin)/2;

    // calc max overlap
    // currently only for one offset
    if(eoffsetmean<0){
      maxoverlap=min(SKIM3_readpool->getRead(rid2).getLenClippedSeq()+eoffsetmean,SKIM3_readpool->getRead(actreadid).getLenClippedSeq());
    }else{
      maxoverlap=min(SKIM3_readpool->getRead(actreadid).getLenClippedSeq()-eoffsetmean,SKIM3_readpool->getRead(rid2).getLenClippedSeq());
    }

    int32 hashesoverlap=hp1max-hp1min;
    int32 perc=100*hashesoverlap/maxoverlap;

    if(perc>=SKIM3_percentrequired) {
      acceptedhits++;

      // increase overlapcounter only for "real" reads,
      //  not for rails
      if(!SKIM3_readpool->getRead(actreadid).isRail()
	 && !SKIM3_readpool->getRead(rid2).isRail()){
	(*SKIM3_overlapcounter)[actreadid]+=1;
	(*SKIM3_overlapcounter)[rid2]+=1;
      }

      matchwithsorter_t tmp;
      tmp.otherid=rid2;
      tmp.eoffset=eoffsetmean;

      tmp.percent_in_overlap=perc;
      tmp.numhashes=numhashes;
      tmp.taken=false;
      tmpmatchwith.push_back(tmp);
      CEBUG("Pushing possible hit with offset: " << tmp.eoffset << endl);
    }

    if(actreadid>=1000000){
      CEBUG(rid2);
      CEBUG("\t" << SKIM3_readpool->getRead(rid2).getLenClippedSeq());
      CEBUG("\t" << hp1min);
      CEBUG("\t" << hp1max);
      //CEBUG("\t" << hp2min);
      //CEBUG("\t" << hp2max);
      CEBUG("\t" << eoffsetmin);
      CEBUG("\t" << eoffsetmax);
      CEBUG("\t" << maxoverlap);
      CEBUG("\t" << hashesoverlap);
      CEBUG("\t" << perc << '%');
      CEBUG('\n');
    }
    if(sI!=SKIM3_readhashmatches.end()) countid=sI->rid2;
  }
  CEBUG("Numhits " << actreadid << "\t" << possiblehits << "\t" << acceptedhits << "\n\n");

  SKIM3_possiblehits+=possiblehits;
  SKIM3_acceptedhits+=acceptedhits;
}
//#define CEBUG(bla)
//#define CEBUGF(bla)
