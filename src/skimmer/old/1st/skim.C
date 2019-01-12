/*
 * Written by Bastien Chevreux (BaCh)
 *
 *
 */
// 	$Id$	

#ifndef lint
static char vcid[] = "$Id$";
#endif /* lint */

#include "skim.H"


void Skim::foolCompiler()
{
#include "stdinc/foolcompiler.C"
}

// Plain vanilla constructor
Skim::Skim(ReadPool & rp)
{
  FUNCSTART("Skim::Skim(ReadPool & rp)");
  //  ERROR("Not implemented yet.");

  prepareSkim(rp);

  FUNCEND();
}

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

ostream & operator<<(ostream &ostr, Skim const &theskim)
{
  FUNCSTART("friend ostream & Skim::operator<<(ostream &ostr, const  &theskim)");
  //  ERROR("Not implemented yet.");

  FUNCEND();
  return ostr;
}

void Skim::discard()
{
  FUNCSTART("Skim::discard()");

  if(SKIM_valid==0){
    //    ERROR("Not implemented yet.");
  }

  FUNCEND();
}


void Skim::prepareSkim(ReadPool & rp)
{
  FUNCSTART("void Skim::prepareSkim(ReadPool & rp)");

  uint32 totalseqlen=0;
  for(uint32 i=0; i<rp.getNOReads(); i++) {
    totalseqlen+=rp.getRead(i).getLenClippedSeq();
  }

  cout << rp.getNOReads() << " sequences to skim totalling " << totalseqlen << " bases." << endl;
  
  for(int32 i=0; i< 65536; i++){
    SKIM_forwardhashcount[i]=0;
    SKIM_reversehashcount[i]=0;
  }

  SKIM_hashperseq.resize(rp.getNOReads(),0);

  SKIM_forwardhashes = new uint16[totalseqlen];
  SKIM_reversehashes = new uint16[totalseqlen];
  uint16 * fhashp= SKIM_forwardhashes;
  uint16 * rhashp= SKIM_reversehashes;

  uint32 totalhashes=0;
  uint32 lasthashes=0;

  for(int32 seqnr=0; seqnr<rp.getNOReads(); seqnr++){
    Read & actread= rp.getRead(seqnr);
    uint32 slen=actread.getLenClippedSeq();
    if (slen < 24) continue;

    fhashp=transformSeqToHash( actread.getClippedSeqAsChar(),
			       slen,
			       actread.getName(),
			       fhashp,
			       SKIM_forwardhashcount);
//    rhashp=transformSeqToHash( actread.getClippedComplementSeqAsChar(),
//				 slen,
//				 actread.getName(),
//				 rhashp,
//				 SKIM_reversehashcount);
    totalhashes=(fhashp-SKIM_forwardhashes);
    SKIM_hashperseq[seqnr]=totalhashes-lasthashes;
    lasthashes=totalhashes;

    cout << "Hashes in read " << seqnr << ": " << SKIM_hashperseq[seqnr] << " (" << totalhashes << ")" << endl;
  }

  cout << "Computed " << totalhashes << " hashes." << endl;

  SKIM_forwardrsph= new uint32[totalhashes];
  SKIM_reversersph= new uint32[totalhashes];

  {
    uint32 hashcount=0;
    for(int i=0; i<65536; i++) {
      SKIM_forwardreadptrfield[i]=hashcount;
      //cout << "Hashptrfield[" << i << "]: " << hashptrfield[i] << endl;
      hashcount+=SKIM_forwardhashcount[i];
    }
    
    fhashp=SKIM_forwardhashes;
    for(int32 seqnr=0; seqnr<rp.getNOReads(); seqnr++){
      for(int i=0; i<SKIM_hashperseq[seqnr]; i++, fhashp++){
	SKIM_forwardrsph[SKIM_forwardreadptrfield[*fhashp]]=seqnr;
	SKIM_forwardreadptrfield[*fhashp]++;
      }
    }

    hashcount=0;
    for(int i=0; i<65536; i++) {
      SKIM_forwardreadptrfield[i]=hashcount;
      //cout << "Hashptrfield[" << i << "]: " << hashptrfield[i] << endl;
      hashcount+=SKIM_forwardhashcount[i];
    }
  }

//  {
//    uint32 hashcount=0;
//    for(int i=0; i<65536; i++) {
//      SKIM_reversereadptrfield[i]=hashcount;
//      //cout << "Hashptrfield[" << i << "]: " << hashptrfield[i] << endl;
//      hashcount+=SKIM_reversehashcount[i];
//    }
//    
//    fhashp=SKIM_reversehashes;
//    for(int32 seqnr=0; seqnr<rp.getNOReads(); seqnr++){
//      for(int i=0; i<SKIM_hashperseq[seqnr]; i++, fhashp++){
//	 SKIM_reversersph[SKIM_reversereadptrfield[*fhashp]]=seqnr;
//	 SKIM_reversereadptrfield[*fhashp]++;
//      }
//    }
//
//    hashcount=0;
//    for(int i=0; i<65536; i++) {
//      SKIM_reversereadptrfield[i]=hashcount;
//      //cout << "Hashptrfield[" << i << "]: " << hashptrfield[i] << endl;
//      hashcount+=SKIM_reversehashcount[i];
//    }
//  }


  FUNCEND();
}


void Skim::go()
{
  cout << "Have a go ...\n";

  vector<uint32> count;
  count.resize(SKIM_hashperseq.size(),0);

  uint16 * hashp=SKIM_forwardhashes;
  for(uint32 i=0; i< SKIM_hashperseq[0];i++, hashp++){
    
    //    cout << " SKIM_forwardreadptrfield[*hashp] " << SKIM_forwardreadptrfield[*hashp];
    uint32 * rsph=SKIM_forwardrsph+SKIM_forwardreadptrfield[*hashp];
    for(uint32 j=0; j<SKIM_forwardhashcount[*hashp]; j++, rsph++) {
      if(*rsph==0) continue;
      cout << "Acthash " << hex << *hashp << dec;
      cout << " ding " << *rsph;
      count[*rsph]++;
    }
    if (SKIM_forwardhashcount[*hashp]>1) cout << endl;
  }

  for(uint32 i=0; i<count.size(); i++){
    cout << "Hit with " << i << ": " << count[i] << endl;
  }
}


uint16 * Skim::transformSeqToHash (const char * seq, uint32 slen, const char * readname, uint16 * hashp, uint32 * hashcount)
{

  enum { SEARCHNEW, SEARCHING, ACGTCONTINUE};
  uint32 state = SEARCHNEW;
  
  uint16 acthash =0;
  uint32 baseok=0;

  for(int32 i=0; i<slen; i++, seq++){
    acthash<<=2;
    baseok++;
    
    switch (*seq) {
    case 'a' :
    case 'A' : break;
    case 'c' :
    case 'C' : {
      acthash+=1;
      break;
    }
    case 'g' :
    case 'G' : {
      acthash+=2;
      break;
    }
    case 't' :
    case 'T' : {
      acthash+=3;
      break;
    }
    case 'n' :
    case 'N' :
    case '*' : {
      acthash=0;
      baseok=0;
	break;
    }
    default : {
      cout << "Unknown base '" << *seq << "' (ASCII " << (uint16) *seq << ") at position " << i << " in _CLIPPED_ sequence " << *readname << endl;
      exit(100);
    }
    }
    
    //      cout << *seq << " ";
    if(baseok>=8) {
      cout << "Hash: " << hex << acthash << dec << endl;
      
      *hashp=acthash;
      hashp++;
      
      hashcount[acthash]++;
    } else {
      //cout << "Missed hash" << endl;
    }
  }
  
  return hashp;
}



