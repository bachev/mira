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

  SKIM_forwardhashok = new uint8[totalseqlen];
  SKIM_reversehashok = new uint8[totalseqlen];
  uint8 * fhashokp= SKIM_forwardhashok;
  uint8 * rhashokp= SKIM_reversehashok;

  uint32 totalhashes=0;

  for(int32 seqnr=0; seqnr<rp.getNOReads(); seqnr++){
    Read & actread= rp.getRead(seqnr);
    uint32 slen=actread.getLenClippedSeq();

    transformSeqToHash( actread.getClippedSeqAsChar(),
			slen,
			actread.getName(),
			fhashp,
			fhashokp,
			SKIM_forwardhashcount);

    fhashp+=slen;
    fhashokp+=slen;

//    rhashp=transformSeqToHash( actread.getClippedComplementSeqAsChar(),
//				 slen,
//				 actread.getName(),
//				 rhashp,
//				 SKIM_reversehashcount);
    totalhashes+=slen;
    SKIM_hashperseq[seqnr]=slen;

    //cout << "Hashes in read " << seqnr << ": " << SKIM_hashperseq[seqnr] << " (" << totalhashes << ")" << endl;
  }

  cout << "Computed " << totalhashes << " hashes." << endl;

  SKIM_forwardrahpsph= new rahp_t[totalhashes];
  SKIM_reverserahpsph= new rahp_t[totalhashes];

  {
    uint32 hashcount=0;
    for(int i=0; i<65536; i++) {
      SKIM_forwardreadptrfield[i]=hashcount;
      //cout << "Hashptrfield[" << i << "]: " << hashptrfield[i] << endl;
      hashcount+=SKIM_forwardhashcount[i];
    }
    
    fhashp=SKIM_forwardhashes;
    fhashokp=SKIM_forwardhashok;
    for(int32 seqnr=0; seqnr<rp.getNOReads(); seqnr++){
      for(int i=0; i<SKIM_hashperseq[seqnr]; i++, fhashp++, fhashokp++){
	if(*fhashokp) {
	  SKIM_forwardrahpsph[SKIM_forwardreadptrfield[*fhashp]].read=seqnr;
	  SKIM_forwardrahpsph[SKIM_forwardreadptrfield[*fhashp]].hashpos=i;
	  SKIM_forwardreadptrfield[*fhashp]++;
	}
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

  uint32 numseq=SKIM_hashperseq.size();
  vector<uint32> count;

  const uint32 maxoffsets=100;

  int16 * offsetarrays = new int16[maxoffsets*numseq];
  vector<uint32> offsetcountused;

  uint16 * hashp=SKIM_forwardhashes;
  uint8  * hashokp=SKIM_forwardhashok;

  uint16 * histogram = new uint16[40000];
  memset(histogram,0,80000);
  uint16 * midhist=histogram+20000;

  // pruefen ob histogramm leer
  {
    uint16 * ptr=histogram;
    for(uint32 i=0; i<40000; i++,ptr++){
      if(*ptr!=0){
	cout << "Arg! Start-Histogramm nicht leer :(\n";
	exit(100);
      }
    }
  }

  uint32 totalpossiblehits=0;
  uint32 totalgoodhits=0;
  uint32 totalexcellenthits=0;

  for(uint32 readnr=0; readnr<numseq; readnr++){
    //for(uint32 readnr=0; readnr<1; readnr++){
    count.resize(0);
    count.resize(numseq,0);
    offsetcountused.resize(0);
    offsetcountused.resize(numseq,0);

    for(int16 readpos=0; readpos< SKIM_hashperseq[readnr];readpos++, hashp++, hashokp++){
      if(*hashokp==0) continue;
      
      rahp_t * rahpsph=SKIM_forwardrahpsph+SKIM_forwardreadptrfield[*hashp];
      for(uint32 j=0, loopto= SKIM_forwardhashcount[*hashp]; j<loopto; j++, rahpsph++) {
	if(rahpsph->read<=readnr) {
	  // hmmm, this could save time on BIG projects (>10k reads)
	  // yep, it does. around 10 to 20% for a 45k read project
	  // but arrays need to be rebuild afterwards if they are to be used again
	  SKIM_forwardreadptrfield[*hashp]++;
	  SKIM_forwardhashcount[*hashp]--;

	  continue;  
	}
	count[rahpsph->read]++;
	if(offsetcountused[rahpsph->read] == maxoffsets) continue;
	offsetarrays[((rahpsph->read)*maxoffsets)+offsetcountused[rahpsph->read]]=readpos-(rahpsph->hashpos);
	offsetcountused[rahpsph->read]++;
      }
    }

    uint32 possiblehit=0;
    uint32 goodhit=0;
    uint32 excellenthit=0;
    for(uint32 seq=0; seq<count.size(); seq++){
      //cout << "Hit with " << i << ": " << count[i] << endl;
      if(count[seq]>=4) {
	uint16 maxhist=0;
	for(uint32 offset=0; offset<offsetcountused[seq]; offset++){
	  midhist[offsetarrays[seq*maxoffsets+offset]]++;
	  if(midhist[offsetarrays[seq*maxoffsets+offset]] > maxhist){
	    maxhist=midhist[offsetarrays[seq*maxoffsets+offset]];
	  }
	}
	possiblehit++;
	if(maxhist>=4 && maxhist <10) {
	  for(uint32 offset=0; offset<offsetcountused[seq]; offset++){
	    int32 localoffset=offsetarrays[seq*maxoffsets+offset];
	    if(midhist[localoffset]>0){
	      int32 localto=localoffset;
	      while(midhist[localto]>0) localto++;
	      int32 localfrom=localoffset;
	      while(midhist[localfrom-1]>0) localfrom--;
	      uint32 maxsum=0;
	      if(localto-localfrom <= 3) {
		for(;localfrom<localto;localfrom++){
		  maxsum+=midhist[localfrom];
		  midhist[localfrom]=0;
		}
	      } else {
		int32 localsub=localfrom;
		int32 sum=midhist[localfrom++];
		sum+=midhist[localfrom++];
		sum+=midhist[localfrom++];
		maxsum=sum;
		for(;localfrom<localto;localfrom++,localsub++) {
		  sum+=midhist[localfrom];
		  sum-=midhist[localsub];
		  midhist[localsub]=0;
		  if(sum>maxsum) {
		    maxsum=sum;
		  }
		}
		midhist[localfrom-1]=0;
		midhist[localfrom-2]=0;
		midhist[localfrom-3]=0;
	      }
	      if(maxsum>=10){
		excellenthit++;
	      } else {
		goodhit++;
	      }
	    }
	  }
	}

	if(maxhist>=10) excellenthit++;


	// histogramm loeschen
	for(uint32 offset=0; offset<offsetcountused[seq]; offset++){
	  midhist[offsetarrays[seq*maxoffsets+offset]]=0;
	}

	// pruefen ob histogramm leer
//	  {
//	    uint16 * ptr=histogram;
//	    for(uint32 i=0; i<4000; i++,ptr++){
//	      if(*ptr!=0){
//		cout << "Arg! End-Histogramm nicht leer :(\n";
//		exit(100);
//	      }
//	    }
//	  }

      }
    }

    cout << "Hits with read " << readnr << ":   excellent(" << excellenthit << ")\tgood(" << goodhit << ")\tpossible(" << possiblehit << ")" <<endl;
    totalpossiblehits+=possiblehit;
    totalgoodhits+=goodhit;
    totalexcellenthits+=excellenthit;
  }

  cout << "\nTotal hits:   excellent(" << totalexcellenthits << ")\tgood(" << totalgoodhits << ")\tpossible(" << totalpossiblehits << ")" <<endl;

  delete [] histogram;
  delete [] offsetarrays;
}


//	    uint16 *from =midhist-1000;
//	    uint16 *addp =from;
//	    uint16 *subp =from;
//	    
//	    uint32 sum= *(addp++);
//	    sum+=*(addp++);
//	    sum+=*(addp++);
//	    uint32 maxsum=sum;
//	    for(int16 i=0; i<2000; i++, addp++, subp++){
//	      sum+=*(addp);
//	      sum-=*(subp);
//	      if(sum>maxsum){
//		maxsum=sum;
//	      }
//	    }
//	    if (maxsum>5){
//	      excellenthit++;
//	    }else{
//	      goodhit++;
//	    }




void Skim::transformSeqToHash (const char * seq, uint32 slen, const char * readname, uint16 * hashp, uint8 * hashokp, uint32 * hashcount)
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
      // just break: behave like an 'a'   (wrong in 75%, but right in 25%
      //  and doesn't break the hash making!)
      // feel free to make addhash+1, +2 or +3 for c, g ot t instead
      
      // or

      // uncomment this to break hash making (which is actually better,
      //  in case of multiple bases with 'n' or '*')
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
      //cout << "Hash: " << hex << acthash << dec << endl;
      
      hashcount[acthash]++;

      *hashp=acthash;
      hashp++;
      *hashokp=1;
      hashokp++;
    } else {
      //cout << "Missed hash" << endl;
      if(i>=8) {
	*hashp=0;
	hashp++;
	*hashokp=0;
	hashokp++;
      }

    }
  }
  
  for(uint32 i=0; i<7; i++, hashp++, hashokp++) {
    *hashp=0;
    *hashokp=0;
  }

}



