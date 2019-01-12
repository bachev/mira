/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2007 and later by Bastien Chevreux
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


// 	$Id$	


/*****************************************************************
 *
 * Example readpool_101
 *
 *
 * Topics covered:
 * - always enclose usage of MIRA classes with a try/catch block
 *   (at least on the main program level)
 * - add new empty reads to a ReadPool
 * - accessing reads in readpool by index or by name
 * - (from class Read: setting some basic values)
 * - (from class Read: output single reads)
 * - dumping the readpool in different formats to a stream
 *
 *****************************************************************/

// include C++ iostreams
#include <iostream>

// include MIRA readpool class
#include "mira/readpool.H"


int main(int argc, char ** argv)
{
  cout << "Example program readpool_101: basic ReadPool functions.\n";

  // ALWAYS use a try/catch block when using MIRA library classes
  //  the way it is shown here. On nasty problems, the classes will
  //  throw exceptions to abort.

  try{
    // for creating a readpool, we need a vector of MIRAParameters
    vector<MIRAParameters> Pv;
    MIRAParameters::setupStdMIRAParameters(Pv);

    // create a readpool
    ReadPool thepool(&Pv);

    // create three new empty reads in the readpool
    thepool.addNewEmptyRead();
    thepool.addNewEmptyRead();
    thepool.addNewEmptyRead();

    // initialise that last read in the pool with some values
    thepool.back().setName("MyRead3");
    thepool.back().setSequenceFromString("ttggtacgtacgttggtac");
    thepool.back().setRQClipoff(thepool.back().getLenSeq()-9);


    // initialise the first read in the pool with some values
    thepool[0].setName("MyRead1");
    thepool[0].setSequenceFromString("gacttgactagctgactgactgacgtacgtac");
    thepool[0].setLSClipoff(7);

    // The second new empty read in the readpool
    // won't be initialised with sequence -> read will be invalid
    thepool[1].setName("MyRead2");


    // dump some pool info to cout stream
    thepool.dumpPoolInfo(cout);

    // fetch the read with index 0 (== the first)
    // and dump to cout as EXP (GAP4DA) format
    Read & aread=thepool[0];
    Read::setCoutType(Read::AS_GAP4DA);
    cout << aread;


    // fetch the read with name "MyRead2" (== the second)
    // and dump to cout as text format (for debugging or so)
    aread=thepool.getRead("MyRead2");
    Read::setCoutType(Read::AS_TEXT);
    cout << aread;


    // dump all the reads (valid and invalid) as CAF to cout
    thepool.dumpAs(cout,Read::AS_CAF,true);

    // dump the valid reads as CAF to cout
    thepool.dumpAs(cout,Read::AS_CAF,false);

    // dump the valid reads with sequencing vector masked to cout
    //  as FASTA
    thepool.dumpAs(cout,Read::AS_SEQVECMASKEDFASTA,false);
    //  as FASTA quality
    thepool.dumpAs(cout,Read::AS_SEQVECMASKEDFASTAQUAL,false);

    // dump the valid reads, only the "good" part without clips
    //  as FASTA
    thepool.dumpAs(cout,Read::AS_CLIPPEDFASTA,false);
    //  as FASTA quality
    thepool.dumpAs(cout,Read::AS_CLIPPEDFASTAQUAL,false);
  }
  // Catch exceptions as thrown by MIRA classes
  catch(Notify n){
    n.handleError("main");
    return 1;
  }
  catch(Flow f){
    cerr << "INTERNAL ERROR: Unexpected exception: Flow()\n";
    return 1;
  }
  // you can define own catches here
  catch(...){
    cerr << "Unknown exception caught, aborting the process.\n\n";
    return 1;
  }


  return 0;
}

