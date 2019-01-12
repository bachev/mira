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
 * Example read_101
 *
 *
 * Topics covered:
 * - always enclose usage of MIRA classes with a try/catch block
 *   (at least on the main program level)
 * - setting name and sequence
 * - setting quality and sequencing vector clips
 * - defining output types
 * - setting quality values
 * - sequence and quality manipulation
 * - setting additional information like template name etc.
 * - adding tags to sequences
 *
 *****************************************************************/

// include C++ iostreams
#include <iostream>

// include MIRA read class
#include "mira/read.H"


int main(int argc, char ** argv)
{
  cout << "Example program read_101: basic Read functions.\n";

  // ALWAYS use a try/catch block when using MIRA library classes
  //  the way it is shown here. On nasty problems, the classes will
  //  throw exceptions to abort.

  try{
    // create a read
    Read r;

    // we want our read at the minimum to have a name and a sequence
    //  for the read
    r.setName("U45pbk14A");
    r.setSequenceFromString("aggtcatgacttctcatgactgcatcatctgctactgatcgtgcatg");

    // have a look at how it looks like ... as FASTA
    // tell the read class to output as FASTA
    Read::setCoutType(Read::AS_FASTA);
    cout << r;

    /**************** Clips and output types ******************/

    // play a bit around: set left and right (quality) clips
    r.setClipoffs(3,r.getLenSeq()-6, false);
    // look at it again (still FASTA)
    cout << r;
    // and now only the clipped part
    Read::setCoutType(Read::AS_CLIPPEDFASTA);
    cout << r;

    // set left sequencing vector clip
    r.setLSClipoff(7);
    // look at it again (still CLIPPEDFASTA)
    cout << r;

    // have a look at it in gap4 EXP file format
    Read::setCoutType(Read::AS_GAP4DA);
    cout << r;

    /******************* Quality values ***********************/
    // the default quality value for each base in a read is 10
    // set it to something else, like 25
    r.setQualities(25);
    cout << r;

    // if you have a vector qith qualities (perhaps loaded from
    //  a file), well, then this also works.
    // here we use a vector with all values at 1
    {
      vector<base_quality_t> quals(r.getLenSeq(),1);
      r.setQualities(quals);
    }
    cout << r;

    /************ Sequence and quality manipulation **********/
    // change, insert and delete a base in the sequence
    // change base at pos 12 to A and give quality 99
    r.changeBaseInSequence('A',99,12);

    // delete base at position 14
    r.deleteBaseFromSequence(14);

    // insert a base 'G' at position 16 with quality 98
    r.insertBaseInSequence('G',98,16, false);

    cout << r;

    /************ Set a few additional information ***********/
    r.setTemplate("ACklp56tft");

    cout << r;

    /***************** Setting tags **************************/
    r.addTag(8,15,"MIRA","An example tag.");

    // oh well, just for fun change the output type to CAF
    // (note that e.g. the cloning vector name is not contained 
    //  in this format)
    Read::setCoutType(Read::AS_CAF);
    cout << r;
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

