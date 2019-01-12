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


//TODO
// an Kommentare Hinweis auf Cut setzen (originalfile, frombase, tobase

#include "io/scf.H"
#include "io/phd.H"


void usage()
{
  cout << "mergeSCFandPHD V1.01    Bastien Chevreux\nUsage:\n";
  cout << "  mergeSCFandPHD scf_infile phd_infile scf_outfile\n";
}

int main(int argc, char **argv)
{
  if(argc!=4) {
    usage();
    return 1;
  }

  try{
    SCF thescf;
    PHD thephd;

    thescf.load(argv[1]);

    thephd.load(argv[2]);

    // cludge, until all classes have been changed to STL and string
    const string & seq=thephd.getSequence();
    uint8 * probs= new uint8[seq.length()];
    uint32 * peakpos= new uint32[seq.length()];
    vector<uint8>::const_iterator probI=thephd.getQualities().begin();
    vector<uint32>::const_iterator peakI=thephd.getPeakIndex().begin();
    for(uint32 i=0; i<seq.length(); i++, probI++, peakI++){
      probs[i]=*probI;
      peakpos[i]=*peakI;
    }

    thescf.setBases(seq.length(), seq.c_str(), probs, peakpos);
    const SCF_Comments * comments =thescf.getComments();
    if(comments != nullptr) {
      uint32 comlen=strlen(comments);
      if(comments[comlen-1] != '\n') thescf.addComment("\n");
    }
    thescf.addComment("TPSW=mergeSCFandPHD V1.01\n");

    thescf.save(argv[3]);
  }
  catch(Notify n){
    n.handleError("main");
    exit(1);
  }


  return 0;
}
