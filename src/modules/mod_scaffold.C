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



#include "modules/mod_scaffold.H"

#include "mira/parameters.H"
#include "mira/contig.H"
#include "mira/readpool_io.H"
#include "mira/maf_parse.H"

#include "version.H"


Scaffolder MiraScaffold::MS_scaffolder;

using namespace std;


vector<MIRAParameters> MS_Pv;


void MiraScaffold::cafmafload_callback(list<Contig> & clist, ReadPool & rp)
{
  for(auto & con : clist){
    MS_scaffolder.storeInfoFreshContig(con);
  }
  clist.clear();
}

int MiraScaffold::mainMiraScaffold(int argc, char ** argv)
{
  //CALLGRIND_STOP_INSTRUMENTATION;

  FUNCSTART("int MiraScaffold::mainMiraScaffold(int argc, char ** argv)");

  string filename("bla.maf");

  try{
    void (*usecrcallback)(list<Contig> &, ReadPool &) = cafmafload_callback;

    ReadPool rpool;
    list<Contig> clist;
    {
      MAFParse mafp(&rpool, &clist, &MS_Pv);
      mafp.setProgressIndicator(true);
      {
	vector<uint32> dummy;
	mafp.load(filename,
		  ReadGroupLib::SEQTYPE_SANGER,
		  1,
		  dummy,
		  false,
		  usecrcallback,
		  nullptr
	  );
      }
    }
    MS_scaffolder.dumpDebug();
  }
  catch(Notify n){
    // Need to close by hand as handleError() will perform a hard exit
    n.handleError(THISFUNC);
  }
  catch(Flow f){
    cerr << "Unexpected exception: Flow()\n";
  }
  catch(...){
    cout.flush();
    cerr.flush();
    cerr << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    abort();
  }

  FUNCEND();
  return 0;
}
