/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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

#include "assembly.H"

#define CEBUG(bla)

using std::cout;
using std::cerr;
using std::endl;

/*************************************************************************
 * works a bit like assemble()
 * after assembly, all reads that are not singlets
 *  OR having a SRMB or WRMB tag are thought of as 'good'
 *
 * Then, all reads marked as good are saved into CAF files by strain,
 *  i.e., each strain gets its own CAF, with repeat tags (SRMx, WRMx)
 *  and allelic variances tags (PALV) and eventually SNP tags (PSNP etc.)
 *
 * Function returns a vector of strings where elements are paired by two:
 *    filename1    strainname1
 *    filename2    strainname2 etc...
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

#if 0
vector<string> Assembly::assembleESTs()
{
  FUNCSTART("Assembly::assembleESTs()");

  //AS_miraparams[0].setPathfinderuseGenomicAlgorithms(true);

  BUGIFTHROW(true,"need redo for PlacedContigReads");

  vector<string> returns;

  assemble();
  saveResults();

  vector<bool> goodreads;
  goodreads.resize(AS_readpool.size(),false);
  {
    CEBUG(AS_contigs.size() << " contigs in list.\n");
    for(const auto & cle : AS_contigs){
      if(cle.getNumReadsInContig()<2) {
	continue;
      }
      for(auto pcrI=cle-getContigReads().begin(); pcrI != cle.getContigReads().end(); ++pcrI){
	if(pcrI.getORPID() >=0 ) goodreads[pcrI.getORPID()]=true;
      }
    }

    for(uint32 ri=0; ri < AS_readpool.size(); ++ri) {
      try {
	if(AS_readpool[ri].hasTag(Read::REA_tagentry_idSROr)) goodreads[ri]=true;
	if(AS_readpool[ri].hasTag(Read::REA_tagentry_idSIOr)) goodreads[ri]=true;
	if(AS_readpool[ri].hasTag(Read::REA_tagentry_idSAOr)) goodreads[ri]=true;
	if(AS_readpool[ri].hasTag(Read::REA_tagentry_idSRMr)) goodreads[ri]=true;
	if(AS_readpool[ri].hasTag(Read::REA_tagentry_idWRMr)) goodreads[ri]=true;
      }
      catch(...){
	CEBUG("Caught: " << AS_readpool[ri].getName() << '\n');
      }
      CEBUG("GR[" << AS_readpool[ri].getName() << "]: " << goodreads[ri] << endl);
    }
  }

  vector<int32> goodstrainids;
  vector<string> goodstrainnames;
  for(uint32 ri=0; ri < AS_readpool.size(); ++ri){
    if(find(goodstrainids.begin(),
	    goodstrainids.end(),
	    AS_readpool[ri].getStrainID()) == goodstrainids.end()) {
      goodstrainids.push_back(AS_readpool[ri].getStrainID());
      if(AS_readpool[ri].getStrainName().size()>0){
	goodstrainnames.push_back(AS_readpool[ri].getStrainName());
      }else{
	goodstrainnames.push_back("default");
      }
    }
  }

  {
    string basename=AS_miraparams[0].getAssemblyParams().as_projectname_out;
    Read::setCoutType(Read::AS_CAF);
    for(uint32 gsi=0; gsi < goodstrainids.size(); ++gsi){
      cout << "Writing CAF for reads for strain " << goodstrainnames[gsi] << endl;
      string filename=basename+"_snpsinSTRAIN_"+goodstrainnames[gsi]+".caf";
      returns.push_back(filename);
      returns.push_back(goodstrainnames[gsi]);
      ofstream cafout(filename.c_str(), ios::out | ios::trunc);
      for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
	if(goodreads[ri]==true && AS_readpool[ri].getStrainID()==goodstrainids[gsi]){
	  AS_readpool[ri].removeGapsFromRead();
	  cafout << AS_readpool[ri];
	}
      }
    }

    // also put reads with no potential SNPs detected into an own "STRAIN"
    {
      string filename=basename+"_nosnps_remain.caf";
      returns.push_back(filename);
      returns.push_back("remain");
      ofstream cafout(filename.c_str(), ios::out | ios::trunc);
      for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
	if(goodreads[ri]==false) {
	  AS_readpool[ri].removeGapsFromRead();
	  cafout << AS_readpool[ri];
	}
      }
    }

  }

  FUNCEND();

  return returns;
}
//#define CEBUG(bla)   {cout << bla; cout.flush(); }

#endif
