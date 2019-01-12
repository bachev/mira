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

#include <list>

#include "errorhandling/errorhandling.H"

#include "mira/gff_save.H"

#include "mira/contig.H"

using std::cout;
using std::cerr;
using std::endl;


// Plain vanilla constructor
GFFSave::GFFSave()
{
  FUNCSTART("GFFSave::GFFSave()");

  zeroVars();

  FUNCEND();
}

void GFFSave::zeroVars()
{
  FUNCSTART("void GFFSave::zeroVars()");

  GFFS_gff3_names.clear();
  GFFS_gff3_paddedcons.clear();
  GFFS_gff3_strainnames.clear();

  GFFS_gff3_strainnames.push_back("AllStrains");
  GFFS_fouts.clear();
  GFFS_fouts.push_back(new std::ofstream);

  FUNCEND();
}


GFFSave::~GFFSave()
{
  FUNCSTART("GFFSave::~GFFSave()");

  discard();

  FUNCEND();
}


void GFFSave::discard()
{
  FUNCSTART("GFFSave::discard()");

  if(GFFS_foutp.is_open()) close();
  zeroVars();

  FUNCEND();
}

void GFFSave::open(const char * filename)
{
  FUNCSTART("void GFFSave::open(const char * filename)");

  BUGIFTHROW(GFFS_foutp.is_open(),"GFFSave object already has a file open? ('"<<GFFS_foutname<<"')\n");

  GFFS_foutname=filename;
  std::string basename = GFFS_foutname;

  GFFS_foutp.open((basename+"_AllStrains_padded.gff3"), std::ios::out);

  GFFS_foutp << "##gff-version 3\n#written by MIRA assembler\n";

  FUNCEND();
}

void GFFSave::close()
{
  FUNCSTART("GFFSave::close()");

  BUGIFTHROW(!is_open(),"GFFSave object not opened ... why the close()?");

  BUGIFTHROW(GFFS_gff3_names.size() != GFFS_gff3_paddedcons.size(),
	     "Inconsisteny in cached data: " << GFFS_gff3_names.size() << " " << GFFS_gff3_paddedcons.size());

  if(!GFFS_gff3_names.empty()){
    GFFS_foutp << "##FASTA\n";
    for(uint32 gffsi=1; gffsi<GFFS_fouts.size(); ++gffsi){
      *(GFFS_fouts[gffsi]) << "##FASTA\n";
    }
    for(uint32 straini=0; straini<GFFS_gff3_strainnames.size(); ++straini){
      auto nI=GFFS_gff3_names.begin();
      auto mpcI=GFFS_gff3_paddedcons.begin();
      for(; nI!=GFFS_gff3_names.end(); ++nI, ++mpcI){
	auto & vs=*mpcI;
	if(straini<vs.size()){
	  bool dooutput=false;
	  if(!vs[straini].empty()){
	    auto sI=vs[straini].begin();
	    for(; sI != vs[straini].end(); ++sI){
	      if(*sI!='X'){
		dooutput=true;
		break;
	      }
	    }
	  }
	  if(dooutput){
	    GFFS_foutp << ">" << *nI << ' ' << GFFS_gff3_strainnames[straini] << '\n';
	    uint32 charcounter=0;
	    auto sI=vs[straini].begin();
	    for(; sI != vs[straini].end(); ++sI, ++charcounter){
	      if(charcounter==1200) {
		charcounter=0;
		GFFS_foutp << '\n';
	      }
	      GFFS_foutp << *sI;
	    }
	    GFFS_foutp << '\n';

	    if(straini>0 && straini < GFFS_fouts.size()){
	      *(GFFS_fouts[straini]) << ">" << *nI << ' ' << GFFS_gff3_strainnames[straini] << '\n';
	      sI=vs[straini].begin();
	      for(; sI != vs[straini].end(); ++sI){
		if(charcounter==1200) {
		  charcounter=0;
		  *(GFFS_fouts[straini]) << '\n';
		}
		if(*sI != '*'){
		  *(GFFS_fouts[straini]) << *sI;
		  ++charcounter;
		}
	      }
	      *(GFFS_fouts[straini]) << '\n';
	    }
	  }
	}
      }
    }
  }

  GFFS_foutp.close();
  for(uint32 gffsi=1; gffsi<GFFS_fouts.size(); ++gffsi){
    GFFS_fouts[gffsi]->close();
    delete GFFS_fouts[gffsi];
  }
  GFFS_fouts.clear();

  zeroVars();

  FUNCEND();
}


void GFFSave::acquireContig(Contig & con, const ReadPool & rp)
{
  FUNCSTART("void GFFSave::acquireContig(Contig & con, const ReadPool & rp)");

  BUGIFTHROW(!is_open(),"GFFSave object not opened?");

  // insert strain names newly seen in this readpool into GFFS_gff3_strainnames if needed
  for(uint32 rpsi=0; rpsi<ReadGroupLib::getNumOfStrains(); ++rpsi){
    if(con.getNumReadsPerStrain(rpsi)==0) continue;
    bool found=false;
    for(uint32 gffsi=0; gffsi < GFFS_gff3_strainnames.size(); ++gffsi){
      if(GFFS_gff3_strainnames[gffsi] == ReadGroupLib::getStrainOfStrainID(rpsi)){
	found=true;
	break;
      }
    }
    if(!found){
      GFFS_gff3_strainnames.push_back(ReadGroupLib::getStrainOfStrainID(rpsi));
      // TODO: for unpadded, also open file here
      GFFS_fouts.push_back(new std::ofstream);
      GFFS_fouts.back()->open((GFFS_foutname+"_"+ReadGroupLib::getStrainOfStrainID(rpsi)+"_unpadded.gff3"), std::ios::out);
      *(GFFS_fouts.back()) << "##gff-version 3\n#written by MIRA assembler\n";
    }
  }

//  {
//    cout << "rp strains: " << ReadGroupLib::getNumOfStrains() << endl;
//    for(uint32 rpsi=0; rpsi<ReadGroupLib::getNumOfStrains(); ++rpsi){
//     cout << rpsi << "\t" << ReadGroupLib::getStrainOfStrainID(rpsi) << endl;
//    }
//
//    cout << "\nGFFS_gff3_strainnames:" << endl;
//    for(uint32 i=0; i<GFFS_gff3_strainnames.size(); ++i){
//      cout << i << "\t" << GFFS_gff3_strainnames[i] << endl;
//    }
//    cout << endl;
//  }

  {
    std::vector<std::string> tmp;
    GFFS_gff3_paddedcons.push_back(tmp);
    GFFS_gff3_paddedcons.back().resize(GFFS_gff3_strainnames.size());
  }
  GFFS_gff3_names.push_back(con.getContigName());

  // now push back the strain sequences
  std::vector<base_quality_t> dummyqual;
  for(uint32 gffsi=0; gffsi<GFFS_gff3_strainnames.size(); ++gffsi){
    int32 rpsi=-1;
    if(gffsi>0){
      for(rpsi=0; rpsi< static_cast<int32>(ReadGroupLib::getNumOfStrains()); ++rpsi){
	if(GFFS_gff3_strainnames[gffsi] == ReadGroupLib::getStrainOfStrainID(rpsi)) break;
      }
    }
    if(rpsi < static_cast<int32>(ReadGroupLib::getNumOfStrains())){
      //cout << "Trying gffsi " << gffsi << " with rpsi " << rpsi << endl;
      con.newConsensusGet(GFFS_gff3_paddedcons.back()[gffsi], dummyqual, rpsi);
    }else{
    }
  }

  // now print out all the tags of the backbone sequence(s)
  // and consensus tags
  con.dumpTagsAsGFF3(GFFS_foutp);

  // TODO: for unpadded seqs
  // start at 1, as strain "0" is the AllStrains "strain"
  for(uint32 gffsi=1; gffsi<GFFS_gff3_paddedcons.back().size(); ++gffsi){
    con.dumpTagsAsGFF3(*(GFFS_fouts[gffsi]), GFFS_gff3_paddedcons.back()[gffsi]);
  }


  FUNCEND();
}



void GFFSave::acquireRead(Read & read)
{
  FUNCSTART("void GFFSave::acquireRead(Read & read)");

  BUGIFTHROW(!is_open(),"GFFSave object not opened?");

  GFFS_gff3_names.push_back(read.getName());

  {
    std::vector<std::string> tmp;
    GFFS_gff3_paddedcons.push_back(tmp);
  }
  GFFS_gff3_paddedcons.back().resize(1);
  read.getSeqAsString(GFFS_gff3_paddedcons.back().back());

  read.dumpTagsAsGFF3(GFFS_foutp);

  FUNCEND();
}



//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//GFFSave::GFFSave(GFFSave const &other)
//{
//  FUNCSTART("GFFSave::GFFSave(GFFSave const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//GFFSave const & GFFSave::operator=(GFFSave const & other)
//{
//  FUNCSTART("GFFSave const & GFFSave::operator=(GFFSave const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, GFFSave const &???)
//{
//  FUNCSTART("friend ostream & GFFSave::operator<<(ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}
