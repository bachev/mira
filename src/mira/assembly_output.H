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


#ifndef _mira_assout_h_
#define _mira_assout_h_


#include "mira/assembly_info.H"
#include "mira/contig.H"
#include "mira/readpool.H"


namespace assout {

  struct shortgbfinfo_t {
    std::string identifier;
    std::string locustag;

    std::string gene;
    std::string function;
    std::string ecnumber;
    std::string product;
    std::string note;

    friend std::ostream & operator<<(std::ostream &ostr, const shortgbfinfo_t & sgi){
      ostr << "identifier: " << sgi.identifier
	   << "\nlocustag  : " << sgi.locustag
	   << "\ngene      : " << sgi.gene
	   << "\nfunction  : " << sgi.function
	   << "\necnumber  : " << sgi.ecnumber
	   << "\nproduct   : " << sgi.product
	   << "\nnote      : " << sgi.note
	   << std::endl;
      return ostr;
    }
  };

  struct snpenv_t {
    std::list<Contig>::iterator cI;
    std::string snpjumpname;

    std::list<gbfsummary_t> gbfslist;
    std::list<shortgbfinfo_t> sgilist;

    uint32 from;
    uint32 to;
    uint32 numsnps;

    friend std::ostream & operator<<(std::ostream &ostr, const snpenv_t & se){
      ostr << "Contig: " << se.cI->getContigName()
	   << "\njumpname: " << se.snpjumpname
	   << "\nfrom      : " << se.from
	   << "\nto        : " << se.to
	   << "\nnumsnps   : " << se.numsnps
	   << std::endl;

      ostr << "gbfslist: ";
      if(se.gbfslist.empty()){
	ostr << "empty\n";
      }else{
	ostr << "\n";
	for(auto & gbfs : se.gbfslist){
	  ostr << gbfs;
	}
      }
      ostr << "sgilist: ";
      if(se.sgilist.empty()){
	ostr << "empty\n";
      }else{
	ostr << "\n";
	std::list<shortgbfinfo_t>::const_iterator sgI=se.sgilist.begin();
	for(auto & sg : se.sgilist){
	  ostr << sg;
	}
      }

      return ostr;
    }

  };


  struct proteinchangesummary_t {
    bool isaffected;

    bool firstcodonisstart;

    bool changedstart;
    bool destroyedstart;
    bool changedstop;
    bool destroyedstop;
    bool prematurestop;

    uint32 mutinintergenic;

    uint32 insertioninlocus;
    uint32 deletioninlocus;
    uint32 silentinlocus;
    uint32 aachangeinlocus;

    uint32 insertionuntranslated;    //
    uint32 deletionuntranslated;     //
    uint32 silentuntranslated;       //
    uint32 aachangeuntranslated;     //

    uint32 proteinlenchanged;
    std::string proteininfeature;

    std::string coveragestatus;
  };


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

  typedef std::map<int32, std::string> strainid2name_t;
  typedef strainid2name_t::value_type s_idname_pair_t;

  void saveContigReadList(std::list<Contig> & clist,
			  const std::string & filename,
			  bool deleteoldfile);
  void saveContigReadList(Contig & con,
			  const std::string & filename,
			  bool deleteoldfile);
  void saveStatistics(std::list<Contig> & clist,
		      const std::string & filename,
		      bool deleteoldfile);
  void saveStatistics(Contig & con,
		      const std::string & filename,
		      bool deleteoldfile);
  void saveAssemblyInfo(AssemblyInfo & asi,
		      const std::string & filename,
		      bool deleteoldfile);
  void saveLargeContigsInfo(AssemblyInfo & asi,
			    const std::string & filename,
			    bool deleteoldfile);
  void saveReadTagList(std::list<Contig> & clist,
		       const std::string & filename,
		       bool deleteoldfile);
  void saveReadTagList(Contig & con,
		       const std::string & filename,
		       bool deleteoldfile);
  void saveConsensusTagList(std::list<Contig> & clist,
			    const std::string & filename,
			    bool deleteoldfile);
  void saveConsensusTagList(Contig & con,
			    const std::string & filename,
			    bool deleteoldfile);
  void getPreviousLocusTag(const std::list<gbfsummary_t> & gbfs,
			   std::list<gbfsummary_t>::const_iterator gbfsI,
			   std::string & locustag,
			   std::string & gene);
  void getNextLocusTag(const std::list<gbfsummary_t> & gbfs,
		       std::list<gbfsummary_t>::const_iterator gbfsI,
		       std::string & locustag,
		       std::string & gene);
  void saveSNPList(std::list<Contig> & clist,
		   const std::string & filename,
		   bool deleteoldfile);
  void saveCoverageInfo(std::list<Contig> & clist,
			const std::string & filename,
			bool deleteoldfile);
  void saveSNPSurroundingAsHTML(std::list<Contig> & clist,
				const std::string & filename,
				bool deleteoldfile);
  void makeAllStrainGenomes(Contig & contigI,
			    base_quality_t minqual,
			    std::string & consseq,
			    std::vector<base_quality_t> & consqual,
			    std::vector<std::string> & strain_consseq,
			    std::vector< std::vector<base_quality_t> > & strain_consqual,
			    strainid2name_t & strainnames_in_contig,
			    bool fillholesinstraingenomes);
  void saveFeatureAnalysis(std::list<Contig> & clist,
			   ReadPool & rpool,
			   const std::string & faname,
			   const std::string & fsname,
			   const std::string & fcname,
			   bool deleteoldfile);
  void adjustCaseOfSequences(std::string & reference,
			     std::vector<std::string> & mutants);

  void saveStrainsAsFASTAQ(std::list<Contig> & clist,
			   const ReadPool & rp,
			   const std::string & paddedfilename,
			   bool asfastq,
			   uint32 mincoverage=0,
			   base_quality_t minqual=0,
			   bool deleteoldfile=false,
			   bool fillholesinstrain=false);
  void saveStrainsAsFASTAQ(Contig & outcon,
			   const ReadPool & rp,
			   const std::string & paddedfilename,
			   bool asfastq,
			   uint32 mincoverage=0,
			   base_quality_t minqual=0,
			   bool deleteoldfile=false,
			   bool fillholesinstrain=false);

  void saveStrainsAsGBF(std::list<Contig> & clist,
			const ReadPool & rp,
			const std::string & paddedfilename,
			base_quality_t minqual,
			bool fillholesinstraingenomes,
			bool deleteoldfile);

  void saveAs_TYPE(std::list<Contig> & clist,
		   const std::string & filename,
		   const uint8 type,
		   bool deleteoldfile);
  void saveAs_TYPE(Contig & con,
		   const std::string & filename,
		   const uint8 type,
		   bool deleteoldfile);
  inline void saveAsTCS(std::list<Contig> & clist,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(clist,filename,Contig::AS_TCS, deleteoldfile);
  };
  inline void saveAsTCS(Contig & con,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(con,filename,Contig::AS_TCS, deleteoldfile);
  };
  inline void saveAsCAF(std::list<Contig> & clist,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(clist,filename,Contig::AS_CAF, deleteoldfile);
  };
  inline void saveAsCAF(Contig & con,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(con,filename,Contig::AS_CAF,deleteoldfile);
  };
  inline void saveAsMAF(std::list<Contig> & clist,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(clist,filename,Contig::AS_MAF, deleteoldfile);
  };
  inline void saveAsMAF(Contig & con,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(con,filename,Contig::AS_MAF,deleteoldfile);
  };
  inline void saveAsTXT(std::list<Contig> & clist,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(clist,filename,Contig::AS_TEXT,deleteoldfile);
  };
  inline void saveAsTXT(Contig & con,
			const std::string & filename,
			bool deleteoldfile){
    saveAs_TYPE(con,filename,Contig::AS_TEXT, deleteoldfile);
  };
  void saveAsGAP4DA(std::list<Contig> & clist,
		    const std::string & dirname,
		    bool deleteoldfile);
  void saveAsGAP4DA(Contig & con,
		    const std::string & dirname,
		    bool deleteoldfile);
  void saveAsFASTA(std::list<Contig> & clist,
		   const std::string & filename,
		   const std::string & paddedfilename,
		   bool deleteoldfile);
  void saveAsFASTA(Contig & con,
		   const std::string & filename,
		   const std::string & paddedfilename,
		   bool deleteoldfile);
  void saveAsACE(std::list<Contig> & clist,
		 const std::string & filename,
		 bool deleteoldfile);
  void saveAsACE(Contig & con,
		 const std::string & filename,
		 bool deleteoldfile);
  void saveAsWiggle(std::list<Contig> & clist,
		    const std::string & filename,
		    bool deleteoldfile,
		    bool gcinsteadcov);
  void saveAsWiggle(Contig & con,
		    const std::string & filename,
		    bool deleteoldfile,
		    bool gcinsteadcov);


  void dumpContigs(std::list<Contig> & clist, std::ostream & fout);
  inline void dumpAsTCS(std::list<Contig> & clist, std::ostream & fout){
    Contig::setCoutType(Contig::AS_TCS);
    dumpContigs(clist,fout);
  };
  inline void dumpAsTCS(Contig & con, std::ostream & fout){
    Contig::setCoutType(Contig::AS_TCS);
    fout << con;
  };
  inline void dumpAsCAF(std::list<Contig> & clist, std::ostream & fout){
    Contig::setCoutType(Contig::AS_CAF);
    dumpContigs(clist,fout);
  };
  inline void dumpAsCAF(Contig & con, std::ostream & fout){
    Contig::setCoutType(Contig::AS_CAF);
    fout << con;
  };
  inline void dumpAsTXT(std::list<Contig> & clist, std::ostream & fout){
    Contig::setCoutType(Contig::AS_TEXT);
    dumpContigs(clist,fout);
  };
  inline void dumpAsTXT(Contig & con, std::ostream & fout){
    Contig::setCoutType(Contig::AS_TEXT);
    fout << con;
  };
  void dumpAsACE(std::list<Contig> & clist, std::ostream & fout);


  void dumpHTMLHeader(const std::string & projectname, std::ostream & htmlout);
  void dumpContigListAsHTML(std::list<Contig> & clist,
			    const std::string & filename,
			    bool deleteoldfile,
			    const std::string & projectname);
  void dumpContigAsHTML(Contig & con,
			const std::string & filename,
			bool deleteoldfile,
			const std::string & projectname);


  // ---------------------- internals
  void dumpSequenceAsGBF_priv(const std::string & seqname,
			      const std::string & strainname,
			      std::string & seq,
			      const std::string & consseq,
			      const std::list<gbfsummary_t> & featurelist,
			      std::ostream & fout);
  void dumpTextAsGBFValueLine_priv(const std::string & key,
				   const std::string & value,
				   std::ostream & fout);

  void saveAsACE_openACE(std::fstream & fio,
			 const std::string & filename,
			 bool deleteoldfile,
			 uint32 & numcontigs,
			 uint32 & numreads);
  void saveAsACE_rewriteHeader(std::fstream & fio, const uint32 numcontigs, const uint32 numreads);

}

#endif
