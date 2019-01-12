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

#ifndef _bas_readpoolio_h_
#define _bas_readpoolio_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "io/generalio.H"
#include "io/fasta.H"
#include "io/fastq-lh.H"
#include "mira/gbf_parse.H"
#include "mira/gff_parse.H"

class MAFParse;
class CAF;

#include "util/progressindic.H"


#include "mira/readpool.H"
#include "mira/contig.H"



class ReadPoolIO
{
private:
  KSEQ_INIT(gzFile, gzread);


  //Variables
private:
  ReadPool * RPIO_rpptr;
  std::list<Contig> * RPIO_clistptr;
  std::vector<MIRAParameters> * RPIO_mpptr;   // parameters for the contigs

   // 0: undefined; 1: not open; 2: open (loading); 3: closed (no more sequences)
  enum loadstatus {LS_UNDEFINED=0,LS_NOTOPEN,LS_OPEN,LS_CLOSED};
  uint32 RPIO_loadstatus;
  std::string RPIO_filetype;
  std::string RPIO_filename1;
  std::string RPIO_optfilename2;
  ReadGroupLib::ReadGroupID RPIO_rgid;
  bool RPIO_countonly;

  enum loadtype {LT_UNDEFINED=0,LT_FASTQ,LT_FASTA,LT_MAF,LT_CAF,LT_GBF,LT_GFF3,LT_FOFNEXP,LT_EXP};
  uint32 RPIO_loadtype;

  ProgressIndicator<int64> * RPIO_progressindic;
  uint64 RPIO_totalreadsloaded;
  std::streamsize RPIO_fsize;

  base_quality_t RPIO_fastq_qualoffset;  // 0 = external correction, q-values not corrected after loading!
  gzFile   RPIO_fastq_fp;
  kseq_t * RPIO_fastq_seq;
  bool RPIO_fastq_transformname;

  std::ifstream RPIO_fasta_fin;
  std::ifstream RPIO_fasta_qin;
  std::streamsize RPIO_fasta_fsize;
  std::string RPIO_missingfastaqual_resolvemsg;
  bool RPIO_fasta_wantsqualfiletoexist;
  bool RPIO_fasta_hasqualfile;
  FASTA RPIO_fasta_ioobj;

  bool RPIO_fastqa_preservecomments;
  bool RPIO_fastqa_checkquals;

  GBF RPIO_gbf_ioobj;
  bool RPIO_gbf_gbfloaded;
  uint32 RPIO_gbf_numtransferred;

  GFFParse RPIO_gff_ioobj;
  bool RPIO_gff_gffloaded;
  uint32 RPIO_gff_numtransferred;

  std::string RPIO_fofnexp_justpath;
  std::vector<std::string> RPIO_fofnexp_names;
  uint32 RPIO_fofnexp_nameiloaded;

  MAFParse * RPIO_maf_parse;

  CAF * RPIO_caf_parse;



public:


  //Functions
private:
  void priv_init();
  void priv_zeroVars();

  void priv_closeFiles();
  void priv_openFiles();
  void priv_openFiles_fastq();
  void priv_openFiles_fasta();
  void priv_openFiles_maf();
  void priv_openFiles_caf();
  void priv_openFiles_gbf();
  void priv_openFiles_gff3();
  void priv_openFiles_fofnexp();
  void priv_openFiles_exp();

  uint64 priv_loadNextSeqs_fastq(uint64 numseqstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_fasta(uint64 numseqstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_maf(uint64 numseqstoload, uint64 numconstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_caf(uint64 numseqstoload, uint64 numconstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_gbf(uint64 numseqstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_gff3(uint64 numseqstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_fofnexp(uint64 numseqstoload, uint64 lenseqstoload);
  uint64 priv_loadNextSeqs_exp();

public:
  ReadPoolIO(ReadPool & rp);
  ~ReadPoolIO();

  ReadPoolIO(ReadPoolIO const &other)=delete;
  ReadPoolIO const & operator=(ReadPoolIO const & other) = delete;

  friend std::ostream & operator<<(std::ostream &ostr, ReadPoolIO const &rpio);

  void setNewReadPool(ReadPool & rp);

  void discard();

  void registerFile(const std::string & filetype,
		    const std::string & filename1,
		    const std::string & optfilename2,
		    const ReadGroupLib::ReadGroupID rgid,
		    bool countonly);
  void registerFile(const std::string & filename1,
		    const std::string & optfilename2,
		    const ReadGroupLib::ReadGroupID rgid,
		    bool countonly);
  uint64 loadNextSeqs(uint64 numseqstoload, uint64 numconstoload=0, uint64 lenseqstoload=-1);

  void setAttributesForContigs(std::list<Contig>  * clist, std::vector<MIRAParameters> * mp);
  void setAttributeReadPool(ReadPool & rp) {RPIO_rpptr=&rp;}
  void setAttributeContigList(std::list<Contig> * clist) {RPIO_clistptr=clist;}
  void setAttributeProgressIndicator(bool b);

  void setAttributeFASTAQualFileWanted(bool b) {RPIO_fasta_wantsqualfiletoexist=b;}
  void setAttributeFASTAMissingQualFileResolveMsg(const std::string & s) {RPIO_missingfastaqual_resolvemsg=s;}

  void setAttributeFASTQQualOffset(base_quality_t q) {RPIO_fastq_qualoffset=q;}
  void setAttributeFASTQTransformName(bool b) {RPIO_fastq_transformname=b;}

  void setAttributeFASTQAPreserveComment(bool b) {RPIO_fastqa_preservecomments=b;}
  void setAttributeFASTQACheckQuals(bool b) {RPIO_fastqa_checkquals=b;}
};


#endif
