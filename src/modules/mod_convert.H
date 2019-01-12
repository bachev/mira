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

#include <fstream>

#include <string>
#include <vector>
#include <list>

#include "mira/parameters.H"
#include "mira/readpool.H"
#include "mira/contig.H"
#include "mira/assembly.H"
#include "mira/sam_collect.H"
#include "mira/gff_save.H"

class ConvPro
{
private:
  static std::vector<MIRAParameters> CP_Pv;

  static std::string CP_fromtype;
  static std::list<std::string> CP_totype;

  static std::list<std::ofstream *> CP_ofs;


  static std::string CP_infile;

  static std::string CP_inbasename;
  static std::string CP_outbasename;

  static std::string CP_renamesequences;
  static std::string CP_renamenamescheme;

  static bool CP_splitcontigs2singlefiles;

  static bool CP_deletestaronlycolumns;
  static bool CP_blinddata;
  static bool CP_fillholesinstraingenomes;
  static bool CP_makecontigs;
  static bool CP_extractreadsinsteadcontigs;
  static bool CP_hardtrim;
  static bool CP_trimnx;

  static std::string CP_namefile;
  static bool CP_keepnamesfromfile;
  static bool CP_sortbyname;

  static bool CP_mustdeletetargetfiles;

  static bool CP_specialtestcode;

  static bool CP_filter2readgroup;
  static bool CP_filter2readgroup_split;

  static base_quality_t CP_minqual;
  static bool CP_needsquality;
  static base_quality_t CP_defaultqual;

  static char CP_recalcconopt;
  static bool CP_recalcfeatureopt_snp;
  static bool CP_recalcfeatureopt_rep;

  static uint32 CP_minbasecoverage;

  static uint32 CP_mincontiglength;
  static bool   CP_minlengthisclipped;

  static uint32 CP_mincontigcoverage;
  static uint32 CP_minnumreads;

  static std::list<Contig> CP_clist;   // needed for CAF & MAF conversion (and GBF)
  static AssemblyInfo CP_assemblyinfo;

  static uint64 CP_readrenamecounter;
  static GFFSave CP_gffsave;
  static SAMCollect CP_samcollect;

  static uint64 CP_numclippedreadsinload; // number of reads which have clips already when loaded (CAF/MAF)
  static bool CP_ulcaseclips; // user wants change case of sequence according to left / right clips?
  static bool CP_mustcaseclips; // mirabait must change case clips for reads not in contigs

  static uint64 CP_yieldmax; // max bases to be written out
  static uint64 CP_yieldwritten; // currently written

private:
  static void usage();
  static void checkTypes(const std::string & fromtype,std::list<std::string> & totype);
  static void putReadsInContigsAndSave(std::vector<MIRAParameters> & Pv, ReadPool & rp);
  static void discardShortReads(std::vector<MIRAParameters> & Pv, ReadPool & rp, uint32 minlength, bool fromclipped);
  static void specialTestCode(std::list<Contig> & clist, ReadPool & rp);

  static void filterToReadGroup(ReadPool & rp);

  static bool contig__nameordercomp(const Contig & a, const Contig & b);
  static void sortContigsByName(std::list<Contig> & clist);
  static void sortPoolByName(ReadPool & rp, std::string & filename);

  static void saveContigList(std::list<Contig> & clist, ReadPool & rp);
  static void saveContigList_helper(std::list<Contig> & clist, ReadPool & rp);

  static void saveReadPool(ReadPool & rp, std::list<std::ofstream *> & ofs);
  static void cafmafload_callback(std::list<Contig> & clist, ReadPool & rp);
  static void readpoolload_callback(ReadPool & rp);
  static std::string createFileNameFromBasePostfixContigAndRead(const std::string & basename,
							   const char * postfix,
							   Contig * actcon = nullptr,
							   Read * actread = nullptr);
  static uint32 openOFSlist(Contig * optcontig, std::list<std::ofstream *> & ofs);
  static void closeOFSList(uint32 howmany, std::list<std::ofstream *> & ofs);

  static void openIFStream(std::ifstream & ifs, std::string fname, std::ios_base::openmode mode);
  static void openOFStream(std::ofstream & ofs, std::string fname, std::ios_base::openmode mode);

public:
  ~ConvPro();

  int mainConvPro(int argc, char ** argv);

  static void closeOpenStreams(std::list<std::ofstream *> & ofsl);

  static bool checkForFromType(const std::string & ftype);
  static bool checkForToType(const std::string & ttype);
  static void guessFromAndToType(const std::string & fnamefrom,
				 std::string & fromtype,
				 std::string * fromstem,
				 const std::string & fnameto,
				 std::list<std::string> & totypes,
				 std::string * tostem);
};

