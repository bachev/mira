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


#include <getopt.h>

#include <iostream>
#include <string>
#include <algorithm>

using namespace std;


const char step1par[] = {
//#include "me_step1.par.H"
};
const char step2par[] = {
//#include "me_step2.par.H"
};
const char step3par[] = {
//#include "me_step3.par.H"
};




void dumpStdMsg()
{
  cout <<
    "To (un-)subscribe the MIRA mailing lists, see:\n"
    "\thttp://www.chevreux.org/mira_mailinglists.html\n\n"
    "After subscribing, mail general questions to the MIRA talk mailing list:\n"
    "\tmira_talk@freelists.org\n\n"
    "\nTo report bugs or ask for features, please use the SourceForge ticketing\nsystem at:\n"
    "\thttp://sourceforge.net/p/mira-assembler/tickets/\n"
    "This ensures that requests do not get lost.\n\n\n";

  bool addnl=false;

  //if(sizeof(size_t) == 4){
  //  cout << "Compiled in 32 bit mode.\n";
  //}else if(sizeof(size_t) == 8){
  //  cout << "Compiled in 64 bit mode.\n";
  //}else{
  //  cout << "Compiled in ??? bit mode.\n";
  //}

  cout << compileinfo;
#ifdef CEBUGFLAG
  cout << "Compiled in debug output mode.\n";
  addnl=true;
#endif
#ifdef TRACEFLAG
  cout << "Compiled with trace mode.\n";
  addnl=true;
#endif
#ifdef BOUNDTRACKFLAG
  cout << "Compiled in boundtracking mode.\n";
  addnl=true;
#endif
#ifdef BUGTRACKFLAG
  cout << "Compiled in bugtracking mode.\n";
  addnl=true;
#endif
#ifdef PARANOIABUGTRACKFLAG
  cout << "Compiled in paranoia bugtracking mode.\n";
  addnl=true;
#endif
#ifdef ENABLE64
  cout << "Compiled with ENABLE64 activated.\n";
  addnl=true;
#else
  cout << "Compiled with ENABLE64 de-activated.\n";
  addnl=true;
#endif
#ifdef MIRAMEMORC
  cout << "Compiled with memory overrun checks, MIRA *will* be slower.\n";
  addnl=true;
#endif

  cout << "Runtime settings (sorry, for debug):"
       << "\n\tSize of size_t  : " << sizeof(size_t)
       << "\n\tSize of uint32  : " << sizeof(uint32)
       << "\n\tSize of uint32_t: " << sizeof(uint32_t)
       << "\n\tSize of uint64  : " << sizeof(uint64)
       << "\n\tSize of uint64_t: " << sizeof(uint64_t)
       << "\nCurrent system: ";
  {
    cout.flush();
    int tmp=system("uname -a");
    // don't complain about unused variable
    (void) tmp;
  }

  if(addnl) cout << endl;
}



void miraESTstep1(vector<MIRAParameters> & Pv, const string & csfile, int argc, char ** argv)
{
  cout << "De-activated atm, sorry." << endl;
  exit(0);

//  // step 1
//  // assemble _all_ reads to squeeze out most of the data
//  // assembled reads will be written per strain to a CAF file
//  //  if they're in a contig >= 2 reads or have SRMr/WRMr tags
//
//  // by setting the AssumeSNPInsteadofRepeats to false, we will have
//  //  SRMB tags set and the assembler will try to resolve what
//  //  it thinks are repeats.
//
//  vector<string> strainfiles;
//
//  //loadparams(Pv,"me_step1.par", &step1par[0]);
//
//  const_cast<assembly_parameters &>(Pv[0].getAssemblyParams()).as_put_asswithmira_tags=true;
//
//  cout << "Starting step 1: assembling all reads of all strains." << endl;
//  cout << "Setting some standard parameters.\n";
//  //as.assemble();
//
//  MIRAParameters::postParsingChanges(Pv);
//  MIRAParameters::dumpAllParams(Pv, cout);
//
//  Assembly as(Pv, false);
//  dumpCommandlineToFile(as, Pv, argc, argv);
//
//  as.loadSequenceData();
//
//  if(Pv[0].getAssemblyParams().as_filecheck_only==true){
//    cout << "File check only was selected, exiting.";
//  }
//  strainfiles=as.assembleESTs();
//  as.discard();
//  cout << "Finished step 1" << endl;
//
//  ofstream fout(csfile, ios::out| ios::trunc);
//
//  for(uint32 i=0; i<strainfiles.size(); i+=2){
//    fout << strainfiles[i] << "\t" << strainfiles[i+1] << endl;
//  }
//  fout.close();
}


uint32 mests2_contigcount=0;
uint32 mests2_straini=0;
ofstream mests2_step2out;
ofstream mests2_step2strout;
vector<string> mests2_strainfiles;

void miraESTstep2_contigBuiltCallback(Contig & con, const ReadPool & rp)
{
  ++mests2_contigcount;

  ostringstream ostr;
  if(con.getNumReadsInContig() > 1){
    ostr << mests2_strainfiles[mests2_straini+1] << "_c" << mests2_contigcount;
  } else {
    ostr << mests2_strainfiles[mests2_straini+1] << "_s" << mests2_contigcount;
  }

  Read dummy;

  string cons;
  vector<base_quality_t> quals;
  con.newConsensusGet(cons, quals);
  dummy.setSequenceFromString(cons);
  dummy.setQualities(quals);

  //cout << "My consensus is: " << cons << endl;
  //cout << "My read is: " << dummy << endl;

  for(uint32 j=0; j<cons.size(); j++){
    if(cons[j]=='*') dummy.changeBaseInSequence('*', 0, j);
  }
  con.transposeReadSRMTagsToContig();
  {
    const vector<Contig::consensustag_t> & ctags=con.getConsensusTags();
    for(uint32 j=0; j<ctags.size(); j++){
      dummy.addTagO(ctags[j]);
    }
  }

  {
    Contig::cccontainer_t & concounts=const_cast<Contig::cccontainer_t&>(con.getConsensusCounts());
    size_t lbound=0;
    size_t rbound=concounts.size();
    auto ccI=concounts.begin();
    for(; ccI!=concounts.end(); lbound++, ccI++){
      if(ccI->total_cov > 1) break;
    }

    ccI=concounts.end();
    for(; ccI!=concounts.begin(); rbound--){
      if((--ccI)->total_cov>1) break;
    }


    //cout << "Clip bounds: " << lbound << "\t" << rbound << endl;

    //Read::setCoutType(Read::AS_TEXTSHORT);
    //cout << "Before clipping:\n" << dummy;

    // single reads will have reversed bounds, but we still
    // want them as they had some marks in the combined assembly,
    // so don't set bound for them as it would get them
    // completely removed from assembly

    //if(lbound>rbound) {
    //  dummy.setClipoffs(0, 1, false);
    //}else{
    //  dummy.setClipoffs(lbound, rbound, false);
    //}

    //Read::setCoutType(Read::AS_TEXTSHORT);
    //cout << "After clipping:\n" << dummy;

  }

  dummy.removeGapsFromRead();
  Read::setCoutType(Read::AS_CAF);
  mests2_step2out << dummy;
  mests2_step2strout << dummy.getName() << "\t"<< mests2_strainfiles[mests2_straini+1] << '\n';

}

void miraESTstep2(vector<MIRAParameters> & Pv, const string & csfile, int argc, char ** argv)
{
  cout << "De-activated atm, sorry." << endl;
  exit(0);

//
//
//  // assemble each strain for itself taking only the good
//  //  reads identified in the previous step
//  // the resulting contigs are transformed into virtual 'reads':
//  //  single coverage at the ends is clipped (exception: reads that
//  //  are completely in single coverage), tags are taken
//  //  and together with virtual base quality the 'read' gets
//  //  written to a CAF file
//
//  //loadparams(Pv,"me_step2.par", &step2par[0]);
//
//  Pv[0].setAssemblyPutAssembledWithMIRATags(false);
//
//  cout << "Starting step 2: assembling each strain for itself" << endl;
//
//  MIRAParameters::postParsingChanges(Pv);
//  MIRAParameters::dumpAllParams(Pv, cout);
//
//  // Load the file of with caf names and strain names
//  ifstream fin(csfile, ios::in|ios::ate);
//  if(!fin){
//    throw Notify(Notify::FATAL, "main", (static_cast<std::string>("File not found: ")+csfile.c_str()).c_str());
//  }
//
//  auto len_fofn=fin.tellg();
//  if(len_fofn==1){
//    throw Notify(Notify::FATAL, "main", (static_cast<std::string>("Zero length file: ")+csfile.c_str()).c_str());
//  }
//  fin.seekg(0, ios::beg);
//
//  string filename, sname;
//  while(GeneralIO::readKeyValue(fin, filename, sname)){
//    mests2_strainfiles.push_back(filename);
//    mests2_strainfiles.push_back(sname);
//  }
//  fin.close();
//
//
//
//  mests2_step2out.open("step2_reads.caf", ios::out | ios::trunc);
//  mests2_step2strout.open("step2_straindata_in.txt", ios::out | ios::trunc);
//
//  mests2_step2strout << "# Automatically generated file" << endl;
//  mests2_step2strout << "# You probably don't want to edit it.\n" << endl;
//  for(mests2_straini=0; mests2_straini<mests2_strainfiles.size(); mests2_straini+=2){
//
//    cout << "Assembly of strain " << mests2_strainfiles[mests2_straini] << "(" << mests2_strainfiles[mests2_straini+1] << ")" << endl;
//
//    Pv[0].generateProjectNames(Pv,"step2_"+mests2_strainfiles[mests2_straini+1]);
//    Pv[0].setAssemblyInfileCAF(const_cast<char *>(mests2_strainfiles[mests2_straini].c_str()));
//    Pv[0].setAssemblyInfileStrainData("step2_straindata_in.txt");
//
//    //P.setAssemblyOutfileCAF((char *)strainfiles[i+1].c_str());
//    //P.setAssemblyOutdirGAP4DA(strainfiles[i+1].c_str());
//    //cout << P;
//
//    Assembly as(Pv, false);
//
//    bool loadok=true;
//    try{
//      as.loadSequenceData();
//    }
//    catch(Notify n){
//      loadok=false;
//      n.setGravity(Notify::WARNING);
//      n.handleError(" miraESTstep2()");
//    }
//    if(!loadok) continue;
//
//    as.setContigBuiltCallback(miraESTstep2_contigBuiltCallback);
//    as.assemble();
//    //as.saveResults();
//
//    cout << "Finished assembly, extracting contigs." << endl;
//
//  }
//  cout << "Closing step2out." << endl;
//  mests2_step2out.close();
//  cout << "Closing step2strout." << endl;
//  mests2_step2strout.close();
//
//  cout << "Done with step 2." << endl;
}


void miraESTstep3(vector<MIRAParameters> & Pv, int argc, char ** argv)
{
  cout << "De-activated atm, sorry." << endl;
  exit(0);

//  // in the last step, the virtual 'reads' get assembled
//  // by setting AssumeSNPInsteadofRepeats to true, PALV (possible
//  //  allellic variation) or PAVS (possible allellic variation
//  //  with SNP) will be set when conflicts occur instead of SRMB
//  // the assembler will therefore work in cluster-mode and not
//  //  break those 'misassemblies'
//
//  cout << "Starting step 3:\n\tclustering contigs\n\tfinding possible allelic variances\n\tfinding possible allelic variances with SNP\n\n";
//
//  //loadparams(Pv,"me_step3.par", &step3par[0]);
//
//  const_cast<assembly_parameters &>(Pv[0].getAssemblyParams()).as_put_asswithmira_tags=true;
//
//  Assembly as(Pv, false);
//
//  as.loadSequenceData();
//  as.assemble();
//  as.saveResults();
}

void miraEST(int argc, char ** argv)
{
  cout << "De-activated atm, sorry." << endl;
  exit(0);

//  cout << "This is miraEST "MIRAVERSION" for EST SNP analysis in strains.\n\n";
//
//  cout << "De-activated atm, step 2&3 need to adapt to new loading system, sorry." << endl;
//  exit(0);
//
//  cout << "Please cite: Chevreux, B., Pfisterer, T., Drescher, B., Driesel, A. J.,\nMueller, W. E., Wetter, T. and Suhai, S. (2004),\nUsing the miraEST Assembler for Reliable and Automated mRNA Transcript\nAssembly and SNP Detection in Sequenced ESTs. Genome Research, 14(6).\n\n";
//
//  //cout << "miraEST has been de-activated in this development version as necessary adaptations there have not been made yet, sorry.\n";
//  //doAbort();
//
//  dumpStdMsg();
//
//  try{
//    vector<MIRAParameters> X;
//    //loadparams(X,"", &step1par[0]);
//    //loadparams(X,"", &step2par[0]);
//    //loadparams(X,"", &step3par[0]);
//  }
//  catch(...){
//    cout << "Internal error: one of the default parameter files caused an error while parsing, aborting.\n";
//    exit(1000);
//  }
//
//  vector<MIRAParameters> Pv;
//  MIRAParameters::setupStdMIRAParameters(Pv);
//  {
//    MIRAParameters::parse(argc, argv, Pv);
//    cout << "\nParameters parsed without error, perfect.\n\n";
//    MIRAParameters::postParsingChanges(Pv);
//  }
//
//  string csfile="step1_res_cafstrainnames.txt";
//
//  uint32 startstep=Pv[0].getSpecialParams().sp_est_startstep;
//
//  switch(startstep) {
//    case 1 : {
//      miraESTstep1(Pv, csfile, argc, argv);
//     break;
//    }
//    case 2 : {
//      miraESTstep2(Pv, csfile, argc, argv);
//      break;
//    }
//    case 3 : {
//      miraESTstep3(Pv, argc, argv);
//      break;
//    }
//  default : {
//    throw Notify(Notify::FATAL, "main", ": miraEST SNP pipeline.step is not 1,2 or 3");
//  }
//  }
//
//  cout << "\n\nEnd of assembly process, thank you for using miraEST." << endl;
//
  return;
}

