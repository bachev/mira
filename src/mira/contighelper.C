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

#include "contighelper.H"

#include "mira/contig.H"



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void dumpContigListAsHTML(std::list<Contig> & clist, const string & projectname, ostream & htmlout)
{

  htmlout << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\">\n\
<html>\n\
<head>\n\
   <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n\
   <meta name=\"GENERATOR\" content=\"MIRA / EdIt (c) Bastien Chevreux & Thomas Pfisterer;\">\n\
   <meta name=\"Author\" content=\"";
  //char * logname= getenv("LOGNAME");
  //htmlout << getlogin() << "\">\n\

  struct passwd * pws = getpwuid(getuid());
  bool noname=true;
  if(pws != nullptr) {
    if(pws->pw_name != nullptr && pws->pw_gecos != nullptr) {
      htmlout << pws->pw_name << " (" << pws->pw_gecos << ")";
      noname=false;
    }
  }
  if(noname) {
    char * namestr=getenv("LOGNAME");
    if(namestr!=nullptr) {
      htmlout << namestr;
      noname=false;
    }
    namestr=getenv("USER");
    if(namestr!=nullptr) {
      if(!noname) htmlout << " / ";
      htmlout << namestr;
      noname=false;
    }
  }

  if(noname){
    htmlout << "unknown";
  }

  htmlout << "\">\n<meta name=\"Description\" content=\"Assembled shotgun project\">\n\
   <title>";
  htmlout << "Project " << projectname << " </title>\n\
  <STYLE TYPE=\"text/css\">\n\
  <!--\n\
  \n\
   .MISM {color:black;  background-color:#20c020;}\n\
   .PRMB {color:black;  background-color:#ff5050;}\n\
   .WRMB {color:black;  background-color:blue;}\n\
   .PROS {color:black;  background-color:#00ced1;}\n\
   .PAOS {color:black;  background-color:#2e8b57;}\n\
   .PIOS {color:black;  background-color:#98fb98;}\n\
   .POLY {color:black;  background-color:#ffff99;}\n\
   .EDxD {color:black;  background-color:#db7093;}\n\
   .EDxI {color:black;  background-color:#db7093;}\n\
   .EDxC {color:black;  background-color:#db7093;}\n\
   .IUPC {color:black;  background-color:#cccccc;}\n\
\n\
  -->\n\
</STYLE>\n\
</head>\n\
<body TEXT=\"#000000\" BGCOLOR=\"#FFFFFF\" LINK=\"#FF0000\" VLINK=\"#551A8B\" ALINK=\"#000088\">\n";

  //   .ALUS {color:black;  background-color:#90ee90;}\n\
  // #66ffff = azure (pi*daumen)

  htmlout << "<h1><center>Tag legend</center></h1>\n";

  htmlout << "<center>\n";
  htmlout << "<table CELLSPACING=0 CELLPADDING=0 NOSAVE >\n";

  htmlout << "<tr><td><tt><SPAN CLASS=\"MISM\">&nbsp;</SPAN> = MISM;</tt> Mismatch (discrepancy) between reads and consensus</td></tr>\n";
  //htmlout << "<tr><td><tt><SPAN CLASS=\"ALUS\">&nbsp;</SPAN> = ALUS;</tt> Repetitive ALU sequence</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"PRMB\">&nbsp;</SPAN> = PRMB;</tt> Probable Repeat Marker Base set by MIRA</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"WRMB\">&nbsp;</SPAN> = WRMB;</tt> Weak Repeat Marker Base set by MIRA</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"PROS\">&nbsp;</SPAN> = PROS;</tt> Possible inteR Organism SNP set by MIRA</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"PAOS\">&nbsp;</SPAN> = PAOS;</tt> Possible intrA Organism SNP set by MIRA</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"PIOS\">&nbsp;</SPAN> = PIOS;</tt> Possible Inter- and intra-Organism SNP set by MIRA</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"POLY\">&nbsp;</SPAN> = POLY;</tt> Poly-A site</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"EDxD\">&nbsp;</SPAN> = EDxD;</tt> Delete operation set by EdIt</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"EDxI\">&nbsp;</SPAN> = EDxI;</tt> Insert operation set by EdIt</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"EDxC\">&nbsp;</SPAN> = EDxC;</tt> Change operation set by EdIt</td></tr>\n";
  htmlout << "<tr><td><tt><SPAN CLASS=\"IUPC\">&nbsp;</SPAN> = IUPAC;</tt> IUPAC base (shows only in HTML output)</td></tr>\n";

  htmlout<< "</table></center>\n";

  htmlout << "<h1><center>Contig List</center></h1>\n";
  auto I=clist.cbegin();
  while(I!=clist.cend()){

    if(I->getContigReads().size() > 1) {
      htmlout << "<a href=\"#" << I->getContigName() << "\">Contig " <<  I->getContigID();
    }else{
      htmlout << "<a href=\"#" << I->getContigName() << "\">Singlet " <<  I->getContigID();
    }
    htmlout << " (" << I->getContigSize() << ")</a>";
    ++I;
    if(I!=clist.cend()){
      htmlout << ", ";
    }else{
      htmlout << "\n<p>\n";
    }
  }

  Contig::setCoutType(Contig::AS_HTML);
  for(I=clist.cbegin();I!=clist.cend();++I){
    htmlout << *I;
  }

  htmlout << "\n</body></html>";

  return;
}
