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

#ifndef _bas_manifest_h_
#define _bas_manifest_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "util/fileanddisk.H"
#include "mira/readgrouplib.H"


class Manifest
{
private:
  //Variables
public:
  struct manifestloadentry_t {
    ReadGroupLib::ReadGroupID rgid;
    std::list<std::string> datanames;
    bool loadasbackbone;
    bool sp_infoonlygiven;   // whether infoonly/exclusion_criterion was given by the user
    bool ts_infoonlygiven;   // whether infoonly/exclusion_criterion was given by the user

    std::list<fnft_t> mainfilesfoundfordata;
    std::list<fnft_t> ancillaryfilesfoundfordata;

    manifestloadentry_t() : loadasbackbone(false), sp_infoonlygiven(false), ts_infoonlygiven(false) {};
  };

  std::vector<manifestloadentry_t> MAN_manifestdata2load;

  //std::list<std::string> MAN_ancillary_datanames;
  //std::list<std::string> MAN_ancillary_datatypes;


private:
  std::string MAN_projectname;
  std::string MAN_job;
  std::string MAN_parameters;

  std::vector<bool> MAN_seentechnology;
  std::string MAN_technologystring;

  //Functions
private:
  void init();
  void zeroVars();

  void slurpInManifest(std::stringstream & mfin, const std::string & origsource, bool resume);

  bool provideFileNames();

public:
  Manifest();
  Manifest(Manifest const &other);
  ~Manifest();

  // use std copy-operator
  //Manifest const & operator=(Manifest const & other);
  friend std::ostream & operator<<(std::ostream &ostr, Manifest const &m);

  void discard();

  void loadManifestFile(const std::string & mfilename, bool resume);
  std::string getFullMIRAParameterString();
  const std::string & getProjectName() const {return MAN_projectname;}
  const std::string & getJob() const {return MAN_job;}
  const std::string & getParameters() const {return MAN_parameters;}
  const std::string & getTechnologies() const {return MAN_technologystring;}

};


#endif
