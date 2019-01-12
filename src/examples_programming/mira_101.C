/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2007 and later by Bastien Chevreux
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

// 	$Id$	


/*****************************************************************
 *
 * Example mira_101
 *
 * Simply copy the behaviour of the "mira" executable (without
 *  miraEST part): perform an assembly as defined by the user on 
 *  the command line
 *
 * Topics covered:
 * - always enclose usage of MIRA classes with a try/catch block
 *   (at least on the main program level)
 * - simple usage of the the MIRAParameter class
 * - simple usage of the Assembly class
 *
 *****************************************************************/


// include C++ iostreams
#include <iostream>

// include the MIRA assembly class
#include "mira/assembly.H"


int main(int argc, char ** argv)
{
  cout << "Example program mira_101: how to make an assembly.\n";

  // ALWAYS use a try/catch block when using MIRA library classes
  //  the way it is shown here. On nasty problems, the classes will
  //  throw exceptions to abort.

  try{

    // MIRAParameters contains all parameters the MIRA classes use,
    //  (especially the ones performing computations) so define one
    // With the advent of different sequencing types, this has changed
    //  to a vector of MIRAParameters

    vector<MIRAParameters> Pv;

    // setup standard values for these parameters
    // you MUST do this

    MIRAParameters::setupStdMIRAParameters(Pv);

    // Parse MIRA parameters from the command line into the class
    // One can also load them from a file via: P.loadParams("somefile")

    MIRAParameters::parse(argc, argv, Pv);

    // dump the MIRAParameters to stdout so that the user knows
    //  with which parameters this assembly will run

    cout << "\nParameters parsed without error, perfect.\n\n";
    MIRAParameters::dumpAllParams(Pv, cout);

    // once the parameters are set, create an instance of the Assembly 
    //  class, passing the parameter vector as reference
    Assembly as(Pv, false);

    // tell class to load all necessary data as instructed by parameters
    as.loadSequenceData();

    // tell class to assemble data as instructed by parameters
    as.assemble();

    // tell class to save final results as instructed by parameters
    as.saveResults();
  }
  // Catch exceptions as thrown by MIRA classes
  catch(Notify n){
    n.handleError("main");
    return 1;
  }
  catch(Flow f){
    cerr << "INTERNAL ERROR: Unexpected exception: Flow()\n";
    return 1;
  }
  // you can define own catches here
  catch(...){
    cerr << "Unknown exception caught, aborting the process.\n\n";
    return 1;
  }

  return 0;
}
