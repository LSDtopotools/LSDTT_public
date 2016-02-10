//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// sediment_routing.cpp
// A program to compute sediment routing
//
// The function take two arguments to main
// The first is the data folder. This is the folder in which the parameter
// files are stored
// The second is the prefix of the parameter files
//
// Developed by:
//  Simon M. Mudd, University of Edinburgh, School of GeoSciences
//  Stuart W.D. Grieve, University of Edinburgh, School of GeoSciences
//  Marie-Alice Harel, University of Edinburgh, School of GeoSciences
//  Martin D. Hurst, British Geological Survey
//
//
// Copyright (C) 2015 Simon M. Mudd 2015
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "../LSDRaster.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDStrahlerLinks.hpp"
#include "../LSDBasin.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDCosmoData.hpp"
#include "../LSDRasterAggregator.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
  // the driver version
  string driver_version = "Driver_version: 0.02";

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "================================================================" << endl;
    cout << "|| Welcome to the routing tool!							   ||" << endl;
    cout << "================================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of the parameter files." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "There must be two parameter files in the named path." << endl;
    cout << "The first MUST have the extension .rasters" << endl;
    cout << " This contains the filenames (including full path) of" << endl;
    cout << " names of the DEMs (without file extensions) to be analysed." << endl;
    cout << "The second file contains the parameter data" << endl;
    cout << "  and MUST have the extension .params" << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "For example if you parameter files are in /home/Sediment_routing/" << endl;
    cout << " in this folder you will have two parameter files, for example" << endl;
    cout << " Corsica_analysis.rasters and Corsica_analysis.param" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./sediment_routing.exe /home/Sediment_routing/ Corsica_analysis" << endl;
    cout << "In windows (the slash directions will change and there is no leading ./)" << endl;
    cout << "sediment_routing.exe c:\\home\\Sediment_routing\\ Corsica_analysis" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string param_name_prefix = argv[2];

  // now load the CRNCosmoData object
  cout << "I'm gonna load some data for you buddy" << endl;
  LSDSedimentRouting SedimentRouter(path_name,param_name_prefix);
  cout << "Got the data" << endl;



}
