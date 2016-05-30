//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// SimpleSnowAndLandslides.cpp
// This creates some simple snow and landsliding raster for use in testing 
// the cosmogenic functions
// 
// Developed by:
//  Simon M. Mudd, University of Edinburgh, School of GeoSciences
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
#include "../LSDSoilHydroRaster.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDStrahlerLinks.hpp"
#include "../LSDBasin.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDCosmoData.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
  // the driver version
  string driver_version = "Driver_version: 0.02";

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Simplified snow shielding tool!      ||" << endl;
    cout << "|| This program makes a snow shielding raster for use  ||" << endl;
    cout << "|| with the CRNBasinwide tools. The output is a snow   ||" << endl;
    cout << "|| shielding raster with annual average SWE values in  ||" << endl;
    cout << "|| units of effective depth (g per cm^2).              ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of sparam file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./SimpleSnowShield .exe home/fieldwork/Chile/CRN/ My_data" << endl;
    cout << "In windows (the slash directions will change and there is no leading ./)" << endl;
    cout << "SimpleSnowShield.exe c:\\fieldwork\\Chile\\CRN\\ My_data" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on cosmo data file format, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << "=========================================================" << endl;
    cout << " You WILL NEED a snow parameter file. " << endl;
    cout << "This file will have the same prefix as the DEM to be printed but the extension .sparam" << endl;
    cout << "It contains a method on the first line:" << endl;
    cout << "Richards OR Bilinear are the options" << endl;
    cout << "and the parameters for these two snow shielding functions" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string DEM_name_prefix = argv[2];
  
  string DEM_f_name = path_name+DEM_name_prefix;
  string DEM_bil_extension = "bil";
  
  string Param_f_name = path_name+DEM_name_prefix+".sparam";
  
  // load the snow parameters
  ifstream snowparamf;
  snowparamf.open(Param_f_name.c_str());
  string method_name;
  float SlopeAscend;
  float SlopeDescend;
  float PeakElevation;
  float PeakSnowpack;
  
  float v;
  float lambda;
  float MaximumSlope;
  
  // get the method
  snowparamf >> method_name;
  
  // get rid of any control characters
  method_name = RemoveControlCharacters(method_name);
  
  if (method_name == "Bilinear")
  {
    snowparamf >> SlopeAscend;
    snowparamf >> SlopeDescend;
    snowparamf >> PeakElevation;
    snowparamf >> PeakSnowpack;
  }
  else if (method_name == "Richards")
  {
    snowparamf >> v;
    snowparamf >>  lambda;
    snowparamf >>  MaximumSlope;  
  }
  else
  {
    cout << "You have not selected a valid method. Options are Bilinear and Richars" << endl;
    cout << "You selected " << method_name << endl;
    cout << "These names are case sensitive, and make sure there are no spaces on the line." << endl;
    exit(EXIT_SUCCESS);
  }
  
  
  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";

  // load the DEM
  LSDRaster topo_test(DEM_f_name, DEM_bil_extension);
  
  // remove the sea (seems to be required if gdal is used in places with nodata)
  topo_test.remove_seas();

  // get the filled file
  float MinSlope = 0.0001;
  cout << "Filling the DEM" << endl;
  LSDRaster filled_topo_test = topo_test.fill(MinSlope);
  
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
  
  // now make a snow raster
  LSDSoilHydroRaster SnowRaster(filled_topo_test);
  
  string snow_ext = "_SnowBL";
  string snow_out_name = path_name+DEM_name_prefix+snow_ext;
  
  // update it with a bilinear snow function
  if (method_name == "Bilinear");
  {
    SnowRaster.SetSnowEffDepthBilinear(SlopeAscend, SlopeDescend, PeakElevation, 
                                     PeakSnowpack, filled_topo_test);
                                 
    SnowRaster.write_raster(snow_out_name,DEM_bil_extension);
  }
  
  if (method_name == "Richards")
  {
    SnowRaster.SetSnowEffDepthRichards(PeakSnowpack, MaximumSlope, v, lambda, filled_topo_test);
    SnowRaster.write_raster(snow_out_name,DEM_bil_extension);
  }
}
  