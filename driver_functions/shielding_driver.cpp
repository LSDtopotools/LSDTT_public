//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//shielding_driver.cpp
//
//Driver file for calculating topographic shielding functions for CRN production
//Testing the implementation of drop shadows/cast shadows
// Code based on:
// Calculation of the cosmogenic nuclide production topographic shielding 
//  scaling factor for large areas using DEMs
//  AT Codilean 
//  Earth Surface Processes and Landforms, 2006 
//
//1Martin D. Hurst  and 2Simon M. Mudd
//1British Geological Survey
//2University of Edinburgh
// Copyright (C) 2015 Martin D. Hurst
//
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDChannel.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDMostLikelyPartitionsFinder.hpp"
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"
using namespace std;

// delcartion for function that makes sure azimuth and inclination have correct
// increment values
void check_azimuth_and_inclination(int& theta_step, int& phi_step);

int main (int nNumberofArgs,char *argv[])
{

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=5)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the totgraphic shielding tool!  ||" << endl;
    cout << "=========================================================" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"  << endl;
    cout << "You MUST provide a DEM in ENVI bil format!"  << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"  << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires four inputs: " << endl;
    cout << "* First the path to the DEM files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of the DEM (that is, without the bil." << endl;
    cout << "   For example if you DEM is Ladakh.bil, you should enter Ladakh" << endl;
    cout << "* Third the increment in azimuth (in degrees) over which you" << endl;
    cout << "   want shielding calculated. Recommended values is 5" << endl;
    cout << "* Fourth the increment in inclination (in degrees) over which you" << endl;
    cout << "   want shielding calculated. Recommended values is 5" << endl;    
    cout << "=========================================================" << endl;
    cout << "This code is ditributed under a GNU lisence by Martin D. Hurst" << endl;
    cout << " at the British Geological Survey and Simon M. Mudd at " << endl;
    cout << " the University of Edinburgh" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string file_name = argv[2];
  string azimuth_str = argv[3];
  string inclination_str = argv[4];
  
  file_name = path_name+file_name;
  
  int azimuth_step =  atoi(azimuth_str.c_str());
  int inclination_step =  atoi(inclination_str.c_str());
  
  // check the parameter values
  check_azimuth_and_inclination(azimuth_step,inclination_step);

  //load dem
  LSDRaster DEM(file_name, "bil");  

  //launch toposhielding
  LSDRaster TopoShielding = DEM.TopographicShielding(azimuth_step,inclination_step);
  TopoShielding.write_raster(file_name+"_TopoShield","bil");
}


// make sure the azimuth and inclinations are factors of 360 and 90, respectively
void check_azimuth_and_inclination(int& theta_step, int& phi_step)
{
  // now check the phi and theta values. These must be a factor of 360 and 90, respectively
  int temp_step;
  if (360%theta_step != 0)
  {
    temp_step = theta_step-1;
    while(360%temp_step != 0)
    {
      temp_step--;
    }
    cout << "Theta was not a factor of 360, changing from: " << theta_step << endl;
    theta_step = temp_step;
    cout << " to: " << theta_step << endl;
  }
  if (90%phi_step != 0)
  {
    temp_step = phi_step-1;
    while(90%temp_step != 0)
    {
      temp_step--;
    }
    cout << "Phi was not a factor of 90, changing from: " << phi_step << endl;
    phi_step = temp_step;
    cout << " to: " << phi_step << endl;
  }
}
