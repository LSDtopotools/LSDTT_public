//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// shielding_driver.cpp
//
// Driver file for calculating topographic shielding functions for CRN production
//
// Code based on Codilean (2006):
// Codilean, A. T. (2006) Calculation of the cosmogenic nuclide production
//  topographic shielding scaling factor for large areas using DEMs
//  Earth Surf. Proc. Land. 31(6), pg 785-794
//  doi: 10.1002/esp.1336
//
// Developed by:
//  Simon M. Mudd, University of Edinburgh, School of GeoSciences
//  Stuart W.D. Grieve, University of Edinburgh, School of GeoSciences
//  Marie-Alice Harel, University of Edinburgh, School of GeoSciences
//  Martin D. Hurst, British Geological Survey
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
  // the driver version
  string driver_version = "Driver_version: 0.01";

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=5)
  {
    cout << endl;
    cout << "=====================================================================" << endl;
    cout << "|| Welcome to the Topographic Shielding tool!                      ||" << endl;
    cout << "|| This program is used to calculate topographic shielding_driver  ||" << endl;
    cout << "|| rasters following the method of Codilean (2006).                ||" << endl;
    cout << "=====================================================================" << endl;
    cout << "This program requires four inputs: " << endl;
    cout << "* First the path to the DEM files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of the DEM (that is, without the .bil)." << endl;
    cout << "   For example if you DEM is Ladakh.bil, you should enter Ladakh" << endl;
    cout << "   (Note that the DEM should be in *.bil format)" << endl;
    cout << "* Third the increment in azimuth (in degrees) over which you" << endl;
    cout << "   want shielding calculated. Recommended values is 5" << endl;
    cout << "* Fourth the increment in inclination (in degrees) over which you" << endl;
    cout << "   want shielding calculated. Recommended values is 5" << endl;
    cout << "=====================================================================" << endl;
    cout << "For more documentation, see readme and online documentation" << endl;
    cout << "=====================================================================" << endl;
    cout << endl;
    exit(EXIT_SUCCESS);
  }

  cout << endl;
  cout << "===============================================================" << endl;
  cout << "Welcome to the Topographic Shielding tool" << endl;
  cout << "This software was developed at the University of Edinburgh," << endl;
  cout << "by the Land Surface Dynamics group. For questions email" << endl;
  cout << "simon.m.mudd _at_ ed.ac.uk" << endl;
  cout << "This software is released under a GNU public license." << endl;
  cout << "You are using " << driver_version << endl;
  cout << "================================================================" << endl;
  cout << "++IMPORTANT++ The DEM must be an ENVI bil format file" << endl;
  cout << "ENVI bil files are required because, unlike asc or flt files, " << endl;
  cout << "they use georeferencing information." << endl;
  cout << "For more information about changing DEM formatting, see: " << endl;
  cout << "http://lsdtopotools.github.io/LSDTT_book/#_gdal_2" << endl;
  cout << "================================================================" << endl;


  //Assign values from input arguments
  string path_name = argv[1];
  string file_name = argv[2];
  string azimuth_str = argv[3];
  string inclination_str = argv[4];

  //Set the DEM filename
  file_name = path_name+file_name;

  //Convert the angle increments to integers
  int azimuth_step =  atoi(azimuth_str.c_str());
  int inclination_step =  atoi(inclination_str.c_str());

  // check the parameter values
  check_azimuth_and_inclination(azimuth_step,inclination_step);

  //load dem
  LSDRaster DEM(file_name, "bil");

  //launch toposhielding
  LSDRaster TopoShielding = DEM.TopographicShielding(azimuth_step,inclination_step);
  TopoShielding.write_raster(file_name+"_TopoShield","bil");

  cout << "Done!" << endl << endl;
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
