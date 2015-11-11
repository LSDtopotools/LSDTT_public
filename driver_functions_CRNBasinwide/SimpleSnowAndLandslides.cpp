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
  if (nNumberofArgs!=4)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Basinwide cosmogenic analysis tool!  ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires three inputs: " << endl;
    cout << "* First the path to the parameter files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of the parameter files." << endl;
    cout << "* The third is a method flag. 0 does things without error analysis or muons," << endl; 
    cout << "  1 is a full analysis without spawned basins." << endl;
    cout << "  2 is a full analysis with spawned basins." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./SimpleSnowAndLandslides.exe home/fieldwork/Chile/CRN/ My_data" << endl;
    cout << "In windows (the slash directions will change and there is no leading ./)" << endl;
    cout << "SimpleSnowAndLandslides.exe c:\\fieldwork\\Chile\\CRN\\ My_data" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on cosmo data file format, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string DEM_name_prefix = argv[2];
  int method_flag = atoi(argv[3]);
  
  string DEM_f_name = path_name+DEM_name_prefix;
  string DEM_bil_extension = "bil";
  
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
  
  // update it with a bilinear snow function
  float SlopeAscend = 0.05;
  float SlopeDescend = -0.1;
  float PeakElevation = 3000;
  float PeakSnowpack = 80;
  SnowRaster.SetSnowEffDepthBilinear(SlopeAscend, SlopeDescend, PeakElevation, 
                                     PeakSnowpack, filled_topo_test);
                                 
  string snow_ext = "_SnowBL";
  string snow_out_name = path_name+DEM_name_prefix+snow_ext;
  SnowRaster.write_raster(snow_out_name,DEM_bil_extension);
  
  float v = 0.01;
  float lambda = 2000;
  float MaximumSlope = 0.1;
  SnowRaster.SetSnowEffDepthRichards(PeakSnowpack, MaximumSlope, v, lambda, filled_topo_test);
 
  string snow_ext2 = "_SnowBL2";
  snow_out_name = path_name+DEM_name_prefix+snow_ext2;
  SnowRaster.write_raster(snow_out_name,DEM_bil_extension);
  
  
  // now we go on to make a rudimetary landsliding raster
  // First get some sources
  int flow_threshold = 300;
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  vector<int> sources = FlowInfo.get_sources_index_threshold(ContributingPixels, flow_threshold);
  LSDRaster NaiveLanslide = FlowInfo.get_upslope_node_mask(sources);
  
  string LS_ext = "_LS";
  string LS_out_name = path_name+DEM_name_prefix+LS_ext;
  NaiveLanslide.write_raster(LS_out_name,DEM_bil_extension);  
}
  