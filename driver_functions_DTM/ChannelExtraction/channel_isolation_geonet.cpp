//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// threshold_quartile_quantile.cpp
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// David T Milodowski
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../TNT/tnt.h"
int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

  string full_name = path_name+f_name;

  ifstream file_info_in;
  file_info_in.open(full_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  string Raster_name;
  string Output_name;
  string q_q_filename_prefix;
  int timesteps;  
  float area_threshold,window_radius;
  string DEM_extension = "flt";
  string temp;
  file_info_in >> temp >> Raster_name
               >> temp >> Output_name
               >> temp >> q_q_filename_prefix 
               >> temp >> timesteps
               >> temp >> window_radius
               >> temp >> area_threshold;
  file_info_in.close();
  // Now create the raster selection vector based on user's selection
  // Elevation
  LSDRaster raster(Raster_name, DEM_extension);
  LSDIndexRaster output_raster = raster.IsolateChannelsGeonet(timesteps, area_threshold, window_radius, q_q_filename_prefix+".txt");
  output_raster.write_raster(Output_name,DEM_extension);
  cout << "DONE" << endl;
}
