//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Chile_test.cpp
//
// This program is used for testing the LSDRaster object
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Simon M. Mudd
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../TNT/tnt.h"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string DEM_name = argv[2];
	string window_radius = argv[3];

	cout << "The path is: " << path_name << "\n The DEM name is: " << DEM_name << "\n The window radius is: " << window_radius << endl;

	string raster_output = "_test";
	string elev_output = "_elev";
  string DEM_extension = "bil";
	vector<int> raster_selection(8,0);  	
	raster_selection[0] = 1;  // get the smoothed elevation
	string temp;
	
	float surface_fitting_window_radius = strtof(window_radius.c_str(),0);
	                                                      
	LSDRaster dem(path_name+DEM_name, DEM_extension);
	vector<LSDRaster> output_rasters;
	output_rasters = dem.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
	
	output_rasters[0].write_raster((path_name+DEM_name+elev_output), DEM_extension);
	dem.write_raster((path_name+DEM_name+raster_output), DEM_extension);
}