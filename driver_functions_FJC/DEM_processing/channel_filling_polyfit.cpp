//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// channel_filling_polyfit.cpp
//
// This program fits a polynomial surface to the DEM with a specified window radius, and
// prints out the smoothed elevation raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb
// University of Edinburgh
// 30/01/17
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
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the DEM name, and the window radius for polyfitting." << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string DEM_name = argv[2];
	string window_radius = argv[3];

	cout << "The path is: " << path_name << "\n The DEM name is: " << DEM_name << "\n The window radius is: " << window_radius << endl;

	string raster_output = "_test";
	string elev_output = "_elev";
	string slope_output = "_slope";
	string HS_output = "_HS";
  string DEM_extension = "bil";
	vector<int> raster_selection(8,0);
	raster_selection[0] = 1;  // get the smoothed elevation
	string temp;

	// load the raster and remove values below 0
	LSDRaster dem(path_name+DEM_name, DEM_extension);
	dem.remove_seas();
	dem.write_raster((path_name+DEM_name+raster_output), DEM_extension);

	float surface_fitting_window_radius = strtof(window_radius.c_str(),0);

	vector<LSDRaster> output_rasters;
	output_rasters = dem.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);

	// smoothed elevation
	output_rasters[0].write_raster((path_name+DEM_name+elev_output), DEM_extension);
	// slope
	output_rasters[1].write_raster((path_name+DEM_name+slope_output), DEM_extension);

	// write smoothed hillshade
	LSDRaster HS = output_rasters[0].hillshade(45, 315, 1);
	HS.write_raster((path_name+DEM_name+elev_output+HS_output), DEM_extension);
}
