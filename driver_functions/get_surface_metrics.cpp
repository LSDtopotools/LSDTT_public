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
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../TNT/tnt.h"
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

	string DEM_name;
	string dem_ext = "_dem";
	vector<string> output_extensions;
	output_extensions.push_back("_elev");
	output_extensions.push_back("_slope");
	output_extensions.push_back("_aspect");
	output_extensions.push_back("_curv");
	output_extensions.push_back("_pl_curv");
	output_extensions.push_back("_pr_curv");
	output_extensions.push_back("_tan_curv");
	output_extensions.push_back("_sp");
  	string DEM_flt_extension = "flt";
	float surface_fitting_window_radius;
	vector<int> raster_selection;  	
	string temp;
	int raster_selector;                                                        
	file_info_in >> temp >> DEM_name
               >> temp >> surface_fitting_window_radius;
	// Now create the raster selection vector based on user's selection
	// Elevation
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Slope
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Aspect
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Curvature
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Planform Curvature
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Profile Curvature
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Tangential Curvature
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	// Classification of stationary points
        file_info_in >> temp >> raster_selector;
	raster_selection.push_back(raster_selector);
	file_info_in.close();

	LSDRaster dem(DEM_name+dem_ext, DEM_flt_extension);
	vector<LSDRaster> output_rasters;
	output_rasters = dem.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
	for(int i = 0; i<raster_selection.size(); ++i)
	{
		if(raster_selection[i]==1)
		{
			string output_file = DEM_name+output_extensions[i];
			output_rasters[i].write_raster(output_file,DEM_flt_extension);
		}
	}
}
