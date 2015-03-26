//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get_drainage_density_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets the drainage density of all basins with a threshold area of
// 0.5 km^2 to compare drainage density with mean basin slope and mean hilltop
// curvature (used as a proxy for erosion rate). It gets drainage density based on a channel
// network methodology described by Clubb et al. (2014). The data are also binned (bin width = 0.01)
//
// Outputs a txt file for the data cloud with the format
// drainage_density mean_hilltop_curvature mean_basin_slope 
//
// Outputs a txt file for the binned data with the format
// binned_CHT CHT_standard_deviation CHT_standard_error binned_DD DD_standard_deviation
// DD_standard_error binned_slope slope_standard_deviation slope_standard_error
//
// Developed by:
//  Fiona J. Clubb
//  Simon M. Mudd
//  Stuart W. D. Grieve
//
// Copyright (C) 2013 Fiona Clubb and Simon M. Mudd 2013
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
#include <iomanip>
#include <math.h>
#include <string.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDBasin.hpp"
#include "../LSDCRNParameters.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
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
  string sources_name; 
	string fill_ext = "_fill";
	file_info_in >> DEM_name >> sources_name;
	float Minimum_Slope;
	file_info_in >> Minimum_Slope;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// Set the no flux boundary conditions
  vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
		// get the filled file
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_flt_extension);
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
  
  //get the sources from raster to vector
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+sources_name),DEM_flt_extension);
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name = "_SO";
	SOArray.write_raster((path_name+DEM_name+SO_name), DEM_flt_extension);
	
	// get all the basins greater than 200,000 m^2
	float AreaThreshold = 500000;
	vector<int> basin_nodes = ChanNetwork.extract_basin_nodes_by_drainage_area(AreaThreshold, FlowInfo);
	vector<int> basin_junctions = ChanNetwork.extract_basin_junctions_from_nodes(basin_nodes, FlowInfo);
	  
  cout << "Writing basin junctions to text file" << endl;
  
  // creating file to write junctions to
  string string_filename;
  string dot = ".";
  string extension = "txt";
  string filename = "_DD_junctions";
  string_filename = DEM_name+filename+dot+extension;
  ofstream DD_junctions;
  DD_junctions.open(string_filename.c_str());
  
  for (int i = 0; i < int(basin_junctions.size()); i++)
  {
    DD_junctions << basin_junctions[i] << endl;
  }
  
  // get flow directions - D8 flowdir for drainage density calculations 
  Array2D<int> FlowDir_array = FlowInfo.get_FlowDirection();
 
  //get the curvature from polynomial fitting and write to a raster
  
  float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
  vector<LSDRaster> surface_fitting;
  string curv_name = "_curv";
  vector<int> raster_selection(8, 0);
  raster_selection[3] = 1;                      // this indicates you want the curvature
  
  surface_fitting = filled_topo_test.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
  LSDRaster Curvature = surface_fitting[3]; 
  //Curvature.write_raster((path_name+DEM_name+curv_name), DEM_flt_extension);
  
  // Get the hilltop curvature
  LSDRaster Hilltops = ChanNetwork.ExtractRidges(FlowInfo);
  LSDRaster CHT_temp = filled_topo_test.get_hilltop_curvature(Curvature, Hilltops);
  LSDRaster CHT = filled_topo_test.remove_positive_hilltop_curvature(CHT_temp);
  string CHT_name = "_CHT";
  CHT.write_raster((path_name+DEM_name+CHT_name), DEM_flt_extension);
  
}