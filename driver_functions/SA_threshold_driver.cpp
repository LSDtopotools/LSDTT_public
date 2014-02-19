//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// SA_threshold_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program takes the coordinates of field mapped channel heads and calculates their 
// drainage area and slope to work out a threshold value AS^2 which is used to predict 
// channel head locations across the rest of the landscape.  It returns two rasters: one of
// the position of the predicted channel heads ("DEM_name_sources") and one of the stream
// network ("DEM_name_SO").
// It takes an input file with the format: 
// X_coordinate Y_coordinate
//
// Developed by:
//  Fiona J. Clubb
//  Simon M. Mudd
//
// Copyright (C) 2013 Fiona J. Clubb and Simon M. Mudd 2013
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
#include "../LSDRasterSpectral.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDBasin.hpp"

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
	file_info_in >> DEM_name;
	float Minimum_Slope;
	file_info_in >> Minimum_Slope;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);

	// Set the no flux boundary conditions
  vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
		// get the filled file
	//cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
		// load the filled DEM
	//LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_flt_extension);
	//filled_topo_test.write_raster((DEM_f_name),DEM_flt_extension);
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	
	//getting the slope raster
  // Calculate polyfit coefficients
  Array2D<float> a;
  Array2D<float> b;
  Array2D<float> c;
  Array2D<float> d;
  Array2D<float> e;
  Array2D<float> f;
  float window_radius = 6;
  filled_topo_test.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster Slope = filled_topo_test.calculate_polyfit_slope(d, e);
	
	//string CP_name = "_CP";
	//ContributingPixels.write_raster((path_name+DEM_name+CP_name), DEM_flt_extension);
  
  //read in the channel head coordinates from file
  string string_filename;
	string filename = "channel_heads_vir";
	string dot = ".";
	string extension = "txt";
  string_filename = filename+dot+extension;
  ifstream coords_list;
  coords_list.open(string_filename.c_str());
  
  float X_coord, Y_coord;
  vector<float> X_coords;
  vector<float> Y_coords;
  
  while (coords_list >> X_coord >> Y_coord)
  {
     cout << "X Coord: " << X_coord << " Y Coord: " << Y_coord << endl;
     X_coords.push_back(X_coord);
     Y_coords.push_back(Y_coord);
  }
  
  //get the threshold from the mapped channel heads
  vector<int> NodeIndices = FlowInfo.get_sources_from_mapped_channel_heads(X_coords, Y_coords);
  
  float total_SA = 0;
  int n_observations = 0;
  //get the average contributing pixels for the channel heads
  for (unsigned int i = 0; i < NodeIndices.size(); i++)
  {
    int node = NodeIndices[i];
    int row,col;
    FlowInfo.retrieve_current_row_and_col(node,row,col);
    int contributing_pixels = FlowInfo.retrieve_contributing_pixels_of_node(node);
    float slope = Slope.get_data_element(row,col);
    float SA_product = contributing_pixels * (slope*slope);
    cout << "SA product: " << SA_product << endl;
    total_SA = total_SA + SA_product;
    ++n_observations;
  }
  
  cout << "total SA: " << total_SA << endl;
  float average_SA = total_SA/n_observations;
  cout << "threshold: " << average_SA << endl;   
  
  // use calculated threshold to get sources
  vector<int> sources;
	sources = FlowInfo.get_sources_slope_area(ContributingPixels, Slope, average_SA);
	cout << "Got the sources" << endl;
	
	string sources_raster_name = "_sources";
	LSDIndexRaster sources_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(sources);
	sources_raster.write_raster((path_name+DEM_name+sources_raster_name),DEM_flt_extension);

	// get the junction network and print out the stream order array
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);        
  cout << "Got the junction network " << endl;
  LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  cout << "got the SO array" << endl;
  
  string SO_name = "_SO";
	SOArray.write_raster((path_name+DEM_name+SO_name),DEM_flt_extension);      
}                                         
