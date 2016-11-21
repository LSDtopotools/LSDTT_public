//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// basin_puncher_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets all the basins of a specified drainage area and prints them
// to a raster.
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
#include <ctime>
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDBasin.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();
	
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
  string CH_name; 
	string fill_ext = "_fill";
	file_info_in >> DEM_name >> CH_name;
	float Minimum_Slope;
	int DrainageArea;
	file_info_in >> Minimum_Slope >> DrainageArea;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_extension = "bil";
	string DEM_asc_extension = "asc";
	string csv_extension = "csv";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_extension);
	cout << "Got the DEM" << endl;
	
	// load the sources
	//LSDRaster sources_raster((path_name+sources_name), DEM_extension);

	// Set the no flux boundary conditions
  vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	// load the filled DEM
	LSDRaster filled_topo_test;
	LSDRaster load_DEM((path_name+DEM_name+fill_ext), DEM_extension);
  filled_topo_test = load_DEM;
	cout << "Got the filled DEM" << endl;
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	
	cout << "\t Loading Sources..." << endl;
	// load the sources
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+CH_name), DEM_extension, 1);
	cout << "Got the sources" << endl;
  
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "Got the channel network" << endl;
  //----------------------------------------------------------------------------------------------------//
  // GET ALL BASINS OF THE SPECIFIED STREAM ORDER
  //----------------------------------------------------------------------------------------------------//
  
  //cout << "Now getting all basins with a drainage area of: " << threshold_area << " m^2" << endl;
	
	//get the vector of basin junctions
	vector<int> BasinNodes = ChanNetwork.extract_basin_nodes_by_drainage_area(DrainageArea, FlowInfo);
	cout << "Got the basin nodes" << endl;
	vector<int> basin_junctions = ChanNetwork.extract_basin_junctions_from_nodes(BasinNodes, FlowInfo);
	cout << "Got the basin junctions" << endl;
	
	//vector<int> basin_junctions = ChanNetwork.extract_basins_order_outlet_junctions(BasinOrder, FlowInfo);
	
	// get raster of basins from the junctin vector
	LSDIndexRaster BasinRaster = ChanNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);
	string basin_ext = "_basins";
	BasinRaster.write_raster((path_name+DEM_name+basin_ext), DEM_extension);
	
	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}