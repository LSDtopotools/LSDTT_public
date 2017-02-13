//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// basin_puncher_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets all the basins of a specified drainage area and prints them
// to a raster.
// It reads in a list of sources from the OS MasterMap Water Network Layer which are
// used to generate the channel network.
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

	string DEM_name, csv_name, mask_name;
	float Minimum_Slope;
	int DrainageArea;
	string fill_ext = "_fill";
	string temp;
	file_info_in >> temp >> DEM_name
							 >> temp >> csv_name
							 >> temp >> mask_name
							 >> temp >> Minimum_Slope
							 >> temp >> DrainageArea;

	file_info_in.close();

	cout << "You are running the basin driver with the following settings:" << endl;
	cout << "\t DEM name: " << DEM_name << endl;
	cout << "\t Source CSV filename: " << csv_name << endl;
	cout << "\t Lochs filename: " << mask_name << endl;
	cout << "\t Minimum slope: " << Minimum_Slope << endl;
	cout <<" \t Drainage area for basin extraction: " << DrainageArea << endl;

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
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_extension);
	//LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
	//filled_topo_test.write_raster((path_name+DEM_name+fill_ext), DEM_extension);
	cout << "Got the filled DEM" << endl;

  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDRaster DrainageArea_raster = FlowInfo.write_DrainageArea_to_LSDRaster();
	string da_ext = "_DA";
	DrainageArea_raster.write_raster((path_name+DEM_name+da_ext), DEM_extension);

	cout << "\t Loading Sources..." << endl;
	// load the sources
  vector<int> sources = FlowInfo.Ingest_Channel_Heads_OS(path_name+csv_name);
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "\t Got channel network!" << endl;

	//write stream order array to a raster
  LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  string SO_name = "_SO";
	SOArray.write_raster((path_name+DEM_name+SO_name), DEM_extension);

	//write junction network to raster
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
	string JI_name = "_JI";
	JIArray.write_raster((path_name+DEM_name+JI_name), DEM_extension);

	//load the Lochs
	LSDRaster MaskRaster((path_name+mask_name), DEM_extension);
	cout << "\t Got the mask raster!" << endl;

  //----------------------------------------------------------------------------------------------------//
  // GET ALL BASINS OF THE SPECIFIED STREAM ORDER
  //----------------------------------------------------------------------------------------------------//

  cout << "Now getting all basins with a drainage area of: " << DrainageArea << " m^2" << endl;

	vector<int> basin_nodes = ChanNetwork.extract_basin_nodes_above_drainage_area_threshold(FlowInfo, DrainageArea);
	cout << "Got the basin nodes!" << endl;

	cout << "Removing lochs..." << endl;
	vector<int> new_basin_nodes = ChanNetwork.modify_basin_nodes_from_mask(basin_nodes, FlowInfo, MaskRaster);

	//checking position of new basin nodes
	string csv_out = "tay_basin_nodes_check";
	ofstream output_file;
	output_file.open((path_name+csv_out+"."+csv_extension).c_str());
	
	string old_basins = "basins_old";
	FlowInfo.print_vector_of_nodeindices_to_csv_file(basin_nodes, path_name+DEM_name+old_basins);

	string new_basins = "basins_new";
	FlowInfo.print_vector_of_nodeindices_to_csv_file(new_basin_nodes, path_name+DEM_name+new_basins);

	// cout << "Number of basins = " << basin_nodes.size() << endl;
	//
	// vector<int> basin_junctions = ChanNetwork.extract_basin_junctions_from_nodes(basin_nodes, FlowInfo);
	// cout << "Got the basin junctions" << endl;
	//
	// // get raster of basins from the junction vector
	// LSDIndexRaster BasinRaster = ChanNetwork.extract_basins_from_junction_vector_nested(basin_junctions, FlowInfo);
	// string basin_ext = "_basins";
	// BasinRaster.write_raster((path_name+DEM_name+basin_ext), DEM_extension);

	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
