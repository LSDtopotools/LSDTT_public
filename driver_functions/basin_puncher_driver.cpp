//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// basin_puncher_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program writes an LSDRaster of a basin of interest from a larger DEM.
// It uses the sources from the DrEICH algorithm to create the  channel network.
// User needs to specify the junction number of the basin required (Line 156)
//  
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
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
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
  string CH_name; 
	string fill_ext = "_Fill";
	string filt_ext = "_filtered";
	file_info_in >> DEM_name >> CH_name;
	float Minimum_Slope;
	int junction_number;
	file_info_in >> Minimum_Slope >> junction_number;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "bil";
	string DEM_asc_extension = "asc";
	string csv_extension = "csv";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
	cout << "Got the DEM" << endl;
	
	// load the sources
	//LSDRaster sources_raster((path_name+sources_name), DEM_flt_extension);

	// Set the no flux boundary conditions
  vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
		// get the filled file
	//cout << "Filling the DEM" << endl;
  //LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
		// load the filled DEM
	LSDRaster filled_topo_test;
	LSDRaster load_DEM((path_name+DEM_name+fill_ext), DEM_flt_extension);
  filled_topo_test = load_DEM;
	cout << "Got the filled DEM" << endl;
	//int NRows = filled_topo_test.get_NRows();
	//cout << "NRows: " << NRows << endl;
	//filled_topo_test.write_raster((DEM_f_name),DEM_flt_extension);
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	
	cout << "\t Loading Sources..." << endl;
	// load the sources
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+CH_name), DEM_flt_extension, 1);
	cout << "Got the sources" << endl;
  
  //----------------------------------------------------------------------------------------------------//
  // GET JUNCTION FOR BASIN EXTRACTION
  //----------------------------------------------------------------------------------------------------//
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
//	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
//	string JI_name = "_JI";	JIArray.write_raster((path_name+DEM_name+JI_name),DEM_flt_extension);
//	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();	
//	string SO_name = "_SO_wiener";
//	SOArray.write_raster((path_name+DEM_name+SO_name),DEM_flt_extension);
//	cout << "Got the junction network" << endl;
	
	//----------------------------------------------------------------------------------------------------//
	
	//----------------------------------------------------------------------------------------------------//
  // EXTRACT BASIN FROM JUNCTION TO GET BASIN RASTER
  //----------------------------------------------------------------------------------------------------//
  
  cout << "Basin junction: " << junction_number << endl;
	LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
	cout << "Generating the basin object" << endl;
  LSDRaster BasinRaster = Basin.write_raster_data_to_LSDRaster(filled_topo_test, FlowInfo);
	cout << "Got the basin raster" << endl;
		
	// Trim raster to get rid of nodata values around the edge
	LSDRaster TrimmedRaster = BasinRaster.RasterTrimmer();
	cout << "Got the trimmed raster" << endl;
  
  //write basin to LSDRaster
  string basin_name = "_basin";
  string jn_name = itoa(junction_number);
	string uscore = "_";
	jn_name = uscore+jn_name;
	TrimmedRaster.write_raster((path_name+DEM_name+basin_name+jn_name), DEM_flt_extension);      
	  	
}