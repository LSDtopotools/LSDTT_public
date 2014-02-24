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
	file_info_in >> DEM_name >> sources_name;
	float Minimum_Slope;
	file_info_in >> Minimum_Slope;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
	
	// load the sources
	LSDRaster sources_raster((path_name+sources_name), DEM_flt_extension);

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
	
	//string NI_name = "_NI";
  //LSDIndexRaster NodeIndex = FlowInfo.write_NodeIndex_to_LSDIndexRaster();
	//NodeIndex.write_raster((path_name+DEM_name+NI_name), DEM_flt_extension);
	
	int NRows = topo_test.get_NRows();
  int NCols = topo_test.get_NCols();
  //float XMinimum = topo_test.get_XMinimum();
  //float YMinimum = topo_test.get_YMinimum();
  //float DataResolution = topo_test.get_DataResolution();
  float NoDataValue = topo_test.get_NoDataValue();
  
	//get the sources from raster to vector
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+sources_name),DEM_flt_extension);
  
  //----------------------------------------------------------------------------------------------------//
  // GET JUNCTION FOR BASIN EXTRACTION
  //----------------------------------------------------------------------------------------------------//
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
	string JI_name = "_JI";
	JIArray.write_raster((path_name+DEM_name+JI_name),DEM_flt_extension);
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();	
	
	//----------------------------------------------------------------------------------------------------//
	
	//----------------------------------------------------------------------------------------------------//
  // EXTRACT BASIN FROM JUNCTION TO GET BASIN RASTER
  //----------------------------------------------------------------------------------------------------//
  
  int junction_number = 129;
  LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
  LSDRaster BasinRaster = Basin.write_raster_data_to_LSDRaster(topo_test, FlowInfo);
  
  //write basin to LSDRaster
  string basin_name = "_basin";
  string jn_name = itoa(junction_number);
	string uscore = "_";
	jn_name = uscore+jn_name;
	BasinRaster.write_raster((path_name+DEM_name+basin_name+jn_name), DEM_flt_extension);      
	  	
}