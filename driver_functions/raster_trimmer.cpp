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
	string DEM_name = argv[2];

	cout << "The path is: " << path_name << " and the DEM name is: " << DEM_name << endl;

	string DEM_flt_extension = "bil";
	// load the DEM
	cout << "\t Loading the raster..." << endl;
	LSDRaster DEM((path_name+DEM_name), DEM_flt_extension);
	cout << "\t Done!" << endl;
	
	// Trim raster to get rid of nodata values around the edge
	cout << "\t Trimming the raster..." << endl;
	LSDRaster TrimmedRaster = DEM.RasterTrimmerSpiral();
	cout << "\t Done!" << endl;
  
  //write basin to LSDRaster
  string trimmed_name = "_trimmed";
	TrimmedRaster.write_raster((path_name+DEM_name+trimmed_name), DEM_flt_extension);      
	  	
}