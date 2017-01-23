//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// channel_extraction_dreich.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program calculates channel heads using the Dreich method, Clubb et al. (2014)
//
// Reference: Clubb, F. J., S. M. Mudd, D. T. Milodowski, M. D. Hurst,
// and L. J. Slater (2014), Objective extraction of channel heads from
// high-resolution topographic data, Water Resour. Res., 50, doi: 10.1002/2013WR015167.
//
// Developed by:
//  Fiona Clubb
//  Simon M. Mudd
//  David T. Milodowski
//  Stuart W.D. Grieve
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
//
// Fiona J. Clubb, Univertsity of Edinburgh
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Stuart W.D. Grieve, University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDRasterSpectral.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../TNT/tnt.h"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
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
	

	
	// read inputs from file
	string DEM_name, output_path;
  float A_0, m_over_n;
  string temp;
	int NJunctions;
  file_info_in >> temp >> DEM_name
							 >> temp >> output_path
	             >> temp >> A_0
	             >> temp >> m_over_n
							 >> temp >> NJunctions;
  file_info_in.close();
	
	// get some file names
	string fill_ext = "_fill";
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_extension = "bil";
	
	cout << "Running the DrEICH channel extraction, parameters are: \n"
			 << "Raster name: " << DEM_name << "\n"
		   << "A_0: " << A_0 << "\n"
		 	 << "m/n: " << m_over_n << "\n";

	// Set the no flux boundary conditions
  vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
  // get the DEM
  LSDRaster filled_topo_test(DEM_f_name, DEM_extension);
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
	
	cout << "\t Got the DEM, ingesting channel heads..." << endl;

  // Load in the sources
	vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_name+"_CH_DrEICH"), DEM_extension);
	
	// get the flow distance
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	
	//get the junction network
	LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);
	
	cout << "\t Writing chi profiles to csv..." << endl;
	
	// Write the chi profiles
  JunctionNetwork.write_valley_hilltop_chi_profiles_to_csv(sources, A_0, m_over_n, FlowInfo, DistanceFromOutlet, filled_topo_test, NJunctions, output_path, DEM_name);			
	cout << "\t Done!" << endl;
}
