//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get_boomerang_plots_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets the drainage density of all the second order basins within the DEM.  It
// also creates three files for each basin which can be used to create a boomerang plot for each,
// using the python visualisation script slope_area_plot_new.py.  The file "junctions.txt" is created
// to allow automation of the visualisation for every basin.
// Drainage density is calculated using the DrEICH algorithm for channel head prediction.  For
// details on its methodology see Clubb et al. (in review).
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
#include "../LSDBasin.hpp"
#include "../LSDChiNetwork.hpp"

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
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
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
  float XMinimum = topo_test.get_XMinimum();
  float YMinimum = topo_test.get_YMinimum();
  float DataResolution = topo_test.get_DataResolution();
  float NoDataValue = topo_test.get_NoDataValue();
  
	//get the sources from raster to vector
	Array2D<float> sources_array = sources_raster.get_RasterData();
	vector<int> sources;
	for (int row = 0; row < NRows; row++)
	{
    for (int col = 0; col < NCols; col++)
    {
      if (sources_array[row][col] != NoDataValue)
      {
        int node = FlowInfo.retrieve_node_from_row_and_column(row,col);
        sources.push_back(node);
      }
    }
  }

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	
	// get the slope and the dinf area for boomerang plots
		// Calculate polyfit coefficients and the tangential curvature
  Array2D<float> a;
  Array2D<float> b;
  Array2D<float> c;
  Array2D<float> d;
  Array2D<float> e;
  Array2D<float> f;
  float window_radius = 6;
  filled_topo_test.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster curv = filled_topo_test.calculate_polyfit_curvature(a,b);
  LSDRaster Slope = filled_topo_test.calculate_polyfit_slope(d, e);
	
	cout << "Getting the junctions of the basins" << endl;
  // get all the junctions of the second order basins
	int BasinOrder = 2;
	vector<int> second_order_junctions = ChanNetwork.extract_basins_order_outlet_junctions(BasinOrder, FlowInfo);
	int no_junctions = second_order_junctions.size();
	
	string string_filename;
	string filename = "junctions";
	string dot = ".";
	string extension = "txt";
	string_filename = filename+dot+extension;
	ofstream junction_list;
	junction_list.open(string_filename.c_str());
	
	for (int i = 0; i < no_junctions; i++)
  {
    junction_list << second_order_junctions[i] << endl;
  }	
	junction_list.close();
	
	// get the drainage density of the basins
	LSDIndexRaster Basins = ChanNetwork.extract_basins_from_junction_vector(second_order_junctions, FlowInfo);
  
  // get flow directions - D8 flowdir for drainage density calculations 
  Array2D<int> FlowDir_array = FlowInfo.get_FlowDirection();
  
  cout << "Calculating drainage density" << endl;
  //get drainage density
  LSDRaster drainage_density = filled_topo_test.DrainageDensity(SOArray, Basins, FlowDir_array);          
  string DD_name = "_DD";
  drainage_density.write_raster((path_name+DEM_name+DD_name),DEM_flt_extension); 
  
  cout << "Calculating mean CHT of basins" << endl;
  // D-infinty flowdirection for HFR   
  Array2D<float> FlowDir = filled_topo_test.D_inf_FlowDir();
  LSDRaster DinfArea = filled_topo_test.D_inf_FlowArea(FlowDir);
  
  //hilltop flow routing
  float critical_slope = 0.4;
  LSDRaster Ridges = ChanNetwork.ExtractRidges(FlowInfo);
  LSDRaster hilltops = ChanNetwork.ExtractHilltops(Ridges, Slope, critical_slope);
  //vector< Array2D<float> > Routed_Hilltop_Data = filled_topo_test.HFR_Driver(hilltops, FlowDir, SOArray, Basins, DEM_name);  
  //write routed_hilltop_data to a series of LSDRasters
  //LSDRaster RoutedHilltops = filled_topo_test.LSDRasterTemplate(Routed_Hilltop_Data[0]);
  
  //calculate CHT and mean slope of basins
  LSDRaster CHT_temp = filled_topo_test.get_hilltop_curvature(curv, hilltops);
  
  //removing pixels where CHT is positive (noise)
  Array2D<float> CHT_array(NRows,NCols,NoDataValue);
  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      float curvature = CHT_temp.get_data_element(row,col);
      if (curvature < 0)
      {
        CHT_array[row][col] = curvature; 
      }
    }
  }
  LSDRaster CHT(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, CHT_array);
  LSDRaster mean_CHT = CHT.BasinAverager(Basins);
  LSDRaster mean_slope = Slope.BasinAverager(Basins);
  
  //write some rasters
  string CHT_name = "_CHT";
  mean_CHT.write_raster((path_name+DEM_name+CHT_name),DEM_flt_extension);
  
  string slope_name = "_slope";
  mean_slope.write_raster((path_name+DEM_name+slope_name),DEM_flt_extension);
  	
	/*cout << "Creating files for boomerang plotting" << endl;
  // get the boomerang plot from each basin
	float log_bin_width = 0.1;
	int SplineResolution = 10000;
	float bin_threshold =  0;
	
  for (int i = 0; i < no_junctions; i++)
	{
    cout << flush << "Junction = " << i+1 << " of " << no_junctions << "\r";
    int junction_number = second_order_junctions[i];
    LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
    Basin.Plot_Boomerang(Slope, DinfArea, FlowInfo, log_bin_width, SplineResolution, bin_threshold, path_name);
  }  */
}
