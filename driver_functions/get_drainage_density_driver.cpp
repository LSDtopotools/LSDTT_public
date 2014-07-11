//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get_drainage_density_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets the drainage density of all basins with a threshold area of
// 100,000 m^2 or above to compare drainage density with mean basin slope and mean hilltop
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
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_flt_extension);
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
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+sources_name),DEM_flt_extension);
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	
	// get all the basins greater than 100,000 m^2
	float AreaThreshold = 100000;
	vector<int> basin_nodes = ChanNetwork.extract_basin_nodes_by_drainage_area(AreaThreshold, FlowInfo);
	vector<int> basin_junctions = ChanNetwork.extract_basin_junctions_from_nodes(basin_nodes, FlowInfo);
	  
  cout << "Extracting the basins" << endl;
  //getting basins to calculate drainage density	
  LSDIndexRaster Basins = ChanNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);
  string basins_name = "_basins";
  Basins.write_raster((path_name+DEM_name+basins_name), DEM_flt_extension);
  
  // get flow directions - D8 flowdir for drainage density calculations 
  Array2D<int> FlowDir_array = FlowInfo.get_FlowDirection();
  
  //cout << "Calculating drainage density" << endl;
  //get drainage density
  //LSDRaster drainage_density = filled_topo_test.DrainageDensity(SOArray, Basins, FlowDir_array);          
  //string DD_name = "_DD";
  //drainage_density.write_raster((path_name+DEM_name+DD_name),DEM_flt_extension);   
  
  //get the slope and the curvature from polynomial fitting
  
  float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
  vector<LSDRaster> surface_fitting;
  string slope_name = "_slope";
  string curv_name = "_curv";
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;                      // this indicates you want the slope
  raster_selection[3] = 1;                      // this indicates you want the curvature
  
  surface_fitting = filled_topo_test.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);

  LSDRaster Slope = surface_fitting[1];
  LSDRaster Curvature = surface_fitting[3];
  
  Slope.write_raster((path_name+DEM_name+slope_name), DEM_flt_extension);
  Curvature.write_raster((path_name+DEM_name+curv_name), DEM_flt_extension);
     
  // Create the text files for writing info to
            
  int no_junctions = basin_junctions.size();  
  cout << "Number of basins: " << no_junctions << endl;
  
  string string_filename;
  string dot = ".";
  string extension = "txt";
  string filename = "_drainage_density_binned";
  string_filename = DEM_name+filename+dot+extension;
  ofstream DD_binned;
  DD_binned.open(string_filename.c_str());
    
  string string_filename2;
  string filename2 = "_drainage_density_cloud";
  string_filename2 = DEM_name+filename2+dot+extension;
  ofstream DD_cloud;
  DD_cloud.open(string_filename2.c_str());
  
  cout << "Created the drainage density file" << endl;
    	
  //get metrics for each basin for plotting
  
  vector<float> DrainageDensities;
  vector<float> Slopes;
  vector<float> CHTs;
  
  for (int i = 0; i < no_junctions; i++)
	{
    cout << flush << "Junction = " << i+1 << " of " << no_junctions << "\r";
    int junction_number = basin_junctions[i];
    int StreamOrder = ChanNetwork.get_StreamOrder_of_Junction(FlowInfo, junction_number);
    
    // Get the hilltop curvature
    int MinOrder = 1;
    int MaxOrder = StreamOrder-1;
    LSDRaster Hilltops = ChanNetwork.ExtractRidges(FlowInfo, MinOrder, MaxOrder);
    LSDRaster CHT_temp = filled_topo_test.get_hilltop_curvature(Curvature, Hilltops);
  
    // Remove hilltop pixels with positive curvature (noise)
  
    Array2D<float> CHT_array = CHT_temp.get_RasterData();
    Array2D<float> CHT_array_final(NRows, NCols, NoDataValue);
    for (int row = 0; row < NRows; row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        if (CHT_array[row][col] < 0)
        {
          CHT_array_final[row][col] = CHT_array[row][col];
        }
      }
    }
    LSDRaster CHT(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,CHT_array_final);

    // set basin parameters
    LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
    Basin.set_FlowLength(SOArray, FlowInfo);
    Basin.set_DrainageDensity();
    Basin.set_SlopeMean(FlowInfo, Slope);
    Basin.set_CHTMean(FlowInfo, CHT);
    //Basin.Plot_Boomerang(Slope, DinfArea, FlowInfo, log_bin_width, SplineResolution, bin_threshold, path_name);
    
    // return basin parameters
    float drainage_density = Basin.get_DrainageDensity();
    cout << "Drainage density: " << drainage_density << endl;
    float basin_slope = Basin.get_SlopeMean();
    float basin_CHT = Basin.get_CHTMean();
    float basin_area = Basin.get_Area();
    if (isnan(basin_CHT) == false)
    {
      Slopes.push_back(basin_slope);
      CHTs.push_back(abs(basin_CHT));
      DrainageDensities.push_back(drainage_density);
    }

    if (drainage_density != NoDataValue || isnan(basin_CHT) == false)
    {
      DD_cloud << drainage_density << " " << basin_CHT << " " << basin_slope << " " << basin_area << endl;
    } 
  }   
  
  float bin_width = 0.01;
  float bin_lower_limit = 0;
  float bin_threshold = 0.02;
  vector<float> CHTX_output;
  vector<float> DDY_output;
  vector<float> midpoints_output;
  vector<float> MedianY_output;
  vector<float> CHT_StDev_X_output;
  vector<float> DD_StDev_Y_output;
  vector<float> CHT_StandardErrorX_output;
  vector<float> DD_StandardErrorY_output;
  vector<int> number_observations_output;
  vector<float> SlopeY_output;
  vector<float> Slope_StDev_Y_output;
  vector<float> Slope_StandardErrorY_output;
  
  cout << "Binning CHT/drainage density" << endl;

  bin_data(CHTs, DrainageDensities, bin_width, CHTX_output, DDY_output, midpoints_output, MedianY_output, 
           CHT_StDev_X_output, DD_StDev_Y_output, CHT_StandardErrorX_output, DD_StandardErrorY_output,
           number_observations_output, bin_lower_limit, NoDataValue); 
           
  RemoveSmallBins(CHTX_output, DDY_output, midpoints_output, CHT_StDev_X_output, DD_StDev_Y_output, CHT_StandardErrorX_output, DD_StandardErrorY_output,
                  number_observations_output, bin_threshold);
           
  cout << "Binning CHT/slope" << endl; 
           
  bin_data(CHTs, Slopes, bin_width, CHTX_output, SlopeY_output, midpoints_output, MedianY_output, 
           CHT_StDev_X_output, Slope_StDev_Y_output, CHT_StandardErrorX_output, Slope_StandardErrorY_output,
           number_observations_output, bin_lower_limit, NoDataValue);        
           
  RemoveSmallBins(CHTX_output, SlopeY_output, midpoints_output, CHT_StDev_X_output, Slope_StDev_Y_output, CHT_StandardErrorX_output, Slope_StandardErrorY_output,
                  number_observations_output, bin_threshold); 
           
  for (int i = 0; i < CHTX_output.size(); i++)
  {
    DD_binned << CHTX_output[i] << " " << CHT_StDev_X_output[i] << " " << CHT_StandardErrorX_output[i] << " " << DDY_output[i] << " " <<
    DD_StDev_Y_output[i] << " " << DD_StandardErrorY_output[i] << " " << SlopeY_output[i] << " " << Slope_StDev_Y_output[i] << " " <<
    Slope_StandardErrorY_output[i] << endl;
  }
}
