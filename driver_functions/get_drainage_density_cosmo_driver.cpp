//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get_drainage_density_cosmo_driver.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets the drainage density of basins from the X and Y coordinates
// of cosmogenic radionuclide sample points to compare drainage density with erosion
// rate. It also outputs files to plot a boomerang for each basin with a cosmogenic
// erosion rate.  It gets drainage density based on a channel network methodology 
// described by Clubb et al. (in review)
// It takes an input file with the format: 
// X_coordinate Y_coordinate Cosmo_Erosion_Rate Cosmo_Point_Error
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
  //float XMinimum = topo_test.get_XMinimum();
  //float YMinimum = topo_test.get_YMinimum();
  //float DataResolution = topo_test.get_DataResolution();
  float NoDataValue = topo_test.get_NoDataValue();
  
	//get the sources from raster to vector
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+sources_name),DEM_flt_extension);
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	
	// get all the junctions of the second order basins
	int BasinOrder = 2;
	vector<int> second_order_junctions = ChanNetwork.extract_basins_order_outlet_junctions(BasinOrder, FlowInfo);
	
  
  //-------------------------------------//
  // DRAINAGE DENSITY
  //-------------------------------------//
  
  cout << "Reading in the cosmo points" << endl;
  // reading in cosmo points - getting upstream junctions and then extracting basins
  string string_filename;
	string filename = "fr1m_cosmo";
	string dot = ".";
	string extension = "txt";
  string_filename = filename+dot+extension;
  ifstream coords_list;
  coords_list.open(string_filename.c_str());
  
  float X_coord, Y_coord, ErosionRate, Error;
  vector<float> X_coords;
  vector<float> Y_coords;
  vector<float> ErosionRates;
  vector<float> Errors;
  while (coords_list >> X_coord >> Y_coord >> ErosionRate >> Error)
  {
     cout << "X Coord: " << X_coord << " Y Coord: " << Y_coord << " Erosion Rate: " << ErosionRate << "Error: " << Error << endl;
     X_coords.push_back(X_coord);
     Y_coords.push_back(Y_coord);
     ErosionRates.push_back(ErosionRate);
     Errors.push_back(Error);
  }
  
 
  cout << "Snapping to channel; getting upstream junctions" << endl;
  int threshold_SO = 1;
  int search_radius = 5;
  vector<int> junction_vector;
  cout << "Line 180" << endl;
  for (unsigned int i = 0; i < X_coords.size(); i++)
  {
    cout << "Line 183" << endl;
    cout << flush << "Snapped = " << i+1 << " of " << X_coords.size() << "\r";
    cout << "Line 185" << endl;
    int NICosmo = ChanNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(X_coords[i], Y_coords[i],
                                                                search_radius, threshold_SO, FlowInfo);
    cout << "Line 187" << endl;
    int JunctionCosmo = ChanNetwork.find_upstream_junction_from_channel_nodeindex(NICosmo, FlowInfo);
    cout << "Line 189" << endl;
    cout << "Junction of cosmo point: " << JunctionCosmo << endl;
    junction_vector.push_back(JunctionCosmo);  
  }                                                     
    
  cout << "Extracting the basins" << endl;
  //getting basins to calculate drainage density	
  LSDIndexRaster Basins = ChanNetwork.extract_basins_from_junction_vector(junction_vector, FlowInfo);
  
  // get flow directions - D8 flowdir for drainage density calculations 
  Array2D<int> FlowDir_array = FlowInfo.get_FlowDirection();
  
  cout << "Calculating drainage density" << endl;
  //get drainage density
  //LSDRaster drainage_density = filled_topo_test.DrainageDensity(SOArray, Basins, FlowDir_array);          
  //string DD_name = "_DD";
  //drainage_density.write_raster((path_name+DEM_name+DD_name),DEM_flt_extension);   
  
  //getting the slope raster
  float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
  vector<LSDRaster> surface_fitting;
  string slope_name = "_slope";
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;                      // this indicates you want the slope
  
  surface_fitting = filled_topo_test.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);

  LSDRaster Slope = surface_fitting[1];  
  Slope.write_raster((path_name+DEM_name+slope_name), DEM_flt_extension);
            
  int no_junctions = junction_vector.size();  
  cout << "Number of basins: " << no_junctions << endl;
  string string_filename2;
  string filename2 = "drainage_density_cosmo_fr1m";
  string_filename2 = filename2+dot+extension;
  ofstream DD_cosmo;
  DD_cosmo.open(string_filename2.c_str());
  
  cout << "Created the drainage density file" << endl;
  
  cout << "Calculating DINF area for boomerang plotting" << endl;
  
  // D-infinty flowdirection   
  Array2D<float> FlowDir = filled_topo_test.D_inf_FlowDir();
  LSDRaster DinfArea = filled_topo_test.D_inf_FlowArea(FlowDir);
  
  //declaring variables for boomerang plotting
  float log_bin_width = 0.1;
	int SplineResolution = 20000;
	float bin_threshold =  0;
	
  //get mean DD of each basin for plotting
  for (int i = 0; i < no_junctions; i++)
	{
    cout << flush << "Junction = " << i+1 << " of " << no_junctions << "\r";
    int junction_number = junction_vector[i];
    float BasinErosionRate = ErosionRates[i];
    float BasinError = Errors[i];
    
    // set basin parameters
    LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
    Basin.set_FlowLength(SOArray, FlowInfo);
    Basin.set_DrainageDensity();
    Basin.set_SlopeMean(FlowInfo, Slope);
    //Basin.Plot_Boomerang(Slope, DinfArea, FlowInfo, log_bin_width, SplineResolution, bin_threshold, path_name);
    
    // return basin parameters
    float drainage_density = Basin.get_DrainageDensity();
    cout << "Drainage density: " << drainage_density << endl;
    float basin_slope = Basin.get_SlopeMean();
    float slope_stdev = Basin.CalculateBasinStdDev(FlowInfo, Slope);
    float slope_sterr = Basin.CalculateBasinStdError(FlowInfo, Slope);

    if (drainage_density != NoDataValue)
    {
      DD_cosmo << drainage_density << " " << basin_slope << " " << slope_stdev << " " << slope_sterr << " " << BasinErosionRate << " " << BasinError << endl;
    } 
  }        
}
