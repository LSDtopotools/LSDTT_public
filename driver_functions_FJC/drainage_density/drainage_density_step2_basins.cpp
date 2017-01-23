//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// drainage_density_step2_basins.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets the drainage density of all the 3rd order basins in the 
// "DEM_name_junctions.txt" file to compare drainage density with mean hilltop
// curvature (used as a proxy for erosion rate). It gets drainage density based on a channel
// network methodology described by Clubb et al. (2014). The data are also binned.
//
// Outputs a txt file for the data cloud with the format
// drainage_density mean_hilltop_curvature drainage_area 
//
// Outputs a txt file for the binned data with the format
// CHT_value CHT_standard_deviation CHT_standard_error DD_value DD_standard_deviation 
// DD_standard_error
//
// Developed by:
//  Fiona J. Clubb
//  Simon M. Mudd
//  Stuart W. D. Grieve
//
// Copyright (C) 2016 Fiona Clubb and Simon M. Mudd 2016
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
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../LSDBasin.hpp"
#include "../../LSDCRNParameters.hpp"

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
  
  // get some file names
    string DEM_name; 
    string sources_name; 
    string fill_ext = "_fill";
	file_info_in >> DEM_name >> sources_name;
	float Minimum_Slope;
	file_info_in >> Minimum_Slope;  
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "bil";
  
    // Get list of junctions into a vector - this file is produced by ./drainage_density_step1_junctions.out
    string ListOfJunctions = DEM_name+"_DD_junctions_bedrock.txt";
    vector<int> junction_list;
    int this_junction;
    string line;
    ifstream junction_input(ListOfJunctions.c_str());
    if (junction_input.is_open())
    {
        while (junction_input >> this_junction)
        {
            junction_list.push_back(this_junction);
        }
        junction_input.close();
    }
 
    int no_junctions = int(junction_list.size());
    cout << "Number of basins: " << no_junctions << endl;
  
    // load in rasters
  
    // get the filled file
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_flt_extension);
  
    // Set the no flux boundary conditions
    vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	//get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
    float NoDataValue = filled_topo_test.get_NoDataValue();
	
	//get the sources from raster to vector
    vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+sources_name),DEM_flt_extension);
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	
	LSDIndexRaster Basin_Raster = ChanNetwork.extract_basins_from_junction_vector(junction_list, FlowInfo);
    string basin_name = "_basins";
    Basin_Raster.write_raster((path_name+DEM_name+basin_name), DEM_flt_extension);
		
    //get the CHT raster
    string CHT_ext = "_CHT";
	LSDRaster CHT((path_name+DEM_name+CHT_ext), DEM_flt_extension);
	
	// get the SO array
	string SO_ext = "_SO";
	LSDIndexRaster SOArray((path_name+DEM_name+SO_ext), DEM_flt_extension);
    
    // get the roughness raster
    string RE_ext = "_s3_1_5";
    LSDRaster Roughness((path_name+DEM_name+RE_ext), DEM_flt_extension);
   	
    //get metrics for each basin for plotting
  
    vector<float> DrainageDensities;
    vector<float> CHTs;
  
    string string_filename2;
    string dot = ".";
    string extension = "txt";
    string filename2 = "_drainage_density_cloud_plus_bedrock";
    string_filename2 = DEM_name+filename2+dot+extension;
    ofstream DD_cloud;
    DD_cloud.open(string_filename2.c_str());
 
    //set roughness threshold for identifying bedrock
    //float threshold = 0.015;
    
    for (int i = 0; i < int(junction_list.size())-1; i++)
	{
        cout << flush << "Junction = " << i+1 << " of " << junction_list.size() << "\r";
        int junction_number = junction_list[i];

        // set basin parameters
        LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
        Basin.set_FlowLength(SOArray, FlowInfo);
        Basin.set_DrainageDensity();
        Basin.set_Perimeter(FlowInfo);
        LSDRaster CHT_basin = Basin.write_raster_data_to_LSDRaster(CHT, FlowInfo);
        LSDRaster CHT_internal = Basin.keep_only_internal_hilltop_curvature(CHT_basin, FlowInfo);
        Basin.set_CHTMean(FlowInfo, CHT_internal);
        
        //get the percentage of bedrock ridgetops
        //LSDRaster basin_elev = Basin.write_raster_data_to_LSDRaster(filled_topo_test, FlowInfo);
        //LSDRaster BasinRoughness = Basin.write_raster_data_to_LSDRaster(Roughness, FlowInfo);
        //float percent_bedrock = basin_elev.get_percentage_bedrock_ridgetops(BasinRoughness, CHT_basin, threshold);
    
        // return basin parameters
        float drainage_density = Basin.get_DrainageDensity();
        float basin_CHT = Basin.get_CHTMean();
        //cout << "Drainage density: " << drainage_density << " Basin CHT: " << basin_CHT << " Percent bedrock: " << percent_bedrock << endl;
        float basin_area = Basin.get_Area(); 
        if (drainage_density != NoDataValue && isnan(basin_CHT) == false)
        {
            CHTs.push_back(abs(basin_CHT));
            DrainageDensities.push_back(drainage_density);
            DD_cloud << drainage_density << " " << basin_CHT << " " << basin_area << endl;
        }
    }   
  
    DD_cloud.close();
  
    cout << "Line 210" << endl;
    // Create the text files for writing info to  
    string string_filename;
    string filename = "_drainage_density_binned_plus_bedrock";
    string_filename = DEM_name+filename+dot+extension;
    ofstream DD_binned;
    DD_binned.open(string_filename.c_str());

  
    cout << "Created the drainage density file" << endl;
  
    float bin_width = 0.005;
    float bin_lower_limit = 0;
    float bin_threshold = 0.05;
    vector<float> CHTX_output;
    vector<float> DDY_output;
    vector<float> midpoints_output;
    vector<float> MedianY_output;
    vector<float> CHT_StDev_X_output;
    vector<float> DD_StDev_Y_output;
    vector<float> CHT_StandardErrorX_output;
    vector<float> DD_StandardErrorY_output;
    vector<int> number_observations_output;
  
    cout << "Binning CHT/drainage density" << endl;

    bin_data(CHTs, DrainageDensities, bin_width, CHTX_output, DDY_output, midpoints_output, MedianY_output, 
           CHT_StDev_X_output, DD_StDev_Y_output, CHT_StandardErrorX_output, DD_StandardErrorY_output,
           number_observations_output, bin_lower_limit, NoDataValue); 
           
    RemoveSmallBins(CHTX_output, DDY_output, midpoints_output, CHT_StDev_X_output, DD_StDev_Y_output,           CHT_StandardErrorX_output, DD_StandardErrorY_output, number_observations_output, bin_threshold);

    for (int i = 0; i < int(CHTX_output.size()); i++)
    {
        DD_binned << CHTX_output[i] << " " << CHT_StDev_X_output[i] << " " << CHT_StandardErrorX_output[i] << " " <<    DDY_output[i] << " " << DD_StDev_Y_output[i] << " " << DD_StandardErrorY_output[i] << endl;
    }
    DD_binned.close();
}
