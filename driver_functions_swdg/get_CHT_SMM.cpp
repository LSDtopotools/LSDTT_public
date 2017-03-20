//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get_CHT.cpp
//
// This program calculates channel heads using an Optimal wiener filter (used by Pelletier,
// 2013) and a quantile-quantile curvature threshold similar to Geonet.  It then uses a
// connected components threshold to create a channel network from the curvature mask.
//
// Following the channel extraction, the driver computes hilltop curvature,
// extracts basins from a list of lat/long coordinates and calculates basin average
// CHT statistics.
//
// References: Pelletier, J.D. (2013) A robust, two-parameter method for the extraction of
// drainage networks from high-resolution digital elevation models (DEMs): Evaluation using
// synthetic and real-world DEMs, Water Resources Research 49(1): 75-89, doi:10.1029/2012WR012452
//
// Passalacqua, P., Do Trung, T., Foufoula-Georgiou, E., Sapiro, G., & Dietrich, W. E.
// (2010). A geometric framework for channel network extraction from lidar: Nonlinear diffusion and
// geodesic paths. Journal of Geophysical Research: Earth Surface (2003�2012), 115(F1).
//
// He, L., Chao, Y., & Suzuki, K. (2008). A run-based two-scan labeling algorithm. Image Processing,
// IEEE Transactions on, 17(5), 749-756.
//
// Developed by:
//  David T. Milodowski
//  Simon M. Mudd
//  Stuart W.D. Grieve
//  Fiona J. Clubb
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
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Stuart W.D. Grieve, University of Edinburgh
// Fiona J. Clubb, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDIndexRaster.hpp"
#include "../TNT/tnt.h"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDBasin.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDParameterParser.hpp"
int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the CHT averaging tool!                  ||" << endl;
    cout << "|| This program has a number of options to compute     ||" << endl;
    cout << "|| basin averaged statistics                           ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Simon M. Mudd and Stuart W.D. Grieve                ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./basin_averaging_tool.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_BasinAvg.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << " http://lsdtopotools.github.io/LSDTT_book/#_chi_analysis_part_3_getting_chi_gradients_for_the_entire_landscape" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);
  
  // for the basin tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // set default float parameters
  int_default_map["basin_order"] = 2;
  int_default_map["threshold_contributing_pixels"] = 1500;
  int_default_map["search_radius_nodes"] = 10;
  int_default_map["connected_components_threshold"] = 100;
  int_default_map["threshold_stream_order"] = 1;
  int_default_map["search_radius_nodes"] = 10;
  

  // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;
  float_default_map["surface_fitting_radius"] = 6;
  float_default_map["roughness_threshold"] = 0.015;
  float_default_map["roughness_radius"] = 5;
  float_default_map["pruning_drainage_area"] = 1000;
  float_default_map["hilltop_slope_threshold"] = 0.4;
  
  
  // set default methods
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["print_basin_raster"] = false;
  bool_default_map["print_stream_order_raster"] = false;

  bool_default_map["check_basin_locations"] = false;
  
  bool_default_map["spawn_basins_from_outlets"] = false;

  bool_default_map["print_channels_to_csv"] = false;  
  bool_default_map["print_junction_csv"] = false;
  bool_default_map["print_sources_to_csv"] = true;
  
  // set default string method
  string_default_map["q_q_filename_prefix"] = "NULL";
  string_default_map["basin_outlet_csv"] = "NULL";
  string_default_map["sample_ID_column_name"] = "sample";
  
  // set default string method
  string_default_map["CHeads_file"] = "NULL";
  
  
  // Use the parameter parser to get the maps of the parameters required for the 
  // analysis

  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();
  
  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  
  cout << "Read filename is:" <<  DATA_DIR+DEM_ID << endl;
  
    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);  

  //Load the elevation data, fill it and generate a FlowInfo object
  LSDRaster DEM(DATA_DIR+DEM_ID, raster_ext);
  LSDRaster FilledDEM = DEM.fill(this_float_map["min_slope_for_fill"]);
  cout << "Filling the DEM. This may take a while..." << endl;
  LSDFlowInfo FlowInfo(boundary_conditions,FilledDEM);

  vector<int> FinalSources;
  // test if you have a channel heads file. If not, extract channels using our in hourse Wiener method.
  if(this_string_map["CHeads_file"] == "NULL")
  {
    cout << "You didn't give me a channel heads file, I am computing the channel heads" << endl;
    // NOTE: this is the Wiener method
    
    //Start of channel extraction
    string QQ_fname = OUT_DIR+OUT_ID+"__qq.txt";
    LSDRasterSpectral raster(DATA_DIR+DEM_ID, raster_ext);
    LSDIndexRaster connected_components = raster.IsolateChannelsWienerQQ(this_float_map["pruning_drainage_area"], 
                            this_float_map["surface_fitting_radius"], QQ_fname);
    
    cout << "Filtering connected components." << endl;
    LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(this_int_map["connected_components_threshold"]);
    LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();
    
    cout << "Making drainage skeleton." << endl;
    LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();
    LSDIndexRaster Ends = skeleton_raster.find_end_points();
    
    cout << "Processing endpoints." << endl;
    Ends.remove_downstream_endpoints(CC_raster, raster);

    //this processes the end points to only keep the upper extent of the channel network
    vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);

    // we need a temp junction network to search for single pixel channels
    LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
    LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

    //Now we have the final channel heads, so we can generate a channel network from them
    FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);
  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    FinalSources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t Got sources!" << endl;
  }
  
  LSDJunctionNetwork JunctionNetwork(FinalSources, FlowInfo);

  //Finally we write the channel heads to file so they can be used in other drivers.
  if (this_bool_map["print_sources_to_csv"])
  {
    cout << "Let me print the sources for you. You can use this later to skip channel extraction." << endl;
    string this_sources_file = OUT_DIR+OUT_ID+"_CH";
    FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(FinalSources,this_sources_file);
  }
  
  vector< int > basin_junctions = JunctionNetwork.ExtractBasinJunctionOrder(this_int_map["basin_order"], FlowInfo);
  LSDIndexRaster Basin_Raster = JunctionNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);

  //surface fitting - this is the second call of this method, this time on the unfiltered topography
  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(0); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(0); //plan curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  cout << "Calculating polyfit." << endl;
  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(this_float_map["surface_fitting_radius"], raster_selection);

  //raster_selection.clear();
  vector<int> raster_selection2;

  raster_selection2.push_back(0);
  raster_selection2.push_back(0);
  raster_selection2.push_back(1);

  cout << "Calculating roughness." << endl;
  vector<LSDRaster> Roughness = FilledDEM.calculate_polyfit_roughness_metrics(this_float_map["surface_fitting_radius"], 
                           this_float_map["roughness_radius"], raster_selection2);

  // extract hilltops
  cout << "Getting the ridgetops." << endl;
  LSDRaster hilltops = JunctionNetwork.ExtractRidges(FlowInfo);

  //get hilltop curvature using filter to remove positive curvatures
  cout << "Mapping the curvature onto the hilltops." << endl;
  LSDRaster cht_raster = FilledDEM.get_hilltop_curvature(Surfaces[3], hilltops);
  LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(cht_raster);

  //perform a gradient filter on the CHT data
 cout << "Hilltops with a gradient filter." << endl;
 LSDRaster CHT_gradient = JunctionNetwork.ExtractHilltops(CHT, Surfaces[1], this_float_map["hilltop_slope_threshold"]);

  //Set up a converter object to transform Lat Long into UTM
  LSDCoordinateConverterLLandUTM Converter;

  // Now get the lat long data
  string basin_outlet_fname = DATA_DIR+this_string_map["basin_outlet_csv"];
  LSDSpatialCSVReader BasinOutlets(RI,basin_outlet_fname);

  // set up some vectors
  vector<float> fUTM_northing;
  vector<float> fUTM_easting;
  vector<string> IDs;

  string column_name = "sample";
  BasinOutlets.get_data_in_raster_for_snapping(this_string_map["sample_ID_column_name"],fUTM_easting,fUTM_northing,IDs);

  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices

  JunctionNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing,
            this_int_map["search_radius_nodes"], this_int_map["threshold_stream_order"], FlowInfo,
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);

  int n_valid_points = int(valid_cosmo_points.size());  //The bumber of points which were within the current DEM
  cout << "I found " << n_valid_points << " in this raster. Here they are:" << endl;
  for (int i = 0; i< n_valid_points; i++)
  {
    int valid_point_index = valid_cosmo_points[i];
    cout << "Point " << i << endl;
    cout << "easting: " << fUTM_easting[valid_point_index] << " northing: " <<  fUTM_northing[valid_point_index] << endl;
    cout << "Snapped junction index: " << snapped_junction_indices[i] << endl << endl;
  }


  // print some stuff. We need to do this before the basin loop since it seems to be crashing there.
  if(this_bool_map["print_junction_csv"])
  {
    cout << "I am printing a junction csv" << endl;
    string JN_fname = OUT_DIR+OUT_ID+"_Junctions.csv";
    vector<int> JunctionList;
    JunctionNetwork.print_junctions_to_csv(FlowInfo, JunctionList,JN_fname);
  }

  if (this_bool_map["print_basin_raster"])
  {
    cout << "I am printing a basin raster." << endl;
    string BR_fname = OUT_DIR+OUT_ID+"_AllBasins.bil";
    LSDIndexRaster BasinRaster = JunctionNetwork.extract_basins_from_junction_vector_nested(basin_junctions, FlowInfo);
  }

  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "I am printing the channel netowk to csv." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_W_CN";
    JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
  }


  //Open a file to write the basin average data to
  ofstream WriteData;
  stringstream ss2;
  ss2 << OUT_DIR+OUT_ID << "_CHT_Data.csv";
  WriteData.open(ss2.str().c_str());

  //write headers
  WriteData << "_ID,min,max,median,mean,range,std_dev,std_err,count,min_gradient,max_gradient,median_gradient,mean_gradient,range_gradient,std_dev_gradient,std_err_gradient,count_gradient,min_internal,max_internal,median_internal,mean_internal,range_internal,std_dev_internal,std_err_internal,count_internal,min_internal_gradient,max_internal_gradient,median_internal_gradient,mean_internal_gradient,range_internal_gradient,std_dev_internal_gradient,std_err_internal_gradient,count_internal_gradient,bedrock_percentage,bedrock_percentage_internal" << endl;

  for(int samp = 0; samp<n_valid_points; samp++)
  {
    cout << "Cout getting the snapped basin: " << snapped_junction_indices[samp] << endl;
    LSDBasin thisBasin(snapped_junction_indices[samp], FlowInfo, JunctionNetwork);

    cout << "Writing CHT rasters." << endl;
    LSDRaster CHT_basin = thisBasin.write_raster_data_to_LSDRaster(CHT, FlowInfo);
    LSDRaster CHT_internal = thisBasin.keep_only_internal_hilltop_curvature(CHT_basin, FlowInfo);
    LSDRaster CHT_basin_gradient = thisBasin.write_raster_data_to_LSDRaster(CHT_gradient, FlowInfo);
    LSDRaster CHT_internal_gradient = thisBasin.keep_only_internal_hilltop_curvature(CHT_basin_gradient, FlowInfo);

    cout << "Calculating bedrock outcrop." << endl;
    float bedrock_full = FilledDEM.get_percentage_bedrock_ridgetops(Roughness[2], CHT_basin, this_float_map["roughness_threshold"]);
    float bedrock_internal = FilledDEM.get_percentage_bedrock_ridgetops(Roughness[2], CHT_internal, this_float_map["roughness_threshold"]);

    // put the 4 permutations of CHT rasters into a vector to loop over
    vector<LSDRaster> CHTs(4, CHT);
    CHTs[1] = CHT_gradient;
    CHTs[2] = CHT_internal;
    CHTs[3] = CHT_internal_gradient;

    cout << "Writing data, the sample is: "  << IDs[valid_cosmo_points[samp]] << endl;
    WriteData << IDs[valid_cosmo_points[samp]];

    for (int a=0; a < 4; ++a)
    {

      float mean = thisBasin.CalculateBasinMean(FlowInfo, CHTs[a]);
      float max = thisBasin.CalculateBasinMax(FlowInfo, CHTs[a]);
      float min = thisBasin.CalculateBasinMin(FlowInfo, CHTs[a]);
      float median = thisBasin.CalculateBasinMedian(FlowInfo, CHTs[a]);
      float stdd = thisBasin.CalculateBasinStdDev(FlowInfo, CHTs[a]);
      float stde = thisBasin.CalculateBasinStdError(FlowInfo, CHTs[a]);
      float range = thisBasin.CalculateBasinRange(FlowInfo, CHTs[a]);
      int count = thisBasin.CalculateNumDataPoints(FlowInfo, CHTs[a]);

      WriteData << "," << min << "," << max << "," << median << "," << mean << "," << range << "," << stdd << "," << stde << "," << count;
    }
    WriteData << "," << bedrock_full << "," << bedrock_internal << endl;

  }

  WriteData.close();

  //Write CHT to a raster and to a csv file
  CHT.write_raster((OUT_DIR+OUT_ID+"_CHT"), raster_ext);
  stringstream ss;
  ss << OUT_DIR+OUT_ID << "_Spatial_CHT.csv";
  bool is_North;
  int UTM_zone;
  FilledDEM.get_UTM_information(UTM_zone, is_North);
  int eId = 22;
  FilledDEM.HilltopsToCSV(CHT, CHT_gradient, Surfaces[1], UTM_zone, is_North, eId, ss.str());



}
