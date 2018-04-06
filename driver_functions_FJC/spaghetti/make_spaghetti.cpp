//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// make spaghetti plots for river profile clustering.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2018 Fiona Clubb
//
// Developer can be contacted by clubb _at_ uni-potsdam.de
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "../../LSDParameterParser.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDRasterInfo.hpp"
#include "../../LSDStatsTools.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDBasin.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDSpatialCSVReader.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the spaghetti making tool, developed by  ||" << endl;
    cout << "|| Fiona Clubb at the University of Potsdam            ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./make_spaghetti.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_spaghetti.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
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

  // parameters for DEM preprocessing
  float_default_map["min_slope_for_fill"] = 0.0001;
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["write hillshade"] = false;

  // parameters for channel extraction
  int_default_map["threshold_contributing_pixels"] = 2000;
  string_default_map["channel heads fname"] = "NULL";
  bool_default_map["print_channels_to_csv"] = false;

  // parameters for basin selection
  int_default_map["basin_order"] = 3;
  bool_default_map["print_basin_raster"] = false;
  bool_default_map["print_junction_csv"] = false;

  // parameters for getting the profile plots
  bool_default_map["print_spaghetti_profiles_to_csv"] = true;

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

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;


  //============================================================================
  // Start gathering necessary rasters
  //============================================================================
  LSDRaster filled_topography;
  // now get the flow info object
  if ( this_bool_map["raster_is_filled"] )
  {
    cout << "You have chosen to use a filled raster." << endl;
    filled_topography = topography_raster;
  }
  else
  {
    cout << "Let me fill that raster for you, the min slope is: "
         << this_float_map["min_slope_for_fill"] << endl;
    filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
  }


  if (this_bool_map["print_filled_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }

  // check to see if you need hillshade
  if (this_bool_map["write hillshade"])
  {
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+OUT_ID+"_hs";
    hs_raster.write_raster(hs_fname,raster_ext);
  }

  //============================================================================
  // Get the flow info
  //============================================================================
  cout << "\t Flow routing..." << endl;
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

  // now deal with the channel network
  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    // calculate the flow accumulation
    cout << "\t Calculating flow accumulation (in pixels)..." << endl;
    LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    // Now get the sources from flow accumulation
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

    cout << "The number of sources is: " << sources.size() << endl;

  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t Got sources!" << endl;
  }

  // Now create the network
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "I am writing the channels to csv." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_W_CN";
    JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_W_CN.geojson";
      LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_W_CN.csv");
      thiscsv.print_data_to_geojson(gjson_name);
    }

  }

  if( this_bool_map["print_junctions_to_csv"])
  {
    cout << "I am writing the junctions to csv." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_JN.csv";
    JunctionNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_JN.geojson";
      LSDSpatialCSVReader thiscsv(channel_csv_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

  // Now get all the basins of the selected order
  vector<int> BasinJunctions = JunctionNetwork.extract_basins_order_outlet_junctions(this_int_map["basin_order"], FlowInfo);

  // get the distance from outlet
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  if( this_bool_map["print_spaghetti_profiles_to_csv"])
  {
    string csv_filename = OUT_DIR+OUT_ID+"_spaghetti_profiles.csv";
    JunctionNetwork.write_river_profiles_to_csv(BasinJunctions, FlowInfo, DistanceFromOutlet, filled_topography, csv_filename);
  }

  if( this_bool_map["print_basin_raster"])
  {
    LSDIndexRaster BasinRaster = JunctionNetwork.extract_basins_from_junction_vector(BasinJunctions, FlowInfo);
    string basins_name = OUT_DIR+OUT_ID+"_basins";
    BasinRaster.write_raster(basins_name, raster_ext);
  }


}
