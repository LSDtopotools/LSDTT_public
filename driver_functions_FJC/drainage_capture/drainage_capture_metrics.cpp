//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// drainage_capture_metrics.cpp
//
// This driver is used to determine the topographic expression of drainage
// captures.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2018 Fiona J. Clubb
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
    cout << "|| Welcome to the drainage capture tool, developed by    ||" << endl;
    cout << "|| Fiona Clubb at the University of Potsdam              ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./drainage_capture_metrics.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_capture.param" << endl;
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

  // set default float parameters
  int_default_map["basin_order"] = 2;
  int_default_map["threshold_contributing_pixels"] = 2000;
  int_default_map["search_radius_nodes"] = 10;

  // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;

  // set default methods
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["print_junctions_to_csv"] = false;
  bool_default_map["print_junction_angles_to_csv"] = false;
  bool_default_map["read_junctions_from_file"] = false;
  bool_default_map["print_junction_statistics_to_csv"] = false;
  bool_default_map["convert_csv_to_geojson"] = false;
  bool_default_map["get_basin_perimeter"] = false;


  bool_default_map["print_basin_raster"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_stream_order_csv"] = false;
  bool_default_map["get_basins_from_outlets"] = false;
  bool_default_map["print_junction_csv"] = false;
  bool_default_map["print_channels_to_csv"] = false;
  bool_default_map["spawn_basins_from_outlets"] = false;
  bool_default_map["spawn_csv_file_from_basin_spawn"] = false;
  bool_default_map["spawn_parameter_files_from_basin_spawn"] = false;
  bool_default_map["write hillshade"] = false;

  // set default string method
  string_default_map["slope method"] = "polynomial";
  string_default_map["averaging_raster_vector"] = "NULL";
  string_default_map["basin_outlet_csv"] = "NULL";
  string_default_map["sample_ID_column_name"] = "IDs";
  string_default_map["parameter_file_for_spawning"] = "NULL";

  // turn on parameter



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

  if (this_bool_map["print_junctions_to_csv"])
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

  // now get the basin perimeters - this is for getting the hypsometry
  // of the perimeter to look for drainage captures
  if (this_bool_map["get_basin_perimeter"])
  {
    cout << "I am getting the basin perimeters" << endl;
    string perimeter_name = OUT_DIR+OUT_ID+"_Perimeters.csv";

    // for now, just get one of the perimeters.
    int JunctionNumber = 0;

    // get the node index of this junction
    int basin_node = JunctionNetwork.get_Node_of_Junction(JunctionNumber);

    // now get the perimeter
    vector<int> perimeter_vec = FlowInfo.basin_edge_extractor(basin_node, topography_raster);
    FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(perimeter_vec, perimeter_name);

    LSDBasin ABasin(JunctionNumber, FlowInfo, JunctionNetwork);
    LSDRaster ThisBasin = ABasin.write_raster_data_to_LSDRaster(filled_topography, FlowInfo);
    ThisBasin.write_raster((OUT_DIR+OUT_ID+"_Perimeters"),"bil");
    ABasin.print_perimeter_hypsometry_to_csv(FlowInfo, perimeter_name, filled_topography);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_Perimeter.geojson";
      LSDSpatialCSVReader thiscsv(perimeter_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

  // get the channel profiles for each basin to look at the hypsometries
  if ( this_bool_map["get_channel_profiles"])
  {

  }

  // get all junction angles to test for captures.
  if( this_bool_map["print_junction_angles_to_csv"])
  {
    cout << "I am testing the junction angle code." << endl;
    string JAngles_csv_name = OUT_DIR+OUT_ID+"_JAngles.csv";
    vector<int> JunctionList;
    JunctionNetwork.print_junction_angles_to_csv(JunctionList, FlowInfo, JAngles_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_JAngles.geojson";
      LSDSpatialCSVReader thiscsv(JAngles_csv_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

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


/*
  // do a few tests
  // first make two vectors
  float rad1 = 0.25*M_PI;
  float rad2 = (0.2)*M_PI;

  float gradient1 = tan(rad1);
  float gradient2 = tan(rad2);

  cout << "Gradient 1 is: " << gradient1 << " and gradient 2 is: " << gradient2 << endl;

  float dx = 0.1;

  vector<float> x1;
  vector<float> y1;
  vector<float> x2;
  vector<float> y2;

  x1.push_back(1.3);
  y1.push_back(2.6);
  x2.push_back(-0.2);
  y2.push_back(0.4);

  for (int i = 0; i<11; i++)
  {
    x1.push_back( x1[i] +dx);
    y1.push_back( y1[i] + dx*gradient1);

    x2.push_back( x2[i] +dx);
    y2.push_back( y2[i] + dx*gradient2);
  }

  // now try the gradient finder
  float R2;
  float i1, g1, i2, g2;
  vector<float> means1 = orthogonal_linear_regression( x1, y1, i1, g1, R2);
  vector<float> means2 = orthogonal_linear_regression( x2, y2, i2, g2, R2);

  cout << "Regressed orthogonal gradients are: " << g1 << " " << g2 << endl;
  cout << "Regressed orthogonal intercepts are: " << i1 << " " << i2 << endl;

  // now normal regression
  least_squares_linear_regression(x1, y1, i1, g1,R2);
  least_squares_linear_regression(x2, y2, i2, g2, R2);

  cout << "Regressed simple gradients are: " << g1 << " " << g2 << endl;
  cout << "Regressed simple intercepts are: " << i1 << " " << i2 << endl;


  // now get angle between two datasets
  bool channels_point_downstream = false;
  float angle = angle_between_two_vector_datasets(x1, y1,x2, y2,channels_point_downstream);
  cout << "Angle is: " << angle << " radians, which is " << deg(angle) << " degrees." << endl;
*/

}
