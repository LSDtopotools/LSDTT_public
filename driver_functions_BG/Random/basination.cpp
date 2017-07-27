//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Basin Averager
//
// A program for looking at different aspects of drainage basins.
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate
// different kinds of chi analysis
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2017 Simon M. Mudd 2017
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDRasterInfo.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDBasin.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDSpatialCSVReader.hpp"
#include "../../LSDParameterParser.hpp"

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Declare the print basins function
// Why don't we put this in junction network?
// Because it uses the LSDBasin object which has a load of dependencies
// so if we put it there we would need to change all the make files.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void print_basins(string basin_raster_name, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Junctions);


int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the basin averaging tool!                ||" << endl;
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
  int_default_map["basin_spawn_padding_pixels"] = 100;
  int_default_map["start_basin_order"] = 1;
  int_default_map["end_basin_order"] = 7;


  // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;

  // set default methods
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["print_basin_raster"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["get_basins_from_outlets"] = false;
  bool_default_map["print_junctions_to_csv"] = false;
  bool_default_map["print_channels_to_csv"] = false;

  // switches for basin spawning
  bool_default_map["spawn_basins_from_outlets"] = false;
  bool_default_map["spawn_csv_file_from_basin_spawn"] = false;
  bool_default_map["spawn_parameter_files_from_basin_spawn"] = false;
  bool_default_map["remove_edge_basins"] = false;

  // simple plotting options
  bool_default_map["write hillshade"] = false;
  bool_default_map["convert_csv_to_geojson"] = false;

  // set default string method
  string_default_map["slope method"] = "polynomial";
  string_default_map["averaging_raster_vector"] = "NULL";
  string_default_map["basin_outlet_csv"] = "NULL";
  string_default_map["sample_ID_column_name"] = "IDs";
  string_default_map["parameter_file_for_spawning"] = "NULL";

  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  string avg_raster_names = "averaging_raster_vector";
  vector<string> vec_of_rasters = LSDPP.parse_string_vector(avg_raster_names);

  cout << endl <<endl << "List of rasters to be averaged" << endl;
  for (int i = 0; i<int(vec_of_rasters.size()); i++)
  {
    cout << "Raster is: " << vec_of_rasters[i] << endl;
  }
  cout << endl;

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

  // debugging
  //string basin_outlet_fname = DATA_DIR+this_string_map["basin_outlet_csv"];
  //cout << "The basin file is: " << basin_outlet_fname << endl;
  //LSDSpatialCSVReader Outlet_CSV_data(RI,basin_outlet_fname);
  //exit(EXIT_FAILURE);

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;
  LSDRaster FillRaster;
  // now get the flow info object
  if ( this_bool_map["raster_is_filled"] )
  {
    cout << "You have chosen to use a filled raster." << endl;
    FillRaster = topography_raster;
  }
  else
  {
    cout << "Let me fill that raster for you, the min slope is: "
         << this_float_map["min_slope_for_fill"] << endl;
    FillRaster = topography_raster.fill(this_float_map["min_slope_for_fill"]);
  }

  // write the fill raster if necessary
  if ( this_bool_map["print_filled_raster"] )
  {
    string fill_ext = "_FILL";
    FillRaster.write_raster((OUT_DIR+OUT_ID+fill_ext), raster_ext);
  }


  // now load the flow info object
  cout << "\t Flow routing..." << endl;
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,FillRaster);

  // now deal with the channel network
  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file != "null")
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

  vector<int> basin_junctions;
  int padding_pixels = this_int_map["basin_spawn_padding_pixels"];
  int sampi = this_int_map["start_basin_order"]; // Order of basins
  int sampf = this_int_map["end_basin_order"]; // Order of basins


  for(int sampt = sampi; sampt<=sampf;sampt ++){
    cout << "I am getting the basins from all basins of order " << sampt << endl;
    if(this_bool_map["remove_edge_basins"])
    {
      basin_junctions = JunctionNetwork.ExtractBasinJunctionOrder(sampt, FlowInfo);
    }
    else
    {
      basin_junctions = JunctionNetwork.ExtractBasinJunctionOrderKeepEdgeBasins(sampt, FlowInfo);
    }
    cout << "The number of basin junctions are: " << basin_junctions.size() << endl;

    cout << "Getting basin" << endl;
    string BRExt = "_Spawned_"+itoa(sampt);;
    // This a good thing, not to be confused with Brexit, which is the most shit thing ever.

    string basin_raster_name = OUT_DIR+OUT_ID+BRExt;
    print_basins(basin_raster_name, FlowInfo, JunctionNetwork,
                   basin_junctions);
  }



}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints the basins and an additional file that has basin centroids
// and labelling information for plotting
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void print_basins(string basin_raster_name, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Junctions)
{
  int N_Juncs = Junctions.size();
  LSDCoordinateConverterLLandUTM Converter;

  // Get some data members for holding basins and the raster
  vector<LSDBasin> AllTheBasins;
  map<int,int> drainage_of_other_basins;
  LSDIndexRaster BasinMasterRaster;

  //cout << "I am trying to print basins, found " << N_BaseLevelJuncs << " base levels." << endl;
  // Loop through the junctions
  for(int BN = 0; BN<N_Juncs; BN++)
  {
    //cout << "Getting basin " << BN << " and the junction is: "  << BaseLevelJunctions[BN] << endl;
    LSDBasin thisBasin(Junctions[BN],FlowInfo, JunctionNetwork);
    //cout << "...got it!" << endl;
    AllTheBasins.push_back(thisBasin);

    // This is required if the basins are nested--test the code which numbers
    // to be overwritten by a smaller basin
    drainage_of_other_basins[Junctions[BN]] = thisBasin.get_NumberOfCells();
  }

  // now loop through everything again getting the raster
  if (N_Juncs > 0)     // this gets the first raster
  {
    BasinMasterRaster = AllTheBasins[0].write_integer_data_to_LSDIndexRaster(Junctions[0], FlowInfo);
  }

  // now add on the subsequent basins
  for(int BN = 1; BN<N_Juncs; BN++)
  {
    AllTheBasins[BN].add_basin_to_LSDIndexRaster(BasinMasterRaster, FlowInfo,
                              drainage_of_other_basins, Junctions[BN]);
  }


  // We need to use bil format since there is georeferencing
  string raster_ext = "bil";
  // print the basin raster
  BasinMasterRaster.write_raster(basin_raster_name, raster_ext);
  cout << "Finished writing the basins!" << endl;

}
