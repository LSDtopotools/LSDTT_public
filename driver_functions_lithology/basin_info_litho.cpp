//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// basin_info_litho
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options to extract basin information
// function to a lithology/geology map. This last can be pre-processed
// from a shapefile following the (future) instructions on LSDMappingTool
// LSD_Geology.py script available on:
// https://github.com/LSDtopotools/LSDMappingTools
// *I am working on this now, this is not ready yet - Boris*
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2016 Simon M. Mudd 2016
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
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDBasin.hpp"
#include "../LSDChiTools.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the chi mapping tool!                    ||" << endl;
    cout << "|| This program has a number of options to make chi    ||" << endl;
    cout << "|| plots and to map out slopes in chi space.           ||" << endl;
    cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./chi_mapping_tool.exe /LSDTopoTools/Topographic_projects/LSDTT_chi_examples/ Xian_example1.driver" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << " https://lsdtopotools.github.io/LSDTopoTools_ChiMudd2014/" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // Basic DEM preprocessing
  float_default_map["minimum_elevation"] = 0.0;
  float_default_map["maximum_elevation"] = 30000;
  float_default_map["min_slope_for_fill"] = 0.0001; // Set up a minimum slope when fillin' the dem to force a flow direction when nature,noise or flat areas create holes
  bool_default_map["raster_is_filled"] = false; // assume base raster is already filled
  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  string_default_map["CHeads_file"] = "NULL"; // channel head location extracted from another methods
  bool_default_map["print_raster_without_seas"] = false;


  // Selecting basins
  int_default_map["threshold_contributing_pixels"] = 1000; // number of pixels in the basins
  int_default_map["minimum_basin_size_pixels"] = 5000;
  int_default_map["maximum_basin_size_pixels"] = 500000;
  bool_default_map["test_drainage_boundaries"] = true;
  bool_default_map["only_take_largest_basin"] = false;
  string_default_map["BaselevelJunctions_file"] = "NULL";
  bool_default_map["extend_channel_to_node_before_receiver_junction"] = true;

  // IMPORTANT: S-A analysis and chi analysis wont work if you have a truncated
  // basin. For this reason the default is to test for edge effects
  bool_default_map["find_complete_basins_in_window"] = true;
  bool_default_map["find_largest_complete_basins"] = false;
  bool_default_map["print_basin_raster"] = false;

  // printing of rasters and data before chi analysis
  bool_default_map["convert_csv_to_geojson"] = false;  // This converts all cv files to geojson (for easier loading in a GIS)

  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_channels_to_csv"] = false;
  bool_default_map["print_junction_index_raster"] = false;
  bool_default_map["print_junctions_to_csv"] = false;
  bool_default_map["print_fill_raster"] = false;
  bool_default_map["print_DrainageArea_raster"] = false;
  bool_default_map["write_hillshade"] = false;
  bool_default_map["print_basic_M_chi_map_to_csv"] = false;

  // This burns a raster value to any csv output of chi data
  // Useful for appending geology data to chi profiles
  bool_default_map["geolithomap_to_csv"] = true;
  string_default_map["geolithomap_prefix"] = "NULL";
  string_default_map["burn_data_csv_column_header"] = "geolitho";

  // these print various basin and source data for visualisation
  bool_default_map["print_source_keys"] = false;
  bool_default_map["print_sources_to_csv"] = false;
  bool_default_map["print_sources_to_raster"] = false;
  bool_default_map["print_baselevel_keys"] = false;
  bool_default_map["find_complete_basins_in_window"] = true;
  bool_default_map["print_basin_raster"] = true;
  bool_default_map["write_hillshade"] = true;
  bool_default_map["print_chi_data_maps"] = true;

  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  cout << "PRINT THE PARAMETERS..." << endl;
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  string BaselevelJunctions_file = LSDPP.get_BaselevelJunctions_file();



  cout << "Read filename is: " <<  DATA_DIR+DEM_ID << endl;
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;


  // -------------------------------------------------Baslevel stuffs----------------------------------------------

  if(BaselevelJunctions_file == "NULL" || BaselevelJunctions_file == "Null" || BaselevelJunctions_file == "null" || BaselevelJunctions_file.empty() == true)
  {
    if(BaselevelJunctions_file.empty())
    {
      cout << "You have a null baselevel junctions file; the string is empty." << endl;
    }
    else
    {
      cout << "You have a null baselevel junctions file, it is: " << BaselevelJunctions_file << endl;
    }
  }
  else
  {
    BaselevelJunctions_file = RemoveControlCharactersFromEndOfString(BaselevelJunctions_file);
    BaselevelJunctions_file = DATA_DIR+BaselevelJunctions_file;
    cout << "You have selected a baselevel junctions file, it is: " << BaselevelJunctions_file << endl;
    cout << "Let me check if it exists..." << endl;

    ifstream test_file;
    test_file.open(BaselevelJunctions_file.c_str());
    if( test_file.fail() )
    {
      cout << "\nWHOOPS the baselevel file: \"" << BaselevelJunctions_file
         << "\" doesn't exist" << endl;
      cout << "I am changing it to a NULL value!" << endl;
      BaselevelJunctions_file = "NULL";
    }
    test_file.close();
  }
 // ----------------------------------------------end of Baslevel stuffs----------------------------------------------

 // ===================== LOADING RASTERS ====================================
    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);

  //============================================================================
  // check to see if the raster for burning exists - lithologic/geologic map
  LSDIndexRaster geolithomap;
  bool geolithomap_exists = false;
  string geolithomap_header = DATA_DIR+this_string_map["geolithomap_prefix"]+".hdr";
  if (this_bool_map["geolithomap_to_csv"])
  {
    cout << "Your lithologic map is: " << endl;
    cout <<  geolithomap_header << endl;
  }
  ifstream burn_head_in;
  burn_head_in.open(geolithomap_header.c_str());
  string burn_fname = DATA_DIR+this_string_map["geolithomap_prefix"];
  if( not burn_head_in.fail() )
  {
    geolithomap_exists = true;

    cout << "The lithologic raster exists. It has a prefix of: " << endl;
    cout <<  burn_fname << endl;
    LSDIndexRaster TempRaster(burn_fname,raster_ext);

    geolithomap = TempRaster;
  }
  else
  {
    cout << "No lithology raster. Please check the prefix is correctly spelled and without the extention" << endl;
    cout << "The file you tried to give me is: " << burn_fname << endl << "But it does not exists" << endl ;
    exit(EXIT_FAILURE);
  }

  // load the  DEM
  LSDRaster topography_raster;
  if (this_bool_map["remove_seas"])
  {
    cout << "I am removing high and low values to get rid of things that should be nodata." << endl;
    LSDRaster start_raster((DATA_DIR+DEM_ID), raster_ext);
    // now get rid of the low and high values
    float lower_threshold = this_float_map["minimum_elevation"];
    float upper_threshold = this_float_map["maximum_elevation"];
    bool belowthresholdisnodata = true;
    LSDRaster Flooded = start_raster.mask_to_nodata_using_threshold(lower_threshold,belowthresholdisnodata);
    belowthresholdisnodata = false;
    topography_raster = Flooded.mask_to_nodata_using_threshold(upper_threshold,belowthresholdisnodata);

    if (this_bool_map["print_raster_without_seas"])
    {
      cout << "I'm replacing your raster with a raster without seas." << endl;
      string this_raster_name = OUT_DIR+OUT_ID;
      topography_raster.write_raster(this_raster_name,raster_ext);
    }
  }
  else
  {
    LSDRaster start_raster((DATA_DIR+DEM_ID), raster_ext);
    topography_raster = start_raster;
  }
  cout << "Got the base dem: " <<  DATA_DIR+DEM_ID << endl;



  cout << "The lithologic raster the base raster and the lithologic raster are loaded " << endl;

  //================= END OF LOADING ============================================================================================================


  //================= Flow - River - Basin extracion from the topographic raster ============================================================================================================

  cout << "I will now divide the topographic raster into basins " << endl;

  // initialise variables to be assigned from .driver file
  // These will all be assigned default values
  cout << "Loading the parameters:" << endl;
  int threshold_contributing_pixels = this_int_map["threshold_contributing_pixels"];
  int minimum_basin_size_pixels = this_int_map["minimum_basin_size_pixels"];
  int maximum_basin_size_pixels = this_int_map["maximum_basin_size_pixels"];
  cout << "threshold_contributing_pixels: " << threshold_contributing_pixels << " - Define the density of the river network" << endl;
  cout << "I'll only select basins containing more than " << minimum_basin_size_pixels << " pixels and less than " << maximum_basin_size_pixels <<" pixels" << endl;
  float reso_dem = topography_raster.get_DataResolution();
  float minibas = minimum_basin_size_pixels * reso_dem* reso_dem;
  float maxibas = maximum_basin_size_pixels * reso_dem* reso_dem;
  cout << "(Between " << minibas <<" square meters and " << maxibas <<" square meters )" << endl;

  // Filling topography to get the right flow direction
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

  if (this_bool_map["print_fill_raster"])
  {
    cout << "Let me print the fill raster for you."  << endl;
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }

  if (this_bool_map["write_hillshade"])
  {
    cout << "Let me print the hillshade for you. " << endl;
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+OUT_ID+"_hs";
    hs_raster.write_raster(hs_fname,raster_ext);
  }


  cout << "\t Flow routing..." << endl;
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

  // calculate the flow accumulation
  cout << "\t Calculating flow accumulation (in pixels)..." << endl;
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  cout << "\t Converting to flow area..." << endl;
  LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

  if (this_bool_map["print_DrainageArea_raster"])
  {
    string DA_raster_name = OUT_DIR+OUT_ID+"_DArea";
    DrainageArea.write_raster(DA_raster_name,raster_ext);
  }

  // calcualte the distance from outlet
  cout << "\t Calculating flow distance..." << endl;
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  cout << "\t Loading Sources..." << endl;
  cout << "\t Source file is... " << CHeads_file << endl;

  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, threshold_contributing_pixels);

    cout << "The number of sources is: " << sources.size() << endl;
  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t Got sources!" << endl;
  }

  // now get the junction network
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

  // Print channels and junctions if you want them.
  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "I am going to print the channel network." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
    JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);

  }

  // Print sources
  if( this_bool_map["print_sources_to_csv"])
  {
    string sources_csv_name = OUT_DIR+OUT_ID+"_ATsources.csv";

    //write channel_heads to a csv file
    FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(sources, sources_csv_name);
    string sources_csv_name_2 = OUT_DIR+OUT_ID+"_ATsources_rowcol.csv";
    FlowInfo.print_vector_of_nodeindices_to_csv_file(sources, sources_csv_name_2);

  }

  // need to get base-level nodes , otherwise these catchments will be missed!
  vector< int > BaseLevelJunctions;
  vector< int > BaseLevelJunctions_Initial;

  //Check to see if a list of junctions for extraction exists
  if (BaselevelJunctions_file == "NULL" || BaselevelJunctions_file == "Null" || BaselevelJunctions_file == "null" || BaselevelJunctions_file.empty() == true)
  {
    cout << "To reiterate, there is no base level junction file. I am going to select basins for you using an algorithm. " << endl;
    // remove basins drainage from edge if that is what the user wants
    if (this_bool_map["find_complete_basins_in_window"])
    {
      cout << "I am going to look for basins in a pixel window that are not influended by nodata." << endl;
      cout << "I am also going to remove any nested basins." << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(FlowInfo, filled_topography, FlowAcc,
                                              this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
    }
    else
    {
      //Get baselevel junction nodes from the whole network
      BaseLevelJunctions_Initial = JunctionNetwork.get_BaseLevelJunctions();

      // now prune these by drainage area
      cout << "Removing basins with fewer than " << minimum_basin_size_pixels << " pixels" << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions_Initial,
                                              FlowInfo, FlowAcc, minimum_basin_size_pixels);
      cout << "Now I have " << BaseLevelJunctions.size() << " baselelvel junctions left. " << endl;

      if (this_bool_map["find_largest_complete_basins"])
      {
        cout << "I am looking for the largest basin not influenced by nodata within all baselevel nodes." << endl;
        BaseLevelJunctions = JunctionNetwork.Prune_To_Largest_Complete_Basins(BaseLevelJunctions,FlowInfo, filled_topography, FlowAcc);
      }
      else
      {
        if (this_bool_map["test_drainage_boundaries"])     // now check for edge effects
        {
          cout << endl << endl << "I am going to remove basins draining to the edge." << endl;
          BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions,FlowInfo, filled_topography);
          //BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge(BaseLevelJunctions,FlowInfo);
        }
      }
    }
  }
  else
  {
    cout << "I am attempting to read base level junctions from a base level junction list." << endl;
    cout << "If this is not a simple text file that only contains itegers there will be problems!" << endl;

    //specify junctions to work on from a list file
    //string JunctionsFile = DATA_DIR+BaselevelJunctions_file;
    cout << "The junctions file is: " << BaselevelJunctions_file << endl;


    vector<int> JunctionsList;
    ifstream infile(BaselevelJunctions_file.c_str());
    if (infile)
    {
      cout << "Junctions File " << BaselevelJunctions_file << " exists" << endl;;
      int n;
      while (infile >> n) BaseLevelJunctions_Initial.push_back(n);
    }
    else
    {
      cout << "Fatal Error: Junctions File " << BaselevelJunctions_file << " does not exist" << endl;
      exit(EXIT_FAILURE);
    }

    // Now make sure none of the basins drain to the edge
    cout << "I am pruning junctions that are influenced by the edge of the DEM!" << endl;
    cout << "This is necessary to have complete basins." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions_Initial, FlowInfo, filled_topography);

  }

  // Correct number of base level junctions
  int N_BaseLevelJuncs = BaseLevelJunctions.size();
  cout << "The number of basins I will analyse is: " << N_BaseLevelJuncs << endl;
  if (N_BaseLevelJuncs == 0)
  {
    cout << "I am stopping here since I don't have any basins to analyse." << endl;
    exit(EXIT_FAILURE);
  }

  // Now getting the basin list - I need a chitool object ??? I think there might be a way to optimise this
  LSDChiTools blasinf(FlowInfo);
  LSDIndexRaster basin_raster = blasinf.get_basin_raster(FlowInfo, JunctionNetwork, BaseLevelJunctions);
  geolithomap.detect_unique_values();
  geolithomap.detect_unique_values();
  basin_raster.write_raster(OUT_DIR+OUT_ID+"_basicheck","bil"); // saving for checking


}
