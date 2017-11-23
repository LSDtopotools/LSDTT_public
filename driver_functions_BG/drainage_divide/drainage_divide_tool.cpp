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
// 
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

    // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;

  int_default_map["threshold_contributing_pixels"] = 2000;
  int_default_map["search_radius_nodes"] = 50;
  int_default_map["padding_pixels"] = 100;
  string_default_map["basin_outlet_csv"] = "NULL";
  bool_default_map["print_individual_basin_raster"] = false;
  bool_default_map["merge_adjacent_raster"] = false;

  // Analysis
  bool_default_map["swath_my_ridge"] = false;
  bool_default_map["extract_my_sources"] = false;




  // set default methods
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["print_basin_raster"] = false;

  bool_default_map["spawn_csv_file_from_basin_spawn"] = false;
  bool_default_map["remove_edge_basins"] = false;

  // simple plotting options
  bool_default_map["write hillshade"] = false;



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

  // debugging
  //string basin_outlet_fname = DATA_DIR+this_string_map["basin_outlet_csv"];
  //cout << "The basin file is: " << basin_outlet_fname << endl;
  //LSDSpatialCSVReader Outlet_CSV_data(RI,basin_outlet_fname);
  //exit(EXIT_FAILURE);

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;

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

 

  // Now we are going to get the basins in question. The junctions of these
  // basins will go into the vector basin_junctions
  // There are several options for this
  // 1) If this_bool_map["get_basins_from_outlets"] it gets basins from
  //     a file containing locations of basin outlets (it snaps to these
  //     outlets)
  // 2) If not true, it defaults to getting basins by basin order
  // 3) There should be a few more options here, one idea is to have a raster filling
  //     function that does some kind of heirarchical sorting to fill the raster.
  vector<int> basin_junctions;
  vector<string> IDs;

  // This uses a list lat-long coordinates and snaps these coordinates to the
  // neares channel.
  if(this_string_map["basin_outlet_csv"] != "NULL")
  {
    cout << "I am getting your basins from a list of nodes. " << endl;

    // first we see if the nodes file exists
    string basin_outlet_fname = DATA_DIR+this_string_map["basin_outlet_csv"];

    cout << "The basin file is: " << basin_outlet_fname << endl;
    LSDSpatialCSVReader Outlet_CSV_data(RI,basin_outlet_fname);

    cout << "I am reading the following samples: " << endl;
    vector<double> latitude = Outlet_CSV_data.get_latitude();
    vector<double> longitude = Outlet_CSV_data.get_longitude();
    int n_samples = int(latitude.size());
    for(int samp = 0; samp<n_samples; samp++)
    {
      cout << " lat: " << latitude[samp] << " long: " << longitude[samp] << endl;
    }

    // now get the junction network
    vector<float> fUTM_easting;
    vector<float> fUTM_northing;
    Outlet_CSV_data.get_x_and_y_from_latlong(fUTM_easting,fUTM_northing);

    int threshold_stream_order = 2;
    vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
    vector<int> snapped_node_indices;       // a vector to hold the valid node indices
    vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
    cout << "The search radius is: " << this_int_map["search_radius_nodes"] << endl;
    JunctionNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing,
                this_int_map["search_radius_nodes"], threshold_stream_order, FlowInfo,
                valid_cosmo_points, snapped_node_indices, snapped_junction_indices);

    cout << "The number of valid points is: " << int(valid_cosmo_points.size()) << endl;

    // Now get the basin rasters
    // HEY HEY| |LOOK HERE: I don't think we really need this since basins are
    // printed at a later stage with the flag "print_basin_raster"
    //string basin_raster_name = OUT_DIR+OUT_ID+"_AllBasins";
    //print_basins(basin_raster_name, FlowInfo, JunctionNetwork,
    //               snapped_junction_indices);
    basin_junctions = snapped_junction_indices;

    int n_valid_points =  int(basin_junctions.size());
    cout << "I am trying to spawn basins, found " << n_valid_points << " valid points." << endl;
    vector<LSDBasin> basin_list;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
      
      cout << "Getting basin" << endl;
      LSDBasin thisBasin(basin_junctions[samp],FlowInfo, JunctionNetwork);
      basin_list.push_back(thisBasin);

      if(this_bool_map["print_individual_basin_raster"])
      {
        cout << "Now getting the basin raster" << endl;
        LSDRaster BasinRaster = thisBasin.TrimPaddedRasterToBasin(this_int_map["padding_pixels"],
                                             FlowInfo,FillRaster);
        cout << "Formatting the path" << endl;

        //string DEM_newpath = ReformatPath(DEMpath);
        string DEMnewname = OUT_DIR+OUT_ID+"_Spawned_"+itoa(samp);
        cout << "Writing a new basin to: " << DEMnewname << endl;
        BasinRaster.write_raster(DEMnewname,"bil");
      }
    }


    vector<int> perimeter;
    if(this_bool_map["merge_adjacent_raster"])
    {
      int n_basins = basin_list.size();
      
      bool basin_check = false;
      // checking if the basins are adjacent
      if(n_basins>1)
      {
        cout << "I am now checking if your basins are adjacents" << endl;
        //bool all_adjacent = true;
        LSDBasin out_basin = basin_list[0];
        for(size_t samp = 0; samp<basin_list.size(); samp++)
        {
          for(size_t sec = 0; sec<basin_list.size(); sec++)
          {
            if(sec != samp && basin_list[samp].is_adjacent(basin_list[sec], FlowInfo))
            {
              basin_check = true;
            }
          }
          if(!basin_check)
          {
            cout<< "Your basin number "<< samp+1 << " is not adjacent with the others. I am aborting" << endl;
            exit(EXIT_FAILURE);
          }
          basin_check = false;
          basin_list[samp].print_perimeter_to_csv(FlowInfo, OUT_DIR+OUT_ID+"_"+itoa(samp)+"_perimeter.csv");
        }
        cout << "Your Basins are adjacent, let me merge their perimeter and remove theIR Adjacent border, It takes a while I need to optimize this part" << endl;
        perimeter = out_basin.merge_perimeter_nodes_adjacent_basins(basin_list, FlowInfo);
      }
      else
      {
        LSDBasin out_basin = basin_list[0];
        perimeter = out_basin.get_Perimeter_nodes();
      }
    }
    else
    {
      LSDBasin out_basin = basin_list[0];
      perimeter = out_basin.get_Perimeter_nodes();
    }

    if(this_bool_map["extract_my_sources"])
    {
      vector<int> selected_sources;
      selected_sources = basin_list[0].get_source_node_from_perimeter(perimeter, FlowInfo, JunctionNetwork, this_int_map["search_radius_nodes"]);
      cout<< perimeter.size() << endl;

      // printing the selected sources
      if(selected_sources.size() == 0){cout<<"No Channel heads detected, you can try to increase you search radius or check your basins, Check as well if your basin is entirely in the DEM" << endl; exit(EXIT_FAILURE);}

      string sources_csv_name = OUT_DIR+OUT_ID+"_ATsourcesDD.csv";

      //write selected channel_heads to a csv file
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(selected_sources, sources_csv_name);
      cout << "I printed the selected sources into " << sources_csv_name << endl;

    //write all channel_heads to a csv file
      sources_csv_name = OUT_DIR+OUT_ID+"_ATsources.csv";
     
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(sources, sources_csv_name);
      cout << "I printed the selected sources into " << sources_csv_name << endl;
    }


    if(this_bool_map["swath_my_ridge"])
    {
      cout<< "TO BUILD" << endl;

    }




  }
  else
  {
    cout << "You need to give me a csv file containing the latitude and longitude of the basin you want to ingest " <<  endl;
  }


  // This function is used to create new rasters that contain only the
  // basin, with a padding (set here to )
//   if (this_bool_map["spawn_basins_from_outlets"])
//   {
//     // now loop through the valid points, spawning a basin each time
//     int n_valid_points =  int(basin_junctions.size());
//     int padding_pixels = this_int_map["basin_spawn_padding_pixels"];

//     // check to see if there is a sample name vector, if not write one with junction indices
//     if (IDs.size() != basin_junctions.size())
//     {
//       vector<string> newID;
//       for (int i = 0; i<n_valid_points; i++)
//       {
//         newID.push_back(itoa(basin_junctions[i]));
//       }
//       IDs = newID;
//     }

//     cout << "I am trying to spawn basins, found " << n_valid_points << " valid points." << endl;
//     for(int samp = 0; samp<n_valid_points; samp++)
//     {

//       cout << "Getting basin" << endl;
//       LSDBasin thisBasin(basin_junctions[samp],FlowInfo, JunctionNetwork);

//       cout << "Now getting the basin raster" << endl;
//       LSDRaster BasinRaster = thisBasin.TrimPaddedRasterToBasin(padding_pixels,
//                                            FlowInfo,FillRaster);
//       cout << "Formatting the path" << endl;

//       //string DEM_newpath = ReformatPath(DEMpath);
//       string DEMnewname = OUT_DIR+OUT_ID+"_Spawned_"+itoa(samp);
//       cout << "Writing a new basin to: " << DEMnewname << endl;
//       BasinRaster.write_raster(DEMnewname,"bil");


//       if(this_bool_map["spawn_csv_file_from_basin_spawn"])
//       {
//         cout << "I am going to spawn csv files with the appropriate basins" << endl;
//         cout << "This sample is: " << IDs[samp] << endl;

//         // first we see if the nodes file exists
//         string basin_outlet_fname = DATA_DIR+this_string_map["basin_outlet_csv"];

//         string csv_spawn_fname = OUT_DIR+OUT_ID+"_CSVSpawned_"+itoa(samp)+".csv";
//         LSDSpatialCSVReader Outlet_CSV_data(RI,basin_outlet_fname);
//         vector<string> data_for_selection;
//         data_for_selection.push_back(IDs[samp]);
//         LSDSpatialCSVReader filtered = Outlet_CSV_data.select_data_to_new_csv_object(this_string_map["sample_ID_column_name"], data_for_selection);
//         filtered.print_data_to_csv(csv_spawn_fname);

//         if(this_bool_map["spawn_parameter_files_from_basin_spawn"])
//         {
//           cout << "Let me update and print a parameter file for you for the spawned basin." << endl;

//           cout << "I am going to parse this file: " << endl;
//           cout << this_string_map["parameter_file_for_spawning"] << endl;
//           // load a parameter parser
//           LSDParameterParser SpawnPP(DATA_DIR,this_string_map["parameter_file_for_spawning"]);

//           // force the parsing of everything in the file
//           SpawnPP.force_parse();

//           // update the file names and paths
//           string parameter_fname = OUT_ID+"_Spawned_"+itoa(samp)+".param";
//           map<string,string> replace_parameters;
//           string new_read_fname = OUT_ID+"_Spawned_"+itoa(samp);
//           string new_write_fname = OUT_ID+"_Spawned_"+itoa(samp);
//           replace_parameters["basin_outlet_csv"] = OUT_ID+"_CSVSpawned_"+itoa(samp)+".csv";

//           cout << "I am loading the file: " << endl << DATA_DIR+this_string_map["parameter_file_for_spawning"] << endl;
//           cout << "I am writing the file: " << endl << parameter_fname << endl;

//           SpawnPP.replace_and_print_parameter_file(parameter_fname,
//                                      DATA_DIR, new_read_fname, OUT_DIR, new_write_fname,
//                                      replace_parameters);
//         }
//       }
//     }
//   }

//   // Print basins to file if wanted
//   if(this_bool_map["print_basin_raster"])
//   {
//     string BRExt = "_AllBasins";
//     // This a good thing, not to be confused with Brexit, which is the most shit thing ever.

//     string basin_raster_name = OUT_DIR+OUT_ID+BRExt;
//     // NOTE: This creates nested basins so if a basin is a subbasin of another basin,
//     // the smaller basin will overwrite the larger basin
//     print_basins(basin_raster_name, FlowInfo, JunctionNetwork,
//                    basin_junctions);
//   }

// }
// //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




// //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// // This function prints the basins and an additional file that has basin centroids
// // and labelling information for plotting
// //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// void print_basins(string basin_raster_name, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
//                                vector<int> Junctions)
// {
//   int N_Juncs = Junctions.size();
//   LSDCoordinateConverterLLandUTM Converter;

//   // Get some data members for holding basins and the raster
//   vector<LSDBasin> AllTheBasins;
//   map<int,int> drainage_of_other_basins;
//   LSDIndexRaster BasinMasterRaster;

//   //cout << "I am trying to print basins, found " << N_BaseLevelJuncs << " base levels." << endl;
//   // Loop through the junctions
//   for(int BN = 0; BN<N_Juncs; BN++)
//   {
//     //cout << "Getting basin " << BN << " and the junction is: "  << BaseLevelJunctions[BN] << endl;
//     LSDBasin thisBasin(Junctions[BN],FlowInfo, JunctionNetwork);
//     //cout << "...got it!" << endl;
//     AllTheBasins.push_back(thisBasin);

//     // This is required if the basins are nested--test the code which numbers
//     // to be overwritten by a smaller basin
//     drainage_of_other_basins[Junctions[BN]] = thisBasin.get_NumberOfCells();
//   }

//   // now loop through everything again getting the raster
//   if (N_Juncs > 0)     // this gets the first raster
//   {
//     BasinMasterRaster = AllTheBasins[0].write_integer_data_to_LSDIndexRaster(Junctions[0], FlowInfo);
//   }

//   // now add on the subsequent basins
//   for(int BN = 1; BN<N_Juncs; BN++)
//   {
//     AllTheBasins[BN].add_basin_to_LSDIndexRaster(BasinMasterRaster, FlowInfo,
//                               drainage_of_other_basins, Junctions[BN]);
//   }


//   // We need to use bil format since there is georeferencing
//   string raster_ext = "bil";
//   // print the basin raster
//   BasinMasterRaster.write_raster(basin_raster_name, raster_ext);
//   cout << "Finished writing the basins!" << endl;

}
