//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_mapping_tool
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate 
// different kinds of chi analysis
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDParameterParser.hpp"

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
  cout << "Finished with exporting basins!" << endl;

}



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
  int_default_map["threshold_contributing_pixels"] = 2000;
  int_default_map["search_radius_nodes"] = 10;
  
  // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;
  
  // set default methods
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["print_basin_raster"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_stream_order_csv"] = false;
  bool_default_map["check_basin_locations"] = false;
  bool_default_map["print_junction_csv"] = false;
  bool_default_map["spawn_basins_from_outlets"] = false;
  
  // set default string method
  string_default_map["slope method"] = "polynomial";
  string_default_map["averaging_raster_vector"] = "NULL";
  string_default_map["basin_outlet_csv"] = "NULL";
  
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
  
  
  // now we check for basins (if you've ask for that)
  vector<int> SnappedJunctions;
  if(this_bool_map["check_basin_locations"])
  {
    // first we see if the nodes file exists
    string basin_outlet_fname = DATA_DIR+this_string_map["basin_outlet_csv"];
    
    cout << "I am checking basins, the outlet points are: " <<  basin_outlet_fname << endl;
  
    ifstream ifs; 
    ifs.open(basin_outlet_fname.c_str());
  
    // check if the parameter file exists
    if( ifs.fail() )
    {
      cout << "\nFATAL ERROR: The basin points file \"" << basin_outlet_fname
           << "\" doesn't exist" << endl;
    }
    else
    {
      // the LL data
      vector<float> latitude;
      vector<float> longitude;
      vector<string> IDs;
      string line_from_file;
      vector<string> this_string_vec;
      vector<string> empty_string_vec;
      
      // now loop through the rest of the lines, getting the data. 
      while( getline(ifs, line_from_file))
      {
        // reset the string vec
        this_string_vec = empty_string_vec;
    
        // create a stringstream
        stringstream ss(line_from_file);
    
        while( ss.good() )
        {
          string substr;
          getline( ss, substr, ',' );
      
          // remove the spaces
          substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
          // remove control characters
          substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
          // add the string to the string vec
          this_string_vec.push_back( substr );
        }
    
        //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
        if ( int(this_string_vec.size()) < 3)
        {    
          cout << "Hey there, I am trying to load your basin outlet data but you seem not to have" << endl;
          cout << "enough columns in your file. I am ignoring a line" << endl;
        }
        else
        {
          // now convert the data
          //cout << "Getting sample name: " <<  this_string_vec[0] << endl;
      
          // let the user know about offending underscores, and replace them
          string s = this_string_vec[0];
          string uscore = "_";
          size_t found = s.find(uscore);
          if (found!=std::string::npos)
          {
            cout << "I found an underscore in the sample name. Replacing with a dash." <<endl;
            replace( s.begin(), s.end(), '_', '-');
            cout << "New sample name is: " << s << endl;
          }

          IDs.push_back( s );
          latitude.push_back( atof( this_string_vec[1].c_str() ) );
          longitude.push_back( atof(this_string_vec[2].c_str() ) );
        }
      }
      ifs.close();

      int N_lats = int(latitude.size());
      cout << "N outlets: " << N_lats << endl;
      for (int i = 0; i< N_lats; i++)
      {
        cout << "lat: " << latitude[i] << " and long: " <<  longitude[i] << endl;
      }

      // now get the junction network
      vector<float> fUTM_easting;
      vector<float> fUTM_northing;

      // find the x any y locations of the points for snapping
      JunctionNetwork.get_x_and_y_from_latlong(latitude, longitude,fUTM_easting,fUTM_northing);

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
      string basin_raster_name = OUT_DIR+OUT_ID+"_AllBasins";
      print_basins(basin_raster_name, FlowInfo, JunctionNetwork,
                   snapped_junction_indices);
      SnappedJunctions = snapped_junction_indices;
    }
  }


  if (this_bool_map["spawn_basins_from_outlets"])
  {
    //========================
    // LOOPING THROUGH BASINS
    //========================
    // now loop through the valid points, getting the cosmo data 
    int n_valid_points =  int(SnappedJunctions.size());
    int padding_pixels = 100;
    cout << "I am trying to spawn basins, found " << n_valid_points << " valid points." << endl;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
       
      cout << "Getting basin" << endl;
      LSDBasin thisBasin(SnappedJunctions[samp],FlowInfo, JunctionNetwork);
      
      cout << "Now getting the basin raster" << endl;  
      LSDRaster BasinRaster = thisBasin.TrimPaddedRasterToBasin(padding_pixels, 
                                           FlowInfo,FillRaster);
      cout << "Formatting the path" << endl;

      //string DEM_newpath = ReformatPath(DEMpath);
      string DEMnewname = OUT_DIR+OUT_ID+"_Spawned_"+itoa(samp);
      cout << "Writing a new basin to: " << DEMnewname << endl;
      BasinRaster.write_raster(DEMnewname,"bil");
    }
  }


  // and get the basins
  cout << "I am getting the basins with basin order: " << this_int_map["basin_order"] << endl;
  vector< int > basin_junctions = JunctionNetwork.ExtractBasinJunctionOrder(this_int_map["basin_order"], FlowInfo);
  LSDIndexRaster Basin_Raster = JunctionNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);

  cout << "I extracted " << basin_junctions.size() << " basins!" << endl;

  LSDBasin A_basin(basin_junctions[0],FlowInfo,JunctionNetwork);
  

  // print stream order to file if wanted
  if(this_bool_map["print_stream_order_raster"])
  {
    LSDIndexRaster SO_array = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
    string SOExt = "_SO";
    SO_array.write_raster(OUT_DIR+OUT_ID+SOExt, raster_ext);
  }

  // print basins to file if wanted
  if(this_bool_map["print_basin_raster"])
  {
    //LSDIndexRaster Basin_Raster = A_basin.write_integer_data_to_LSDIndexRaster(basin_junctions[0], FlowInfo);
    string BRExt = "_BR";
    Basin_Raster.write_raster(OUT_DIR+OUT_ID+BRExt, raster_ext);
  }
  
  // print stream order CSV to file if wanted
  if(this_bool_map["print_stream_order_csv"])
  {
    string CN_fname = OUT_DIR+OUT_ID+"_ChanNet";
    JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, CN_fname);
  }


  if(this_bool_map["print_junction_csv"])
  {
    string JN_fname = OUT_DIR+OUT_ID+"_Junctions.csv";
    vector<int> JunctionList;
    JunctionNetwork.print_junctions_to_csv(FlowInfo, JunctionList,JN_fname);
  }
  

}
