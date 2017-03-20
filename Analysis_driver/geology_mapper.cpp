//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Geology mapper
//
// This program takes two arguments, the path name and the driver name
// It maps geological data onto csv files
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
#include "../LSDSpatialCSVReader.hpp"



int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the geology mapping tool!                ||" << endl;
    cout << "|| This program has a number of options to burn        ||" << endl;
    cout << "|| geological data onto csv files                      ||" << endl;
    cout << "|| Simon M. Mudd                                       ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./geology_mapping_tool.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_Geol.param" << endl;
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
  
  // The name of the csv file
  string_default_map["csv_for_burning"] = "NULL";
  string_default_map["updated_csv_name"] = "updated.csv";
  string_default_map["geojson_name"] = "new_json.geojson";
  string_default_map["x_column_name"] = "easting";
  string_default_map["y_column_name"] = "northing";
  
  // soe switches that determine the output
  bool_default_map["print_updated_csv"] = true;
  bool_default_map["print_geojson"] = false;
  
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
  LSDIndexRaster geology_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;

  if(this_string_map["csv_for_burning"]=="NULL" || this_string_map["csv_for_burning"]=="Null" || this_string_map["csv_for_burning"]=="null")
  {
    cout << "I don't have a csv file for burning. Please enter one with the parameter flag csv_for_burning" << endl;
    exit(EXIT_SUCCESS);
  }
  else
  {
    // Try to open the file
    LSDSpatialCSVReader CSVFile(RI,DATA_DIR+this_string_map["csv_for_burning"]);
    
    // now check if it has lat and long
    bool have_latlong =  CSVFile.check_if_latitude_and_longitude_exist();
    if (have_latlong)
    {
      cout << "This csv has lat long coordinates" << endl;
    }
    else
    {
      cout << "I didn't find lat-long coordinates" << endl;
      cout << "I am going to get lat-long from x and y coordinates" << endl;
      CSVFile.get_latlong_from_x_and_y(this_string_map["x_column_name"], this_string_map["y_column_name"]);
    }
    
    cout << "Let me check all the data columns" << endl;
    bool have_all_columns = CSVFile.check_if_all_data_columns_same_length();
    if (have_all_columns)
    {
      cout << "Data columns seem to be in order" << endl;
    }
    else
    {
      cout << "You seem to be missing some data in your columns." << endl;
      exit(EXIT_SUCCESS);
    }
    
    
    cout << "I am burning the raster to the csv file." << endl;
    string column_name = "geology";
    CSVFile.burn_raster_data_to_csv(geology_raster,column_name);
    
    
    if(this_bool_map["print_updated_csv"])
    {
      string out_csv_name = OUT_DIR+this_string_map["updated_csv_name"];
      CSVFile.print_data_to_csv(out_csv_name);
    }
    
    if(this_bool_map["print_geojson"])
    {
      cout << "I am printing a geojson file for you." << endl;
      string out_json_name = OUT_DIR+this_string_map["geojson_name"];
      CSVFile.print_data_to_geojson(out_json_name);
    }

  }
}
