//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// DEM_preprocessing
//
// This preprocesses DEMs, checking their 
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDParameterParser.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDBasin.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDShapeTools.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the DEM_preprocessing tool!              ||" << endl;
    cout << "|| This program tried to fix common problems in your   ||" << endl;
    cout << "|| DEM such as funky nodata, wrong projection, etc.    ||" << endl;
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
    cout << "./DEM_preprocessing.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_preprocess.param" << endl;
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
  
  // set default in parameter
  float_default_map["minimum_elevation"] = 0.0;
  float_default_map["maximum_elevation"] = 30000;
  float_default_map["filling_window_radius"] = 50;
  float_default_map["relief_radius"] = 100;
  float_default_map["relief_threshold"] = 50;

  string_default_map["filling_raster_fname"] = "NULL";

  // set default methods
  bool_default_map["fill_nodata"] = false;
  bool_default_map["remove_low_relief"] = false;
  bool_default_map["write_relief_raster"] = false;
  bool_default_map["find_holes"] = false;
  bool_default_map["fill_holes_with_coarse_raster"] = false;
  
  int_default_map["hole_filling_steps"] = 500;
  int_default_map["hole_filling_sweeps"] = 50;

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
  
  // now get rid of the low and high values
  float lower_threshold = this_float_map["minimum_elevation"];
  float upper_threshold = this_float_map["maximum_elevation"];
  bool belowthresholdisnodata = true;
  LSDRaster Flooded = topography_raster.mask_to_nodata_using_threshold(lower_threshold,belowthresholdisnodata);
  belowthresholdisnodata = false;
  LSDRaster FinalRaster = Flooded.mask_to_nodata_using_threshold(upper_threshold,belowthresholdisnodata);
  
  float window_radius = this_float_map["filling_window_radius"];
  if(this_bool_map["fill_nodata"])
  {
    cout << "I am going to fill internal nodata in your raster. This might take a little while." << endl;
    LSDRaster FilledDEM = FinalRaster.alternating_direction_nodata_fill_irregular_raster(window_radius);
    FinalRaster = FilledDEM;
  }
  
  // This will remove low relief areas and replace with NoData
  if(this_bool_map["remove_low_relief"])
  {
    float relief_radius = this_float_map["relief_radius"];
    float relief_threshold = this_float_map["relief_threshold"];
    
    cout << "Let me calculate the relief for you. Radius: " << relief_radius << " thresh: " <<  relief_threshold << endl;
    
    int relief_method = 0;    // THis means a square kernal is used. 
                              // A curcyular kernal is more computationally intensive
                              // and this is a rudimentary calculation so doesn't need
                              // fancy circular windows. 
    LSDRaster Relief =  FinalRaster.calculate_relief(relief_radius, relief_method);
    
    // Now change to nodata
    bool belowthresholdisnodata = true;
    LSDRaster Thresholded = FinalRaster.mask_to_nodata_using_threshold_using_other_raster(relief_threshold,belowthresholdisnodata, Relief);
    FinalRaster = Thresholded; 

    if(this_bool_map["fill_nodata"])
    {
      cout << "I am going to fill internal nodata once more, because you have removed flat areas." << endl;
      LSDRaster FilledDEM = FinalRaster.alternating_direction_nodata_fill_irregular_raster(window_radius);
      FinalRaster = FilledDEM;
    }      
      
    if(this_bool_map["write_relief_raster"])
    {
      string relief_str = "_REL";
      Relief.write_raster(OUT_DIR+OUT_ID+relief_str,raster_ext);
    }
  }
  
  //cout << "Find holes is: " << this_bool_map["find_holes"] << endl;
  if(this_bool_map["find_holes"])
  {
    cout << "I am finding some holes for you. " << endl;
    LSDIndexRaster LookForHoles = FinalRaster.create_binary_isdata_raster();
    
    string LFH_str = "_LFH";
    LookForHoles.write_raster(OUT_DIR+OUT_ID+LFH_str,raster_ext);
  
    LSDIndexRaster Holer = LookForHoles.find_holes_with_nodata_bots(this_int_map["hole_filling_steps"], this_int_map["hole_filling_sweeps"]);
    
    string Holer_str = "_HOLES";
    Holer.write_raster(OUT_DIR+OUT_ID+Holer_str,raster_ext);
  
  }
  
  
  if(this_bool_map["fill_holes_with_coarse_raster"])
  {
    cout << "Getting the locations of the holes. " << endl;
    LSDIndexRaster LookForHoles = FinalRaster.create_binary_isdata_raster();
    
    vector<float> UTME;
    vector<float> UTMN;
    vector<int> rows_of_nodes;
    vector<int> cols_of_nodes;
    LookForHoles.get_points_in_holes_for_interpolation(this_int_map["hole_filling_steps"], 
                                  this_int_map["hole_filling_sweeps"],
                                  UTME,UTMN,rows_of_nodes,cols_of_nodes);
    
    cout << "Number of hole points is: " << UTME.size() << endl;
    
    
    string filling_raster_fname = DATA_DIR+this_string_map["filling_raster_fname"];
    cout << "I am getting the filling raster for interpolation. " << endl;
    LSDRaster filling_raster(filling_raster_fname,raster_ext);
    
    vector<float> filled_data = filling_raster.interpolate_points_bilinear(UTME, UTMN);
    
    LSDRaster HolesFilled = FinalRaster.fill_with_interpolated_data(rows_of_nodes,cols_of_nodes, filled_data);
    
    string HF_str = "_HF";
    HolesFilled.write_raster(OUT_DIR+OUT_ID+HF_str,raster_ext);
  }
  

  // Now write the final raster
  string preprocess_str = "_PP";
  FinalRaster.write_raster(OUT_DIR+OUT_ID+preprocess_str,raster_ext);
  
}
