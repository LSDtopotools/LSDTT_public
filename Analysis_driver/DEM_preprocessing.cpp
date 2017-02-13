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

  // some parameters for snapping to basins
  int_default_map["threshold_contributing_pixels"] = 1000;
  int_default_map["search_radius_nodes"] = 10;
  string_default_map["basin_outlet_csv"] = "NULL";

  // set default methods
  bool_default_map["fill_nodata"] = false;
  bool_default_map["remove_low_relief"] = false;
  bool_default_map["write_relief_raster"] = false;
  bool_default_map["check_basin_locations"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  
  
  
  
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
  
  // Now print the new file
  string preprocess_str = "_PP";
  FinalRaster.write_raster(OUT_DIR+OUT_ID+preprocess_str,raster_ext);
  
  // now we check for basins (if you've ask for that)
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
      // read the file


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

      // now get a flow info
  
      // Fill this raster
      float fill_slope = 0.0001;
      cout << "Filling the raster" << endl;
      LSDRaster filled_raster = FinalRaster.fill(fill_slope);
  
      // a boundary condition for the flow info object
      vector<string> boundary_conditionst(4);
      boundary_conditionst[0] = "n";
      boundary_conditionst[1] = "n";
      boundary_conditionst[2] = "n";
      boundary_conditionst[3] = "n";
      boundary_conditions = boundary_conditionst;
      
      // get the flow info
      LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);

      // get contributing pixels (needed for junction network)
      LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  
      // get the sources
      vector<int> sources;
      sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["threshold_contributing_pixels"]);
  
      // now get the junction network
      vector<float> fUTM_easting;
      vector<float> fUTM_northing;
      LSDJunctionNetwork JNetwork(sources, FlowInfo);
      
      if (this_bool_map["print_stream_order_raster"])
      {
        LSDIndexRaster SOArray = JNetwork.StreamOrderArray_to_LSDIndexRaster();
        string SO_raster_name = OUT_DIR+OUT_ID+"_SO";
        SOArray.write_raster(SO_raster_name,raster_ext);
      }
      
      JNetwork.get_x_and_y_from_latlong(latitude, longitude,fUTM_easting,fUTM_northing);

      int threshold_stream_order = 2;
      vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
      vector<int> snapped_node_indices;       // a vector to hold the valid node indices
      vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
      cout << "The search radius is: " << this_int_map["search_radius_nodes"] << endl;
      JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
                this_int_map["search_radius_nodes"], threshold_stream_order, FlowInfo, 
                valid_cosmo_points, snapped_node_indices, snapped_junction_indices);
      
      cout << "The number of valid points is: " << int(valid_cosmo_points.size()) << endl;
      
      // Now get the basin rasters
      string basin_raster_name = OUT_DIR+OUT_ID+"_AllBasins";
      print_basins(basin_raster_name, FlowInfo, JNetwork,
                   snapped_junction_indices);
    }
  }
  
  
  
}
