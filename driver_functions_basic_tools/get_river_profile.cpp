//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// get_river_profile.cpp
// This function calcualtes chi across a landscape. Run without
// arguments for instructions
//
// Developed by:
//  Simon M. Mudd
//
// Copyright (C) 2013 Simon M. Mudd 2013
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#include <iostream>
#include <sstream>
#include <cstdlib>
//#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDJunctionNetwork.hpp"
using namespace std;


int main(int argc, char *argv[])
{

  if (argc!=3)
  {
    cout << "===========================================================" << endl;
    cout << "|| Welcome to the river profile tool!                    ||" << endl;
    cout << "|| I am used to calculate the chi coordinate from a DEM. ||" << endl;
    cout << "===========================================================" << endl;
    cout << "This program requires two inputs: " << endl << endl;
    cout << "* First the path to the DEM." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "   (Either \\ or / depending on your operating system.)" << endl << endl;
    cout << "* Second the name of your parameter file " << endl;
    cout << "============================================================" << endl;
    cout << "An example call is: " << endl;
    cout << "./get_river_profile.out ~home/basins/Chile/ LSDTT_river_profile.param" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

	// maps for setting default parameters
	map<string,int> int_default_map;
	map<string,float> float_default_map;
	map<string,bool> bool_default_map;
	map<string,string> string_default_map;

  // set default int parameters
  int_default_map["search_radius"] = 10;
  int_default_map["Threshold_SO"] = 1;

	// set default float parameters
	float_default_map["threshold_contributing_pixels"] = 400;

  // set default bool parameters
  bool_default_map["Ingest_Channel_Heads"] = false;

  // set default string parameters
  string_default_map["CHeads_file"] = "NULL";
	string_default_map["coords_csv_file"] = "NULL";

	// Use the parameter parser to get the maps of the parameters required for the
	// analysis
	// load parameter parser object
	LSDParameterParser LSDPP(path_name,f_name);
	LSDPP.force_bil_extension();

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
  string DEM_extension =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();

  // some error checking
  if (this_string_map["coords_csv_file"] == "NULL")
  {
    cout << "FATAL ERROR: I can't read the coordinates of your profile points. Please include the CSV filename in your parameter file!" << endl;
    exit(EXIT_SUCCESS);
  }

  // load the raster
  LSDRaster new_raster(DATA_DIR+DEM_ID,DEM_extension);

  // remove the sea (seems to be required if gdal is used in places with nodata)
  new_raster.remove_seas();

  // get the filled raster
  float min_slope = 0.0001;
  LSDRaster filled_raster = new_raster.fill(min_slope);

  cout << "Filled raster, getting flow info" << endl;
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  cout << "Got flow info" << endl;

  vector<int> sources;
  // getting the channel network
  if (CHeads_file == "NULL")
  {
    cout << "I'm using an area threshold to get your channel network..." << endl;
    // get some relevant rasters
    LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
    LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["threshold_contributing_pixels"]);
  }
  else
  {
    cout << "I'm reading in your channel heads from the csv file..." << endl;
    sources = FlowInfo.Ingest_Channel_Heads(DATA_DIR+CHeads_file,0);
  }
  cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;

  // reading in the csv file with the lat long points
  cout << "\t Reading in the csv file" << endl;
  LSDSpatialCSVReader ProfilePoints(filled_raster, DATA_DIR+this_string_map["coords_csv_file"]);
  vector<float> UTME;
  vector<float> UTMN;
  string column_name = "Point";
  vector<string> point_locations;
  ProfilePoints.get_data_in_raster_for_snapping(column_name, UTME, UTMN, point_locations);
  cout << "\t Got the x and y locations" << endl;
  string csv_outname = "_UTM_check.csv";
  ProfilePoints.print_UTM_coords_to_csv(UTME, UTMN, (DATA_DIR+DEM_ID+csv_outname));

  // snap to nearest channel
	vector<int> valid_indices;
	vector<int> snapped_nodes;
	vector<int> snapped_JNs;
	ChanNetwork.snap_point_locations_to_channels(UTME, UTMN, this_int_map["search_radius"], this_int_map["Threshold_SO"], FlowInfo, valid_indices, snapped_nodes, snapped_JNs);

  cout << "\t Getting the upstream and downstream points of your channel..." << endl;
  int upstream_node, downstream_node, downstream_jn;
  int found_points = 0;
  for (int i = 0; i < int(point_locations.size()); i++)
  {
    if (point_locations[i] == "upstream")
    {
      upstream_node = snapped_nodes[i];
      found_points++;
    }
    if (point_locations[i] == "downstream")
    {
      downstream_node = snapped_nodes[i];
      downstream_jn = snapped_JNs[i];
      found_points++;
    }
  }

  if (int(valid_indices.size()) == 2 || found_points == 2)
  {
    cout << "I've read in both the upstream and downstream points! Now getting your river profile..." << endl;
  }
  else
  {
    cout << "I couldn't read in all your points. Please check your latitude and longitude coordinates..." << endl;
  }

  LSDIndexChannel ThisChannel(upstream_node, downstream_node, FlowInfo);
  string jn_name = itoa(downstream_jn);
  string out_name = "_profile_"+jn_name;
  ThisChannel.write_channel_to_csv(DATA_DIR,DEM_ID+out_name);
  cout << "Done!" << endl;
}
