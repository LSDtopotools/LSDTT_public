//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiTools
// Land Surface Dynamics ChiTools object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for performing various analyses in chi space
//
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDChiTools.cpp
// LSDChiTools object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDChiTools_CPP
#define LSDChiTools_CPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDChiTools.hpp"
#include "LSDBasin.hpp"
#include "LSDChiNetwork.hpp"
using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDRaster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDRaster& ThisRaster)
{
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDRaster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDIndexRaster& ThisRaster)
{
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDFlowInfo
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDFlowInfo& ThisFI)
{
  NRows = ThisFI.get_NRows();
  NCols = ThisFI.get_NCols();
  XMinimum = ThisFI.get_XMinimum();
  YMinimum = ThisFI.get_YMinimum();
  DataResolution = ThisFI.get_DataResolution();
  NoDataValue = ThisFI.get_NoDataValue();
  GeoReferencingStrings = ThisFI.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDFlowInfo
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDJunctionNetwork& ThisJN)
{
  NRows = ThisJN.get_NRows();
  NCols = ThisJN.get_NCols();
  XMinimum = ThisJN.get_XMinimum();
  YMinimum = ThisJN.get_YMinimum();
  DataResolution = ThisJN.get_DataResolution();
  NoDataValue = ThisJN.get_NoDataValue();
  GeoReferencingStrings = ThisJN.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This resets all the data maps
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::reset_data_maps()
{
  map<int,float> empty_map;
  vector<int> empty_vec;

  M_chi_data_map = empty_map;
  b_chi_data_map = empty_map;
  elev_data_map = empty_map;
  chi_data_map = empty_map;
  flow_distance_data_map = empty_map;
  drainage_area_data_map = empty_map;
  node_sequence = empty_vec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc)
{

  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;

  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
// Same as above but with floats
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc)
{

  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;

  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to convert a node position with a row and column to a lat
// and long coordinate
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_lat_and_long_locations(int row, int col, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter)
{
  // get the x and y locations of the node
  double x_loc,y_loc;
  get_x_and_y_locations(row, col, x_loc, y_loc);

  // get the UTM zone of the node
  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);
  //cout << endl << endl << "Line 1034, UTM zone is: " << UTM_zone << endl;


  if(UTM_zone == NoDataValue)
  {
    lat = NoDataValue;
    longitude = NoDataValue;
  }
  else
  {
    // set the default ellipsoid to WGS84
    int eId = 22;

    double xld = double(x_loc);
    double yld = double(y_loc);

    // use the converter to convert to lat and long
    double Lat,Long;
    Converter.UTMtoLL(eId, yld, xld, UTM_zone, is_North, Lat, Long);


    lat = Lat;
    longitude = Long;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_UTM_information(int& UTM_zone, bool& is_North)
{

  // set up strings and iterators
  map<string,string>::iterator iter;

  //check to see if there is already a map info string
  string mi_key = "ENVI_map_info";
  iter = GeoReferencingStrings.find(mi_key);
  if (iter != GeoReferencingStrings.end() )
  {
    string info_str = GeoReferencingStrings[mi_key] ;

    // now parse the string
    vector<string> mapinfo_strings;
    istringstream iss(info_str);
    while( iss.good() )
    {
      string substr;
      getline( iss, substr, ',' );
      mapinfo_strings.push_back( substr );
    }
    UTM_zone = atoi(mapinfo_strings[7].c_str());
    //cout << "Line 1041, UTM zone: " << UTM_zone << endl;
    //cout << "LINE 1042 LSDRaster, N or S: " << mapinfo_strings[7] << endl;

    // find if the zone is in the north
    string n_str = "n";
    string N_str = "N";
    is_North = false;
    size_t found = mapinfo_strings[8].find(N_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    //cout << "is_North is: " << is_North << endl;

  }
  else
  {
    UTM_zone = NoDataValue;
    is_North = false;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname,
                                 float A_0, float m_over_n, float area_threshold)
{

  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());

  chi_map_csv_out.precision(9);

  float chi_coord;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  chi_map_csv_out << "latitude,longitude,chi" << endl;

  LSDRaster Chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0, area_threshold);

  float NDV = Chi.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      chi_coord =  Chi.get_data_element(row,col);

      if (chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out << latitude << "," << longitude  << "," << chi_coord << endl;
      }
    }
  }

  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2. You feed it the chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname,
                                 LSDRaster& chi_coord)
{

  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());



  float this_chi_coord;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  chi_map_csv_out << "latitude,longitude,chi" << endl;

  float NDV = chi_coord.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      this_chi_coord = chi_coord.get_data_element(row,col);

      if (this_chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out.precision(9);
        chi_map_csv_out << latitude << "," << longitude  << ",";
        chi_map_csv_out.precision(5);
        chi_map_csv_out << this_chi_coord << endl;
      }
    }
  }

  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2. You feed it the chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname,
                                 LSDRaster& chi_coord, LSDIndexRaster& basin_raster)
{

  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());



  float this_chi_coord;
  int this_basin_number;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  chi_map_csv_out << "latitude,longitude,chi,basin_junction" << endl;

  float NDV = chi_coord.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      this_chi_coord = chi_coord.get_data_element(row,col);
      this_basin_number = basin_raster.get_data_element(row,col);


      if (this_chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out.precision(9);
        chi_map_csv_out << latitude << "," << longitude  << ",";
        chi_map_csv_out.precision(5);
        chi_map_csv_out << this_chi_coord << ",";
        chi_map_csv_out << this_basin_number << endl;
      }
    }
  }

  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating segments from all sources in a DEM
// The sources and their outlets are supplied by the source and outlet nodes
// vectors. These are generated from the LSDJunctionNetwork function
// get_overlapping_channels
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int target_nodes,
                                    int n_iterations, int skip,
                                    int minimum_segment_length, float sigma)
{

  // IMPORTANT THESE PARAMETERS ARE NOT USED BECAUSE CHI IS CALUCALTED SEPARATELY
  // However we need to give something to pass to the Monte carlo functions
  // even through they are not used (they are inherited)
  float A_0 = 1;
  float m_over_n = 0.5;

  // These elements access the chi data
  vector< vector<float> > chi_m_means;
  vector< vector<float> > chi_b_means;
  vector< vector<float> > chi_coordinates;
  vector< vector<int> > chi_node_indices;

  // these are for the individual channels
  vector<float> these_chi_m_means;
  vector<float> these_chi_b_means;
  vector<float> these_chi_coordinates;
  vector<int> these_chi_node_indices;

  // these are maps that will store the data
  map<int,float> m_means_map;
  map<int,float> b_means_map;
  map<int,float> chi_coord_map;
  map<int,float> elev_map;
  map<int,float> area_map;
  map<int,float> flow_distance_map;
  vector<int> node_sequence_vec;

  // these are vectors that will store information about the individual nodes
  // that allow us to map the nodes to specific channels during data visualisation

  // These two maps have each node in the channel (the index)
  // linked to a key (either the baselevel key or source key)
  map<int,int> these_source_keys;
  map<int,int> these_baselevel_keys;

  // These two maps link a keys, which are incrmented by one, to the
  // junction or node of the baselevel or source
  map<int,int> this_key_to_source_map;
  map<int,int> this_key_to_baselevel_map;

  // these are for working with the FlowInfo object
  int this_node,row,col;
  int this_base_level, this_source_node;

  // get the number of channels
  int source_node_tracker = -1;
  int baselevel_tracker = -1;
  int n_channels = int(source_nodes.size());
  for(int chan = 0; chan<n_channels; chan++)
  {
    cout << "Sampling channel " << chan+1 << " of " << n_channels << endl;

    // get the base level
    this_base_level = FlowInfo.retrieve_base_level_node(source_nodes[chan]);
    //cout << "Got the base level" << endl;

    // If a key to this base level does not exist, add one.
    if ( this_key_to_baselevel_map.find(this_base_level) == this_key_to_baselevel_map.end() )
    {
      baselevel_tracker++;
      cout << "Found a new baselevel. The node is: " << this_base_level << " and key is: " << baselevel_tracker << endl;
      this_key_to_baselevel_map[this_base_level] = baselevel_tracker;
    }

    // now add the source tracker
    source_node_tracker++;
    this_source_node = source_nodes[chan];
    this_key_to_source_map[this_source_node] = source_node_tracker;

    cout << "The source key is: " << source_node_tracker << " and basin key is: " << baselevel_tracker << endl;

    // get this particular channel (it is a chi network with only one channel)
    LSDChiNetwork ThisChiChannel(FlowInfo, source_nodes[chan], outlet_nodes[chan],
                                Elevation, FlowDistance, DrainageArea,chi_coordinate);

    // split the channel
    //cout << "Splitting channels" << endl;
    ThisChiChannel.split_all_channels(A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma);

    // monte carlo sample all channels
    //cout << "Entering the monte carlo sampling" << endl;
    ThisChiChannel.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations, skip, minimum_segment_length, sigma);

    // okay the ChiNetwork has all the data about the m vales at this stage.
    // Get these vales and print them to a raster
    chi_m_means = ThisChiChannel.get_m_means();
    chi_b_means = ThisChiChannel.get_b_means();
    chi_coordinates = ThisChiChannel.get_chis();
    chi_node_indices = ThisChiChannel.get_node_indices();

    // now get the number of channels. This should be 1!
    int n_channels = int(chi_m_means.size());
    if (n_channels != 1)
    {
      cout << "Whoa there, I am trying to make a chi map but something seems to have gone wrong with the channel extraction."  << endl;
      cout << "I should only have one channel per look but I have " << n_channels << " channels." << endl;
    }

    // now get the m_means out
    these_chi_m_means = chi_m_means[0];
    these_chi_b_means = chi_b_means[0];
    these_chi_coordinates = chi_coordinates[0];
    these_chi_node_indices = chi_node_indices[0];

    //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;


    int n_nodes_in_channel = int(these_chi_m_means.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {

      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;

      // only take the nodes that have not been found
      if (m_means_map.find(this_node) == m_means_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);

        //cout << "This is a new node; " << this_node << endl;
        m_means_map[this_node] = these_chi_m_means[node];
        b_means_map[this_node] = these_chi_b_means[node];
        chi_coord_map[this_node] = these_chi_coordinates[node];
        elev_map[this_node] = Elevation.get_data_element(row,col);
        area_map[this_node] = DrainageArea.get_data_element(row,col);
        flow_distance_map[this_node] = FlowDistance.get_data_element(row,col);
        node_sequence_vec.push_back(this_node);

        these_source_keys[this_node] = source_node_tracker;
        these_baselevel_keys[this_node] = baselevel_tracker;

      }
      else
      {
        //cout << "I already have node: " << this_node << endl;
      }
    }
  }

  cout << "I am all finished segmenting the channels!" << endl;

  // set the object data members
  M_chi_data_map =m_means_map;
  b_chi_data_map = b_means_map;
  elev_data_map = elev_map;
  chi_data_map = chi_coord_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_sequence_vec;

  source_keys_map = these_source_keys;
  baselevel_keys_map = these_baselevel_keys;
  key_to_source_map = this_key_to_source_map;
  key_to_baselevel_map = this_key_to_baselevel_map;

  // get the fitted elevations
  calculate_segmented_elevation(FlowInfo);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is a much more rudimentary version that mimics the
// channel steepness caluclations.
// chi needs tobe calculated outside of the function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator_rudimentary(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int regression_nodes)
{

  // the data is stored in maps, for easier testing if a node has been
  // visited.
  // You might consider having these as data elements in the object so you don't
  // have to pass them
  map<int,float> gradient_data_map;
  map<int,float> intercept_data_map;
  map<int,float> R2_data_map;
  map<int,float> chi_coordinate_data_map;
  map<int,float> elevation_data_map;
  map<int,float> flow_distance_map;
  map<int,float> area_map;
  vector<int> node_order;

  // check if the number of nodes are odd .If not add 1
  if (regression_nodes % 2 == 0)
  {
    cout << "Hello user. You need an odd number of regression nodes." << endl;
    regression_nodes = regression_nodes+1;
    cout << " Changing your regression nodes to " << regression_nodes << endl;
  }

  // now get the midpoint
  int mp_nodes = (regression_nodes-1)/2;

  //cout << "The number of mp nodes is: " << mp_nodes << endl;

  // these keep track of the beginning and ending nodes of a given channel
  int channel_start_node;
  int channel_end_node;
  float channel_end_elevation;

  // vectors for holding the chi elevation data
  vector<float> chi_vec;
  vector<float> elev_vec;
  vector<float> empty_vec;

  // these are extracted from the channel segments using linear regression
  float intercept,gradient,R_squared;

  //float this_chi;
  //float this_elev;
  int this_mp_node;
  int this_end_node;
  int this_start_node;

  // these are for getting information out of the FlowInfo object
  int row,col, this_node;
  int r_node, r_row,r_col;          // reciever row and column.

  // The way this works is that it starts at the top of a channel. It then works
  // its way down and find the node that is the midpoint and the node that is the
  // end point. The midpoint node is where the data will be recorded.
  // It then puts the data from the start node to the end node into a vector
  // and performs a linear regression of this vector. The regression data from these
  // vectors are recorded at the nodes.
  // We then want to cover all the nodes with data so what happens if some nodes
  // do not become midpoints?
  // We could start at the top and get the first midpoint.
  // From there we can work our way down checking if the top of the regression segment
  // is more than one node down from the end point...

  // get the number of channels
  int n_channels = int(source_nodes.size());
  // now loop through the channels
  for(int chan = 0; chan<n_channels; chan++)
  {
    channel_start_node = source_nodes[chan];
    channel_end_node = outlet_nodes[chan];

    // Get the elevation of the end node as a secondary check of the ending of the channel
    // segment
    FlowInfo.retrieve_current_row_and_col(channel_end_node,row,col);
    channel_end_elevation = Elevation.get_data_element(row,col);

    // reset the flag for ending the channel
    bool is_end_of_channel = false;

    // set the segment start node to the channel start node
    this_start_node = channel_start_node;

    // now retrieve the midpoint node
    this_node = channel_start_node;
    for(int n = 0; n<mp_nodes; n++)
    {
      FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
      this_node = r_node;
    }
    this_mp_node = this_node;
    this_node = r_node;

    // now go down one step
    FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
    this_node = r_node;

    // now get the end node
    for(int n = 0; n<mp_nodes; n++)
    {
      FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
      this_node = r_node;
    }
    this_end_node = this_node;

    //================================================
    // This loop is for bug checking
    //this_node = this_start_node;
    //do
    //{
    //  // get the elevation and chi vectors by following the flow
    //  cout << "This node is: " << this_node << endl;
    //  FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    //  FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
    //  this_node = r_node;
    //}
    //while(this_node != this_end_node);
    //
    //cout << "And the midpoint node was: " << this_mp_node << endl;
    //================================================

    // we search down the channel, collecting linear regressions at the
    // midpoint of the intervals
    while (not is_end_of_channel)
    {
      // get a vector of chi and elevation from the start node to the end node
      chi_vec = empty_vec;
      elev_vec = empty_vec;

      // copy the data elements into the vecotrs. This is a little stupid
      // because one might just use a deque to pop the first element
      // and push the last, but the linear regression takes vectors,
      // not deques so you would have to copy the deque element-wise anyway
      // If you wanted, you could speed this up by implementing a linear regression
      // of deques, but that will need to wait for another day.
      this_node = this_start_node;
      do
      {
        // get the elevation and chi vectors by following the flow
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        chi_vec.push_back(chi_coordinate.get_data_element(row,col));
        elev_vec.push_back(Elevation.get_data_element(row,col));

        FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
        this_node = r_node;
      } while(this_node != this_end_node);

      // do a linear regression on the segment
      least_squares_linear_regression(chi_vec,elev_vec, intercept, gradient, R_squared);

      // now add the intercept and gradient data to the correct node
      // only take data that has not been calculated before
      // The channels are in order of descending length so data from
      // longer channels take precidence.
      if (gradient_data_map.find(this_mp_node) == gradient_data_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_mp_node,row,col);
        gradient_data_map[this_mp_node] = gradient;
        intercept_data_map[this_mp_node] = intercept;
        R2_data_map[this_mp_node] = R_squared;
        chi_coordinate_data_map[this_mp_node] = chi_coordinate.get_data_element(row,col);
        elevation_data_map[this_mp_node] = Elevation.get_data_element(row,col);
        flow_distance_map[this_mp_node] = FlowDistance.get_data_element(row,col);
        area_map[this_mp_node] = DrainageArea.get_data_element(row,col);
        node_order.push_back(this_mp_node);
      }
      else
      {
        is_end_of_channel = true;
      }

      // now move all the nodes down one
      FlowInfo.retrieve_receiver_information(this_start_node,r_node,r_row,r_col);
      this_start_node = r_node;

      FlowInfo.retrieve_receiver_information(this_mp_node,r_node,r_row,r_col);
      this_mp_node = r_node;

      FlowInfo.retrieve_receiver_information(this_end_node,r_node,r_row,r_col);
      this_end_node = r_node;

      // check if we are at the end of the channel
      if (this_end_node == channel_end_node)
      {
        is_end_of_channel = true;
      }
      // also check if the end node is lower elevation than the end node,
      // just to try and stop the channel passing the end node
      FlowInfo.retrieve_current_row_and_col(this_end_node,row,col);
      if (channel_end_elevation > Elevation.get_data_element(row,col))
      {
        is_end_of_channel = true;
      }
    }          // This finishes the regression segment loop
  }            // This finishes the channel and resets channel start and end nodes

  // set the data objects
  M_chi_data_map = gradient_data_map;
  b_chi_data_map = intercept_data_map;
  elev_data_map = elevation_data_map;
  chi_data_map = chi_coordinate_data_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_order;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the basins and an additional file that has basin centroids
// and labelling information for plotting
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDChiTools::get_basin_raster(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
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

  return BasinMasterRaster;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to tag channels with a segment number
// It decides on segments if the M_Chi value has changed so should only be used
// with chi networks that have used a skip of 0 and a monte carlo itertions of 1
// This data is used by other routines to look at the spatial distribution of
// hillslope-channel coupling.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::segment_counter(LSDFlowInfo& FlowInfo)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node;
  int segment_counter = 0;
  map<int,int> this_segment_counter_map;
  float last_M_chi, this_M_chi;

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    this_node = node_sequence[0];
    last_M_chi =  M_chi_data_map[this_node];

    for (int n = 0; n< n_nodes; n++)
    {

      // Get the M_chi from the current node
      this_node = node_sequence[n];
      this_M_chi = M_chi_data_map[this_node];

      // If the M_chi has changed, increment the segment counter
      if (this_M_chi != last_M_chi)
      {
        segment_counter++;
        last_M_chi = this_M_chi;
      }

      // Print the segment counter to the data map
      this_segment_counter_map[this_node]  = segment_counter;
    }
  }
  segment_counter_map = this_segment_counter_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to tag channels with a segment number
// It decides on segments if the M_Chi value has changed so should only be used
// with chi networks that have used a skip of 0 and a monte carlo itertions of 1
// This data is used by other routines to look at the spatial distribution of
// hillslope-channel coupling.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::segment_counter_knickpoint(LSDFlowInfo& FlowInfo, float threshold_knickpoint, float threshold_knickpoint_length)
{
  // these are for extracting element-wise data from the channel profiles.
  int abs_threshhold_knickpoint = abs (threshold_knickpoint);
  int this_node = 0;
  int segment_counter_knickpoint = 0; // count the number of knickpoints
  int segment_counter = 0; // count the number of segments
  map<int,float> this_segment_counter_knickpoint_map;
  map<int,int> this_segment_counter_map;
  map<int,int> this_segment_knickpoint_sign_map;
  map<int,int> this_segment_length_map;
  float last_M_chi, this_M_chi;
  float delta_m = 0; // difference between last and new m_chi
  int knickpoint_sign = 0; // sign of the knickpoint: + =1 and - = -1
  float temp_delta_m = 0; // debugging stuff
  float this_segment_length = 0;
  int last_node = 0;
  int n_nodes_segment = 0;
  float x1_temp =0;
  float y1_temp =0;
  float x2_temp =0;
  float y2_temp =0;
  bool same_channel = true;
  float max_knickpoint_value =0;
  int new_knickpoint_counter =0;


  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    this_node = node_sequence[0];
    last_M_chi =  M_chi_data_map[this_node];

    for (int n = 0; n< n_nodes; n++)
    {

      // set the nodes number and keep information about the previous one
      if(n>0)
      {
        last_node = this_node;
      }
      this_node = node_sequence[n];
      // Get the M_chi from the current node
      this_M_chi = M_chi_data_map[this_node];
      // increment the segment node counter
      n_nodes_segment++;

      // If the M_chi has changed, do stuffs
      if (this_M_chi != last_M_chi)
      {
        segment_counter++; // increment the  segment counter
        delta_m= last_M_chi/this_M_chi; // Ratio between last and new chi steepness
        if(delta_m<1){knickpoint_sign = -1;} else {knickpoint_sign = 1;} // Assign the knickpoint sign value
        //delta_m = abs(delta_m); // required to detect all the knickpoints
        //if(delta_m > temp_delta_m) {temp_delta_m = delta_m;} // debugging stuff
        // now checking if the difference of m_chi between two segment is not due to a channel change
          // first retrieving xy coordinates
        FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
        FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
          // Then check if the distance betweenthe two is more than 2 nodes (distance between two points via pytagore or thing like this)
        if (sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp)) > (2*FlowInfo.get_DataResolution()))
        {
          same_channel = false;
        }
          // done
        // Check if the threshold is (over)reached
        //if(delta_m > abs_threshhold_knickpoint && same_channel) if we are are using the threshold value
        if(true) // useless thing
        {
          segment_counter_knickpoint++; // number of knickpoints
          this_segment_counter_knickpoint_map[this_node] = delta_m; // adding the knickpoint value
          if(delta_m>max_knickpoint_value){max_knickpoint_value=delta_m; cout << max_knickpoint_value <<"||"<<this_segment_counter_knickpoint_map[this_node] << endl; } // assign the new knickpoint max value
        }
        //this_segment_length = n_nodes_segment * FlowInfo.get_DataResolution(); // getting the length of the segment using the resolution * number of nodes
        same_channel = true; // Set back the same channel parameter to true
        // now assign the segment lenght to all the point of the segment
        for(int i = n_nodes_segment ; i > 0 ; i--)
        {
          this_segment_length_map[node_sequence[n - i]] =  this_segment_length;
        }

        last_M_chi = this_M_chi;
        this_segment_knickpoint_sign_map[this_node] = knickpoint_sign; // assign the segment polarity
        n_nodes_segment = 0; // reinitialyse the # of node for the next segment
        this_segment_length = 0; // same
      }
      else
      {
        // incrementing the segment length
        if (n>0)
        {
          // Now calculating the distance between the two last nodes

          // Retrieving the x/y coordinates of the last two nodes
          FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
          FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
          // calculate and increment the distance from the last node
          this_segment_length += sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
          //cout << "distance : " << this_segment_length;

        }
      }


      // Print the segment counter to the data map
      this_segment_counter_map[this_node]  = segment_counter;
    }

    if(threshold_knickpoint_length>0)
    {


      // calculating equivalent between nodes_of_knickpoints and nodes to be able to properly navigate inside the map
      int nodes_of_knickpoints [segment_counter_knickpoint];
      int knickpoint_id = 0;
      int last_knickpoint_id = 0;
      for(int i = 0; i< n_nodes; i++)
      {
        this_node = node_sequence[i];
        if(this_segment_counter_knickpoint_map.count(this_node)) // supposed to test if this map has a value assigned for this key
        {
          nodes_of_knickpoints[knickpoint_id] = this_node;
          last_knickpoint_id = knickpoint_id; // stock the last knickpoint_id
          knickpoint_id++;
        }
      }
      // setting up the calculation for the length threshold
      bool still_processing = true;
      bool still_processing_total = true;
      float old_max_knickpoint_value = max_knickpoint_value;
      float distance_to_process_down = threshold_knickpoint_length/2; // represents the threshold down the knickpoint
      float distance_to_process_up = threshold_knickpoint_length/2; // represents the threshold up the knickpoint
      int number_of_erase = 0;

      while (still_processing_total)
      {
        still_processing = true;
        vector <int> knickpoint_to_delete;
        knickpoint_id = 0;

        map<int, float>::iterator it = this_segment_counter_knickpoint_map.begin();
        while(still_processing && it != this_segment_counter_knickpoint_map.end())
        {
          //cout<<it->first<<" :: "<<it->second<<endl;
          if(it->second == max_knickpoint_value)
          {
            distance_to_process_down = threshold_knickpoint_length/2;
            distance_to_process_up = threshold_knickpoint_length/2;
            // let's know test the distance before and beyond this nodes
            for(int g = 1; distance_to_process_down>0 || distance_to_process_up >0 ;g++)
            {

              //Calculate the case down
              if(this_segment_counter_knickpoint_map.count(nodes_of_knickpoints[knickpoint_id - g]) && distance_to_process_down>0)
              {
                //cout << "down exists" << endl;
                this_node = nodes_of_knickpoints[knickpoint_id - g];
                last_node = nodes_of_knickpoints[knickpoint_id - g+1];
                FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
                float temp_distance = sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                if(temp_distance < distance_to_process_down)
                {
                  distance_to_process_down -= temp_distance;
                  knickpoint_to_delete.push_back(this_node);
                } else{distance_to_process_down = 0;}
              }
              else{distance_to_process_down = 0;}
              // calculate the case up
              if(this_segment_counter_knickpoint_map.count(nodes_of_knickpoints[knickpoint_id + g]) && distance_to_process_up> 0)
              {
                //cout << "up exists" << endl;
                this_node = nodes_of_knickpoints[knickpoint_id + g];
                last_node = nodes_of_knickpoints[knickpoint_id + g-1];
                FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
                float temp_distance = sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                if(temp_distance < distance_to_process_down)
                {
                  distance_to_process_down-=temp_distance;
                  knickpoint_to_delete.push_back(this_node);
                } else{distance_to_process_down = 0;}
              }
              else {distance_to_process_up = 0;}
            }
            // Now I have to erase everything I planned to
            for (std::vector<int>::iterator it2 = knickpoint_to_delete.begin() ; it2 != knickpoint_to_delete.end(); it2++)
            {
              this_segment_counter_knickpoint_map.erase(*it2);
              //cout << *it2 <<endl;
              number_of_erase++;
            }

            it = this_segment_counter_knickpoint_map.begin();
            knickpoint_to_delete.clear();

            // calculate the new maximum
            old_max_knickpoint_value = max_knickpoint_value;
            max_knickpoint_value = 0;
            for(int j = 0; j< segment_counter_knickpoint; j++)
            {
              //cout << "Am I reaching this point?? j: "<< j << endl;
              this_node = nodes_of_knickpoints[j];
              if(this_segment_counter_knickpoint_map.count(this_node))
              {
                if(this_segment_counter_knickpoint_map[this_node] > max_knickpoint_value && this_segment_counter_knickpoint_map[this_node] < old_max_knickpoint_value && j!=knickpoint_id )
                  {
                    max_knickpoint_value = this_segment_counter_knickpoint_map[this_node];
                    //if(this_segment_counter_knickpoint_map[this_node] == max_knickpoint_value) {cout<< "biiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiite"<< endl;}
                  }
              }
            }
            // reset the internal loop
            still_processing = false;


          }
          else
          {

            it++;
            knickpoint_id++;
            cout <<"application of the length threshold: " << it->first << "||" << number_of_erase <<"||"<<old_max_knickpoint_value<<"||" <<max_knickpoint_value <<endl;
            if(it->first == nodes_of_knickpoints[last_knickpoint_id] || max_knickpoint_value == 0)
            {
              still_processing_total = false;
            }
          }
        }
      }
    }
    // alright, trying a new stuff for the length threshold this time











    /*


    if(false)
    {
      // now sorting and calculating the knickpoint values by length (hopefully)
      int nodes_of_knickpoints [segment_counter_knickpoint];
      map<int,vector <float> > temp_knickpoint_map;
      //map <int,float,float>
      int knickpoint_id = 0;
      for(int i = 0; i< n_nodes; i++)
      {
        this_node = node_sequence[i];
        if(this_segment_counter_knickpoint_map.count(this_node)) // supposed to test if this map has a value assigned for this key
        {
          nodes_of_knickpoints[knickpoint_id] = this_node;
          FlowInfo.get_x_and_y_from_current_node(this_node, x1_temp, y1_temp);
          float this_node_float = this_node;
          vector<float> temp_vecta (3);
          temp_vecta.push_back(this_node);
          temp_vecta.push_back(x1_temp);
          temp_vecta.push_back(y1_temp);
          temp_knickpoint_map[knickpoint_id] = temp_vecta;
          knickpoint_id++;
          cout << temp_vecta[1] << endl;

          //cout << this_node << " || " << nodes_of_knickpoints[knickpoint_id] << endl;
        }
      }
    }




    // this part begin the calculation for the length threshold to erase too dense areas and set new variables
    float distance_to_process_down = threshold_knickpoint_length/2; // represents the threshold
    float distance_to_process_up = threshold_knickpoint_length/2; // represents the threshold
    float distance_to_substract = 0; // temp variable that remove the distance already processed from the threshold
    bool still_processing = true; // tell the for loop when to reloop fromn the beginning
    bool still_processing_total = false; // tell the main loop when to stop
    int old_max = 0; // used to stock the previous maximum data
    new_knickpoint_counter = segment_counter_knickpoint; // will return the number of knickpoints after deletion of the oldest
    //number_of_nodes_to_investigate_length = threshold_knickpoint_length;
    // beginning the calculation
    cout << "beginning the length calculation stuffs" << endl;




    // old method, will erase if I found another way
    while(still_processing_total) // This will be shutted down when the last node will have been processed
    {
      still_processing = true; // Required to trigger next loop
      for(int i = 0; i< segment_counter_knickpoint && still_processing; i++)
      {
        this_node = nodes_of_knickpoints[i]; // go through the knickpoints
        if(this_segment_counter_knickpoint_map.count(this_node)) // check if this node still exists
        {
          if(  this_segment_counter_knickpoint_map[this_node]==max_knickpoint_value) // Check if it is the current maximum value
          {
            for(int g = 1; distance_to_process_down > 0 || distance_to_process_up > 0 ; g++) // If so, test the adjacent nodes in order to delete the required ones
            {
              // now getting the coordinates of the wanted nodes
              cout<<"aue coute"<<endl;
              if(this_segment_length_map.count(this_node-g) && this_segment_length_map.count(this_node-g+1))
              {

                FlowInfo.get_x_and_y_from_current_node(this_node-g, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node-g+1, x2_temp, y2_temp);


                // Check if it exists and if it is on the same river
                if(this_segment_counter_knickpoint_map.count(this_node-g) && distance_to_process_down >0)
                {
                  if(sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp))<= distance_to_process_down)
                  {
                    this_segment_counter_knickpoint_map.erase(this_node-g);
                    this_segment_knickpoint_sign_map.erase(this_node-g);
                    new_knickpoint_counter--;
                    cout << "something to test blablabla" << endl;
                    distance_to_process_down -= sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                  }
                  else
                  {
                    distance_to_process_down = 0;
                  }
                }
              }
              else
              {
                distance_to_process_down = 0;
              }
              if(this_segment_length_map.count(this_node+g) && this_segment_length_map.count(this_node+g-1))
              {
                FlowInfo.get_x_and_y_from_current_node(this_node+g, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node+g-1, x2_temp, y2_temp);
                if(this_segment_counter_knickpoint_map.count(this_node+g) && distance_to_process_up>0)
                {
                  if(sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp))<= distance_to_process_up)
                  {
                    cout << "before erase :" << this_segment_knickpoint_sign_map[this_node+g] << endl;
                    this_segment_counter_knickpoint_map.erase(this_node+g);
                    this_segment_knickpoint_sign_map.erase(this_node+g);
                    new_knickpoint_counter--;
                    cout << "after erase :" << this_segment_knickpoint_sign_map[this_node+g] << endl;
                    distance_to_process_up -= sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                  }
                  else
                  {
                    distance_to_process_up=0;
                  }
                }
              }
              else
              {
                distance_to_process_up=0;
              }
            }

            old_max = max_knickpoint_value;
            max_knickpoint_value = 0;
            for(int j = 0; j< segment_counter_knickpoint; j++)
            {
              //cout << "Am I reaching this point?? j: "<< j << endl;
              this_node = nodes_of_knickpoints[j];
              if(this_segment_counter_knickpoint_map.count(this_node))
              {
                if(this_segment_counter_knickpoint_map[this_node] > max_knickpoint_value && this_segment_counter_knickpoint_map[this_node] <= old_max)
                  {
                    max_knickpoint_value =this_segment_counter_knickpoint_map[this_node];
                    //if(this_segment_counter_knickpoint_map[this_node] == max_knickpoint_value) {cout<< "biiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiite"<< endl;}
                  }
              }
            }
            still_processing = false;
            distance_to_process_up = threshold_knickpoint_length/2;
            distance_to_process_down = threshold_knickpoint_length/2;
          }
        }
        if(i == segment_counter_knickpoint-1)
        {
          still_processing = false;
          still_processing_total = false;
        }
        else{still_processing = true;}
      }
    }*/
  }


  cout << "segment_counter_knickpoint is   " << new_knickpoint_counter << "/" << segment_counter << " delta max is " << temp_delta_m << endl;
  // print everything in the public/protected variables
  segment_counter_knickpoint_map = this_segment_counter_knickpoint_map;
  segment_knickpoint_sign_map = this_segment_knickpoint_sign_map;
  segment_length_map = this_segment_length_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the fitted elevations: It uses m_chi and b_chi
// data to get the fitted elevation of the channel points.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_segmented_elevation(LSDFlowInfo& FlowInfo)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node;
  map<int,float> this_segmented_elevation_map;
  float this_M_chi, this_b_chi, this_chi, this_segemented_elevation;

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {

      // Get the M_chi and b_chi from the current node
      this_node = node_sequence[n];
      this_M_chi = M_chi_data_map[this_node];
      this_b_chi = b_chi_data_map[this_node];
      this_chi = chi_data_map[this_node];

      // calculate elevations simply based on the fact that we are fitting segments
      // with the equation z = M_chi*chi+b_chi
      this_segemented_elevation = this_M_chi*this_chi+this_b_chi;

      // Print the segment counter to the data map
      this_segmented_elevation_map[this_node]  = this_segemented_elevation;
    }
  }
  segmented_elevation_map = this_segmented_elevation_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_full(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // test to see if there is segment numbering
  bool have_segments = false;
  if( segment_counter_map.size() == node_sequence.size())
  {
    have_segments = true;
  }

  // test to see if the fitted elevations have been calculated
  bool have_segmented_elevation = false;
  if( segmented_elevation_map.size() == node_sequence.size())
  {
    have_segmented_elevation = true;
  }


  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,chi,elevation,flow distance,drainage area,m_chi,b_chi,source_key,basin_key";
  if(have_segmented_elevation)
  {
    chi_data_out << ",segmented_elevation";
  }
  if (have_segments)
  {
    chi_data_out << ",segment_number";
  }
  chi_data_out << endl;




  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << ","
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      if(have_segmented_elevation)
      {
        chi_data_out << "," << segmented_elevation_map[this_node];
      }
      if (have_segments)
      {
        chi_data_out << "," << segment_counter_map[this_node];
      }
      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Development function to Print data maps to file including knickpoints
// BG
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_full_knickpoints(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // test to see if there is segment numbering
  bool have_segments = false;
  if( segment_counter_map.size() == node_sequence.size())
  {
    have_segments = true;
  }

  // test to see if the fitted elevations have been calculated
  bool have_segmented_elevation = false;
  if( segmented_elevation_map.size() == node_sequence.size())
  {
    have_segmented_elevation = true;
  }


  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,chi,elevation,flow distance,drainage area,m_chi,b_chi,source_key,basin_key";
  if(have_segmented_elevation)
  {
    chi_data_out << ",segmented_elevation";
  }
  if (have_segments)
  {
    chi_data_out << ",segment_number";
  }
  chi_data_out << ",knickpoints,knickpoint_sign,segment_length"; // add the knickpoint col
  chi_data_out << endl;

  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << ","
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      if(have_segmented_elevation)
      {
        chi_data_out << "," << segmented_elevation_map[this_node];
      }
      if (have_segments)
      {
        chi_data_out << "," << segment_counter_map[this_node];
      }

      chi_data_out << "," << segment_counter_knickpoint_map[this_node];
      chi_data_out << "," << segment_knickpoint_sign_map[this_node];
      chi_data_out << "," << segment_length_map[this_node];
      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_source_keys(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col, key;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  map<int, int>::iterator it;

  // open the data file
  ofstream  source_keys_out;
  source_keys_out.open(filename.c_str());
  source_keys_out << "latitude,longitude,source_node,source_key" << endl;

  // loop through the source key map
  for ( it = key_to_source_map.begin(); it != key_to_source_map.end(); it++ )
  {
    key = it->second;
    this_node = it->first;
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    source_keys_out.precision(9);
    source_keys_out << latitude << ","
                   << longitude << "," << this_node << ",";
    source_keys_out.precision(5);
    source_keys_out << key << endl;
  }

  source_keys_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_baselevel_keys(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, this_junc,row,col, key;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  map<int, int>::iterator it;

  // open the data file
  ofstream  baselevel_keys_out;
  baselevel_keys_out.open(filename.c_str());
  baselevel_keys_out << "latitude,longitude,baselevel_node,baselevel_junction,baselevel_key" << endl;

  // loop through the source
  for ( it = key_to_baselevel_map.begin(); it != key_to_baselevel_map.end(); it++ )
  {
    key = it->second;
    this_node = it->first;
    this_junc = JunctionNetwork.get_Junction_of_Node(this_node, FlowInfo);
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    baselevel_keys_out.precision(9);
    baselevel_keys_out << latitude << ","
                   << longitude << "," << this_node << ",";
    baselevel_keys_out.precision(5);
    baselevel_keys_out << this_junc <<"," << key << endl;
  }

  baselevel_keys_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints the basins and an additional file that has basin centroids
// and labelling information for plotting
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_basins(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Junctions, string base_filename)
{
  int N_Juncs = Junctions.size();
  LSDCoordinateConverterLLandUTM Converter;


  // Get some data members for holding basins and the raster
  vector<LSDBasin> AllTheBasins;
  map<int,int> drainage_of_other_basins;
  LSDIndexRaster BasinMasterRaster;

  string basin_raster_name = base_filename+"_AllBasins";
  string basin_info_name = base_filename+"_AllBasinsInfo.csv";

  ofstream basin_info_out;
  basin_info_out.open(basin_info_name.c_str());
  basin_info_out << "latitude,longitude,outlet_latitude,outlet_longitude,outlet_junction" << endl;

  // Make sure the full lat-long information is printed
  basin_info_out.precision(9);

  // These store row and column information for converting the outlet and centroid to
  // latitude and longitude
  int centroid_i, centroid_j, outlet_i, outlet_j;
  double out_lat,out_long, cen_lat, cen_long;

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


    // get the centroid and outlet locations
    centroid_i = thisBasin.get_Centroid_i();
    centroid_j = thisBasin.get_Centroid_j();

    outlet_i = thisBasin.get_Outlet_i();
    outlet_j = thisBasin.get_Outlet_j();

    // Find the latitude and longitude of the outlet and centroid
    get_lat_and_long_locations(centroid_i, centroid_j, cen_lat, cen_long, Converter);
    get_lat_and_long_locations(outlet_i, outlet_j, out_lat, out_long, Converter);

    basin_info_out << cen_lat << "," << cen_long << "," << out_lat << "," << out_long << "," << Junctions[BN] << endl;
  }
  basin_info_out.close();

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





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_basic(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,m_chi,b_chi" << endl;

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(6);
      chi_data_out << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << "," << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




#endif
