//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// channel_extraction_pelletier.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program calculates channel heads using Pelletier (2013)
//
// Reference: Pelletier, J.D. (2013) A robust, two-parameter method for the extraction of
// drainage networks from high-resolution digital elevation models (DEMs): Evaluation using
// synthetic and real-world DEMs, Water Resources Research 49(1): 75-89, doi:10.1029/2012WR012452
//
// Developed by:
//  Fiona Clubb
//  Simon M. Mudd
//  David T. Milodowski
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp" 
#include "../LSDRasterInfo.hpp"
#include "../LSDParameterParser.hpp"

int main (int nNumberofArgs,char *argv[])
{
  //start the clock
  clock_t begin = clock();

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the channel extraction tool!             ||" << endl;
    cout << "|| This program has a number of options to extract     ||" << endl;
    cout << "|| channel networks.                                   ||" << endl;
    cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
    cout << "|| Fiona J. Clubb and Simon M. Mudd"                   ||
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./channel_extraction_tool.exe /LSDTopoTools/Topographic_projects/test_data channel_analysis.param" << endl;
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

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // set default float parameters
  int_default_map["threshold_contributing_pixels"] = 100;
  int_default_map["connected_components_threshold"] = 100;
  
  
  int_default_map["n_iterations"] = 20;
  int_default_map["minimum_segment_length"] = 10;
  int_default_map["n_nodes_to_visit"] = 10;
  int_default_map["target_nodes"] = 80;
  int_default_map["skip"] = 2;
  int_default_map["threshold_pixels_for_chi"] = 1000;
  int_default_map["minimum_basin_size_pixels"] = 1000;
  
  
  // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;
  float_default_map["surface_fitting_radius"] = 6;
  float_default_map["pruning_drainage_area"] = 1000;


  float_default_map["A_0"] = 1;
  float_default_map["m_over_n"] = 0.5;
  float_default_map["sigma"] = 20;
  
  // set default methods
  bool_default_map["print_area_threshold_channels"] = true;
  bool_default_map["print_dreich_channels"] = false;
  bool_default_map["print_pelletier_channels"] = false;
  bool_default_map["print_wiener_channels"] = false;
  
  bool_default_map["print_stream_order_raster"] = true;
  bool_default_map["print_sources_to_raster"] = false;
  bool_default_map["print_fill_raster"] = false;
  bool_default_map["print_DrainageArea_raster"] = false;
  bool_default_map["write hillshade"] = false;
  
  bool_default_map["print_sources_to_csv"] = true;
  bool_default_map["print_channels_to_csv"] = true;
  

  // set default string method
  string_default_map["CHeads_file"] = "NULL";
  
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
  
  cout << "Read filename is: " <<  DATA_DIR+DEM_ID << endl;
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;
  
  // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);
  
  // load the base raster
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;
  
  
  //============================================================================
  // Start gathering necessary rasters
  //============================================================================
  // Get the fill
  cout << "Filling topography." << endl;
  LSDRaster filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
  
  if (this_bool_map["print_fill_raster"])
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
  
  // get the flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);
  
  //===============================================================
  // AREA THRESHOLD
  //===============================================================
  if (this_bool_map["print_area_threshold_channels"])
  {
    cout << "I am calculating channels using an area threshold." << endl;
    cout << "Only use this if you aren't that bothered about where the channel heads actually are!" << endl;
    
    // get some relevant rasters
    LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
    LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
    
    //get the sources: note: this is only to select basins!
    vector<int> sources;
    sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["threshold_contributing_pixels"]);
  
    // now get the junction network
    LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
    
    
    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
      string sources_csv_name = OUT_DIR+OUT_ID+"_ATsources.csv";
      
      //write channel_heads to a csv file
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(sources, sources_csv_name);
    }

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_ATsources";
      
      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(sources);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_AT_SO";
      
      //write stream order array to a raster
      LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }
  
    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_AT_CN";
      ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
    }
  }

  //===============================================================
  // DREICH
  //===============================================================
  if (this_bool_map["print_dreich_channels"])
  {
    cout << "I am calculating channels using the dreich algorighm (DOI: 10.1002/2013WR015167)." << endl;
  
  }

  //===============================================================
  // PELLETIER
  //===============================================================
  if (this_bool_map["print_pelletier_channels"])
  {
    cout << "I am calculating channels using the pelletier algorighm (doi:10.1029/2012WR012452)." << endl;
  }
  
  //===============================================================
  // WIENER
  //===============================================================
  if (this_bool_map["print_weiner_channels"])
  {
    cout << "I am calculating channels using the weiner algorithm (doi:10.1029/2012WR012452)." << endl;
    cout << "This algorithm was used by Clubb et al. (2016, DOI: 10.1002/2015JF003747) " << endl;
    cout << "and Grieve et al. (2016, doi:10.5194/esurf-4-627-2016) " << endl;
    cout << "and combines elements of the Pelletier and Passalacqua et al  methods: " << endl;
    cout << " doi:10.1029/2012WR012452 and doi:10.1029/2009JF001254" << endl;

    // initiate the spectral raster
    LSDRasterSpectral Spec_raster(topography_raster);
    
    string QQ_fname = OUT_DIR+OUT_ID+"__qq.txt";
    
    cout << "I am am getting the connected components using a weiner QQ filter." << endl;
    LSDIndexRaster connected_components = Spec_raster.IsolateChannelsWienerQQ(this_float_map["pruning_area"], 
                                                       this_float_map["surface_fitting_radius"], QQ_fname);
    
    cout << "I am filtering by connected components" << endl;
    LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(this_int_map["connected_components_threshold"]);
    LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();

    cout << "I am thinning the network to a skeleton." << endl;
    LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();
    
    cout << "I am finding the finding end points" << endl;
    LSDIndexRaster Ends = skeleton_raster.find_end_points();
    Ends.remove_downstream_endpoints(CC_raster, Spec_raster);

    //this processes the end points to only keep the upper extent of the channel network
    cout << "getting channel heads" << endl;
    vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);

    // we need a temp junction network to search for single pixel channels
    LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
    LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

    cout << "removing single px channels" << endl;
    vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);

    //Now we have the final channel heads, so we can generate a channel network from them
    LSDJunctionNetwork ChanNetwork(FinalSources, FlowInfo);

    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
      string sources_csv_name = OUT_DIR+OUT_ID+"_Wsources.csv";
      
      //write channel_heads to a csv file
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(FinalSources, sources_csv_name);
    }

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_Wsources";
      
      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(FinalSources);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_W_SO";
      
      //write stream order array to a raster
      LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }
  
    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_W_CN";
      ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
    }

  }  
  



  // Done, check how long it took
  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  cout << "DONE! Have fun with your channels! Time taken (secs): " << elapsed_secs << endl;

}
