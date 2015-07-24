//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// channel_isolation_wiener_qq.cpp
// make with channel_isolation_wiener_qq.make
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// David T Milodowski
// Stuart W.D Grieve
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDRasterSpectral.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../TNT/tnt.h"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

  string full_name = path_name+f_name;

  ifstream file_info_in;
  file_info_in.open(full_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  string Raster_name;
  string Output_name;
  string q_q_filename_prefix;
  float area_threshold,window_radius;
  string DEM_extension = "flt";
  string temp;
  int connected_components_threshold;
  file_info_in >> temp >> Raster_name
               >> temp >> Output_name
               >> temp >> q_q_filename_prefix 
               >> temp >> window_radius
               >> temp >> area_threshold
	       >> temp >> connected_components_threshold;
  file_info_in.close();
  // Now create the raster selection vector based on user's selection
  // Elevation
  LSDRasterSpectral raster(Raster_name, DEM_extension);
  LSDIndexRaster output_raster = raster.IsolateChannelsWienerQQ(area_threshold, window_radius, q_q_filename_prefix+".txt");
  cout << "filter by connected components" << endl;
  //LSDIndexRaster output_raster(Output_name,DEM_extension);
  output_raster = output_raster.filter_by_connected_components(connected_components_threshold);
  //LSDIndexRaster output_raster(Output_name,DEM_extension);
  cout << "thin network to skeleton" << endl;
  LSDIndexRaster skeleton_raster = output_raster.thin_to_skeleton();
  cout << "finding end points" << endl;
  LSDIndexRaster Ends = skeleton_raster.find_end_points();
  //output_raster = output_raster.ConnectedComponents();
  output_raster.write_raster(Output_name+"_cc", DEM_extension);
  skeleton_raster.write_raster(Output_name+"_skeleton",DEM_extension);
  Ends.write_raster(Output_name+"_end_points",DEM_extension);
  
  
  //Now we can process the end points to get only the channel heads - SWDG
  
  cout << "Starting channel head processing" << endl;
  
  //First we need to load the elevation data, fill it and generate a FlowInfo object
  LSDRaster DEM(Raster_name, DEM_extension);
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);
  vector<string> BoundaryConditions(4, "No Flux");
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);
  
  //this processes the end points to only keep the upper extent of the channel network
  cout << "getting channel heads" << endl;
  vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);
      
  // we need a temp junction network to search for single pixel channels
  LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
  LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  
  cout << "removing single px channels" << endl;
  vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);
  
  //Now we have the final channel heads, so we can generate a channel network from them
  LSDJunctionNetwork JunctionNetwork(FinalSources, FlowInfo);
  LSDIndexRaster StreamNetwork = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  
  //Finally we write the channel network and channel heads to files so they can be used in other drivers. 
  LSDIndexRaster Heads = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(FinalSources);
  Heads.write_raster((Output_name+"_CH"),DEM_extension);
  StreamNetwork.write_raster((Output_name+"_st_final"),DEM_extension);
  
  cout << "DONE" << endl;
}
