//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// channel_heads_wiener_qq_driver.cpp
// make with channel_heads_wiener_qq.make
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// David T Milodowski
// Stuart W.D Grieve
// Fiona J. Clubb
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
  float area_threshold,window_radius, A_0, m_over_n;
  string DEM_extension = "flt";
  string temp;
  int connected_components_threshold;
  file_info_in >> temp >> Raster_name
               >> temp >> Output_name
               >> temp >> q_q_filename_prefix 
               >> temp >> window_radius
               >> temp >> area_threshold
	             >> temp >> connected_components_threshold
	             >> temp >> A_0
	             >> temp >> m_over_n;
  file_info_in.close();
  // Now create the raster selection vector based on user's selection
  // Elevation
  LSDRasterSpectral raster(Raster_name, DEM_extension);
  LSDIndexRaster connected_components = raster.IsolateChannelsWienerQQ(area_threshold, window_radius, q_q_filename_prefix+".txt");
  cout << "filter by connected components" << endl;
  //LSDIndexRaster output_raster(Output_name,DEM_extension);
  LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(connected_components_threshold);
  LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();
  //LSDIndexRaster output_raster(Output_name,DEM_extension);
  cout << "thin network to skeleton" << endl;
  LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();
  cout << "finding end points" << endl;
  LSDIndexRaster Ends = skeleton_raster.find_end_points();
  Ends.remove_downstream_endpoints(CC_raster, raster);

  //write some rasters
  connected_components_filtered.write_raster(Output_name+"_cc", DEM_extension);
  skeleton_raster.write_raster(Output_name+"_skeleton",DEM_extension);
  Ends.write_raster(Output_name+"_end_points",DEM_extension);
  
  
  //Now we can process the end points to get only the channel heads - SWDG
  
  cout << "Starting channel head processing" << endl;
  
  //First we need to load the elevation data, fill it and generate a FlowInfo object
  LSDRaster DEM(Raster_name, DEM_extension);
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);
  string fill_name = "_fill";
  FilledDEM.write_raster((Output_name+fill_name), DEM_extension);
  
  vector<string> BoundaryConditions(4, "No Flux");
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);
  
  string HS_name = "_HS";
	LSDRaster HS = FilledDEM.hillshade(45, 315, 1);
	HS.write_raster((Output_name+HS_name),DEM_extension);
  
  //this processes the end points to only keep the upper extent of the channel network
  cout << "getting channel heads" << endl;
  vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);
  cout << "processed all end points" << endl;
      
  // we need a temp junction network to search for single pixel channels
  LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
  LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  
  cout << "removing single px channels" << endl;
  vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);
  
  // using these sources as the input to run the DrEICH algorithm  - FJC
   
  //Generate a channel netowrk from the sources
  LSDJunctionNetwork JunctionNetwork(FinalSources, FlowInfo);
	LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
	//string JN_name = "_JI";
	//JIArray.write_raster(Output_name+JN_name, DEM_extension);
  
  LSDIndexRaster StreamNetwork = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  string SO_name = "_SO";
  StreamNetwork.write_raster(Output_name+SO_name, DEM_extension);
    
  //Get the outlet junctions of each of the valleys
  vector<int> valley_nodes = JunctionNetwork.get_outlet_nodes_from_sources(FlowInfo, FinalSources);
  cout << "Got valley nodes, proceeding to chi analysis" << endl;
  
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
  
	// Calculate the channel head nodes
  int MinSegLength = 30;
  
  vector<int> ChannelHeadNodes = JunctionNetwork.GetChannelHeadsChiMethodFromValleys(valley_nodes, MinSegLength, A_0, m_over_n,
									                    FlowInfo, DistanceFromOutlet, FilledDEM);
                                                 
  //write channel heads to a raster
  string CH_name = "_CH_DrEICH";
  LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(ChannelHeadNodes);
  Channel_heads_raster.write_raster((Output_name+CH_name),DEM_extension);

  //create a channel network based on these channel heads
  LSDJunctionNetwork NewChanNetwork(ChannelHeadNodes, FlowInfo);
  //int n_junctions = NewChanNetwork.get_Number_of_Junctions();
  LSDIndexRaster SOArrayNew = NewChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  string SO_name_new = "_SO_from_CH_DrEICH_test";

  SOArrayNew.write_raster((Output_name+SO_name_new),DEM_extension);                                               
}
