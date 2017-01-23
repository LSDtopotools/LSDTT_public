//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// basin_metrics_driver_resubmission_version.cpp
// David Milodowski July 2014
// Analysis for Ecology paper
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../TNT/tnt.h"
int main (int nNumberofArgs,char *argv[])
{
  // ADD ALL PATH RELATED STUFF HERE!
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


	string raster_ID;
	string data_path;
	string dem_ext = "_dem";
	string ch_ext = "_CH_dreich";
  string flt_ext = "flt";

  string temp;                
	int threshold;
	float Minimum_Slope;
	float A_0;
	float m_over_n;
	int no_connecting_nodes;
	float curv_threshold;
	
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the header file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}	
  file_info_in >> temp >> data_path;                                     
	file_info_in >> temp >> raster_ID;
	file_info_in >> temp >> threshold;
	file_info_in >> temp >> Minimum_Slope;
	file_info_in >> temp >> A_0;
	file_info_in >> temp >> m_over_n;
	file_info_in >> temp >> no_connecting_nodes;
	file_info_in >> temp >> curv_threshold;
	file_info_in.close();

	LSDRaster dem(data_path+raster_ID+dem_ext, flt_ext);
	
  int NRows = dem.get_NRows();
	int NCols = dem.get_NCols();
	float NoDataValue = dem.get_NoDataValue();
  
  cout << "filling DEM" << endl;
  LSDRaster filled_dem = dem.fill(Minimum_Slope);
//   filled_dem.write_raster("test","flt");
  cout << "creating FlowInfo" << endl;
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
  
	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_dem);
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	
		//get the sources: note: this is only to select basins!
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

  float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
  vector<LSDRaster> surface_fitting;
  LSDRaster tan_curvature;
  string curv_ext = "_tan_curv";
  vector<int> raster_selection(8, 0);
  raster_selection[6] = 1;                      // this indicates you only want the tangential curvature
  surface_fitting = filled_dem.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);

  // now print the tangential curvature raster to file.
  for(int i = 0; i<int(raster_selection.size()); ++i)
	{
		if(raster_selection[i]==1)
		{
      tan_curvature = surface_fitting[i];
      tan_curvature.write_raster((path_name+raster_ID+curv_ext), flt_ext);
    }
  }
  // Find the valley junctions
  Array2D<float> tan_curv_array = tan_curvature.get_RasterData();
  cout << "got tan curvature array" << endl;
	Array2D<int> valley_junctions = ChanNetwork.find_valleys(FlowInfo, tan_curv_array, sources, no_connecting_nodes);
	cout << "Got valley junctions, proceeding to chi analysis" << endl;

	// Calculate the channel head nodes
  int MinSegLength = 10;
  vector<int> ChannelHeadNodes = ChanNetwork.GetChannelHeadsChiMethodFromValleys(valley_junctions, MinSegLength,
                                                 A_0, m_over_n, FlowInfo, DistanceFromOutlet, filled_dem);

  //write channel_heads to a csv file
  FlowInfo.print_vector_of_nodeindices_to_csv_file(ChannelHeadNodes, data_path+raster_ID+ch_ext);

	//create a channel network based on these channel heads
	LSDJunctionNetwork NewChanNetwork(ChannelHeadNodes, FlowInfo);
  LSDIndexRaster SOArrayNew = NewChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name_new = "_SO_from_CH";


}
