//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// get_floodplain_mask
//
// This program takes two arguments, the path name and the driver name
// It gets a mask of the floodplain from a DEM. User has the option to filter the DEM
// using a Perona-Malik filter before getting the floodplain. After the inital floodplain mask
// is found, patches identified that are not connected to the channel network are removed. 
// Small holes are then removed from the floodplain mask.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 19/10/16
// Fiona J. Clubb
// University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../../LSDStatsTools.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDFloodplain.hpp"


int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();
	
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

	string DEM_ID;
	string RASTER_NAME;
	string CH_name;
	string input_path;
  string dem_ext = "_dem";
  //string sources_ext = "_CH_wiener";
  string DEM_extension = "bil";
  string csv_extension = "csv";
  string txt_extension = ".txt";
  
  // initialise variables to be assigned from .driver file
  int threshold_SO, FilterTopo, window_radius, lower_percentile, upper_percentile, minimum_patch_size, junction_number, search_distance;
	float Minimum_Slope, surface_fitting_window_radius, threshold_condition;
  string temp;
  
  // read in the parameters                                                          
	file_info_in >> temp >> DEM_ID
               >> temp >> RASTER_NAME
							 >> temp >> CH_name
               >> temp >> input_path
               >> temp >> Minimum_Slope
               >> temp >> threshold_SO
               >> temp >> FilterTopo
               >> temp >> surface_fitting_window_radius
               >> temp >> window_radius
							 >> temp >> threshold_condition
							 >> temp >> lower_percentile
		           >> temp >> upper_percentile
							 >> temp >> minimum_patch_size
							 >> temp >> search_distance
							 >> temp >> junction_number;
                   
	file_info_in.close();

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	LSDRaster filled_topo_test;
   
  if(FilterTopo == 1)
  {
     // load the DEM
		 cout << "Loading the DEM..." << endl;
     LSDRaster topo_test((input_path+DEM_ID), DEM_extension);   
     
     // filter using Perona Malik
     int timesteps = 50;
     float percentile_for_lambda = 90;
     float dt = 0.1;
     topo_test = topo_test.PeronaMalikFilter(timesteps, percentile_for_lambda, dt); 
     
     // fill
     filled_topo_test = topo_test.fill(Minimum_Slope);
     string fill_name = "_filtered";
     filled_topo_test.write_raster((path_name+DEM_ID+fill_name), DEM_extension);   
  }
  else
  {
    //previously done the filtering and filling, just load the filled DEM
    LSDRaster load_DEM((input_path+DEM_ID+"_filtered"), DEM_extension);
    filled_topo_test = load_DEM;
  }
  

  cout << "\t Flow routing..." << endl;
	// get a flow info object
 	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
	
  // calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	cout << "\t Loading Sources..." << endl;
	// load the sources
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_ID+CH_name), csv_extension, 2);
  cout << "\t Got sources!" << endl;
	
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;
	
	//print out the junction network
	string JI_name = "_JI";
  LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
  JIArray.write_raster((path_name+DEM_ID+JI_name), DEM_extension);
	
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name = "_SO";
	SOArray.write_raster((path_name+DEM_ID+SO_name), DEM_extension);
    
   
  //calculate the channel relief
  cout << "\t Getting relief relative to channel" << endl;
  cout << "\t Threshold stream order = " << threshold_SO << endl;
  LSDRaster ChannelRelief = ChanNetwork.calculate_relief_from_channel(filled_topo_test, FlowInfo, threshold_SO);
  string relief_name = "_channel_relief";
  ChannelRelief.write_raster((input_path+DEM_ID+relief_name), DEM_extension);
  cout << "\t Got the relief!" << endl;
     
  //get the slope
  cout << "\t Calculating slope..." << endl;
  vector<LSDRaster> surface_fitting;
  LSDRaster Slope;
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;             // this means you want the slope
  surface_fitting = filled_topo_test.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
  Slope = surface_fitting[1];
  cout << "\t Done!" << endl;
  string slope_name = "_slope";
  Slope.write_raster((input_path+DEM_ID+slope_name), DEM_extension);

  // get the channel relief and slope threshold using quantile-quantile plots
  cout << "Getting channel relief threshold from QQ plots" << endl;
  string qq_fname = path_name+DEM_ID+"_qq_relief.txt";
  float relief_threshold_from_qq = ChannelRelief.get_threshold_for_floodplain_QQ(qq_fname, threshold_condition, lower_percentile, upper_percentile);
  
  cout << "Getting slope threshold from QQ plots" << endl;
  string qq_slope = path_name+DEM_ID+"_qq_slope.txt";
  float slope_threshold_from_qq = Slope.get_threshold_for_floodplain_QQ(qq_slope, threshold_condition, lower_percentile, upper_percentile);
	
	cout << "Relief threshold: " << relief_threshold_from_qq << " Slope threshold: " << slope_threshold_from_qq << endl;
	
	// get the distance from outlet
	LSDRaster DistFromOutlet = FlowInfo.distance_from_outlet();
	
	// TESTING FLOODPLAIN OBJECT
	
	// get the floodplain object
	LSDFloodplain Floodplain(ChannelRelief, Slope, ChanNetwork, FlowInfo, relief_threshold_from_qq, slope_threshold_from_qq, minimum_patch_size, threshold_SO);
	
	// get the relief relative to nearest channel
	Floodplain.Get_Relief_of_Nearest_Channel(ChanNetwork, FlowInfo, filled_topo_test, DistFromOutlet, threshold_SO, search_distance);
	LSDRaster relief_final = Floodplain.print_ChannelRelief_to_Raster();
	string relief_ext = "_relief_final";
	relief_final.write_raster((input_path+DEM_ID+relief_ext), DEM_extension);
	
	LSDIndexRaster BinaryRaster = Floodplain.print_BinaryRaster();
	string bin_ext = "_FP";
	BinaryRaster.write_raster((input_path+DEM_ID+bin_ext), DEM_extension);
	
//	cout << "This junction number is: " << junction_number << endl;
//	Floodplain.get_distance_upstream_along_main_stem(junction_number, ChanNetwork, FlowInfo, DistFromOutlet);
//	LSDRaster UpstreamDistance = Floodplain.print_UpstreamDistance_to_Raster();
//	string dist_ext = "_upstream_dist";
//	UpstreamDistance.write_raster((input_path+DEM_ID+dist_ext), DEM_extension);
//	
//	// write to text file
//	string filename = "_floodplain_data.txt";
//	Floodplain.print_ChannelRelief_to_File(input_path+DEM_ID+filename);
//	
//	string filename_binned = "_floodplain_data_binned.txt";
//	float bin_width = 50;
//	float bin_lower_limit = 0;
//	float bin_threshold = 0;
//	Floodplain.print_Binned_ChannelRelief_to_File(input_path+DEM_ID+filename_binned, bin_width, bin_lower_limit, bin_threshold);
		
	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
