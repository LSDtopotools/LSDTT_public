//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// extract_terraces_driver
//
// This program takes two arguments, the path name and the driver name
// It extracts all the terraces from a DEM. User has the option to filter the DEM
// using a Perona-Malik filter before selecting the terraces.
// For each terrace it calculates the relief relative to the nearest channel.  User 
// must specify the junction number of the main stem channel - creates a text file for an XY
// plot of terrace elevations vs distance upstream.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 28/10/16
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
#include "../../LSDTerrace.hpp"


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
  int threshold_SO, FilterTopo, window_radius, lower_percentile_relief, upper_percentile_relief, lower_percentile_slope, upper_percentile_slope, minimum_patch_size, junction_number, search_distance, remove_channel_patches, MainStemAnalysis;
	float Minimum_Slope, surface_fitting_window_radius, threshold_condition, bin_width;
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
							 >> temp >> lower_percentile_relief
		           >> temp >> upper_percentile_relief
							 >> temp >> lower_percentile_slope
							 >> temp >> upper_percentile_slope
							 >> temp >> minimum_patch_size
							 >> temp >> search_distance
						   >> temp >> MainStemAnalysis
							 >> temp >> junction_number
							 >> temp >> bin_width
							 >> temp >> remove_channel_patches;
                   
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
  float relief_threshold_from_qq = ChannelRelief.get_threshold_for_floodplain_QQ(qq_fname, threshold_condition, lower_percentile_relief, upper_percentile_relief);
  
  cout << "Getting slope threshold from QQ plots" << endl;
  string qq_slope = path_name+DEM_ID+"_qq_slope.txt";
  float slope_threshold_from_qq = Slope.get_threshold_for_floodplain_QQ(qq_slope, threshold_condition, lower_percentile_slope, upper_percentile_slope);
	
	cout << "Relief threshold: " << relief_threshold_from_qq << " Slope threshold: " << slope_threshold_from_qq << endl;
	
	// get the distance from outlet
	LSDRaster DistFromOutlet = FlowInfo.distance_from_outlet();
	
	// TESTING TERRACE OBJECT
	
	//float relief_thresh = 100;
	//float slope_thresh = 0.1;
	// get the terrace object
	cout << "Remove channels? " << remove_channel_patches << endl;
	LSDTerrace Terraces(ChannelRelief, Slope, ChanNetwork, FlowInfo, relief_threshold_from_qq, slope_threshold_from_qq, minimum_patch_size, threshold_SO, remove_channel_patches);
	
	LSDIndexRaster TerraceLocations = Terraces.print_ConnectedComponents_to_Raster();
	string CC_ext = "_CC";
	TerraceLocations.write_raster((input_path+DEM_ID+CC_ext), DEM_extension);
	
	// get the relief relative to nearest channel
	Terraces.Get_Relief_of_Nearest_Channel(ChanNetwork, FlowInfo, filled_topo_test, DistFromOutlet, threshold_SO, search_distance);
	LSDRaster relief_final = Terraces.print_ChannelRelief_to_Raster();
	string relief_ext = "_terrace_relief_final";
	relief_final.write_raster((input_path+DEM_ID+relief_ext), DEM_extension);
	
	if (MainStemAnalysis == 1)
	{
		cout << "This junction number is: " << junction_number << endl;
		Terraces.get_terraces_along_main_stem(junction_number, ChanNetwork, FlowInfo, DistFromOutlet);
		LSDRaster UpstreamDistance = Terraces.print_UpstreamDistance_to_Raster();
		string dist_ext = "_upstream_dist";
		UpstreamDistance.write_raster((input_path+DEM_ID+dist_ext), DEM_extension);

		// print main stem relief and distance

		LSDRaster MainStemRelief = Terraces.print_ChannelRelief_to_Raster_MainStem();
		string ms_relief_ext = "_relief_MS";
		MainStemRelief.write_raster((input_path+DEM_ID+ms_relief_ext), DEM_extension);

		LSDRaster MainStemDist = Terraces.print_UpstreamDistance_to_Raster_MainStem();
		string ms_dist_ext = "_dist_MS";
		MainStemDist.write_raster((input_path+DEM_ID+ms_dist_ext), DEM_extension);

		// write to text file
		string filename = "_terraces_data.txt";
		Terraces.print_ChannelRelief_to_File(input_path+DEM_ID+filename);

		string filename_binned = "_terraces_data_binned.txt";
		float bin_lower_limit = 0;
		float bin_threshold = 0;
		Terraces.print_Binned_ChannelRelief_to_File(input_path+DEM_ID+filename_binned, bin_width, bin_lower_limit, bin_threshold);
	}
		
	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
