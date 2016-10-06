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
// 26/10/15
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
  string flt_extension = "bil";
  string csv_extension = "csv";
  string txt_extension = ".txt";
  
  // initialise variables to be assigned from .driver file
  int threshold_SO, FilterTopo, window_radius, lower_percentile, upper_percentile, minimum_patch_size, search_distance, junction_number;
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
     LSDRaster topo_test((input_path+DEM_ID), flt_extension);   
     
     // filter using Perona Malik
     int timesteps = 50;
     float percentile_for_lambda = 90;
     float dt = 0.1;
     topo_test = topo_test.PeronaMalikFilter(timesteps, percentile_for_lambda, dt); 
     
     // fill
     filled_topo_test = topo_test.fill(Minimum_Slope);
     string fill_name = "_filtered";
     filled_topo_test.write_raster((path_name+DEM_ID+fill_name), flt_extension);   
  }
  else
  {
    //previously done the filtering and filling, just load the filled DEM
    LSDRaster load_DEM((input_path+DEM_ID+"_filtered"), flt_extension);
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
	
//	cout << "N sources: " << sources.size() << endl;
//	for (i =0; i < int(sources.size(); i++)
//	{
//		cout << sources[i] << endl;
//  }
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;
	
	//print out the junction network
	string JI_name = "_JI";
  LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
  JIArray.write_raster((path_name+DEM_ID+JI_name), flt_extension);
	
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name = "_SO";
	SOArray.write_raster((path_name+DEM_ID+SO_name), flt_extension);
    
   
  //calculate the channel relief
  cout << "\t Getting relief relative to channel" << endl;
  cout << "\t Threshold stream order = " << threshold_SO << endl;
  LSDRaster ChannelRelief = ChanNetwork.calculate_relief_from_channel(filled_topo_test, FlowInfo, threshold_SO);
  string relief_name = "_channel_relief";
  ChannelRelief.write_raster((input_path+DEM_ID+relief_name), flt_extension);
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
  Slope.write_raster((input_path+DEM_ID+slope_name), flt_extension);
    
  //get the drainage area using Dinf
  //Array2D<float> dinf = filled_topo_test.D_inf_FlowDir();
  //LSDRaster dinf_rast = filled_topo_test.LSDRasterTemplate(dinf);
  //LSDRaster DinfArea = filled_topo_test.D_inf_units();
  
  //get the topographic index raster
  //LSDRaster TopoIndex = filled_topo_test.calculate_topographic_index(DinfArea, Slope);
  //string TopoIndex_name = "_TI_Dinf";
  //TopoIndex.write_raster((input_path+DEM_ID+TopoIndex_name), flt_extension);

  // get the channel relief and slope threshold using quantile-quantile plots
  cout << "Getting channel relief threshold from QQ plots" << endl;
  string qq_fname = path_name+DEM_ID+"_qq_relief.txt";
  float relief_threshold_from_qq = ChannelRelief.get_threshold_for_floodplain_QQ(qq_fname, threshold_condition, lower_percentile, upper_percentile);
  
  cout << "Getting slope threshold from QQ plots" << endl;
  string qq_slope = path_name+DEM_ID+"_qq_slope.txt";
  float slope_threshold_from_qq = Slope.get_threshold_for_floodplain_QQ(qq_slope, threshold_condition, lower_percentile, upper_percentile);
	
	cout << "Relief threshold: " << relief_threshold_from_qq << " Slope threshold: " << slope_threshold_from_qq << endl;

  //get the potential floodplain mask
  cout << "\t Getting the floodplain mask" << endl;
  LSDIndexRaster FloodplainRaster_temp = filled_topo_test.get_potential_floodplain_patches(ChannelRelief, Slope, relief_threshold_from_qq, slope_threshold_from_qq);
  
  cout << "\t Connected components" << endl;
  LSDIndexRaster ConnectedComponents = FloodplainRaster_temp.ConnectedComponents();

//  
//  //remove holes in the connected components raster
//  cout << "\t Removing holes" << endl;
//  LSDIndexRaster ConnectedComponents_final = ChannelPatches.remove_holes_in_patches_connected_components(window_radius);
  //string CC_new_name = "_CC_no_holes";
  //ConnectedComponents_final.write_raster((input_path+DEM_ID+CC_new_name), flt_extension); 
	
	// remove patches smaller than a certain number of pixels
	LSDIndexRaster ConnectedComponents_final = ConnectedComponents.RemoveSmallPatches(minimum_patch_size);
	string CC_name = "_CC_filt";
  ConnectedComponents_final.write_raster((input_path+DEM_ID+CC_name), flt_extension); 
	
	//remove patches of identified floodplain that are not connected to the channel network
//  cout << "\t Separating into floodplain and terraces" << endl;
//  LSDIndexRaster ChannelPatches, TerracePatches;
//	ChanNetwork.separate_floodplain_and_terrace_patches(ConnectedComponents_final, ChannelPatches, TerracePatches, threshold_SO);
//	
//  //get a binary raster of floodplain pixels
//  int value = 1;
//  int ndv = -9999;
//  LSDIndexRaster FloodplainMask = ChannelPatches.ConvertToBinary(value, ndv);
//  string mask_name = "_FP";
//  FloodplainMask.write_raster((input_path+DEM_ID+mask_name), flt_extension); 
//	
//	cout << "\t Getting channel relief for the floodplain" << endl;
//	LSDRaster FP_Relief_Masked = ChannelRelief.ExtractByMask(FloodplainMask);
//	string FP_relief_masked = "_FP_relief_masked";
//	FP_Relief_Masked.write_raster((input_path+DEM_ID+FP_relief_masked), flt_extension); 
//	

	
	// get the channel relief for the patches - this is to stop one terrace mapping to 2 different channels.  DO THIS SEPARATELY FOR FLOODPLAIN AND TERRACES OR IT ALL GOES TO SHIT
	
	// get the distance from outlet
	LSDRaster DistFromOutlet = FlowInfo.distance_from_outlet();
	
	
	// get the raster of channel relief masked by the terrace mask
//	cout << "\t Getting channel relief for the terraces" << endl;
//	LSDRaster TerraceReliefMasked = ChanNetwork.calculate_relief_from_channel_connected_components(filled_topo_test, TerracePatches, DistFromOutlet, FlowInfo, threshold_SO, search_distance);
//	string Terrace_relief_masked = "_terrace_relief_masked";
//	TerraceReliefMasked.write_raster((input_path+DEM_ID+Terrace_relief_masked), flt_extension); 
//	
//	// get the channel relief data to a text file for plotting
//	
//	//merge the rasters together
//	LSDRaster MergedRaster = FP_Relief_Masked.MergeRasters(TerraceReliefMasked);
//	string merged_ext = "_relief_merged";
//	MergedRaster.write_raster((input_path+DEM_ID+merged_ext), flt_extension); 
//	
//	string out_file_ext = "_channel_relief";
//	MergedRaster.write_RasterData_to_text_file(input_path+DEM_ID+out_file_ext);
	
	//test the new relief and distance code
	LSDRaster MainStemRelief, UpstreamDistance;
	cout << "The junction number is: " << junction_number << endl;
	
	// get the main stem channel from the junction number
	LSDIndexChannel MainStem = ChanNetwork.generate_longest_index_channel_from_junction(junction_number, FlowInfo, DistFromOutlet);
	LSDIndexRaster ChannelRaster = MainStem.print_index_channel_to_index_raster();
	string csv_fname = "_main_stem";
	ChannelRaster.FlattenToCSV(path_name+DEM_ID+csv_fname);
	ChanNetwork.get_information_about_nearest_main_stem_channel_connected_components(ConnectedComponents_final, filled_topo_test, DistFromOutlet, FlowInfo, MainStem, junction_number, MainStemRelief, UpstreamDistance);
	
	string relief_ext = "_relief_MS";
	string dist_ext = "_dist_MS";
	
	MainStemRelief.write_raster((input_path+DEM_ID+relief_ext), flt_extension); 
	UpstreamDistance.write_raster((input_path+DEM_ID+dist_ext), flt_extension); 
	
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
