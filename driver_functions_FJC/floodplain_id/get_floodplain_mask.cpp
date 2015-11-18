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
// Fiona J Clubb
// University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../../LSDStatsTools.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../LSDRaster.hpp"
//#include "../../LSDRasterSpectral.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDShapeTools.hpp"


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

	string DEM_ID;
	string RASTER_NAME;
	string input_path;
  string dem_ext = "_dem";
  string sources_ext = "_CH";
  string flt_extension = "flt";
  string txt_extension = ".txt";
  string Swath_ext = "_swath_trans";
  string Long_Swath_ext = "_swath_long";
  
  // initialise variables to be assigned from .driver file
  int threshold_SO, FilterTopo, peak_distance;
	float Minimum_Slope,bin_width_relief, bin_width_slope, peak_threshold;
  string temp;
                                                            
	file_info_in >> temp >> DEM_ID
               >> temp >> RASTER_NAME
               >> temp >> input_path
               >> temp >> Minimum_Slope
               >> temp >> FilterTopo
               >> temp >> threshold_SO
               >> temp >> bin_width_relief
               >> temp >> bin_width_slope
               >> temp >> peak_threshold
               >> temp >> peak_distance;
                   
	file_info_in.close();

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// load the DEM
	LSDRaster topo_test((input_path+DEM_ID), flt_extension);  
  
  if(FilterTopo == 1)
  {
    // filter using Perona Malik
    int timesteps = 50;
    float percentile_for_lambda = 90;
    float dt = 0.1;
    topo_test = topo_test.PeronaMalikFilter(timesteps, percentile_for_lambda, dt);    
  }
  // fill
  LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
  cout << "\t Flow routing..." << endl;
	// get a flow info object
 	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
	
  // calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	cout << "\t Loading Sources..." << endl;
	// load the sources
	string CH_name = "_CH_DrEICH";
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_ID+CH_name),flt_extension);
  cout << "\t Got sources!" << endl;
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;
  //string SO_ext = "_SO_DrEICH";
  string JI_ext = "_JI_DrEICH";
  //LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
 	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();

	//SOArray.write_raster((output_path+DEM_ID+SO_ext),flt_extension);
  JIArray.write_raster((input_path+DEM_ID+JI_ext),flt_extension);
  
  //load in the slope and relief rasters
   //LSDRaster Slope((input_path+DEM_ID+"_slope"), flt_extension);
   //LSDRaster Relief((input_path+DEM_ID+"_relief"), flt_extension);
   
   //calculate the channel relief
   cout << "\t Getting relief relative to channel" << endl;
   LSDRaster ChannelRelief = ChanNetwork.calculate_relief_from_channel(filled_topo_test, FlowInfo, threshold_SO);
   string relief_name = "_channel_relief";
   ChannelRelief.write_raster((input_path+DEM_ID+relief_name), flt_extension);
   cout << "\t Got the relief!" << endl;
     
  //get the slope
  cout << "\t Calculating slope..." << endl;
  float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
  vector<LSDRaster> surface_fitting;
  LSDRaster Slope;
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;             // this means you want the slope
  surface_fitting = filled_topo_test.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
  Slope = surface_fitting[1];
  cout << "\t Done!" << endl;
  string slope_name = "_slope";
  Slope.write_raster((input_path+DEM_ID+slope_name), flt_extension);

  //get the relief threshold
  string relief_fname = DEM_ID+"_histogram_relief.txt";
  float relief_threshold = ChannelRelief.get_threshold_for_floodplain(bin_width_relief, peak_threshold, peak_distance, relief_fname);
  
  //get the slope threshold
  string slope_fname = DEM_ID+"_histogram_slope.txt";
  float slope_threshold = Slope.get_threshold_for_floodplain(bin_width_slope, peak_threshold, peak_distance, slope_fname);

  //get the potential floodplain mask
  cout << "\t Getting the floodplain mask" << endl;
  LSDIndexRaster FloodplainRaster_temp = topo_test.get_potential_floodplain_patches(ChannelRelief, Slope, relief_threshold, slope_threshold);
  
  cout << "\t Connected components" << endl;
  LSDIndexRaster ConnectedComponents = FloodplainRaster_temp.ConnectedComponents();
  string CC_name = "_CC_filt";
  ConnectedComponents.write_raster((input_path+DEM_ID+CC_name), flt_extension); 
  
  cout << "\t Removing hillslope patches" << endl;
  LSDIndexRaster FloodplainRaster_modified = ChanNetwork.remove_hillslope_patches_from_floodplain_mask(ConnectedComponents);
  
  cout << "\t Removing holes" << endl;
  int window_radius = 60;
  LSDIndexRaster FloodplainRaster_noholes = FloodplainRaster_modified.remove_holes_in_patches(window_radius);
  LSDIndexRaster FloodplainRaster_final = FloodplainRaster_noholes.remove_holes_in_patches(window_radius);
  LSDIndexRaster FloodplainRaster_superfinal = FloodplainRaster_final.remove_checkerboard_pattern();
    
  cout << "\t Done!" << endl;
  
  string FP_name2 = "_FP_final";
  FloodplainRaster_superfinal.write_raster((input_path+DEM_ID+FP_name2), flt_extension);
}
