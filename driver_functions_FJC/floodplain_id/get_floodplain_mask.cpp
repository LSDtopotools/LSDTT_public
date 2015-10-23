//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_map_sensitivity_analysis
//
// This program takes two arguments, the path name and the driver name
// It then runs the chi analysis for the parameter space defined in the driver
// file, allowing the user to examine the effects of changing m/n value,
// number of target nodes, minimum segment length, sigma and target skip value.
// At present it just spits out an output file for each iteration.
// Now also kicks out a table with the parameters and best fit m/n.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 18/03/2014
// David Milodowski
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
	string output_path;
  string dem_ext = "_dem";
  string sources_ext = "_CH";
  string flt_extension = "flt";
  string txt_extension = ".txt";
  string Swath_ext = "_swath_trans";
  string Long_Swath_ext = "_swath_long";
  
  // initialise variables to be assigned from .driver file
  int junction_number,NormaliseTransProfile,NormaliseLongProfile,FilterTopo;
	float pruning_threshold = 1;
	float HalfWidth,BinWidth,Minimum_Slope,histogram_bin_width,knickpoint_position;
  string temp;
                                                            
	file_info_in >> temp >> DEM_ID
               >> temp >> RASTER_NAME
               >> temp >> input_path
               >> temp >> output_path
               >> temp >> junction_number
               >> temp >> Minimum_Slope
               >> temp >> HalfWidth
               >> temp >> BinWidth
               >> temp >> NormaliseTransProfile
               >> temp >> NormaliseLongProfile
               >> temp >> FilterTopo
               >> temp >> histogram_bin_width
               >> temp >> knickpoint_position;
                   
	string jn_name = itoa(junction_number);
	string uscore = "_";
	jn_name = uscore+jn_name;
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
  
  string SO_ext = "_SO_DrEICH";
  string JI_ext = "_JI_DrEICH";
  LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
 	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();

	SOArray.write_raster((output_path+DEM_ID+SO_ext),flt_extension);
  JIArray.write_raster((output_path+DEM_ID+JI_ext),flt_extension);
  
  //load in the slope and relief rasters
  LSDRaster Slope((input_path+DEM_ID+"_slope"), flt_extension);
  LSDRaster Relief((input_path+DEM_ID+"_relief"), flt_extension);
  
//   //get the relief
//   cout << "\t Calculating relief..." << endl;
//   float kernelWidth = 100;
//   int kernelType = 1;
//   LSDRaster Relief = topo_test.calculate_relief(kernelWidth, kernelType);
//   cout << "\t Done!" << endl;
//   string relief_name = "_relief";
//   Relief.write_raster((output_path+DEM_ID+relief_name), flt_extension);
//   
//   //get the slope
//   cout << "\t Calculating slope..." << endl;
//   float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
//   vector<LSDRaster> surface_fitting;
//   LSDRaster Slope;
//   vector<int> raster_selection(8, 0);
//   raster_selection[1] = 1;             // this means you want the slope
//   surface_fitting = filled_topo_test.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
//   Slope = surface_fitting[1];
//   cout << "\t Done!" << endl;
//   string slope_name = "_slope";
//   Slope.write_raster((output_path+DEM_ID+slope_name), flt_extension);
  
  //get the potential floodplain mask
  cout << "\t Getting the floodplain mask" << endl;
  float relief_threshold = 5;
  float slope_threshold = 0.5;
  LSDIndexRaster FloodplainRaster_temp = topo_test.get_potential_floodplain_patches(Relief, Slope, relief_threshold, slope_threshold);
  
  cout << "\t Connected components" << endl;
  LSDIndexRaster ConnectedComponents = FloodplainRaster_temp.ConnectedComponents();
  string CC_name = "_CC_filt";
  ConnectedComponents.write_raster((output_path+DEM_ID+CC_name), flt_extension); 
  
  cout << "\t Removing hillslope patches" << endl;
  LSDIndexRaster FloodplainRaster_modified = ChanNetwork.remove_hillslope_patches_from_floodplain_mask(ConnectedComponents);
  
  cout << "\t Removing holes" << endl;
  int window_radius = 50;
  LSDIndexRaster FloodplainRaster_noholes = FloodplainRaster_modified.remove_holes_in_patches(window_radius);
  
  int smaller_radius = 10;
  LSDIndexRaster FloodplainRaster_final = FloodplainRaster_noholes.remove_holes_in_patches(smaller_radius);
    
  cout << "\t Done!" << endl;
  string FP_name = "_FP_filt";
  FloodplainRaster_final.write_raster((output_path+DEM_ID+FP_name), flt_extension); 
}
