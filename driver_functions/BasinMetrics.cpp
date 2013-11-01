//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Basin_Metrics.cpp
//
// Driver created to automate extraction of hillslope and basin metrics from DEM 
// files.
//
// Needs a parameter file with the following data, one value per line:
// DEM_name 
// DEM_file_extension 
// MinSlope 
// window_radius 
// critical_slope 
// log_bin_width
// SplineResolution
//
// And a file with a list of junction numbers of basins to be studied. If no basins are to be
// analysed pass in a single *valid* junction number to allow the rest of the code to run
// in future revisions this will be better dealt with.
//
// The HilltopData file generated is ~ 1.5x the size of the DEM file, so make 
// sure there is enough space to store this file.  
// 
// Version 1.1 - now gives flowarea from boomerang plots in spatial units instead of pixels.
//
// Version 1.2 - Modified to get and write the 2 LH measures from the Boomerang plotting function.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Stuart W.D. Grieve
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDChannel.hpp"
#include "../LSDChannelNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDMostLikelyPartitionsFinder.hpp"


int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the names of the paramter and junction files." << endl;
		exit(EXIT_SUCCESS);
	}

  //get input file names
  string param_file_name = argv[1];
  string junction_file_name = argv[2];

  //load parameters
  string DEM_name;
  string DEM_file_extension;
  double MinSlope;
  double critical_slope;
  double log_bin_width;
  double window_radius;
  double SplineResolution;
	
  ifstream param_file_in;
	param_file_in.open(param_file_name.c_str());
	if( param_file_in.fail() )
	{
		cout << "\nFATAL ERROR: the parameters file \"" << param_file_name << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
  
  param_file_in >> DEM_name >> DEM_file_extension >> MinSlope >> window_radius >> critical_slope >> log_bin_width >> SplineResolution;
  param_file_in.close();


  //load list of junctions  
  ifstream read_junction;
  read_junction.open(junction_file_name.c_str());
  
  if( read_junction.fail() )
	{
		cout << "\nFATAL ERROR: the junctions file \"" << junction_file_name << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
  
  vector<int> BasinJunctions;
  int tmp;
  while(read_junction >> tmp){
    BasinJunctions.push_back(tmp);
  }
  
  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  // coefficent matrices for polyfit routine
  Array2D<double> a;
  Array2D<double> b;
  Array2D<double> c;
  Array2D<double> d;
  Array2D<double> e;
  Array2D<double> f;
  
  //Start processessing topographic data here
  
  //Load and fill the DEM
	LSDRaster DEM(DEM_name, DEM_file_extension);
	LSDRaster FilledDEM = DEM.fill(MinSlope);
 	
	// get a flow info object
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  //------------------------------------------------------------------------------------------
  //This gets the channel network and will be updated with Fiona's channel head stuff 
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	int threshold = 800;
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);
	//end of channel head stuff
  //------------------------------------------------------------------------------------------
  
  //get basins from list of junctions
  LSDIndexRaster Basins = ChanNetwork.extract_basins_from_junction_vector(BasinJunctions, FlowInfo);
    
  //Create slope, curvature, aspect and elevation rasters 
  FilledDEM.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster curv = FilledDEM.calculate_polyfit_curvature(a,b);
  LSDRaster elevation = FilledDEM.calculate_polyfit_elevation(f);  //smoothing the elevation may give better boom plots, but would then need to re-generate slope and area: will this not just amplify noise?
  LSDRaster aspect = FilledDEM.calculate_polyfit_aspect(d,e);  
  LSDRaster slope = FilledDEM.calculate_polyfit_slope(d,e);
  
  //so maybe a smoothed_slope raster generated here for boom plots??
  //elevation.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  //LSDRaster smoothed_slope = FilledDEM.calculate_polyfit_slope(d,e);
    
  //get d infinity flowarea - in spatial units
  LSDRaster dinf_area = FilledDEM.D_inf_units();
  //LSDRaster dinf_area_smoothed = elevation.D_inf();   //could use the smoothed elevation here though? will need to test this with real data
    
  // get flow directions - D8 flowdir for drainage density calculations 
  Array2D<int> FlowDir_array = FlowInfo.get_FlowDirection();
   
  //get stream network
  LSDIndexRaster StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  
  //get basin area
  LSDRaster Area = FilledDEM.BasinArea(Basins);  
  
  //get drainage density
  LSDRaster drainage_density = FilledDEM.DrainageDensity(StreamNetwork, Basins, FlowDir_array);
  
  // extract ridges and then hilltops based on critical slope
  LSDRaster Ridges = ChanNetwork.ExtractRidges(FlowInfo);  
  LSDRaster hilltops = ChanNetwork.ExtractHilltops(Ridges, slope, critical_slope);
 
  // D-infinty flowdirection for HFR 
  Array2D<double> FlowDir_inf = FilledDEM.D_inf_FlowDir();
  
  //hilltop flow routing
  vector< Array2D<double> > Routed_Hilltop_Data = FilledDEM.HFR_Driver(hilltops, FlowDir_inf, StreamNetwork, Basins, DEM_name);
  
  //write routed_hilltop_data to a series of LSDRasters
  LSDRaster RoutedHilltops = FilledDEM.LSDRasterTemplate(Routed_Hilltop_Data[0]);
  LSDRaster RoutedHilltopLengths = FilledDEM.LSDRasterTemplate(Routed_Hilltop_Data[1]);
  LSDRaster RoutedHilltopSlopes = FilledDEM.LSDRasterTemplate(Routed_Hilltop_Data[2]);
  LSDRaster RoutedHilltopAspects = FilledDEM.LSDRasterTemplate(Routed_Hilltop_Data[3]);
  LSDRaster RoutedHilltopReliefs = FilledDEM.LSDRasterTemplate(Routed_Hilltop_Data[4]);
    
  //get hilltop curvature for all routed hilltop pixels
  LSDRaster CHT = FilledDEM.get_hilltop_curvature(curv, RoutedHilltops); // this is the routed hilltops

  //basin puncher - generate sub tiles for basins of interest for boomerang plotting  
  vector<LSDRaster> Basin_Slopes = slope.BasinPuncher(BasinJunctions, Basins);
  vector<LSDRaster> Basin_Drainage_Areas = dinf_area.BasinPuncher(BasinJunctions, Basins);
  
  //empty array to store the 2 LH measurements for each basin generated by the boomerang plots
  Array2D<double> HillslopeLengths(int(BasinJunctions.size()), 3);
   
  //boomerang plotting 
  for (int q = 0; q < int(Basin_Slopes.size()); ++q){
    //convert junction number here into a string to pass in for boomerang output filenames  
    stringstream ss;
    ss << BasinJunctions[q];
    string junction_string = ss.str();    
    pair<double,double> LH_Pair = Basin_Slopes[q].Boomerang(Basin_Slopes[q], Basin_Drainage_Areas[q], junction_string, log_bin_width, SplineResolution); 
    
    //write LH values and Basin IDs to the array
    HillslopeLengths[q][0] = BasinJunctions[q];
    HillslopeLengths[q][1] = LH_Pair.first;
    HillslopeLengths[q][2] = LH_Pair.second;          
  }
                                                               
  //Collect mean Basin Metrics into a text file 
  FilledDEM.CollectBasinMetrics(Basins, slope, elevation, aspect, 
                              Area, drainage_density, CHT, RoutedHilltopLengths,
                              RoutedHilltopSlopes, RoutedHilltopReliefs, RoutedHilltopAspects, HillslopeLengths, critical_slope, DEM_name);
  
}
