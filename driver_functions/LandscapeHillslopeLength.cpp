//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LandscapeHillslopeLength.cpp
//
// Driver created to automate extraction of hillslope lengths husing hilltop flow routing 
// from a full DEM tile.
//
// Needs a parameter file with the following data, one value per line:
// DEM_name 
// DEM_file_extension 
// MinSlope
// Drainage Threshold 
// window_radius 
// critical_slope 
//
// Run with the arguments LandscapeHillslopeLength.out <input data path> <output data path> <name of parameter file> 
//
// Outputs a text file (into the output directory called <DEM_filename> + "_HilltopData" ) with the format:
// "hilltop_i hilltop_j hilltop_easting hilltop_northing stream_i stream_j stream_easting stream_northing stream_id basin_id relief lh aspect slope"
//
// from which the hilltops which reach a defined order's lengths can be filtered out using python.
//
// The HilltopData file generated is ~ 1.5x the size of the DEM file, so make 
// sure there is enough space to store this file.  
// 
// Version 0.1 - Hack of CollectBasinMetrics to be run on a full DEM tile.
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
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDMostLikelyPartitionsFinder.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs != 4)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path to the input data, the path to write the outputs and the name of the paramter file." << endl;
		exit(EXIT_SUCCESS);
	}

  //get the input and output paths
  string InputPath = argv[1]; 
  string OutputPath = argv[2];
  
  //get input file names
  string param_file_name = argv[3];

  //load parameters
  string DEM_name;
  string DEM_file_extension;
  float MinSlope;
  int threshold;
  float critical_slope;
  float window_radius;
  
  ifstream param_file_in;
	param_file_in.open((InputPath+param_file_name).c_str());
	if( param_file_in.fail() )
	{
		cout << "\nFATAL ERROR: the parameters file \"" << param_file_name << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
  
  param_file_in >> DEM_name >> DEM_file_extension >> MinSlope >> threshold >> window_radius >> critical_slope;
  param_file_in.close();
  
  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  // coefficent matrices for polyfit routine
  Array2D<float> a;
  Array2D<float> b;
  Array2D<float> c;
  Array2D<float> d;
  Array2D<float> e;
  Array2D<float> f;

  //Start processessing topographic data here
  
  //Load and fill the DEM
	LSDRaster DEM((InputPath+DEM_name), DEM_file_extension);
	LSDRaster FilledDEM = DEM.fill(MinSlope);
 	
	// get a flow info object
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  //get channel network
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	vector<int> sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  
  //Create slope, curvature, aspect and elevation rasters 
  FilledDEM.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);  
  LSDRaster slope = FilledDEM.calculate_polyfit_slope(d,e);  
   
  //get stream network
  LSDIndexRaster StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
  
  //create blank Basins index raster using the streamnetwork as a template - has the effect
  // of listing every hilltop px in a basin with the junction ID of NDV
  Array2D<int> FakeBasins(StreamNetwork.get_NRows(), StreamNetwork.get_NCols(), StreamNetwork.get_NoDataValue());
  LSDIndexRaster Basins = StreamNetwork.LSDRasterTemplate(FakeBasins);
      
  // extract ridges and then hilltops based on critical slope
  LSDRaster Ridges = ChanNetwork.ExtractRidges(FlowInfo);  
  LSDRaster hilltops = ChanNetwork.ExtractHilltops(Ridges, slope, critical_slope);
 
  // D-infinty flowdirection for HFR 
  Array2D<float> FlowDir_inf = FilledDEM.D_inf_FlowDir();
  
  //hilltop flow routing
  vector< Array2D<float> > Routed_Hilltop_Data = FilledDEM.HFR_Driver(hilltops, FlowDir_inf, StreamNetwork, Basins, (OutputPath+DEM_name));
  
  //write routed_hilltop_data to an LSDRaster
  LSDRaster RoutedHilltops = FilledDEM.LSDRasterTemplate(Routed_Hilltop_Data[0]);
      
}
