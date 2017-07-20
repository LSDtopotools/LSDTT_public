//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LH_Driver.cpp
//
// Driver created to automate extraction of hillslope profiles
//
// Driver expects an unfilled DEM in the given directory in flt format with the name format <prefix>.bil,
// a channel heads file in the same directory with the format <prefix>_CH.flt and a floodplain raster
// <prefix>_FloodPlain.flt if floodplains are to be excluded from the analysis.
//
// Run driver with the following arguments:
//
// path to the input files with a trailing slash
// filename prefix without an underscore
// window radius value in spatial units for surface fitting
// basin order the strahler order of basins to be extracted
// switch to exclude floodplains (1) or not exclude them (0) 
// switch to write rasters 0 == do not write rasters and 1 == write rasters
//
// A usage example is:
//nice ./LH_Driver.out /home/s0675405/DataStore/Final_Paper_Data/NC/ NC 7 2 1 0
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Martin D. Hurst
// University of Glasgow
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
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{
	
	//Test for correct input arguments
	if (nNumberofArgs!=7)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the filename prefix, window radius, "; 
		cout << "basin order, a switch to use or exclude floodplains and a switch to write rasters if desired." << endl;
		exit(EXIT_SUCCESS);
	}
	
	//load input arguments
	string path = argv[1];
	string filename = argv[2];
	
	// load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);
  
  // make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // setup default parameters
  string_default_map["CHeads_file"] = "NULL";
  string_default_map["ChannelSegments_file"] = "NULL";
	string_default_map["Floodplain_file"] = "NULL";
	 
	float_default_map["MinSlopeForFill"] 0.0001;
	float_default_map["WindowRadius"] = 12.;
	int_default_map["threshold_contributing_pixels"] = 1000;
	int_default_map["BasinOrder"] = 2;
	this_int_map["StreamNetworkPadding"] = 0;
	
	bool_default_map["RemovePositiveHilltops"] = true;
	bool_default_map["RemoveSteepHilltops"] = true;
	float_default_map["Threshold_Hilltop_Gradient"] = 0.4;
	
	bool_default_map["print_fill_raster"] = false;
	bool_default_map["write_hillshade"] = false;
	bool_default_map["write_hillshade"] = false;
  bool_default_map["write_slope"] = false;
  bool_default_map["write_aspect"] = false;
  bool_default_map["write_curvature"] = false;
  bool_default_map["write_planform_curvature"] = false;
  bool_default_map["write_stream_network"] = false;
  bool_default_map["write_basins"] = false;
  bool_default_map["write_hilltops"] = false;
  bool_default_map["write_hillslope_length"] = false;
  bool_default_map["write_hillslope_gradient"] = false;
  bool_default_map["write_hillslope_relief"] = false;
  
	// these params do not need changed during normal use of the HFR algorithm
	bool_default_map["print_hillslope_traces"] = false;
	int_default_map["hillslope_trace_thinning"] = 1;
	string_default_map["hillslope_traces_file"] = "";
	bool_default_map["hillslope_traces_basin_filter"] = false;													
	
	// Use the parameter parser to get the maps of the parameters required for the analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();
  
	// Now print the parameters for bug checking
  cout << "PRINT THE PARAMETERS..." << endl;
  LSDPP.print_parameters();
  
  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string RasterExt =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  string BaselevelJunctions_file = LSDPP.get_BaselevelJunctions_file();
    	
	//load the DEM
	cout << "\n\t Loading the DEM...";
	LSDRaster DEM((OUT_DIR+OUT_ID), RasterExt);	
	cout << " done."
		
	//Fill the DEM
	cout << "\n\tFilling DEM...";
	LSDRaster FilledDEM = DEM.fill(this_float_map["MinSlopeForFill"]);
	cout << " done."
	
	//Surface fitting to get slope, aspect, curvature and planform curvature
	static const int Arr[] = {0,1,1,1,1,0,0,0}
	vector<int> RasterSelection (Arr, Arr + sizeof(Arr) / sizeof(Arr[0]));
	cout << "\n\tCalculating surface metrics...";
	vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(WindowRadius, RasterSelection);
	cout << " done."
	
	//LSDRaster slope = Surfaces[1];
	//LSDRaster aspect = Surfaces[2];	 

	// get a flow info object	
	cout << "\n\tExtracting the drainage network...";
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);
	cout << " done."
	
	// ge the channel heads either from file or from a threshold drainage area
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    cout << "\n\tGetting sources from a threshold of "<< threshold_contributing_pixels << " pixels...";
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);
    cout << " done.";  
  }
  else
  {
    cout << "\n\tLoading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << " done.";  
  }
  
	//set up the stream network to get index raster for channel mapping
	LSDIndexRaster StreamNetwork;
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	
	if (ChannelSegments_file == "NULL" || ChannelSegments_file == "Null" || ChannelSegments_file == "null")
  {
    cout << "\n\tGetting stream network from source nodes..."
	  StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
    cout << " done.";
  }
  else
  {
    cout << "\n\tGetting stream network from Channel segments file..." << DATA_DIR+ChannelSegments_file;
    StreamNetwork = LSDIndexRaster(ChannelSegments_file, RasterExt);
    cout << " done.";
  }

  //pad the stream network by a number of pixels
  if (this_int_map["StreamNetworkPadding"] > 0)
  {
    cout << "\n\tPadding the stream network by" << this_int_map["StreamNetworkPadding"] << "pixels";
    StreamNetwork.PadRaster(this_int_map["StreamNetworkPadding"]);
    cout << " done.";
  }

  //load floodplain and merge with the channel network if required, otherwise the 
 	//floodplain mask will only contain the channel data
  if (Floodplain_file == "NULL" || Floodplain_file == "Null" || Floodplain_file == "null") {}
  else
  {
    cout << "\n\tCombining the channelnetwork and floodplain masks...";
		LSDIndexRaster Floodplains((Floodplain_file), ext);
		StreamNetwork.MergeIndexRasters(Floodplains);
		cout << " done.";
	}
	
	//Extract basins based on input stream order
	//Or just use a list of basins?
//	cout << "\n\tExtracting basins for a strahler stream order of " << BasinOrder << "...";
//	vector<int> BasinsJunctions = ChanNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
//	LSDIndexRaster BasinsRaster = ChanNetwork.extract_basins_from_junction_vector(BasinsJunctions, FlowInfo);
//	cout << " done.";
	
	// Extract hilltops
	cout << "\n\tExtracting hilltop network...";
	LSDRaster Hilltops = ChanNetwork.ExtractRidges(FlowInfo); 
  cout << " done.";
  
  //get hilltop curvature using filters to remove positive curvatures and steep slopes			 
  cout << "\n\tGetting hilltop curvature...";
	LSDRaster CHT = FilledDEM.get_hilltop_curvature(Surfaces[3], Hilltops);
	cout << " done.";
	
	if (this_bool_map["RemovePositiveHilltops"])
	{
    cout << "\n\tFiltering out positive hilltop curvature...";	
	  LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(CHT);
	  cout << " done.";
	}
	
	if (this_bool_map["RemoveSteepHilltops"])
	{
	  cout << "\n\tFiltering out hilltops with slopes greater than" << this_float_map["Threshold_Hilltop_Gradient"] << "...";	
    LSDRaster CHT = ChanNetwork.ExtractHilltops(CHT, Surfaces[1], this_float_map["Threshold_Hilltop_Gradient"]);
    cout << " done.";
  }
  
	// Run hilltop flow routing 
	vector<int> Target_Basin_Vector;
	cout << "\n\tRunning hillslope flow routing...";	
 	vector< Array2D<float> > HFR_Arrays = FlowInfo.HilltopFlowRouting(FilledDEM, Hilltops, Surfaces[1], StreamNetwork, Surfaces[2], (DATA_DIR+DEM_ID), BasinsRaster, Surfaces[4], this_bool_map["print_hillslope_traces"], this_int_map["hillslope_trace_thinning"], this_string_map["hillslope_traces_file"], this_bool_map["hillslope_traces_basin_filter"], Target_Basin_Vector);
 	cout << "done.";

  // Write rasters to file

  //fill
  if (this_bool_map["print_fill_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    cout << "\n\tWriting filled topography to "  << filled_raster_name <<RasterExt;
    FilledDEM.write_raster(filled_raster_name,RasterExt);
  }

  //hillshade
  if (this_bool_map["write_hillshade"])
  {
    string hillshade_raster_name = OUT_DIR+OUT_ID+"_Hillshade";
    cout << "Writing hillshade to " << hillshade_raster_name << RasterExt;
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster HillshadeRaster = DEM.hillshade(hs_altitude,hs_azimuth,hs_z_factor);
    HillshadeRaster.write_raster(hillshade_raster_name,RasterExt);
  }

  //slope
  if (this_bool_map["write_slope"])
  {
    string slope_raster_name = OUT_DIR+OUT_ID+"_Slope";
    cout << "\n\tWriting slope to " << slope_raster_name << RasterExt << "...";
    Surfaces[1].write_raster((OUT_DIR+OUT_ID+"_Slope"),RasterExt);
    cout << "done.";
  }
  //aspect
  if (this_bool_map["write_aspect"])
  {
    string aspect_raster_name = OUT_DIR+OUT_ID+"_Aspect";
    cout << "\n\tWriting aspect to " << aspect_raster_name << RasterExt << "...";
    Surfaces[2].write_raster(aspect_raster_name,RasterExt);
    cout << "done.";
  }
  //curvature
  if (this_bool_map["write_curvature"])
  {
    string curvature_raster_name = OUT_DIR+OUT_ID+"_Curvature";
    cout << "\n\tWriting curvature to " << curvature_raster_name << RasterExt << "...";
    Surfaces[3].write_raster(curvature_raster_name,RasterExt);
    cout << "done.";
  }
  //planform curvature
  if (this_bool_map["write_planform_curvature"])
  {
    string plancurv_raster_name = OUT_DIR+OUT_ID+"_PlanCurv";
    cout << "\n\tWriting planform curvature to " << plancurv_raster_name << RasterExt << "...";
    Surfaces[4].write_raster(plancurv_raster_name, RasterExt);
    cout << "done.";
  }
  //stream network
  if (this_bool_map["write_stream_network"])
  {
    string stnet_raster_name = OUT_DIR+OUT_ID+"_MaskedStreamNetwork";
    cout << "\n\tWriting stream newtork mask to " << stnet_raster_name << RasterExt << "...";
    StreamNetwork.write_raster(stnet_raster_name, RasterExt);
    cout << "done.";
  }
  //basins
  if (this_bool_map["write_basins"])
  {
    string basins_raster_name = OUT_DIR+OUT_ID+"_Basins";
    cout << "\n\tWriting basins mask to " << basins_raster_name << RasterExt << "...";
    BasinsRaster.write_raster(basins_raster_name, RasterExt);
    cout << "done.";
  }
  //hilltops
  if (this_bool_map["write_hilltops"])
  {
    string hilltops_raster_name = OUT_DIR+OUT_ID+"_Hilltops";
    cout << "\n\tWriting hilltops to " << hilltops_raster_name << RasterExt << "...";
    CHT.write_raster(hilltops_raster_name, RasterExt);
    cout << "done.";
  }
  //hillslope length raster
  if (this_bool_map["write_hillslope_length"])
  {
    string hillslope_length_raster_name = OUT_DIR+OUT_ID+"_HFR_LH";
    cout << "\n\tWriting hillslope_lengths to " << hillslope_length_raster_name << RasterExt << "...";
    LSDRaster HFR_LH = CHT.LSDRasterTemplate(HFR_Arrays[1]);
    HFR_LH.write_raster(hillslope_length_raster_name,RasterExt); 
    cout << "done.";
  }
  //hillslope gradient raster
  if (this_bool_map["write_hillslope_gradient"])
  {
    string hillslope_gradient_raster_name = OUT_DIR+OUT_ID+"_HFR_SLP";
    cout << "\n\tWriting hillslope gradients to " << hillslope_gradient_raster_name << RasterExt << "...";
    LSDRaster HFR_Slope = CHT.LSDRasterTemplate(HFR_Arrays[2]);
    HFR_Slope.write_raster(hillslope_gradient_raster_name,RasterExt);
    cout << "done.";
  }
  //hillslope relief raster
  if (this_bool_map["write_hillslope_relief"])
  {
    string hillslope_relief_raster_name = OUT_DIR+OUT_ID+"_HFR_Relief";
    cout << "\n\tWriting hillslope relief to " << hillslope_relief_raster_name << RasterExt << "...";
    LSDRaster relief = CHT.LSDRasterTemplate(HFR_Arrays[3]);
    relief.write_raster(hillslope_relief_raster_name,RasterExt);
    cout << "done.";
  }
}
