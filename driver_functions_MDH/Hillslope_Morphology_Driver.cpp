//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Hillslope_Morphology_Driver.cpp
//
// Driver file for automating the extraction of hillslope traces and hillslope morphology
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
#include "../LSDParameterParser.hpp"

int main (int nNumberofArgs,char *argv[])
{
	cout << endl;
	cout << "======================================================" << endl;
	cout << "|| Welcome to the hillslope morphology tool!        ||" << endl;
	cout << "|| This program maps hilltops and extracts          ||" << endl;
	cout << "|| hillslope morphlogy by running individual traces ||" << endl;
	cout << "|| from the hilltops to channels.                   ||" << endl;
	cout << "|| This program was developed by the Land Surface   ||" << endl;
	cout << "|| Dynamics group at University of Edinburgh        ||" << endl;
	cout << "======================================================" << endl;	
  cout << endl;
  
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
    cout << "Error: This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "  ./Hillslope_Morphology_Tool.exe /ProjectFolder/ project_hillslopes_driver.param" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << "  see readme and online documentation." << endl;
    cout << "  web address goes here" << endl;
    cout << "------------------------------------------------------" << endl;
    exit(EXIT_SUCCESS);
	}
	
	//load input arguments
	string path = argv[1];
	string filename = argv[2];
	
	// load parameter parser object
  LSDParameterParser LSDPP(path,filename);
  
  // make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // setup default parameters
  // input files
  string_default_map["CHeads_file"] = "NULL";
  string_default_map["ChannelSegments_file"] = "NULL";
	string_default_map["Floodplain_file"] = "NULL";
	
	// surface metrics 
	float_default_map["MinSlopeForFill"] = 0.0001;
	bool_default_map["raster_is_filled"] = false;
	float_default_map["WindowRadius"] = 12.;
  bool_default_map["use_Dinf_flow"] = false;
  
   //Defining hilltops
	int_default_map["StreamNetworkPadding"] = 0;
	int_default_map["min_stream_order_to_extract_basins"] = 0;
	int_default_map["max_stream_order_to_extract_basins"] = 100;
	bool_default_map["RemovePositiveHilltops"] = true;
	bool_default_map["RemoveSteepHilltops"] = true;
	float_default_map["Threshold_Hilltop_Gradient"] = 0.4;
	bool_default_map["MaskHilltopstoBasins"] = true;
	
	// Selecting basins
  int_default_map["threshold_contributing_pixels"] = 1000;
  int_default_map["minimum_basin_size_pixels"] = 5000;
  int_default_map["maximum_basin_size_pixels"] = 500000;
  bool_default_map["test_drainage_boundaries"] = true;
  bool_default_map["only_take_largest_basin"] = false;
  string_default_map["BaselevelJunctions_file"] = "NULL";
  
  // basin. For this reason the default is to test for edge effects
  bool_default_map["find_complete_basins_in_window"] = false;
  bool_default_map["find_largest_complete_basins"] = false;
  bool_default_map["print_basin_raster"] = false;
  
	//writing results	
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
  bool_default_map["write_hilltop_curvature"] = false;
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
  LSDPP.print_parameters();
  
  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string RasterExt =  LSDPP.get_dem_read_extension();
  vector<string> BoundaryConditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  string BaselevelJunctions_file = LSDPP.get_BaselevelJunctions_file();
  string ChannelSegments_file = LSDPP.get_ChannelSegments_file();
  string Floodplain_file = LSDPP.get_Floodplain_file();
  
  //load the DEM
	cout << "\tLoading the DEM..." << endl;
	LSDRaster DEM((DATA_DIR+DEM_ID), RasterExt);	
		
	//Fill the DEM
	LSDRaster FilledDEM;
	if ( this_bool_map["raster_is_filled"] )
  {
    cout << "\tDEM already filled." << endl;
    FilledDEM = DEM;
  }
  else
  {
    cout << "\tFilling DEM..." << endl;
	  FilledDEM = DEM.fill(this_float_map["MinSlopeForFill"]);
	}
	
	//Surface fitting to get slope, aspect, curvature and planform curvature
	static const int Arr[] = {0,1,1,1,1,0,0,0};
	vector<int> RasterSelection (Arr, Arr + sizeof(Arr) / sizeof(Arr[0]));
	cout << "\tCalculating surface metrics..." << endl;
	vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(this_float_map["WindowRadius"], RasterSelection);
  LSDRaster Aspect = Surfaces[2];
  
	// get a flow info object	
	cout << "\tExtracting the drainage network..." << endl;
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);
	
	if (this_bool_map["use_Dinf_flow"])
	{
	  //get d infinity flowdirection
    cout << "\tGetting Dinf flow directions" << endl;
    Array2D<float> dinf = FilledDEM.D_inf_FlowDir();
    Aspect = FilledDEM.LSDRasterTemplate(dinf);
	}
	
	// calculate the flow accumulation
  cout << "\tCalculating flow accumulation (in pixels)..." << endl;
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  
	// get the channel heads either from file or from a threshold drainage area
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    cout << "\tGetting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels..." << endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);
  }
  else
  {
    cout << "\tLoading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
  }
  
	//set up the stream network to get index raster for channel mapping
	LSDIndexRaster StreamNetwork;
	LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);
	
	if (ChannelSegments_file == "NULL" || ChannelSegments_file == "Null" || ChannelSegments_file == "null")
  {
    cout << "\tGetting stream network from source nodes..." << endl;
	  StreamNetwork = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  }
  else
  {
    cout << "\tGetting stream network from Channel segments file... " << DATA_DIR+ChannelSegments_file << endl;
    StreamNetwork = LSDIndexRaster((DATA_DIR+ChannelSegments_file), RasterExt);
  }

  //pad the stream network by a number of pixels
  if (this_int_map["StreamNetworkPadding"] > 0)
  {
    cout << "\tPadding the stream network by " << this_int_map["StreamNetworkPadding"] << " pixels..." << endl;
    StreamNetwork.PadRaster(this_int_map["StreamNetworkPadding"]);
  }

  //load floodplain and merge with the channel network if required, otherwise the 
 	//floodplain mask will only contain the channel data
  if (Floodplain_file == "NULL" || Floodplain_file == "Null" || Floodplain_file == "null") {}
  else
  {
    cout << "\tCombining the channelnetwork and floodplain masks..." << endl;
		LSDIndexRaster Floodplains((Floodplain_file), RasterExt);
		StreamNetwork.MergeIndexRasters(Floodplains);
	}
	
  //Check to see if a list of junctions for extraction exists
  // need to get base-level nodes , otherwise these catchments will be missed!
  vector< int > BaseLevelJunctions;
  vector< int > BaseLevelJunctions_Initial;
  if (BaselevelJunctions_file == "NULL" || BaselevelJunctions_file == "Null" || BaselevelJunctions_file == "null" || BaselevelJunctions_file.empty() == true)
  {
    cout << "\tSelecting basins..." << endl;
    // remove basins drainage from edge if that is what the user wants
    if (this_bool_map["find_complete_basins_in_window"])
    {
      cout << "\tFinding basins not influended by nodata and removing nested basins..." << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(FlowInfo, FilledDEM, FlowAcc,
                                              this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
    }
    else
    {
      //Get baselevel junction nodes from the whole network
      BaseLevelJunctions_Initial = JunctionNetwork.get_BaseLevelJunctions();

      // now prune these by drainage area
      cout << "\tRemoving basins with fewer than " << this_int_map["minimum_basin_size_pixels"] << " pixels" << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions_Initial,
                                              FlowInfo, FlowAcc, this_int_map["minimum_basin_size_pixels"]);
      cout << "\tNow I have " << BaseLevelJunctions.size() << " baselelvel junctions." << endl;

      if (this_bool_map["find_largest_complete_basins"])
      {
        cout << "\tFinding largest basin not influenced by nodata..." << endl;
        BaseLevelJunctions = JunctionNetwork.Prune_To_Largest_Complete_Basins(BaseLevelJunctions,FlowInfo, FilledDEM, FlowAcc);
      }
      else
      {
        if (this_bool_map["test_drainage_boundaries"])     // now check for edge effects
        {
          cout << endl << endl << "\tRemoving basins draining to the edge." << endl;
          BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions,FlowInfo, FilledDEM);
          //BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge(BaseLevelJunctions,FlowInfo);
        }
      }
    }
  }
  else if (this_bool_map["extract_basins_by_stream_order"])
  {
    //Or just use a list of basins?
  	cout << "\tExtracting baselevel nodes for a strahler stream order of " << this_int_map["stream_order_to_extract_basins"] << "..." << endl;
  	BaseLevelJunctions = JunctionNetwork.ExtractBasinJunctionOrder(this_int_map["stream_order_to_extract_basins"], FlowInfo);

  }
  else
  {
    //specify junctions to work on from a list file
    string JunctionsFile = DATA_DIR+DEM_ID+"_junctions.list";

    cout << "\tReading junctions from a junction list... " << JunctionsFile << endl;

    vector<int> JunctionsList;
    ifstream infile(JunctionsFile.c_str());
    if (infile)
    {
      int n;
      while (infile >> n) BaseLevelJunctions_Initial.push_back(n);
    }
    else
    {
      cout << "Fatal Error: Junctions File " << JunctionsFile << " does not exist" << endl;
      exit(EXIT_FAILURE);
    }
    
    // Now make sure none of the basins drain to the edge
    cout << "\tPruning junctions that drain to the edge of the DEM..." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions_Initial, FlowInfo, FilledDEM);
  }
		
	//get the basins
	cout << "\tExtracting basins..." << endl;
	LSDIndexRaster BasinsRaster = JunctionNetwork.extract_basins_from_junction_vector(BaseLevelJunctions, FlowInfo);

	// Extract hilltops
	cout << "\tExtracting hilltop network..." << endl;
	LSDRaster Hilltops = JunctionNetwork.ExtractRidges(FlowInfo, this_int_map["min_stream_order_to_extract_basins"], this_int_map["max_stream_order_to_extract_basins"]); 
  
  // Mask to only use hilltops inside basins being analysed
  if (this_bool_map["MaskHilltopstoBasins"])
  {
    cout << "\tMasking hilltops to basins being analysed..." << endl;
    Hilltops.MaskRaster(BasinsRaster);
  }
   
	if (this_bool_map["RemoveSteepHilltops"])
	{
	  cout << "\tFiltering out hilltops with slopes greater than " << this_float_map["Threshold_Hilltop_Gradient"] << "..." << endl;	
	  string Condition = "<";
    LSDIndexRaster SteepHilltopsMask = Surfaces[1].Create_Mask(Condition,this_float_map["Threshold_Hilltop_Gradient"]);
    Hilltops.MaskRaster(SteepHilltopsMask);
  }

  if (this_bool_map["RemovePositiveHilltops"])
	{
	  cout << "\tFiltering out hilltops with positive curvature..." << endl;	
	  float Zero = 0;
	  string Condition = "<";
    LSDIndexRaster NegativeCurvatureMask = Surfaces[3].Create_Mask(Condition,Zero);
    Hilltops.MaskRaster(NegativeCurvatureMask);
	}
	
  //get hilltop curvature using filters to remove positive curvatures and steep slopes			 
  cout << "\tGetting hilltop curvature..." << endl;
	LSDRaster CHT = FilledDEM.get_hilltop_curvature(Surfaces[3], Hilltops);
		
	// Run hilltop flow routing 
	vector<int> Target_Basin_Vector;
	cout << "\tRunning hillslope flow routing..." << endl;	
	string HillslopeTracesFolder = DATA_DIR;
 	vector< Array2D<float> > HFR_Arrays = FlowInfo.HilltopFlowRouting(FilledDEM, Hilltops, Surfaces[1], StreamNetwork, Surfaces[2], (DATA_DIR+DEM_ID), BasinsRaster, Surfaces[4], this_bool_map["print_hillslope_traces"], this_int_map["hillslope_trace_thinning"], HillslopeTracesFolder, this_bool_map["hillslope_traces_basin_filter"], Target_Basin_Vector);
 	cout << "\tDone" << endl;

  // Write rasters to file

  //fill
  if (this_bool_map["print_fill_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    cout << "\tWriting filled topography to "  << filled_raster_name <<RasterExt << endl;
    FilledDEM.write_raster(filled_raster_name,RasterExt);
  }

  //hillshade
  if (this_bool_map["write_hillshade"])
  {
    string hillshade_raster_name = OUT_DIR+OUT_ID+"_Hillshade";
    cout << "Writing hillshade to " << hillshade_raster_name << RasterExt << endl;
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
    cout << "\tWriting slope to " << slope_raster_name << RasterExt << "..." << endl;
    Surfaces[1].write_raster((OUT_DIR+OUT_ID+"_Slope"),RasterExt);
    cout << "done.";
  }
  //aspect
  if (this_bool_map["write_aspect"])
  {
    string aspect_raster_name = OUT_DIR+OUT_ID+"_Aspect";
    cout << "\tWriting aspect to " << aspect_raster_name << RasterExt << "..." << endl;
    Surfaces[2].write_raster(aspect_raster_name,RasterExt);
    cout << "done.";
  }
  //curvature
  if (this_bool_map["write_curvature"])
  {
    string curvature_raster_name = OUT_DIR+OUT_ID+"_Curvature";
    cout << "\tWriting curvature to " << curvature_raster_name << RasterExt << "..." << endl;
    Surfaces[3].write_raster(curvature_raster_name,RasterExt);
    cout << "done.";
  }
  //planform curvature
  if (this_bool_map["write_planform_curvature"])
  {
    string plancurv_raster_name = OUT_DIR+OUT_ID+"_PlanCurv";
    cout << "\tWriting planform curvature to " << plancurv_raster_name << RasterExt << "..." << endl;
    Surfaces[4].write_raster(plancurv_raster_name, RasterExt);
    cout << "done.";
  }
  //stream network
  if (this_bool_map["write_stream_network"])
  {
    string stnet_raster_name = OUT_DIR+OUT_ID+"_MaskedStreamNetwork";
    cout << "\tWriting stream newtork mask to " << stnet_raster_name << RasterExt << "..." << endl;
    StreamNetwork.write_raster(stnet_raster_name, RasterExt);
    cout << "done.";
  }
  //basins
  if (this_bool_map["write_basins"])
  {
    string basins_raster_name = OUT_DIR+OUT_ID+"_Basins";
    cout << "\tWriting basins mask to " << basins_raster_name << RasterExt << "..." << endl;
    BasinsRaster.write_raster(basins_raster_name, RasterExt);
    cout << "done.";
  }
  //hilltops
  if (this_bool_map["write_hilltops"])
  {
    string hilltops_raster_name = OUT_DIR+OUT_ID+"_Hilltops";
    cout << "\tWriting hilltops to " << hilltops_raster_name << RasterExt << "..." << endl;
    Hilltops.write_raster(hilltops_raster_name, RasterExt);
    cout << "done.";
  }
  //hilltop curvature
  if (this_bool_map["write_hilltop_curvature"])
  {
    string hilltop_curvature_raster_name = OUT_DIR+OUT_ID+"_HFR_CHT";
    cout << "\tWriting hilltop curvature to " << hilltop_curvature_raster_name << RasterExt << "..." << endl;
    Hilltops.write_raster(hilltop_curvature_raster_name, RasterExt);
    cout << "done.";
  }
  //hillslope length raster
  if (this_bool_map["write_hillslope_length"])
  {
    string hillslope_length_raster_name = OUT_DIR+OUT_ID+"_HFR_LH";
    cout << "\tWriting hillslope_lengths to " << hillslope_length_raster_name << RasterExt << "..." << endl;
    LSDRaster HFR_LH = CHT.LSDRasterTemplate(HFR_Arrays[1]);
    HFR_LH.write_raster(hillslope_length_raster_name,RasterExt); 
    cout << "done.";
  }
  //hillslope gradient raster
  if (this_bool_map["write_hillslope_gradient"])
  {
    string hillslope_gradient_raster_name = OUT_DIR+OUT_ID+"_HFR_SLP";
    cout << "\tWriting hillslope gradients to " << hillslope_gradient_raster_name << RasterExt << "..." << endl;
    LSDRaster HFR_Slope = CHT.LSDRasterTemplate(HFR_Arrays[2]);
    HFR_Slope.write_raster(hillslope_gradient_raster_name,RasterExt);
    cout << "done.";
  }
  //hillslope relief raster
  if (this_bool_map["write_hillslope_relief"])
  {
    string hillslope_relief_raster_name = OUT_DIR+OUT_ID+"_HFR_Relief";
    cout << "\tWriting hillslope relief to " << hillslope_relief_raster_name << RasterExt << "..." << endl;
    LSDRaster relief = CHT.LSDRasterTemplate(HFR_Arrays[3]);
    relief.write_raster(hillslope_relief_raster_name,RasterExt);
    cout << "done.";
  }
}
