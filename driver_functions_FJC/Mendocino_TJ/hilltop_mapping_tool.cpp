//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// hilltop_mapping_tool.cpp
// make with make -f hilltop_mapping_tool.make
// FJC 20/11/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDChannel.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannel.hpp"
#include "../../LSDMostLikelyPartitionsFinder.hpp"
#include "../../LSDBasin.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDParameterParser.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDSpatialCSVReader.hpp"

int main (int nNumberofArgs,char *argv[])
{

	//Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the hilltops tool!        	              ||" << endl;
    cout << "|| This program takes the latitude and longitude       ||" << endl;
    cout << "|| of a point which defines the outlet of a basin.     ||" << endl;
		cout << "|| It then gets a load of info about the hilltops		  ||" << endl;
    cout << "|| upstream of that point.                             ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Fiona J. Clubb while at SAFL.                       ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./hilltop_mapping_tool.out /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_hilltops.param" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }
  string path_name = argv[1];
  string f_name = argv[2];

	// maps for setting default parameters
	map<string,int> int_default_map;
	map<string,float> float_default_map;
	map<string,bool> bool_default_map;
	map<string,string> string_default_map;
  map<string,double> double_default_map;

	// set default int parameters
	int_default_map["Chan area threshold"] = 1000;
  int_default_map["HalfWidth"] = 1000;
  int_default_map["search_radius"] = 25;
  int_default_map["Threshold_SO"] = 3;

  // set default float parameters
  float_default_map["surface_fitting_window_radius"] = 6;
  float_default_map["Min slope filling"] = 0.0001;

	// set default string parameters
	string_default_map["CHeads_file"] = "NULL";
  string_default_map["coords_csv_file"] = "NULL";

  // params for printing
	bool_default_map["write_hillshade"] = false;
	bool_default_map["write_slope"] = true;
	bool_default_map["write_aspect"] = true;
	bool_default_map["write_curvature"] = true;
	bool_default_map["write_hilltops"] = true;
	bool_default_map["write_cht"] = true;
	bool_default_map["write_hilltops"] = true;
	bool_default_map["load_previous_rasters"] = false;
	bool_default_map["print_stream_order_raster"] = false;
	bool_default_map["use_filtered_topography"] = false;

  // params for choosing which analysis you want to do
  bool_default_map["analyse_hilltops_along_swath"] = false;

	// Use the parameter parser to get the maps of the parameters required for the
	// analysis
	// load parameter parser object
	LSDParameterParser LSDPP(path_name,f_name);
	LSDPP.force_bil_extension();

	LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map,double_default_map);
	map<string,float> this_float_map = LSDPP.get_float_parameters();
	map<string,int> this_int_map = LSDPP.get_int_parameters();
	map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
	map<string,string> this_string_map = LSDPP.get_string_parameters();
  map<string,double>  this_double_map = LSDPP.get_double_parameters();

	// Now print the parameters for bug checking
	LSDPP.print_parameters();

	// location of the files
	string DATA_DIR =  LSDPP.get_read_path();
	string DEM_ID =  LSDPP.get_read_fname();
	string OUT_DIR = LSDPP.get_write_path();
	string OUT_ID = LSDPP.get_write_fname();
	string DEM_extension =  LSDPP.get_dem_read_extension();
	vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
	string CHeads_file = LSDPP.get_CHeads_file();

	LSDRaster FilledDEM;
	if (this_bool_map["load_previous_rasters"])
	{
		LSDRaster load_DEM((DATA_DIR+DEM_ID+"_Fill"), DEM_extension);
		FilledDEM = load_DEM;
	}
	else
	{
		if (this_bool_map["use_filtered_topography"])
		{
			//load dem
			LSDRaster DEM((DATA_DIR+DEM_ID+"_filtered"), DEM_extension);
		  FilledDEM = DEM.fill(this_float_map["Min slope filling"]);
		}
		else
		{
			//load dem
			LSDRaster DEM((DATA_DIR+DEM_ID), DEM_extension);
			FilledDEM = DEM.fill(this_float_map["Min slope filling"]);
		}
	}

  //surface fitting
  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(1); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(1); //plan curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

	cout << "\nCalculating surface metrics\n" << endl;
  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(this_float_map["surface_fitting_window_radius"], raster_selection);
  LSDRaster slope = Surfaces[1];
  LSDRaster aspect = Surfaces[2];
	LSDRaster curvature = Surfaces[3];

	if (this_bool_map["write_slope"])
	{
		string slope_ext = "_slope";
		slope.write_raster((DATA_DIR+DEM_ID+slope_ext), DEM_extension);
	}
	if (this_bool_map["write_aspect"])
	{
		string aspect_ext = "_aspect";
		aspect.write_raster((DATA_DIR+DEM_ID+aspect_ext), DEM_extension);
	}
	if (this_bool_map["write_curvature"])
	{
		string curv_ext = "_curv";
		curvature.write_raster((DATA_DIR+DEM_ID+curv_ext), DEM_extension);
	}

  cout << "\nGetting drainage network and basins\n" << endl;

  // get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,FilledDEM);

	// some error checking
	vector<int> sources;
	if (CHeads_file == "NULL")
	{
		cout << "I can't find your channel heads file so I'm going to use an area threshold to extract the sources" << endl;
		LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
		sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["Chan area threshold"]);
	}
	else
	{
		cout << "\t Loading the sources" << endl;
		// load the sources
		vector<int> sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv", 2);
		cout << "\t Got sources!" << endl;
	}

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "\t Got the channel network" << endl;

	if (this_bool_map["print_stream_order_raster"])
	{
		LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
		string SO_raster_name = DATA_DIR+DEM_ID+"_SO";
		SOArray.write_raster(SO_raster_name,DEM_extension);
	}

  cout << "\nExtracting hilltops and hilltop curvature" << endl;

  // extract hilltops - no critical slope filtering is performed here
  LSDRaster hilltops = ChanNetwork.ExtractRidges(FlowInfo);
	if (this_bool_map["write_hilltops"])
	{
		string ht_ext = "_hilltops";
		hilltops.write_raster((DATA_DIR+DEM_ID+ht_ext), DEM_extension);
	}

  //get hilltop curvature using filter to remove positive curvatures
  LSDRaster cht_raster = FilledDEM.get_hilltop_curvature(Surfaces[3], hilltops);
  LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(cht_raster);

	if (this_bool_map["write_cht"])
	{
		string CHT_ext = "_cht";
		CHT.write_raster((DATA_DIR+DEM_ID+CHT_ext), DEM_extension);
	}

  if (this_bool_map["analyse_hilltops_along_swath"])
  {
    // reading in the csv file with the lat long points
  	cout << "\t Reading in the csv file, filename is: " << this_string_map["coords_csv_file"] << endl;
  	LSDSpatialCSVReader SwathPoints(FilledDEM, DATA_DIR+this_string_map["coords_csv_file"]);
  	vector<float> UTME;
  	vector<float> UTMN;
  	SwathPoints.get_x_and_y_from_latlong(UTME, UTMN);
  	cout << "\t Got the x and y locations" << endl;
  	string csv_outname = "_UTM_check.csv";
  	SwathPoints.print_UTM_coords_to_csv(UTME, UTMN, (DATA_DIR+DEM_ID+csv_outname));

  	// snap to nearest channel
  	vector<int> valid_indices;
  	vector<int> snapped_nodes;
  	vector<int> snapped_JNs;
  	ChanNetwork.snap_point_locations_to_channels(UTME, UTMN, this_int_map["search_radius"], this_int_map["Threshold_SO"], FlowInfo, valid_indices, snapped_nodes, snapped_JNs);

  	cout << "The number of valid points is: " << int(valid_indices.size()) << endl;

  	if (int(valid_indices.size()) == 2)
  	{
  		// get the channel between these points
  		cout << "Got channel nodes: " << snapped_nodes[0] << ", " << snapped_nodes[1] << endl;
  		LSDIndexChannel BaselineChannel(snapped_nodes[0], snapped_nodes[1], FlowInfo);
  		vector<double> X_coords;
  		vector<double> Y_coords;
  		BaselineChannel.get_coordinates_of_channel_nodes(X_coords, Y_coords);
      // get the point data from the BaselineChannel
      PointData BaselinePoints = get_point_data_from_coordinates(X_coords, Y_coords);

      // get the swath
      cout << "\t creating swath template" << endl;
      LSDSwath TestSwath(BaselinePoints, FilledDEM, this_int_map["HalfWidth"]);

      // now get the raster values along the swath
      TestSwath.write_RasterValues_along_swath_to_csv(CHT, 0, DATA_DIR+DEM_ID+"_swath_cht_data.csv");

      string baseline_filename = "_baseline.csv";
      TestSwath.print_baseline_to_csv(FilledDEM, DATA_DIR+DEM_ID+baseline_filename);
    }
  }
}
