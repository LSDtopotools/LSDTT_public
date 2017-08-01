//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// get_swath_from_latlong.cpp
// Takes in a csv file with latitude and longitude and gets a swath profile
// along a channel between the points.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../../LSDRaster.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDSpatialCSVReader.hpp"
#include "../../LSDParameterParser.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();

	//Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the swath tool!  					              ||" << endl;
    cout << "|| This program takes in a csv file with latitude and  ||" << endl;
		cout << "|| longitude and gets a swath profile along            ||" << endl;
		cout << "|| a channel between the points												||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Fiona J. Clubb												              ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./swath_profiler.out /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_Swath.param" << endl;
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

	// set default int parameters
	int_default_map["search_radius"] = 50;
	int_default_map["HalfWidth"] = 500;
	int_default_map["NormaliseToBaseline"] = 0;
	int_default_map["threshold_contributing_pixels"] = 1000;

	// set default string parameters
	string_default_map["csv_file"] = "Merapi_swath_points.csv";
	string_default_map["CHeads_format"] = "csv";

	bool_default_map["trim_raster"] = true;

	// set default float parameters
	float_default_map["Minimum_Slope"] = 0.0001;

	// Use the parameter parser to get the maps of the parameters required for the
	// analysis
	// load parameter parser object
	LSDParameterParser LSDPP(path_name,f_name);
	LSDPP.force_bil_extension();

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
	string raster_ext =  LSDPP.get_dem_read_extension();
	vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
	string CHeads_file = LSDPP.get_CHeads_file();

  //string Swath_ext = "_swath_trans";
  //string Long_Swath_ext = "_swath_long";
  //string BV_ext = "_baseline_values";
  cout << "starting the test run... here we go!" << endl;

	// load the DEM
	cout << "Loading the DEM..." << endl;
	LSDRaster RasterTemplate((DATA_DIR+DEM_ID), raster_ext);

	cout << "Filling..." << endl;
	LSDRaster FilledRaster = RasterTemplate.fill(this_float_map["Minimum_Slope"]);

	cout << "\t Flow routing..." << endl;
	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions, FilledRaster);

	// reading in the csv file with the lat long points
	cout << "Reading in the csv file" << endl;
	LSDSpatialCSVReader SwathPoints(FilledRaster, DATA_DIR+this_string_map["csv_file"]);
	vector<float> UTME;
	vector<float> UTMN;
	SwathPoints.get_x_and_y_from_latlong(UTME, UTMN);
	cout << "Got the x and y locations" << endl;
	string csv_outname = "_UTM_check.csv";
	SwathPoints.print_UTM_coords_to_csv(UTME, UTMN, (DATA_DIR+DEM_ID+csv_outname));

	// get the channel network
	vector<int> sources;
	if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {

		int threshold_contributing_pixels = this_int_map["threshold_contributing_pixels"];
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
		// calculate the flow accumulation
	  cout << "\t Calculating flow accumulation (in pixels)..." << endl;

	  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, threshold_contributing_pixels);

    cout << "The number of sources is: " << sources.size() << endl;
  }
  else
  {
		cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
		sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), this_string_map["CHeads_format"], 2);
		cout << "\t Got sources!" << endl;

  }


	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	// snap to nearest channel
	int threshold_SO = 2;
	vector<int> valid_indices;
	vector<int> snapped_nodes;
	vector<int> snapped_JNs;
	ChanNetwork.snap_point_locations_to_channels(UTME, UTMN, this_int_map["search_radius"], threshold_SO, FlowInfo, valid_indices, snapped_nodes, snapped_JNs);

	cout << "The number of valid points is: " << int(valid_indices.size()) << endl;

	if (int(valid_indices.size()) == 2)
	{
		// get the channel between these points
		cout << "Got channel nodes: " << snapped_nodes[0] << ", " << snapped_nodes[1] << endl;
		LSDIndexChannel BaselineChannel(snapped_nodes[0], snapped_nodes[1], FlowInfo);
		//BaselineChannel.write_channel_to_csv(DATA_DIR,DEM_ID+"_chan_check");
		vector<double> X_coords;
		vector<double> Y_coords;
		BaselineChannel.get_coordinates_of_channel_nodes(X_coords, Y_coords);
		// get the point data from the BaselineChannel
		PointData BaselinePoints = get_point_data_from_coordinates(X_coords, Y_coords);

		// get the swath
		cout << "\t creating swath template" << endl;
		LSDSwath TestSwath(BaselinePoints, FilledRaster, this_int_map["HalfWidth"]);

		cout << "\n\t Getting raster from swath" << endl;
		LSDRaster SwathRaster = TestSwath.get_raster_from_swath_profile(FilledRaster, this_int_map["NormaliseToBaseline"]);
		string swath_ext = "_swath_raster";
		if(this_bool_map["trim_raster"])
		{
			// Reducing the Raster
      cout << "I am now trying to reduce your raster, by removing part the NoDataValue" << endl;
      LSDRaster LightRaster = SwathRaster.RasterTrimmerPadded(50);
			LightRaster.write_raster((DATA_DIR+DEM_ID+swath_ext), raster_ext);
		}
		else
		{
			SwathRaster.write_raster((DATA_DIR+DEM_ID+swath_ext), raster_ext);
		}

		// get the raster values along the swath
		vector <vector <float> > ElevationValues = TestSwath.get_RasterValues_along_swath(RasterTemplate, this_int_map["NormaliseToBaseline"]);

		// push back results to file for plotting
		ofstream output_file;
		string output_fname = "_swath_elevations.csv";
		output_file.open((path_name+DEM_ID+output_fname).c_str());
		output_file << "Distance,Mean,Min,Max" << endl;
		for (int i = 0; i < int(ElevationValues[0].size()); ++i)
		{
			output_file << ElevationValues[0][i] << "," << ElevationValues[1][i] << "," << ElevationValues[2][i] << "," << ElevationValues[3][i] << endl;
		}
		output_file.close();

	}
	else
	{
		cout << "I was unable to find a channel between those coordinates! Check your coordinates or increase the search radius." << endl;
	}

	// // get the swath
  // cout << "\t creating swath template" << endl;
  // LSDSwath TestSwath(BaselinePoints, RasterTemplate, HalfWidth);
	//
	// cout << "\n\t Getting raster from swath" << endl;
	// LSDRaster SwathRaster = TestSwath.get_raster_from_swath_profile(RasterTemplate, NormaliseToBaseline);
	// string swath_ext = "_swath_raster";
	// SwathRaster.write_raster((DATA_DIR+DEM_ID+swath_ext), raster_ext);
	//
  // // get the raster values along the swath
	// vector <vector <float> > ElevationValues = TestSwath.get_RasterValues_along_swath(RasterTemplate, NormaliseToBaseline);
	//
	// // push back results to file for plotting
	// ofstream output_file;
	// string output_fname = "_swath_elevations.csv";
	// output_file.open((DATA_DIR+DEM_ID+output_fname).c_str());
	// output_file << "Distance,Mean,Min,Max" << endl;
	// for (int i = 0; i < int(ElevationValues[0].size()); ++i)
	// {
	// 	output_file << ElevationValues[0][i] << "," << ElevationValues[1][i] << "," << ElevationValues[2][i] << "," << ElevationValues[3][i] << endl;
	// }
	// output_file.close();
	//
	// // Done, check how long it took
	// clock_t end = clock();
	// float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	// cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
