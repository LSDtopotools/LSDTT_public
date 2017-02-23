//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// run_typology_model.cpp
//
// This program runs the river typology model.  Firstly it fills the DEM and then smooths
// it based on the polynomial surface fitting.  It then generates a channel network based
// on an area threshold and splits the channel network into segments with a unique ID for
// each segment.  The segment length is calculated as:
// Length = MinReachLength * sqrt(DrainageArea (km))
// where the MinReachLength is set by the user.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb
// University of Edinburgh
// 06/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <ctime>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDSpatialCSVReader.hpp"
#include "../../TNT/tnt.h"

int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();

	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the and the parameter file name." << endl;
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

	string DEM_name, csv_filename;
	float surface_fitting_window_radius, Minimum_Slope;
	int AreaThreshold, MinReachLength;
	string temp;
	file_info_in >> temp >> DEM_name
							 >> temp >> Minimum_Slope
							 >> temp >> AreaThreshold
							 >> temp >> surface_fitting_window_radius
							 >> temp >> MinReachLength
							 >> temp >> csv_filename;

	file_info_in.close();

	cout << "You are running the typology model analysis with the following settings:" << endl;
	cout << "\t DEM name: " << DEM_name << endl;
	cout << "\t Minimum slope for filling: " << Minimum_Slope << endl;
	cout << "\t Area threshold for channel extraction: " << AreaThreshold << endl;
	cout << "\t Surface fitting radius: " << surface_fitting_window_radius << endl;
	cout <<	"\t Minimum reach length: " << MinReachLength << endl;
	cout << "\t Baseline CSV filename: " << csv_filename << endl;

	string raster_output = "_test";
	string elev_output = "_elev";
	string slope_output = "_slope";
	string HS_output = "_HS";
  string DEM_extension = "bil";
	string fill_ext = "_fill";
	string discharge_ext = "_Qmed";
	vector<int> raster_selection(8,0);
	raster_selection[0] = 1;  // get the smoothed elevation
	raster_selection[1] = 1;  // get the slope

	// load the raster and remove values below 0
	// LSDRaster dem(path_name+DEM_name, DEM_extension);
	// dem.remove_seas();
	// //dem.write_raster((path_name+DEM_name+raster_output), DEM_extension);
	//
	// //====================================================================//
	// // 							              	SMOOTHING                             //
	// //====================================================================//
	//
	// cout << "\t Running polyfitting..." << endl;
	//
	// vector<LSDRaster> output_rasters;
	// output_rasters = dem.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
	//
	// // smoothed elevation
	// output_rasters[0].write_raster((path_name+DEM_name+elev_output), DEM_extension);
	// // slope
	// output_rasters[1].write_raster((path_name+DEM_name+slope_output), DEM_extension);
	//
	// // write smoothed hillshade
	// LSDRaster HS = output_rasters[0].hillshade(45, 315, 1);
	// HS.write_raster((path_name+DEM_name+elev_output+HS_output), DEM_extension);
	//
	// cout << "\t Filling the DEM..." << endl;
	// // fill
	// LSDRaster filled_DEM = output_rasters[0].fill(Minimum_Slope);
	// filled_DEM.write_raster((path_name+DEM_name+fill_ext), DEM_extension);
	// cout << "Got the filled DEM" << endl;
	LSDRaster filled_DEM((path_name+DEM_name+fill_ext), DEM_extension);

	// load the discharge raster
	LSDRaster Discharge((path_name+DEM_name+discharge_ext), DEM_extension);

	// buffer the discharge raster
	float buffer_radius = 150;
	LSDRaster BufferedQ = Discharge.BufferRasterData(buffer_radius);
	string buffer_ext = "_Qmed_buffer";
	BufferedQ.write_raster((path_name+DEM_name+buffer_ext), DEM_extension);

	//====================================================================//
	// 							      	CHANNEL NETWORK EXTRACTION                    //
	//====================================================================//

	cout << "\t Getting channel network..." << endl;
	// Set the no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	//get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_DEM);
	// get the sources
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	vector<int> sources = FlowInfo.get_sources_index_threshold(ContributingPixels, AreaThreshold);

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

	// print the channel network to raster
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name = "_SO";
	SOArray.write_raster((path_name+DEM_name+SO_name), DEM_extension);

	cout << "\t Reading in the baseline channel network" << endl;
	vector<int> NIs;
	vector<float> X_coords;
	vector<float> Y_coords;
	//FlowInfo.get_nodeindices_from_csv((path_name+csv_filename), NIs, X_coords, Y_coords);

	LSDSpatialCSVReader BaselineSources(filled_DEM, (path_name+csv_filename));
	BaselineSources.get_nodeindices_from_x_and_y_coords(FlowInfo, X_coords, Y_coords, NIs);

	// get the catch_ID and hydro_code from the baseline sources
	string catch_name = "CATCH_";
	string hydro_code_name = "HYDRO_CODE";
	vector<int> CatchIDs = BaselineSources.data_column_to_int(catch_name);
	vector<int> HydroCodes = BaselineSources.data_column_to_int(hydro_code_name);

	int search_radius_nodes = 25;
	int threshold_SO = 1;
	vector<int> valid_points;
	vector<int> snapped_NIs;
	vector<int> snapped_JNs;
	ChanNetwork.snap_point_locations_to_channels(X_coords, Y_coords, search_radius_nodes, threshold_SO, FlowInfo, valid_points, snapped_NIs, snapped_JNs);

	string snapped_file = "_snapped";
	FlowInfo.print_vector_of_nodeindices_to_csv_file(snapped_NIs, (path_name+DEM_name+snapped_file));

	cout << "\t Now splitting the channel into segments..." << endl;
	// now test the split channel function
	LSDIndexRaster ChannelSegments;
	vector < vector<int> > SegmentInfoInts;
	vector < vector<float> > SegmentInfoFloats;
	int search_radius = 25;
	ChanNetwork.TypologyModel(FlowInfo, sources, snapped_NIs, MinReachLength, search_radius, filled_DEM, BufferedQ, ChannelSegments, SegmentInfoInts, SegmentInfoFloats);
	string segment_ext = "_segments";
	ChannelSegments.write_raster((path_name+DEM_name+segment_ext), DEM_extension);

	// removing segments not downstream of baseline sources
	//cout << "\t Removing segments not downstream of baseline sources..." << endl;
	//ChanNetwork.remove_tributary_segments(FlowInfo, snapped_NIs, SegmentInfoInts, SegmentInfoFloats);

	//print to csv file
	cout << "\t Writing csv file of channel segments..." << endl;
	ChanNetwork.print_channel_segments_to_csv(FlowInfo, SegmentInfoInts, SegmentInfoFloats, (path_name+DEM_name));

	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end-begin) / CLOCKS_PER_SEC;
	cout << "Done, time taken (secs): " << elapsed_secs << endl;
}
