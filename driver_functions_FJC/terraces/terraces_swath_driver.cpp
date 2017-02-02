//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// terraces_swath_driver.cpp
// Extract information about terraces using a shapefile of the main stem channel.
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
#include "../../LSDRaster.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDTerrace.hpp"

int main (int nNumberofArgs,char *argv[])
{
	if (nNumberofArgs != 3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }
  string path_name = argv[1];
  string f_name = argv[2];
  cout << "The path name is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if (file_info_in.fail())
	{
		cout << "\nFATAL ERROR: the header file\"" << full_name
				 << "\" doesn't exist" << endl;
	}

	string DEM_ID, CH_name, Baseline_file, temp;
	int threshold_SO, FilterTopo, lower_percentile_relief, upper_percentile_relief, lower_percentile_slope, upper_percentile_slope, minimum_patch_size, search_distance;
	float surface_fitting_window_radius, Minimum_Slope, threshold_condition, RemoveChannelThreshold, HalfWidth;
	string CC_ext = "_CC";
	string DEM_extension = "bil";
	string csv_extension = "csv";

	// read in the parameters
	file_info_in >> temp >> DEM_ID
							 >> temp >> CH_name
						   >> temp >> Baseline_file
							 >> temp >> HalfWidth
							 >> temp >> threshold_SO
							 >> temp >> FilterTopo
							 >> temp >> Minimum_Slope
							 >> temp >> surface_fitting_window_radius
							 >> temp >> threshold_condition
							 >> temp >> lower_percentile_relief
 		           >> temp >> upper_percentile_relief
 							 >> temp >> lower_percentile_slope
 							 >> temp >> upper_percentile_slope
 							 >> temp >> RemoveChannelThreshold
 							 >> temp >> minimum_patch_size
 							 >> temp >> search_distance;

	file_info_in.close();

  //string Swath_ext = "_swath_trans";
  //string Long_Swath_ext = "_swath_long";
  //string BV_ext = "_baseline_values";
  cout << "starting the test run... here we go!" << endl;

	LSDRaster RasterTemplate;

	if(FilterTopo == 1)
	{
		 // load the DEM
		 cout << "Loading the DEM..." << endl;
		 LSDRaster RasterTemplate((path_name+DEM_ID), DEM_extension);

		 // filter using Perona Malik
		 int timesteps = 50;
		 float percentile_for_lambda = 90;
		 float dt = 0.1;
		 RasterTemplate = RasterTemplate.PeronaMalikFilter(timesteps, percentile_for_lambda, dt);

		 // fill
		 RasterTemplate = RasterTemplate.fill(Minimum_Slope);
		 string fill_name = "_filtered";
		 RasterTemplate.write_raster((path_name+DEM_ID+fill_name), DEM_extension);
	}
	else
	{
		//previously done the filtering and filling, just load the filled DEM
		LSDRaster load_DEM((path_name+DEM_ID+"_filtered"), DEM_extension);
		RasterTemplate = load_DEM;
	}

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	cout << "\t Flow routing..." << endl;
	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions, RasterTemplate);

	cout << "\t Loading the sources" << endl;
	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	// load the sources
	vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_ID+CH_name), csv_extension, 2);
	cout << "\t Got sources!" << endl;

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;

  cout << "\t loading baseline points" << endl;
  PointData BaselinePoints = LoadShapefile(path_name+Baseline_file.c_str());

  cout << "\t creating swath template" << endl;
  LSDSwath TestSwath(BaselinePoints, RasterTemplate, HalfWidth);
  vector<float> percentiles;
  percentiles.push_back(0);
  percentiles.push_back(25);
  percentiles.push_back(50);
  percentiles.push_back(75);
  percentiles.push_back(100);
  int NormaliseTransProfile = 1;
  int NormaliseLongProfile = 1;

	cout << "\t Getting raster from swath" << endl;
	int NormaliseToBaseline = 1;
	LSDRaster SwathRaster = TestSwath.get_raster_from_swath_profile(RasterTemplate, NormaliseToBaseline);

  // get the slope
	cout << "\t Getting the slope" << endl;
  vector<LSDRaster> surface_fitting;
  LSDRaster Slope;
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;             // this means you want the slope
  surface_fitting = RasterTemplate.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
  Slope = surface_fitting[1];

	float mask_threshold = 1.0;
	bool below = 0;
	// remove any stupid slope values
	LSDRaster Slope_new = Slope.mask_to_nodata_using_threshold(mask_threshold, below);

	// get the channel relief and slope threshold using quantile-quantile plots
	cout << "Getting channel relief threshold from QQ plots" << endl;
	string qq_fname = path_name+DEM_ID+"_qq_relief.txt";
	float relief_threshold_from_qq = SwathRaster.get_threshold_for_floodplain_QQ(qq_fname, threshold_condition, lower_percentile_relief, upper_percentile_relief);

	cout << "Getting slope threshold from QQ plots" << endl;
	string qq_slope = path_name+DEM_ID+"_qq_slope.txt";
	float slope_threshold_from_qq = Slope_new.get_threshold_for_floodplain_QQ(qq_slope, threshold_condition, lower_percentile_slope, upper_percentile_slope);

	cout << "Relief threshold: " << relief_threshold_from_qq << " Slope threshold: " << slope_threshold_from_qq << endl;

	cout << "Removing pixels within " << RemoveChannelThreshold << " m of the modern channel" << endl;
	// get the terrace pixels
	LSDTerrace Terraces(SwathRaster, Slope_new, ChanNetwork, FlowInfo, relief_threshold_from_qq, slope_threshold_from_qq, minimum_patch_size, threshold_SO, RemoveChannelThreshold);
	LSDIndexRaster ConnectedComponents = Terraces.print_ConnectedComponents_to_Raster();
	ConnectedComponents.write_raster((path_name+DEM_ID+CC_ext), DEM_extension);

	cout << "\t Testing connected components" << endl;
	vector <vector <float> > CC_vector = TestSwath.get_connected_components_along_swath(ConnectedComponents, RasterTemplate, NormaliseLongProfile);

	// push back results to file for plotting
	ofstream output_file_CC;
	string output_fname = "_terrace_swath_plots.txt";
	output_file_CC.open((path_name+DEM_ID+output_fname).c_str());
	for (int i = 0; i < int(CC_vector[0].size()); ++i)
	{
		output_file_CC << CC_vector[0][i] << " " << CC_vector[1][i] << " " << CC_vector[2][i] << endl;
	}
	output_file_CC.close();
}
