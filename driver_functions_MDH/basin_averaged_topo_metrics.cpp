//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// basin_averaged_metrics.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program gets the drainage density and topographic metricsof all catchments in the
// "DEM_name_junctions.txt" file to compare drainage density with mean hilltop
// curvature (used as a proxy for erosion rate). It gets drainage density based on a channel
// network methodology described by Pelletier (2013) and implemented by Clubb et al. (2014). 
//
// Outputs a txt file for the data cloud with the format
// catchment_ID drainage_density drainage_area mean_hilltop_curvature mean_hillslope_length mean_hillslope_gradient 
//
// Outputs a txt file for the binned data with the format
// CHT_value CHT_standard_deviation CHT_standard_error DD_value DD_standard_deviation
// DD_standard_error
//
// Developed by:
//  Fiona J. Clubb
//  Simon M. Mudd
//  Stuart W. D. Grieve
//	 Martin D. Hurst
//
// Copyright (C) 2017 Fiona Clubb and Simon M. Mudd 2016
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDBasin.hpp"
#include "../LSDCRNParameters.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the DEM name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string DEM_name = argv[2];

	cout << "The path is: " << path_name << " and the DEM name is: " << DEM_name << endl;

  // get some file names
	string fill_ext = "_Fill";
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string FileExtension = "bil";

  // Get list of junctions into a vector - this file is produced by ./drainage_density_step1_junctions.out
  string ListOfJunctions = DEM_name+"_junctions.list";
  vector<int> junction_list;
  int this_junction;
  string line;
  ifstream junction_input(ListOfJunctions.c_str());
  if (junction_input.is_open())
  {
      while (junction_input >> this_junction)
      {
          junction_list.push_back(this_junction);
      }
      junction_input.close();
  }

  int no_junctions = int(junction_list.size());
  cout << "Number of basins: " << no_junctions << endl;

	// load in rasters

	// get the filled file
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), FileExtension);

	// Set the no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	//get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);

	//create stream network from mapped sources
	cout << "\t Loading Sources..." << endl;
	vector<int> sources;
	string CHeads_file = path_name+DEM_name+"_CH_wiener_nodeindices_for_Arc";
	string CHeads_file_ext = CHeads_file + ".csv";
	if (ifstream(CHeads_file_ext.c_str()))
	{
		cout << "Loading channel heads from the file: " << CHeads_file << endl;
		sources = FlowInfo.Ingest_Channel_Heads(CHeads_file, "csv",2);
		cout << "\t Got sources!" << endl;
	}
	else
	{
	 	cout << endl << endl << endl << "==================================" << endl;
	 	cout << "The channel head file is null. " << endl;
		exit(EXIT_FAILURE);
	}

	// now get the junction network
	cout << "Initializing Channel network" << endl;
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "got channel_network" << endl;

	//Check to see if a list of junctions for extraction exists
	string JunctionsFile = path_name+DEM_name+"_junctions.list";
	vector<int> JunctionsList;
	ifstream infile(JunctionsFile.c_str());
	if (infile)
	{
		cout << "Junctions File " << JunctionsFile << " exists" << endl;;
		int n;
		while (infile >> n) JunctionsList.push_back(n);
	}
	else 
	{
		cout << "Junctions File " << JunctionsFile << " does not exist" << endl;
		exit(EXIT_FAILURE);
	}

	// Get a raster with all basins
	LSDIndexRaster Basin_Raster = ChanNetwork.extract_basins_from_junction_vector(junction_list, FlowInfo);
	string basin_name = "_basins";
	Basin_Raster.write_raster((path_name+DEM_name+basin_name), FileExtension);

	//Load the rasters from the hillslope analysis that we want to analyse
  	//get the CHT raster
  	string CHT_ext = "_CHT";
	LSDRaster CHT((path_name+DEM_name+CHT_ext), FileExtension);
	
	//get the LH raster
	string LH_ext = "_HFR_LH";
	LSDRaster LH((path_name+DEM_name+LH_ext), FileExtension);

	//get the S raster
	string SLP_ext = "_HFR_SLP";
	LSDRaster SLP((path_name+DEM_name+SLP_ext), FileExtension);

	//get the MChi raster
	string MChi_ext = "_MChi";
	LSDRaster MChi((path_name+DEM_name+MChi_ext), FileExtension);

	// get the SO array
	string SO_ext = "_SO";
	LSDIndexRaster SOArray((path_name+DEM_name+SO_ext), FileExtension);

	//get metrics for each basin for plotting
	float DD, Mean, Median, SD, SE, P5th, P16th, P84th, P95th, BinWidth;
	int Percentile, junction_number;
	vector<float> DataVector;
	string OutputFile, Hist_Filename;

	OutputFile = DEM_name+"_basin_topo_metrics.txt";
  	ofstream OutputStream;
	OutputStream.open(OutputFile.c_str());
	OutputStream 	<< "ID DD CHT_Mean CHT_Median CHT_SD CHT_SE CHT_P5th CHT_P16th CHT_P84th CHT_P95th "
						<< "LH_Mean LH_Median LH_SD LH_SE LH_P5th LH_P16th LH_P84th LH_P95th "
						<< "SLP_Mean SLP_Median SLP_SD SLP_SE SLP_P5th SLP_P16th SLP_P84th SLP_P95th "
						<< "MChi_Mean MChi_Median MChi_SD MChi_SE MChi_P5th MChi_P16th MChi_P84th MChi_P95th" << endl;

	for (int i = 0; i < int(junction_list.size()); i++)
	{
		junction_number = junction_list[i];
		cout << flush << "Junction is " << junction_number << ";  " << i+1 << "/" << junction_list.size() << "\r";
		
		// set basin parameters
		LSDBasin Basin(junction_number, FlowInfo, ChanNetwork);
		Basin.set_FlowLength(SOArray, FlowInfo);
		Basin.set_DrainageDensity();
		Basin.set_Perimeter(FlowInfo);
		
		// get the drainage density of the catchment
		DD = Basin.get_DrainageDensity();
		
		//write to file
		OutputStream << junction_number << " " << DD << " ";

		// hilltop curvature
		LSDRaster CHT_Basin = Basin.write_raster_data_to_LSDRaster(CHT, FlowInfo);
		Mean = Basin.CalculateBasinMean(FlowInfo,CHT_Basin);
		Median = Basin.CalculateBasinMedian(FlowInfo,CHT_Basin);
		SD = Basin.CalculateBasinStdDev(FlowInfo,CHT_Basin);
		SE = Basin.CalculateBasinStdError(FlowInfo,CHT_Basin);
		Percentile = 5;
		P5th = Basin.CalculateBasinPercentile(FlowInfo,CHT_Basin,Percentile);
		Percentile = 16;
		P16th = Basin.CalculateBasinPercentile(FlowInfo,CHT_Basin,Percentile);
		Percentile = 84;
		P84th = Basin.CalculateBasinPercentile(FlowInfo,CHT_Basin,Percentile);
		Percentile = 95;
		P95th = Basin.CalculateBasinPercentile(FlowInfo,CHT_Basin,Percentile);

		//write to file
		OutputStream << Mean << " " << Median << " " << SD << " " << SE << " " << P5th << " " << P16th << " " << P84th << P95th << " ";

		//Generate histogram
		DataVector = CHT_Basin.get_RasterData_vector_No_NDVs();
		BinWidth = 0.001;
		Hist_Filename = "Basin_"+string(itoa(junction_number))+"_CHT_hist.txt";
		print_histogram(DataVector, BinWidth, Hist_Filename);

		//Hillslope length
		LSDRaster LH_Basin = Basin.write_raster_data_to_LSDRaster(LH, FlowInfo);
		Mean = Basin.CalculateBasinMean(FlowInfo,LH_Basin);
		Median = Basin.CalculateBasinMedian(FlowInfo,LH_Basin);
		SD = Basin.CalculateBasinStdDev(FlowInfo,LH_Basin);
		SE = Basin.CalculateBasinStdError(FlowInfo,LH_Basin);
		Percentile = 5;
		P5th = Basin.CalculateBasinPercentile(FlowInfo,LH_Basin,Percentile);
		Percentile = 15;
		P16th = Basin.CalculateBasinPercentile(FlowInfo,LH_Basin,Percentile);
		Percentile = 84;
		P84th = Basin.CalculateBasinPercentile(FlowInfo,LH_Basin,Percentile);
		Percentile = 95;
		P95th = Basin.CalculateBasinPercentile(FlowInfo,LH_Basin,Percentile);

		//write to file
		OutputStream << Mean << " " << Median << " " << SD << " " << SE << " " << P5th << " " << P16th << " " << P84th << P95th << " ";

		//Generate histogram
		DataVector = LH_Basin.get_RasterData_vector_No_NDVs();
		BinWidth = 1.;
		Hist_Filename = "Basin_"+string(itoa(junction_number))+"_LH_hist.txt";
		print_histogram(DataVector, BinWidth, Hist_Filename);

		// Slope
		LSDRaster SLP_Basin = Basin.write_raster_data_to_LSDRaster(SLP, FlowInfo);
		Mean = Basin.CalculateBasinMean(FlowInfo,SLP_Basin);
		Median = Basin.CalculateBasinMedian(FlowInfo,SLP_Basin);
		SD = Basin.CalculateBasinStdDev(FlowInfo,SLP_Basin);
		SE = Basin.CalculateBasinStdError(FlowInfo,SLP_Basin);
		Percentile = 5;
		P5th = Basin.CalculateBasinPercentile(FlowInfo,SLP_Basin,Percentile);
		Percentile = 16;
		P16th = Basin.CalculateBasinPercentile(FlowInfo,SLP_Basin,Percentile);
		Percentile = 84;
		P84th = Basin.CalculateBasinPercentile(FlowInfo,SLP_Basin,Percentile);
		Percentile = 95;
		P95th = Basin.CalculateBasinPercentile(FlowInfo,SLP_Basin,Percentile);

		//write to file
		OutputStream << Mean << " " << Median << " " << SD << " " << SE << " " << P5th << " " << P16th << " " << P84th << P95th << " ";

		//Generate histogram
		DataVector = SLP_Basin.get_RasterData_vector_No_NDVs();
		BinWidth = 0.01;
		Hist_Filename = "Basin_"+string(itoa(junction_number))+"_SLP_hist.txt";
		print_histogram(DataVector, BinWidth, Hist_Filename);

		// Chi steepness
		LSDRaster MChi_Basin = Basin.write_raster_data_to_LSDRaster(MChi, FlowInfo);
		Mean = Basin.CalculateBasinMean(FlowInfo,MChi_Basin);
		Median = Basin.CalculateBasinMedian(FlowInfo,MChi_Basin);
		SD = Basin.CalculateBasinStdDev(FlowInfo,MChi_Basin);
		SE = Basin.CalculateBasinStdError(FlowInfo,MChi_Basin);
		Percentile = 5;
		P5th = Basin.CalculateBasinPercentile(FlowInfo,MChi_Basin,Percentile);
		Percentile = 16;
		P16th = Basin.CalculateBasinPercentile(FlowInfo,MChi_Basin,Percentile);
		Percentile = 84;
		P84th = Basin.CalculateBasinPercentile(FlowInfo,MChi_Basin,Percentile);
		Percentile = 95;
		P95th = Basin.CalculateBasinPercentile(FlowInfo,MChi_Basin,Percentile);

		//write to file
		OutputStream << Mean << " " << Median << " " << SD << " " << SE << " " << P5th << " " << P16th << " " << P84th << P95th << endl;
		
		//Generate histogram
		DataVector = MChi_Basin.get_RasterData_vector_No_NDVs();
		BinWidth = 0.01;
		Hist_Filename = "Basin_"+string(itoa(junction_number))+"_MChi_hist.txt";
		print_histogram(DataVector, BinWidth, Hist_Filename);

	}    
}
