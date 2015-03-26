//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Longitudinal_channel_swath_erosion_profiler.cpp
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
//
// This driver file extracts the channel network and then creates a longitudinal
// swath profile of the erosion in the mainstem channel from the erosion raster 
// (from LSDCatchmentModel or CAESAR-lisflood model(s). 
//
// /// (Not yet) This version makes use of the dreich algorithm for finding the sources of the channels (Not yet)
//
// Declan Valters 
// March 2015
//
// 0) Check for both the elevation (topography) DEM and the erosion DEM
// 1) Extract the channel network from the topography DEM
// 2) Write the channel network to a points file in the LSDSwathProfile format
// 3) Use this swath template to extract the longitudinal swath profile of the erosion
// 	  variation along the channel profile.
// 4) Append this to the channel file for further plotting/analysis
//
//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDSwathProfile.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "\nFATAL ERROR: wrong number inputs. The program needs the PATH name, the PARAMETER FILE name" << endl;
		exit(EXIT_FAILURE);
	}

	string path_name = argv[1];
	
	// make sure there is a slash on the end of the file
	string lchar = path_name.substr(path_name.length()-2,1);
	string slash = "/";      
	if (lchar != slash)
	{ 
	cout << "\nYou forgot the frontslash at the end of the path. Appending." << endl; 
	path_name = path_name+slash;
	} 		
	
	// parse the parameter file name from the command line arguments
	string f_name = argv[2];

	// Printing the path and filename(s) to screen
	cout << "\nThe path is: " << path_name << endl << " and the parameter file name is: " << f_name << endl;

	string full_name = path_name + f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the parameter file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	
	string DEM_name;
	string EROSIONDEM_name;
	string raster_extension;
	string fill_ext = "_fill";
	string sources_ext = "_CH";    // for the dreich channel heads
	
	file_info_in >> DEM_name >> EROSIONDEM_name >> raster_extension;
	
	float MinSlope; 
	file_info_in >> MinSlope;
	
	int SourceThreshold;
	file_info_in >> SourceThreshold;
	
	float SwathHalfWidth;
	file_info_in >> SwathHalfWidth;
	
	float SwathBinWidth;
	file_info_in >> SwathBinWidth;
	
	file_info_in.close();
	
	cout << "\nLONGITUDINAL CHANNEL EROSION SWATH PROFILE driver function." << endl;
	cout << "This driver file extracts the channel network and then creates a longitudinal \n \
swath profile of the erosion in the mainstem channel from the erosion raster \n \
(from LSDCatchmentModel or CAESAR-lisflood model(s)." << endl;
	cout << "Author: DAV - March 2015" << endl;
	cout << "\nThe Elevation DEM is: " << DEM_name << endl;
	cout << "The Erosion amount DEM is: " << EROSIONDEM_name << endl;
	cout << "The Swath Width for averaging the erosion in the channel is: " << SwathHalfWidth << endl;
	cout << "The Source Threshold pixel number is: " << SourceThreshold << endl;

	string DEM_fill_name = path_name+DEM_name+fill_ext;

	// Load the terrain DEM
	LSDRaster topo_test((path_name+DEM_name), raster_extension);

	// Get the filled DEM
	cout << "\nFilling the DEM..." << endl;
	LSDRaster filled_topo_test = topo_test.fill(MinSlope);
	filled_topo_test.write_raster((DEM_fill_name),raster_extension);

	// Set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	// Get a flow info object for computing the channel network later on
	cout << "\nGetting FlowInfo object... " << endl;
	LSDFlowInfo FlowInfo(boundary_conditions, filled_topo_test);
    cout << "\ngot FlowInfo object." << endl;
	
	// Find the number of contributing pixels within the flowinfo object
	string CP_name =  path_name+DEM_name+"_CP";
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	ContributingPixels.write_raster(CP_name,raster_extension);	
	
	// Calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	// Get the channel head sources
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, SourceThreshold);
	cout << "\ngot the channel head sources using the SIMPLE THRESHOLD method" << endl;
	cout << "The contributing pixel threshold for channel sources was: " << SourceThreshold << endl;

	// Get the junction network
	cout << "\nInitializing Channel network..." << endl;
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "got channel_network." << endl;
	
  // Now get a junction and look for the longest channel upstream (the mainstem)
  // DAV - I test it here with a hard coded junction number since the raster is already per
  // 		clipped to a basin, so the first junction should be the outlet more or less (???)
  int junction_number = 0;
  cout << "\nCreating main stem channel..." << endl;
  LSDIndexChannel main_stem = ChanNetwork.generate_longest_index_channel_in_basin(junction_number, FlowInfo, DistanceFromOutlet);
  cout << "got main stem channel, with n_nodes " << main_stem.get_n_nodes_in_channel() <<  endl;	

	// Now create the Channel Tree, this contains basic channel profile data and will be converted in PointData later
	cout << "\nCreating IndexChannelTree..." << endl;
	int organization_switch = 1;
	int pruning_switch = 1;
	LSDIndexChannelTree ChannelTree(FlowInfo, ChanNetwork, junction_number, organization_switch,
                                        DistanceFromOutlet, pruning_switch, SourceThreshold);   // pruning_threshold = SourceThreshold (is good yes?)
  cout << "\ngot index channel tree." << endl;
    
  // Interim step of writing the channel tree to a file for the next part
	string Chan_fname = "_ChanNet";
	string Chan_ext = ".chan";
	string Chan_for_swath_ingestion_fname = path_name + DEM_name + Chan_fname + Chan_ext;
	ChannelTree.print_LSDChannels_for_chi_network_ingestion(FlowInfo,
                             filled_topo_test, DistanceFromOutlet, Chan_for_swath_ingestion_fname);
    
  cout << "\nWrote index channel tree to file: " << Chan_for_swath_ingestion_fname << endl;   
     
  // Now we are going to use a nifty tool from LSDShapeTools to write the index channel tree to a PointData object
  cout << "\nCreating PointData object from the tree file just written..." << endl;
  PointData BaselinePoints = LoadChannelTree(Chan_for_swath_ingestion_fname, 0);   // Zero argument tells it to just pick out the mainstem
    
	//-=-=-=-=-=-=-=-=-=-=-=-=//
	// NOW FOR THE COOL STUFF //
	//-=-=-=-=-=-=-=-=-=-=-=-=//
	
	// Load the raster template
  string RasterTemplate_file = DEM_name;
  string Long_Swath_ext = "_swath_long";
  string BV_ext = "_baseline_values";
  cout << "STARTING SWATHING... here we go!" << endl;
	//RasterTemplate_file = RasterTemplate_file;   // Huh?
	
  cout << "\t Loading template raster..." << endl;
  LSDRaster RasterTemplate(RasterTemplate_file.c_str(),raster_extension);
  cout << "\nRaster template loaded using: " << RasterTemplate_file << endl;

	// Load the baseline point data (created above)
	
	// Create Swath Template
  cout << "\nCreating swath template..." << endl;
  LSDSwath TestSwath(BaselinePoints, RasterTemplate, SwathHalfWidth);
  vector<float> percentiles;
  percentiles.push_back(0);
  percentiles.push_back(25);
  percentiles.push_back(50);
  percentiles.push_back(75);
  percentiles.push_back(100);
  int NormaliseTransProfile = 1;
  int NormaliseLongProfile = 0;
  
  //cout << "\n\t writing output \n\t\t - transverse profile" << endl;
  //TestSwath.write_transverse_profile_to_file(RasterTemplate, percentiles, BinWidth, RasterTemplate_file.c_str(),NormaliseTransProfile);
  
  cout << "\n Writing longitudinal profile..." << endl;
  TestSwath.write_longitudinal_profile_to_file(RasterTemplate, percentiles, SwathBinWidth, RasterTemplate_file.c_str(),NormaliseLongProfile);
  cout << "long profile swath written to file." << endl;
  
  //LSDRaster Swath(RasterTemplate.get_NRows(), RasterTemplate.get_NCols(), RasterTemplate.get_XMinimum(), RasterTemplate.get_YMinimum(),
  //                RasterTemplate.get_DataResolution(), RasterTemplate.get_NoDataValue(), TestSwath.get_DistanceToBaselineArray()); 
                   
  LSDRaster Long_Swath(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_DistanceAlongBaselineArray());
                  
  LSDRaster BaselineValues(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_BaselineValueArray());
                  
  //string output_file = RasterTemplate_file + Swath_ext;
  //Swath.write_raster(output_file.c_str(),flt_ext);
  
  string output_file2 = RasterTemplate_file + Long_Swath_ext;
  Long_Swath.write_raster(output_file2.c_str(),raster_extension);
	cout << "S.I.G." << endl;	
	
}
