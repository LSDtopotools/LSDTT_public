//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// compare_to_field_data
//
// This program takes three arguments, the path name, the file name for the XY coords, and // the DEM name WITHOUT EXTENSION.
// It takes a list of XY coordinates for mapped floodplain initiation points.  It snaps the // XY coordinate to the nearest channel and outputs the mean snapping distance.  It then
// finds the nearest connected components floodplain patch to the field mapped FIP, and 
// finds the furthest upstream pixel from this point.  It then calculates the distance
// between this point and the mapped FIP.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 08/09/16
// Fiona J. Clubb
// University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../../LSDStatsTools.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDShapeTools.hpp"


int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();
	
	//Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the coordinate file name, and the DEM name WITHOUT EXTENSION" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];
	string DEM_name = argv[3];

	cout << "The path is: " << path_name << ", the filename is: " << f_name << ", and the DEM name is: " << DEM_name << endl;

	string full_name = path_name+f_name;

	// Open the file with the coordianates
	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the coordinates file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	
	vector<string> FIP_IDs;
	vector<float> X_coords;
	vector<float> Y_coords;
	string FIP_ID;
	float X_coord, Y_coord;
	
	// Read in the coordinates
	if (file_info_in.is_open())
	{
		while (file_info_in >> FIP_ID >> X_coord >> Y_coord)
		{
			cout << "FIP ID is: " << FIP_ID << " X: " << X_coord << " Y: " << Y_coord << endl;
			FIP_IDs.push_back(FIP_ID);
			X_coords.push_back(X_coord);
			Y_coords.push_back(Y_coord);
		}
	}
	file_info_in.close();
	
	// declare some useful strings
	string input_path;
  string DEM_extension = "bil";
	string fill_ext = "_fill";
	string sources_ext = "_CH_wiener";

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	// load the filled DEM
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_extension);
	// get the flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
														 
	// get the channel network
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_name+sources_ext),DEM_extension);
  
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	
	// Load in the floodplain raster
	string CC_name = "_FP";
	LSDRaster ConnectedComponents((path_name+DEM_name+CC_name), DEM_extension);
	
	// Snap each point to the nearest channel
	cout << "Snapping points to the nearest channel" << endl;
														 
	vector<float> NodeIndices;													 
  int threshold_SO = 2;
  int search_radius = 20;
	int search_distance_FIP = 200;
	
	// output file
	ofstream output_file;
	string output_ext = "_error_distances.txt";
	string filename = path_name+DEM_name+output_ext;
	cout << "Filename is: " << filename << endl;
	output_file.open(filename.c_str());
	
  for (int i = 0; i < int(X_coords.size()); i++)
  {
    //cout << flush << "Snapped = " << i+1 << " of " << X_coords.size() << "\r";
    int NI = ChanNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(X_coords[i], Y_coords[i], search_radius, threshold_SO, FlowInfo);
		NodeIndices.push_back(NI);
		cout << "ID: " << FIP_IDs[i] << endl;
		//get the distance to the predicted floodplain initiation point
		float distance = ChanNetwork.find_distance_to_nearest_floodplain_pixel(NI, search_distance_FIP, ConnectedComponents, FlowInfo);
		cout << "Distance: " << distance << endl;
		//write to file
		output_file << FIP_IDs[i] << " " << distance << endl;
  }     
	output_file.close();
	
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
