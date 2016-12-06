//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// write_junctions_driver.cpp
//
// This program is used to write a junction file
// The junction file can be used to select a junction for a channel tree
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Simon M. Mudd
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
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"


int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
    exit(EXIT_SUCCESS);
  }

  // make sure there is a slash on the end of the file
  string path_name = argv[1];
  path_name = FixPath(path_name);
	
	
	string f_name = argv[2];

	cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the parameter file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	

	string DEM_name;
	string fill_ext = "_fill";
	file_info_in >> DEM_name;
	float MinSlope;
	int threshold_contributing_pixels;
	string CHeads_file = "Null";
	file_info_in >> MinSlope >> threshold_contributing_pixels;
	file_info_in.close();
	
	cout << "\nYou are running the write junctions driver." << endl
       <<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
       <<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
       <<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;

	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_bil_extension = "bil";

  // load the DEM
  LSDRaster topo_test((path_name+DEM_name), DEM_bil_extension);
  
  // remove the sea (seems to be required if gdal is used in places with nodata)
  topo_test.remove_seas();

	// get the filled file
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(MinSlope);
	filled_topo_test.write_raster((DEM_f_name),DEM_bil_extension);

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// make the hillshade (this is faster than doing it in arc
	string HS_name = "_HS";
	LSDRaster HS = filled_topo_test.hillshade(45, 315, 1);
	HS.write_raster((path_name+DEM_name+HS_name),DEM_bil_extension);
	
	cout << "NRows: " << filled_topo_test.get_NRows() << endl;
	

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
	cout << "got FlowInfo" << endl;


	string CP_name =  path_name+DEM_name+"_CP";
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	ContributingPixels.write_raster(CP_name,DEM_bil_extension);

	// Write flow info to file
	string FI_fname = "_flowinfo";
	FlowInfo.pickle((path_name+DEM_name+FI_fname));

	// calcualte the distance from outlet
	cout << "\t Calculating flow distance..." << endl;
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();


	//create stream network from mapped sources
	cout << "\t Loading Sources..." << endl;
	vector<int> sources;
	CHeads_file = path_name+DEM_name+"_CH_wiener_nodeindices_for_Arc";

	if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  	{
   	cout << endl << endl << endl << "==================================" << endl;
   	cout << "The channel head file is null. " << endl;
		cout << "Getting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
   	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold_contributing_pixels);
   
		cout << "The number of sources is: " << sources.size() << endl;
	}
	else
	{
		cout << "Loading channel heads from the file: " << path_name+CHeads_file << endl;
		sources = FlowInfo.Ingest_Channel_Heads(path_name+CHeads_file, "csv",2);
		cout << "\t Got sources!" << endl;
	}

	// now get the junction network
	cout << "Initializing Channel network" << endl;
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "got channel_network" << endl;

	// get the stream orders and the junctions
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();

	string SO_name = "_SO";
	string JI_name = "_JI";

	SOArray.write_raster((path_name+DEM_name+SO_name),DEM_bil_extension);
	JIArray.write_raster((path_name+DEM_name+JI_name),DEM_bil_extension);

}
