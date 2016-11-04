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
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"


int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	
  // make sure there is a slash on the end of the file
  string lchar = path_name.substr(path_name.length()-2,1);
  string slash = "/";      
  if (lchar != slash)
  { 
    cout << "You forgot the frontslash at the end of the path. Appending." << endl; 
    path_name = path_name+slash;
  } 		
	
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
	string sources_name;
	string fill_ext = "_fill";
	file_info_in >> DEM_name >> sources_name;
	file_info_in.close();
	
	cout << "\nYou are running the write junctions driver." << endl
       <<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
       <<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
       <<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;

	string DEM_f_name = DEM_name+fill_ext;
	string DEM_bil_extension = "bil";

  // load the DEM
  LSDRaster filled_topo_test((path_name+DEM_f_name), DEM_bil_extension);

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);

  cout << "got FlowInfo" << endl;

	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	cout << "\t Loading Sources..." << endl;
	// load the sources
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_name+sources_name), DEM_bil_extension, 1);
  cout << "\t Got sources!" << endl;

	// now get the junction network
	cout << "Initializing Channel network" << endl;
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "got channel_network" << endl;

	// get the stream orders and the junctions
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
	string JI_name = "_JI";
	JIArray.write_raster((path_name+DEM_name+JI_name),DEM_bil_extension);

}
