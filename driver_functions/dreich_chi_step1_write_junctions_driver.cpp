//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Modification of the get junctions to use the DrEICH algorithm to create chi 
// profiles rather than using an area threshold.
//
// NOTE: if using the same driver file as for the original channel writing driver 
//(chi_step2_write_channel_file_driver.cpp) then you MUST delete the line with the threshold
// area (line 3 in previous parameter file)
// NOTE#2: I have moved the declaration of the input raster type into the driver file - my DEMs
// are in .flt format rather than .bil.  The format of the DEM is now declared in line 2 of the
// driver file without the first "." e.g. "flt" or "bil"
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
	string raster_extension;
	string fill_ext = "_fill";
	string sources_ext = "_CH";
	file_info_in >> DEM_name >> raster_extension;
	float MinSlope;
	file_info_in >> MinSlope;
	file_info_in.close();
	
	cout << "\nYou are running the write junctions driver." << endl;
       //<<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
       //<<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
       //<<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;

	string DEM_f_name = path_name+DEM_name+fill_ext;
	//string DEM_bil_extension = "bil";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), raster_extension);

	// get the filled file
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(MinSlope);
	filled_topo_test.write_raster((DEM_f_name),raster_extension);

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// make the hillshade (this is faster than doing it in arc
	string HS_name = "_HS";
	LSDRaster HS = filled_topo_test.hillshade(45, 315, 1);
	HS.write_raster((path_name+DEM_name+HS_name),raster_extension);
	
	cout << "NRows: " << filled_topo_test.get_NRows() << endl;
	

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);

  cout << "got FlowInfo" << endl;


  string CP_name =  path_name+DEM_name+"_CP";
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	ContributingPixels.write_raster(CP_name,raster_extension);

	//string FI_fname = "_flowinfo";
	//FlowInfo.pickle((path_name+DEM_name+FI_fname));



	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	// get the sources
	vector<int> sources;
	sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_name+sources_ext),raster_extension);

	// now get the junction network
	cout << "Initializing Channel network" << endl;
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	cout << "got channel_network" << endl;

	// get the stream orders and the junctions
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();

	string SO_name = "_SO_dreich";
	string JI_name = "_JI_dreich";

	SOArray.write_raster((path_name+DEM_name+SO_name),raster_extension);
	JIArray.write_raster((path_name+DEM_name+JI_name),raster_extension);


}
