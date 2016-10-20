//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// change_raster_value_driver.cpp
//
// This program should be able to to change a specific value for a raster
// 
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
    //open the file and perform i/o operation with file I guess
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
    string path_DEM;
	file_info_in >> DEM_name;
    file_info_in >> path_DEM;
	float X;
	float Y;
    float value;
	file_info_in >> X >> Y >> value;
	file_info_in.close();
	
	cout << "\nYou are running the write junctions driver." << endl
       <<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
       <<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
       <<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;

	string DEM_f_name = path_DEM+DEM_name+fill_ext;
	string DEM_bil_extension = "bil";
    
    path_name = path_DEM; 
    
  // load the DEM
  LSDRaster topo_test((path_DEM+DEM_name), DEM_bil_extension);

  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";
    
  // load the Index Raster
    LSDFlowInfo topo_index(topo_test);
  
  // remove the sea (seems to be required if gdal is used in places with nodata)
  topo_test.remove_seas();

 
    // Find the corresponding Nodes with the X and Y values  
    int nbNode;
    nbNode = topo_index.get_node_index_of_coordinate_point (X, Y);
    int currentRow = 0;
    int currentCol = 0;
    topo_index.retrieve_current_row_and_col (nbNode, currentRow, currentCol);
    
    // Set the new value
    topo_test.set_data_element (currentRow, currentCol, value);
    
    //write the output raster
    string output_name = "output";
    topo_test.write_raster(path_DEM + output_name, "bil");



}
