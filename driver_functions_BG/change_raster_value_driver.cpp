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
  
  
  cout << "\nYou are running the write junctions driver." << endl
     <<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
     <<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
     <<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;
  
  
  //open the file and perform i/o operation 
	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
  
    //check the file
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the parameter file \"" << full_name
		   << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}
	
    //Load parameters  
	string DEM_name;
	string path_DEM;
  string output_name;
	file_info_in >> DEM_name;
  file_info_in >> path_DEM;
  file_info_in >> output_name;
  	
  
  
  // initiate variables
  int row, col;
  float value;
  vector<int> rows;
  vector<int> cols;
  vector<float> values;

    // read the file
  while(file_info_in >> row >> col >> value)
  {
      rows.push_back(row);
      cols.push_back(col);
      values.push_back(value);
  }
  file_info_in.close();

  // now print the results
  int n_data_points = int(rows.size());
  for (int i = 0; i < n_data_points; i++)
  {
    cout << "i: " << i << " row: " << rows[i] << " col: " << cols[i] << " value: " << values[i] << endl;
  }
  cout << " There are " << n_data_points << " changes to process" << endl;
  	
	string DEM_bil_extension = "bil";
    
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

 
  // Find the corresponding Nodes (pixel) with the X and Y coordinates, then with the rows and column to finally change the raster value  
  int nbNode;
  int currentRow = 0;
  int currentCol = 0;
  int count_temp = 0;
  
  
  
  
  
  for (int i=0; i< n_data_points; i++){
    nbNode = topo_index.get_node_index_of_coordinate_point (rows[i], cols[i]);
    topo_index.retrieve_current_row_and_col (nbNode, currentRow, currentCol);
    topo_test.set_data_element (currentRow, currentCol, values[i]);
    count_temp = i+1;
    cout << "Operating change " << count_temp << endl;
  }

  //write the output raster
  topo_test.write_raster(path_DEM + output_name, "bil");


}
