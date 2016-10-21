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
	file_info_in >> DEM_name;
    file_info_in >> path_DEM;
    int max_param = 50;
	   // This array stocks the changes param: [i][0]: X, [i][1]: Y and [i][2]: value
    float param[max_param][2];
        //count the number of changes to operate
    int number_of_changes = 0;
    string test_end_file;
    
        //Reads the XYvalue params, but first test if it is the end of the file or not
    for (int i=0; i<max_param; i++){
        file_info_in >> test_end_file;
        if (test_end_file.compare("end") != 0){
            param[i][0] = atof(test_end_file.c_str());
            file_info_in >> param[i][1] >> param[i][2];
            cout << " Changes are X: " << param [i][0] << " Y: " << param [i][1] << " New value: " << param[i][2] << endl;
            number_of_changes++;
        }
        else {
            cout << "End statement detected, " << number_of_changes << " changes should be processed" << endl << endl;
            file_info_in.close();
            break;
        }
    }

	
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
    
    for (int i=0; i<number_of_changes; i++){
        nbNode = topo_index.get_node_index_of_coordinate_point (param[i][0], param[i][1]);
        topo_index.retrieve_current_row_and_col (nbNode, currentRow, currentCol);
        topo_test.set_data_element (currentRow, currentCol, param[i][2]);
        cout << "Operating change " << i << endl;
    }

    //write the output raster
    string output_name = "output2";
    topo_test.write_raster(path_DEM + output_name, "bil");


}
