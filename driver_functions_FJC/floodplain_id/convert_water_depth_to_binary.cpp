//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// convert_water_depth_to_binary
//
// This program takes two arguments, the path name and the driver name
// It takes in a water depth raster created by LSDCatchmentModel (in .asc format)
// and outputs a binary raster of flood inundation using a small positive threshold which
// is set by the user in the driver file.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 30/05/16
// Fiona J. Clubb
// University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
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
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];

	cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the header file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}

	string DEM_name;
  string DEM_extension = "asc";
	string bil_extension = "bil";
	string temp;
	float threshold;	
  
	// read in the parameters
	file_info_in >> temp >> DEM_name
							 >> temp >> threshold;
	
	cout << "The parameters are: DEM name = " << DEM_name << endl;
	cout << "Threshold = " << threshold << endl;
	
	// get the water depth raster
	LSDRaster WaterDepth((path_name+DEM_name), DEM_extension);
	cout << "Got the water depth raster" << endl;
	
	// remove water depths below the threshold
	LSDRaster WaterDepthFinal = WaterDepth.RemoveBelow(threshold);
	cout << "Removed values below the threshold" << endl;
	
	// convert to binary 
	int Value = 1;
	int ndv = 0;
	LSDIndexRaster BinaryRaster = WaterDepthFinal.ConvertToBinary(Value, ndv);
  string mask_name = "_binary";
  BinaryRaster.write_raster((path_name+DEM_name+mask_name), DEM_extension); 
}
