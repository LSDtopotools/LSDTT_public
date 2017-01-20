//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// compare_terrace_areas
//
// Compare the percentage of terrace pixels between the mapped and predicted terraces.
// Both must be in bil format - can translate shapefiles to raster using gdal_rasterize
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 18/01/17
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

	string predicted_raster;
	string actual_raster;
  string DEM_extension = "bil";
	string temp;
	string input_path;	
  
	// read in the parameters
	file_info_in >> temp >> predicted_raster
							 >> temp >> actual_raster;
	
	// load the predicted raster
	LSDIndexRaster PredictedRaster((path_name+predicted_raster), DEM_extension);
	cout << "Got the predicted raster" << endl;
	// load the actual raster
	LSDIndexRaster ActualRaster((path_name+actual_raster), DEM_extension);
	cout << "Got the actual raster" << endl;
	
	//Get the number of pixels in the DEM
	int NRows = PredictedRaster.get_NRows();
	int NCols = PredictedRaster.get_NCols();
	int NPixels = NRows*NCols;
	cout << "N pixels: " << NPixels << endl;
	
	// calculate the TPs, FPs, TNs, and FNs
	float PercentageDiff = PredictedRaster.GetAreaDifference(ActualRaster);
	cout << "Area: " << PercentageDiff << endl;
	
}
