//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// reliability_sensitivity
//
// This program takes two arguments, the path name and the driver name
// It calculates the reliability and sensitivity indices of two rasters: the raster of predicted
// values (e.g. floodplain pixels), and the raster of values to be tested against.
// In the case of floodplain initiation this will be the mask predicted by the geomorphic
// algorithm, and the flood hazard maps.
// Flood hazard maps: Zone of no flooding has a raster value of 2, high-risk zones (100 year flood)
// have values of 1 and 3.
//
// The reliability and sensitivity are calculated following Orlandini et al. (2011)
// Reliability (R) = sum(TP) / (sum(TP) + sum(FP))
// Sensitivity (S) = sum(TP) / (sum(TP) + sum (FN))
// where TP = true positive, FP = false positive, FN = false negative
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 04/05/16
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
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string predicted_raster = argv[2];
	string actual_raster = argv[3];

	cout << "The path is: " << path_name << endl;
	cout << "The predicted raster is: " << predicted_raster << endl;
	cout << "The actual raster is: " << actual_raster << endl;

  string DEM_extension = "bil";

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
	vector<float> quality_results = PredictedRaster.AnalysisOfQuality(ActualRaster);

	// print
	for (int i = 0; i < int(quality_results.size()); i++)
	{
		cout << float(quality_results[i]) << endl;
	}

}
