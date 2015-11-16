//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// PolyFitWindowSize.cpp
//
// Driver created to automate estimation of the corect window size to use for the
// polyfit routine.
//
// Driver expects an unfilled DEM in the given directory of the given name and format.
//
// Run driver with the following arguments:
//
// path to the DEM file with a trailing slash
// DEM filename without the file extension or dot
// file extension without the dot

// A usage example is:
// nice ./PolyFitWindowSize.out /home/s0675405/DataStore/Final_Paper_Data/NC/ NC_DEM flt
//
// Output data is stored in the input directory in a file called <DEM_Name>_Window_Size_Data.txt
// and can be plotted using WindowSize.py  found at https://github.com/sgrieve/GeneralAnalysis/
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// June 2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDShapeTools.hpp"

int main(int nNumberofArgs, char *argv[])
{
  //Test for correct input arguments
	if (nNumberofArgs!=6)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the filename prefix, the DEM file format, the window size in spatial units and the curvature type switch.";
		exit(EXIT_FAILURE);
	}

  //get input args
  string path = argv[1];
  string Prefix = argv[2];
  string DEM_Format = argv[3];
  int WindowSize = atoi(argv[4]);
	int Curvature_Type = atoi(argv[5]);

  //surface fitting
  vector<int> raster_selection;
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

	if (Curvature_Type == 3){
    raster_selection.push_back(1); //curvature
  }
	else if (Curvature_Type == 4){
		raster_selection.push_back(0);
    raster_selection.push_back(1); //planform curvature
  }
	else if (Curvature_Type == 5){
		raster_selection.push_back(0);
		raster_selection.push_back(0);
		raster_selection.push_back(1); //profile curvature
	}
	else if (Curvature_Type == 6){
		raster_selection.push_back(0);
		raster_selection.push_back(0);
		raster_selection.push_back(0);
		raster_selection.push_back(1); //tangential curvature
	}
  else{
		cout << "FATAL ERROR: The Curvature type selected was invalid. Use a value of 3, 4, 5 or 6 to get curvature, plan curvature, profile curvature or tan curvature.";
		exit(EXIT_FAILURE);
	}

  //set up a writer to write the output data
  ofstream WriteData;

  //create an output filename based on the dem name
  stringstream ss;

  ss << path << Prefix << "_CurvatureResData_" << Curvature_Type << ".txt";
  WriteData.open(ss.str().c_str());

  //write headers
  WriteData << "resoulution 2pc 25pc median mean 75pc 98pc minimum maximum" << endl;

  //array of resolutions to load
  int Resolutions[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 60};

  // vectors to hold the stats about the fitted surface
  vector<float> Curv_vec;

  for (int a = 0; a < 26; ++a){

	  cout << "Processing DEM " << a+1 << " of " << "26" << endl;

		//load the DEM
    //build the string of the filename to load
		stringstream ss2;
		ss2 << path << Prefix << "_" << Resolutions[a] << "_DEM";
		LSDRaster DEM(ss2.str(), DEM_Format);

	  int CurrentWindowSize = WindowSize * Resolutions[a];
	  vector<LSDRaster> Surfaces = DEM.calculate_polyfit_surface_metrics(CurrentWindowSize, raster_selection);
	  LSDRaster curvature = Surfaces[Curvature_Type];

	  //go through the landscape and get every curvature value into a 1D vector
    Curv_vec = Flatten_Without_Nodata(curvature.get_RasterData(), curvature.get_NoDataValue());

    vector<float> Boxplot = BoxPlot(Curv_vec);

	  //write the values to the output file
	  WriteData << Resolutions[a] << " " << Boxplot[0] << " " << Boxplot[1] << " " << Boxplot[2] << " " << Boxplot[3] << " " << Boxplot[4] << " " << Boxplot[5] << " " << Boxplot[6] << " " << Boxplot[7];

		WriteData << endl;

  }
  WriteData.close();

}
