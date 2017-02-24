//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// get_elevation_profiles.cpp
// Takes in a baseline shapefile and gets the mean, max, and min elevation at each point
// along the shapefile.
// Writes the values to a csv file for visualisation.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../../LSDRaster.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDJunctionNetwork.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();

	if (nNumberofArgs != 3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }
  string path_name = argv[1];
  string f_name = argv[2];
  cout << "The path name is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if (file_info_in.fail())
	{
		cout << "\nFATAL ERROR: the header file\"" << full_name
				 << "\" doesn't exist" << endl;
	}

	string DEM_ID, Baseline_file, temp;
	int NormaliseToBaseline;
	float HalfWidth;
	string DEM_extension = "bil";
	string csv_extension = "csv";

	// read in the parameters. need to have an input string before each parameter in the input file. I will update this with the parameter parser at some point when I have time...
	file_info_in >> temp >> DEM_ID
						   >> temp >> Baseline_file
							 >> temp >> HalfWidth
							 >> temp >> NormaliseToBaseline;

	file_info_in.close();

  //string Swath_ext = "_swath_trans";
  //string Long_Swath_ext = "_swath_long";
  //string BV_ext = "_baseline_values";
  cout << "starting the test run... here we go!" << endl;

	// load the DEM
	cout << "Loading the DEM..." << endl;
	LSDRaster RasterTemplate((path_name+DEM_ID), DEM_extension);

  cout << "\t loading baseline points" << endl;
  PointData BaselinePoints = LoadShapefile(path_name+Baseline_file.c_str());

	// get the swath
  cout << "\t creating swath template" << endl;
  LSDSwath TestSwath(BaselinePoints, RasterTemplate, HalfWidth);

	cout << "\n\t Getting raster from swath" << endl;
	LSDRaster SwathRaster = TestSwath.get_raster_from_swath_profile(RasterTemplate, NormaliseToBaseline);
	string swath_ext = "_swath_raster";
	SwathRaster.write_raster((path_name+DEM_ID+swath_ext), DEM_extension);

  // get the raster values along the swath
	vector <vector <float> > ElevationValues = TestSwath.get_RasterValues_along_swath(RasterTemplate, NormaliseToBaseline);

	// push back results to file for plotting
	ofstream output_file;
	string output_fname = "_swath_elevations.csv";
	output_file.open((path_name+DEM_ID+output_fname).c_str());
	output_file << "Distance,Mean,Min,Max" << endl;
	for (int i = 0; i < int(ElevationValues[0].size()); ++i)
	{
		output_file << ElevationValues[0][i] << "," << ElevationValues[1][i] << "," << ElevationValues[2][i] << "," << ElevationValues[3][i] << endl;
	}
	output_file.close();

	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
