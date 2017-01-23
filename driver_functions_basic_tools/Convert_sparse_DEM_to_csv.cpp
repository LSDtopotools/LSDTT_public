//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// basin_grabber.cpp
//
// This function grabs basins and prints a basin file
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
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=4)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Convert_sparse_DEM_to_csv tool!      ||" << endl;
    cout << "|| I am used to change some raster data into a csv.    ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl << endl;
    cout << "* First the path to the parameter files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "   (Either \\ or / depending on your operating system.)" << endl << endl;
    cout << "* Second the prefix of the DEM." << endl;
    cout << "   For example, if the DEM is called Spain.bil the filename is Spain." << endl << endl;
    cout << "* Third, the extension of the DEM. " << endl;
    cout << "  Options are: bil, flt and asc." << endl << endl;
    cout << "=========================================================" << endl;
    cout << "An example call is: " << endl;
    cout << "Convert_sparse_to_csv ~home/basins/Chile/ Chile_test bil" << endl;
    exit(EXIT_SUCCESS);
  }

  //load input arguments
  string path = argv[1];
  string filename = argv[2];
  string extension = argv[3];

  //load dem
  LSDRaster DEM((path+filename), extension);

  // remove seas
  DEM.remove_seas();

  // now write to csv
  DEM.FlattenToCSV(path+filename);


}
