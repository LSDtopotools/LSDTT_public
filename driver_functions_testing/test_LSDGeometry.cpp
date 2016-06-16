//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Test_node_locations.cpp
//
// This is a function for testing the location and printing of points from
// and to rasters
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
#include <ctime>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDGeometry.hpp"
#include "../LSDFileTools.hpp"

int main (int nNumberofArgs,char *argv[])
{

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=4)
  {
    cout << "==============================================================" << endl;
    cout << "|| Welcome to test_LSDGeometry tool!                        ||" << endl;
    cout << "|| I am used to test the writing and reading of geometry data. ||" << endl;
    cout << "==============================================================" << endl;
    cout << "This program requires two inputs: " << endl << endl;
    cout << "* First the path to the DEM." << endl;
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
  
  LSDRasterInfo RI((path+filename), extension);


  // get the min and max x and y locations
  vector<float> XYMinMax = DEM.get_XY_MinMax();
  cout << "Got the Min and Max values" << endl;
  cout << "Xmin: " << XYMinMax[0] << endl;
  cout << "YMin: " << XYMinMax[1] << endl;
  cout << "Xmax: " << XYMinMax[2] << endl;
  cout << "YMax: " << XYMinMax[3] << endl;
  
  // these get the ranges of X and Y for random location generation
  float XRange = XYMinMax[2]-XYMinMax[0];
  float YRange = XYMinMax[3]-XYMinMax[1];
  long seed = time(NULL);
  
  // these are the vectors for holding the random points
  vector<float> x_locs;
  vector<float> y_locs;
  
  float xl,yl;
  
  // get some random points to generate a geometry object
  for (int i = 0; i< 2; i++)
  {
    // get random x,y locations
    //cout << "Ran:" << ran3(&seed) << endl;
    
    xl = XYMinMax[0]+XRange*ran3(&seed);
    yl = XYMinMax[1]+YRange*ran3(&seed);
    
    // add to vectors
    x_locs.push_back(xl);
    y_locs.push_back(yl);
  }
  
  // get the UTM info from the raster
  int UTM_zone;
  bool is_North;
  DEM.get_UTM_information(UTM_zone, is_North);
  
  
  // now make a geometry object
  LSDGeometry SomePoints(x_locs,y_locs, UTM_zone, is_North);
  
  // get the latlong
  SomePoints.convert_points_to_LatLong();
  
  bool isNorth = SomePoints.get_isNorth();
  int UTMZone = SomePoints.get_UTMZone();
  cout << "Is north? " << isNorth << " and UTMZone: " << UTMZone << endl;

  // now print to csv
  string fname_prefix = "Points_Test";
  SomePoints.print_points_to_csv(path,fname_prefix);
  
  LSDPolyline ALine(x_locs,y_locs, UTM_zone, is_North);
  ALine.make_simple_polyline();
  
  vector<int> RowOfNodes;
  vector<int> ColOfNodes;
  ALine.find_row_and_col_of_points(RI, RowOfNodes, ColOfNodes);
  
  int n_nodes = int(RowOfNodes.size());
  for(int i = 0; i<n_nodes; i++)
  {
    cout << "i:" << i << " row: " << RowOfNodes[i]<< " col: " << ColOfNodes[i] << endl;
  } 

  vector<int> affected_rows;
  vector<int> affected_cols;
  int start_node = 0;
  int end_node = 1;
  ALine.get_affected_pixels_in_line_segment_brute_force(RI,affected_rows, affected_cols, 
                                    start_node, end_node); 

}
