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

int main (int nNumberofArgs,char *argv[])
{

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=4)
  {
    cout << "==============================================================" << endl;
    cout << "|| Welcome to theTest_node_locations tool!                  ||" << endl;
    cout << "|| I am used to test the writing and reading of point data. ||" << endl;
    cout << "|| I make some random points within a DEM and write them    ||" << endl;
    cout << "|| and read them, making sure everything is in the correct  ||" << endl;
    cout << "|| place. You can check in a GIS.                           ||" << endl;
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

  // remove seas
  DEM.remove_seas();
  
  // now get a flow info object
  LSDFlowInfo FlowInfo(DEM);
  cout << "Got flow info" << endl;
  
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
  
  // get some random points to print to a csv file 
  for (int i = 0; i< 50; i++)
  {
    // get random x,y locations
    //cout << "Ran:" << ran3(&seed) << endl;
    
    xl = XYMinMax[0]+XRange*ran3(&seed);
    yl = XYMinMax[1]+YRange*ran3(&seed);
    
    // add to vectors
    x_locs.push_back(xl);
    y_locs.push_back(yl);
  }
  
  // now write the node index raster
  cout << "Writing node info to LSDIndexRaster" << endl;
  LSDIndexRaster NodeIndex = FlowInfo.write_NodeIndex_to_LSDIndexRaster();
  cout << "Wrote node info to LSDIndexRaster" << endl;
  
  // now get the node indices of the points and print to csv
  ofstream random_points;
  string cvfname = path+"random_points.csv";
  random_points.open(cvfname.c_str());
  random_points << "x,y,id" << endl;
  random_points.precision(8);
  
  
  int n_nodes = x_locs.size();
  int this_ni;
  for (int i = 0; i<n_nodes; i++)
  {
    this_ni = FlowInfo.get_node_index_of_coordinate_point(x_locs[i], y_locs[i]);
    
    // print to csv
    random_points << x_locs[i] << "," << y_locs[i] << "," << this_ni << endl;
    
  }

  // now get some sources
  vector<int> sources;
  int threshold = 250;
  LSDIndexRaster FlowPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  sources = FlowInfo.get_sources_index_threshold(FlowPixels, threshold);
  
  
  
  // now write the channel network to raster
  LSDJunctionNetwork JN(sources,FlowInfo);
  LSDIndexRaster SO = JN.StreamOrderArray_to_LSDIndexRaster();
  
  // now print the sources to a csv file
  string outfilename = path+"Raster_sources";
  FlowInfo.print_vector_of_nodeindices_to_csv_file(sources, outfilename);
  outfilename = outfilename+"_nodeindices_for_Arc";
  
  // now read the sources and get the basins
  int read_csv_flag = 2;      // this is for reading x,y data
  cout << "Reading sources" << endl;
  vector<int> new_sources = FlowInfo.Ingest_Channel_Heads(outfilename,"csv",read_csv_flag);
  cout << "Got sources" << endl;
  
  // check if the sources have the same nodeindex
  int n_sources = new_sources.size();
  cout << endl << endl << "Here are the new and old sources!" << endl;
  for(int i = 0; i<n_sources; i++)
  {
    cout << sources[i] << "\t" << new_sources[i] << endl;
  }
  
  // use these sources to get the basins
  cout << "Getting network" << endl;
  LSDJunctionNetwork NewChanNetwork(new_sources, FlowInfo);
  cout << "Got network" << endl;
  
  cout << "Getting channel head junctions" << endl;
  vector< int > basin_junctions = NewChanNetwork.get_Junctions_of_Sources(FlowInfo);
  cout << "Got channel head junctions" << endl;
  
  cout << "Getting basin raster" << endl;
  LSDIndexRaster Basin_Raster = NewChanNetwork.extract_basins_from_junctions_rudimentary(basin_junctions, FlowInfo);
  cout << "Got basin raster" << endl;
  
  // Now print the SO as csv
  string SO_prefix =path+"Stream_order_flattened";
  SO.FlattenToCSV(SO_prefix);
  
  // now write the DEMS
  string bil_ext = "bil";
  NodeIndex.write_raster(path+filename+"_NI",bil_ext);
  SO.write_raster(path+filename+"_SO",bil_ext);
  Basin_Raster.write_raster(path+filename+"_Basins",bil_ext);


}
