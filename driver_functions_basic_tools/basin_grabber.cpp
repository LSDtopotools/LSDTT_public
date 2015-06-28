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
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDChannel.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDMostLikelyPartitionsFinder.hpp"
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=7)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Basingrabber tool!                   ||" << endl;
    cout << "|| I am used to collect basins from a DEM.             ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires six inputs: " << endl << endl;
    cout << "* First the path to the parameter files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "   (Either \\ or / depending on your operating system.)" << endl << endl;
    cout << "* Second the prefix of the DEM. The DEM must be a .bil file" << endl;
    cout << "   For example, if the DEM is called Spain.bil the filename is Spain." << endl << endl;
    cout << "* Third, a flag for reading a source file. " << endl;
    cout << "  0 for false, 1 for true using nodeindices, 2 for true using x,y locations." << endl << endl;
    cout << "* Fourth, either the name of the sources file or the" << endl
         << "   drainage area threshold in pixels for channels." << endl << endl;
    cout << "* Fifth, the window size, in metres, of the polynomial window." << endl;
    cout << "   The window will be changed to sqrt 2*DataResolution if the" << endl;
    cout << "   provided value is smaller than this radius." << endl << endl;
    cout << "* Sixth, the order of the basins to be extracted." << endl;
    cout << "=========================================================" << endl;
    cout << "An example call is: " << endl;
    cout << "basin_grabber ~home/basins/Chile/ Chile_test 1 macon_lidem_utm_bl_er_CH_nodeindices_for_Arc 5 3" << endl;
    exit(EXIT_SUCCESS);
  }

  //load input arguments
  string path = argv[1];
  string filename = argv[2];
  int SourceReadFlag = atoi(argv[3]);

  int threshold_pixels = 10;
  string sources_filename;
  bool channel_network_from_file = false;
  int read_csv_flag = 0;    // a flag for reading the csv file

  cout << endl <<"===============================" << endl;
  cout << "Are you reading the sources from file?" << endl;


  if (SourceReadFlag == 0)
  {
    cout << "No, I am using a threshold" << endl;
    threshold_pixels = atoi(argv[4]);
  }
  else
  {
    sources_filename = argv[4];
    cout << "Yes, I am using this csv file: " << sources_filename << endl;
    channel_network_from_file = true;
    if (SourceReadFlag == 1)
    {
      cout << "I am going to use the rows and columns read the csv file" << endl;
      read_csv_flag = SourceReadFlag;
    }
    else if (SourceReadFlag == 2)
    {
      cout << "I am going to use the x and y locations to read the csv file" << endl;
      read_csv_flag = SourceReadFlag;
    }
    else
    {
      cout << "I am goung to use the nodeindices to read the csv file" << endl;
    }
    
    
  }
  cout <<"===============================" << endl;
  float window_radius = atof(argv[5]);
  int BasinOrder = atoi(argv[6]);

  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  //load dem
  LSDRaster DEM((path+filename), "bil");

  //Fill
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  //surface fitting
  vector<int> raster_selection;

  float CriticalSlope = 1.2; //this needs modified to generate proper E*R* data


  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(1); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(window_radius, raster_selection);
  LSDRaster slope = Surfaces[1];
  LSDRaster aspect = Surfaces[2];

  cout << "\nGetting drainage network and basins\n" << endl;

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  // a sources vector
  vector<int> sources;

  if(channel_network_from_file)
  {
    //get stream net from channel heads
    cout << "IMPORTANT: your sources file should be a csv!!" << endl;
    sources = FlowInfo.Ingest_Channel_Heads((path+sources_filename), "csv",read_csv_flag); //swap to csv?
  }
  else
  {
    // get contributing pixels (needed for junction network)
    LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    // get the sources
    sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold_pixels);
  }

  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  LSDIndexRaster StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();

  LSDIndexRaster Basin_Raster;

  //Extract basins based on input stream order
  cout << "\nGetting basins\n";
  if (BasinOrder == 0)
  {
    vector< int > basin_junctions = ChanNetwork.get_Junctions_of_Sources(FlowInfo);
    Basin_Raster = ChanNetwork.extract_basins_from_junctions_rudimentary(basin_junctions, FlowInfo);
  }
  else
  {
    vector< int > basin_junctions = ChanNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
    Basin_Raster = ChanNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);
  }

  cout << "\nExtracting hilltops and hilltop curvature" << endl;

  // extract ridges and then hilltops based on critical slope
  LSDRaster Ridges = ChanNetwork.ExtractRidges(FlowInfo);
  LSDRaster Hilltops = ChanNetwork.ExtractHilltops(Ridges, slope, CriticalSlope);
  
  // get the drainage area (D8)
  LSDRaster DrainageAgea = FlowInfo.write_DrainageArea_to_LSDRaster();
  

  //if the user requests the raster to be written, write the rasters
  int WriteRasters = 1;       // just print for now, deal with logic later
  if (WriteRasters == 1){
    cout << "Writing Rasters\n" << endl;
    FilledDEM.write_raster((path+filename+"_Fill"), "bil");
    Surfaces[1].write_raster((path+filename+"_Slope"),"bil");
    Surfaces[2].write_raster((path+filename+"_Aspect"),"bil");
    Surfaces[3].write_raster((path+filename+"_Curvature"),"bil");
    StreamNetwork.write_raster((path+filename+"_STNET"), "bil");
    Hilltops.write_raster((path+filename+"_Hilltops"), "bil");
    Basin_Raster.write_raster((path+filename+"_Basins"), "bil");
    DrainageAgea.write_raster((path+filename+"_D8DA"), "bil");
  }
}
