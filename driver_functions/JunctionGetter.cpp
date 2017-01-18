//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// JunctionGetter.cpp
//
// Simple driver to get the junctions of all non-beheaded basins in a DEM of a
// given stream order. Pass in the path to the data, the dem name without the extension,
// and the desired stream order.
//
// Writes junctions to a file in the same directory as the input DEM called 
// Junctions.txt which can be passed straight into the BasinMetrics.cpp driver.
//
// Must use a *.flt raster.
//
// Version 1.0
//
// 24/10/13
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Stuart W.D. Grieve
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

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path, the name of the DEM file and the stream order of basin junctions to be extracted." << endl;
		exit(EXIT_SUCCESS);
	}

  //get input path and file names
  string DataPath = argv[1];
  string param_file_name = argv[2];
  string DEM_name = DataPath + param_file_name;
  string OutputJunctionFile = DataPath + "Junctions.txt";
  
  //get basin order as an integer
  int BasinOrder = atoi(argv[3]);
    
  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  string DEM_file_extension = "flt"; 
  
  //Load and fill the DEM
	LSDRaster DEM(DEM_name, DEM_file_extension);
  float MinSlope = 0.0001; 
  LSDRaster FilledDEM = DEM.fill(MinSlope);
 	
	// get a flow info object
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  //------------------------------------------------------------------------------------------
  //This gets the channel network and will be updated with Fiona's channel head stuff 
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	int threshold = 800;
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	//end of channel head stuff
  //------------------------------------------------------------------------------------------
  
  vector<int> Junctions = ChanNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
    
  ofstream output(OutputJunctionFile.c_str()); //is this right??
  
  for (int q = 0; q < int(Junctions.size()); ++q){
    output << Junctions[q] << endl;
  }
   
  output.close(); 
  cout << "\n\nList of junctions saved in: " << OutputJunctionFile << endl; 
}
