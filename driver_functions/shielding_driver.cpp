/*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	shielding_driver.cpp

	Driver file for testing topographic shielding functions for CRN production
	Testing the implementation of drop shadows/cast shadows
	
	Martin D. Hurst
	British Geological Survey
	University of Edinburgh

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
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
  
//  //Test for correct input arguments
//	if (nNumberofArgs!=7)
//	{
//		cout << "FATAL ERROR: wrong number inputs. The program needs the path (with trailing slash), the filename prefix, window radius, basin order, critical slope and a switch to write rasters if desired." << endl;
//		exit(EXIT_SUCCESS);
//	}
  
  //load input arguments
  string filename = "bh";
  
  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  //load dem
  LSDRaster DEM((filename+"_dem"), "flt");  

  //launch toposhielding
  LSDRaster TopoShielding = DEM.TopographicShielding(5,5);
  TopoShielding.write_raster(filename+"_TopoShield","flt");
}
