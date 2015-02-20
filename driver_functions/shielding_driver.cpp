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

  //run hillshade
  for (int Azimuth = 25.; Azimuth < 35; Azimuth += 45)
  {
    stringstream ss1;
		ss1 << Azimuth;
		string AzimuthString = ss1.str();
    for (int Zenith = 35; Zenith <= 35; Zenith+=10)
    {
    	LSDIndexRaster Shadows = DEM.CastShadows((float)Azimuth,(float)Zenith);
    	stringstream ss2;
		  ss2 << Zenith;
		  string ZenithString = ss2.str();
    	//string AzString = string(Zenith);
    	Shadows.write_raster((filename+"_Shadows_"+AzimuthString+"_"+ZenithString),"flt");
    }
  }
}
