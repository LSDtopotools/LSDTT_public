//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LH_Driver.cpp
//
// Driver created to automate extraction of hillslope profiles
//
// Driver expects an unfilled DEM in the given directory in flt format with the name format <prefix>.bil,
// a channel heads file in the same directory with the format <prefix>_CH.flt and a floodplain raster
// <prefix>_FloodPlain.flt if floodplains are to be excluded from the analysis.
//
// Run driver with the following arguments:
//
// path to the input files with a trailing slash
// filename prefix without an underscore
// window radius value in spatial units for surface fitting
// basin order the strahler order of basins to be extracted
// switch to exclude floodplains (1) or not exclude them (0) 
// switch to write rasters 0 == do not write rasters and 1 == write rasters
//
// A usage example is:
//nice ./LH_Driver.out /home/s0675405/DataStore/Final_Paper_Data/NC/ NC 7 2 1 0
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Martin D. Hurst
// University of Glasgow
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
	
	//Test for correct input arguments
	if (nNumberofArgs!=2)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash)." << endl;
		exit(EXIT_SUCCESS);
	}
	
	//load input arguments
	string path = argv[1];
	
  string ChannelSegments_file = "bolinas_Segments";
  string RasterExt = "bil";
  LSDIndexRaster StreamNetwork(path+ChannelSegments_file, RasterExt);
  cout << " done.";

  int NPixels = 2;  
  cout << "\n\tPadding the stream network by " << NPixels << " pixels" << endl;
  StreamNetwork.PadRaster(NPixels);
  string padded_raster_name = path+"padded";
  StreamNetwork.write_raster(padded_raster_name,RasterExt);
  cout << " done.";
}
