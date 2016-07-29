//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// fill_no_data.cpp
//
// Driver to fill in no data holes in the raster. It is // first trimmed to remove no data around the edges
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Fiona J. Clubb
// University of Edinburgh
// 29/07/16
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
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the DEM filename, and tghe DEM file format."; 
		exit(EXIT_FAILURE);
	}
	
  //get input args
  string path = argv[1];
  string DEM_Name = argv[2];
  string DEM_Format = argv[3];

  //load the DEM 
  LSDRaster DEM((path+DEM_Name), DEM_Format);
	
	// fill in the holes
	int window_size = 10;
	LSDRaster FilledDEM = DEM.alternating_direction_nodata_fill_with_trimmer(window_size);
	
	string Filled_name = "_noholes";
	FilledDEM.write_raster((path+DEM_Name+Filled_name), DEM_Format);
	

}
