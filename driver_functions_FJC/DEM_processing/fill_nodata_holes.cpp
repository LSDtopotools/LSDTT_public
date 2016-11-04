//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// fill_nodata_holes.cpp
// make with fill_nodata_holes.make
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../../LSDRaster.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
  if (nNumberofArgs!=2)
  {
      cout << "FATAL ERROR: wrong number inputs. The program needs the path name." << endl;
      exit(EXIT_SUCCESS);
  }

	string path_name = argv[1];
  string DEM_ID = argv[2];
	string DEM_extension = "bil";

	//load the DEM
  cout << "\t Loading the DEM" << endl;
	LSDRaster DEM(path_name+DEM_ID, DEM_extension);

	cout << "\t Filling no data holes" << endl;
	// no data filling
	int window_radius = 50;
	LSDRaster FilledDEM = DEM.alternating_direction_nodata_fill_with_trimmer(window_radius);
	
	string fill_ext = "_NoHoles";
	FilledDEM.write_raster((path_name+DEM_ID+fill_ext), DEM_extension);
}
