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
  if (nNumberofArgs!=3)
  {
      cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the DEM name without extension." << endl;
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
	int window_radius = 250;
	//LSDRaster FilledDEM = DEM.nodata_fill_irregular_raster(window_radius);

	LSDRaster FilledDEM = DEM.alternating_direction_nodata_fill_irregular_raster(window_radius);
	
	string fill_ext = "_NoHoles";
	FilledDEM.write_raster((path_name+DEM_ID+fill_ext), DEM_extension);
	
	// write a hillshade of the filled DEM for checking
	string HS_ext = "_HS";
	LSDRaster HS = FilledDEM.hillshade(45, 315, 1);
	HS.write_raster((path_name+DEM_ID+HS_ext),DEM_extension);
	
}
