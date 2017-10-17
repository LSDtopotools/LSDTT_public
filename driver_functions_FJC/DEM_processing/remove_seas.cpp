//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// remove_seas.cpp
// make with remove_seas.make
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

  string fname = path_name+DEM_ID+"."+DEM_extension;
  ifstream file_info_in;
  file_info_in.open(fname.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << fname
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

	LSDRaster DEM(path_name+DEM_ID, DEM_extension);

  // remove the seas
	cout << "\t Removing seas" << endl;
  DEM.remove_seas();

  //removing some weird values
  // float threshold = 0.5;
  // bool below_threshold = true;
  // DEM.mask_to_nodata_using_threshold(threshold, below_threshold);

  // trim the raster_array
  cout << "\t Trimming the raster" << endl;
  //DEM.RasterTrimmer();
	DEM.write_raster((path_name+DEM_ID), DEM_extension);
}
