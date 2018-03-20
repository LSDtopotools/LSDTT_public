//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// change_junctions_to_xy
// make with change_junctions_to_xy.make
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
  if (nNumberofArgs!=4)
  {
      cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the DEM name without extension, and the junction file name." << endl;
      exit(EXIT_SUCCESS);
  }

	string path_name = argv[1];
  string DEM_ID = argv[2];
  string junction_fname = argv[3];
	string DEM_extension = "bil";

	//load the DEM
  cout << "\t Loading the DEM" << endl;
	LSDRaster DEM(path_name+DEM_ID, DEM_extension);

	cout << "Getting the junction list" << endl;

}
