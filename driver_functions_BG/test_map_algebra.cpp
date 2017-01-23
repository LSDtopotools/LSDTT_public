//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// change_raster_value_driver.cpp
//
// This program should be able to to change a specific value for a raster
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
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
 

int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
  cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
  exit(EXIT_SUCCESS);
  }

  // make sure there is a slash on the end of the file
  string path_name = argv[1]; 
  path_name = FixPath(path_name);
	

  string f_name = argv[2];
  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;
  string full_name = path_name+f_name;
  
  
  cout << "\nYou are testing a map algebra function." << endl
     <<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
     <<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
     <<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;
  
  string DEM_bil_extension = "bil";
    
  // load the DEM
  LSDRaster topo_test(full_name, DEM_bil_extension);
  
  
  float Minimum_Slope = 0.0001;
  cout << "I'm filling your raster" << endl;
  LSDRaster fill_raster = topo_test.fill(Minimum_Slope);
  
  LSDRaster difference_raster = fill_raster.MapAlgebra_subtract(topo_test);
  
  string fill_ext = "_FILL";
  string diff_ext = "_DIFF";
  string fill_name = full_name+fill_ext;
  string diff_name = full_name+diff_ext;
  
  
  //write the output raster
  fill_raster.write_raster(fill_name, "bil");
  difference_raster.write_raster(diff_name, "bil");
}
