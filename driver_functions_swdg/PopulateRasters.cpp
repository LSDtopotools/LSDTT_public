#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../TNT/tnt.h"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=4)
  {
    cout << "FATAL ERROR: wrong number of inputs. The program needs the path name, the DEM name and the DEM format." << endl;
    exit(EXIT_FAILURE);
  }

  string Path = argv[1];
  string DEMName = argv[2];
  string RasterFormat = argv[3];

  LSDRaster DEM((Path + DEMName), RasterFormat);

  LSDRaster RoverT = DEM.PoupulateRasterGaussian(0.0005, 0.0003);
  LSDRaster Cs = DEM.PoupulateRasterGaussian(0, 4000);
  LSDRaster Cr = DEM.PoupulateRasterGaussian(3, 12500);
  LSDRaster rhos = DEM.PoupulateRasterGaussian(2000, 2325);
  LSDRaster phi = DEM.PoupulateRasterGaussian(30, 37.5);
  LSDRaster D = DEM.PoupulateRasterGaussian(0.1, 1.05);

  RoverT.write_raster((Path+"RoverT"), RasterFormat);
  Cs.write_raster((Path+"Cs"), RasterFormat);
  Cr.write_raster((Path+"Cr"), RasterFormat);
  rhos.write_raster((Path+"rhos"), RasterFormat);
  phi.write_raster((Path+"phi"), RasterFormat);
  D.write_raster((Path+"D"), RasterFormat);

}
