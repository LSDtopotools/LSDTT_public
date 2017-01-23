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
  if (nNumberofArgs!=3)
  {
    cout << "FATAL ERROR: wrong number of inputs. The program needs the path name and the driver file name" << endl;
    exit(EXIT_FAILURE);
  }

  string Path = argv[1];
  string ParamFile = argv[2];

  string FullName = Path + ParamFile;

  ifstream file_info_in;
  file_info_in.open(FullName.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the parameter file:\n\n " << FullName
         << "\n\ndoesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  string RasterFormat;
  string DEMName;
  string CrName;
  string CsName;
  string DName;
  string rhosName;
  string phiName;
  string RoverTName;
  string OutPath;
  string temp;

  file_info_in >> temp >> RasterFormat
               >> temp >> DEMName
               >> temp >> CrName
               >> temp >> CsName
	             >> temp >> DName
	             >> temp >> rhosName
	             >> temp >> phiName
	             >> temp >> RoverTName
	             >> temp >> OutPath;
  file_info_in.close();

  //constants: Gravity and water density
  float g = 9.81;
  float rhow = 1000;

  LSDRaster DEM((Path + DEMName), RasterFormat);
  float MinSlope = 0.0001;

  LSDRaster FilledDEM = DEM.fill(MinSlope);

  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(0); //aspect
  raster_selection.push_back(0); //curvature
  raster_selection.push_back(0); //plan curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  int WindowSize = 3;

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(WindowSize, raster_selection);
  LSDRaster Slope = Surfaces[1];

  LSDRaster SlopeDeg = Slope.calculate_slope_angles();

  cout << "got the slope" << endl;
  LSDRaster DrainageArea = FilledDEM.D_inf_units();

cout << "area" << endl;
  LSDRaster Cr((Path + CrName), RasterFormat);
  LSDRaster Cs((Path + CsName), RasterFormat);
  LSDRaster D((Path + DName), RasterFormat);

cout << "3 rasters" << endl;

  LSDRaster rhos((Path + rhosName), RasterFormat);
  LSDRaster phi((Path + phiName), RasterFormat);
  LSDRaster RoverT((Path + RoverTName), RasterFormat);

cout << "3 rasters" << endl;

  LSDRaster h = D.Calculate_h(SlopeDeg);
  cout << "h" << endl;
  LSDRaster w = DrainageArea.Calculate_w(RoverT, SlopeDeg);
  cout << "RoverT" << endl;
  LSDRaster r = rhos.Calculate_r(rhow);
  cout << "r" << endl;
  LSDRaster C = Cr.Calculate_C(Cs, h, rhos, g);

cout << "C" << endl;
  LSDRaster Fs = SlopeDeg.Calculate_sinmap_Fs(C, w, r, phi);

cout << "Fs" << endl;

  Fs.write_raster((OutPath+"Fs"), RasterFormat);

}
