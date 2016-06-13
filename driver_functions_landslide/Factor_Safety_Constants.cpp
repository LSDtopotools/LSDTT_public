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
#include "../LSDSoilHydroRaster.hpp"
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

  //Load parameter values

  string RasterFormat;
  string DEMName;
  float CrValue;
  float CsValue;
  float DValue;
  float rhosValue;
  float phiValue;
  float RoverTValue;
  int WindowSize;
  float g;
  float rhow;
  string OutPath;
  string temp;

  file_info_in >> temp >> RasterFormat
               >> temp >> DEMName
               >> temp >> CrValue
               >> temp >> CsValue
	             >> temp >> DValue
	             >> temp >> rhosValue
	             >> temp >> phiValue
	             >> temp >> RoverTValue
               >> temp >> g
               >> temp >> rhow
               >> temp >> WindowSize
	             >> temp >> OutPath;
  file_info_in.close();

  // Load the DEM and fill it
  LSDRaster DEM((Path + DEMName), RasterFormat);
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  //Get the slope using the surface fitting routines
  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(0); //aspect
  raster_selection.push_back(0); //curvature
  raster_selection.push_back(0); //plan curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(WindowSize, raster_selection);

  LSDRaster SlopeDeg = Surfaces[1].calculate_slope_angles();

  //Get the D-infinity drainage area
  LSDRaster DrainageArea = FilledDEM.D_inf_units();

  //populate SoilHydroRasters with parameter values
  LSDSoilHydroRaster RoverT(FilledDEM, RoverTValue);
  LSDSoilHydroRaster Cs(FilledDEM, CsValue);
  LSDSoilHydroRaster Cr(FilledDEM, CrValue);
  LSDSoilHydroRaster rhos(FilledDEM, rhosValue);
  LSDSoilHydroRaster phi(FilledDEM, phiValue);
  LSDSoilHydroRaster D(FilledDEM, DValue);

  //Calculate dimensionless parameters used in final Fs calculation
  LSDSoilHydroRaster h = D.Calculate_h(SlopeDeg);
  LSDSoilHydroRaster w = RoverT.Calculate_w(SlopeDeg, DrainageArea);
  LSDSoilHydroRaster r = rhos.Calculate_r(rhow);
  LSDSoilHydroRaster C = Cr.Calculate_C(Cs, h, rhos, g);

  //Calculate factor of safety
  LSDSoilHydroRaster Fs = C.Calculate_sinmap_Fs(SlopeDeg, w, r, phi);

  //Write factor of safety to file
  Fs.write_raster((OutPath+"Fs"), RasterFormat);

}
