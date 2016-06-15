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
  float low_CrValue;
  float hi_CrValue;
  float low_CsValue;
  float hi_CsValue;
  float low_DValue;
  float hi_DValue;
  float low_phiValue;
  float hi_phiValue;
  float low_RoverTValue;
  float hi_RoverTValue;
  float low_rhosValue;
  float hi_rhosValue;
  float rhow;
  float g;
  int WindowSize;
  string OutPath;
  string temp;

  file_info_in >> temp >> RasterFormat
               >> temp >> DEMName
               >> temp >> low_CrValue >> hi_CrValue
               >> temp >> low_CsValue >> hi_CsValue
	             >> temp >> low_DValue >> hi_DValue
	             >> temp >> low_phiValue >> hi_phiValue
	             >> temp >> low_RoverTValue >> hi_RoverTValue
               >> temp >> low_rhosValue >> hi_rhosValue
               >> temp >> rhow
               >> temp >> g
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
  LSDRaster Slope = Surfaces[1];
  LSDRaster SlopeDeg = Surfaces[1].calculate_slope_angles();

  //Get the D-infinity drainage area
  LSDRaster DrainageArea = FilledDEM.D_inf_units();

  //populate SoilHydroRasters with lower bounds parameter values
  LSDSoilHydroRaster low_RoverT(FilledDEM, low_RoverTValue);
  LSDSoilHydroRaster low_Cs(FilledDEM, low_CsValue);
  LSDSoilHydroRaster low_Cr(FilledDEM, low_CrValue);
  LSDSoilHydroRaster low_phi(FilledDEM, low_phiValue);
  LSDSoilHydroRaster low_D(FilledDEM, low_DValue);
  LSDSoilHydroRaster low_rhos(FilledDEM, low_rhosValue);

  //populate SoilHydroRasters with upper bounds parameter values
  LSDSoilHydroRaster hi_RoverT(FilledDEM, hi_RoverTValue);
  LSDSoilHydroRaster hi_Cs(FilledDEM, hi_CsValue);
  LSDSoilHydroRaster hi_Cr(FilledDEM, hi_CrValue);
  LSDSoilHydroRaster hi_phi(FilledDEM, hi_phiValue);
  LSDSoilHydroRaster hi_D(FilledDEM, hi_DValue);
  LSDSoilHydroRaster hi_rhos(FilledDEM, hi_rhosValue);

  //Calculate dimensionless parameters used in final Fs calculations
  LSDSoilHydroRaster hi_h = hi_D.Calculate_h(SlopeDeg);
  LSDSoilHydroRaster hi_w = hi_RoverT.Calculate_w(SlopeDeg, DrainageArea);
  LSDSoilHydroRaster hi_C = hi_Cr.Calculate_C(hi_Cs, hi_h, low_rhos, g);
  LSDSoilHydroRaster hi_r = hi_rhos.Calculate_r(rhow);

  LSDSoilHydroRaster low_h = low_D.Calculate_h(SlopeDeg);
  LSDSoilHydroRaster low_w = low_RoverT.Calculate_w(SlopeDeg, DrainageArea);
  LSDSoilHydroRaster low_C = low_Cr.Calculate_C(low_Cs, low_h, hi_rhos, g);
  LSDSoilHydroRaster low_r = low_rhos.Calculate_r(rhow);

  //Calculate factor of safety
  LSDSoilHydroRaster Fs_min = hi_C.Calculate_sinmap_Fs(SlopeDeg, low_w, low_r, hi_phi);
  LSDSoilHydroRaster Fs_max = low_C.Calculate_sinmap_Fs(SlopeDeg, hi_w, hi_r, low_phi);

  LSDSoilHydroRaster SI = low_r.Calculate_sinmap_SI(Slope, DrainageArea, low_C, hi_C, low_phi, hi_phi, low_RoverT, hi_RoverT, low_r, hi_r, Fs_min,  Fs_max);

  //Write factor of safety to file
  Fs_min.write_raster((OutPath+"Fs_min"), RasterFormat);
  Fs_max.write_raster((OutPath+"Fs_max"), RasterFormat);

  //write Stability index to a file
  SI.write_raster((OutPath+"SI"), RasterFormat);

}
