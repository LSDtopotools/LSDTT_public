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
#include "../LSDShapeTools.hpp"
#include "../LSDGeometry.hpp"
#include "../LSDRasterInfo.hpp"
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
  string low_CrInput;
  string hi_CrInput;
  string low_CsInput;
  string hi_CsInput;
  string low_DInput;
  string hi_DInput;
  string low_phiInput;
  string hi_phiInput;
  string low_RoverTInput;
  string hi_RoverTInput;
  string low_rhosInput;
  string hi_rhosInput;
  float rhow;
  float g;
  int WindowSize;
  string OutPath;
  string temp;
  int is_float;

  file_info_in >> temp >> RasterFormat
               >> temp >> DEMName
               >> temp >> low_CrInput >> hi_CrInput
               >> temp >> low_CsInput >> hi_CsInput
	             >> temp >> low_DInput >> hi_DInput
	             >> temp >> low_phiInput >> hi_phiInput
	             >> temp >> low_RoverTInput >> hi_RoverTInput
               >> temp >> low_rhosInput >> hi_rhosInput
               >> temp >> rhow
               >> temp >> g
               >> temp >> WindowSize
	             >> temp >> OutPath;
  file_info_in.close();

  // Load the DEM and fill it
  LSDRaster DEM((Path+DEMName), RasterFormat);
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  LSDRasterInfo DEMInfo(DEM);

  vector<string> BoundaryConditions(4, "No Flux");
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  vector<PointData> TestPoly = LoadPolyline("/home/sgrieve/LSData/fakeroad_segments.shp");

  vector< vector<int> > UpslopeNodes;

  for (int q = 0; q < int(TestPoly.size()); ++q){

    LSDPolyline test(TestPoly[q].X, TestPoly[q].Y, 17);

    vector<int> nodes = test.get_flowinfo_nodes_of_line(DEMInfo, FlowInfo);
    vector<int> tmp1;
    vector<int> tmp2;

    for (int w = 0; w < int(nodes.size()); ++w){
      tmp2 = FlowInfo.get_upslope_nodes(nodes[w]);
      tmp1.insert(tmp1.end(), tmp2.begin(), tmp2.end());

    }

    UpslopeNodes.push_back(tmp1); //may want to run each tmp1 through the unique method to remove dupes?
    tmp1.clear();
    tmp2.clear();

  }

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

  // Crunch through all the input parameters and load them into SoilHydroRasters

  //root_cohesion
  is_float = isFloat(low_CrInput);
  //Constructs a placeholder LSDRaster in cases where we have no need of one as the parameter has been passed as a float
  LSDRaster low_CrRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + low_CrInput), RasterFormat);
  LSDSoilHydroRaster low_Cr = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(low_CrInput.c_str())) : LSDSoilHydroRaster(FilledDEM, low_CrRaster, 0);
  is_float = isFloat(hi_CrInput);
  LSDRaster hi_CrRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + hi_CrInput), RasterFormat);
  LSDSoilHydroRaster hi_Cr = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(hi_CrInput.c_str())) : LSDSoilHydroRaster(FilledDEM, hi_CrRaster, 0);

  //soil_cohesion
  is_float = isFloat(low_CsInput);
  LSDRaster low_CsRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + low_CsInput), RasterFormat);
  LSDSoilHydroRaster low_Cs = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(low_CsInput.c_str())) : LSDSoilHydroRaster(FilledDEM, low_CsRaster, 0);
  is_float = isFloat(hi_CsInput);
  LSDRaster hi_CsRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + hi_CsInput), RasterFormat);
  LSDSoilHydroRaster hi_Cs = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(hi_CsInput.c_str())) : LSDSoilHydroRaster(FilledDEM, hi_CsRaster, 0);

  //soil_depth
  is_float = isFloat(low_DInput);
  LSDRaster low_DRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + low_DInput), RasterFormat);
  LSDSoilHydroRaster low_D = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(low_DInput.c_str())) : LSDSoilHydroRaster(FilledDEM, low_DRaster, 0);
  is_float = isFloat(hi_DInput);
  LSDRaster hi_DRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + hi_DInput), RasterFormat);
  LSDSoilHydroRaster hi_D = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(hi_DInput.c_str())) : LSDSoilHydroRaster(FilledDEM, hi_DRaster, 0);

  //soil_friction_angle
  is_float = isFloat(low_phiInput);
  LSDRaster low_phiRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + low_phiInput), RasterFormat);
  LSDSoilHydroRaster low_phi = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(low_phiInput.c_str())) : LSDSoilHydroRaster(FilledDEM, low_phiRaster, 0);
  is_float = isFloat(hi_phiInput);
  LSDRaster hi_phiRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + hi_phiInput), RasterFormat);
  LSDSoilHydroRaster hi_phi = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(hi_phiInput.c_str())) : LSDSoilHydroRaster(FilledDEM, hi_phiRaster, 0);

  //recharge_over_transmissivity
  is_float = isFloat(low_RoverTInput);
  LSDRaster low_RoverTRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + low_RoverTInput), RasterFormat);
  LSDSoilHydroRaster low_RoverT = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(low_RoverTInput.c_str())) : LSDSoilHydroRaster(FilledDEM, low_RoverTRaster, 0);
  is_float = isFloat(hi_RoverTInput);
  LSDRaster hi_RoverTRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + hi_RoverTInput), RasterFormat);
  LSDSoilHydroRaster hi_RoverT = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(hi_RoverTInput.c_str())) : LSDSoilHydroRaster(FilledDEM, hi_RoverTRaster, 0);

  //soil_density
  is_float = isFloat(low_rhosInput);
  LSDRaster low_rhosRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + low_rhosInput), RasterFormat);
  LSDSoilHydroRaster low_rhos = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(low_rhosInput.c_str())) : LSDSoilHydroRaster(FilledDEM, low_rhosRaster, 0);
  is_float = isFloat(hi_rhosInput);
  LSDRaster hi_rhosRaster = (is_float == 1) ? LSDRaster() : LSDRaster((Path + hi_rhosInput), RasterFormat);
  LSDSoilHydroRaster hi_rhos = (is_float == 1) ? LSDSoilHydroRaster(FilledDEM, atof(hi_rhosInput.c_str())) : LSDSoilHydroRaster(FilledDEM, hi_rhosRaster, 0);

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

  //Link Stability index to polyline data
  vector<float> avg_si = SI.AverageSIs(UpslopeNodes, FlowInfo);
  LSDSoilHydroRaster AVGs =  SI.WriteAvgSIs(avg_si, UpslopeNodes, FlowInfo);

  //Write factor of safety to file
  Fs_min.write_raster((OutPath+"Fs_min_NC"), RasterFormat);
  Fs_max.write_raster((OutPath+"Fs_max_NC"), RasterFormat);

  //write Stability index to a file
  SI.write_raster((OutPath+"SI_NC"), RasterFormat);
  AVGs.write_raster((OutPath+"SI_AVG"), RasterFormat);

}
