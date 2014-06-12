//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// CloudBaseTest.cpp
// A test module for the LSDCloudBase module
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// David T. Milodowski
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDRaster.hpp"
#include "../LSDSwathProfile.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{
	if (nNumberofArgs != 5)
  {
    cout << "FATAL ERROR: You must provide 4 input arguements\n\t - an input file with the baseline coordinates\n\t - a template raster\n\t - the swath half width\n\t - the profile resolution\n\t" << endl;
    exit(EXIT_SUCCESS);
  }
  string Baseline_file = argv[1];
  string RasterTemplate_file = argv[2];
  float HalfWidth = atof(argv[3]);
  float BinWidth = atof(argv[4]);
  
  string flt_ext = "flt";  
  string Swath_ext = "_swath_trans";
  string Long_Swath_ext = "_swath_long";
  cout << "starting the test run... here we go!" << endl;
	RasterTemplate_file = RasterTemplate_file;
	
  cout << "\t loading template raster" << endl;
  LSDRaster RasterTemplate(RasterTemplate_file.c_str(),flt_ext);
  
  cout << "\t loading baseline points" << endl;
  PointData BaselinePoints = LoadShapefile(Baseline_file.c_str());
  
  cout << "\t creating swath template" << endl;
  LSDSwath TestSwath(BaselinePoints, RasterTemplate, HalfWidth);
  vector<float> percentiles;
  percentiles.push_back(0);
  percentiles.push_back(25);
  percentiles.push_back(50);
  percentiles.push_back(75);
  percentiles.push_back(100);
  cout << "\t writing output" << endl;
  TestSwath.write_transverse_profile_to_file(RasterTemplate, percentiles, BinWidth, RasterTemplate_file.c_str());
  TestSwath.write_longitudinal_profile_to_file(RasterTemplate, percentiles, BinWidth, RasterTemplate_file.c_str());
  LSDRaster Swath(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_DistanceToBaselineArray());  
  LSDRaster Long_Swath(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_DistanceAlongBaselineArray());
  string output_file = RasterTemplate_file+Swath_ext;
  Swath.write_raster(output_file.c_str(),flt_ext);
  string output_file2 = RasterTemplate_file+Long_Swath_ext;
  Long_Swath.write_raster(output_file2.c_str(),flt_ext);
	cout << "S.I.G." << endl;
}
