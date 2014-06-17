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
#include "../../LSDRaster.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDShapeTools.hpp"

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
  string BV_ext = "_baseline_values";
  string SE_ext = "_supraelevation";
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
  int NormaliseTransProfile = 1;
  int NormaliseLongProfile = 0;
  cout << "\n\t writing output \n\t\t - transverse profile" << endl;
  TestSwath.write_transverse_profile_to_file(RasterTemplate, percentiles, BinWidth, RasterTemplate_file.c_str(),NormaliseTransProfile);
  cout << "\t - longitudinal profile" << endl;
  TestSwath.write_longitudinal_profile_to_file(RasterTemplate, percentiles, BinWidth, RasterTemplate_file.c_str(),NormaliseLongProfile);
  cout << "\t - profile templates" << endl;
  LSDRaster Swath(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_DistanceToBaselineArray());  
  LSDRaster Long_Swath(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_DistanceAlongBaselineArray());
  LSDRaster BaselineValues(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),TestSwath.get_BaselineValueArray());
                  
  Array2D<float> ChannelElevation = TestSwath.get_BaselineValueArray();
  Array2D<float> SupraElevationArray(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_NoDataValue());
  
  for(int i = 0; i<RasterTemplate.get_NRows(); ++i)
  {
    for(int j = 0; j<RasterTemplate.get_NCols(); ++j)
    {
      if(ChannelElevation[i][j] != RasterTemplate.get_NoDataValue()) SupraElevationArray[i][j] = RasterTemplate.get_data_element(i,j) - ChannelElevation[i][j];
    }
  }
  
  LSDRaster SupraElevation(RasterTemplate.get_NRows(),RasterTemplate.get_NCols(),RasterTemplate.get_XMinimum(),RasterTemplate.get_YMinimum(),
                  RasterTemplate.get_DataResolution(),RasterTemplate.get_NoDataValue(),SupraElevationArray);
  
                  
  string output_file = RasterTemplate_file+Swath_ext;
  Swath.write_raster(output_file.c_str(),flt_ext);
  string output_file2 = RasterTemplate_file+Long_Swath_ext;
  Long_Swath.write_raster(output_file2.c_str(),flt_ext);
  string output_file3 = RasterTemplate_file+BV_ext;
  BaselineValues.write_raster(output_file3.c_str(),flt_ext);
	string output_file4 = RasterTemplate_file+SE_ext;
  SupraElevation.write_raster(output_file4.c_str(),flt_ext);
	cout << "S.I.G." << endl;
}
