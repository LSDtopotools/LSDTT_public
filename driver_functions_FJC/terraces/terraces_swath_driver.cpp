//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// terraces_swath_driver.cpp
// Extract information about terraces using a shapefile of the main stem channel.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb
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
	if (nNumberofArgs != 3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }
  string path_name = argv[1];
  string f_name = argv[2];
  cout << "The path name is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if (file_info_in.fail())
	{
		cout << "\nFATAL ERROR: the header file\"" << full_name
				 << "\" doesn't exist" << endl;
	}

	string DEM_ID, Baseline_file, temp;
	string CC_ext = "_CC";
	string DEM_extension = "bil";
	float HalfWidth;

	// read in the parameters
	file_info_in >> temp >> DEM_ID
						   >> temp >> Baseline_file
							 >> temp >> HalfWidth;

	file_info_in.close();

  string Swath_ext = "_swath_trans";
  string Long_Swath_ext = "_swath_long";
  string BV_ext = "_baseline_values";
  cout << "starting the test run... here we go!" << endl;

  cout << "\t Loading the DEM" << endl;
  LSDRaster Elevation((path_name+DEM_ID), DEM_extension);

	cout << "\t Loading the terraces" << endl;
	LSDRaster ConnectedComponents((path_name+DEM_ID+CC_ext), DEM_extension);

  cout << "\t loading baseline points" << endl;
  PointData BaselinePoints = LoadShapefile(path_name+Baseline_file.c_str());

  cout << "\t creating swath template" << endl;
  LSDSwath TestSwath(BaselinePoints, ConnectedComponents, HalfWidth);
  vector<float> percentiles;
  percentiles.push_back(0);
  percentiles.push_back(25);
  percentiles.push_back(50);
  percentiles.push_back(75);
  percentiles.push_back(100);
  int NormaliseTransProfile = 1;
  int NormaliseLongProfile = 1;
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
  string output_file = RasterTemplate_file+Swath_ext;
  Swath.write_raster(output_file.c_str(),flt_ext);
  string output_file2 = RasterTemplate_file+Long_Swath_ext;
  Long_Swath.write_raster(output_file2.c_str(),flt_ext);
	cout << "S.I.G." << endl;
}
