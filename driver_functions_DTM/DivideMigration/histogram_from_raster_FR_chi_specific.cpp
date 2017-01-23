//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// histogram_from_raster.cpp
//
// This program makes a histogram of raster values
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// David T Milodowski
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../TNT/tnt.h"
int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

  string full_name = path_name+f_name;

  ifstream file_info_in;
  file_info_in.open(full_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  string raster1_name,raster2_name,raster3_name;
  string extension;
  float upper_limit,lower_limit;
  string temp;
  float bin_width;
  file_info_in >> temp >> raster1_name;
  file_info_in >> temp >> raster2_name;
  file_info_in >> temp >> raster3_name;
  file_info_in >> temp >> extension;
  file_info_in >> temp >> upper_limit;
  file_info_in >> temp >> lower_limit;
  file_info_in >> temp >> bin_width;
  file_info_in.close();
  
  LSDRaster raster1(raster1_name, extension);
  LSDRaster raster2(raster2_name, extension);
  LSDRaster raster3(raster3_name, extension);
  int NRows = raster1.get_NRows();
  int NCols = raster1.get_NCols();
  float NoDataValue = raster1.get_NoDataValue();
  
  int NRows_10m = raster3.get_NRows();
  int NCols_10m = raster3.get_NCols();
  
  // Loop through raster, collecting values into a vector
  vector<float> raster_values;
  for(int i = 0; i < NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      if(raster1.get_data_element(i,j)!=NoDataValue) raster_values.push_back(raster1.get_data_element(i,j));
    }
  }
  for(int i = 0; i < raster2.get_NRows(); ++i)
  {
    for(int j = 0; j < raster2.get_NCols(); ++j)
    {
      if(raster2.get_data_element(i,j)!=NoDataValue) raster_values.push_back(raster2.get_data_element(i,j));
    }
  }
  vector<float> raster_values_10m;
  for(int i = 0; i < NRows_10m; ++i)
  {
    for(int j = 0; j < NCols_10m; ++j)
    {
      if(raster3.get_data_element(i,j)!=NoDataValue) raster_values_10m.push_back(raster3.get_data_element(i,j));
    }
  }
  cout << raster_values.size() << endl;
  cout << raster_values_10m.size() << endl;
  vector<float> Midpoints,LLims,ULims,PD;
  vector<int> Count;
  vector<float> Midpoints_10m,LLims_10m,ULims_10m,PD_10m;
  vector<int> Count_10m;
  calculate_histogram_fixed_limits(raster_values, bin_width, lower_limit, upper_limit, Midpoints,LLims,ULims,Count,PD);
  calculate_histogram_fixed_limits(raster_values_10m, bin_width, lower_limit, upper_limit, Midpoints_10m,LLims_10m,ULims_10m,Count_10m,PD_10m);

  // Write ouput file
  cout << "writing output files" << endl;
  string output_filename = raster1_name + "_histogram.txt";
  ofstream ofs;
  ofs.open(output_filename.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Midpoint LowerLim UpperLim PD\n";
  for(int i = 0; i<Midpoints.size();++i)
  {
    ofs << Midpoints[i] << " " << LLims[i] << " " << ULims[i] << " " << PD[i] <<  " " << Count[i] << "\n";
  }
  ofs.close();

  string output_filename_10m = raster3_name + "_histogram.txt";
 
  ofs.open(output_filename_10m.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Midpoint LowerLim UpperLim PD\n";
  for(int i = 0; i<Midpoints_10m.size();++i)
  {
    ofs << Midpoints_10m[i] << " " << LLims_10m[i] << " " << ULims_10m[i] << " " << PD_10m[i] << " " << Count_10m[i] << "\n";
  }
  ofs.close();
}
