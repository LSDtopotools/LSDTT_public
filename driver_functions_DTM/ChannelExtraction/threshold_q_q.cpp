//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// threshold_quartile_quantile.cpp
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
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../TNT/tnt.h"
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

	string Raster_name;
	string Output_name;
	string q_q_filename;
	string DEM_extension = "flt";
	float standard_normal_variate;	
	string temp;
	int raster_selector;                                                        
	file_info_in >> temp >> Raster_name
	             >> temp >> Output_name
	             >> temp >> q_q_filename
               >> temp >> standard_normal_variate;
	// Now create the raster selection vector based on user's selection
	// Elevation
  LSDRaster raster(Raster_name, DEM_extension);
  vector<float> values;
  for(int i = 0; i < raster.get_NRows(); ++i)
  {
    for(int j = 0; j < raster.get_NCols(); ++j)
    {
      if(raster.get_data_element(i,j) != raster.get_NoDataValue())
      {
        values.push_back(raster.get_data_element(i,j));
      }
    }
  }
//   float threshold = get_quartile_quantile_value(values, standard_normal_variate);

  
  // write a q-q file for plotting - going from -3 <= standard normal variate < +3 (i.e. the 99.7% CI)
  vector<float> quantile_values,normal_variates,mn_values;
  int N_points = 1000;
  quantile_quantile_analysis(values, quantile_values, normal_variates, mn_values, N_points);
  ofstream ofs;
  ofs.open(q_q_filename.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "normal_variate value\n";
  for(int i = 0; i<quantile_values.size();++i)
  {
    ofs << normal_variates[i] << " " << quantile_values[i] << " " << mn_values[i] << "\n";
  }
  ofs.close();
  
  // Find q-q threshold
  cout << "\t finding deviation from Gaussian distribution to define q-q threshold" << endl;
  vector<int> indices;
  int flag = 0;
  float threshold_condition=0.99;
  float curvature_threshold = 0.1;
  for(int i = 0; i<quantile_values.size(); ++i)
  {
    if((normal_variates[i] >= 0) && (flag == 0))
    {
      if(mn_values[i]<threshold_condition*quantile_values[i])
      {
        flag = 1;
        curvature_threshold = quantile_values[i];
      }
    }
  }
  
  cout << curvature_threshold << endl;
  
  cout << "\t Creating channel mask based on curvature threshold" << endl;
  Array2D<float> binary_raster(raster.get_NRows(),raster.get_NCols(),raster.get_NoDataValue());
  for(int i = 0; i < raster.get_NRows(); ++i)
  {
    for(int j = 0; j < raster.get_NCols(); ++j)
    {
      if(raster.get_data_element(i,j) != raster.get_NoDataValue())
      {
        if(raster.get_data_element(i,j) >= curvature_threshold) binary_raster[i][j]=1;
        else binary_raster[i][j]=0;
      }
    }
  }
  LSDRaster output_raster(raster.get_NRows(),raster.get_NCols(),raster.get_XMinimum(), raster.get_YMinimum(), raster.get_DataResolution(), raster.get_NoDataValue(),binary_raster);
  output_raster.write_raster(Output_name,DEM_extension);
  cout << "DONE" << endl;
}
