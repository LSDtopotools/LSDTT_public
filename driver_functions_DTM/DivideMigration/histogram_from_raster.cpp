//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// sample_roughness_lengthscales.cpp
//
// This program is used for testing the LSDRaster object
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

	string DEM_name;
	string dem_ext = "_dem";
	string class_ext = "_class";
	string chan_MD_ext = "_channel_mask";
  vector<string> output_extensions;
	output_extensions.push_back("_s1");
	output_extensions.push_back("_s2");
	output_extensions.push_back("_s3");
  string DEM_flt_extension = "flt";
	
  float surface_fitting_window_radius,normal_variability_radius;
	int rock_class,soil_class,veg_class;
	vector<int> raster_selection;  	
	string temp;
	int raster_selector;
	float bin_width;
	file_info_in >> temp >> DEM_name;
  file_info_in >> temp >> surface_fitting_window_radius;
  file_info_in >> temp >> normal_variability_radius;
  file_info_in >> temp >> rock_class;
  file_info_in >> temp >> soil_class;
  file_info_in >> temp >> veg_class;
  file_info_in >> temp >> bin_width;
           
  file_info_in.close();

  stringstream ss;
  ss << surface_fitting_window_radius << ' ' << normal_variability_radius;
  string lengthscale_1, lengthscale_2; 
  ss >> lengthscale_1 >> lengthscale_2;
  string lengthscale_ext = "_" + lengthscale_1 + "_" + lengthscale_2;
  
  
  cout << "\tRock class = " << rock_class << "\n\tSoil class = " << soil_class << "\n\tVeg class = " << veg_class << endl;
	LSDRaster dem(DEM_name+dem_ext, DEM_flt_extension);
	LSDIndexRaster classification((DEM_name+class_ext),DEM_flt_extension);
	LSDRaster MDChannels((DEM_name+chan_MD_ext),DEM_flt_extension);
  int NRows = dem.get_NRows();
	int NCols = dem.get_NCols();
	float NoDataValue = dem.get_NoDataValue();
  
  string s1_file = DEM_name+output_extensions[0]+lengthscale_ext;
  string s2_file = DEM_name+output_extensions[1]+lengthscale_ext;
  string s3_file = DEM_name+output_extensions[2]+lengthscale_ext;
  LSDRaster s1(s1_file, DEM_flt_extension);
  LSDRaster s2(s2_file, DEM_flt_extension);
  LSDRaster s3(s3_file, DEM_flt_extension);

  
  // Loop through DEM, avoiding edges, buffered by proximity to margin.  Loading
  // Store row and column
  vector<int> RockRowMajorIndex;
  vector<int> SoilRowMajorIndex;
  vector<int> RowIndex(NRows*NCols,int(NoDataValue));
  vector<int> ColIndex(NRows*NCols,int(NoDataValue));
  cout << "\t retrieving cell indices for each class" << endl;
  for(int i = (0+normal_variability_radius+1); i < (NRows-(normal_variability_radius+1)); ++i)
  {
    for(int j = (0+normal_variability_radius+1); j < (NCols-(normal_variability_radius+1)); ++j)
    {
      if(dem.get_data_element(i,j)!=NoDataValue && classification.get_data_element(i,j)==rock_class && MDChannels.get_data_element(i,j)==0) RockRowMajorIndex.push_back(NCols*i + j);
      if(dem.get_data_element(i,j)!=NoDataValue && classification.get_data_element(i,j)==soil_class && MDChannels.get_data_element(i,j)==0) SoilRowMajorIndex.push_back(NCols*i + j);
      if(dem.get_data_element(i,j)!=NoDataValue && classification.get_data_element(i,j)==veg_class && MDChannels.get_data_element(i,j)==0) SoilRowMajorIndex.push_back(NCols*i + j);
    }
  }
  for(int i = 0; i < NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      RowIndex[NCols*i + j]=i;
      ColIndex[NCols*i + j]=j;
    }
  }
  int N_rock = RockRowMajorIndex.size();
  int N_soil = SoilRowMajorIndex.size();
  cout << "\t sampling data...\n\t\t N rock-class = " << N_rock << "\n\t\t N soil-class = " << N_soil << endl;  
  // Sample RowMajorIndex
  vector<int> RockRowMajorIndexSample = RockRowMajorIndex;//sample_without_replacement(RockRowMajorIndex, N);
  vector<int> SoilRowMajorIndexSample = SoilRowMajorIndex;//sample_without_replacement(SoilRowMajorIndex, N);
  
  // EXTRACT ROUGHNESS METRICS -> ROCK EXPOSURE  
  cout << "SAMPLING SURFACE ROUGHNESS METRICS FOR ROCK EXPOSURE" << endl;
  vector<float> s1_sample_rock(N_rock,NoDataValue);
  vector<float> s2_sample_rock(N_rock,NoDataValue);
  vector<float> s3_sample_rock(N_rock,NoDataValue);

  // Loop through the sample indices, and store the values in a vector for 
  // each eigenvalue
  for(int i = 0; i < N_rock; ++i)
  {
    // unwrap row and col
    int row = RowIndex[RockRowMajorIndexSample[i]];
    int col = ColIndex[RockRowMajorIndexSample[i]]; 
    s1_sample_rock[i]=s1.get_data_element(row,col);
    s2_sample_rock[i]=s2.get_data_element(row,col);
    s3_sample_rock[i]=s3.get_data_element(row,col);
   if(s1_sample_rock[i]<0 || s2_sample_rock[i]<0 || s3_sample_rock[i]<0) cout << RockRowMajorIndexSample[i] << " " << row << " " << col << " " << s1_sample_rock[i] << " " << s2_sample_rock[i] << " " << s3_sample_rock[i] << endl;
  }                           
  vector<float> Midpoints,LLims,ULims,PD_s1_soil,PD_s2_soil,PD_s3_soil,PD_s1_rock,PD_s2_rock,PD_s3_rock;
  vector<int> Count;
  
  calculate_histogram_fixed_limits(s1_sample_rock, bin_width, 0, 1, Midpoints,LLims,ULims,Count,PD_s1_rock);
  calculate_histogram_fixed_limits(s2_sample_rock, bin_width, 0, 1, Midpoints,LLims,ULims,Count,PD_s2_rock);
  calculate_histogram_fixed_limits(s3_sample_rock, bin_width, 0, 1, Midpoints,LLims,ULims,Count,PD_s3_rock);
    
  // EXTRACT ROUGHNESS METRICS -> SOIL
  // Now repeat for bedrock exposure
  cout << "SAMPLING SURFACE ROUGHNESS METRICS FOR SOIL" << endl;
  vector<float> s1_sample_soil(N_soil,NoDataValue);
  vector<float> s2_sample_soil(N_soil,NoDataValue);
  vector<float> s3_sample_soil(N_soil,NoDataValue);

  // Loop through the sample indices, and store the values in a vector for 
  // each eigenvalue
  for(int i = 0; i < N_soil; ++i)
  {
    // unwrap row and col
    int row = RowIndex[SoilRowMajorIndexSample[i]];
    int col = ColIndex[SoilRowMajorIndexSample[i]]; 
    s1_sample_soil[i]=s1.get_data_element(row,col);
    s2_sample_soil[i]=s2.get_data_element(row,col);
    s3_sample_soil[i]=s3.get_data_element(row,col);
   if(s1_sample_soil[i]<0 || s2_sample_soil[i]<0 || s3_sample_soil[i]<0) cout << SoilRowMajorIndexSample[i] << " " << row << " " << col << " " << s1_sample_soil[i] << " " << s2_sample_soil[i] << " " << s3_sample_soil[i] << endl;
  }
  calculate_histogram_fixed_limits(s1_sample_soil, bin_width, 0, 1, Midpoints,LLims,ULims,Count,PD_s1_soil);
  calculate_histogram_fixed_limits(s2_sample_soil, bin_width, 0, 1, Midpoints,LLims,ULims,Count,PD_s2_soil);
  calculate_histogram_fixed_limits(s3_sample_soil, bin_width, 0, 1, Midpoints,LLims,ULims,Count,PD_s3_soil);
  
  // Write ouput file
  string output_filename = "co_" + lengthscale_ext+"_roughness_lengthscales.txt";
  ofstream ofs;
  ofs.open(output_filename.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "Midpoint LowerLim UpperLim PD_s1_rock PD_s2_rock PD_s3_rock PD_s1_soil PD_s2_soil PD_s3_soil\n";
  for(int i = 0; i<Midpoints.size();++i)
  {
    ofs << Midpoints[i] << " " << LLims[i] << " " << ULims[i] << " " << PD_s1_rock[i] << " " << PD_s2_rock[i] << " " << PD_s3_rock[i] << " " << PD_s1_soil[i] << " " << PD_s2_soil[i] << " " << PD_s3_soil[i] << "\n";
  }
  ofs.close();
}
