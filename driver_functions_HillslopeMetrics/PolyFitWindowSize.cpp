//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// PolyFitWindowSize.cpp
//
// Driver created to automate estimation of the corect window size to use for the 
// polyfit routine.
//
// Driver expects an unfilled DEM in the given directory of the given name and format.
//
// Run driver with the following arguments:
//
// path to the DEM file with a trailing slash
// DEM filename without the file extension or dot
// file extension without the dot

// A usage example is:
// nice ./PolyFitWindowSize.out /home/s0675405/DataStore/Final_Paper_Data/NC/ NC_DEM flt
//
// Output data is stored in the input directory in a file called <DEM_Name>_Window_Size_Data.txt
// and can be plotted using WindowSize.py  found at https://github.com/sgrieve/GeneralAnalysis/
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// June 2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDShapeTools.hpp"

int main(int nNumberofArgs, char *argv[])
{
  //Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the DEM filename, and tghe DEM file format."; 
		exit(EXIT_FAILURE);
	}
	
  //get input args
  string path = argv[1];
  string DEM_Name = argv[2];
  string DEM_Format = argv[3];
    
  //surface fitting
  vector<int> raster_selection;
  raster_selection.push_back(0);
  raster_selection.push_back(0);  
  raster_selection.push_back(0); 
  raster_selection.push_back(1); //curvature

  //load the DEM and fill it
  LSDRaster DEM((path+DEM_Name), DEM_Format);
	float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  //set up a writer to write the output data
  ofstream WriteData;
  
  //create an output filename based on the dem name
  stringstream ss;
  ss << path << DEM_Name<< "_Window_Size_Data.txt";
  
  WriteData.open(ss.str().c_str());

  //write headers
  WriteData << "Length_scale Curv_mean Curv_stddev Curv_iqr" << endl;
  
  //array of window sizes to test 
  //this is a little arbritrary but reflects the sizes used by Roering et al. 2010
  int WindowSizes[] = {1, 2, 3, 4, 5, 7, 8, 10,15,20,25,50,70,90,100};
  
  // floats to hold the stats about the fitted surface
  float Curv_mean;
  float Curv_stddev;
  float Curv_median;
  float Curv_UpperQuartile;
  float Curv_LowerQuartile;
  float Curv_MaxValue;
  float Curv_iqr;
  vector<float> Curv_vec;
  
  for (int w = 0; w < 15; ++w){
    
    cout << "Processing surface " << w+1 << " of " << "15" << endl;
                                            
    vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(WindowSizes[w], raster_selection);
    LSDRaster curvature = Surfaces[3];   
  
    //reset values for next run
    Curv_mean = 0;
    Curv_stddev = 0;
    Curv_median = 0;
    Curv_UpperQuartile = 0;
    Curv_LowerQuartile = 0;
    Curv_MaxValue = 0;
    Curv_iqr = 0;
    Curv_vec.clear();  
  
    //go through the landscape and get every curvature value into a 1D vector    
    for (int i = 0; i < int(curvature.get_NRows()); ++i){
      for (int j = 0; j < int(curvature.get_NCols()); ++j){
        if (curvature.get_data_element(i,j) != curvature.get_NoDataValue()){
          Curv_vec.push_back(curvature.get_data_element(i,j));
        }
      }
    }  
  
    //calculate the std dev and iqr 
    get_distribution_stats(Curv_vec, Curv_mean, Curv_median, Curv_UpperQuartile, Curv_LowerQuartile, Curv_MaxValue);
    Curv_stddev = get_standard_deviation(Curv_vec, Curv_mean);
    Curv_iqr = Curv_UpperQuartile - Curv_LowerQuartile;

    //write the values to the output file  
    WriteData << WindowSizes[w] << " " << Curv_mean << " " << Curv_stddev << " " << Curv_iqr << endl; 
    
    }
                                          
  WriteData.close();

}
