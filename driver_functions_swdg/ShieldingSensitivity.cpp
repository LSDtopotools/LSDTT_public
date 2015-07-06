//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// ShieldingSensitivity.cpp
//
// Driver created to automate estimation of the corect theta,phi steps to use for the 
// topgraphic shielding calculations.
//
// Driver expects an unfilled DEM in the given directory of the given name and format.
//
// Run driver with the following arguments:
//
// path to the DEM file with a trailing slash
// DEM filename without the file extension or dot
// file extension without the dot
//
// A usage example is:
// nice ./ShieldingSensitivity.out /home/s0675405/DataStore/Final_Paper_Data/NC/ NC_DEM flt
//
// Output data is stored in the input directory in a series of data files called <DEM_Name>_<theta step>_<phi step>_shield.txt
// and can be plotted using ShieldingSensitivity.py found at https://github.com/sgrieve/GeneralAnalysis/
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// June 18 2015
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
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the DEM filename, and the DEM file format."; 
		exit(EXIT_FAILURE);
	}
	
  //get input args
  string path = argv[1];
  string DEM_Name = argv[2];
  string DEM_Format = argv[3];
    
  //load the DEM and fill it
  LSDRaster DEM((path+DEM_Name), DEM_Format);


  LSDRaster Trimmed = DEM.RasterTrimmer();
  
  //set up a writer to write the output data
  ofstream WriteData;
 
  //arraya of theta and phi steps to test 
  int PhiSteps[] = {1, 2, 3, 5, 8, 10, 15, 20, 30, 45, 60, 90};
  int ThetaSteps[] = {1, 2, 3, 5, 8, 10, 15, 20, 30, 45, 60, 90, 180, 360};
   
  for (int q = 11; q >= 0; --q){
    for (int w = 13; w >= 0; --w){
    
      cout << "Processing theta, phi pair of " << ThetaSteps[w] << " " << PhiSteps[q] << endl;
      cout << "Currently processing Theta " << (q+1) << " of " << "14 & Phi " << (w+1) << " of " << "12" << endl; 
     
      //create an output filename based on the dem name
      stringstream ss;
      ss << path << DEM_Name<< "_" << ThetaSteps[w] << "_" << PhiSteps[q] << "_shield.txt"; 
                                            
      LSDRaster Shielded = Trimmed.TopographicShielding(ThetaSteps[w], PhiSteps[q]);
      Shielded.FlattenToFile(ss.str().c_str());
   
    }
  }
                                          


}
