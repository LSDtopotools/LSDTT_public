//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// GenerateShielding.cpp
//
// Driver created to generate shielding rasters for cosmo basins.
//
// Driver expects an unfilled DEM in the given directory of the given name and format and
// a text file with the name of each basin **without the extension** to be processed on its own line.
//
// Run driver with the following arguments:
//
// path to the DEM file with a trailing slash
// DEM filename without the file extension or dot
// file extension without the dot
// basin list filename
// theta step
// phi step
//
// A usage example is:
// nice ./GenerateShielding.out /home/s0675405/DataStore/NC/ NC_DEM flt Basins.txt 5 5
//
// Output data is stored in the input directory in a series of *.flt files called <Basin_Name>_SH.flt
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// July 06 2015
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
	if (nNumberofArgs!=7)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the DEM filename, the DEM file format, the basin list filename, and a theta and phi step."; 
		exit(EXIT_FAILURE);
	}

  //get input args
  string path = argv[1];
  string DEM_Name = argv[2];
  string DEM_Format = argv[3];
  string basins = argv[4];
  int theta_step = atoi(argv[5]);
  int phi_step = atoi(argv[6]);
    
  //load the DEM
  LSDRaster DEM((path+DEM_Name), DEM_Format);
    
  ifstream infile;  
                     
  stringstream ss;
  ss << path << basins;                
  infile.open(ss.str().c_str());
  
  //load the list of basin filenames 
  string temp;
  vector<string> BasinList;
  while (infile >> temp){
    BasinList.push_back(temp);
  }
  
  //loop over each basin, punch it out, trim it and generate the shielding raster.
  
  for (int q = 0; q < int(BasinList.size()); ++q){
    
    cout << "Processing " << q+1 << " of " << BasinList.size() << endl;
    
    //load basin raster
    LSDRaster Basin((path+BasinList[q]), DEM_Format);
    
    //process basin into a dem for shielding
    LSDRaster CutDEM = DEM.CookieCutRaster(Basin);
    LSDRaster Trimmed = CutDEM.RasterTrimmer();
    LSDRaster Shielded = Trimmed.TopographicShielding(theta_step, phi_step);
        
    //write the shielding raster to the working directory
    Shielded.write_raster((path+BasinList[q]+"_SH2"),DEM_Format);
  
  }
        
}
