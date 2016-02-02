//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Compile_Dreich_lengths.cpp
//
// Driver created to generate channel length measurements data at a range of input resolutions.
//
// Driver expects unfilled DEMs at a range of grid resolutions in the input directory 
//
// Run the driver with the following arguments:
//
// path to the DEM files with a trailing slash
// output path with trailing slash
// DEM filename Prefix
// file extension without the dot
//
// DEM naming convention should be <prefix>_<resolution>_DEM
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Stuart W.D. Grieve
// University of Edinburgh
// November 2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDFlowInfo.hpp"

int main(int nNumberofArgs, char *argv[])
{
  //Test for correct input arguments
	if (nNumberofArgs!=5)
	{
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path to the dems, the output path (with trailing slashes), the filename prefix and the DEM file format.";
		exit(EXIT_FAILURE);
	}

  //get input args
  string DEMpath = argv[1];
  string path = argv[2];
  string Prefix = argv[3];
  string DEM_Format = argv[4];

  //array of resolutions to load
  int Resolutions[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 60};

  // vectors to hold the length strings
  vector<string> Length_Data;

  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");


  for (int a = 0; a < 26; ++a){

	cout << "Processing DEM " << a+1 << " of " << "26" << endl;

	//load the DEM
        //build the string of the filename to load
	stringstream ss2;
	ss2 << DEMpath << Prefix << "_" << Resolutions[a] << "_DEM";
	LSDRaster DEM(ss2.str(), DEM_Format);
     
  //Fill 
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  Array2D<int> FlowDir = FlowInfo.get_FlowDirection();

  //load the channel network raster
  //build the string of the filename to load        
	stringstream ss3;
	ss3 << path << Prefix << "_" << Resolutions[a] << "_SO_DrEICH";
	LSDIndexRaster Chans(ss3.str(), DEM_Format);

  string tmp = FilledDEM.ChannelLengthByOrder(Chans, FlowDir);
  Length_Data.push_back(tmp);

  }


  //set up a writer to write the output data
  ofstream WriteData;

  //create an output filename based on the dem name
  stringstream ss;

  ss << path << Prefix << "_Dreich_Length.txt";
  WriteData.open(ss.str().c_str());

  for (int i = 0; i < int(Length_Data.size()); ++i){

    WriteData << Length_Data[i] << endl;

  }

  WriteData.close();

}
