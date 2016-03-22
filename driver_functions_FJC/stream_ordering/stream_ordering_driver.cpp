//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// compare_slopes_at_channel_head.cpp
// make with comparea_slopes_at_channel_head.make
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb
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
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDStrahlerLinks.hpp"
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

	//read in the info from the parameter file  
	string DEM_ID;
	string Sources_name;
	float MinSlope;
	string DEM_extension = "bil";
	string temp;
	file_info_in >> temp >> DEM_ID
							 >> temp >> Sources_name
							 >> temp >> MinSlope;

	file_info_in.close();

	//load the DEM
//   cout << "\t Filling the DEM" << endl;
//    LSDRaster DEM(path_name+DEM_ID, DEM_extension);
//    float MinSlope = 0.0001;
//    LSDRaster FilledDEM = DEM.fill(MinSlope);
//    FilledDEM.write_raster((path_name+DEM_ID+"_fill"), DEM_extension);
	LSDRaster FilledDEM((path_name+DEM_ID+"_fill"), DEM_extension);

	//get a FlowInfo object
	cout << "\t Flow routing..." << endl;
	vector<string> BoundaryConditions(4, "No Flux");
	// get a flow info object
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

	//get the channel heads
	cout << "\t Ingesting channel heads" << endl;
	vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_ID+Sources_name), DEM_extension);
	
	//get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;

	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();

	string JI_name = "_JI";
	JIArray.write_raster((path_name+DEM_ID+JI_name),DEM_extension);
	
	//get the number of streams of each order and print to screen

	int max_SO = ChanNetwork.get_maximum_stream_order();
	cout << "Max stream order = " << max_SO << endl;
	for (int i = 1; i < int(max_SO+1); i++)
	{
		int n_streams = ChanNetwork.get_number_of_streams(FlowInfo, i);
		cout << "Stream order: " << i << " Number streams: " << n_streams << endl;
	}
  
  //TESTING STRAHLER LINKS CODE
  
  LSDStrahlerLinks Links(ChanNetwork, FlowInfo);
  vector<int> number_streams = Links.get_number_of_streams();
	
    
}
