//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Elevation_Profiles.cpp
//
// Driver created to automate extraction of hillslope elevation profiles.
//
// Driver expects an unfilled DEM in the given directory in flt format with the name format <prefix>_DEM.flt and
// a channel heads file in the same directory with the format <prefix>_CH.flt
//
// Run driver with the following arguments:
//
// path to the input files with a trailing slash
// path to write the trace files with a trailing slash *folder must exist*
// filename prefix without an underscore
// window radius value in spatial units for surface fitting
// basin order the strahler order of basins to be extracted
// switch to exclude floodplains (1) or not exclude them (0) 
// switch to write rasters 0 == do not write rasters and 1 == write rasters
//
// A usage example is:
// nice ./Elevation_Profiles.out /home/s0675405/DataStore/Final_Paper_Data/NC/ /home/s0675405/DataStore/Final_Paper_Data/NC/traces/ NC 7 2 0
//                                /home/s0675405/DataStore/Rory_Data/                  /home/s0675405/DataStore/Rory_Data/traces/
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Stuart W.D. Grieve
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDChannel.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannel.hpp"
#include "../LSDMostLikelyPartitionsFinder.hpp"
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{
  
  //Test for correct input arguments
	if (nNumberofArgs!=7)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path to the input data (with trailing slash), ";
    cout << "the filename prefix, the path to write the trace data to (with trailing slash), a window radius, "; 
    cout << "basin order, a switch to use or exclude floodplains and a switch to write rasters if desired." << endl;
		exit(EXIT_SUCCESS);
	}
  
  //load input arguments
  string path = argv[1];
	string trace_path = argv[2];
  string filename = argv[3];
  float window_radius = atof(argv[4]);
  int BasinOrder = atoi(argv[5]);  //2
  int WriteRasters = atoi(argv[6]);  //0 (do not write rasters) or 1 (write rasters)
   
  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  //load dem
  LSDRaster DEM((path+filename+"_DEM"), "flt");  
  
  //Fill 
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);
  
  //surface fitting
  vector<int> raster_selection;
  
  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope 
  raster_selection.push_back(0); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(0); //plan curvature
  raster_selection.push_back(0); 
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(window_radius, raster_selection);
  LSDRaster slope = Surfaces[1];  
  
  cout << "\nGetting drainage network and basins\n" << endl;

  // get a flow info object
	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);
  
  //get stream net from channel heads
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path+filename+"_DEM_CH"), "flt"); //swap to csv? 
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  LSDIndexRaster StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
                                                         
  //Extract basins based on input stream order
  vector< int > basin_junctions = ChanNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
  LSDIndexRaster Basin_Raster = ChanNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);
  
  cout << "\nExtracting hilltops" << endl;
  // extract hilltops - no critical slope filtering is performed here
  LSDRaster hilltops = ChanNetwork.ExtractRidges(FlowInfo); 
  
  //Array2D<float> Fake_HT(FilledDEM.get_NRows(),FilledDEM.get_NCols(),-9999);
  
  //int fake_i = 2500;
  //int fake_j = 2632;
  //Fake_HT[fake_i][fake_j] = 1;
  
  //LSDRaster hilltops = FilledDEM.LSDRasterTemplate(Fake_HT);   
     
  //get d infinity flowdirection and flow area
  Array2D<float> dinf = FilledDEM.D_inf_FlowDir();
  LSDRaster dinf_rast = FilledDEM.LSDRasterTemplate(dinf);
    
  cout << "Starting hilltop flow routing\n" << endl;
  //start of Hilltop flow routing
  string prefix = (path+filename);  //set a path to write the hillslope length data to, based on the input path and filename given by the user
  
  // these params do not need changed during normal use of the HFR algorithm
  bool print_paths_switch = true;
  int thinning = 1;
  bool basin_filter_switch = false;                          
  vector<int> Target_Basin_Vector;

  //run HFR    
  vector< Array2D<float> > HFR_Arrays = FlowInfo.HilltopFlowRouting_Profile(FilledDEM, hilltops, slope, StreamNetwork, dinf_rast, prefix, Basin_Raster, print_paths_switch, thinning, trace_path, basin_filter_switch, Target_Basin_Vector);

  //end of HFR 

  

  //if the user requests the raster to be written, write the rasters
  if (WriteRasters == 1){
    cout << "Writing Rasters\n" << endl;                                   
    FilledDEM.write_raster((path+filename+"_Fill"), "flt");
    Surfaces[1].write_raster((path+filename+"_Slope"),"flt");
    Surfaces[2].write_raster((path+filename+"_Aspect"),"flt");
    Surfaces[3].write_raster((path+filename+"_Curvature"),"flt");
    StreamNetwork.write_raster((path+filename+"_STNET"), "flt"); 
    Basin_Raster.write_raster((path+filename+"_Basins"), "flt"); 
         
    //perform a hillshade
    LSDRaster Hillshade = FilledDEM.hillshade(45.0,315.0,1.0);
    Hillshade.write_raster((path+filename+"_HS"),"flt");
  
  }
}
