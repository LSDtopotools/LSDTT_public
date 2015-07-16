//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LH_Driver.cpp
//
// Driver created to automate extraction of hillslope lengths and basin metrics from DEM 
// files.
//
// Driver expects an unfilled DEM in the given directory in flt format with the name format <prefix>_DEM.flt,
// a channel heads file in the same directory with the format <prefix>_CH.flt and a floodplain raster
// <prefix>_FloodPlain.flt if floodplains are to be excluded from the analysis.
//
// Run driver with the following arguments:
//
// path to the input files with a trailing slash
// filename prefix without an underscore
// window radius value in spatial units for surface fitting
// basin order the strahler order of basins to be extracted
// switch to exclude floodplains (1) or not exclude them (0) 
// switch to write rasters 0 == do not write rasters and 1 == write rasters
//
// A usage example is:
//nice ./LH_Driver.out /home/s0675405/DataStore/Final_Paper_Data/NC/ NC 7 2 1 0
//nice ./LH_Driver.out /home/s0675405/DataStore/Final_Paper_Data/PA/ PA 5 2 1 0
//nice ./LH_Driver.out /home/s0675405/DataStore/Final_Paper_Data/CA/ CA 12 2 1 0
//nice ./LH_Driver.out /home/s0675405/DataStore/Final_Paper_Data/OR/ OR 7.5 2 1 0
// 
//
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
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDChannel.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannel.hpp"
#include "../../LSDMostLikelyPartitionsFinder.hpp"
#include "../../LSDBasin.hpp"
#include "../../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{
  
  //Test for correct input arguments
  if (nNumberofArgs!=7)
  {
    cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the filename prefix, window radius, "; 
    cout << "basin order, a switch to use or exclude floodplains and a switch to write rasters if desired." << endl;
    exit(EXIT_SUCCESS);
  }
  
  //load input arguments
  string path = argv[1];
  string filename = argv[2];
  float window_radius = atof(argv[3]); //6
  int BasinOrder = atoi(argv[4]);  //2
  //int FloodplainSwitch = atoi(argv[5]);
  int MultipixelSwitch = atoi(argv[5]);
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
  raster_selection.push_back(1); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(1); //plan curvature
  raster_selection.push_back(0); 
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(window_radius, raster_selection);
  LSDRaster slope = Surfaces[1];
  LSDRaster aspect = Surfaces[2];   
  
  cout << "\nGetting drainage network and basins\n" << endl;

  // get a flow info object
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);
  
  //get stream net from channel heads
//   vector<int> sources = FlowInfo.Ingest_Channel_Heads((path+filename+"_dem_CH"), "flt"); //swap to csv?   
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path+ChannelHeads_name), "csv",2);
  LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  LSDIndexRaster StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();

  //load floodplain and merge with the channel network if required, otherwise the 
  //floodplain mask will only contain the channel data
//   LSDIndexRaster ChannelAndFloodplain;
// 
//   if (FloodplainSwitch == 1){
//     LSDIndexRaster Floodplains((path+filename+"_FloodPlain"), "flt");
//     ChannelAndFloodplain = StreamNetwork.MergeChannelWithFloodplain(Floodplains);
//   }
//   else{
//     ChannelAndFloodplain = StreamNetwork;
//   }
//   
  LSDIndexRaster MultipixelNetwork;
  if (MultipixelSwitch == 1){
    LSDIndexRaster MultipixelChannels((path+filename+"_channel_mask"), "flt");
    MultipixelNetwork = StreamNetwork.MergeChannelWithFloodplain(MultipixelChannels);
  }
  else{
    MultipixelNetwork = StreamNetwork;
  }

  // load rock exposure raster
  LSDRaster RockExposure((path+filename+"_rock_exposure"), "flt");
  
  //Extract basins based on input stream order
  vector< int > basin_junctions = ChanNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
  LSDIndexRaster Basin_Raster = ChanNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);
  
  cout << "\nExtracting hilltops and hilltop curvature" << endl;
  
  // extract hilltops - no critical slope filtering is performed here
  LSDRaster hilltops = ChanNetwork.ExtractRidges(FlowInfo);   
     
  //get hilltop curvature using filter to remove positive curvatures         
  LSDRaster cht_raster = FilledDEM.get_hilltop_curvature(Surfaces[3], hilltops);
  LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(cht_raster);  
  
  //get d infinity flowdirection and flow area
  Array2D<float> dinf = FilledDEM.D_inf_FlowDir();
  LSDRaster dinf_rast = FilledDEM.LSDRasterTemplate(dinf);
  LSDRaster DinfArea = FilledDEM.D_inf_units();
  
  cout << "Starting hilltop flow routing\n" << endl;
  
  //start of Hilltop flow routing
  string prefix = (path+filename);  //set a path to write the hillslope length data to, based on the input path and filename given by the user
  
  // these params do not need changed during normal use of the HFR algorithm
  bool print_paths_switch = false;
  int thinning = 1;
  string trace_path = "";
  bool basin_filter_switch = false;
  vector<int> Target_Basin_Vector;

  //run HFR    
  vector< Array2D<float> > HFR_Arrays = FlowInfo.HilltopFlowRoutingBedrock(FilledDEM, hilltops, slope, MultipixelNetwork, aspect, prefix, Basin_Raster, Surfaces[4], print_paths_switch, thinning, trace_path, basin_filter_switch, Target_Basin_Vector, RockExposure);
   
  LSDRaster HFR_LH = hilltops.LSDRasterTemplate(HFR_Arrays[1]);
  LSDRaster HFR_Slope = hilltops.LSDRasterTemplate(HFR_Arrays[2]);
  LSDRaster relief = hilltops.LSDRasterTemplate(HFR_Arrays[3]);
  //end of HFR 

  //create lsdbasin objects in a loop over each junction number
  vector< LSDBasin > Basins;

  //slope area plotting parameters - these defaults are usually fine
  float log_bin_width = 0.1;
  int SplineResolution = 10000;
  int bin_threshold = 0;
  float CriticalSlope = 1.2; //this needs modified to generate proper E*R* data
  
  cout << "\nCreating each LSDBasin" << endl;
  
  //loop over each basin, generating an LSDBasin object which contains that basin's measurements
  for (int w = 0; w < int(basin_junctions.size()); ++w){
    
    cout << (w+1) << " / " << basin_junctions.size() << endl;
    
    LSDBasin Basin(basin_junctions[w], FlowInfo, ChanNetwork);
    Basin.set_FlowLength(StreamNetwork, FlowInfo);
    Basin.set_DrainageDensity();
    Basin.set_all_HillslopeLengths(FlowInfo, HFR_LH, slope, DinfArea, log_bin_width, SplineResolution, bin_threshold);
    Basin.set_SlopeMean(FlowInfo, slope);
    Basin.set_AspectMean(FlowInfo, aspect);
    Basin.set_ElevationMean(FlowInfo, FilledDEM);
    Basin.set_ReliefMean(FlowInfo, relief);
    Basin.set_CHTMean(FlowInfo, CHT);
    Basin.set_EStar_RStar(CriticalSlope);
    Basin.set_BedrockFraction(FlowInfo, RockExposure);
    
    Basins.push_back(Basin);
                             
  }
 
  //create a filestream to write the output data
  // use the input arguments to generate a path and filename for the output file
  ofstream WriteData;                 
  stringstream ss;
  ss << path << filename << "_PaperData.txt";                
  WriteData.open(ss.str().c_str());

  //write headers
  WriteData << "BasinID HFR_mean HFR_median HFR_stddev HFR_stderr HFR_Nvalues HFR_range HFR_min HFR_max SA_binned_LH SA_Spline_LH LH_Density Area Basin_Slope_mean Basin_Slope_median Basin_Slope_stddev Basin_Slope_stderr Basin_Slope_Nvalues Basin_Slope_range Basin_Slope_min Basin_Slope_max Basin_elev_mean Basin_elev_median Basin_elev_stddev Basin_elev_stderr Basin_elev_Nvalues Basin_elev_Range Basin_elev_min Basin_elev_max Aspect_mean CHT_mean CHT_median CHT_stddev CHT_stderr CHT_Nvalues CHT_range CHT_min CHT_max EStar RStar HT_Slope_mean HT_Slope_median HT_Slope_stddev HT_Slope_stderr HT_Slope_Nvalues HT_Slope_range HT_Slope_min HT_Slope_max HT_relief_mean HT_relief_median HT_relief_stddev HT_relief_stderr HT_relief_Nvalues HT_relief_range HT_relief_min HT_relief_max BedrockFraction" << endl;

  cout << "\nWriting data to file\n" << endl;

  //write all data to the opened file, ensuring that there are data points to be written in each basin                                         
  for (int q = 0; q < int(Basins.size()); ++q){
    // only work where we have data points
    if (Basins[q].CalculateNumDataPoints(FlowInfo, HFR_LH) != 0 && Basins[q].CalculateNumDataPoints(FlowInfo, slope) != 0 && Basins[q].CalculateNumDataPoints(FlowInfo, FilledDEM) != 0 && Basins[q].CalculateNumDataPoints(FlowInfo, CHT) != 0 && Basins[q].CalculateNumDataPoints(FlowInfo, HFR_Slope) != 0 && Basins[q].CalculateNumDataPoints(FlowInfo, relief) != 0){
      
      // BasinID
      WriteData << Basins[q].get_Junction()<< " ";
      
      //HFR
      WriteData << Basins[q].get_HillslopeLength_HFR() << " " << Basins[q].CalculateBasinMedian(FlowInfo, HFR_LH) << " " << Basins[q].CalculateBasinStdDev(FlowInfo, HFR_LH) << " " << Basins[q].CalculateBasinStdError(FlowInfo, HFR_LH) << " " << Basins[q].CalculateNumDataPoints(FlowInfo, HFR_LH) << " " << Basins[q].CalculateBasinRange(FlowInfo, HFR_LH) << " " << Basins[q].CalculateBasinMin(FlowInfo, HFR_LH) << " " << Basins[q].CalculateBasinMax(FlowInfo, HFR_LH)<< " ";         
      
      //SA_Bins
      WriteData << Basins[q].get_HillslopeLength_Binned()<< " ";
      
      //SA_Spline
      WriteData << Basins[q].get_HillslopeLength_Spline()<< " ";
      
      //Density
      WriteData << Basins[q].get_HillslopeLength_Density()<< " ";
      
      //Area
      WriteData << Basins[q].get_Area()<< " ";
      
      //Slope_Basin
      WriteData << Basins[q].get_SlopeMean() << " " << Basins[q].CalculateBasinMedian(FlowInfo, slope) << " " << Basins[q].CalculateBasinStdDev(FlowInfo, slope) << " " << Basins[q].CalculateBasinStdError(FlowInfo, slope) << " " << Basins[q].CalculateNumDataPoints(FlowInfo, slope) << " " << Basins[q].CalculateBasinRange(FlowInfo, slope) << " " << Basins[q].CalculateBasinMin(FlowInfo, slope) << " " << Basins[q].CalculateBasinMax(FlowInfo, slope)<< " ";
      
      //Elev_Basin
      WriteData << Basins[q].get_ElevationMean() << " " << Basins[q].CalculateBasinMedian(FlowInfo, FilledDEM) << " " << Basins[q].CalculateBasinStdDev(FlowInfo, FilledDEM) << " " << Basins[q].CalculateBasinStdError(FlowInfo, FilledDEM) << " " << Basins[q].CalculateNumDataPoints(FlowInfo, FilledDEM) << " " << Basins[q].CalculateBasinRange(FlowInfo, FilledDEM) << " " << Basins[q].CalculateBasinMin(FlowInfo, FilledDEM) << " " << Basins[q].CalculateBasinMax(FlowInfo, FilledDEM)<< " ";
      
      //Aspect_Basin
      WriteData << Basins[q].get_AspectMean()<< " ";
      
      //CHT
      WriteData << Basins[q].get_CHTMean() << " " << Basins[q].CalculateBasinMedian(FlowInfo, CHT) << " " << Basins[q].CalculateBasinStdDev(FlowInfo, CHT) << " " << Basins[q].CalculateBasinStdError(FlowInfo, CHT) << " " << Basins[q].CalculateNumDataPoints(FlowInfo, CHT) << " " << Basins[q].CalculateBasinRange(FlowInfo, CHT) << " " << Basins[q].CalculateBasinMin(FlowInfo, CHT) << " " << Basins[q].CalculateBasinMax(FlowInfo, CHT) << " ";
      
      //EStar
      WriteData << Basins[q].get_EStar()<< " ";
      
      //RStar
      WriteData << Basins[q].get_RStar()<< " ";
      
      //Slope_mean
      WriteData << Basins[q].CalculateBasinMean(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateBasinMedian(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateBasinStdDev(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateBasinStdError(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateNumDataPoints(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateBasinRange(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateBasinMin(FlowInfo, HFR_Slope) << " " << Basins[q].CalculateBasinMax(FlowInfo, HFR_Slope)<< " ";
      
      //Relief_mean
      WriteData << Basins[q].get_ReliefMean() << " " << Basins[q].CalculateBasinMedian(FlowInfo, relief) << " " << Basins[q].CalculateBasinStdDev(FlowInfo, relief) << " " << Basins[q].CalculateBasinStdError(FlowInfo, relief) << " " << Basins[q].CalculateNumDataPoints(FlowInfo, relief) << " " << Basins[q].CalculateBasinRange(FlowInfo, relief) << " " << Basins[q].CalculateBasinMin(FlowInfo, relief) << " " << Basins[q].CalculateBasinMax(FlowInfo, relief)<< " ";
  
      //BedrockFraction
      WriteData << Basins[q].get_BedrockFraction() << "\n";
    }    
  } 
   
  // close the output file
  WriteData.close();

  //if the user requests the raster to be written, write the rasters
  if (WriteRasters == 1){
    cout << "Writing Rasters\n" << endl;                                   
    FilledDEM.write_raster((path+filename+"_Fill"), "flt");
    Surfaces[1].write_raster((path+filename+"_Slope"),"flt");
    Surfaces[2].write_raster((path+filename+"_Aspect"),"flt");
    Surfaces[3].write_raster((path+filename+"_Curvature"),"flt");
    StreamNetwork.write_raster((path+filename+"_STNET"), "flt"); 
    Basin_Raster.write_raster((path+filename+"_Basins"), "flt"); 
    CHT.write_raster((path+filename+"_CHT"),"flt");
    HFR_LH.write_raster((path+filename+"_HFR_LH"),"flt"); 
    HFR_Slope.write_raster((path+filename+"_HFR_SLP"),"flt");
    relief.write_raster((path+filename+"_Relief"),"flt");
    
    //perform a hillshade
    LSDRaster Hillshade = FilledDEM.hillshade(45.0,315.0,1.0);
    Hillshade.write_raster((path+filename+"_HS"),"flt");
  
  }
}
