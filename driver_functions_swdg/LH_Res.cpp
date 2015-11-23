//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LH_Res.cpp
//
// Driver created generate boxplot data for hillslope length and relief measurements
// taken at a range of grid resolutions. Filters the data to exclude any value
// below 2.0, as discussed in the hillslope length manuscript.
//
// Driver expects a drirectory of unfilled DEMs with the name format <prefix>_<Resolution>_DEM.<Format>,
// a channel heads file in the same directory with the format <prefix>_DEM_CH.csv.
//
// Run driver with the following arguments:
//
// path to the input files with a trailing slash
// filename prefix without an underscore
// DEM file format: bil, flt or asc
// window radius value in spatial units for surface fitting
// basin order the strahler order of basins to be extracted
// switch to write rasters 0 == do not write rasters and 1 == write rasters
//
// A usage example is:
// nice ./LH_Res.out /home/s0675405/Data/SC/ SC bil 4 2 1
//
// Outputs a pair of text files and a csv file for every resolution step.
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
		cout << "FATAL ERROR: wrong number of inputs. The program needs the path (with trailing slash), the filename prefix, the DEM file format, the window radius, ";
    cout << "basin order, and a switch to write rasters if desired." << endl;
		exit(EXIT_SUCCESS);
	}

  //load input arguments
  string Path = argv[1];
	string Prefix = argv[2];
  string DEM_Format = argv[3];
  float window_radius = atof(argv[4]); //6
  int BasinOrder = atoi(argv[5]);  //2
  int WriteRasters = atoi(argv[6]);  //0 (do not write rasters) or 1 (write rasters)

  //set up writers to write the output data
  stringstream ssLH;
  stringstream ssR;
  ssLH << Path << Prefix << "_LHResData.txt";
  ssR << Path << Prefix << "_RResData.txt";

  ofstream WriteLHData;
  WriteLHData.open(ssLH.str().c_str());
  ofstream WriteRData;
  WriteRData.open(ssR.str().c_str());

  //write headers
  WriteLHData << "resolution 2pc 25pc median mean 75pc 98pc minimum maximum" << endl;
  WriteRData << "resolution 2pc 25pc median mean 75pc 98pc minimum maximum" << endl;

  //set boundary conditions
  vector<string> BoundaryConditions(4, "No Flux");

  //array of resolutions to load
  int Resolutions[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 60};

  for (int a = 0; a < 26; ++a){
	  cout << "Processing DEM " << a+1 << " of " << "26" << endl;

    //create an output filename based on the dem name and the resolution
    stringstream ss;
    ss << Path << Prefix << "_" << Resolutions[a];
    string Filename = ss.str();

    //load dem
    LSDRaster DEM((Filename + "_DEM"), DEM_Format);

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

    int CurrentWindowSize = window_radius * Resolutions[a];

    vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(CurrentWindowSize, raster_selection);
    LSDRaster slope = Surfaces[1];
    LSDRaster aspect = Surfaces[2];

    cout << "\nGetting drainage network and basins\n" << endl;

    // get a flow info object
  	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

    //get stream net from channel heads
    vector<int> sources = FlowInfo.Ingest_Channel_Heads((Filename+"_DEM_CH"), "csv", 2); //this filename may be wrong
    LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
    LSDIndexRaster StreamNetwork = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();

    //Extract basins based on input stream order
    vector< int > basin_junctions = ChanNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
    LSDIndexRaster Basin_Raster = ChanNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);

    cout << "\nExtracting hilltops and hilltop curvature" << endl;

    // extract hilltops - no critical slope filtering is performed here
    LSDRaster hilltops = ChanNetwork.ExtractRidges(FlowInfo);

    //get hilltop curvature using filter to remove positive curvatures
    LSDRaster cht_raster = FilledDEM.get_hilltop_curvature(Surfaces[3], hilltops);
    LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(cht_raster);

    cout << "Starting hilltop flow routing\n" << endl;

    // these params do not need changed during normal use of the HFR algorithm
    bool print_paths_switch = false;
    int thinning = 1;
    string trace_path = "";
    bool basin_filter_switch = false;
    vector<int> Target_Basin_Vector;

    //run HFR
    vector< Array2D<float> > HFR_Arrays = FlowInfo.HilltopFlowRouting(FilledDEM, hilltops, slope, StreamNetwork, aspect, Filename, Basin_Raster, Surfaces[4], print_paths_switch, thinning, trace_path, basin_filter_switch, Target_Basin_Vector);

    LSDRaster HFR_LH = hilltops.LSDRasterTemplate(HFR_Arrays[1]);
    LSDRaster relief = hilltops.LSDRasterTemplate(HFR_Arrays[3]);

    //Filter Relief and LH data to remove any values < 2, as in Grieve et al (2015)
    LSDRaster LH = HFR_LH.RemoveBelow(2.0);
    LSDRaster Relief = relief.RemoveBelow(2.0);

    //go through the lh raster and get every value into a 1D vector
    vector<float> LH_vec = Flatten_Without_Nodata(LH.get_RasterData(), LH.get_NoDataValue());
    vector<float> Boxplot = BoxPlot(LH_vec);

    //write the values to the output file
    WriteLHData << Resolutions[a] << " " << Boxplot[0] << " " << Boxplot[1] << " " << Boxplot[2] << " " << Boxplot[3] << " " << Boxplot[4] << " " << Boxplot[5] << " " << Boxplot[6] << " " << Boxplot[7];

    WriteLHData << endl;

    //go through the relief raster and get every value into a 1D vector
    vector<float> R_vec = Flatten_Without_Nodata(Relief.get_RasterData(), Relief.get_NoDataValue());
    vector<float> Boxplot_R = BoxPlot(R_vec);

    //write the values to the output file
    WriteRData << Resolutions[a] << " " << Boxplot_R[0] << " " << Boxplot_R[1] << " " << Boxplot_R[2] << " " << Boxplot_R[3] << " " << Boxplot_R[4] << " " << Boxplot_R[5] << " " << Boxplot_R[6] << " " << Boxplot_R[7];

    WriteRData << endl;

    //if the user requests the rasters to be written, write the rasters
    if (WriteRasters == 1){
      cout << "Writing Rasters\n" << endl;

      Surfaces[1].write_raster((Filename + "_Slope"), DEM_Format);
      CHT.write_raster((Filename + "_CHT"), DEM_Format);
      LH.write_raster((Filename + "_HFR_LH"), DEM_Format);
      Relief.write_raster((Filename + "_Relief"), DEM_Format);

    }
  }
    WriteLHData.close();
    WriteRData.close();
}
