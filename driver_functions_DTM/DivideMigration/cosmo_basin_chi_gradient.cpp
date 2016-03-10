//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_map_sensitivity_analysis
//
// This program takes two arguments, the path name and the driver name
// It then runs the chi analysis for the parameter space defined in the driver
// file, allowing the user to examine the effects of changing m/n value,
// number of target nodes, minimum segment length, sigma and target skip value.
// At present it just spits out an output file for each iteration.
// Now also kicks out a table with the parameters and best fit m/n.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// 18/03/2014
// David Milodowski
// University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../../LSDStatsTools.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDRasterSpectral.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDBasin.hpp"


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

  string DEM_ID,raster_ext,ChiGradientArray_name,data_file_name;
  string dem_ext = "_dem";
  string fill_ext = "_fill";
  string acc_ext = "_ACC";
  string txt_extension = ".txt";
  
  // initialise variables to be assigned from .driver file
  string DATA_DIR, OUTPUT_DIR, CHeads_file;
  int threshold; 
  float minimum_stream_order;
  float A_0;
  float movern;
  float Minimum_Slope;
  int n_iterations;
  string temp;
  int minimum_segment_length;
  float sigma;
  int target_nodes;
  int skip;
  
  file_info_in >> temp >> DATA_DIR
               >> temp >> OUTPUT_DIR
               >> temp >> DEM_ID
               >> temp >> CHeads_file
               >> temp >> ChiGradientArray_name
               >> temp >> data_file_name
               >> temp >> raster_ext
               >> temp >> Minimum_Slope;
  file_info_in.close();

  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";

  // load the DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID+dem_ext), raster_ext);

  int NRows = topography_raster.get_NRows();
  int NCols = topography_raster.get_NCols();
  float NoData = topography_raster.get_NoDataValue();
  float XMin =topography_raster.get_XMinimum(); 
  float YMin = topography_raster.get_YMinimum();
  float Resolution = topography_raster.get_DataResolution();
  
  LSDRaster filled_topography = topography_raster.fill(Minimum_Slope);
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);
  // calcualte the distance from outlet
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);;
  cout << "\t Got sources!" << endl;
  
  // now get the junction network
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);
  
  // Load in Chi Gradient Array
  LSDRaster ChiGradientArray(OUTPUT_DIR+ChiGradientArray_name,raster_ext);
  
  // Read in cosmo basin summary file
  vector<float> CosmoErate,CosmoError;
  float Erate,Error;
  vector<string> SampleNames;
  int junc;
  vector<int> BasinJunctions;
  string name;
  float rate,error,sd;
  int count = 0;

  ifstream data_file_in;
  data_file_in.open(data_file_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  data_file_in >> temp >> temp >> temp >> temp;
  while( data_file_in >> name >> Erate >> Error >> junc) 
  {
    ++count;
    cout << flush << "\t\t line " << count << "\r";
    CosmoErate.push_back(Erate);
    CosmoError.push_back(Error);
    SampleNames.push_back(name);
    BasinJunctions.push_back(junc);
  }
  data_file_in.close();
  cout << "\t\t Data read in complete" << endl;  

  int N_basins = BasinJunctions.size();
  vector<float> MeanChiGradient,SDChiGradient,SErrChiGradient,CatchmentArea;
  for(int i = 0; i<N_basins; ++i)
  {  
    vector<float> ChiGradients;
    int basin_junction = BasinJunctions[i];
    // get basin nodes  
//     int outlet_node = JunctionNetwork.get_Node_of_Junction(basin_junction);
    int outlet_node = JunctionNetwork.get_penultimate_node_from_stream_link(basin_junction, FlowInfo);
    vector<int> donor_nodes = FlowInfo.get_upslope_nodes(outlet_node);
    int N_nodes = donor_nodes.size();
    cout << SampleNames[i] << "; number of upstream nodes: " << N_nodes << endl;
    Array2D<float> basin(NRows,NCols,NoData);
    for(int i_node = 0; i_node < N_nodes; ++i_node)
    {
      int row, col;
      FlowInfo.retrieve_current_row_and_col(donor_nodes[i_node], row, col);
      basin[row][col] = 1;
      if(ChiGradientArray.get_data_element(row,col) != NoData)
      {
        ChiGradients.push_back(ChiGradientArray.get_data_element(row,col));
      }
    }
    MeanChiGradient.push_back(get_mean(ChiGradients));
    SDChiGradient.push_back(get_standard_deviation(ChiGradients,get_mean(ChiGradients)));
    SErrChiGradient.push_back(get_standard_error(ChiGradients,get_mean(ChiGradients)));
    CatchmentArea.push_back(float(N_nodes)*Resolution*Resolution);
    LSDRaster basin_raster(NRows,NCols,XMin,YMin,Resolution,NoData,basin);
    basin_raster.write_raster(OUTPUT_DIR+DEM_ID+"_"+SampleNames[i],"flt");
  } 

  // Now write output file
  cout << "writing output file" << endl;
  ofstream ofs;
  string outfile = OUTPUT_DIR+DEM_ID+"_cosmo_basin_chi_gradient.txt";
  ofs.open(outfile.c_str());
  if(ofs.fail()){
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "SampleID CatchmentArea CosmoErate CosmoError MeanChiGradient SD SErr\n";
  for(int i = 0; i < N_basins; ++i)
  {
    ofs << SampleNames[i] << " " << CatchmentArea[i] << " " << CosmoErate[i] << " " << CosmoError[i] << " " << MeanChiGradient[i] << " " << SDChiGradient[i] << " " << SErrChiGradient[i] << "\n"; 
  }
  ofs.close();
}
