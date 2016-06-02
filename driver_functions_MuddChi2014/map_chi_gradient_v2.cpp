//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_map_sensitivity_analysis version 2
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
#include "../LSDStatsTools.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDBasin.hpp"
#include "../LSDChiTools.hpp"
#include "../LSDShapeTools.hpp"


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
    cout << "\nFATAL ERROR: The parameter \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  string DEM_ID,raster_ext;
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
  int threshold_pixels_for_chi;
  int threshold_contributing_pixels;
  int minimum_basin_size_pixels;
  bool test_drainage_boundaries = true;

  file_info_in >> temp >> DATA_DIR
               >> temp >> OUTPUT_DIR
               >> temp >> DEM_ID
               >> temp >> CHeads_file
               >> temp >> raster_ext
               >> temp >> minimum_stream_order
               >> temp >> Minimum_Slope
               >> temp >> A_0
               >> temp >> movern
               >> temp >> n_iterations
               >> temp >> target_nodes
               >> temp >> minimum_segment_length
               >> temp >> sigma
               >> temp >> skip
               >> temp >> threshold_pixels_for_chi
               >> temp >> threshold_contributing_pixels
               >> temp >> minimum_basin_size_pixels
               >> temp >> test_drainage_boundaries;


  file_info_in.close();

  // make sure the threshold contributing pixel data is not nonsense
  if(threshold_contributing_pixels<=0)
  {
    cout << "You gave me a nonsense threshold contributing pixel value: " << threshold_contributing_pixels << endl;
    cout << "Defaulting to 50" << endl;
    threshold_contributing_pixels = 50;
  }


  cout << "PARAMETERS FOR Chi mapping\n\t DEM_ID = " << DEM_ID
               << "\n\t DATA_DIR " << DATA_DIR
               << "\n\t OUTPUT_DIR " << OUTPUT_DIR
               << "\n\t CHeads_file " << CHeads_file
               << "\n\t raster_ext " << raster_ext
               << "\n\t Basin Stream Order " << minimum_stream_order
               << "\n\t Minimum Slope (for fill function) " << Minimum_Slope
               << "\n\t A_0 " <<  A_0
               << "\n\t m/n: " <<  movern
               << "\n\t number of MC iterations " << n_iterations
               << "\n\t min segment length: " <<  minimum_segment_length
               << "\n\t target nodes: " <<  target_nodes
               << "\n\t sigma: " <<  sigma
               << "\n\t Skip: " <<  skip 
               << "\n\t threshold_pixels_for_chi: " << threshold_pixels_for_chi 
               << "\n\t test_drainage_boundaries: " << test_drainage_boundaries 
               << "\n\t minimum_basin_size_pixels: " << minimum_basin_size_pixels
               << "\n\t threshold_contributing_pixels: " << threshold_contributing_pixels
               << endl << endl;

  // Additional parameters for chi analysis
  int organization_switch = 1;
  int pruning_switch = 3;
  float pruning_threshold;// = target_stream_order - 1;
  
  cout << "Setting boundary conditions" << endl;

  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";
  
  // check the DEM
  cout << "I am checking to make sure the DEM exists." << endl;

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;
  
  int NRows = topography_raster.get_NRows();
  int NCols = topography_raster.get_NCols();
  float NoData = topography_raster.get_NoDataValue();
  float XMin =topography_raster.get_XMinimum();
  float YMin = topography_raster.get_YMinimum();
  float Resolution = topography_raster.get_DataResolution();
  map<string,string> GRS = topography_raster.get_GeoReferencingStrings();

  cout << "Filling topography." << endl;
  LSDRaster filled_topography = topography_raster.fill(Minimum_Slope);
  
  string filled_raster_name = DATA_DIR+DEM_ID+"_Fill";
  
  //filled_topography.write_raster(filled_raster_name,raster_ext);
  
  cout << "\t Flow routing..." << endl;
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

  // calcualte the distance from outlet
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  // calculate the flow accumulation
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();
  
  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file != "null")
  {
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, threshold_contributing_pixels);
    
    cout << "The number of sources is: " << sources.size() << endl;
    
  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t Got sources!" << endl;
  }

  // now get the junction network
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

  string SO_ext = "_SO";
  string JI_ext = "_JI";
  LSDIndexRaster SOArray = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
  SOArray.write_raster((OUTPUT_DIR+DEM_ID+SO_ext),raster_ext);
  JIArray.write_raster((OUTPUT_DIR+DEM_ID+JI_ext),raster_ext);
  int max_stream_order = JunctionNetwork.get_maximum_stream_order();
  Array2D<float> ChiGradientArray(NRows,NCols,NoData);
  
  cout << "The maximum stream order is: " << max_stream_order << endl;

  // need to get base-level nodes , otherwise these catchments will be missed!
  vector< int > BaseLevelJunctions_Initial = JunctionNetwork.get_BaseLevelJunctions();
  vector< int > BaseLevelJunctions;
  int N_BaseLevelJuncs = BaseLevelJunctions_Initial.size();
  cout << "The number of base level junctions is: " << N_BaseLevelJuncs << endl; 

  // remove basins drainage from edge if that is what the user wants
  
  cout << "Right, let me check the drainage basins. " << endl;
  if (test_drainage_boundaries)
  {
    cout << "Test_drainage_bondaries: " << test_drainage_boundaries << endl;
  
    cout << endl << endl << "I am going to remove any basins draining to the edge." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge(BaseLevelJunctions_Initial,FlowInfo); 
  }
  else
  {
    cout << "I'm not going to remove any drainage basins drainage to the edge." << endl;
    BaseLevelJunctions = BaseLevelJunctions_Initial;
  }
  
  // remove base level junctions for which catchment is too small
  BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions, 
                                              FlowInfo, FlowAcc, minimum_basin_size_pixels);

  // Correct number of base level junctions
  N_BaseLevelJuncs = BaseLevelJunctions.size();

  // calculate chi for the entire DEM
  float threshold_area_for_chi = Resolution*Resolution*float(threshold_pixels_for_chi);
  LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,threshold_area_for_chi);
  
  string chi_coord_string = OUTPUT_DIR+DEM_ID+"_chi_coord";
  chi_coordinate.write_raster(chi_coord_string,raster_ext);

  
  // now use a ChiTool object to print the chi tree to csv
  LSDChiTools ChiTool(FlowInfo);
  //ChiTool.chi_map_to_csv(FlowInfo, chi_csv_fname, A_0, movern, threshold_area_for_chi);
  
  // the precision needs to be high for the lat long and UTM coordinates to acutally
  // give a reasonable number of significant digits. 
  //chitree_out.precision(12);
  //chitree_out.open(csv_fname.c_str());
  //chitree_out << "id,x,y,lat,long,chi,mean_chi_gradient" << endl;
  
    // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;
  
  // now get the overlapping channels from the junction network file
  vector<int> source_nodes;
  vector<int> outlet_nodes;
  int n_nodes_to_visit = 10;
  JunctionNetwork.get_overlapping_channels(FlowInfo, BaseLevelJunctions, DistanceFromOutlet, 
                                    source_nodes,outlet_nodes,n_nodes_to_visit);


  string sources_fname = OUTPUT_DIR+DEM_ID+"_SourcesList.csv";
  FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(source_nodes, sources_fname);

  string outlet_fname = OUTPUT_DIR+DEM_ID+"_OutletList.csv";
  FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(outlet_nodes, outlet_fname);
  
  
  ChiTool.chi_map_automator(FlowInfo, source_nodes, outlet_nodes,
                            topography_raster, DistanceFromOutlet, 
                            DrainageArea, chi_coordinate, target_nodes, 
                            n_iterations, skip, minimum_segment_length, sigma);
  

  string csv_fname = OUTPUT_DIR+DEM_ID+"_ChiSegmented_tree.csv";
  string csv_full_fname = OUTPUT_DIR+DEM_ID+"_ChiFullSegmented_tree.csv";
  string chi_csv_fname =  OUTPUT_DIR+DEM_ID+"_Chi.csv";
  ChiTool.print_data_maps_to_file_basic(FlowInfo, csv_fname);
  ChiTool.print_data_maps_to_file_full(FlowInfo, csv_full_fname);


}
