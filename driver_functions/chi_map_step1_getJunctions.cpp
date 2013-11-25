//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Chile_test.cpp
//
// This program is used for testing the LSDRaster object
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Simon M. Mudd
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
#include "../../LSDChannelNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../TNT/tnt.h"
int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
	if (nNumberofArgs!=2)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the target junction number" << endl;
		exit(EXIT_SUCCESS);
	}

	int TRUNK_CHANNEL_JUNCTION = atoi(argv[1]); 
	
  string DEM_name = "../cr_data/cr_dem";
	string DEM_flt_extension = "flt";
	string DEM_asc_extension = "asc";

	LSDRaster topo_test(DEM_name, DEM_flt_extension);
	//cout << "LINE 41" << endl;

	double MinSlope = 0.0001;
	int NRows = topo_test.get_NRows();
	int NCols = topo_test.get_NCols();
	double NoDataValue = topo_test.get_NoDataValue();
	double XMinimum = topo_test.get_XMinimum();
	double YMinimum = topo_test.get_YMinimum();
	double DataResolution = topo_test.get_DataResolution();
	
	
  
  LSDRaster filled_topo_test = topo_test.fill(MinSlope);
// 
// 	string DEM_outname = "../cr_data/cr_dem_fill";
//   filled_topo_test.write_raster(DEM_outname,DEM_flt_extension);

	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	string SO_name = "../cr_data/cr_dem_SO";
// 	string J_name = "cr_J";
 	string JI_name = "../cr_data/cr_dem_JI";
// 	string FDist_name = "cr_FDist";

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
  	//string FI_fname = "_flowinfo";
	//FlowInfo.pickle((path_name+DEM_name+FI_fname));
	// calculate the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
// 	DistanceFromOutlet.write_raster(FDist_name,DEM_flt_extension);

	// LOAD SOURCES!!!
	LSDIndexRaster Sources("../cr_data/cr_dem_CH",DEM_flt_extension);
  //int threshold = 10000;
 	vector<int> sources;
 	for(int i = 0; i < NRows; ++i)
 	{
    for(int j=0; j< NCols; ++j)
    {
      if(Sources.get_data_element(i,j)==1) // source here
      {
        sources.push_back(FlowInfo.retrieve_node_from_row_and_column(i, j));
      }
    } 
  }
 	
//   LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
//   sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);
  

  // GENERATE CHANNEL NETWORK
	cout << "Creating channel network from sources" << endl;
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);

// 	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
//  	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
// 
// 	SOArray.write_raster(SO_name,DEM_flt_extension);
//   JIArray.write_raster(JI_name,DEM_flt_extension);

  // LOCATE THE JUNCTIONS FOR THE TARGET TRIBUTARIES
	// now get a junction and look for the longest channel upstream
	cout << "Getting tributary junctions" << endl;
	int pruning_switch = 0;
	double pruning_threshold = 25000;
  vector<int> target_tributary_juncs = ChanNetwork.get_pruned_tributaries_from_main_stem(FlowInfo, ChanNetwork, TRUNK_CHANNEL_JUNCTION, DistanceFromOutlet, pruning_switch, pruning_threshold);
//   vector<int> target_tributary_juncs, node_of_tributaries;
//   LSDIndexChannel main_stem = ChanNetwork.generate_longest_index_channel_in_basin(TRUNK_CHANNEL_JUNCTION,FlowInfo, DistanceFromOutlet); 
//   ChanNetwork.extract_tributary_junctions_to_main_stem(main_stem, FlowInfo, target_tributary_juncs, node_of_tributaries);
  // Check that basins are not truncated by edge of study area
  cout << "testing basins for truncations" << endl; 
  vector<int> target_tributary_juncs_filt, target_tributary_juncs_rows, target_tributary_juncs_cols; 
  Array2D<int> TargetJuncs(NRows,NCols, NoDataValue);
  int node, row, col;
  for(int i = 0; i < int(target_tributary_juncs.size()); ++i)
  {
    //cout << "a" << endl;
    node = 	ChanNetwork.get_Node_of_Junction(target_tributary_juncs[i]);
    //cout << "b" << endl;
    bool IsTruncated = ChanNetwork.node_tester(FlowInfo, target_tributary_juncs[i]);
    
    if(IsTruncated == false)
    {
      FlowInfo.retrieve_current_row_and_col(node,row,col);
      target_tributary_juncs_rows.push_back(row);
      target_tributary_juncs_cols.push_back(col);
      target_tributary_juncs_filt.push_back( target_tributary_juncs[i] );

      TargetJuncs[row][col] = target_tributary_juncs[i];
    }
  }
  cout << "\t selected " << target_tributary_juncs_filt.size() << " tributaries for chi-analysis" << endl;

  
  
	// WRITE OUTPUT FILE OF JUNCTIONS TO BE ANALYSED
  stringstream ss;
  ss << TRUNK_CHANNEL_JUNCTION;
  string trunk_channel_junc = ss.str();
  LSDIndexRaster TargetJunctions(NRows,NCols,XMinimum,YMinimum,DataResolution,int(NoDataValue),TargetJuncs);
	TargetJunctions.write_raster(("cr_TargetJunctions"+trunk_channel_junc),"flt");
  string file_prefix = "tributary_junctions_for_chi_map_step2.";
	string file_extension = ".txt";  	
  string fname = file_prefix+trunk_channel_junc+file_extension;
  
  cout << "Printing junctions to: " << fname << endl;
  
  ofstream ofs;
  ofs.open(fname.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
		exit(EXIT_FAILURE);
  }
  
  for(int i = 0; i < int(target_tributary_juncs_filt.size()); ++i)
  {
    ofs << target_tributary_juncs_filt[i] << " " << target_tributary_juncs_rows[i] << " " << target_tributary_juncs_cols[i] << "\n";
  }
  
  ofs.close();

}
