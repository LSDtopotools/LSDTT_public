//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_get_profiles
//
// This program takes two arguments, the path name and the driver name
// It then reads the .chan file with the prefix listed in the first row of
// the driver file and creates a number of chi profiles that are forced
// with the m/n ratios determied by parameters in teh driver file on rows
// 9-11.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//
// Simon M. Mudd
// Declan Valters     (implemented m/n loop)
// University of Edinburgh
//
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
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name, and the junction number" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];
	string trunk_channel_junction = argv[3];

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

	string DEM_name;
	//string dem_ext = "_dem";
  string DEM_flt_extension = "flt";
  string fill_ext = "_fill";
  string sources_ext = "_CH";
	int temp;
	double pruning_threshold;
	int threshold;
	double A_0;
	int minimum_segment_length;
	double sigma;
	double start_movern;
	double d_movern;
	double Minimum_Slope;
	int n_movern;
	int target_nodes;
	int n_iterations;
	double fraction_dchi_for_variation;
	double vertical_interval;
	double horizontal_interval;
	double area_thin_frac;
	int target_skip;
                                                          // change junction_number to temp
	file_info_in >> DEM_name >> Minimum_Slope >> threshold >> temp
				>> pruning_threshold >> A_0 >> minimum_segment_length >> sigma >> start_movern
				>> d_movern >> n_movern >> target_nodes >> n_iterations >> fraction_dchi_for_variation
				>> vertical_interval >> horizontal_interval >> area_thin_frac >> target_skip;


	cout << "Paramters of this run: " << endl
		   << "DEM ID: " << DEM_name << endl
	     << "trunk channel junction number: " << trunk_channel_junction << endl
	     << "pruning_threshold: " << pruning_threshold << endl
	     << "threshold: " << threshold << endl
	     << "A_0: " << A_0 << endl
	     << "minimum_segment_length: " << minimum_segment_length << endl
	     << "sigma: " << sigma << endl
	     << "start_movern " << start_movern << endl
	     << "d_movern: " << d_movern << endl
	     << "n_movern: " << n_movern << endl
	     << "target_nodes: " << target_nodes << endl
	     << "n_iterarions: " << n_iterations << endl
	     << "fraction_dchi_for_variation: " << fraction_dchi_for_variation << endl
	     << "target_skip is: " << target_skip << endl;


	file_info_in.close();
	
  // Get list of junctions into a vector - this file is produced by ./chi_map_step1_getJunctions.out
  string ListOfJunctions = "tributary_junctions_for_chi_map_step2." + trunk_channel_junction + ".txt";
  vector<int> junction_list, junction_rows, junction_cols;
  int this_junction, junction_row,junction_col;
  string line;
  ifstream junction_input(ListOfJunctions.c_str());
  if (junction_input.is_open())
  {
    while (junction_input >> this_junction >> junction_row >> junction_col)
    {
      junction_list.push_back(this_junction);
      junction_rows.push_back(junction_row);
      junction_cols.push_back(junction_col);
    }
    junction_input.close();
  }
   
  int NumberOfJunctions = junction_list.size();
  
  // set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// load the filled DEM
  LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
//	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_flt_extension);
  double MinSlope = 0.0001;
  LSDRaster filled_topo_test = topo_test.fill(MinSlope);

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
  
	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	// load the sources
	LSDIndexRaster Sources((path_name+DEM_name+sources_ext), DEM_flt_extension);
  vector<int> sources;
  int NRows = Sources.get_NRows();
  int NCols = Sources.get_NCols();
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
 	
	// now get the junction network
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);
//   LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
//   JIArray.write_raster("testJI","flt");
//   LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
//   SOArray.write_raster("testSO","flt");
  for(int i_junction = 0; i_junction < NumberOfJunctions; ++ i_junction)
  {
  	cout << "RUNNING CHI ANALYSIS from junction " << i_junction+1 << " of " << NumberOfJunctions << ": " << junction_list[i_junction] << endl;
    //int junction_number = junction_list[i_junction];
    int junction_number = ChanNetwork.retrieve_junction_number_at_row_and_column(junction_rows[i_junction],junction_cols[i_junction]);
    string jn_name = itoa(junction_number);
	  string uscore = "_";
	  jn_name = uscore+jn_name;
    // now get a junction and look for the longest channel upstream
  	cout << "creating main stem" << endl;
  	LSDIndexChannel main_stem = ChanNetwork.generate_longest_index_channel_in_basin(junction_number, FlowInfo, DistanceFromOutlet);
    cout << "got main stem channel, with n_nodes " << main_stem.get_n_nodes_in_channel() <<  endl;
  
//  	string Basin_name = "_basin";
//  	LSDIndexRaster BasinArray = ChanNetwork.extract_basin_from_junction(junction_number,junction_number,FlowInfo);
//   	BasinArray.write_raster((path_name+DEM_name+Basin_name+jn_name),DEM_flt_extension);
  
  	// now get the best fit m over n for all the tributaries
  	int organization_switch = 1;                                                  // Decide how I want to thin the dataset
  	int pruning_switch = 0;
  	LSDIndexChannelTree ChannelTree(FlowInfo, ChanNetwork, junction_number, organization_switch,
                                          DistanceFromOutlet, pruning_switch, pruning_threshold);
  
  	// print a file that can be ingested bt the chi fitting algorithm
   	string Chan_fname = "_ChanNet";
   	string Chan_ext = ".chan";
  	string Chan_for_chi_ingestion_fname = path_name+DEM_name+Chan_fname+jn_name+Chan_ext;
   	ChannelTree.print_LSDChannels_for_chi_network_ingestion(FlowInfo, filled_topo_test, DistanceFromOutlet, Chan_for_chi_ingestion_fname);
                               
  	// get a string with some paramter values
  	string sigma_str;
  	string skip_str;
  	string msl_str;
  	string tn_str;
  	string param_str;
  	sigma_str = static_cast<ostringstream*>( &(ostringstream() << sigma) )->str();
  	skip_str = static_cast<ostringstream*>( &(ostringstream() << target_skip) )->str();
  	msl_str = static_cast<ostringstream*>( &(ostringstream() << minimum_segment_length) )->str();
  	tn_str = static_cast<ostringstream*>( &(ostringstream() << target_nodes) )->str();
  
  	param_str = uscore+sigma_str+uscore+skip_str+uscore+msl_str+uscore+tn_str;
  
  // Create chi network automatically here, rather than reading in a file
  	// create the chi network
  	LSDChiNetwork ChiNetwork(Chan_for_chi_ingestion_fname);
  	LSDChiNetwork ChiNetwork_extended(Chan_for_chi_ingestion_fname);
  	ChiNetwork_extended.extend_tributaries_to_outlet();
  
  
  // 	//=-=-=-=-=-=-=-=-=-=-=-=-
  // 	// now caluculate the slope area data
  // 	//=-=-=-=-=-=-=-=-=-=-=-=-
  // 	string SA_vert_fname = "SAVert";
  // 	string SA_ext = ".SAdata";
  // 	ChiNetwork.slope_area_extraction_vertical_intervals(vertical_interval, area_thin_frac,
  // 	                                              (path_name+DEM_name+SA_vert_fname+jn_name+SA_ext));
  // 
  // 	string SA_horiz_fname = "SAHoriz";
  // 	ChiNetwork.slope_area_extraction_horizontal_intervals(horizontal_interval, area_thin_frac,
  // 	                                              (path_name+DEM_name+SA_horiz_fname+jn_name+SA_ext));
  
  	//=-=-=-=-=-=-=-=-=
  	///////////////////
  	// Now do the segment finding
  	///////////////////
  	//=-=-=-=-=-=-=-=-=
  	//
  	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	// loop through MULTIPLE M/N vlaues
    // produces separate profile files for each value of m/n
  
    for(int i = 0; i<n_movern; i++)
    {
  	  double this_movern = start_movern + double(i)*d_movern;
  
  	  cout << "This m/n is: " << this_movern << endl;
  
  	  string fpt_ext = ".tree";
  
  	  //label each profile file with the value of movn used
  	  //string prefix_movn = std::to_string(this_movern);
  
  	  // convert the m/n ratio to a string for the output filename
  	  string prefix_movn = static_cast<ostringstream*>( &(ostringstream() << this_movern) )->str();
  
  		// get the breaks of all the channels
     	cout << "splitting extended chi network" << endl;
      ChiNetwork_extended.split_all_channels(A_0, this_movern, n_iterations,
  	  					target_skip, target_nodes, minimum_segment_length, sigma);
      cout << "Here we go... starting Monte Carlo sampling..." << endl;
  	  // monte carlo sample all channels
  	  ChiNetwork_extended.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, this_movern, n_iterations,
  							target_skip, minimum_segment_length, sigma);
  
  	  string fpt_mc = "_fullProfileMC_forced_" + prefix_movn+param_str;
  
  	  ChiNetwork_extended.print_channel_details_to_file_full_fitted((path_name+DEM_name+fpt_mc+jn_name+fpt_ext));
    }
  }
}
