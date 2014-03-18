//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_map_sensitivity_analysis
//
// This program takes two arguments, the path name and the driver name
// It then runs the chi analysis for the parameter space defined in the driver
// file, allowing the user to examine the effects of changing m/n value,
// number of target nodes, minimum segment length, sigma and target skip value.
// At present it just spits out an output file for each iteration. In future, it
// will kick out the files in a nicer format :-).
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
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"


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

	string DEM_ID;
	string dem_ext = "_dem";
  string fill_ext = "_fill";
  string sources_ext = "_CH";
  string DEM_flt_extension = "flt";
  
  // initialise variables to be assigned from .driver file
  int junction_number;
	float pruning_threshold;
	int threshold;
	float A_0;
	float start_movern;
	float d_movern;
	float Minimum_Slope;
	int n_movern;
	int n_iterations;
  string temp;
  int start_min_segment_length, d_min_segment_length, n_min_segment_length;
  float start_Sigma, d_Sigma;
  int n_Sigma;
  int start_TargetNodes, d_TargetNodes, n_TargetNodes;
  int start_Skip, d_Skip, n_Skip;
                                                          
	file_info_in >> temp >> DEM_ID
               >> temp >> junction_number
               >> temp >> Minimum_Slope
               >> temp >> pruning_threshold 
               >> temp >> A_0              
               >> temp >> start_movern
               >> temp >> d_movern 
               >> temp >> n_movern
               >> temp >> n_iterations
               >> temp >> start_min_segment_length
               >> temp >> d_min_segment_length
               >> temp >> n_min_segment_length
               >> temp >> start_TargetNodes
               >> temp >> d_TargetNodes
               >> temp >> n_TargetNodes
               >> temp >> start_Sigma
               >> temp >> d_Sigma
               >> temp >> n_Sigma   
               >> temp >> start_Skip
               >> temp >> d_Skip
               >> temp >> n_Skip;
                   
	string jn_name = itoa(junction_number);
	string uscore = "_";
	jn_name = uscore+jn_name;
	file_info_in.close();

  cout << "PARAMETERS FOR SENSITIVITY ANALYSIS\n\t DEM_ID = " << DEM_ID
               << "\n\t Junction Number " << junction_number
               << "\n\t Minimum Slope (for fill function) " << Minimum_Slope
               << "\n\t Pruning Threshold (area) " <<  pruning_threshold 
               << "\n\t A_0 " <<  A_0              
               << "\n\t m/n: start " <<  start_movern << "; increment " << d_movern << "; number " << n_movern
               << "\n\t number of MC iterations " << n_iterations
               << "\n\t min segment length: start " <<  start_min_segment_length << "; increment " << d_min_segment_length << "; number " << n_min_segment_length
               << "\n\t target nodes: start " <<  start_TargetNodes << "; increment " << d_TargetNodes << "; number " << n_TargetNodes
               << "\n\t sigma: start " <<  start_Sigma << "; increment " << d_Sigma << "; number " << n_Sigma
               << "\n\t Skip: start " <<  start_Skip << "; increment " << d_Skip << "; number " << n_Skip << endl;

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// load the  DEM
	LSDRaster topo_test((path_name+DEM_ID+dem_ext), DEM_flt_extension);
  LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
  cout << "\t Flow routing..." << endl;
	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);

	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  cout << "\t Loading Sources..." << endl;
	// load the sources
	LSDIndexRaster Sources((path_name+DEM_ID+sources_ext), DEM_flt_extension);
  vector<int> sources;
	for(int i = 0; i < topo_test.get_NRows(); ++i)
	{
    for(int j = 0; j < topo_test.get_NCols(); ++j)
    {
      if(Sources.get_data_element(i,j) != Sources.get_NoDataValue()) sources.push_back(Sources.get_data_element(i,j));
    }
  }
  
//   sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);
  cout << "\t Got sources!" << endl;
	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

	// now get a junction and look for the longest channel upstream
	cout << "creating main stem" << endl;
	LSDIndexChannel main_stem = ChanNetwork.generate_longest_index_channel_in_basin(junction_number, FlowInfo, DistanceFromOutlet);
        cout << "got main stem channel, with n_nodes " << main_stem.get_n_nodes_in_channel() <<  endl;

	// now get the best fit m over n for all the tributaries
	int organization_switch = 1;                                                  // Decide how I want to thin the dataset
	int pruning_switch = 0;
	LSDIndexChannelTree ChannelTree(FlowInfo, ChanNetwork, junction_number, organization_switch,
                                        DistanceFromOutlet, pruning_switch, pruning_threshold);

	// print a file that can be ingested bt the chi fitting algorithm
	string Chan_fname = "_ChanNet";
	string Chan_ext = ".chan";
	string Chan_for_chi_ingestion_fname = path_name+DEM_ID+Chan_fname+jn_name+Chan_ext;
	ChannelTree.print_LSDChannels_for_chi_network_ingestion(FlowInfo, filled_topo_test, DistanceFromOutlet, Chan_for_chi_ingestion_fname);
	
// Create chi network automatically here, rather than reading in a file
	// create the chi network
	LSDChiNetwork ChiNetwork(Chan_for_chi_ingestion_fname);
	LSDChiNetwork ChiNetwork_extended(Chan_for_chi_ingestion_fname);
	ChiNetwork_extended.extend_tributaries_to_outlet();
	
  // Loop through range of target nodes to be tested
  for(int i_target_nodes = 0; i_target_nodes < n_TargetNodes; ++ i_target_nodes)
  {
    int target_nodes = start_TargetNodes + i_target_nodes*d_TargetNodes;
    // Loop through range of minimum segment lengths
    for(int i_min_seg = 0; i_min_seg < n_min_segment_length; ++ i_min_seg)
    {
      int minimum_segment_length = start_min_segment_length + i_min_seg*d_min_segment_length;
      // Loop through range of Sigma values
      for(int i_sigma = 0; i_sigma < n_Sigma; ++i_sigma)
      {
        float sigma = start_Sigma + float(i_sigma)*d_Sigma;
        //Loop through Skip values
        for(int i_skip = 0; i_skip < n_Skip; ++i_skip)
        {
          int target_skip = start_Skip + i_skip*d_Skip;
          
          // Now run the chi-analysis to contruct the best fit m/n value for each
          // scenario, and kick out chi-profiles for each iteration.
          cout << "Params for this iteration:\n\t TargetNodes = " << target_nodes << "\n\t MinSegmentLength = " << minimum_segment_length
                << "\n\t Sigma = " << sigma << "\n\t Skip = " << target_skip << endl;                 
        	// first get a string with some paramter values; needed for file names
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
        
          // Loop through range of m/n values
          for(int i = 0; i<n_movern; i++)
          {
        	  float this_movern = start_movern + float(i)*d_movern;
        
        	  cout << "This m/n is: " << this_movern << endl;
        
        	  string fpt_ext = ".tree";
        
        	  //label each profile file with the value of movn used
        	  //string prefix_movn = std::to_string(this_movern);
        
        	  // convert the m/n ratio to a string for the output filename
        	  string prefix_movn = static_cast<ostringstream*>( &(ostringstream() << this_movern) )->str();
        
        		// get the breaks of all the channels
           	ChiNetwork_extended.split_all_channels(A_0, this_movern, n_iterations,
        	  					target_skip, target_nodes, minimum_segment_length, sigma);
        
        	  // monte carlo sample all channels
        	  ChiNetwork_extended.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, this_movern, n_iterations,
        							target_skip, minimum_segment_length, sigma);
        
        	  string fpt_mc = "_fullProfileMC_forced_" + prefix_movn+param_str;
        
        	  ChiNetwork_extended.print_channel_details_to_file_full_fitted((path_name+DEM_ID+fpt_mc+param_str+jn_name+fpt_ext));
        	} 
                                                                               
          //--------------------------------------------------------------------
      	  // GET BEST FIT m/n value for basin using these parameters
      	  // get the best fit m over n ratio for this basin
        	float bf_cum_movn_ms, bf_colinear_movn_breaks;
        	string bfchi_fname = "_BFmovern";
        	string fpt2_ext = ".movern";
        
        	vector<float> m_over_n_values;
        	vector<float> AICc_mean_breaks;
        	vector<float> AICc_stdd_breaks;
        	vector< vector<float> > AICc_vals;
        	vector< vector<float> > AICc_stddev;
        
        	// get the best fit m over n in two different ways
        	int Monte_Carlo_switch = 1;
        	bf_colinear_movn_breaks = ChiNetwork.search_for_best_fit_m_over_n_colinearity_test_with_breaks(A_0,
        						 n_movern, d_movern,start_movern,
        						 minimum_segment_length, sigma, target_skip, target_nodes, n_iterations,
        						 m_over_n_values, AICc_mean_breaks, AICc_stdd_breaks,
        						 Monte_Carlo_switch);
        
        	bf_cum_movn_ms = ChiNetwork_extended.search_for_best_fit_m_over_n_individual_channels_with_breaks_monte_carlo(A_0,
        	                     n_movern, d_movern,start_movern,
        						 minimum_segment_length, sigma, target_skip, target_nodes, n_iterations,
        						 m_over_n_values, AICc_vals, AICc_stddev);
        
        	// now print the results from the m_over_n analysis
        	ofstream mn_analysis_out;
        	string mn_fo_fname = (path_name+DEM_ID+bfchi_fname+param_str+jn_name+fpt2_ext);
        	mn_analysis_out.open(mn_fo_fname.c_str());
        	mn_analysis_out << "-99 ";
        	for (int i = 0; i< n_movern; i++)
        	{
        		mn_analysis_out << m_over_n_values[i] << " ";
        	}
        	mn_analysis_out << endl << "-99 ";
        	for (int i = 0; i< n_movern; i++)
        	{
        		mn_analysis_out << AICc_mean_breaks[i] << " ";
        	}
        	mn_analysis_out << endl << "-99 ";
        	for (int i = 0; i< n_movern; i++)
        	{
        		mn_analysis_out << AICc_stdd_breaks[i] << " ";
        	}
        	mn_analysis_out << endl;
        
        
        	for(int chan = 0; chan<  ChiNetwork.get_n_channels(); chan++)
        	{
        		mn_analysis_out << chan << " ";
        		for (int i = 0; i< n_movern; i++)
        		{
        			mn_analysis_out << AICc_vals[chan][i] << " ";
        		}
        		mn_analysis_out << endl;
        		mn_analysis_out << chan << " ";
        		for (int i = 0; i< n_movern; i++)
        		{
        			mn_analysis_out << AICc_stddev[chan][i] << " ";
        		}
        		mn_analysis_out << endl;
        	}
        	mn_analysis_out.close();                                            
          //--------------------------------------------------------------------
        	cout << "best fit m over n collinear: " << bf_colinear_movn_breaks
               << " and mainstem: " << bf_cum_movn_ms << endl;
          //-------------------------------------------------------------------- 
        }
      }
    }
  }
}
