//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexChannelTree.cpp
// source code for the LSDIndexChannelTree object
// this object spawns vectors of LSDIndexChannels
// they can be indexed by the LSDCahnnel network, but can also be independant of the
// channel network, storing longest channels from sources, for example
// This object is designed to be flexible, it can be used either with the
// LSDFlowInfo or LSDChannelNetwork object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// <your name here>
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		30/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include <string>
#include <fstream>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDChannel.hpp"
#include "LSDIndexChannelTree.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDChannelNetwork.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannelTree_CPP
#define LSDIndexChannelTree_CPP


void LSDIndexChannelTree::create(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork, int starting_junction)
{
  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  DataResolution = FlowInfo.get_DataResolution();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  NoDataValue = ChannelNetwork.get_NoDataValue();

  	organization_switch = 0;

	// get the upstream junctions
	upstream_junction_list = ChannelNetwork.get_upslope_junctions(starting_junction);

	// get the number of upstream junctions
	int n_us_junctions = upstream_junction_list.size();

	int this_junction;
	outlet_junction = starting_junction;
	outlet_node = ChannelNetwork.get_Node_of_Junction(starting_junction);

	for(int i = 0; i< n_us_junctions; i++)
	{
		this_junction = upstream_junction_list[i];
		IndexChannelVector.push_back( ChannelNetwork.generate_link_index_channel_from_junction(this_junction,FlowInfo) );
	}

}


void LSDIndexChannelTree::create(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
					      	int starting_junction, int org_switch, LSDRaster& DistanceFromOutlet)
{
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	DataResolution = FlowInfo.get_DataResolution();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	NoDataValue = ChannelNetwork.get_NoDataValue();

  	organization_switch = org_switch;

	if (organization_switch == 0)
	{
		// get the upstream junctions
		upstream_junction_list = ChannelNetwork.get_upslope_junctions(starting_junction);

		// get the number of upstream junctions
		int n_us_junctions = upstream_junction_list.size();

		int this_junction;
		outlet_junction = starting_junction;
		outlet_node = ChannelNetwork.get_Node_of_Junction(starting_junction);

		for(int i = 0; i< n_us_junctions; i++)
		{
			this_junction = upstream_junction_list[i];
			IndexChannelVector.push_back( ChannelNetwork.generate_link_index_channel_from_junction(this_junction,FlowInfo) );
		}
	}
	else if (organization_switch == 1)
	{
		int nodes_in_channel;

		// initiate the main stem
		LSDIndexChannel main_stem = ChannelNetwork.generate_longest_index_channel_in_basin(starting_junction,
		    					FlowInfo, DistanceFromOutlet);
		nodes_in_channel = main_stem.get_n_nodes_in_channel();
		IndexChannelVector.push_back(main_stem);
		receiver_channel.push_back(0);
		node_on_receiver_channel.push_back(nodes_in_channel-1);

		// now get the tributaries for the main stem
		vector<int> tributary_junctions;
		vector<int> node_of_tributaries;
		ChannelNetwork.extract_tributary_junctions_to_main_stem(main_stem, FlowInfo, tributary_junctions, node_of_tributaries);

		// get the number of main stem tributaries
		int n_main_stem_tributaries = tributary_junctions.size();
		//cout << "number main_stem trib: " << n_main_stem_tributaries << endl;
		// now loop through the tributaries, getting the longest channel for each, and linking them to the main stem
		for (int trib = 0; trib<n_main_stem_tributaries; trib++)
		{
			//cout << "trib junction: " << tributary_junctions[trib] << " node on main stem: "
            //                  << node_of_tributaries[trib] << endl;
			LSDIndexChannel main_stem_tributary= ChannelNetwork.generate_longest_index_channel_in_basin(tributary_junctions[trib],
																							FlowInfo, DistanceFromOutlet);
			IndexChannelVector.push_back(main_stem_tributary);

			// in this case the receiver is the main stem (index 0) so push back zero
			receiver_channel.push_back(0);
			node_on_receiver_channel.push_back(node_of_tributaries[trib]);

		}
	}
	else
	{
		cout << "LINE 1xx LSDIndexChannelTree Invalid ortganization switch!" << endl;
		exit(EXIT_FAILURE);
	}

}


// this is the same create function as above except it includes a pruning method that rejects
// some of the subbasins based on various criteria set by an integer called pruning switch
// pruning_switch == 0  channels are only added if they exceed a threshold drainage area
// pruning_switch == 1  channels are only added if the ratio between them and the mainstem
//						exceeds a certain value (pruning_threshold)
// pruning_switch == 2	channels are only added if the ratio between them and the area of the
//						mainstem _at the junction_ exceeds a certain value
void LSDIndexChannelTree::create(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
							  	int starting_junction, int org_switch, LSDRaster& DistanceFromOutlet,
									int pruning_switch, double pruning_threshold)
{
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	DataResolution = FlowInfo.get_DataResolution();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	NoDataValue = ChannelNetwork.get_NoDataValue();
	double pixel_area = DataResolution*DataResolution;

  	organization_switch = org_switch;

	if (organization_switch == 0)
	{
		// get the upstream junctions
		upstream_junction_list = ChannelNetwork.get_upslope_junctions(starting_junction);

		// get the number of upstream junctions
		int n_us_junctions = upstream_junction_list.size();

		int this_junction;
		outlet_junction = starting_junction;
		outlet_node = ChannelNetwork.get_Node_of_Junction(starting_junction);

		for(int i = 0; i< n_us_junctions; i++)
		{
			this_junction = upstream_junction_list[i];
			IndexChannelVector.push_back( ChannelNetwork.generate_link_index_channel_from_junction(this_junction,FlowInfo) );
		}
	}
	else if (organization_switch == 1)
	{
		int nodes_in_channel;

		// initiate the main stem
		LSDIndexChannel main_stem = ChannelNetwork.generate_longest_index_channel_in_basin(starting_junction,
																							FlowInfo, DistanceFromOutlet);
		nodes_in_channel = main_stem.get_n_nodes_in_channel();
		IndexChannelVector.push_back(main_stem);
		receiver_channel.push_back(0);
		node_on_receiver_channel.push_back(nodes_in_channel-1);
		double main_stem_drainage_area = double(main_stem.get_contributing_pixels_at_outlet(FlowInfo))*pixel_area;
		double contributing_pixel_area;

		// now get the tributaries for the main stem
		vector<int> tributary_junctions;
		vector<int> node_of_tributaries;
		ChannelNetwork.extract_tributary_junctions_to_main_stem(main_stem, FlowInfo, tributary_junctions, node_of_tributaries);

		// get the number of main stem tributaries
		int n_main_stem_tributaries = tributary_junctions.size();
		//cout << "number main_stem trib: " << n_main_stem_tributaries << endl;
		// now loop through the tributaries, getting the longest channel for each, and linking them to the main stem
		for (int trib = 0; trib<n_main_stem_tributaries; trib++)
		{
			LSDIndexChannel main_stem_tributary= ChannelNetwork.generate_longest_index_channel_in_basin(tributary_junctions[trib],
																								FlowInfo, DistanceFromOutlet);
			contributing_pixel_area = pixel_area*double(main_stem_tributary.get_contributing_pixels_at_penultimate_node(FlowInfo));
			// enter pruning logic
			if (pruning_switch == 0)
			{
				if (contributing_pixel_area >pruning_threshold)
				{
					IndexChannelVector.push_back(main_stem_tributary);
					// in this case the receiver is the main stem (index 0) so push back zero
					receiver_channel.push_back(0);
					node_on_receiver_channel.push_back(node_of_tributaries[trib]);
				}
			}
			else if (pruning_switch == 1)
			{
				//cout << "tributary number: " << trib << " ratio: " << contributing_pixel_area/main_stem_drainage_area << endl;

				if (contributing_pixel_area/main_stem_drainage_area >pruning_threshold)
				{
					IndexChannelVector.push_back(main_stem_tributary);
					// in this case the receiver is the main stem (index 0) so push back zero
					receiver_channel.push_back(0);
					node_on_receiver_channel.push_back(node_of_tributaries[trib]);
				}
			}
			else if (pruning_switch == 2)
			{
				double main_stem_junc_area = double(main_stem.get_contributing_pixels_at_node(node_of_tributaries[trib], FlowInfo))*
										pixel_area;
				if ( contributing_pixel_area/main_stem_junc_area > pruning_threshold)
				{
					IndexChannelVector.push_back(main_stem_tributary);
					// in this case the receiver is the main stem (index 0) so push back zero
					receiver_channel.push_back(0);
					node_on_receiver_channel.push_back(node_of_tributaries[trib]);
				}
			}
			else
			{
				IndexChannelVector.push_back(main_stem_tributary);
				// in this case the receiver is the main stem (index 0) so push back zero
				receiver_channel.push_back(0);
				node_on_receiver_channel.push_back(node_of_tributaries[trib]);
			}

		}
	}
	else
	{
		cout << "LINE 1xx LSDIndexChannelTree Invalid ortganization switch!" << endl;
		exit(EXIT_FAILURE);
	}

}



// this function calcualtes the chi value starting from the bottom node of the channel tree and working its way up
// note that junctions are the top of the channel
vector< vector<double> > LSDIndexChannelTree::calculate_chi_from_channel_tree(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
                                                     double m_over_n, double A_0)
{
	if(organization_switch != 0)
	{
		cout << "LSDIndexChannelTree you can't run LSDIndexChannelTree::calculate_chi_from_channel_tree with this channel organization" << endl;
		exit(EXIT_FAILURE);
	}

	double root2 = 1.41421356;
	double diag_length = root2*DataResolution;
	double dx;
	double pixel_area = DataResolution*DataResolution;
	int curr_node;
	//cout << "diag length:" << diag_length << " and data resolution: " << DataResolution << endl;

	vector< vector<double> > chi_vectors;
	int n_channels = IndexChannelVector.size();

	// read through the tree
	cout << "the number of channels is: " << n_channels << endl;
	if (n_channels == 0)
	{
		cout << "LINE 70 LSDIndexChannelTree, empty channel network" << endl;
		exit(EXIT_FAILURE);
	}

	int n_nodes_in_link;
	//int current_node_index;

	// the junctions are at the top of the channel segment.
	// this will get the channel downslope of the starting junction
	// on the first channel you need to assign a chi value of 0
	n_nodes_in_link = IndexChannelVector[0].get_n_nodes_in_channel();

	// the bottom node is at chi of zero
	// initiate the chi vector
	vector<double> empty_vec;
	vector<double> chi_temp(n_nodes_in_link,0.0);
	// now loop up through the channel, adding chi values
	// note, the channel index are arranges with upstream element first, so you need to go through the channel
	// in reverse order
	//cout << "YOYOMA nodes in channel: " << n_nodes_in_link << endl;
	for (int ChIndex = n_nodes_in_link-2; ChIndex>=0; ChIndex--)
	{
	 	//cout << "ChIndex is: " << ChIndex << endl;
		curr_node = IndexChannelVector[0].get_node_in_channel(ChIndex);
		if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
		{
			dx = diag_length;
		}
		else
		{
			dx = DataResolution;
		}
		//cout << "dx is: " << dx << endl;

		chi_temp[ChIndex] = dx*(pow( (A_0/ (double(
			                    FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area) ),
			                    m_over_n))
		                       + chi_temp[ChIndex+1];

	}
	// add the chi_vector to the vector of vecotrs
	chi_vectors.push_back(chi_temp);

	int junction_of_this_link;
	int receiver_junction_of_this_link;
	int receiver_channel_index;

	// the chi iterator
  	vector< vector<double> >::iterator chi_iter;

	// loop through the rest of the channels, calcualting chi as you go
	for (int i = 1; i<n_channels; i++)
	{
		n_nodes_in_link = IndexChannelVector[i].get_n_nodes_in_channel();

		// reset chi_temp
		chi_temp = empty_vec;
		chi_temp.resize(n_nodes_in_link,-99);
		//cout << endl << "link number: " <<  " this link has " << n_nodes_in_link << " nodes" <<endl;

		// within this loop all channels are upstream of the source junction.
		// the end node of a channel (its downstream node) is the start node (its upstream node)
		// if its receiver channel. So you need to find the reciever vector. To do this you need to
		// know the index of the channel in the IndexChannelVector that is the downstream junction
		junction_of_this_link = upstream_junction_list[i];

		// get index in the upslope_junctions_list of the reciever junction
		// because of the fastscape ordering, this should already have a chi vector assigned
		receiver_junction_of_this_link = ChannelNetwork.get_Receiver_of_Junction(junction_of_this_link);
		receiver_channel_index = ChannelNetwork.map_junction_to_upslope_junction_list(upstream_junction_list, receiver_junction_of_this_link);

    	//cout << " the reciever_junc is " << receiver_junction_of_this_link
        // << " which is indexed at " << receiver_channel_index << endl;

		// now get the most uplsope chi of the reciever channel
		chi_iter = chi_vectors.begin();
		for (int j = 0; j<receiver_channel_index; j++)
		{
      		chi_iter++;
    	}
    	chi_temp[n_nodes_in_link-1] = (*chi_iter)[0];
    	cout << "the downstream chi is " <<   chi_temp[n_nodes_in_link-1] << endl;

    	// now loop upslope through the channel. The chi values increment
    	// from the channel pixel immediately downslope.
		for (int ChIndex = n_nodes_in_link-2; ChIndex>=0; ChIndex--)
		{
			curr_node = IndexChannelVector[i].get_node_in_channel(ChIndex);
			if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
			{
				dx = diag_length;
			}
			else
			{
				dx = DataResolution;
			}

			chi_temp[ChIndex] = dx*(pow( (A_0/ (double(
									FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area) ),
									m_over_n))
								   + chi_temp[ChIndex+1];
			//cout << "link " << i<<", node " << curr_node << " and chi: " << chi_temp[ChIndex] << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
		}

    	chi_vectors.push_back(chi_temp);

	}



  return chi_vectors;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function prints chi values. It is used on the channel tree when channels are organized by links
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannelTree::print_chi_vs_elevation_from_channel_tree(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
                                                     double m_over_n, double A_0, string chi_vs_elev_fname)
{

	if(organization_switch != 0)
	{
		cout << "LSDIndexChannelTree you can't run LSDIndexChannelTree::print_chi_vs_elevation_from_channel_tree with this channel organization" << endl;
		exit(EXIT_FAILURE);
	}

	int n_channels = IndexChannelVector.size();
	int n_nodes_in_link,current_node_index;
	vector< vector<double> > chi_vectors = calculate_chi_from_channel_tree(FlowInfo, ChannelNetwork,
                                                     m_over_n, A_0);

	// the chi iterator
  	vector< vector<double> >::iterator chi_iter;
    chi_iter =  chi_vectors.begin();
	double current_chi;
	double current_elev;
	int curr_row,curr_col;

	ofstream chi_elev_out;
	chi_elev_out.open(chi_vs_elev_fname.c_str());

	for (int i = 0; i<n_channels; i++)
	{
		n_nodes_in_link = IndexChannelVector[i].get_n_nodes_in_channel();

	  	for (int this_node = 0; this_node<n_nodes_in_link; this_node++)
		{
			current_node_index = IndexChannelVector[i].get_node_in_channel(this_node);
			current_chi  =  (*chi_iter)[this_node];
			FlowInfo.retrieve_current_row_and_col(current_node_index,curr_row,curr_col);
			current_elev = Elevation.get_data_element(curr_row,curr_col);

			//cout << "current_node: " << current_node_index << " and chi: " << current_chi << endl;
			//chi_elev_out << current_chi << " " << current_elev << endl;
		}
		chi_iter++;
	}
	chi_elev_out.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function minimizes the R^2 of the main stem channel assuming it is in steady state
// that is assuming the entire main stem is undergoing the same uplift
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double LSDIndexChannelTree::fit_m_over_n_mainstem(vector<double>& m_over_n_values, vector<double>& R_squared,
												double A_0, LSDFlowInfo& FlowInfo, LSDRaster& Elevation_Raster,
												double start_movn, double increment_movn, int n_movn)
{
	double downslope_chi = 0.0;

	// we loop through m_over_n values getting R_2 for each
	vector<double> R_squared_local;
	vector<double> m_over_n_values_local;
	double this_movn;
	double max_r2 = 0;
	double max_movn = 0;

	//cout << "LSDIndexChannelTree::fit_m_over_n_mainstem n_movn is: " << n_movn << endl;

	for(int i = 0; i<n_movn; i++)
	{
		this_movn = start_movn+ double(i)*increment_movn;
		m_over_n_values_local.push_back(this_movn);
		LSDChannel MainStem(downslope_chi, this_movn, A_0, IndexChannelVector[0], FlowInfo, Elevation_Raster);

		vector<double> chi = MainStem.get_Chi();
		vector<double> elev = MainStem.get_Elevation();
		vector<double> residuals;

		// do least squares regression
		vector<double> results = simple_linear_regression(chi, elev, residuals);

		//cout << "LSDIndexChannelTree::fit_m_over_n_mainstem m_ov_n: " << this_movn << " and R2: " << results[2] << endl;

		// replace max_r2 if it is the maximum r2 value
		if(results[2] > max_r2)
		{
			max_r2 = results[2];
			max_movn = this_movn;
		}
		R_squared_local.push_back(results[2]);
	}

	//cout << "The best fit m over n ratio is:" <<max_movn << endl;
	return max_movn;
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function works on organiation_switch == 1, where you have a main stem channel and the longest
// tributaries to this channel.
// It takes these Index channels and then turns them into channels with full drainage area,
// chi and elevation data
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<LSDChannel> LSDIndexChannelTree::retrieve_LSDChannels_from_tree(double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
{
	if (organization_switch != 1)
	{
		cout << "LSDIndexChannelTree you can't run LSDIndexChannelTree::retrieve_LSDChannels_from_tree" << endl;
		cout << "with this channel organization, organization switch: " << organization_switch << endl;
		exit(EXIT_FAILURE);
	}

	vector<LSDChannel> vector_of_channels;
	int receiver_chann;

	// get the main stem chi
	double downslope_chi = 0.0;
	LSDChannel MainStem(downslope_chi, m_over_n, A_0, IndexChannelVector[0], FlowInfo, Elevation_Raster);

	vector_of_channels.push_back(MainStem);

	// now loop through the tributaries
	int n_tributaries = IndexChannelVector.size();

	if (n_tributaries > 1)
	{
		for(int trib = 1; trib<n_tributaries; trib++)
		{
			// figure out what channel this drains to
			receiver_chann = receiver_channel[trib];

			// get the chi value at the tributary
			downslope_chi = vector_of_channels[receiver_chann].retrieve_chi_at_channel_node(node_on_receiver_channel[trib]);

			// make the new channel
			LSDChannel TempTrib(downslope_chi, m_over_n, A_0, IndexChannelVector[trib], FlowInfo, Elevation_Raster);

			// add it to the vector
			vector_of_channels.push_back(TempTrib);
		}
	}

	return vector_of_channels;
}


// this function uses the segment fitting tool to look for the best fit values of m over n
double LSDIndexChannelTree::search_for_best_fit_m_over_n(double A_0, int n_movern, double d_movern,double start_movern,
						       int minimum_segment_length, double sigma, int target_nodes,
                                  LSDFlowInfo& FlowInfo,  LSDRaster& Elevation_Raster, string fname)
{
  double m_over_n;
  int n_channels = IndexChannelVector.size();

  // data structures to hold information about segments from all channels
  // these are the best fit data for individual channels
  vector< vector<double> > b_vecvec(n_channels);
  vector< vector<double> > m_vecvec(n_channels);
  vector< vector<double> > DW_vecvec(n_channels);
  vector< vector<double> > r2_vecvec(n_channels);
  vector< vector<double> > thinned_chi_vecvec(n_channels);
  vector< vector<double> > thinned_elev_vecvec(n_channels);
  vector< vector<double> > fitted_elev_vecvec(n_channels);
  vector< vector<int> > node_ref_thinned_vecvec(n_channels);
  vector<vector<int> > these_segment_lengths_vecvec(n_channels);
  vector<double> MLE_vec(n_channels);
  vector<int> n_segments_vec(n_channels);
  vector<int> n_data_nodes_vec(n_channels);
  vector<double> AICc_vec(n_channels, 9999);
  vector<double> best_m_over_n(n_channels);

  // data structures to hold information about segments from all channels
  // these are the best fit data for the cumulative channels
  vector< vector<double> > cum_b_vecvec(n_channels);
  vector< vector<double> > cum_m_vecvec(n_channels);
  vector< vector<double> > cum_DW_vecvec(n_channels);
  vector< vector<double> > cum_r2_vecvec(n_channels);
  vector< vector<double> > cum_thinned_chi_vecvec(n_channels);
  vector< vector<double> > cum_thinned_elev_vecvec(n_channels);
  vector< vector<double> > cum_fitted_elev_vecvec(n_channels);
  vector< vector<int> > cum_node_ref_thinned_vecvec(n_channels);
  vector<vector<int> > cum_these_segment_lengths_vecvec(n_channels);

  // these data are for the cumulative AICs
  vector<double> AICc_combined_vec(n_movern);
  vector<double> m_over_n_vec(n_movern);

  // these are from the individual channels, which are replaced each time a new channel is analyzed
  vector<double> m_vec;
  vector<double> b_vec;
  vector<double> r2_vec;
  vector<double> DW_vec;
  vector<double> fitted_y;
  int n_data_nodes;
  int this_n_segments;
  double this_MLE, this_AIC, this_AICc;
  vector<int> these_segment_lengths;
  vector<double> chi_thinned;
  vector<double> elev_thinned;
  vector<double> elev_fitted;
  vector<int> node_ref_thinned;

  // some data about the AICc

  // loop through m_over_n values
  for(int movn = 0; movn< n_movern; movn++)
    {
      m_over_n = double(movn)*d_movern+start_movern;

      // get the vector of channels for this m_over_n
      vector<LSDChannel> vector_of_channels = retrieve_LSDChannels_from_tree(m_over_n, A_0, FlowInfo,
									 Elevation_Raster);

      vector<double> MLEs_thischan(n_channels);
      vector<int> n_segs_thischan(n_channels);
      vector<int> n_datanodes_thischan(n_channels);


      // now loop through channels
      for (int chan = 0; chan<n_channels; chan++)
	{
	  // get the channels for this m over n ratio
	  vector_of_channels[chan].find_most_likeley_segments(minimum_segment_length, sigma, target_nodes,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
							      this_AIC, this_AICc );
	  // check to see if the AICc value is the smallest
	  // if so add the data to the best fit data elements
	  if (this_AICc < AICc_vec[chan])
	    {
	       b_vecvec[chan] = b_vec;
	       m_vecvec[chan] = m_vec;
	       DW_vecvec[chan] = DW_vec;
	       r2_vecvec[chan] = r2_vec;
	       thinned_chi_vecvec[chan] = chi_thinned;
	       thinned_elev_vecvec[chan] = elev_thinned;
	       fitted_elev_vecvec[chan] = elev_fitted;
	       node_ref_thinned_vecvec[chan] = node_ref_thinned;
	       these_segment_lengths_vecvec[chan] = these_segment_lengths;
	       MLE_vec[chan] = this_MLE;
	       n_segments_vec[chan] = this_n_segments;
	       n_data_nodes_vec[chan] = n_data_nodes;
	       AICc_vec[chan] = this_AICc;
	       best_m_over_n[chan] = m_over_n;
	    }

	  // add the data from this channel to the vectors that will be used to calcualte cumulative AICc
	  MLEs_thischan[chan] = this_MLE;
	  n_segs_thischan[chan] = this_n_segments;
	  n_datanodes_thischan[chan] = n_data_nodes;

	}

      //now calculate the cumulative AICc for this m over n
      double thismn_AIC;
      double thismn_AICc;

      int n_total_segments = 0;
      int n_total_nodes = 0;
      double cumulative_MLE = 1;

      for (int chan = 0; chan<n_channels; chan++)
	{
	  //cout << "JUMBO m_ov_n node: " << movn << " chan: " << chan << " MLE: " << MLEs_thischan[chan] << endl;
	  n_total_segments += n_segs_thischan[chan];
	  n_total_nodes += n_datanodes_thischan[chan];
	  cumulative_MLE = MLEs_thischan[chan]*cumulative_MLE;
	}


      // these AIC and AICc values are cumualtive for a given m_over_n
      thismn_AIC = 4*n_total_segments-2*log(cumulative_MLE);		// the 4 comes from the fact that
			                                                       // for each segment there are 2 parameters
      thismn_AICc =  thismn_AIC + 2*n_total_segments*(n_total_segments+1)/(n_total_nodes-n_total_segments-1);
      AICc_combined_vec[movn] = thismn_AICc;
      m_over_n_vec[movn] = m_over_n;

      //cout << endl << endl <<"m_over_n: " << m_over_n << " and combined AICc: " << thismn_AICc << endl;
      //cout << "this cumulative MLE: " << cumulative_MLE << " n_segs: " << n_total_segments << " and n_nodes: " << n_total_nodes << endl;

    }



    //cout << "and the cumulative m_over n values"<< endl;
    double min_cum_AICc = 9999;
    double bf_cum_movn = start_movern;
    for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
      {
		//cout << "m over n: " << m_over_n_vec[mn] << " and AICc: " << AICc_combined_vec[mn] << endl;
		// if this is the minimum, store the m over n value
		if(AICc_combined_vec[mn] < min_cum_AICc)
	  	{
	   	 	min_cum_AICc = AICc_combined_vec[mn];
	    	bf_cum_movn = m_over_n_vec[mn];
	  	}
      }

    // now get the cumulative best fit channels
    vector<LSDChannel> cum_bf_vector_of_channels = retrieve_LSDChannels_from_tree(bf_cum_movn, A_0, FlowInfo,
									 Elevation_Raster);
    // now loop through channels
    for (int chan = 0; chan<n_channels; chan++)
      {
	// get the channels for this m over n ratio
	cum_bf_vector_of_channels[chan].find_most_likeley_segments(minimum_segment_length, sigma, target_nodes,
					    b_vec, m_vec, r2_vec,DW_vec,chi_thinned, elev_thinned,
					    elev_fitted, node_ref_thinned,these_segment_lengths,
					    this_MLE, this_n_segments, n_data_nodes,
                                            this_AIC, this_AICc );
	 cum_b_vecvec[chan] = b_vec;
	 cum_m_vecvec[chan] = m_vec;
         cum_DW_vecvec[chan] = DW_vec;
         cum_r2_vecvec[chan] = r2_vec;
         cum_thinned_chi_vecvec[chan] = chi_thinned;
         cum_thinned_elev_vecvec[chan] = elev_thinned;
         cum_fitted_elev_vecvec[chan] = elev_fitted;
         cum_node_ref_thinned_vecvec[chan] = node_ref_thinned;
         cum_these_segment_lengths_vecvec[chan] = these_segment_lengths;
      }


    // write a file
    ofstream best_fit_info;
	best_fit_info.open(fname.c_str());

  	best_fit_info << "N_channels: " << n_channels << endl;
  	best_fit_info << "m_over_n_for_channels ";
    for (int ch = 0; ch<n_channels; ch++)
    {
		best_fit_info << "  " << best_m_over_n[ch];
    }
  	best_fit_info << endl << "m_over_n_values ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << m_over_n_vec[mn];
	}
	best_fit_info << endl << "cumulative_AICc: ";
  	for (int mn = 0; mn< int(m_over_n_vec.size()); mn++)
    {
		best_fit_info << " " << AICc_combined_vec[mn];
	}
	best_fit_info << endl;
	cout << "The_best_fit_cumulative_m_over_n_is: " << bf_cum_movn << endl;
	for (int chan = 0; chan<n_channels; chan++)
	{
		vector<int> seglength = cum_these_segment_lengths_vecvec[chan];
		vector<double> m_val = cum_m_vecvec[chan];
		vector<double> b_val = cum_b_vecvec[chan];
		vector<double> DW_val = cum_DW_vecvec[chan];
		vector<double> r2_val = cum_r2_vecvec[chan];
		int n_segs_this_channel = seglength.size();

		best_fit_info << "Channel " << chan << " segment_length";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << seglength[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_gradient";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << m_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_intercept";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << b_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_DW_stat";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << DW_val[i];
		}
		best_fit_info << endl;
		best_fit_info << "Channel " << chan << " segment_r2";
		for (int i = 0; i<n_segs_this_channel; i++)
		{
			best_fit_info << " " << r2_val[i];
		}
		best_fit_info << endl;
	}


    for (int chan = 0; chan<n_channels; chan++)
    {
		chi_thinned = thinned_chi_vecvec[chan];
		elev_thinned = thinned_elev_vecvec[chan];
		elev_fitted = fitted_elev_vecvec[chan];

		// print the cumulative best fit profiles
        int thin_n_nodes = chi_thinned.size();

		for(int i = 0; i<thin_n_nodes; i++)
	  	{
	    	best_fit_info << chan << " " << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
	  	}
     }
     return bf_cum_movn;
}








//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function prints chi and elevation, along with flow distance and the number of the tributary
// it all goes to one file
//
// the file format is
// channel_number node_index row column flow_dist chi elevation drainage_area
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannelTree::print_LSDChannels_from_tree(double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster, LSDRaster& FlowDistance, string fname)
{
	if (organization_switch != 1)
	{
		cout << "LSDIndexChannelTree you can't run LSDIndexChannelTree::retrieve_LSDChannels_from_tree" << endl;
		cout << "with this channel organization, organization switch: " << organization_switch << endl;
		exit(EXIT_FAILURE);
	}

	// open the outfile
	ofstream channelfile_out;
	channelfile_out.open(fname.c_str());

	// get the vector of channels
	vector<LSDChannel> vector_of_channels = retrieve_LSDChannels_from_tree(m_over_n, A_0, FlowInfo,Elevation_Raster);

	int n_channels = vector_of_channels.size();
	int n_nodes_in_channel;
	int node,row,col;
	double elev,chi,drain_area,flow_dist;
	//loop through the channels
	for (int i = 0; i< n_channels; i++)
	{
		// get the number of nodes in the channel
		n_nodes_in_channel =IndexChannelVector[i].get_n_nodes_in_channel();

		// now loop through the channel, printing out the data.
		for(int ch_node= 0; ch_node<n_nodes_in_channel; ch_node++)
		{
			IndexChannelVector[i].get_node_row_col_in_channel(ch_node, node, row, col);
			vector_of_channels[i].retrieve_node_information(ch_node, elev, chi, drain_area);
			flow_dist = FlowDistance.get_data_element(row,col);

			// print data to file
			channelfile_out << i << " " << node << " " << row << " " << col << " " << flow_dist << " "
			                << chi << " " << elev << " " << drain_area << endl;
		}
	}

	channelfile_out.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function prints chi and elevation, along with flow distance and the number of the tributary
// it all goes to one file
//
// the file format is
// channel_number node_index row column flow_dist chi elevation drainage_area
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannelTree::print_LSDChannels_for_chi_network_ingestion(LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster, LSDRaster& FlowDistance, string fname)
{
	if (organization_switch != 1)
	{
		cout << "LSDIndexChannelTree you can't run LSDIndexChannelTree::retrieve_LSDChannels_from_tree" << endl;
		cout << "with this channel organization, organization switch: " << organization_switch << endl;
		exit(EXIT_FAILURE);
	}

	// open the outfile
	ofstream channelfile_out;
	channelfile_out.open(fname.c_str());

	double m_over_n = 0.5;
	double A_0 = 1;

	// get the vector of channels
	vector<LSDChannel> vector_of_channels = retrieve_LSDChannels_from_tree(m_over_n, A_0, FlowInfo,Elevation_Raster);

	int n_channels = vector_of_channels.size();
	int n_nodes_in_channel;
	int node,row,col;
	double elev,chi,drain_area,flow_dist;

	// first print out some data about the dem
	channelfile_out << get_NRows() << endl;
	channelfile_out << get_NCols() << endl;
	channelfile_out << get_XMinimum() << endl;
	channelfile_out << get_YMinimum() << endl;
	channelfile_out << get_DataResolution() << endl;
	channelfile_out << get_NoDataValue() << endl;

	//loop through the channels
	for (int i = 0; i< n_channels; i++)
	{
		// get the number of nodes in the channel
		n_nodes_in_channel =IndexChannelVector[i].get_n_nodes_in_channel();

		// now loop through the channel, printing out the data.
		for(int ch_node= 0; ch_node<n_nodes_in_channel; ch_node++)
		{
			IndexChannelVector[i].get_node_row_col_in_channel(ch_node, node, row, col);
			vector_of_channels[i].retrieve_node_information(ch_node, elev, chi, drain_area);
			flow_dist = FlowDistance.get_data_element(row,col);

			// print data to file
			channelfile_out << i << " " << receiver_channel[i] << " " << node_on_receiver_channel[i] << " "
			                << node << " " << row << " " << col << " " << flow_dist << " "
			                << " " << elev << " " << drain_area << endl;
		}
	}

	channelfile_out.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes the channel tree and prints it to an LSDIndexRaster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDIndexChannelTree::TributaryArray_to_LSDIndexRaster()
{
	if(organization_switch != 1)
	{
		cout << "LSDIndexChannelTree you can't run LSDIndexChannelTree::create_tributary_array with this channel organization" << endl;
		cout << "organization_switch: " << organization_switch << endl;
		exit(EXIT_FAILURE);
	}

	Array2D<int> Channel_array(NRows,NCols,NoDataValue);
	int n_channels = IndexChannelVector.size();
	int node,row,col;

	cout << "n channels: " << n_channels  << endl;

	for(int chan = 0; chan<n_channels; chan++)
	{
		cout << "channel number: " << chan << endl;
		int n_nodes_in_channel = IndexChannelVector[chan].get_n_nodes_in_channel();


		for(int i = 0; i<n_nodes_in_channel-1; i++)
		{
			IndexChannelVector[chan].get_node_row_col_in_channel(i, node, row, col);
			//cout << "row: " << RowSequence[i] << " col: " << ColSequence[i] << endl;
			Channel_array[row][col]= chan;
		}

		// the last node, which is the downstream junction, will default to the
		// receiver channel
		IndexChannelVector[chan].get_node_row_col_in_channel(n_nodes_in_channel-1, node, row, col);
		Channel_array[row][col]=receiver_channel[chan];
	}

	LSDIndexRaster Channel_loc(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Channel_array);
	return Channel_loc;
}

#endif
