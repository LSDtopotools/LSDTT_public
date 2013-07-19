//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexChannelTree.hpp
// header file for the LSDIndexChannelTree object
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

#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDChannel.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDChannelNetwork.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannelTree_H
#define LSDIndexChannelTree_H


class LSDIndexChannelTree
{
	public:
	LSDIndexChannelTree(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork, int starting_junction)
							{ create(FlowInfo, ChannelNetwork, starting_junction); }
	LSDIndexChannelTree(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
	                    int starting_junction, int org_switch, LSDRaster& DistanceFromOutlet)
					       	{ create(FlowInfo, ChannelNetwork, starting_junction, org_switch, DistanceFromOutlet); }
	LSDIndexChannelTree(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
	                    int starting_junction, int org_switch, LSDRaster& DistanceFromOutlet,
	                    int pruning_switch, double pruning_threshold)
					       	{ create(FlowInfo, ChannelNetwork, starting_junction, org_switch, DistanceFromOutlet,
							         pruning_switch, pruning_threshold); }

	// acess function
	vector< LSDIndexChannel > get_LSDIndexChannelVector()	{ return IndexChannelVector; }


	// this calcualtes chi values along the channel
	vector< vector<double> > calculate_chi_from_channel_tree(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
											double m_over_n, double A_0);

	void print_chi_vs_elevation_from_channel_tree(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
                                                     double m_over_n, double A_0, string chi_vs_elev_fname);

	// this calcualtes the best fit m over n on the main stem channel.
	double fit_m_over_n_mainstem(vector<double>& m_over_n_values, vector<double>& R_squared,
				    	double A_0, LSDFlowInfo& FlowInfo, LSDRaster& Elevation_Raster,
				      	double start_movn, double increment_movn, int n_movn);

	// this returns a tributary array
	LSDIndexRaster TributaryArray_to_LSDIndexRaster();

	// this creates a vector of LSDChannels, they contain area and chi information
	vector<LSDChannel> retrieve_LSDChannels_from_tree(double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster);

    double search_for_best_fit_m_over_n(double A_0, int n_movern, double d_movern,double start_movern,
                                          int minimum_segment_length, double sigma, int target_nodes,
					                      LSDFlowInfo& FlowInfo,  LSDRaster& Elevation_Raster,  string fname);

	// this prints a file that contiains all the channel information. It can be used to plot
	// and analyze the channel profiles.
	void print_LSDChannels_from_tree(double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster, LSDRaster& FlowDistance, string fname);

	// this prints all the channels for ingestion into the chi analysis object
	// data extracted form this file can be used in a standalone chi analysis function
	void print_LSDChannels_for_chi_network_ingestion(LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster, LSDRaster& FlowDistance, string fname);


	// get functions
	// these get data elements
	int get_NRows() const				{ return NRows; }
	int get_NCols() const				{ return NCols; }
	double get_XMinimum() const			{ return XMinimum; }
	double get_YMinimum() const			{ return YMinimum; }
	double get_DataResolution() const	{ return DataResolution; }
	int get_NoDataValue() const			{ return NoDataValue; }

	protected:

	// data for georeferencing
	int NRows;			// number of rows
	int NCols;			// number of columns
	double XMinimum;
	double YMinimum;

	// metadata
	double DataResolution;
	int NoDataValue;

	int outlet_junction;
	int outlet_node;
	int organization_switch;						// there are a number of ways to organize this
													// data and this switch tells the object
													// how its data are organized.
													// it will reject member function operations if the data type is incorrect
	vector<int> upstream_junction_list;
	vector< LSDIndexChannel > IndexChannelVector;	// a vector containing all the index channel nodes

	vector<int> receiver_channel;
	vector<int> node_on_receiver_channel;

	private:
	void create(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork, int starting_junction);
	void create(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
									int starting_junction, int org_switch, LSDRaster& DistanceFromOutlet);
	void create(LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChannelNetwork,
									int starting_junction, int org_switch, LSDRaster& DistanceFromOutlet,
									int pruning_switch, double pruning_threshold);
};

#endif
