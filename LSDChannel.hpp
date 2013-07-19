//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChannelNetwork.hpp
// header file for the LSDChannelNetwork object
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
// Version 0.0.1		30/08/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexChannel.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDChannel_H
#define LSDChannel_H


class LSDChannel: public LSDIndexChannel
{
	public:
	LSDChannel(LSDIndexChannel& InChann)
							{  create_LSDC(InChann); }
	LSDChannel(int StartNode, int EndNode, LSDFlowInfo& FlowInfo)
							{ create_LSDC(StartNode, EndNode, FlowInfo); }
	LSDChannel(int StartJunction, int StartNode,
	                int EndJunction, int EndNode, LSDFlowInfo& FlowInfo)
							{ create_LSDC(StartJunction,StartNode,
							  EndJunction,EndNode, FlowInfo); }
	LSDChannel(int StartNode, int EndNode, double downslope_chi,
                             double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
							{ create_LSDC(StartNode, EndNode, downslope_chi,
                              m_over_n, A_0, FlowInfo, Elevation_Raster); }
	LSDChannel(double downslope_chi, double m_over_n, double A_0,
							LSDIndexChannel& InChann, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
							{ create_LSDC(downslope_chi, m_over_n, A_0,
								InChann, FlowInfo, Elevation_Raster); }


	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=INHERITED=-=-=-=-=-
	// the following are incherited from LSDIndexChannel
	//  // get functions
	//  // these get data elements
	//  int get_StartJunction() const			{ return StartJunction; }
	//  int get_EndJunction() const				{ return EndJunction; }
	//  int get_StartNode() const				{ return StartNode; }
	//  int get_EndNode() const					{ return EndNode; }
	//
	//  int get_NRows() const				{ return NRows; }
	//  int get_NCols() const				{ return NCols; }
	//  double get_XMinimum() const			{ return XMinimum; }
	//  double get_YMinimum() const			{ return YMinimum; }
	//  double get_DataResolution() const	{ return DataResolution; }
	//  int get_NoDataValue() const		{ return NoDataValue; }
	//
	//  vector<int> get_RowSequence() const		{ return RowSequence; }
	//  vector<int> get_ColSequence() const		{ return ColSequence; }
	//  vector<int> get_NodeSequence() const	{ return NodeSequence; }
	//
	//  int get_n_nodes_in_channel() const		{return int(NodeSequence.size()); }
	//
	//  int get_node_in_channel(int n_node);
	//  void get_node_row_col_in_channel(int n_node, int& node, int& row, int& col);
	//
	//  // this prints the channel to an LSDIndexRaster in order to see where the channel is
	//  LSDIndexRaster print_index_channel_to_index_raster();
	//=-=-=-=-=-=-=INHERITED=-=-=-=-=-
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	// access the data
	// these are primarily used for fitting of the profiles
	vector<double> get_Chi()		{ return Chi; }
	vector<double> get_Elevation()	{ return Elevation; }

	// this calcualtes chi
	void calculate_chi(double downslope_chi, double m_over_n, double A_0, LSDFlowInfo& FlowInfo );

	// fnctions for getting data out of the channel
	double retrieve_chi_at_channel_node(int ch_node) { return Chi[ch_node]; }
	void retrieve_node_information(int ch_node, double& elev, double& chi, double& drainarea)
					{ elev = Elevation[ch_node]; chi = Chi[ch_node]; drainarea = DrainageArea[ch_node]; }

	// this function looks for the most likeley segments
        void find_most_likeley_segments(int minimum_segment_length, double sigma, int target_nodes,
                                             vector<double>& b_vec, vector<double>&  m_vec,
					     vector<double>& 	r2_vec,vector<double>&  DW_vec, 
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
                                             vector<double>& fitted_elev, vector<int>& node_ref_thinned,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes, 
					     double& this_AIC, double& this_AICc );
        // this function looks for the best fit of a channel for a range of m_over_n values where the channel has segments
        void find_best_fit_m_over_n_with_segments(int n_movern, double d_movern,double start_movern, 
						  double downslope_chi, double A_0, LSDFlowInfo& FlowInfo, 
					     int minimum_segment_length, double sigma, double target_nodes );

        // this function loops through m_over_n looking for best fit segments
        void find_best_fit_m_over_n_with_segments();


	protected:

	// This is an inherited class so
	// NOTE that there are DATA ELEMENTS INHERITED FROM LSDIndexChannel

	vector<double> Elevation;
	vector<double> Chi;
	vector<double> DrainageArea;


	private:
	void create_LSDC(LSDIndexChannel& IndexChannel);
	void create_LSDC(int StartNode, int EndNode, LSDFlowInfo& FlowInfo);
	void create_LSDC(int StartJunction, int StartNode,
	            int EndJunction, int EndNode, LSDFlowInfo& FlowInfo);
	void create_LSDC(int SJN, int EJN, double downslope_chi,
                             double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster);
	void create_LSDC(double downslope_chi, double m_over_n, double A_0,
						LSDIndexChannel& InChann, LSDFlowInfo& FlowInfo,
                        LSDRaster& Elevation_Raster);
};


#endif
