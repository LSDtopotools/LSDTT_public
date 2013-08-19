/** @file LSDChannel.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 0.0.1
@brief This object inherets from LSDIndexChannel and is used for chi analysis. 

@date 30/08/2012
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexChannel.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDChannel_H
#define LSDChannel_H

///@brief This object inherets from LSDIndexChannel and is used for chi analysis. 
class LSDChannel: public LSDIndexChannel
{
	public:
	
  /// @brief Creates an LSDChannel by copying from an IndexChannel.
  ///
  /// @details The starting node is upstream and the ending node is downstream.
  /// In this create function the junction indices are left blank (this can
  /// describe a channel between two arbitraty points.
  /// @param InChann LSDIndexChannel object.
	LSDChannel(LSDIndexChannel& InChann)
							{  create_LSDC(InChann); }

  /// @brief Creates an index channel with just the node index of the starting and ending nodes.
  /// @details The starting node is upstream and the ending node is downstream.
  /// In this create function the junction indices are left blank (this can
  /// describe a channel between two arbitraty points.
  /// @param StartNode Starting node.
  /// @param EndNode Ending node.
  /// @param FlowInfo LSDFlowInfo object. 							
	LSDChannel(int StartNode, int EndNode, LSDFlowInfo& FlowInfo)
							{ create_LSDC(StartNode, EndNode, FlowInfo); }
	
  /// @brief Creates an index channel with just the node index of the starting and ending nodes also includes junction information.
  /// @details The starting node is upstream and the ending node is downstream.
  /// In this create function the junction indices are left blank (this can
  /// describe a channel between two arbitraty points.
  /// @param StartJunction Starting junction.
  /// @param StartNode Starting node.
  /// @param EndJunction Ending junction.
  /// @param EndNode Ending node.
  /// @param FlowInfo LSDFlowInfo object. 
  LSDChannel(int StartJunction, int StartNode,
	                int EndJunction, int EndNode, LSDFlowInfo& FlowInfo)
							{ create_LSDC(StartJunction,StartNode,
							  EndJunction,EndNode, FlowInfo); }
	
	
  /// @brief This calculates all the channel areas, elevations and chi parameters based on for a starting node index and ending node index.
	/// @param StartNode Starting node.
  /// @param EndNode Ending node.
  /// @param downslope_chi Downslope Chi value.
	/// @param m_over_n m over n ratio.
	/// @param A_0 A_0 value.		
	/// @param FlowInfo LSDFlowInfo object. 
	/// @param Elevation_Raster Elevation LSDRaster object.
  LSDChannel(int StartNode, int EndNode, double downslope_chi,
                             double m_over_n, double A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
							{ create_LSDC(StartNode, EndNode, downslope_chi,
                              m_over_n, A_0, FlowInfo, Elevation_Raster); }
	
	
	
	
  /// @brief This calculates all the channel areas, elevations and chi parameters based on for a given LSDChannelIndex.
	/// @param downslope_chi Downslope Chi value.
	/// @param m_over_n m over n ratio.
	/// @param A_0 A_0 value.
  /// @param InChann LSDIndexChannel object.	
	/// @param FlowInfo LSDFlowInfo object. 
	/// @param Elevation_Raster Elevation LSDRaster object.
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
	
  /// @return Vector of chi values.
  vector<double> get_Chi()		{ return Chi; }
  /// @return Vector of elevation values.
	vector<double> get_Elevation()	{ return Elevation; }

	/// @brief This function uses a flow info object to calculate the chi values in the channel.
	/// @param downslope_chi Downslope Chi value.
	/// @param m_over_n m over n ratio.
	/// @param A_0 A_0 value.	
	/// @param FlowInfo LSDFlowInfo object. 
	void calculate_chi(double downslope_chi, double m_over_n, double A_0, LSDFlowInfo& FlowInfo );

	/// @brief Get chi value at channel node.
	/// @param ch_node Integer node index.
	/// @return chi value at channel node.
  double retrieve_chi_at_channel_node(int ch_node) { return Chi[ch_node]; }
	
	/// @brief Get node chi value, elevation and drainage area.
	/// @param ch_node Integer node index.
	/// @param elev Elevation data.
	/// @param chi Chi Value
	/// @param drainarea Drainage area.	
	void retrieve_node_information(int ch_node, double& elev, double& chi, double& drainarea)
					{ elev = Elevation[ch_node]; chi = Chi[ch_node]; drainarea = DrainageArea[ch_node]; }

	/// @brief This function looks for the most likeley segments.
	/// @param minimum_segment_length
  /// @param sigma Sigma value.
  /// @param target_nodes
  /// @param b_vec Vector of b values.
  /// @param m_vec Vector of m values.
  /// @param r2_vec Vector of r-squared values.
  /// @param DW_vec Vector of Durbin-Watson values.
  /// @param thinned_chi 
  /// @param thinned_elev
  /// @param fitted_elev
  /// @param node_ref_thinned
  /// @param these_segment_lengths
  /// @param this_MLE
  /// @param this_n_segments
  /// @param n_data_nodes
  /// @param this_AIC
  /// @param this_AICc
  void find_most_likeley_segments(int minimum_segment_length, double sigma, int target_nodes,
                                  vector<double>& b_vec, vector<double>&  m_vec,
					     vector<double>& 	r2_vec,vector<double>&  DW_vec, 
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
                                             vector<double>& fitted_elev, vector<int>& node_ref_thinned,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes, 
					     double& this_AIC, double& this_AICc );
  
  /// @brief This function looks for the best fit of a channel for a range of m_over_n values where the channel has segments.
  /// @param n_movern
  /// @param d_movern
  /// @param start_movern
	/// @param downslope_chi Downslope Chi value.
	/// @param A_0 A_0 value.	
	/// @param FlowInfo LSDFlowInfo object. 
  /// @param minimum_segment_length
  /// @param sigma Sigma value.
  /// @param target_nodes
  void find_best_fit_m_over_n_with_segments(int n_movern, double d_movern,double start_movern, 
						  double downslope_chi, double A_0, LSDFlowInfo& FlowInfo, 
					     int minimum_segment_length, double sigma, double target_nodes );

  /// @brief This function loops through m_over_n looking for best fit segments.
  void find_best_fit_m_over_n_with_segments();
  
	protected:

	// This is an inherited class so
	// NOTE that there are DATA ELEMENTS INHERITED FROM LSDIndexChannel

  /// @brief Elevation vector.
	vector<double> Elevation;
	/// @brief Chi vector.
  vector<double> Chi;
	/// @brief Drainage area vector.
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
