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
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannel_H
#define LSDIndexChannel_H


class LSDIndexChannel
{
	public:
	LSDIndexChannel()		{ create(); }
	LSDIndexChannel(int StartNode, int EndNode, LSDFlowInfo& FlowInfo)
							{ create(StartNode, EndNode, FlowInfo); }
	LSDIndexChannel(int StartJunction, int StartNode,
	                int EndJunction, int EndNode, LSDFlowInfo& FlowInfo)
							{ create(StartJunction,StartNode,
							  EndJunction,EndNode, FlowInfo); }

	// get functions
	// these get data elements
	int get_StartJunction() const			{ return StartJunction; }
	int get_EndJunction() const				{ return EndJunction; }
	int get_StartNode() const				{ return StartNode; }
	int get_EndNode() const					{ return EndNode; }

	int get_NRows() const				{ return NRows; }
	int get_NCols() const				{ return NCols; }
	double get_XMinimum() const			{ return XMinimum; }
	double get_YMinimum() const			{ return YMinimum; }
	double get_DataResolution() const	{ return DataResolution; }
	int get_NoDataValue() const		{ return NoDataValue; }

	vector<int> get_RowSequence() const		{ return RowSequence; }
	vector<int> get_ColSequence() const		{ return ColSequence; }
	vector<int> get_NodeSequence() const	{ return NodeSequence; }

	// this tells how many nodes are in the channel
	int get_n_nodes_in_channel() const		{return int(NodeSequence.size()); }

	// this gets the node index at a given node in the index channel
	int get_node_in_channel(int n_node);

	// get the number of contributing pixels at a given node in the channel
	int get_contributing_pixels_at_node(int n_node, LSDFlowInfo& FlowInfo);

	// this gets the node, row, and column index
	void get_node_row_col_in_channel(int n_node, int& node, int& row, int& col);

	// this gets the contriubting pixels at the outlet of the channel
	int get_contributing_pixels_at_outlet(LSDFlowInfo& FlowInfo)
								{ return FlowInfo.retrieve_contributing_pixels_of_node(EndNode); }
	// this is like the previous function, but gets the pixels at the penultimate node
	// this is useful when calculating area of basins where the tributary junction is at the EndNode
	int get_contributing_pixels_at_penultimate_node(LSDFlowInfo& FlowInfo);

	// this prints the channel to an LSDIndexRaster in order to see where the channel is
	LSDIndexRaster print_index_channel_to_index_raster();


	protected:

	// data for georeferencing
	int NRows;			// number of rows
	int NCols;			// number of columns
	double XMinimum;
	double YMinimum;

	// metadata
	double DataResolution;
	int NoDataValue;


	// data for georeferencing
	int StartJunction;				// the starting junction (numbered withing LSDChannelNetwork object
	int StartNode;			// the node index of the starting Junction (numbered in the FlowInfo object
	int EndJunction;
	int EndNode;

	vector<int> RowSequence;
	vector<int> ColSequence;
	vector<int> NodeSequence;

	private:
	void create();
	void create(int StartNode, int EndNode, LSDFlowInfo& FlowInfo);
	void create(int StartJunction, int StartNode,
	            int EndJunction, int EndNode, LSDFlowInfo& FlowInfo);


};



#endif
