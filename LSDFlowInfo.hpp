//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo.hpp
// header file for the LSDFlowInfo object
// LSD stands for Land Surface Dynamics
// This is a data object which generates and then
// stores information about flow routing
// It is the object that is used to generate
// contributing area, etc
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
// Version 0.0.1		29/08/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// TODO: add a create function that lets the user load a LSDFlowInfo object from files
//			// then you need a binary 'pickle' function that prints all the binary information
// 			// into one file.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDFlowInfo_H
#define LSDFlowInfo_H


class LSDFlowInfo
{
	public:

	// declare the channel network object to be a friend
	friend class LSDChannelNetwork;

	LSDFlowInfo()										{ create(); }
	LSDFlowInfo(vector<string> b_conditions, LSDRaster& TopoRaster)
									{ create(b_conditions, TopoRaster); }
	void add_to_stack(int lm_index, int& j_index, int bl_node);

	// some functions that print out indeces to rasters
	LSDIndexRaster write_NodeIndex_to_LSDIndexRaster();
	LSDIndexRaster write_FlowDirection_to_LSDIndexRaster();
	LSDIndexRaster write_FlowLengthCode_to_LSDIndexRaster();
	LSDIndexRaster write_FlowDriection_to_LSDIndexRaster_Arcformat();

	// this prints the flow information to file
	void print_flow_info_vectors(string filename);

	// functions for getting flow, discharge, sediment flux, etc
	LSDIndexRaster calculate_n_pixels_contributing_from_upslope();
	LSDIndexRaster calculate_n_pixels_contributing_from_upslope_include_self();

	// this calculates area and makes an index into the s vector for efficient
	// calculation of the basin upslope of a given node.
	void calculate_upslope_reference_indices();
	LSDIndexRaster calculate_n_pixels_contributing_from_upslope_include_self_from_indices();

	// algorithms for basin collection
	vector<int> get_upslope_nodes(int node_number_outlet);

	// flow length functions
	LSDRaster distance_from_outlet();

	// functions for getting sources
	vector<int> get_sources_index_threshold(LSDIndexRaster& FlowPixels, int threshold);

	protected:

	// data for georeferencing
	int NRows;			// number of rows
	int NCols;			// number of columns
	int n_data_nodes;
	double XMinimum;
	double YMinimum;

	// metadata
	double DataResolution;
	int NoDataValue;

	// data about the flow network
	int NDataNodes;

	// the indexing data
	Array2D<int> NodeIndex;				// an array that says what node number
										// is at a given row and column
	Array2D<int> FlowDirection;			// a raster of flow direction information
										// in this format:
										// 7  0 1
										// 6 -1 2
										// 5  4 3
										// nodes with flow direction of -1 drain
										// to themselvs and are base level/sink nodes
	Array2D<int> FlowLengthCode;		// a code to denote the flow lengthfrom the node
										// to its reciever node (note each node
										// has one and only one receiver)
										// 0 == no receiver/self receiver (base level)
										// 1 == cardinal direction, flow length = DataResolution
										// 2 == diagonal, flow length = DataResolution*(1/sqrt(2));
	vector<int> RowIndex;				// This stores the row of a node in the vectorized
										// node index. It, combined with ColIndex, is the
										// inverse of NodeIndex
	vector<int> ColIndex;				// This stores the column of a node in the vectorized
										// node index. It, combined with RowIndex, is the
										// inverse of NodeIndex
	vector<int> BaseLevelNodeList;		// a list of base level nodes
	vector<int> NDonorsVector;			// stores the number of donors to each node
	vector<int> ReceiverVector;			// stores the node index of the receiving node
	vector<int> DeltaVector;			// stores the delta vector which proivides is used
										// to index into the donor stack and order contributing nodes
										// see Braun and Willett
	vector<int> DonorStackVector;		// This is a vector that stores the donor nodes of of the nodes
										// and is indexed by the DeltaVector
	vector<int> SVector;				// this vector is used to caluculate flow accumulation; for each
										// base level node it progresses from a hilltop to a confluence and then
										// jumps to the next hilltop so that by cascading down through
										// the node indices in this list one can quickly calcualte
										// drainage area, discharge, sediment flux, etc.
	vector<int> BLBasinVector;			// this stores the base level node for all of the nodes in the
										// DEM

	vector<int> SVectorIndex;			// this points to the starting point in the S vector of each node
	vector<int> NContributingNodes;		// the number of contributing nodes !!INCULDING SELF!! to a current
										// pixel. It is used in conjunction with the SVectorIndex to build
										// basins upslope of any and all nodes in the node list
	vector<string> BoundaryConditions;	// stores the boundary conditions in a vector of four strings
										// the conditions are N[0] E[1] S[2] W[3]

	private:
	void create();
	void create(vector<string> temp_BoundaryConditions, LSDRaster& TopoRaster);
};

#endif
