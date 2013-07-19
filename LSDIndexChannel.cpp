// source code for the LSDChannelIndex object
// this object contains the node indeces as well as the
// row and col indices for individual channel segments
// these could be arranged arbitrailiy accoridng to channel
// junctions or simply nodes downstream of a given node and upstream
// of another arbitrary node EndNode

// there should also be an LSDChannel object which contains
// the actual data of the channel including things like
// elevation, distance downstream, area, etc.
// however this needs to be thrugh out carefully since it
// seems different applications will need different kinds
// and amounts of data in the LSDChannel object. One perhaps
// could use a list of vectors, each vector containing the
// double data along the channel but this list is expandable
// so one could have channels with just elevation and distance data
// but could also have more detailed channels.


#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexChannel.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannel_CPP
#define LSDIndexChannel_CPP


void LSDIndexChannel::create()
{
	//cout << "LSDIndexChannel You need to initialize with some parameters" << endl;
	//exit(EXIT_FAILURE);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// first create routine
// creates an index channel with just the node index of the starting and ending nodes
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(int SJN, int EJN, LSDFlowInfo& FlowInfo)
{

	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();

	StartJunction = -1;
	EndJunction = -1;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;


		}
		else
		{
			curr_node = receive_node;
		}
	}

	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// second create routine
// creates an index channel with just the node index of the starting and ending nodes
// also includes junction information
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(int SJ, int SJN, int EJ, int EJN, LSDFlowInfo& FlowInfo)
{
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();

	StartJunction = SJ;
	EndJunction = EJ;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		//cout << "receive_node: " << receive_node << " and Endnode: " << EndNode << endl;

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;

		}
		else
		{
			curr_node = receive_node;
		}
	}
	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the n contributing pixels at one before the final pixel.
// useful for when you want to get the basin area just above the tributary
// junction, which is usually at EndNode
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_contributing_pixels_at_penultimate_node(LSDFlowInfo& FlowInfo)
{
	if(StartNode==EndNode)
	{
		return FlowInfo.retrieve_contributing_pixels_of_node(EndNode);
	}
	else
	{
		int n_nodes = RowSequence.size();
		return FlowInfo.retrieve_contributing_pixels_of_node( NodeSequence[n_nodes-2]);
	}
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the n contributing pixels at the node in the channel (not the node index)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_contributing_pixels_at_node(int n_node, LSDFlowInfo& FlowInfo)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	if (n_node < 0)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	return FlowInfo.retrieve_contributing_pixels_of_node( NodeSequence[n_node] );
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the node index (from flowInfo object
// of a node in the channel
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_node_in_channel(int n_node)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	if (n_node < 0)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	return NodeSequence[n_node];
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the node index (from flowInfo object) and row and column indices
// into an LSDRaster or IndexRaster
// of a node in the channel
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_node_row_col_in_channel(int n_node, int& node, int& row, int& col)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	if (n_node < 0)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	node = NodeSequence[n_node];
	row = RowSequence[n_node];
	col = ColSequence[n_node];
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the channel onto an index raster. Used to test where the channel is.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDIndexChannel::print_index_channel_to_index_raster()
{
	int n_nodes_in_channel = int(NodeSequence.size());
	cout << "NRows: " << NRows << " NCols: " << NCols << endl;

	Array2D<int> Channel_array(NRows,NCols,NoDataValue);
	for(int i = 0; i<n_nodes_in_channel; i++)
	{
		//cout << "row: " << RowSequence[i] << " col: " << ColSequence[i] << endl;
		Channel_array[RowSequence[i]][ColSequence[i]]= 1;
	}

	LSDIndexRaster Channel_loc(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Channel_array);
	return Channel_loc;
}



#endif
