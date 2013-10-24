//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// calculate_channel_heads.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program calculates channel heads using a chi method described in
// Clubb et al. (manuscript in prep)
//
// Developed by:
//  Fiona Clubb
//  Simon M. Mudd
//
// Copyright (C) 2013 Fiona Clubb and Simon M. Mudd 2013
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb, Univertsity of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Modifications:
// Simon M. Mudd, 25/9/2013: Attempted to speed up code
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "LSDStatsTools.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDChannelNetwork.hpp"
#include "LSDIndexChannelTree.hpp"
#include "LSDChiNetwork.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
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

	string DEM_name;
	string fill_ext = "_fill";
	file_info_in >> DEM_name;
	//int junction_number;
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

	file_info_in >> Minimum_Slope >> threshold
				>> pruning_threshold >> A_0 >> minimum_segment_length >> sigma >> start_movern
				>> d_movern >> n_movern >> target_nodes >> n_iterations >> fraction_dchi_for_variation
				>> vertical_interval >> horizontal_interval >> area_thin_frac >> target_skip;


	cout << "Paramters of this run: " << endl
	     //<< "junction number: " << junction_number << endl
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
	     << "vertical interval: " << vertical_interval << endl
	     << "horizontal interval: " << horizontal_interval << endl
	     << "area thinning fraction for SA analysis: " << area_thin_frac << endl
	     << "target_skip is: " << target_skip << endl;


	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);

	// get the filled file
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
	filled_topo_test.write_raster((DEM_f_name),DEM_flt_extension);

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);


	// calculate the distance from outlet and contributing pixels
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  	string dist_raster_fname = "_FD";
  	DistanceFromOutlet.write_raster((path_name+DEM_name+dist_raster_fname),DEM_flt_extension);

  	string NI_fname = "_NI";
	//LSDIndexRaster NodeIndexRaster =  FlowInfo.write_NodeIndex_to_LSDIndexRaster();
  //	NodeIndexRaster.write_raster((path_name+DEM_name+NI_fname),DEM_flt_extension);



	// get the sources: note: this is only to select basins!
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);

	// now get the junction network
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);

	// get the stream orders and the junctions
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();

	string SO_name = "_SO";
	string JI_name = "_JI";

	SOArray.write_raster((path_name+DEM_name+SO_name),DEM_flt_extension);
	JIArray.write_raster((path_name+DEM_name+JI_name),DEM_flt_extension);


	//int NRows = FlowInfo.get_NRows();
	//int NCols = FlowInfo.get_NCols();
	//double NoDataValue = FlowInfo.get_NoDataValue();
	//double XMinimum = FlowInfo.get_XMinimum();
	//double YMinimum = FlowInfo.get_YMinimum();
	//double DataResolution = FlowInfo.get_DataResolution();

  	int basin_order = 3;
  	double m_over_n = 0.3;
  	int MinSegLength = 10;

  // get the channel heads
	vector<int> Channel_heads = ChanNetwork.GetChannelHeadsChiMethodBasinOrder(basin_order,
	                            MinSegLength, A_0, m_over_n,FlowInfo, DistanceFromOutlet,
									            filled_topo_test);

	//print channel heads to a raster
	string CH_name = "_CH";
	LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(Channel_heads);
	Channel_heads_raster.write_raster((path_name+DEM_name+CH_name),DEM_flt_extension);
	
	//create a channel network based on these channel heads
	LSDChannelNetwork NewChanNetwork(Channel_heads, FlowInfo);
	LSDIndexRaster SOArrayNew = NewChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name_new = "_SO_from_CH";
	
	SOArrayNew.write_raster((path_name+DEM_name+SO_name_new),DEM_flt_extension);

}
