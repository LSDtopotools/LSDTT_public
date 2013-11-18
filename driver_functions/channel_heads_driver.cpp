//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// channel_heads_part2_driver.cpp
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
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include "LSDStatsTools.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDChannelNetwork.hpp"
#include "LSDIndexChannelTree.hpp"
#include "LSDChiNetwork.hpp"
#include "fftw-3.3.1/api/fftw3.h"
#include "TNT/jama_lu.h"
#include "TNT/tnt.h"

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
	int threshold;
	double Minimum_Slope;
	double A_0;
	double m_over_n;
	int no_connecting_nodes;

	file_info_in >> Minimum_Slope >> threshold >> A_0 >> m_over_n >> no_connecting_nodes;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
	
  //Filter the DEM using optimal wiener filtering
  LSDRasterSpectral SpectralRaster(topo_test);                                    
	LSDRaster topo_data_filtered = SpectralRaster.fftw2D_wiener();
	// Write raster of filtered DEM
	//topo_data_filtered.write_raster(DEM_outname,DEM_flt_extension);	

	// Set the no flux boundary conditions
  vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
		// get the filled file
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
	filled_topo_test.write_raster((DEM_f_name),DEM_flt_extension);
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	
	//string NI_name = "_NI";
  //LSDIndexRaster NodeIndex = FlowInfo.write_NodeIndex_to_LSDIndexRaster();
	//NodeIndex.write_raster((path_name+DEM_name+NI_name), DEM_flt_extension);
	
	//get the sources: note: this is only to select basins!
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold); 

	// now get the junction network
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);

	//string tan_curvature_name = "ind_curv";
	//string chan_heads_name = "ind_chan_heads";
	
	// Get the valleys using the contour curvature
	
	// Calculate polyfit coefficients and the tangential curvature
  Array2D<double> a;
  Array2D<double> b;
  Array2D<double> c;
  Array2D<double> d;
  Array2D<double> e;
  Array2D<double> f;
  //int NRows = topo_test.get_NRows();
  //int NCols = topo_test.get_NCols();
  //double XMinimum = topo_test.get_XMinimum();
  //double YMinimum = topo_test.get_YMinimum();
  //double DataResolution = topo_test.get_DataResolution();
  //double NoDataValue = topo_test.get_NoDataValue();
  double window_radius=6;
  
  topo_data_filtered.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster tan_curvature = topo_data_filtered.calculate_polyfit_tangential_curvature(a ,b ,c ,d, e);
  
  // Find the valley junctions
  Array2D<double> tan_curv_array = tan_curvature.get_RasterData();
  cout << "got tan curvature array" << endl;
	Array2D<int> valley_junctions = ChanNetwork.find_valleys(FlowInfo, tan_curv_array, sources, no_connecting_nodes);
	
	// Write the valley junctions to an LSDIndexRaster
  //string VJ_name = "_VJ";
  //LSDIndexRaster ValleyJunctionsRaster (NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,valley_junctions);
	//ValleyJunctionsRaster.write_raster((path_name+DEM_name+VJ_name),DEM_flt_extension);
	
	cout << "Got valley junctions, proceeding to chi analysis" << endl;
	
	// Calculate the channel head nodes
  int MinSegLength = 10;
  vector<int> ChannelHeadNodes = ChanNetwork.GetChannelHeadsChiMethodFromValleys(valley_junctions, MinSegLength,
                                                 A_0, m_over_n, FlowInfo, DistanceFromOutlet, filled_topo_test);
									                    
	//write channel heads to a raster
	string CH_name = "_CH";
	LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(ChannelHeadNodes);
	Channel_heads_raster.write_raster((path_name+DEM_name+CH_name),DEM_flt_extension);
	
	//create a channel network based on these channel heads
	LSDChannelNetwork NewChanNetwork(ChannelHeadNodes, FlowInfo);
	LSDIndexRaster SOArrayNew = NewChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	string SO_name_new = "_SO_from_CH";
	
	SOArrayNew.write_raster((path_name+DEM_name+SO_name_new),DEM_flt_extension);	
                                
}
