//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// channel_heads_part1_driver.cpp
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
#include "LSDRasterSpectral.hpp"
#include "fftw-3.3.1/api/fftw3.h"

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
	int basin_order;

	file_info_in >> Minimum_Slope >> threshold >> A_0 >> m_over_n >> basin_order;

	// get some file names
	string DEM_f_name = path_name+DEM_name+fill_ext;
	string DEM_flt_extension = "flt";

	// load the DEM
	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
	
	//LSDRasterSpectral SpectralRaster(topo_test);                                    
	//char* file_id = "test_spectral";                                                
	//double LogBinWidth = 0.1;                                                       
	//SpectralRaster.fftw2D_spectral_analysis(file_id, LogBinWidth);

	//LSDRaster topo_data_filtered = SpectralRaster.fftw2D_wiener();
	
	// get the filled file
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
	//filled_topo_test.write_raster((DEM_f_name),DEM_flt_extension);

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);   
	
	int NRows = FlowInfo.get_NRows();
	int NCols = FlowInfo.get_NCols();
	double NoDataValue = FlowInfo.get_NoDataValue();
	double XMinimum = FlowInfo.get_XMinimum();
	double YMinimum = FlowInfo.get_YMinimum();
	double DataResolution = FlowInfo.get_DataResolution();

	// calculate the distance from outlet and contributing pixels
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  //string dist_raster_fname = "_FD";
  //DistanceFromOutlet.write_raster((path_name+DEM_name+dist_raster_fname),DEM_flt_extension);

  //string NI_fname = "_NI";
	//LSDIndexRaster NodeIndexRaster =  FlowInfo.write_NodeIndex_to_LSDIndexRaster();
  //NodeIndexRaster.write_raster((path_name+DEM_name+NI_fname),DEM_flt_extension);

  //get the sources: note: this is only to select basins!
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold); 

	// now get the junction network
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);

  int junction_number = 0;
  double bin_width = 10;
  vector<int> junction_list = ChanNetwork.extract_basins_order_outlet_junctions(basin_order, FlowInfo);
  int max_junctions = junction_list.size();
  cout << "No of third order junctions: " << max_junctions << endl;
  //Array2D<int> ChannelHeadPixelsAll(NRows,NCols,0);
  Array2D<int> source_raster_all(NRows,NCols,NoDataValue);
  string string_filename;
	string filename = "source_list";
	string dot = ".";
	string extension = "txt";
  string_filename = filename+dot+extension;
  ofstream source_list;
  source_list.open(string_filename.c_str());
  int contributing_pixels = 0;
  
  for (int i=0; i<max_junctions; i++)
  {
    junction_number = junction_list[i];
    vector<int> SourceNodes = ChanNetwork.GetSourceNodesChiMethodAllPixels(junction_number, A_0, m_over_n, bin_width,
                                                                         FlowInfo, filled_topo_test);    
    
    int row, col;
    Array2D<int> source_raster(NRows,NCols,NoDataValue);
    for (int i=0; i < SourceNodes.size(); i++)
    {
       FlowInfo.retrieve_current_row_and_col(SourceNodes[i], row, col);
       source_raster[row][col] = SourceNodes[i]; 
       contributing_pixels = FlowInfo.retrieve_contributing_pixels_of_node(SourceNodes[i]);
       source_list << SourceNodes[i] << " " << contributing_pixels << endl;   
    }
    
    for (int row = 0; row < NRows; row++)
	  {
      for (int col = 0; col < NCols; col++)
      {
      if (source_raster[row][col] != NoDataValue) 
        {
            source_raster_all[row][col] = source_raster[row][col];
        }
      }
    }
  }  
  
  source_list.close();  
  //LSDIndexRaster source_pixels_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,source_raster_all);
  //string chan_heads_name = "_CH";
  //source_pixels_raster.write_raster((path_name+DEM_name+chan_heads_name),DEM_flt_extension); 
}