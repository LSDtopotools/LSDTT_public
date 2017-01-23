//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//  pelletier_channel_heads_driver.cpp
//
// This program is used to predict channel head locations based on the method proposed by
// Pelletier (2013).  The method involves using a wiener filter in order to remove high freuqency
// topographic noise (modelled as white nosie in his paper).  It then creates a contour curvature 
// map and identifies channel heads as pixels greater than a user defined contour curvature threshold
// value.  The threshold curvature can also be defined as a multiple of the standard deviation of the
// curvature.
//
// Reference: Pelletier (2013) A robust, two-parameter method for the extraction of drainage
// networks from high-resolution digital elevation models (DEMs): Evaluation using synthetic and real-world
// DEMs, Water Resources Research 49: 1-15
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb
// Simon M. Mudd 
// David T. Milodowski
// 
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../fftw-3.3.1/api/fftw3.h"
#include "../TNT/jama_lu.h"
#include "../TNT/tnt.h"

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
	float Minimum_Slope;
	file_info_in >> Minimum_Slope;

  // STEP 1: Spectral filtering
	string DEM_flt_extension = "flt";
	string DEM_asc_extension = "asc";
	string DEM_outname = "fr1m_filt";

	LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
  cout << "LINE 47" << endl;
	
  LSDRasterSpectral SpectralRaster(topo_test);                                    
	//char* file_id = "test_spectral";                                                
	//float LogBinWidth = 0.1;                                                       
	//SpectralRaster.fftw2D_spectral_analysis(file_id, LogBinWidth);

	LSDRaster topo_data_filtered = SpectralRaster.fftw2D_wiener();
	// Write raster of filtered DEM
	topo_data_filtered.write_raster((path_name+DEM_outname),DEM_flt_extension);	

	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	float threshold = 50;
	// get the filled file
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(Minimum_Slope);
  
  //get a FlowInfo object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test); 
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
	
	//get the sources: note: this is only to select basins!
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold); 

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

	string tan_curvature_name = "fr1m_curv";
	string chan_heads_name = "fr1m_chan_heads";
	
	// STEP 2: Creating a contour curvature map
	
	// Calculate polyfit coefficients and the tangential curvature
  Array2D<float> a;
  Array2D<float> b;
  Array2D<float> c;
  Array2D<float> d;
  Array2D<float> e;
  Array2D<float> f;
  float window_radius=6;
  float tan_curv_threshold = 0.1;
  
  topo_data_filtered.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster tan_curvature = topo_data_filtered.calculate_polyfit_tangential_curvature(a ,b ,c ,d, e);
  // Write curvature raster
  tan_curvature.write_raster(tan_curvature_name, DEM_flt_extension);
  Array2D<float> tan_curv_array = tan_curvature.get_RasterData();
  
	// Call on the member function
  LSDIndexRaster source_nodes = topo_data_filtered.calculate_pelletier_channel_heads(window_radius, tan_curv_threshold,
                                                                                      tan_curv_array); 
	source_nodes.write_raster((path_name+chan_heads_name),DEM_flt_extension);	
	
	Array2D<int> pelletier_sources = source_nodes.get_RasterData();
  vector<int> source_vector;
  int n_sources = 0;
  int NRows = FlowInfo.get_NRows();
  int NCols = FlowInfo.get_NCols();
  int NoDataValue = FlowInfo.get_NoDataValue();
  
  for (int row = 0; row < NRows; row++)
	{
    for (int col = 0; col < NCols; col++)
    {
       if (pelletier_sources[row][col] == 1)
       {
          int current_node = FlowInfo.retrieve_node_from_row_and_column(row, col);
          if (current_node != NoDataValue)
          {
            source_vector.push_back(current_node);
            ++n_sources;
          }          
       }
    }
  }
  
  LSDJunctionNetwork ChanNetworkPel(source_vector, FlowInfo);
  LSDIndexRaster SOArrayNew = ChanNetworkPel.StreamOrderArray_to_LSDIndexRaster();
	string SO_name_new = "fr1m_SO_from_PEL";	
	SOArrayNew.write_raster((path_name+SO_name_new),DEM_flt_extension);	
}
