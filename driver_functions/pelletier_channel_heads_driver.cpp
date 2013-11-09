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
#include "LSDStatsTools.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDChannelNetwork.hpp"
#include "fftw-3.3.1/api/fftw3.h"
#include "TNT/jama_lu.h"
#include "TNT/tnt.h"

int main ()
{

  // STEP 1: Spectral filtering
  string DEM_name = "ind_17n";
	string DEM_flt_extension = "flt";
	string DEM_asc_extension = "asc";
	string DEM_outname = "ind_filt";

	LSDRaster topo_test(DEM_name, DEM_flt_extension);
  cout << "LINE 47" << endl;
	
  LSDRasterSpectral SpectralRaster(topo_test);                                    
	//char* file_id = "test_spectral";                                                
	//double LogBinWidth = 0.1;                                                       
	//SpectralRaster.fftw2D_spectral_analysis(file_id, LogBinWidth);

	LSDRaster topo_data_filtered = SpectralRaster.fftw2D_wiener();
	// Write raster of filtered DEM
	topo_data_filtered.write_raster(DEM_outname,DEM_flt_extension);	

	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";
	
	double Minimum_Slope = 0.0001;	
	double threshold = 50;
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
	LSDChannelNetwork ChanNetwork(sources, FlowInfo);

	string tan_curvature_name = "ind_curv";
	string chan_heads_name = "ind_chan_heads";
	
	// STEP 2: Creating a contour curvature map
	
	// Calculate polyfit coefficients and the tangential curvature
  Array2D<double> a;
  Array2D<double> b;
  Array2D<double> c;
  Array2D<double> d;
  Array2D<double> e;
  Array2D<double> f;
  double window_radius=6;
  double tan_curv_threshold = 0.1;
  
  topo_data_filtered.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster tan_curvature = topo_data_filtered.calculate_polyfit_tangential_curvature(a ,b ,c ,d, e);
  // Write curvature raster
  tan_curvature.write_raster(tan_curvature_name, DEM_flt_extension);
  Array2D<double> tan_curv_array = tan_curvature.get_RasterData();
  
	// Call on the member function
  LSDIndexRaster source_nodes = topo_data_filtered.calculate_pelletier_channel_heads(window_radius, tan_curv_threshold,
                                                                                      tan_curv_array); 
	source_nodes.write_raster(chan_heads_name,DEM_flt_extension);	
	
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
  
  LSDChannelNetwork ChanNetworkPel(source_vector, FlowInfo);
  LSDIndexRaster SOArrayNew = ChanNetworkPel.StreamOrderArray_to_LSDIndexRaster();
	string SO_name_new = "ind_17n_SO_from_PEL";	
	SOArrayNew.write_raster(SO_name_new,DEM_flt_extension);	
}
