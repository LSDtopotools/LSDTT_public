//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get_CHT.cpp
//
// This program calculates channel heads using an Optimal wiener filter (used by Pelletier,
// 2013) and a quantile-quantile curvature threshold similar to Geonet.  It then uses a
// connected components threshold to create a channel network from the curvature mask.
//
// Following the channel extraction, the driver computes hilltop curvature,
// extracts basins from a list of lat/long coordinates and calculates basin average
// CHT statistics.
//
// References: Pelletier, J.D. (2013) A robust, two-parameter method for the extraction of
// drainage networks from high-resolution digital elevation models (DEMs): Evaluation using
// synthetic and real-world DEMs, Water Resources Research 49(1): 75-89, doi:10.1029/2012WR012452
//
// Passalacqua, P., Do Trung, T., Foufoula-Georgiou, E., Sapiro, G., & Dietrich, W. E.
// (2010). A geometric framework for channel network extraction from lidar: Nonlinear diffusion and
// geodesic paths. Journal of Geophysical Research: Earth Surface (2003�2012), 115(F1).
//
// He, L., Chao, Y., & Suzuki, K. (2008). A run-based two-scan labeling algorithm. Image Processing,
// IEEE Transactions on, 17(5), 749-756.
//
// Developed by:
//  David T. Milodowski
//  Simon M. Mudd
//  Stuart W.D. Grieve
//  Fiona J. Clubb
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
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Stuart W.D. Grieve, University of Edinburgh
// Fiona J. Clubb, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDIndexRaster.hpp"
#include "../TNT/tnt.h"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDBasin.hpp"
int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "FATAL ERROR: wrong number of inputs. The program needs the path name, the driver file name" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  string full_name = path_name+f_name;

  ifstream file_info_in;
  file_info_in.open(full_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  string Raster_name;
  string Output_name;
  string q_q_filename_prefix;
  float area_threshold;
  float window_radius;
  string DEM_extension = "bil";
  string temp;
  int BasinOrder;
  float slope_threshold;
  string ll_filename;
  int connected_components_threshold;
  int threshold_stream_order;
  int search_radius_nodes;
  int UTM_zone;
  int eId;
  bool isNorth;
  float roughness_threshold;
  float roughness_radius;

  file_info_in >> temp >> Raster_name
               >> temp >> Output_name
               >> temp >> q_q_filename_prefix
               >> temp >> window_radius
               >> temp >> area_threshold
	             >> temp >> connected_components_threshold
               >> temp >> BasinOrder
               >> temp >> ll_filename
               >> temp >> slope_threshold
               >> temp >> threshold_stream_order
               >> temp >> search_radius_nodes
               >> temp >> UTM_zone
               >> temp >> eId
               >> temp >> isNorth
               >> temp >> roughness_threshold
               >> temp >> roughness_radius;
  file_info_in.close();

  // Now get the lat long data
  ifstream ll_data;
  ll_data.open(ll_filename.c_str());

  // variables for ingesting data
  string temp_string;
  float Lat;
  float Long;
  int id;

  // the LL data
  vector<float> latitude;
  vector<float> longitude;
  vector<int> IDs;

  // skip the header
  ll_data >> temp_string;
  char delimiter;

  // now get the data
  while ((ll_data >> id >> delimiter >> Lat >> delimiter >> Long) && (delimiter == ',')){
    latitude.push_back(Lat);
    longitude.push_back(Long);
    IDs.push_back(id);
  }
  ll_data.close();

  //Start of channel extraction
  LSDRasterSpectral raster(Raster_name, DEM_extension);
  LSDIndexRaster connected_components = raster.IsolateChannelsWienerQQ(area_threshold, window_radius, q_q_filename_prefix+".txt");
  LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(connected_components_threshold);
  LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();
  LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();
  LSDIndexRaster Ends = skeleton_raster.find_end_points();
  Ends.remove_downstream_endpoints(CC_raster, raster);

  //Load the elevation data, fill it and generate a FlowInfo object
  LSDRaster DEM(Raster_name, DEM_extension);
  float MinSlope = 0.0001;
  LSDRaster FilledDEM = DEM.fill(MinSlope);
  vector<string> BoundaryConditions(4, "No Flux");
  LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

  //this processes the end points to only keep the upper extent of the channel network
  vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);

  // we need a temp junction network to search for single pixel channels
  LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
  LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

  //Now we have the final channel heads, so we can generate a channel network from them
  vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);
  LSDJunctionNetwork JunctionNetwork(FinalSources, FlowInfo);

  //Finally we write the channel heads to file so they can be used in other drivers.
  FlowInfo.print_vector_of_nodeindices_to_csv_file(FinalSources,(Output_name+"_CH"));

  vector< int > basin_junctions = JunctionNetwork.ExtractBasinJunctionOrder(BasinOrder, FlowInfo);
  LSDIndexRaster Basin_Raster = JunctionNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);

  //surface fitting - this is the second call of this method, this time on the unfiltered topography
  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(0); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(0); //plan curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  vector<LSDRaster> Surfaces = FilledDEM.calculate_polyfit_surface_metrics(window_radius, raster_selection);

  //raster_selection.clear();
  vector<int> raster_selection2;

  raster_selection2.push_back(0);
  raster_selection2.push_back(0);
  raster_selection2.push_back(1);

  vector<LSDRaster> Roughness = FilledDEM.calculate_polyfit_roughness_metrics(window_radius, roughness_radius, raster_selection2);

  // extract hilltops
  LSDRaster hilltops = JunctionNetwork.ExtractRidges(FlowInfo);

  //get hilltop curvature using filter to remove positive curvatures
  LSDRaster cht_raster = FilledDEM.get_hilltop_curvature(Surfaces[3], hilltops);
  LSDRaster CHT = FilledDEM.remove_positive_hilltop_curvature(cht_raster);

  //perform a gradient filter on the CHT data
  LSDRaster CHT_gradient = JunctionNetwork.ExtractHilltops(CHT, Surfaces[1], slope_threshold);

  //Set up a converter object to transform Lat Long into UTM
  LSDCoordinateConverterLLandUTM Converter;

  int N_samples = int(latitude.size());  // number of coordinate pairs loaded

  // set up some vectors
  vector<float> fUTM_northing(N_samples,0);
  vector<float> fUTM_easting(N_samples,0);

  double this_Northing;
  double this_Easting;

  // loop through the coordinates, converting to UTM
  for(int i = 0; i<N_samples; i++)
  {
    Converter.LLtoUTM_ForceZone(eId, latitude[i], longitude[i], this_Northing, this_Easting, UTM_zone);
    fUTM_easting[i] = float(this_Easting);
    fUTM_northing[i] = float(this_Northing);
  }

  // snap to hilltop function here
  vector<int> SnappedNodes;
  vector<int> Valid_node_IDs;
  FlowInfo.snap_to_hilltops(fUTM_easting, fUTM_northing, search_radius_nodes, CHT, SnappedNodes, Valid_node_IDs);

  int n_valid_points = int(SnappedNodes.size());  //The number of points which were within the current DEM

  if (n_valid_points != N_samples){
    cout << "Not every point was located within the DEM" << endl;
  }

  //Open a file to write the basin average data to
  ofstream WriteData;
  stringstream ss2;
  ss2 << Output_name << "_CHT_Data_Ridge.csv";
  WriteData.open(ss2.str().c_str());

  //write headers
  WriteData << "_ID,min,max,median,mean,range,std_dev,std_err,min_gradient,max_gradient,median_gradient,mean_gradient,range_gradient,std_dev_gradient,std_err_gradient,bedrock_percentage" << endl;

  int threshold = 10;

  for(int samp = 0; samp<n_valid_points; samp++)
  {
    //vector of nodeindexes from new snapping fn
    //convert nodeindexes to i,j to feed into sampler

    int a;
    int b;

    FlowInfo.retrieve_current_row_and_col(SnappedNodes[samp], a, b);

    vector< vector<float> > Samples = FilledDEM.Sample_Along_Ridge(CHT, CHT_gradient, Roughness[2], a, b, threshold);

    WriteData << IDs[Valid_node_IDs[samp]];

    for (int a=0; a < 2; ++a){

      float mean = get_mean(Samples[a]);
      float max = *max_element(Samples[a].begin(), Samples[a].end());
      float min = *min_element(Samples[a].begin(), Samples[a].end());
      float median = get_percentile(Samples[a], 50.0);
      float stdd = get_standard_deviation(Samples[a], mean);
      float stde = get_standard_error(Samples[a], stdd);
      float range = get_range_from_vector(Samples[a], CHT.get_NoDataValue());

      WriteData << "," << min << "," << max << "," << median << "," << mean << "," << range << "," << stdd << "," << stde;
    }

    int bedrock_count = 0;

    for (int w=0; w < int(Samples[3].size()); ++w){
      if (Samples[3][w] > roughness_threshold){
        bedrock_count++;
      }
    }

    float bedrock_full = (bedrock_count / int(Samples[3].size())) * 100.0;
    WriteData << "," << bedrock_full << endl;
  }

  WriteData.close();

}
