//beginning of the LSDFloodplain class

#ifndef LSDFloodplain_H
#define LSDFloodplain_H

#include <vector>
#include <string>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

///@brief Class to store information about floodplains and terraces...
class LSDFloodplain
{
  public:

  /// @brief create function
  ///
  /// @details This populates the binary array and connected components array for the floodplain
	/// given rasters of channel relief and slope and thresholds for both. Each pixel
  /// must be below the slope and channel relief threshold to be classified as floodplain.
  /// @author FJC
	/// 18/10/16
  LSDFloodplain(LSDRaster& ChannelRelief, LSDRaster& Slope, float relief_threshold, float slope_threshold, int min_patch_size)
					{ create(ChannelRelief, Slope, relief_threshold, slope_threshold, min_patch_size); }
	
	/// @return Number of rows as an integer.
  int get_NRows() const        { return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const        { return NCols; }
  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const        { return XMinimum; }
  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const        { return YMinimum; }
  /// @return Data resolution as an integer.
  float get_DataResolution() const        { return DataResolution; }
  /// @return No Data Value as an integer.
  int get_NoDataValue() const        { return NoDataValue; }
	/// @return Georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }  
	
	/// @brief This function gets the information about all the floodplain pixels connected to the main stem channel from a junction
	/// @details Takes a junction number and generates the main stem channel from this point. THe information about each floodplain or terrace pixel is then calculated relative to the main channel.
	/// @param junction_number junction number of interest
	/// @param ChanNetwork LSDJunctionNetwork object
	/// @param FlowInfo LSDFlow info object
	/// @param DistFromOutlet LSDRaster of distances from the outlet
	/// @param ElevationRaster LSDRaster of elevations
	/// @author FJC
	/// @date 18/10/16
	void get_main_stem_information(int junction_number, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet, LSDRaster& ElevationRaster);
	
	/// @brief This function separates te floodplain into floodplain patches and terrace patches depending 
	/// on whether they are connected to the channel network.
	/// @param ChanNetwork LSDJunctionNetwork object
	/// @param threshold_SO threshold SO to count as the channel network
	/// @param FloodplainPatches Empty LSDIndexRaster to store floodplain patches
	/// @param TerracePatches Empty LSDIndexRaster to store terrace patches
	/// @author FJC
	/// @date 18/10/16
	void separate_floodplain_and_terrace_patches(LSDJunctionNetwork& ChanNetwork, float threshold_SO, LSDIndexRaster& FloodplainPatches, LSDIndexRaster& TerracePatches);
	
	/// FUNCTIONS TO GENERATE RASTERS
	
	/// @brief This function prints the channel relief compared to the main stem to a raster
	/// @return ChannelRelief LSDRaster of channel relief
	/// @author FJC
	/// @date 18/10/16
	LSDRaster print_ChannelRelief_to_Raster();
	
	/// @brief This function prints the upstream distance compared of the nearest main stem to a raster
	/// @return UpstreamDist LSDRaster of upstream distance
	/// @author FJC
	/// @date 19/10/16
	LSDRaster print_UpstreamDistance_to_Raster();
	
	/// @brief This function prints the flow lengths to the nearest main stem to a raster
	/// @return Flow Lengths LSDRaster of flow length
	/// @author FJC
	/// @date 19/10/16
	LSDRaster print_FlowLengths_to_Raster();

  protected:
	
	/// Number of rows
	int NRows;
	/// Number of columns
	int NCols;
	/// X minimum
	float XMinimum;
	/// Y minimum
	float YMinimum;
	/// Data resolution
	float DataResolution;
	/// No data value
	int NoDataValue;
	/// A map of strings for holding georeferencing information
  map<string,string> GeoReferencingStrings;
		
	/// Relief threshold
	float relief_threshold;
	/// Slope threshold
	float slope_threshold;
	
  /// The binary array of floodplain data
	Array2D<int> BinaryArray;
	/// The array of connected components
	Array2D<int> ConnectedComponents_Array;
	
	/// Junction number - used to calculate information about the patches compared to the main stem channel
	int JunctionNumber;
	
	/// These members are set by the get_mainstem_information function for a specific junction
	
	/// vector of nodes upslope of the junction
	vector<int> UpslopeNodes;
	/// array of the nearest NI on the main stem to each pixel
	Array2D<int> MainStemNIs;
	/// array of relief relative to main stem
	Array2D<float> ChannelRelief_array;
	/// array of distance upstream along main stem
	Array2D<float> UpstreamDistance_array;
	/// array of flow lengths from the main stem
	Array2D<float> FlowLength_array;
	
  private:
	void create(LSDRaster& ChannelRelief, LSDRaster& Slope, float relief_threshold, float slope_threshold, int min_patch_size);

};

#endif
