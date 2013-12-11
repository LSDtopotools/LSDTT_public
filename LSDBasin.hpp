//beginning of the LSDBasin object

#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDChannelNetwork.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDBasin_H
#define LSDBasin_H

class LSDBasin
{

  public:
  
  LSDBasin(int Junction, LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChanNet)
											{ create(Junction, FlowInfo, ChanNet); }

  // getters for the basin attributes go here
  
  /// @return Number of rows as an integer.
	int get_NRows() const				{ return NRows; }
	/// @return Number of columns as an integer.
  int get_NCols() const				{ return NCols; }
  /// @return Minimum X coordinate as an integer.
	double get_XMinimum() const			{ return XMinimum; }
	/// @return Minimum Y coordinate as an integer.
	double get_YMinimum() const			{ return YMinimum; }
	/// @return Data resolution as an integer.
	double get_DataResolution() const	{ return DataResolution; }
	/// @return No Data Value as an integer.
	int get_NoDataValue() const			{ return NoDataValue; }
											
  
  protected:
  
  ///Number of rows.
  int NRows;
  ///Number of columns.
	int NCols;
	///Minimum X coordinate.
  double XMinimum;
	///Minimum Y coordinate.
	double YMinimum;

	///Data resolution.
	double DataResolution;
	///No data value.
	int NoDataValue;
  
  
  /// Junction Number of the basin, serves as a unique ID of a basin.
  int Junction;
  
  ///Vector of all nodes in basin.
  vector<int> BasinNodes;
  
  /// Number of DEM cells. 
  int NumberOfCells;
  
  /// Area in spatial units of the basin.
  double Area;
  
  /// Boolean to show if the basin is beheaded or not.
  bool Beheaded;
  
  /// i index of the outlet pixel
  int Outlet_i;
  
  /// j index of the outlet pixel
  int Outlet_j;
  
  /// The Strahler order of the basin
  int BasinOrder;
  
  ///The i index of the centroid of the basin
  int Centroid_i;
  
  ///The j index of the centroid of the basin
  int Centroid_j;

  /// Mean basin slope.
  double SlopeMean;
  /// Mean basin elevation.
  double ElevationMean;
  /// Mean basin aspect.
  double AspectMean;
  /// Mean basin relief.
  double ReliefMean;
  /// Mean basin planform curvature.
  double PlanCurvMean;
  /// Mean basin profile curvature.
  double ProfileCurvMean;
  /// Mean basin total curvature.
  double TotalCurvMean;
  /// Max basin planform curvature.
  double PlanCurvMax;
  /// Max basin profile curvature.
  double ProfileCurvMax;
  /// Max basin total curvature.
  double TotalCurvMax;
  /// Mean hillslope length from hilltop flow routing.
  double HillslopeLength_HFR;
  /// Mean hillslope length from binned boomerang plots.
  double HillslopeLength_Binned;
  /// Mean hillslope length from spline curves on boomerang plots.
  double HillslopeLength_Spline;
  /// Mean hillslope length from drainage density.
  double HillslopeLength_Density;
  /// Basin flowlength.
  double FlowLength;
  /// Basin drainage density.
  double DrainageDensity;
  /// Basin perimiter.
  double Perimiter;
  /// Cosmo erosion rate.
  double CosmoErosionRate;
  /// Other erosion rate.
  double OtherErosionRate;
  /// Mean basin hilltop curvature.
  double CHTMean;
  /// Basin E* value.
  double EStar;
  /// Basin R* value.
  double RStar;


  											
  private:
	void create(int Junction, LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChanNet);

};

#endif



