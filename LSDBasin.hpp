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
	
	/// @return Nodes of the basin.
	vector<int> get_BasinNodes() const { return BasinNodes; }
	
	/// @return Number of cells in the basin.
	int get_NumberOfCells() const { return NumberOfCells; }
	
  /// @return Area of basin in spataial units.
  double get_Area() const { return Area; }
  
  /// @return Junction Number of the basin.
  int get_Junction() const { return Junction; }
  
  /// @return Stream order of basin.
  int get_BasinOrder() const { return BasinOrder; }
  
  /// @return Boolean value of whether a basin is beheaded or not.										
  bool get_Beheaded() const { return Beheaded; }
  
  /// @return i index of outlet pixel.
  int get_Outlet_i() const { return Outlet_i; }

  /// @return j index of outlet pixel.
  int get_Outlet_j() const { return Outlet_j; }
  
  ///@return i index of Centroid pixel.
  int get_Centroid_i() const { return Centroid_i; }

  ///@return j index of Centroid pixel.
  int get_Centroid_j() const { return Centroid_j; }
  
  //Getters of basin parameters
  
  /// @return Mean slope.
  double get_SlopeMean() const { return SlopeMean; }
  /// @return Mean elevation.
  double get_ElevationMean() const { return ElevationMean; }
  /// @return Mean aspect.
  double get_AspectMean() const { return AspectMean; }
  /// @return Mean relief.
  double get_ReliefMean() const { return ReliefMean; }
  /// @return Mean plan curvature.
  double get_PlanCurvMean() const { return PlanCurvMean; }
  /// @return Mean profile curvature.
  double get_ProfileCurvMean() const { return ProfileCurvMean; }
  /// @return Mean total curvature.
  double get_TotalCurvMean() const { return TotalCurvMean; }
  /// @return Max plan curvature.
  double get_PlanCurvMax() const { return PlanCurvMax; }
  /// @return Max profile curvature.
  double get_ProfileCurvMax() const { return ProfileCurvMax; }
  /// @return Max total curvature.
  double get_TotalCurvMax() const { return TotalCurvMax; }
  /// @return Hillslope length from hilltop flow routing.
  double get_HillslopeLength_HFR() const { return HillslopeLength_HFR; }
  /// @return Hillslope length from boomerang bins.
  double get_HillslopeLength_Binned() const { return HillslopeLength_Binned; }
  /// @return Hillslope length from boomerang splines.
  double get_HillslopeLength_Spline() const { return HillslopeLength_Spline; }
  /// @return Hillslope length from drainage density.
  double get_HillslopeLength_Density() const { return HillslopeLength_Density; }
  /// @return Flow length.
  double get_FlowLength() const { return FlowLength; }
  /// @return Drainage Density.
  double get_DrainageDensity() const { return DrainageDensity; }
  /// @return Basin perimiter.
  double get_Perimiter() const { return Perimiter; }
  /// @return Cosmo erosion rate.
  double get_CosmoErosionRate() const { return CosmoErosionRate; }
  /// @return Other eroision rate.
  double get_OtherErosionRate() const { return OtherErosionRate; }
  /// @return Mean hilltop curvature.
  double get_CHTMean() const { return CHTMean; }
  /// @return E* value.
  double get_EStar() const { return EStar; }
  /// @return R* value.
  double get_RStar() const { return RStar; }  
  
      
  protected:
  
  //These instance variables are set at initialisation
  
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
  
  
  //These instance variables are used to store calculated basin parameters 
  
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



