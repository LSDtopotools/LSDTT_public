#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDChannelNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDBasin.hpp"

using namespace std;
using namespace TNT;

#ifndef LSDBasin_CPP
#define LSDBasin_CPP


void LSDBasin::create(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDChannelNetwork& ChanNet){

  //NO BOUNDS CHECKING ON JunctionNumber

  //setting all of the instance variables for the given junction

  NRows = ChanNet.NRows;
	NCols = ChanNet.NCols;
	XMinimum = ChanNet.XMinimum;
	YMinimum = ChanNet.YMinimum;
	DataResolution = ChanNet.DataResolution;
	NoDataValue = ChanNet.NoDataValue;

  Junction = JunctionNumber;
  
  LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction, ChanNet.JunctionVector[Junction],
                                                     ChanNet.ReceiverVector[Junction], ChanNet.JunctionVector[ChanNet.ReceiverVector[Junction]], FlowInfo);

  int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  int basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  BasinNodes = FlowInfo.get_upslope_nodes(basin_outlet);
  
  NumberOfCells = BasinNodes.size();
  Area = NumberOfCells * (DataResolution*DataResolution);
  
  Beheaded = ChanNet.node_tester(FlowInfo, Junction);

  FlowInfo.retrieve_current_row_and_col(ChanNet.get_Node_of_Junction(Junction), Outlet_i, Outlet_j);
  
  
  BasinOrder = ChanNet.StreamOrderVector[Junction];


  int i_max = 0;
  int i_min = 9999999; //a very large number
  int j_max = 0;
  int j_min = 9999999; //a very large number
  
  int i = 0;
  int j = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
    if (i > i_max){i_max = i;}
    else if (i < i_min){i_min = i;}
    if (j > j_max){j_max = j;}
    else if (j < j_min){j_min = j;}
    
  }
  
  Centroid_i = (i_max - i_min)/2;
  Centroid_j = (j_max - j_min)/2;   //how do these handle 0.5s ??



  //finished setting all the instance variables
  
  
  // now we set up empty variables to store properties of the basin
  // these are populated as they are required using the set methods in LSDBasin
  
  
  SlopeMean = NoDataValue;
  ElevationMean = NoDataValue;
  AspectMean = NoDataValue;
  ReliefMean = NoDataValue;
  PlanCurvMean = NoDataValue;
  ProfileCurvMean = NoDataValue;
  TotalCurvMean = NoDataValue;
  PlanCurvMax = NoDataValue;
  ProfileCurvMax = NoDataValue;
  TotalCurvMax = NoDataValue;
  HillslopeLength_HFR = NoDataValue;
  HillslopeLength_Binned = NoDataValue;
  HillslopeLength_Spline = NoDataValue;
  HillslopeLength_Density = NoDataValue;
  FlowLength = NoDataValue;
  DrainageDensity = NoDataValue;
  Perimiter = NoDataValue;
  CosmoErosionRate = NoDataValue;
  OtherErosionRate = NoDataValue;
  CHTMean = NoDataValue;
  EStar = NoDataValue;
  RStar = NoDataValue;
   
  //finished creating empty variables 

}


#endif