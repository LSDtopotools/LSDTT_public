     //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStrahlerLinks
// Land Surface Dynamics StrahlerLinks
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  that keeps track of Strahler ordered stream links and computes
//  various statistics on them
//
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2014 Simon M. Mudd 2014
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStrahlerLinks.hpp
// LSDStrahlerLinks object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.1.0		26/10/2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDStrahlerLinks_HPP
#define LSDStrahlerLinks_HPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDStrahlerLinks.hpp"
using namespace std;

///@brief Object to create a channel network from an LSDFlowInfo object.
class LSDStrahlerLinks
{
	public:
    /// @brief This defines a Strahler links object, is empty
    /// @author SMM
    /// @date 26/10/14
	  LSDStrahlerLinks()   { create(); }
   
  /// @brief This creates the LSDStrahlerLinks object
  /// @param JNetwork a junction network object
  /// @param FlowInfo LSDFlowInfo object.
  /// @param Sources vector of source nodes.
  /// @author SMM
  /// @date 26/10/14
	LSDStrahlerLinks(LSDJunctionNetwork& JNetwork, LSDFlowInfo& FlowInfo)
											{ create(JNetwork, FlowInfo); }
	
	protected:
    ///Number of rows.
    int NRows;
    ///Number of columns.
    int NCols;
    ///Minimum X coordinate.
    float XMinimum;
    ///Minimum Y coordinate.
    float YMinimum;

    ///Data resolution.
    float DataResolution;
    ///No data value.
    int NoDataValue;

    ///A map of strings for holding georeferencing information
    map<string,string> GeoReferencingStrings;
    
	  /// a vec vec containing the sources of all the Strahler links
    vector< vector<int> > SourceJunctions;
    
    /// a vec vec containing the end junctions of te Strahler links
    vector< vector<int> > ReceiverJunctions;
    
    	
	private:
	  void create();
	  void create(LSDJunctionNetwork& JN, LSDFlowInfo& FI);
};

#endif