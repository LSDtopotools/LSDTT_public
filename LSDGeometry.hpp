//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDGeometry
// Land Surface Dynamics LSDGeometry objects
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for dealing with geometric data
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
// Copyright (C) 2016 Simon M. Mudd 2016
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
// either version 3 of the License, or (at your option) any later version.
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
// LSDGeometry.cpp
// LSDGeometry object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDGeometry_HPP
#define LSDGeometry_HPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
using namespace std;
using namespace TNT;


/// @brief This object packages a number of tools for chi analysis
class LSDGeometry
{
  public:

    /// @brief Empty create function. Leads to some empty vectors
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry()  { create(); }

    /// @brief Create with two vectors
    /// @detail The X should be either Easting or Longitude
    ///  the Y vector should be either Northing or latitude
    ///  The function checks which one it is. 
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<double> x, vector<double> y)  { create(x,y); }
 
    /// @brief Create with two vectors
    /// @detail The X should be either Easting or Longitude
    ///  the Y vector should be either Northing or latitude
    ///  The function checks which one it is. 
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<float> x, vector<float> y)  { create(x,y); }
 
    /// @brief Create with two vectors
    /// @detail The X should be either Easting or Longitude
    ///  the Y vector should be either Northing or latitude
    ///  The function checks which one it is. 
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<double> x, vector<double> y, string CoordSys)  { create(x,y,CoordSys); }
    

  protected:

    vector<double> UTMPoints_Easting;
    vector<double> UTMPoints_Northing;
    
    vector<double> WGS84Points_latitude;
    vector<double> WGS84Points_longitude;

  private:
    void create();
    void create(vector<double> x, vector<double> y);
    void create(vector<float> x, vector<float> y);
    void create(vector<double> x, vector<double> y,string coord_sys);

};

#endif