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

    /// @brief Create with two vectors: no UTM provided so assumed lat/long
    /// @detail The X should be Longitude
    ///  the Y vector should be latitude
    /// @param x The longitdue data in a vector
    /// @param y the latitude data in a vector
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<double> x, vector<double> y)  { create(x,y); }
 
    /// @brief Create with two vectors: no UTM provided so assumed lat/long
    /// @detail The X should be Longitude
    ///  the Y vector should be latitude
    /// @param x The longitdue data in a vector
    /// @param y the latitude data in a vector
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<float> x, vector<float> y)  { create(x,y); }
 
    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    ///  This version assumes northern hemisphere 
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<double> x, vector<double> y, int UTMZone)  { create(x,y,UTMZone); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing 
    ///  This version assumes northern hemisphere 
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<float> x, vector<float> y, int UTMZone)  { create(x,y,UTMZone); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing 
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<double> x, vector<double> y, int UTMZone, bool isNorth)  { create(x,y,UTMZone, isNorth); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing 
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<float> x, vector<float> y, int UTMZone, bool isNorth)  { create(x,y,UTMZone, isNorth); }

    /// @brief This function converts points from Lat/Long to UTM.
    ///  The UTM zone is set as the zone of the first data point
    ///  If there is no data in the Lat/Long data vectors no action is taken
    /// @author SMM
    /// @date 10/06/2016
    void convert_points_to_UTM();

    /// @brief This function converts points from UTM to Lat/Long.
    ///  The UTM zone is set as the zone of the first data point
    ///  If there is no data in the Lat/Long data vectors no action is taken
    /// @author SMM
    /// @date 13/06/2016
    void convert_points_to_LatLong();

    /// @brief This prints the points to a csv file. It will contain both UTM and
    ///  lat-long coordinates. The UTM zone is the zone of the first point, 
    ///  The lat long coordinates are in WGS84
    /// @param  path The path to the outfile. Needs the trailing slash
    /// @param file_prefix The prefix of the file **before extension**. That is, 
    ///  this function will add the .csv to the end of the filename
    /// @author SMM
    /// @date 13/06/2016
    void print_points_to_csv(string path, string file_prefix);

    // the getter functions
    int get_UTMZone() { return UTMZone; }
    bool get_isNorth() { return isNorth; }
    vector<double> get_UTMPoints_Easting() { return UTMPoints_Easting; }
    vector<double> get_UTMPoints_Northing() { return UTMPoints_Northing; }
    vector<double> get_WGS84Points_latitude() { return WGS84Points_latitude; }
    vector<double> get_WGS84Points_longitude() { return WGS84Points_longitude; }


  protected:
    
    int UTMZone;
    bool isNorth;

    vector<double> UTMPoints_Easting;
    vector<double> UTMPoints_Northing;
    
    vector<double> WGS84Points_latitude;
    vector<double> WGS84Points_longitude;
    

    

  private:
    void create();
    void create(vector<double> x, vector<double> y);
    void create(vector<float> x, vector<float> y);
    void create(vector<double> x, vector<double> y, int UTMZone);
    void create(vector<float> x, vector<float> y, int UTMZone);
    void create(vector<double> x, vector<double> y, int UTMZone, bool isNorth);
    void create(vector<float> x, vector<float> y, int UTMZone, bool isNorth);

};

#endif