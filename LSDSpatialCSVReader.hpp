//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDSpatialCSVReader.hpp
// Land Surface Dynamics SpatialCSVReader
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for reading csv data. The data needs to have latitude and longitude 
//  in WGS84 coordinates.
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2017 Simon M. Mudd 2017
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
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "LSDRaster.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
using namespace std;

#ifndef LSDSpatialCSVReader_HPP
#define LSDSpatialCSVReader_HPP

class LSDSpatialCSVReader
{
  public:

    /// @brief Create an LSDSpatialCSVReader from a raster and csv filenname
    /// @param ThisRaster An LSDRaster object
    /// @param csv_fname The name of the csv file including extension  and path
    /// @author SMM
    /// @date 16/02/2017
    LSDSpatialCSVReader(LSDRaster& ThisRaster, string csv_fname)  { create(ThisRaster,csv_fname); }

    /// @brief Create an LSDChiTools from a raster.
    /// @param ThisRaster An LSDIndexRaster object
    /// @param csv_fname The name of the csv file including extension  and path
    /// @author SMM
    /// @date 16/02/2017
    LSDSpatialCSVReader(LSDRasterInfo& ThisRaster, string csv_fname)  { create(ThisRaster,csv_fname); }

    /// @brief This loads a csv file, grabbing the latitude and longitude, 
    ///  and putting the rest of the data into data maps
    /// @param filename The name of the csv file including path and extension
    /// @author SMM
    /// @date 16/02/2017
    void load_csv_data(string filename);

    /// @brief this prints the latitude and longitude to screen
    /// @author SMM
    /// @date 17/02/2017
    void print_lat_long_to_screen();


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
    
    /// A vector of the latitude (in WGS84)
    vector<double> latitude;
    
    /// A vector of the longitude (in WGS84)
    vector<double> longitude;
    
    /// The map to hold all the read data
    map<string, vector<string> > data_map;
  
  private:
  
    void create(LSDRaster&,string);
    
    void create(LSDRasterInfo&,string);
    
  
    

};

#endif