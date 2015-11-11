//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDSoilHydroRaster
// Land Surface Dynamics Raster for manipulating soil and hydroligcal data
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing raster data, with a particular focus on topography
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDSoilHydroRaster.cpp
// cpp file for the LSDSoilHydroRaster object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
// Stuart Grieve, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 1.0.0		12/11/2015
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDSoilHydroRaster_CPP
#define LSDSoilHydroRaster_CPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include "LSDRaster.hpp"
#include "LSDSoilHydroRaster.hpp"
#include "LSDStatsTools.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDShapeTools.hpp"
using namespace std;


// empty default create function
void LSDSoilHydroRaster::create()
{}

// Create function that just copies a raster into the hydro raster
void LSDSoilHydroRaster::create(LSDRaster& OtherRaster)
{
  NRows = OtherRaster.get_NRows();
  NCols = OtherRaster.get_NCols();
  XMinimum = OtherRaster.get_XMinimum();
  YMinimum = OtherRaster.get_YMinimum();
  DataResolution = OtherRaster.get_DataResolution();
  NoDataValue = OtherRaster.get_NoDataValue();
  GeoReferencingStrings = OtherRaster.get_GeoReferencingStrings();
  RasterData = OtherRaster.get_RasterData();
}

// Create function that takes the dimensions and georeferencing of a raster
// but then sets all data to value
void LSDSoilHydroRaster::create(LSDRaster& OtherRaster, float value)
{
  NRows = OtherRaster.get_NRows();
  NCols = OtherRaster.get_NCols();
  XMinimum = OtherRaster.get_XMinimum();
  YMinimum = OtherRaster.get_YMinimum();
  DataResolution = OtherRaster.get_DataResolution();
  NoDataValue = OtherRaster.get_NoDataValue();
  GeoReferencingStrings = OtherRaster.get_GeoReferencingStrings();

  // set the raster data to be a certain value
  Array2D<float> data(NRows,NCols,value); 
  RasterData = data.copy();

}

// Creates a raster from raw data
void LSDSoilHydroRaster::create(int ncols, int nrows, float xmin, float ymin,
                float cellsize, float ndv, Array2D<float> data)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;

  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
}

// Creates a raster from raw data, this time with the georeferencing strings
void LSDSoilHydroRaster::create(int ncols, int nrows, float xmin, float ymin,
            float cellsize, float ndv, Array2D<float> data, map<string,string> GRS)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;
  GeoReferencingStrings = GRS;

  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
}             

#endif