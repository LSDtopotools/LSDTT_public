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
#include <cmath>
#include "LSDRaster.hpp"
#include "LSDSoilHydroRaster.hpp"
#include "LSDStatsTools.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDShapeTools.hpp"
#include "LSDFlowInfo.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// empty default create function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::create()
{}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Create function that just copies a raster into the hydro raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Create function that takes the dimensions and georeferencing of a raster
// but then sets all data to value, setting the NoDataValues to 
// the NoData of the raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
  Array2D<float> data(NRows,NCols,NoDataValue); 
  
  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (OtherRaster.get_data_element(row,col) != NoDataValue)
      {
        data[row][col] = value;
      }
    }
  }
  
  RasterData = data.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates a raster from raw data
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates a raster from raw data, this time with the georeferencing strings
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets all non nodata pixels to value
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::SetHomogenousValues(float value)
{
  // set the raster data to be a certain value
  Array2D<float> data(NRows,NCols,NoDataValue); 
  
  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        RasterData[row][col] = value;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates a snow thickenss (effective, in g cm^-2 for cosmogenic 
// applications)  based on a bilinear model such as that of P Kirchner: http://escholarship.org/uc/item/9zn1c1mk#page-8
// The paper is here: http://www.hydrol-earth-syst-sci.net/18/4261/2014/hess-18-4261-2014.html
// This paper also agrees withy this general trend:
// http://www.the-cryosphere.net/8/2381/2014/tc-8-2381-2014.pdf
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::SetSnowEffDepthBilinear(float SlopeAscend, float SlopeDescend, 
                          float PeakElevation, float PeakSnowpack, LSDRaster& Elevation)
{
  float LocalElevation;
  float ascendEffDepth;
  float descendEffDepth;
  float thisEffDepth = 0;
  
  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        LocalElevation = Elevation.get_data_element(row,col);
        
        if (LocalElevation != NoDataValue)
        {
          // get the effective depth on both the ascending and descending limb
          ascendEffDepth = SlopeAscend*(LocalElevation-PeakElevation)+PeakSnowpack;
          descendEffDepth = SlopeDescend*(LocalElevation-PeakElevation)+PeakSnowpack;
          
          // the correct depth is the lesser of the two
          if (ascendEffDepth < descendEffDepth)
          {
            thisEffDepth =  ascendEffDepth;
          }
          else
          {
            thisEffDepth = descendEffDepth;
          }
          
          // if the depth is less than zero, then set to zero
          if(thisEffDepth <0)
          {
            thisEffDepth = 0;
          }
          
          RasterData[row][col] =  thisEffDepth;

        }
        else        // if there ins't any elevation data, set the snow data to NoData
        {
          RasterData[row][col] = NoDataValue;
        }
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates a snow thickenss (effective, in g cm^-2 for cosmogenic 
// applications)  based on a richard's equsion sigmoidal growth model
// It was propoesd to represent peak SWE so we cruedly apply it to average annual SWE
// see 
// http://onlinelibrary.wiley.com/doi/10.1002/2015GL063413/epdf
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::SetSnowEffDepthRichards(float MaximumEffDepth, float MaximumSlope, float v, 
                          float lambda, LSDRaster& Elevation)
{
  // Don't let V be less than or equal to zero
  if (v <= 0)
  {
    v = 0.001;
  }

  // some variables to speed up compuation
  float exp_term;
  float thisEffDepth = 0;
  float elev_mulitplier = (MaximumSlope/MaximumEffDepth)*pow((1+v),1+(1/v));
  float LocalElevation;
  
  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        LocalElevation = Elevation.get_data_element(row,col);
        
        if (LocalElevation != NoDataValue)
        {
          // get the effective depth using the richards sigmoidal gorth function
          exp_term = 1+v*exp(elev_mulitplier*(lambda-LocalElevation));
          thisEffDepth = MaximumEffDepth*pow(exp_term,-(1/v));
          
          // if the depth is less than zero, then set to zero
          if(thisEffDepth <0)
          {
            thisEffDepth = 0;
          }
          
          // update the data
          RasterData[row][col] =  thisEffDepth;

        }
        else        // if there ins't any elevation data, set the snow data to NoData
        {
          RasterData[row][col] = NoDataValue;
        }
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is an incredibly rudimentary function used to modify landslide raster
// It takes a few rasters from the 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::NaiveLandslide(LSDRaster& FilledElevation, int initiationPixels,
                                      int MinPixels, float landslide_thickness)
{
  // Get a flow info object
  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";
  
  
  // some values from the rasters
  float local_elev;
  float local_mask;
  
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,FilledElevation);  
  
  // get the contributing pixels
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  vector<int> sources = FlowInfo.get_sources_index_threshold(ContributingPixels, initiationPixels);
  
  // get a value vector for the landslides
  vector<float> landslide_thicknesses;
  for (int i = 0; i< int(sources.size()); i++)
  {
    landslide_thicknesses.push_back(landslide_thickness);
  }
  
  
  // get the mask
  LSDRaster Mask = FlowInfo.get_upslope_node_mask(sources,landslide_thicknesses);
  
  // now set all points that have elevation data but not landslide data to 
  // the value of the landslide thickness, removing data that is below the minium
  // pixel area
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      local_elev =  FilledElevation.get_data_element(row,col);
      local_mask =  Mask.get_data_element(row,col);
      
      RasterData[row][col] = local_mask;
      
      // Turn nodata points into 0s
      if( local_mask == NoDataValue)
      {
        RasterData[row][col] = 0.0;
      }
      
      // remove data where there is no topographic information
      if( local_elev == NoDataValue)
      {
        RasterData[row][col] = NoDataValue;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#endif