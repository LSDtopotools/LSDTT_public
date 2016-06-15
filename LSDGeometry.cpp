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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDGeometry.cpp
// LSDGeometry object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDGeometry_CPP
#define LSDGeometry_CPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDGeometry.hpp"
using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry from an LSDRaster 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create()
{
  // this does nothing
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry object. This needs to be lat-long data (since it)
// doesn't come with UTM information
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create(vector<double> x, vector<double> y)
{
  // This reads in x and y vectors. 
  // It judges whether or not these are UTM coordinates by checking if
  // any of the data elements are over the maximum or minimum of latitude
  // or longitude
  int n_x_points = int(x.size());
  int n_y_points = int(y.size());
  
  bool is_UTM = false;
  
  if (n_x_points != n_y_points)
  {
    cout << "X and Y vectors not the same size! No data loaded. " << endl;
  }
  else if (n_x_points == 0)
  {
    cout << "Your data vectors are empty! " << endl;
  }
  else
  {
    int i = 0;
    do
    {
      // Check if latitude/Northing outside of 90: if it is this is UTM
      if (y[i] > 90 || y[i] <-90)
      {
        is_UTM = true;
      }
      
      // Check if longitude/easting outside of 180: if it is this is UTM
      if((x[i] > 180 || x[i] <-180))
      {
        is_UTM = true;
      }
      i++;
    } while (not is_UTM || i < n_x_points);
  }
  
  if(is_UTM)
  {
    cout << "This appears to be UTM data but you have not provided a zone!" << endl;
  }
  else
  {
    WGS84Points_longitude = x; 
    WGS84Points_latitude = y; 
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry object. This needs to be lat-long data (since it)
// doesn't come with UTM information
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create(vector<float> x, vector<float> y)
{
  int n_x_points = int(x.size());
  int n_y_points = int(y.size());
  vector<double> x_double;
  vector<double> y_double;
  if (n_x_points != n_y_points)
  {
    cout << "X and Y vectors not the same size! No data loaded. " << endl;
  }
  for(int i = 0; i<n_x_points; i++)
  {
    x_double.push_back( double(x[i]));
    y_double.push_back( double(y[i]));
  }
  
  create(x_double, y_double);
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry object. This is in UTM format and so has UTM information 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create(vector<double> x, vector<double> y, int thisUTMzone)
{
  // This reads in x and y vectors. 
  // It judges whether or not these are UTM coordinates by checking if
  // any of the data elements are over the maximum or minimum of latitude
  // or longitude
  int n_x_points = int(x.size());
  int n_y_points = int(y.size());
  
  if (n_x_points != n_y_points)
  {
    cout << "X and Y vectors not the same size! No data loaded. " << endl;
  }
  else if (n_x_points == 0)
  {
    cout << "Your data vectors are empty! " << endl;
  }
  else
  {
    UTMPoints_Easting = x; 
    UTMPoints_Northing = y;
    UTMZone = thisUTMzone;
    isNorth = true; 
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry object. This is in UTM format and so has UTM information 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create(vector<float> x, vector<float> y, int thisUTMzone)
{
  int n_x_points = int(x.size());
  int n_y_points = int(y.size());
  vector<double> x_double;
  vector<double> y_double;
  if (n_x_points != n_y_points)
  {
    cout << "X and Y vectors not the same size! No data loaded. " << endl;
  }
  for(int i = 0; i<n_x_points; i++)
  {
    x_double.push_back( double(x[i]));
    y_double.push_back( double(y[i]));
  }
  
  create(x_double, y_double, thisUTMzone);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry object. This is in UTM format and so has UTM information
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create(vector<double> x, vector<double> y, int thisUTMzone, bool ThisisNorth)
{
  // This reads in x and y vectors. 
  // It judges whether or not these are UTM coordinates by checking if
  // any of the data elements are over the maximum or minimum of latitude
  // or longitude
  int n_x_points = int(x.size());
  int n_y_points = int(y.size());
  
  if (n_x_points != n_y_points)
  {
    cout << "X and Y vectors not the same size! No data loaded. " << endl;
  }
  else if (n_x_points == 0)
  {
    cout << "Your data vectors are empty! " << endl;
  }
  else
  {
    UTMPoints_Easting = x; 
    UTMPoints_Northing = y;
    UTMZone = thisUTMzone;
    isNorth = ThisisNorth; 
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry object. This is in UTM format and so has UTM information 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::create(vector<float> x, vector<float> y, int thisUTMzone, bool ThisisNorth)
{
  int n_x_points = int(x.size());
  int n_y_points = int(y.size());
  vector<double> x_double;
  vector<double> y_double;
  if (n_x_points != n_y_points)
  {
    cout << "X and Y vectors not the same size! No data loaded. " << endl;
  }
  for(int i = 0; i<n_x_points; i++)
  {
    x_double.push_back( double(x[i]));
    y_double.push_back( double(y[i]));
  }
  
  create(x_double, y_double, thisUTMzone,ThisisNorth);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This converts points from Lat/Long to UTM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::convert_points_to_UTM()
{

  if( WGS84Points_latitude.size() == 0)
  {
    cout << "Trying to convert from Latitude and Longitude but you don't have any data." << endl;
  }
  else
  {
    LSDCoordinateConverterLLandUTM Converter;
  
    // set the default ellipsoid to WGS84
    int eId = 22;
  
    double Northing,Easting;
    double Lat,Long;
  
    int thisZone = 1;     // This gets replaced in the LLtoUTM operation
  
    vector<double> new_UTM_Northing;
    vector<double> new_UTM_Easting;
    // get UTMZone of the first coordinate and set that as the object zone
    Converter.LLtoUTM(eId, WGS84Points_latitude[0], WGS84Points_longitude[0],  Northing, Easting, thisZone);
    
    // get the zone and the isNorth parameters
    UTMZone = thisZone;
    if(WGS84Points_latitude[0] > 0)
    {
      isNorth = true;
    }
    else
    {
      isNorth = false;
    }

    int n_nodes = int(WGS84Points_longitude.size());
  
    for(int i = 0; i<n_nodes; i++)
    {
      Lat =WGS84Points_latitude[i]; 
      Long =WGS84Points_longitude[i];
     
     // We use force zone here since all points will be forced into the UTM zone 
     // of the first data element. 
      Converter.LLtoUTM_ForceZone(eId, Lat, Long,  Northing, Easting, UTMZone);
    
      new_UTM_Northing.push_back(Northing);
      new_UTM_Easting.push_back(Easting);
    }
    UTMPoints_Easting = new_UTM_Easting;
    UTMPoints_Northing = new_UTM_Northing;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDGeometry from UTM to Lat/Long
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::convert_points_to_LatLong()
{

  
  if( UTMPoints_Northing.size() == 0)
  {
    cout << "Trying to convert from UTM but you don't have any data." << endl;
  }
  else
  {
    LSDCoordinateConverterLLandUTM Converter;
  
    // set the default ellipsoid to WGS84
    int eId = 22;
  
    double Northing,Easting;
    double Lat,Long;
  
    vector<double> new_WGS84Points_latitude;
    vector<double> new_WGS84Points_longitude;

    int n_nodes = int(UTMPoints_Easting.size());
  
    for(int i = 0; i<n_nodes; i++)
    {
      Northing =UTMPoints_Northing[i]; 
      Easting =UTMPoints_Easting[i];
     
      Converter.UTMtoLL(eId, Northing, Easting, UTMZone, isNorth, Lat,Long);
    
      new_WGS84Points_latitude.push_back(Lat);
      new_WGS84Points_longitude.push_back(Long);
    }
    WGS84Points_latitude = new_WGS84Points_latitude;
    WGS84Points_longitude = new_WGS84Points_longitude;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the underlying point data to csv
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::print_points_to_csv(string path, string file_prefix)
{

  string fname = path+file_prefix+".csv";

  int n_UTM_nodes = int(UTMPoints_Northing.size());
  int n_WGS_nodes = int(WGS84Points_latitude.size());
  
  bool have_UTM = false;
  bool have_WGS = false;
  
  if (n_UTM_nodes > 0)
  {
    have_UTM = true;
  }
  if (n_WGS_nodes > 0)
  {
    have_WGS = true;
  }
  
  if (not have_UTM && not have_WGS)
  {
    cout << "Trying to print points but there is no data!" << endl;
  }
  else
  {
    if (have_UTM && not have_WGS)
    {
      convert_points_to_LatLong();
      n_WGS_nodes = int(WGS84Points_latitude.size());
    }
    if (not have_UTM && have_WGS)
    {
      convert_points_to_UTM();
      n_UTM_nodes = int(UTMPoints_Northing.size());
    }
    
    ofstream csv_out;
    csv_out.open(fname.c_str());
    csv_out << "latitude,longitude,Northing,Easting" << endl;
    for (int i = 0; i< n_UTM_nodes; i++)
    {
      csv_out.precision(7);
      csv_out << WGS84Points_latitude[i] << "," << WGS84Points_longitude[i] << ",";
      csv_out.precision(9);
      csv_out << UTMPoints_Northing[i] << "," << UTMPoints_Easting[i] << endl;
    }
    csv_out.close();
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets the row and column of all points in the raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::find_row_and_col_of_points(LSDRasterInfo& RI, vector<int>& RowOfNodes, vector<int>& ColOfNodes)
{
  // empty vectors for adding the rows and columns
  vector<int> row_vec;
  vector<int> col_vec;

  // first check to see if there is any UTM data
  if (UTMPoints_Northing.size() == 0)
  {
    // If you don't have UTM points, convert WGS to UTM
    convert_points_to_UTM();
  }
  
  // Get the X and Y minimums
  float XMinimum = RI.get_XMinimum();
  float YMinimum = RI.get_YMinimum();
  int NoDataValue = int(RI.get_NoDataValue());
  float DataResolution = RI.get_DataResolution();
  int NRows = RI.get_NRows();
  int NCols = RI.get_NCols();

  int this_row;
  int this_col;
  
  // Shift origin to that of dataset
  float X_coordinate_shifted_origin;
  float Y_coordinate_shifted_origin;
  
  // the actual coordinates
  float X_coordinate;
  float Y_coordinate;
  
  // now loop through nodes
  int N_nodes =  int(UTMPoints_Northing.size());
  for(int i = 0; i<N_nodes; i++)
  {
    this_row = NoDataValue;
    this_col = NoDataValue;

    X_coordinate = UTMPoints_Easting[i];
    Y_coordinate = UTMPoints_Northing[i];
  
    // Shift origin to that of dataset
    X_coordinate_shifted_origin = X_coordinate - XMinimum;
    Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

    // Get row and column of point
    int col_point = int(X_coordinate_shifted_origin/DataResolution);
    int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

    //cout << "Getting row and col, " << row_point << " " << col_point << endl;
    if(col_point > 0 && col_point < NCols-1)
    {
      this_col = col_point;
    }
    if(row_point > 0 && row_point < NRows -1)
    {
      this_row = row_point;
    }
    
    row_vec.push_back(this_row);
    col_vec.push_back(this_col);
  }

  RowOfNodes = row_vec;
  ColOfNodes = col_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets the row and column of a single point. 
// It includes out of bounds points (i.e., with rows and columsn either less than zero
//  or greater than NRows or NCols)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::find_row_and_col_of_point_inc_out_of_bounds(LSDRasterInfo& RI, 
                      int point_index, int& RowOfNode, int& ColOfNode, bool IsOutOfBounds)
{

  // Get the X and Y minimums
  float XMinimum = RI.get_XMinimum();
  float YMinimum = RI.get_YMinimum();
  int NoDataValue = int(RI.get_NoDataValue());
  float DataResolution = RI.get_DataResolution();
  int NRows = RI.get_NRows();
  int NCols = RI.get_NCols();

  int this_row = NoDataValue;
  int this_col = NoDataValue;
  bool OutOfBounds = true;

  // Shift origin to that of dataset
  float X_coordinate_shifted_origin;
  float Y_coordinate_shifted_origin;
  
  // the actual coordinates
  float X_coordinate;
  float Y_coordinate;

  // first check to see if there is any UTM data
  if (UTMPoints_Northing.size() == 0)
  {
    // If you don't have UTM points, convert WGS to UTM
    convert_points_to_UTM();
  }
  
  // Test to see if node index is in range
  int N_nodes =  int(UTMPoints_Northing.size());
  if (point_index < 0 || point_index >= N_nodes)
  {
    cout << "Your node index is out of the bounds of the point data." << endl;
  
  }
  else
  {
  
    X_coordinate = UTMPoints_Easting[point_index];
    Y_coordinate = UTMPoints_Northing[point_index];
  
    // Shift origin to that of dataset
    X_coordinate_shifted_origin = X_coordinate - XMinimum;
    Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

    // Get row and column of point
    int col_point = int(X_coordinate_shifted_origin/DataResolution);
    int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

    //cout << "Getting row and col, " << row_point << " " << col_point << endl;
    if(col_point > 0 && col_point < NCols-1)
    {
      OutOfBounds = false;
    }
    if(row_point > 0 && row_point < NRows -1)
    {
      OutOfBounds = false;
    }
    
    this_row = row_point;
    this_col = col_point;
  }

  RowOfNode = this_row;
  ColOfNode = this_col;
  IsOutOfBounds = OutOfBounds;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks if there is UTM data and if not updates it
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDGeometry::check_and_update_UTM()
{
  int n_nodes = int(UTMPoints_Northing.size());
  if (n_nodes == 0)
  {
    convert_points_to_UTM();
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// These functions get the maximum and minimum values for Northing and easting
// In general these are for scaling svg output. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDGeometry::get_max_UTM_Northing()
{
  check_and_update_UTM();
  int n_nodes = int(UTMPoints_Northing.size());
  double max_Northing = UTMPoints_Northing[0];
  
  for(int n = 0; n<n_nodes; n++)
  {
    if(UTMPoints_Northing[n] > max_Northing)
    {
      max_Northing = UTMPoints_Northing[n]; 
    }
  }
  return max_Northing;

}
double LSDGeometry::get_min_UTM_Northing()
{
  check_and_update_UTM();
  int n_nodes = int(UTMPoints_Northing.size());
  double min_Northing = UTMPoints_Northing[0];
  
  for(int n = 0; n<n_nodes; n++)
  {
    if(UTMPoints_Northing[n] < min_Northing)
    {
      min_Northing = UTMPoints_Northing[n]; 
    }
  }
  return min_Northing;
}
double LSDGeometry::get_max_UTM_Easting()
{
  check_and_update_UTM();
  int n_nodes = int(UTMPoints_Easting.size());
  double max_Easting = UTMPoints_Easting[0];
  
  for(int n = 0; n<n_nodes; n++)
  {
    if(UTMPoints_Easting[n] > max_Easting)
    {
      max_Easting = UTMPoints_Easting[n]; 
    }
  }
  return max_Easting;

}
double LSDGeometry::get_min_UTM_Easting()
{
  check_and_update_UTM();
  int n_nodes = int(UTMPoints_Easting.size());
  double min_Easting = UTMPoints_Easting[0];
  
  for(int n = 0; n<n_nodes; n++)
  {
    if(UTMPoints_Easting[n] < min_Easting)
    {
      min_Easting = UTMPoints_Easting[n]; 
    }
  }
  return min_Easting;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks to see if the NodeOrder exists and if not simply links the poitns in order
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPolyline::make_simple_polyline()
{
  vector<int> new_node_order;
  
  int n_order = (node_order.size());
  if(n_order == 0)
  {
    int n_nodes;
    n_nodes = int(UTMPoints_Easting.size());
    if( n_nodes == 0)
    {
      n_nodes = int(WGS84Points_latitude.size());
    }
    
    for(int i = 0; i<n_nodes; i++)
    {
      new_node_order.push_back(i);
    }
    node_order =  new_node_order;
    
  }
  else
  {
    cout << "This polyline already has a node order: I won't do anything" << endl;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks to see if the NodeOrder exists and if not simply links the points in order
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPolyline::force_simple_polyline()
{
  vector<int> new_node_order;
  
  int n_nodes;
  n_nodes = int(UTMPoints_Easting.size());
  if( n_nodes == 0)
  {
    n_nodes = int(WGS84Points_latitude.size());
  }
    
  for(int i = 0; i<n_nodes; i++)
  {
    new_node_order.push_back(i);
  }
  node_order = new_node_order;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDPolyline::get_affected_pixels_in_line_segment(LSDRasterInfo& RI,
                                                      vector<int>& affected_rows, vector<int>& affected_cols, 
                                                      int start_node, int end_node)
{
  // get some info from the RasterInfo object
  float XMinimum = RI.get_XMinimum();
  float YMinimum = RI.get_YMinimum();
  int NoDataValue = int(RI.get_NoDataValue());
  float DataResolution = RI.get_DataResolution();
  int NRows = RI.get_NRows();
  int NCols = RI.get_NCols();


  // make empty vectors that will contain the nodes
  vector<int> node_row;
  vector<int> node_col;
  
  // check to make sure there are UTM coords
  check_and_update_UTM();
  
  int n_nodes = int(UTMPoints_Easting.size());
  
  bool valid_nodes = true;
  
  // get the starting node location in UTM
  if (start_node < 0 || start_node >= n_nodes)
  {
    cout << "Start node does not have a valid index" << endl;
    valid_nodes = false;
  }
  if (end_node < 0 || end_node >= n_nodes)
  {
    cout << "End node does not have a valid index" << endl;
    valid_nodes = false;
  }
  
  // The start and end nodes are valid, enter the analysis
  if (valid_nodes)
  {
    double UTM_East_start = UTMPoints_Easting[start_node];
    double UTM_North_start = UTMPoints_Northing[start_node]; 

    double UTM_East_end = UTMPoints_Easting[end_node];
    double UTM_North_end = UTMPoints_Northing[end_node]; 
    
    // get the row and column of both the start and end node
    int start_row,start_col;
    int end_row,end_col;
    bool OOB_start, OOB_end;
    find_row_and_col_of_point_inc_out_of_bounds(RI, start_node, start_row, start_col, OOB_start);
    find_row_and_col_of_point_inc_out_of_bounds(RI, end_node, end_row, end_col, OOB_end);
    
    
    int current_col;
    int current_row;
    
    bool IsIncreasingCol;
    bool IsIncreasingRow;
    
    
    double denominator = UTM_East_end-UTM_East_start;
    double numerator = UTM_North_end-UTM_North_start;
    if ( (end_row-start_row) == 0)
    {
      // loop through all the nodes between the start and end
      if(start_col < 0)
      {
        if (end_col < 0)
        {
          current_col =NoDataValue;
        }
        else
        {
          current_col = 0;
        }
        
        
        
        
      }
    }
    else if (numerator == 0)
    {
      // This means the row always stays the same
    }
    else
    {
      double slope = numerator/denominator;
    }

  }

}



#endif