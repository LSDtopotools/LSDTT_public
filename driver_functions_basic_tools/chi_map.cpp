//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_map.cpp
// This function calcualtes chi across a landscape. Run without
// arguments for instructions
//
// Developed by:
//  Simon M. Mudd
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#include <iostream>
#include <sstream>
#include <cstdlib>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDFlowInfo.hpp"
using namespace std;


int main(int argc, char *argv[])
{

  if (argc!=6)
  {
    cout << "===========================================================" << endl;
    cout << "|| Welcome to the chi_map tool!                          ||" << endl;
    cout << "|| I am used to calculate the chi coordinate from a DEM. ||" << endl;
    cout << "===========================================================" << endl;
    cout << "This program requires six inputs: " << endl << endl;
    cout << "* First the path to the DEM." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "   (Either \\ or / depending on your operating system.)" << endl << endl;
    cout << "* Second the prefix of the DEM. " << endl;
    cout << "   For example, if the DEM is called Spain.bil the filename is Spain." << endl << endl;
    cout << "* Third, the DEM format. Must be either bil, flt or asc. " << endl;
    cout << "   If you write something the computer doesn't understand it will default to bil" << endl << endl;
    cout << "* Fourth, a threshold area for a channel (in m^2). " << endl << endl;
    cout << "* Fifth, a flag to tell the code to print out the data as a csv (1 for true)." << endl;
    cout << "   WARNING you should only do this if you have a small DEM or " << endl
         << "   set a reasonably large drainage area threshold or you'll get a huge file" << endl;
    cout << "============================================================" << endl;
    cout << "An example call is: " << endl;
    cout << "./chi_map.exe ~home/basins/Chile/ Chile_test bil 50 1" << endl;
    exit(EXIT_SUCCESS);
  }

  string pathname = argv[1];
  string lchar = pathname.substr(pathname.length()-2,1);
  string slash = "/";
  cout << "lchar is " << lchar << " and slash is " << slash << endl;
      
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;  
    pathname = pathname+slash;
  } 
  cout << "The pathname is: " << pathname << endl;


  string dem_name = argv[2];
  cout << "The DEM filename is: " << dem_name << endl;
  string full_DEM_name = pathname+dem_name;
  cout << "The full path is: " << full_DEM_name << endl;
    
  string DEM_ext = argv[3];

  float threshold_area = atof(argv[4]);
  cout << "The threshold area is: " << threshold_area << endl;
  
  int temp_csv_flag = atoi(argv[5]);
  bool write_csv = false;
  if(temp_csv_flag == 1)
  {
    cout << endl << endl << endl << "---------" << endl << "I'll print a csv" << endl;
    write_csv = true;
  }

  // load the raster
  LSDRaster new_raster(full_DEM_name,DEM_ext);
  
  // remove the sea (seems to be required if gdal is used in places with nodata)
  new_raster.remove_seas();
  
  // get the filled raster
  float min_slope = 0.0001;
  LSDRaster filled_raster = new_raster.fill(min_slope);
  
  // get the flow info 
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "Periodic";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "Periodic";
		
  cout << "Filled raster, getting flow info" << endl;
  LSDFlowInfo FI(boundary_conditions, filled_raster);
  cout << "Got flow info" << endl;
	
  float m_over_n = 0.5;
  float A_0 = 1000;
	
  //vector<int> bl_nodes(2);
  //bl_nodes[0] = 156860657;
  //bl_nodes[2] = 6177171;
  LSDRaster chi_map = FI.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0, 
                                                threshold_area);
  cout << "got chi map" << endl;                                              
  //LSDRaster chi_map = FI.get_upslope_chi_from_multiple_starting_nodes(bl_nodes,
  //                                              m_over_n, A_0, 
  //                                              threshold_area);
  
  string hs_name = "_HS";
  hs_name = full_DEM_name+hs_name;
  LSDRaster HS = filled_raster.hillshade(45,315,1); 
  HS.write_raster(hs_name,DEM_ext);


 
  float hundredtimesmn = m_over_n*100;
  int m_n = int(hundredtimesmn);
  int A_0int = int(A_0);
  int TAint = int(threshold_area);

  string m_nstr = itoa(m_n);
  m_nstr = "_mn"+m_nstr;

  string A0str = itoa(A_0int);
  A0str = "_A0"+A0str;

  string TAstr = itoa(TAint);
  TAstr = "_TA"+TAstr;

  string chistr = "_chi"+m_nstr+A0str+TAstr;
  
  if(write_csv)
  {
    string pointsname = full_DEM_name+chistr+"_points";
    cout << "Writing csv, the prefix is: " << pointsname << endl;
    chi_map.FlattenToCSV(pointsname);
  }
  
                     
  string cm_fname = full_DEM_name+chistr;                                             
  chi_map.write_raster(cm_fname, DEM_ext);                                             
}
