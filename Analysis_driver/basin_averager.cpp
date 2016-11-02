//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_mapping_tool
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate 
// different kinds of chi analysis
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDBasin.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDParameterParser.hpp"

int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the basin averaging tool!                ||" << endl;
    cout << "|| This program has a number of options to compute     ||" << endl;
    cout << "|| basin averaged statistics                           ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Simon M. Mudd and Stuart W.D. Grieve                ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./basin_averaging_tool.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_BasinAvg.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << " http://lsdtopotools.github.io/LSDTT_book/#_chi_analysis_part_3_getting_chi_gradients_for_the_entire_landscape" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  LSDParameterParser LSDPP(path_name,f_name);

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  
  int_default_map["basin_order"] = 2;
  float_default_map["min_slope_for_fill"] = 0.0001;
  
  // Use the parameter parser to get the maps of the parameters required for the 
  // analysis
  map<string,float> these_float_parameters = LSDPP.set_float_parameters(float_default_map);
  map<string,int> these_int_parameters = LSDPP.set_int_parameters(int_default_map);

}
