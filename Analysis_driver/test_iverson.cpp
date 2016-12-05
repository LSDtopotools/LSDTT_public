//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// test_iverson
//
// This tests the porewater pressure routines
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDPorewaterParams.hpp"
#include "../LSDPorewaterColumn.hpp"
#include "../LSDStatsTools.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{

  string path_name;
  
  // load porewater parameter object
  LSDPorewaterParams LSDPP;
  
  //Test for correct input arguments
  if (nNumberofArgs == 1 )
  {
    LSDPorewaterParams TempPP;
    LSDPP = TempPP;
    
    path_name = "./";
  }
  else if (nNumberofArgs == 3 )
  {
    path_name = argv[1];
    string f_name = argv[2];

    // Make sure the path has an extension
    path_name = FixPath(path_name);

    LSDPorewaterParams TempPP(path_name,f_name);
    LSDPP = TempPP;
  }
  else
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the test_iverson program                 ||" << endl;
    cout << "|| This program is for testing the porewater column    ||" << endl;
    cout << "|| in LSDTopoTools. These are used for hydrology and   ||" << endl;
    cout << "|| slope stability calculations.                       ||" << endl;
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
    cout << "./test_iverson.out /LSDTopoTools/Topographic_projects/Test_data/ SoilColumn.param" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  // get a steady state column    
  //cout << "I just loaded my data!!!!"<< endl;
  //LSDPP.print_parameters_to_screen();
  
  LSDPorewaterColumn LSD_PC(LSDPP);
  LSDPP.print_parameters_to_screen();
  
  cout << "THe K_sat is: " << LSDPP.get_K_sat() << endl;
  
  
  cout << "Loading some rainfall data" << endl;
  string rainfall_fname = "MidasDoverSmall.csv";
  vector<float> intensities;
  vector<int> days;
  LSDPP.parse_MIDAS_rainfall_file(path_name, rainfall_fname,days,intensities);
  
  // lets make a duration intensity record to test the model
  vector<float> intensities_new;
  vector<float> durations_weeks;
  
  intensities_new.push_back(0.5);
  intensities_new.push_back(1.0);
  intensities_new.push_back(0.5);
  
  durations_weeks.push_back(5);
  durations_weeks.push_back(6);
  durations_weeks.push_back(4);
  
  vector<float> durations_seconds = LSDPP.weeks_to_seconds(durations_weeks);
  
  cout << "Durations: " <<durations_seconds[0] << "," << durations_seconds[1] << "," << durations_seconds[2] << endl;
  
  // now run the model
  float week_of_pressure = 11;
  float sec_of_pressure = LSDPP.weeks_to_seconds(week_of_pressure);
  cout << "Sec of pressure = " <<  sec_of_pressure << endl;
  
  LSD_PC.CalculatePsiFromTimeSeries(durations_seconds, intensities_new,LSDPP, sec_of_pressure);
  
  vector<float> FS = LSD_PC.FS(LSDPP);
  
  // now get a vector of weeks
  vector<float> times;
  float this_weeks;
  for(int i = 0; i<15; i++)
  {
    this_weeks = float(i);
    times.push_back(LSDPP.weeks_to_seconds(this_weeks));
  }
  
  float minimum_depth = 0.2;
  LSD_PC.ScanTimeseriesForFailure(durations_seconds, intensities_new,LSDPP, 
                                    minimum_depth, times);
  


}
