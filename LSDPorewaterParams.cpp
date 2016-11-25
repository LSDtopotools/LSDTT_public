//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDPorewaterParams
// Land Surface Dynamics PorewterParams object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object interfaces with teh porewater column object
//  In a landsacpe each pixel will have its own pore ressures but 
//  the parameters will be constant (or have similar statistical properties)
//  across a landscape. This structure tries to minimize memory requirements
//
// Developed by:
//  Simon M. Mudd
//  Stuart W.D. Grieve
//
// Copyright (C) 2016 Simon M. Mudd 2013 6
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


/** @file LSDPorewaterColumn.hpp
@author Simon M. Mudd, University of Edinburgh
@author Stuart W. D. Grieve, University of Edinburgh

**/



#ifndef LSDPorewaterParams_CPP
#define LSDPorewaterParams_CPP

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "LSDPorewaterParams.hpp"
#include "LSDParameterParser.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Empty create function
// Starts with some defaults. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::create()
{
  D_0 = 0.00001;
  d = 2;
  alpha = 0.1;
  Iz_over_K_steady = 0.2;
  
  vector<float> this_z;
  for(int i = 0; i<31; i++)
  {
    this_z.push_back(float(i)*0.1);
  }
  Depths = this_z;

  calculate_beta();
  calculate_D_hat();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function reads from a parameter file. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::create(string paramfile_path, string paramfile_name)
{
  LSDParameterParser LSDPP(paramfile_path, paramfile_name);
  
  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // set default float parameters
  float_default_map["D_0"] = 0.00001;
  float_default_map["K_sat"] = 5e-8;
  float_default_map["d"] = 2;
  float_default_map["alpha"] = 0.1;
  float_default_map["depth_spacing"] = 0.1;
  float_default_map["Iz_over_K_steady"] = 0.2;
  
  // set default in parameter
  int_default_map["n_depths"] = 31;
  
  bool_default_map["use_depth_vector"] = false;
  
  string depth_vector_key = "depth_vector";
  string_default_map["depth_vector"] = "";
  
  // Get the parameters
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  //map<string,string> this_string_map = LSDPP.get_string_parameters();
  
  D_0 = this_float_map["D_0"];
  d = this_float_map["d"];
  alpha = this_float_map["alpha"];
  Iz_over_K_steady = this_float_map["Iz_over_K_steady"];
  K_sat = this_float_map["K_sat"];

  // parameters for getting the depths
  float depth_spacing = this_float_map["depth_spacing"];
  int n_depths = this_int_map["n_depths"];
  
  vector<float> these_extracted_depths = LSDPP.parse_float_vector(depth_vector_key);
  vector<float> these_calculated_depths;
  for(int i = 0; i<n_depths; i++)
  {
    these_calculated_depths.push_back(float(i)*depth_spacing);
  }
  
  if (these_extracted_depths.size() == 0)
  {
    these_extracted_depths = these_calculated_depths;
  }
  
  // If we are going to use the depth vector, read it from the 
  // comma separated list
  if( this_bool_map["use_depth_vector"]  )
  {
    Depths = these_extracted_depths;
  }
  else
  {
    Depths = these_calculated_depths;
  }
  
  calculate_beta();
  calculate_D_hat();
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates beta
// Comes from Iverson's eq 27
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::calculate_beta()
{

  beta = cos(alpha)*cos(alpha)-Iz_over_K_steady;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates D_hat
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::calculate_D_hat()
{

  D_hat = 4*D_0*(cos(alpha)*cos(alpha));
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This uses the parameters to get a steady state pressure profile
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterParams::calculate_steady_psi()
{
  vector<float> Psi;
  for(int i = 0 ; i < int(Depths.size()); i++ )
  {
    Psi.push_back(beta*(Depths[i]-d));
  }
  return Psi;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses a rainfall file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_rainfall_file(string path, string filename, vector<float>& intensities)
{
  
  string fname = FixPath(path)+ filename;
  

  
  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  float this_rain; 
  vector<float> rain_vec;
  
  vector<string> HeaderInfo = ReadCSVHeader(path, filename);
  
  // now find the rainfall column
  string rain_string = "rainfall_rate";
  string this_string;
  int rain_column = 0; 
  for(int i = 0; i< int(HeaderInfo.size()); i++)
  {
    cout << "Header["<<i<<"]: " << HeaderInfo[i] << endl;
    this_string = HeaderInfo[i];
    if (this_string.compare(rain_string) == 0)
    {
      cout << "I found the rain rate, it is column " << i << endl;
      rain_column = i;
    }
  }
  
  // now we work through the file. 
  // make sure the filename works
  ifstream ifs(fname.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv cosmo data file, but the file" << filename
         << "doesn't exist; LINE 245 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  // get the first line  and discard
  getline(ifs, line_from_file);

  // now loop through the rest of the lines, getting the data. 
  while( getline(ifs, line_from_file))
  {
    // reset the string vec
    this_string_vec = empty_string_vec;
    
    // create a stringstream
    stringstream ss(line_from_file);
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    // Now extract the rain rate
    this_rain =  atof(this_string_vec[rain_column].c_str());
    rain_vec.push_back(this_rain);
    
  }
  intensities = rain_vec;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses a rainfall file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_MIDAS_rainfall_file(string path, string filename,vector<int>& days, vector<float>& intensities )
{
  
  string fname = FixPath(path)+ filename;
  
  // These are 
  vector<int> day_since_1900;

  
  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  float this_rain; 
  int this_date;
  vector<float> rain_vec;
  
  vector<string> HeaderInfo = ReadCSVHeader(path, filename);
  
  // now find the data columns column
  string rain_string = "prcp_amt";
  string this_string;
  int rain_column = 0; 
  string date_string = "days_since_1900";
  int date_column = 0; 
  for(int i = 0; i< int(HeaderInfo.size()); i++)
  {
    cout << "Header["<<i<<"]: " << HeaderInfo[i] << endl;
    this_string = HeaderInfo[i];
    if (this_string.compare(rain_string) == 0)
    {
      cout << "I found the rain rate, it is column " << i << endl;
      rain_column = i;
    }
    if (this_string.compare(date_string) == 0)
    {
      cout << "I found the date, it is column " << i << endl;
      date_column = i;
    }
  }


  // now we work through the file. 
  // make sure the filename works
  ifstream ifs(fname.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv cosmo data file, but the file" << filename
         << "doesn't exist; LINE 245 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  // get the first line  and discard
  getline(ifs, line_from_file);

  // now loop through the rest of the lines, getting the data. 
  while( getline(ifs, line_from_file))
  {
    // reset the string vec
    this_string_vec = empty_string_vec;
    
    // create a stringstream
    stringstream ss(line_from_file);
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    // Now extract the rain rate
    this_rain =  atof(this_string_vec[rain_column].c_str());
    this_date = atoi(this_string_vec[date_column].c_str());
    rain_vec.push_back(this_rain);
    day_since_1900.push_back(this_date);
  }
  intensities = rain_vec;
  days = day_since_1900;
  
  vector<float> new_intensities;
  vector<float> durations;
  parse_MIDAS_duration_intensities(days, intensities, durations);
  
  // bug checking
  for(int i = 0; i<int(durations.size()); i++)
  {
    cout << "duration["<<i<<"]: " << durations[i] << " intensity: "  << intensities[i] << endl;
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses MIDAS intensity and day data to give multiday durations if
// there are empty data slots or repeated rainfall amounts
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_MIDAS_duration_intensities(vector<int>& days, vector<float>& intensities, vector<float>& durations)
{
  vector<float> new_intensities;
  vector<float> new_durations;

  // we first do a sweep getting intermediate days
  vector<float> intermediate_intensities;
  intermediate_intensities.push_back(intensities[0]);
  int n_records = int(intensities.size());
  for(int i = 1; i<n_records; i++)
  {
    // If we have incremented by a single day, just push back the data
    if(days[i] == days[i-1]+1)
    {
      intermediate_intensities.push_back(intensities[i]);
    }
    else
    {
      // here we need to push back zeros for the missing days
      int day_dif = days[i]-days[i-1]-1;
      for(int j = 0; j<day_dif; j++)
      {
        intermediate_intensities.push_back(0.0);
      }
      // now push back the next intensity
      intermediate_intensities.push_back(intensities[i]);
    }
  }
  
  // now print all these for bug checking
  int n_int_intensities = int(intermediate_intensities.size());
  int this_day = days[0];
  //for(int i = 0; i<n_int_intensities; i++)
  //{
  //  cout << "Day: " << this_day << " i: " << intermediate_intensities[i] << endl;
  //  this_day++;
  //}
  
  // Now go back over the intensities and set the durations and intensities
  int this_consecutive = 0;
  float last_intensity = intermediate_intensities[0];
  for(int i = 1; i<n_int_intensities; i++)
  {
    // add a consecutive day
    this_consecutive++;
    
    // check to see if the intensity has changed
    // if the intensity has not changed, you don't do anything.
    if (intermediate_intensities[i] != last_intensity)
    {
      // the intensity has changed. Add the last one to the record and update 
      // the last changed intensity
      new_intensities.push_back(last_intensity);
      new_durations.push_back(float(this_consecutive));
      last_intensity = intermediate_intensities[i];
      this_consecutive = 0;
    }
  }
  // now you need to deal with the last point. There is always data there
  if(this_consecutive == 0)
  {
    // This means that the last data point was a new intensity, so there is only
    // one day left
    new_intensities.push_back(last_intensity);
    new_durations.push_back(1.0);
  }
  else
  {
    // this means that the last data point was repeated, so we need to use
    // different logic
    new_intensities.push_back(last_intensity);
    new_durations.push_back(float(this_consecutive));
  }
  
  
  // check to see if it worked
  //int n_intensities = int(new_intensities.size());
  //for(int i = 0; i<n_intensities; i++)
  //{
  //  cout << "Duration: " << new_durations[i] << " i: " << new_intensities[i] << endl;
  //}
  
  // now reset the vectors
  intensities = new_intensities;
  durations = new_durations; 
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis prints the parameters to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterParams::print_parameters_to_screen()
{
  cout << "Your parameters are: " << endl;
  cout << "D_0: " << D_0 << endl;   
  cout << "D_hat: " << D_hat << endl;
  cout << "alpha: " << alpha << endl;
  cout << "d: " << d << endl;
  cout << "beta: " << beta << endl;
  cout << "Iz_over_K_steady: " << Iz_over_K_steady << endl;
  cout << "And the depths are: " << endl;
  for (int i = 0; i< int(Depths.size()); i++)
  {
    cout << Depths[i] << ",";
  }
  cout << endl;
}

#endif