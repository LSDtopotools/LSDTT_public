//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDPorewaterColumn
// Land Surface Dynamics PorewterColumn object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object calculates porewater pressure based on the Iverson 2000 WRR
// model
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



#ifndef LSDPorewaterColumn_CPP
#define LSDPorewaterColumn_CPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "LSDPorewaterColumn.hpp"
#include "LSDPorewaterParams.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Empty create function
// Starts with some defaults. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterColumn::create()
{
  row = 0;
  col = 0;
  node_index = 0;
  cout << "I am an empty LSDPorewaterColumn object." <<endl;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function just uses an initial Psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterColumn::create(vector<float> Initial_Psi)
{
  row = 0;
  col = 0;
  node_index = 0;
  Psi = Initial_Psi;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function takes a porewater parameter object and uses the
// steady infiltration rate to set psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterColumn::create(LSDPorewaterParams LSDPP)
{
  row = 0;
  col = 0;
  node_index = 0;
  Psi = LSDPP.calculate_steady_psi();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the Psi values to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterColumn::PrintPsiToScreen()
{
  for(int i = 0; i<int(Psi.size()); i++)
  {
    cout << "Psi["<<i<<"]: " << Psi[i] << endl;
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This caluclates the response function
// THis comes from iverson's equation 27e
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDPorewaterColumn::CalculateResponseFunction(float t_star)
{
    float R;
    
    float sqrt_term = sqrt(t_star/M_PI);
    float exp_term = exp(-1/t_star);
    
    float multiple_bit = sqrt_term*exp_term;
    
    if (t_star != 0)
    {
      R = multiple_bit- erfcf(1/ (sqrt(t_star)));
    }
    else   // If t_star is 0, then 1/sqrt(t_star) is infinity, meaning erfc is 0)
    {
      R = 0;
    }

    return R;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This caluclates the Psi value based on iverson's equation 27
// Has only the transient component of psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::CalcualtePsiTransient(LSDPorewaterParams& LSDPP, float Iz_over_Kz, 
                           float t_star, float T_star)
{
    vector<float> Depths = LSDPP.get_Depths();
    vector<float> transient_Psi(Depths.size());
    
    float R;
    if (t_star < T_star)
    {
      R = CalculateResponseFunction(t_star);
    }
    else
    {
      R = CalculateResponseFunction(t_star-T_star);
    }
    
    // This solves the equation, based on the response function (R_fn),
    // which is equation 27e
    for (int i = 0; i< int(Depths.size()) ; i++)
    {
      transient_Psi[i] = Depths[i]*Iz_over_Kz*R;
    }

    return transient_Psi;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Compute psi from equation 27a and b, but using dimensional time
// A bit slow since I haven't vectorised the calculations.
// Only calculates the transient component of psi for use with 
// time series of rainfall
// times need to be in seconds
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::CalculatePsiDimensionalTimeTransient(LSDPorewaterParams& LSDPP, float t, float T, float Iz_over_Kz)
{

  vector<float> Depths = LSDPP.get_Depths();
  vector<float> transient_Psi(Depths.size());
  float t_star;
  float T_star;
  float zsquare;
  
  float D_hat = LSDPP.get_D_hat();
  float R;
  
  // loop through depths: each depth has a different t_star and T_star since
  // these depend on depth
  for(int i = 0; i< int(Depths.size()) ; i++)
  {
    // first get the nondimensional time. Note that according to
    // equations 27c,d the dimensionless time is a function of depth,
    // so each point below the surface has a different t_star and T_star
    zsquare = Depths[i]*Depths[i];
    t_star = t * D_hat / zsquare;
    T_star = T * D_hat / zsquare;
    
    if (t_star < T_star)
    {
      R = CalculateResponseFunction(t_star);
    }
    else
    {
      R = CalculateResponseFunction(t_star-T_star);
    }
    transient_Psi[i] =Depths[i]*Iz_over_Kz*R;
  }
  return transient_Psi;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This caluclates the Psi value based on iverson's equation 27
// Has only the transient component of psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterColumn::CalculatePsiFromTimeSeries(vector<float> durations, vector<float> intensities, 
                                LSDPorewaterParams& LSDPP, float t)
{
  // Get the steady state time
  vector<float> steady_psi = LSDPP.calculate_steady_psi();
  vector<float> cumulative_psi = steady_psi;

  // Now we try to construct the transient pressure. 
  // loop through the record getting cumulative times
  vector<float> starting_times;
  starting_times.push_back(0);
  float cumulative_time = 0;
  int count = 0; 
  bool end_count_found = false;
  int end_count;
  
  for (int i = 0; i< int(durations.size()); i++)
  {
    cumulative_time += durations[i];
    
    // the cumulative time is the time at the end of this timestep. 
    // if the cumulative time  is less than the time of simulation, 
    // then we need to acount for this pulse of rainfall        
    if (t < cumulative_time)
    {
      if (end_count_found == false)
      {
        end_count_found = true;
        end_count = count;
      }
        
      count++;
      starting_times.push_back(cumulative_time);
    }
  }
  
  // we don't need the last element
  starting_times.pop_back();

  // If we didn't find the end count it means the rainfall records have ended and we need
  // all of the data        
  if (end_count_found == false)
  {
    // The minus one is needed since we have counted past the end of the index
    end_count = count-1;
  }

  // okay, now get the transients from superposition 
  // First we need to figure out how many of these we will need
  float eff_t, this_intensity, this_duration;
  vector<float> this_transient_Psi;
  for(int i = 0; i< int(starting_times.size()); i++)
  {
    if(i<= end_count)
    {
      eff_t = t-starting_times[i];
      this_intensity = intensities[i];
      this_duration = durations[i];
      
      // get this steps Psi value
      this_transient_Psi = CalculatePsiDimensionalTimeTransient(LSDPP, eff_t, this_duration, this_intensity);

      // add this step's transient Psi values.
      for(int i = 0; i<int(cumulative_psi.size()); i++)
      {
        cumulative_psi[i]+=this_transient_Psi[i];
      }
    }
  }
}




#endif