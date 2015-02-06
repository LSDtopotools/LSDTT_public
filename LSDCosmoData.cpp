//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDCosmoData.cpp
// Land Surface Dynamics CosmoData
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for keeping track of cosmogenic data
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
#include <fstream>
#include <math.h>
#include <iostream>
#include <map>
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDCosmoData.hpp"
#include "LSDRaster.hpp"
using namespace std;

#ifndef LSDCosmoData_CPP
#define LSDCosmoData_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The default create function. Doesn't do anything
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::create()
{}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function that actually loads the file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::create(string filename)
{
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: The file" << filename
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  // now populate the standardisation maps
  standards_Be10["07KNSTD"] = 1.0;
  standards_Be10["KNSTD"] = 0.9042;
  standards_Be10["NIST_Certified"] = 1.0425;
  standards_Be10["LLNL31000"] = 0.8761;
  standards_Be10["LLNL10000"] = 0.9042;
  standards_Be10["LLNL3000"] = 0.8644;
  standards_Be10["LLNL1000"] = 0.9313;
  standards_Be10["LLNL300"] = 0.8562;
  standards_Be10["NIST_30000"] = 0.9313;
  standards_Be10["NIST_30200"] = 0.9251;
  standards_Be10["NIST_30300"] = 0.9221;
  standards_Be10["NIST_30600"] = 0.9130;
  standards_Be10["NIST_27900"] = 1.0;
  standards_Be10["S555"] = 0.9124;
  standards_Be10["S2007"] = 0.9124;
  standards_Be10["BEST433"] = 0.9124;
  standards_Be10["BEST433N"] = 1.0;
  standards_Be10["S555N"] = 1.0;
  standards_Be10["S2007N"] = 1.0;

  standards_Al26["KNSTD"] = 1.0;
  standards_Al26["ZAL94"] = 0.9134;
  standards_Al26["SMAL11"] = 1.021;
  standards_Al26["0"] = 1.0;
  standards_Al26["ZAL94N"] = 1.0;
  standards_Al26["ASTER"] = 1.021;
  standards_Al26["Z92-0222"] = 1.0;

  // now load the file. 
  double this_lat;
  double this_long;
  double this_conc;
  double this_uncert;
  
  string this_sample_name;
  string this_standard;
  string this_nuclide;
  
  while( ifs >> this_sample_name >> this_nuclide  >> this_lat >> this_long 
             >> this_conc >> this_uncert >> this_standard)
  {
    sample_name.push_back(this_sample_name);
    latitude.push_back(this_lat);
    longitude.push_back(this_long);
    nuclide.push_back(this_nuclide);
    Concentration_unstandardised.push_back(this_conc);
    Concentration_uncertainty_unstandardised.push_back(this_uncert);
    standardisation.push_back(this_standard);
  }
  
  // now loop through the data, getting the standardised concentrations
  N_samples = int(sample_name.size());
  
  for(int i = 0; i<N_samples; i++)
  {
    // read the nuclide
    if(nuclide[i] == "Be10")
    {
      // check to see if standard exists
      if(standards_Be10.find(standardisation[i]) == standards_Be10.end())
      {
        cout << "Standardisation not found, assuming 07KNSTD" << endl;
        standardisation[i] = "07KNSTD";
      }
      
      // adjust the 10Be concentration
      Concentration.push_back(Concentration_unstandardised[i]*
                                standards_Be10[ standardisation[i] ]);
      Concentration_uncertainty.push_back(Concentration_uncertainty_unstandardised[i]*
                                standards_Be10[ standardisation[i] ]); 
    }
    else if (nuclide[i] == "Al26")
    {
      // check to see if standard exists
      if(standards_Al26.find(standardisation[i]) == standards_Al26.end())
      {
        cout << "Standardisation not found, assuming KNSTD" << endl;
        standardisation[i] = "KNSTD";
      }
      
      // adjust the 26Al concentration
      Concentration.push_back(Concentration_unstandardised[i]*
                                standards_Al26[ standardisation[i] ]);
      Concentration_uncertainty.push_back(Concentration_uncertainty_unstandardised[i]*
                                standards_Al26[ standardisation[i] ]);   
    }
    else
    {
      cout << "You haven't selected a valid nuclide, choices are Be10 and Al26" << endl;
      cout << "You chose: " << nuclide[i] << ", defaulting to Be10 and 07KNSTD" << endl;
      nuclide[i] = "Be10";
      Concentration.push_back(Concentration_unstandardised[i]);
      Concentration_uncertainty.push_back(Concentration_uncertainty_unstandardised[i]);
    }
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function converts data to UTM coordinates. You have to tell it 
// the UTM zone
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::convert_to_UTM(int UTM_zone)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;
  
  // set up some temporary vectors
  vector<double> this_UTMN(N_samples,0);
  vector<double> this_UTME(N_samples,0);
  
  double this_Northing;
  double this_Easting;
  
  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    Converter.LLtoUTM(eId, latitude[i], longitude[i], 
                      this_Northing, this_Easting, UTM_zone);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
  }
  
  UTM_easting = this_UTME;
  UTM_northing = this_UTMN;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function converts data to UTM coordinates. It determines the UTM from
// a raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::convert_to_UTM(LSDRaster& Raster)
{
  // get the UTM zone from the raster
  int UTM_zone;
  bool is_North;
  Raster.get_UTM_information(UTM_zone, is_North);
  
  // convert the data
  convert_to_UTM(UTM_zone);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#endif