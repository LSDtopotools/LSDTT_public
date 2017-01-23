//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// CRONUS_emulator.cpp
// This program takes data from the CRN calculator and converts them to 
// denudation rates as calculated by the CRONUS2.2 calculator. 
// The code is a port of the CRONUS2.2 matlab scripts
//
// Developed by:
//  Simon M. Mudd, University of Edinburgh, School of GeoSciences
//  Stuart W.D. Grieve, University of Edinburgh, School of GeoSciences
//  Marie-Alice Harel, University of Edinburgh, School of GeoSciences
//  Martin D. Hurst, British Geological Survey
//
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "../LSDRaster.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDStrahlerLinks.hpp"
#include "../LSDBasin.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDCosmoData.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
  // the driver version
  string driver_version = "Driver_version: 0.02";

  cout << "================================================================" << endl;
  cout << "|| Welcome to the CRONUS emulator tool!                       ||" << endl;
  cout << "|| This tool takes output from the CRN basinwide calculator   ||" << endl;
  cout << "|| and then calculates the erosion rates using an emulator of ||" << endl;
  cout << "|| the CRONUS2.2 online calculator (it replicates the         ||" << endl;
  cout << "|| underlying matlab functions). It takes 2 arguments:        ||" << endl;
  cout << "|| 1. The path to the data files 2. The filename prefix       ||" << endl;
  cout << "================================================================" << endl;

  string path_name = argv[1];
  string f_name = argv[2];
  string filename = path_name+f_name+"_CRNResults.csv";
  string line_from_file;
  
  string file_out_name =  path_name+f_name+"_CRONUSEmulator.csv";
  
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  
  vector<string> sample_name;
  vector<string> sample_nuclide;
  vector<double> sample_Concentration_standardised;
  vector<double> sample_Concentration_uncertainty_standardised;
  vector<double> sample_latitude;
  vector<double> sample_longitude;
  vector<double> sample_eff_pressure;
  vector<double> sample_shielding;
  
  // the atmospheric data is in the folder with the driver_functions
  string path_to_atmospheric_data = "./";

  // print the muon comparison file
  LSDCRNParameters LSDCRNP;

  // get the information out of the file
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv cosmo data file, but the file" << filename
         << "doesn't exist; LINE 245 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    // discard the first four lines
    getline(ifs, line_from_file);
    getline(ifs, line_from_file);
    getline(ifs, line_from_file);
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
      
      //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
      if ( int(this_string_vec.size()) < 7)
      {
        cout << "Hey there, I am trying to load your cosmo data but you seem not to have" << endl;
        cout << "enough columns in your file. I am ignoring a line" << endl;
      }
      else
      {
        // now convert the data
        //cout << "Getting sample name: " <<  this_string_vec[0] << endl;
        
        // let the user know about offending underscores, and replace them
        string s = this_string_vec[1];
        string uscore = "_";
        size_t found = s.find(uscore);
        if (found!=std::string::npos)
        {
          cout << "I found an underscore in the sample name. Replacing with a dash." <<endl;
          replace( s.begin(), s.end(), '_', '-');
          cout << "New sample name is: " << s << endl;
        }

        sample_name.push_back( s );
        sample_latitude.push_back( atof( this_string_vec[3].c_str() ) );
        sample_longitude.push_back( atof(this_string_vec[4].c_str() ) );
        sample_nuclide.push_back( this_string_vec[2] );
        sample_Concentration_standardised.push_back( atof(this_string_vec[5].c_str() ) );
        sample_Concentration_uncertainty_standardised.push_back( atof(this_string_vec[6].c_str() ) );
        sample_eff_pressure.push_back( atof(this_string_vec[24].c_str() ) );
        sample_shielding.push_back( atof(this_string_vec[13].c_str() ) );
        
        
      }  // Ended logic for getting data
    }    // Ended logic for reading lines
  }      // Ended logic for file input
  
  // now loop through these samples
  ofstream CRONUS_out;
  CRONUS_out.open(file_out_name.c_str());
  
  int n_samp = sample_name.size();
  cout << "I found " << n_samp << " samples!" << endl;
  for(int samp = 0; samp < n_samp; samp++)
  {
    cout << "Sample details: Name: " << sample_name[samp] << " lat: " << sample_latitude[samp]
         << " Conc: " <<  sample_Concentration_standardised[samp] << " eff p:" << sample_eff_pressure[samp] << endl;
  
    // initiate a particle. We'll just repeatedly call this particle
    // for the sample.
    int startType = 0;
    double Xloc = 0;
    double Yloc = 0;
    double  startdLoc = 0.0;
    double  start_effdloc = 0.0;
    double startzLoc = 0.0;

    // create a particle at zero depth
    LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                               startdLoc, start_effdloc, startzLoc);

    // now get the cronus emulator
    double pressure = sample_eff_pressure[samp];
    double site_lat = sample_latitude[samp];
    double rho = 2650;
    double topo_scale =sample_shielding[samp];
    double snow_scale = 1;
    double N_26Al = 0;
    double sample_del26 =  0;
    double N_10Be_test = sample_Concentration_standardised[samp];
    double sample_del10 = sample_Concentration_uncertainty_standardised[samp];
    vector<double> erateinfo = eroded_particle.CRONUS_get_Al_Be_erosion(LSDCRNP, pressure,
                      site_lat, rho, N_10Be_test, N_26Al,sample_del10, sample_del26,
                      topo_scale,  snow_scale);
                      
    CRONUS_out << sample_name[samp] << "," << topo_scale << ",-99," << erateinfo[2]
               << "," << erateinfo[0] << "," << erateinfo[1] << ","
               << erateinfo[3] << "," <<  erateinfo[4] << endl;
  }
  CRONUS_out.close();

}
