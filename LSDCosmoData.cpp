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
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <ctype.h>
#include <sstream>
#include <algorithm> 
#include <vector>
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDCosmoData.hpp"
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDBasin.hpp"
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
void LSDCosmoData::create(string filename, string filetype)
{
  if(filetype != "csv" && filetype != "txt")
  {
    cout << "LSDCosmoData line 77 You have not selected a valid filetype!" << endl;
    cout << "Options are csv and txt" << endl;
    cout << "Defaulting to txt" << endl;
    filetype = "txt";
  }
  
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

  if(filetype == "csv")
  {
    load_csv_cosmo_data(filename);
  }
  else
  {
    load_txt_cosmo_data(filename);
  }
  
  // now loop through the data, getting the standardised concentrations
  N_samples = int(sample_name.size());
  
  // create the vec vec for holding sample results
  vector<double> empty_vec;
  vector< vector<double> > result_vecvec;
  
  for(int i = 0; i<N_samples; i++)
  {
    // populate the results vector
    result_vecvec.push_back(empty_vec);
  
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
  erosion_rate_results = result_vecvec;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a csv file of cosmogenic data
// The first row is a header
// the following rows contain the data. 
// The data columns are:
//  column[0]: sample_name (NO SPACES OR COMMAS!!)
//  column[1]: latitude (decimal degrees)
//  column[2]: longitude (decimal degrees)
//  column[3]: Nuclide (Be10 or Al26)
//  column[4]: Nuclide concentration (atoms per gram)
//  column[6]: Nuclide uncertainty (atoms per gram)
//  column[7]: standardisation
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_csv_cosmo_data(string filename)
{
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv data file, but the file" << filename
         << "doesn't exist; LINE 219 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  // initiate temporary vectors
  vector<string> temp_sample_name;
  vector<double> temp_latitude;
  vector<double> temp_longitude;
  vector<string> temp_nuclide;
  vector<double> temp_Concentration_unstandardised;
  vector<double> temp_Concentration_uncertainty_unstandardised;
  vector<string> temp_standardisation;

  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  
  // discard the first line
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
      
      // remove constrol characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    // now convert the data
    temp_sample_name.push_back( this_string_vec[0] );
    temp_latitude.push_back( atof( this_string_vec[1].c_str() ) );
    temp_longitude.push_back( atof(this_string_vec[2].c_str() ) );
    temp_nuclide.push_back( this_string_vec[3] );
    temp_Concentration_unstandardised.push_back( atof(this_string_vec[4].c_str() ) );
    temp_Concentration_uncertainty_unstandardised.push_back( atof(this_string_vec[5].c_str() ) );
    temp_standardisation.push_back( this_string_vec[6] );
  
  }
  
  // now update the data members
  sample_name = temp_sample_name;
  latitude = temp_latitude;
  longitude = temp_longitude;
  nuclide = temp_nuclide;
  Concentration_unstandardised = temp_Concentration_unstandardised;
  Concentration_uncertainty_unstandardised = 
                         temp_Concentration_uncertainty_unstandardised;
  standardisation = temp_standardisation;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads the filenames 
// for the DEMs or single value parameters.
// It reads the DEM, either the snow shiled raster or a single value
// the self shield raster name or a single value, and the topo shield raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_DEM_and_shielding_filenames_csv(string filename)
{
  // this vecvec holds data for determining the dem, the snow shielding
  // the self shielding and the topographic sheilding
  vector< vector<string> > temp_DEM_names_vecvec;
  
  // a string for null values
  string null_str = "NULL";

  // this vecvec holds information about snow, self and toposhielding. 
  vector< vector<double> > temp_snow_self_topo_shielding_params; 

  
  // a vector of strings for holding the DEM names.
  // Elements without a DEM get a null value
  vector<string> null_string_vec(4,null_str);
  vector<string> this_snow_self_shield_names;
  
  // a vector for holding parameter values. 
  vector<double> empty_snow_self(2,0);
  vector<double> this_snow_self;
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv filenames file, but the file" << filename
         << "doesn't exist; LINE 219 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  
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
      
      // remove constrol characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    // reset the vectors for this line
    this_snow_self_shield_names = null_string_vec;
    this_snow_self = empty_snow_self;
    
    // now we need to see how many data elements we have
    int n_strings_in_line =  int(this_string_vec.size());
    
    // the can hold 1, 2, 3 or 4 elements. 
    // If it holds 1, it only has the name of the DEM
    // If it has 2, it has name of DEM and snow shielding
    // If it has 3, it has name of DEM and self shielding
    // If it has 4, the last element must be the name of the toposhield raster
    //  which has been precalculated
    if(n_strings_in_line == 1)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
    }
    else if(n_strings_in_line == 2)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
      
      // test if the second element is a number
      string second = this_string_vec[1];
      if(isnum(second[0]))
      {
        this_snow_self[0] = atof(this_string_vec[1].c_str());
      }
      else
      {
        this_snow_self_shield_names[1] = this_string_vec[1];
      }
      
    }
    else if(n_strings_in_line == 3)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
      
      // test if the second element is a number
      string second = this_string_vec[1];
      if(isnum(second[0]))
      {
        this_snow_self[0] = atof(this_string_vec[1].c_str());
      }
      else
      {
        this_snow_self_shield_names[1] = this_string_vec[1];
      }

      // test if the third element is a number
      string third = this_string_vec[2];
      if(isnum(third[0]))
      {
        this_snow_self[1] = atof(this_string_vec[2].c_str());
      }
      else
      {
        this_snow_self_shield_names[2] = this_string_vec[2];
      }            
    }    
    else if(n_strings_in_line == 4)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
      
      // test if the second element is a number
      string second = this_string_vec[1];
      if(isnum(second[0]))
      {
        this_snow_self[0] = atof(this_string_vec[1].c_str());
      }
      else
      {
        this_snow_self_shield_names[1] = this_string_vec[1];
      }

      // test if the third element is a number
      string third = this_string_vec[2];
      if(isnum(third[0]))
      {
        this_snow_self[1] = atof(this_string_vec[2].c_str());
      }
      else
      {
        this_snow_self_shield_names[2] = this_string_vec[2];
      }                  
      
      this_snow_self_shield_names[3] = this_string_vec[3];
    }
    else 
    {
      cout << "LSDCosmoData line 383, your cosmo data file names has the" << endl;
      cout << "wrong number of elements on this line. " << endl;
      cout << "I am only going to keep the DEM name" << endl;
      this_snow_self_shield_names[0] = this_string_vec[0];
    } 

    // push the parameters back into the vecvecs
    temp_DEM_names_vecvec.push_back(this_snow_self_shield_names);
    temp_snow_self_topo_shielding_params(this_snow_self);
  } 
  
  DEM_names_vecvec = temp_DEM_names_vecvec;
  snow_self_topo_shielding_params = temp_snow_self_topo_shielding_params;    
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a text file of cosmogenic data
// The first row is a header
// the following rows contain the data. 
// The data columns are:
//  column[0]: sample_name (NO SPACES OR COMMAS!!)
//  column[1]: latitude (decimal degrees)
//  column[2]: longitude (decimal degrees)
//  column[3]: Nuclide (Be10 or Al26)
//  column[4]: Nuclide concentration (atoms per gram)
//  column[6]: Nuclide uncertainty (atoms per gram)
//  column[7]: standardisation
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_txt_cosmo_data(string filename)
{
  cout << "Opening text file: " << filename << endl;
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: The file" << filename
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  // initiate temporary vectors
  vector<string> temp_sample_name;
  vector<double> temp_latitude;
  vector<double> temp_longitude;
  vector<string> temp_nuclide;
  vector<double> temp_Concentration_unstandardised;
  vector<double> temp_Concentration_uncertainty_unstandardised;
  vector<string> temp_standardisation;

  // initiate the string to hold the file
  string line_from_file;

  // now load the file. 
  double this_lat;
  double this_long;
  double this_conc;
  double this_uncert;
  
  string this_sample_name;
  string this_standard;
  string this_nuclide;
  
  // discard the first line
  ifs >> this_sample_name >> this_sample_name  >> this_sample_name >> this_sample_name 
             >> this_sample_name >> this_sample_name >> this_sample_name;
  
  while( ifs >> this_sample_name >> this_lat >> this_long >> this_nuclide 
             >> this_conc >> this_uncert >> this_standard)
  {
    temp_sample_name.push_back(this_sample_name);
    temp_latitude.push_back(this_lat);
    temp_longitude.push_back(this_long);
    temp_nuclide.push_back(this_nuclide);
    temp_Concentration_unstandardised.push_back(this_conc);
    temp_Concentration_uncertainty_unstandardised.push_back(this_uncert);
    temp_standardisation.push_back(this_standard);
  }

  // now update the data members
  sample_name = temp_sample_name;
  latitude = temp_latitude;
  longitude = temp_longitude;
  nuclide = temp_nuclide;
  Concentration_unstandardised = temp_Concentration_unstandardised;
  Concentration_uncertainty_unstandardised = 
                         temp_Concentration_uncertainty_unstandardised;
  standardisation = temp_standardisation;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function prints the data to screen
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_data_to_screen()
{
  cout << endl;
  cout << "==========================================================" << endl;
  cout << "PRINTING COSMO DATA HELD IN LSDCOSMODATA OBJECT" << endl;
  cout << "Sample_name\tLatitude\tLongitude\tNuclide\tConcentration\tUncertainty\tscaling\n";
  for(int i = 0; i<N_samples; i++)
  {
    cout << sample_name[i] << "\t"  << latitude[i] << "\t"  << longitude[i] << "\t"
         << nuclide[i] << "\t"  << Concentration_unstandardised[i] << "\t"
         << Concentration_uncertainty_unstandardised[i] << "\t"
         << standardisation[i] << "\n";
  }
  cout << "==========================================================" << endl;
  cout << endl << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the results to screen
// simple is for external, muon and production uncertainty only
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_simple_results_to_screen(double rho)
{
    
    cout << "=======================================================" << endl;
    cout << "Printing results of the analysis" << endl;
    for (int i = 0; i<N_samples; i++)
    {
      
      
      // don't print the results unless they exist
      if (int(erosion_rate_results[i].size()) > 0)
      {
        // get the results from this sample
        vector<double> erate_analysis = erosion_rate_results[i];
        
        cout << "-----------------------------------------------------" << endl;
        cout << "Sample " << sample_name[i] << " , a " << nuclide[i] << " sample" << endl;
        cout << "latitude:\t" << latitude[i] << "\tlongitude:\t" << longitude[i] << endl;
        cout << "Concentration: " << Concentration[i] << " +/- " 
             << Concentration_uncertainty[i] << " atoms/g" << endl;
        cout << "Erate is: " 
             <<  erate_analysis[0] << " g/cm^2/yr" << endl;
        cout << "The erosion rate for rho = "<< rho << " is: " << endl
             << erate_analysis[0]*10/rho << " m/yr and " 
             << erate_analysis[0]*1e6/rho << " cm/kyr" << endl;
        cout << "Ext uncert: " << erate_analysis[1] << " muon uncert: " << erate_analysis[2]
             << " production uncert: " <<  erate_analysis[3] << " g/cm^2/yr"  << endl;
        cout << "Total uncertainty in g/cm^2/yr: " << erate_analysis[4] << endl;
      }
    }
    cout << "=======================================================" << endl;

}


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



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the determination of cosmogenic erosion rates
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::basic_cosmogenic_analysis(int search_radius_nodes, 
                            int threshold_stream_order, LSDRaster& Elevations,
                            LSDRaster& TopoShield,
                            LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JNetwork)
{
  // the atmospheric data is in the folder with the driver_functions
  string path_to_atmospheric_data = "./";
  
  // first, convert the data into this UTM zone
  convert_to_UTM(Elevations);
  
  // now find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);

  // after this operation the three int vectors valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertainties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }
  
  // some environment variables
  double prod_uncert_factor = 1;          // this is a legacy parameter.
  string Muon_scaling = "Braucher";       // defaul muon scaling
  //bool data_from_outlet_only = true;      // this needs to be turned off later!!!
  
  
  // some temporary doubles to hold the nuclide concentrations
  double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
  double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin
  double test_N, test_dN;       // concentration and uncertainty of the nuclide in basin
  
  // now loop through the valid points, getting the cosmo data 
  for(int samp = 0; samp<n_valid_points; samp++)
  {
    if( valid_nuclide_names[samp] == "Be10")
    {
      test_N10 = valid_concentrations[samp];
      test_dN10 = valid_concentration_uncertainties[samp];
      test_N26 = 1e9;
      test_dN26 = 0;
      test_N = test_N10;
      test_dN = test_dN10;
    }
    else if( valid_nuclide_names[samp] == "Al26")
    {
      test_N10 = 1e9;
      test_dN10 = 0;
      test_N26 = valid_concentrations[samp];
      test_dN26 = valid_concentration_uncertainties[samp];
      test_N = test_N26;
      test_dN = test_dN26;
    }
    else
    {
      cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
      cout << "Defaulting to Be10" << endl;
      valid_nuclide_names[samp] = "Be10";
      test_N10 = valid_concentrations[samp];
      test_dN10 = valid_concentration_uncertainties[samp];
      test_N26 = 1e9;
      test_dN26 = 0;
      test_N = test_N10;
      test_dN = test_dN10;
    }
    
    cout << "Valid point is: " << valid_cosmo_points[samp] << " X: " 
         << UTM_easting[valid_cosmo_points[samp]] << " Y: "
         << UTM_northing[valid_cosmo_points[samp]] << endl;
    cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
         << snapped_junction_indices[samp] << endl;
    LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                            test_N10,test_dN10, test_N26,test_dN26);
    
    // populate the scaling vectors
    thisBasin.populate_scaling_vectors(FlowInfo, Elevations, TopoShield,
                                       path_to_atmospheric_data);
    
    // get the atmospheric pressure for bug checking. THis will print to screen
    //thisBasin.get_atmospheric_pressure(FlowInfo, Elevations, path_to_atmospheric_data);
    
    // now do the analysis
    vector<double> erate_analysis = thisBasin.full_CRN_erosion_analysis(test_N, 
                                        valid_nuclide_names[samp], test_dN, 
                                        prod_uncert_factor, Muon_scaling);
    
    erosion_rate_results[ valid_cosmo_points[samp] ] = erate_analysis;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-






#endif