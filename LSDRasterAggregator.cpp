//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterAggregator.cpp
// Land Surface Dynamics RasterAggregator
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
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDRasterAggregator.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDRasterAggregator_CPP
#define LSDRasterAggregator_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The default create function. Doesn't do anything
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::create()
{}

void LSDRasterAggregator::create(string path_name, string param_name_prefix)
{
  // set the names of the data members
  path = path_name;
  param_name = param_name_prefix;

  // Set the parameters
  // The default slope parameter for filling. Do not change. 
  min_slope = 0.0001;

  // a boundary condition for the flow info object
  vector<string> boundary_conditionst(4);
  boundary_conditionst[0] = "n";
  boundary_conditionst[1] = "n";
  boundary_conditionst[2] = "n";
  boundary_conditionst[3] = "n";
  boundary_conditions = boundary_conditionst;

  // remove control characters from these strings
  path_name.erase(remove_if(path_name.begin(), path_name.end(), ::iscntrl), path_name.end());
  param_name_prefix.erase(remove_if(param_name_prefix.begin(), param_name_prefix.end(), ::iscntrl), param_name_prefix.end());
  
  string combined_filename = path_name+param_name_prefix;
  
  // load the file that contains the path to both the cosmo data and the 
  // DEMs
  
  // set up extensions to the files
  string rasters_ext = ".rasters";
  string params_ext = ".param";
  
  // get the filenames to open
  string rasters_fname =  path_name+param_name_prefix+rasters_ext;
  string parameters_fname = path_name+param_name_prefix+params_ext;

  // load the data. 
  cout << "Loading the names of the rasters." << endl;
  load_raster_filenames(rasters_fname);

  cout << "Loading parameters" << endl;
  load_parameters(parameters_fname);

  // check the rasters
  cout << "Hold on while I check the rasters. I don't want to eat dirty rasters!" << endl;
  check_rasters();
  cout << "The rasters seem okay. Mmmm, rasters." << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads parameters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::load_parameters(string filename)
{
  map<string,string> temp_map;
  
  ifstream infile;
  infile.open(filename.c_str());
  string parameter, value, lower;

  if( infile.fail() )
  {
    cout << "Parameter file: " << filename << " does not exist." << endl;
    cout << "Using default parameters." << endl;
  }

  // now ingest parameters
  while (infile.good())
  {
    parse_line(infile, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
    {
      lower[i] = tolower(parameter[i]);
    }

    cout << "parameter is: " << lower << " and value is: " << value << endl;

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);
    
    // add the parameter to the map
    if(temp_map.find(parameter) == temp_map.end())
    {
      temp_map[parameter] = value;
    }
    else
    {
      cout << "WARNING; you have two parameters with the same name." << endl;
      cout << "Using the most recent parameter value" << endl;
      temp_map[parameter] = value;
    }
  }
  infile.close();
  parameter_map = temp_map;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads the filenames 
// for the DEMs or single value parameters.
// It reads the DEM, either the snow shield raster or a single value
// the self shield raster name or a single value, and the topo shield raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::load_raster_filenames(string filename)
{
  //cout << "Getting filenames" << endl;
  
  // this vecvec holds data for determining the dem, the snow shielding
  // the self shielding and the topographic sheilding
  vector< string > temp_raster_names;
  
  // a string for null values
  string null_str = "NULL";
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv filenames file, but the file" << filename
         << "doesn't exist; LINE 348 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  
  string line_from_file;
  // now loop through the rest of the lines, getting the data. 
  while( getline(ifs, line_from_file))
  {
    line_from_file = RemoveControlCharactersFromEndOfString(line_from_file);
    temp_raster_names.push_back(line_from_file);
  }

  raster_filenames = temp_raster_names;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks the rasters for georeferencing and scaling
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::check_rasters()
{
  // check to see if the read extension is set
  
  string bil_ext = "bil";
  string null_str = "NULL";
  string raster_ext;
  if(parameter_map.find("dem read extension") == parameter_map.end())
  {
    cout << "You did not set the read extension. Defaulting to bil" << endl;
    raster_ext = bil_ext;
  }
  else
  {
    raster_ext = parameter_map["dem read extension"];
  }
  
  // check to make sure the read extension make sense
  if (raster_ext != "asc" && raster_ext != "flt" && raster_ext != "bil")
  {
    cout << "You didn't choose a valid raster format." << endl;
    cout << "Options are asc, bil and flt." << endl;
    cout << "You chose: " << raster_ext << endl;
    cout << "defaulting to bil, but your DEMs might not load" << endl;
    
    raster_ext = bil_ext;
  }

  // loop through the lines in the files, checking to see if the 
  // georeferencing is equivalent
  int N_rasters = int(raster_filenames.size());
  
  if(N_rasters < 1)
  {
    cout << "Sorry, I didn't find any rasters." << endl;
  }
  else if (N_rasters == 1)
  {
    cout << "You've only got one raster. Are you sure your raster file is correct?" << endl;
  }
  else
  {
    // get the information from the DEM
    // get the info from the DEM
    LSDRasterInfo DEM_info(raster_filenames[0],raster_ext);
    
    for(int iRaster = 1; iRaster<N_rasters; iRaster++)
    {
      // now compare with the other DEMs
      // first snow shielding
      if(raster_filenames[iRaster] != null_str)
      {
        LSDRasterInfo ThisRaster_info(raster_filenames[iRaster],bil_ext);
        if( ThisRaster_info!=DEM_info)
        {
          cout << "A raster " <<raster_filenames[iRaster] << " is not same shape and size as the DEM!" << endl;
          cout << "Setting it to NULL" << endl;
          raster_filenames[iRaster] = null_str;
        }
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-










void LSDSedimentRouting::create(string path_name, string param_name_prefix)
{
  cout << "Hello friends. I am creating a LSDSedimentRouting object" << endl;

  // set the names of the data members
  path = path_name;
  param_name = param_name_prefix;

  // Set the parameters
  // The default slope parameter for filling. Do not change. 
  min_slope = 0.0001;

  // a boundary condition for the flow info object
  vector<string> boundary_conditionst(4);
  boundary_conditionst[0] = "n";
  boundary_conditionst[1] = "n";
  boundary_conditionst[2] = "n";
  boundary_conditionst[3] = "n";
  boundary_conditions = boundary_conditionst;

  // remove control characters from these strings
  path_name.erase(remove_if(path_name.begin(), path_name.end(), ::iscntrl), path_name.end());
  param_name_prefix.erase(remove_if(param_name_prefix.begin(), param_name_prefix.end(), ::iscntrl), param_name_prefix.end());
  
  string combined_filename = path_name+param_name_prefix;
  
  // load the file that contains the path to both the cosmo data and the 
  // DEMs
  
  // set up extensions to the files
  string rasters_ext = ".rasters";
  string params_ext = ".param";
  
  // get the filenames to open
  string rasters_fname =  path_name+param_name_prefix+rasters_ext;
  string parameters_fname = path_name+param_name_prefix+params_ext;

  // load the data. 
  cout << "Loading the names of the rasters." << endl;
  load_raster_filenames(rasters_fname);

  cout << "Loading parameters" << endl;
  load_parameters(parameters_fname);
  
  cout << "I am going to check the parameter values now. " << endl;
  check_parameter_values();
  print_parameter_values_to_screen();

  // check the rasters
  cout << "Hold on while I check the rasters. I don't want to eat dirty rasters!" << endl;
  check_rasters();
  cout << "The rasters seem okay. Mmmm, rasters." << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks the parameter values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSedimentRouting::check_parameter_values()
{
  string this_string; 
  
  if(parameter_map.find("n_lithologies") == parameter_map.end())
  {
    cout << "You are missing the number of lithologies!!" << endl;
  }
  else
  {
    N_lithologies = atoi(parameter_map["n_lithologies"].c_str());
  }
  
  
  // check on the erodibility_coefficients
  if(parameter_map.find("erodibility_coefficients") == parameter_map.end())
  {
    cout << "You are missing the erodibility_coefficients!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["erodibility_coefficients"]);
    
    // a temporary vector for holding the data
    vector<float> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }
    erodibility_coefficients = temp_vec;
  }
  
  
  // check on the fertility_coefficients
  if(parameter_map.find("fertility_coefficients") == parameter_map.end())
  {
    cout << "You are missing the fertility_coefficients!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["fertility_coefficients"]);
    
    // a temporary vector for holding the data
    vector<float> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }
    fertility_coefficients = temp_vec;
  }

  // check on the source_1mm
  if(parameter_map.find("source_1mm") == parameter_map.end())
  {
    cout << "You are missing the source_1mm!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["source_1mm"]);
    
    // a temporary vector for holding the data
    vector<float> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }
    source_1mm = temp_vec;
  }
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks the parameter values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSedimentRouting::print_parameter_values_to_screen()
{
  cout << "===================================================" << endl;
  cout << "Here are your parameter values"  << endl;
  cout << "===================================================" << endl;
  cout << "N_lithologies: " << N_lithologies << endl;
  cout << "erodibility_coefficients:" << endl;
  for (int i = 0; i< int(erodibility_coefficients.size()); i++)
  {
    cout << " " <<erodibility_coefficients[i];
  }
  cout << endl;
  cout << "fertility_coefficients:" << endl;
  for (int i = 0; i< int(erodibility_coefficients.size()); i++)
  {
    cout << " " <<fertility_coefficients[i];
  }
  cout << endl;
  cout << "source_1mm:" << endl;
  for (int i = 0; i< int(erodibility_coefficients.size()); i++)
  {
    cout << " " <<source_1mm[i];
  }
  cout << endl;
  
  cout << "===================================================" << endl;
}


#endif