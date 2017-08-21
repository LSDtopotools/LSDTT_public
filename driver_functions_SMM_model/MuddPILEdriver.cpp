//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// MuddPILEdriver
// The main driver function for the MuddPILE model
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2017 Simon M. Mudd 2017
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



#include <sstream>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "../LSDRasterModel.hpp"
#include "../LSDParticleColumn.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDParameterParser.hpp"
using namespace std;


int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "================================================================" << endl;
    cout << "|| Welcome to the MuddPILEdirver!                             ||" << endl;
    cout << "|| This program drives the MuddPILE lanscape evolution model. ||" << endl;
    cout << "|| One day this model will have documentation.                ||" << endl;
    cout << "|| This program was developed by Simon M. Mudd                ||" << endl;
    cout << "||  at the University of Edinburgh                            ||" << endl;
    cout << "================================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./MuddPILEdriver.exe /LSDTopoTools/Topographic_projects/Model_runs/ MuddPILE_test.driver" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // Parameters for initiating the model domain
  bool_default_map["read_initial_raster"] = false;
  int_default_map["NRows"] = 200;
  int_default_map["NCols"] = 400;
  float_default_map["DataResolution"] = 30;
  
  // Parameters for the initial surface
  bool_default_map["use_diamond_square_initial"] = true;
  float_default_map["diamond_square_relief"] = 10;
  int_default_map["diamond_square_feature_order"] = 8;
  bool_default_map["superimpose_parabola"] = true;
  float_default_map["parabola_relief"] = 10;
  bool_default_map["roughen_surface"] = true;
  float_default_map["roughness_relief"] = 0.25;
  bool_default_map["taper_edges"] = true;
  int_default_map["taper_rows"] = 10;
  
  // Parameters for spinning up the simulation
  float_default_map["spinup_K"] = 0.001;
  float_default_map["spinup_time"] = 50000;
  bool_default_map["staged_spinup"] = false;
  
  
  // control of how long you want the simulations to run
  
  
  // paramters for controlling model output
  int_default_map["print_interval"] = 10;
  bool_default_map["write_hillshade"] = false;

  // Use the parameter parser to get the maps of the parameters required for the 
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  
  
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;


  // Initiate the model object
  LSDRasterModel mod;


  // print parameters to screen
  mod.print_parameters();

  // see if we want to load a prior DEM
  bool create_initial_surface = false;
  if(this_bool_map["read_initial_raster"])
  {
    cout << "I am going to try to read an intial raster for you." << endl;
    cout << "The read filename is: " <<  DATA_DIR+DEM_ID << endl;
    cout << "I am going to IGNORE initial surface instructions!" << endl;
    string header =  DATA_DIR+DEM_ID+".hdr";
    ifstream file_info_in;
    file_info_in.open(header.c_str());
    // check if the parameter file exists
    if( file_info_in.fail() )
    {
      cout << "I found the header. I am loading this initial file. " << endl;
      LSDRaster temp_raster(DATA_DIR+DEM_ID,"bil");
      LSDRasterModel temp_mod(temp_raster);
      mod = temp_mod;
    }
    else
    {
      cout << "Header doesn't exist. I am switching to a fixed rows and columns." << endl;
      cout << "I will also create an initial surface for you." << endl;
      create_initial_surface = true;
    }
  }
  else
  {
    cout << "You have chosen not to read an initial raster so I will create an initial surface for you." << endl;
    create_initial_surface = true;
  }
  
  // create an initial surface if needed
  if(create_initial_surface)
  {
    cout << endl << endl << "=============================================" << endl;
    cout << "Creating and initial surface. " << endl;
    cout << "NRows: " << this_int_map["NRows"]<< ", NCols: " << this_int_map["NCols"] 
         << ", DataResolution: " << this_float_map["DataResolution"] << endl;
    mod.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"]);
    
    if (this_bool_map["use_diamond_square_initial"])
    {
      cout << "   I am starting with a fractal surface created by the diamond square algorithm." << endl;
      cout << "   It has a relief of " << this_int_map["diamond_square_feature_order"] << endl;
      cout << "   Let me check on the largest possible scale of features in the pseudo fractal surface." << endl;
      
      // we need to check the feature order to make sure it is not bigger than the largest dimension of the raster
      int largest_dimension;
      if( this_int_map["NRows"] >= this_int_map["NCols"]) 
      {
        largest_dimension = this_int_map["NRows"];
      }      
      else
      {
        largest_dimension = this_int_map["NCols"];
      }
      // now get the largest possible dimension
      float lps = floor( log2( float(largest_dimension) ) );
      int largest_possible_scale = int(lps);
      if( this_int_map["diamond_square_feature_order"] > largest_possible_scale )
      {
        this_int_map["diamond_square_feature_order"] = largest_possible_scale;
      }
      cout << "    The largest dimnesion is: " << largest_dimension << " pixels." << endl;
      cout << "    So the biggest scale is: " << pow(2,largest_possible_scale) << " pixels or 2 to the " << largest_possible_scale << " power" << endl;
      cout << "    I am setting the scale to: " << this_int_map["diamond_square_feature_order"] << endl; 

      // Now actually build the fractal surface
      mod.intialise_diamond_square_fractal_surface(this_int_map["diamond_square_feature_order"], this_float_map["diamond_square_relief"]);
    }
    
    
  }
  


}
