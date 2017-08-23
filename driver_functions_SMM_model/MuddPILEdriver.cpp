//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// MuddPILEdriver
// The main driver function for the MuddPILE model
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
  int_default_map["NRows"] = 200;
  int_default_map["NCols"] = 400;
  float_default_map["DataResolution"] = 30;

  // Parameters that are used if you load a raster
  bool_default_map["read_initial_raster"] = false;
  float_default_map["minimum_elevation"] = 0.0;
  float_default_map["maximum_elevation"] = 30000;
  float_default_map["min_slope_for_fill"] = 0.0001;
  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  bool_default_map["print_raster_without_seas"] = false;

  // paramters for controlling model output
  int_default_map["print_interval"] = 10;
  bool_default_map["write_hillshade"] = true;

  // Parameters for the initial surface
  bool_default_map["use_diamond_square_initial"] = true;
  float_default_map["diamond_square_relief"] = 16;
  int_default_map["diamond_square_feature_order"] = 8;
  bool_default_map["taper_edges"] = true;
  int_default_map["taper_rows"] = 10;
  bool_default_map["superimpose_parabola"] = true;
  float_default_map["parabola_relief"] = 6;
  bool_default_map["roughen_surface"] = true;
  bool_default_map["fill_then_roughen_surface"] = true;
  float_default_map["roughness_relief"] = 0.25;

  // Parameters for spinning up the simulation
  bool_default_map["spinup"] = false;
  float_default_map["spinup_K"] = 0.001;
  float_default_map["spinup_U"] = 0.001;
  float_default_map["spinup_dt"] = 250;
  float_default_map["spinup_time"] = 20000;
  bool_default_map["staged_spinup"] = true;

  // control of m and n, and paramters for chi
  float_default_map["A_0"] = 1;
  float_default_map["m"] = 0.5;
  float_default_map["n"] = 1;
  float_default_map["dt"] = 250;
  int_default_map["print_interval"] = 10;

  // control of snapping to steady state
  bool_default_map["snap_to_steady"] = true;
  float_default_map["snapped_to_steady_uplift"] = 0.0001;
  float_default_map["snapped_to_steady_relief"] = 400;

  // Some parameters for very rudimentary steady forcing
  bool_default_map["run_steady_forcing"] = true;
  float_default_map["rudimentary_steady_forcing_time"] = 100000;
  float_default_map["rudimentary_steady_forcing_uplift"] = 0.0005;

  // Parameters for hillslopes
  bool_default_map["hillslopes_on"] = false;

  // Some parameters for transient forcing
  bool_default_map["run_transient_forcing"] = false;
  float_default_map["transient_forcing_time"] = 500000;
  float_default_map["transient_forcing_uplift"] = 0.001;

  // Use the parameter parser to get the maps of the parameters required for the analysis
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

  // Initiate the model object
  LSDRasterModel mod;

  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;
  mod.add_path_to_names(OUT_DIR);
  mod.set_name(OUT_DIR+OUT_ID);

  // set the print intervals and other parameters
  mod.set_m(this_float_map["m"]);
  mod.set_n(this_float_map["n"]);
  mod.set_print_hillshade(this_bool_map["write_hillshade"]);
  mod.set_timeStep( this_float_map["dt"] );
  mod.set_print_interval(this_int_map["print_interval"]);

  // print parameters to screen
  mod.print_parameters();

  // need this to keep track of the end time
  float current_end_time = 0;

  //============================================================================
  // Logic for reading an initial surface
  // It will look for a raster but if it doesn't find one it will default back to
  // fixed rows and columns. If you do load a raster it will override any
  // instructions about NRows, NCols, and DataResolution
  //============================================================================
  bool create_initial_surface = false;
  if(this_bool_map["read_initial_raster"])
  {
    cout << "I am going to try to read an intial raster for you." << endl;
    cout << "The read filename is: " <<  DATA_DIR+DEM_ID << endl;
    cout << "I am going to IGNORE initial surface instructions!" << endl;
    string header =  DATA_DIR+DEM_ID+".hdr";
    cout << "The full read path is: " << header << endl;
    ifstream file_info_in;
    file_info_in.open(header.c_str());
    // check if the parameter file exists
    if( not file_info_in.fail() )
    {
      cout << "I found the header. I am loading this initial file. " << endl;
      LSDRaster temp_raster(DATA_DIR+DEM_ID,"bil");
      if (this_bool_map["remove_seas"])
      {
        cout << "I am removing high and low values to get rid of things that should be nodata." << endl;
        float lower_threshold = this_float_map["minimum_elevation"];
        float upper_threshold = this_float_map["maximum_elevation"];
        bool belowthresholdisnodata = true;
        LSDRaster Flooded = temp_raster.mask_to_nodata_using_threshold(lower_threshold,belowthresholdisnodata);
        belowthresholdisnodata = false;
        temp_raster = Flooded.mask_to_nodata_using_threshold(upper_threshold,belowthresholdisnodata);
        if (this_bool_map["print_raster_without_seas"])
        {
          cout << "I'm replacing your input raster with a raster without seas." << endl;
          string this_raster_name = OUT_DIR+OUT_ID;
          temp_raster.write_raster(this_raster_name,raster_ext);
        }
      }

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

  //============================================================================
  // Logic for creating an initial surface
  // There are a number of options here, but the defaults are the ones we have
  // found to create the nicest initial surfaces.
  // The diamond square routine creates a pseudo-fractal surface. The algorithm
  // comes from Minecraft's creator, believe it or not.
  // If you use only the diamond square, you will get pits in the middle of the
  //  intial surface which will fill, leaving areas of the DEM with very
  //  straight channels. To mitigate this you can superimpose a parabola
  //  on the surface, and you can also run a fill function and then add roughness
  //  to the filled DEM (these are the defaults).
  // The greater the relief of the parabola relative to the relief of the
  //  diamond square fractal, the straighter your channels and basins will be
  //  after dissection.
  //============================================================================
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
      cout << "   It has a relief of " << this_int_map["diamond_square_relief"] << endl;
      cout << "   Let me check on the largest possible scale of features in the pseudo fractal surface." << endl;

      // we need to check the feature order to make sure it is not bigger than the largest dimension of the raster
      int smallest_dimension;
      if( this_int_map["NRows"] <= this_int_map["NCols"])
      {
        smallest_dimension = this_int_map["NRows"];
      }
      else
      {
        smallest_dimension = this_int_map["NCols"];
      }
      // now get the largest possible dimension
      float lps = floor( log2( float(smallest_dimension) ) );
      int largest_possible_scale = int(lps);
      if( this_int_map["diamond_square_feature_order"] > largest_possible_scale )
      {
        this_int_map["diamond_square_feature_order"] = largest_possible_scale;
      }
      cout << "    The largest dimnesion is: " << smallest_dimension << " pixels." << endl;
      cout << "    So the biggest scale is: " << pow(2,largest_possible_scale) << " pixels or 2 to the " << largest_possible_scale << " power" << endl;
      cout << "    I am setting the scale to: " << this_int_map["diamond_square_feature_order"] << endl;

      // Now actually build the fractal surface
      mod.intialise_diamond_square_fractal_surface(this_int_map["diamond_square_feature_order"], this_float_map["diamond_square_relief"]);

      if (this_bool_map["taper_edges"])
      {
        cout << "   Let me taper the edges of this fractal surface for you so that everything drains to the edge." << endl;
        cout << "   I am tapering along the " << this_int_map["taper_rows"] << " row closes to the N and S boundaries" << endl;
        mod.initialise_taper_edges_and_raise_raster(this_int_map["taper_rows"]);
      }
    }
    if (this_bool_map["superimpose_parabola"])
    {
      cout << "   I am superimposing a parabola with a relief of "  << this_float_map["parabola_relief"] << " metres" << endl;
      mod.superimpose_parabolic_surface(this_float_map["parabola_relief"]);
    }
    if(this_bool_map["roughen_surface"])
    {
      cout << "   I am going to roughen the surface for you. Roughness elements will " << endl;
      cout << "   have a maximum amplitude of " << this_float_map["roughness_relief"]<< " metres." << endl;
      mod.set_noise(this_float_map["roughness_relief"]);
      mod.random_surface_noise();
    }
    if(this_bool_map["fill_then_roughen_surface"])
    {
      cout << "   I am going to fill and then roughen the surface for you. Roughness elements will " << endl;
      cout << "   have a maximum amplitude of " << this_float_map["roughness_relief"]<< " metres." << endl;
      mod.fill(this_float_map["min_slope_for_fill"]);
      mod.set_noise(this_float_map["roughness_relief"]);
      mod.random_surface_noise();
    }
    cout << "Finished creating an initial surface. " << endl;
    cout << "=============================================" << endl << endl << endl;
  }

  //============================================================================
  // Logic for a spinning up the model.
  // Here spinup is used to dissect the landscape. Usually you create a
  // fractal surface first and then dissect it to get an initial channel and
  // basin geometry. You can then snap it to steady state with the snapping
  // functions.
  // This uses fluvial only! You need to run hillslope diffusion after this
  // if you want a steady landscape with hillslopes
  //============================================================================
  if(this_bool_map["spinup"])
  {
    cout << "I am going to spin the model up for you." << endl;
    cout << "This will rapidly develop a drainage network. It uses only fluvial incision." << endl;
    cout << "The typical spinup method is to dissect the landscape and then bring it to " << endl;
    cout << "steady state using the snap_to_steady flag." << endl;

    // turn the hillslope diffusion off
    mod.set_hillslope(false);

    current_end_time = this_float_map["spinup_time"];

    mod.set_endTime(this_float_map["spinup_time"]);
    mod.set_timeStep( this_float_map["spinup_dt"] );
    mod.set_K(this_float_map["spinup_K"]);
    mod.set_baseline_uplift(this_float_map["spinup_U"]);
    mod.set_current_frame(1);
    mod.run_components_combined();

    if(this_bool_map["staged_spinup"])
    {
      current_end_time = this_float_map["spinup_time"]+2*this_float_map["spinup_time"];
      mod.set_endTime(current_end_time);
      mod.set_K(this_float_map["spinup_K"]*0.1);
      mod.set_baseline_uplift(this_float_map["spinup_U"]*0.1);
      mod.run_components_combined();
    }
  }

  //============================================================================
  // Logic for snapping to steady state using equation 4a from Mudd et al 2014
  // Every pixel is considered to be a channel. If you want hillslopes
  // you will need to run the model after this with the hillslopes turned on
  //============================================================================
  if(this_bool_map["snap_to_steady"])
  {
    cout << "I am going to snap the landscape to steady state. " << endl;
    cout << "The way this works is that chi is calculated and then the steady state" << endl;
    cout << " elevation is calculated using equation 4a from Mudd et al 2014 JGR-ES" << endl;
    cout << "This method assumes only fluvial erosion. If you want hillslopes you need" << endl;
    cout << "To turn them on and then let the model run to steady state. " << endl;

    float new_K = mod.fluvial_snap_to_steady_state_tune_K_for_relief(this_float_map["snapped_to_steady_uplift"], this_float_map["snapped_to_steady_relief"]);
    cout << "Getting a steady solution for a landscape with relief of " << this_float_map["snapped_to_steady_relief"]
         << " metres and uplift of " << this_float_map["snapped_to_steady_uplift"]*1000 << " mm per year." << endl;
    cout << "The new K is: " << new_K << endl;
    mod.set_K(new_K);
  }

  //============================================================================
  // Logic for a rudimentary steady forcing of uplift
  //============================================================================
  if(this_bool_map["run_steady_forcing"])
  {
    if( not this_bool_map["hillslopes_on"] )
    {
      cout << "I'm turning hillslope diffusion off." << endl;
      mod.set_hillslope(false);
    }
    cout << "Let me run some steady forcing for you. " << endl;
    current_end_time = current_end_time+float_default_map["rudimentary_steady_forcing_time"];
    mod.set_endTime(current_end_time);
    mod.set_baseline_uplift(float_default_map["rudimentary_steady_forcing_uplift"]);
    mod.run_components_combined();
  }

  //============================================================================
  // Logic for a transient forcing. This at the moment is set with block
  // uplift (I hope that is right...)
  // FJC 23/08/17
  //============================================================================
  if(this_bool_map["run_transient_forcing"])
  {
    if( not this_bool_map["hillslopes_on"] )
    {
      cout << "I'm turning hillslope diffusion off." << endl;
      mod.set_hillslope(false);
    }
    cout << "Let me run some transient forcing for you. " << endl;
    current_end_time = current_end_time+float_default_map["transient_forcing_time"];
    mod.set_endTime(current_end_time);
    //mod.set_baseline_uplift(float_default_map["rudimentary_steady_forcing_uplift"]);
    mod.set_uplift(0, float_default_map["transient_forcing_uplift"]);
    mod.run_components_combined();
  }

}
