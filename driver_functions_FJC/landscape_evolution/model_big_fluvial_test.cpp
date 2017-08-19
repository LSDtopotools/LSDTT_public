//----------------------------------------------------------------------------//
// steady_state_fluvial_model.cpp
//
// This script creates a steady state landscape for a specific m/n landscape.
// The model is run in 3 stages:
//
// 1. The initial surface is created. Firstly, we generate a fractal surface
// using a diamond square algorithm with a desired relief.  We then taper the
// edges of the domain to sea level.  After the fractal surface is generated
// we add a parabola to the topography in order to stop lots of random depressions
// in the landscape. The user needs to specify the relief of the parabola - at
// the moment this is set to fractal_relief/2.
//
// 2. The initial surface is fully dissected by setting a high K value and uplift
// rate. We then run the model for 50,000 years.
//
// 3. The appropriate parameters are set using the chi values. We calculate chi
// for every point in the model, and use this to calculate the appropriate K value
// to get a desired relief of the final model domain.  We then use this K value
// and uplift rate to calculate the elevation of the surface, and snap to steady
// state.
//
// Authors: Simon M. Mudd and Fiona J. Clubb
// University of Edinburgh
//----------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include "../../LSDRasterModel.hpp"
#include "../../LSDParticleColumn.hpp"
#include "../../LSDParticle.hpp"
#include "../../LSDCRNParameters.hpp"
using namespace std;

bool file_check(string name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char *argv[])
{
  LSDRasterModel mod;
  stringstream ss;
  string pathname;
  string start_ext;
  string full_start_name;
  bool is_from_initial = false;

  if (argc <= 1)
  {
    cout << "ERROR: you need to enter a data folder and, \n"
         << "optionally, a parameter name and a run name" << endl;
    exit(0);
  }
  if (argc > 1)
  {
    pathname = argv[1];
    string lchar = pathname.substr(pathname.length()-2,1);
    string slash = "/";
    cout << "lchar is " << lchar << " and slash is " << slash << endl;

    if (lchar != slash)
    {
      cout << "You forgot the frontslash at the end of the path. Appending." << endl;
      pathname = pathname+slash;
    }
    cout << "The pathname is: " << pathname << endl;
  }
  if (argc >2 )
  {
    string param_name = argv[2];
    cout << "The parameter filename is: " << param_name << endl;
    string full_param_name = pathname+param_name;
    cout << "The full path is: " << full_param_name << endl;
    mod.initialize_model(full_param_name);

    // add the path to the default filenames
    mod.add_path_to_names( pathname);
  }
  else
  {
    cout << "\n###################################################" << endl;
    cout << "No parameter file supplied" << endl;
    cout << "Creating a template parameter file (template_param)" << endl;
    cout << "###################################################" << endl;
    // first change the default dimensions
    int newrows = 540;
    int newcols = 1080;
    float datares = 30;
    mod.resize_and_reset(newrows,newcols,datares);

    // set the end time, print interval, etc
    mod.set_K(0.00001);
    mod.set_endTime(50000);
    mod.set_print_interval(5);
    mod.set_print_hillshade(true);

    string template_param = "template_param.param";
    string full_template_name = pathname+template_param;
    mod.make_template_param_file(full_template_name);

    // add the path to the default filenames
    mod.add_path_to_names( pathname);

    // add random asperities to the surface of default model
    mod.random_surface_noise();
  }

  if (argc > 3)
  {
    string start_name = argv[3];
    full_start_name = pathname+start_name;
    // parameter two, run name
    //
    start_ext = "_fromstart";

    mod.set_name(full_start_name+start_ext);
    is_from_initial = true;

    cout << "This starts from an inital condition. I am turning the wash function off!" << endl;
    mod.set_threshold_drainage(-1);
    mod.set_current_frame(0);
  }
  else if (argc == 3)
  {
    if (file_check(mod.get_name() + "_report"))
    {
      string ans;
      cerr << "A run with the name '" << mod.get_name() << "' already exsist, do you wish to overwrite it? (y/n) ";
      cin >> ans;
      if (ans != "y")
      {
        cout << "You will need to choose another run name, exiting" << endl;
        exit(0);
      }
      else
      cout << "\nOverwriting" << endl;
    }
    ss << "rm " << mod.get_name() << "[0-9]*.asc";
    cout << ss.str() << endl;
    system(ss.str().c_str());
    ss.str("");
    ss << "rm " << mod.get_name() << "*_sa";
    cout << ss.str() << endl;
    system(ss.str().c_str());
  }

  // print parameters to screen
  mod.print_parameters();


  if (argc <= 3)
  {
    cout << "No initial topography loaded, running to a steady condition" << endl;
    // This is the desired relief of the fractal surface. If this is large, it takes longer
    // for the channels to fully dissect the landscape, but if it is too low
    // I think you will get very straight channels.
    float desired_relief = 10;

    // run to steady state using a large fluvial incision rate
    // (I don't know why but this seems to encourage complete dissection)
    // and turn the hillslopes off, just to save some time
    mod.set_hillslope(false);

    // lets use Declans fourier thing
    bool use_fourier = true;
    if (use_fourier)
    {
      cout << "I am going to start with a fractal surface" << endl;
      //float fractal_D = 2.7;
      //mod.intialise_fourier_fractal_surface(fractal_D);

      //mod.intialise_fourier_fractal_surface_v2(fractal_D,desired_relief);
      // The feature order sets the largest wavelength of repeating features.
      // The number of pixels of the largest feature is this 2^n where n is the
      // feature order
      int feature_order = 9;
      mod.intialise_diamond_square_fractal_surface(feature_order, desired_relief);

      // Make sure the edges are tapered
      int rows_to_taper = 10;
      mod.initialise_taper_edges_and_raise_raster(rows_to_taper);

      // add a bit of noise to get rid of straight channels on edge
      mod.set_noise(0.02);
      mod.random_surface_noise();

      // add a parabola to the topography
      // It is set with a relief. If this is similar to the relief of the fractal
      // surface then you are less likely to get straight channels in the middle
      // of the landscape where the model has had to fill low points created by the
      // fractal algorithm, but you will be more likely to get parallel rather
      // than dendritic channels.
      float parabolic_relief = desired_relief/2;
      mod.superimpose_parabolic_surface(parabolic_relief);
      cout << "Finished with the fractal surface " << endl;

      int frame = 555;
      mod.print_rasters(frame);
    }
    else
    {
      cout << "I will start with a low relief, randomly perturbed parabola" << endl;
      // initiate a very low relief parabolic surface
      // if the relief is too great you tend to get a network that
      // has many closely spaced parallel channels
      float max_elev = 0.2;
      float new_noise = max_elev/2;
      mod.set_noise(new_noise);
      float edge_offset = 0.0;
      mod.initialise_parabolic_surface(max_elev, edge_offset);
    }


    // now run to steady state
    // Set a relatively high K value and a high uplift rate to dissect the landscape
    // This happens quite fast at these settings! You might not even need
    // 50000 years and could possible shorten the end time.
    cout << "Dissecting landscape." << endl;
    mod.set_endTime(50000);
    mod.set_timeStep( 250 );
    mod.set_K(0.01);
    float new_baseline_uplift = 0.0025;
    mod.set_baseline_uplift(new_baseline_uplift);
    mod.set_current_frame(1);
    mod.run_components_combined();

    //mod.reach_steady_state();

    int frame = 777;
    //mod.print_rasters(frame);

    // a bit more dissection
    // I guess I just do this to make sure full dissection occurs.
    // You might reduce the end time. Note that the time is cumulative so this has
    // to be greater than the previous end time.
    mod.set_endTime(100000);
    mod.set_K(0.001);
    float U = 0.0001;    // a tenth of a mm per year
    mod.set_baseline_uplift(U);
    mod.run_components_combined();


    // Now set to steady state
    // This gives you the correct relief for a given U by back calculating K.
    // It is calculated directly from the chi plots.
    // Note that if the landscape is not completely dissected you will get weird
    // looking landscapes at this point.

    // Note that if you want to change m or n use
    // mod.set_m(0.5);
    // mod.set_n(1.5);
    // the m/n ratio will be calculated from these.
    mod.set_baseline_uplift(U);
    desired_relief = 1000;
    float new_K = mod.fluvial_snap_to_steady_state_tune_K_for_relief(U, desired_relief);
    cout << "Getting a steady solution for a landscape with relief of " << desired_relief
         << " metres and uplift of " << U*1000 << " mm per year." << endl;
    cout << "The new K is: " << new_K << endl;

    cout << "Now let me print the raster for you" << endl;
    mod.set_print_hillshade(true);

    // The frames are added to the filenames, so these set specific filenames
    // for, in this case, the steady state landscape.
    frame = 999;
    mod.print_rasters(frame);

    // add some time
    mod.set_endTime(150000);
    mod.run_components_combined();

    //float new_end_time = mod.get_current_time();
    //float this_end_time = mod.get_endTime();
    //float run_end_time = this_end_time;
    //float brief_end_time = run_end_time*3;
    //new_end_time = new_end_time+brief_end_time;
    //mod.set_endTime(new_end_time);
    //mod.run_components_combined();


  }
  else if(argc > 3)
  {
    cout << "you have chosen to load the file: " <<   full_start_name << endl;
    cout << "I am assuming this is an asc file";

    string asc_ext = "asc";
    mod.read_raster(full_start_name, asc_ext);
  }

  // Now run for a little while to diffuse out the hillslopes
  // this will run until the end time specified by
  // the data member
  // It will also print to file
  //cout << "Running now with hillslopes, until "
  //     << mod.get_endTime() << " years" << endl;
  //mod.run_components_combined();

  // now set the tilting going
  //int new_mode = 1;
  //float new_baseline_uplift = mod.get_max_uplift();
  //float new_max_uplift = (mod.get_max_uplift())*2;

  //mod.set_uplift(new_mode,new_max_uplift );
  //mod.set_baseline_uplift(new_baseline_uplift);

  // add some time
  //float new_end_time = mod.get_current_time();
  //float this_end_time = mod.get_endTime();
  //float run_end_time = this_end_time;
  //float brief_end_time = run_end_time*3;
  //new_end_time = new_end_time+brief_end_time;
  //mod.set_endTime(new_end_time);
  //mod.set_print_hillshade(true);

  /*
  // now reduce the fluvial efficiency and run again
  float new_K = mod.get_K();
  new_K = new_K/4;
  mod.set_K(new_K);

  float new_end_time = mod.get_current_time();
  float this_end_time = mod.get_endTime();
  new_end_time = new_end_time+this_end_time;

  mod.set_endTime(new_end_time);
  mod.run_components();

	//cout << "\ay";
	//mod.slope_area();
	//mod.show();
	//ss.str("");
	//ss << "./graph.py " << mod.get_name();
	//system(ss.str().c_str());
  */

	return 0;
}
