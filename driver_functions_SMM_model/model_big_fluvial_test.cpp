#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include "../LSDRasterModel.hpp"
#include "../LSDParticleColumn.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
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
    int newrows = 500;
    int newcols = 1000;
    float datares = 30;
    mod.resize_and_reset(newrows,newcols,datares);

    // set the end time, print interval, etc
    mod.set_K(0.00001);
    mod.set_endTime(50000);  
    mod.set_print_interval(100);
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
    float desired_relief = 25;
    
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
      int feature_order = 8;
      mod.intialise_diamond_square_fractal_surface(feature_order, desired_relief);
      
      // Make sure the edges are tapered
      int rows_to_taper = 10;
      mod.initialise_taper_edges_and_raise_raster(rows_to_taper);
      cout << "Finished with the fractal surface " << endl;
      
      int frame = 9997;
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
    cout << "Dissecting landscape." << endl;
    mod.set_K(0.01);
    float new_baseline_uplift = 0.0025;
    mod.set_baseline_uplift(new_baseline_uplift);
    mod.reach_steady_state();
        
    int frame = 777;
    mod.print_rasters(frame);
    
    // a bit more dissection
    mod.set_K(0.001);
    float U = 0.0001;    // a tenth of a mm per year
    mod.set_baseline_uplift(U);
    mod.reach_steady_state();
    
    frame = 778;
    mod.print_rasters(frame);
    
    // Now set to steady state
    mod.set_baseline_uplift(U);
    desired_relief = 1000;
    float new_K = mod.fluvial_snap_to_steady_state_tune_K_for_relief(U, desired_relief);
    cout << "Getting a steady solution for a landscape with relief of " << desired_relief 
         << " metres and uplift of " << U*1000 << " mm per year." << endl;
    cout << "The new K is: " << new_K << endl;
    
    cout << "Now let me print the raster for you" << endl;
    mod.set_print_hillshade(true);
    
    frame = 999;
    mod.print_rasters(frame);

    mod.run_components_combined();

  
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
