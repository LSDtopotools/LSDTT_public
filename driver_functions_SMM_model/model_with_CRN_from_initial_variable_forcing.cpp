#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include <python2.7/Python.h>
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
  //LSDRasterModel mod(0,0);
  LSDRasterModel mod;
  stringstream ss;
  string pathname;
  string start_ext;
  string full_start_name;
  bool is_from_initial = false;
  
  bool is_variable_U = false;
  bool is_variable_D = false;
  bool is_variable_K = false;
  
  if (argc <= 1)
  {
    cout << "=============================================================\n" 
         << "Welcome to a raster model driver; this is for running models that have time varying" << endl
         << "uplift and sediment transport and fluvial erosion coefficients"
         << "This driver takes 2, 3 or 4 arguments to the command line" << endl
         << "Argument 1: The path name (with a slash at the end)" << endl
         << "Argument 2: Switch to turn variable uplift on (0 == off, 1 == on)" << endl
         << "Argument 3: Switch to turn variable D (for hillslopes) on (0 == off, 1 == on)" << endl
         << "Argument 4: Switch to turn variable K (for channels) on (0 == off, 1 == on)" << endl
         << "Argument 5: The parameter file name" << endl
         << "Argument 6: The name of the inital condition raster (without the extension)" << endl
         << "*For example:" << endl
         << "./model_with_CRN_from_initial.out /home/smudd/SMMDataStore/analysis_for_papers/SOS_paper/model_cosmo_tracker/" << endl
         << "0 1 1 Variable_BL_Padova.param InitialForCRN" << endl
         << "=============================================================\n";
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
  if (argc > 2)
  {
    int us = atoi(argv[2]);
    if (us == 1)
    {
      is_variable_U = true;
      cout << "You have chosen variable uplift" << endl;
    }
    else
    {
      cout << "Uplift will remain constant in this run" << endl;
    }
    
  }
  if (argc > 3)
  {
    int ds = atoi(argv[3]);
    if (ds == 1)
    {
      is_variable_D = true;
      cout << "You have chosen variable diffusivity" << endl;
    }
    else
    {
      cout << "Diffusivity will remain constant in this run" << endl;
    }
  }
  if (argc > 4)
  {
    int ks = atoi(argv[4]);
    if (ks == 1)
    {
      is_variable_K = true;
      cout << "You have chosen variable K" << endl;
    }
    else
    {
      cout << "Kwill remain constant in this run" << endl;
    }
  }
  if (argc >5 )
  {
    string param_name = argv[5];
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
    int newrows = 150;
    int newcols = 300;
    mod.resize_and_reset(newrows,newcols);

    // set the end time, print interval, etc
    mod.set_K(0.0001);
    mod.set_endTime(50000);  
    mod.set_print_interval(25);

    string template_param = "template_param.param";
    string full_template_name = pathname+template_param;		
    mod.make_template_param_file(full_template_name);

    // add the path to the default filenames
    mod.add_path_to_names( pathname);

    // add random asperities to the surface of default model
    mod.random_surface_noise();
  }
  if (argc > 6)
  {
    string start_name = argv[6];
    full_start_name = pathname+start_name;
    // parameter two, run name
    //
    //start_ext = "_fromstart";  
    
    //mod.set_name(full_start_name+start_ext);
    is_from_initial = true;
    
    cout << "This starts from an inital condition. I am turning the wash function off!" << endl;
    cout << "The full starting condition is called: " << full_start_name << endl;
    mod.set_threshold_drainage(-1);
    mod.set_current_frame(0);
  }
  else if (argc == 6)
  {
    if (file_check(mod.get_name() + "_report"))
    {
      string ans;
      cerr << "A run with the name '" << mod.get_name() 
           << "' already exsist, do you wish to overwrite it? (y/n) ";
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


  cout << "argc is: " << argc << endl;

  if (argc <= 6)
  {
    cout << "No initial topography loaded, running to a steady condition" << endl;
    // initiate a very low relief parabolic surface
    // if the relief is too great you tend to get a network that 
    // has many closely spaced parallel channels  
    float max_elev = 0.2;
    float new_noise = max_elev/2;
    mod.set_noise(new_noise);
    float edge_offset = 0.0;
    mod.initialise_parabolic_surface(max_elev, edge_offset);
    
    // run to steady state using a large fluvial incision rate 
    // (I don't know why but this seems to encourage complete dissection)
    // and turn the hillslopes off, just to save some time
    mod.set_hillslope(false);
    
    // we want to develop the channel network but in a low relief landscape, 
    // so we temporarily crank up the fluvial erodibility
    float overall_K = mod.get_K();
    float high_K = 1.4*overall_K;
    mod.set_K(high_K);
    cout << "Setting K to : " << high_K << " was: " << overall_K << endl;
   
    // now run to steady state
    cout << "Running fluvial only under steady forcing" << endl;
    mod.reach_steady_state();
  
    // for some reason this doesn't seem to get it completely to steady state
    // so run some more
    mod.run_components();	

    // reset the fluvial component
    mod.set_K(overall_K);
    cout << "Reset K to original value" << endl;
   
    // set the hillslope parameters to true 
    mod.set_hillslope(true);
    mod.set_nonlinear(true);
  
    // you'll need to extend the end time
    float new_end_time = mod.get_current_time();
    float this_end_time = mod.get_endTime();
    float run_end_time = this_end_time;
    float brief_end_time = run_end_time*3;
    new_end_time = new_end_time+brief_end_time; 
    mod.set_endTime(new_end_time);
  }
  else if(argc > 6)
  {
    cout << "you have chosen to load the file: " <<   full_start_name << endl;
    cout << "I am assuming this is an asc file";
    
    string asc_ext = "asc";
    mod.read_raster(full_start_name, asc_ext);  
  }
  
  
  int startType = 0;
  int startDepth = 5;
  double particle_spacing = 0.2;
  int column_spacing = 3;
  LSDCRNParameters CRNParam;
  CRNParam.set_Neutron_only_parameters();     // set to neutron only production
  double rho_r = 2000;
  double this_U = mod.get_max_uplift();         // this is in m/yr
  //double eff_U = rho_r*this_U/10;             // uplift rate in g/cm^2/yr
  
  vector<int> CRNc_rows;
  vector<int> CRNc_cols;
  
  // initiate the CRNParticleColumns
  vector<LSDParticleColumn> ErodedParticles_vec; 
  vector<LSDParticleColumn> CRNParticleColumns_vec = 
     mod.initiate_steady_CRN_columns(column_spacing,CRNc_rows, CRNc_cols, rho_r,
                  this_U, startType,startDepth,particle_spacing, CRNParam);
  
  // force steady state
  mod.force_initial_steady_state();
  
  // Now run for a little while to diffuse out the hillslopes
  // this will run until the end time specified by
  // the data member
  // It will also print to file
  cout << "Running now with hillslopes, until "
       << mod.get_endTime() << " years" << endl; 
  // now start to vary diffusitivy
  //
  //mod.set_D_mode(new_D_mode);             
  
  int new_D_mode = 1;
  int new_K_mode = 1;
  double uplift_amplitude_fraction = 1.0;
  
  string sU, sD, sK; // strings for appending the run name
    
  if (is_variable_U)
  {
    cout << "This model will be run with variable uplift." << endl;
    mod.set_periodic_uplift(uplift_amplitude_fraction); 
    sU = "1";
  }
  else
  {
    sU = "0";
  }
  if (is_variable_D)
  {
     cout << "This model will be run with variable hillslope diffusivity." << endl; 
     mod.set_D_mode(new_D_mode);
     sD = "1";
  }
  else
  {
    sD = "0";
  }
  if (is_variable_K)
  {
     cout << "This model will be run with variable K." << endl; 
     mod.set_K_mode(new_K_mode);
     sK = "1";
  }
  else
  {
    sK = "0";
  }
  
  // now set appended run name:
  cout << endl <<endl << "=================" << endl << "Appending run name" << endl;
  string append_run_name = "_"+sU+"_"+sD+"_"+sK+"var";
  mod.append_run_name(append_run_name);

  

  mod.run_components_combined_cell_tracker( CRNParticleColumns_vec,
                      ErodedParticles_vec, startType, startDepth, particle_spacing, 
                      CRNParam);

  //mod.run_components_combined();


/*
  // add some time
  float new_end_time = mod.get_current_time();
  float this_end_time = mod.get_endTime();
  float run_end_time = this_end_time;
  float brief_end_time = run_end_time;
  new_end_time = new_end_time+brief_end_time; 
  mod.set_endTime(new_end_time);
  
  // now start to vary diffusitivy
  int new_D_mode = 1;
  mod.set_D_mode(new_D_mode);

  cout << "Running now with variable D, until "
       << mod.get_endTime() << " years" << endl; 
  mod.run_components_combined();
*/
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
