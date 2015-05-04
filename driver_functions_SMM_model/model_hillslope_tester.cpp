#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include <python2.7/Python.h>
#include "../LSDRasterModel.hpp"

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
  }
	else
	{
		cout << "\n###################################################" << endl;
		cout << "No parameter file supplied" << endl;
		cout << "Creating a template parameter file (template_param)" << endl;
		cout << "###################################################" << endl;
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
	  string run_name = argv[3];
	  string full_run_name = pathname+run_name;
		// parameter two, run name
		mod.set_name(full_run_name);
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

	mod.print_parameters();

  // first change the size of the domain
  int newrows = 30;
  int newcols = 10;
  float new_resolution = 2.5;
  mod.resize_and_reset(newrows,newcols, new_resolution);

  // and turn the fluvial off, just to save some time
  mod.set_fluvial(false);
  
  // give an uplift rate for the surface
  float uplift_rate = 0.0005;     // this is in mm/yr
  
  float max_elev = 2.0;
  float edge_offset = 0.0;
  

  
  // get a steady state hillslope
  mod.initialise_nonlinear_SS(uplift_rate);
  
  // print this surface. This is for checking the steady state solution
  string DEM_name = "hillslope_test_SS";
  DEM_name = pathname+DEM_name;
  string asc_ext = "asc";
  mod.write_raster(DEM_name,asc_ext);
  
  // Now change to parabola
  mod.initialise_parabolic_surface(max_elev, edge_offset);

  // print this surface
  DEM_name = "hillslope_test_parabola";
  DEM_name = pathname+DEM_name;
  mod.write_raster(DEM_name,asc_ext);  
  
	// make an uplift field for the tester
	int uplift_mode = 0;
	float timestep = 50.0;
	mod.set_timeStep(timestep);
	
	// need to divide upflift rate by ts because this member function returns absoulte
	// elevation gain in a timestep rather than the uplift rate
	Array2D<float> uplift_field = mod.generate_uplift_field( uplift_mode, uplift_rate/timestep );
	
	cout << "Uplift field is: " <<  uplift_field << endl;
	
	int nrows = uplift_field.dim1();
	int ncols = uplift_field.dim2();
	Array2D<float> zero_array(nrows,ncols,0.0);
	
	// now run a few timesteps.
	float t_ime = 0;
	float iteration_tolerance = 0.000001;
	
	// initiate the assembler matrix
	mod.MuddPILE_initiate_assembler_matrix();
	
	cout << "entering nl loop" << endl;
	for(int ts = 1; ts<2001; ts++)
	{
	  
    mod.MuddPILE_nonlinear_creep_timestep(uplift_field, zero_array, iteration_tolerance);
    t_ime += timestep;
    if (ts%50 == 0)
    {
      cout << "time is: " << t_ime << endl;
    }
  }
  
  // print the file at the end
  string DEM_name2 = "hillslope_test_end";
  DEM_name2 = pathname+DEM_name2;
  mod.write_raster(DEM_name2,asc_ext);  
	
	
	//mod.run_model_from_steady_state();
	//mod.run_model();
	//cout << "\a";
	//mod.slope_area();
	//mod.show();
	//ss.str("");
	//ss << "./graph.py " << mod.get_name();
	//system(ss.str().c_str());

	return 0;
	
}
