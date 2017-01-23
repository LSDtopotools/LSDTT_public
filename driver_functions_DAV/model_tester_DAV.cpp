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
	if (argc > 1)
		// parameter one, param file name
		mod.initialize_model(argv[1]);
	else
	{
		cout << "\n###################################################" << endl;
		cout << "No parameter file supplied" << endl;
		cout << "Creating a template parameter file (template_param)" << endl;
		cout << "###################################################" << endl;
		mod.make_template_param_file("template_param");
		
		// add random asperities to the surface of default model
		mod.random_surface_noise();
	}

	if (argc > 2)
		// parameter two, run name
		mod.set_name(argv[2]);
	else if (argc == 2)
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
	
	cout << "Initialising Fractal Landscape Surface..." << endl;
	mod.intialise_fourier_fractal_surface(0.5);
	cout << "Fractal landscape succesfully created!" << endl;
	
	mod.print_parameters();

	//mod.reach_steady_state();
	//mod.run_model_from_steady_state();
	mod.run_model();
	cout << "\a";
	//mod.slope_area();
	mod.show();
	//ss.str("");
	//ss << "./graph.py " << mod.get_name();
	//system(ss.str().c_str());

	return 0;
}
