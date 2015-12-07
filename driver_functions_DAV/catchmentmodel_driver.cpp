#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include <python2.7/Python.h>
#include "../LSDCatchmentModel.hpp"

bool file_check(string name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char *argv[])
{
	std::cerr << "argc = " << argc << '\n';
    if ((argc > 0) && (argv[1] == NULL)) std::cerr << "Urk! 1\n";
    if ((argc > 1) && (argv[2] == NULL)) std::cerr << "Urk! 2\n";
    if ((argc > 2) && (argv[3] == NULL)) std::cerr << "Urk! 3\n";
    std::cerr << "Past Urks\n";
    
	std::stringstream ss;
	
	if (argc < 3)
	{
		std::cout << "\n###################################################" << std::endl;
		std::cout << "No parameter file supplied" << std::endl;
		std::cout << "You must supply a path and parameter file!" << std::endl;
		std::cout << "see http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/" << std::endl;
		std::cout << "LSDCatchmentModel docs for assistance." << std::endl;
		std::cout << "###################################################" << std::endl;
		
		exit(0);  // Game over, try again.
	}

	if (argc == 3)
	// Bingo!
	// DAV - Copied this out of "main_loop()" in origial CAESAR (after initialising 
	{
		
		std::string pname(argv[1]);
		std::string pfname(argv[2]);
		// The path name and the parameter file name, respectively.
		// Remember: argc[0] is the program name that you just typed in to the terminal.
		std::cout << "The pathname is: " << pname << " and the paramter file is: " << pfname << std::endl;
		
		LSDCatchmentModel modelrun(pname, pfname); // This create function calls 
    // initialise_variables() as well as creating an instance of LSDCatchmentModel called 'modelrun'
		std::cout << "The user-defined parameters have been ingested from the param file." << std::endl;
    
    modelrun.initialise_model_domain_extents();
    std::cout << "The model domain has been set from reading the elevation DEM header." << std::endl;
    
    modelrun.initialise_arrays(); // could be part of create() ??
    // also calls zero_values()
    std::cout << "The model field arrays have been initialised and zeroed where necessary." << std::endl;
    
		modelrun.load_data(); // Loads data from external files (DEM, hydroindex etc.)
    std::cout << "The terrain array and supplementary input data has been loaded." << std::endl;
    
    //modelrun.check_DEM_edge_condition();
    //std::cout << "Checking edge cells for suitable catchment outlet point..." << std::endl;
    // This function checks that there is at least one pixel on the edge of the
    // model domain for water/sediment to flow out of. (i.e. DEM cannot 
    // be surrounded entrely by NoData values.
    
    //modelrun.print_rainfall_data(); // Prints the rainfall data for checking

		modelrun.run_components();
	}
	
	else
	{
		std::cout << "\n###################################################" << std::endl;
		std::cout << "Either you have supplied too many arguments" << std::endl;
		std::cout << "or, somehow, you have supplied a negative number" << std::endl;
		std::cout << "of arguments. Specify the path name and the parameter file name only." << std::endl;
		std::cout << "There are no other options with this version yet." << std::endl;
		std::cout << "###################################################" << std::endl;	
		
		exit(0);	
	}

	std::cout << "The model ran successfully!" << std::endl;
	std::cout << "SCIENCE IS DONE" << std::endl;
	return 0;
}
