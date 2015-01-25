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
	{
		
		std::string pname = argv[1];
		std::string pfname = argv[2];
		// The path name and the parameter file name, respectively.
		// Remember: argc[0] is the program name that you just typed in to the terminal.
		std::cout << "The pathname is: " << pname << " and the paramter file is: " << pfname << std::endl;
		
		LSDCatchmentModel modelrun(pname, pfname);
		std::cout << "The parameters have been successfully ingested and the model initialised" << std::endl;
		//modelrun.run_components();
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
