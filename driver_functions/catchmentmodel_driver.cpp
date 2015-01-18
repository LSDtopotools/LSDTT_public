#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include <python2.7/Python.h>
#include "../LSDCatchmentModel.hpp"

using namespace std;

bool file_check(string name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char *argv[])
{
	stringstream ss;
	
	if (argc < 2)
	{
		cout << "\n###################################################" << endl;
		cout << "No parameter file supplied" << endl;
		cout << "You must supply a path and parameter file!" << endl;
		cout << "see http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/" << endl;
		cout << "LSDCatchmentModel docs for assistance." << endl;
		cout << "###################################################" << endl;
		
		exit(0);  // Game over, try again.
	}

	if (argc == 2)
	// Bingo!
	{
		string pname = argv[1];
		string pfname = argv[2];
		// The path name and the parameter file name, respectively.
		// Remember: argc[0] is the program name that you just typed in to the terminal.
		LSDCatchmentModel modelrun(pname, pfname);
		
		modelrun.run_components();
	}
	
	else
	{
		cout << "\n###################################################" << endl;
		cout << "Either you have supplied too many arguments" << endl;
		cout << "or, somehow, you have supplied a negative number" << endl;
		cout << "of arguments. Specify the path name and the parameter file name only." << endl;
		cout << "There are no other options with this version yet." << endl;
		cout << "###################################################" << endl;		
	}

	cout << "The model ran successfully!" << endl;
	cout << "SCIENCE IS DONE" << endl;
	return 0;
}
