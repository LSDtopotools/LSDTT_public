#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../../LSDStatsTools.hpp"


int main (int nNumberofArgs,char *argv[])
{
	cout << "this is a small program to test statistical functions without having to run heavy analysis" << endl;


	// creatintg a vector of float
	vector<float> mouss(10000), mouss_TVD(10000);
	// generating random values
	for (size_t it = 0; it <mouss.size(); it++)
	{
		float thisnumb  = ((double) rand() / (RAND_MAX)) + ((double) rand() / (RAND_MAX)) + 1 + ((double) rand() / (RAND_MAX)) *100  + ((double) rand() / (RAND_MAX)) *2000  ;
		// cout << thisnumb << endl;
		mouss[it] = thisnumb;
	}


	// TVD
	const double lambda = 12.6;
	mouss_TVD = TV1D_denoise_v2(mouss, lambda);

	// checking the results
	for (size_t it = 0; it <mouss.size(); it++)
	{
		cout << mouss[it] << " -----> " << mouss_TVD[it] << endl;

	}

	
}