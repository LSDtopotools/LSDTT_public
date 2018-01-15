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
		float thisnumb  =  ((double) rand() / (RAND_MAX)) *10  ;
		// cout << thisnumb << endl;
		mouss[it] = thisnumb;
	}


	// TVD
	float lambda = 0;
	mouss_TVD = TV1D_denoise_v2(mouss, lambda);

	ofstream  chi_data_out;
	string filename = "test.csv";
    chi_data_out.open(filename.c_str());
    chi_data_out << "ID,original,TVD" << endl;
	// checking the results
	for (int it = 0; it <mouss.size(); it++)
	{
		//cout << mouss[it] << " -----> " << mouss_TVD[it] << endl;
		chi_data_out << it << "," << mouss[it] << "," << mouss_TVD[it] << endl;

	}
    chi_data_out.close();
	
}