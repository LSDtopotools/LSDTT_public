//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Lengthscale_correlations.cpp
// David Milodowski September 2014
// Creates spectrally filtered rasters at different lengthscales and compares the
// correlation of the two datasets.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../TNT/tnt.h"
int main (int nNumberofArgs,char *argv[])
{
  // ADD ALL PATH RELATED STUFF HERE!
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];

	cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;

	string raster_name;
	string data_path;
	string output_id;
	int window_option,N_iterations;
	float log_bin_width;
	string flt_ext = "flt";
	string ENVI_ext = "bil";

  string temp;                
	
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the header file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}	
  file_info_in >> temp >> data_path;                                     
	file_info_in >> temp >> raster_name;                                   
	file_info_in >> temp >> output_id;
	file_info_in >> temp >> N_iterations;
	file_info_in >> temp >> window_option;
	file_info_in >> temp >> log_bin_width;
	
	
	file_info_in.close();

  // Load in data
	LSDRaster raster(data_path+raster_name, flt_ext);
   
  // Perform full spectral analysis  
  LSDRasterSpectral SpectralRaster(raster);
  SpectralRaster.full_spectral_analysis(log_bin_width,N_iterations,window_option);
  SpectralRaster.print_radial_spectrum(output_id);
  SpectralRaster.print_binned_spectrum(output_id, log_bin_width);

}                                                                                 