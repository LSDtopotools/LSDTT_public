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
#include "../LSDRasterInfo.hpp"
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
  string asc_ext = "asc";
  string raster_format;
  string raster_ext;

  string temp;

  file_info_in.open(full_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: the header file \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  file_info_in >> temp >> data_path;
  data_path = RemoveControlCharactersFromEndOfString(data_path);

  file_info_in >> temp >> raster_name;
  raster_name = RemoveControlCharactersFromEndOfString(raster_name);

  file_info_in >> temp >> output_id;
  output_id = RemoveControlCharactersFromEndOfString(output_id);

  file_info_in >> temp >> N_iterations;

  file_info_in >> temp >> window_option;

  file_info_in >> temp >> log_bin_width;

  file_info_in >> temp >> raster_format;
  raster_format = RemoveControlCharactersFromEndOfString(raster_format);

  file_info_in.close();

  // print the spectral data to the data folder.
  output_id = data_path+output_id;


  // now check the raster format
  string lower = raster_format;
  for (unsigned int i=0; i<raster_format.length(); ++i)
  {
    lower[i] = tolower(raster_format[i]);
  }

  if (lower == "bil" || lower == "envi")
  {
    raster_ext = ENVI_ext;
    cout << "You have selected the ENVI bil format" << endl;
  }
  else if (lower == "flt" || lower == "float")
  {
    raster_ext = flt_ext;
    cout << "You have selected the flt format. Note that no georeferencing information" << endl
         << "will be preserved. Perhaps you should try the ENVI bil format?" << endl;
  }
  else if (lower == "asc" || lower == "ascii")
  {
    raster_ext = asc_ext;
    cout << "You have selected the ascii format. Note that no georeferencing information" << endl
         << "will be preserved. Also you are dealing with very large files," << endl
         << "Perhaps you should try the ENVI bil format?" << endl;
  }
  else
  {
    cout << "You did not choose a raster extension, defaulting to ENVI bil format" << endl;
    raster_ext = ENVI_ext;
  }

  // Here is the filename
  string DEM_f_name = data_path+raster_name;
  
  // check to see if the data exists
  LSDRasterInfo RInfo(DEM_f_name, raster_ext);

  // Load in data
  
  cout << "You are attempting to load the file: " <<  DEM_f_name<<"."<<ENVI_ext<< endl;
  LSDRaster raw_raster(DEM_f_name, raster_ext);
  cout << "I've loaded the raster" << endl;
  


  // convert to float by using the polyfit function
  //float window_radius = 61.0;
  //vector<int> raster_selection(8,0);
  //raster_selection[0] = 1;                    // set to return the smoothed surface
  //vector<LSDRaster> polyfit_rasters =
  //    raw_raster.calculate_polyfit_surface_metrics(window_radius, raster_selection);

  // remove the seas
  float sea_threshold = 1.0;
  raw_raster.mask_to_nodata_below_threshold(sea_threshold);

  // print the raw raster for bug checking
  raw_raster.write_raster(data_path+"raw",raster_ext);

  // Now try to get rid of internal nodata
  //cout << "I am trying to get rid of internal nodata values. " << endl;
  //int window_width = 5;
  //raw_raster.alternating_direction_nodata_fill(window_width);

  // now trim the raster
  cout << "Now I am trimming the raster." << endl;
  LSDRaster trimmed = raw_raster.RasterTrimmerSpiral();

  cout << "Trimmed raster, printing" << endl;

  string Trim_name = "_TRIM";
  trimmed.write_raster((DEM_f_name+Trim_name),raster_ext);

  cout << "Now performing spectral analysis" << endl;

  // Perform full spectral analysis
  LSDRasterSpectral SpectralRaster(trimmed);
  cout << "Loaded the SpectralRaster" << endl;
  
  SpectralRaster.full_spectral_analysis(log_bin_width,N_iterations,window_option);
  SpectralRaster.print_radial_spectrum(output_id);
  SpectralRaster.print_binned_spectrum(output_id, log_bin_width);

}