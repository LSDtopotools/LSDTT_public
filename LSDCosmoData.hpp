//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDCosmoData.hpp
// Land Surface Dynamics CosmoData
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for keeping track of cosmogenic data
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd                                                    
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <fstream>
#include <math.h>
#include <iostream>
#include <map>
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
using namespace std;

#ifndef LSDCosmoData_HPP
#define LSDCosmoData_HPP

class LSDCosmoData
{
  public:
  
    /// @brief the default constructor. This doesn't do anything. 
    LSDCosmoData()                                { create(); }
    
    /// @brief This constructor requires a filename
    /// @detail The file format is: 
    /// sample_name nuclide latitude longitude concentration uncertainty standard
    ///  the sample name cannot have a space 
    ///  the nuclide name must be Be10 or Al26
    ///  lat and long are in decimal degrees
    ///  concetration is in atoms per gram
    ///  uncertainty is the AMS error in atoms per gram
    ///  standard is the laboratory standardisation. The options can be found
    ///  in the code of this function; the numbers come from Balco et al's (2009)
    ///  CRONUS calculator. 
    /// @param filename the string of the file within which the cosmo data
    ///  is stored. This INCLUDES PATH AND EXTENSION.
    /// @param ftype the type of file, options are csv and txt. If not one of these
    ///  options the code assumes txt
    /// @author SMM
    /// @date 06/02/2015
    LSDCosmoData( string path, string file_prefix)  { create(path,file_prefix); }
    
    /// @brief  THis function loads the crn data rom either a csv or text file
    /// @param filename the name of the file
    /// @param filetype this is either csv or txt
    /// @author SMM
    /// @date 26/02/2015
    void load_cosmogenic_data(string filename, string filetype);
    
    /// @brief This function load a csv file containing names of DEMs and 
    ///  (possibly) sheilding rasters or shielding parameters
    /// @detail The file should have one DEM name per row. Each row
    ///  can have between 1 and 4 elements. 
    /// The elements are:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    /// If elements are missing, they are considered null arguments. 
    /// @param filename the name of the csv file
    /// @author SMM
    /// @date 26/02/2015
    void load_DEM_and_shielding_filenames_csv(string filename);
    
    /// @brief This loads parameters for the comsogenic calculations
    ///  it uses the parse_files frinction in LSDStatsTools, so the format
    ///  is the paramterer name, followed by ":", followed by a space and then
    ///  the parameter value
    /// @param filename a string of the full filename
    /// @author SMM
    /// @date 02/03/2015
    void load_parameters(string filename);
    
    /// @brief this gets the names of the DEMs to be used in the analysis
    /// @detail only returns the DEM, not snow shielding, topo shielding, etc
    ///  rasters
    /// @return DEM_fnames a vector of fname strings
    /// @author SMM 
    /// @date 19/03/2015
    vector<string> get_DEM_fnames();
    
    /// @brief this gets the names of the snow shielding rasters
    //  to be used in the analysis
    /// @detail returns only snow shielding names. If name does not exist,  
    ///  returns NULL.
    /// @return Snow_fnames a vector of fname strings
    /// @author SMM 
    /// @date 07/07/2015
    vector<string> get_Snow_fnames();

    /// @brief this gets the names of the self shielding rasters
    //  to be used in the analysis
    /// @detail returns only snow shielding names. If name does not exist,  
    ///  returns NULL.
    /// @return Self_fnames a vector of fname strings
    /// @author SMM 
    /// @date 07/07/2015
    vector<string> get_Self_fnames();

    /// @brief This function checks to make sure parameter values are
    ///  valid for the cosmo data
    /// @author SMM
    /// @date 03/03/2015
    void check_parameter_values();
    
    /// @brief this function checks the existence and georeferencing of 
    ///  the rasters outlined in the file list
    /// @author SMM
    /// @date 03/03/2015
    void check_rasters();
    
    /// @detail This function takes a DEM and then spawns dems that are clipped
    ///  to the basin boundaries, with padding to incorporate surrounding high
    ///  topography and CRN samples that might lie outside the basin due to 
    ///  GPS errors
    /// @param DEM_fname the name of the DEM with FULL PATH but excluding extension
    /// @param basin_padding_px the number of pixels that should be padded from
    ///  the basin
    /// @author SMM
    /// @date 18/03/2015
    vector<string> spawn_clipped_basins(string DEM_fname, int basin_padding_px);

    /// @brief This dirves the spawning of basins
    /// @param path This is a string containing the path to the data files (needs / at the end)
    /// @param prefix the prefix of the data files
    /// @param padding_pixels the number of pixels with which to pad the basins
    /// @author SMM
    /// @date 10/07/2015
    void BasinSpawnerMaster(string path, string prefix, int padding_pixels);

    /// @brief This calculates topographic shielding for basins listed in the 
    ///   _CRNRasters.csv file
    /// @detail Shielding rasters are printed to the same folder as the DEM
    /// @param path This is a string containing the path to the data files (needs / at the end)
    /// @param prefix the prefix of the data files
    /// @author SMM
    /// @date 15/07/2015
    void RunShielding(string path, string prefix);

    /// @brief this function calculates the UTM coordinates of all the sample
    ///  points for a given UTM zone. 
    /// @param UTM_zone the UTM zone
    /// @author SMM
    /// @date 06/02/2015
    void convert_to_UTM(int UTM_zone);
    
    /// @brief this function calculates the UTM coordinates of all the sample
    ///  points for a given UTM zone. It determines the UTM zone from a raster
    /// @param Raster the LSDRaster from which the UTM zone is determined 
    /// @author SMM
    /// @date 06/02/2015
    void convert_to_UTM(LSDRaster& Raster);    

    /// @brief This function calculates the CRN erosion rate from the data
    ///  stored within the data elements of the object
    /// @detail This analysis does not calculate snow shielding or landsliding
    /// @param DEM_prefix a string holding the prefix of the DEM (without the bil)
    ///  NOTE the DEM must be in bil format.
    /// @author SMM
    /// @date 09/02/2015 
    void basic_cosmogenic_analysis(string DEM_prefix);

    /// @brief This function computes erosion rates and uncertainties for 
    ///  a given DEM. It is wrapped by a function that goes through
    ///  the list of DEM, providing this function with the raster names
    ///  and the parameters for the model run
    /// @param Raster_names a vector of strings with 4 elements:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    ///  It there is no DEM then this is set to "NULL"
    ///  @param CRN_params this contains the single shielding depths for snow
    ///   and self shielding if the rasters are not supplied. 
    /// @author SMM
    /// @date 28/02/2015
    void full_shielding_cosmogenic_analysis(vector<string> Raster_names,
                            vector<double> CRN_params);

    /// @brief This function wraps the cosmogenic rate calculators.
    /// @detail Looks throught the vecvecs listing file locations and then
    ///   finds valid CRN data, and runs the erosion rate routine for these
    ///   data.
    /// @param method_flag: an int that sets the method
    ///  default: use snow and self shielding
    ///  method_flag == 0, basic analysis (snow and self set to 0, calculate toposhield)
    ///  method_flag == 1, snow and self sheilding, can load toposhield if provided
    /// @author SMM
    /// @date 28/02/2015
    void calculate_erosion_rates(int method_flag);

    /// @brief This function prints sevear rasters to file:
    ///  1) Pixel-by-pixel production scaling
    ///  2) Pixel-by-pixel combined scaling (production plus combined shielding)
    ///  3) Pixel-by-pixel combined shielding
    ///  4) The cosmo concetration coming from each pixel at the proscribed erosion rate
    /// @param Raster_names a vector containing the names of DEM, snow, self and toposhield rasters
    /// @param CRN_params a vector of parameter values
    /// @author SMM
    /// @date 23/03/2015
    void full_shielding_raster_printer(vector<string> Raster_names,
                                        vector<double> CRN_params);

    /// @brief this function prints the data held in the the data members
    ///  to screen. Is used for bug checking. 
    /// @detail Note the function does not print the standardised values, only raw values.
    /// @author SMM
    /// @date 09/02/2015
    void print_data_to_screen(); 

    /// @brief This prints the cosmo data to a new file 
    /// @param path is the pathname. Needs a "/"  at the end
    /// @param prefix: the name of the file before _CRNData.csv
    /// @author SMM
    /// @date 14/07/2015
    void print_renamed_cosmo_data(string path, string prefix);

    /// @brief This prints the parameter data to a new file 
    /// @param path is the pathname. Needs a "/"  at the end
    /// @param prefix: the name of the file before _CRNData.csv
    /// @author SMM
    /// @date 14/07/2015
    void print_renamed_parameter_data(string path, string prefix);

    /// @brief this prints the file structure data to screen
    ///  it is a list of the DEMs, snow shielding rasters, self shelding rasters
    ///  and topo shielding rasters used in the analysis
    ///  the later three rasters can be NULL values, and the snow and self shielding
    ///  rasters can be replaced by single values
    /// @author SMM
    /// @date 26/02/2015
    void print_file_structures_to_screen();

    /// @brief this prints all the data, parameters and file structures to screen
    ///  It is used to keep a record of CRN erosion rate computations so
    ///  analyses can be reproduced
    /// @param outfilename the name of the outfile including full path
    /// @author SMM
    /// @date 02/03/2015
    void print_all_data_parameters_and_filestructures(string outfilename);

    /// @brief Prints the simple results to the screen
    /// @detail The 'simple' is because it only looks at external, 
    ///  muon, and production uncertainties. 
    /// @param rho the rock denisty in kg/m^3
    /// @author SMM
    /// @date 10/02/2015
    void print_simple_results_to_screen(double rho);

    /// @brief This prints the results to a csv file and to a file for
    ///  passing to CRONUS
    /// @detail The columns in the CSV file are:
    ///  sample_name
    ///  nuclide
    ///  latitude
    ///  longitude
    ///  concentration (atms/g)
    ///  concentration_uncert (atoms/g)
    ///  erosion rate g_percm2_peryr
    ///  erosion rate AMS_uncert g_percm2_peryr
    ///  muon_uncert g_percm2_peryr
    ///  production_uncert g_percm2_peryr
    ///  total_uncert g_percm2_peryr
    ///  AvgProdScaling dimensionless
    ///  AverageTopoShielding dimensionless
    ///  AverageSelfShielding dimensionless
    ///  AverageSnowShielding dimensionless
    ///  AverageCombinedScaling dimensionless (this is averaged production scaling times toposhielding)
    ///  outlet_latitude
    ///  OutletPressure hPa
    ///  OutletEffPressure hPa (pressure needed to get basin averaged production scaling)
    ///  centroid_latitude
    ///  CentroidPressure hPa
    ///  CentroidEffPressure (pressure needed to get basin averaged production scaling)
    ///  ErosionRate_in_cmperkyr (to check against cosmocalc, assumes 2650 kg/m^2)
    ///  ErosionRate_COSMOCALC_in_g_percm2_peryr (assumes 2650 kg/m^2): The erosion
    ///   rate you would get if you took production weighted scaling and used
    ///   cosmocalc. 
    ///  ErosionRate_COSMOCALC_cmperkyr (assumes 2650 kg/m^2): The erosion
    ///   rate you would get if you took production weighted scaling and used
    ///   cosmocalc. 
    /// @author SMM
    /// @date 12/03/2015
    void print_results();

    /// @brief This function prints several rasters to file:
    ///  1) Pixel-by-pixel production scaling
    ///  2) Pixel-by-pixel combined scaling (production plus combined shielding)
    ///  3) Pixel-by-pixel combined shielding
    ///  4) The cosmo concetration coming from each pixel at the proscribed erosion rate
    /// @author SMM
    /// @date 23/03/2015
    void print_rasters();
    
  protected:
    
    /// the number of samples
    int N_samples;
    
    /// the path to the cosmo data. Also used to print results
    string path;
    
    /// the prefix of the parameter files
    string param_name;
    
    /// A vector of the sample names
    vector<string> sample_name;
    
    /// a vector holding the latitude of the samples
    vector<double> latitude;
    
    /// a vector holding the longitude of the samples
    vector<double> longitude;
    
    /// a vector holding the UTM of the samples
    vector<double> UTM_easting;
    
    /// a vector holding the UTM northing of the samples
    vector<double> UTM_northing;
    
    /// the nuclide. Only options are Be10 and 26Al. 
    vector<string> nuclide;
    
    /// a string holding the standardisation
    vector<string> standardisation;

    /// a vector holding the concetration, in atoms per gram, of the sample
    vector<double> Concentration;
    
    /// a vector holding the AMS uncertainty of the sample, in atoms per gram
    vector<double> Concentration_uncertainty; 
    
    /// a vector holding the concetration, in atoms per gram, of the sample
    /// this is for data before standardisation
    vector<double> Concentration_unstandardised;
    
    /// a vector holding the AMS uncertainty of the sample, in atoms per gram
    /// this is for data before standardisation
    vector<double> Concentration_uncertainty_unstandardised; 
    
    /// a vector of vectors holding the results of the cosmogenic analysis
    vector< vector<double> > erosion_rate_results;
    
    /// a standardisation map for Be10
    map<string,double> standards_Be10;
    
    /// a standardisation map for Al26
    map<string,double> standards_Al26;
    
    /// a string vecvec for holding the DEM names involved in the analysis
    ///  elements are
    ///  [0] = DEM name
    ///  [1] = snow_shield raster
    ///  [2] = self shield raster
    ///  [3] = topo shield raster
    ///  When there is no raster, this is "NULL"
    vector< vector<string> > DEM_names_vecvec;
    
    /// A vector holding the parameters for snow and self shielding if it is
    ///  a single value.
    vector< vector<double> > snow_self_topo_shielding_params;
    
    /// The minimum slope for the fill function
    float min_slope;
  
    /// The number of pixels for a channel
    int source_threshold;
    
    /// The number of pixels over which to search for a channel
    int search_radius_nodes;
  
    /// The minimum stream order to be considered a channel for cosmo sampling
    int threshold_stream_order;

    /// The azimuth step for topographic shielding calculations
    int theta_step;
  
    /// The inclination step for topographic sheilding calculations
    int phi_step;

    /// an uncertainty parameter which was superceded by the new
    /// error analyses but I have been too lazy to remove it. Does nothing 
    double prod_uncert_factor;
    
    /// The muon production scaling. Options are "Braucher", "Granger" and "Schaller"
    string Muon_scaling;       

    /// the atmospheric data is in the folder with the driver_functions, 
    /// but can be changed if necessary.
    string path_to_atmospheric_data;
  
    /// The boundary conditions for the flow info object
    vector<string> boundary_conditions;
    
    //---------------Flags for writing files---------------------
    /// Write toposheild rasters if they don't exist
    bool write_TopoShield_raster;
    
    /// Write a LSDIndexRaster with the basins
    bool write_basin_index_raster;
    
    /// Write the shielding and scaling rasters
    bool write_full_scaling_rasters;
    
    //-----------------Information used in cosmogenic calculators---------------
    /// This contains data with all sorts of scaling parameters
    /// for calculation of erosion rates using other calculators
    map<string, map<int,double> > MapOfProdAndScaling;

    //---------------basin metrics---------------------
    /// The Basin Relief
    vector<double> MBS;

  private:
  
    /// @brief the empty create function
    void create();
    
    /// @brief the create function used when calling a file
    /// @param the path to the file. Must have the "/" at the end
    /// @param  file_prefice the prefix (without extension) of the parameter files
    /// @author SMM
    /// @date 02/02/2015
    void create(string path, string file_prefix);
    
    /// @brief This loads data from a text file
    /// @detail The data columns are:
    ///  column[0]: sample_name (NO SPACES OR COMMAS!!)
    ///  column[1]: latitude (decimal degrees)
    ///  column[2]: longitude (decimal degrees)
    ///  column[3]: Nuclide (Be10 or Al26)
    ///  column[4]: Nuclide concentration (atoms per gram)
    ///  column[6]: Nuclide uncertainty (atoms per gram)
    ///  column[7]: standardisation
    /// @param filename the name of the file WITH PATH AND EXTENSION
    /// @author SMM
    /// @date 09/02/2015
    void load_txt_cosmo_data(string filename);
    
    /// @brief This loads data from a csv file
    /// @detail The data columns are:
    ///  column[0]: sample_name (NO SPACES OR COMMAS!!)
    ///  column[1]: latitude (decimal degrees)
    ///  column[2]: longitude (decimal degrees)
    ///  column[3]: Nuclide (Be10 or Al26)
    ///  column[4]: Nuclide concentration (atoms per gram)
    ///  column[6]: Nuclide uncertainty (atoms per gram)
    ///  column[7]: standardisation
    /// @param filename the name of the file WITH PATH AND EXTENSION
    /// @author SMM
    /// @date 09/02/2015
    void load_csv_cosmo_data(string filename);
};


#endif