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
    LSDCosmoData( string fname, string ftype)  { create(fname,ftype); }
    
    /// @brief this function prints the data held in the the data members
    ///  to screen. Is used for bug checking. 
    /// @detail Note the function does not print the standardised values, only raw values.
    /// @author SMM
    /// @date 09/02/2015
    void print_data_to_screen(); 
    
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
    
  protected:
    
    /// the number of samples
    int N_samples;
    
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
    
    /// a standardisation map for Be10
    map<string,double> standards_Be10;
    
    /// a standardisation map for Al26
    map<string,double> standards_Al26;
  
  private:
  
    /// @brief the empty create function
    void create();
    
    /// @brief the create function used when calling a file
    /// @param filename the name of the file WITH PATH AND EXTENSION
    /// @param  filtype the type of the file (either csv or txt)
    /// @author SMM
    /// @date 02/02/2015
    void create(string filename, string filetype);
    
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