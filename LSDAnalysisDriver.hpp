//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDAnalysisDriver
// Land Surface Dynamics Analysis Driver
//
// This object parses parameter files and drives analysis. Its purpose is
// to stop having to write a bunch of .cpp driver functions and instead
// be able to write a parameter files without compiling
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2014 Simon M. Mudd 2014
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
#include <iostream>
#include <vector>
#include <list>
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDStatsTools.hpp"
using namespace std;

#ifndef LSDAnalysisDriver_H
#define LSDAnalysisDriver_H


/// @brief This is a class for a particle that can be tracked through simulations
/// and retains data about position and chemical content
class LSDAnalysisDriver
{
	public:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Constructors
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    /// @brief The default constructor. It is the only possible constructor
	  LSDAnalysisDriver()			{ create(); }
	

 
    
	protected:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // DATA MEMBERS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    /// the path to the datafiles
    string pathname;
    
    /// the name of the parameter file
    string param_fname;
    
    /// the extension of the dem. Must be either asc or flt. There is a function
    /// for checking this. 
		string dem_extension;
		
		/// Extension for reading DEMs. Correspondence with other extensions is checked
		string dem_read_extension;
    
    /// Extension for writing DEMs. Correspondence with other extensions is checked
    string dem_write_extension;
    
    /// Path to files being written. Default is pathname
    string write_path;
    
    /// file prefix of files to be written. Default is the param name prefix 
    string write_fname;
    
    /// This tells the driver whether or not to write the fill raster
    bool write_fill;
		
    /// This tells the driver whether or not to write the node index raster
    bool write_nodeindex;	
    
    
  private:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // CREATE FUNCTIONS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
    /// @brief default create function
    /// @details this asks for a pathname and a filename of the parameter file
    /// It then opens the paramter file and ingests the information
    /// @author SMM
    /// @date 29/07/2014
    void create();	

    /// @brief This is the main function for parsing the parameter file
    /// @param pathname the path to the paramter file
    /// @param param_fname the name of the parameter file
    /// @author SMM
    /// @date 29/07/2014
    void ingest_data(string pname, string p_fname);
    
    /// @brief this adds a slash to the end of the pathname
    /// @author SMM
    /// @date 29/07/2014
    void check_pathname_for_slash();	

    /// @brief this returns the string before the last dot in a string. 
    /// so for example if you gave it paramfile.param it would return paramfile
    /// @param this_string the string you want truncated
    /// @return the truncated string
    /// @author SMM
    /// @date 29/07/2014
    string get_string_before_dot(string this_string);

};

#endif