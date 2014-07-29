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
#include <string>
#include <map>
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDStatsTools.hpp"
using namespace std;

#ifndef LSDAnalysisDriver_H
#define LSDAnalysisDriver_H


/// @brief This is a class to manage running LSDTopoTools. It parses a parameter
/// file and then manages running of analyses. 
/// @details The intention of this object is to run analyses via parameter
/// files and not through numerous compiled driver functions. We eventually
/// will want some kind of 'recorder' so that any time this object
/// runs an analysis it gives a full report of what analyses were run so that
/// results are reproducable
class LSDAnalysisDriver
{
	public:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Constructors
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    /// @brief The default constructor. 
    /// @details this asks for a pathname and a filename of the parameter file
    /// It then opens the paramter file and ingests the information
    /// @author SMM
    /// @date 29/07/2014    
	  LSDAnalysisDriver()			{ create(); }
	
    /// @brief this constructor just reads the param file given by the path and
    /// filename. You must give the parameter file extension!
    /// @param pname the pathname to the parameter file
    /// @param fname the filename of the parameter file !!INCLUDING EXTENSION!!
    /// @author SMM
    /// @date 29/07/2014	
    LSDAnalysisDriver(string pname, string pfname)			{ create(pname, pfname); }
 
    
	protected:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // DATA MEMBERS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    /// This vector holds various rasters computed during the trun
    vector<LSDRaster> vector_of_LSDRasters;
    
    /// This holds the flow info object
    LSDFlowInfo FlowInfo;
    
    /// index into the vector of rasters for the nodeindex raster
    int nodeindex_index;
        
    /// the path to the datafiles
    string pathname;
    
    /// the name of the parameter file
    string param_fname;
    		
		/// Extension for reading DEMs. Correspondence with write extensions is checked
		string dem_read_extension;
    
    /// Extension for writing DEMs. Correspondence with read extensions is checked
    string dem_write_extension;
    
    /// Path to files being written. Default is pathname
    string write_path;
    
    /// file prefix of files to be written. Default is the param name prefix 
    string write_fname;

    /// Path to files being read. Default is pathname
    string read_path;
    
    /// file prefix of files to be written. Default is the param name prefix 
    string read_fname;       

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Parameters for various analyses
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    /// the minimum slope for the fill function
    float min_slope_for_fill;
    
    /// the four boundary conditions on the raster for the flow info object 
    vector<string> boundary_conditions;

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Switches for running different analyses
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-                                
    /// this tells the program if it will need the fill raster
    bool need_fill;
    
    /// this tells the program if it will need a flow info object
    bool need_flow_info;

    /// This tells the driver whether or not to write the fill raster
    bool write_fill;
		
    /// This tells the driver whether or not to write the node index raster
    bool write_nodeindex;	
    
    /// This map holds all the possible analyses
    map<string,bool> analyses_switches; 
    
    /// This is a map that determines if various rasters are needed for the
    /// analysis. This ensures things like the fill raster are only calculated
    /// once
    map<string,bool> raster_switches;

    /// This is a map that tell where the indices into the raster vecs are    
    map<string,int> raster_indices;
    

    
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

    /// @brief create function
    /// @param pname the pathname to the parameter file
    /// @param fname the filename of the parameter file !!INCLUDING EXTENSION!!
    /// @author SMM
    /// @date 29/07/2014
    void create(string pname, string fname);	


    /// @brief This is the main function for parsing the parameter file
    /// @param pathname the path to the paramter file
    /// @param param_fname the name of the parameter file
    /// @author SMM
    /// @date 29/07/2014
    void ingest_data(string pname, string p_fname);
    
    /// @brief This wraps all the analyses, it first looks in the 
    /// raster_switces map to see what data it needs, and then
    /// calculates the necessary rasters. It then goes through and writes the
    /// rasters requested by the user. 
    /// @author SMM
    /// @date 29/07/2014
    void run_analyses();
    
    /// @brief This simply loads the base raster into the vector of rasters
    /// @author SMM
    /// @date 29/07/2014
    void read_base_raster();       
    
    /// @brief This calculates the fill raster
    /// @author SMM
    /// @date 29/07/2014
    void fill_raster();

    /// @brief This checks to see if boundary condtions have been assigned and 
    /// if not defaults to no flux boundaries
    /// @author SMM
    /// @date 29/07/2014
    void check_boundary_conditions();

    /// @brief This check to see if the filenames, paths and extensions have
    /// been assigned. If not it changes these to defaults
    /// @author SMM
    /// @date 29/07/2014
    void check_file_extensions_and_paths();
    
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