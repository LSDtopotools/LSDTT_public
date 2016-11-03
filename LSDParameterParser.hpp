//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDParameterParser.hpp
//
// Land Surface Dynamics Parameter Parser Object
//
// An object for keeping track of parameters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for calculating concntration of environmental tracers, CRNs, TCN, fallout
//  nuclides
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2016 Simon M. Mudd 2016
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
#include <string>
#include <vector>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
using namespace std;
using namespace TNT;


#ifndef LSDParameterParser_HPP
#define LSDParameterParser_HPP

/// @brief Object to perform flow routing.
class LSDParameterParser
{
  public:
    /// @brief The create function. This is default and throws an warning.
    /// @author SMM
    /// @date 11/02/16
    LSDParameterParser()        { create(); }

    /// @brief Creates a LSDParamParser object using a path and filename
    /// @param PathName the name of the path
    /// @param FileName the name of the parameter file
    /// @author SMM
    /// @date 3/07/2015
    LSDParameterParser(string PathName, string FileName)
                     { create(PathName,FileName); }  

    /// @brief This is the function for reading parameters
    /// @ detail It is identical to the parse_line function in LSDStatsTools but
    ///   has a bigger buffer to cope with long path names
    /// @param infile The file being read
    /// @param pramater a string of the parameter name: returned from the parser
    ///   this is always changed to lower case only!!
    /// @param value A string holding a value. It will need to be further parsed
    ///   depending on the name of the parameter (i.e., to a float, int, string, etc)
    /// @author SMM
    /// @date 02/11/16
    void LSDPP_parse_line(ifstream &infile, string &parameter, string &value);

    /// @brief This takes a parameter file and parses all the parameters into
    ///  a map with string keys and values. Further functions are used to convert these
    ///  into the appropriate data types
    /// @author SMM
    /// @date 03/11/2016
    void LSDPP_parse_file_into_parameter_map(string FullName);
    
    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are floats
    /// @author SMM
    /// @date 03/11/2016
    void parse_float_parameters(map<string,float> default_map);

    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are int
    /// @author SMM
    /// @date 03/11/2016
    void parse_int_parameters(map<string,int> default_map);
    
    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are bool
    /// @author SMM
    /// @date 03/11/2016
    void parse_bool_parameters(map<string,bool> default_map);
    
    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are strings
    /// @author SMM
    /// @date 03/11/2016
    void parse_string_parameters(map<string,string> default_map);        



    /// @brief This forces the read and wirte extensions to bil
    /// @author SMM
    /// @date 02/11/16
    void force_bil_extension(); 
    
    /// @return read_extension
    string get_dem_read_extension() const        { return dem_read_extension; }
    /// @return write_extension
    string get_dem_write_extension() const        { return dem_write_extension; }
    /// @return write_path
    string get_write_path() const      { return write_path; }
    /// @return read_path
    string get_read_path() const      { return read_path; }
    /// @return write_fname
    string get_write_fname() const  { return write_fname; }
    /// @return read_fname
    string get_read_fname() const      { return read_fname; }
    /// @return CHeads_file
    string get_CHeads_file() const      { return CHeads_file; }
    /// @return the boundary conditions
    vector<string> get_boundary_conditions() const     { return boundary_conditions; }

    /// @brief This is the main function for parsing the parameter file
    /// @param FullName The full name, including path, of the parameter file
    /// @author SMM
    /// @date 02/11/2016
    void ingest_data(string FullName);

    /// @brief This checks to see if boundary condtions have been assigned and 
    /// if not defaults to no flux boundaries
    /// @author SMM
    /// @date 02/11/2016
    void check_boundary_conditions();

    /// @brief This check to see if the filenames, paths and extensions have
    /// been assigned. If not it changes these to defaults
    /// @author SMM
    /// @date 02/11/2016
    void check_file_extensions_and_paths();

    /// @brief this returns the string before the last dot in a string. 
    /// so for example if you gave it paramfile.param it would return paramfile
    /// @param this_string the string you want truncated
    /// @return the truncated string
    /// @author SMM
    /// @date 29/07/2014
    string get_string_before_dot(string this_string);
    
    /// @brief This takes a map of default float parameters and updates them
    ///  with any parameters that have been fed to the parser
    /// @param float_default_map A map containing the defualt parameters
    /// @return a map containing the parameters for this run
    /// @author SMM
    /// @date 02/11/2016
    map<string,float> set_float_parameters(map<string,float> float_default_map);
    
    /// @brief This takes a map of default int parameters and updates them
    ///  with any parameters that have been fed to the parser
    /// @param int_default_map A map containing the defualt parameters
    /// @return a map containing the parameters for this run
    /// @author SMM
    /// @date 02/11/2016
    map<string,int> set_int_parameters(map<string,int> int_default_map);
    
    /// @brief This takes a map of default bool parameters and updates them
    ///  with any parameters that have been fed to the parser
    /// @param bool_default_map A map containing the defualt parameters
    /// @return a map containing the parameters for this run
    /// @author SMM
    /// @date 02/11/2016
    map<string,bool> set_bool_parameters(map<string,bool> bool_default_map);

    
  protected:

    /// The path to the parameter file (used if no read or wirte path is supplied)
    string param_file_path;
    
    /// The path to the parameter file (used if no read or wirte path is supplied)
    string param_fname;

    /// The exteniosn of the rasters; either bil, asc or flt
    string dem_read_extension;
    
    /// The exteniosn of the rasters; either bil, asc or flt
    string dem_write_extension;

    /// Path to which files will be written
    string write_path;
    
    /// Path from which files will be written
    string read_path;
    
    /// Prefix of files to be written (i.e., no path, no extension)
    string write_fname;
    
    /// Prefix of files to be read (i.e., no path, no extension)
    string read_fname;
    
    /// The prefix of the channelheads filename
    string CHeads_file;

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    /// the four boundary conditions on the raster for the flow info object 
    vector<string> boundary_conditions;

    /// This is the map of all the paramters read from file, in string format
    /// you use the default parameter maps to convert the strings to the
    /// parameter values you want. 
    map<string,string> parameter_map;
    
    /// This is used for bug checking: all defaults used are then passed to this
    /// map, which can be printed and used to check the spelling of the default
    /// parameters
    map<string,string> defaults_used_map;

    /// This holds float parameters
    map<string,float> float_parameters;
    
    /// This holds integer parameters
    map<string,int> int_parameters;

    /// This holds string parameters
    map<string,string> string_parameters;
    
    /// This holds bool parameters
    map<string,bool> bool_parameters;
    
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Switches for running different analyses
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    /// This map holds all the possible analyses
    map<string,bool> analyses_switches; 

    /// This vector holds various rasters computed during the run
    map<string,string> map_of_LSDRasters;
    
    /// as above, but these are index rasters
    map<string,LSDIndexRaster> map_of_LSDIndexRasters;
  
    /// This is a map  container that determines if various rasters are needed for the
    /// analysis. This ensures things like the fill raster are only calculated
    /// once
    map<string,bool> raster_switches;

    /// This is a map that tell where the indices into the raster vecs are    
    map<string,int> raster_indices;
    

    
    /// This holds integer vectors. Can be used to get sources
    map<string, vector<int> > integer_vector_map;
    
    /// This holds names of supporting files, for example files that contain
    /// node of junction indices to be loaded. 
    map<string,string> support_file_names;

    /// This holds names of methods. For example, if the key is drainage_area_method, the string is
    /// the method which is used to calculate drainage area
    map<string,string> method_map;

  private:
    void create();
    void create(string PathName, string FileName);
};

#endif