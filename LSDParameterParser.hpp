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

    /// @brief Creates a LSDParamParser object with a filename. Must include path
    /// @param fname String of the binary flowinfo data file to be read.
    /// @author SMM
    /// @date 11/02/16
    LSDParameterParser(string FullName)    { create(FullName); }

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
    /// @date 11/02/16
    void LSDPP_parse_line(ifstream &infile, string &parameter, string &value);
    
    /// @brief This is the main function for parsing the parameter file
    /// @param FullName The full name, including path, of the parameter file
    /// @author SMM
    /// @date 02/11/2016
    void ingest_data(string FullName);

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
    
  protected:

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
    
    /// This holds float parameters
    map<string,float> float_parameters;
    
    /// This holds integer parameters
    map<string,int> int_parameters;
    
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
    void create(string FullName);
    void create(string PathName, string FileName);
};

#endif