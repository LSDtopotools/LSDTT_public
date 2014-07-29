//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDAnalysisDriver.cpp
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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
#include <vector>
#include "LSDStatsTools.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDParticle.hpp"
#include "LSDParticleColumn.hpp"
#include "LSDAnalysisDriver.hpp"
using namespace std;

#ifndef LSDAnalysisDriver_CPP
#define LSDAnalysisDriver_CPP

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// The default constructor. This asks the user for a pathname and
// param filename
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::create()
{
  cout << "I need a parameter file to run. Please enter the path: " << endl;
  cin >> pathname;
  check_pathname_for_slash();

  cout << "Now I need a parameter filename: " << endl;
  cin >> param_fname;
  
  ingest_data(pathname, param_fname);
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function gets all the data from a parameter file
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::ingest_data(string pname, string p_fname)
{
  // the full name of the file
  string full_name = pname+p_fname;

  ifstream infile;
	infile.open(full_name.c_str());
	string parameter, value, lower;
	
	// now ingest parameters 
	while (infile.good())
	{
		parse_line(infile, parameter, value);
		lower = parameter;
		if (parameter == "NULL")
			continue;
		for (unsigned int i=0; i<parameter.length(); ++i)
			lower[i] = tolower(parameter[i]);

		if 	(lower == "dem read extension")		dem_extension = value;
		else if (lower == "dem read extension")		dem_read_extension = value;
    else if (lower == "dem write extension")		dem_write_extension = value;
    else if (lower == "write path")		write_path = value;
    else if (lower == "write fname")		write_fname = value;
    else if (lower == "write fill")		write_fill 	= (value == "true") ? true : false;
		else if (lower == "write nodeindex")		write_nodeindex 	= (value == "true") ? true : false;		

		else
    {
    	cout << "Line " << __LINE__ << ": No parameter '" 
           << parameter << "' expected.\n\t> Check spelling." << endl;
    }
	}
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a private function that makes sure the path has a slash at the end
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_pathname_for_slash()
{
  string lchar = pathname.substr(pathname.length()-2,1);
  string slash = "/";
  //cout << "lchar is " << lchar << " and slash is " << slash << endl;
     
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;  
    pathname = pathname+slash;
  } 
  cout << "The pathname is: " << pathname << endl;  
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function strips the text after the final dot in a string
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string LSDAnalysisDriver::get_string_before_dot(string this_string)
{
  string cut_string;
  unsigned found = this_string.find_last_of("."); 
  cut_string = this_string.substr(0,found); 
  return cut_string;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#endif
