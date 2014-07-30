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
#include <string>
#include <map>
#include "LSDStatsTools.hpp"
#include "LSDRaster.hpp"
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
// This constructor runs with the path and filename
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::create(string pname, string fname)
{
  pathname = pname;
  check_pathname_for_slash();
  param_fname = fname;
  
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
  string parameter, value, lower, lower_val;
  string bc;
	
  cout << "Parameter filename is: " << full_name << endl;
	
  // now ingest parameters 
  while (infile.good())
  {
    parse_line(infile, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
    {
      lower[i] = tolower(parameter[i]);
    }
	
    cout << "parameter is: " << lower << " and value is: " << value << endl;

    if 	(lower == "dem read extension")		dem_read_extension = value;
    else if (lower == "dem write extension")		dem_write_extension = value;
    else if (lower == "write path")		write_path = value;
    else if (lower == "write fname")		write_fname = value;
    else if (lower == "read path")		read_path = value;
    else if (lower == "read fname")		read_fname = value; 
 
    // parameter for fill
    else if (lower == "min_slope_for_fill")	
    {
      float_parameters["min_slope_for_fill"] = atof(value.c_str());
    }
    
    // parameters for flow info
    else if (lower == "boundary conditions")
    {
      // get the boundary value in lowercase
      lower_val = value;
		  for (unsigned int i=0; i<value.length(); ++i)
			  lower_val[i] = tolower(value[i]);
      vector<string> temp_bc(4);
      bc 	= lower_val;
      
      // now loop through collecting boundary conidtions
      for (int i = 0; i<4; i++)
      {
    	  string this_bc = bc.substr(i,1);
    	  cout << "Component " << i << " of the bc string: " << this_bc << endl;
        if (this_bc.find("p") != 0 && this_bc.find("p") != 0 && this_bc.find("n") != 0)
        {
          cout << "boundary condition not periodic, baselevel or noflux!" << endl;
          cout << "defaulting to no flux" << endl;
          temp_bc[i] = "n"; 
        }
        else
        {
          temp_bc[i] = this_bc;
        }
      }
      boundary_conditions = temp_bc;
    } 
       
    // parameters for chi analysis
    else if (lower == "nodeindex fname for chi map")		
    {
      support_file_names["nodeindex_fname_for_chi_map"] = atof(value.c_str());
    }
    else if (lower == "A_0")	
    {
      float_parameters["A_0"] = atof(value.c_str());
    }
    else if (lower == "m_over_n")	
    {
      float_parameters["m_over_n"] = atof(value.c_str());
    }  
    else if (lower == "threshold_area_for_chi")	
    {
      float_parameters["threshold_area_for_chi"] = atof(value.c_str());
    }
    
    // parameters for slope fitting
    else if (lower == "polyfit_window_radius")	
    {
      float_parameters["polyfit_window_radius"] = atof(value.c_str());
    }
            
    // now the various functions
    else if (lower == "write fill")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_fill"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
    } 
    else if (lower == "write nodeindex")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_nodeindex"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_flowinfo"] = temp_bool;
      raster_switches["need_nodeindex"] = temp_bool;
    } 		    
    else if (lower == "write chi map")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_chi_map"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_flowinfo"] = temp_bool;
      raster_switches["need_chi_map"] = temp_bool;
    } 	       
    else
    {
    	cout << "Line " << __LINE__ << ": No parameter '" 
           << parameter << "' expected.\n\t> Check spelling." << endl;
    }
  }
  infile.close();
	
  check_file_extensions_and_paths();
	
  compute_rasters_from_raster_switches();
  write_rasters_from_analysis_switches();

	
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a wrapper function that loops through the maps of raster switches
// gets the desired rasters, and then prints where necessary
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::compute_rasters_from_raster_switches()
{
  // read the base raster
  if(raster_switches.find("need_base_raster") == raster_switches.end())
  {
    cout<<"I need to read the base raster!!!"<<endl;
    read_base_raster();
  }
  
  // get the fill raster  
  if(raster_switches.find("need_fill") == raster_switches.end())
  {
    cout<<"I need to compute fill!!!"<<endl;
    
    // check to see if the base raster is loaded
    if(raster_indices.find("base_raster") == raster_indices.end())
    {
      cout << "Base raster hasn't been loaded. Loading it now." << endl;
      read_base_raster();
      
      // now the base raster is in the map
      fill_raster();  
    }
    else
    {
      fill_raster();  
    }
  } 
  
  // get the flow info
  if(raster_switches.find("need_flowinfo") == raster_switches.end())
  { 
    // only calculate flow info if it has not already been calculated
    if (not got_flowinfo)
    {  
      calculate_flowinfo();
    }
  }  
  
  // check to see if you need the nodeindex raster
  if(raster_switches.find("need_nodeindex") == raster_switches.end())
  {
    // check to see if it has already been calculated
    if(raster_indices.find("nodeindex") == raster_indices.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_nodeindex();
    }
  }  

  // check to see if you need the chi raster
  if(raster_switches.find("need_chi_map") == raster_switches.end())
  {
    // check to see if it has already been calculated
    if(raster_indices.find("chi_map") == raster_indices.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_chi_map();
    }
  }  
  
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Reads the base raster for analysis
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::write_rasters_from_analysis_switches()
{
  if(analyses_switches.find("write_fill") == raster_switches.end())
  {
    if(raster_indices.find("fill") == raster_indices.end()) 
    {
      cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string fill_seperator = "_fill";
    string fill_fname = write_path+write_fname+fill_seperator;

    int fill_index = raster_indices["fill"];
    vector_of_LSDRasters[fill_index].write_raster(fill_fname,dem_write_extension);
  }

  if(analyses_switches.find("write_nodeindex") == raster_switches.end())
  {
    if(raster_indices.find("nodeindex") == raster_indices.end()) 
    {
      cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string NI_seperator = "_NI";
    string NI_fname = write_path+write_fname+NI_seperator;

    int NI_index = raster_indices["nodeindex"];
    vector_of_LSDRasters[NI_index].write_raster(NI_fname,dem_write_extension);
  }

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    

  


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Reads the base raster for analysis
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::read_base_raster()
{
  // check to see if you've already got the base raster
  if(raster_indices.find("base_raster") == raster_indices.end())
  {
    // you don't have it. Calculate it here. 
    string full_raster_name = read_path+read_fname;
    cout <<"Reading the raster: " << endl;
    cout << full_raster_name + "." + dem_read_extension << endl;
    LSDRaster BaseRaster(full_raster_name,dem_read_extension);
    
    int RV_size = vector_of_LSDRasters.size();
    int base_raster_index = RV_size;
    vector_of_LSDRasters.push_back(BaseRaster);
    raster_indices["base_raster"] =  base_raster_index;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the fill DEM
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::fill_raster()
{
  // first check to make sure the base raster exists
  if(raster_indices.find("base_raster") == raster_indices.end())
  {
    cout << "Base raster doesn't exist! Reading it now." << endl;
    read_base_raster();
  }
  
  // see if the min_slope_for_fill has been initialised
  // If not,  set to default. 
  if(float_parameters.find("min_slope_for_fill") == float_parameters.end())
  {
    float_parameters["min_slope_for_fill"] = 0.0001;  
  }
  
  // now run the fill
  // first check to see if it has already been calculated
  if(raster_indices.find("fill") == raster_indices.end())
  {
    int base_raster_index =  raster_indices["base_raster"];
    LSDRaster temp_fill = vector_of_LSDRasters[base_raster_index].fill( float_parameters["min_slope_for_fill"] );
  
    int n_rasters = int(vector_of_LSDRasters.size());
    int fill_raster_index = n_rasters;
    vector_of_LSDRasters.push_back(temp_fill);
    raster_indices["fill"] =  fill_raster_index;    
  }
  else
  {
    cout << "Fill raster exists, its index is " << raster_indices["fill"] << endl;
  }
  
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the flow info object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_flowinfo()
{
  // first check to see if this already exists
  if (not got_flowinfo)
  {
    // it doens't exist. Calculate it here. 
  
    // this requires the fill raster. See if it exists
    if(raster_indices.find("fill") == raster_indices.end())
    {
      // it doesn't exit. Calculate it. 
      fill_raster();  
    }
    
    // now the fill exists. Check to see if boundary conditions exist
    check_boundary_conditions();
    
    // okay, everything should be ready for flow info calculation
    int fill_index =  raster_indices["fill"];
    LSDFlowInfo temp_FI(boundary_conditions, vector_of_LSDRasters[fill_index]);
    
    // set the data members
    FlowInfo = temp_FI;
    got_flowinfo = true;
  } 
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the nodeindex object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_nodeindex()
{  
  // see if you've already got this raster
  if(raster_indices.find("nodeindex") == raster_indices.end())
  {
    // you don't have it. Calculate it here.   
    // this requires the flow info object. See if it exists
    if(not got_flowinfo)
    {
      // it doesn't exit. Calculate it. 
      calculate_flowinfo(); 
    }
  
    // now you have the LSDFlowInfo object. Spit out the nodeindex raster
    LSDIndexRaster temp_NI = FlowInfo.write_NodeIndex_to_LSDIndexRaster();
    
    int n_IR = int(vector_of_LSDIndexRasters.size());
    vector_of_LSDIndexRasters.push_back(temp_NI);
    raster_indices["nodeindex"] =  n_IR;  
  }  
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates chi map
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_chi_map()
{  
  // see if you've already got this raster
  if(raster_indices.find("chi_map") == raster_indices.end())
  {
    // you don't have it. Calculate it here.   
    // this requires the flow info object. See if it exists
    if(not got_flowinfo)
    {
      // it doesn't exit. Calculate it. 
      calculate_flowinfo(); 
    }
    
    // now see if you've got the parameters
    if(float_parameters.find("A_0") == float_parameters.end())
    {
      // you don't have A_0. Replace with default
      float_parameters["A_0"] = 1000;
    }
    if(float_parameters.find("m_over_n") == float_parameters.end())
    {
      // you don't have m_over_n. Replace with default
      cout << "You haven't assigned m/n. Defaulting to 0.45" << endl;
      float_parameters["m_over_n"] = 0.45;
    }    
    if(float_parameters.find("threshold_area_for_chi_map") == float_parameters.end())
    {
      // you don't have m_over_n. Replace with default
      cout << "You haven't assigned the threshold_area_for_chi_map."
           << "Defaulting to 0" << endl;
      float_parameters["threshold_area_for_chi"] = 0;
    }   
      
    // now you have the LSDFlowInfo object. Spit out the chi raster
    LSDRaster temp_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes( 
                                 float_parameters["m_over_n"], 
                                 float_parameters["m_over_n"],
                                 float_parameters["threshold_area_for_chi"]);
    
    int n_chi = int(vector_of_LSDIndexRasters.size());
    vector_of_LSDRasters.push_back(temp_chi);
    raster_indices["chi_map"] =  n_chi;  
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
// this is a private function that makes sure the path has a slash at the end
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_boundary_conditions()
{
  if( int(boundary_conditions.size()) != 4)
  {
    cout << "Boundary conditions not assigned! Defaulting to no flux."  << endl;
    vector<string> temp_bc(4);
    for (int i = 0; i< 4; i++)
    {
      temp_bc[i] = "n";
    }
    boundary_conditions = temp_bc;
  }
   
  for (int i =0; i< 4; i++)
  {
    cout << "Boundary["<<i<<"]: "<<boundary_conditions[i]<< endl;
  }
  
  
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function checks the file extensions for reading and writing DEMs
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_file_extensions_and_paths()
{

  // first check the extensions
  if (dem_read_extension != "asc"  && dem_read_extension != "flt" && 
      dem_write_extension != "asc"  && dem_write_extension != "flt")
  {
    cout << "Raster file extension not assigned! Defaulting to flt format." << endl;
    dem_read_extension = "flt";
    dem_write_extension = "flt";
  }
  else
  {
    if (dem_read_extension != "asc"  && dem_read_extension != "flt")
    {
      //cout << "DEM read extension not assigned, defaulting to write extension." << endl;
      dem_read_extension = dem_write_extension;
    }
    else
    {
      //cout << "DEM write extension not assigned, defaulting to read extension." << endl;
      dem_write_extension = dem_read_extension;         
    }   
  }
  
  // now check the paths
  //cout << "Write path length is: " << write_path.length() << endl;
  if (write_path.length() == 0)
  {
    write_path = pathname;
    if (read_path.length() != 0)
    {
      write_path = read_path;
    }    
  }
  if (write_fname.length() == 0)
  {
    if (read_fname.length() != 0)
    {
      write_fname = read_fname;
    }
    write_fname = get_string_before_dot(param_fname);
    //cout << "Write fname not assigned, defaulting to name of parameter file." << endl;
    //cout << "The write fname is: " << write_fname << endl;
  }

  // now check the path
  //cout << "Read path length is: " << read_path.length() << endl;
  if (read_path.length() == 0)
  {
    read_path = write_path;
  }
  if (read_fname.length() == 0)
  {
    read_fname = get_string_before_dot(param_fname);
    //cout << "Read fname not assigned, defaulting to name of parameter file." << endl;
    //cout << "The read fname is: " << read_fname << endl;
  }
  
  cout << "The full read fname is:\n " << read_path+read_fname << endl;
  cout << "The full write fname is:\n " << write_path+read_fname << endl;
  cout << "The read and write extensions are: " << dem_read_extension 
       << " " << dem_write_extension << endl;
    
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
