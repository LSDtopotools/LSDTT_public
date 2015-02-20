//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Basinwide_cosmogenic_analysis.cpp
// A program to compute basinwide cosmogenic erosion rates.
//
// The function take two arguments to main
// The first is the data folder. This is the folder in which the parameter
// files are stored
// The second is the prefix of the parameter files
// 
// Developed by:
//  Simon M. Mudd, University of Edinburgh, School of GeoSciences
//  Stuart W.D. Grieve, University of Edinburgh, School of GeoSciences
//  Marie-Alice Harel, University of Edinburgh, School of GeoSciences
//  Martin D. Hurst, British Geological Survey
// 
//
// Copyright (C) 2015 Simon M. Mudd 2015
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDRaster.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDStrahlerLinks.hpp"
#include "../LSDBasin.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDCosmoData.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
  // the driver version
  string driver_version = "Driver_version: 0.01";

  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Basinwide cosmogenic analysis tool!  ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of the parameter files." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "There must be two parameter files in the named path." << endl;
    cout << "The first MUST have the extension .CRNfiles" << endl;
    cout << " This contains the filenames (including full path) of" << endl;
    cout << " the data file containing the cosmogenic data and all of the " << endl;
    cout << " names of the DEMs (without file extensions) to be analysed" << endl;
    cout << "The second file contains parameters for the analyses" << endl;
    cout << "  and MUST have the extension CRNparam" << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "For example if you parameter files are in home/fieldwork/Chile/CRN/" << endl;
    cout << " in this folder you will have two parameter files, for example" << endl;
    cout << " My_data.CRNfiles and My_data.CRNparam" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./Basinwide_CRN.exe home/fieldwork/Chile/CRN/ My_data" << endl;
    cout << "In windows (the slash directions will change and there is no leading ./)" << endl;
    cout << "Basinwide_CRN.exe c:\\fieldwork\\Chile\\CRN\\ My_data" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on cosmo data file format, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string param_name_prefix = argv[2];
  
  // remove control characters from these strings
  path_name.erase(remove_if(path_name.begin(), path_name.end(), ::iscntrl), path_name.end());
  param_name_prefix.erase(remove_if(param_name_prefix.begin(), param_name_prefix.end(), ::iscntrl), param_name_prefix.end());
  
  // load the file that contains the path to both the cosmo data and the 
  // DEMs
  string file_list_ext = ".CRNfiles";
  string file_list_fname = path_name+param_name_prefix+file_list_ext;

  // open the file and get both the cosmo data name and the list of DEMs
  ifstream file_list_in;
  file_list_in.open(file_list_fname.c_str());
 
  // load the name of the cosmo data
  string cosmo_data_name;
  file_list_in >> cosmo_data_name;
  
  //  eliminate control characters and find out if it is a csv file
  cosmo_data_name.erase(remove_if(cosmo_data_name.begin(), cosmo_data_name.end(), 
                                              ::iscntrl), cosmo_data_name.end());
  string csv_ext = ".csv";
  size_t found = cosmo_data_name.find(csv_ext);
  string filetype = "txt";
  if (found != string::npos)
  {
    filetype = "csv";
  }

  // now load the names of the DEMs
  vector<string> DEM_file_list;
  string this_file;
  while(file_list_in >> this_file)
  {
    // get rid of control characters
    this_file.erase(remove_if(this_file.begin(), this_file.end(), 
                                            ::iscntrl), this_file.end());
    DEM_file_list.push_back(this_file);
  }
  int N_DEMS = int(DEM_file_list.size());

  // now load the CRNCosmoData object
  LSDCosmoData CosmoData(cosmo_data_name,filetype); 

  cout << "===========================================================" << endl;
  cout << "Welcome to the Basinwide cosmogenic analysis tool" << endl;
  cout << "This software was developed at the University of Edinburgh," << endl;
  cout << "by the Land Surface Dynamics group. For questions email" << endl;
  cout << "simon.m.mudd _at_ ed.ac.uk" << endl;
  cout << "This software is released under a GNU public license." << endl;
  cout << "You are using " << driver_version << endl;
  cout << "===========================================================" << endl;
  cout << "Your cosmogenic data is below." << endl; 
  cout << "If this looks wrong, check your filenames." << endl;
  CosmoData.print_data_to_screen();
  cout << "You are about to look through these DEMs for cosmo basins:" << endl;
  cout << "..........................................................." << endl;
  for (int dem = 0; dem < N_DEMS; dem++)
  {
    cout << DEM_file_list[dem] << endl;
  }
  cout << "..........................................................." << endl;
  cout << "++IMPORTANT++ There must be ENVI bil files with these names." << endl;
  cout << "ENVI bil files are required because, unlike asc or flt files, " << endl;
  cout << "they use georeferencing information, which is used in the analyses." << endl;
  cout << "For more information about changing DEM formatting, see: " << endl;
  cout << "http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl;
  cout << "===========================================================" << endl;

  // load the parameter file
  string param_ext = ".CRNparam";
  string param_fname = path_name+param_name_prefix+param_ext;

  // a boundary condition for the flow info object
  // get the flow info 
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "n";
  boundary_conditions[1] = "n";
  boundary_conditions[2] = "n";
  boundary_conditions[3] = "n";
 
  // The default slope parameter for filling. Do not change. 
  float min_slope = 0.0001;
  
  // parameters for making stream networks and looking for channels
  int source_threshold = 12;
  int search_radius_nodes = 1;
  int threshold_stream_order = 1;

  // the values of theta and phi step are based on testing by S. Grieve 
  // Note that Codilian reccomends 5,5 but 10,15 leads to minimal errors
  int theta_step = 30;
  int phi_step = 30;
  
  // density of rock. Most calculations do this in shielding depth but 
  // printing results to screen includes a transaformation to length per time
  // units
  double rho = 2650; 

  // extensions for filenames
  string fill_ext = "_fill";
  string DEM_bil_extension = "bil";
  
  // now go through the list of DEMs, extracting the data
  // One of the main motivations for doing this in a loop is to not store the 
  // data for each analysis
  for (int dem = 0; dem <N_DEMS; dem++)
  {

    // load the DEM
    LSDRaster topo_test(DEM_file_list[dem], DEM_bil_extension);
    
    // fill the raster
    LSDRaster filled_raster = topo_test.fill(min_slope);
    
    // get the flow info
    LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);

    // get contributing pixels (needed for junction network)
    LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    // get the sources
    vector<int> sources;
    sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);

    // now get the junction network
    LSDJunctionNetwork JNetwork(sources, FlowInfo);
  
    // get the topographic shielding
    cout << "Starting topogrpahic shielding" << endl;
    LSDRaster T_shield = filled_raster.TopographicShielding(theta_step, phi_step);
    
    // transform the cosmo data to this DEM UTM zone
    CosmoData.convert_to_UTM(filled_raster);
    
    // get the cosmogenic erosion rates
    CosmoData.basic_cosmogenic_analysis(search_radius_nodes, 
                            threshold_stream_order, filled_raster,
                            T_shield, FlowInfo, JNetwork);
    
    CosmoData.print_simple_results_to_screen(rho);       
    
  }
  
  

  
}
  