//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// CRN_stochastic_river.cpp
//
// This program calculates CRN accumulation as particles make their way
// down a river. Uses a brute force stochastic approach
//
// Developed by:
//  Simon M. Mudd
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "../LSDRaster.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDStrahlerLinks.hpp"
#include "../LSDBasin.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDCosmoData.hpp"
using namespace std;



// velocity is in metres per year  
double EvolveParticle(double start_concentration, double velocity,
                    vector<double> Latitude, vector<double> Longitude, 
                    vector<double> Elevation, vector<double> FlowDistance)
{
  // initiate a particle. We'll just repeatedly call this particle
  // for the sample. 
  int startType = 0; 
  double Xloc = 0;
  double  startdLoc = 0.0;
  double  start_effdloc = 0.0;
  double startzLoc = 0.0;
  
  int n_nodes = int(Latitude.size());
  
  double ThisFlowDistance = 0;
  double EndFlowDistance = FlowDistance[n_nodes-1];
  cout << "The ending flow distance is " << EndFlowDistance << endl;

  double null_conc = 0;
  
  // create a particle at zero depth
  LSDCRNParticle this_particle(startType, Xloc, startzLoc,
                               startdLoc, start_effdloc, null_conc,null_conc,start_concentration);

  cout << "The starting concentration is: " << this_particle.getConc_21Ne() << endl;

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;
  string path_to_atmospheric_data = "./";
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  
  // Update to the latest production rates
  LSDCRNP.set_newCRONUS_parameters();
  
  int current_node = 0;
  double Next_flow_distance = FlowDistance[1];
  bool at_end = false;
  double this_lat,this_long, this_Elev;
  double this_pressure;
  double Fsp = 1;  // default for stone scaling
  double erate = 0;
  double dt = 1;
  double this_time = 0;

  double production;
  
  while(not at_end)
  {
    // calculate production 
    this_lat = Latitude[current_node];
    this_long = Longitude[current_node];
    this_Elev = Elevation[current_node];
    this_pressure = LSDCRNP.NCEPatm_2(this_lat, this_long,this_Elev);
  
    production = LSDCRNP.stone2000sp(this_lat,this_pressure, Fsp);
    
    // update the scaling
    double no_shield = 1;
    //cout << "production is: " << production << endl;
    LSDCRNP.set_scaling(production, no_shield, no_shield);
    
    // update the concentration
    this_particle.update_21Ne_conc(dt, erate, LSDCRNP);
    
    //cout << "Concentration is: " << this_particle.getConc_21Ne() << endl;
    
    // now move the particle 
    ThisFlowDistance += velocity*dt;
    this_time+=dt;
    
    //cout << "New flow distance is: " << ThisFlowDistance << " and next flow distance is: " << Next_flow_distance << endl;
    
    if (ThisFlowDistance > Next_flow_distance)
    {
      current_node++;
      
      if (current_node == n_nodes)
      {
        at_end = true;
      }
      else
      {
        Next_flow_distance = FlowDistance[current_node];
      }
    }
  }
  
  // get the concentration
  double final_conc = this_particle.getConc_21Ne();
  cout << "The final concentration is: " << final_conc << " after: " << this_time << endl;
  return final_conc;

}

void ReadRiver(LSDSpatialCSVReader& CSVFile, string Elevation_name, string FlowDistance_name, bool NeedToReverse,
               vector<double>& Latitude, vector<double>& Longitude, 
               vector<double>& Elevation, vector<double>& FlowDistance)
{
  cout << "Got the data!" << endl;
  
  vector<double> latitude = CSVFile.get_latitude();
  vector<double> longitude = CSVFile.get_longitude();
  vector<float> elevation = CSVFile.data_column_to_float(Elevation_name);
  vector<float> flow_distance = CSVFile.data_column_to_float(FlowDistance_name);
  
  
  //cout << "Start FD = " << flow_distance[0] << " and later " << flow_distance[10] << endl;
  
  
  // the 
  if(NeedToReverse)
  {
    cout << "Let me reverse those vectors for you." << endl;
    reverse(latitude.begin(),latitude.end());
    reverse(longitude.begin(),longitude.end());
    reverse(elevation.begin(),elevation.end());
    reverse(flow_distance.begin(),flow_distance.end());
  }
  
  //cout << "Start FD = " << flow_distance[0] << " and later " << flow_distance[10] << endl;
  
  // prepare flow distance
  int n_fd = int(flow_distance.size());
  double final_FD =   flow_distance[0];
  for(int i = 0; i< n_fd; i++)
  {
    
    flow_distance[i] = final_FD-flow_distance[i];
  }
  
  //cout << "Start FD = " << flow_distance[0] << " and later " << flow_distance[10] << endl;
  
  
  // convert to doubles
  vector<double> FD;
  vector<double> E;
  for(int i = 0; i< n_fd; i++)
  {
    FD.push_back(double(flow_distance[i]));
    E.push_back(double(elevation[i]));
  }
  
  
  Latitude = latitude;
  Longitude = longitude;
  Elevation = E;
  FlowDistance = FD;
}


int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the stochastic river tool.               ||" << endl;
    cout << "|| It assumes constant velocity of a particle down a   ||" << endl;
    cout << "|| river channel and accumulates 21Ne.                 ||" << endl;
    cout << "|| An extremely basic accumulation that is only meant  ||" << endl;
    cout << "|| to give basic insight into particl transit times    ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Simon M. Mudd                                       ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./basin_averaging_tool.exe /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_BasinAvg.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << " http://lsdtopotools.github.io/LSDTT_book/#_chi_analysis_part_3_getting_chi_gradients_for_the_entire_landscape" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);
  
  // for the basin tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // set default float parameters
  int_default_map["UTM_zone"] = 13;

  
  // set default in parameter
  float_default_map["start_21Ne_conc"] = 8.82E6;
  
  // set default methods
  bool_default_map["is_North"] = true;
  bool_default_map["NeedToReverse"] = false;

  
  // set default string method
  string_default_map["channel_file"] = "test_channel.csv";
  string_default_map["Easting_column_name"] = "X";
  string_default_map["Northing_column_name"] = "Y";
  string_default_map["Elevation_column_name"] = "Elevation";
  string_default_map["FlowDistance_column_name"] = "FlowDistance";

  // Use the parameter parser to get the maps of the parameters required for the 
  // analysis

  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();
  
  // okay, we are going to load the channel first
  LSDSpatialCSVReader CSVFile((DATA_DIR+this_string_map["channel_file"]));
  //CSVFile.print_lat_long_to_screen();
  CSVFile.print_data_map_keys_to_screen();
  
  // set the UTM. 
  CSVFile.set_UTM_information(this_int_map["UTM_zone"], this_bool_map["is_North"]);
  
  // Convert X and Y to lat-long
  CSVFile.get_latlong_from_x_and_y(this_string_map["Easting_column_name"], this_string_map["Northing_column_name"]);

  vector<double> Latitude;
  vector<double> Longitude; 
  vector<double> Elevation;
  vector<double> FlowDistance;

  // read the river data
  ReadRiver(CSVFile, this_string_map["Elevation_column_name"], this_string_map["FlowDistance_column_name"], this_bool_map["NeedToReverse"],
            Latitude, Longitude, Elevation, FlowDistance);


  double Ne21_Conc = double(this_float_map["start_21Ne_conc"]);
  double velocity = 3;      // metres per year
  double de = 0.1;
  double start_vel = 0.1;
  
  string data_out_fname = OUT_DIR+OUT_ID+"_excess21Ne.csv";
  ofstream data_out;
  data_out.open(data_out_fname.c_str());
  data_out << "velocity (m/yr),excess 21Ne (atoms/gram)" << endl;
  
  
  for( int i = 0; i<100; i++)
  {
    velocity = double(i)*de+start_vel;
    cout << "velocity is: " << velocity << endl;
    
    double FinalConc =  EvolveParticle(Ne21_Conc, velocity,Latitude, Longitude, 
                                     Elevation, FlowDistance);
                                     
    data_out << velocity << "," << FinalConc- Ne21_Conc << endl;

  }
  data_out.close();

}
