//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Production_comparison.cpp
// A program to compare the results from various production schemes
//
// There are no arguments to main
// 
// Developed by:
//  Simon M. Mudd, University of Edinburgh, School of GeoSciences
//  Stuart W.D. Grieve, University of Edinburgh, School of GeoSciences
//  Marie-Alice Harel, University of Edinburgh, School of GeoSciences
//  Martin D. Hurst, British Geological Survey
// 
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
#include "../LSDShapeTools.hpp"
#include "../LSDCosmoData.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
  // the driver version
  string driver_version = "Driver_version: 0.02";

  cout << "================================================================" << endl;
  cout << "|| Welcome to the production comparison tool!                 ||" << endl;
  cout << "|| This tool prints some files so you can compare production  ||" << endl;
  cout << "|| rates of the cosmogenic nuclide 10Be from muons.           ||" << endl;
  cout << "================================================================" << endl;

  // filename for the muon production comparison
  string filename = "Muon_production_comparison.csv";
  
  // the atmospheric data is in the folder with the driver_functions
  string path_to_atmospheric_data = "./";
  
  // print the muon comparison file
  LSDCRNParameters LSDCRNP;
  LSDCRNP.Print_10Beproduction_csv(filename, path_to_atmospheric_data);
  

  
    // initiate a particle. We'll just repeatedly call this particle
  // for the sample. 
  int startType = 0; 
  double Xloc = 0;
  double Yloc = 0;
  double  startdLoc = 0.0;
  double  start_effdloc = 0.0;
  double startzLoc = 0.0;
  
  // create a particle at zero depth
  LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                               startdLoc, start_effdloc, startzLoc);
  
  // These are parameters for the CRONUS emulator
  // get some parameters for Stone production
  vector<double> Prefs_st = LSDCRNP.get_Stone_Pref();
  
  // we are going to print these at a fixed location
  double site_lat = 70;
  double site_lon = 0;
  double site_elev = 0;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();
  double pressure = LSDCRNP.NCEPatm_2(site_lat, site_lon, site_elev);  
  
  // get the Stone scaling
  double Fsp = 1.0;     // for the initial guess we don't adjust Fsp (as in CRONUS)
  double stoneP = LSDCRNP.stone2000sp(site_lat, pressure, Fsp);
  //cout << "LINE 1718 Lat: " << lat << " pressure: " << pressure << " stone: " << stoneP << endl;
  
  // retrieve the pre-scaling factors
  double P_ref_St_10 = Prefs_st[0];
  double P_ref_St_26 = Prefs_st[1];
  //cout << "P ref stone 10: " <<  P_ref_St_10 <<  " P_ref_St_26: " << P_ref_St_26 << endl;
  
  // precalculate the P_mu vectors
  vector<double> z_mu;
  vector<double> P_mu_z_10Be;
  vector<double> P_mu_z_26Al;
  double effective_dLoc = 0;      // we look at a particle at the surface. 
  LSDCRNP.get_CRONUS_P_mu_vectors(pressure, effective_dLoc, z_mu, 
                                  P_mu_z_10Be,P_mu_z_26Al);
                                  
  // get the scaled production rates
  double P_sp_10Be = P_ref_St_10*stoneP;
  double P_sp_26Al = P_ref_St_26*stoneP;
  double thickSF = 1;     // no self shielding, the scaling factor is 1
  double N_Be10, N_Al26;
  
  double N_Schaller, N_Granger, N_Braucher, N_newCRONUS, N_newCRONUSShielded;
  
  // now create an erosion gradient, with particles. 
  // These are in g/cm^2/yr
  int N_Erate = 101;
  double erate_begin = 0.001;
  double erate_end = 0.2;
  double derate = (erate_end-erate_begin)/N_Erate;
  vector<double> erate_vec;
  double this_erate; 
  
  // this tells the scaling routine to on;y scale for 10Be;
  vector<bool> nuclide_scaling_switches(4,false);
  nuclide_scaling_switches[0]= true;
  
  ofstream N_out;
  N_out.open("Nuclide_Conc_vs_Erate.csv");
  N_out << "Nuclide concetrations as a function of erosion rates for different production schemes" << endl;
  N_out << "All concentrations are in atoms/g" << endl;
  N_out << "Erate(g/cm^2/yr),CRONUS,Schaller,Granger,Braucher,newCRONUS," 
        << "logN_CRONUS,logN_Schaller,logN_Granger,logN_Braucher,logN_newCRONUS,N_newCronus_0p9shielded" << endl;
  
  // Loop through the erosion rates getting the nuclude concentrations
  for (int i = 0; i<N_Erate; i++)
  {
    this_erate =  double(i)*derate+erate_begin;
    erate_vec.push_back(this_erate);
    
    // first the CRONUS concentration
    eroded_particle.CRONUS_calculate_N_forward(this_erate, LSDCRNP, z_mu, P_mu_z_10Be, 
                               P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                               N_Be10,  N_Al26);

    LSDCRNP.set_Schaller_parameters();
    LSDCRNP.scale_F_values(stoneP,nuclide_scaling_switches);
    eroded_particle.update_10Be_SSfull(this_erate,LSDCRNP);
    N_Schaller = eroded_particle.getConc_10Be();
    
    LSDCRNP.set_Granger_parameters();
    LSDCRNP.scale_F_values(stoneP,nuclide_scaling_switches);
    eroded_particle.update_10Be_SSfull(this_erate,LSDCRNP);
    N_Granger = eroded_particle.getConc_10Be();
    
    LSDCRNP.set_Braucher_parameters();
    LSDCRNP.scale_F_values(stoneP,nuclide_scaling_switches);
    eroded_particle.update_10Be_SSfull(this_erate,LSDCRNP);
    N_Braucher = eroded_particle.getConc_10Be();
    
    LSDCRNP.set_newCRONUS_parameters();
    LSDCRNP.scale_F_values(stoneP,nuclide_scaling_switches);
    eroded_particle.update_10Be_SSfull(this_erate,LSDCRNP);
    N_newCRONUS = eroded_particle.getConc_10Be();
    
    LSDCRNP.set_newCRONUS_parameters();
    LSDCRNP.scale_F_values(stoneP*0.8,nuclide_scaling_switches);
    eroded_particle.update_10Be_SSfull(this_erate,LSDCRNP);
    N_newCRONUSShielded = eroded_particle.getConc_10Be();
    
        
    N_out << this_erate << "," << N_Be10 << "," << N_Schaller << "," << N_Granger
          << "," << N_Braucher << "," << N_newCRONUS <<","
          << log(N_Be10) << "," << log(N_Schaller) << "," << log(N_Granger) 
          << "," << log(N_Braucher) << "," << log(N_newCRONUS) << "," << N_newCRONUSShielded <<  endl;  
  }

  N_out.close();
  

  // error checking
  this_erate = 0.001;
  eroded_particle.CRONUS_calculate_N_forward(this_erate, LSDCRNP, z_mu, P_mu_z_10Be, 
                               P_mu_z_26Al, thickSF, P_sp_10Be, P_sp_26Al,
                               N_Be10,  N_Al26);
  cout << "Erate is: " << this_erate << " and N_Be10 is: " << N_Be10 << endl;


  // now get the cronus emulator
  double rho = 2650;
  double topo_scale =1;
  double snow_scale = 1;
  double N_26Al = 2690000;
  double sample_del26 =  88000 ;
  double N_10Be_test = 176015;
  double sample_del10 = 13000;
  vector<double> erateinfo = eroded_particle.CRONUS_get_Al_Be_erosion(LSDCRNP, pressure,
                      site_lat, rho, N_10Be_test, N_26Al,sample_del10, sample_del26,
                      topo_scale,  snow_scale);
                      
  cout << "Effective erate is: " << erateinfo[0] << endl;
}
  