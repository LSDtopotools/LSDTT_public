#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>
#include <fstream>
#include "../LSDParticleColumn.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
using namespace std;

bool file_check(string name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char *argv[])
{

  cout << "argc is: " << argc << endl;

  if (argc <= 3)
  {
    cout << "=============================================================\n"
         << "Welcome to a test function for LSDCRN columns." << endl
         << "These columns are compased of LSDCRN particles and can be used" << endl
         << "to investigate how CRN concentrations in detrital sediments" << endl
         << "might vary under different erosion scenarios." << endl
         << "The program takes 3 arguments:" << endl
         << "The first argument is the particles "
         << "=============================================================\n";
    exit(0);
  }

  double periodicity = atof(argv[1]);
  cout << "periodicity is: " << periodicity << endl;
  
  double uplift = atof(argv[2]);
  cout << "uplift_rate is: " << uplift << endl;
  
  double amplitude_fraction = atof(argv[3]);
  double erosion_amp = amplitude_fraction*2*uplift;
  cout << "amplitude fraction is: " << amplitude_fraction << " and erosion amplitude is: "
       << erosion_amp << endl;

  // Set up the LSDCRN parameters
  int startType = 0;
  int startDepth = 3;
  double particle_spacing = 0.1;
  LSDCRNParameters CRNParam;
  CRNParam.set_Neutron_only_parameters();     // set to neutron only production

  // initiate an empty column
  LSDParticleColumn temp_col;

  // set a default elevation
  double zeta_at_cell = 100;
  
  double this_x = 0; 
  double this_y = 0; 
  
  // the uplift rate in m/yr
  // double uplift = 0.00001;
  

  // the density
  double rho = 2000;
  
  // the effective uplift rate
  double effective_uplift = uplift*rho*0.1;
  
  // the erosion rate, and surface locations
  double erosion_rate;
  double zeta_old;
  
  // make a steady state column
  //cout << "LINE 3032, initiating SS column" << endl;
  temp_col.initiate_SS_cosmo_column_3CRN(startType, this_x, this_y, startDepth, particle_spacing, 
                                             zeta_at_cell, effective_uplift, CRNParam); 
  //cout << "LINE 3035, Initiated column" << endl;
      
  // make sure this is a rock only simulation          
  double ST = 0;          
  temp_col.set_SoilThickness(ST);
  
  
  // now open the files for printing
  ofstream cosmo_data_out;
  cosmo_data_out.open("Test_cosmo_column.csv");
  cosmo_data_out << "time,erosion_rate,uplift_rate,app_erosion_10Be,app_erosion_14C,app_erosion_21Ne" << endl;
  
  // set up the time 
  double dt = 100;
  double t_ime = 0;
  double  end_time = 200000;
  
  double startxLoc = 0;
  double startyLoc = 0;
  
  // now run forward in time
  while(t_ime < end_time)
  {
    // replace zeta_old
    zeta_old = zeta_at_cell;
    
    // get the erosion rate (which is periodic)
    erosion_rate = erosion_amp*sin( (t_ime)*2*M_PI/periodicity )+uplift;
    
    // get the new zeta
    zeta_at_cell = zeta_at_cell-erosion_rate*dt+dt*uplift;
    
    // now update the column.
    LSDParticleColumn this_eroded_column = 
                 temp_col.update_CRN_list_rock_only_eros_limit_3CRN(
                       dt, uplift,
                       startType, startDepth,
                       startxLoc,  startyLoc,
                       zeta_old,zeta_at_cell,
                       particle_spacing, CRNParam);

    // Now get the apparent erosion rates
    vector<double> this_app_erosion = 
       temp_col.calculate_app_erosion_3CRN_neutron_rock_only(CRNParam);
    
    cosmo_data_out << t_ime << "," << erosion_rate << "," << uplift << "," 
                   << this_app_erosion[0] << "," << this_app_erosion[1] << ","
                   << this_app_erosion[2] << endl;
    
    // increment time
    t_ime += dt; 
    
    
  
  }
  cosmo_data_out.close();
  
  return 0;
}
