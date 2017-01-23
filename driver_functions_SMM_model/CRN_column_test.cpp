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

  if (argc <= 4)
  {
    cout << "=============================================================\n"
         << "Welcome to a test function for LSDCRN columns." << endl
         << "These columns are compased of LSDCRN particles and can be used" << endl
         << "to investigate how CRN concentrations in detrital sediments" << endl
         << "might vary under different erosion scenarios." << endl
         << "The program takes 6 arguments:" << endl
         << "The first argument is the periodicity " << endl
         << "The second argument is the uplift rate " << endl
         << "The third argumnt is the amplitude fraction, which is the ratio" << endl
         << "  of the erosion amplitude to the uplift rate" << endl
         << "The fourth argument is the timestep" << endl
         << "The fifth agument is the particle spacing (suggest 0.1)" << endl
         << "The sixth argument is the start depth (suggest 3) " << endl
         << "The seventh argument is a switch: " << endl
         << "    0 == sinusoidal" << endl
         << "    1 == square wave" << endl
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

  double dt = atof(argv[4]);
  cout << "dt: " << dt << endl;
  
  double particle_spacing = atof(argv[5]);
  cout << "particle spacing: " << particle_spacing << endl;
  
  double startDepth = atof(argv[6]);
  cout << "startDepth: " << startDepth << endl;

  int Vary_Type = atoi(argv[7]);
  cout << "Vary_Type: " << Vary_Type << endl;

  // Set up the LSDCRN parameters
  int startType = 0;
  //double startDepth = 3;
  //double particle_spacing = 0.1;
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
  //double dt = 100;
  double t_ime = 0;
  double  end_time = periodicity*4;
  double print_interval = periodicity/20;
  int pi = int(print_interval/dt);
  
  double startxLoc = 0;
  double startyLoc = 0;
  
  // the minium and maximum erosion rates and apparent erosion rates
  double min_erate = 1000;
  double max_erate = -1000;
  double min_app_erate = 1000;
  double max_app_erate = -1000;
  
  // now run forward in time
  int step = 0;
  while(t_ime < end_time)
  {
    // replace zeta_old
    zeta_old = zeta_at_cell;
    
    if (Vary_Type == 0)
    {
      // get the erosion rate (which is periodic)
      erosion_rate = erosion_amp*sin( (t_ime)*2*M_PI/periodicity )+uplift;
    }
    else if (Vary_Type == 1)
    {
      // a stupid inefficient way to get a square wave
      double test_er = erosion_amp*sin( (t_ime)*2*M_PI/periodicity );
      if(test_er >= 0)
      {
        erosion_rate = uplift+2*uplift*amplitude_fraction;
      }
      else
      {
        erosion_rate = uplift-2*uplift*amplitude_fraction;
      }
    }
    else
    {
      // get the erosion rate (which is periodic)
      erosion_rate = erosion_amp*sin( (t_ime)*2*M_PI/periodicity )+uplift;
    }
    
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
    
    // get the minimum and maximum erates   
    if(erosion_rate > max_erate)
    {
      max_erate = erosion_rate;
    }
    if(erosion_rate < min_erate)
    {
      min_erate = erosion_rate;
    }
    if(this_app_erosion[0] > max_app_erate)
    {
      max_app_erate = this_app_erosion[0];
    }
    if(this_app_erosion[0] < min_app_erate)
    {
      min_app_erate = this_app_erosion[0];
    }    
    
       
    if(step%pi == 0)
    {
      cosmo_data_out << t_ime << "," << erosion_rate << "," << uplift << "," 
                     << this_app_erosion[0] << "," << this_app_erosion[1] << ","
                     << this_app_erosion[2] << endl; 
    }
    
    // increment time
    t_ime += dt;
    step++; 
    
  }
  cosmo_data_out.close();
  
  cout << endl << endl << "------------------------------------------------" << endl;
  cout << "dt is: " << dt << " and U in mm/yr is: " << uplift*1000 << endl;
  cout << "min erate diff in mm/yr: " << 1000*(min_app_erate - min_erate) << endl
       << " and err: " << 1000*(min_app_erate - min_erate)/(uplift*1000) << endl;
  cout << "max erate diff in mm/yr: " << 1000*(max_app_erate - max_erate) << endl
       << " and err: " << 1000*(max_app_erate - max_erate)/(uplift*1000) << endl;
  cout << "------------------------------------------------" << endl;
  
  return 0;
}
