//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDCRNParameters.hpp
//
// Land Surface Dynamics Cosmogenic Radionuclide Parameters Object
//
// This keeps track of paramters used to calculate the evolution of 
// in situ cosmogenic nuclides. 
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
#include <map>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
#include "LSDCRNParameters.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDCRNParameters_CPP
#define LSDCRNParameters_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// the LSDCRNParameters object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// function to set CRN parameters
void LSDCRNParameters::create()
{
  S_t = 1;

  // from Vermeesh 2007
  lambda_10Be = 456e-9;		// in yr-1
  lambda_26Al = 980e-9;		// in yr-1
  lambda_14C = 121e-6;		// in yr-1
  lambda_36Cl = 230e-8;		// in yr-1

  // from Vermeesh 2007
  P0_10Be = 5.11;					// in a/g/yr
  P0_26Al = 30.31;				// in a/g/yr
  P0_14C = 5.86;					// in a/g/yr
  P0_36Cl = 55.45;				// in a/g/yr
  P0_21Ne = 20.29;				// in a/g/yr
  P0_3He = 97.40;					// in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 0.9724;
  F_10Be[1] = 0.0186;
  F_10Be[2] = 0.004;
  F_10Be[3] = 0.005;

  // dimensionless
  F_26Al[0] = 0.9655;
  F_26Al[1] = 0.0233;
  F_26Al[2] = 0.005;
  F_26Al[3] = 0.0062;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.0691;
  F_14C[2] = 0.0809;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.903;
  F_36Cl[1] = 0.0447;
  F_36Cl[2] = 0.05023;
  F_36Cl[3] = 0.0;
}

// this function gets the parameters used to convert elevation to 
// pressure
void LSDCRNParameters::load_parameters_for_atmospheric_scaling(string path_to_data)
{
  cout.precision(8);
  
  // first load the levels
  levels.push_back(1000);
  levels.push_back(925);
  levels.push_back(850);
  levels.push_back(700);
  levels.push_back(600);
  levels.push_back(500);
  levels.push_back(400);
  levels.push_back(300);
  
  // the dimensions of the data
  int n_levels = 8;
  int NRows = 73;
  int NCols = 145;
  Array2D<double> new_slp(NRows,NCols,0.0);
  Array2D<double> new_meant(NRows,NCols,0.0);
    
  // now load the mean sea level pressure
  string filename = "NCEP2.bin";
  filename = path_to_data+filename;
  cout << "Loading mean sea level, file is: " << endl << filename << endl;

  ifstream ifs_data(filename.c_str(), ios::in | ios::binary);
  if( ifs_data.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  
  

  double temp;
  cout << "The size of a double is: " << sizeof(temp) << endl;
  for (int i=0; i<NCols; ++i)
  {
    for (int j=0; j<NRows; ++j)
    {
      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      new_slp[j][i] = temp;
      //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    }
  }
  
  for (int i=0; i<NCols; ++i)
  {
    for (int j=0; j<NRows; ++j)
    {
      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      new_meant[j][i] = temp;
      //cout << "new_meant100["<<j+1<<"]["<<i+1<<"]: " << new_meant[j][i] << endl;
    }
  }  
  
  // now get the indices
  vector<double> temp_lat(NRows,0.0);
  for (int i=0; i<NRows; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    temp_lat[i] = temp;
    //cout << "Lat["<<i+1<<"]: " << temp_lat[i] << endl;
  }
  vector<double> temp_long(NCols,0.0);
  for (int i=0; i<NCols; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    temp_long[i] = temp;
    //cout << "Long["<<i+1<<"]: " << temp_long[i] << endl;
  }  
  
  ifs_data.close();
  
  
  // now the data with levels
  filename = "NCEP_hgt.bin";
  filename = path_to_data+filename;
  cout << "Loading hgt, file is: " << endl << filename << endl;

  ifstream ifs_data2(filename.c_str(), ios::in | ios::binary);
  if( ifs_data2.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  
  

  // get the gm heights
  vector< Array2D<double> > vec_hgt_gm_array;
  for (int lvl = 0; lvl < n_levels; lvl++)
  {
    Array2D<double> current_hgt_array(NRows,NCols,0.0);
    for (int i=0; i<NCols; ++i)
    {
      for (int j=0; j<NRows; ++j)
      {
        ifs_data2.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        current_hgt_array[j][i] = temp;
        //cout << "new_slp["<<lvl<<"]["<<j+1<<"]["<<i+1<<"]: " << current_hgt_array[j][i] << endl;
      }
    }
    //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    vec_hgt_gm_array.push_back(current_hgt_array.copy());
  }

  // now the gp heights
  vector< Array2D<double> > vec_hgt_gp_array;
  for (int lvl = 0; lvl < n_levels; lvl++)
  {
    Array2D<double> current_hgt_array(NRows,NCols,0.0);
    for (int i=0; i<NCols; ++i)
    {
      for (int j=0; j<NRows; ++j)            
      {
        ifs_data2.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        current_hgt_array[j][i] = temp;
        //cout << "new_slp["<<lvl<<"]["<<j+1<<"]["<<i+1<<"]: " << current_hgt_array[j][i] << endl;
      }
    }
    //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    vec_hgt_gp_array.push_back(current_hgt_array.copy());
  }
  ifs_data2.close();
  
  // now update the data elements
  NCEPlat = temp_lat;
  NCEPlon = temp_long;
  
  meanslp = new_slp.copy();
  meant1000 = new_meant.copy();
  gp_hgt = vec_hgt_gp_array;
  gm_hgt = vec_hgt_gp_array;  
  
  cout << "Size lat: " << NCEPlat.size() << " size long: " << NCEPlon.size() << endl;
  cout << "size slp:" << meanslp.dim1() << " " << meanslp.dim2() << endl;
  cout << "size t1000: " << meant1000.dim1() << " " << meant1000.dim2() << endl;
  

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function sets a number of paramters that are used for acalucaltion
// of scaling and error propigation
//
// They are constants used in the CRONUS caluclator, and have been ported
// from make_al_be_consts_v22.m written by Greg Balco
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_CRONUS_data_maps()
{
  map<string,double> temp_map;
  
  // 10Be decay: this is specific to the CRONUS calculator
  temp_map["l10"] = -log(0.5)/1.387e6; // Chmeleff/Korschinek value
  // Al-26 decay constant -- value compatible with Nishiizumi standards
  temp_map["l26"] = 9.83e-7;
  

  // Note that the uncertainty is not used in exposure-age or erosion-rate
  // calculators. Here only for development purposes
  double dldt = -log(0.5)*(1/(1.387e6*1.387e6));
  temp_map["dell10"] = sqrt((dldt*1/(0.012e6*0.012e6)));
  temp_map["dell26"] = 2.5e-8;

  //Effective attenuation length for spallation in rock
  // Commonly accepted value: 160 g/cm2
  // For discussion see Gosse and Phillips (2000)
  temp_map["Lsp"] = 160.0;

  // Fsp - fraction of total production by spallation rather than muons
  // For use with Lal/Stone scaling scheme in exposure age calculation
  // For details, see Stone (2000)
  // This aspect of Stone(2000) is de-emphasized in version 2. These constants
  // are only in use for historical comparisons and quick initial guesses for 
  // the exposure age and erosion rate solvers. 

  // Note that they are probably incorrect WRT Be-10 restandardization. Don't
  // use these for anything important. 
  temp_map["Fsp10"] = 0.978;
  temp_map["Fsp26"] = 0.974;

  // Be-10 standardization info.
  // Standards comparison/conversion lookup table
  temp_map["Be_std_07KNSTD"] = 1.0000;
  temp_map["Be_std_KNSTD"] = 0.9042;
  temp_map["Be_std_NIST_Certified"] = 1.0425;
  temp_map["Be_std_LLNL31000"] =  0.8761;
  temp_map["Be_std_LLNL10000"] = 0.9042;
  temp_map["Be_std_LLNL3000"] = 0.8644;
  temp_map["Be_std_LLNL1000"] = 0.9313;
  temp_map["Be_std_LLNL300"] = 0.8562;
  temp_map["Be_std_NIST_30000"] =  0.9313;
  temp_map["Be_std_NIST_30200"] = 0.9251;
  temp_map["Be_std_NIST_30300"] = 0.9221;
  temp_map["Be_std_NIST_30600"] = 0.9130;
  temp_map["Be_std_NIST_27900"] = 1;
  temp_map["Be_std_S555"] = 0.9124;
  temp_map["Be_std_S2007"] = 0.9124;
  temp_map["Be_std_BEST433"] = 0.9124;
  temp_map["Be_std_BEST433N"] = 1;
  temp_map["Be_std_S555N"] = 1;
  temp_map["Be_std_S2007N"] = 1;;

  // Same for Al-26. A zero placeholder is
  // also allowed. 
  temp_map["Al_std_KNSTD"] = 1.0;
  temp_map["Al_std_ZAL94"] = 0.9134;
  temp_map["Al_std_SMAL11"] = 1.021;
  temp_map["Al_std_0"] = 1.0;
  temp_map["Al_std_ZAL94N"] = 1.0;
  temp_map["Al_std_ASTER"] = 1.021;
  temp_map["Al_std_Z92-0222"] = 1;

  // Reference production rates at SLHL for spallation according to various
  // scaling schemes. Letter codes De, Du, Li, and St refer to
  // scaling schemes of Desilets et al. (2006), Dunai (2001), Lifton (2006),
  // and Stone (2000) respectively. The final letter code Lm refers to the
  // paleomagnetically corrected version of the Lal 1991/Stone 2000 scaling
  // factors.

  // These reference production rates are derived using the current version of
  // get_al_be_age.m. 

  // Al-26 and Be-10 production rates are linked by the production ratio,
  // not determined independently. This is because there are not very
  // many geological-calibration sites for Al-26. See the calibration data
  // sets and the accompanying paper for more information. 

  // This reflects the v 2.1 air pressure upgrade.
  // Also the update to 07KNSTD in version 2.2. 

  // Be-10 production rates. Brent's values. Greg agrees. 
  temp_map["P10_ref_St"] = 4.49;
  temp_map["delP10_ref_St"] = 0.39;

  temp_map["P10_ref_Lm"] = 4.39; 
  temp_map["delP10_ref_Lm"] = 0.37;

  temp_map["P10_ref_De"] = 4.41;
  temp_map["delP10_ref_De"] = 0.52;

  temp_map["P10_ref_Du"] = 4.43;
  temp_map["delP10_ref_Du"] = 0.52;

  temp_map["P10_ref_Li"] = 4.87;
  temp_map["delP10_ref_Li"] = 0.48;

  // Al-26 production rates are derived from Be-10 production rates 
  double R2610 = 6.1*1.106; // Update assumed production ratio
  temp_map["P26_ref_St"] = temp_map["P10_ref_St"]*R2610;
  temp_map[".delP26_ref_St"] = temp_map["delP10_ref_St"]*R2610;
  temp_map["P26_ref_Lm"] = temp_map["P10_ref_Lm"]*R2610; 
  temp_map["delP26_ref_Lm"] = temp_map["delP10_ref_Lm"]*R2610;
  temp_map["P26_ref_De"] = temp_map["P10_ref_De"]*R2610;
  temp_map["delP26_ref_De"] = temp_map["delP10_ref_De"]*R2610;
  temp_map["P26_ref_Du"] = temp_map["P10_ref_Du"]*R2610;
  temp_map["delP26_ref_Du"] = temp_map["delP10_ref_Du"]*R2610;
  temp_map["P26_ref_Li"] = temp_map["P10_ref_Li"]*R2610;
  temp_map["delP26_ref_Li"] = temp_map["delP10_ref_Li"]*R2610;

  // Muon interaction cross-sections. All follow Heisinger (2002a,b).
  // Note that the energy-dependence-of-muon-interaction-cross-section
  // exponent alpha is treated as model-dependent -- it's internal to 
  // P_mu_total.m and can't be passed.  

  temp_map["Natoms10"] = 2.006e22;
  temp_map["Natoms26"] = 1.003e22;

  // Be-10 interaction cross-sections
  // Restandardized by ETH-07KNSTD factor for version 2.2.1. 
  temp_map["k_neg10"] = (0.704 * 0.1828 * 0.0043)/1.096;
  temp_map["delk_neg10"] = (0.704 * 0.1828 * 0.0003)/1.096;
  temp_map["sigma190_10"] = (0.094e-27)/1.096;
  temp_map["delsigma190_10"] = (0.013e-27)/1.096;

  // Al-26 interaction cross-sections
  temp_map["k_neg26"] = 0.296 * 0.6559 * 0.022;
  temp_map["delk_neg26"] = 0.296 * 0.6559 * 0.002;
  temp_map["sigma190_26"] = 1.41e-27;
  temp_map["delsigma190_26"] = 0.17e-27;

  // Paleomagnetic records for use in time-dependent production rate schemes
  // Derived from Nat Lifton's compilation of paleomagnetic data from
  // various sources. See Lifton et al. (2006) and Pigati and Lifton (2005).

  // Load the magnetic field data
  // load PMag_Mar07
  // Relative dipole moment and time vector
  // al_be_consts.M = MM0; 
  // al_be_consts.t_M = t_M; 
  // These start at 7500 yr -- time slices are 7500,8500,9500,10500,11500
  // in order to use data from Yang et al; subsequent time slices are 
  // 12000:1000:800000 for SINT800 data; final two time points are 801000
  // and Inf. 

  // Cutoff rigidity blocks for past 6900 yr. 
  // TTRc and IHRC are lon x lat x time blocks of Rc values for the past 
  // 6900 years.
  // Both are derived by Nat Lifton from the magnetic field reconstructions of
  // Korte and Constable. 
  // TTRC has cutoff rigidity obtained by trajectory tracing -- these are for
  // the Lifton and Desilets scaling factors. IHRc has cutoff rigidity
  // obtained by finding magnetic inclination and horizontal field strength
  // from the field model, then applying Equation 2 of Dunai(2001). 
  // al_be_consts.TTRc = TTRc; % data block
  // al_be_consts.IHRc = IHRc; % data block
  // al_be_consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
  // al_be_consts.lon_Rc = lon_Rc;
  // al_be_consts.t_Rc = t_Rc; % time vector for Rc data block

  // Effective pole positions and field strengths inferred from K and C field
  // reconstructions for last 7000 yr. These are used in the
  // paleomagnetically-corrected implementation of the Lal SF. They are for
  // the same times as the RC slices in the data block above. Again,
  // generated by Nat Lifton -- hence KCL = Korte-Constable-Lifton. 

  // al_be_consts.MM0_KCL = MM0_KCL;
  // al_be_consts.lat_pp_KCL = lat_pp_KCL;
  // al_be_consts.lon_pp_KCL = lon_pp_KCL;

  // Solar variability from Lifton et al. 2005
  // Has been averaged and resampled to the same time slices as everything
  // else. 

  // al_be_consts.S = S; 
  // al_be_consts.SInf = 0.95; % Long-term mean S value;

  // set the data member map
  CRONUS_data_map = temp_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  
// This function sets CRONUS muon production paramteters
//
// Calculates the production rate of Al-26 or Be-10 by muons
// as a function of depth below the surface z (g/cm2) and
// site atmospheric pressure h (hPa).
//
// out.phi_vert_slhl muons/cm2/s/sr
// out.R_vert_slhl muons/g/s/sr
// out.R_vert_site muons/g/s/sr
// out.phi_vert_site muons/cm2/s/sr
// out.phi muons/cm2/yr
// out.R muons/g/yr
// out.P_fast atoms/g/yr
// out.P_neg atoms/g/yr
// out.Beta nondimensional
// out.Ebar GeV
// out.H g/cm2
// out.LZ g/cm2
//
// This uses the scheme in Heisinger and others (2002, 2 papers). The
// vertically traveling muon flux is scaled to the site elevation using
// energy-dependent attenuation lengths from Boezio et al. (2000). See the 
// hard-copy documentation for detailed citations and a full discussion of
// the calculation. 
//
// Note that some constants are internal to the function. The only ones that
// get passed from upstream are the ones that a) are nuclide-specific, or b) 
// actually have quoted uncertainties in Heisinger's papers. 
// The fraction of muons that are negative is internal; so is the
// energy-dependence exponent alpha.
//
// Original Written by Greg Balco -- UW Cosmogenic Nuclide Lab
// balcs@u.washington.edu
// March, 2006
// Part of the CRONUS-Earth online calculators: 
//      http://hess.ess.washington.edu/math
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::P_mu_total(double z,double h)
{
  map<string,double> temp_data_map;

  // calculator the atmospheric depth in g/cm2
  double H = (1013.25 - h)*1.019716;

  // find the vertical flux at SLHL
  double a = 258.5*(pow(100,2.66));
  double b = 75*(pow(100,1.66));
  double phi_vert_slhl = (a/((z+21000.0)*((pow((z+1000),1.66)) + b)))
                            *exp(-5.5e-6* z);

  // The above expression is only good to 2e5 g/cm2. We don't ever consider production
  // below that depth. The full-depth scheme appears in the comments below.
  // ------ begin full-depth flux equations -------
  //phiz_1 = (a./((z+21000).*(((z+1000).^1.66) + b))).*exp(-5.5e-6 .* z);
  //phiz_2 = 1.82e-6.*((121100./z).^2).*exp(-z./121100) + 2.84e-13;
  //out(find(z<200000)) = phiz_1(find(z<200000));
  //out(find(z>=200000)) = phiz_2(find(z>=200000));
  // ------ end full-depth flux equations -------

  // find the stopping rate of vertical muons at SLHL
  // this is done in a subfunction Rv0, because it gets integrated later.
  double R_vert_slhl = Rv0(z);

  // find the stopping rate of vertical muons at site
  double R_vert_site = R_vert_slhl*exp(H/LZ(z));

  // find the flux of vertical muons at site
  // integrate
  // ends at 200,001 g/cm2 to avoid being asked for an zero
  // range of integration -- 
  // get integration tolerance -- want relative tolerance around
  // 1 part in 10^4. 
  double tol = phi_vert_slhl*1e-4;
  double temp = tol;
  
  cout << "YO YOU HAVE NOT FINISHED YOU NEED TO NUMERICALLY INTEGRATE ON LINE 538 OF LSDCRNPARAMETERS" << endl;
  //  [temp,fcnt] = quad(@(x) Rv0(x).*exp(H./LZ(x)),z(a),(2e5+1),tol);
  //  % second variable assignment here to preserve fcnt if needed
  double phi_vert_site = temp;
   
  // invariant flux at 2e5 g/cm2 depth - constant of integration
  // calculated using commented-out formula above
  double phi_200k = (a/((2.0e5+21000.0)*((pow((2.0e5+1000.0),1.66)) + b)))
                      *exp(-5.5e-6 * 2.0e5);
  phi_vert_site = phi_vert_site + phi_200k;

  // find the total flux of muons at site
  // angular distribution exponent
  double nofz = 3.21 - 0.297*log((z+H)/100.0 + 42.0) + 1.21e-5*(z+H);
  // derivative of same
  double dndz = (-0.297/100.0)/((z+H)/100.0 + 42.0) + 1.21e-5;

  // caluculate phi in muons/cm2/s
  double phi_temp = (phi_vert_site*2* M_PI) / (nofz+1.0);

  // convert to muons/cm2/yr
  double phi = phi_temp*60.0*60.0*24.0*365.0;

  // find the total stopping rate of muons at site in muons/g/s
  double R_temp = (2*M_PI/(nofz+1.0))*R_vert_site 
                  - phi_vert_site*(-2*M_PI*(1/((nofz+1.0)*(nofz+1.0))))*dndz;
    
  // convert to negative muons/g/yr
  double R = R_temp*0.44*60.0*60.0*24.0*365.0;

  // Now calculate the production rates. 
  // Depth-dependent parts of the fast muon reaction cross-section
  double Beta = 0.846 - 0.015 * log((z/100.0)+1.0) 
                      + 0.003139 * (log((z/100.0)+1.0)*log((z/100.0)+1.0));
  double Ebar = 7.6 + 321.7*(1 - exp(-8.059e-6*z)) 
                    + 50.7*(1-exp(-5.05e-7*z));

  // internally defined constants
  double aalpha = 0.75;
  
  // this needs some logic for the isotope type
  double sigma0_Be10 = CRONUS_data_map["sigma190_10"]/(pow(190.0,aalpha));
  double sigma0_Al26 = CRONUS_data_map["sigma190_26"]/(pow(190.0,aalpha));
  
  // fast muon production
  double P_fast_Be10 = phi*Beta*(pow(Ebar,aalpha))
                          *sigma0_Be10*CRONUS_data_map["Natoms10"];
  double P_fast_Al26 = phi*Beta*(pow(Ebar,aalpha))
                          *sigma0_Al26*CRONUS_data_map["Natoms26"];
  
  // negative muon capture
  double P_neg_Be10 = R*CRONUS_data_map["k_neg10"];
  double P_neg_Al26 = R*CRONUS_data_map["k_neg26"];

  temp_data_map["phi_vert_slhl"] = phi_vert_slhl;
  temp_data_map["R_vert_slhl"] = R_vert_slhl;
  temp_data_map["phi_vert_site"] = phi_vert_site;
  temp_data_map["R_vert_site"] = R_vert_site;
  temp_data_map["phi"] = phi;
  temp_data_map["R"] = R;
  temp_data_map["Beta"] = Beta;
  temp_data_map["Ebar"] = Ebar;
  temp_data_map["P_fast_10Be"] = P_fast_Be10;
  temp_data_map["P_fast_26Al"] = P_fast_Al26;
  temp_data_map["P_neg_10Be"] = P_neg_Be10;
  temp_data_map["P_neg_26Al"] = P_neg_Al26;
  temp_data_map["H"] = H;
  temp_data_map["LZ"] = LZ(z);

  CRONUS_muon_data = temp_data_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this subfunction returns the stopping rate of vertically traveling muons
// as a function of depth z at sea level and high latitude.
// Modified from Greg Balco's CRONUS calculator
// z is the depth below the surface in g/cm^2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::Rv0(double z)
{
  double a = exp(-5.5e-6*z);
  double b = z + 21000.0;
  double c = pow((z + 1000.0),1.66) + 1.567e5;
  double dadz = -5.5e-6 * exp(-5.5e-6*z);
  double dbdz = 1.0;
  double dcdz = 1.66*(pow((z + 1000),0.66));
  
  double out = -5.401e7*(b*c*dadz-a*(c*dbdz+b*dcdz))/(b*b*c*c);
  return out;

  // full depth calculation appears in comments below
  // testing indicates this isn't really necessary
  //R_1 = -5.401e7 .* (b.*c.*dadz - a.*(c.*dbdz + b.*dcdz))./(b.^2 .* c.^2);
  //f = (121100./z).^2;
  //g = exp(-z./121100);
  //dfdz = (-2.*(121100.^2))./(z.^3);
  //dgdz = -exp(-z./121100)./121100;
  //R_2 = -1.82e-6.*(g.*dfdz + f.*dgdz);
  //out(find(z<200000)) = R_1(find(z<200000));
  //out(find(z>=200000)) = R_2(find(z>=200000));
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// this subfunction returns the effective atmospheric attenuation length for
// muons of range Z
// z is the depth in g/cm^2
//
// Original by Greg Balco as part of the CRONUS calculator
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::LZ(double z)
{

  //define range/momentum relation
  // table for muons in standard rock in Groom and others 2001
  // units are range in g cm-2 (column 2)
  //momentum in MeV/c (column 1)
  vector<double> data_for_LZ_range;
  vector<double> data_for_LZ_momentum;
  data_for_LZ_momentum.push_back(4.704e1);
  data_for_LZ_range.push_back(8.516e-1);
  data_for_LZ_momentum.push_back(5.616e1); 
  data_for_LZ_range.push_back(1.542e0);
  data_for_LZ_momentum.push_back(6.802e1); 
  data_for_LZ_range.push_back(2.866e0);
  data_for_LZ_momentum.push_back(8.509e1);
  data_for_LZ_range.push_back(5.698e0);
  data_for_LZ_momentum.push_back(1.003e2);
  data_for_LZ_range.push_back(9.145e0);
  data_for_LZ_momentum.push_back(1.527e2);
  data_for_LZ_range.push_back(2.676e1);
  data_for_LZ_momentum.push_back(1.764e2); 
  data_for_LZ_range.push_back(3.696e1);
  data_for_LZ_momentum.push_back(2.218e2);
  data_for_LZ_range.push_back(5.879e1);
  data_for_LZ_momentum.push_back(2.868e2);
  data_for_LZ_range.push_back(9.332e1);
  data_for_LZ_momentum.push_back(3.917e2);
  data_for_LZ_range.push_back(1.524e2);
  data_for_LZ_momentum.push_back(0.945e2);
  data_for_LZ_range.push_back(2.115e2);
  data_for_LZ_momentum.push_back(8.995e2);
  data_for_LZ_range.push_back(4.418e2);
  data_for_LZ_momentum.push_back(1.101e3);
  data_for_LZ_range.push_back(5.534e2);
  data_for_LZ_momentum.push_back(1.502e3);
  data_for_LZ_range.push_back(7.712e2);
  data_for_LZ_momentum.push_back(2.103e3);
  data_for_LZ_range.push_back(1.088e3);
  data_for_LZ_momentum.push_back(3.104e3);
  data_for_LZ_range.push_back(1.599e3);
  data_for_LZ_momentum.push_back(4.104e3);
  data_for_LZ_range.push_back(2.095e3);
  data_for_LZ_momentum.push_back(8.105e3);
  data_for_LZ_range.push_back(3.998e3);
  data_for_LZ_momentum.push_back(1.011e4);
  data_for_LZ_range.push_back(4.920e3);
  data_for_LZ_momentum.push_back(1.411e4);
  data_for_LZ_range.push_back(6.724e3);
  data_for_LZ_momentum.push_back(2.011e4);
  data_for_LZ_range.push_back(9.360e3);
  data_for_LZ_momentum.push_back(3.011e4);
  data_for_LZ_range.push_back(1.362e4);
  data_for_LZ_momentum.push_back(4.011e4);
  data_for_LZ_range.push_back(1.776e4);
  data_for_LZ_momentum.push_back(8.011e4);
  data_for_LZ_range.push_back(3.343e4);
  data_for_LZ_momentum.push_back(1.001e5);
  data_for_LZ_range.push_back(4.084e4);
  data_for_LZ_momentum.push_back(1.401e5);
  data_for_LZ_range.push_back(5.495e4);
  data_for_LZ_momentum.push_back(2.001e5);
  data_for_LZ_range.push_back(7.459e4);
  data_for_LZ_momentum.push_back(3.001e5); 
  data_for_LZ_range.push_back(1.040e5);
  data_for_LZ_momentum.push_back(4.001e5);
  data_for_LZ_range.push_back(1.302e5);
  data_for_LZ_momentum.push_back(8.001e5);
  data_for_LZ_range.push_back(2.129e5);

  // deal with zero situation
  if(z < 1)
  {
    z = 1.0;
  }

  // obtain momenta
  // use log-linear interpolation
  int n_momentum_dpoints = int(data_for_LZ_momentum.size());
  vector<double> log_momentum;
  vector<double> log_range;
  for(int i = 0; i<n_momentum_dpoints; i++)
  {
    log_momentum.push_back(log(data_for_LZ_momentum[i]));
    log_range.push_back(log(data_for_LZ_range[i]));
  }
  double P_MeVc = exp(interp1D_ordered(log_range,log_momentum,z));

  // obtain attenuation lengths
  double out = 263.0 + 150*(P_MeVc/1000.0);
  return out;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets the parameters to those used by Granger (need reference year!)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_Granger_parameters()
{
  S_t = 1;

  // from Vermeesh 2007
  lambda_10Be = 456e-9;		// in yr-1
  lambda_26Al = 980e-9;		// in yr-1
  lambda_14C = 121e-6;		// in yr-1
  lambda_36Cl = 230e-8;		// in yr-1

  // from Vermeesh 2007
  P0_10Be = 5.11;					// in a/g/yr
  P0_26Al = 30.31;				// in a/g/yr
  P0_14C = 5.86;					// in a/g/yr
  P0_36Cl = 55.45;				// in a/g/yr
  P0_21Ne = 20.29;				// in a/g/yr
  P0_3He = 97.40;					// in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 0.9724;
  F_10Be[1] = 0.0186;
  F_10Be[2] = 0.004;
  F_10Be[3] = 0.005;

  // dimensionless
  F_26Al[0] = 0.9655;
  F_26Al[1] = 0.0233;
  F_26Al[2] = 0.005;
  F_26Al[3] = 0.0062;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.0691;
  F_14C[2] = 0.0809;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.903;
  F_36Cl[1] = 0.0447;
  F_36Cl[2] = 0.05023;
  F_36Cl[3] = 0.0;
}

// function to set CRN parameters
// based on the Vermeesch approximation of the Schaller et al (2000)
// formulation
void LSDCRNParameters::set_Schaller_parameters()
{
  S_t =1;

  // from Vermeesh 2007
  lambda_10Be = 456e-9;		// in yr-1
  lambda_26Al = 980e-9;		// in yr-1
  lambda_14C = 121e-6;		// in yr-1
  lambda_36Cl = 230e-8;		// in yr-1

  // from Vermeesh 2007
  P0_10Be = 5.11;					// in a/g/yr
  P0_26Al = 30.31;				// in a/g/yr
  P0_14C = 5.86;					// in a/g/yr
  P0_36Cl = 55.45;				// in a/g/yr
  P0_21Ne = 20.29;				// in a/g/yr
  P0_3He = 97.40;					// in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 0.964;
  F_10Be[1] = 0.0266;
  F_10Be[2] = -0.0074;
  F_10Be[3] = 0.0168;

  // dimensionless
  F_26Al[0] = 0.9575;
  F_26Al[1] = 0.0315;
  F_26Al[2] = -0.009;
  F_26Al[3] = 0.02;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.1363;
  F_14C[2] = 0.0137;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.903;
  F_36Cl[1] = 0.0793;
  F_36Cl[2] = 0.0177;
  F_36Cl[3] = 0.0;
}

// this forces a neutron only calculation
void LSDCRNParameters::set_Neutron_only_parameters()
{
  S_t =1;

  // from Vermeesh 2007
  lambda_10Be = 456e-9;		// in yr-1
  lambda_26Al = 980e-9;		// in yr-1
  lambda_14C = 121e-6;		// in yr-1
  lambda_36Cl = 230e-8;		// in yr-1

  // from Vermeesh 2007
  P0_10Be = 5.11;					// in a/g/yr
  P0_26Al = 30.31;				// in a/g/yr
  P0_14C = 5.86;					// in a/g/yr
  P0_36Cl = 55.45;				// in a/g/yr
  P0_21Ne = 20.29;				// in a/g/yr
  P0_3He = 97.40;					// in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 1;
  F_10Be[1] = 0;
  F_10Be[2] = 0;
  F_10Be[3] = 0;

  // dimensionless
  F_26Al[0] = 1;
  F_26Al[1] = 0;
  F_26Al[2] = 0;
  F_26Al[3] = 0;

  // dimensionless
  F_14C[0] = 1;
  F_14C[1] = 0;
  F_14C[2] = 0;
  F_14C[3] = 0;

  // dimensionless
  F_36Cl[0] = 1;
  F_36Cl[1] = 0;
  F_36Cl[2] = 0;
  F_36Cl[3] = 0;
}


// this function takes a single scaling factor for
// elevation scaling, self shielding, snow shielding,
// and latitude scaling and produces scaling factors
// for each production mechamism.
// the scaling follows the approach of vermeesch 2008
// it uses a 'virstual' shielding depth to calcualte
// the updated scaling factors
void LSDCRNParameters::scale_F_values(double single_scaling)
{
  double tol = 1e-7;
  double x = 0;
  double new_x = 0;
  double test_scaling = 1e8;
  double dx =-10;

  // first do 10Be
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_10Be[0]+
                   exp(-new_x/Gamma[1])*F_10Be[1]+
                   exp(-new_x/Gamma[2])*F_10Be[2]+
                   exp(-new_x/Gamma[3])*F_10Be[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_10Be[0] = exp(-new_x/Gamma[0])*F_10Be[0];
  F_10Be[1] = exp(-new_x/Gamma[1])*F_10Be[1];
  F_10Be[2] = exp(-new_x/Gamma[2])*F_10Be[2];
  F_10Be[3] = exp(-new_x/Gamma[3])*F_10Be[3];

  //cout << "FINISHED 10Be x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_10Be[0] << endl << F_10Be[1] << endl << F_10Be[2] << endl << F_10Be[3] << endl;

  // now do 26Al
  x = 0;
  new_x = 0;
  test_scaling = 1e8;
  dx =-10;
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_26Al[0]+
                   exp(-new_x/Gamma[1])*F_26Al[1]+
                   exp(-new_x/Gamma[2])*F_26Al[2]+
                   exp(-new_x/Gamma[3])*F_26Al[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_26Al[0] = exp(-new_x/Gamma[0])*F_26Al[0];
  F_26Al[1] = exp(-new_x/Gamma[1])*F_26Al[1];
  F_26Al[2] = exp(-new_x/Gamma[2])*F_26Al[2];
  F_26Al[3] = exp(-new_x/Gamma[3])*F_26Al[3];

  //cout << "FINISHED 26Al x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_26Al[0] << endl << F_26Al[1] << endl << F_26Al[2] << endl << F_26Al[3] << endl;

  // now do 36Cl
  x = 0;
  new_x = 0;
  test_scaling = 1e8;
  dx =-10;
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_36Cl[0]+
                   exp(-new_x/Gamma[1])*F_36Cl[1]+
                   exp(-new_x/Gamma[2])*F_36Cl[2]+
                   exp(-new_x/Gamma[3])*F_36Cl[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_36Cl[0] = exp(-new_x/Gamma[0])*F_36Cl[0];
  F_36Cl[1] = exp(-new_x/Gamma[1])*F_36Cl[1];
  F_36Cl[2] = exp(-new_x/Gamma[2])*F_36Cl[2];
  F_36Cl[3] = exp(-new_x/Gamma[3])*F_36Cl[3];

  //cout << "FINISHED 36Cl x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_36Cl[0] << endl << F_36Cl[1] << endl << F_36Cl[2] << endl << F_36Cl[3] << endl;

  // now do 14C
  x = 0;
  new_x = 0;
  test_scaling = 1e8;
  dx =-10;
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_14C[0]+
                   exp(-new_x/Gamma[1])*F_14C[1]+
                   exp(-new_x/Gamma[2])*F_14C[2]+
                   exp(-new_x/Gamma[3])*F_14C[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_14C[0] = exp(-new_x/Gamma[0])*F_14C[0];
  F_14C[1] = exp(-new_x/Gamma[1])*F_14C[1];
  F_14C[2] = exp(-new_x/Gamma[2])*F_14C[2];
  F_14C[3] = exp(-new_x/Gamma[3])*F_14C[3];

  //cout << "FINISHED 14C x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_14C[0] << endl << F_14C[1] << endl << F_14C[2] << endl << F_14C[3] << endl;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function converts elevation to atmospheric pressure
// It follows the cronus calculator NCEPatm_2 function
// written by Greg Balco
// SMM
// 4/12/2014
//
// Looks up surface pressure and 1000 mb temp from NCEP reanalysis
// and calculates site atmospheric pressures using these as inputs to the
// standard atmosphere equation. 
//
// Syntax: pressure = NCEPatm_2(site_lat,site_lon,site_elv);
// 
// Requires:
///       site_lat: latitude (DD). Southern hemisphere is negative.
//       site_lon: longitude (DD). Western hemisphere is negative.
//           Tries to deal with 0-360 longitudes gracefully.
//       site_elv: elevation (m).
//
// Returns site pressure in hPa.
//
// Vectorized. Send vectors of equal length.
//
// Note: this must load the data file NCEP2.mat whenever called. 
// Repeated calls to this function will be slow for this reason. 
//
// Also: This function is OK but not great for Antarctica.
// Use antatm.m instead. 
//
// Remember: it is always better to estimate the average pressure at your 
// site using a pressure-altitude relation obtained from nearby station
// data.
//
// Original m code Written by Greg Balco -- UW Cosmogenic Nuclide Lab
// balcs@u.washington.edu
// October, 2007
// Part of the CRONUS-Earth online calculators: 
//      http://hess.ess.washington.edu/math
//
// Copyright 2001-2007, University of Washington
// All rights reserved
// Developed in part with funding from the National Science Foundation.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::NCEPatm_2(double site_lat, double site_lon, double site_elev)
{
  // deal with negative loingitudes
  if(site_lon < 0)
  {
    site_lon = site_lon+360.0;
  }
  
  // check to see if data is loaded:
  if (int(gm_hgt.size()) != 8)
  {
    string path_to_data;
    cout << "You didn't load the NCEP data. Doing that now. " << endl;
    cout << "Enter path to data files: " << endl;
    cin >> path_to_data;
    load_parameters_for_atmospheric_scaling(path_to_data);
  }
  
  // now, interpolate sea level pressure and temperature
  double site_slp = interp2D_bilinear(NCEPlat, NCEPlon, meanslp, 
                        site_lat, site_lon);
  double site_T = interp2D_bilinear(NCEPlat, NCEPlon, meant1000, 
                        site_lat, site_lon);
                        
  
  
  double site_T_degK = site_T + 273.15;

  // Some More parameters
  double gmr = -0.03417; // Assorted constants (this has come from Greg Balco's code)
  double dtdz = 0.0065;  // Lapse rate from standard atmosphere 
  
  // Calculate site pressure using the site-specific SLP and T1000 with the
  // standard atmosphere equation.

  //cout << site_T_degK << endl;
  //cout << "Log term: " << log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) << endl;
  //cout << "Exp term: " << exp( (gmr/dtdz)*( log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) ) ) << endl;

  double out = site_slp*exp( (gmr/dtdz)*( log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) ) );
  
  cout << endl;
  //cout << "Site sea level pressure: " << site_slp << " and site Temp: "<< site_T 
  //     << " and pressure: " << out << endl << endl <<endl;
  return out;


}



#endif
