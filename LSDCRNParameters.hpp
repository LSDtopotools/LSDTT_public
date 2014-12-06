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
#include <iostream>
#include <vector>
#include <map>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDCRNParameters_H
#define LSDCRNParameters_H

/// @brief This class contains parameters used in cosmogenic nuclide calculations
/// It sits seperately from the particle object since it applies to an
/// entire environment and not just an individual particle. 
/// Seperating the object in this way reduces memory redundancy 
class LSDCRNParameters
{
  public:
  
  /// @brief The default constructor. It is the only possible constructor
  LSDCRNParameters()			{ create(); }
  
  /// This is a friend class so that it can be called from the particle 
  friend class LSDCRNParticle;

  /// @brief function for loading parameters that allow pressure calculation
  /// from elevation
  /// @author SMM
  /// @date 02/12/2014
  void load_parameters_for_atmospheric_scaling(string path_to_params);
  
  /// @brief This function sets a numer of parameters that are used
  /// to replicate the CRONUS calculator.
  /// @details the original parameters are derived from the 
  /// make_al_be_consts_v22
  /// Written by Greg Balco -- Berkeley Geochronology Center
  ///  balcs@u.washington.edu -- balcs@bgc.org
  ///  February, 2008
  ///  Part of the CRONUS-Earth online calculators: 
  ///     http://hess.ess.washington.edu/math
  /// @author SMM
  /// @date 06/12/2014
  void set_CRONUS_data_maps();

  /// @brief alculates the production rate of Al-26 or Be-10 by muons
  /// @detail This uses the scheme in Heisinger and others (2002, 2 papers). The
  ///  vertically traveling muon flux is scaled to the site elevation using
  ///  energy-dependent attenuation lengths from Boezio et al. (2000). See the 
  ///  hard-copy documentation for detailed citations and a full discussion of
  ///  the calculation. 
  ///  Note that some constants are internal to the function. The only ones that
  ///  get passed from upstream are the ones that a) are nuclide-specific, or b) 
  ///  actually have quoted uncertainties in Heisinger's papers. 
  ///  The fraction of muons that are negative is internal; so is the
  ///  energy-dependence exponent alpha.
  ///  Original Written by Greg Balco -- UW Cosmogenic Nuclide Lab
  ///  balcs@u.washington.edu
  ///  March, 2006
  ///  Part of the CRONUS-Earth online calculators: 
  ///      http://hess.ess.washington.edu/math
  /// @param z depth below the surface z (g/cm2)
  /// @param h atmospheric pressure (hPa)
  /// @author SMM
  /// @date 06/12/2014
  void P_mu_total(double z,double h);
 
  /// @brief this subfunction returns the stopping rate of vertically traveling muons
  /// as a function of depth z at sea level and high latitude.
  /// @detail Modified from Greg Balco's CRONUS calculator
  /// @param z is the depth below the surface in g/cm^2
  /// @return Rv0 the muon stopping rate
  /// @author SMM
  /// @date 06/12/2014
  double Rv0(double z); 

  /// @brief this subfunction returns the effective atmospheric attenuation length for
  /// muons of range Z
  /// @detail Original by Greg Balco as part of the CRONUS calculator
  /// @param z is the depth in g/cm^2
  /// @return effective atmospheric attenuation length in g/cm^2
  /// @author SMM
  /// @date 06/12/2014
  double LZ(double z);
  
  // functions for altering the parameter values
  
  /// @brief This resets the F, Gamma and P0 values so that they conform to 
  /// Granger scaling. Adopted from from Vermeesh 2007
  /// @author SMM
  /// @date 01/01/2010
  void set_Granger_parameters();
  
  /// @brief This resets the F, Gamma and P0 values so that they conform to 
  /// Schaller (2009) scaling. Adopted from from Vermeesh 2007
  /// @author SMM
  /// @date 01/01/2010	
  void set_Schaller_parameters();

  /// @brief This sets the F values to use neutron only production
  /// @details F0 == 1, all other F values == 0
  /// @author SMM
  /// @date 14/07/2014	
  void set_Neutron_only_parameters();
  
  /// @brief this function takes a single scaling factor for
  /// elevation scaling, self shielding, snow shielding,
  /// and latitude scaling and produces scaling factors
  /// for each production mechamism.
  /// the scaling follows the approach of vermeesch 2008
  /// it uses a 'virtual' shielding depth to calculate
  /// the updated scaling factors
  /// @param single_scaling a lumped scaling factor	
  /// @author SMM
  /// @date 01/01/2010
  void scale_F_values(double single_scaling);	// parameter values
  
  /// @brief this changes the 10Be decay. It is here because
  /// 10Be decay rates reported in the literature have changed
  /// @author SMM
  /// @date 01/01/2010
  void update_10Be_decay(double new_decay)	{ lambda_10Be = new_decay; }
  
  /// @brief this changes the 10Be P0 value. It is here because
  /// 10Be decay rates reported in the literature have changed
  /// @author SMM
  /// @date 01/01/2010	
  void update_10Be_P0(double new_P0)			{ P0_10Be = new_P0; }

  /// @brief This calcualtes the atmospheric pressure given latidude, longitude
  /// and elevation
  /// @details Looks up surface pressure and 1000 mb temp from NCEP reanalysis
  /// and calculates site atmospheric pressures using these as inputs to the
  /// standard atmosphere equation. 
  /// Also: This function is OK but not great for Antarctica.
  /// Use antatm.m instead. 
  /// Remember: it is always better to estimate the average pressure at your 
  /// site using a pressure-altitude relation obtained from nearby station
  /// data.
  ///
  /// Original m code Written by Greg Balco -- UW Cosmogenic Nuclide Lab
  /// @param site_lat latitude (DD). Southern hemisphere is negative.
  /// @param site_lon longitude (DD). Western hemisphere is negative.
  ///       Tries to deal with 0-360 longitudes gracefully.
  /// @param site_elv elevation (m).
  /// @return site pressure in hPa.
  /// @author SMM
  /// @date 04/12/2014
  double NCEPatm_2(double site_lat, double site_lon, double site_elev);
  
  
  
  
  private:	
  /// @brief This is called by the default constructor. 
  /// It is the only possible constructor
  void create();

  /// The attenudation lengths in g/cm^2. Each element in the array
  /// refers to a different production mechanism
  double Gamma[4];			// attenuation legths in g/cm^2

  /// F values allocating production to spallation and muon production for 10Be
  double F_10Be[4];
  
  /// F values allocating production to spallation and muon production for 14C
  double F_14C[4];
  
  /// F values allocating production to spallation and muon production for 26Al
  double F_26Al[4];
  
  /// F values allocating production to spallation and muon production for 36Cl
  double F_36Cl[4];
  
  /// F values allocating production to spallation and muon production for 21Ne
  double F_21Ne[4];
  
  /// F values allocating production to spallation and muon production for 3He
  double F_3He[4];

  /// topographic shielding
  /// other shielding and scaling calucalted using the scale_F_values function
  double S_t;

  /// decay rate for 10Be in yr-1
  double lambda_10Be;	
  
  /// decay rate for 26Al in yr-1	
  double lambda_26Al;	
  
  /// decay rate for 14C in yr-1	
  double lambda_14C;	
  
  /// decay rate for 36Cl in yr-1	
  double lambda_36Cl;		

  /// production rate for 10Be in a/g/yr
  double P0_10Be;	
  
  /// production rate for 26Al in a/g/yr		
  double P0_26Al;		
  
  /// production rate for 14C in a/g/yr
  double P0_14C;			
  
  /// production rate for 36Cl in a/g/yr
  double P0_36Cl;		
  
  /// production rate for 21Ne in a/g/yr
  double P0_21Ne;			
  
  /// production rate for 3He in a/g/yr
  double P0_3He;			
  
  /// This is a data map used for storing CRONUS calcluator parameters
  map<string,double> CRONUS_data_map;
  
  /// This is a data map used for storing CRONUS muon parameters
  map<string,double> CRONUS_muon_data;
  
  /// levels: the levels for the atmospheric scaling of pressure
  vector<double> levels;
  
  /// This is an index for the latitudes for atmospheric scaling
  vector<double> NCEPlat;
  
  /// This is an index for the longitudes for atmospheric scaling
  vector<double> NCEPlon;
  
  /// This is an array holding sea level perssures
  Array2D<double> meanslp;
  
  /// This is an array healing mean temperatures
  Array2D<double> meant1000;
  
  /// This is a vector of arrays holding something called gp_hgt;
  vector< Array2D<double> > gp_hgt;
  
  /// This is a vector of arrays holding something called gp_hgt;
  vector< Array2D<double> > gm_hgt;  
  
  
};

#endif
