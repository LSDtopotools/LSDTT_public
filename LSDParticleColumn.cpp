//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDParticleColumn.cpp
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for holding particles in a column
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
using namespace std;

#ifndef LSDParticleColumn_CPP
#define LSDParticleColumn_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  create, the default constructor
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::create()
{
   Row = 0;
   Col = 0;
   NodeIndex = 0;
   
   RockDensity = 2000;
   SoilDensity = 1300;
   
   UseDenstyProfile = false;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates a column with particles at constant depth spacing and
// wit SS cosmo concentration
// SMM 25 July 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::initiate_SS_cosmo_column_3CRN(int start_type, 
          double startxLoc, double startyLoc,
		      double start_depth, double particle_spacing, 
		      double zeta, double eff_eros_rate,
		      LSDCRNParameters& CRN_param)
{
  // get number of particles
  int N_particles = double(start_depth/particle_spacing)+1;

  // create the list
  list<LSDCRNParticle> CRN_plist;

  double this_depth;

  // now loop over the depths, inserting particles and setting them to steady state
  double top_depth = start_depth-double(N_particles-1)*particle_spacing;
  cout << "Inserting particles, top depth is: " << top_depth << endl;
  for (int p = 0; p< N_particles; p++)
  {
    // first the depth
    this_depth = top_depth+double(p)*particle_spacing;

    // now insert the particle into the list
    insert_particle_into_column(start_type, startxLoc, startyLoc,
			      this_depth,zeta);
  }

  // now loop through the particles, setting to steady
  list<LSDCRNParticle>::iterator part_iter;
  part_iter = CRN_plist.begin();
  while(part_iter != CRN_plist.end())
  {
    // update the CRN_concntrations
    ( *part_iter ).update_10Be_SSfull(eff_eros_rate,CRN_param);
    ( *part_iter ).update_14C_SSfull(eff_eros_rate,CRN_param);
    ( *part_iter ).update_21Ne_SSfull(eff_eros_rate,CRN_param);
 
    part_iter++;   
  }

  CRNParticleList = CRN_plist;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
	      

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// update_list
// this is the function that inserts a particle into a list
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::insert_particle_into_column(int start_type, 
   double startxLoc, double startyLoc,
	 double start_depth, double zeta)
{
	double d = start_depth;
	double eff_d;                   // effective depth in g/cm^2
	double z_p = zeta-start_depth;  // starting elevation of the particle
	
	// first check to see if we use depth column
	if (UseDenstyProfile)
	{
    cout << "LSDParticleColumn::insert_part_into_column, you asked to use density\n";
    cout << "profile but this is not implemented yet. Fatal error";
    exit(EXIT_FAILURE);
  }
  else
  {
    if (d > SoilThickness)
    {
       // get effective depth if the particle is below the soil layer
       eff_d = SoilThickness*SoilDensity*0.1+ (d-SoilThickness)*RockDensity*0.1;
    }
    else
    {
      // gets the effective depth if the particle is in the soil
      eff_d =  SoilThickness*SoilDensity*0.1;
    }
  }
		
	LSDCRNParticle CRN_tp(start_type,startxLoc,startyLoc,d,eff_d,z_p);
	CRNParticleList.push_back(CRN_tp);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#endif
