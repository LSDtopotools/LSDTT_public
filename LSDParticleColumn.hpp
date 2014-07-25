//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDParticleColumn
// Land Surface Dynamics Particle Column
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for tracing particles and retaining geochemical information
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
#include <list>
#include "LSDCRNParameters.hpp"
#include "LSDParticle.hpp"
using namespace std;

#ifndef LSDParticleColumn_H
#define LSDParticleColumn_H


/// @brief This is a class for a particle that can be tracked through simulations
/// and retains data about position and chemical content
class LSDParticleColumn
{
	public:
    /// @brief The default constructor. It is the only possible constructor
	  LSDParticleColumn()			{ create(); }
	


	protected:
  
    /// The row of the particle column. 
    /// This allows the obect to be linked to an LSDRaster	
    int Row;
    
    /// The row of the particle column. 
    /// This allows the obect to be linked to an LSDRaster	    
    int Col;
    
    /// The row of the particle column. 
    /// This allows the obect to be linked to an LSDFlowInfo object	    
    int NodeIndex;
    
    /// the soil density in kg/m^3
    double SoilDensity;
    
    /// the rock density in kh/m^3
    double RockDensity;
    
    /// the soil thickness. It is a float since LSDRasters are floats
    float SoilThickness;
    
    /// the cellsize of the raster if the column is linked to a raster.
    /// It is a float since LSDRasters are floats
    /// used to track if the particle is withing the cell
    float DataResolution;
    
    /// This tells the model if it needs to use the density profile
    bool UseDenstyProfile;
    
    /// This is the density profile. The depths are the points where
    /// density has been sampled
    vector<double> DensityDepths;
    
    /// This is the density profile. The densities (in kg/m^3) are the density 
    /// measuremetns at depths set by DenstiyDepths    
    vector<double> DensityDensities;
    
    /// This is the list for holding the particles
    list<LSDCRNParticle> CRNParticleList;
    
    /// @brief This function initiates a list with particles evenly spaced
    /// @param start_type the starting type of the particle
    /// @param startxLoc the starting x location of the particle
    /// @param startxLoc the starting y location of the particle
    /// @param start_depth the starting depth of the particle in metres
    /// @param particle_spacing the spacing in metres between particles (in depth)
    /// @param zeta the surface elevation (m)
    /// @param eff_eros_rate the erosion rate in g/cm^2/yr
    /// @param CRN_param a CRNParameters object
    /// @author SMM
    /// @date 25/07/2014
    void initiate_SS_cosmo_column_3CRN(int start_type, 
          double startxLoc, double startyLoc,
		      double start_depth, double particle_spacing, 
		      double zeta, double eff_eros_rate,
		      LSDCRNParameters& CRN_param);

    /// @brief This inserts particles into the column
    /// @param start_type the starting type of the particle
    /// @param startxLoc the starting x location of the particle
    /// @param startxLoc the starting y location of the particle
    /// @param start_depth the starting depth of the particle in metres
    /// @param zeta the elevation of the surface
    void insert_particle_into_column(int start_type, 
       double startxLoc, double startyLoc,
	     double start_depth, double zeta);
    
    
  private:
    void create();		
};

#endif