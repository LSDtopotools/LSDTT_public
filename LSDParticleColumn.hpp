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


/// This is a class for a particle that can be tracked through simulations
/// and retains data about position and chemical content
class LSDParticleColumn
{
	public:
    /// @brief The default constructor. It is the only possible constructor
	  LSDParticleColumn()			{ create(); }
	


	protected:
  
    /// The row of the particlecolumn. 
    /// This allows the obect to be linked to an LSDRaster	
    int Row;
    
    /// The row of the particlecolumn. 
    /// This allows the obect to be linked to an LSDRaster	    
    int Col;
    
    /// The row of the particlecolumn. 
    /// This allows the obect to be linked to an LSDFlowInfo object	    
    int NodeIndex;
    
    /// the soil density in kg/m^3
    double SoilDensity;
    
    /// the rock density in kh/m^3
    double RockDensity;
    
    /// the soil thickness. It is a float since LSDRasters are floats
    float SoilThickness;
    
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
    
    /// 
    
  private:
    void create();		
};

#endif