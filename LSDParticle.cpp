//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDParticle
// Land Surface Dynamics Particle
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
#include <fstream>
#include <math.h>
#include <iostream>
#include "LSDStatsTools.hpp"
#include "LSDParticle.hpp"
using namespace std;

const double one_min_exp_neg_2 = 1-exp(-2);
const double one_min_exp_neg_5 = 1-exp(-5);

// default constructor
void LSDParticle::create()
 {
  Type = 0;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
 }

// start with a specific type
void LSDParticle::create( int StartType )
 {
  Type = StartType;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
 }


void LSDParticle::create( int StartType, double StartAge )
 {
  Type = StartType;
  CellIndex = -1;
  Age = StartAge;
 }

void LSDParticle::create( int StartType, int StartCI,  double StartAge, double StartOSLage, double StartxLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = StartCI;
  Age = StartAge;
  OSLage = StartOSLage;
  xLoc = StartxLoc;
  dLoc = StartdLoc;
 }

void LSDParticle::create(int StartType, double StartxLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
  xLoc = StartxLoc;
  dLoc = StartdLoc;
 }

LSDParticle& LSDParticle::operator=(const LSDParticle& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.getType(),rhs.getCellIndex(),rhs.getAge(),rhs.getOSLage(), rhs.getxLoc(),rhs.getdLoc());
   }
  return *this;
 }

std::ostream& operator<<(std::ostream& os, const LSDParticle& tP)
 {
  os <<   tP.getType() << " " << tP.getCellIndex() << " " << tP.getAge() << " "
       << tP.getOSLage() << " " << tP.getxLoc() << " " << tP.getdLoc();
  return os;
 }

void LSDParticle::incrementAge(double dt)
 {
  if (Age < 0)
   Age = dt;
  else
   Age += dt;
  if (OSLage > 0)
   OSLage += dt;
 }

void LSDParticle::setCellIndex(int tempCI)
{
	CellIndex = tempCI;
}

void LSDParticle::OSLexpose()
 {
	 OSLage = 0;
 }

 void LSDParticle::SoilAgeExpose()
  {
 	 Age = 0;
 }

void LSDParticle::update_xLoc(double new_xLoc)
{
	xLoc = new_xLoc;
}

void LSDParticle::displaceReflect(double dx,double dd,double h, double dt)
 {
  xLoc += dx;

  //std::cout << "tPart.cpp LINE 77 dx  is: " << dd << std::endl;
  double dOld = dLoc;
  double dNew = dLoc+dd;
  if (dNew > h)
   {
    //std::cout << "tPart.cpp LINE 84 dNew  is: " << dNew << " and dd is: " << dd
    //          << " and dOld is: "<< dOld << std::endl;
    dLoc = 2*h - dd - dOld;
    //std::cout << "tPart.cpp LINE 86 dLoc  is: " << dLoc << std::endl;
   }
  else if (dNew <= 0)
   {
    dLoc = -(dd+dOld);
    OSLage = (dd+dOld)*dt/dd;
   }
  else
   dLoc = dNew;

  //if (dLoc > h)
  // std::cout << "tPart.cpp LINE 86 dLoc  is: " << dLoc << std::endl;
 }

int  LSDParticle::test_domain(double lambda)
 {
  int td;
  if (xLoc >= 0 && xLoc <= lambda)
   td = 1;
  else
   {
    td = -1;
    OSLage = -9999;
    Age = -9999;
    xLoc = -9999;
    dLoc = -9999;
   }
  return td;
 }


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDCRNParticle object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// functions for the CRN loaded tracer particle
void LSDCRNParticle::create()
 {
	double rho_r = 2650;		// in kg/m^3

	Type = 0;
	CellIndex = -1;
	Age = 0;
	OSLage = 0;
	xLoc = 0;			// in meters
	dLoc = 5;			// in meters
	zetaLoc = 100;
    Conc_10Be = 0.0;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = rho_r*dLoc*0.1;	// the 0.1
	 									// converts between kg/m^2
	 									// to g/cm^2
 }

void LSDCRNParticle::create(int startType, double startxLoc,
	            double startzLoc)
{
	double rho_r = 2650;		// in kg/m^3

	Type = startType;
	CellIndex = -1;
	Age = 0;
	OSLage = 0;
	xLoc = startxLoc;			// in meters
	dLoc = 0;			// in meters
	zetaLoc = startzLoc;
    Conc_10Be = 0.0;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = rho_r*dLoc*0.1;	// the 0.1
	 									// converts between kg/m^2
	 									// to g/cm^2
 }

void LSDCRNParticle::create(int startType, double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc)
{
	Type = startType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	zetaLoc = startzLoc;
    Conc_10Be = 0.0;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = start_effdloc;
}

void LSDCRNParticle::create(int startType, double startxLoc,double startzeta_Loc,
 					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C)
{
	Type = startType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = start_effdLoc;
	zetaLoc = startzeta_Loc;
    Conc_10Be = start_C10Be;
	Conc_26Al = 0;
	Conc_36Cl = 0;
	Conc_14C = start_C14C;
	Conc_21Ne = 0;
	Conc_3He = 0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
}

void LSDCRNParticle::create(int startType, double startxLoc,double startzeta_Loc,
 					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C, double start_21Ne)
{
	Type = startType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = start_effdLoc;
	zetaLoc = startzeta_Loc;
    Conc_10Be = start_C10Be;
	Conc_26Al = 0;
	Conc_36Cl = 0;
	Conc_14C = start_C14C;
	Conc_21Ne = start_21Ne;
	Conc_3He = 0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
}

void LSDCRNParticle::create(int startType, int startCellIndex, double startAge, double startOSLAge,
	            double startxLoc,double startdLoc, double startefdLoc,
				double startzLoc, double start_C10Be, double start_C26Al,
				double start_C36Cl, double start_C14C,
				double start_C21Ne, double start_C3He,
				double start_Cf7Be, double start_Cf10Be,
				double start_Cf210Pb, double start_Cf137Cs)
{
	Type = startType;
	CellIndex = startCellIndex;
	Age = startAge;
	OSLage = startOSLAge;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = startefdLoc;
	zetaLoc = startzLoc;
    Conc_10Be = start_C10Be;
	Conc_26Al = start_C26Al;
	Conc_36Cl = start_C36Cl;
	Conc_14C = start_C14C;
	Conc_21Ne = start_C21Ne;
	Conc_3He = start_C3He;
	Conc_f7Be = start_Cf7Be;
	Conc_f10Be = start_Cf10Be;
	Conc_f210Pb = start_Cf210Pb;
	Conc_f137Cs = start_Cf137Cs;
}


LSDCRNParticle& LSDCRNParticle::operator=(const LSDCRNParticle& rhs)
{
	if (&rhs != this)
    {
		create(rhs.getType(), rhs.getCellIndex(), rhs.getAge(),rhs.getOSLage(),rhs.getxLoc(),rhs.getdLoc(),
    	         rhs.geteffective_dLoc(),rhs.get_zetaLoc(),rhs.getConc_10Be(),
    	         rhs.getConc_26Al(), rhs.getConc_36Cl(), rhs.getConc_14C(),
    	         rhs.getConc_21Ne(), rhs.getConc_3He(),
    	         rhs.getConc_f7Be(), rhs.getConc_f10Be(),
    	         rhs.getConc_f210Pb(), rhs.getConc_f137Cs());
    }
    return *this;
}

LSDCRNParticle& LSDCRNParticle::operator=(LSDCRNParticle& rhs)
{
	if (&rhs != this)
    {
		create(rhs.getType(), rhs.getCellIndex(), rhs.getAge(),rhs.getOSLage(),rhs.getxLoc(),rhs.getdLoc(),
     	         rhs.geteffective_dLoc(),rhs.get_zetaLoc(),rhs.getConc_10Be(),
     	         rhs.getConc_26Al(), rhs.getConc_36Cl(), rhs.getConc_14C(),
    	         rhs.getConc_21Ne(), rhs.getConc_3He(),
    	         rhs.getConc_f7Be(), rhs.getConc_f10Be(),
    	         rhs.getConc_f210Pb(), rhs.getConc_f137Cs());
    }
    return *this;
}

// this function just sets the concentration of cosmogenic in situ nuclides
void LSDCRNParticle::update_cosmo_conc_const(double C_10Be, double C_26Al, double C_36Cl,
								 double C_14C, double C_21Ne, double C_3He)
{
	Conc_10Be = C_10Be;
	Conc_26Al = C_26Al;
	Conc_36Cl = C_36Cl;
	Conc_14C  = C_14C;
	Conc_21Ne = C_21Ne;
	Conc_3He  = C_3He;
}


// this function updates the concentration of CRNs in a particle
// the model assumes that during the timestep the change in the
// 'depth' of the particle occurs ofver a constant rate.
// The depth in this case is an equvalent depth...it is linearly
// proportional to depth if overlying density is constant, but
// really represetnt the mass per area above a point in the subsurface
// thus the erosion rate represent a mass removal rate.
void LSDCRNParticle::update_10Be_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	// fist berillium
	double Be_exp = exp(-dt*CRNp.lambda_10Be);
	//double Be_depth_exp = exp(-d_0/Gamma_Be1);

    //cout << "LINE 236 LSDParticle.cpp " << endl;
    //cout << lambda_10Be << " " << Be_exp << endl;
    //cout << "starting conc: " << Conc_10Be << endl;
    //cout << "starting depth: " << effective_dLoc << endl;
    //cout << "erosion_rate: "<< erosion_rate << endl;
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_10Be[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_10Be))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
		//cout << "F_10Be["<<i<<"]: " << CRNp.F_10Be[i] << " and sum_term: " << sum_term <<endl;
	}

    //cout << "and sum term is: " << sum_term << endl;

	Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*CRNp.P0_10Be*sum_term;
	//cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

// this function updates the 10Be concentration if there is a linear
// increase (or decrease) in erosion rate.
void LSDCRNParticle::update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, LSDCRNParameters& CRNp)
{
	// fist berillium
	double Be_exp = exp(-dt*CRNp.lambda_10Be);
	//double Be_depth_exp = exp(-d_0/Gamma_Be1);

    //cout << "LINE 236 LSDParticle.cpp " << endl;
    //cout << lambda_10Be << " " << Be_exp << endl;
    //cout << "starting conc: " << Conc_10Be << endl;
    //cout << "starting depth: " << effective_dLoc << endl;
    //cout << "erosion_rate: "<< erosion_rate << endl;
	double sum_term = 0;
	double L_j, M_j, A_j, B_j, D_j;
	double erfi_arg1, erfi_arg2;
	for (int j = 0; j<4; j++)
	{
		A_j = sqrt(0.5*CRNp.Gamma[j]/alpha);
		B_j = sqrt(2*alpha*CRNp.Gamma[j]);
		D_j = sqrt(0.5*alpha/CRNp.Gamma[j]);
		L_j = exp(-( (erosion_rate+CRNp.Gamma[j]*CRNp.lambda_10Be)*(erosion_rate+CRNp.Gamma[j]*CRNp.lambda_10Be)+
					2*alpha*(effective_dLoc+dt*CRNp.Gamma[j]*CRNp.lambda_10Be) )/(B_j*B_j) );
		erfi_arg1 = erosion_rate/B_j + A_j*CRNp.lambda_10Be;
		erfi_arg2 = D_j*dt + erfi_arg1;
		M_j = erfi(erfi_arg2) - erfi(erfi_arg1);

		//std::cout << "j is: " << j << endl;
		//std::cout << "erfi_arg1 is: " << erfi_arg1 << endl;
		//std::cout << "erfi_arg2 is: " << erfi_arg2 << endl;
		//std::cout << "M_j term is: " << M_j << endl;
		//std::cout << "L_j term is: " << L_j << endl;


		sum_term+= (CRNp.F_10Be[j]*A_j*L_j*M_j);
	}

    //cout << "and sum term is: " << sum_term << endl;

	Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*CRNp.P0_10Be*sum_term*sqrt(M_PI);
	//cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

void LSDCRNParticle::update_10Be_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{

	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double Be_exp = exp(-dt*CRNp.lambda_10Be);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_10Be))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_10Be);

	Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*Be_exp*CRNp.P0_10Be*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This updates the 10Be concetration if erosion is steady using the full
// 4 attenuation depth model of Vermeesch (2007)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParticle::update_10Be_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_10Be[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
	}

    //cout << "and sum term is: " << sum_term << endl;

	Conc_10Be = CRNp.S_t*CRNp.P0_10Be*sum_term;
	//cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_10Be_neutron_only(double rho, LSDCRNParameters& CRNp)
{
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eff_eros;                    // in g/cm2/yr
	double app_eros;                        // in m/yr
	
	double exp_term = exp(-effective_dLoc/Gamma_neutron);

	app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_10Be
                                /Conc_10Be-CRNp.lambda_10Be);
  app_eros = app_eff_eros*10/rho;                              
	return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


void LSDCRNParticle::update_26Al_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Al_exp = exp(-dt*CRNp.lambda_26Al);

	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_26Al[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_26Al))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
	}

	Conc_26Al = Conc_26Al*Al_exp +  CRNp.S_t*Al_exp*CRNp.P0_26Al*sum_term;
}



void LSDCRNParticle::update_26Al_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2
	double Al_exp = exp(-dt*CRNp.lambda_26Al);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_26Al))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_26Al);

	Conc_26Al = Conc_26Al*Al_exp +  CRNp.S_t*Al_exp*CRNp.P0_26Al*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_26Al_neutron_only(double rho, LSDCRNParameters& CRNp)
{
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eff_eros;                    // in g/cm2/yr
	double app_eros;                        // in m/yr
	
	double exp_term = exp(-effective_dLoc/Gamma_neutron);

	app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_26Al
                                /Conc_26Al-CRNp.lambda_26Al);
  app_eros = app_eff_eros*10/rho;                              
	return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDCRNParticle::update_26Al_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_26Al[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
	}

	Conc_26Al = CRNp.S_t*CRNp.P0_26Al*sum_term;
}

void LSDCRNParticle::update_14C_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double C_exp = exp(-dt*CRNp.lambda_14C);

	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_14C[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_14C))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_14C);
	}

	Conc_14C = Conc_14C*C_exp +  CRNp.S_t*C_exp*CRNp.P0_14C*sum_term;
}

void LSDCRNParticle::update_14C_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

	double C_exp = exp(-dt*CRNp.lambda_14C);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_14C))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_14C);

	Conc_14C = Conc_14C*C_exp +  CRNp.S_t*C_exp*CRNp.P0_14C*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_14C_neutron_only(double rho, LSDCRNParameters& CRNp)
{
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eff_eros;                    // in g/cm2/yr
	double app_eros;                        // in m/yr
	
	double exp_term = exp(-effective_dLoc/Gamma_neutron);

	app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_14C
                                /Conc_14C-CRNp.lambda_14C);
  app_eros = app_eff_eros*10/rho;                              
	return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDCRNParticle::update_14C_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_14C[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_14C);
	}

	Conc_14C = CRNp.S_t*CRNp.P0_14C*sum_term;
}

void LSDCRNParticle::update_36Cl_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Cl_exp = exp(-dt*CRNp.lambda_36Cl);

	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_36Cl[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_36Cl))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_36Cl);
	}

	Conc_36Cl = Conc_36Cl*Cl_exp +  CRNp.S_t*Cl_exp*CRNp.P0_36Cl*sum_term;
}

void LSDCRNParticle::update_36Cl_conc_neutron_only(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

	double Cl_exp = exp(-dt*CRNp.lambda_36Cl);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_36Cl))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_36Cl);


	Conc_36Cl = Conc_36Cl*Cl_exp +  CRNp.S_t*Cl_exp*CRNp.P0_36Cl*sum_term;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_36Cl_neutron_only(double rho, LSDCRNParameters& CRNp)
{
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eff_eros;                    // in g/cm2/yr
	double app_eros;                        // in m/yr
	
	double exp_term = exp(-effective_dLoc/Gamma_neutron);

	app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_36Cl
                                /Conc_36Cl-CRNp.lambda_36Cl);
  app_eros = app_eff_eros*10/rho;                              
	return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




void LSDCRNParticle::update_36Cl_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_36Cl[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_36Cl);
	}

	Conc_36Cl = CRNp.S_t*CRNp.P0_36Cl*sum_term;
}

void LSDCRNParticle::update_21Ne_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron= CRNp.Gamma[0];					// in g/cm^2

	Conc_21Ne = Conc_21Ne +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_21Ne*
	             (exp(dt*erosion_rate/Gamma_neutron) - 1)/erosion_rate;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_21Ne(double rho, LSDCRNParameters& CRNp)
{
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eff_eros;                    // in g/cm2/yr
	double app_eros;                        // in m/yr
	
	double exp_term = exp(-effective_dLoc/Gamma_neutron);

	app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_21Ne
                                /Conc_21Ne);
  app_eros = app_eff_eros*10/rho;                              
	return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


void LSDCRNParticle::update_21Ne_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

	//cout << endl << endl <<"BEFORE Conc_21Ne: " << Conc_21Ne << endl;
	if (erosion_rate*erosion_rate < 0.0000000001)
	{
		Conc_21Ne = 0;
	}
	else
	{
		Conc_21Ne = CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_21Ne/erosion_rate;
	}
	//cout << "AFTER Conc_21Ne: " << Conc_21Ne << endl;
}

void LSDCRNParticle::update_3He_conc(double dt,double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2
	Conc_3He = Conc_3He +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_3He*
	             (exp(dt*erosion_rate/Gamma_neutron) - 1)/erosion_rate;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the apparent erosion rate. 
// if rho is in kg/m^3 this returns erosion in m/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParticle::apparent_erosion_3He(double rho, LSDCRNParameters& CRNp)
{
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eff_eros;                    // in g/cm2/yr
	double app_eros;                        // in m/yr
	
	double exp_term = exp(-effective_dLoc/Gamma_neutron);

	app_eff_eros = Gamma_neutron*(exp_term*CRNp.S_t*CRNp.P0_3He
                                /Conc_3He);
  app_eros = app_eff_eros*10/rho;                              
	return app_eros;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



void LSDCRNParticle::update_3He_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2


	if (erosion_rate*erosion_rate < 0.0000000001)
	{
		Conc_3He = 0;
	}
	else
	{
		Conc_3He = Conc_3He +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_3He/erosion_rate;
	}
}

void LSDCRNParticle::update_all_CRN(double dt, double erosion_rate, LSDCRNParameters& CRNp)
{
	//cout << "LINE 445 LSDCRNParticle.cpp updating CRN" << endl;
	update_10Be_conc(dt,erosion_rate, CRNp);
	update_26Al_conc(dt,erosion_rate, CRNp);
	update_36Cl_conc(dt,erosion_rate, CRNp);
	update_14C_conc(dt,erosion_rate, CRNp);
	update_21Ne_conc(dt,erosion_rate, CRNp);
	update_3He_conc(dt,erosion_rate, CRNp);
}

void LSDCRNParticle::update_all_CRN_SSfull(double erosion_rate, LSDCRNParameters& CRNp)
{
	//cout << "LINE 445 LSDCRNParticle.cpp updating CRN" << endl;
	update_10Be_SSfull(erosion_rate, CRNp);
	update_26Al_SSfull(erosion_rate, CRNp);
	update_36Cl_SSfull(erosion_rate, CRNp);
	update_14C_SSfull(erosion_rate, CRNp);
	update_21Ne_SSfull(erosion_rate, CRNp);
	update_3He_SSfull(erosion_rate, CRNp);
}

void LSDCRNParticle::update_all_CRN_neutron_only(double dt, double erosion_rate, LSDCRNParameters& CRNp)
{
	update_10Be_conc_neutron_only(dt,erosion_rate, CRNp);
	update_26Al_conc_neutron_only(dt,erosion_rate, CRNp);
	update_36Cl_conc_neutron_only(dt,erosion_rate, CRNp);
	update_14C_conc_neutron_only(dt,erosion_rate, CRNp);
	update_21Ne_conc(dt,erosion_rate, CRNp);
	update_3He_conc(dt,erosion_rate, CRNp);
}

// this updates fallout radionuclides
// NOTE!!!!
// the untis of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// we convert rho into g/cm^3
// adn deltad into cm
// withing the LSDCRNParticle object
// the units of k_f10Be here are cm^2/g
void LSDCRNParticle::update_fallout10Be_simple_density(double dt, double M_supply_surface,
					double rho_skg, double k_f10Be, double deltad_m, LSDCRNParameters& CRNp)
{
	// first find which depth interval the particle is in
	int depth_interval = int(dLoc/deltad_m);
	double d_top = double(depth_interval)*deltad_m*100;
	double d_bottom = double(depth_interval+1)*deltad_m*100;
	double deltad = deltad_m*100;
	// the factor of 100 is to convert to cm

	// convert density to g/cm^3
	double rho_s = rho_skg/1000;

	// get the cutoff depth
	double cutoff_depth = 5/(rho_s*k_f10Be);

	if (dLoc*100 > cutoff_depth)
	{
		Conc_f10Be += - Conc_f10Be*CRNp.lambda_10Be;
	}
	else
	{
		Conc_f10Be += dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top) -exp(-k_f10Be*rho_s*d_bottom) )/
		              (deltad*rho_s*one_min_exp_neg_5) - Conc_f10Be*CRNp.lambda_10Be;
		//		cout << " added conc: " <<  dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top)
		//						-exp(-k_f10Be*rho_s*d_bottom) )/
		//		              (deltad*rho_s*one_min_exp_neg_5) << endl;
	}
}

// this updates fallout radionuclides
// NOTE!!!!
// the units of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// we convert rho into g/cm^3
// adn deltad into cm
// within the LSDCRNParticle object
// the units of k_f10Be here are cm^2/g
// chi_f10Be is the fraction of shallow fallout 10Be
void LSDCRNParticle::update_fallout10Be_simple_density_2exp(double dt, double M_supply_surface,
					double rho_skg, double k1_f10Be, double k2_f10Be, double chi_f10Be, double deltad_m, LSDCRNParameters& CRNp)
{
	// first find which depth interval the particle is in
	int depth_interval = int(dLoc/deltad_m);
	double d_top = double(depth_interval)*deltad_m*100;
	double d_bottom = double(depth_interval+1)*deltad_m*100;
	double deltad = deltad_m*100;
	// the factor of 100 is to convert to cm

	// convert density to g/cm^3
	double rho_s = rho_skg/1000;

	// get the cutoff depth for k1
	double cutoff_depth1 = 5/(rho_s*k1_f10Be);
	double cutoff_depth2 = 5/(rho_s*k2_f10Be);

	if (dLoc*100 > cutoff_depth2)
	{
		Conc_f10Be += - Conc_f10Be*CRNp.lambda_10Be;
	}
	if (dLoc*100 < cutoff_depth2 && dLoc*100 > cutoff_depth1)
	{
		Conc_f10Be += dt*M_supply_surface*(1-chi_f10Be)*( exp(-k2_f10Be*rho_s*d_top) -exp(-k2_f10Be*rho_s*d_bottom) )/
		              (deltad*rho_s*one_min_exp_neg_5) - Conc_f10Be*CRNp.lambda_10Be;
	}
	else
	{
		Conc_f10Be += dt*M_supply_surface*(1-chi_f10Be)*( exp(-k2_f10Be*rho_s*d_top) -exp(-k2_f10Be*rho_s*d_bottom) )/
		               (deltad*rho_s*one_min_exp_neg_5) +
		              dt*M_supply_surface*chi_f10Be*( exp(-k1_f10Be*rho_s*d_top) -exp(-k1_f10Be*rho_s*d_bottom) )/
		               (deltad*rho_s*one_min_exp_neg_5) -
		              Conc_f10Be*CRNp.lambda_10Be;
		//		cout << " added conc: " <<  dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top)
		//						-exp(-k_f10Be*rho_s*d_bottom) )/
		//		              (deltad*rho_s*one_min_exp_neg_5) << endl;
	}
}

// this function resets the depth and effective depth
void LSDCRNParticle::update_depths(double d, double ed)
{
	dLoc=d;
	effective_dLoc=ed;
}

// this function changes the effective depth (i.e., the shielding depth)
// for a constant rate of erosion
void LSDCRNParticle::erode_mass_only(double dt, double mass_erosion_rate)
{
	effective_dLoc-= mass_erosion_rate*dt;
}

// this function changes the effective depth (i.e., the shielding depth)
// for erosion that is changing linearly in time
void LSDCRNParticle::erode_mass_only_linear_increase(double dt, double mass_erosion_rate, double alpha)
{
	effective_dLoc-= 0.5*dt*(alpha*dt+2*mass_erosion_rate);
}

// this function simply resets zeta
void LSDCRNParticle::update_zetaLoc(double new_zeta)
{
	zetaLoc = new_zeta;
}

// this function resets the zeta locations using an updated
// surface elevation that preserves the d and effective d locations
// it is only for use with a model that has no soil
void LSDCRNParticle::update_zetaLoc_with_new_surface(double new_zeta)
{
	zetaLoc = new_zeta-dLoc;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





