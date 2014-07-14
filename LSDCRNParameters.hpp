// LSDCRNParameters.hpp
// header file for the discrete particle object
// this version of the dicrete particle keeps track
// of two integers, the time and the position

#include <iostream>
#include <vector>
using namespace std;

#ifndef LSDCRNParameters_H
#define LSDCRNParameters_H

class LSDCRNParameters
{
	public:
	LSDCRNParameters()			{ create(); }
	friend class LSDCRNParticle;

	// functions for altering the parameter values
	void set_Granger_parameters();
	void set_Schaller_parameters();
	void scale_F_values(double single_scaling);	// parameter values
	void update_10Be_decay(double new_decay)	{ lambda_10Be = new_decay; }
	void update_10Be_P0(double new_P0)			{ P0_10Be = new_P0; }

	private:
	void create();

	double Gamma[4];			// attenuation legths in g/cm^2

	// these F values allocate production to spallation and muon-derived production
	double F_10Be[4];
	double F_14C[4];
	double F_26Al[4];
	double F_36Cl[4];
	double F_21Ne[4];
	double F_3He[4];

	// topographic shielding
	// other shielding and scaling calucalted using the scale_F_values function
	double S_t;

	// decay rates
	double lambda_10Be;		// in yr-1
	double lambda_26Al;		// in yr-1
	double lambda_14C;		// in yr-1
	double lambda_36Cl;		// in yr-1

	// production rates
	double P0_10Be;			// in a/g/yr
	double P0_26Al;			// in a/g/yr
	double P0_14C;			// in a/g/yr
	double P0_36Cl;			// in a/g/yr
	double P0_21Ne;			// in a/g/yr
	double P0_3He;			// in a/g/yr
};

#endif
