// LSDCRNParameters
// This is an object that holds parameters for the cosmogenic 
// particles. It calculates and sets parameters within a function
// so that these do not have to be calculated for each particle. 
#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include "LSDStatsTools.hpp"
#include "LSDCRNParameters.hpp"
using namespace std;

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
