// LSDParticle.hpp
// implementation of the LSDParticle.cpp object
#include <iostream>
#include <vector>
#include "LSDCRNParameters.hpp"
using namespace std;

#ifndef LSDParticle_H
#define LSDParticle_H

// empty class definition so that we can use friend functions
class LSDCRNParameters;

class LSDParticle
{
 public:
    LSDParticle()					{ create(); }
    LSDParticle( int StartType)			{ create(StartType); }
    LSDParticle( int StartType, double StartAge)	{ create(StartType, StartAge); }
    LSDParticle( int StartType, int StartCI, double StartAge, double StartOSLage, double StartxLoc, double StartdLoc)
    						{ create(StartType, StartCI, StartAge, StartOSLage, StartxLoc, StartdLoc); }
    LSDParticle( int StartType, double StartxLoc, double StartdLoc)
    						{ create(StartType, StartxLoc, StartdLoc); }

    int    getType() const			{ return Type; }
    int    getCellIndex() const		{ return CellIndex; }
    double getAge() const			{ return Age; }
    double getOSLage() const		{ return OSLage; }
    double getxLoc() const			{ return xLoc; }
    double getdLoc() const			{ return dLoc; }

    LSDParticle(const LSDParticle& tP)
    	{ create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc()); }
    LSDParticle(LSDParticle& tP)
    	{ create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc()); }

    LSDParticle& operator=(const LSDParticle& tP);

    std::ostream& operator<< (std::ostream&);

    void incrementAge(double dt);
    void setCellIndex(int CI);
    void OSLexpose();
    void SoilAgeExpose();
	  void update_xLoc(double new_xLoc);
    void displaceReflect(double dx,double dd,double h,double dt);
    int  test_domain(double lambda);

 protected:
    int Type;
    int CellIndex;
    double Age;
    double OSLage;
    double xLoc;
    double dLoc;

 private:
 
    /// @brief creates a default particle of type 0, 
    /// cell index of -1, age of 0 and OSLage of -9999
    /// @author SMM
    /// @date 01/01/2008
    void create();
    
    /// @brief creates a default particle 
    /// cell index of -1, age of 0 and OSLage of -9999
    /// @param StartType the starting type of the particle
    /// @author SMM
    /// @date 01/01/2008
    void create(int);
    
    /// @brief creates a default particle of 
    /// cell index of -1, and OSLage of -9999
    /// @param StartType the starting type of the particle
    /// @param StartAge the starting age of the particle
    /// @author SMM
    /// @date 01/01/2008    
    void create(int, double);
    
    void create(int, int, double, double, double, double);
    void create(int, double, double);
};

// the CRN tracer particle object
// this object tracks the CRN concnentration in a number
// of particles
class LSDCRNParticle: public LSDParticle
{
	public:
	LSDCRNParticle()			{ create(); }
	LSDCRNParticle(int startType, double startxLoc,double startdLoc,
					double start_effdloc, double startzloc)
							{ create(startType,startxLoc,startdLoc,
							  start_effdloc,startzloc); }
	LSDCRNParticle(int startType, double startxLoc,double startzeta_Loc)
							{ create(startType,startxLoc,startzeta_Loc); }
	LSDCRNParticle(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C)
							{ create(startType,startxLoc,startzeta_Loc,
								startdLoc, start_effdLoc,
								start_C10Be, start_C14C); }
	LSDCRNParticle(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C, double start_21Ne)
							{ create(startType,startxLoc,startzeta_Loc,
								startdLoc, start_effdLoc,
								start_C10Be, start_C14C, start_21Ne); }
	LSDCRNParticle(int startType, int startCellIndex, double startAge, double startOSLAge,
					double startxLoc,double startdLoc, double startefdLoc,
					double startzLoc, double start_C10Be, double start_C26Al,
					double start_C36Cl, double start_C14C,
					double start_C21Ne, double start_C3He,
					double start_Cf7Be, double start_Cf10Be,
					double start_Cf210Pb, double start_Cf137Cs)
							{ create(startType, startCellIndex, startAge, startOSLAge,
								startxLoc, startdLoc, startefdLoc,
								startzLoc, start_C10Be, start_C26Al,
								start_C36Cl, start_C14C,
								start_C21Ne, start_C3He,
								start_Cf7Be, start_Cf10Be,
								start_Cf210Pb, start_Cf137Cs); }


    LSDCRNParticle(const LSDCRNParticle& tP)
    	{ create(tP.getType(), tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc(),
    	         tP.geteffective_dLoc(),tP.get_zetaLoc(),tP.getConc_10Be(),
    	         tP.getConc_26Al(), tP.getConc_36Cl(), tP.getConc_14C(),
    	         tP.getConc_21Ne(), tP.getConc_3He(),
    	         tP.getConc_f7Be(), tP.getConc_f10Be(),
    	         tP.getConc_f210Pb(), tP.getConc_f137Cs()); }
    LSDCRNParticle(LSDCRNParticle& tP)
    	{ create(tP.getType(), tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc(),
    	         tP.geteffective_dLoc(),tP.get_zetaLoc(),tP.getConc_10Be(),
    	         tP.getConc_26Al(), tP.getConc_36Cl(), tP.getConc_14C(),
    	         tP.getConc_21Ne(), tP.getConc_3He(),
    	         tP.getConc_f7Be(), tP.getConc_f10Be(),
    	         tP.getConc_f210Pb(), tP.getConc_f137Cs());   }

	LSDCRNParticle& operator=(const LSDCRNParticle& tP);
	LSDCRNParticle& operator=(LSDCRNParticle& tP);



	// accessing function to get at the data elements
	double getConc_10Be() const			{ return Conc_10Be; }
	double getConc_26Al() const			{ return Conc_26Al; }
	double getConc_36Cl() const			{ return Conc_36Cl; }
	double getConc_14C() const			{ return Conc_14C; }
	double getConc_21Ne() const			{ return Conc_21Ne; }
	double getConc_3He() const			{ return Conc_3He; }
	double getConc_f7Be() const			{ return Conc_f7Be; }
	double getConc_f10Be() const		{ return Conc_f10Be; }
	double getConc_f210Pb() const		{ return Conc_f210Pb; }
	double getConc_f137Cs() const		{ return Conc_f137Cs; }
	double geteffective_dLoc() const	{ return effective_dLoc; }
	double get_zetaLoc() const			{ return zetaLoc; }

	// update nuclide concentrations in the event of constant erosion
	// erosion is in g/cm^2/yr
	void update_10Be_conc(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_26Al_conc(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_14C_conc(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_36Cl_conc(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_21Ne_conc(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_3He_conc(double dt,double erosion, LSDCRNParameters& CRNp);

	// update nuclide concentrations if erosion is increasing (or decreasing) linearly
	// alpha is the change in erosion rate--do not set alpha to zero
	void update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, LSDCRNParameters& CRNp);

	// update nuclide concentrations using only spallation
	void update_10Be_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_26Al_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_14C_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);
	void update_36Cl_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);

	// update nuclide concentrations with constant values
	void update_cosmo_conc_const(double C_10Be, double C_26Al, double C_36Cl,
								 double C_14C, double C_21Ne, double C_3He);

	// update the concentration of nuclides if there is a steady state
	// profile
	void update_10Be_SSfull(double erosion, LSDCRNParameters& CRNp);
	void update_26Al_SSfull(double erosion, LSDCRNParameters& CRNp);
	void update_14C_SSfull(double erosion, LSDCRNParameters& CRNp);
	void update_36Cl_SSfull(double erosion, LSDCRNParameters& CRNp);
	void update_21Ne_SSfull(double erosion, LSDCRNParameters& CRNp);
	void update_3He_SSfull(double erosion, LSDCRNParameters& CRNp);

	// functions collecting the updating functions: these update
	// all nuclides at once
	void update_all_CRN(double dt, double erosion, LSDCRNParameters& CRNp);
	void update_all_CRN_neutron_only(double dt, double erosion, LSDCRNParameters& CRNp);
	void update_all_CRN_SSfull(double erosion_rate, LSDCRNParameters& CRNp);

	// caluclate the apparent erosion from nuclide concentrations
	double apparent_erosion_10Be_neutron_only(double rho, LSDCRNParameters& CRNp);
	double apparent_erosion_26Al_neutron_only(double rho, LSDCRNParameters& CRNp);
	double apparent_erosion_14C_neutron_only(double rho, LSDCRNParameters& CRNp);
	double apparent_erosion_36Cl_neutron_only(double rho, LSDCRNParameters& CRNp);
	double apparent_erosion_21Ne(double rho, LSDCRNParameters& CRNp);
	double apparent_erosion_3He(double rho, LSDCRNParameters& CRNp);

	// functions for dealing with fallout numclides
	void update_fallout10Be_simple_density(double dt, double M_supply_surface,
					double rho_s, double k_f10Be, double deltad, LSDCRNParameters& CRNp);
	void update_fallout10Be_simple_density_2exp(double dt, double M_supply_surface,
					double rho_skg, double k1_f10Be, double k2_f10Be, double chi_f10Be,
					double deltad_m, LSDCRNParameters& CRNp);

	// functions for managing the shielding depth, depth, and elevation of
	// particles
	void update_depths(double delta_d, double delta_ed);
	void update_zetaLoc(double new_zeta);
	void update_zetaLoc_with_new_surface(double new_zeta);
	void erode_mass_only(double dt, double mass_erosion_rate);
	void erode_mass_only_linear_increase(double dt, double mass_erosion_rate, double alpha);

  protected:
  double effective_dLoc;		// the effective depth: in g/cm^2
	double zetaLoc;				// the elevation relative to arbitrary
								// datum (m)
	double Conc_10Be;			// concnetration of 10Be in atoms/g
	double Conc_26Al;			// a/g
	double Conc_36Cl;			// a/g
	double Conc_14C;			// a/g
	double Conc_21Ne;			// a/g
	double Conc_3He;			// a/g
	double Conc_f7Be;			// fallout units tba
	double Conc_f10Be; 			// fallout, units a/g
	double Conc_f210Pb;			// fallout, units tba
	double Conc_f137Cs;			// fallout, units tba



	private:
	// functions for creating particles
	void create();
	void create(int startType, double startxLoc,
	            double startzLoc);
	void create(int startType, double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzloc);
	void create(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C);
	void create(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C, double start_21Ne);
	void create(int startType, int startCellIndex, double startAge, double startOSLAge,
	            double startxLoc,double startdLoc, double startefdLoc,
				double startzLoc, double start_C10Be, double start_C26Al,
				double start_C36Cl, double start_C14C,
				double start_C21Ne, double start_C3He,
				double start_Cf7Be, double start_Cf10Be,
				double start_Cf210Pb, double start_Cf137Cs);


};


#endif
