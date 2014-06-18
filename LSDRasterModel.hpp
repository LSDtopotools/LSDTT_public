///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDRasterModel.cpp
/// cpp file for the LSDRasterModel object
/// LSD stands for Land Surface Dynamics
/// This object provides an environment for landscape evolution modelling, which can then
/// be integrated with the topographic analysis tools to efficiently analyse model runs.
///                                                                                     
/// The landscape evolution model uses implicit methods to provide stability with
/// relatively long timesteps.  Fluvial erosion is solved following Braun and Willet (2013)
/// using the fastscape algorithm, whilst hillslope sediment transport is modelled as a
/// non-linear diffusive sediment flux, following the implicit scheme developed for
/// MuDDPile.
///
/// The aim is to have two complimentary models:
/// i) a simple coupled hillslope-channel model in which large scale landscape dynamics
/// can be modelled
/// ii) a more complex treatment of hillslopes explicitly incorporating the role of
/// vegetation in driving sediment production and transport, and that copes with the  
/// with the transition from soil mantled-bedrock hillslopes at high erosion rates.
///
/// In order to run the model, one needs a parameter file that should be read by the 
/// driver function.
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
/// @author David T. Milodowski, University of Edinburgh
/// @author Martin D. Hurst, British Geological Survey
/// @author Fiona Clubb, University of Edinburgh
/// @author Stuart Grieve, University of Edinburgh
/// @author James Jenkinson, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1		24/07/2013
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TNT/tnt.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "LSDRaster.hpp"
#include "LSDRasterSpectral.hpp"
#include "LSDJunctionNetwork.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterModel_H
#define LSDRasterModel_H

///@brief Create model objects to use LSDRaster methods on synthetic landscapes.
class LSDRasterModel: public LSDRasterSpectral
{
	public:
	/// @brief Create a deafult LSDRasterModel (100x100)
	/// @return An instance of LSDRasterModel
	LSDRasterModel() { create(); }

	/// @brief Create a LSDRasterModel from a parameter file
	/// @return An instance of LSDRasterModel
	/// @param master_param A filenam for the master parameter file
	LSDRasterModel( string master_param )  { create(master_param); }
	/// @brief Create an LSDRasterModel from a file.
	/// Uses a filename and file extension
	/// @return LSDRasterModel
	/// @param filename A String, the file to be loaded.
	/// @param extension A String, the file extension to be loaded.
	LSDRasterModel(string filename, string extension)	{ create(filename, extension); default_parameters();}
	/// @brief Create an LSDRasterModel from memory.
	/// @return LSDRasterModel
	/// @param nrows An integer of the number of rows.
	/// @param ncols An integer of the number of columns.
	/// @param xmin A float of the minimum X coordinate.
	/// @param ymin A float of the minimum Y coordinate.
	/// @param cellsize A float of the cellsize.
	/// @param ndv An integer of the no data value.
	/// @param data An Array2D of floats in the shape nrows*ncols,
	///containing the data to be written.
	LSDRasterModel(int nrows, int ncols, float xmin, float ymin,
	          float cellsize, float ndv, Array2D<float> data)
								{ default_parameters(); create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }
	/// @brief Create an LSDRasterModel from an LSDRaster.
	/// @return LSDRasterModel
	/// @param An_LSDRaster LSDRaster object.							
	LSDRasterModel(LSDRaster& An_LSDRaster) 				{ create(An_LSDRaster); default_parameters();}

	/// @brief Create a blank raster nodel
	/// @return LSDRasterModel
	/// @param NCols Height of raster
	/// @param NRows Width of raster
	LSDRasterModel(int NRows, int NCols);

	/// operator
	LSDRasterModel& operator=(const LSDRasterModel& LSDR);

	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// INITIALISATION MODULE
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// This module initialises the model runs, calling the required function from
	/// the initial topography and loads the parameters from the parameter file.
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	///  Initialise_model_with_parabolic_surface
	///----------------------------------------------------------------------------
	void initialize_model(
	  string& parameter_file, string& run_name, float& dt, float& EndTime, float& PrintInterval,
	  float& k_w, float& b, float& m, float& n, float& K, float& ErosionThreshold, 
	  float& K_nl, float& S_c, float& UpliftRate, float& PrecipitationRate,
	  float& NorthBoundaryElevation, float& SouthBoundaryElevation,
	  Array2D<float>& PrecipitationFlux, Array2D<float>& SlopesBetweenRows,
	  Array2D<float>& SlopesBetweenColumns, Array2D<float>& ErosionRate);

	void initialize_model( string parameter_file );

	/// --------------------------------------------------------------------------
	/// Check steady state
	/// --------------------------------------------------------------------------
	void check_steady_state( void );

	void check_recording( void );

	bool check_end_condition( void );

	void check_periodicity_switch( void );

	bool check_if_hung( void );
	     
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// BUFFER SURFACE
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// 2 functions to buffer the raster surface.
	///----------------------------------------------------------------------------
	/// This first function is used as a simple way to implement boundary conditions,
	/// particularly no flux and periodic boundary conditions.
	/// The buffered surface has NRows+2 rows and NCols+2 columns.
	/// The integer b_type sets the type of boundary conditions, but currently there
	/// is only one implementation: no flux across N and S; periodic for E and W. 
	///----------------------------------------------------------------------------
	LSDRasterModel create_buffered_surf(int b_type);
	///----------------------------------------------------------------------------
	/// This second version has periodic boundaries at E and W boundaries, and 
	/// Neumann boundary conditions (prescribed elevations) at the N and S
	/// boundaries.
	///------------------------------------------------------------------------------
	LSDRasterModel create_buffered_surf(float South_boundary_elevation,float North_boundary_elevation);
	
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// CALCULATE EROSION RATES
	/// Simple function that creates an array with the erosion rates for a given 
	/// timestep
	///----------------------------------------------------------------------------
	Array2D<float> calculate_erosion_rates( void );
	float get_erosion_at_cell(int row, int col);
	
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// UPLIFT SURFACE
	/// Apply uplift field to the raster.  Overloaded function so that the first
	/// simply considers uniform uplift, the second allows user to use a prescribed
	/// uplift fields of greater complexity, for example taking account of fault
	/// geometry.
	///----------------------------------------------------------------------------
	LSDRasterModel uplift_surface(float UpliftRate, float dt);
	///------------------------------------------------------------------------------
	/// Specified uplift field 
	/// uplift field should be specified as an array with the same dimensions as the
	/// elevation raster, permitting non-uniform uplift fields to be applied in the 
	/// model.
	///----------------------------------------------------------------------------
	LSDRasterModel uplift_surface(Array2D<float> UpliftRate, float dt);
	
	/// ---------------------------------------------------------------------------
	/// Intrinsice method of uplifting the Raster
	/// Uplift field attribute is incremented onto RasterData itself
	/// ---------------------------------------------------------------------------
	void uplift_surface( void );

	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	/// CREATE PRECIPITION FLUX ARRAY
	/// Produces precipitation array from provided precipitation rate.
	///---------------------------------------------------------------------------
	Array2D<float> precip_array_from_precip_rate(float precip_rate);
	
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// TOPOGRAPHIC DERIVATIVES
	/// Specifically, this function gets the topographic slopes, as required for
	/// the sediment flux calculations.  The slopes are stored as two matrices, one
	/// that stores slopes between rows, the other which for slopes between 
	/// columns.  Note that this is a finite volume model that utilises cubic model
	/// voxels. Sediment fluxes are only permitted through the faces.
	///
	/// For slopes between columns, the entry at S[row][col] refers to the slope
	/// between zeta at node [row][col] and at node [row][col+1].  Likewise for the
	/// slopes between rows.  In short, the center points of the slopes are offset
	/// by 1/2 a node spacing in the positive direction.
	/// Note that there are NCols +1 and NRows +1 columns and rows respectively
	///----------------------------------------------------------------------------
	void get_slopes(Array2D<float>& SlopesBetweenRows, Array2D<float>& SlopesBetweenCols);
	///----------------------------------------------------------------------------
	/// get_topographic_divergence
	/// gets the topographic divergence at each point in the model domain.  Use
	/// buffered topography
	///----------------------------------------------------------------------------
	Array2D<float> get_topographic_divergence();

	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// HYDROLOGICAL TOOLS
	///----------------------------------------------------------------------------
	/// calculate_channel_width_wolman
	/// This function calculates channel width using the wolman method.
	/// NOTE: typically Q_w will be in m^3/s.
	/// EXAMPLE: in Salmon River, Idaho (Emmett, 1975 cited in Knighton 1988):
	///          k_w = 2.77 and b = 0.56. b is often assumed to be 0.5
	///----------------------------------------------------------------------------
	float calculate_channel_width_wolman(float Q_w, float k_w, float b);
	///----------------------------------------------------------------------------
	/// array_channel_width_wolman
	/// this function calcualtes channel width in a stand alone module so the widths
	/// can be tested
	///----------------------------------------------------------------------------
	Array2D<float> array_channel_width_wolman(Array2D<float>& Q_w, float& k_w, float& b);
	
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// EROSION RATES/SEDIMENT FLUXES
	///
	///----------------------------------------------------------------------------
	/// this caluclates the fluvial erosion rate at each point
	///----------------------------------------------------------------------------
	Array2D<float> calculate_fluvial_erosion_rate(Array2D<float> ChannelWidth, Array2D<float> Q_w,
							Array2D<float> TopoDivergence, float K, float n, float m, float eros_thresh);                                  

	
	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// IMPLICIT MODEL COMPONENTS
	///------------------------------------------------------------------------------
	/// Implicit schemes for combination of hillslope sediment transport using
	/// non-linear hillslope transport law, and fluvial erosion.  This is essentially
	/// the implicit implementation of MuddPILE, but has been modified so that now
	/// fluvial erosion is undertaken using FASTSCAPE (Braun and Willet, 2013), which
	/// greatly increases computational efficiency.
	///------------------------------------------------------------------------------
	/// calculate_k_values_for_assembly_matrix/mtl_initiate_assembler_matrix
	/// this function creates vectors of integers that refer to the k values, that is
	/// the index into the vectorized matrix of zeta values, that is used in the assembly matrix
	/// the number of elements in the k vectors is N_rows*N_cols
	///------------------------------------------------------------------------------
	void calculate_k_values_for_assembly_matrix(int NRows, int NCols, vector<int>& k_value_i_j,
	                    vector<int>& k_value_ip1_j,	vector<int>& k_value_im1_j, vector<int>& k_value_i_jp1, 
	                    vector<int>& k_value_i_jm1);
	                    
	/// mtl_initiate_assembler_matrix                      
	void mtl_initiate_assembler_matrix(int& problem_dimension,				     
			     float& inv_dx_S_c_squared, float& inv_dy_S_c_squared, float& dx_front_term, 
	             float& dy_front_term, vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,
			     vector<int>& vec_k_value_im1_j, vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1);                    					
	//------------------------------------------------------------------------------
	/// mtl_assemble_matrix
	/// this function assembles the solution matrix for nonlinear creep transport
	//------------------------------------------------------------------------------
	void mtl_assemble_matrix(Array2D<float>& zeta_last_iter, Array2D<float>& zeta_last_timestep,
						 Array2D<float>& zeta_this_iter, Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate,
	           mtl::compressed2D<float>& mtl_Assembly_matrix, mtl::dense_vector<float>& mtl_b_vector,
						 float dt, int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared, 
					   float dx_front_term, float dy_front_term,
	           float South_boundary_elevation, float North_boundary_elevation,
	           vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,vector<int>& vec_k_value_im1_j,
						 vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1);
	///------------------------------------------------------------------------------
	/// mtl_solve_assembler_matrix
	/// this function assembles the solution matrix
	///------------------------------------------------------------------------------
	void mtl_solve_assembler_matrix(Array2D<float>& zeta_last_iter, Array2D<float>& zeta_last_timestep,
						 Array2D<float>& zeta_this_iter, Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate,
						 float dt, int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
					   float dx_front_term, float dy_front_term,
	           vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,
	           vector<int>& vec_k_value_i_jp1, std::vector<int>& vec_k_value_i_jm1,
	           float South_boundary_elevation, float North_boundary_elevation);
	///------------------------------------------------------------------------------
	/// nonlinear_creep_timestep
	/// do a creep timestep.  This function houses the above two functions to
	/// undertake model timestep using implicit implementation of the nonlinear
	/// transport law.
	/// NOTE you need to run mtl_initiate_assembler_matrix before you run this function
	///------------------------------------------------------------------------------
	void nonlinear_creep_timestep(Array2D<float>& fluvial_erosion_rate, float iteration_tolerance,
				int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared,
				float dx_front_term, float dy_front_term, vector<int>& vec_k_value_i_j,
				vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,
				vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1,
				float South_boundary_elevation, float North_boundary_elevation);

	/// -----------------------------------------------------------------------------
	/// Soil diffusion method
	/// Container for all the finite volume components
	/// -----------------------------------------------------------------------------
	void soil_diffusion_fv( void );
	
	/// -----------------------------------------------------------------------------
	/// Finite difference matrix
	/// -----------------------------------------------------------------------------
	mtl::compressed2D<float> generate_fd_matrix( int dimension, int size, bool periodic );
	mtl::dense_vector <float> build_fd_vector( int dimension, int size );
	//void repack_fd_vector(mtl::dense_vector <float> &data_vector, int dimension); 

	mtl::compressed2D<float> generate_fv_matrix( int dimension, int size, bool periodic );
	mtl::dense_vector <float> build_fv_vector( int dimension, int size );
	void repack_vector(mtl::dense_vector <float> &data_vector, int dimension); 

	/// -----------------------------------------------------------------------------
	/// Soil diffusion using linear flux model
	/// Solved using finite difference
	/// -----------------------------------------------------------------------------
	void soil_diffusion_fd_linear( void );

	void soil_diffusion_fv_nonlinear( void );

	///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	/// RUN MODEL
	///------------------------------------------------------------------------------
	/// A series of wrapper functions that implement the numerical model
	///------------------------------------------------------------------------------
	/// implicit_hillslope_and_fluvial
	/// This function sets up a landscape evolution model run incorporating fluvial
	/// erosion and hillslope erosion via non-linear creep.  It calls the implicit
	/// implementation and returns the topography after the final timestep.
	/// The user should provide the parameter file which sets out the details of the 
	/// model run.   
	///------------------------------------------------------------------------------
	LSDRasterModel run_model_implicit_hillslope_and_fluvial(string param_file);
		
  /// @brief This wrapper just calls the run_components method. Parameters used are those
  /// stored as data members
	/// @author JAJ
	/// @date 01/01/2014	
	void run_model( void );
	
	void run_model_from_steady_state( void );
	
  /// @brief This wrapper just calls the run_components method. Parameters used are those
  /// stored as data members. This one actually calls the erosion laws
	/// @author JAJ
	/// @date 01/01/2014		
	void run_components( void );

	/// This method forces the landscape into its steady state profile, by using periodic forcing. 
	/// This is much more efficient than using static forcing (as in run model), but doesn't give
	/// a nice animation of an evolving landscape
	/// Swings and roundabouts
	void reach_steady_state( void );

	/// -----------------------------------------------------------------------
	/// Reset model - reset erosion values to 0 after a complete model run
	/// -----------------------------------------------------------------------
	void reset_model( void );

	/// -----------------------------------------------------------------------
	/// Fastscape, implicit finite difference solver for stream power equations
	/// O(n)
	/// Method takes a series of paramaters and solves the stream power equation
	/// At a future timestep in linear time
	/// -----------------------------------------------------------------------
	void fluvial_incision( void );
  
  /// @brief This assumes that all sediment transported from rivers into 
  /// channels is removed. It checks the raster to see where the channels are
  /// which at this point is determined by a threshold drainage area, 
  /// and then removes all the sediment to those pixels
	/// @author JAJ
	/// @date 01/01/2014	
	void wash_out( void );
	
	/// fluvial_erosion_rate
	LSDRaster fluvial_erosion_rate(float timestep, float K, float m, float n, vector <string> boundary);

	/// --------------------------------------------------------------------
	/// Runs flexural isostatic calculations
	/// Uses fourier filtering method
	/// Pelletier (2008)
	/// --------------------------------------------------------------------
	LSDRasterModel run_isostatic_correction( void );

	/// --------------------------------------------------------------------
  /// @brief Adds random noise to each pixel in range [min, max]
	/// @param minium random addition
	/// @param maximum random addition
	/// @author JAJ
	/// @date 01/01/2014
	void random_surface_noise( float min, float max );

  /// --------------------------------------------------------------------
  /// @brief Adds random noise to each pixel using the noise data member
	/// @author SMM
	/// @date 17/06/2014
	void random_surface_noise();

	/// --------------------------------------------------------------------
	/// Creates uplift field from a set of templates
	/// int mode specifies a mode of uplift:
	/// 	(0) - block uplift
	///	1   - tilt block
	///	2   - gaussian
	///	3   - quadratic
	/// --------------------------------------------------------------------
	Array2D <float> generate_uplift_field( int mode, float max_uplift); 


	/// -------------------------------------------------------------------
	/// Gets the uplift value at a given cell
	/// this method is implemented as a memory saving measure, rather than
	/// storing the uplift field in memory
	///
	/// Some methods still implemented still use this uplift field
	/// It's advisable this is changed, otherwise the size of rasters that
	/// can be modelled will be severely reduced
	/// -------------------------------------------------------------------
	float get_uplift_at_cell(int i, int j);

	/// -------------------------------------------------------------------
	/// Correct for isostasy using Airy model 
	/// -------------------------------------------------------------------
	
	void Airy_isostasy( void );

	/// -------------------------------------------------------------------
	/// Correct for isostasy using flexural model 
	/// -------------------------------------------------------------------
	
	void flexural_isostasy( float alpha );

	void flexural_isostasy_alt( void );

	/// -------------------------------------------------------------------
	/// Calculates depth of topographic root using FFT methods inherited from LSDRasterSpectral
	/// -------------------------------------------------------------------
	Array2D <float> calculate_root( void );
	Array2D <float> calculate_airy( void );

	/// -------------------------------------------------------------------
	/// Check whether current node is a base level node
	/// -------------------------------------------------------------------
	bool is_base_level(int i, int j);

	void interpret_boundary(short &dimension, bool &periodic, int &size);

	/// -------------------------------------------------------------------
	/// Setter methods 
	/// -------------------------------------------------------------------
	
	void set_boundary_conditions(vector <string> bc) 	{ for (int i=0; i<4; ++i) {bc[i][0] = tolower(bc[i][0]);} boundary_conditions = bc; }
	void set_timeStep( float dt )				{ timeStep = dt; }
	void set_endTime( float time )				{ endTime = time; }
	void set_num_runs( int num )				{ num_runs = num; }
	void set_uplift( Array2D <float> uplift )		{ uplift_field = uplift; }
	void set_uplift( int mode, float max_rate )		{ uplift_field = generate_uplift_field( mode, max_rate ); this->max_uplift = max_rate; }
	void set_steady_state_tolerance( float tol )		{ steady_state_tolerance = tol; }
	void set_K( float K )					{ this->K_fluv = K; }
	void set_D( float D )					{ this->K_soil = D; }
	void set_rigidity( float D )				{ this->rigidity = D; }
	void set_m( float m )					{ this->m = m; }
	void set_n( float n )					{ this->n = n; }
	void set_threshold_drainage( float area )		{ this->threshold_drainage = area; }
	void set_S_c( float degrees )				{ S_c = tan(degrees*3.14159265358/180); }
	void set_periodicity( float time )			{ periodicity = time; }
	void set_periodicity_2( float time )			{ periodicity_2 = time; }
	void snap_periodicity( void );
	void set_print_interval( int num_steps )		{ print_interval = num_steps; }
	void set_K_mode( short mode )				{ K_mode = mode; }
	void set_D_mode( short mode )				{ D_mode = mode; }
	void set_period_mode( short mode )			{ period_mode = mode; }

	void set_name( string name )				{this->name = name;}
	void set_report_name( string name )			{report_name = name;}

	// Setters for turning on/off model components
	void set_fluvial( bool on_status )			{ fluvial = on_status; }
	void set_hillslope( bool on_status )			{ hillslope = on_status; }
	void set_isostasy( bool on_status )			{ isostasy = on_status; }
	void set_flexure( bool on_status )			{ flexure = on_status; }
	void set_quiet( bool on_status )			{ quiet = on_status; }

	/// -------------------------------------------------------------------
	/// Getter methods 
	/// -------------------------------------------------------------------
	
	string get_name( void )					{ return name; }
	float get_K( void );
	float get_D( void );
	float get_max_uplift( void );			// NOTE; while this is currently a very trivial, and arguably unecessary method, it should be used
							// and developed if someone wants to integrate some sort of changing uplift field
	float find_max_boundary( int boundary_number );
	
	/// ------------------------------------------------------------------
	/// Display list of parameters
	/// ------------------------------------------------------------------
	void print_parameters( void );

	/// ------------------------------------------------------------------
	/// Write model reports
	/// ------------------------------------------------------------------
	void write_report( void );
	void cycle_report( float, float, float);

	void final_report( void );

	void print_rasters( int frame_num );

	/// ------------------------------------------------------------------
	/// Print slope area data
	/// Probably fits better into LSDRaster, but requires LSDFlowInfo
	/// ------------------------------------------------------------------
	void slope_area_data( string name );

	/// ------------------------------------------------------------------
	/// Produce a template of a parameter file to be supplied to the model
	/// ------------------------------------------------------------------
	void make_template_param_file(string filename);

	/// ------------------------------------------------------------------
	/// Display method
	/// Wrapper for python script that animates the model output
	/// ------------------------------------------------------------------
	void show( void );

	/// ----------------------------------------------------------------
	/// Dave's stuff
	/// ----------------------------------------------------------------
	void DAVE_initiate_assembler_matrix(int& problem_dimension, float& inv_dx_S_c_squared, float& inv_dy_S_c_squared, float& dx_front_term, 
               float& dy_front_term, vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,
					     vector<int>& vec_k_value_im1_j, vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1);
	void DAVE_calculate_k_values_for_assembly_matrix(int NRows, int NCols, vector<int>& k_value_i_j,
                      vector<int>& k_value_ip1_j,	vector<int>& k_value_im1_j, vector<int>& k_value_i_jp1, vector<int>& k_value_i_jm1);											
	void DAVE_nonlinear_creep_timestep(Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate, float iteration_tolerance, int problem_dimension,
				float inv_dx_S_c_squared, float inv_dy_S_c_squared, float dx_front_term, float dy_front_term, vector<int>& vec_k_value_i_j,
				vector<int>& vec_k_value_ip1_j, vector<int>& vec_k_value_im1_j,	vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1,
					  float South_boundary_elevation, float North_boundary_elevation);
	void DAVE_solve_assembler_matrix(Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate, int problem_dimension, float inv_dx_S_c_squared,
			float inv_dy_S_c_squared, float dx_front_term, float dy_front_term, vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,
			vector<int>& vec_k_value_im1_j, vector<int>& vec_k_value_i_jp1, std::vector<int>& vec_k_value_i_jm1, float South_boundary_elevation, float North_boundary_elevation);
	void DAVE_assemble_matrix(Array2D<float>& uplift_rate, Array2D<float>& fluvial_erosion_rate, mtl::compressed2D<float>& mtl_Assembly_matrix, mtl::dense_vector<float>& mtl_b_vector,
						 int problem_dimension, float inv_dx_S_c_squared, float inv_dy_S_c_squared, float dx_front_term, float dy_front_term,
		     float South_boundary_elevation, float North_boundary_elevation, vector<int>& vec_k_value_i_j, vector<int>& vec_k_value_ip1_j,vector<int>& vec_k_value_im1_j,
						 vector<int>& vec_k_value_i_jp1, vector<int>& vec_k_value_i_jm1);
	void DAVE_wrapper( void );



	void write_root(string name, string ext);


		
	protected:
	// Various parameters used in throughout the model run
	// These are set to the default values with default_parameters( void )
	// and then a parameter file is read to overwrite any of these explicitly set
	bool			quiet;				// Supress output messages
	bool			initialized;			// Check whether initialize_model has been run
	bool			steady_state;			// Check whether initial steady state has been arrived at
	bool			initial_steady_state;		// Switch for first steady state arrival, used for activating periodic forcing parameters
	bool			cycle_steady_check;		// Check for steady state using periodic forcing (not sure about this one)
	bool			recording;			// Switch for start of recording data
	bool			reporting;			// Whether or not to write the run report
	vector <string>		boundary_conditions;		// Boundary conditions of model NESW
	string			name;				// Name of the model run
	string 			report_name;
	float 			current_time;			// Current time
	float			time_delay;			// Time at which initial steady state was reached
	float 			timeStep;			// Time in between each calculation
	float			endTime;			// Time at end of model run
	short			endTime_mode;
	int			num_runs;
	Array2D <float> 	uplift_field;			// Field of uplift for each cell
	int 			uplift_mode;
	float			max_uplift;
	float			steady_state_tolerance;
	float			steady_state_limit;
	float 			m, n;				// Dimensionless exponents for SPL
	float			K_fluv, K_soil;			// Fluvial and soil coefficients
	float			threshold_drainage;		// Drainage area above which soil will be flushed from the system
	float			S_c;				// Critical slope (for non-linear soil creep)
	float			rigidity;			// Flexural rigidity of plate
	int			print_interval;			// Interval at which output is written
	bool			print_elevation;		
	bool			print_erosion;			// Whether or not to print the erosion field at each print interval
	bool			print_erosion_cycle;
	bool			print_hillshade;		// Whether or not to print a hillshade raster instead of topographic	
	bool			print_slope_area;		// Whether or not to produce a slope area report 
	Array2D <float>		root_depth;			// Depth of topographic root
	//	Measures of landscape response
	float			erosion;
	float			erosion_last_step;
	vector<float>		erosion_cycle_record;
	float			total_erosion;			// Offset from uplift from each cell
	float			min_erosion;
	float 			max_erosion;
	float 			response;			// maximum response over a single run
	float			total_response;			// response over all model runs (to be divided by num_runs)
	float			noise;
	float			report_delay;

	Array2D <float>		zeta_old;
	Array2D <float>		steady_state_data;
	Array2D <float>		erosion_cycle_field;

	// Parameters for periodic forcing components
	short K_mode;			// Whether or not K is periodic
	short D_mode;			// Whether or not D is periodic
	short period_mode;		// Whether or not there is only one periodicity or two
					// 1 (default) one periodicity used without
					// 2 Two periodicities that switch at a given interval
					// 3 Two periodicities used as a compound sin wave
					// 4 Same as three, but weightings switch at a given interval (as in 2)

	float K_amplitude;		// Amplitude of K wave
	float D_amplitude;		// Amplitude of D wave
	float periodicity;		// Periodicty
	float periodicity_2;		// 2nd periodicity (if period_mode = 2)
	int   cycle_number;
	float p_weight;			// Ratio for weight of periodicity in period_mode 3 or 4
	float switch_time;		// Time at which switch happens (time mode is same as endTime_mode)
	float switch_delay;		// Similar to time delay (but for switching periodicities)

	// Components switches
	bool			fluvial;
	bool			hillslope;
	bool			nonlinear;			// Whether nonlinear mode is used
	bool			isostasy;
	bool			flexure;			// Whether flexural isostasy will be used

  Array2D<float> zeta_last_iter;
  Array2D<float> zeta_last_timestep;
  Array2D<float> zeta_this_iter;


	private:
	void create();
	void create(string master_param);
	void create(string filename, string extension);
	void create(int ncols, int nrows, float xmin, float ymin,
		    float cellsize, float ndv, Array2D<float> data);
	void create(LSDRaster& An_LSDRaster);
	void default_parameters( void );

	float periodic_parameter( float base_param, float amplitude );
	float square_wave_parameter( float base_param, float amplitude );
	float stream_K_fluv(void);		// I don't like that these two are split up, but I was rushed and couldn't figure out how to pass an ifstream
	float stream_K_soil(void);
	void _run_components( void );
};

#endif


//  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  /// INITIAL TOPOGRAPHY MODULE
//  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  /// There are a number of options that will be built into this to give 
//  /// flexibility when it comes to setting up the LEM.
//  /// This includes:
//  /// i) Artifial parabolic surfaces - other surfaces may be added later as 
//  /// needed
//  /// ii) Existing topographic datasets (for example a real DEM) - this is simply
//  /// done by loading the data as LSDRasterModel rather than an LSDRaster
//  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  /// Create_parabolic_surface
//  /// This function generates a parabolic surface with elevations on the North
//  /// and South edges at zero, and in the middle at 'PeakElevation'
//  LSDRasterModel create_parabolic_surface(int NRows, int NCols, float dx, 
//    float PeakElevation, float RandomAmplitude, float EdgeOffset);
//  
//  ///----------------------------------------------------------------------------
//  /// Create Nonlinear Steady State Hillslope
//  /// Does what it says on the tin
//  LSDRasterModel create_nonlinear_SS_hillslope(float K_nl, float S_c, float U);
	   
//   //------------------------------------------------------------------------------
//   // nonlinear_diffusion - MIGHT STILL BE A WORK IN PROGRESS - CHECK WITH SIMON
//   // calculate fluxes using the nonlinear soil flux model - use the buffered
//   // topographic raster
//   //------------------------------------------------------------------------------
//   LSDRasterModel nonlinear_diffusion(Array2D<float>& SlopesBetweenRows, Array2D<float>& SlopesBetweenColumns, 
//                 float K_nl, float S_c, float dt);
//                 
//   ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//   /// EXPLICIT MODEL COMPONENTS
//   ///------------------------------------------------------------------------------
//   /// A series of functions that carry out parts of the number crunching in the
//   /// explicit version of the model
//   ///------------------------------------------------------------------------------
//   /// creep_and_fluvial_timestep
//   /// For a given timestep in the model, calculates the hillslope and channel
//   /// sediment fluxes and adjusts topography accordingly
//   ///------------------------------------------------------------------------------
//   LSDRasterModel creep_and_fluvial_timestep(float& t_ime, float dt, float uplift_rate,
// 								float South_boundary_elevation, float North_boundary_elevation,
// 								float D_nl, float S_c, float k_w, float b, float K, float n, float m, 
//                 float erosion_threshold, LSDRasterModel& ZetaOld, LSDRasterModel& ZetaRasterBuff,
// 								Array2D<float>& SlopesBetweenRows, Array2D<float>& SlopesBetweenColumns,
// 								Array2D<float>& ErosionRateArray, Array2D<float>& precip_flux,
//                 Array2D<float>& Q_w, Array2D<float>& ChannelWidthArray, 
//                 Array2D<float>& TopoDivergence, Array2D<float>& FluvialErosionRateArray);   

	       
