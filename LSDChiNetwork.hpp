//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiNetwork.hpp
// header file for the LSDChiNetwork object
// this object is used to examine a network of channels in chi space
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		28/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <vector>
#include <string>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;


#ifndef LSDChiNetwork_H
#define LSDChiNetwork_H

class LSDChiNetwork
{
	public:
		LSDChiNetwork(string channel_network_fname)
					{ create( channel_network_fname ); }

		int get_n_channels()	{ return int(node_indices.size()); }

		// get functions. These are used for interfacing with
		// the LSDRaster object (not in the standalone version)
		int get_NRows() const				{ return NRows; }
		int get_NCols() const				{ return NCols; }
		double get_XMinimum() const			{ return XMinimum; }
		double get_YMinimum() const			{ return YMinimum; }
		double get_DataResolution() const	{ return DataResolution; }
		double get_NoDataValue() const		{ return NoDataValue; }

		// printing routines for bug checking
		void print_channel_details_to_screen(int channel_number);
		void print_channel_details_to_file(string fname, double A_0, double m_over_n);
		void print_channel_details_to_file_full_fitted(string fname);
		void print_channel_details_to_file_full_fitted(string fname, int target_nodes,
		                                               int minimum_segment_length);

		// this extends the tributary channels all the way to the outlet
		// In its current version this only works if the tributaries all drain
		// to the mainstem
		void extend_tributaries_to_outlet();

		// routine for returning calculated data to an array
		Array2D<double> data_to_array(int data_member);

		// routines for getting slope-area data
		// these print to file at the moment.
		void slope_area_extraction_vertical_intervals(double interval, double area_thin_fraction,
															string fname);
		void slope_area_extraction_horizontal_intervals(double interval, double area_thin_fraction,
															string fname);

		// routines for calculating chi and maniplating chi
		void calculate_chi(double A_0, double m_over_n);
		double calculate_optimal_chi_spacing(int target_nodes);
		int calculate_skip(int target_nodes);
		int calculate_skip(int target_nodes, vector<double>& sorted_chis);
		int calculate_skip(int target_nodes, int channel_number);

		// routines for calcualting the most likeley segments.
		void find_most_likeley_segments(int channel,int minimum_segment_length,
						 double sigma, int N, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc );
		void find_most_likeley_segments_dchi(int channel,int minimum_segment_length,
						 double sigma, double dchi, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc );
		void find_most_likeley_segments_monte_carlo(int channel, int minimum_segment_length,
						 double sigma, int mean_skip, int skip_range, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc );
		void find_most_likeley_segments_monte_carlo_dchi(int channel, int minimum_segment_length,
						 double sigma, double mean_dchi, double variation_dchi, vector<double>& b_vec,
					     vector<double>& m_vec, vector<double>& r2_vec,vector<double>& DW_vec,
					     vector<double>& thinned_chi, vector<double>& thinned_elev,
					     vector<double>& fitted_elev, vector<int>& node_reference,
					     vector<int>& these_segment_lengths,
					     double& this_MLE, int& this_n_segments, int& n_data_nodes,
					     double& this_AIC, double& this_AICc );

		// the master routine for calculating the best fit m over n values for a channel network
		// this is based on a fixed value of dchi
		double search_for_best_fit_m_over_n_dchi(double A_0, int n_movern, double d_movern, double start_movern,
						       int minimum_segment_length, double sigma, int target_nodes_mainstem, string fname);

		// the master routine for calculating the best fit m over n values for a channel network
		double search_for_best_fit_m_over_n(double A_0, int n_movern, double d_movern, double start_movern,
						       int minimum_segment_length, double sigma, int target_nodes_mainstem, string fname);

		// similar to above, but calculates the mainstem and the tributaries seperately.
		double search_for_best_fit_m_over_n_seperate_ms_and_tribs(double A_0, int n_movern, double d_movern, double start_movern,
						       int minimum_segment_length, double sigma, int target_nodes_mainstem, string fname);

		// this function randomly selects data from the mainstem and the channels and creates
		// a composite channel. This composite channel is then tested for linear segments.
		// the algorithm therefere tries to maximize both colinearity and linearity of the channel
		// network.
		double search_for_best_fit_m_over_n_colinearity_test(double A_0, int n_movern, double d_movern,
								        double start_movern, int minimum_segment_length, double sigma,
						      			int target_nodes, int n_iterations,
						      			vector<double>& m_over_n_values,
						      			vector<double>& AICc_mean, vector<double>& AICc_sdtd);

		// this function is similar to the one below above but allows breaks in the channel
		double search_for_best_fit_m_over_n_colinearity_test_with_breaks(double A_0, int n_movern, double d_movern,
						       double start_movern, int minimum_segment_length, double sigma,
						       int target_skip, int target_nodes, int n_iterations,
						       vector<double>& m_over_n_values, vector<double>& AICc_mean, vector<double>& AICc_sdtd,
						       int Monte_Carlo_switch);


		// this gets the best fit m over n values of all the individual tributaries
		double search_for_best_fit_m_over_n_individual_channels_with_breaks(double A_0, int n_movern, double d_movern,
								        double start_movern, int minimum_segment_length, double sigma,
						      			int target_skip, int target_nodes, int n_iterations,
						      			vector<double>& m_over_n_values, vector< vector<double> >& AICc_vals);

		// this gets the best fit m over n values of all the individual tributaries
		// it uses a monte carlo appraoach so all tributaries have botht the mean and the variability
		// of the AICc values reported
		double search_for_best_fit_m_over_n_individual_channels_with_breaks_monte_carlo(double A_0, int n_movern,
							   double d_movern,double start_movern, int minimum_segment_length, double sigma,
						       int target_skip, int target_nodes, int n_iterations,
						       vector<double>& m_over_n_values,
						       vector< vector<double> >& AICc_means, vector< vector<double> >& AICc_stddev);


		// this routine uses a monte carlo approach to repeatedly sampling all the data in the
		// channel network using a reduced number of data elements and then poulating each channel node with
		// a distribution of m, b and fitted elevation values. These then can be averaged and details of their variation
		// calculated.
		// this is based on a fixed value of dchi
		void monte_carlo_sample_river_network_for_best_fit_dchi(double A_0, double m_over_n, int n_iterations,
															double fraction_dchi_for_variation,
															int minimum_segment_length, double sigma,
															int target_nodes_mainstem);

		// this routine uses a monte carlo approach to repeatedly sampling all the data in the
		// channel network using a reduced number of data elements and then poulating each channel node with
		// a distribution of m, b and fitted elevation values. These then can be averaged and details of their variation
		// calculated.
		void monte_carlo_sample_river_network_for_best_fit(double A_0, double m_over_n, int n_iterations,
															int mean_skip, int skip_range,
															int minimum_segment_length, double sigma);


		// this routine uses a monte carlo approach to repeatedly sampling all the data in the
		// channel network using a reduced number of data elements and then poulating each channel node with
		// a distribution of m, b and fitted elevation values. These then can be averaged and details of their variation
		// calculated.
		void monte_carlo_sample_river_network_for_best_fit_after_breaks(double A_0, double m_over_n, int n_iterations,
				int skip, int minimum_segment_length, double sigma);

		// this function uses a monte carlo sampling approach to try and split channels into
		// so that the channel is sampled at the target skipping interval
		void monte_carlo_split_channel(double A_0, double m_over_n, int n_iterations,
				int target_skip, int target_nodes,
				int minimum_segment_length, double sigma, int chan, vector<int>& break_nodes);

		// this function uses a monte carlo sampling approach to try and split channels into
		// so that the channel is sampled at the target skipping interval. It does it with a
		// colinear dataset
		void monte_carlo_split_channel_colinear(double A_0, double m_over_n, int n_iterations,
				int target_skip, int target_nodes,
				int minimum_segment_length, double sigma,
				vector<double> reverse_Chi, vector<double> reverse_Elevation, vector<int>& break_nodes);

		// this function splits all the channels in one go
		void split_all_channels(double A_0, double m_over_n, int n_iterations,
				int target_skip, int target_nodes, int minimum_segment_length, double sigma);

		// this function gets the AICc after breaking the channel
		double calculate_AICc_after_breaks(double A_0, double m_over_n,
				int skip, int minimum_segment_length, double sigma, int chan, vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE);

		// this function gets the AICc after breaking the channel
		// it does this for n_iterations and returns a vector with all of
		// the AICc values for each iteration reported. This can then be used
		// to calculate the statistics of the AICc to tell if the minimum
		// AICc is significantly different from the other AICc values for different
		// values of m/n
		vector<double> calculate_AICc_after_breaks_monte_carlo(double A_0, double m_over_n,
				int target_skip, int minimum_segment_length, double sigma, int chan, vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE,
				int n_iterations);

		// this function gets the AICc after breaking the channel with a colinear dataset
		// the reverse_chi and reverse_elevation data has to be provided
		double calculate_AICc_after_breaks_colinear(double A_0, double m_over_n,
						int skip, int minimum_segment_length, double sigma,
						vector<double> reverse_chi, vector<double> reverse_elevation,
						vector<int> break_nodes,
						int& n_total_segments, int& n_total_nodes, double& cumulative_MLE);

		// this function gets the AICc after breaking the channel with a colinear dataset
		// the reverse_chi and reverse_elevation data has to be provided
		//
		// it uses a monte carlo scheme and returns a vector with all of the AICc values calcluated
		// from the analyses
		vector<double> calculate_AICc_after_breaks_colinear_monte_carlo(double A_0, double m_over_n,
				int skip, int minimum_segment_length, double sigma,
				vector<double> reverse_Chi, vector<double> reverse_Elevation,
				vector<int> break_nodes,
				int& n_total_segments, int& n_total_nodes, double& cumulative_MLE,
				int n_iterations);

		// this routine test to see if channels are long enough to get a decent fitting from
		// the segment finding algorithms
		void is_channel_long_enough_test(int minimum_segment_length,int N);

	protected:
		// data for georeferencing
		int NRows;			// number of rows
		int NCols;			// number of columns
		double XMinimum;
		double YMinimum;

		// metadata
		double DataResolution;
		double NoDataValue;


		vector< vector<int> > node_indices;			// node indices: used in conjunction with
													// other LSD topographic tool objects and not
													// necessary for standalone program
		vector< vector<int> > row_indices;			// row indices: used in conjunction with
													// other LSD topographic tool objects and not
													// necessary for standalone program
		vector< vector<int> > col_indices;			// column indices: used in conjunction with
													// other LSD topographic tool objects and not
													// necessary for standalone program

		vector< vector<double> > elevations;		// the elevations along the channels
		vector< vector<double> > flow_distances;	// flow distances along channels. Used to integrate
													// to arrive at chi
		vector< vector<double> > drainage_areas;

		vector< vector<double> > chis;				// the chi values for the channels. This data will
													// be overwritten as m_over_n changes

		vector<int> node_on_receiver_channel;		// this is the node on the reciever channel where the
													// tributary enters the channel. Used to find the downstream
													// chi value of a channel
		vector<int> receiver_channel;				// this is the channel that the tributary enters.


		// the following data elements are fitted ci slopes and intercepts, as well
		// as best fit elevation that are generated after monte-carlo sampling
		// they are only filled with data after calling the
		// monte_carlo_sample_river_network_for_best_fit function
		double m_over_n_for_fitted_data;			// this stores the m over n value use to generate the
													// means and standard deviations of the network properties
		double A_0_for_fitted_data;
		vector<int> is_tributary_long_enough;		// this vector is the same size as
													// the number of channels and is 1 if the channel
													// analysis n_nodes > 3* minimum_segment_length
		vector< vector<double> > chi_m_means;
		vector< vector<double> > chi_m_standard_deviations;
		vector< vector<double> > chi_m_standard_errors;
		vector< vector<double> > chi_b_means;
		vector< vector<double> > chi_b_standard_deviations;
		vector< vector<double> > chi_b_standard_errors;
		vector< vector<double> > all_fitted_elev_means;
		vector< vector<double> > all_fitted_elev_standard_deviations;
		vector< vector<double> > all_fitted_elev_standard_errors;
		vector< vector<double> > chi_DW_means;
		vector< vector<double> > chi_DW_standard_deviations;
		vector< vector<double> > chi_DW_standard_errors;
		vector< vector<int> > n_data_points_used_in_stats;		// the parameters are generated using a
																// monte carlo approach and not all nodes will
																// have the same number of data points, so the number
																// of data points is stored

		vector< vector<int> > break_nodes_vecvec;				// this vector holds the vectors containing
																// the node locations of breaks in the segments

	private:
		void create(string channel_network_fname);
};

#endif
