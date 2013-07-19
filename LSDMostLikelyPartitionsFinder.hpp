//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDMostLikelyPartitionsFinder.hpp
// header file for the LSDMostLikelysPartitionFinder object
// this object looks for the most likeley partitions or segments
// of 2D data, and is principally used to identify segments of
// differing channel steepness in chi-zeta space
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// <your name here>
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		03/01/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <vector>
#include "LSDStatsTools.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;


#ifndef LSDMostLikelyPartitionsFinder_H
#define LSDMostLikelyPartitionsFinder_H

class LSDMostLikelyPartitionsFinder
{
	public:
		LSDMostLikelyPartitionsFinder(int this_min_seg_length, vector<double> this_x_data, vector<double> this_y_data)
					{ create( this_min_seg_length, this_x_data, this_y_data); }

		// some get functions
        int get_n_nodes()                    {return int(x_data.size()); }
		vector<double> get_MLE_of_segments()	{return MLE_of_segments;}
        vector<double> get_x_data()             {return x_data; }
        vector<double> get_y_data()             {return y_data; }

		// functions for thinning the data
		void reset_derived_data_members();

		void thin_data_target_dx_preserve_data(double dx);
        void thin_data_target_dx_preserve_data(double dx, vector<int>& node_ref);
		LSDMostLikelyPartitionsFinder spawn_thinned_data_target_dx_preserve_data(double dx);

		void thin_data_skip(int N, vector<int>& node_ref);

		void thin_data_target_dx_linear_interpolation(double dx);
		LSDMostLikelyPartitionsFinder spawn_thinned_data_target_dx_linear_interpolation(double dx);


        void thin_data_monte_carlo_skip(int Mean_skip,int skip_range, vector<int>& node_ref);


		void thin_data_monte_carlo_dchi(double mean_dchi, double variation_dchi, vector<int>& node_ref);


		// function for lookg at the x and y data
		void print_x_y_data_to_screen();

		// this function drives the whole shebang
		void best_fit_driver_AIC_for_linear_segments(vector<double> sigma_values);
		void best_fit_driver_AIC_for_linear_segments(double sigma_values);
		void get_data_from_best_fit_lines(int node, vector<double> sigma_values,
					      	vector<double>& b_values, vector<double>& m_values,
					      	vector<double>& r2_values, vector<double>& DW_values,
						vector<double>& fitted_y,vector<int>& seg_lengths,
						double& this_MLE, int& this_n_segments, int& this_n_nodes,
                                                double& this_AIC, double& this_AICc);
		void get_start_and_end_x_for_segments(vector<double>& start_x, vector<double>& end_x, vector<int> seg_lengths);


		// these functions populate the arrays used for calcualting the best fit segments
		void calculate_segment_matrices(double sigma);
		void populate_segment_matrix(int start_node, int end_node, double no_data_value,double sigma);

		// get the maximum likilihood of all the possible number of segments
		void find_max_like_of_segments();

		// this drives the partitioning algorithm
		void partition_driver_to_vecvecvec(int k);

		// functions that perform components of the partioning
		void partitions_with_minimum_length(int n, int k, int t, vector<int>& p);
		void partition_assign(int t, vector<int>& p);
		int LSDpartitions_min( int x, int y);

		// functions for working with liklihoods. Can transform likelihoods between different sigma values
		Array2D<double> normalize_like_matrix_to_sigma_one(double sigma);
		vector<double> normalize_like_vector_to_sigma_one(double sigma, vector<double> like_vector);
		void change_normalized_like_matrix_to_new_sigma(double sigma, Array2D<double>& sig1_like_array);
		vector<double> change_normalized_like_vector_to_new_sigma(double sigma, vector<double> sig1_like_vector);
		vector<double> transform_like_from_sigma1_to_sigma2(double sigma1,
													vector<double> sig1_like_vector, double sigma2);

		// these get the best fit number of segments for a variety of sigma values.
		// because of the way the AIC and AICc algorithms work, the minimum AIC for different number of segments
		// can vary depending on sigma, which we don't know. Therefore, we include this function to scan for best fits across
		// different values fo sigma
		void get_n_segments_for_various_sigma(vector<double> sigma_values);
		void calculate_AIC_of_segments_with_variable_sigma(double sigma,
										vector<double>& AIC_of_segments,
										vector<double>& AICc_of_segments);

		// this function extracts the m, b, r^2 and DW statistic of the most likeley segments
		void get_properties_of_best_fit_segments(int bestfit_segments_node,
										 vector<double>& m_values, vector<double>& b_values,
										 vector<double>& r2_values, vector<double>& DW_values);

		// some functions for printing results to screen
		void print_to_screen_most_likeley_segment_lengths();
		void print_AIC_and_AICc_to_screen(vector<double> sigma_values);

	protected:
		int minimum_segment_length;
		vector<double> x_data;
		vector<double> y_data;

		double base_sigma;					// the base sigma value from which the MLE of the segments
											// is calcluated

		// arrays containing the properties of the segments
		// the arrays are indexed so the row is the starting node and the column is the ending node
		Array2D<double> like_array;					// liklihood
		Array2D<double> m_array;					// slope
		Array2D<double> b_array;					// intercept
		Array2D<double> rsquared_array;				// r^2
		Array2D<double> DW_array;					// Durbin-Watson statistic to test if the residuals
													// are autocorrelated, used to determine if the
													// segment is truly linear. Values less than 1
													// indicate that the segment is probably not linear
													// values < 1.5 should arouse suspicion.

		vector<double> MLE_of_segments;				// the maximum likelihood of the different number of
													// segments
		vector< vector<int> > segments_for_each_n_segments;
													// each element of this vector contains the
													// most likeley segments for that number of segments
													// so for example segments_for_each_n_segments[3]
													// is a vector containing the lengths of the most likeley
													// segments for 4 segments (note 0 indexing)
		vector< vector < vector<int> > > partitions;
													// this vecvecvec: top index references the number of
													// segments. Second layer loopps through the possible
													// partitions of that number of segments. Third layer
													// is the individual segment lenghts.

		vector<int> best_fit_AIC;
		vector<int> best_fit_AICc;
		vector< vector<double> > AIC_for_each_n_segments;
		vector< vector<double> > AICc_for_each_n_segments;

	private:
		void create(int this_min_seg_length, vector<double> this_x_data, vector<double> this_y_data);
};

#endif
