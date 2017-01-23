//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// m_chi_from_cosmo_basins_driver.cpp
//
// This program gets the mean M chi value for each basin with a catchment-averaged erosion
// rate. It takes an input parameter file with the m/n value, A_0 and 
// It also needs a text file with the details of the cosmo points to be read in:
// X_coordinate Y_coordinate Erosion_Rate Error
// This file must be called "DEM_name_cosmo.txt"
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp"
#include "../TNT/tnt.h"


int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the driver file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];

	cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the header file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}

	string DEM_name;
	//string dem_ext = "_dem";
  string DEM_flt_extension = "flt";
  string fill_ext = "_fill";
  string sources_ext = "_CH";
	float pruning_threshold;
	float A_0;
	int minimum_segment_length;
	float sigma;
	float m_over_n;
	float Minimum_Slope;
	int target_nodes;
	int n_iterations;
	float fraction_dchi_for_variation;
	float vertical_interval;
	float horizontal_interval;
	float area_thin_frac;
	int target_skip;
                                                        
	file_info_in >> DEM_name >> Minimum_Slope >> pruning_threshold >> A_0 >> minimum_segment_length
        >> sigma >> m_over_n >> target_nodes >> n_iterations >> fraction_dchi_for_variation
				>> vertical_interval >> horizontal_interval >> area_thin_frac >> target_skip;


	cout << "Paramters of this run: " << endl
		   << "DEM ID: " << DEM_name << endl
	     << "pruning_threshold: " << pruning_threshold << endl
	     << "A_0: " << A_0 << endl
	     << "minimum_segment_length: " << minimum_segment_length << endl
	     << "sigma: " << sigma << endl
	     << "m over n: " << m_over_n << endl
	     << "target_nodes: " << target_nodes << endl
	     << "n_iterarions: " << n_iterations << endl
	     << "fraction_dchi_for_variation: " << fraction_dchi_for_variation << endl
	     << "target_skip is: " << target_skip << endl;


	file_info_in.close();
	
	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// load the filled DEM
  LSDRaster topo_test((path_name+DEM_name), DEM_flt_extension);
//	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), DEM_flt_extension);
  float MinSlope = 0.0005;
  LSDRaster filled_topo_test = topo_test.fill(MinSlope);
  string fill_name = "_fill";
  filled_topo_test.write_raster((path_name+DEM_name+fill_name), DEM_flt_extension);

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);
  
	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	//get the sources from raster to vector
  vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_name+sources_ext),DEM_flt_extension);
  
  // now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	LSDIndexRaster JIArray = ChanNetwork.JunctionIndexArray_to_LSDIndexRaster();
	LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
	
	//string SO_name = "_SO_test";
	//string JI_name = "_JI_test";
	
	//SOArray.write_raster((path_name+DEM_name+SO_name),DEM_flt_extension);
	//JIArray.write_raster((path_name+DEM_name+JI_name),DEM_flt_extension);
	  
    cout << "Reading in the cosmo points" << endl;
  // reading in cosmo points - getting upstream junctions and then extracting basins
  string string_filename;
	string filename = "_cosmo";
	string dot = ".";
	string extension = "txt";
  string_filename = DEM_name+filename+dot+extension;
  ifstream coords_list;
  coords_list.open(string_filename.c_str());
  
  float X_coord, Y_coord, ErosionRate, Error;
  vector<float> X_coords;
  vector<float> Y_coords;
  vector<float> ErosionRates_temp;
  vector<float> Errors_temp;
  while (coords_list >> X_coord >> Y_coord >> ErosionRate >> Error)
  {
     cout << "X Coord: " << X_coord << " Y Coord: " << Y_coord << " Erosion Rate: " << ErosionRate << " Error: " << Error << endl;
     X_coords.push_back(X_coord);
     Y_coords.push_back(Y_coord);
     ErosionRates_temp.push_back(ErosionRate);
     Errors_temp.push_back(Error);
  }
  
  float NoDataValue = topo_test.get_NoDataValue();
  float DataResolution = topo_test.get_DataResolution();
  cout << "Snapping to channel; getting upstream junctions" << endl;
  int threshold_SO = 2;
  int search_radius = 20;
  vector<int> junction_vector;
  vector<float> ErosionRates;
  vector<float> Errors;
  for (unsigned int i = 0; i < X_coords.size(); i++)
  {
    cout << flush << "Snapped = " << i+1 << " of " << X_coords.size() << "\r";
    int NICosmo = ChanNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(X_coords[i], Y_coords[i],
                                                                search_radius, threshold_SO, FlowInfo);

    if (NICosmo != NoDataValue)
    {
      int JunctionCosmo = ChanNetwork.find_upstream_junction_from_channel_nodeindex(NICosmo, FlowInfo);
      cout << "Junction of cosmo point: " << JunctionCosmo << endl;
      junction_vector.push_back(JunctionCosmo); 
      ErosionRates.push_back(ErosionRates_temp[i]);
      Errors.push_back(Errors_temp[i]); 
    }
  }   
   
  int NumberOfJunctions = junction_vector.size();
  
  string string_filename2;
  string filename2 = "_m_chi";
  string_filename2 = DEM_name+filename2+dot+extension;
  ofstream m_chi;
  m_chi.open(string_filename2.c_str());

  for(int i_junction = 0; i_junction < NumberOfJunctions; ++ i_junction)
  {
  	cout << "RUNNING CHI ANALYSIS from junction " << i_junction+1 << " of " << NumberOfJunctions << ": " << junction_vector[i_junction] << endl;
    int junction_number = junction_vector[i_junction];
    string jn_name = itoa(junction_number);
	  string uscore = "_";
	  jn_name = uscore+jn_name;
    // now get a junction and look for the longest channel upstream
  	cout << "creating main stem" << endl;
  	LSDIndexChannel main_stem = ChanNetwork.generate_longest_index_channel_in_basin(junction_number, FlowInfo, DistanceFromOutlet);
    cout << "got main stem channel, with n_nodes " << main_stem.get_n_nodes_in_channel() <<  endl;
  
//  	string Basin_name = "_basin";
//  	LSDIndexRaster BasinArray = ChanNetwork.extract_basin_from_junction(junction_number,junction_number,FlowInfo);
//   	BasinArray.write_raster((path_name+DEM_name+Basin_name+jn_name),DEM_flt_extension);
  
  	// now get the best fit m over n for all the tributaries
  	int organization_switch = 1;                                                  // Decide how I want to thin the dataset
  	int pruning_switch = 0;
  	LSDIndexChannelTree ChannelTree(FlowInfo, ChanNetwork, junction_number, organization_switch,
                                          DistanceFromOutlet, pruning_switch, pruning_threshold);
  
  	// print a file that can be ingested bt the chi fitting algorithm
   	string Chan_fname = "_ChanNet";
   	string Chan_ext = ".chan";
  	string Chan_for_chi_ingestion_fname = path_name+DEM_name+Chan_fname+jn_name+Chan_ext;
   	ChannelTree.print_LSDChannels_for_chi_network_ingestion(FlowInfo, filled_topo_test, DistanceFromOutlet, Chan_for_chi_ingestion_fname);
   	ChannelTree.convert_chan_file_for_ArcMap_ingestion(Chan_for_chi_ingestion_fname);
                               
  	// get a string with some paramter values
  	string sigma_str;
  	string skip_str;
  	string msl_str;
  	string tn_str;
  	string param_str;
  	sigma_str = static_cast<ostringstream*>( &(ostringstream() << sigma) )->str();
  	skip_str = static_cast<ostringstream*>( &(ostringstream() << target_skip) )->str();
  	msl_str = static_cast<ostringstream*>( &(ostringstream() << minimum_segment_length) )->str();
  	tn_str = static_cast<ostringstream*>( &(ostringstream() << target_nodes) )->str();
  
  	param_str = uscore+sigma_str+uscore+skip_str+uscore+msl_str+uscore+tn_str;
  
  // Create chi network automatically here, rather than reading in a file
  	// create the chi network
  	LSDChiNetwork ChiNetwork(Chan_for_chi_ingestion_fname);
  	LSDChiNetwork ChiNetwork_extended(Chan_for_chi_ingestion_fname);
  	ChiNetwork_extended.extend_tributaries_to_outlet();
  
  
  // 	//=-=-=-=-=-=-=-=-=-=-=-=-
  // 	// now caluculate the slope area data
  // 	//=-=-=-=-=-=-=-=-=-=-=-=-
  // 	string SA_vert_fname = "SAVert";
  // 	string SA_ext = ".SAdata";
  // 	ChiNetwork.slope_area_extraction_vertical_intervals(vertical_interval, area_thin_frac,
  // 	                                              (path_name+DEM_name+SA_vert_fname+jn_name+SA_ext));
  // 
  // 	string SA_horiz_fname = "SAHoriz";
  // 	ChiNetwork.slope_area_extraction_horizontal_intervals(horizontal_interval, area_thin_frac,
  // 	                                              (path_name+DEM_name+SA_horiz_fname+jn_name+SA_ext));
  
  	//=-=-=-=-=-=-=-=-=
  	///////////////////
  	// Now do the segment finding
  	///////////////////
  	//=-=-=-=-=-=-=-=-=
  
    cout << "Creating the chi profile" << endl;
    
    string fpt_ext = ".tree";
    string prefix_movn = static_cast<ostringstream*>( &(ostringstream() << m_over_n) )->str();
    
    // get the breaks of all the channels
    cout << "splitting extended chi network" << endl;
    ChiNetwork_extended.split_all_channels(A_0, m_over_n, n_iterations,
  	  					target_skip, target_nodes, minimum_segment_length, sigma);
    cout << "Here we go... starting Monte Carlo sampling..." << endl;
  	// monte carlo sample all channels
  	ChiNetwork_extended.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations,
  							target_skip, minimum_segment_length, sigma);
  							
    // get the mean m chi, standard deviation and standard error
  	vector< vector<float> > m_means = ChiNetwork_extended.get_m_means(); 	
  	int n_channels = m_means.size();
  	
  	float total_m = 0;
  	int n_observations = 0;
  	vector<float> m_chi_means;
    for (int i = 0; i < n_channels; i++)
  	{
      vector<float> m_mean = m_means[i];
        for (int j = 0; j < m_mean.size(); j++)
        {
          m_chi_means.push_back(m_mean[j]);
        }
    }
    float m_chi_mean = get_mean(m_chi_means);
    float m_chi_stdev = get_standard_deviation(m_chi_means, m_chi_mean);
    float m_chi_sterr = get_standard_error(m_chi_means, m_chi_stdev);

    cout << "Mean Mchi of basin: " << m_chi_mean << " Standard deviation: " << m_chi_stdev << " Standard error: " << m_chi_sterr << endl;
    
    int node = ChanNetwork.get_Node_of_Junction(junction_number);
    int ContributingPixels = FlowInfo.retrieve_contributing_pixels_of_node(node);
    float DrainageArea = ContributingPixels * DataResolution * DataResolution;
       
    m_chi << ErosionRates[i_junction] << " " << Errors[i_junction] << " " << m_chi_mean << " " << m_chi_stdev << " " << m_chi_sterr << " " << DrainageArea << endl;
    
    string fpt_mc = "_fullProfile" + param_str;
    string arc_name = "_fullProfile_for_Arc" + param_str;
    
    ChiNetwork_extended.print_channel_details_to_file_full_fitted((path_name+DEM_name+fpt_mc+jn_name+fpt_ext));
    ChiNetwork_extended.print_channel_details_to_file_full_fitted_for_ArcMap((path_name+DEM_name+arc_name+jn_name+fpt_ext));
  }
}
