//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Modification of the channel writing driver to use the DrEICH algorithm to create chi 
// profiles rather than using an area threshold.
//
// NOTE: if using the same driver file as for the original channel writing driver 
//(chi_step2_write_channel_file_driver.cpp) then you MUST delete the line with the threshold
// area (line 3 in previous parameter file)
// NOTE#2: I have moved the declaration of the input raster type into the driver file - my DEMs
// are in .flt format rather than .bil.  The format of the DEM is now declared in line 2 of the
// driver file without the first "." e.g. "flt" or "bil"
//
// Updated to take a list of floodplain initiation points (X Y coordinates), snap them to the nearest
// channel, get upstream junction and then write channel file for each junction.
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb
// Simon M. Mudd
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

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	
  // make sure there is a slash on the end of the file
  string lchar = path_name.substr(path_name.length()-2,1);
  string slash = "/";      
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl; 
    path_name = path_name+slash;
  } 		
	
	string f_name = argv[2];

  cout << "\nYou are running the channel writing driver." << endl;
      // <<"IMPORTANT: this has been updated to load an ENVI DEM, whith extension .bil" << endl
      // <<"You can convert your DEM to this file format using gdal_translate, with -of ENVI" << endl
      // <<"See documentation at: http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/gdal_notes.html" << endl << endl;
       
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
	string sources_ext = "_CH";
	string fill_ext = "_fill";
	string raster_extension;
	file_info_in >> DEM_name >> raster_extension;
	float pruning_threshold;
	//int threshold;
	float A_0;
	int minimum_segment_length;
	float sigma;
	float start_movern;
	float d_movern;
	float Minimum_Slope;
	int n_movern;
	int target_nodes;
	int n_iterations;
	float fraction_dchi_for_variation;
	float vertical_interval;
	float horizontal_interval;
	float area_thin_frac;
	int target_skip;

	file_info_in >> Minimum_Slope >> pruning_threshold >> A_0 >> minimum_segment_length >> sigma >> start_movern
				>> d_movern >> n_movern >> target_nodes >> n_iterations >> fraction_dchi_for_variation
				>> vertical_interval >> horizontal_interval >> area_thin_frac >> target_skip;


	cout << "Paramters of this run: " << endl
	     //<< "junction number: " << junction_number << endl
	     << "pruning_threshold: " << pruning_threshold << endl
	     << "A_0: " << A_0 << endl
	     << "minimum_segment_length: " << minimum_segment_length << endl
	     << "sigma: " << sigma << endl
	     << "start_movern " << start_movern << endl
	     << "d_movern: " << d_movern << endl
	     << "n_movern: " << n_movern << endl
	     << "target_nodes: " << target_nodes << endl
	     << "n_iterarions: " << n_iterations << endl
	     << "fraction_dchi_for_variation: " << fraction_dchi_for_variation << endl
	     << "vertical interval: " << vertical_interval << endl
	     << "horizontal interval: " << horizontal_interval << endl
	     << "area thinning fraction for SA analysis: " << area_thin_frac << endl
	     << "target_skip is: " << target_skip << endl;


// 	string jn_name = itoa(junction_number);
// 	string uscore = "_";
// 	jn_name = uscore+jn_name;
	file_info_in.close();

	string DEM_f_name = path_name+DEM_name+fill_ext;

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";


	// load the filled DEM
	LSDRaster filled_topo_test((path_name+DEM_name+fill_ext), raster_extension);

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);

	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

	// get the sources
	vector<int> sources;
	sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_name+sources_ext),raster_extension);

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
	
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//                  Get the floodplain initiation points from file
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=	
  
  cout << "Reading in the floodplain initiation points" << endl;
  // reading in the floodplain points - getting downstream junctions and then extracting basins
  string string_filename;
	string filename = "_floodplain";
	string dot = ".";
	string extension = "txt";
  string_filename = DEM_name+filename+dot+extension;
  ifstream coords_list;
  coords_list.open(string_filename.c_str());
  
  float X_coord, Y_coord;
  vector<float> X_coords;
  vector<float> Y_coords;
  while (coords_list >> X_coord >> Y_coord)
  {
     cout << "X Coord: " << X_coord << " Y Coord: " << Y_coord << endl;
     X_coords.push_back(X_coord);
     Y_coords.push_back(Y_coord);
  }
  
  // Snap the floodplain points to the nearest channel
  cout << "Snapping FIPs to channel; getting junctions" << endl;
  int threshold_SO = 1;
  int search_radius = 10;
  vector<int> snapped_junction_indices;
  for (unsigned int i = 0; i < X_coords.size(); i++)
  {
    cout << flush << "Snapped = " << i+1 << " of " << X_coords.size() << "\r";
    int NICosmo = ChanNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(X_coords[i], Y_coords[i],
                                                                search_radius, threshold_SO, FlowInfo);

    int JunctionCosmo = ChanNetwork.find_upstream_junction_from_channel_nodeindex(NICosmo, FlowInfo);
    snapped_junction_indices.push_back(JunctionCosmo);  
  }  
  
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//                     Write channel file for each junction
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=	  
            
  for (int i=0; i < snapped_junction_indices.size(); i++)
  {
    int junction_number = snapped_junction_indices[i];
    cout << "Junction of FIP: " << junction_number << endl;
    
    string jn_name = itoa(junction_number);
   	string uscore = "_";
  	jn_name = uscore+jn_name;
  
  	// now get a junction and look for the longest channel upstream
  	cout << "creating main stem" << endl;
  	LSDIndexChannel main_stem = ChanNetwork.generate_longest_index_channel_in_basin(junction_number, FlowInfo, DistanceFromOutlet);
          cout << "got main stem channel, with n_nodes " << main_stem.get_n_nodes_in_channel() <<  endl;
  
  	string Basin_name = "_basin_dreich";
  	LSDIndexRaster BasinArray = ChanNetwork.extract_basin_from_junction(junction_number,junction_number,FlowInfo);
  	BasinArray.write_raster((path_name+DEM_name+Basin_name+jn_name),raster_extension);
  
  	// now get the best fit m over n for all the tributaries
  	int organization_switch = 1;
  	int pruning_switch = 1;
  	LSDIndexChannelTree ChannelTree(FlowInfo, ChanNetwork, junction_number, organization_switch,
                                          DistanceFromOutlet, pruning_switch, pruning_threshold);
  
  	// print a file that can be ingested bt the chi fitting algorithm
  	string Chan_fname = "_ChanNet_dreich";
  	string Chan_ext = ".chan";
  	string Chan_for_chi_ingestion_fname = path_name+DEM_name+Chan_fname+jn_name+Chan_ext;
  	ChannelTree.print_LSDChannels_for_chi_network_ingestion(FlowInfo,
                               filled_topo_test, DistanceFromOutlet, Chan_for_chi_ingestion_fname);
  	ChannelTree.convert_chan_file_for_ArcMap_ingestion(Chan_for_chi_ingestion_fname);
  }
}
