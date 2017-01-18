//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Chile_test.cpp
//
// This program is used for testing the LSDRaster object
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Simon M. Mudd
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDFlowInfo.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDIndexChannelTree.hpp"
#include "../../LSDChiNetwork.hpp"
#include "../../TNT/tnt.h"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the driver file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string driver_file = argv[2];
  
	cout << "The path is: " << path_name << " and the filename is: " << driver_file << endl;

	string full_name = path_name+driver_file;

	ifstream file_info_in;
	file_info_in.open(full_name.c_str());
	if( file_info_in.fail() )
	{
		cout << "\nFATAL ERROR: the header file \"" << full_name
		     << "\" doesn't exist" << endl;
		exit(EXIT_FAILURE);
	}

	string DEM_name,DEM_extension,ListOfChiProfiles;
	file_info_in >> DEM_name;
	file_info_in >> DEM_extension;
	file_info_in >> ListOfChiProfiles;
  file_info_in.close();

	cout << "Run details: " << endl
	     << "DEM: " << DEM_name << endl 
	     << "Extension: " << DEM_extension << endl
	     << "ChannelTree list: " << ListOfChiProfiles << endl;

  LSDRaster DEM(DEM_name,DEM_extension);
  int NCols = DEM.get_NCols();
  int NRows = DEM.get_NRows();
  float NoDataValue = DEM.get_NoDataValue();
  float XMinimum = DEM.get_XMinimum();
  float YMinimum = DEM.get_YMinimum();
  float DataResolution = DEM.get_DataResolution();
  
  // Generated arrays for outputs.
  Array2D<float> chi_array(NRows,NCols,NoDataValue);
  Array2D<float> chi_gradient_array(NRows,NCols,NoDataValue);
  Array2D<float> chi_gradient_sd_array(NRows,NCols,NoDataValue);
  cout << "Reading list of channel tree files" << endl;
  // Read list of ChannelTree_files
  vector<string> ChiProfile_list;
  string line;
  ifstream FileList(ListOfChiProfiles.c_str());
  if (FileList.is_open())
  {
    while (getline (FileList,line))
    {
      ChiProfile_list.push_back(line);
    }
    FileList.close();
  }
  cout << "Loop through list of ChannelTree files, loading data into the above arrays" << endl;
  // Loop through list of ChannelTree files, loading data into the above arrays
  for(int i_file = 0; i_file < int(ChiProfile_list.size()); ++i_file)
  {
    // open profile
    string ChiProfile = ChiProfile_list[i_file];
    ifstream ChiProfile_in;
  	ChiProfile_in.open(ChiProfile.c_str());

    if( ChiProfile_in.fail() )
  	{
  		cout << "\nFATAL ERROR: the file \"" << ChiProfile << "\" doesn't exist" << endl;
  		exit(EXIT_FAILURE);
  	}  	
    
    // read in data
		float A_0, m_ov_n;
    ChiProfile_in >> A_0 >> m_ov_n;
		
    int node, row, col, channel_number, receiver_channel, node_on_receiver_channel, n_data_points_uic;
		float elevation,flow_distance,drainage_area,chi,
            m_mean,m_standard_deviation,m_standard_error,
            b_mean,b_standard_deviation,b_standard_error,
            DW_mean,DW_standard_deviation,DW_standard_error,
            fitted_elev_mean,fitted_elev_standard_deviation,fitted_elev_standard_error;
    
    while (ChiProfile_in.eof()==false)  
    {
      ChiProfile_in >> channel_number
                    >> receiver_channel >> node_on_receiver_channel
                    >> node >> row >> col
                    >> flow_distance >> chi >> elevation >> drainage_area 
                    >> n_data_points_uic
                    >> m_mean >> m_standard_deviation >> m_standard_error
						        >> b_mean >> b_standard_deviation >> b_standard_error
						        >> DW_mean >> DW_standard_deviation >> DW_standard_error
                    >> fitted_elev_mean >> fitted_elev_standard_deviation >> fitted_elev_standard_error;
      
      chi_array[row][col] = chi;
      chi_gradient_array[row][col] = m_mean;
      chi_gradient_sd_array[row][col] = m_standard_deviation;
    }
    ChiProfile_in.close();
  }
  
  // Write output rasters  
  LSDRaster ChiRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, chi_array);
  ChiRaster.write_raster(DEM_name+"_chi",DEM_extension);
  LSDRaster ChiMValueRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, chi_gradient_array);
  ChiMValueRaster.write_raster(DEM_name+"_chi_m_values",DEM_extension);
  LSDRaster ChiMSDRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, chi_gradient_sd_array);
  ChiMSDRaster.write_raster(DEM_name+"_chi_m_sd",DEM_extension);

}
