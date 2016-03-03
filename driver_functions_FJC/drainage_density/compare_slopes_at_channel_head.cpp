//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// compare_slopes_at_channel_head.cpp
// make with compare_slopes_at_channel_head.make
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "../../LSDStatsTools.hpp"
#include "../../LSDRaster.hpp"
#include "../../LSDFlowInfo.hpp"
int main (int nNumberofArgs,char *argv[])
{
    //Test for correct input arguments
    if (nNumberofArgs!=3)
    {
        cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the driver file name" << endl;
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

    //read in the info from the parameter file  
    string DEM_ID;
    string Sources_name;
    string Slope_name;
    string DEM_extension = "bil";
    string temp;
    file_info_in >> temp >> DEM_ID
                 >> temp >> Sources_name
                 >> temp >> Slope_name;

    file_info_in.close();
    
    //load the DEM
    cout << "\t Filling the DEM" << endl;
    LSDRaster DEM(path_name+DEM_ID, DEM_extension);
    float MinSlope = 0.0001;
    LSDRaster FilledDEM = DEM.fill(MinSlope);
    FilledDEM.write_raster((path_name+DEM_ID+"_fill"), DEM_extension);

    
    //get a FlowInfo object
    cout << "\t Flow routing..." << endl;
    vector<string> BoundaryConditions(4, "No Flux");
	// get a flow info object
 	LSDFlowInfo FlowInfo(BoundaryConditions,FilledDEM);

    //get the channel heads
    cout << "\t Ingesting channel heads" << endl;
    vector<int> sources = FlowInfo.Ingest_Channel_Heads((path_name+DEM_ID+Sources_name), DEM_extension);
    
    //get downslope nodes from channel heads
    vector<int> DownslopeSources;
    vector<int> FinalHeads_Lower;
    float MoveDist = 2;
    FlowInfo.MoveChannelHeadDown(sources, MoveDist, DownslopeSources, FinalHeads_Lower);
    
    //get upslope nodes from channel heads
    vector<int> UpslopeSources;
    vector<int> FinalHeads_Upper;
    FlowInfo.MoveChannelHeadUp(sources, MoveDist, FilledDEM, UpslopeSources,  FinalHeads_Upper);
    
    cout << "CHeads: " << sources.size() << " Lower: " << DownslopeSources.size() << " Upper: " << UpslopeSources.size() << endl;
    
    //get the slope raster 
    float surface_fitting_window_radius = 6;      // the radius of the fitting window in metres
    vector<LSDRaster> surface_fitting;
    vector<int> raster_selection(8, 0);
    raster_selection[1] = 1;                      // this indicates you want the slope
  
    surface_fitting = FilledDEM.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
    LSDRaster Slope = surface_fitting[1]; 
    Slope.write_raster((path_name+DEM_ID+Slope_name), DEM_extension);

    //get the nodata value
    float NoDataValue = FilledDEM.get_NoDataValue();
    
    vector<float> PercentDiff;
    //get the difference between the upslope and downslope nodes
    for (int i =0; i < int(sources.size()); i++)
    {
        if (UpslopeSources[i] != NoDataValue && DownslopeSources[i] != NoDataValue)
        {
            //upslope sources
            int upper_row, upper_col;
            FlowInfo.retrieve_current_row_and_col(UpslopeSources[i], upper_row, upper_col);
            float UpperSlope = Slope.get_data_element(upper_row, upper_col);
            
            //downslope sources
            int lower_row, lower_col;
            FlowInfo.retrieve_current_row_and_col(DownslopeSources[i], lower_row, lower_col);
            float LowerSlope = Slope.get_data_element(lower_row, lower_col);
            
            //slope of channel head
            int row, col;
            FlowInfo.retrieve_current_row_and_col(sources[i], row, col);
            float CHSlope = Slope.get_data_element(row, col);
            
            //get the percentage difference compared to the channel head slope
            // (Upper slope - lower slope) / channel head slope
            float percent_diff = abs(UpperSlope - LowerSlope)/CHSlope;
            //cout << percent_diff << endl;
            if (percent_diff < 2) PercentDiff.push_back(percent_diff);    
        }
    }
    
    float mean_diff = get_mean(PercentDiff);
    cout << "Mean difference: " << mean_diff << endl;
}
