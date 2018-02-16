//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// get_elevation_profiles.cpp
// Takes in a baseline shapefile and gets the mean, max, and min elevation at each point
// along the shapefile.
// Writes the values to a csv file for visualisation.
//
// Notes: You need the shapefile to be a points file. 
// The point file needs to be in UTM coordinates. 
// It prints a raster and a csv
// I am not sure what the csv does. I get almost all nodata
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
#include <ctime>
#include "../../LSDRaster.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDParameterParser.hpp"

int main (int nNumberofArgs,char *argv[])
{
  //start the clock
  clock_t begin = clock();

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the swath tool!  					              ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Fiona J. Clubb												              ||" << endl;
    cout << "|| and Simon M Mudd                                    ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./swath_profiler.out /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_Swath.param" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }
  string path_name = argv[1];
  string f_name = argv[2];

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // set default parameters
  float_default_map["HalfWidth"] = 500;
  int_default_map["NormaliseToBaseline"] = 0;
  
  // Some parameters for the profiles
  float_default_map["swath_bin_spacing"] = 500;

  // set default string parameters
  string_default_map["Baseline_file"] = "Swath_points.shp";
  string_default_map["CHeads_format"] = "csv";
  
  // The analyses you want
  bool_default_map["print_DistanceToBaseline_raster"] = false;
  bool_default_map["print_DistanceAlongBaseline_raster"] = false;
  bool_default_map["print_BaselineValue_raster"] = false;
  bool_default_map["print_longitudinal_swath_profile"] = false;
  bool_default_map["print_baseline_csv"] = false;

  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);
  LSDPP.force_bil_extension();

  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();


  //string Swath_ext = "_swath_trans";
  //string Long_Swath_ext = "_swath_long";
  //string BV_ext = "_baseline_values";
  cout << "starting the test run... here we go!" << endl;

  // load the DEM
  cout << "Loading the DEM..." << endl;
  LSDRaster RasterTemplate((DATA_DIR+DEM_ID), raster_ext);

  cout << "\t loading baseline points" << endl;
  PointData BaselinePoints = LoadShapefile(DATA_DIR+this_string_map["Baseline_file"].c_str());

  // get the swath
  cout << "\t creating swath template" << endl;
  LSDSwath TestSwath(BaselinePoints, RasterTemplate, this_float_map["HalfWidth"]);


  if(  this_bool_map["print_DistanceToBaseline_raster"])
  {
    cout << "I am printing the distance to baseline raster." << endl;
    string of_name = DATA_DIR+DEM_ID+"_DistToBL";
    LSDRaster SR = TestSwath.get_raster_DistanceToBaselineArray(RasterTemplate);
    SR.write_raster(of_name,raster_ext);
  }
  if(  this_bool_map["print_DistanceAlongBaseline_raster"])
  {
    cout << "I am printing the distance along baseline raster." << endl;
    string of_name = DATA_DIR+DEM_ID+"_DistAlBL";
    LSDRaster SR = TestSwath.get_raster_DistanceAlongBaselineArray(RasterTemplate);
    SR.write_raster(of_name,raster_ext);
  }
  if(  this_bool_map["print_BaselineValue_raster"])
  {
    cout << "I am printing the baseline value raster." << endl;
    string of_name = DATA_DIR+DEM_ID+"BLV";
    LSDRaster SR = TestSwath.get_raster_BaselineValueArray(RasterTemplate);
    SR.write_raster(of_name,raster_ext);
  }



  if(this_bool_map["print_longitudinal_swath_profile"])
  {
    cout << "I am printing the longitudianl swath profile." << endl;
    // try to get the longitudunal swath profile
    vector<float> desired_percentiles;
    desired_percentiles.push_back(0);
    desired_percentiles.push_back(25);
    desired_percentiles.push_back(50);
    desired_percentiles.push_back(75);
    desired_percentiles.push_back(100);
    float BinWidth = 500;
    vector<float> mid_points;
    vector<float> mean_profile;
    vector<float> sd_profile;
    vector< vector<float> > output_percentile_profiles;
  
    // this gets the longitudinal swath
    TestSwath.get_longitudinal_swath_profile(RasterTemplate, desired_percentiles, this_float_map["swath_bin_spacing"],
                                            mid_points, mean_profile, sd_profile, 
                                            output_percentile_profiles,this_int_map["NormaliseToBaseline"]);
  
    // push back results to file for plotting
    ofstream output_file;
    string output_fname = "_swath_elevations.csv";
    output_file.open((DATA_DIR+DEM_ID+output_fname).c_str());
    output_file << "Distance,Mean,Min,Max" << endl;
    for (int i = 0; i < int(mid_points.size()); ++i)
    {
      output_file << mid_points[i] << "," << mean_profile[i] << "," << output_percentile_profiles[0][i] << "," << output_percentile_profiles[4][i] << endl;
    }
    output_file.close();
  }

  if(this_bool_map["print_longitudinal_swath_profile"])
  {
    string output_fname3 = "_baseline.csv";
    output_fname3 = DATA_DIR+DEM_ID+output_fname3;  
    TestSwath.print_baseline_to_csv(RasterTemplate,output_fname3);
  }

  // Done, check how long it took
  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
