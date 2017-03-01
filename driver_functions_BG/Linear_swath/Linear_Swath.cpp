//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Linear_swath.cpp
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to produce
// linear swath profiles between a list of two lat/long points
//
// B.G.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2016 Simon M. Mudd 2016
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../../LSDRaster.hpp"
#include "../../LSDRasterInfo.hpp"
#include "../../LSDIndexRaster.hpp"
#include "../../LSDJunctionNetwork.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDParameterParser.hpp"
#include "../../LSDSwathProfile.hpp"
#include "../../LSDShapeTools.hpp"
#include "../../LSDSpatialCSVReader.hpp"

// This funcion converts Lat/long coordinates into UTM
void get_x_and_y_from_latlong(vector<float>& UTME,vector<float>& UTMN, int& Number_of_lines, int UTM_zone_func)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;
  // number of samples
  int N_samples =  Number_of_lines;

  // set up some temporary vectors
  vector<float> this_UTMN(N_samples,0);
  vector<float> this_UTME(N_samples,0);

  double this_Northing;
  double this_Easting;

  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    //cout << "Converting point " << i << " to UTM." << endl;
    Converter.LLtoUTM_ForceZone(eId, UTMN[i], UTME[i],
                      this_Northing, this_Easting, UTM_zone_func);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
    //cout << "Easting: " << this_Easting << " and northing: " << this_Northing << endl;
  }

  UTME = this_UTME;
  UTMN = this_UTMN;
  // done
}

int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Linear Swath Tool!                   ||" << endl;
    cout << "|| This program has a number of options to produce     ||" << endl;
    cout << "|| Swath profile(s) between pair of lat/long points    ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Boris Gailleton, Fiona Clubb and Simon Mudd         ||" << endl;
    cout << "|| at the University of Edinburgh                      ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./Linear_swath.out /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_LS.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << " http://lsdtopotools.github.io/LSDTT_book/#_chi_analysis_part_3_getting_chi_gradients_for_the_entire_landscape" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the basin tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters

  map<string,string> string_default_map;
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;


  // set default  parameters

  float_default_map["HalfWidth"] = 1000; // HalfWidth of the profile in meters
  float HalfWidth = 50000; //TEMPORARY MANUAL SETTING OF THE HalfWidth (I need to sort the parameter reader bug)
  float_default_map["d_space"] = 250; // spacing between the sampling points
  float d_space = 100; //TEMPORARY MANUAL SETTING OF THE d_space (I need to sort the parameter reader bug)

  string_default_map["coordinate_csv_file"] = "example_coordinate.csv"; // csv file theat host Coordinates


  // Use the parameter parser to get the maps of the parameters required for the
  // analysis

  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,string> this_string_map = LSDPP.get_string_parameters();
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();


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

  cout << "Read filename is:" <<  DATA_DIR+DEM_ID << endl;

    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);

  string coordinate_csv_fname = DATA_DIR+"example_coordinate.csv";  // SAME REASON THAN ABOVE (I need to sort the parameter reader bug)



  // Loading the csv file and everything
  ifstream ifs;
  ifs.open(coordinate_csv_fname.c_str());

  // check if the parameter file exists
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: The coordinate points file \"" << coordinate_csv_fname
         << "\" doesn't exist" << endl;
  }
  else
  {
    // Host the coordinates of the "A" and "B" points
    vector<float> latitudeA;
    vector<float> longitudeA;
    vector<float> latitudeB;
    vector<float> longitudeB;

    // Other paraneters
    vector<string> this_string_vec;
    vector<string> empty_string_vec;
    string line_from_file;
    vector<string> IDs;
    int n_lines = 0; // number of data lines readed

    // now loop through the rest of the lines, getting the data.
    while( getline(ifs, line_from_file))
    {
      // reset the string vec
      this_string_vec = empty_string_vec;

      // create a stringstream
      stringstream ss(line_from_file);

      while( ss.good() )
      {
        string substr;
        getline( ss, substr, ',' );

        // remove the spaces
        substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

        // remove control characters
        substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

        // add the string to the string vec
        this_string_vec.push_back( substr );
      }

      if ( int(this_string_vec.size()) < 3)
      {
        cout << "Hey there, I am trying to load your basin outlet data but you seem not to have" << endl;
        cout << "enough columns in your file. I am ignoring a line" << endl;
      }
      else
      {
        // now convert the data

        string s = this_string_vec[0];
        IDs.push_back( s );
        latitudeA.push_back( atof( this_string_vec[0].c_str() ) );
        longitudeA.push_back( atof(this_string_vec[1].c_str() ) );
        latitudeB.push_back( atof( this_string_vec[2].c_str() ) );
        longitudeB.push_back( atof(this_string_vec[3].c_str() ) );
        n_lines++;
      }
    }
    ifs.close();

    // Hosts the points
    vector<float> points_before_processing;

    // host the UTm cordinates
    vector<float> XA = longitudeA;
    vector<float> XB = longitudeB;
    vector<float> YA = latitudeA;
    vector<float> YB = latitudeB;

    //conversion of lat/long into UTM
    int UTM_zone;
    bool is_North;
    topography_raster.get_UTM_information(UTM_zone, is_North);
    get_x_and_y_from_latlong(XA,YA,n_lines,UTM_zone);
    get_x_and_y_from_latlong(XB,YB,n_lines,UTM_zone);
    // done


    for(int i=0; i<n_lines;i++)
    {
      // get the swath
      cout << "\t creating swath template for row " << i << endl;
      points_before_processing.push_back(YA[i]);
      points_before_processing.push_back(XA[i]);
      points_before_processing.push_back(YB[i]);
      points_before_processing.push_back(XB[i]);
      cout << "I am creating the swath profile " << i<< " between: " << YA[i] << "/" << XA[i]<< " and " << YB[i]<< "/" << XB[i]<<endl;
      LSDSwath TestSwath(points_before_processing, topography_raster, HalfWidth, d_space);
      points_before_processing.clear();
      // preparing the output parameters
      ostringstream oss;
      ostringstream oss1;
      bool NormaliseToBaseline = false;
    	LSDRaster SwathRaster = TestSwath.get_raster_from_swath_profile(topography_raster, NormaliseToBaseline);
      oss << "_swath_raster_" << i;
    	//string swath_ext = ("_swath_raster"+i.str());
      // Writing the DEM
      string DEM_extension = "bil";
    	SwathRaster.write_raster((path_name+DEM_ID+oss.str()), DEM_extension);

      // get the raster values along the swath
    	vector <vector <float> > ElevationValues = TestSwath.get_RasterValues_along_swath(topography_raster, NormaliseToBaseline);

    	// push back results to file for plotting
    	ofstream output_file;
      oss1 << "_swath_elevations_" << i;
    	//string output_fname = "_swath_elevations_" + i.str();
      string CSV_ext = ".csv";
    	output_file.open((path_name+DEM_ID+oss1.str()+CSV_ext).c_str());
      // writing the csv
    	output_file << "Distance,Mean,Min,Max" << endl;
    	for (int j = 0; j < int(ElevationValues[0].size()); j++)
    	{
    		output_file << ElevationValues[0][j] << "," << ElevationValues[1][j] << "," << ElevationValues[2][j] << "," << ElevationValues[3][j] << endl;
    	}
    	output_file.close();

    }
  }
cout << "End of the program"<< endl ;
}
