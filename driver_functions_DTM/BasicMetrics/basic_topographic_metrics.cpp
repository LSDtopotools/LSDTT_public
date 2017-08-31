//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// basic_topographic_metrics.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program extracts some basic topographic metrics from a DEM:
// - slope
// - aspect
// - topographic wetness
// 
// Routine written by D.T. Milodowski, leaning heavily on earlier scripts by
// S.M. Mudd 
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
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

int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=4)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name, the file name and the file extension" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];
  string f_ext = argv[3];

  string lchar = path_name.substr(path_name.length()-2,1);
  string slash = "/";
  cout << "lchar is " << lchar << " and slash is " << slash << endl;

  // make sure there is a slash at the end of the pathname      
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;  
    path_name = path_name+slash;
  } 
  cout << "The pathname is: " << path_name << endl;

  // make sure the extension is correct
  if(f_ext != "asc" && f_ext != "flt")
  {
    cout << "You did not choose a valid file option, options are asc and flt" << endl;
    cout << "you selected: " << f_ext << endl;
    exit(EXIT_SUCCESS);
  }

  string full_name = path_name+f_name;
  cout << "The full path is: " << full_name << endl;
  //string DEM_name;
  string fill_ext = "_fill";

  // get some file names
  string DEM_f_name = path_name+f_name;
  string DEM_ext = f_ext;

  // load the DEM
  LSDRaster topo(DEM_f_name, f_ext);

  // get the properties of the DEM
  //int NRows = topo_test.get_NRows();
  //int NCols = topo_test.get_NCols();
  //float NoDataValue = topo_test.get_NoDataValue();
  //float XMinimum = topo_test.get_XMinimum();
  //float YMinimum = topo_test.get_YMinimum();
  //float DataResolution = topo.get_DataResolution();

  // get the filled file
  float MinSlope = 0.00001;
  cout << "Filling the DEM" << endl;
  LSDRaster filled_topo = topo.fill(MinSlope);
  
  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";

  float surface_fitting_window_radius = 30;      // the radius of the fitting window in metres
  vector<LSDRaster> surface_fitting;
  LSDRaster slope;
  LSDRaster aspect;
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;             // this means you want the slope
  raster_selection[2] = 1;             // this means you want the aspect
  surface_fitting = filled_topo.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);

  // now print the derivative rasters to file.
  aspect = surface_fitting[6];
  aspect.write_raster(tan_curvature_name, DEM_ext);
  slope = surface_fitting[1];
  slope.write_raster(slope_name, DEM_ext);

  // now get flow area and topographic index
  LSDRaster D_inf_area = filled_topo.D_inf_units();
  LSDRaster TopoIndex = filled_topo.calculate_topographic_index(D_inf_area,slope);
  string area_name = "_DinfDA";
  string TI_name = "_topo_index";
  area_name = path_name+f_name+area_name;
  TI_name = path_name+f_name+TI_name;

  D_inf_area.write_raster(area_name,DEM_ext);
  TopoIndex.write_raster(TW_name,DEM_ext);

}
