//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// extract_typical_metrics.cpp
// A driver function for use with the Land Surace Dynamics Topo Toolbox
// This program calculates channel heads using a chi method described in
// Clubb et al. (manuscript in prep)
//
// Developed by:
//  Simon M. Mudd
//
// Copyright (C) 2013 Simon M. Mudd 2013
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
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp"

int main (int nNumberofArgs,char *argv[])
{


	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;
		exit(EXIT_SUCCESS);
	}

	string path_name = argv[1];
	string f_name = argv[2];

	cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

	string full_name = path_name+f_name;


	//string DEM_name;
	string fill_ext = "_fill";

	// get some file names
	string DEM_f_name = path_name+f_name;
	string DEM_flt_extension = "flt";

	// load the DEM
	LSDRaster topo_test(DEM_f_name, DEM_flt_extension);

	// get the properties of the DEM
	//int NRows = topo_test.get_NRows();
	//int NCols = topo_test.get_NCols();
	//float NoDataValue = topo_test.get_NoDataValue();
	//float XMinimum = topo_test.get_XMinimum();
	//float YMinimum = topo_test.get_YMinimum();
	float DataResolution = topo_test.get_DataResolution();

	// get the filled file
	float MinSlope = 0.00001;
	cout << "Filling the DEM" << endl;
	LSDRaster filled_topo_test = topo_test.fill(MinSlope);
	filled_topo_test.write_raster((DEM_f_name+fill_ext),DEM_flt_extension);

	// set no flux boundary conditions
	vector<string> boundary_conditions(4);
	boundary_conditions[0] = "No";
	boundary_conditions[1] = "no flux";
	boundary_conditions[2] = "no flux";
	boundary_conditions[3] = "No flux";

	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions,filled_topo_test);

	// get the contributing pixels
	LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

	// calculate parameters for drainage extraction based on DEM resolution
	int threshold;
	float window_radius;

	ifstream param_file;
	string param_filename = "topo_metrics.param";
	param_filename = param_filename+path_name;

	param_file.open(param_filename.c_str());

	if( param_file.fail() )
	{

		if (DataResolution < 1.5)
		{
			window_radius = 7.0;
			threshold = 500;
		}
		else if (DataResolution >= 1.5 && DataResolution < 5.5)
		{
			window_radius = 15.1;
			threshold = 250;
		}
		else if (DataResolution >= 5.5 && DataResolution < 10.5)
		{
			window_radius = 20.1;
			threshold = 150;
		}
		else if (DataResolution >= 10.5 && DataResolution < 30.5)
		{
			window_radius = 60.1;
			threshold = 100;
		}
		else if (DataResolution >= 30.5)
		{
			window_radius = 180.1;
			threshold = 25;
		}
		cout << "No parameter files, setting parameters to defaults, "
		     << " resolution: " << DataResolution << " window_radius: " << window_radius << endl
		     << " and threshold pixels for channel: " << threshold << endl;

	}
	else
	{
		param_file >> threshold >> window_radius;
	}
	param_file.close();




	// get the sources: note: this is only to select basins!
	vector<int> sources;
	sources = FlowInfo.get_sources_index_threshold(ContributingPixels, threshold);

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

	string curvature_name = "_curv";
	string tan_curvature_name = "_tan_curv";
	string slope_name = "_slope";
	string pf_curvature_name = "_planform_curv";
	string prf_curvature_name = "_profile_curv";
	string hs_name = "_hs";
	curvature_name = path_name+f_name+curvature_name;
	tan_curvature_name = path_name+f_name+tan_curvature_name;
	pf_curvature_name = path_name+f_name+pf_curvature_name;
	prf_curvature_name = path_name+f_name+prf_curvature_name;
	slope_name = path_name+f_name+slope_name;
	hs_name = path_name+f_name+hs_name;


	// Calculate polyfit coefficients
	Array2D<float> a;
  	Array2D<float> b;
  	Array2D<float> c;
  	Array2D<float> d;
  	Array2D<float> e;
  	Array2D<float> f;
  	topo_test.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);

	// get various derived rasters
	LSDRaster hs_raster = topo_test.hillshade(45.0, 315.0, 1.0);
	LSDRaster slope_raster = topo_test.calculate_polyfit_slope(d, e);
	LSDRaster curvature = topo_test.calculate_polyfit_curvature(a ,b);
	LSDRaster tan_curvature = topo_test.calculate_polyfit_tangential_curvature(a ,b ,c ,d, e);
	LSDRaster pf_curvature = topo_test.calculate_polyfit_curvature(a ,b);
	LSDRaster prf_curvature = topo_test.calculate_polyfit_planform_curvature(a ,b ,c ,d, e);

	// Write rasters
	hs_raster.write_raster(hs_name, DEM_flt_extension);
	curvature.write_raster(curvature_name, DEM_flt_extension);
	tan_curvature.write_raster(tan_curvature_name, DEM_flt_extension);
	slope_raster.write_raster(slope_name, DEM_flt_extension);
	pf_curvature.write_raster(pf_curvature_name, DEM_flt_extension);
	prf_curvature.write_raster(prf_curvature_name, DEM_flt_extension);


	// now get flow area and distance
	LSDRaster D_inf_area = filled_topo_test.D_inf_units();
	LSDRaster Dist_from_outlet = FlowInfo.distance_from_outlet();

	string CP_name = "_CP";
	string area_name = "_DinfDA";
	string FD_name = "_FDist";
	CP_name = path_name+f_name+CP_name;
	area_name = path_name+f_name+area_name;
	FD_name = path_name+f_name+FD_name;

	D_inf_area.write_raster(area_name,DEM_flt_extension);
	Dist_from_outlet.write_raster(FD_name,DEM_flt_extension);
	ContributingPixels.write_raster(CP_name,DEM_flt_extension);


}
