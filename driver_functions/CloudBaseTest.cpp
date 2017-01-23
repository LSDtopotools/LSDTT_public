//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// CloudBaseTest.cpp
// A test module for the LSDCloudBase module
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// David T. Milodowski
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDCloudBase.hpp"
// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
#include <liblas/liblas.hpp>

int main (int nNumberofArgs,char *argv[])
{
	if (nNumberofArgs != 5){
    cout << "FATAL ERROR: You must provide 4 input arguements\n\t - x coordinate of search point\n\t - y coordinate of search point\n\t - a search radius in metres\n\t - a .las file with the return heights\n\t" << endl;
    exit(EXIT_SUCCESS);
  }
  float searchPoint_x = atof(argv[1]);
  float searchPoint_y = atof(argv[2]);
  float searchRadius = atof(argv[3]);
  string las_file = argv[4];  
  cout << "starting the test run... here we go!" << endl;	
	int filter = 0;
	double res = 1;
// 	string las_file = "/geos/u64/shared/LSDdata/Idaho_PointCloud/id_659500_5038000.las";	
	LSDCloud TestCloud(las_file, filter,res);
	cout << TestCloud.get_NPts() << endl;
  
//   float searchPoint_x = 100;
//   float searchPoint_y = 100;
//   float searchRadius = 10;
  vector<float> pointValues;
  vector<int> pointIndex;
  vector<float> pointSquaredDistance;

 	TestCloud.RadiusSearch2D(searchPoint_x, searchPoint_y, searchRadius, pointValues, pointIndex, pointSquaredDistance);
  cout << pointValues.size() << endl;
	cout << "S.I.G." << endl;
}
