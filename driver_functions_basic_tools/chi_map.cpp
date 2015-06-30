#include <iostream>
#include <sstream>
#include <cstdlib>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDFlowInfo.hpp"
using namespace std;


int main(int argc, char *argv[])
{

  if (argc != 4)
  {
    cout << "ERROR: you need to enter a data folder, \n" 
         << "a DEM name (without extension) and the extension name" << endl;
    exit(0);
  }

  string pathname = argv[1];
  string lchar = pathname.substr(pathname.length()-2,1);
  string slash = "/";
  cout << "lchar is " << lchar << " and slash is " << slash << endl;
      
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;  
    pathname = pathname+slash;
  } 
  cout << "The pathname is: " << pathname << endl;


  string dem_name = argv[2];
  cout << "The DEM filename is: " << dem_name << endl;
  string full_DEM_name = pathname+dem_name;
  cout << "The full path is: " << full_DEM_name << endl;
    
  string DEM_ext = argv[3];

  LSDRaster new_raster(full_DEM_name,DEM_ext);
  
  float min_slope = 0.0001;
  LSDRaster filled_raster = new_raster.fill(min_slope);
  
  // get the flow info 
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "Periodic";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "Periodic";
		
  cout << "Filled raster, getting flow info" << endl;
  LSDFlowInfo FI(boundary_conditions, filled_raster);
  cout << "Got flow info" << endl;
	
  float m_over_n = 0.5;
  float A_0 = 1000;
  float threshold_area = 200;
	
  //vector<int> bl_nodes(2);
  //bl_nodes[0] = 156860657;
  //bl_nodes[2] = 6177171;
  LSDRaster chi_map = FI.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0, 
                                                threshold_area);
  cout << "got chi map" << endl;                                              
  //LSDRaster chi_map = FI.get_upslope_chi_from_multiple_starting_nodes(bl_nodes,
  //                                              m_over_n, A_0, 
  //                                              threshold_area);
  
  string hs_name = "_HS";
  hs_name = full_DEM_name+hs_name;
  LSDRaster HS = filled_raster.hillshade(45,315,1); 
  HS.write_raster(hs_name,DEM_ext);


 
  float hundredtimesmn = m_over_n*100;
  int m_n = int(hundredtimesmn);
  int A_0int = int(A_0);
  int TAint = int(threshold_area);

  string m_nstr = itoa(m_n);
  m_nstr = "_mn"+m_nstr;

  string A0str = itoa(A_0int);
  A0str = "_A0"+A0str;

  string TAstr = itoa(TAint);
  TAstr = "_TA"+TAstr;

  string chistr = "_chi"+m_nstr+A0str+TAstr;
  
                     
  string cm_fname = full_DEM_name+chistr;                                             
  chi_map.write_raster(cm_fname, DEM_ext);                                             
}
