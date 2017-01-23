#include <iostream>
#include <fstream>
#include "../LSDRaster.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=5)
  {
    cout << "========================================================"  << endl
         << "Welcome to the difference DEM tool. " << endl
         << "This program needs the path name," << endl
         << "the extension name, and the name of two DEMs" << endl
         << "For example, a command line could be: " << endl
         << "LINUX: "
         << "./DifferenceDEM.exe /home/myDEMS/ asc DEM_one DEM_two" << endl
         << "WINDOWS (cygwin):" << endl
         << "DifferenceDEM.exe c:/home/myDEMS/ asc DEM_one DEM_two" << endl << endl
         << "NOTE you need a slash at the end of the path name." << endl
         << "--------------------------------------"  << endl
         << "IMPORTANT: The DEMs need to have the same dimensions!" << endl
         << "========================================================"  << endl;
         
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string dem_ext = argv[2];
  string f_name1 = argv[3];
  string f_name2 = argv[4];
  
  if (dem_ext != "asc" && dem_ext != "bil" && dem_ext != "flt")
  {
    cout << "You have not chosen a valid DEM extenison." << endl;
    cout << "Valid extensions are bil, asc and flt" << endl;
    exit(EXIT_SUCCESS);
  }
  else
  {

    // get some file names
    string DEM_f_name1 = path_name+f_name1;
    string DEM_f_name2 = path_name+f_name2;

    // load the DEMs
    LSDRaster topo_test1(DEM_f_name1, dem_ext);
    LSDRaster topo_test2(DEM_f_name2, dem_ext);
    
    // now compare the DEMs
    float mean_difference = topo_test1.difference_rasters(topo_test2);
    
    cout << "The mean difference between rasters is: " << mean_difference << endl;
  }
  
  return 0;
}
    
    