#include <iostream>
#include <sstream>
#include <cstdlib>
#include "../LSDRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDAnalysisDriver.hpp"
using namespace std;


int main(int argc, char *argv[])
{
  string path_name;
  string f_name;
  
  //Test for correct input arguments
  if (argc!=3)
  {
    cout << "FATAL ERROR: wrong number inputs. The program needs the path name and the file name" << endl;  
    exit(EXIT_SUCCESS);
  }
  else
  {
    path_name = argv[1];
    f_name = argv[2];
  }

  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

  //string full_name = path_name+f_name;
  LSDAnalysisDriver ADriver(path_name,f_name);                                      
}
