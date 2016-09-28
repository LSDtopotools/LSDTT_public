#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <string>
#include "../LSDRaster.hpp"
#include "../LSDStatsTools.hpp"
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

  path_name = FixPath(path_name);
  string full_file = path_name+f_name;

  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;
  cout << "I am going to check if that file actually exists" << endl;
  ifstream param_in(full_file.c_str());
  if( param_in.fail() )
  {
    cout << "\nFATAL ERROR: Your parameter file: \"" << full_file
         << "\" doesn't exist" << std::endl;
    cout << "You need to check the directory and filename! These files are case sensitive!" << endl;
    exit(EXIT_FAILURE);
  }
  cout << "The file exists, hopefully these are the first three lines:" << endl;
  string this_line;
  getline(param_in,this_line);
  cout << this_line << endl; 
  getline(param_in,this_line);
  cout << this_line << endl;   
  getline(param_in,this_line);
  cout << this_line << endl; 
  param_in.close();
  
  
  //cout << "Hi there. I am going to test maps" << endl;
  //map<string,string> method_map;
  //method_map["yo"] = "mama";
  //cout << "Map of yo: " << method_map["yo"] << endl;
  //string joe = "joe";
  //string momma = "momma";
  //pair<string,string> p(joe,momma);
  //method_map.insert(p);
  //cout << "mm joe: " << method_map[joe] <<" and string version: " << method_map["joe"] << endl;

  //string full_name = path_name+f_name;
  LSDAnalysisDriver ADriver(path_name,f_name);                                      
}
