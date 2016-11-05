// LSDNetCDFTools.cpp

/// Note: you need to compile with -lnetcdf and -lnetcdf_cxx4 flags
/// You will need the need the netCDF-4 C++ libraries installed
/// Most Linux package managers will install this for you using
/// yum, dnf, or apt-get.

/// A tool for reading netCDF-4 files into LSDTopoTools
///
/// NetCDF (The network common data format is a file format for storing
/// multivariate climate, atmospheric, etc. data. (Although it can be used to store
/// any type of gridded data. You could use it to store states of topography through
/// time, for example.)

/// Full documentation of the netCDF C++ API can be found at:
/// http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-cxx
///
/// @author DAV
/// @date 2016-11-05

#include <iostream>
#include <netcdf>
#include "TNT/tnt.h"

using namespace netCDF;
using namespace netCDF::exceptions;

// Return this in event of a problem.
static const int NC_ERR = 2;

// Originally had it as a class but not much point at this stage
//class LSDNetCDFTools
//{
//public:

/// @brief Reads in a netCDF field and writes it into
/// a TNT array 2D.
/// @return An error code, (0 for success, NC_ERR for failure)
/// @details Template allows to specify different TNT array data types
/// e.g. double, float, int etc.
template<typename T>
int read_tnt_array2d(TNT::Array2D<T>*const tnt_array);  
// The pointer to the tnt array is constant, but not the actual array itself
// i.e. we can change the array thtough the pointer, but ensure the pointer
// remains constant and cannot be changed to point to something else, which would be 
// a bad idea...

//};

template<typename T>
int read_tnt_array2d(TNT::Array2D<T>*const tnt_array)
{
   int NX = tnt_array->dim1();
   int NY = tnt_array->dim2();

   std::cout << "Rows: " << NX << ", " << "Cols: " << NY << std::endl;

   try
   {
   // This is a C-style array.
   // We tend not to use these in LSDTopoTools
   // but typically used in netCDF.
   // int dataIn[NX][NY];

   // This is a TNT (Template Numerical Toolkit)
   // 2D array (allocated at runtime)
   //tnt_array(NX,NY,0);

   // Open the file for read access
   NcFile dataFile("simple_xy.nc", NcFile::read);

   // Retrieve the variable named "data"
   NcVar data = dataFile.getVar("data");
   if(data.isNull()) return NC_ERR;

   // Ingest the data into the array
   // We pass a reference to the first element of
   // the TNT array [0][0]. getVar writes in the values.
   data.getVar(tnt_array[0][0]);

   // Note: with a C-style array you would do this:
   // data.getVar(dataIn)
   // Since C arrays are already pointers to memory.

   // Check the values.  // DV - this does not work
   // Something strange to do with pointers...
//   for (int i = 0; i < NX; i++)
//   {
//      for (int j = 0; j < NY; j++)
//      {
//         // Check that we haven't gone out of bounds...
//         if (*tnt_array[i][j] != i * NY + j) return NC_ERR;
////         else
////         {
////           std::cout << *tnt_array[i][j] << " ";
////         }

//      }
//      //std::cout << std::endl;
//   }
   
   // The netCDF file is automatically closed by the NcFile destructor
   std::cout << "*** SUCCESS reading example file simple_xy.nc!" << std::endl;

   return 0;
   }
   catch(NcException& e)
   {
     e.what();
     std::cout<<"FAILURE*************************************"<< std::endl;
     return NC_ERR;
   }
}

// Test.
int main()
{
  TNT::Array2D<int> myArray2D(6,12, 0);

  read_tnt_array2d(&myArray2D);

  for (int i = 0; i<6; i++)
  {
    for (int j = 0; j<12; j++)
    {
      //myArray2D[i][j] += 7;
      std::cout << myArray2D[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
