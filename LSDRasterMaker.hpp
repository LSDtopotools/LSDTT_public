///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDRasterMaker.cpp
/// header for the RasterMaker object
/// The raster maker is a series of simple functions to make some rasters
/// with different properties. 
/// The initial use is mainly to make rasters for use in the raster model
/// for uplift and K
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1    01/09/2017
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterMaker_H
#define LSDRasterMaker_H

///@brief Create model objects to use LSDRaster methods on synthetic landscapes.
class LSDRasterMaker: public LSDRaster
{
  public:

    /// @brief Constructor. Create an LSDRasterMaker from an LSDRaster.
    /// @return LSDRasterMaker
    LSDRasterMaker()
    {
      create();
    }


    /// @brief Constructor. Create an LSDRasterMaker from an LSDRaster.
    /// @return LSDRasterMaker
    /// @param An_LSDRaster LSDRaster object.
    LSDRasterMaker(LSDRaster& An_LSDRaster)
    {
      create(An_LSDRaster);
    }

    /// @brief Constructor. Create an LSDRasterMaker from a file.
    /// Uses a filename and file extension
    /// @return LSDRasterMaker
    /// @param filename A String, the file to be loaded.
    /// @param extension A String, the file extension to be loaded.
    LSDRasterMaker(string filename, string extension)
    {
      create(filename, extension);
    }

    /// @brief Constructor. Create a blank raster nodel
    /// @return LSDRasterModel
    /// @param NCols Height of raster
    /// @param NRows Width of raster
    LSDRasterMaker(int NRows, int NCols)
    {
      create(NRows,NCols);
    }

    /// @brief This just returns the raster model object data as a raster
    /// @return A raster with the data from the LSDRasterModel
    /// @author SMM
    /// @date 01/09/2017
    LSDRaster return_as_raster();
    
    /// @brief This resizes the LSDRasterModel, resetting some flags in the process,
    /// as well as setting many of the Array2D data members to be empty arrays
    /// The raster data in the end is a random surface (determined by the noise
    /// data member)
    /// This overloaded version also changes the data resolution
    /// @param new_rows the new number of rows
    /// @param new_cols the new number of columns
    /// @param new_resolution the new data resolution
    /// @param new_value the new value of the raster
    /// @author SMM
    /// @date 01/09/2017
    void resize_and_reset( int new_rows, int new_cols, float new_resolution, float new_value );


    //void random_horizontal_strips(int minimum_strip_size, int maximum_strip_size, float minimum_value, float maximum_value);
    
    /// @brief This makes square blobs randomly in the DEM
    /// @param minimum_blob_size minimum blob size
    /// @param maximum_blob_size maximum blob size
    /// @param minimum_value minimum value on the raster
    /// @param maximum_value maximum value on the raster
    /// @parameter n_blobs the number of square blobs
    /// @author SMM
    /// @date 01/09/2017
    void random_square_blobs(int minimum_blob_size, int maximum_blob_size, float minimum_value, float maximum_value, int n_blobs);

    /// @brief This function adds some sine waves together to get patterns
    /// @param x_coefficients Coefficients of the sin waves in the x direction.
    /// @param y_coefficients Coefficients of the sin waves in the y direction.
    /// @author SMM
    /// @date 01/09/2017
    void sine_waves(vector<float> x_coefficients, vector<float> y_coefficients);

  protected:


    
  private:
    /// @brief Make a 100x100 raster
    void create();
    
    /// @brief Load a raster
    void create(string filename, string extension);
    
    /// @brief Make a rastermaker from another raster
    void create(LSDRaster& An_LSDRaster);
    
    /// @brief create a raster from NRows and NCols information
    void create(int NRows, int NCols);

};

#endif
