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

    /// @brief Class destructor
    ~LSDRasterMaker( void );

    /// @brief operator assignment
    LSDRasterMaker& operator=(const LSDRasterMaker& LSDR);

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
