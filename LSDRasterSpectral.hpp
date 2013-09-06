//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterSpectral
// Land Surface Dynamics StatsTools
//
// An object for manipulating rasters developed for the University of Edinburgh
//  Land Surface Dynamics group topographic toolbox. This is a derivative class
// from LSDRaster, for use specifically with spectral analysis.
//
// These tools have been seperated from the LSDRaster class mainly because
//  they require the FFTW library and are therefore less portable than
//  the standard LSDRaster object.
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** @file LSDRasterSpectral.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 0.0.1
@brief This object performs spectral analysis.
@details It is seperate from LSDRaster simply because it requires the FFTW package so this can be removed from compilation to retain portability.

@date 02/04/2013
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterSpectral_H
#define LSDRasterSpectral_H

/// @brief This object performs spectral analysis.
class LSDRasterSpectral: public LSDRaster
{
	public:
	//LSDRasterSpectral()										{ create(); }

  /// @brief Create an LSDRasterSpectral from a file.
  /// Uses a filename and file extension
  /// @return LSDRasterSpectral
  /// @param filename A String, the file to be loaded.
  /// @param extension A String, the file extension to be loaded.
  /// @author SMM
  /// @date 18/12/2012
  LSDRasterSpectral(string filename, string extension)	{ create(filename, extension); }
	/// @brief Create an LSDRasterSpectral from memory.
  /// @return LSDRasterSpectral
  /// @param nrows An integer of the number of rows.
  /// @param ncols An integer of the number of columns.
  /// @param xmin A double of the minimum X coordinate.
  /// @param ymin A double of the minimum Y coordinate.
  /// @param cellsize A double of the cellsize.
  /// @param ndv An integer of the no data value.
  /// @param data An Array2D of doubles in the shape nrows*ncols,
  ///containing the data to be written.
  /// @author SMM
  /// @date 18/12/2012
  LSDRasterSpectral(int nrows, int ncols, double xmin, double ymin,
	          double cellsize, double ndv, Array2D<double> data)
								{ create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }
	/// @brief Create an LSDRasterSpectral from an LSDRaster object.
  /// @param An_LSDRaster LSDRaster object.
  /// @return LSDRasterSpectral
  /// @author SMM
  /// @date 18/12/2012
  LSDRasterSpectral(LSDRaster& An_LSDRaster) 				{ create(An_LSDRaster); }

	/// Assignment operator.
	LSDRasterSpectral& operator=(const LSDRasterSpectral& LSDR);


	  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // FAST FOURIER TRANSFORM MODULE
    //------------------------------------------------------------------------------

    /// @brief Computes the forward fast fourier transform of a 2D discrete dataset.
    /// @param InputArray = zeta_padded (padded DEM).
    /// @param transform_direction = -1.
    /// @param OutputArrayReal = Real 2D spectrum.
    /// @param OutputArrayImaginary = Imaginary 2D spectrum.
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
   	/// @author David Milodowski
    /// @date 18/12/2012
    void dfftw2D_fwd(Array2D<double>& InputArray, Array2D<double>& OutputArrayReal, Array2D<double>& OutputArrayImaginary,
	                 int transform_direction, int Ly, int Lx);

    /// @brief Computes the inverse fast fourier transform of a 2D discrete dataset.
    /// @param InputArrayReal = Real component of 2D spectrum.
    /// @param InputArrayImaginary = Imaginary component of 2D spectrum.
    /// @param OutputArray = reconstructed DEM.
    /// @param transform_direction = -1.
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
  	/// @author David Milodowski
    /// @date 18/12/2012
    void dfftw2D_inv(Array2D<double>& InputArrayReal, Array2D<double>& InputArrayImaginary,
  	                 Array2D<double>& OutputArray, int transform_direction, int Ly, int Lx);

    /// @brief Detrend Data.
    ///
    /// @details Fit plane by least squares regression and use coefficients to determine local slope ax + by + c = z.
    /// @param zeta Input elevation data.
    /// @param zeta_detrend Output detrended elevation data.
    /// @param trend_plane Output array of trend plane.
    /// @param nrows Number of rows.
    /// @param ncols Number of columns.
    /// @param ndv No data value.
  	/// @author David Milodowski
    /// @date 18/12/2012
    void detrend2D(Array2D<double>& zeta, Array2D<double>& zeta_detrend,
  	               Array2D<double>& trend_plane, int nrows, int ncols, double ndv);

    /// @brief Hann Window Module.
    ///
    /// @details Use 2D elliptical Hann (raised cosine) window on data matrix, to reduce spectral leakage and retain good frequency resolution.
    /// @param zeta_detrend Detrended elevation data
    /// @param zeta_Hann2D Output windowed data.
    /// @param Hann2D Output Hann window.
    /// @param WSS Summed square of the weighting coefficients.
    /// @param nrows Number of rows.
    /// @param ncols Number of columns.
    /// @param ndv No data value.
  	/// @author David Milodowski
    /// @date 18/12/2012
    void window_data_Hann2D(Array2D<double>& zeta_detrend, Array2D<double>& zeta_Hann2D,
  	                        Array2D<double>& Hann2D, double& WSS, int nrows, int ncols, int ndv);

    /// @brief SHIFT ORIGIN OF SPECTRUM IN FOURIER DOMAIN.
    ///
    /// @details The output of the DFT algorithm must be rearranged to place the zero wavenumber element near the center of the array.
    /// @param spectrum_real
    /// @param spectrum_imaginary
    /// @param spectrum_real_shift
    /// @param spectrum_imaginary_shift
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
  	/// @author David Milodowski
    /// @date 18/12/2012
    void shift_spectrum(Array2D<double>& spectrum_real,  Array2D<double>& spectrum_imaginary,
  	                    Array2D<double>& spectrum_real_shift, Array2D<double>& spectrum_imaginary_shift, int Ly, int Lx);

    /// @brief DE-SHIFT ORIGIN OF SPECTRUM.
    ///
    /// @details Inverse process of shift_spectrum() to return filtered spectrum to original format required for the inverse fourier transform algorithm.
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param FilteredSpectrumReal_deshift
    /// @param FilteredSpectrumImaginary_deshift
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
  	/// @author David Milodowski
    /// @date 18/12/2012
    void shift_spectrum_inv(Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                        Array2D<double>& FilteredSpectrumReal_deshift, Array2D<double>& FilteredSpectrumImaginary_deshift,
  	                        int Ly, int Lx);

    /// @brief CALCULATE THE DFT PERIODOGRAM.
    ///
    /// @details Multiply fourier analysis output by complex conjugate and normalises.
    /// @param spectrum_real_shift
    /// @param spectrum_imaginary_shift
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
    /// @param WSS Summed square of the weighting coefficients.
    /// @return 2D array of DFT Periodogram.
   	/// @author David Milodowski
    /// @date 18/12/2012
    Array2D<double> calculate_2D_PSD(Array2D<double>& spectrum_real_shift, Array2D<double>& spectrum_imaginary_shift,
  	                        int Lx, int Ly, double WSS);

    /// @brief GET RADIAL POWER SPECTRUM.
    ///
    /// @details Collapse 2D PSD into a radial PSD.
    /// @param P_DFT
    /// @param RadialPSD_output
    /// @param RadialFrequency_output
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
    /// @param WSS Summed square of the weighting coefficients.
    /// @param dem_res DEM resolution.
   	/// @author David Milodowski
    /// @date 18/12/2012
    void calculate_radial_PSD(Array2D<double>& P_DFT, vector<double>& RadialPSD_output, vector<double>& RadialFrequency_output,
  	                        int Lx, int Ly, double WSS, double dem_res);

    /// @brief COMPUTE DISCRETE FAST FOURIER TRANSFORM OF A REAL, 2-DIMENSIONAL DATASET.
    ///
    /// @details Computes the 2D and radial power spectra of a 2D array.
    /// @param file_id File identifier to prefix output files
    /// @param LogBinWidth Wwidth of the logarithmically spaced bins. For topography, suggest this is 0.1 to start.
  	/// @author David Milodowski
    /// @date 18/12/2012
  	void fftw2D_spectral_analysis(char* file_id, double LogBinWidth);

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	// FUNCTIONS TO ADD WEIGHTS TO FOURIER SPECTRA (FOR USE IN SPECTRA FILTERS)
  	//------------------------------------------------------------------------------

    /// @brief BANDPASS FILTER.
    ///
    /// @details Filter array to band between frequency bands f1 and f2.  The bandpass filter
    /// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
    /// @param dem_res DEM resolution.
    /// @param f1
    /// @param f2
   	/// @author David Milodowski
    /// @date 18/12/2012
    void bandpass_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                     Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double f1, double f2);

    /// @brief LOWPASS FILTER.
    ///
    /// @details Filter array to retain frequencies below f1.  The filter edge is a radial gaussian function with a SD of |f2-f1|/3.
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
    /// @param dem_res DEM resolution.
    /// @param f1
    /// @param f2
   	/// @author David Milodowski
    /// @date 18/12/2012
    void lowpass_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                    Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double f1, double f2);

    /// @brief HIGHPASS FILTER.
    ///
    /// @details Filter array to retain frequencies above f1.  The filter edge is a radial gaussian function with a SD of |f2-f1|/3.
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
    /// @param dem_res DEM resolution.
    /// @param f1
    /// @param f2
  	/// @author David Milodowski
    /// @date 18/12/2012
    void highpass_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                     Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double f1, double f2);

    /// @brief WIENER FILTER.
    ///
    /// @details The Wiener filter is a spectral filter that removes noise from an image or DEM.
    ///  Weights in filter given by amplitudes of noise and signal: \n\n
    ///        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param Ly Array y dimension.
    /// @param Lx Array x dimension.
    /// @param dem_res DEM resolution.
    /// @param WSS Summed square of the weighting coefficients.
   	/// @author David Milodowski
    /// @date 18/12/2012
  	void wiener_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                   Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double WSS);

  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	// MAIN FUNCTIONS USING SPECTRAL FILTERS
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  	/// @brief FAST FOURIER TRANSFORM FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
  	///
  	/// Note that FLow <= FHigh
  	///
  	/// There are three types of filters depending on the intentions of the user
  	///
  	/// BANDPASS FILTER (FilterType = 1) \n
  	/// Filter array to band between frequency bands f1 and f2.  The bandpass filter
  	/// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
  	/// \n\n
  	/// LOWPASS FILTER (FilterType = 2)\n
  	/// Filter array to retain frequencies below f1.  The filter edge is a radial
  	/// gaussian function with a SD of |f2-f1|/3.  f1 is the frequency below which
  	/// the filter starts to taper; f2 is the frequency at which the filter tapers to
  	/// zero. If f1 = f2, the edge is effectively a step function.
  	/// \n\n
  	/// HIGHPASS FILTER (FilterType = 3) \n
  	/// Filter array to retain frequencies above f2.  The filter edge is a radial
  	/// gaussian function with a SD of |f2-f1|/3.  f2 is the frequency below which
  	/// the filter starts to taper; f1 is the frequency at which the filter tapers to
  	/// zero. If f1 = f2, the edge is effectively a step function.
  	/// \n\n
  	/// A second type of bandpass filter is possible by combining the highpass and
  	/// lowpass filters.
  	///
  	/// @param FilterType
  	/// @param FLow
  	/// @param FHigh
  	/// @author David Milodowski
    /// @date 18/12/2012
  	LSDRaster fftw2D_filter(int FilterType, double FLow, double FHigh);

  	/// @brief WIENER FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
  	///
  	/// The Wiener filter is a spectral filter that removes noise from an image or
  	/// DEM.  Essentially, it works on the principle that the observed spectrum
  	/// contains the superposition of the real signal and an additional noise signal,
  	/// which we want to remove.  If we know, or can make a reasonable guess at the
  	/// noise, N(f), and signal, S(f), parts of the spectrum then we can remove the
  	/// noise using the filter:
  	/// \n\n
  	///        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
  	/// \n\n
  	/// For topography; at long wavelengths the topographic signal obeys an
  	/// approximate power law relationship between amplitude and frequency,
  	/// decreasing as the frequency increases (and wavelength decreases).  Noise
  	/// typically dominates the high frequency part of the spectrum.  Thus at high
  	/// frequencies the spectrum is dominated by noise, and the filter weight goes to
  	/// zero.  In contrast, at low frequencies, the signal dominates and the filter
  	/// weight goes to 1.
  	/// \n\n
  	/// The optimal wiener filter is described in more detail in Numerical Recipes,
  	/// 13.3, p149.
  	/// \n\n
  	/// The exact structure of the noise is worth thinking about.  White noise, which
  	/// is random, has equal power across all wavelengths.  In the instance of
  	/// topography, noise can be created by a whole range of sources, from rock
  	/// exposure, to pit and mound topography, to unfiltered vegetation etc.  It is
  	/// likely that these sources will not produce purely white noise, but rather
  	/// will show an element of structure.  This program makes two assumptions about
  	/// the noise: i) it dominates the signal at high frequencies (close to the
  	/// Nquist frequency) and ii) we can reasonably model this using a linear fit in
  	/// log-log space - i.e. it obeys some form of power law function between
  	/// frequency and amplitude.  Note that if the noise in the signal is really
  	/// white noise, then the power law function for the noise would simply have an
  	/// exponent of zero.  I prefer this formulation because it permits the
  	/// characterisation of the noise model without assuming that the noise has a
  	/// particular structure (white noise, pink noise etc.)
  	/// @author David Milodowski
    /// @date 18/12/2012
  	LSDRaster fftw2D_wiener();

	private:
	void create();
	void create(string filename, string extension);
	void create(int ncols, int nrows, double xmin, double ymin,
	            double cellsize, double ndv, Array2D<double> data);
	void create(LSDRaster& An_LSDRaster);
};

#endif
