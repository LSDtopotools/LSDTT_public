//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterSpectral.hpp
// header file for the LSDRasterSpectral object
// LSD stands for Land Surface Dynamics
// This object perform spectral analysis and is seperate from LSDRaster
// simply because it requires the FFTW package so this can be removed
// from compilation to retain portability
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
// Stuart Grieve, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		02/04/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterSpectral_H
#define LSDRasterSpectral_H

class LSDRasterSpectral: public LSDRaster
{
	public:
	//LSDRasterSpectral()										{ create(); }
	LSDRasterSpectral(string filename, string extension)	{ create(filename, extension); }
	LSDRasterSpectral(int nrows, int ncols, double xmin, double ymin,
	          double cellsize, double ndv, Array2D<double> data)
								{ create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }
	LSDRasterSpectral(LSDRaster& An_LSDRaster) 				{ create(An_LSDRaster); }

	// operator
	LSDRasterSpectral& operator=(const LSDRasterSpectral& LSDR);


	  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // FAST FOURIER TRANSFORM MODULE
    //------------------------------------------------------------------------------
    // Computes both the forward and inverse fast fourier transforms of a 2D
    // discrete dataset.
    //------------------------------------------------------------------------------
    // FOR FORWARD TRANSFORM:
    //    - InputArray = zeta_padded (padded DEM)
    //    - transform_direction = -1
    //    - OutputArray = 2D spectrum
    void dfftw2D_fwd(Array2D<double>& InputArray, Array2D<double>& OutputArrayReal, Array2D<double>& OutputArrayImaginary,
	                 int transform_direction, int Ly, int Lx);
  	//------------------------------------------------------------------------------
    // FOR INVERSE TRANSFORM:
    //    - InputArrays = Real and Imaginary components of 2D spectrum
    //    - transform_direction = 1
    //    - OutputArray = reconstructed DEM
    void dfftw2D_inv(Array2D<double>& InputArrayReal, Array2D<double>& InputArrayImaginary,
  	                 Array2D<double>& OutputArray, int transform_direction, int Ly, int Lx);
  	//------------------------------------------------------------------------------
    // DETREND DATA
    // Fit plane by least squares regression and use coefficients to determine
    // local slope ax + by + c = z
    // Returns detrended array, and an array with the trend plane
    void detrend2D(Array2D<double>& zeta, Array2D<double>& zeta_detrend,
  	               Array2D<double>& trend_plane, int nrows, int ncols, double ndv);
  	//------------------------------------------------------------------------------
    // HANN WINDOW MODULE
    // Use 2D elliptical Hann (raised cosine) window on data matrix, to reduce
    // spectral leakage and retain good frequency resolution.
    // Returns windowed data, the Hann window and also the summed square of the
    // weighting coefficients, WSS.
    void window_data_Hann2D(Array2D<double>& zeta_detrend, Array2D<double>& zeta_Hann2D,
  	                        Array2D<double>& Hann2D, double& WSS, int nrows, int ncols, int ndv);
  	//------------------------------------------------------------------------------
    // SHIFT ORIGIN OF SPECTRUM IN FOURIER DOMAIN
    // The output of the DFT algorithm must be rearranged to place the zero
    // wavenumber element near the center of the array.
    void shift_spectrum(Array2D<double>& spectrum_real,  Array2D<double>& spectrum_imaginary,
  	                    Array2D<double>& spectrum_real_shift, Array2D<double>& spectrum_imaginary_shift, int Ly, int Lx);
  	//------------------------------------------------------------------------------
    // DE-SHIFT ORIGIN OF SPECTRUM
    // Inverse process of above to return filtered spectrum to original format
    // required for the inverse fourier transform algorithm.
    void shift_spectrum_inv(Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                        Array2D<double>& FilteredSpectrumReal_deshift, Array2D<double>& FilteredSpectrumImaginary_deshift,
  	                        int Ly, int Lx);
  	//------------------------------------------------------------------------------
    // CALCULATE THE DFT PERIODOGRAM
    // Multiply fourier analysis output by complex conjugate and normalises.
    // Returns 2D PSD as only output
    Array2D<double> calculate_2D_PSD(Array2D<double>& spectrum_real_shift, Array2D<double>& spectrum_imaginary_shift,
  	                        int Lx, int Ly, double WSS);
  	//------------------------------------------------------------------------------
    // GET RADIAL POWER SPECTRUM
    // Collapse 2D PSD into a radial PSD
    void calculate_radial_PSD(Array2D<double>& P_DFT, vector<double>& RadialPSD_output, vector<double>& RadialFrequency_output,
  	                        int Lx, int Ly, double WSS, double dem_res);
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // COMPUTE DISCRETE FAST FOURIER TRANSFORM OF A REAL, 2-DIMENSIONAL DATASET.
    // Computes the 2D and radial power spectra of a 2D array.
    // Input arguements are:
    //                      i)  file identifier to prefix output files
    //                      ii) the width of the logarithmically spaced bins.
    //                          For topography, suggest this is 0.1 to start
  	void fftw2D_spectral_analysis(char* file_id, double LogBinWidth);


    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	// FUNCTIONS TO ADD WEIGHTS TO FOURIER SPECTRA (FOR USE IN SPECTRA FILTERS)
  	//------------------------------------------------------------------------------
    // BANDPASS FILTER
    // Filter array to band between frequency bands f1 and f2.  The bandpass filter
    // is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
    void bandpass_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                     Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double f1, double f2);
  	//------------------------------------------------------------------------------
    // LOWPASS FILTER
    // Filter array to retain frequencies below f1.  The filter edge is a radial
    // gaussian function with a SD of |f2-f1|/3.
    void lowpass_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                    Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double f1, double f2);
  	//------------------------------------------------------------------------------
    // HIGHPASS FILTER
    // Filter array to retain frequencies above f1.  The filter edge is a radial
    // gaussian function with a SD of |f2-f1|/3.
    void highpass_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                     Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double f1, double f2);
    //------------------------------------------------------------------------------
    // WIENER FILTER
    // The Wiener filter is a spectral filter that removes noise from an image or
    // DEM.  Weights in filter given by amplitudes of noise and signal:
    //        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
  	void wiener_filter(Array2D<double>& RawSpectrumReal, Array2D<double>& RawSpectrumImaginary,
  	                   Array2D<double>& FilteredSpectrumReal, Array2D<double>& FilteredSpectrumImaginary,
  	                      int Lx, int Ly, double dem_res, double WSS);



  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	// MAIN FUNCTIONS USING SPECTRAL FILTERS
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	// FAST FOURIER TRANSFORM FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
  	//
  	// Note that FLow <= FHigh
  	//
  	// There are three types of filters depending on the intentions of the user
  	//
  	// BANDPASS FILTER (FilterType = 1)
  	// Filter array to band between frequency bands f1 and f2.  The bandpass filter
  	// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
  	//
  	// LOWPASS FILTER (FilterType = 2)
  	// Filter array to retain frequencies below f1.  The filter edge is a radial
  	// gaussian function with a SD of |f2-f1|/3.  f1 is the frequency below which
  	// the filter starts to taper; f2 is the frequency at which the filter tapers to
  	// zero. If f1 = f2, the edge is effectively a step function.
  	// HIGHPASS FILTER (FilterType = 3)
  	//
  	// Filter array to retain frequencies above f2.  The filter edge is a radial
  	// gaussian function with a SD of |f2-f1|/3.  f2 is the frequency below which
  	// the filter starts to taper; f1 is the frequency at which the filter tapers to
  	// zero. If f1 = f2, the edge is effectively a step function.
  	//
  	// A second type of bandpass filter is possible by combining the highpass and
  	// lowpass filters.
  	//------------------------------------------------------------------------------
  	// David Milodowski, 18/12/2012
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	LSDRaster fftw2D_filter(int FilterType, double FLow, double FHigh);
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	// WIENER FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	//
  	// The Wiener filter is a spectral filter that removes noise from an image or
  	// DEM.  Essentially, it works on the principle that the observed spectrum
  	// contains the superposition of the real signal and an additional noise signal,
  	// which we want to remove.  If we know, or can make a reasonable guess at the
  	// noise, N(f), and signal, S(f), parts of the spectrum then we can remove the
  	// noise using the filter:
  	//
  	//        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
  	//
  	// For topography; at long wavelengths the topographic signal obeys an
  	// approximate power law relationship between amplitude and frequency,
  	// decreasing as the frequency increases (and wavelength decreases).  Noise
  	// typically dominates the high frequency part of the spectrum.  Thus at high
  	// frequencies the spectrum is dominated by noise, and the filter weight goes to
  	// zero.  In contrast, at low frequencies, the signal dominates and the filter
  	// weight goes to 1.
  	//
  	// The optimal wiener filter is described in more detail in Numerical Recipes,
  	// 13.3, p149.
  	//
  	// The exact structure of the noise is worth thinking about.  White noise, which
  	// is random, has equal power across all wavelengths.  In the instance of
  	// topography, noise can be created by a whole range of sources, from rock
  	// exposure, to pit and mound topography, to unfiltered vegetation etc.  It is
  	// likely that these sources will not produce purely white noise, but rather
  	// will show an element of structure.  This program makes two assumptions about
  	// the noise: i) it dominates the signal at high frequencies (close to the
  	// Nquist frequency) and ii) we can reasonably model this using a linear fit in
  	// log-log space - i.e. it obeys some form of power law function between
  	// frequency and amplitude.  Note that if the noise in the signal is really
  	// white noise, then the power law function for the noise would simply have an
  	// exponent of zero.  I prefer this formulation because it permits the
  	// characterisation of the noise model without assuming that the noise has a
  	// particular structure (white noise, pink noise etc.)
  	//
  	//------------------------------------------------------------------------------
  	// David Milodowski, 18/12/2012
  	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  	LSDRaster fftw2D_wiener();


	private:
	void create();
	void create(string filename, string extension);
	void create(int ncols, int nrows, double xmin, double ymin,
	            double cellsize, double ndv, Array2D<double> data);
	void create(LSDRaster& An_LSDRaster);
};

#endif
