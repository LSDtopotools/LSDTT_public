cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(LSDTopoTools)

# Dependencies
# ============

# TODO: issue FFTW3LibraryDepends.cmake missing in most fftw3 installs
# https://github.com/FFTW/fftw3/issues/130
#find_package(fftw3 REQUIRED)
#message(STATUS "Found fttw3: ${FFTW3_INCLUDE_DIRS}/fftw3")


# Build
# =====
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")


include_directories(/TNT)
add_executable(LSDTT_BasicMetrics.exe Analysis_driver/LSDTT_BasicMetrics.cpp LSDRaster.cpp LSDIndexRaster.cpp LSDFlowInfo.cpp LSDStatsTools.cpp LSDShapeTools.cpp LSDJunctionNetwork.cpp LSDIndexChannel.cpp LSDChannel.cpp LSDMostLikelyPartitionsFinder.cpp LSDParameterParser.cpp LSDSpatialCSVReader.cpp LSDRasterSpectral.cpp LSDRasterInfo.cpp LSDIndexChannelTree.cpp LSDChiTools.cpp LSDChiNetwork.cpp LSDBasin.cpp LSDParticle.cpp LSDCRNParameters.cpp)
target_link_libraries(LSDTT_BasicMetrics.exe fftw3)
