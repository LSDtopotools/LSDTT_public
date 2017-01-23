cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")
project(CosmoChiGradient)
include_directories(/home/s1143956/Analysis/LSDTopotools/ /home/s1143956/Analysis/LSDTopotools/TNT/)
link_directories(/home/s1143956/Analysis/LSDTopotools/)
add_executable(../histogram_from_raster_FR_chi_specific.out ../histogram_from_raster_FR_chi_specific.cpp ../../../LSDMostLikelyPartitionsFinder.cpp ../../../LSDIndexRaster.cpp ../../../LSDRaster.cpp ../../../LSDFlowInfo.cpp ../../../LSDJunctionNetwork.cpp ../../../LSDChiNetwork.cpp ../../../LSDIndexChannel.cpp ../../../LSDChannel.cpp ../../../LSDIndexChannelTree.cpp ../../../LSDStatsTools.cpp ../../../LSDShapeTools.cpp ../../../LSDBasin.cpp ../../../LSDParticle.cpp  ../../../LSDCRNParameters.cpp ../../../LSDCosmoData.cpp ../../../LSDRasterInfo.cpp )
target_link_libraries(../histogram_from_raster_FR_chi_specific.out m stdc++ fftw3)

