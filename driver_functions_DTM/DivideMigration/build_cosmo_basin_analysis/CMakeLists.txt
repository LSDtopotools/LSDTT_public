cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")
project(CosmoChiGradient)
include_directories(/home/s1143956/Analysis/LSDTopotools/ /home/s1143956/Analysis/LSDTopotools/TNT/)
link_directories(/home/s1143956/Analysis/LSDTopotools/)
add_executable(../cosmo_basin_chi_gradient.out ../cosmo_basin_chi_gradient.cpp ../../../LSDMostLikelyPartitionsFinder.cpp ../../../LSDIndexRaster.cpp ../../../LSDRaster.cpp ../../../LSDFlowInfo.cpp ../../../LSDJunctionNetwork.cpp ../../../LSDChiNetwork.cpp ../../../LSDIndexChannel.cpp ../../../LSDChannel.cpp ../../../LSDIndexChannelTree.cpp ../../../LSDStatsTools.cpp ../../../LSDShapeTools.cpp ../../../LSDBasin.cpp ../../../LSDParticle.cpp  ../../../LSDCRNParameters.cpp ../../../LSDCosmoData.cpp ../../../LSDRasterInfo.cpp )
target_link_libraries(../cosmo_basin_chi_gradient.out m stdc++ fftw3)

