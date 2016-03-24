cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")
project(ChiMap)
include_directories(../../../ ../../../TNT/)
link_directories(../../../)
add_executable(../map_chi_gradient.out ../map_chi_gradient.cpp ../../../LSDMostLikelyPartitionsFinder.cpp ../../../LSDIndexRaster.cpp ../../../LSDRaster.cpp ../../../LSDFlowInfo.cpp ../../../LSDJunctionNetwork.cpp ../../../LSDChiNetwork.cpp ../../../LSDIndexChannel.cpp ../../../LSDChannel.cpp ../../../LSDIndexChannelTree.cpp ../../../LSDStatsTools.cpp ../../../LSDShapeTools.cpp ../../../LSDBasin.cpp ../../../LSDParticle.cpp  ../../../LSDCRNParameters.cpp ../../../LSDCosmoData.cpp ../../../LSDRasterInfo.cpp )
target_link_libraries(../map_chi_gradient.out m stdc++ fftw3)

