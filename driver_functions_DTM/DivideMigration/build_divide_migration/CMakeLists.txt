cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")
project(DivideMigration)
include_directories(/home/s1143956/Analysis/LSDTopotools/ /home/s1143956/Analysis/LSDTopotools/TNT/)
link_directories(/home/s1143956/Analysis/LSDTopotools/)
add_executable(../map_divide_migration.out ../map_divide_migration.cpp ../../../LSDMostLikelyPartitionsFinder.cpp ../../../LSDIndexRaster.cpp ../../../LSDRaster.cpp ../../../LSDFlowInfo.cpp ../../../LSDJunctionNetwork.cpp ../../../LSDChiNetwork.cpp ../../../LSDIndexChannel.cpp ../../../LSDChannel.cpp ../../../LSDIndexChannelTree.cpp ../../../LSDStatsTools.cpp ../../../LSDShapeTools.cpp ../../../LSDBasin.cpp ../../../LSDParticle.cpp  ../../../LSDCRNParameters.cpp ../../../LSDCosmoData.cpp ../../../LSDRasterInfo.cpp )
target_link_libraries(../map_divide_migration.out m stdc++ fftw3)

