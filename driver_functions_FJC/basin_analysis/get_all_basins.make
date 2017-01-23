# make with make -f get_all_basins.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=get_all_basins.cpp ../../LSDIndexRaster.cpp ../../LSDRaster.cpp ../../LSDFlowInfo.cpp ../../LSDIndexChannel.cpp ../../LSDStatsTools.cpp ../../LSDJunctionNetwork.cpp ../../LSDChannel.cpp ../../LSDMostLikelyPartitionsFinder.cpp ../../LSDBasin.cpp ../../LSDParticle.cpp ../../LSDCRNParameters.cpp ../../LSDShapeTools.cpp

OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=get_all_basins.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
