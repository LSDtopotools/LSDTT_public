# make with make -f drainage_density_step2_basins.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=drainage_density_step2_basins.cpp \
    ../../LSDMostLikelyPartitionsFinder.cpp \
    ../../LSDIndexRaster.cpp \
    ../../LSDRaster.cpp \
    ../../LSDFlowInfo.cpp \
    ../../LSDJunctionNetwork.cpp \
    ../../LSDIndexChannel.cpp \
    ../../LSDChannel.cpp \
    ../../LSDIndexChannelTree.cpp \
    ../../LSDStatsTools.cpp \
    ../../LSDShapeTools.cpp \
    ../../LSDBasin.cpp \
    ../../LSDParticle.cpp \
    ../../LSDCRNParameters.cpp
LIBS   = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=drainage_density_step2_basins.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
