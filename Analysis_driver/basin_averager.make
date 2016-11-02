# make with make -f basin_averager.make

CC=g++
CFLAGS=-c -Wall -O3 -g 
OFLAGS = -Wall -O3 -g 
LDFLAGS= -Wall
SOURCES=basin_averager.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDRasterInfo.cpp \
    ../LSDBasin.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDStatsTools.cpp \
    ../LSDJunctionNetwork.cpp \
    ../LSDIndexChannel.cpp \
    ../LSDChannel.cpp \
    ../LSDMostLikelyPartitionsFinder.cpp \
    ../LSDShapeTools.cpp \
    ../LSDAnalysisDriver.cpp \
    ../LSDCRNParameters.cpp \
    ../LSDParticle.cpp \
    ../LSDParameterParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=basin_averager.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
