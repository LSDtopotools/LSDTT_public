# make with make -f geology_mapper.make

CC=g++
CFLAGS=-c -Wall -O3 -g 
OFLAGS = -Wall -O3 -g 
LDFLAGS= -Wall
SOURCES=geology_mapper.cpp \
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
    ../LSDSpatialCSVReader.cpp \
    ../LSDShapeTools.cpp \
    ../LSDAnalysisDriver.cpp \
    ../LSDCRNParameters.cpp \
    ../LSDParticle.cpp \
    ../LSDParameterParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=geology_mapper.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
