# make with make -f DEM_preprocessing.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=DEM_preprocessing.cpp \
             ../LSDMostLikelyPartitionsFinder.cpp \
             ../LSDIndexRaster.cpp \
             ../LSDRaster.cpp \
             ../LSDRasterInfo.cpp \
             ../LSDFlowInfo.cpp \
             ../LSDJunctionNetwork.cpp \
             ../LSDIndexChannel.cpp \
             ../LSDChannel.cpp \
             ../LSDIndexChannelTree.cpp \
             ../LSDStatsTools.cpp \
             ../LSDShapeTools.cpp \
             ../LSDChiNetwork.cpp \
             ../LSDBasin.cpp \
             ../LSDParticle.cpp \
             ../LSDChiTools.cpp \
             ../LSDParameterParser.cpp \
             ../LSDCRNParameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=DEM_preprocessing.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
