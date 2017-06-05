# make with make -f knickpoint_v2.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=knickpoint_v2.cpp \
             ../../LSDMostLikelyPartitionsFinder.cpp \
             ../../LSDIndexRaster.cpp \
             ../../LSDRaster.cpp \
             ../../LSDRasterInfo.cpp \
             ../../LSDFlowInfo.cpp \
             ../../LSDJunctionNetwork.cpp \
             ../../LSDIndexChannel.cpp \
             ../../LSDChannel.cpp \
             ../../LSDIndexChannelTree.cpp \
             ../../LSDStatsTools.cpp \
             ../../LSDShapeTools.cpp \
             ../../LSDChiNetwork.cpp \
             ../../LSDBasin.cpp \
             ../../LSDParticle.cpp \
             ../../LSDChiTools.cpp \
             ../../LSDParameterParser.cpp \
             ../../LSDSpatialCSVReader.cpp \
             ../../LSDCRNParameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=knickpoint_v2.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
