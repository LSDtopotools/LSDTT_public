# make with make -f test_junction_angles.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=test_junction_angles.cpp \
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
EXECUTABLE=test_junction_angles.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
