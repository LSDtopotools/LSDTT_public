# make with make -f LSDTT_BasicMetrics.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=LSDTT_BasicMetrics.cpp \
         ../LSDIndexRaster.cpp \
         ../LSDRaster.cpp \
         ../LSDChiTools.cpp \
         ../LSDChiNetwork.cpp \
         ../LSDBasin.cpp \
         ../LSDParticle.cpp \
         ../LSDCRNParameters.cpp \
         ../LSDFlowInfo.cpp \
         ../LSDIndexChannel.cpp \
         ../LSDStatsTools.cpp \
         ../LSDJunctionNetwork.cpp \
         ../LSDRasterInfo.cpp \
         ../LSDParameterParser.cpp \
         ../LSDSpatialCSVReader.cpp \
         ../LSDChannel.cpp \
         ../LSDMostLikelyPartitionsFinder.cpp \
         ../LSDShapeTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=LSDTT_BasicMetrics.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
