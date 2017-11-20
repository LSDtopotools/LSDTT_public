CC=g++
CFLAGS=-c -Wall -O3 -pg -g
OFLAGS = -Wall -O3 -pg -g
LDFLAGS= -Wall
SOURCES= hilltop_mapping_tool.cpp \
        ../../LSDMostLikelyPartitionsFinder.cpp \
        ../../LSDIndexRaster.cpp \
        ../../LSDRaster.cpp \
        ../../LSDFlowInfo.cpp \
        ../../LSDJunctionNetwork.cpp \
        ../../LSDIndexChannel.cpp \
        ../../LSDChannel.cpp \
        ../../LSDStatsTools.cpp \
        ../../LSDBasin.cpp \
        ../../LSDShapeTools.cpp \
        ../../LSDParticle.cpp \
        ../../LSDCRNParameters.cpp \
        ../../LSDParameterParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE= hilltop_mapping_tool.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
