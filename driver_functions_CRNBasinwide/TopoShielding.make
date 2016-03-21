CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3

SOURCES= TopographicShielding.cpp \
        ../LSDMostLikelyPartitionsFinder.cpp \
        ../LSDChiNetwork.cpp \
        ../LSDIndexRaster.cpp \
        ../LSDRaster.cpp \
        ../LSDFlowInfo.cpp \
        ../LSDJunctionNetwork.cpp \
        ../LSDIndexChannel.cpp \
        ../LSDChannel.cpp \
        ../LSDIndexChannelTree.cpp \
        ../LSDStatsTools.cpp \
        ../LSDBasin.cpp \
        ../LSDShapeTools.cpp \
        ../LSDParticle.cpp \
        ../LSDCRNParameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE= TopoShielding.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
