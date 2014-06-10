CC=g++
CFLAGS=-c -Wall -O3 -pg -g
OFLAGS = -Wall -O3 -pg -g
LDFLAGS= -Wall
SOURCES= LH_Resolution.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDChiNetwork.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDJunctionNetwork.cpp ../LSDIndexChannel.cpp ../LSDChannel.cpp ../LSDIndexChannelTree.cpp ../LSDStatsTools.cpp ../LSDBasin.cpp ../LSDShapeTools.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE= LH_Resolution.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
