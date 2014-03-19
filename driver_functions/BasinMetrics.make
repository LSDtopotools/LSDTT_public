# make with make -f BasinMetrics.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES= BasinMetrics.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDChiNetwork.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDJunctionNetwork.cpp ../LSDIndexChannel.cpp ../LSDChannel.cpp ../LSDIndexChannelTree.cpp ../LSDStatsTools.cpp ../LSDShapeTools.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE= BasinMetrics.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
