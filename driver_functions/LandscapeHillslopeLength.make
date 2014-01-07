# make with make -f LandscapeHillslopeLength.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES= LandscapeHillslopeLength.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDChiNetwork.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDJunctionNetwork.cpp ../LSDIndexChannel.cpp ../LSDChannel.cpp ../LSDIndexChannelTree.cpp ../LSDStatsTools.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE= LandscapeHillslopeLength.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
