# make with make -f extract_typical_metrics.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=extract_typical_metrics.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDChiNetwork.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDChannelNetwork.cpp ../LSDIndexChannel.cpp ../LSDChannel.cpp ../LSDIndexChannelTree.cpp ../LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=extract_typical_metrics.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
