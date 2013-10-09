# make with make -f BasinMetrics.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
#SOURCES=basin_metrics.cpp LSDMostLikelyPartitionsFinder.cpp LSDChiNetwork.cpp LSDIndexRaster.cpp LSDRaster.cpp LSDRasterSpectral.cpp LSDFlowInfo.cpp LSDChannelNetwork.cpp LSDIndexChannel.cpp LSDChannel.cpp LSDIndexChannelTree.cpp LSDStatsTools.cpp
SOURCES=BasinMetrics.cpp LSDMostLikelyPartitionsFinder.cpp LSDChiNetwork.cpp LSDIndexRaster.cpp LSDRaster.cpp LSDFlowInfo.cpp LSDChannelNetwork.cpp LSDIndexChannel.cpp LSDChannel.cpp LSDIndexChannelTree.cpp LSDStatsTools.cpp
#LIBS   = -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
# EXECUTABLE=basin_metrics.exe
EXECUTABLE=BasinMetrics.out

all: $(SOURCES) $(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS)
#	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
