# make with make -f channel_heads_part2.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=channel_heads_part2_driver.cpp LSDMostLikelyPartitionsFinder.cpp LSDIndexRaster.cpp LSDRaster.cpp LSDFlowInfo.cpp LSDChannelNetwork.cpp LSDIndexChannel.cpp LSDChannel.cpp LSDIndexChannelTree.cpp LSDStatsTools.cpp LSDRasterSpectral.cpp
LIBS= -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=channel_heads_part2.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
