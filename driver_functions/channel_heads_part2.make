# make with make -f channel_heads_part2.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=channel_heads_part2_driver.cpp LSDIndexRaster.cpp LSDRaster.cpp LSDFlowInfo.cpp LSDIndexChannel.cpp LSDStatsTools.cpp LSDRasterSpectral.cpp LSDChannelNetwork.cpp LSDChannel.cpp LSDMostLikelyPartitionsFinder.cpp
LIBS   = -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=channel_heads_part2.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
