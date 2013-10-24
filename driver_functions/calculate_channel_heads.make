# make with make -f chi_step2_write_channel_file.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=calculate_channel_heads_driver.cpp LSDMostLikelyPartitionsFinder.cpp LSDIndexRaster.cpp LSDRaster.cpp LSDFlowInfo.cpp LSDChannelNetwork.cpp LSDIndexChannel.cpp LSDChannel.cpp LSDIndexChannelTree.cpp LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=channel_heads.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
