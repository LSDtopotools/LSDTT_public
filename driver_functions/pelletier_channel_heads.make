# make with make -f pelletier_channel_heads.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=pelletier_channel_heads_driver.cpp LSDIndexRaster.cpp LSDRaster.cpp LSDFlowInfo.cpp LSDIndexChannel.cpp LSDStatsTools.cpp LSDRasterSpectral.cpp LSDChannelNetwork.cpp LSDChannel.cpp LSDMostLikelyPartitionsFinder.cpp
LIBS   = -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=wiener_filter.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
