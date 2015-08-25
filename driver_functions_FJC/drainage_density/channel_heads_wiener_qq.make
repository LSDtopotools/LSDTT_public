# make with make -f channel_heads_wiener_qq.make

CC=g++
CFLAGS=-c -Wall -O3 
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=channel_heads_wiener_qq_driver.cpp ../../LSDMostLikelyPartitionsFinder.cpp ../../LSDIndexRaster.cpp ../../LSDRaster.cpp ../../LSDRasterSpectral.cpp ../../LSDFlowInfo.cpp ../../LSDJunctionNetwork.cpp ../../LSDIndexChannel.cpp ../../LSDChannel.cpp ../../LSDIndexChannelTree.cpp ../../LSDStatsTools.cpp ../../LSDChiNetwork.cpp ../../LSDShapeTools.cpp
LIBS= -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=channel_heads_wiener_qq.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
