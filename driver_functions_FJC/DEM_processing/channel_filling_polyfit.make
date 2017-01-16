# make with make -f channel_filling_polyfit.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=channel_filling_polyfit.cpp ../../LSDIndexRaster.cpp ../../LSDRaster.cpp ../../LSDFlowInfo.cpp ../../LSDIndexChannel.cpp ../../LSDStatsTools.cpp ../../LSDJunctionNetwork.cpp ../../LSDChannel.cpp ../../LSDMostLikelyPartitionsFinder.cpp ../../LSDShapeTools.cpp
LIBS   = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=channel_filling_polyfit.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
