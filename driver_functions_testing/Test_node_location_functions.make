# make with make -f basin_grabber.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=Test_node_location_functions.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDStatsTools.cpp \
    ../LSDShapeTools.cpp \
    ../LSDJunctionNetwork.cpp \
    ../LSDIndexChannel.cpp \
    ../LSDChannel.cpp \
    ../LSDMostLikelyPartitionsFinder.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Test_node_location_functions.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
