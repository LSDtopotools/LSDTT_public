# make with make -f write_dreich_junctions.make

CC=g++
CFLAGS=-c -Wall -O3 
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=write_dreich_junctions.cpp \
             ../../LSDMostLikelyPartitionsFinder.cpp \
             ../../LSDIndexRaster.cpp \
             ../../LSDRaster.cpp \
             ../../LSDFlowInfo.cpp \
             ../../LSDJunctionNetwork.cpp \
             ../../LSDIndexChannel.cpp \
             ../../LSDChannel.cpp \
             ../../LSDIndexChannelTree.cpp \
             ../../LSDStatsTools.cpp \
             ../../LSDShapeTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=write_dreich_junctions.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
