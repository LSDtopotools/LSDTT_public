# make with make -f change_raster_value.make

CC=g++
CFLAGS=-c -Wall -O3 
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=change_raster_value_driver.cpp \
             ../LSDMostLikelyPartitionsFinder.cpp \
             ../LSDIndexRaster.cpp \
             ../LSDRaster.cpp \
             ../LSDFlowInfo.cpp \
             ../LSDJunctionNetwork.cpp \
             ../LSDIndexChannel.cpp \
             ../LSDChannel.cpp \
             ../LSDIndexChannelTree.cpp \
             ../LSDStatsTools.cpp \
             ../LSDShapeTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=change_raster_value.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
