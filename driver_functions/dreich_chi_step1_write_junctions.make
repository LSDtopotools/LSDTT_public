# make with make -f dreich_chi_step1_write_junctions.make

CC=g++
CFLAGS=-c -Wall -O3 
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=dreich_chi_step1_write_junctions_driver.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDJunctionNetwork.cpp ../LSDIndexChannel.cpp ../LSDChannel.cpp ../LSDIndexChannelTree.cpp ../LSDStatsTools.cpp ../LSDShapeTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=dreich_chi1_write_junctions.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
