# make with make -f basin_puncher.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=basin_puncher_driver.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDIndexChannel.cpp ../LSDStatsTools.cpp ../LSDJunctionNetwork.cpp ../LSDChannel.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDBasin.cpp
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=basin_puncher.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
