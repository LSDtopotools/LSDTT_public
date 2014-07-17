# make with make -f get_estar_rstar.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=get_estar_rstar_driver.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDIndexChannel.cpp ../LSDStatsTools.cpp ../LSDJunctionNetwork.cpp ../LSDChannel.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDBasin.cpp ../LSDShapeTools.cpp
LIBS   = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=get_estar_rstar.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
