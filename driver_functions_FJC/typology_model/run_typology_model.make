# make with make -f run_typology_model.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=run_typology_model.cpp ../../LSDIndexRaster.cpp ../../LSDRaster.cpp ../../LSDFlowInfo.cpp ../../LSDIndexChannel.cpp ../../LSDStatsTools.cpp ../../LSDJunctionNetwork.cpp ../../LSDChannel.cpp ../../LSDMostLikelyPartitionsFinder.cpp ../../LSDShapeTools.cpp ../../LSDSpatialCSVReader.cpp
LIBS   = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=run_typology_model.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
