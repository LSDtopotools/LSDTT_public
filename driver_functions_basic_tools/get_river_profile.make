# make with make -f chi_map.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=get_river_profile.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDFlowInfo.cpp \
		../LSDIndexChannel.cpp \
		../LSDParameterParser.cpp \
		../LSDSpatialCSVReader.cpp \
		../LSDJunctionNetwork.cpp \
		../LSDStatsTools.cpp \
		../LSDChannel.cpp \
		../LSDMostLikelyPartitionsFinder.cpp \
    ../LSDShapeTools.cpp
LIBS= -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=get_river_profile.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
