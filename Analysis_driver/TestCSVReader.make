# make with make -f TestCSVReader.make

CC=g++
CFLAGS=-c -std=c++11 -Wall -O3
OFLAGS = -std=c++11 -Wall -O3
LDFLAGS= -std=c++11 -Wall
SOURCES=TestCSVReader.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDRasterInfo.cpp \
    ../LSDStatsTools.cpp \
    ../LSDShapeTools.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDSpatialCSVReader.cpp \
    ../LSDParameterParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=TestCSVReader.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
