# make with make -f TestCSVReader.make

CC=g++
CFLAGS=-c -Wall -O3 -g 
OFLAGS = -Wall -O3 -g 
LDFLAGS= -Wall
SOURCES=TestCSVReader.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDRasterInfo.cpp \
    ../LSDStatsTools.cpp \
    ../LSDShapeTools.cpp \
    ../LSDSpatialCSVReader.cpp \
    ../LSDParameterParser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=basin_averager.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
