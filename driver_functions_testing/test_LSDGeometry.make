# make with make -f basin_grabber.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=test_LSDGeometry.cpp \
    ../LSDGeometry.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDRasterInfo.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDStatsTools.cpp \
    ../LSDShapeTools.cpp


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test_LSDGeometry.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
