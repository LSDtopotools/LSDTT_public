# make with make -f fill_nodata_holes.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=change_feet_to_metres.cpp \
    ../../LSDIndexRaster.cpp \
    ../../LSDRaster.cpp \
    ../../LSDShapeTools.cpp \
    ../../LSDStatsTools.cpp
LIBS = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=change_feet_to_metres.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
