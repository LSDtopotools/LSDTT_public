# make with make -f compare_slopes_at_channel_head.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=compare_slopes_at_channel_head.cpp \
    ../../LSDIndexRaster.cpp \
    ../../LSDRaster.cpp \
    ../../LSDFlowInfo.cpp \
    ../../LSDShapeTools.cpp \
    ../../LSDStatsTools.cpp
LIBS = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=compare_slopes_at_channel_head.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
