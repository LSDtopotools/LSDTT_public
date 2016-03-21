# make with make -f stream_ordering.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=stream_ordering_driver.cpp \
    ../../LSDIndexRaster.cpp \
    ../../LSDRaster.cpp \
    ../../LSDFlowInfo.cpp \
    ../../LSDShapeTools.cpp \
		../../LSDJunctionNetwork.cpp \
		../../LSDIndexChannel.cpp \
		../../LSDChannel.cpp \
		../../LSDMostLikelyPartitionsFinder.cpp \
    ../../LSDStrahlerLinks.cpp \
    ../../LSDStatsTools.cpp
LIBS = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=stream_ordering.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
