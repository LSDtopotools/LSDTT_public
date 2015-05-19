# make with make -f child_drainage_density.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=child_drainage_density.cpp \
        ../../LSDIndexRaster.cpp \
        ../../LSDRaster.cpp \
        ../../LSDFlowInfo.cpp \
        ../../LSDIndexChannel.cpp \
        ../../LSDStatsTools.cpp \
        ../../LSDJunctionNetwork.cpp \
        ../../LSDChannel.cpp \
        ../../LSDMostLikelyPartitionsFinder.cpp \
        ../../LSDBasin.cpp \
        ../../LSDShapeTools.cpp \
        ../../LSDCRNParameters.cpp \
        ../../LSDParticle.cpp
LIBS   = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=child_drainage_density.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
