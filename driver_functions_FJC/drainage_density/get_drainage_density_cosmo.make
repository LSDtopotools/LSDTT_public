# make with make -f get_drainage_density_cosmo.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=get_drainage_density_cosmo_driver.cpp \
    ../../LSDMostLikelyPartitionsFinder.cpp \
    ../../LSDIndexRaster.cpp \
    ../../LSDRaster.cpp \
    ../../LSDFlowInfo.cpp \
    ../../LSDJunctionNetwork.cpp \
    ../../LSDIndexChannel.cpp \
    ../../LSDChannel.cpp \
    ../../LSDIndexChannelTree.cpp \
    ../../LSDStatsTools.cpp \
    ../../LSDShapeTools.cpp \
    ../../LSDBasin.cpp \
    ../../LSDParticle.cpp \
    ../../LSDCRNParameters.cpp
LIBS   = -lm -lstdc++
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=get_drainage_density_cosmo.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
