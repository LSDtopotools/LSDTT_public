# make with make -f get_CHT.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=get_CHT.cpp \
        ../LSDMostLikelyPartitionsFinder.cpp \
        ../LSDIndexRaster.cpp \
        ../LSDRaster.cpp \
        ../LSDRasterSpectral.cpp \
        ../LSDFlowInfo.cpp \
        ../LSDJunctionNetwork.cpp \
        ../LSDIndexChannel.cpp \
        ../LSDChannel.cpp \
        ../LSDIndexChannelTree.cpp \
        ../LSDStatsTools.cpp \
        ../LSDChiNetwork.cpp \
        ../LSDShapeTools.cpp \
        ../LSDBasin.cpp \
        ../LSDCosmoData.cpp \
        ../LSDCRNParameters.cpp \
        ../LSDRasterInfo.cpp \
        ../LSDParticle.cpp
LIBS= -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=get_CHT.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
