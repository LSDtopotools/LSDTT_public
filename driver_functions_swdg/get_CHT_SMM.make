# make with make -f get_CHT_SMM.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=get_CHT_SMM.cpp \
        ../LSDMostLikelyPartitionsFinder.cpp \
        ../LSDStatsTools.cpp \
        ../LSDRasterInfo.cpp \
        ../LSDParameterParser.cpp \
        ../LSDIndexRaster.cpp \
        ../LSDRaster.cpp \
        ../LSDRasterSpectral.cpp \
        ../LSDFlowInfo.cpp \
        ../LSDJunctionNetwork.cpp \
        ../LSDIndexChannel.cpp \
        ../LSDChannel.cpp \
        ../LSDIndexChannelTree.cpp \
        ../LSDChiNetwork.cpp \
        ../LSDShapeTools.cpp \
        ../LSDBasin.cpp \
        ../LSDCosmoData.cpp \
        ../LSDCRNParameters.cpp \
        ../LSDParticle.cpp
LIBS= -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=get_CHT.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
