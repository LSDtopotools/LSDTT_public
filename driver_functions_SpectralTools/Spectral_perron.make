# make with make -f Spectral_perron.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=Spectral_analysis_Perron_driver.cpp \
        ../LSDIndexRaster.cpp \
        ../LSDRaster.cpp \
        ../LSDRasterInfo.cpp \
        ../LSDFlowInfo.cpp \
        ../LSDIndexChannel.cpp \
        ../LSDStatsTools.cpp \
        ../LSDRasterSpectral.cpp \
        ../LSDJunctionNetwork.cpp \
        ../LSDMostLikelyPartitionsFinder.cpp \
        ../LSDChannel.cpp \
        ../LSDShapeTools.cpp
LIBS   = -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Spectral_perron.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
