# make with make -f chi_mapping_tool.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=basin_info_litho.cpp \
             ../LSDMostLikelyPartitionsFinder.cpp \
             ../LSDIndexRaster.cpp \
             ../LSDRaster.cpp \
             ../LSDRasterInfo.cpp \
             ../LSDFlowInfo.cpp \
             ../LSDJunctionNetwork.cpp \
             ../LSDIndexChannel.cpp \
             ../LSDChannel.cpp \
             ../LSDIndexChannelTree.cpp \
             ../LSDStatsTools.cpp \
             ../LSDShapeTools.cpp \
             ../LSDChiNetwork.cpp \
             ../LSDBasin.cpp \
             ../LSDParticle.cpp \
             ../LSDChiTools.cpp \
             ../LSDParameterParser.cpp \
             ../LSDSpatialCSVReader.cpp \
             ../LSDCRNParameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=lithobasinfo.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ../*.o *.o *.out *.exe
