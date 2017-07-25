CC=g++
CFLAGS=-c -Wall -O3 -pg -g
OFLAGS = -Wall -O3 -pg -g
LDFLAGS= -Wall
SOURCES= Hillslope_Morphology_Driver.cpp \
        ../LSDMostLikelyPartitionsFinder.cpp \
        ../LSDIndexRaster.cpp \
        ../LSDRaster.cpp \
        ../LSDFlowInfo.cpp \
        ../LSDJunctionNetwork.cpp \
        ../LSDIndexChannel.cpp \
        ../LSDChannel.cpp \
        ../LSDStatsTools.cpp \
        ../LSDBasin.cpp \
        ../LSDShapeTools.cpp \
        ../LSDParticle.cpp \
        ../LSDParameterParser.cpp \
        ../LSDCRNParameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE= Hillslope_Morphology_Tool.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -f ../*.o *.o *.out *.exe
