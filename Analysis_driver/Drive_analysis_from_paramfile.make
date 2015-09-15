# make with make -f Drive_analysis_from_paramfile.make

CC=g++
CFLAGS=-c -Wall -O3 -g 
OFLAGS = -Wall -O3 -g 
LDFLAGS= -Wall
SOURCES=Drive_analysis_from_paramfile.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDStatsTools.cpp \
    ../LSDJunctionNetwork.cpp \
    ../LSDIndexChannel.cpp \
    ../LSDChannel.cpp \
    ../LSDMostLikelyPartitionsFinder.cpp \
    ../LSDShapeTools.cpp \
    ../LSDAnalysisDriver.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=LSDTT_analysis_from_paramfile.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
