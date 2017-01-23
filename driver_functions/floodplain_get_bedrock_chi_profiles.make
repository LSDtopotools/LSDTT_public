# make with make -f floodplain_get_bedrock_chi_profiles.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=floodplain_get_bedrock_chi_profiles_driver.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDChiNetwork.cpp ../LSDStatsTools.cpp ../LSDShapeTools.cpp ../LSDJunctionNetwork.cpp ../LSDFlowInfo.cpp ../LSDRaster.cpp ../LSDIndexRaster.cpp ../LSDIndexChannelTree.cpp ../LSDChannel.cpp ../LSDIndexChannel.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=floodplain_get_bedrock_chi_profiles.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
