# make with make -f longitudinal_channel_swath_erosion_profile.make

CC=g++
CFLAGS= -c -std=c++11 -Wall -O3
OFLAGS = -Wall -lboost_system -O3 
LDFLAGS= -Wall -lboost_system
SOURCES=longitudinal_channel_swath_erosion_profile.cpp \
../LSDMostLikelyPartitionsFinder.cpp \
../LSDIndexRaster.cpp \
../LSDRaster.cpp \
../LSDFlowInfo.cpp \
../LSDJunctionNetwork.cpp \
../LSDIndexChannel.cpp \
../LSDChannel.cpp \
../LSDIndexChannelTree.cpp \
../LSDStatsTools.cpp \
../LSDShapeTools.cpp \
../LSDSwathProfile.cpp \
../LSDCloudBase.cpp
LIBS=-L/usr/local/lib/liblas.so -L/usr/lib64/ImageMagick-6.8.8/modules-Q16/coders/pcl.so
INC=-I/usr/local/include/liblas -I/usr/include/pcl-1.7


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=long_profile_erosion_swath.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(INC) $(LIBS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $(LIBS) $< -o $@
