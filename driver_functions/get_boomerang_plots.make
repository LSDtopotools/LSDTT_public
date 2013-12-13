# make with make -f get_boomerang_plots.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=get_boomerang_plots_driver.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDIndexChannel.cpp ../LSDStatsTools.cpp ../LSDRasterSpectral.cpp ../LSDChannelNetwork.cpp ../LSDChannel.cpp ../LSDMostLikelyPartitionsFinder.cpp ../LSDBasin.cpp
LIBS   = -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=Chile_test3.exe
EXECUTABLE=get_boomerang_plots.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
