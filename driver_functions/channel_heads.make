# make with make -f chi_map.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=chi_map.cpp ../LSDIndexRaster.cpp ../LSDRaster.cpp ../LSDFlowInfo.cpp ../LSDStatsTools.cpp
LIBS= -lm -lstdc++ -lfftw3
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=chi_map.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
