# make with make -f basin_averager.make

CC=g++
CFLAGS=-c -Wall -O3 
OFLAGS = -Wall -O3 
LDFLAGS= -Wall
SOURCES=statistest.cpp \
    ../../LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=testistics.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
