# make with make -f difference_DEMS.make
# Prints the difference between DEMS to screen 

CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
SOURCES = difference_DEMS.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDShapeTools.cpp \
    ../LSDRaster.cpp \
    ../LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=difference_DEMs.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
