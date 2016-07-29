CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
SOURCES = fill_no_data.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDShapeTools.cpp \
    ../LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fill_no_data.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
