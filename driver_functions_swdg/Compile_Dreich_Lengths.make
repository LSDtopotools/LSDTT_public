CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
SOURCES = Compile_Dreich_Lengths.cpp \
    ../LSDIndexRaster.cpp \
    ../LSDRaster.cpp \
    ../LSDFlowInfo.cpp \
    ../LSDShapeTools.cpp \
    ../LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Compile_Dreich_Lengths.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
