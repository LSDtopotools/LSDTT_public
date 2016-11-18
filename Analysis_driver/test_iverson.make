# make with make -f test_iverson.make

CC=g++
CFLAGS=-c -Wall -O3 -g 
OFLAGS = -Wall -O3 -g 
LDFLAGS= -Wall
SOURCES=test_iverson.cpp \
    ../LSDPorewaterColumn.cpp \
    ../LSDPorewaterParams.cpp \
    ../LSDParameterParser.cpp \
    ../LSDStatsTools.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test_iverson.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
