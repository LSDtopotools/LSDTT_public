# make with make -f CRN_column_test.make
# This makes a testing function for a CRN column 

CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES = CRN_column_test.cpp \
		../LSDIndexRaster.cpp \
		../LSDShapeTools.cpp \
		../LSDRaster.cpp \
		../LSDStatsTools.cpp \
		../LSDFlowInfo.cpp \
		../LSDParticle.cpp \
		../LSDParticleColumn.cpp \
		../LSDCRNParameters.cpp
OBJ = $(SOURCES:.cpp=.o)
EXEC = CRN_column_test.exe

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
