CC = g++
CFLAGS= -c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES = model_with_CRN_from_initial2.cpp \
		../LSDRasterSpectral.cpp \
		../LSDIndexRaster.cpp \
		../LSDShapeTools.cpp \
		../LSDRaster.cpp \
		../LSDRasterModel.cpp \
		../LSDStatsTools.cpp \
		../LSDFlowInfo.cpp \
		../LSDParticle.cpp \
		../LSDParticleColumn.cpp \
		../LSDCRNParameters.cpp
SCRIPTS = animate.py
OBJ = $(SOURCES:.cpp=.o)
#LIBS = -lfftw3 -lpython2.7 -g -O0 -D_GLIBCXX_DEBUG
LIBS = -lfftw3 -lpython2.7 -Wwrite-strings
EXEC = model_with_CRN_from_initial2.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
