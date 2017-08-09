CC = g++
CFLAGS= -c -I/LSDTopoTools/boost_mtl_minimal -Wall -O3
OFLAGS = -I/LSDTopoTools/boost_mtl_minimal -Wall -O3
LDFLAGS= -Wall
SOURCES = model_divide_migration.cpp \
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
#LIBS = -lfftw3 -g -O0 -D_GLIBCXX_DEBUG
LIBS = -lfftw3 -Wwrite-strings
EXEC = model_divide_migration.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
