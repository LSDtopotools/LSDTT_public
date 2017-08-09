CC = g++
CFLAGS= -c -I/home/smudd/libs/boost_1_64_0 -Wall -O3
OFLAGS = -I/home/smudd/libs/boost_1_64_0 -Wall -O3
LDFLAGS= -Wall
SOURCES = model_transient_fastscape.cpp \
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
EXEC = model_transient_fastscape.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
