CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -I/home/s0923330/libs/boost_1_64_0 -Wall -O3
LDFLAGS= -Wall
SOURCES = stream_power_model.cpp \
		../../LSDRasterSpectral.cpp \
		../../LSDIndexRaster.cpp \
		../../LSDShapeTools.cpp \
		../../LSDRaster.cpp \
		../../LSDRasterModel.cpp \
		../../LSDStatsTools.cpp \
		../../LSDParticleColumn.cpp \
		../../LSDParticle.cpp \
		../../LSDCRNParameters.cpp \
		../../LSDFlowInfo.cpp
SCRIPTS = animate.py
OBJ = $(SOURCES:.cpp=.o)
#LIBS = -lfftw3 -lpython2.7 -g -O0 -D_GLIBCXX_DEBUG
LIBS = -lfftw3 -lpython2.7 -Wwrite-strings
EXEC = stream_power_model.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
