# make with make -f model_with_CRN_from_initial_variable.make
# You will need libraries for this: Boost and FFTW. 
# I am afraid these are not the easiest to install (sorry).
# You can go here for details:
# http://www.geos.ed.ac.uk/~smudd/LSDTT_docs/html/external_instructions.html
# Also your life will be easier if you have a linux operating system. 
# You can install one using a virtual machine on your windows or fruit-based computer:
# http://www.geos.ed.ac.uk/~smudd/TopoTutorials/html/outside_edin.html
# Some degree of patience is required: we have tried to make this software
# as portable as possible but for some computationally intensive tasks libraries
# are needed. 

CC = g++
CFLAGS= -c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES = model_with_CRN_from_initial_variable_forcing.cpp \
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
EXEC = model_with_CRN_from_initial_variable.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
