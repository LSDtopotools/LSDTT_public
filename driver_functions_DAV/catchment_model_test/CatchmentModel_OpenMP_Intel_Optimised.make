CC = icc 
CFLAGS= -c -Wall -O2 -std=c++11 -openmp
OFLAGS = -Wall -O2 -std=c++11 -openmp
LDFLAGS= -Wreturn-type
SOURCES = ../catchmentmodel_driver.cpp \
			../../LSDCatchmentModel.cpp \
			../../LSDRaster.cpp \
			../../LSDIndexRaster.cpp \
			../../LSDStatsTools.cpp \
			../../LSDShapeTools.cpp
SCRIPTS = 
OBJ = $(SOURCES:.cpp=.o)
#LIBS = -lfftw3 -lpython2.7 -g -O0 -D_GLIBCXX_DEBUG
#LIBS = -lfftw3 -lpython2.7 -Wwrite-strings
LIBS = 
EXEC = CatchmentModel_OpenMP_Intel_Optimised.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@


clean:
	rm -f ../../$(OBJ) ../catchmentmodel_driver.o CatchmentModel_OpenMP_Intel_Optimised.out
