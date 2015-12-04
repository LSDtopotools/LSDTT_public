CC = g++
CFLAGS= -c -Wreturn-type -O0 -g3 -std=gnu++11 -D_GLIBCXX_DEBUG
OFLAGS = -Wreturn-type -O0 -g3 -std=gnu++11 -D_GLIBCXX_DEBUG
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
EXEC = CatchmentModel.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@


clean:
	rm -f ../../$(OBJ) ../catchmentmodel_driver.o CatchmentModel.out
