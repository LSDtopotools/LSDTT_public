CC = g++
CFLAGS= -c -Wall -O3 -pg -std=gnu++11
OFLAGS = -Wall -O3 -pg -std=gnu++11
LDFLAGS= -Wall
SOURCES = catchmentmodel_driver.cpp \
			../LSDCatchmentModel.cpp \
			../LSDRaster.cpp \
			../LSDIndexRaster.cpp \
			../LSDStatsTools.cpp \
			../LSDShapeTools.cpp
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
