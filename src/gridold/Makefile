# Makefile for GridDock
# Version 0.1a | Summer 2009 | Stephen J. Bush
#  > single class, basic functionality

# Version 0.1 | Winter 2010 (?) | Stephen J. Bush
#  > expansion of bDocker, bGrid, and bPoints

# Version 0.1.1 | Winter 2010 (?) | Stephen J. Bush
#  > expansion of bStamp

# Version 0.1.2 | March 2010 | Stephen J. Bush
#  > addition of bSort

# Version 0.2 | April (2010) | Stephen J. Bush
#  > addition of bDelTess
#  > addition of bHex
#  > inclusion of namespaces

# Helpful Hints
# gcc -MM *.cpp
# http://www.redantigua.com/c-make-automatic-variables.html

# Compiler Variables
CXX=g++
CXX_FLAGS=-fopenmp -Wall -g -O0
LD_FLAGS=
# -O2 -pg

# Files
OBJECTS_NODEP = bStamp.o bSort.o bHex.o bSnapp.o bAA.o bSys.o bCalc.o
OBJECTS_COREDEP = bPoints.o bSimplex.o 
OBJECTS_ALLDEP = bList.o bGrid.o bDelTess.o bCentroid.o bMC.o bDocker.o
OBJECTS_DT = bSort.o bSys.o bHex.o bAA.o bPoints.o bSnapp.o bSimplex.o bDelTess.o bCentroid.o
OBJECTS = $(OBJECTS_NODEP) $(OBJECTS_COREDEP) $(OBJECTS_ALLDEP)
HEADERS = $(OBJECTS:.o=.h)
SOURCES = $(OBJECTS:.o=.cpp)
EXTRA_FILES = main.cpp Makefile dtess.cpp
DOCK = griddock.e
SCORE = score.e
DTESS = dtess

all: $(DOCK)

score: $(SCORE)

dtess: $(DTESS)

$(DOCK): main.o $(OBJECTS)
	$(CXX) $(CXX_FLAGS) main.o $(OBJECTS) -o $(DOCK)
	@echo ' '

main.o: main.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c $<
	@echo ' '

$(DTESS): dtess.o $(OBJECTS_DT)
	$(CXX) $(CXX_FLAGS) dtess.o $(OBJECTS_DT) -o $(DTESS)
	@echo ' '

dtess.o: dtess.cpp $(OBJECTS_DT)
	$(CXX) $(CXX_FLAGS) -c $<
	@echo ' '

$(SCORE): score.o $(OBJECTS)
	$(CXX) $(CXX_FLAGS) score.o $(OBJECTS) -o $(SCORE)
	@echo ' '

score.o: score.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c $<
	@echo ' '

# not exactly efficient, but it works
$(OBJECTS_NODEP): %.o: %.cpp %.h
	$(CXX) $(CXX_FLAGS) -c $<
	@echo ' '

$(OBJECTS_COREDEP): %.o: %.cpp %.h
	$(CXX) $(CXX_FLAGS) -c $<
	@echo ' '

$(OBJECTS_ALLDEP): %.o: %.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c $<
	@echo ' '

clean:
	rm -f $(DOCK) $(SCORE) $(OBJECTS)

backup:
	@-if [ ! -e .backup ]; then mkdir .backup; fi;
	@zip .backup/backup_`date +%d-%m-%y_%H.%M`.zip $(SOURCES) $(HEADERS) $(EXTRA_FILES)
