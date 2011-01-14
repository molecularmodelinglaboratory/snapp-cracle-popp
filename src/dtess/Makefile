
# Variable definition
OBJS =
HDRS = $(OBJS:.o=.h)
SRCS = $(OBJS:.o=.cpp)
DRVR = dtess
XTRA = $(DRVR).cpp Makefile

CXX=g++
CXX_FLAGS=-Wall -g -O3 -fopenmp

EXE = dtess
LIBCMMN=libbcommon.so
LIBDIR = $(GRID_DIR)lib/
SRCDIR = $(GRID_DIR)src/
BINDIR = $(GRID_DIR)bin/
INCLUDES = -I../gridcommon/
LDFLAGS = -L$(LIBDIR) -Wl,-rpath,$(LIBDIR)

all: $(EXE)

$(EXE): $(DRVR).o $(OBJS)
	@echo ' '
	$(CXX) $(CXX_FLAGS) $(LDFLAGS) $(DRVR).o $(OBJS) -lgridcommon -o $(EXE).e
	@mv ./$(EXE).e ../
	@ln -sf $(SRCDIR)$(EXE).e $(BINDIR)$(EXE)
	@echo ' '

$(DRVR).o: $(DRVR).cpp $(OBJS)
	$(CXX) $(CXX_FLAGS) $(INCLUDES) -c $<
	@echo ' '

$(OBJS): %.o: %.cpp %.h
	$(CXX) $(CXX_FLAGS) $(INCLUDES) -c $<
	@echo ' '

clean:
	rm -f $(OBJS) $(DRVR).o *~

backup:
	@-if [ ! -e .backup ]; then mkdir .backup; fi;
	@zip .backup/backup_`date +%d-%m-%y_%H.%M`.zip $(SRCS) $(HDRS) $(XTRA)
