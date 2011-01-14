# Variable definition
OBJS_ND = bVBase.o bSort.o bHex.o bSnapp.o bAA.o bSys.o bCalc.o
OBJS_CD = bPoints.o bSimplex.o
OBJS_AD = bDelTess.o bCentroid.o
OBJS = $(OBJS_ND) $(OBJS_CD) $(OBJS_AD)
HDRS_ND = $(OBJS_ND:.o=.h)
HDRS = $(OBJS:.o=.h)
SRCS = $(OBJS:.o=.cpp)
XTRA = Makefile

LIBCOMMON = libgridcommon.so
INSTALLDIR = $(GRID_DIR)/lib/

CXX=g++
CXX_FLAGS=-fPIC -fopenmp -Wall -c -g -O3
LD_FLAGS=-shared -Wl,-soname,$(LIBCOMMON).1

default: $(LIBCOMMON)

$(LIBCOMMON): $(OBJS)
	@echo ' '
	$(CXX) $(LD_FLAGS) -o $(LIBCOMMON).1.0 $(OBJS)
	@mv $(LIBCOMMON).1.0 $(INSTALLDIR)
	@ln -sf $(INSTALLDIR)$(LIBCOMMON).1.0 $(INSTALLDIR)$(LIBCOMMON)
	@ln -sf $(INSTALLDIR)$(LIBCOMMON).1.0 $(INSTALLDIR)$(LIBCOMMON).1
	@echo ' '

$(OBJS_ND): %.o: %.cpp %.h
	$(CXX) $(CXX_FLAGS) $<
	@echo ' '

$(OBJS_CD): %.o: %.cpp %.h $(HDRS_ND)
	$(CXX) $(CXX_FLAGS) $<
	@echo ' '

$(OBJS_AD): %.o: %.cpp $(HDRS)
	$(CXX) $(CXX_FLAGS) $<
	@echo ' '

#~ install:
#~ cp $(LIBCOMMON) $(INSTALLDIR)

clean:
	rm -f $(OBJS) $(LIBCOMMON) *~

backup:
	@-if [ ! -e .backup ]; then mkdir .backup; fi;
	@zip .backup/backup_`date +%d-%m-%y_%H.%M`.zip $(SRCS) $(HDRS) $(XTRA)

