
GRID_DIR := /data/Dropbox/griddock/
export GRID_DIR

DOCKDIR := griddock
DTSSDIR := dtess
CMMNDIR := gridcommon
SUBDIRS := $(CMMNDIR) $(DTSSDIR) $(DOCKDIR)

# Default Target

.PHONY: all dock libs dtess clean backup

all:
	@echo
	@echo "#######################################"
	@echo "### BUILDING ALL TARGETS ###"
	@echo "#######################################"
	@echo
	for i in $(SUBDIRS) ; do \
	( cd $$i ; make ) ; \
	done

# Rules to build individual targets

libs:
	@cd $(CMMNDIR) ; make -B

dock:
	@cd $(CMMNDIR) ; make
	@cd $(DOCKDIR) ; make

dtess:
	@cd $(CMMNDIR) ; make
	@cd $(DTSSDIR) ; make

clean:
	rm -f *~
	for i in $(SUBDIRS) ; do \
	( cd $$i ; make clean) ; \
	done

backup:
	for i in $(SUBDIRS) ; do \
	( cd $$i ; make backup ) ; \
	done

