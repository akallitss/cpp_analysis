#------------------------------------------------------------------------------
# MAKEFILE : Read LEcroy software
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Defining the directory structure and some flags.
#------------------------------------------------------------------------------
# REST_LIB = $(REST_PATH_SED)/lib
# REST_INC = $(REST_PATH_SED)/inc

CFLAGS = $(shell root-config --cflags)
LIBS   = $(shell root-config --libs)
GLIBS  = $(shell root-config --glibs)

# Programmes to be linked.
PROGS = analyseTree
PROGS_AUX = $(PROGS:=.out)

#------------------------------------------------------------------------------
# Defining the compiling targets.
#------------------------------------------------------------------------------
all: Wellcome $(PROGS_AUX)

#------------------------------------------------------------------------------
# Wellcome message.
#------------------------------------------------------------------------------
Wellcome:
	@echo -e "\033[40m\033[0;1;32m ------------------------------ \033[0m"
	@echo -e "\033[40m\033[0;1;32m  Lecroy Read Makefile          \033[0m"
	@echo -e "\033[40m\033[0;1;32m  New Version: 31st May 2018    \033[0m"
	@echo -e "\033[40m\033[0;1;32m ------------------------------ \033[0m"
# 	@echo -e "\033[47m\033[0;1;33m ------------------------------ \033[0m"
# 	@echo -e "\033[47m\033[0;1;33m  T. Papaevangelou              \033[0m"
# 	@echo -e "\033[47m\033[0;1;33m ------------------------------ \033[0m"
	@echo

#------------------------------------------------------------------------------
# Generation of the executable.
#------------------------------------------------------------------------------
%.out:
	@echo -e "\033[40m\033[0;32m  Linking" $(@:.out=) "........\033[0m"
## 	@g++ $(CFLAGS) $(LIBS) $(GLIBS) -o bin/$(@:.out=) -I$(ROOTSYS)/include -I$(REST_INC) -L$(ROOTSYS)/lib -L$(REST_LIB) $(@:.out=.cxx) -lCore -lCint -lMatrix -lrest -lGui -lHist -lTree -lGraf -lMinuit -lPhysics -lGpad -lGeom -ldl -lGraf3d
# 	@echo -e g++ $(CFLAGS) $(LIBS) $(GLIBS) -o ./$(@:.out=) -I$(ROOTSYS)/include -L$(ROOTSYS)/lib $(@:.out=.cxx) -lCore -lTree -ldl 
	@g++ $(CFLAGS) $(LIBS) $(GLIBS) -g -o ./$(@:.out=) -I$(ROOTSYS)/include -L$(ROOTSYS)/lib $(@:.out=.cxx) -lCore -lTree -ldl  -lRIO -lNet -lHist 
# 	@g++ -o ./$(@:.out=)  $(@:.out=.cxx)

#------------------------------------------------------------------------------
# Cleaning.
#------------------------------------------------------------------------------
clean: 
	@echo -e "\033[40m\033[0;34m  Removing programmes ........\033[0m"
# 	@rm bin/*
#------------------------------------------------------------------------------
