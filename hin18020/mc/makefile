
#ifneq ($(MAKECMSGOALS), clean)
#		include $(SRCS:.C=.)

CC = g++
ROOT = -pthread -stdlib=libc++ -std=c++11 -m64 -I/Applications/root_v6.14.04/include -L/Applications/root_v6.14.04/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lpthread -stdlib=libc++ -lm -ldl -Wl,-rpath,/Applications/root_v6.14.04/lib 

CFLAGS = -fpic -g $(ROOT)
LDFLAGS = -shared $(ROOT)
DEPS = $(shell echo *.h)
FBASE = /Users/tabris/frameLite/
LPATH =$(FBASE)/lib/ 
INCLUDE= -I $(FBASE)/plugin/ -I $(FBASE)/core/interface/
INCLUDE=-I/Users/tabris/frameLite/plugin -I/Users/tabris/frameLite/core/interface
PYTHIA8PATH=/Applications/physPlayGround/pythia8235
MOREPYTHIA8=-I$(PYTHIA8PATH)/include -O3 -pedantic -W -Wall -Wshadow -L$(PYTHIA8PATH)/lib -Wl,-rpath,$(PYTHIA8PATH)/lib -lpythia8 -ldl

FASTJETPATH = /Applications/physPlayGround/fastjet-install
FASTJET=$(shell $(FASTJETPATH)/bin/fastjet-config --cxxflags --libs --plugins)

SRCS=$(wildcard *.C)
OBJS=$(SRCS:%.C=%.o)
BINS=$(SRCS:%.C=%)

all: $(OBJS) $(BINS)

$(OBJS): $(SRCS) 
	@echo $(SRCS)
	$(CC) $(ROOT) $(MOREPYTHIA8) $(INCLUDE) $(FASTJET) -c $< -o $@ 

$(BINS): $(OBJS) $(FBASE)/lib/libedm.so
	$(CC) -I $(LPATH) -L $(LPATH) -ledm $(ROOT) $(MOREPYTHIA8) $(FASTJET) $(INCLUDE) $< -o $@

.PHONY: clean

clean:
	@echo $(shell $(FASTJETPATH)/bin/fastjet-config)
	rm -rf *.o *.*~ ../plugin/*.*~
