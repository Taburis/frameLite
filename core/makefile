
CC = g++
ROOT = -pthread -stdlib=libc++ -std=c++11 -m64 -I/Applications/root_v6.14.04/include -L/Applications/root_v6.14.04/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lpthread -stdlib=libc++ -lm -ldl -Wl,-rpath,/Applications/root_v6.14.04/lib 
CFLAGS = -fpic -g $(ROOT)
LDFLAGS = -shared $(ROOT)
DEPS = $(shell echo *.h)
LPATH = $(shell echo /Users/tabris/xflow/core/src/)
INCLUDE=-I /Users/tabris/cmsProjects/frameLite/core/src/ -I /Users/tabris/cmsProjects/frameLite/plugin/
SRC=src
BIN=bin
#CFLAGSFORDYLIB = -c -shared -fpic $(ROOT)

all: bin/rootEDM.o libedm.so

#%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $(CFLAGS)


libedm.so: bin/rootEDM.o
	$(CC) -shared bin/rootEDM.o -o ../lib/libedm.so $(ROOT)

bin/rootEDM.o: src/rootEDM.cc src/rootEDM.h
	$(CC) -o bin/rootEDM.o src/rootEDM.cc -c -fpic $(ROOT)


.PHONY: clean

clean:
	rm -rf bin/*.o ./*.*~ src/*.*~
