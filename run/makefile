
CC = g++
ROOT = -pthread -stdlib=libc++ -std=c++11 -m64 -I/Applications/root_v6.14.04/include -L/Applications/root_v6.14.04/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lpthread -stdlib=libc++ -lm -ldl -Wl,-rpath,/Applications/root_v6.14.04/lib 
CFLAGS = -fpic -g $(ROOT)
LDFLAGS = -shared $(ROOT)
DEPS = $(shell echo *.h)
LPATH = $(shell echo /Users/tabris/cmsProjects/frameLite/lib/)
INCLUDE=-I /Users/tabris/cmsProjects/frameLite/core/src/ -I /Users/tabris/cmsProjects/frameLite/plugin/
#CFLAGSFORDYLIB = -c -shared -fpic $(ROOT)

all: main.o main

#%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $(CFLAGS)

#main.o: main.cc 
#	$(CC) -o main.o -c main.cc $(ROOT) $(INCLUDE)

main.o: wtaAxis.C 
	$(CC) -o main.o -c wtaAxis.C $(ROOT) $(INCLUDE)
main: main.o /Users/tabris/cmsProjects/frameLite/lib/libedm.so
	$(CC) -o runObj main.o -I $(LPATH) -L $(LPATH) -ledm $(ROOT) $(INCLUDE) 


.PHONY: clean

clean:
	rm -rf *.o 
