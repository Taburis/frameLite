
ROOT = -pthread -stdlib=libc++ -std=c++11 -m64 -I/Applications/root_v6.14.04/include -L/Applications/root_v6.14.04/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lpthread -stdlib=libc++ -lm -ldl -Wl,-rpath,/Applications/root_v6.14.04/lib 
#CFLAGS = -fpic -g $(ROOT)
#LDFLAGS = -shared $(ROOT)
#LPATH = $(shell echo /Users/tabris/xflow/core/src/)
CC = g++
INCLUDE= -I $(PWD)/interface/ -I $(PWD)/../plugin/
SRCS:= $(wildcard src/*.cc)
BINS:= $(SRCS:src/%.cc=bin/%.o)
#CFLAGSFORDYLIB = -c -shared -fpic $(ROOT)
all: $(BINS) libedm.so

bin/%.o: src/%.cc
	@echo "=====  Creating object " $@ "=========="
	$(CC) -c $< -o $@ -fpic $(INCLUDE) $(ROOT)

libedm.so: $(wildcard bin/*.o)
	$(CC) -shared bin/rootEDM.o -o ../lib/libedm.so $(ROOT) $(INCLUDE)


.PHONY: clean
.PHONY: test


clean:
	rm -rf bin/*.o ./*.*~ src/*.*~
