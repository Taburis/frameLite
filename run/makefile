
#ifneq ($(MAKECMSGOALS), clean)
#		include $(SRCS:.C=.)

CC = g++
ROOT = -pthread -stdlib=libc++ -std=c++11 -m64 -I/Applications/root_v6.14.04/include -L/Applications/root_v6.14.04/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lpthread -stdlib=libc++ -lm -ldl -Wl,-rpath,/Applications/root_v6.14.04/lib 

CFLAGS = -fpic -g $(ROOT)
LDFLAGS = -shared $(ROOT)
DEPS = $(shell echo *.h)
LPATH =$(PWD)/../lib/ 
INCLUDE= -I $(PWD)/../plugin/ -I $(PWD)/../core/interface/
INCLUDE=-I/Users/tabris/frameLite/plugin -I/Users/tabris/frameLite/core/interface

SRCS=$(wildcard *.C)
OBJS=$(SRCS:%.C=%.o)
BINS=$(SRCS:%.C=%)

#all: $(OBJS) $(BINS)

$(OBJS): $(SRCS) 
	@echo $(SRCS)
	$(CC) $(ROOT) $(INCLUDE) -c $< -o $@ 

$(BINS): $(OBJS) $(PWD)/../lib/libedm.so
	$(CC) -I $(LPATH) -L $(LPATH) -ledm $(ROOT) $(INCLUDE) $< -o $@

.PHONY: clean
.PHONY: test.C

test.C:
	@echo $(OBJS)
	@echo $(BINS)


clean:
	rm -rf *.o *.*~ ../plugin/*.*~
