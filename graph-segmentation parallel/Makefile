INCDIR = -I.
DBG    = -g
OPT    = -O3
CPP    = g++
OPENMP = -fopenmp
CFLAGS = $(DBG) $(OPT) $(INCDIR) $(OPENMP)
LINK   = -lm 

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@

all: segment

segment: segment.cpp segment-image.h segment-graph.h disjoint-set.h convolve.h
	$(CPP) $(CFLAGS) -o segment segment.cpp cycletimer.c $(LINK)

clean:
	/bin/rm -f segment *.o

clean-all: clean
	/bin/rm -f *~ 



