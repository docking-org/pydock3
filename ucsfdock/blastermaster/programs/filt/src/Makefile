# makefile for DOCK 3.5 chemgrid directory
#
MAKE = make

F77 = pgf77
FFLAGS = -byteswapio -tp p6 -fast -O2 -I.. -Bstatic 

all: filt

FILT = ../bin/filt

clean:
	/bin/rm -f $(FILT)

filt: 
	$(F77) $(FFLAGS) filter1.f -o $(FILT)

