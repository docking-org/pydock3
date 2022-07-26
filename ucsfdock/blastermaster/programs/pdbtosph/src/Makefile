# makefile for DOCK 3.5 chemgrid directory
#
MAKE = make

#F77 = pgf77
F77 = gfortran
FFLAGS = -byteswapio -tp p6 -fast -O2 -I.. -Bstatic 
#FFLAGS = -fconvert=swap -fall-intrinsics -static
FFLAGS = -fconvert=swap -static

#FFLAGS = -O -u
PDBTOSPH = ../bin/pdbtosph
all: pdbtosph

clean:
	/bin/rm -f *.o $(PDBTOSPH)

pdbtosph: 
	$(F77) $(FFLAGS) pdbtosph.f -o $(PDBTOSPH)

