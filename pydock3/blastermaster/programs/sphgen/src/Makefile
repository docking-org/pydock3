# makefile for DOCK 3.5 chemgrid directory
#
MAKE = make
BIN = ../

F77 = pgf77
FFLAGS = -byteswapio -tp p6 -fast -O2 -I.. -Bstatic 

#FFLAGS = -O -u
OBJS = aesthetics.o chemgrid.o dist.o ddist.o dconst.o grdout.o parmrec.o 
OBJ2 = getxpdb.o conect2.o
EXES = chemgrid aesthetics getfiltered getmols getsdf linreg getxpdb conect2 splitmol filt.exe

all: chemgrid getsdf splitmol getmols getxpdb filt

install: all
	mv chemgrid big

clean:
	/bin/rm -f *.o $(EXES)

chemgrid: $(OBJS)
	$(F77)  $(FFLAGS) $(OBJS) -o chemgrid

#chemgrid.o: chemgrid.h
#ddist.o: chemgrid.h
#dconst.o:  chemgrid.h
#grdout.o:  chemgrid.h
#parmrec.o: parmrec.h


aesthetics: aesthetics.o
	$(F77) $(FFLAGS) -o aesthetics aesthetics.o

getfiltered: getfiltered.o
	$(F77) $(FFLAGS) -o getfiltered getfiltered.o

getmols: getmols.o aesthetics.o
	$(F77) $(FFLAGS) -o getmols getmols.o aesthetics.o

getsdf: getsdf.o aesthetics.o
	$(F77) $(FFLAGS) -o getsdf getsdf.o aesthetics.o

getxpdb: getxpdb.o aesthetics.o
	$(F77) $(FFLAGS) -o getxpdb getxpdb.o aesthetics.o

linreg: linreg.o
	$(F77) $(FFLAGS) -o linreg linreg.o

rmsd: rmsd.o
	$(F77) $(FFLAGS) -o rmsd rmsd.o

splitat: splitat.o
	$(F77) $(FFLAGS) -o splitat splitat.o

splitmol: splitmol.o aesthetics.o
	$(F77) $(FFLAGS) -o splitmol splitmol.o aesthetics.o

conect2: conect2.o
	$(F77) $(FFLAGS) -o conect2 conect2.o

filt: 
	$(F77) $(FFLAGS) ../filter1.f -o filt.exe

sphgen: 
	$(F77) $(FFLAGS) ../sphgen.f -o sphgen

%.o: ../%.f
	$(F77) -c $(FFLAGS) -o $*.o ../$*.f

