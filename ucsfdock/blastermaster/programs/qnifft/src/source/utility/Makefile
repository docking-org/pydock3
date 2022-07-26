# gnu
# debug
#FFLAGS = -g -O -Wuninitialized -ffixed-line-length-132
# optimize
FFLAGS = -fno-ugly -O2 -Wall -ffixed-line-length-132
FCOMP = g77
#
BIN = $(DELDIR)/bin
#-----------------
# rules- general
.f.o:
	$(FCOMP) -c $(FFLAGS) $*.f
#-----------------
ionintm_rad: ionintm_rad.o
	$(FCOMP) $(FFLAGS) -o ionintm_rad ionintm_rad.o
	mv ionintm_rad $(BIN)
ionintm: ionintm.o
	$(FCOMP) $(FFLAGS) -o ionintm ionintm.o
	mv ionintm $(BIN)
sumfld: sumfld.o
	$(FCOMP) $(FFLAGS) -o sumfld sumfld.o
	mv sumfld $(BIN)
chkchg: chkchg.o 
	$(FCOMP) $(FFLAGS) -o chkchg chkchg.o 
	mv chkchg $(BIN)
addmap: addmap.o 
	$(FCOMP) $(FFLAGS) -o addmap addmap.o 
	mv addmap $(BIN)
