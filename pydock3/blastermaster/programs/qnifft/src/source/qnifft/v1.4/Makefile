#
# sgi
#
# debug
#FFLAGS = -g -C -col120
#FFLAGS = -s -O2 -col120
#fc = f77
#
# lahey lf95 flags for debugging
#FFLAGS = --verbose --sav --trace --trap diou --warn --wide --chk
# lahey lf95 flags for optimal code
#FFLAGS = --verbose --sav --trace --trap diou --warn --wide -O
#fc = lf95
#
# gnu f77 flags for debugging
#
#FFLAGS = -g -O -Wuninitialized -ffixed-line-length-132
FFLAGS = -fno-ugly -O2 -Wall -ffixed-line-length-132
fc = f77
#
# compaq debug, opt
#FFLAGS = -g -old_f77 -check_bounds -extend_source -check overflow -check underflow 
#FFLAGS = -old_f77 -extend_source  -O -fast 
#fc = f77
#
LFLAGS = 
#LIBS =  /home/sharp/lib/numrec.a /home/sharp/lib/utility.a
DEST = $(DELDIR)/bin
#
#-----------------------
# rules- iris:
.f.o:
	$(fc) -c $(FFLAGS) $*.f
#-----------------------
#
OBJ = qnifft14.o elb.o up.o charge.o radius.o getpdb.o \
 phintp.o scaler.o chrgit.o setbc.o expand6.o wrteps.o cputime.o  mkmaps.o \
 smtheps.o phierg.o inewt.o getpar.o qhelp.o chrgal.o dtstmp.o
OBJ1 = qnifft14.o elb.o up.o charge.o radius.o getpdb.o \
 phintp.o scaler.o chrgit.o setbc.o expand6.o wrteps.o cputime.o  mkmaps.o \
 smtheps.o phierg2.o inewt2.o getpar.o qhelp.o chrgal.o dtstmp.o
FOBJ = qnifft14.f elb.f up.f charge.f radius.f getpdb.f \
 phintp.f scaler.f chrgit.f setbc.f expand6.f wrteps.f cputime.f  mkmaps.f \
 smtheps.f phierg.f inewt2.f getpar.f qhelp.f chrgal.f dtstmp.f
#
qnifft: $(OBJ1)
	echo "      data lstmod / ' `date`'/" > lstmod.h
	f77 $(FFLAGS) -c dtstmp.f
	f77 $(FFLAGS) -o qnifft $(OBJ1)
	mv qnifft $(DEST)/qnifft14
$(OBJ1): qdiffpar.h
qnifft_old: $(OBJ)
	echo "      data lstmod / ' `date`'/" > lstmod.h
	touch dtstmp.f
	f77 $(FFLAGS) -o qnifft $(OBJ)
	mv qnifft $(DEST)/qnifft14
#$(FOBJ): qdiffpar.h
#
