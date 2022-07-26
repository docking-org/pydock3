#makefile for portland group compiler. Ryan Coleman 2012.
PGPATH = /raid3/software/pgi/9.0.4/linux86-64/9.0-4/bin
#FFLAGS = -O3 -byteswapio -Mextend -tp px-32 -fastsse -Minline -Mipa=fast,inline:10,libinline,libopt,vestigial -Munroll=c:8,m:4,n:8 -Mfprelaxed -Mvect=sse,assoc,altcode,short -Mcache_align
FFLAGS = -O3 -byteswapio -Mextend -tp px-32 -fastsse -Minline -Mipa=fast,inline:10,libinline,libopt,vestigial -Munroll=c:8,m:4,n:8 -Mfprelaxed -Mvect=sse,assoc,altcode,short -Mcache_align
FCOMP = $(PGPATH)/pgf95

#
BIN = $(DELDIR)/bin
#-----------------
# rules- general
.f.o:
	$(FCOMP) -c $(FFLAGS) $*.f
#-----------------

#
FOBJ = qnifft22.f elb.f up.f charge.f radius.f getpdb.f \
 phintp.f scaler.f chrgit.f setbc.f wrteps.f cputime.f  \
 phierg.f inewt.f getpar.f qhelp.f chrgal.f dtstmp.f dmpeps.f phirxn.f phifrc.f \
 sasgen.f mapgen.f girdgen.f utility.f pntchk.f mkepsd.f
#
OBJ1 = qcall.o qnifft22_sub.o qcen.o elb.o up.o charge.o radius.o getpdb.o \
 phintp.o scaler.o chrgit.o setbc.o wrteps.o cputime.o  \
 phierg.o inewt.o getpar.o qhelp.o chrgal.o dtstmp.o dmpeps.o phirxn.o phifrc.o \
 mapit.o sasgen.o mapgen.o girdgen.o utility.o pntchk.o mkepsd.o
#
qnifft_huge: 
	/bin/cp qdiffpar193.h qdiffpar.h
	echo "      data lstmod / ' `date`'/" > lstmod.h
	$(FCOMP) $(FFLAGS) -c dtstmp.f
	$(FCOMP) $(FFLAGS) -Bstatic $(FOBJ) 
	mv a.out qnifft22_193_pgf_32
#
#
qnifft_small: 
	/bin/cp qdiffpar65.h qdiffpar.h
	echo "      data lstmod / ' `date`'/" > lstmod.h
	$(FCOMP) $(FFLAGS) -c dtstmp.f
	$(FCOMP) $(FFLAGS) $(FOBJ) 
	mv a.out qnifft22_65
#
qnifft_medium: 
	/bin/cp qdiffpar129.h qdiffpar.h
	echo "      data lstmod / ' `date`'/" > lstmod.h
	$(FCOMP) $(FFLAGS) -c dtstmp.f
	$(FCOMP) $(FFLAGS) $(FOBJ) 
	mv a.out qnifft22_129
#
qnifft_large: 
	/bin/cp qdiffpar161.h qdiffpar.h
	echo "      data lstmod / ' `date`'/" > lstmod.h
	$(FCOMP) $(FFLAGS) -c dtstmp.f
	$(FCOMP) $(FFLAGS) $(FOBJ) 
	mv a.out qnifft22_161
clean:
	/bin/rm *.o *.oo qnifft22_193_pgf_32

qnifft_sub: $(OBJ1)
	echo "      data lstmod / ' `date`'/" > lstmod.h
	$(FCOMP) $(FFLAGS) -c dtstmp.f
	$(FCOMP) $(FFLAGS) -o qcall $(OBJ1) 
	mv qcall qcall22a
$(OBJ1): qdiffpar.h

install:  qnifft_huge
	cp qnifft22_193_pgf_32 ../../../../bin
