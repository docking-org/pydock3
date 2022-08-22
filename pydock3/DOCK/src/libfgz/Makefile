# Set this to the target binary to be created. 32b or 64b
DEFAULT_SIZE = 64
#DEFAULT_COMPILER = pgf
DEFAULT_COMPILER = gnu


#Set default compiler
ifndef COMPILER
COMPILER = $(DEFAULT_COMPILER)
endif

# To build 32b version, "make SIZE=32"
ifndef SIZE
SIZE = $(DEFAULT_SIZE)
endif

# To build GNU version using gfortran, "make GNU=1"

ifeq ($(COMPILER),gnu)
#GNU SPECIFIC

FC = gfortran
CC = gcc

ARCH64 = -mtune=core2
ARCH32 = 

#-O3 is typically buggy in GNU still??
OPT = -g -O2 

FFLAGS = $(ARCH) $(OPT)
CFLAGS = $(ARCH) $(OPT)
endif

ifeq  ($(COMPILER),ifort)
# Path needs to be to 64bit compiler to compile 64b targets 
IFORTPATH = /opt/intel/Compiler/11.1/064/bin/intel64/
FC = $(IFORTPATH)/ifort
CC = gcc

ARCH64 =
ARCH32 = 

OPT = -g -O3  -inline -ip 
FFLAGS = $(ARCH) $(OPT)
CFLAGS = -O2 

endif

ifeq ($(COMPILER),pgf)

# PORTALND GROUP SPECIFIC

#Path to compiler capable of compiling 64bit targets.  This needs to be
#capable of 64b targets so that either 32b or 64b targets can be compiled.
#Thus, use the compiler in .../linux86-x64/... and not /linux86/
#PGPATH = /raid3/software/pgi/9.0.4/linux86-64/9.0-4/bin
PGPATH = /nfs/soft/pgi/current/linux86-64/12.10/bin
#/nfs/soft/pgi/v12_10/linux86-64/12.10/bin/pgf95
#PGPATH = /mnt/nasapps/production/PGI/19.5/linux86-64-llvm/19.5/bin
FC = $(PGPATH)/pgf95 
CC = $(PGPATH)/pgcc 

#Set these to the lowest common 64b and 32b architectures to optimize for.
#ARCH64 = -tp penryn-64
#ARCH32 = -tp piv
ARCH64 = -tp px-64
ARCH32 = -tp px-32

OPT = -O3 -Minline -Mipa=fast,inline:10,libinline -Munroll=c:8,m:4,n:8 -Mfprelaxed -Mvect=sse,assoc,altcode,short -Mcache_align 
FFLAGS = $(ARCH) $(OPT)
CFLAGS = $(ARCH) $(OPT)

endif

# Intel optimizing
#FC = ifort -O3
#CC = icc -O3

ifeq ($(SIZE),64)
ARCH = $(ARCH64)
else 
ARCH = $(ARCH32)
endif


SOURCES = gzopen.F gzclose.F gzread.F gzwrite.F gzbread.F gzbwrite.F \
			gztell.F gzseek.F gzflush.F ftn_gzio.c

libfgz.a: $(SOURCES)
	$(FC) $(FFLAGS) -c *.F
	$(CC) $(CFLAGS) -c ftn_gzio.c
	ar -r libfgz$(SIZE).a *.o

test: test.f libfgz.a
	$(FC) $(FFLAGS) -o test test.f -L. -lfgz$(SIZE) -lz

clean:
	rm -f *.o *.oo libfgz64.a libfgz32.a test *~
