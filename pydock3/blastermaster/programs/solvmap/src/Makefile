
#FC=gfortran
FC=pgf77
CC=/usr/bin/gcc

# debugging on for now
#FLAGS = -g -fcray-pointer -DLINUX
#FLAGS = -tp piv -g -Mbounds -Mdclchk -DLINUX
#FLAGS = -tp piv -g -Mdclchk -DLINUX
#FLAGS = -g -DLINUX
FLAGS = -tp piv -O3 -fast -DLINUX
#CFLAGS = $(FLAGS)
#CFLAGS = -m32 -g -DLINUX
CFLAGS = -m32 -O3 -DLINUX
LFLAGS = $(FLAGS) -Bstatic

SRCF = sas.f cube.f
SRCC = memalloc.c

OBJF=$(SRCF:.f=.o)
OBJC=$(SRCC:.c=.o)

.f.o:
	$(FC) $(FLAGS) -c  $*.f

.c.o:
	$(CC) $(CFLAGS) -c  $*.c

sevsolv: $(OBJF) $(OBJC) $(SRCF) $(SRCC) solvmap.f solvmap.o
	$(FC) $(LFLAGS) -o $@ $(OBJF) $(OBJC) solvmap.o

$(OBJF): dist.h pointer.h

clean:
	rm -f *.o
	rm -f core*
