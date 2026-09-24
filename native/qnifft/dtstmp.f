	subroutine dtstmp()
c A call to dtstmp() from any program, 
c in conjunction with Makefile of form:
c
c********************************************
c*** # FFLAGS = -O1 
c*** OBJ = dtstmp.o sub1.o sub2.o main.o 
c*** 
c*** .f.o:
c***	fc $(FFLAGS) -c $*.f
c*** 	fc $(FLAGS)  -c $*.f
c*** main: $(OBJ) 
c***	echo "	data lstmod / ' `date`'/" > lstmod.h
c***	fc -c dtstmp.f
c*** 	fc $(FFLAGS) -o main $(OBJ)
c*** dtstmp.o: lstmod.h
c********************************************
c
c will print date, time, and time of program compilation
c
      character*8 hour
      character*9 day
	character*30 lstmod
	include 'lstmod.h'
c-------------------------------------------------
c
c print title, description
c
c      call date(day)
c      call time(hour) ! kas 1-feb-01 time not supported under lf95
	write(6,*)'  '
	write(6,*)' last modified: ',lstmod
c      write(6,*)' run on ',day,' at ',hour ! kas 1-feb-01 time not supported under lf95
      write(6,*)' run on ',day ! kas 1-feb-01 time not supported under lf95
	write(6,*)'  '
	return
	end
