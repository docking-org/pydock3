c---------------------------------------------------------------------
c                          PROGRAM CHEMGRID
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
c  --calculates the grids used in DOCK 3.0 for force field type scoring,
c  given a receptor pdb file with hydrogens on the polar atoms, and the
c  appropriate parameter tables.  4 grids are created and stored in 3
c  files:  a bump grid/file, an electrostatic potential grid/file and
c  van der Waals repulsion and attraction grids stored in the same 
c  file.
c  ECM  9/93  --version to read input box from SHOWBOX; no output box
c
c  --informative output files:
c    OUTCHEM    restatement of input parameters; messages pertaining
c               to calculation of the grids
c    OUTPARM    messages pertaining to parameterization of receptor
c               atoms; net charge on the receptor molecule (including
c               any ions or waters in the receptor pdb file); created
c               in subroutine parmrec
c    PDBPARM    shows which parameters have been associated with each
c               atom in the receptor pdb file; created in subroutine
c               parmrec, and used in subroutines dconst (constant
c               dielectric) and ddist (distance-dependent dielectric)
c---------------------------------------------------------------------
      include 'chemgrid.h'
c
      character*80 vdwfil
      real scrd(3), sumcrd(3), com(3)
c  scrd()--coordinates of a sphere center
c  sumcrd()--sum of sphere center coordinates so far
c  com()--grid center coordinates; user input or sphere cluster center
c   of mass
      character*80 recfil, sphfil, boxfil, grdfil
c  recfil--pdb-format receptor file (input file)
c  sphfil--file containing one or more sphere clusters, read when the
c   user wishes to center the grids on a sphere cluster center of mass
c   (input file)
c  boxfil--pdb-format file for displaying the grid boundaries (output
c   file)
c  grdfil--prefix name for grid files (output)
      character*80 table, dumlin
c  table--table containing receptor atom parameters
      character*1 ctrtyp
c  ctrtyp--character flag for grid center option; 'u' or 'U' for user
c   input, otherwise a sphere cluster center of mass
      integer i, j, n
c  variables for reading sphere coordinates:
      integer cntemp, nstemp, clnum, nsph
      logical done
c  done--whether or not sphere cluster center of mass has been calculated
c
      open (unit=1, file='INCHEM', status='old')
      open (unit=2, file='OUTCHEM', status='new')
c
      read (1, 1000) recfil
 1000 format (A80)
      write (2, *) 'receptor pdb file:'
      write (2, 1000) recfil
      read (1, 1000) table
      write (2, *) 'receptor parameters will be read from:'
      write (2, 1000) table
      read (1, 1000) vdwfil
      write (2, *) 'van der Waals parameter file:'
      write (2, 1000) vdwfil
c  
      call parmrec(recfil, table, vdwfil, 2)
c  
      read (1, 1000) boxfil
      write (2, *) 'input box file defining grid location:'
      write (2, 1000) boxfil
c
      open (unit=3, file=boxfil, status='old')
      read (3, 1000) dumlin
      read (3, 1001) (com(i), i=1,3)
 1001 format (25x, 3f8.3)
      read (3, 1002) (boxdim(i), i=1,3)
 1002 format (29x, 3f8.3)
      close (3)
c
      write (2, *) 'box center coordinates [x y z]:'
      write (2, *) (com(i), i=1,3)
      write (2, *) 'box x-dimension = ', boxdim(1)
      write (2, *) 'box y-dimension = ', boxdim(2)
      write (2, *) 'box z-dimension = ', boxdim(3)
c
c  --set offset to xmin, ymin, zmin of box
c
      do 65 i=1,3
        offset(i)=com(i) - boxdim(i)/2.0
   65 continue
c
      read (1, *) grddiv
      write (2, *) 'grid spacing in angstroms'
      write (2, *) grddiv
      npts=1
c
c  --convert box dimensions to grid units, rounding upwards
c  --note that points per side .ne. side length in grid units,
c    because lowest indices are (1,1,1) and not (0,0,0)
c
      do 70 i=1,3
        grddim(i)=int(boxdim(i)/grddiv + 1.0)
        grdpts(i)=grddim(i) + 1
        npts=npts*grdpts(i)
   70 continue
      if (npts .gt. maxpts) then
        write (2, *) 'maximum number of grid points exceeded--'
        write (2, *) 'decrease box size, increase grid spacing, or'
        write (2, *) 'increase parameter maxpts'
        write (2, *) 'program stops'
        stop
      endif
      write (2, *) 'grid points per side [x y z]:'
      write (2, *) (grdpts(i), i=1,3)
      write (2, *) 'total number of grid points = ', npts
      read (1, *) estype
      if (estype .ne. 0) then
        estype=1
        write (2, *) 'a distance-dependent dielectric will be used'
      else
        write (2, *) 'a constant dielectric will be used'
      endif
      read (1, *) esfact
      write (2, '(A31, A15, F6.2)') 'the dielectric function will be',
     &' multiplied by ', esfact
   75 continue
      read (1, *) cutoff
      write (2, *) 'cutoff distance for energy calculations:'
      write (2, *) cutoff
      cutsq=cutoff*cutoff
c
c  --convert cutoff to grid units, rounding up (only add 1 rather
c    than 2, because differences in indices rather than the 
c    absolute indices are required)
c
      grdcut=int(cutoff/grddiv + 1.0)
      read (1, *) pcon, ccon
      write (2, *) 'distances defining bumps with receptor atoms:'
      write (2, '(A21, F5.2)') 'receptor polar atoms ', pcon
      write (2, '(A22, F5.2)') 'receptor carbon atoms ', ccon
      pconsq=pcon*pcon
      cconsq=ccon*ccon
      read (1, 1000) grdfil
      write (2, *) 'output grid prefix name:'
      write (2, 1000) grdfil
      close (1)
c
c  --initialize grid
c
      do 80 n=1, maxpts
        aval(n)=0.0
        bval(n)=0.0
        esval(n)=0.0
        bump(n)='F'
   80 continue
c
      if (estype .eq. 0) then
        call dconst(3, grdcut, grddiv, grdpts, esfact, offset)
      else
        call ddist(3, grdcut, grddiv, grdpts, esfact, offset)
      endif
c
      call grdout(grdfil, 3, npts, grddiv, grdpts, offset)
c
      close (2)
      end
c---------------------------------------------------------------------
