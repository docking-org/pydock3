c-------------------------------------------------------------------------
c                        distribution version
c      header file for CHEMGRID                             ECMeng  4/91
c-------------------------------------------------------------------------
      integer maxpts
      parameter (maxpts=20000000)
c  maxpts--maximum number of grid points
      integer npts
c  npts--number of grid points
      real aval(maxpts), bval(maxpts), esval(maxpts)
      character*1 bump(maxpts)
c  aval(), bval(), esval(), bump()--values stored "at" grid points
      real rsra, rsrb, rcrg, rcrd(3)
      integer natm, vdwn
c  rsra, rsrb, rcrg, rcrd(), natm, vdwn--values for current receptor atom
      integer nearpt(3)
c  nearpt()--3D indices of grid point closest to current receptor atom
      real gcrd(3)
c  gcrd()--coordinates in angstroms of current grid point
      real grddiv
c  grddiv--spacing of grid points in angstroms
      real boxdim(3)
c  boxdim()--box dimensions in angstroms (x,y,z)
      real offset(3)
c  offset()--box xmin, ymin, zmin in angstroms
      integer grddim(3)
c  grddim()--box dimensions in grid units (x,y,z)
      integer grdpts(3)
c  grdpts()--number of grid points along box dimensions (x,y,z)
c    NOTE: grdpts(i)=griddim(i) + 1   (lowest indices are (1,1,1))
      integer estype
c  estype--type of electrostatic calculation desired:
c    0 = use constant dielectric function
c    1 = use distance-dependent dielectric function
c    2 = use previously generated (DelPhi) electrostatic potential map
      real esfact
c  esfact--factor to multiply dielectric by when estype = 0 or estype = 1;
c    not read or used when estype = 2
c    examples:  D = 1    estype = 0, esfact = 1
c               D = 4    estype = 0, esfact = 4
c               D = r    estype = 1, esfact = 1
c               D = 4r   estype = 1, esfact = 4
      real cutoff, cutsq, pcon, ccon, pconsq, cconsq
c  cutoff--cutoff distance for energy calculations
c  cutsq--cutoff distance squared
c  pcon (ccon)--distance defining a bump with a polar atom (a carbon)
c    of the receptor
c  pconsq, cconsq--the squares of pcon and ccon, respectively
      integer grdcut
c  grdcut--cutoff, in grid units
      real dist2
c  dist2--function to calculate distance squared
      integer indx1
c  indx1--function to convert the 3-dimensional (virtual) indices of a
c    grid point to the actual index in a 1-dimensional array

c
      common
     &/rmaps/ aval, bval, esval
     &/cmap/ bump
     &/vals/ cutsq, pconsq, cconsq
c-------------------------------------------------------------------------
