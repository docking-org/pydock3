c----------------------------------------------------------------------
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
      subroutine dconst(unitno, grdcut, grddiv, grdpts, esfact, offset)
c
c  --called from CHEMGRID
c  --increments vdw and electrostatics values at grid points, using
c    a constant dielectric function                 ECMeng    4/91
c----------------------------------------------------------------------
c
      include 'chemgrid.h'
c
      real mincon, minsq
      parameter (mincon=0.0001)
      integer unitno, i, j, k, n
      real r2, r6
c
      minsq=mincon*mincon
c
c  --open parameterized receptor file (from subroutine parmrec)
c
      open (unit=unitno, file='PDBPARM', status='old')
c
  100 read (unitno, 1006, end=500) natm, vdwn, rsra, rsrb, rcrg,
     &(rcrd(i), i=1,3)
 1006 format (2I5, 2(1x, F8.2), 1x, F8.3, 1x, 3F8.3)
      if (vdwn .le. 0) go to 100
c
c  --subtract offset from receptor atom coordinates, find the 3D indices
c    of the nearest grid point (adding 1 because the lowest indices
c    are (1,1,1) rather than (0,0,0)); ignore receptor atoms farther
c    from the grid than the cutoff distance
c
      do 110 i=1,3
        rcrd(i)=rcrd(i) - offset(i)
        nearpt(i)=nint(rcrd(i)/grddiv) + 1
        if (nearpt(i) .gt. (grdpts(i) + grdcut)) go to 100 
        if (nearpt(i) .lt. (1 - grdcut)) go to 100 
  110 continue
c
c --loop through grid points within the cutoff cube (not sphere) of
c   the current receptor atom, but only increment values if the grid
c   point is within the cutoff sphere for the atom
c
      do 400 i=max(1,(nearpt(1)-grdcut)),
     &min(grdpts(1),(nearpt(1)+grdcut))
        gcrd(1)=float(i-1)*grddiv
        do 300 j=max(1,(nearpt(2)-grdcut)),
     &  min(grdpts(2),(nearpt(2)+grdcut))
          gcrd(2)=float(j-1)*grddiv
          do 200 k=max(1,(nearpt(3)-grdcut)),
     &    min(grdpts(3),(nearpt(3)+grdcut))
            gcrd(3)=float(k-1)*grddiv
            n = indx1(i,j,k,grdpts)
            r2 = dist2(rcrd,gcrd)
            if (r2 .gt. cutsq) go to 120
            if (r2 .lt. minsq) then
              bump(n)='X'
              r2 = minsq
            else if(((r2 .lt. cconsq .and. vdwn .le. 5) .or. (r2 .lt.
     &      pconsq .and. vdwn .ge. 8)) .and. bump(n) .eq. 'F') then
              bump(n)='T'
            endif
            r6 = r2*r2*r2
            aval(n)=aval(n) + rsra/(r6*r6)
            bval(n)=bval(n) + rsrb/r6
            esval(n)=esval(n) + 332.0*rcrg/(esfact*sqrt(r2))
  120       continue
  200     continue
  300   continue
  400 continue
      go to 100
  500 continue
      close (unitno)
      return
      end
c----------------------------------------------------------------------
