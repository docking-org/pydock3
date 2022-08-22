c
c       rotates coords in x1 to give least-
c       squares fit to coords in x0
c       uses hermans and ferro program.
c
c-----------------------------------------------------------------------

      subroutine orient(n, x0, x1, cg0, cg1, rot)

      implicit none

      integer n
c        n:  number of points to be matched.
      real x0(3,*)
c        x0:  stationary coordinates.
      real x1(3,*)
c        x1:  coordinates to be least squares fit to x0.
      real cg0(3), cg1(3)
c        cg0:  center of mass of x0.
c        cg1:  center of mass of x1.
      real rot(3,3) !rotation matrix

      real tol
      parameter (tol = 0.001)
      !jklyu, 2019.03.31, just for debug
      !parameter (tol = 0.0001)
c        tol:  tolerance level check for whether iteration
c              can be stopped.
c
c     variables --
      integer i, j, k
c        i:  do loop index.
c        j:  do loop index.
c        k:  do loop index.
      integer ict
c        ict:  counter for number of iterations.  must do at
c              least 3 iterations.
      integer ix, iy, iz
c        ix:  pointer used in iterative least squares.  may
c             have the value 1, 2, or 3.  the value changes on
c             each iteration.
c        iy:  pointer used in iterative least squares.  may
c             have the value 1, 2, or 3.  the value changes on
c             each iteration.
c        iz:  pointer used in iterative least squares.  may
c             have the value 1, 2, or 3.  the value changes on
c             each iteration.
      integer iflag
c        iflag:  indicator variable for whether iterations can be
c                stopped.  must do at least 3 iterations and iflag
c                must equal 0 to be finished.
      real an
c        an:  real number equal to the inverse of the number
c             of points to be matched.
      real sig, gam, sg, bb, cc
c        sig:  temporary variable used in calculation of rotation
c              matrix.
c        gam:  temporary variable used in calculation of rotation
c              matrix.
c        sg:  temporary variable used in calculation of rotation
c              matrix.
c        bb:  temporary variable used in calculation of rotation
c              matrix.
c        cc:  temporary variable used in calculation of rotation
c              matrix.
c     arrays --
      real aa(3,3)
c        aa:  correlation matrix used in calculation of rotation
c             matrix.

c     center of gravity calcn
      an = 1. / float(n)
      iflag = 0
      ix = 0
      do i = 1, 3
        cg0(i) = 0.0
        cg1(i) = 0.0
      enddo
      do i = 1, n
        do j = 1, 3
          cg0(j) = cg0(j) + x0(j,i)
          cg1(j) = cg1(j) + x1(j,i)
        enddo
      enddo
      do i = 1, 3
        cg0(i) = cg0(i)*an
        cg1(i) = cg1(i)*an
      enddo

      do i = 1, 3
        do j = 1, 3
          aa(i,j) = 0.0
        enddo
      enddo
c     calcn of correlation matrix
      do k = 1, n
        aa(1,1) = aa(1,1) + (x1(1,k)-cg1(1)) * (x0(1,k)-cg0(1))
        aa(2,1) = aa(2,1) + (x1(2,k)-cg1(2)) * (x0(1,k)-cg0(1))
        aa(3,1) = aa(3,1) + (x1(3,k)-cg1(3)) * (x0(1,k)-cg0(1))
        aa(1,2) = aa(1,2) + (x1(1,k)-cg1(1)) * (x0(2,k)-cg0(2))
        aa(2,2) = aa(2,2) + (x1(2,k)-cg1(2)) * (x0(2,k)-cg0(2))
        aa(3,2) = aa(3,2) + (x1(3,k)-cg1(3)) * (x0(2,k)-cg0(2))
        aa(1,3) = aa(1,3) + (x1(1,k)-cg1(1)) * (x0(3,k)-cg0(3))
        aa(2,3) = aa(2,3) + (x1(2,k)-cg1(2)) * (x0(3,k)-cg0(3))
        aa(3,3) = aa(3,3) + (x1(3,k)-cg1(3)) * (x0(3,k)-cg0(3))
      enddo
      do i = 1, 3
        do j = 1, 3
          rot(i, j) = 0.0
        enddo
        rot(i, i) = 1.0
      enddo
c       from here to 70, iterative rotation scheme
      ict = 0
      goto 51
   50 continue
      ix = ix+1
      if (ix .lt. 4) goto 52
      if (iflag .eq. 0) goto 70
   51 continue
      iflag = 0
      ix = 1
   52 continue
      ict = ict+1
      if (ict .gt. 100) goto 70
      !jklyu, 2019.03.31, just for debug
      !if (ict .gt. 10000) goto 70
      iy = ix+1
      if (iy .eq. 4) iy = 1
      iz = 6-ix-iy
      sig = aa(iz,iy) - aa(iy,iz)
      gam = aa(iy,iy) + aa(iz,iz)
      sg = sqrt(sig*sig + gam*gam)
      if (sg .eq. 0.0) goto 50
      sg = 1./sg
      if (abs(sig) .le. tol*abs(gam)) goto 50
      do k = 1, 3
        bb = gam*aa(iy,k) + sig*aa(iz,k)
        cc = gam*aa(iz,k) - sig*aa(iy,k)
        aa(iy,k) = bb*sg
        aa(iz,k) = cc*sg
        bb = gam*rot(iy,k) + sig*rot(iz,k)
        cc = gam*rot(iz,k) - sig*rot(iy,k)
        rot(iy,k) = bb*sg
        rot(iz,k) = cc*sg
      enddo
      iflag = 1
      goto 50
   70 continue

      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
