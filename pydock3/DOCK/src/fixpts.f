c this subroutine fixes the points that orient will use in defining
c the rotation translation matrix.
c BKS, 07/89
c remove global variables RGC, 09/2011
c----------------------------------------------------------------------

      subroutine fixpts(numnod, spcorr, cl, iasign, xos, xor, maxpts, 
     &            maxatmpos, maxctr)

      implicit none

c max sizes (GFortran needs these first)
      integer maxpts
      integer maxatmpos
      integer maxctr

      integer numnod !how many sphere pairs are being used
      real    spcorr(3, maxpts)
c       spcorr:  receptor sphere center coordinates.
      real    cl(3, maxatmpos)
c       cl:  ligand atom center coordinates.  
      integer iasign(maxctr, 2)
      real xos(3, maxctr)
      real xor(3, maxctr)

      integer count
c     just for debugging, jklyu, 2019.03.28
c      character (len=*), parameter :: PDBFORM =
c     &  '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,11X,A2,A2)'

      !just for debugging, jklyu, 2019.03.31
c      write(6,fmt="(A7)",advance="no") "REMARK "
c      do count = 1, numnod
c        !just for debugging, jklyu, 2019.03.30
c        write(6,fmt="(I3,A1,I3,A1)",advance="no") iasign(count, 2),"-"
c     &       ,iasign(count, 1)," "
c      enddo
      !just for debugging, jklyu, 2019.03.30
      !write(6,*)
      !write(6,fmt="(A11)") "match pairs"
      do count = 1, numnod
        !copy receptor spheres into xor
        xor(1, count) = spcorr(1, iasign(count, 2))
        xor(2, count) = spcorr(2, iasign(count, 2))
        xor(3, count) = spcorr(3, iasign(count, 2))
        !copy ligand spheres into xos
        xos(1, count) = cl(1, iasign(count, 1))
        xos(2, count) = cl(2, iasign(count, 1))
        xos(3, count) = cl(3, iasign(count, 1))
        !just for debugging, jklyu, 2019.03.30
c        write(6,PDBFORM) "ATOM  ",iasign(count, 2),
c     &       "C","","REC", "A", 1,"",
c     &       xor(1, count),
c     &       xor(2, count),
c     &       xor(3, count),
c     &       0.00, 0.00,
c     &       "C","" 
c        write(6,PDBFORM) "ATOM  ",iasign(count, 1),
c     &       "C","","LIG", "A", 1,"",
c     &       xos(1, count),
c     &       xos(2, count),
c     &       xos(3, count),
c     &       0.00, 0.00,
c     &       "C","" 
      enddo
      !just for debugging, jklyu, 2019.03.30
      !write(6,*)
      !write(6,fmt="(A11)") "match pairs"
      !write(6,"(A6)") "ENDMDL"

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
c
