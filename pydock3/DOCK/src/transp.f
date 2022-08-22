      subroutine transp(nr, rcols, spheres0)

      use errorformats
      use filenums
      use spheres

      implicit none

      integer nr
      integer rcols(*)
      type(spherest), intent(inout) :: spheres0

      integer r, newnum
      integer nonzer

      nonzer = 0
      do r = 1, nr
      if (rcols(r) .le. 0) cycle
      newnum = spheres0%rectrn(rcols(r))
      if (newnum .le. 0) then
        write(OUTDOCK, HALT3)
     &       'tranlg: invalid ligand color number: ', rcols(r)
        write(OUTDOCK, HALT0)
        stop
      endif
      rcols(r) = newnum
      nonzer = nonzer + 1
      enddo
      write(OUTDOCK, 125) nonzer
 125  format(i5,' non-zero colors in receptor sphere cluster')
      return
      end
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c

