      subroutine addcol(nam, spheres0)

      use errorformats
      use filenums
      use spheres

      implicit none

      !character*30 nam ! OLD
      character (len=30) nam
      type(spherest), intent(inout) :: spheres0

      if (spheres0%numcol .ge. MAXCOL) then
        write(OUTDOCK, HALT3) 'Too many colors, maximum is: ', MAXCOL
        write(OUTDOCK, HALT0)
        stop
      endif
      spheres0%numcol = spheres0%numcol + 1
      spheres0%colnam(spheres0%numcol) = nam

      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c

