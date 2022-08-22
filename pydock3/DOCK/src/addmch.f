      subroutine addmch(lgnam, spnam, spheres0)

      use errorformats
      use filenums
      use spheres

      implicit none

      character (len=30) lgnam
      character (len=30) spnam
      type(spherest), intent(inout) :: spheres0

      spheres0%nlrmat = spheres0%nlrmat + 1
      if (spheres0%nlrmat .gt. MAXCOLMATCH) then
        write(OUTDOCK, HALT1) 'too many color matches'
        write(OUTDOCK, HALT0)
        stop
      endif
      spheres0%premat(1, spheres0%nlrmat) = lgnam
      spheres0%premat(2, spheres0%nlrmat) = spnam
      
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
