      function lkorad(nam, options0, spheres0)

      use optionstype
      use spheres

      implicit none

      character (len=30) nam
      type(options), intent(in) :: options0
      type(spherest), intent(inout) :: spheres0

      integer n
      integer lkorad, lkpcol !functions      

      n = lkpcol(nam, options0, spheres0)
      if (n .gt. 0) then
        lkorad = n
        return
      endif
      call addcol(nam, spheres0)
      n = lkpcol(nam, options0, spheres0)
      lkorad = n
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

