      function lkpcol(nam, options0, spheres0)
    
      use optionstype
      use spheres

      implicit none

      character (len=30) nam
      type(options), intent(in) :: options0
      type(spherest), intent(inout) :: spheres0

      character (len=30) nam1, nam2
      integer n
      integer lkpcol

      if (spheres0%numcol .le. 0) then
        lkpcol = 0
        return
      endif
      do n = 1, spheres0%numcol
        nam1 = nam
        nam2 = spheres0%colnam(n)
c       if not case-sensitive, convert both to lower case
        if (.not. options0%casen) then
          call tolow(nam1)
          call tolow(nam2)
        endif
        if (nam1 .eq. nam2) then
          lkpcol = n
          return
        endif
      enddo
      lkpcol = 0
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
