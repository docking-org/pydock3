      integer function nxtsp(line,j)

      implicit none
      character*(*) line
      integer i,j,l

      l = len(line)
      do i = j, l
        nxtsp = i
        if (line(i:i) .eq. ' ') exit
      enddo
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
