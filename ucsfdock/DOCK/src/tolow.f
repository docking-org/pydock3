      subroutine tolow(word)

      implicit none
      character*(*) word

      character u
      character (len=26) lower
      character (len=26) upper
      integer i, j, n

      lower = 'abcdefghijklmnopqrstuvwxyz'
      upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      n = len(word)

      do j = 1, n
        u = word(j:j)
        i = index(upper, u)
        if (i .le. 0) then
          cycle
        endif
        word(j:j) = lower(i:i)
      enddo

      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
