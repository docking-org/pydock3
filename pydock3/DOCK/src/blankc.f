
      subroutine blankc(fline, lline, MAXLNS, lines)

      implicit none

      integer fline, lline
c where INDOCK is stored...
      integer, intent(in) :: MAXLNS
      character (len=300), dimension(MAXLNS), intent(inout) :: lines

      integer i, j, k

      do i = fline, lline
        do j = 1, 132
          if (lines(i)(j:j) .eq. '#') then
            do k = j, 132
              lines(i)(k:k) = ' '
            enddo
            cycle
          endif
        enddo
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
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c
